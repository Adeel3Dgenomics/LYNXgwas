import java.io.*;
import java.nio.file.*;
import java.util.*;

/**
 * Standalone regression test (no external test framework — this project has none) for the
 * distance-based pruning added to GwasQc's lambda_GC computation (previously computed over every
 * SNP with no independence filter at all).
 */
public class GwasQcTest {

    public static void main(String[] args) throws Exception {
        int failures = 0;
        failures += testPruningReducesDenseBlock();
        failures += testMissingChrPosFallsBackGracefully();

        if (failures == 0) {
            System.out.println("PASS: all GwasQc tests passed");
        } else {
            System.out.println("FAIL: " + failures + " test(s) failed");
            System.exit(1);
        }
    }

    private static int testPruningReducesDenseBlock() throws Exception {
        // 100 SNPs packed into a single 10kb block (one dense LD block) on chr1, all with the same
        // strong chi2, plus 3 well-separated SNPs elsewhere at neutral chi2. Unpruned, the dense
        // block dominates the median. Pruned (one per 250kb), the block should collapse to 1 SNP.
        Path tmp = Files.createTempFile("qc-test", ".tsv");
        try (PrintWriter pw = new PrintWriter(new FileWriter(tmp.toFile()))) {
            pw.println("chrom\tpos\tp");
            for (int i = 0; i < 100; i++) {
                pw.println("1\t" + (1_000_000 + i * 100) + "\t1e-8"); // dense block, strong signal
            }
            pw.println("1\t10000000\t0.5");
            pw.println("2\t5000000\t0.5");
            pw.println("3\t5000000\t0.5");
        }

        Config cfg = new Config();
        cfg.gwasFile = tmp.toString();
        cfg.colChr = "chrom"; cfg.colPos = "pos"; cfg.colPvalue = "p";

        GwasQc.Result r = GwasQc.compute(cfg);
        Files.deleteIfExists(tmp);

        int failures = 0;
        if (!r.ok) { System.out.println("FAIL: compute() returned not-ok: " + r.error); return 1; }
        if (r.nSnpsScanned != 103) {
            System.out.println("FAIL: expected 103 scanned SNPs, got " + r.nSnpsScanned);
            failures++;
        }
        // Pruned: 1 from the dense block + 1 more on chr1 (well-separated) + 1 on chr2 + 1 on chr3 = 4.
        if (r.nSnps != 4) {
            System.out.println("FAIL: expected pruning to collapse the 100-SNP dense block to 1 " +
                "(4 total kept SNPs), got " + r.nSnps);
            failures++;
        }
        if (failures == 0) {
            System.out.println("PASS: distance-based pruning collapses a dense 100-SNP LD block " +
                "to 1 representative SNP (103 scanned -> 4 kept)");
        }
        return failures;
    }

    private static int testMissingChrPosFallsBackGracefully() throws Exception {
        // If chr/pos columns aren't mapped, pruning must not crash — it should fall back to using
        // every SNP, matching the pre-fix behavior rather than failing the whole QC check.
        Path tmp = Files.createTempFile("qc-test-nochrpos", ".tsv");
        try (PrintWriter pw = new PrintWriter(new FileWriter(tmp.toFile()))) {
            pw.println("p");
            for (int i = 0; i < 10; i++) pw.println("0.1");
        }

        Config cfg = new Config();
        cfg.gwasFile = tmp.toString();
        cfg.colChr = "chrom_not_present"; cfg.colPos = "pos_not_present"; cfg.colPvalue = "p";

        GwasQc.Result r = GwasQc.compute(cfg);
        Files.deleteIfExists(tmp);

        if (!r.ok || r.nSnps != 10 || r.nSnpsScanned != 10) {
            System.out.println("FAIL: expected graceful no-pruning fallback (10 scanned, 10 kept), got ok=" +
                r.ok + " nSnps=" + r.nSnps + " nSnpsScanned=" + r.nSnpsScanned);
            return 1;
        }
        System.out.println("PASS: missing chr/pos columns fall back to unpruned computation without crashing");
        return 0;
    }
}
