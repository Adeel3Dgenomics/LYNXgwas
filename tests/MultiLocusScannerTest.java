import java.io.*;
import java.nio.file.*;
import java.util.*;
import loci.LociIdentifier;

/**
 * Standalone regression test (no external test framework — this project has none) for the
 * MultiLocusScanner bug found during the six-disease validation: a dataset whose raw GWAS file
 * is not sorted by position, or whose chromosome column uses a different label convention
 * ("chr1" vs "1") than the loci list it's being scanned against, previously showed near-zero
 * coverage across the whole Locus Matrix comparison. Both conditions are exercised here together.
 *
 * Run with: java -cp bin;lib/* MultiLocusScannerTest   (after compiling this file into bin/ too)
 * Exits 0 and prints PASS on success, exits 1 and prints FAIL with details otherwise.
 */
public class MultiLocusScannerTest {

    public static void main(String[] args) throws Exception {
        int failures = 0;
        failures += testUnsortedFileAndMixedChrLabels();
        failures += testNormalizeChr();

        if (failures == 0) {
            System.out.println("PASS: all MultiLocusScanner tests passed");
        } else {
            System.out.println("FAIL: " + failures + " test(s) failed");
            System.exit(1);
        }
    }

    private static int testUnsortedFileAndMixedChrLabels() throws Exception {
        // Two loci on the SAME chromosome, one labeled "1" (as if discovered from a dataset using
        // bare chromosome numbers) matching what a differently-labeled dataset's own file ("chr1")
        // must still be matched against.
        LociIdentifier.IdentifiedLocus locusA = new LociIdentifier.IdentifiedLocus();
        locusA.chr = "1"; locusA.start = 1_000_000; locusA.end = 1_200_000;
        LociIdentifier.IdentifiedLocus locusB = new LociIdentifier.IdentifiedLocus();
        locusB.chr = "1"; locusB.start = 5_000_000; locusB.end = 5_200_000;
        List<LociIdentifier.IdentifiedLocus> loci = Arrays.asList(locusA, locusB);

        // Deliberately: (1) uses "chr1" instead of "1", (2) rows are NOT sorted by position —
        // locus B's best SNP appears in the file before locus A's, and locus A's rows are also
        // internally out of order.
        Path tmp = Files.createTempFile("mls-test", ".tsv");
        try (PrintWriter pw = new PrintWriter(new FileWriter(tmp.toFile()))) {
            pw.println("chrom\tpos\tp\tea\tnea");
            pw.println("chr1\t5100000\t1e-9\tA\tG");   // locus B best SNP, appears first
            pw.println("chr1\t1100000\t0.5\tA\tG");    // locus A weak SNP, out of position order
            pw.println("chr1\t1050000\t1e-6\tA\tG");   // locus A best SNP, appears AFTER a later position
            pw.println("chr1\t9000000\t0.9\tA\tG");    // outside both loci, must be ignored
            pw.println("chr1\t5150000\t0.2\tA\tG");    // locus B weak SNP
        }

        Config cfg = new Config();
        cfg.gwasFile = tmp.toString();
        cfg.colChr = "chrom"; cfg.colPos = "pos"; cfg.colPvalue = "p";
        cfg.colEa = "ea"; cfg.colNea = "nea";
        cfg.outputDir = "output";

        MultiLocusProgress progress = new MultiLocusProgress();
        MultiLocusResult result = MultiLocusScanner.scan(
            Collections.singletonList(cfg), Collections.singletonList("test-ds"), loci, null, progress);
        Files.deleteIfExists(tmp);

        int failures = 0;
        if (result.loci.size() != 2) {
            System.out.println("FAIL: expected 2 loci rows, got " + result.loci.size());
            return 1;
        }
        // Loci are re-sorted by chr/start internally, so index 0 = locus A (start 1,000,000).
        MultiLocusResult.LocusRow rowA = result.loci.get(0);
        MultiLocusResult.LocusRow rowB = result.loci.get(1);
        MultiLocusResult.DatasetLocusStat statA = rowA.cells.get("test-ds");
        MultiLocusResult.DatasetLocusStat statB = rowB.cells.get("test-ds");

        if (statA == null || Double.isNaN(statA.bestP) || statA.bestP != 1e-6) {
            System.out.println("FAIL: locus A best p expected 1e-6, got " +
                (statA == null ? "null" : statA.bestP));
            failures++;
        }
        if (statB == null || Double.isNaN(statB.bestP) || statB.bestP != 1e-9) {
            System.out.println("FAIL: locus B best p expected 1e-9, got " +
                (statB == null ? "null" : statB.bestP));
            failures++;
        }
        if (failures == 0) {
            System.out.println("PASS: unsorted file + mixed chr-label ('chr1' vs '1') both matched correctly");
        }
        return failures;
    }

    private static int testNormalizeChr() {
        int failures = 0;
        String[][] cases = {
            {"1", "1"}, {"chr1", "1"}, {"Chr1", "1"}, {"CHR1", "1"},
            {"X", "X"}, {"chrX", "X"}, {"MT", "MT"}, {"chrMT", "MT"}, {"22", "22"}, {"chr22", "22"},
        };
        for (String[] c : cases) {
            String got = MultiLocusScanner.normalizeChr(c[0]);
            if (!got.equals(c[1])) {
                System.out.println("FAIL: normalizeChr(\"" + c[0] + "\") expected \"" + c[1] + "\", got \"" + got + "\"");
                failures++;
            }
        }
        if (failures == 0) System.out.println("PASS: normalizeChr handles chr-prefix/case variants");
        return failures;
    }
}
