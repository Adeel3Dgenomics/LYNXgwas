import java.io.*;
import java.nio.file.*;
import java.util.*;
import java.util.zip.GZIPOutputStream;

/**
 * Regression test for RegulatoryEnrichmentAnalyzer, using a small hand-built fixture (a single
 * peak per mark, a handful of foreground/background SNP positions) chosen so the 2x2 contingency
 * table, odds ratio, and Fisher's exact p-value can all be worked out independently by hand.
 *
 * Fixture: one peak on chr1 [1000,2000] for each of the 3 marks (H3K27ac/H3K4me1/H3K4me3), for a
 * synthetic EID "TESTE".
 *
 * Foreground SNPs (5): 1200, 1500, 1900 (inside the peak -> 3 "in"), 500, 2500 (outside -> 2 "not").
 * Background SNPs (10): 1100, 1800 (inside -> 2 "in"), 100,200,300,400,3000,3100,3200,3300 (8 "not").
 *
 * By-hand 2x2 table (foreground-in-peak, foreground-not, background-in-peak, background-not):
 *   a=3, b=2, c=2, d=8
 * By-hand odds ratio: (a*d)/(b*c) = (3*8)/(2*2) = 24/4 = 6.0 (no zero cells here, plain formula).
 * The p-value is NOT re-derived by hand (Fisher's exact hypergeometric tail-sum is already
 * independently verified against the classic tea-tasting example in EnrichmentAnalyzerTest) — this
 * test instead cross-checks that RegulatoryEnrichmentAnalyzer.run(...) produces the exact same
 * p-value as calling EnrichmentAnalyzer.fisherExactTwoSided(3,2,2,8) directly on the same margins,
 * which really confirms the table (a,b,c,d) was built correctly from the interval-overlap
 * classification, not that Fisher's exact itself is correct (already covered elsewhere).
 */
public class RegulatoryEnrichmentAnalyzerTest {

    public static void main(String[] args) throws Exception {
        int failures = 0;
        Path root = Files.createTempDirectory("regenrich_test");
        try {
            String eid = "TESTE";
            for (String mark : RegulatoryPeakIndex.MARKS) {
                writePeakFile(root, eid, mark, "chr1\t1000\t2000\tPeak1\t500\t.\t10.0\t20.0\t18.0\t50");
            }
            RegulatoryPeakIndex idx = new RegulatoryPeakIndex(root.toString());

            List<RegulatoryEnrichmentAnalyzer.SnpPos> fg = Arrays.asList(
                new RegulatoryEnrichmentAnalyzer.SnpPos("chr1", 1200), // in peak
                new RegulatoryEnrichmentAnalyzer.SnpPos("chr1", 1500), // in peak
                new RegulatoryEnrichmentAnalyzer.SnpPos("chr1", 1900), // in peak
                new RegulatoryEnrichmentAnalyzer.SnpPos("chr1", 500),  // not in peak
                new RegulatoryEnrichmentAnalyzer.SnpPos("chr1", 2500)  // not in peak
            );
            List<RegulatoryEnrichmentAnalyzer.SnpPos> bg = Arrays.asList(
                new RegulatoryEnrichmentAnalyzer.SnpPos("chr1", 1100), // in peak
                new RegulatoryEnrichmentAnalyzer.SnpPos("chr1", 1800), // in peak
                new RegulatoryEnrichmentAnalyzer.SnpPos("chr1", 100),
                new RegulatoryEnrichmentAnalyzer.SnpPos("chr1", 200),
                new RegulatoryEnrichmentAnalyzer.SnpPos("chr1", 300),
                new RegulatoryEnrichmentAnalyzer.SnpPos("chr1", 400),
                new RegulatoryEnrichmentAnalyzer.SnpPos("chr1", 3000),
                new RegulatoryEnrichmentAnalyzer.SnpPos("chr1", 3100),
                new RegulatoryEnrichmentAnalyzer.SnpPos("chr1", 3200),
                new RegulatoryEnrichmentAnalyzer.SnpPos("chr1", 3300)
            );

            RegulatoryEnrichmentAnalyzer.Result result =
                RegulatoryEnrichmentAnalyzer.run(eid, "Test Tissue", 5e-8, fg, bg, idx);

            failures += check("result.eid preserved", "TESTE".equals(result.eid));
            failures += check("result has 3 marks", result.marks.size() == 3);

            double expectedP = EnrichmentAnalyzer.fisherExactTwoSided(3, 2, 2, 8);
            for (RegulatoryEnrichmentAnalyzer.MarkResult mr : result.marks) {
                String tag = mr.mark + ": ";
                failures += check(tag + "n_foreground_in_peak == 3", mr.nForegroundInPeak == 3);
                failures += check(tag + "n_foreground_total == 5", mr.nForegroundTotal == 5);
                failures += check(tag + "n_background_in_peak == 2", mr.nBackgroundInPeak == 2);
                failures += check(tag + "n_background_total == 10", mr.nBackgroundTotal == 10);
                failures += checkClose(tag + "odds ratio == 6.0 (hand-computed 3*8/(2*2))", mr.oddsRatio, 6.0, 1e-9);
                failures += checkClose(tag + "p matches EnrichmentAnalyzer.fisherExactTwoSided(3,2,2,8) directly",
                    mr.pValue, expectedP, 1e-12);
            }

            // ── Unmapped disease (eid == null) must return an empty mark list, not an error ──
            RegulatoryEnrichmentAnalyzer.Result unmapped =
                RegulatoryEnrichmentAnalyzer.run(null, null, 5e-8, fg, bg, idx);
            failures += check("unmapped eid returns null eid in result", unmapped.eid == null);
            failures += check("unmapped eid returns empty marks list (not an error)", unmapped.marks.isEmpty());

            // ── Odds ratio zero-cell handling: Haldane-Anscombe correction (+0.5 to every cell) ──
            // a=0 -> by hand: ((0+0.5)*(8+0.5)) / ((5+0.5)*(2+0.5)) = (0.5*8.5)/(5.5*2.5) = 4.25/13.75
            double expectedCorrected = (0.5 * 8.5) / (5.5 * 2.5);
            double corrected = RegulatoryEnrichmentAnalyzer.oddsRatio(0, 5, 2, 8);
            failures += checkClose("odds ratio zero-cell Haldane-Anscombe correction", corrected, expectedCorrected, 1e-9);
            failures += check("odds ratio zero-cell correction is finite (not Infinity/NaN)",
                Double.isFinite(corrected));

            // ── Plain (no zero cells) odds ratio uses the direct formula ──
            failures += checkClose("odds ratio no zero cells uses plain a*d/(b*c)",
                RegulatoryEnrichmentAnalyzer.oddsRatio(3, 2, 2, 8), 6.0, 1e-12);

            if (failures == 0) {
                System.out.println("PASS: all RegulatoryEnrichmentAnalyzer tests passed");
            } else {
                System.out.println("FAIL: " + failures + " test(s) failed");
                System.exit(1);
            }
        } finally {
            deleteRecursive(root.toFile());
        }
    }

    private static void writePeakFile(Path root, String eid, String mark, String line) throws IOException {
        Path dir = root.resolve(eid);
        Files.createDirectories(dir);
        File gz = dir.resolve(eid + "-" + mark + ".narrowPeak.gz").toFile();
        try (GZIPOutputStream out = new GZIPOutputStream(new FileOutputStream(gz))) {
            out.write((line + "\n").getBytes("UTF-8"));
        }
    }

    private static void deleteRecursive(File f) {
        File[] children = f.listFiles();
        if (children != null) for (File c : children) deleteRecursive(c);
        f.delete();
    }

    private static int check(String label, boolean cond) {
        if (!cond) {
            System.out.println("FAIL: " + label);
            return 1;
        }
        System.out.println("PASS: " + label);
        return 0;
    }

    private static int checkClose(String label, double actual, double expected, double tol) {
        if (Double.isNaN(actual) || Math.abs(actual - expected) > Math.abs(tol) + 1e-15) {
            System.out.println("FAIL: " + label + " — expected " + expected + ", got " + actual);
            return 1;
        }
        System.out.println("PASS: " + label + " (" + actual + ")");
        return 0;
    }
}
