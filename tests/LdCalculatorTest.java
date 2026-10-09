/**
 * Standalone regression test (no external test framework — this project has none) for the
 * LD-window memory-safety cap added after the MHC-region JVM native-memory crash discovered
 * during the six-disease validation (see manuscript Section 3 / DECISIONS.md).
 */
public class LdCalculatorTest {

    public static void main(String[] args) {
        int failures = 0;
        failures += check("typical 250kb-padded locus stays uncapped",
            LdCalculator.computeWindowKb(1_000_000, 1_450_000), 452);
        failures += check("1.5Mb-span locus stays uncapped (below cap)",
            LdCalculator.computeWindowKb(0, 1_500_000), 1502);
        failures += check("8.8Mb MHC-sized locus is capped at MAX_LD_WINDOW_KB",
            LdCalculator.computeWindowKb(25_500_000, 34_300_000), LdCalculator.MAX_LD_WINDOW_KB);
        failures += check("cap value itself is a sane, non-degenerate size",
            LdCalculator.MAX_LD_WINDOW_KB > 0 && LdCalculator.MAX_LD_WINDOW_KB <= 5000 ? 1 : 0, 1);

        // PLINK --r2 square output whose first SNP is monomorphic: its row starts with "nan". That row is data,
        // not a header; reading it as a header shifted every row of the LD triangle up by one.
        try {
            java.io.File f = java.io.File.createTempFile("ldsq", ".ld");
            try (java.io.PrintWriter pw = new java.io.PrintWriter(f)) {
                pw.println("nan	nan	nan");
                pw.println("nan	1	0.25");
                pw.println("nan	0.25	1");
            }
            double[][] m = LdCalculator.parsePairwiseMatrix(f, 3);
            f.delete();
            failures += check("monomorphic first SNP: its pairs stay missing", Double.isNaN(m[0][1]) && Double.isNaN(m[0][2]) ? 1 : 0, 1);
            failures += check("monomorphic first SNP: rows are not shifted (r2 of SNPs 2 and 3)", Math.round(m[1][2] * 100), 25);
            failures += check("monomorphic first SNP: matrix stays symmetric", m[2][1] == m[1][2] ? 1 : 0, 1);
        } catch (java.io.IOException e) {
            failures += check("square LD file could be written: " + e.getMessage(), 0, 1);
        }

        if (failures == 0) {
            System.out.println("PASS: all LdCalculator tests passed");
        } else {
            System.out.println("FAIL: " + failures + " test(s) failed");
            System.exit(1);
        }
    }

    private static int check(String label, long actual, long expected) {
        if (actual != expected) {
            System.out.println("FAIL: " + label + " — expected " + expected + ", got " + actual);
            return 1;
        }
        System.out.println("PASS: " + label);
        return 0;
    }
}
