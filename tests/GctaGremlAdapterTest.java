import java.io.*;
import java.nio.file.*;
import java.util.*;

/**
 * Standalone regression test (no external test framework — this project has none) for
 * GctaGremlAdapter: .hsq parsing against a hand-computed fixture, the binary-not-found path
 * (via the shared GctaBinaryResolver also used by CojoAdapter), and the honest
 * refusal-rather-than-fabrication behavior when no phenotype file is available — the
 * correctness property this adapter is required to uphold per DECISIONS_PHASE2.md's stance on
 * never silently degrading to a mocked result.
 */
public class GctaGremlAdapterTest {

    public static void main(String[] args) throws Exception {
        int failures = 0;
        failures += testParseHsq();
        failures += testGctaBinaryNotFound();
        failures += testRemlRefusesWithoutPhenotypeFile();
        failures += testRemlRefusesWhenPhenotypeFileMissing();

        if (failures == 0) {
            System.out.println("PASS: all GctaGremlAdapter tests passed");
        } else {
            System.out.println("FAIL: " + failures + " test(s) failed");
            System.exit(1);
        }
    }

    /**
     * .hsq fixture (GCTA's real output format: tab-separated Source/Variance/SE, with SE absent
     * on some rows). Values below are internally consistent by construction, not independently
     * re-derived statistics (that's GCTA's job, not this parser's) — this test verifies field
     * extraction is faithful to the file, including the two rows (logL, n) that have no third
     * (SE) column at all, which must come back as SE=NaN rather than crashing or misreading the
     * next row's Source as an SE value.
     *   V(G)/Vp   = 0.350000  (this IS the row parseHsq/reml() must surface as h²)
     *   V(G)/Vp SE = 0.045000
     *   n (last row, Variance-only) = 5000, SE absent -> NaN
     */
    private static int testParseHsq() throws Exception {
        Path tmp = Files.createTempFile("gcta-reml", ".hsq");
        try (PrintWriter pw = new PrintWriter(new FileWriter(tmp.toFile()))) {
            pw.println("Source\tVariance\tSE");
            pw.println("V(G)\t0.350000\t0.050000");
            pw.println("V(e)\t0.650000\t0.048000");
            pw.println("Vp\t1.000000\t0.030000");
            pw.println("V(G)/Vp\t0.350000\t0.045000");
            pw.println("logL\t-1234.560000");
            pw.println("LRT\t12.340000");
            pw.println("df\t1.000000");
            pw.println("Pval\t0.000220");
            pw.println("n\t5000.000000");
        }

        List<GctaGremlAdapter.HsqResult> rows = GctaGremlAdapter.parseHsq(tmp.toFile());
        Files.deleteIfExists(tmp);

        int failures = 0;
        failures += checkEq("row count (9 data rows)", rows.size(), 9);

        GctaGremlAdapter.HsqResult vg = findSource(rows, "V(G)");
        failures += checkNotNull("V(G) row present", vg);
        if (vg != null) {
            failures += checkClose("V(G) variance", vg.variance, 0.35, 1e-9);
            failures += checkClose("V(G) SE", vg.se, 0.05, 1e-9);
        }

        GctaGremlAdapter.HsqResult h2Row = findSource(rows, "V(G)/Vp");
        failures += checkNotNull("V(G)/Vp row present (this is h²)", h2Row);
        if (h2Row != null) {
            failures += checkClose("V(G)/Vp (h²) variance", h2Row.variance, 0.35, 1e-9);
            failures += checkClose("V(G)/Vp (h²) SE", h2Row.se, 0.045, 1e-9);
        }

        GctaGremlAdapter.HsqResult nRow = findSource(rows, "n");
        failures += checkNotNull("n row present", nRow);
        if (nRow != null) {
            failures += checkClose("n variance-column value", nRow.variance, 5000.0, 1e-9);
            failures += check("n row has NO SE (no third column in the fixture) -> NaN, not misread from next line",
                Double.isNaN(nRow.se));
        }

        return failures;
    }

    /**
     * makeGrm()/reml() both resolve GCTA via the shared GctaBinaryResolver (also used by
     * CojoAdapter after its refactor) — confirms the exact CojoAdapter-style message text is
     * preserved for a nonexistent path.
     */
    private static int testGctaBinaryNotFound() throws Exception {
        Path dir = Files.createTempDirectory("gcta-notfound-test");
        File matchedDir = new File(dir.toFile(), "matched");
        File runDir = new File(dir.toFile(), "run");
        matchedDir.mkdirs();
        // matched_ref.bim deliberately absent — but the binary check must fire first regardless.
        String bogusGcta = new File(dir.toFile(), "no_such_gcta64.exe").getAbsolutePath();

        boolean threw = false;
        String message = null;
        try {
            GctaGremlAdapter.makeGrm(matchedDir, runDir, bogusGcta);
        } catch (IOException e) {
            threw = true;
            message = e.getMessage();
        }

        int failures = 0;
        failures += check("makeGrm() throws IOException for a nonexistent GCTA binary", threw);
        if (threw) {
            boolean matchesCojoStyle = message != null
                && message.startsWith("GCTA binary not found at: ")
                && message.contains("Place gcta64.exe in the bin/ folder.");
            failures += check("message matches CojoAdapter's exact wording: " + message, matchesCojoStyle);
        }

        deleteRecursive(dir.toFile());
        return failures;
    }

    /**
     * reml() must refuse — with a clear, specific explanation — rather than fabricate a
     * phenotype file or return a mocked h² when the caller passes no phenotype file at all.
     * Uses a fake-but-existing "GCTA binary" (a plain temp file — verify() only checks
     * existence, it never executes anything) and a fake-but-existing GRM prefix so the failure
     * being tested is unambiguously the phenotype-file check, not an earlier one.
     */
    private static int testRemlRefusesWithoutPhenotypeFile() throws Exception {
        Path dir = Files.createTempDirectory("gcta-noPheno-test");
        File fakeGcta = new File(dir.toFile(), "fake_gcta64.exe");
        fakeGcta.createNewFile();
        String grmPrefix = new File(dir.toFile(), "fake_grm").getAbsolutePath();
        new File(grmPrefix + ".grm.bin").createNewFile();
        File runDir = new File(dir.toFile(), "run");

        boolean threw = false;
        String message = null;
        try {
            GctaGremlAdapter.reml(grmPrefix, null, runDir, fakeGcta.getAbsolutePath());
        } catch (IOException e) {
            threw = true;
            message = e.getMessage();
        }

        int failures = 0;
        failures += check("reml() throws IOException when phenoFile is null", threw);
        if (threw) {
            boolean explainsHonestly = message != null
                && message.toLowerCase().contains("phenotype")
                && message.toLowerCase().contains("does not currently ingest");
            failures += check("message honestly explains the missing individual-level phenotype data "
                + "(no fabricated h²): " + message, explainsHonestly);
        }

        deleteRecursive(dir.toFile());
        return failures;
    }

    /** Same honest-refusal contract when a phenotype path IS given but doesn't exist on disk. */
    private static int testRemlRefusesWhenPhenotypeFileMissing() throws Exception {
        Path dir = Files.createTempDirectory("gcta-missingPhenoFile-test");
        File fakeGcta = new File(dir.toFile(), "fake_gcta64.exe");
        fakeGcta.createNewFile();
        String grmPrefix = new File(dir.toFile(), "fake_grm").getAbsolutePath();
        new File(grmPrefix + ".grm.bin").createNewFile();
        File runDir = new File(dir.toFile(), "run");
        File missingPheno = new File(dir.toFile(), "does_not_exist.pheno");

        boolean threw = false;
        String message = null;
        try {
            GctaGremlAdapter.reml(grmPrefix, missingPheno, runDir, fakeGcta.getAbsolutePath());
        } catch (IOException e) {
            threw = true;
            message = e.getMessage();
        }

        int failures = 0;
        failures += check("reml() throws IOException when the given phenotype file doesn't exist", threw);
        if (threw) {
            boolean namesPath = message != null && message.contains(missingPheno.getAbsolutePath());
            failures += check("message names the missing phenotype file path: " + message, namesPath);
        }

        deleteRecursive(dir.toFile());
        return failures;
    }

    private static GctaGremlAdapter.HsqResult findSource(List<GctaGremlAdapter.HsqResult> rows, String source) {
        for (GctaGremlAdapter.HsqResult r : rows) if (source.equals(r.source)) return r;
        return null;
    }

    private static void deleteRecursive(File f) {
        File[] children = f.listFiles();
        if (children != null) for (File c : children) deleteRecursive(c);
        f.delete();
    }

    private static int checkEq(String label, Object actual, Object expected) {
        if (!Objects.equals(actual, expected)) {
            System.out.println("FAIL: " + label + " — expected " + expected + ", got " + actual);
            return 1;
        }
        System.out.println("PASS: " + label + " (" + actual + ")");
        return 0;
    }

    private static int checkClose(String label, double actual, double expected, double tol) {
        if (Math.abs(actual - expected) > tol) {
            System.out.println("FAIL: " + label + " — expected " + expected + ", got " + actual);
            return 1;
        }
        System.out.println("PASS: " + label + " (" + actual + ")");
        return 0;
    }

    private static int check(String label, boolean cond) {
        if (!cond) {
            System.out.println("FAIL: " + label);
            return 1;
        }
        System.out.println("PASS: " + label);
        return 0;
    }

    private static int checkNotNull(String label, Object o) {
        if (o == null) {
            System.out.println("FAIL: " + label + " — was null");
            return 1;
        }
        System.out.println("PASS: " + label);
        return 0;
    }
}
