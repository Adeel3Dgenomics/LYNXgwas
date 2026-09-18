import java.io.*;
import java.nio.file.*;
import java.util.*;

/**
 * Standalone regression test (no external test framework — this project has none) for
 * SusieAdapter. R is not installed on this machine, so susie_run.R itself is never executed —
 * only the Java-side input-generation logic is tested: susie_input.tsv formatting, ref-panel-ID
 * substitution via posToRefId (readBimIds), numeric filtering of bad beta/se rows, and the
 * required-sampleN guard. SusieAdapter has no separate result-parsing method (result.tsv /
 * susie_diag.json are only ever produced by susie_run.R itself, and PluginEngine's generic
 * mapOutput() reads them declaratively — there is nothing SusieAdapter-specific to unit test
 * there), so there is no parsing test here; see the class comment on ColocAdapter for the
 * analogous "generic mapOutput()" pattern.
 *
 * Note on ID fallback: unlike CojoAdapter/MagmaAdapter (which must never fall back to the GWAS's
 * own SNP id — see MagmaAdapterTest/CojoAdapterTest), SusieAdapter DOES fall back to the GWAS's
 * own id when a SNP has no ref-panel BIM match (`String useId = refId != null ? refId : f[0];`).
 * This is intentional, not a bug: susie_run.R does its own `intersect(df$SNP, bim_win$SNP)`
 * against the real PLINK bim IDs afterward (see the R script body), so a fallback id that
 * doesn't match any real bim SNP is simply filtered out downstream in R rather than silently
 * matched to the wrong genotype column the way an unmatched COJO/MAGMA id could be. This test
 * documents that current fallback behavior rather than treating it as a defect.
 */
public class SusieAdapterTest {

    public static void main(String[] args) throws Exception {
        int failures = 0;
        failures += testPrepareRunWritesInputTsv();
        failures += testPrepareRunFallsBackToOwnIdWhenBimMissingEntirely();
        failures += testPrepareRunThrowsWhenSampleNZero();

        if (failures == 0) {
            System.out.println("PASS: all SusieAdapter tests passed");
        } else {
            System.out.println("FAIL: " + failures + " test(s) failed");
            System.exit(1);
        }
    }

    /**
     * matched_ref.bim (chr, snpId, cM, pos, A1, A2) has entries for chr1:1000000 (refid_1) and
     * chr1:2000000 (refid_2). harmonized_gwas.tsv rows, in file order:
     *   gwasA  chr1:1000000  beta=0.5   se=0.05  P=0.001  -> matches refid_1 -> id substituted
     *   gwasB  chr1:2000000  beta=-0.3  se=0.04  P=0.002  -> matches refid_2 -> id substituted
     *   gwasC  chr1:3000000  beta=0.2   se=0.03  P=0.003  -> NO bim match -> id falls back to
     *          its own GWAS id "gwasC" (see class comment — intentional, not the Fix-2 rule)
     *   gwasD  chr1:4000000  beta=NA    se=0.02           -> dropped (beta.equals("NA"))
     *   gwasE  chr1:5000000  beta=abc   se=0.02           -> dropped (NumberFormatException)
     *   gwasF  chr1:6000000  beta=0.1   se=0              -> dropped (se<=0)
     *   gwasG  chr1:7000000  beta=0.1   se=-0.01          -> dropped (se<=0, negative)
     *   gwasH  chr1:8000000  beta=NaN   se=0.02           -> Double.parseDouble("NaN") parses
     *          successfully to Double.NaN (does NOT throw), so this exercises the SEPARATE
     *          Double.isNaN(b) guard, not the catch block — must be dropped, not crash, and is
     *          distinct from gwasE's plain-garbage-string case.
     * Expected: 3 kept rows (gwasA->refid_1, gwasB->refid_2, gwasC->gwasC), 5 dropped.
     */
    private static int testPrepareRunWritesInputTsv() throws Exception {
        Path dir = Files.createTempDirectory("susie-happy-test");
        File harmonizedDir = new File(dir.toFile(), "harmonized");
        File matchedDir = new File(dir.toFile(), "matched");
        File runDir = new File(dir.toFile(), "run");
        harmonizedDir.mkdirs();
        matchedDir.mkdirs();

        try (PrintWriter pw = new PrintWriter(new FileWriter(new File(matchedDir, "matched_ref.bim")))) {
            pw.println("1\trefid_1\t0\t1000000\tA\tG");
            pw.println("1\trefid_2\t0\t2000000\tC\tT");
        }

        try (PrintWriter pw = new PrintWriter(new FileWriter(new File(harmonizedDir, "harmonized_gwas.tsv")))) {
            pw.println("snp_id\tchr\tpos\tea\tnea\tpvalue\tbeta\tse");
            pw.println("gwasA\t1\t1000000\tA\tG\t0.001\t0.5\t0.05");
            pw.println("gwasB\t1\t2000000\tC\tT\t0.002\t-0.3\t0.04");
            pw.println("gwasC\t1\t3000000\tA\tG\t0.003\t0.2\t0.03");
            pw.println("gwasD\t1\t4000000\tA\tG\t0.004\tNA\t0.02");
            pw.println("gwasE\t1\t5000000\tA\tG\t0.005\tabc\t0.02");
            pw.println("gwasF\t1\t6000000\tA\tG\t0.006\t0.1\t0");
            pw.println("gwasG\t1\t7000000\tA\tG\t0.007\t0.1\t-0.01");
            pw.println("gwasH\t1\t8000000\tA\tG\t0.008\tNaN\t0.02");
        }

        SusieAdapter.prepareRun(harmonizedDir, matchedDir, runDir, /*sampleN*/5000, /*maxCausal*/10,
            /*coverage*/0.95, /*ldShrink*/0.010, /*windowKb*/500);

        int failures = 0;
        File inputTsv = new File(runDir, "susie_input.tsv");
        failures += check("susie_input.tsv exists", inputTsv.exists());
        List<String> lines = Files.readAllLines(inputTsv.toPath());
        failures += checkEq("susie_input.tsv line count (header + 3 kept rows)", lines.size(), 4);
        failures += checkEq("header", lines.get(0), "SNP\tBP\tA1\tA2\tBETA\tSE\tP");

        Map<String, String[]> rows = new LinkedHashMap<>();
        for (int i = 1; i < lines.size(); i++) {
            String[] f = lines.get(i).split("\t", -1);
            rows.put(f[0], f);
        }
        failures += check("kept exactly {refid_1, refid_2, gwasC}: " + rows.keySet(),
            rows.keySet().equals(new LinkedHashSet<>(Arrays.asList("refid_1", "refid_2", "gwasC"))));

        String[] rA = rows.get("refid_1");
        if (rA != null) {
            failures += checkEq("refid_1 (gwasA) BP", rA[1], "1000000");
            failures += checkEq("refid_1 (gwasA) A1", rA[2], "A");
            failures += checkEq("refid_1 (gwasA) A2", rA[3], "G");
            failures += checkClose("refid_1 (gwasA) BETA", Double.parseDouble(rA[4]), 0.5, 1e-9);
            failures += checkClose("refid_1 (gwasA) SE", Double.parseDouble(rA[5]), 0.05, 1e-9);
            failures += checkClose("refid_1 (gwasA) P", Double.parseDouble(rA[6]), 0.001, 1e-9);
        }
        String[] rB = rows.get("refid_2");
        if (rB != null) {
            failures += checkClose("refid_2 (gwasB) BETA (negative, preserved verbatim)",
                Double.parseDouble(rB[4]), -0.3, 1e-9);
        }
        String[] rC = rows.get("gwasC");
        failures += check("gwasC row uses its OWN snp_id as fallback (no ref-panel match at chr1:3000000)",
            rC != null);

        boolean anyBadRow = rows.keySet().stream().anyMatch(k ->
            Arrays.asList("gwasD", "gwasE", "gwasF", "gwasG", "gwasH").contains(k));
        failures += check("none of the 5 bad-stat rows (NA/non-numeric/se<=0/NaN) were written",
            !anyBadRow);

        File rScript = new File(runDir, "susie_run.R");
        failures += check("susie_run.R generated", rScript.exists());
        String rContent = new String(Files.readAllBytes(rScript.toPath()), "UTF-8");
        failures += check("R script substitutes n_samples=5000", rContent.contains("n_samples <- 5000"));
        failures += check("R script substitutes L=10", rContent.contains("L <- 10"));
        failures += check("R script substitutes coverage=0.95", rContent.contains("coverage <- 0.95"));
        failures += check("R script substitutes ld_shrink=0.010", rContent.contains("ld_shrink <- 0.010"));
        failures += check("R script substitutes window_kb=500", rContent.contains("window_kb <- 500"));
        failures += check("R script uses estimate_residual_variance=FALSE (ref-panel LD, not in-sample)",
            rContent.contains("estimate_residual_variance=FALSE"));

        deleteRecursive(dir.toFile());
        return failures;
    }

    /** matched_ref.bim absent entirely -> readBimIds() returns an empty map gracefully (no crash), every kept row falls back to its own GWAS id. */
    private static int testPrepareRunFallsBackToOwnIdWhenBimMissingEntirely() throws Exception {
        Path dir = Files.createTempDirectory("susie-nobim-test");
        File harmonizedDir = new File(dir.toFile(), "harmonized");
        File matchedDir = new File(dir.toFile(), "matched"); // matched_ref.bim deliberately never written
        File runDir = new File(dir.toFile(), "run");
        harmonizedDir.mkdirs();
        matchedDir.mkdirs();

        try (PrintWriter pw = new PrintWriter(new FileWriter(new File(harmonizedDir, "harmonized_gwas.tsv")))) {
            pw.println("snp_id\tchr\tpos\tea\tnea\tpvalue\tbeta\tse");
            pw.println("gwasOnly\t1\t1000000\tA\tG\t0.001\t0.4\t0.05");
        }

        SusieAdapter.prepareRun(harmonizedDir, matchedDir, runDir, 5000, 10, 0.95, 0.01, 500);

        int failures = 0;
        List<String> lines = Files.readAllLines(new File(runDir, "susie_input.tsv").toPath());
        failures += checkEq("line count (header + 1 row)", lines.size(), 2);
        failures += check("row uses GWAS's own id when bim file is entirely absent: " + lines.get(1),
            lines.get(1).startsWith("gwasOnly\t"));

        deleteRecursive(dir.toFile());
        return failures;
    }

    private static int testPrepareRunThrowsWhenSampleNZero() throws Exception {
        Path dir = Files.createTempDirectory("susie-noN-test");
        File harmonizedDir = new File(dir.toFile(), "harmonized");
        File matchedDir = new File(dir.toFile(), "matched");
        File runDir = new File(dir.toFile(), "run");

        boolean threw = false;
        String message = null;
        try {
            SusieAdapter.prepareRun(harmonizedDir, matchedDir, runDir, 0, 10, 0.95, 0.01, 500);
        } catch (IOException e) {
            threw = true;
            message = e.getMessage();
        }

        int failures = 0;
        failures += check("prepareRun() throws IOException when sampleN<=0", threw);
        if (threw) {
            failures += check("message names SuSiE and N requirement: " + message,
                message != null && message.contains("Sample size (N) is required for SuSiE"));
        }
        deleteRecursive(dir.toFile());
        return failures;
    }

    private static void deleteRecursive(File f) {
        File[] children = f.listFiles();
        if (children != null) for (File c : children) deleteRecursive(c);
        f.delete();
    }

    private static int checkEq(String label, Object actual, Object expected) {
        if (!Objects.equals(actual, expected)) {
            System.out.println("FAIL: " + label + " — expected [" + expected + "], got [" + actual + "]");
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
}
