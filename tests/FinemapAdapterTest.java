import java.io.*;
import java.nio.file.*;
import java.util.*;

/**
 * Standalone regression test (no external test framework — this project has none) for
 * FinemapAdapter. FINEMAP has no official Windows build, so run_finemap.R / the FINEMAP binary
 * itself are never executed — only the Java-side z-file/ld/master-file generation is tested.
 *
 * Primary focus, per this task's brief, is the allele-orientation logic in orientedRefFreq
 * (frequency/beta sign must track the same allele as the reference panel's own A1, exactly the
 * class of bug this project has already been bitten by once — see CojoAdapter's Fix 1):
 *   - GWAS EA matches ref A1 directly       -> ref freq used as-is
 *   - GWAS EA matches ref A2 (flipped)      -> ref freq must be 1 - refFreq
 *   - GWAS EA matches NEITHER ref allele    -> orientedRefFreq returns null, don't guess (freq
 *                                              stays at the 0.5 placeholder default rather than
 *                                              being silently misassigned)
 * This orientation logic is only exercised when the GWAS's own MAF column is NA/out-of-range,
 * since GWAS-provided MAF is always preferred when present and valid.
 */
public class FinemapAdapterTest {

    public static void main(String[] args) throws Exception {
        int failures = 0;
        failures += testPrepareRunOrientsRefFrequencyAndFiltersRows();
        failures += testPrepareRunThrowsWhenSampleNZero();
        failures += testZFileOrderFollowsGwasFileIterationOrderNotLdSnpOrderIndex();

        if (failures == 0) {
            System.out.println("PASS: all FinemapAdapter tests passed");
        } else {
            System.out.println("FAIL: " + failures + " test(s) failed");
            System.exit(1);
        }
    }

    /**
     * Fixture:
     *   matched_ref.bim (chr, snpId, cM, pos, A1, A2):
     *     1  refid_1  0  1000000  A  G   (A1=A, A2=G)
     *     1  refid_2  0  2000000  C  T   (A1=C, A2=T)
     *     1  refid_3  0  3000000  A  C   (A1=A, A2=C)
     *   ref_freq.tsv (header + snp_id, A1, freq — freq is for the panel's OWN A1):
     *     refid_1  A  0.30
     *     refid_2  C  0.65
     *     refid_3  A  0.10
     *   ld_snp_order.txt (chr:pos:a1:a2, the real LdMatrixComputer format — FinemapAdapter only
     *   reads parts[0]+":"+parts[1], i.e. "chr:pos"): all three positions present, so none of the
     *   three rows below are dropped by the posLookup gate.
     *
     *   harmonized_gwas.tsv rows:
     *     snpX chr1:1000000 EA=A NEA=G maf=0.28 (valid, non-NA)
     *       -> GWAS's own MAF wins outright; freq = 0.28 exactly, ref panel never consulted.
     *     snpY chr1:2000000 EA=T NEA=C maf=NA
     *       -> falls back to ref: EA "T" == ref A2 ("T") -> FLIPPED. Hand-computed:
     *          freq = 1 - refFreq(refid_2) = 1 - 0.65 = 0.35.
     *     snpZ chr1:3000000 EA=T NEA=G maf=NA
     *       -> falls back to ref: EA "T" matches neither ref A1 ("A") nor A2 ("C") for refid_3
     *          -> orientedRefFreq returns null -> freq stays at the 0.5 placeholder default
     *          (documents current "don't guess" behavior — not miscomputed, just not corrected).
     */
    private static int testPrepareRunOrientsRefFrequencyAndFiltersRows() throws Exception {
        Path dir = Files.createTempDirectory("finemap-orient-test");
        File harmonizedDir = new File(dir.toFile(), "harmonized");
        File ldDir = new File(dir.toFile(), "ld");
        File matchedDir = new File(dir.toFile(), "matched");
        File runDir = new File(dir.toFile(), "run");
        harmonizedDir.mkdirs();
        ldDir.mkdirs();
        matchedDir.mkdirs();

        try (PrintWriter pw = new PrintWriter(new FileWriter(new File(matchedDir, "matched_ref.bim")))) {
            pw.println("1\trefid_1\t0\t1000000\tA\tG");
            pw.println("1\trefid_2\t0\t2000000\tC\tT");
            pw.println("1\trefid_3\t0\t3000000\tA\tC");
        }
        try (PrintWriter pw = new PrintWriter(new FileWriter(new File(matchedDir, "ref_freq.tsv")))) {
            pw.println("snp_id\tA1\tfreq");
            pw.println("refid_1\tA\t0.30");
            pw.println("refid_2\tC\t0.65");
            pw.println("refid_3\tA\t0.10");
        }
        try (PrintWriter pw = new PrintWriter(new FileWriter(new File(ldDir, "ld_snp_order.txt")))) {
            pw.println("1:1000000:A:G");
            pw.println("1:2000000:C:T");
            pw.println("1:3000000:A:C");
        }
        // 3x3 placeholder LD matrix (values are arbitrary — this test does not check LD content
        // correctness, only that the file is copied through as space-delimited).
        try (PrintWriter pw = new PrintWriter(new FileWriter(new File(ldDir, "ld_r.matrix")))) {
            pw.println("1.0\t0.2\t0.1");
            pw.println("0.2\t1.0\t0.3");
            pw.println("0.1\t0.3\t1.0");
        }
        try (PrintWriter pw = new PrintWriter(new FileWriter(new File(harmonizedDir, "harmonized_gwas.tsv")))) {
            pw.println("snp_id\tchr\tpos\tea\tnea\tpvalue\tbeta\tse\tor\tn\tmaf\tinfo\trsid\tvarid");
            pw.println("snpX\t1\t1000000\tA\tG\t0.001\t0.40\t0.05\t1.5\t1000\t0.28\t0.9\trsX\tvX");
            pw.println("snpY\t1\t2000000\tT\tC\t0.002\t-0.20\t0.04\t0.8\t1000\tNA\t0.9\trsY\tvY");
            pw.println("snpZ\t1\t3000000\tT\tG\t0.003\t0.10\t0.02\t1.1\t1000\tNA\t0.9\trsZ\tvZ");
        }

        FinemapAdapter.prepareRun(harmonizedDir, ldDir, matchedDir, runDir, /*sampleN*/1000, /*maxCausal*/5);

        int failures = 0;
        File zFile = new File(runDir, "finemap.z");
        failures += check("finemap.z exists", zFile.exists());
        List<String> lines = Files.readAllLines(zFile.toPath());
        failures += checkEq("finemap.z line count (header + 3 rows)", lines.size(), 4);
        failures += checkEq("header", lines.get(0), "rsid chromosome position allele1 allele2 maf beta se");

        Map<String, String[]> rows = new LinkedHashMap<>();
        for (int i = 1; i < lines.size(); i++) {
            String[] f = lines.get(i).split(" ", -1);
            rows.put(f[0], f);
        }

        String[] rX = rows.get("refid_1");
        failures += check("snpX resolved to ref panel id refid_1", rX != null);
        if (rX != null) {
            failures += checkClose("snpX freq = own GWAS MAF (0.28), ref panel never consulted",
                Double.parseDouble(rX[5]), 0.28, 1e-9);
        }

        String[] rY = rows.get("refid_2");
        failures += check("snpY resolved to ref panel id refid_2", rY != null);
        if (rY != null) {
            failures += checkClose(
                "snpY freq FLIPPED: EA=T matches ref A2, so freq = 1 - refFreq(0.65) = 0.35",
                Double.parseDouble(rY[5]), 0.35, 1e-9);
            failures += checkClose("snpY beta unchanged (-0.20) — only frequency is oriented, not beta sign, "
                + "because beta already refers to the harmonized EA in harmonized_gwas.tsv",
                Double.parseDouble(rY[6]), -0.20, 1e-9);
        }

        String[] rZ = rows.get("refid_3");
        failures += check("snpZ resolved to ref panel id refid_3", rZ != null);
        if (rZ != null) {
            failures += checkClose(
                "snpZ freq stays at the 0.5 default: EA=T matches neither ref A1(A) nor A2(C) "
                + "for refid_3, orientedRefFreq returns null rather than guessing",
                Double.parseDouble(rZ[5]), 0.5, 1e-9);
        }

        File ldOut = new File(runDir, "finemap.ld");
        failures += check("finemap.ld exists (LD matrix copied through)", ldOut.exists());
        List<String> ldLines = Files.readAllLines(ldOut.toPath());
        failures += checkEq("finemap.ld line count", ldLines.size(), 3);
        failures += check("finemap.ld is space-delimited, not tab-delimited: " + ldLines.get(0),
            ldLines.get(0).equals("1.0 0.2 0.1"));

        File masterFile = new File(runDir, "finemap.master");
        failures += check("finemap.master exists", masterFile.exists());
        List<String> masterLines = Files.readAllLines(masterFile.toPath());
        failures += checkEq("finemap.master line count (header + 1 row)", masterLines.size(), 2);
        failures += check("finemap.master's n_samples field is 1000: " + masterLines.get(1),
            masterLines.get(1).endsWith(";1000"));

        deleteRecursive(dir.toFile());
        return failures;
    }

    private static int testPrepareRunThrowsWhenSampleNZero() throws Exception {
        Path dir = Files.createTempDirectory("finemap-noN-test");
        File harmonizedDir = new File(dir.toFile(), "harmonized");
        File ldDir = new File(dir.toFile(), "ld");
        File matchedDir = new File(dir.toFile(), "matched");
        File runDir = new File(dir.toFile(), "run");

        boolean threw = false;
        String message = null;
        try {
            FinemapAdapter.prepareRun(harmonizedDir, ldDir, matchedDir, runDir, 0, 5);
        } catch (IOException e) {
            threw = true;
            message = e.getMessage();
        }

        int failures = 0;
        failures += check("prepareRun() throws IOException when sampleN<=0", threw);
        if (threw) {
            failures += check("message names FINEMAP and N requirement: " + message,
                message != null && message.contains("Sample size (N) is required for FINEMAP"));
        }
        deleteRecursive(dir.toFile());
        return failures;
    }

    /**
     * DOCUMENTS CURRENT BEHAVIOR — flagged as an uncertain correctness risk in the final report,
     * not fixed here (per this task's instruction to document rather than guess when unsure).
     *
     * finemap.z rows are written by iterating harmonized_gwas.tsv in ITS OWN file order and
     * simply skipping rows whose chr:pos isn't in ld_snp_order.txt (a set-membership check via
     * posLookup). Nothing re-sorts the kept rows into ld_snp_order.txt's own order before writing
     * them to finemap.z. FINEMAP (and this adapter's own run_finemap.R ABF fallback) both assume
     * finemap.z row i lines up positionally with finemap.ld row/column i. In this pipeline the
     * two files are expected to already agree (both are ultimately position-sorted upstream —
     * ld_snp_order.txt from a position-sorted BIM window, harmonized_gwas.tsv from a
     * position-filtered GWAS stream) — but nothing in FinemapAdapter itself enforces or verifies
     * that agreement. This test constructs a harmonized_gwas.tsv deliberately out of position
     * order (as could happen with a not-perfectly-sorted input GWAS file) and shows that
     * finemap.z comes out in GWAS-file order, not ld_snp_order.txt order — i.e. confirms the
     * ordering is NOT independently guaranteed by this method.
     */
    private static int testZFileOrderFollowsGwasFileIterationOrderNotLdSnpOrderIndex() throws Exception {
        Path dir = Files.createTempDirectory("finemap-order-test");
        File harmonizedDir = new File(dir.toFile(), "harmonized");
        File ldDir = new File(dir.toFile(), "ld");
        File matchedDir = new File(dir.toFile(), "matched");
        File runDir = new File(dir.toFile(), "run");
        harmonizedDir.mkdirs();
        ldDir.mkdirs();
        matchedDir.mkdirs();

        try (PrintWriter pw = new PrintWriter(new FileWriter(new File(matchedDir, "matched_ref.bim")))) {
            pw.println("1\trefid_1\t0\t1000000\tA\tG");
            pw.println("1\trefid_2\t0\t2000000\tA\tG");
        }
        try (PrintWriter pw = new PrintWriter(new FileWriter(new File(matchedDir, "ref_freq.tsv")))) {
            pw.println("snp_id\tA1\tfreq");
            pw.println("refid_1\tA\t0.30");
            pw.println("refid_2\tA\t0.40");
        }
        // ld_snp_order.txt in ASCENDING position order (position 1000000 first, 2000000 second) —
        // this is the order FINEMAP's LD matrix rows/columns will be in.
        try (PrintWriter pw = new PrintWriter(new FileWriter(new File(ldDir, "ld_snp_order.txt")))) {
            pw.println("1:1000000:A:G");
            pw.println("1:2000000:A:G");
        }
        try (PrintWriter pw = new PrintWriter(new FileWriter(new File(ldDir, "ld_r.matrix")))) {
            pw.println("1.0\t0.5");
            pw.println("0.5\t1.0");
        }
        // harmonized_gwas.tsv deliberately lists the HIGHER position (2000000) BEFORE the lower
        // one (1000000) — the reverse of ld_snp_order.txt.
        try (PrintWriter pw = new PrintWriter(new FileWriter(new File(harmonizedDir, "harmonized_gwas.tsv")))) {
            pw.println("snp_id\tchr\tpos\tea\tnea\tpvalue\tbeta\tse\tor\tn\tmaf\tinfo\trsid\tvarid");
            pw.println("snpHigh\t1\t2000000\tA\tG\t0.001\t0.5\t0.05\t1.6\t1000\t0.3\t0.9\trsH\tvH");
            pw.println("snpLow\t1\t1000000\tA\tG\t0.002\t0.2\t0.04\t1.2\t1000\t0.3\t0.9\trsL\tvL");
        }

        FinemapAdapter.prepareRun(harmonizedDir, ldDir, matchedDir, runDir, 1000, 5);

        List<String> lines = Files.readAllLines(new File(runDir, "finemap.z").toPath());
        int failures = 0;
        failures += checkEq("finemap.z line count (header + 2 rows)", lines.size(), 3);
        boolean firstRowIsHighPos = lines.get(1).startsWith("refid_2 "); // refid_2 = pos 2000000
        failures += check("CURRENT BEHAVIOR (documented, not asserted-safe): finemap.z's first "
            + "data row is refid_2 (position 2000000), matching harmonized_gwas.tsv's own row "
            + "order — NOT re-sorted to ld_snp_order.txt's position order, where position "
            + "1000000 (refid_1) comes first. If finemap.ld's row/column order is taken from "
            + "ld_snp_order.txt (position-ascending) while finemap.z is NOT re-sorted to match, "
            + "FINEMAP would silently pair each SNP's z-score with the WRONG row/column of the "
            + "LD matrix whenever the two input files disagree on order: " + lines.get(1),
            firstRowIsHighPos);

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
