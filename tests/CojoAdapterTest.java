import java.io.*;
import java.nio.file.*;
import java.util.*;

/**
 * Standalone regression test (no external test framework — this project has none) for
 * CojoAdapter v2.1's Java-side input generation: the GCTA-binary-not-found guard (shared
 * GctaBinaryResolver, same as GctaGremlAdapterTest), the required ref_freq.tsv guard, the
 * sample-size guard (global N or per-SNP N), and the per-SNP audit trail covering all six
 * documented fixes' Java-side halves:
 *   Fix 1 — allele/beta/freq orientation verified per SNP (kept vs freq_flipped vs dropped)
 *   Fix 2 — SNPs without a ref panel BIM id are DROPPED, never fall back to the GWAS's own id
 *   Fix 3 — per-SNP N accepted without requiring global N (and used per-row when present)
 *   Fix 4 — LD-GWAS consistency diagnostic (consistency_summary.json) substituted into run_cojo.R
 *   Fix 5 — freq-discrepancy-warning parsing is generated INTO run_cojo.R (real GCTA .jma.cojo
 *           freq/freq_geom columns only exist after an actual GCTA run, so this fix's own logic
 *           lives in the generated R text, not in Java — this test verifies that R code is
 *           actually wired into the script rather than executing it)
 *   Fix 6 — not exercised here (requires parsing an actual empty .jma.cojo from a real GCTA run)
 */
public class CojoAdapterTest {

    public static void main(String[] args) throws Exception {
        int failures = 0;
        failures += testPrepareRunThrowsWhenGctaBinaryMissing();
        failures += testPrepareRunThrowsWhenRefFreqMissing();
        failures += testPrepareRunThrowsWhenNoSampleSizeAvailable();
        failures += testPrepareRunWritesMaFileWithOrientationAndDrops();
        failures += testPrepareRunWiresLdConsistencyIntoRScript();
        failures += testPrepareRunWiresFreqDiscrepancyParsingIntoRScript();

        if (failures == 0) {
            System.out.println("PASS: all CojoAdapter tests passed");
        } else {
            System.out.println("FAIL: " + failures + " test(s) failed");
            System.exit(1);
        }
    }

    private static int testPrepareRunThrowsWhenGctaBinaryMissing() throws Exception {
        Path dir = Files.createTempDirectory("cojo-nogcta-test");
        File harmonizedDir = new File(dir.toFile(), "harmonized");
        File matchedDir = new File(dir.toFile(), "matched");
        File runDir = new File(dir.toFile(), "run");
        String bogusGcta = new File(dir.toFile(), "no_such_gcta64.exe").getAbsolutePath();

        boolean threw = false;
        String message = null;
        try {
            CojoAdapter.prepareRun(harmonizedDir, matchedDir, runDir, 1000, 5e-8, 0.9, bogusGcta);
        } catch (IOException e) {
            threw = true;
            message = e.getMessage();
        }

        int failures = 0;
        failures += check("prepareRun() throws IOException when GCTA binary is missing", threw);
        if (threw) {
            failures += check("message matches GctaBinaryResolver's exact wording: " + message,
                message != null && message.startsWith("GCTA binary not found at: ")
                    && message.contains("Place gcta64.exe in the bin/ folder."));
        }
        deleteRecursive(dir.toFile());
        return failures;
    }

    private static int testPrepareRunThrowsWhenRefFreqMissing() throws Exception {
        Path dir = Files.createTempDirectory("cojo-norefFreq-test");
        File harmonizedDir = new File(dir.toFile(), "harmonized");
        File matchedDir = new File(dir.toFile(), "matched");
        File runDir = new File(dir.toFile(), "run");
        harmonizedDir.mkdirs();
        matchedDir.mkdirs(); // ref_freq.tsv deliberately never written
        File fakeGcta = new File(dir.toFile(), "fake_gcta64.exe");
        fakeGcta.createNewFile();

        try (PrintWriter pw = new PrintWriter(new FileWriter(new File(harmonizedDir, "harmonized_gwas.tsv")))) {
            pw.println("snp_id\tchr\tpos\tea\tnea\tpvalue\tbeta\tse\tn");
            pw.println("snpA\t1\t1000000\tA\tG\t0.001\t0.3\t0.05\t1000");
        }

        boolean threw = false;
        String message = null;
        try {
            CojoAdapter.prepareRun(harmonizedDir, matchedDir, runDir, 1000, 5e-8, 0.9, fakeGcta.getAbsolutePath());
        } catch (IOException e) {
            threw = true;
            message = e.getMessage();
        }

        int failures = 0;
        failures += check("prepareRun() throws IOException when ref_freq.tsv is missing", threw);
        if (threw) {
            failures += check("message names ref_freq.tsv: " + message,
                message != null && message.contains("ref_freq.tsv not found in matched directory"));
        }
        deleteRecursive(dir.toFile());
        return failures;
    }

    /** Fix 3: sampleN<=0 AND no usable per-SNP N anywhere in the file -> must refuse (never guess an N). */
    private static int testPrepareRunThrowsWhenNoSampleSizeAvailable() throws Exception {
        Path dir = Files.createTempDirectory("cojo-noN-test");
        File harmonizedDir = new File(dir.toFile(), "harmonized");
        File matchedDir = new File(dir.toFile(), "matched");
        File runDir = new File(dir.toFile(), "run");
        harmonizedDir.mkdirs();
        matchedDir.mkdirs();
        File fakeGcta = new File(dir.toFile(), "fake_gcta64.exe");
        fakeGcta.createNewFile();

        try (PrintWriter pw = new PrintWriter(new FileWriter(new File(matchedDir, "ref_freq.tsv")))) {
            pw.println("snp_id\tA1\tfreq");
            pw.println("refid_1\tA\t0.4");
        }
        // No "n" column at all in the header -> anySnpHasN stays false.
        try (PrintWriter pw = new PrintWriter(new FileWriter(new File(harmonizedDir, "harmonized_gwas.tsv")))) {
            pw.println("snp_id\tchr\tpos\tea\tnea\tpvalue\tbeta\tse");
            pw.println("snpA\t1\t1000000\tA\tG\t0.001\t0.3\t0.05");
        }

        boolean threw = false;
        String message = null;
        try {
            CojoAdapter.prepareRun(harmonizedDir, matchedDir, runDir, 0, 5e-8, 0.9, fakeGcta.getAbsolutePath());
        } catch (IOException e) {
            threw = true;
            message = e.getMessage();
        }

        int failures = 0;
        failures += check("prepareRun() throws IOException when sampleN<=0 and no per-SNP N column exists", threw);
        if (threw) {
            failures += check("message names COJO and N requirement: " + message,
                message != null && message.contains("Sample size (N) is required for COJO"));
        }
        deleteRecursive(dir.toFile());
        return failures;
    }

    /**
     * Full fixture exercising Fix 1 (orientation), Fix 2 (no-fallback drop), and Fix 3 (per-SNP N).
     *
     * matched_ref.bim (chr, snpId, cM, pos, A1, A2):
     *   1  refid_1  0  1000000  A  G   (A1=A, A2=G)
     *   1  refid_2  0  2000000  C  T   (A1=C, A2=T)
     *   1  refid_4  0  4000000  A  C   (A1=A, A2=C)   [refid_3 deliberately absent -> tests Fix 2]
     *   1  refid_5  0  5000000  A  G   (A1=A, A2=G)   [present in bim, absent from ref_freq.tsv]
     *
     * ref_freq.tsv (freq of the ref panel's OWN A1):
     *   refid_1  A  0.40
     *   refid_2  C  0.70
     *   refid_4  A  0.55
     *   (refid_5 has no entry -> dropped_no_freq)
     *
     * harmonized_gwas.tsv rows (header includes "n" so Fix-3's per-SNP N path is exercised):
     *   snpKeep    chr1:1000000 EA=A NEA=G beta=0.30  se=0.05 N=500
     *     -> EA==refA1(A) -> freqEA = refFreq(refid_1) = 0.40 exactly ("kept"); per-SNP N=500 used.
     *   snpFlip    chr1:2000000 EA=T NEA=C beta=-0.25 se=0.04 N=NA
     *     -> EA==refA2(T) -> FLIPPED: freqEA = 1 - refFreq(refid_2) = 1 - 0.70 = 0.30
     *        ("freq_flipped"); N is NA in this row -> falls back to the global sampleN=1000.
     *   snpNoRef       chr1:3000000 -> no chr:pos entry anywhere in matched_ref.bim -> Fix 2:
     *        dropped_no_ref_id (never falls back to using "snpNoRef" itself as the .ma SNP id).
     *   snpNoFreq      chr1:5000000 -> refid_5 exists in the bim but has no ref_freq.tsv entry
     *        -> dropped_no_freq.
     *   snpOrientBad   chr1:4000000 EA=G NEA=T -> refid_4's A1=A, A2=C; EA "G" matches neither
     *        -> dropped_orientation.
     *   snpNaNBeta     chr1:6000000 beta="NaN" -> Double.parseDouble("NaN") parses successfully
     *        to Double.NaN (no exception) -> caught by the SEPARATE Double.isNaN(betaVal) check
     *        -> dropped_stats (must be skipped, not crash).
     *   snpGarbageBeta chr1:7000000 beta="abc" -> NumberFormatException -> dropped_stats (the
     *        other dropped_stats code path, distinct from the NaN-string case above).
     *
     * Expected: written=2 (snpKeep, snpFlip); droppedNoRefId=1; droppedNoFreq=1;
     * droppedOrientation=1; droppedStats=2; freqFlipped=1; return value == written == 2.
     */
    private static int testPrepareRunWritesMaFileWithOrientationAndDrops() throws Exception {
        Path dir = Files.createTempDirectory("cojo-happy-test");
        File harmonizedDir = new File(dir.toFile(), "harmonized");
        File matchedDir = new File(dir.toFile(), "matched");
        File runDir = new File(dir.toFile(), "run");
        harmonizedDir.mkdirs();
        matchedDir.mkdirs();
        File fakeGcta = new File(dir.toFile(), "fake_gcta64.exe");
        fakeGcta.createNewFile();

        try (PrintWriter pw = new PrintWriter(new FileWriter(new File(matchedDir, "matched_ref.bim")))) {
            pw.println("1\trefid_1\t0\t1000000\tA\tG");
            pw.println("1\trefid_2\t0\t2000000\tC\tT");
            pw.println("1\trefid_4\t0\t4000000\tA\tC");
            pw.println("1\trefid_5\t0\t5000000\tA\tG");
        }
        try (PrintWriter pw = new PrintWriter(new FileWriter(new File(matchedDir, "ref_freq.tsv")))) {
            pw.println("snp_id\tA1\tfreq");
            pw.println("refid_1\tA\t0.40");
            pw.println("refid_2\tC\t0.70");
            pw.println("refid_4\tA\t0.55");
        }
        try (PrintWriter pw = new PrintWriter(new FileWriter(new File(harmonizedDir, "harmonized_gwas.tsv")))) {
            pw.println("snp_id\tchr\tpos\tea\tnea\tpvalue\tbeta\tse\tn");
            pw.println("snpKeep\t1\t1000000\tA\tG\t0.001\t0.30\t0.05\t500");
            pw.println("snpFlip\t1\t2000000\tT\tC\t0.002\t-0.25\t0.04\tNA");
            pw.println("snpNoRef\t1\t3000000\tA\tG\t0.003\t0.10\t0.02\t300");
            pw.println("snpNoFreq\t1\t5000000\tA\tG\t0.004\t0.15\t0.03\t300");
            pw.println("snpOrientBad\t1\t4000000\tG\tT\t0.005\t0.20\t0.02\t300");
            pw.println("snpNaNBeta\t1\t6000000\tA\tG\t0.006\tNaN\t0.02\t300");
            pw.println("snpGarbageBeta\t1\t7000000\tA\tG\t0.007\tabc\t0.02\t300");
        }

        int written = CojoAdapter.prepareRun(harmonizedDir, matchedDir, runDir, /*sampleN*/1000,
            /*pCutoff*/5e-8, /*collinear*/0.9, fakeGcta.getAbsolutePath());

        int failures = 0;
        failures += checkEq("return value == 2 kept SNPs", written, 2);

        List<String> maLines = Files.readAllLines(new File(runDir, "cojo_input.ma").toPath());
        failures += checkEq(".ma line count (header + 2 kept rows)", maLines.size(), 3);
        failures += checkEq(".ma header", maLines.get(0), "SNP\tA1\tA2\tfreq\tb\tse\tp\tN");

        Map<String, String[]> maRows = new LinkedHashMap<>();
        for (int i = 1; i < maLines.size(); i++) {
            String[] f = maLines.get(i).split("\t", -1);
            maRows.put(f[0], f);
        }
        failures += check("kept exactly {refid_1, refid_2}: " + maRows.keySet(),
            maRows.keySet().equals(new LinkedHashSet<>(Arrays.asList("refid_1", "refid_2"))));

        String[] kept = maRows.get("refid_1");
        if (kept != null) {
            failures += checkEq("refid_1 A1 (GWAS's own EA, unchanged)", kept[1], "A");
            failures += checkEq("refid_1 A2 (GWAS's own NEA, unchanged)", kept[2], "G");
            failures += checkClose("refid_1 freq = refFreq(A1) as-is = 0.40 (EA matches ref A1)",
                Double.parseDouble(kept[3]), 0.40, 1e-6);
            failures += checkClose("refid_1 b = 0.30 (unchanged)", Double.parseDouble(kept[4]), 0.30, 1e-6);
            failures += checkEq("refid_1 N = 500 (per-SNP N used, Fix 3)", kept[7], "500");
        }
        String[] flipped = maRows.get("refid_2");
        if (flipped != null) {
            failures += checkEq("refid_2 A1 (GWAS's own EA \"T\", NOT remapped to ref's A1 \"C\")", flipped[1], "T");
            failures += checkClose("refid_2 freq FLIPPED = 1 - refFreq(0.70) = 0.30 (EA matches ref A2)",
                Double.parseDouble(flipped[3]), 0.30, 1e-6);
            failures += checkClose("refid_2 b = -0.25 (beta sign unchanged — only freq is oriented)",
                Double.parseDouble(flipped[4]), -0.25, 1e-6);
            failures += checkEq("refid_2 N = 1000 (row's own N was NA -> falls back to global sampleN)",
                flipped[7], "1000");
        }

        List<String> audit = Files.readAllLines(new File(runDir, "ma_audit.tsv").toPath());
        String auditText = String.join("\n", audit);
        failures += check("audit records snpNoRef as dropped_no_ref_id (Fix 2: never falls back "
            + "to using the GWAS's own id \"snpNoRef\" in place of a ref panel id): " + auditText,
            auditText.contains("snpNoRef\t\tdropped_no_ref_id"));
        failures += check("audit records snpNoFreq as dropped_no_freq",
            auditText.contains("snpNoFreq\trefid_5\tdropped_no_freq"));
        failures += check("audit records snpOrientBad as dropped_orientation (EA=G matches neither "
            + "refid_4's A1=A nor A2=C)",
            auditText.contains("snpOrientBad\trefid_4\tdropped_orientation"));
        boolean nanBetaMentioned = auditText.contains("snpNaNBeta");
        boolean garbageBetaMentioned = auditText.contains("snpGarbageBeta");
        failures += check("bad-stat rows (NaN-string beta, non-numeric beta) are silently skipped "
            + "(dropped_stats has no audit trail column of its own in the current code — verified "
            + "instead via the .ma file NOT containing either row) rather than crashing",
            !maRows.containsKey("snpNaNBeta") && !maRows.containsKey("snpGarbageBeta")
                && written == 2);

        File rScript = new File(runDir, "run_cojo.R");
        failures += check("run_cojo.R generated", rScript.exists());

        deleteRecursive(dir.toFile());
        return failures;
    }

    /** Fix 4: an existing ../ld/consistency_summary.json must have its verdict/flagged-count substituted verbatim into run_cojo.R. */
    private static int testPrepareRunWiresLdConsistencyIntoRScript() throws Exception {
        Path dir = Files.createTempDirectory("cojo-ldconsist-test");
        File analysisRoot = dir.toFile();
        File harmonizedDir = new File(analysisRoot, "harmonized");
        File ldDir = new File(analysisRoot, "ld");
        File matchedDir = new File(analysisRoot, "matched");
        File runDir = new File(analysisRoot, "run");
        harmonizedDir.mkdirs();
        ldDir.mkdirs();
        matchedDir.mkdirs();
        File fakeGcta = new File(analysisRoot, "fake_gcta64.exe");
        fakeGcta.createNewFile();

        try (PrintWriter pw = new PrintWriter(new FileWriter(new File(matchedDir, "matched_ref.bim")))) {
            pw.println("1\trefid_1\t0\t1000000\tA\tG");
        }
        try (PrintWriter pw = new PrintWriter(new FileWriter(new File(matchedDir, "ref_freq.tsv")))) {
            pw.println("snp_id\tA1\tfreq");
            pw.println("refid_1\tA\t0.40");
        }
        try (PrintWriter pw = new PrintWriter(new FileWriter(new File(harmonizedDir, "harmonized_gwas.tsv")))) {
            pw.println("snp_id\tchr\tpos\tea\tnea\tpvalue\tbeta\tse\tn");
            pw.println("snpKeep\t1\t1000000\tA\tG\t0.001\t0.30\t0.05\t500");
        }
        // consistency_summary.json fixture: verdict="high_warn", flagged_snps=3 (hand-picked
        // values, chosen only to be distinguishable from any default/zero value in the code).
        try (PrintWriter pw = new PrintWriter(new FileWriter(new File(ldDir, "consistency_summary.json")))) {
            pw.println("{\"verdict\":\"high_warn\",\"flagged_snps\":3,\"other_field\":\"ignored\"}");
        }

        CojoAdapter.prepareRun(harmonizedDir, matchedDir, runDir, 1000, 5e-8, 0.9, fakeGcta.getAbsolutePath());

        String rContent = new String(Files.readAllBytes(new File(runDir, "run_cojo.R").toPath()), "UTF-8");
        int failures = 0;
        failures += check("run_cojo.R substitutes ld_consistency_verdict <- 'high_warn' from consistency_summary.json",
            rContent.contains("ld_consistency_verdict <- 'high_warn'"));
        failures += check("run_cojo.R substitutes ld_consistency_flagged <- 3 from consistency_summary.json",
            rContent.contains("ld_consistency_flagged <- 3"));

        deleteRecursive(dir.toFile());
        return failures;
    }

    /**
     * Fix 5: the freq-discrepancy-warning parsing logic (comparing GCTA's .jma.cojo freq vs
     * freq_geno columns) is R code generated by this Java method — there is no separate Java
     * parser to unit-test directly (it only ever runs after a real GCTA invocation produces a
     * real .jma.cojo file). This test verifies the generated R script actually contains that
     * Fix-5 logic, rather than actually running GCTA.
     */
    private static int testPrepareRunWiresFreqDiscrepancyParsingIntoRScript() throws Exception {
        Path dir = Files.createTempDirectory("cojo-freqdisc-test");
        File harmonizedDir = new File(dir.toFile(), "harmonized");
        File matchedDir = new File(dir.toFile(), "matched");
        File runDir = new File(dir.toFile(), "run");
        harmonizedDir.mkdirs();
        matchedDir.mkdirs();
        File fakeGcta = new File(dir.toFile(), "fake_gcta64.exe");
        fakeGcta.createNewFile();

        try (PrintWriter pw = new PrintWriter(new FileWriter(new File(matchedDir, "matched_ref.bim")))) {
            pw.println("1\trefid_1\t0\t1000000\tA\tG");
        }
        try (PrintWriter pw = new PrintWriter(new FileWriter(new File(matchedDir, "ref_freq.tsv")))) {
            pw.println("snp_id\tA1\tfreq");
            pw.println("refid_1\tA\t0.40");
        }
        try (PrintWriter pw = new PrintWriter(new FileWriter(new File(harmonizedDir, "harmonized_gwas.tsv")))) {
            pw.println("snp_id\tchr\tpos\tea\tnea\tpvalue\tbeta\tse\tn");
            pw.println("snpKeep\t1\t1000000\tA\tG\t0.001\t0.30\t0.05\t500");
        }

        CojoAdapter.prepareRun(harmonizedDir, matchedDir, runDir, 1000, 5e-8, 0.9, fakeGcta.getAbsolutePath());

        String rContent = new String(Files.readAllBytes(new File(runDir, "run_cojo.R").toPath()), "UTF-8");
        int failures = 0;
        failures += check("run_cojo.R computes jma$freq_diff <- abs(jma$freq - jma$freq_geno) (Fix 5)",
            rContent.contains("jma$freq_diff <- abs(jma$freq - jma$freq_geno)"));
        failures += check("run_cojo.R flags rows with freq_diff > 0.15 as freq_discrepancies (Fix 5)",
            rContent.contains("jma[jma$freq_diff > 0.15, , drop=FALSE]"));
        failures += check("run_cojo.R parses GCTA's own log for 'large difference of allele frequency' warnings",
            rContent.contains("large difference of allele frequency"));

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
