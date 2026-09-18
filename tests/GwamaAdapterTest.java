import java.io.*;
import java.nio.file.*;
import java.util.*;

/**
 * Standalone regression test (no external test framework — this project has none) for
 * GwamaAdapter — the smallest of the six analysis adapters (56 lines). It has no allele
 * orientation logic (it is a documented pass-through single-cohort GWAMA input writer, per its
 * own class comment), so this test focuses on: column extraction/reordering into GWAMA's
 * MARKERNAME/EA/NEA/BETA/SE/N/EAF/STRAND layout, NA-beta/se row dropping, the N/EAF fallback
 * defaults ("0" / "NA") when those optional columns are themselves NA or absent, the
 * too-short-row skip (f.length<8), the two failure paths (missing harmonized_gwas.tsv, and zero
 * usable rows), and the gwama_input.txt cohort-list file.
 */
public class GwamaAdapterTest {

    public static void main(String[] args) throws Exception {
        int failures = 0;
        failures += testPrepareRunWritesCohortAndListFiles();
        failures += testPrepareRunThrowsWhenGwasFileMissing();
        failures += testPrepareRunThrowsWhenNoUsableSnps();

        if (failures == 0) {
            System.out.println("PASS: all GwamaAdapter tests passed");
        } else {
            System.out.println("FAIL: " + failures + " test(s) failed");
            System.exit(1);
        }
    }

    /**
     * harmonized_gwas.tsv fixture using the project's real 14-column header (see
     * LocusGwasExtractor.java:133 — snp_id chr pos ea nea pvalue beta se or n maf info rsid varid,
     * indices 0-13):
     *
     *   Row 1 (snpKeep): beta=0.5 se=0.1 n=5000 maf=0.3 -> all present, kept as-is.
     *   Row 2 (snpBadBeta): beta=NA -> dropped (GwamaAdapter checks beta.equals("NA")).
     *   Row 3 (snpFallback): n=NA, maf=NA -> N column falls back to the literal string "0"
     *          (code: `(f.length > 9 && !f[9].equals("NA") && !f[9].isEmpty()) ? f[9] : "0"`),
     *          EAF column falls back to "NA" (code: `(f.length > 10 && !f[10].equals("NA")) ?
     *          f[10] : "NA"`).
     *   Row 4 (snpShort): only 5 tab-separated fields -> f.length<8 -> skipped entirely, must
     *          not appear anywhere in cohort1.txt and must not throw.
     *
     * Expected cohort1.txt (header + 2 data rows, in gwas-file iteration order):
     *   MARKERNAME  EA  NEA  BETA  SE    N     EAF  STRAND
     *   snpKeep     A   G    0.5   0.1   5000  0.3  +
     *   snpFallback C   T    -0.2  0.05  0     NA   +
     * count returned by the write loop = 2 (verified indirectly via row count in cohort1.txt,
     * since prepareRun itself returns void — its only externally observable "count" is the
     * printed summary line and the file it wrote).
     */
    private static int testPrepareRunWritesCohortAndListFiles() throws Exception {
        Path dir = Files.createTempDirectory("gwama-happy-test");
        File harmonizedDir = new File(dir.toFile(), "harmonized");
        File runDir = new File(dir.toFile(), "run");
        harmonizedDir.mkdirs();

        try (PrintWriter pw = new PrintWriter(new FileWriter(new File(harmonizedDir, "harmonized_gwas.tsv")))) {
            pw.println("snp_id\tchr\tpos\tea\tnea\tpvalue\tbeta\tse\tor\tn\tmaf\tinfo\trsid\tvarid");
            pw.println("snpKeep\t1\t1000000\tA\tG\t0.001\t0.5\t0.1\t1.65\t5000\t0.3\t0.9\trsA\tvA");
            pw.println("snpBadBeta\t1\t2000000\tA\tG\t0.002\tNA\t0.05\tNA\t4000\t0.4\t0.9\trsB\tvB");
            pw.println("snpFallback\t1\t3000000\tC\tT\t0.003\t-0.2\t0.05\t0.82\tNA\tNA\t0.95\trsC\tvC");
            pw.println("snpShort\t1\t4000000\tA");
        }

        GwamaAdapter.prepareRun(harmonizedDir, runDir);

        int failures = 0;
        File cohortFile = new File(runDir, "cohort1.txt");
        failures += check("cohort1.txt exists", cohortFile.exists());
        List<String> lines = Files.readAllLines(cohortFile.toPath());
        failures += checkEq("cohort1.txt line count (header + 2 kept rows)", lines.size(), 3);
        failures += checkEq("cohort1.txt header", lines.get(0),
            "MARKERNAME\tEA\tNEA\tBETA\tSE\tN\tEAF\tSTRAND");
        failures += checkEq("row 1 (snpKeep, all fields present)", lines.get(1),
            "snpKeep\tA\tG\t0.5\t0.1\t5000\t0.3\t+");
        failures += checkEq("row 2 (snpFallback: N NA->\"0\", EAF NA->\"NA\")", lines.get(2),
            "snpFallback\tC\tT\t-0.2\t0.05\t0\tNA\t+");
        boolean mentionsBadBeta = lines.stream().anyMatch(l -> l.contains("snpBadBeta"));
        failures += check("snpBadBeta (beta=NA) does not appear in cohort1.txt", !mentionsBadBeta);
        boolean mentionsShort = lines.stream().anyMatch(l -> l.contains("snpShort"));
        failures += check("snpShort (too few columns) does not appear in cohort1.txt", !mentionsShort);

        File listFile = new File(runDir, "gwama_input.txt");
        failures += check("gwama_input.txt exists", listFile.exists());
        List<String> listLines = Files.readAllLines(listFile.toPath());
        failures += checkEq("gwama_input.txt has exactly one cohort path", listLines.size(), 1);
        failures += check("gwama_input.txt's single line is cohort1.txt's absolute path: " + listLines.get(0),
            listLines.get(0).equals(cohortFile.getAbsolutePath().replace("\\", "/")));

        deleteRecursive(dir.toFile());
        return failures;
    }

    private static int testPrepareRunThrowsWhenGwasFileMissing() throws Exception {
        Path dir = Files.createTempDirectory("gwama-nogwas-test");
        File harmonizedDir = new File(dir.toFile(), "harmonized"); // deliberately not created
        File runDir = new File(dir.toFile(), "run");

        boolean threw = false;
        String message = null;
        try {
            GwamaAdapter.prepareRun(harmonizedDir, runDir);
        } catch (IOException e) {
            threw = true;
            message = e.getMessage();
        }

        int failures = 0;
        failures += check("prepareRun() throws IOException when harmonized_gwas.tsv is missing", threw);
        if (threw) {
            failures += check("message mentions harmonized_gwas.tsv: " + message,
                message != null && message.contains("harmonized_gwas.tsv"));
        }
        deleteRecursive(dir.toFile());
        return failures;
    }

    /** Every row has beta=NA -> zero usable SNPs -> prepareRun() must refuse rather than write an empty/useless cohort file. */
    private static int testPrepareRunThrowsWhenNoUsableSnps() throws Exception {
        Path dir = Files.createTempDirectory("gwama-empty-test");
        File harmonizedDir = new File(dir.toFile(), "harmonized");
        File runDir = new File(dir.toFile(), "run");
        harmonizedDir.mkdirs();

        try (PrintWriter pw = new PrintWriter(new FileWriter(new File(harmonizedDir, "harmonized_gwas.tsv")))) {
            pw.println("snp_id\tchr\tpos\tea\tnea\tpvalue\tbeta\tse\tor\tn\tmaf\tinfo\trsid\tvarid");
            pw.println("snpA\t1\t1000000\tA\tG\t0.001\tNA\t0.1\tNA\t5000\t0.3\t0.9\trsA\tvA");
        }

        boolean threw = false;
        String message = null;
        try {
            GwamaAdapter.prepareRun(harmonizedDir, runDir);
        } catch (IOException e) {
            threw = true;
            message = e.getMessage();
        }

        int failures = 0;
        failures += check("prepareRun() throws IOException when no row has usable beta/se", threw);
        if (threw) {
            failures += check("message explains no usable SNPs: " + message,
                message != null && message.toLowerCase().contains("no usable"));
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

    private static int check(String label, boolean cond) {
        if (!cond) {
            System.out.println("FAIL: " + label);
            return 1;
        }
        System.out.println("PASS: " + label);
        return 0;
    }
}
