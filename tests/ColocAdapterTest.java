import java.io.*;
import java.nio.file.*;
import java.util.*;

/**
 * Standalone regression test (no external test framework — this project has none) for
 * ColocAdapter's Java-side trait1/trait2 matching, allele harmonization (sign-flip on swapped
 * alleles, drop on mismatched alleles), the four input-validation guards, and — per this task's
 * explicit brief — whether the MHC/extended-HLA region is excluded from colocalization.
 *
 * *** MHC EXCLUSION FINDING (see testMhcRegionIsNotExcluded below) ***
 * DECISIONS_PHASE2.md's phase-recap line records "exclude MHC from fine-mapping/coloc" as a
 * precaution that was "actively applied" this session. That is NOT what the actual code does.
 * BaseStepPipeline.overlapsMhc(Locus, String genomeBuild) exists and IS wired into the base
 * pipeline step — but only to print an informational WARNING line into that step's log; it
 * explicitly does NOT block or exclude anything (see BaseStepPipeline.java's own comment: "We
 * don't block analysis ... just surface a clear warning"). ColocAdapter.prepareRun() goes
 * further than "warn but don't exclude": it never references BaseStepPipeline.overlapsMhc, MHC,
 * or the locus's own chr/position at all — grep confirms the only occurrence of the substring
 * "locus." in the whole file is inside an unrelated error message string ("...in this locus."),
 * not a field access. The `Locus locus` parameter ColocAdapter.prepareRun() takes is otherwise
 * unused. So for coloc specifically there is no warning AND no exclusion — an MHC-region SNP is
 * silently colocalized exactly like any other. Per this task's instruction ("if it doesn't
 * [exclude MHC], that's a real bug to flag, not fix without asking"), this is reported as a real,
 * unfixed gap rather than patched here — colocalization statistical-logic changes are easy to
 * get subtly wrong and were explicitly called out as needing sign-off rather than a same-day fix.
 * testMhcRegionIsNotExcluded() below documents this current (undesired) behavior with a
 * hand-built MHC-region fixture rather than silently assuming it.
 */
public class ColocAdapterTest {

    public static void main(String[] args) throws Exception {
        int failures = 0;
        failures += testPrepareRunMatchesAndHarmonizesAlleles();
        failures += testPrepareRunThrowsWhenTrait2FileMissing();
        failures += testPrepareRunThrowsWhenTrait2PathDoesNotExist();
        failures += testPrepareRunThrowsWhenSampleNZero();
        failures += testPrepareRunThrowsWhenNoOverlap();
        failures += testMhcRegionIsNotExcluded();

        if (failures == 0) {
            System.out.println("PASS: all ColocAdapter tests passed");
        } else {
            System.out.println("FAIL: " + failures + " test(s) failed");
            System.exit(1);
        }
    }

    /**
     * harmonized_gwas.tsv (trait 1, header per LocusGwasExtractor's 14-column format):
     *   snpA chr1:1000000 EA=A NEA=G beta=0.40  se=0.05 maf=0.30
     *   snpB chr1:2000000 EA=C NEA=T beta=-0.25 se=0.04 maf=0.20
     *   snpE chr1:5000000 EA=A NEA=G beta=0.10  se=0.02 maf=0.40
     *
     * trait2.tsv (space-delimited, header "CHR POS A1 A2 BETA SE EAF"):
     *   1 1000000 A G 0.35  0.06 0.28  -> alleles ALIGNED with snpA (A/G == A/G) -> no flip
     *   1 2000000 T C -0.30 0.05 0.75  -> alleles SWAPPED vs snpB (T/C == nea/ea reversed of C/T)
     *                                     -> beta2 must be sign-flipped: -(-0.30) = +0.30
     *   1 5000000 A C 0.20  0.03 0.50  -> allele MISMATCH vs snpE (A/C matches neither A/G nor
     *                                     G/A) -> row must be dropped, not guessed
     *   1 9999999 A G 0.10  0.02 0.40  -> no trait-1 SNP at this position at all -> dropped
     *
     * Expected coloc_input.tsv: exactly 2 rows (snpA, snpB); snpE and the unmatched position
     * must not appear.
     */
    private static int testPrepareRunMatchesAndHarmonizesAlleles() throws Exception {
        Path dir = Files.createTempDirectory("coloc-happy-test");
        File harmonizedDir = new File(dir.toFile(), "harmonized");
        File runDir = new File(dir.toFile(), "run");
        harmonizedDir.mkdirs();

        try (PrintWriter pw = new PrintWriter(new FileWriter(new File(harmonizedDir, "harmonized_gwas.tsv")))) {
            pw.println("snp_id\tchr\tpos\tea\tnea\tpvalue\tbeta\tse\tor\tn\tmaf\tinfo\trsid\tvarid");
            pw.println("snpA\t1\t1000000\tA\tG\t0.001\t0.40\t0.05\t1.5\t1000\t0.30\t0.9\trsA\tvA");
            pw.println("snpB\t1\t2000000\tC\tT\t0.002\t-0.25\t0.04\t0.9\t1200\t0.20\t0.9\trsB\tvB");
            pw.println("snpE\t1\t5000000\tA\tG\t0.010\t0.10\t0.02\t1.1\t900\t0.40\t0.9\trsE\tvE");
        }

        File trait2File = new File(dir.toFile(), "trait2.tsv");
        try (PrintWriter pw = new PrintWriter(new FileWriter(trait2File))) {
            pw.println("CHR POS A1 A2 BETA SE EAF");
            pw.println("1 1000000 A G 0.35 0.06 0.28");
            pw.println("1 2000000 T C -0.30 0.05 0.75");
            pw.println("1 5000000 A C 0.20 0.03 0.50");
            pw.println("1 9999999 A G 0.10 0.02 0.40");
        }

        Locus locus = new Locus(0, "1", 900000, 6000000, 0);
        Config config = new Config();
        config.sampleN = 1000;

        ColocAdapter.prepareRun(harmonizedDir, runDir, locus, config, trait2File.getAbsolutePath(),
            "quant", /*trait2N*/2000, /*trait2NCases*/0, /*p1*/1e-4, /*p2*/1e-4, /*p12*/1e-5,
            /*analysisRoot*/null, /*restrictToFinemapped*/false, /*pipThreshold*/0.5);

        int failures = 0;
        File inputTsv = new File(runDir, "coloc_input.tsv");
        failures += check("coloc_input.tsv exists", inputTsv.exists());
        List<String> lines = Files.readAllLines(inputTsv.toPath());
        failures += checkEq("coloc_input.tsv line count (header + 2 matched rows)", lines.size(), 3);
        failures += checkEq("header", lines.get(0), "snp_id\tchr\tpos\tbeta1\tse1\tmaf1\tbeta2\tse2\tmaf2");

        Map<String, String[]> rows = new LinkedHashMap<>();
        for (int i = 1; i < lines.size(); i++) {
            String[] f = lines.get(i).split("\t", -1);
            rows.put(f[0], f);
        }
        failures += check("kept exactly {snpA, snpB}, snpE (mismatch) and the unmatched pos excluded: "
            + rows.keySet(), rows.keySet().equals(new LinkedHashSet<>(Arrays.asList("snpA", "snpB"))));

        String[] rA = rows.get("snpA");
        if (rA != null) {
            failures += checkClose("snpA beta1 = 0.40 (trait 1, unchanged)", Double.parseDouble(rA[3]), 0.40, 1e-9);
            failures += checkClose("snpA se1 = 0.05", Double.parseDouble(rA[4]), 0.05, 1e-9);
            failures += checkClose("snpA maf1 = 0.30 (trait 1's own MAF column)", Double.parseDouble(rA[5]), 0.30, 1e-9);
            failures += checkClose("snpA beta2 = 0.35 (alleles A/G aligned with trait 1 -> no flip)",
                Double.parseDouble(rA[6]), 0.35, 1e-9);
            failures += checkClose("snpA se2 = 0.06", Double.parseDouble(rA[7]), 0.06, 1e-9);
            failures += checkClose("snpA maf2 = 0.28 (trait 2's own EAF column)", Double.parseDouble(rA[8]), 0.28, 1e-9);
        }
        String[] rB = rows.get("snpB");
        if (rB != null) {
            failures += checkClose("snpB beta1 = -0.25 (trait 1, unchanged)", Double.parseDouble(rB[3]), -0.25, 1e-9);
            failures += checkClose(
                "snpB beta2 SIGN-FLIPPED: trait2's A1/A2=T/C is trait1's NEA/EA reversed (C/T) "
                + "-> beta2 = -(-0.30) = 0.30", Double.parseDouble(rB[6]), 0.30, 1e-9);
        }

        deleteRecursive(dir.toFile());
        return failures;
    }

    private static int testPrepareRunThrowsWhenTrait2FileMissing() throws Exception {
        Path dir = Files.createTempDirectory("coloc-notrait2-test");
        File harmonizedDir = new File(dir.toFile(), "harmonized");
        File runDir = new File(dir.toFile(), "run");
        Locus locus = new Locus(0, "1", 900000, 6000000, 0);
        Config config = new Config();
        config.sampleN = 1000;

        boolean threw = false;
        String message = null;
        try {
            ColocAdapter.prepareRun(harmonizedDir, runDir, locus, config, null, "quant", 2000, 0,
                1e-4, 1e-4, 1e-5, null, false, 0.5);
        } catch (IOException e) {
            threw = true;
            message = e.getMessage();
        }

        int failures = 0;
        failures += check("prepareRun() throws IOException when trait2File is null", threw);
        if (threw) {
            failures += check("message names trait 2 requirement: " + message,
                message != null && message.contains("Trait 2 summary statistics file is required"));
        }
        deleteRecursive(dir.toFile());
        return failures;
    }

    private static int testPrepareRunThrowsWhenTrait2PathDoesNotExist() throws Exception {
        Path dir = Files.createTempDirectory("coloc-badtrait2path-test");
        File harmonizedDir = new File(dir.toFile(), "harmonized");
        File runDir = new File(dir.toFile(), "run");
        Locus locus = new Locus(0, "1", 900000, 6000000, 0);
        Config config = new Config();
        config.sampleN = 1000;
        String missingPath = new File(dir.toFile(), "does_not_exist.tsv").getAbsolutePath();

        boolean threw = false;
        String message = null;
        try {
            ColocAdapter.prepareRun(harmonizedDir, runDir, locus, config, missingPath, "quant", 2000, 0,
                1e-4, 1e-4, 1e-5, null, false, 0.5);
        } catch (IOException e) {
            threw = true;
            message = e.getMessage();
        }

        int failures = 0;
        failures += check("prepareRun() throws IOException when trait2File path does not exist on disk", threw);
        if (threw) {
            failures += check("message names the missing file: " + message,
                message != null && message.contains("Trait 2 file not found"));
        }
        deleteRecursive(dir.toFile());
        return failures;
    }

    private static int testPrepareRunThrowsWhenSampleNZero() throws Exception {
        Path dir = Files.createTempDirectory("coloc-noN-test");
        File harmonizedDir = new File(dir.toFile(), "harmonized");
        File runDir = new File(dir.toFile(), "run");
        harmonizedDir.mkdirs();
        File trait2File = new File(dir.toFile(), "trait2.tsv");
        try (PrintWriter pw = new PrintWriter(new FileWriter(trait2File))) {
            pw.println("CHR POS A1 A2 BETA SE EAF");
            pw.println("1 1000000 A G 0.1 0.02 0.3");
        }
        try (PrintWriter pw = new PrintWriter(new FileWriter(new File(harmonizedDir, "harmonized_gwas.tsv")))) {
            pw.println("snp_id\tchr\tpos\tea\tnea\tpvalue\tbeta\tse\tor\tn\tmaf\tinfo\trsid\tvarid");
            pw.println("snpA\t1\t1000000\tA\tG\t0.001\t0.4\t0.05\t1.5\t1000\t0.3\t0.9\trsA\tvA");
        }
        Locus locus = new Locus(0, "1", 900000, 6000000, 0);
        Config config = new Config(); // sampleN defaults to 0

        boolean threw = false;
        String message = null;
        try {
            ColocAdapter.prepareRun(harmonizedDir, runDir, locus, config, trait2File.getAbsolutePath(),
                "quant", 2000, 0, 1e-4, 1e-4, 1e-5, null, false, 0.5);
        } catch (IOException e) {
            threw = true;
            message = e.getMessage();
        }

        int failures = 0;
        failures += check("prepareRun() throws IOException when config.sampleN<=0", threw);
        if (threw) {
            failures += check("message names coloc trait-1 N requirement: " + message,
                message != null && message.contains("Sample size (N) is required for coloc"));
        }
        deleteRecursive(dir.toFile());
        return failures;
    }

    private static int testPrepareRunThrowsWhenNoOverlap() throws Exception {
        Path dir = Files.createTempDirectory("coloc-nooverlap-test");
        File harmonizedDir = new File(dir.toFile(), "harmonized");
        File runDir = new File(dir.toFile(), "run");
        harmonizedDir.mkdirs();
        try (PrintWriter pw = new PrintWriter(new FileWriter(new File(harmonizedDir, "harmonized_gwas.tsv")))) {
            pw.println("snp_id\tchr\tpos\tea\tnea\tpvalue\tbeta\tse\tor\tn\tmaf\tinfo\trsid\tvarid");
            pw.println("snpA\t1\t1000000\tA\tG\t0.001\t0.4\t0.05\t1.5\t1000\t0.3\t0.9\trsA\tvA");
        }
        File trait2File = new File(dir.toFile(), "trait2.tsv");
        try (PrintWriter pw = new PrintWriter(new FileWriter(trait2File))) {
            pw.println("CHR POS A1 A2 BETA SE EAF");
            pw.println("1 7777777 A G 0.1 0.02 0.3"); // different position entirely
        }
        Locus locus = new Locus(0, "1", 900000, 8000000, 0);
        Config config = new Config();
        config.sampleN = 1000;

        boolean threw = false;
        String message = null;
        try {
            ColocAdapter.prepareRun(harmonizedDir, runDir, locus, config, trait2File.getAbsolutePath(),
                "quant", 2000, 0, 1e-4, 1e-4, 1e-5, null, false, 0.5);
        } catch (IOException e) {
            threw = true;
            message = e.getMessage();
        }

        int failures = 0;
        failures += check("prepareRun() throws IOException when trait1/trait2 share no SNPs", threw);
        if (threw) {
            failures += check("message explains no overlap: " + message,
                message != null && message.contains("No overlapping, allele-consistent SNPs"));
        }
        deleteRecursive(dir.toFile());
        return failures;
    }

    /**
     * *** DOCUMENTS A REAL, UNFIXED GAP — see the class-level comment above for full context. ***
     * A SNP squarely inside the GRCh37 MHC region (chr6:28,477,897-33,448,354, per
     * BaseStepPipeline's own MHC_GRCH37_START/END constants) is placed at chr6:30,000,000 — the
     * midpoint of that range, unambiguously inside it. This test first sanity-checks the fixture
     * itself really does overlap the MHC using BaseStepPipeline.overlapsMhc() (the project's own
     * existing MHC-detection logic), then runs ColocAdapter.prepareRun() and asserts — as CURRENT,
     * NOT DESIRED, behavior — that the MHC SNP is NOT excluded from coloc_input.tsv. If a future
     * fix adds MHC exclusion to ColocAdapter, this specific assertion will need to be inverted;
     * until then it documents the gap so it cannot silently regress further or be assumed fixed.
     */
    private static int testMhcRegionIsNotExcluded() throws Exception {
        Path dir = Files.createTempDirectory("coloc-mhc-test");
        File harmonizedDir = new File(dir.toFile(), "harmonized");
        File runDir = new File(dir.toFile(), "run");
        harmonizedDir.mkdirs();

        long mhcPos = 30_000_000L; // inside [28,477,897 .. 33,448,354] for GRCh37
        try (PrintWriter pw = new PrintWriter(new FileWriter(new File(harmonizedDir, "harmonized_gwas.tsv")))) {
            pw.println("snp_id\tchr\tpos\tea\tnea\tpvalue\tbeta\tse\tor\tn\tmaf\tinfo\trsid\tvarid");
            pw.println("snpMHC\t6\t" + mhcPos + "\tA\tG\t0.001\t0.30\t0.05\t1.2\t1000\t0.25\t0.9\trsMHC\tvMHC");
        }
        File trait2File = new File(dir.toFile(), "trait2.tsv");
        try (PrintWriter pw = new PrintWriter(new FileWriter(trait2File))) {
            pw.println("CHR POS A1 A2 BETA SE EAF");
            pw.println("6 " + mhcPos + " A G 0.28 0.05 0.22");
        }

        Locus locus = new Locus(0, "6", mhcPos - 1000, mhcPos + 1000, 0);
        Config config = new Config();
        config.sampleN = 1000;
        config.genomeBuild = "GRCh37";

        int failures = 0;
        failures += check("fixture sanity check: BaseStepPipeline.overlapsMhc() confirms this locus "
            + "really is inside the GRCh37 MHC region (otherwise this test would prove nothing)",
            BaseStepPipeline.overlapsMhc(locus, config.genomeBuild));

        ColocAdapter.prepareRun(harmonizedDir, runDir, locus, config, trait2File.getAbsolutePath(),
            "quant", 2000, 0, 1e-4, 1e-4, 1e-5, null, false, 0.5);

        List<String> lines = Files.readAllLines(new File(runDir, "coloc_input.tsv").toPath());
        boolean mhcSnpPresent = lines.stream().anyMatch(l -> l.startsWith("snpMHC\t"));
        failures += check("GAP (flagged, not fixed): the MHC-region SNP 'snpMHC' IS still present "
            + "in coloc_input.tsv — ColocAdapter.prepareRun() performs no MHC exclusion or even a "
            + "warning (unlike the base pipeline step, which at least warns via overlapsMhc()), "
            + "despite DECISIONS_PHASE2.md recording MHC exclusion from coloc as an applied "
            + "precaution. See this file's class-level comment for the full write-up.",
            mhcSnpPresent);

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
