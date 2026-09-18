import java.io.*;
import java.nio.file.*;
import java.util.*;

/**
 * Standalone regression test (no external test framework — this project has none) for
 * MagmaAdapter: .genes.out / .gsa.out parsing, findMagmaBinary() resolution failure, and the
 * ref-panel-ID-matching robustness (MagmaAdapter.buildPvalFile) that plays the same correctness
 * role here as CojoAdapter's allele-orientation check does for beta/freq. MAGMA itself has no
 * orientation logic to extract: its gene-based test consumes unsigned P-values, which have no
 * allele direction to flip — see the class-level comment on MagmaAdapter.
 *
 * Test classes here have no package statement, matching every source file under src/ (flat
 * default package), so package-private static helpers such as parseGenesOut/parseGsaOut/
 * buildPvalFile are callable directly (see LdCalculatorTest calling LdCalculator's
 * package-private computeWindowKb() the same way) — no reflection needed.
 */
public class MagmaAdapterTest {

    public static void main(String[] args) throws Exception {
        int failures = 0;
        failures += testParseGenesOut();
        failures += testParseGsaOut();
        failures += testFindMagmaBinaryNotFound();
        failures += testBuildPvalFileDropsSnpWithNoRefPanelMatch();

        if (failures == 0) {
            System.out.println("PASS: all MagmaAdapter tests passed");
        } else {
            System.out.println("FAIL: " + failures + " test(s) failed");
            System.exit(1);
        }
    }

    /**
     * .genes.out fixture with two hand-written rows:
     *   Row 1: ZSTAT=3.5 -> its P (2.32629e-04) is the well-known standard-normal upper-tail
     *          reference value P(Z>3.5) = 1 - Phi(3.5) ~ 0.0002326 (a commonly tabulated
     *          constant, independently checkable against any standard normal table), included
     *          here as a sanity check that the fixture row itself looks like a realistic MAGMA
     *          gene-test output row, not because the parser recomputes it.
     *   Row 2: ZSTAT=0.0 -> P=0.5 exactly, by symmetry of the standard normal distribution
     *          around zero (P(Z>0) = 0.5), a trivial hand-computable value with no dependency
     *          on any table lookup.
     * The test itself checks that parseGenesOut() extracts every column verbatim from the file
     * — i.e. it is a parser-fidelity test, not a re-derivation of MAGMA's internal statistics.
     */
    private static int testParseGenesOut() throws Exception {
        Path tmp = Files.createTempFile("magma-genes", ".genes.out");
        try (PrintWriter pw = new PrintWriter(new FileWriter(tmp.toFile()))) {
            pw.println("# TOTAL_GENES = 2");
            pw.println("GENE       CHR START   STOP    NSNPS NPARAM N     ZSTAT   P");
            pw.println("ENSG00001  1   1000000 1010000 25    5      50000 3.5000  2.32629e-04");
            pw.println("ENSG00002  1   1050000 1060000 10    3      50000 0.0000  5.00000e-01");
        }

        List<MagmaAdapter.GeneResult> genes = MagmaAdapter.parseGenesOut(tmp.toFile());
        Files.deleteIfExists(tmp);

        int failures = 0;
        if (genes.size() != 2) {
            System.out.println("FAIL: expected 2 parsed gene rows, got " + genes.size());
            return 1;
        }

        MagmaAdapter.GeneResult g1 = genes.get(0);
        failures += checkEq("gene1 GENE", g1.gene, "ENSG00001");
        failures += checkEq("gene1 CHR", g1.chr, "1");
        failures += checkEq("gene1 START", g1.start, 1000000L);
        failures += checkEq("gene1 STOP", g1.stop, 1010000L);
        failures += checkEq("gene1 NSNPS", g1.nsnps, 25);
        failures += checkEq("gene1 NPARAM", g1.nparam, 5);
        failures += checkEq("gene1 N", g1.n, 50000);
        failures += checkClose("gene1 ZSTAT", g1.zstat, 3.5, 1e-9);
        failures += checkClose("gene1 P (== standard-normal upper-tail reference for Z=3.5)", g1.p, 2.32629e-04, 1e-9);

        MagmaAdapter.GeneResult g2 = genes.get(1);
        failures += checkEq("gene2 GENE", g2.gene, "ENSG00002");
        failures += checkClose("gene2 ZSTAT", g2.zstat, 0.0, 1e-9);
        failures += checkClose("gene2 P (== 0.5 by normal-distribution symmetry at Z=0)", g2.p, 0.5, 1e-9);

        return failures;
    }

    /**
     * .gsa.out fixture with two hand-written rows for a competitive gene-set test:
     * VARIABLE TYPE NGENES BETA BETA_STD SE P. Values are arbitrary but plausible; this test
     * checks parser fidelity, exactly as testParseGenesOut() does.
     */
    private static int testParseGsaOut() throws Exception {
        Path tmp = Files.createTempFile("magma-gsa", ".gsa.out");
        try (PrintWriter pw = new PrintWriter(new FileWriter(tmp.toFile()))) {
            pw.println("VARIABLE   TYPE  NGENES BETA     BETA_STD SE       P");
            pw.println("PATHWAY_A  SET   42     0.150000 0.080000 0.050000 6.20000e-03");
            pw.println("PATHWAY_B  SET   17     0.010000 0.005000 0.060000 8.60000e-01");
        }

        List<MagmaAdapter.GeneSetResult> sets = MagmaAdapter.parseGsaOut(tmp.toFile());
        Files.deleteIfExists(tmp);

        int failures = 0;
        if (sets.size() != 2) {
            System.out.println("FAIL: expected 2 parsed gene-set rows, got " + sets.size());
            return 1;
        }
        MagmaAdapter.GeneSetResult s1 = sets.get(0);
        failures += checkEq("set1 VARIABLE", s1.setId, "PATHWAY_A");
        failures += checkEq("set1 NGENES", s1.nGenes, 42);
        failures += checkClose("set1 BETA", s1.beta, 0.15, 1e-9);
        failures += checkClose("set1 BETA_STD", s1.betaStd, 0.08, 1e-9);
        failures += checkClose("set1 SE", s1.se, 0.05, 1e-9);
        failures += checkClose("set1 P", s1.p, 6.20000e-03, 1e-9);

        MagmaAdapter.GeneSetResult s2 = sets.get(1);
        failures += checkEq("set2 VARIABLE", s2.setId, "PATHWAY_B");
        failures += checkEq("set2 NGENES", s2.nGenes, 17);

        return failures;
    }

    /**
     * findMagmaBinary() must throw a clear IOException naming the checked bin/ path when magma
     * is not present there. appRoot is a fresh empty temp directory, so its bin/ subfolder never
     * exists; the PATH fallback is exercised against the machine's real PATH, which is not
     * expected to contain a "magma"/"magma.exe" binary in this environment.
     */
    private static int testFindMagmaBinaryNotFound() throws Exception {
        Path tmpRoot = Files.createTempDirectory("magma-approot");
        boolean threw = false;
        String message = null;
        try {
            MagmaAdapter.findMagmaBinary(tmpRoot.toFile());
        } catch (IOException e) {
            threw = true;
            message = e.getMessage();
        }
        Files.deleteIfExists(tmpRoot);

        int failures = 0;
        failures += check("findMagmaBinary throws IOException when not found", threw);
        if (threw) {
            boolean mentionsBin = message != null && message.contains("bin");
            boolean mentionsMagma = message != null && message.toLowerCase().contains("magma");
            failures += check("IOException message names the checked bin/ location: " + message, mentionsBin);
            failures += check("IOException message names the magma binary: " + message, mentionsMagma);
        } else {
            System.out.println("NOTE: this only fails to demonstrate not-found behavior if this "
                + "machine happens to have a 'magma'/'magma.exe' binary on its real PATH.");
        }
        return failures;
    }

    /**
     * MAGMA's --pval SNP column must match the reference panel's own BIM IDs used in --snp-loc
     * (same correctness reason CojoAdapter never falls back to the GWAS's own SNP id — Fix 2).
     * Builds a harmonized_gwas.tsv with two SNPs, only one of which has a matching chr:pos entry
     * in matched_ref.bim, and confirms buildPvalFile() writes exactly the matched one (under its
     * ref-panel ID, not its GWAS id) while dropping and counting the unmatched one.
     */
    private static int testBuildPvalFileDropsSnpWithNoRefPanelMatch() throws Exception {
        Path dir = Files.createTempDirectory("magma-refid-test");
        File harmonizedDir = new File(dir.toFile(), "harmonized");
        File matchedDir = new File(dir.toFile(), "matched");
        File runDir = new File(dir.toFile(), "run");
        harmonizedDir.mkdirs();
        matchedDir.mkdirs();

        // Ref panel BIM: only "refid_1" at chr1:1000000 is present; the GWAS's chr1:2000000 SNP
        // (gwas_snp_2) has no corresponding ref panel entry at all.
        try (PrintWriter pw = new PrintWriter(new FileWriter(new File(matchedDir, "matched_ref.bim")))) {
            pw.println("1\trefid_1\t0\t1000000\tA\tG");
        }

        // harmonized_gwas.tsv: snp_id chr pos ea nea pvalue beta se (CojoAdapter/FinemapAdapter's convention)
        try (PrintWriter pw = new PrintWriter(new FileWriter(new File(harmonizedDir, "harmonized_gwas.tsv")))) {
            pw.println("snp_id\tchr\tpos\tea\tnea\tpvalue\tbeta\tse");
            pw.println("gwas_snp_1\t1\t1000000\tA\tG\t0.001\t0.2\t0.05"); // matches refid_1 -> kept
            pw.println("gwas_snp_2\t1\t2000000\tA\tG\t0.002\t0.3\t0.06"); // no ref panel match -> dropped
        }

        MagmaAdapter.PvalBuildResult result = MagmaAdapter.buildPvalFile(harmonizedDir, matchedDir, runDir, 10000);

        int failures = 0;
        failures += checkEq("written count (1 ref-matched SNP kept)", result.written, 1);
        failures += checkEq("dropped-no-ref-id count (1 unmatched SNP dropped)", result.droppedNoRefId, 1);
        failures += checkEq("dropped-stats count (no bad-stat rows in this fixture)", result.droppedStats, 0);
        failures += check("hasPerSnpN is false (fixture has no N column, global N=10000 was supplied)",
            !result.hasPerSnpN);

        List<String> lines = Files.readAllLines(result.pvalFile.toPath());
        failures += checkEq(".pval file line count (header + 1 kept SNP)", lines.size(), 2);
        boolean rowUsesRefId = lines.size() > 1 && lines.get(1).startsWith("refid_1\t");
        failures += check(".pval file's kept row is written under the REF PANEL id (refid_1), "
            + "not the GWAS's own id (gwas_snp_1) — matches CojoAdapter's Fix 2 rule: " + (lines.size() > 1 ? lines.get(1) : "<none>"),
            rowUsesRefId);
        boolean mentionsUnmatched = lines.stream().anyMatch(l -> l.contains("gwas_snp_2"));
        failures += check("the unmatched SNP (gwas_snp_2) does not appear anywhere in the .pval file",
            !mentionsUnmatched);

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
}
