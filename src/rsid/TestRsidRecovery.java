package rsid;

import htsjdk.samtools.util.*;
import htsjdk.tribble.index.*;
import htsjdk.tribble.index.tabix.*;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.variantcontext.writer.*;
import htsjdk.variant.vcf.*;
import java.io.*;
import java.nio.file.*;
import java.util.*;

/**
 * Test harness for rsID recovery module.
 * 1. Generates synthetic bgzipped+tabix-indexed test VCF using HTSJDK
 * 2. Runs assertion-based tests covering every matching branch
 * 3. Optionally runs against a real locus from the actual dbSNP database
 *
 * Run:  java -cp "bin;lib/*" rsid.TestRsidRecovery
 */
public class TestRsidRecovery {

    static int passed = 0, failed = 0;
    static Path testDir;

    public static void main(String[] args) throws Exception {
        System.out.println("=== rsID Recovery Test Harness ===\n");

        testDir = Files.createTempDirectory("rsid_test_");

        try {
            // Phase 1: Synthetic tests
            createSyntheticTestData();
            runSyntheticTests();

            // Phase 2: Real data test (single small locus)
            runRealDataTest();

        } finally {
            deleteDir(testDir.toFile());
        }

        System.out.printf("%n=== Results: %d passed, %d failed ===%n", passed, failed);
        if (failed > 0) System.exit(1);
    }

    // ── Synthetic test data ──────────────────────────────────────

    static void createSyntheticTestData() throws Exception {
        System.out.println("── Creating synthetic test VCF ──");

        Path vcfPath = testDir.resolve("1.vcf.gz");

        // Define test VCF header with minimal sequence dictionary
        Set<VCFHeaderLine> headerLines = new LinkedHashSet<>();
        headerLines.add(new VCFInfoHeaderLine("dbSNPBuildID", 1,
            VCFHeaderLineType.Integer, "dbSNP Build ID"));
        headerLines.add(new VCFContigHeaderLine(
            Collections.singletonMap("ID", "1"), 0));
        VCFHeader header = new VCFHeader(headerLines, Collections.emptyList());

        htsjdk.samtools.SAMSequenceDictionary dict = new htsjdk.samtools.SAMSequenceDictionary(
            Collections.singletonList(new htsjdk.samtools.SAMSequenceRecord("1", 300000000)));
        header.setSequenceDictionary(dict);

        // Write bgzipped VCF using HTSJDK
        VariantContextWriterBuilder builder = new VariantContextWriterBuilder()
            .setOutputPath(vcfPath)
            .setReferenceDictionary(dict)
            .unsetOption(Options.INDEX_ON_THE_FLY)
            .setOption(Options.ALLOW_MISSING_FIELDS_IN_HEADER);

        try (VariantContextWriter writer = builder.build()) {
            writer.writeHeader(header);

            // Test case 1: Forward match (nea==ref, ea==alt)
            writer.add(makeVC("1", 1000, "rs1001", "A", Arrays.asList("G"), 155));

            // Test case 2: Reverse match (ea==ref, nea==alt)
            writer.add(makeVC("1", 2000, "rs2001", "C", Arrays.asList("T"), 155));

            // Test case 3: Multi-allelic (ALT=A,G,T) — GWAS alt is one of several
            writer.add(makeVC("1", 3000, "rs3001", "C", Arrays.asList("A", "G", "T"), 155));

            // Test case 4: Two rsID versions at same position — pick newest
            writer.add(makeVC("1", 4000, "rs4001", "A", Arrays.asList("G"), 130));
            writer.add(makeVC("1", 4000, "rs4002", "A", Arrays.asList("G"), 155));

            // Test case 5: Position present but alleles don't match
            writer.add(makeVC("1", 5000, "rs5001", "A", Arrays.asList("C"), 155));

            // Test case 6: No record at position 6000 (gap — tests no_pos)
        }

        // Create tabix index using HTSJDK
        TabixIndex index = IndexFactory.createTabixIndex(
            vcfPath, new VCFCodec(), TabixFormat.VCF, null);
        index.write(vcfPath.resolveSibling(vcfPath.getFileName() + ".tbi"));

        System.out.printf("  Created %s + .tbi%n", vcfPath.getFileName());
    }

    static VariantContext makeVC(String chr, int pos, String id, String ref, List<String> alts, int build) {
        Allele refAllele = Allele.create(ref, true);
        List<Allele> alleles = new ArrayList<>();
        alleles.add(refAllele);
        for (String a : alts) alleles.add(Allele.create(a, false));

        VariantContextBuilder vcb = new VariantContextBuilder()
            .chr(chr).start(pos).stop(pos + ref.length() - 1)
            .alleles(alleles).id(id);
        vcb.attribute("dbSNPBuildID", build);
        return vcb.make();
    }

    // ── Synthetic assertions ─────────────────────────────────────

    static void runSyntheticTests() throws Exception {
        System.out.println("\n── Running synthetic tests ──");

        RsidRecovery recovery = new RsidRecovery(testDir.toString(), "hg19");
        Map<Integer, List<DbSnpRecord>> posIndex = recovery.buildPositionIndex("1", 1, 10000);

        check("Index loaded positions", posIndex.size() >= 5);

        // Test 1: Forward match — GWAS nea==ref(A), ea==alt(G)
        MatchResult r1 = RsidMatcher.match("1", 1000, "G", "A", posIndex);
        check("Forward match: assigned rs1001", "rs1001".equals(r1.assignedRsid));
        check("Forward match: reason", "matched_forward".equals(r1.matchReason));

        // Test 2: Reverse match — GWAS ea==ref(C), nea==alt(T)
        MatchResult r2 = RsidMatcher.match("1", 2000, "C", "T", posIndex);
        check("Reverse match: assigned rs2001", "rs2001".equals(r2.assignedRsid));
        check("Reverse match: reason", "matched_reverse".equals(r2.matchReason));

        // Test 3: Multi-allelic — GWAS ea=G, nea=C (ref=C, alts=A,G,T → forward)
        MatchResult r3 = RsidMatcher.match("1", 3000, "G", "C", posIndex);
        check("Multi-allelic: assigned rs3001", "rs3001".equals(r3.assignedRsid));
        check("Multi-allelic: reason forward", "matched_forward".equals(r3.matchReason));

        // Test 3b: Multi-allelic with different alt
        MatchResult r3b = RsidMatcher.match("1", 3000, "T", "C", posIndex);
        check("Multi-allelic alt T: assigned rs3001", "rs3001".equals(r3b.assignedRsid));

        // Test 4: Two versions at same position — newest (build 155) wins
        MatchResult r4 = RsidMatcher.match("1", 4000, "G", "A", posIndex);
        check("Version selection: assigned rs4002 (newest)", "rs4002".equals(r4.assignedRsid));
        check("Version selection: 2 candidates", r4.nCandidatesAtPos == 2);

        // Test 5: Position present but alleles mismatch (ref=A, alt=C; GWAS has T,G)
        MatchResult r5 = RsidMatcher.match("1", 5000, "T", "G", posIndex);
        check("Allele mismatch: no assigned rsid", r5.assignedRsid == null);
        check("Allele mismatch: reason", "pos_only_allele_mismatch".equals(r5.matchReason));
        check("Allele mismatch: pos_only_rsid set", "rs5001".equals(r5.posOnlyRsid));

        // Test 6: Position absent
        MatchResult r6 = RsidMatcher.match("1", 6000, "A", "G", posIndex);
        check("No position: no assigned rsid", r6.assignedRsid == null);
        check("No position: reason", "no_pos".equals(r6.matchReason));
        check("No position: 0 candidates", r6.nCandidatesAtPos == 0);

        // Test 7: Case-insensitive alleles (lowercase GWAS vs uppercase VCF)
        MatchResult r7 = RsidMatcher.match("1", 1000, "g", "a", posIndex);
        check("Case insensitive: assigned rs1001", "rs1001".equals(r7.assignedRsid));
        check("Case insensitive: reason forward", "matched_forward".equals(r7.matchReason));

        // Print CSV output
        System.out.println("\n── Synthetic CSV output ──");
        System.out.println(MatchResult.csvHeader());
        for (MatchResult r : Arrays.asList(r1, r2, r3, r3b, r4, r5, r6, r7)) {
            System.out.println(r.toCsv());
        }

        RsidRecovery.printSummary(Arrays.asList(r1, r2, r3, r3b, r4, r5, r6, r7));
    }

    // ── Real data test ───────────────────────────────────────────

    static void runRealDataTest() throws Exception {
        System.out.println("\n── Real data test (single locus) ──");

        String snpdbDir = "C:\\Users\\alsammana\\Documents\\DataResources\\SNPdb\\hg19";
        String gwasFile = "input/Asianq2EURFIN.tsv";
        String lociFile = "input/loci.txt";

        if (!Files.exists(Paths.get(snpdbDir, "1.vcf.gz"))) {
            System.out.println("  [SKIP] dbSNP files not found at: " + snpdbDir);
            return;
        }
        if (!Files.exists(Paths.get(gwasFile))) {
            System.out.println("  [SKIP] GWAS file not found: " + gwasFile);
            return;
        }
        if (!Files.exists(Paths.get(lociFile))) {
            System.out.println("  [SKIP] Loci file not found: " + lociFile);
            return;
        }

        // Read first locus only for quick test
        String firstChr;
        int firstStart, firstEnd;
        try (BufferedReader br = new BufferedReader(new FileReader(lociFile))) {
            String header = br.readLine();
            String[] hcols = header.trim().split("\t");
            int iChr = -1, iStart = -1, iEnd = -1;
            for (int i = 0; i < hcols.length; i++) {
                String h = hcols[i].trim().toLowerCase();
                if (h.equals("meta_chr"))   iChr = i;
                if (h.equals("meta_start")) iStart = i;
                if (h.equals("meta_end"))   iEnd = i;
            }
            String line = br.readLine();
            String[] f = line.trim().split("\t");
            firstChr   = f[iChr].trim();
            firstStart = Integer.parseInt(f[iStart].trim());
            firstEnd   = Integer.parseInt(f[iEnd].trim());
        }

        System.out.printf("  Locus 1: chr%s:%,d-%,d%n", firstChr, firstStart, firstEnd);

        // Build position index from real dbSNP
        RsidRecovery recovery = new RsidRecovery(snpdbDir, "hg19");
        long t0 = System.currentTimeMillis();
        Map<Integer, List<DbSnpRecord>> posIndex =
            recovery.buildPositionIndex(firstChr, firstStart, firstEnd);
        long indexMs = System.currentTimeMillis() - t0;
        System.out.printf("  dbSNP index: %,d positions loaded in %,d ms%n",
            posIndex.size(), indexMs);

        // Read GWAS SNPs in this locus
        List<String[]> gwasSnps = new ArrayList<>();
        try (BufferedReader br = new BufferedReader(new FileReader(gwasFile), 1024 * 1024)) {
            String header = br.readLine();
            String[] hcols = header.trim().split("\t");
            int iChr2 = colIdx(hcols, "chrom"), iPos = colIdx(hcols, "pos");
            int iEa = colIdx(hcols, "ea"), iNea = colIdx(hcols, "nea");

            String line;
            while ((line = br.readLine()) != null) {
                if (line.isEmpty()) continue;
                String[] f = line.split("\t");
                String chr = f[iChr2].trim();
                if (!chr.equals(firstChr)) continue;
                int pos;
                try { pos = Integer.parseInt(f[iPos].trim()); }
                catch (NumberFormatException e) { continue; }
                if (pos < firstStart || pos > firstEnd) continue;
                gwasSnps.add(new String[]{chr, f[iPos].trim(), f[iEa].trim(), f[iNea].trim()});
            }
        }
        System.out.printf("  GWAS SNPs in locus: %,d%n", gwasSnps.size());

        // Match
        t0 = System.currentTimeMillis();
        List<MatchResult> results = RsidRecovery.matchDirect(firstChr,
            gwasSnps.toArray(new String[0][]), posIndex);
        long matchMs = System.currentTimeMillis() - t0;

        System.out.printf("  Matching completed in %,d ms%n", matchMs);

        RsidRecovery.printSummary(results);

        // Write CSV
        Path csvOut = Paths.get("output/rsid_recovery_locus1.csv");
        try (PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(csvOut.toString())))) {
            pw.println(MatchResult.csvHeader());
            for (MatchResult r : results) pw.println(r.toCsv());
        }
        System.out.printf("  CSV written to %s%n", csvOut);

        int matched = (int) results.stream()
            .filter(r -> r.matchReason.startsWith("matched_")).count();
        check("Real data: recovery rate > 0%", matched > 0);
        check("Real data: processed all SNPs", results.size() == gwasSnps.size());
    }

    // ── Helpers ──────────────────────────────────────────────────

    static void check(String name, boolean condition) {
        if (condition) { System.out.println("  PASS  " + name); passed++; }
        else           { System.out.println("  FAIL  " + name); failed++; }
    }

    static int colIdx(String[] cols, String name) {
        for (int i = 0; i < cols.length; i++)
            if (cols[i].trim().equalsIgnoreCase(name)) return i;
        return -1;
    }

    static void deleteDir(File dir) {
        if (dir == null || !dir.exists()) return;
        File[] files = dir.listFiles();
        if (files != null) for (File f : files) {
            if (f.isDirectory()) deleteDir(f);
            else f.delete();
        }
        dir.delete();
    }
}
