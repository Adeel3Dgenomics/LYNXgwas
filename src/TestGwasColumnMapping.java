import java.io.*;
import java.nio.file.*;
import java.util.*;

/**
 * Tests for Step 4: GwasParser configurable column mapping.
 * Run:  java -cp bin TestGwasColumnMapping
 */
public class TestGwasColumnMapping {

    static int passed = 0, failed = 0;

    public static void main(String[] args) throws Exception {
        System.out.println("=== Step 4 GWAS Column Mapping Tests ===\n");

        testStandardMapping();
        testDifferentHeaders();
        testMissingRequiredColumn();
        testOptionalColumnConfiguredButMissing();
        testOptionalColumnsCarriedThrough();
        testNoOptionalFieldsInvented();

        System.out.printf("%n=== Results: %d passed, %d failed ===%n", passed, failed);
        if (failed > 0) System.exit(1);
    }

    // TEST 4.1: GWAS file with headers CHR,BP,P,A1,A2,SNP,OR and matching config.
    //           Assert: parses correctly, OR column carried through; no optional fields silently invented.
    static void testStandardMapping() throws Exception {
        Path tmp = Files.createTempDirectory("lynx_gwas41_");
        Path gwasFile = tmp.resolve("gwas.tsv");
        Files.writeString(gwasFile,
            "CHR\tBP\tP\tA1\tA2\tSNP\tOR\n" +
            "1\t1000\t5e-8\tA\tG\trs123\t1.25\n" +
            "1\t2000\t0.01\tC\tT\trs456\t0.85\n" +
            "1\t3000\t0.5\tG\tA\trs789\t1.01\n");

        Config config = new Config();
        config.gwasFile  = gwasFile.toString();
        config.colChr    = "CHR";
        config.colPos    = "BP";
        config.colPvalue = "P";
        config.colEa     = "A1";
        config.colNea    = "A2";
        config.colRsid   = "SNP";
        config.colOr     = "OR";
        config.colVarid  = "";
        config.maxSnpsPerLocus = 5000;

        Locus locus = new Locus(1, "1", 500, 3500, 500);
        List<Locus> loci = new ArrayList<>(Collections.singletonList(locus));

        GwasParser.parse(config, loci);

        check("4.1 Parsed 3 SNPs", locus.snps.size() == 3);
        check("4.1 First SNP id=rs123", "rs123".equals(locus.snps.get(0).id));
        check("4.1 First SNP pos=1000", locus.snps.get(0).pos == 1000);
        check("4.1 First SNP ea=A",     "A".equals(locus.snps.get(0).ea));
        check("4.1 First SNP nea=G",    "G".equals(locus.snps.get(0).nea));

        // OR column carried through
        check("4.1 First SNP OR=1.25",  locus.snps.get(0).oddsRatio == 1.25);
        check("4.1 Second SNP OR=0.85", locus.snps.get(1).oddsRatio == 0.85);

        // No optional fields invented (beta, se, etc. should all be NaN)
        check("4.1 beta not invented", Double.isNaN(locus.snps.get(0).beta));
        check("4.1 se not invented",   Double.isNaN(locus.snps.get(0).se));
        check("4.1 n not invented",    Double.isNaN(locus.snps.get(0).sampleN));
        check("4.1 maf not invented",  Double.isNaN(locus.snps.get(0).maf));
        check("4.1 info not invented", Double.isNaN(locus.snps.get(0).infoScore));

        deleteDir(tmp.toFile());
    }

    // TEST 4.2: A second GWAS file uses different headers (Chromosome, Position, Pval, ...)
    //           with its own project's config.properties mapping.
    static void testDifferentHeaders() throws Exception {
        Path tmp = Files.createTempDirectory("lynx_gwas42_");

        // File A: standard headers
        Path gwasA = tmp.resolve("gwas_a.tsv");
        Files.writeString(gwasA,
            "CHR\tBP\tP\tA1\tA2\n" +
            "1\t1500\t1e-6\tA\tG\n");

        // File B: completely different headers
        Path gwasB = tmp.resolve("gwas_b.tsv");
        Files.writeString(gwasB,
            "Chromosome\tPosition\tPval\tAllele1\tAllele2\n" +
            "1\t1500\t2e-6\tC\tT\n");

        // Config A
        Config configA = new Config();
        configA.gwasFile = gwasA.toString();
        configA.colChr = "CHR"; configA.colPos = "BP"; configA.colPvalue = "P";
        configA.colEa = "A1"; configA.colNea = "A2";
        configA.colVarid = ""; configA.colRsid = "";
        configA.maxSnpsPerLocus = 5000;

        // Config B
        Config configB = new Config();
        configB.gwasFile = gwasB.toString();
        configB.colChr = "Chromosome"; configB.colPos = "Position"; configB.colPvalue = "Pval";
        configB.colEa = "Allele1"; configB.colNea = "Allele2";
        configB.colVarid = ""; configB.colRsid = "";
        configB.maxSnpsPerLocus = 5000;

        Locus locusA = new Locus(1, "1", 1000, 2000, 1000);
        List<Locus> lociA = new ArrayList<>(Collections.singletonList(locusA));
        GwasParser.parse(configA, lociA);

        Locus locusB = new Locus(1, "1", 1000, 2000, 1000);
        List<Locus> lociB = new ArrayList<>(Collections.singletonList(locusB));
        GwasParser.parse(configB, lociB);

        check("4.2 Project A parsed 1 SNP",    locusA.snps.size() == 1);
        check("4.2 Project B parsed 1 SNP",    locusB.snps.size() == 1);
        check("4.2 Project A p=1e-6",          locusA.snps.get(0).pvalue == 1e-6);
        check("4.2 Project B p=2e-6",          locusB.snps.get(0).pvalue == 2e-6);
        check("4.2 Project A ea=A",            "A".equals(locusA.snps.get(0).ea));
        check("4.2 Project B ea=C",            "C".equals(locusB.snps.get(0).ea));

        // Verify no hardcoded column names — parse with config B's different names
        check("4.2 Not hardcoded to CHR",
            locusB.snps.get(0).chr.equals("1")); // parsed via "Chromosome" not "CHR"

        deleteDir(tmp.toFile());
    }

    // TEST 4.3: config.properties is missing col.pvalue (uses default "p"),
    //           but the file header has no "p" column.
    //           Assert: clear error naming "col.pvalue" and "p-value".
    static void testMissingRequiredColumn() throws Exception {
        Path tmp = Files.createTempDirectory("lynx_gwas43_");
        Path gwasFile = tmp.resolve("gwas.tsv");
        Files.writeString(gwasFile,
            "CHR\tBP\tPVALUE\tA1\tA2\n" +
            "1\t1000\t0.05\tA\tG\n");

        // Config uses default colPvalue="p" but file header has "PVALUE" not "p"
        Config config = new Config();
        config.gwasFile  = gwasFile.toString();
        config.colChr    = "CHR";
        config.colPos    = "BP";
        // colPvalue remains "p" (default) — won't match "PVALUE"
        config.colEa     = "A1";
        config.colNea    = "A2";
        config.colVarid  = "";
        config.colRsid   = "";
        config.maxSnpsPerLocus = 5000;

        Locus locus = new Locus(1, "1", 500, 1500, 500);
        List<Locus> loci = new ArrayList<>(Collections.singletonList(locus));

        String errorMsg = null;
        try {
            GwasParser.parse(config, loci);
        } catch (IOException e) {
            errorMsg = e.getMessage();
        }

        check("4.3 Throws an error", errorMsg != null);
        check("4.3 Error mentions 'col.pvalue'",
            errorMsg != null && errorMsg.contains("col.pvalue"));
        check("4.3 Error mentions 'p-value'",
            errorMsg != null && errorMsg.contains("p-value"));
        check("4.3 Error mentions configured value 'p'",
            errorMsg != null && errorMsg.contains("='p'"));
        check("4.3 Error lists available columns",
            errorMsg != null && errorMsg.contains("PVALUE"));

        deleteDir(tmp.toFile());
    }

    // TEST 4.4: config maps col.or=ODDS_RATIO but file header has no such column.
    //           Assert: clear error naming the missing column, not a silent skip.
    static void testOptionalColumnConfiguredButMissing() throws Exception {
        Path tmp = Files.createTempDirectory("lynx_gwas44_");
        Path gwasFile = tmp.resolve("gwas.tsv");
        Files.writeString(gwasFile,
            "CHR\tBP\tP\tA1\tA2\tOR\n" +
            "1\t1000\t0.05\tA\tG\t1.2\n");

        Config config = new Config();
        config.gwasFile  = gwasFile.toString();
        config.colChr    = "CHR";
        config.colPos    = "BP";
        config.colPvalue = "P";
        config.colEa     = "A1";
        config.colNea    = "A2";
        config.colVarid  = "";
        config.colRsid   = "";
        config.colOr     = "ODDS_RATIO"; // configured but doesn't exist in header
        config.maxSnpsPerLocus = 5000;

        Locus locus = new Locus(1, "1", 500, 1500, 500);
        List<Locus> loci = new ArrayList<>(Collections.singletonList(locus));

        String errorMsg = null;
        try {
            GwasParser.parse(config, loci);
        } catch (IOException e) {
            errorMsg = e.getMessage();
        }

        check("4.4 Throws an error", errorMsg != null);
        check("4.4 Error mentions 'col.or'",
            errorMsg != null && errorMsg.contains("col.or"));
        check("4.4 Error mentions 'ODDS_RATIO'",
            errorMsg != null && errorMsg.contains("ODDS_RATIO"));
        check("4.4 Error mentions 'odds ratio'",
            errorMsg != null && errorMsg.contains("odds ratio"));
        check("4.4 Error lists available columns",
            errorMsg != null && errorMsg.contains("OR"));

        deleteDir(tmp.toFile());
    }

    // Full optional columns: beta, OR, SE, N, MAF, INFO all present and carried through
    static void testOptionalColumnsCarriedThrough() throws Exception {
        Path tmp = Files.createTempDirectory("lynx_gwas4opt_");
        Path gwasFile = tmp.resolve("gwas.tsv");
        Files.writeString(gwasFile,
            "CHR\tBP\tP\tA1\tA2\tBETA\tOR\tSE\tN\tMAF\tINFO\n" +
            "1\t1000\t1e-8\tA\tG\t0.35\t1.42\t0.06\t50000\t0.12\t0.98\n" +
            "1\t2000\t0.01\tC\tT\tNA\t.\t\t.\tnan\t0.75\n");

        Config config = new Config();
        config.gwasFile  = gwasFile.toString();
        config.colChr = "CHR"; config.colPos = "BP"; config.colPvalue = "P";
        config.colEa = "A1"; config.colNea = "A2";
        config.colBeta = "BETA"; config.colOr = "OR"; config.colSe = "SE";
        config.colN = "N"; config.colMaf = "MAF"; config.colInfo = "INFO";
        config.colVarid = ""; config.colRsid = "";
        config.maxSnpsPerLocus = 5000;

        Locus locus = new Locus(1, "1", 500, 2500, 500);
        List<Locus> loci = new ArrayList<>(Collections.singletonList(locus));
        GwasParser.parse(config, loci);

        check("opt: 2 SNPs parsed", locus.snps.size() == 2);

        Snp s1 = locus.snps.get(0);
        check("opt: s1 beta=0.35",   s1.beta == 0.35);
        check("opt: s1 or=1.42",     s1.oddsRatio == 1.42);
        check("opt: s1 se=0.06",     s1.se == 0.06);
        check("opt: s1 n=50000",     s1.sampleN == 50000.0);
        check("opt: s1 maf=0.12",    s1.maf == 0.12);
        check("opt: s1 info=0.98",   s1.infoScore == 0.98);

        // Second SNP has NA/./empty/nan values → all should be NaN
        Snp s2 = locus.snps.get(1);
        check("opt: s2 beta=NaN (NA)",   Double.isNaN(s2.beta));
        check("opt: s2 or=NaN (.)",      Double.isNaN(s2.oddsRatio));
        check("opt: s2 se=NaN (empty)",  Double.isNaN(s2.se));
        check("opt: s2 n=NaN (.)",       Double.isNaN(s2.sampleN));
        check("opt: s2 maf=NaN (nan)",   Double.isNaN(s2.maf));
        check("opt: s2 info=0.75",       s2.infoScore == 0.75);

        deleteDir(tmp.toFile());
    }

    // When optional columns are NOT configured, no fields should be invented
    static void testNoOptionalFieldsInvented() throws Exception {
        Path tmp = Files.createTempDirectory("lynx_gwas4none_");
        Path gwasFile = tmp.resolve("gwas.tsv");
        Files.writeString(gwasFile,
            "CHR\tBP\tP\tA1\tA2\tBETA\tOR\n" +
            "1\t1000\t1e-8\tA\tG\t0.35\t1.42\n");

        // Config does NOT map beta or OR
        Config config = new Config();
        config.gwasFile  = gwasFile.toString();
        config.colChr = "CHR"; config.colPos = "BP"; config.colPvalue = "P";
        config.colEa = "A1"; config.colNea = "A2";
        config.colVarid = ""; config.colRsid = "";
        // colBeta, colOr, etc. all remain "" (default)
        config.maxSnpsPerLocus = 5000;

        Locus locus = new Locus(1, "1", 500, 1500, 500);
        List<Locus> loci = new ArrayList<>(Collections.singletonList(locus));
        GwasParser.parse(config, loci);

        Snp snp = locus.snps.get(0);
        check("no-invent: beta=NaN", Double.isNaN(snp.beta));
        check("no-invent: or=NaN",   Double.isNaN(snp.oddsRatio));
        check("no-invent: se=NaN",   Double.isNaN(snp.se));
        check("no-invent: n=NaN",    Double.isNaN(snp.sampleN));
        check("no-invent: maf=NaN",  Double.isNaN(snp.maf));
        check("no-invent: info=NaN", Double.isNaN(snp.infoScore));

        deleteDir(tmp.toFile());
    }

    // ── Helpers ──────────────────────────────────────────────────────────

    static void check(String name, boolean condition) {
        if (condition) { System.out.println("  PASS  " + name); passed++; }
        else           { System.out.println("  FAIL  " + name); failed++; }
    }

    static void deleteDir(File dir) {
        File[] files = dir.listFiles();
        if (files != null) for (File f : files) {
            if (f.isDirectory()) deleteDir(f);
            else f.delete();
        }
        dir.delete();
    }
}
