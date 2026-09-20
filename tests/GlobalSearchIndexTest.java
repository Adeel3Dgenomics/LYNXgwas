import java.io.*;
import java.nio.charset.StandardCharsets;
import java.nio.file.*;
import java.util.*;

/**
 * Builds 3 tiny synthetic fake "projects" (hand-written config.properties + data/manifest.json +
 * data/locus_N.json, matching the real schema read by GlobalSearchIndex) in a temp directory, then
 * exercises every case called out in the task:
 *   - a gene-contains match (SNP significant AND inside a gene's own [start,end])
 *   - a gene-nearest-fallback match (SNP significant but outside every gene's boundary in that locus)
 *   - an rsID exact match
 *   - a same-build position match
 *   - a GRCh37->GRCh38 position match via real liftover (skipped, with a clear printed message, if
 *     resources/hg19ToHg38.over.chain isn't present on this machine)
 *   - the unsupported GRCh38->GRCh37 direction returning the explicit "not compared" flag
 *   - threshold filtering (a SNP below genome-wide significance must never appear in results)
 *
 * Fixture layout:
 *   proj_a (GRCh37, effect.type blank): locus 1 has genes GENE_IN[1000,2000] and GENE_OUT[5000,6000],
 *     nearest_genes=[NEARBY1,NEARBY2], and 3 SNPs:
 *       rs1001 chr1:1500 p=1e-9   -> inside GENE_IN            -> "contains" match
 *       rs1002 chr1:9000 p=2e-8   -> outside every gene        -> "nearest" match (NEARBY1, NEARBY2)
 *       rs1003 chr1:1600 p=0.5    -> not significant           -> must never appear
 *   proj_b (GRCh37, effect.type=beta): locus 1 has one real, independently-documented anchor SNP,
 *     rs429358 at its well-known hg19 position chr19:45,411,941 (p=1e-10) — the same anchor
 *     GenomeLiftover.java's own javadoc validates against, lifting to exactly chr19:44,908,684.
 *   proj_c (GRCh38, effect.type=OR): locus 1 has rs777 at chr2:100000 (p=5e-9) — used to prove a
 *     GRCh37 position query never silently compares against this project's un-liftable GRCh38 data.
 */
public class GlobalSearchIndexTest {

    static int failures = 0;

    public static void main(String[] args) throws Exception {
        Path root = Files.createTempDirectory("gsi_test");
        try {
            writeProjectA(root);
            writeProjectB(root);
            writeProjectC(root);

            GlobalSearchIndex idx = GlobalSearchIndex.build(root.toFile(), 5e-8);

            check("index summary: 3 projects indexed", idx.projectCount, 3);
            check("index summary: 2 loci indexed (proj_a, proj_b — proj_c also has 1 = 3 total)", idx.lociCount, 3);
            // significant SNPs only: proj_a has 2 (rs1001, rs1002 — rs1003 excluded), proj_b has 1, proj_c has 1
            check("index summary: only significant SNPs counted", idx.snpCount, 4);
            // Regression: summaryJson() once emitted a double comma before "built_at" (invalid JSON) —
            // caught by manual curl verification against real project data, fixed by passing first=true
            // to the kv() call since a comma had already been appended just before it.
            check("summaryJson() has no double commas (valid JSON)", idx.summaryJson().contains(",,"), false);

            // ── Gene-contains match ──────────────────────────────────────
            GlobalSearchIndex.Result geneIn = idx.searchGene("GENE_IN");
            check("GENE_IN: exactly 1 contains-match", geneIn.hits.size(), 1);
            if (!geneIn.hits.isEmpty()) {
                GlobalSearchIndex.Hit h = geneIn.hits.get(0);
                check("GENE_IN hit: match_type=contains", h.matchType, "contains");
                check("GENE_IN hit: snp_id=rs1001", h.snpId, "rs1001");
                check("GENE_IN hit: project_id=proj_a", h.projectId, "proj_a");
                check("GENE_IN hit: project_name resolved via ProjectMetadata", h.projectName, "Project A Test");
                check("GENE_IN hit: effect_type blank handled gracefully (not null, not guessed)", h.effectType, "");
                // Must not crash / must produce valid-looking JSON even with a blank effect type.
                String json = h.toJson();
                check("GENE_IN hit JSON contains effect_type key", json.contains("\"effect_type\":\"\""), true);
            }

            // GENE_OUT contains no significant SNP at all -> zero hits, never a false contains-match.
            GlobalSearchIndex.Result geneOut = idx.searchGene("GENE_OUT");
            check("GENE_OUT: zero matches (no significant SNP falls inside it)", geneOut.hits.size(), 0);

            // ── Gene-nearest-fallback match ──────────────────────────────
            GlobalSearchIndex.Result nearby1 = idx.searchGene("NEARBY1");
            check("NEARBY1: exactly 1 nearest-fallback match", nearby1.hits.size(), 1);
            if (!nearby1.hits.isEmpty()) {
                check("NEARBY1 hit: match_type=nearest", nearby1.hits.get(0).matchType, "nearest");
                check("NEARBY1 hit: snp_id=rs1002", nearby1.hits.get(0).snpId, "rs1002");
            }
            GlobalSearchIndex.Result nearby2 = idx.searchGene("NEARBY2");
            check("NEARBY2: exactly 1 nearest-fallback match (same SNP, second nearest gene)", nearby2.hits.size(), 1);

            // gene search is case-insensitive
            check("gene search is case-insensitive", idx.searchGene("gene_in").hits.size(), 1);

            // ── rsID exact match ──────────────────────────────────────────
            GlobalSearchIndex.Result rsHit = idx.searchSnp("rs1001");
            check("rs1001: exactly 1 hit", rsHit.hits.size(), 1);
            if (!rsHit.hits.isEmpty()) {
                check("rs1001 hit: chr=1", rsHit.hits.get(0).chr, "1");
                check("rs1001 hit: pos=1500", rsHit.hits.get(0).pos, 1500L);
            }

            // ── Threshold filtering ─────────────────────────────────────
            check("rs1003 (p=0.5, below genome-wide significance) never appears", idx.searchSnp("rs1003").hits.size(), 0);
            check("rs1003's position also never matches (excluded at index-build time)",
                idx.searchPosition("1", 1600, "GRCh37").hits.size(), 0);

            // ── Same-build position match ─────────────────────────────────
            GlobalSearchIndex.Result samePos = idx.searchPosition("1", 1500, "GRCh37");
            check("same-build GRCh37 position match: 1 hit", samePos.hits.size(), 1);
            if (!samePos.hits.isEmpty()) check("same-build position hit is rs1001", samePos.hits.get(0).snpId, "rs1001");

            // ── GRCh37->GRCh38 position match via real liftover ────────────
            File chainFile = new File("resources/hg19ToHg38.over.chain");
            if (!chainFile.isFile()) {
                System.out.println("SKIP: GRCh37->GRCh38 liftover sub-test — resources/hg19ToHg38.over.chain "
                    + "not present on this machine. This is a missing external resource, not a code failure; "
                    + "the rest of the suite still runs and is graded normally.");
            } else {
                GlobalSearchIndex.Result lifted = idx.searchPosition("19", 44908684, "GRCh38");
                check("GRCh37->GRCh38 liftover: rs429358 found at its lifted hg38 position", lifted.hits.size(), 1);
                if (!lifted.hits.isEmpty()) check("lifted hit is rs429358", lifted.hits.get(0).snpId, "rs429358");
            }

            // ── Unsupported GRCh38->GRCh37 direction ───────────────────────
            GlobalSearchIndex.Result unsupported = idx.searchPosition("2", 100000, "GRCh37");
            check("GRCh38 project's SNP never appears as a false hit under a GRCh37 query",
                containsProject(unsupported.hits, "proj_c"), false);
            boolean flagged = false;
            for (GlobalSearchIndex.NotCompared nc : unsupported.notCompared) {
                if ("proj_c".equals(nc.projectId)) {
                    flagged = true;
                    check("proj_c not-compared reason mentions liftover being unavailable",
                        nc.reason != null && nc.reason.toLowerCase(Locale.ROOT).contains("liftover unavailable"), true);
                }
            }
            check("proj_c (GRCh38) is explicitly flagged as not-compared for a GRCh37 query", flagged, true);

            if (failures == 0) {
                System.out.println("PASS: all GlobalSearchIndex tests passed");
            } else {
                System.out.println("FAIL: " + failures + " test(s) failed");
                System.exit(1);
            }
        } finally {
            deleteRecursive(root.toFile());
        }
    }

    private static boolean containsProject(List<GlobalSearchIndex.Hit> hits, String projectId) {
        for (GlobalSearchIndex.Hit h : hits) if (projectId.equals(h.projectId)) return true;
        return false;
    }

    // ── Fixture writers ──────────────────────────────────────────────────

    private static void writeProjectA(Path root) throws IOException {
        File dir = new File(root.toFile(), "proj_a");
        new File(dir, "data").mkdirs();
        writeConfig(dir, "GRCh37", "");
        writeFile(new File(dir, "project.json"), "{\"id\":\"proj_a\",\"name\":\"Project A Test\"}");
        writeFile(new File(dir, "data/manifest.json"),
            "{\"total_loci\":1,\"loci\":["
          + "{\"id\":\"a1\",\"index\":1,\"chr\":\"1\",\"start\":900,\"end\":9200,\"locus_name\":\"Locus 1\","
          + "\"top_snp\":\"rs1001\",\"top_snp_pval\":1e-9,\"nearest_genes\":[\"NEARBY1\",\"NEARBY2\"]}"
          + "]}");
        writeFile(new File(dir, "data/locus_1.json"),
            "{\"id\":\"a1\",\"locus_index\":1,\"locus_name\":\"Locus 1\",\"chr\":\"1\",\"start\":900,\"end\":9200,"
          + "\"padded_start\":700,\"padded_end\":9400,\"ref_panel\":\"\","
          + "\"top_snp\":{\"id\":\"rs1001\",\"chr\":\"1\",\"pos\":1500,\"pvalue\":1e-9},"
          + "\"nearest_genes\":[\"NEARBY1\",\"NEARBY2\"],"
          + "\"gwas_snps\":["
          +   "{\"id\":\"rs1001\",\"chr\":\"1\",\"pos\":1500,\"pvalue\":1e-9,\"beta\":0.2,\"se\":0.05,\"ea\":\"a\",\"nea\":\"g\"},"
          +   "{\"id\":\"rs1002\",\"chr\":\"1\",\"pos\":9000,\"pvalue\":2e-8,\"beta\":0.15,\"se\":0.04,\"ea\":\"c\",\"nea\":\"t\"},"
          +   "{\"id\":\"rs1003\",\"chr\":\"1\",\"pos\":1600,\"pvalue\":0.5,\"beta\":0.01,\"se\":0.02,\"ea\":\"a\",\"nea\":\"t\"}"
          + "],"
          + "\"genes\":["
          +   "{\"gene_name\":\"GENE_IN\",\"gene_id\":\"ENSG0001\",\"strand\":\"+\",\"start\":1000,\"end\":2000,\"transcripts\":[]},"
          +   "{\"gene_name\":\"GENE_OUT\",\"gene_id\":\"ENSG0002\",\"strand\":\"-\",\"start\":5000,\"end\":6000,\"transcripts\":[]}"
          + "],"
          + "\"ld_triangle\":null,\"locus_context\":null}");
    }

    private static void writeProjectB(Path root) throws IOException {
        File dir = new File(root.toFile(), "proj_b");
        new File(dir, "data").mkdirs();
        writeConfig(dir, "GRCh37", "beta");
        writeFile(new File(dir, "project.json"), "{\"id\":\"proj_b\",\"name\":\"Project B Test\"}");
        writeFile(new File(dir, "data/manifest.json"),
            "{\"total_loci\":1,\"loci\":["
          + "{\"id\":\"b1\",\"index\":1,\"chr\":\"19\",\"start\":45400000,\"end\":45420000,\"locus_name\":\"Locus 1\","
          + "\"top_snp\":\"rs429358\",\"top_snp_pval\":1e-10,\"nearest_genes\":[\"APOE\"]}"
          + "]}");
        writeFile(new File(dir, "data/locus_1.json"),
            "{\"id\":\"b1\",\"locus_index\":1,\"locus_name\":\"Locus 1\",\"chr\":\"19\",\"start\":45400000,\"end\":45420000,"
          + "\"padded_start\":45200000,\"padded_end\":45620000,\"ref_panel\":\"\","
          + "\"top_snp\":{\"id\":\"rs429358\",\"chr\":\"19\",\"pos\":45411941,\"pvalue\":1e-10},"
          + "\"nearest_genes\":[\"APOE\"],"
          + "\"gwas_snps\":["
          +   "{\"id\":\"rs429358\",\"chr\":\"19\",\"pos\":45411941,\"pvalue\":1e-10,\"beta\":0.3,\"se\":0.02,\"ea\":\"c\",\"nea\":\"t\"}"
          + "],"
          + "\"genes\":[{\"gene_name\":\"APOE\",\"gene_id\":\"ENSG00000130203\",\"strand\":\"+\",\"start\":45409039,\"end\":45412650,\"transcripts\":[]}],"
          + "\"ld_triangle\":null,\"locus_context\":null}");
    }

    private static void writeProjectC(Path root) throws IOException {
        File dir = new File(root.toFile(), "proj_c");
        new File(dir, "data").mkdirs();
        writeConfig(dir, "GRCh38", "OR");
        writeFile(new File(dir, "project.json"), "{\"id\":\"proj_c\",\"name\":\"Project C Test\"}");
        writeFile(new File(dir, "data/manifest.json"),
            "{\"total_loci\":1,\"loci\":["
          + "{\"id\":\"c1\",\"index\":1,\"chr\":\"2\",\"start\":99000,\"end\":101000,\"locus_name\":\"Locus 1\","
          + "\"top_snp\":\"rs777\",\"top_snp_pval\":5e-9,\"nearest_genes\":[\"SOMEGENE\"]}"
          + "]}");
        writeFile(new File(dir, "data/locus_1.json"),
            "{\"id\":\"c1\",\"locus_index\":1,\"locus_name\":\"Locus 1\",\"chr\":\"2\",\"start\":99000,\"end\":101000,"
          + "\"padded_start\":98800,\"padded_end\":101200,\"ref_panel\":\"\","
          + "\"top_snp\":{\"id\":\"rs777\",\"chr\":\"2\",\"pos\":100000,\"pvalue\":5e-9},"
          + "\"nearest_genes\":[\"SOMEGENE\"],"
          + "\"gwas_snps\":["
          +   "{\"id\":\"rs777\",\"chr\":\"2\",\"pos\":100000,\"pvalue\":5e-9,\"beta\":1.2,\"se\":0.1,\"ea\":\"g\",\"nea\":\"a\"}"
          + "],"
          + "\"genes\":[{\"gene_name\":\"SOMEGENE\",\"gene_id\":\"ENSGX\",\"strand\":\"+\",\"start\":99500,\"end\":99900,\"transcripts\":[]}],"
          + "\"ld_triangle\":null,\"locus_context\":null}");
    }

    private static void writeConfig(File dir, String genomeBuild, String effectType) throws IOException {
        StringBuilder sb = new StringBuilder();
        sb.append("gwas.file=dummy_gwas.txt\n");
        sb.append("loci.file=dummy_loci.txt\n");
        sb.append("gff3.file=dummy.gff3\n");
        sb.append("ref.panel.path=\n");
        sb.append("col.chr=chr\n");
        sb.append("col.pos=pos\n");
        sb.append("col.pvalue=p\n");
        sb.append("col.ea=ea\n");
        sb.append("col.nea=nea\n");
        sb.append("effect.type=").append(effectType).append("\n");
        sb.append("genome.build=").append(genomeBuild).append("\n");
        writeFile(new File(dir, "config.properties"), sb.toString());
    }

    private static void writeFile(File f, String content) throws IOException {
        f.getParentFile().mkdirs();
        try (OutputStreamWriter w = new OutputStreamWriter(new FileOutputStream(f), StandardCharsets.UTF_8)) {
            w.write(content);
        }
    }

    private static void deleteRecursive(File f) {
        File[] children = f.listFiles();
        if (children != null) for (File c : children) deleteRecursive(c);
        f.delete();
    }

    private static int check(String label, Object actual, Object expected) {
        if (!Objects.equals(actual, expected)) {
            System.out.println("FAIL: " + label + " — expected " + expected + ", got " + actual);
            failures++;
            return 1;
        }
        System.out.println("PASS: " + label + " (" + actual + ")");
        return 0;
    }
}
