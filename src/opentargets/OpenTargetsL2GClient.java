package opentargets;

import rsid.RateLimiter;

import java.io.*;
import java.net.URI;
import java.net.http.*;
import java.nio.charset.StandardCharsets;
import java.nio.file.*;
import java.time.Duration;
import java.util.*;

/**
 * Locus-to-gene (L2G) causal-gene prediction via the live Open Targets Platform GraphQL API
 * (api.platform.opentargets.org/api/v4/graphql). Unlike the GWAS Catalog associations, L2G scores are
 * NOT available as a plain bulk file — the Platform only publishes them as ~554MB of Apache Parquet,
 * which would need a real Parquet-decoding library (a first for this codebase; every other data source
 * here is hand-parsed text). The live API avoids that dependency entirely, and — unlike the Catalog's
 * REST API — nests everything in one query: search(rsid) resolves to a variant ID, then
 * variant -> credibleSets -> l2GPredictions comes back in a single round trip (verified live, ~0.3s for
 * a real, heavily-studied locus). So this stays two HTTP calls per *locus* (not per SNP, not per hit),
 * cached to disk, rather than the N+1 explosion that made the Catalog's live REST API unusable.
 *
 * L2G only exists for loci Open Targets has already ingested from a specific published study — this
 * returns an empty gene list (not an error) for anything it hasn't seen, including any genuinely novel
 * locus by construction.
 */
public class OpenTargetsL2GClient {

    private static final String GRAPHQL_URL = "https://api.platform.opentargets.org/api/v4/graphql";
    private static final RateLimiter LIMITER = new RateLimiter(3.0);
    private final HttpClient http = HttpClient.newBuilder()
        .connectTimeout(Duration.ofSeconds(10))
        .version(HttpClient.Version.HTTP_1_1) // same GOAWAY-after-one-request issue seen on the GWAS Catalog API
        .build();
    private final Path cacheDir;

    public OpenTargetsL2GClient(String projectOutputDir) {
        this.cacheDir = Paths.get(projectOutputDir, "cache", "opentargets_l2g");
    }

    public static class GenePrediction {
        public String gene;
        public double maxScore;
        public int supportingStudies;
    }

    /** One tissue/cell-type's enhancer-activity evidence for a gene (ENCODE E2G / ABC-model-based). */
    public static class EnhancerEvidence {
        public String biosample;
        public double score;
        public int distanceToTss;
        public String pmid;
    }

    /** One gene's aggregated enhancer-activity evidence — narrower and patchier than L2G (it depends on
     *  which cell types have been experimentally profiled, so a real, well-studied locus can still come
     *  back empty), but adds *which tissue* the link is active in, which the L2G score alone doesn't say. */
    public static class EnhancerGenePrediction {
        public String gene;
        public double maxScore;
        public List<EnhancerEvidence> evidence = new ArrayList<>(); // one per biosample/study, sorted by score desc
    }

    public static class Result {
        public String rsid;
        public String variantId;      // resolved Open Targets variant ID, or null if unmatched
        public int credibleSetCount;  // how many prior studies' credible sets overlap this variant
        public List<GenePrediction> genes = new ArrayList<>();          // sorted by maxScore descending
        public List<EnhancerGenePrediction> enhancerGenes = new ArrayList<>(); // sorted by maxScore descending
    }

    /** Looks up the most likely causal gene(s) for this rsID's locus. Empty gene list, not null/error,
     *  when Open Targets has no data for it (including any genuinely novel locus). */
    public Result lookup(String rsid) throws IOException, InterruptedException {
        String cached = readCache(rsid);
        if (cached != null) return parseCache(rsid, cached);

        Result r = fetchLive(rsid);
        writeCache(rsid, r);
        return r;
    }

    private Result fetchLive(String rsid) throws IOException, InterruptedException {
        Result r = new Result();
        r.rsid = rsid;

        String searchBody = "{\"query\":\"query($q: String!) { search(queryString: $q, entityNames: [\\\"variant\\\"]) { hits { id } } }\",\"variables\":{\"q\":\"" + escJson(rsid) + "\"}}";
        String searchResp = post(searchBody);
        String variantId = firstStringField(extractArrayBlock(searchResp, "hits"), "id");
        if (variantId == null) return r; // no matching variant at all — nothing more to do

        r.variantId = variantId;
        // enhancerToGenes rides along in this SAME query — Open Targets exposes it directly off
        // `variant`, so adding it costs zero extra HTTP round trips over the L2G-only version.
        String variantQuery = "{\"query\":\"query($id: String!) { variant(variantId: $id) { " +
            "credibleSets { count rows { l2GPredictions { rows { score target { approvedSymbol } } } } } " +
            "enhancerToGenes { rows { score target { approvedSymbol } biosampleName distanceToTss pmid } } " +
            "} }\",\"variables\":{\"id\":\"" + escJson(variantId) + "\"}}";
        String variantResp = post(variantQuery);
        String variantBlock = extractObjectBlock(variantResp, "variant");

        String credibleSetsBlock = extractObjectBlock(variantBlock, "credibleSets");
        if (credibleSetsBlock != null) {
            String countStr = firstStringOrNumberField(credibleSetsBlock, "count");
            try { r.credibleSetCount = countStr != null ? Integer.parseInt(countStr) : 0; } catch (NumberFormatException ignored) {}

            Map<String, double[]> geneAgg = new LinkedHashMap<>(); // gene -> {maxScore, studyCount}
            for (String csRow : splitJsonObjects(extractArrayBlock(credibleSetsBlock, "rows"))) {
                String l2gBlock = extractObjectBlock(csRow, "l2GPredictions");
                if (l2gBlock == null) continue;
                for (String predRow : splitJsonObjects(extractArrayBlock(l2gBlock, "rows"))) {
                    String targetBlock = extractObjectBlock(predRow, "target");
                    String gene = targetBlock != null ? firstStringField(targetBlock, "approvedSymbol") : null;
                    String scoreStr = firstStringOrNumberField(predRow, "score");
                    if (gene == null || scoreStr == null) continue;
                    double score;
                    try { score = Double.parseDouble(scoreStr); } catch (NumberFormatException e) { continue; }
                    double[] agg = geneAgg.computeIfAbsent(gene, k -> new double[]{0, 0});
                    if (score > agg[0]) agg[0] = score;
                    agg[1] += 1;
                }
            }
            for (Map.Entry<String, double[]> e : geneAgg.entrySet()) {
                GenePrediction gp = new GenePrediction();
                gp.gene = e.getKey();
                gp.maxScore = e.getValue()[0];
                gp.supportingStudies = (int) e.getValue()[1];
                r.genes.add(gp);
            }
            r.genes.sort((a, b) -> Double.compare(b.maxScore, a.maxScore));
        }

        String enhancerBlock = extractObjectBlock(variantBlock, "enhancerToGenes");
        if (enhancerBlock != null) {
            Map<String, EnhancerGenePrediction> enhByGene = new LinkedHashMap<>();
            for (String row : splitJsonObjects(extractArrayBlock(enhancerBlock, "rows"))) {
                String targetBlock = extractObjectBlock(row, "target");
                String gene = targetBlock != null ? firstStringField(targetBlock, "approvedSymbol") : null;
                String scoreStr = firstStringOrNumberField(row, "score");
                if (gene == null || scoreStr == null) continue;
                double score;
                try { score = Double.parseDouble(scoreStr); } catch (NumberFormatException e) { continue; }

                EnhancerEvidence ev = new EnhancerEvidence();
                ev.biosample = firstStringField(row, "biosampleName");
                ev.score = score;
                String distStr = firstStringOrNumberField(row, "distanceToTss");
                try { ev.distanceToTss = distStr != null ? Integer.parseInt(distStr) : 0; } catch (NumberFormatException ignored) {}
                ev.pmid = firstStringField(row, "pmid");

                EnhancerGenePrediction egp = enhByGene.computeIfAbsent(gene, k -> {
                    EnhancerGenePrediction g = new EnhancerGenePrediction();
                    g.gene = k;
                    return g;
                });
                if (score > egp.maxScore) egp.maxScore = score;
                egp.evidence.add(ev);
            }
            r.enhancerGenes = new ArrayList<>(enhByGene.values());
            for (EnhancerGenePrediction egp : r.enhancerGenes)
                egp.evidence.sort((a, b) -> Double.compare(b.score, a.score));
            r.enhancerGenes.sort((a, b) -> Double.compare(b.maxScore, a.maxScore));
        }

        return r;
    }

    // ---- HTTP -------------------------------------------------------------------------------

    private String post(String jsonBody) throws IOException, InterruptedException {
        LIMITER.acquire();
        HttpRequest req = HttpRequest.newBuilder(URI.create(GRAPHQL_URL))
            .header("Content-Type", "application/json")
            .timeout(Duration.ofSeconds(15))
            .POST(HttpRequest.BodyPublishers.ofString(jsonBody, StandardCharsets.UTF_8))
            .build();
        HttpResponse<String> resp = http.send(req, HttpResponse.BodyHandlers.ofString(StandardCharsets.UTF_8));
        if (resp.statusCode() != 200)
            throw new IOException("Open Targets returned HTTP " + resp.statusCode() + " for " + jsonBody);
        return resp.body();
    }

    // ---- disk cache (one small file per rsID) ------------------------------------------------

    private String readCache(String rsid) {
        Path f = cacheDir.resolve(safe(rsid) + ".tsv");
        if (!Files.exists(f)) return null;
        try { return Files.readString(f, StandardCharsets.UTF_8); }
        catch (IOException e) { return null; }
    }

    private static final String ENHANCER_SECTION_MARKER = "---ENHANCER---";

    private void writeCache(String rsid, Result r) {
        try {
            Files.createDirectories(cacheDir);
            StringBuilder sb = new StringBuilder();
            sb.append(r.variantId == null ? "" : r.variantId).append('\t').append(r.credibleSetCount).append('\n');
            for (GenePrediction gp : r.genes)
                sb.append(gp.gene).append('\t').append(gp.maxScore).append('\t').append(gp.supportingStudies).append('\n');
            sb.append(ENHANCER_SECTION_MARKER).append('\n');
            // One line per (gene, tissue) evidence row — not one line per gene — since each row also
            // carries its own score/distance/PMID, not just a bag of tissue names.
            for (EnhancerGenePrediction egp : r.enhancerGenes)
                for (EnhancerEvidence ev : egp.evidence)
                    sb.append(egp.gene).append('\t').append(ev.score).append('\t')
                      .append(nzTsv(ev.biosample)).append('\t').append(ev.distanceToTss).append('\t')
                      .append(nzTsv(ev.pmid)).append('\n');
            Files.writeString(cacheDir.resolve(safe(rsid) + ".tsv"), sb.toString(), StandardCharsets.UTF_8);
        } catch (IOException e) {
            System.err.println("[OpenTargetsL2GClient] cache write failed: " + e.getMessage());
        }
    }

    private static String nzTsv(String s) { return s == null ? "" : s.replace("\t", " ").replace("\n", " "); }

    private Result parseCache(String rsid, String cached) {
        Result r = new Result();
        r.rsid = rsid;
        String[] lines = cached.split("\n", -1);
        if (lines.length > 0 && !lines[0].isEmpty()) {
            String[] head = lines[0].split("\t", -1);
            r.variantId = head.length > 0 && !head[0].isEmpty() ? head[0] : null;
            try { r.credibleSetCount = head.length > 1 ? Integer.parseInt(head[1]) : 0; } catch (NumberFormatException ignored) {}
        }
        boolean inEnhancerSection = false;
        Map<String, EnhancerGenePrediction> enhByGene = new LinkedHashMap<>();
        for (int i = 1; i < lines.length; i++) {
            if (lines[i].equals(ENHANCER_SECTION_MARKER)) { inEnhancerSection = true; continue; }
            if (lines[i].isEmpty()) continue;
            if (inEnhancerSection) {
                String[] parts = lines[i].split("\t", -1);
                if (parts.length < 5) continue;
                EnhancerEvidence ev = new EnhancerEvidence();
                try { ev.score = Double.parseDouble(parts[1]); } catch (NumberFormatException ignored) {}
                ev.biosample = parts[2].isEmpty() ? null : parts[2];
                try { ev.distanceToTss = Integer.parseInt(parts[3]); } catch (NumberFormatException ignored) {}
                ev.pmid = parts[4].isEmpty() ? null : parts[4];
                String gene = parts[0];
                EnhancerGenePrediction egp = enhByGene.computeIfAbsent(gene, k -> {
                    EnhancerGenePrediction g = new EnhancerGenePrediction();
                    g.gene = k;
                    return g;
                });
                if (ev.score > egp.maxScore) egp.maxScore = ev.score;
                egp.evidence.add(ev);
                continue;
            }
            String[] parts = lines[i].split("\t", -1);
            if (parts.length < 3) continue;
            GenePrediction gp = new GenePrediction();
            gp.gene = parts[0];
            try { gp.maxScore = Double.parseDouble(parts[1]); } catch (NumberFormatException ignored) {}
            try { gp.supportingStudies = Integer.parseInt(parts[2]); } catch (NumberFormatException ignored) {}
            r.genes.add(gp);
        }
        r.enhancerGenes = new ArrayList<>(enhByGene.values());
        for (EnhancerGenePrediction egp : r.enhancerGenes) egp.evidence.sort((a, b) -> Double.compare(b.score, a.score));
        r.enhancerGenes.sort((a, b) -> Double.compare(b.maxScore, a.maxScore));
        return r;
    }

    private static String safe(String key) { return key.replaceAll("[^A-Za-z0-9_.-]", "_"); }
    private static String escJson(String s) { return s.replace("\\", "\\\\").replace("\"", "\\\""); }

    // ---- minimal hand-rolled JSON extraction (same no-library convention as the rest of this codebase) --

    private static String extractArrayBlock(String json, String key) {
        if (json == null) return "[]";
        int i = json.indexOf("\"" + key + "\"");
        if (i < 0) return "[]";
        int start = json.indexOf('[', i);
        if (start < 0) return "[]";
        int depth = 0;
        for (int j = start; j < json.length(); j++) {
            char c = json.charAt(j);
            if (c == '[') depth++;
            else if (c == ']') { depth--; if (depth == 0) return json.substring(start, j + 1); }
        }
        return "[]";
    }

    private static String extractObjectBlock(String json, String key) {
        if (json == null) return null;
        int i = json.indexOf("\"" + key + "\"");
        if (i < 0) return null;
        int start = json.indexOf('{', i);
        if (start < 0) return null;
        int depth = 0;
        for (int j = start; j < json.length(); j++) {
            char c = json.charAt(j);
            if (c == '{') depth++;
            else if (c == '}') { depth--; if (depth == 0) return json.substring(start, j + 1); }
        }
        return null;
    }

    /** Splits a top-level JSON array of objects into individual object strings (brace-depth aware). */
    private static List<String> splitJsonObjects(String arrayJson) {
        List<String> out = new ArrayList<>();
        int depth = 0, start = -1;
        for (int i = 0; i < arrayJson.length(); i++) {
            char c = arrayJson.charAt(i);
            if (c == '{') { if (depth == 0) start = i; depth++; }
            else if (c == '}') { depth--; if (depth == 0 && start >= 0) out.add(arrayJson.substring(start, i + 1)); }
        }
        return out;
    }

    private static String firstStringField(String json, String key) {
        if (json == null) return null;
        int i = json.indexOf("\"" + key + "\"");
        if (i < 0) return null;
        int colon = json.indexOf(':', i);
        int q1 = json.indexOf('"', colon + 1);
        if (q1 < 0) return null;
        int q2 = json.indexOf('"', q1 + 1);
        while (q2 > 0 && json.charAt(q2 - 1) == '\\') q2 = json.indexOf('"', q2 + 1);
        if (q2 < 0) return null;
        return json.substring(q1 + 1, q2).replace("\\\"", "\"");
    }

    private static String firstStringOrNumberField(String json, String key) {
        if (json == null) return null;
        int i = json.indexOf("\"" + key + "\"");
        if (i < 0) return null;
        int colon = json.indexOf(':', i);
        int j = colon + 1;
        while (j < json.length() && Character.isWhitespace(json.charAt(j))) j++;
        if (j < json.length() && json.charAt(j) == '"') return firstStringField(json, key);
        int end = j;
        while (end < json.length() && "-+.0123456789eE".indexOf(json.charAt(end)) >= 0) end++;
        if (end == j) return null;
        return json.substring(j, end);
    }
}
