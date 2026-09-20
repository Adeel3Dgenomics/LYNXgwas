package catalog;

import rsid.RateLimiter;

import java.io.*;
import java.net.URI;
import java.net.http.*;
import java.nio.charset.StandardCharsets;
import java.nio.file.*;
import java.time.Duration;
import java.util.*;

/**
 * Cross-phenotype lookup against the public EBI GWAS Catalog REST API (www.ebi.ac.uk/gwas/rest/api).
 * Given a lead SNP's rsID, returns every trait association already reported for that variant —
 * the "what else is this locus known for" panel.
 *
 * API shape (verified live 2026-09-19): the SNP associations endpoint returns each association's
 * p-value inline but NOT its trait name or study — those need one follow-up GET per association
 * (/associations/{id}/efoTraits, /associations/{id}/study). No projection embeds all three at once.
 * That N+1 shape is why every lookup here is disk-cached per rsID and rate-limited, the same pattern
 * already used for the NCBI/gnomAD rsID-recovery providers in the rsid/ package.
 */
public class GwasCatalogClient {

    private static final String BASE = "https://www.ebi.ac.uk/gwas/rest/api";
    private static final RateLimiter LIMITER = new RateLimiter(3.0); // polite default, matches NCBI provider's default
    private final HttpClient http = HttpClient.newBuilder()
        .connectTimeout(Duration.ofSeconds(10))
        .version(HttpClient.Version.HTTP_1_1) // the Catalog's HTTP/2 endpoint GOAWAYs this client after 1 request
        .build();
    private final Path cacheDir;

    public GwasCatalogClient(String projectOutputDir) {
        this.cacheDir = Paths.get(projectOutputDir, "cache", "gwas_catalog");
    }

    public static class Hit {
        public String trait;
        public String efoId;
        public double pvalue;
        public String pubmedId;
        public String studyAccession;
        public String initialSampleSize;
        public String reportedTrait;

        @Override public String toString() {
            return String.format("%-45s p=%.2e  study=%s  PMID=%s  n=%s",
                trait, pvalue, studyAccession, pubmedId, initialSampleSize);
        }
    }

    /** Returns every reported association for this rsID, most significant first. Empty list, not null, if none found. */
    public List<Hit> lookup(String rsid) throws IOException, InterruptedException {
        String cached = readCache(rsid);
        String assocJson = cached != null ? cached : fetchAssociationsRaw(rsid);
        if (cached == null) writeCache(rsid, assocJson);

        List<Hit> hits = new ArrayList<>();
        for (String assocBlock : splitJsonObjects(extractArrayBlock(assocJson, "associations"))) {
            Hit h = new Hit();
            h.pvalue = parsePvalue(assocBlock);
            String assocHref = extractHref(assocBlock, "self");
            if (assocHref == null) continue;
            String assocId = assocHref.substring(assocHref.lastIndexOf('/') + 1);

            String efoJson = fetchWithCache(rsid + "_efo_" + assocId,
                BASE + "/associations/" + assocId + "/efoTraits");
            h.trait = firstStringField(efoJson, "trait");
            h.efoId = firstStringField(efoJson, "shortForm");

            String studyJson = fetchWithCache(rsid + "_study_" + assocId,
                BASE + "/associations/" + assocId + "/study");
            h.pubmedId = extractNestedString(studyJson, "publicationInfo", "pubmedId");
            h.studyAccession = firstStringField(studyJson, "accessionId");
            h.initialSampleSize = firstStringField(studyJson, "initialSampleSize");
            h.reportedTrait = extractNestedString(studyJson, "diseaseTrait", "trait");

            if (h.trait == null) h.trait = h.reportedTrait != null ? h.reportedTrait : "(trait unavailable)";
            hits.add(h);
        }
        hits.sort(Comparator.comparingDouble(h -> h.pvalue));
        return hits;
    }

    // ---- HTTP + tiny disk cache ------------------------------------------------

    private String fetchAssociationsRaw(String rsid) throws IOException, InterruptedException {
        return get(BASE + "/singleNucleotidePolymorphisms/" + rsid + "/associations");
    }

    private String fetchWithCache(String cacheKey, String url) throws IOException, InterruptedException {
        String cached = readCache(cacheKey);
        if (cached != null) return cached;
        String body = get(url);
        writeCache(cacheKey, body);
        return body;
    }

    private String get(String url) throws IOException, InterruptedException {
        LIMITER.acquire();
        HttpRequest req = HttpRequest.newBuilder(URI.create(url))
            .header("Accept", "application/json")
            .timeout(Duration.ofSeconds(15))
            .GET().build();
        HttpResponse<String> resp = http.send(req, HttpResponse.BodyHandlers.ofString(StandardCharsets.UTF_8));
        if (resp.statusCode() == 404) return "{}";
        if (resp.statusCode() != 200)
            throw new IOException("GWAS Catalog returned HTTP " + resp.statusCode() + " for " + url);
        return resp.body();
    }

    private String readCache(String key) {
        Path f = cacheDir.resolve(safe(key) + ".json");
        if (!Files.exists(f)) return null;
        try { return Files.readString(f, StandardCharsets.UTF_8); }
        catch (IOException e) { return null; }
    }

    private void writeCache(String key, String body) {
        try {
            Files.createDirectories(cacheDir);
            Files.writeString(cacheDir.resolve(safe(key) + ".json"), body, StandardCharsets.UTF_8);
        } catch (IOException e) {
            System.err.println("[GwasCatalogClient] cache write failed: " + e.getMessage());
        }
    }

    private static String safe(String key) { return key.replaceAll("[^A-Za-z0-9_.-]", "_"); }

    // ---- minimal hand-rolled JSON field extraction (same no-library convention as JsonExporter/LocalServer) --

    private static double parsePvalue(String block) {
        String pStr = firstStringOrNumberField(block, "pvalue");
        if (pStr != null) { try { return Double.parseDouble(pStr); } catch (NumberFormatException ignored) {} }
        String mantissa = firstStringOrNumberField(block, "pvalueMantissa");
        String exponent = firstStringOrNumberField(block, "pvalueExponent");
        if (mantissa != null && exponent != null) {
            try { return Double.parseDouble(mantissa) * Math.pow(10, Double.parseDouble(exponent)); }
            catch (NumberFormatException ignored) {}
        }
        return Double.NaN;
    }

    private static String extractArrayBlock(String json, String key) {
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

    private static String extractHref(String block, String linkName) {
        int linksIdx = block.indexOf("\"_links\"");
        if (linksIdx < 0) return null;
        int nameIdx = block.indexOf("\"" + linkName + "\"", linksIdx);
        if (nameIdx < 0) return null;
        return firstStringField(block.substring(nameIdx), "href");
    }

    private static String firstStringField(String json, String key) {
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

    private static String extractNestedString(String json, String objectKey, String fieldKey) {
        int i = json.indexOf("\"" + objectKey + "\"");
        if (i < 0) return null;
        int braceStart = json.indexOf('{', i);
        if (braceStart < 0) return null;
        int depth = 0, braceEnd = braceStart;
        for (int j = braceStart; j < json.length(); j++) {
            char c = json.charAt(j);
            if (c == '{') depth++;
            else if (c == '}') { depth--; if (depth == 0) { braceEnd = j; break; } }
        }
        return firstStringField(json.substring(braceStart, braceEnd + 1), fieldKey);
    }
}
