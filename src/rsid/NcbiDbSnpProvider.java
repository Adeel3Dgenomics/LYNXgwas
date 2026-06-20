package rsid;

import java.io.*;
import java.net.*;
import java.nio.charset.StandardCharsets;

/**
 * NCBI dbSNP Variation Services provider.
 * Queries by genomic position, verifies alleles (forward/reverse).
 */
public class NcbiDbSnpProvider implements RsidApiProvider {

    private final String baseUrl;
    private final String apiKey;
    private final RateLimiter limiter;
    private final int timeoutMs;
    private final int maxRetries;

    public NcbiDbSnpProvider(String baseUrl, String apiKey, double ratePerSec,
                             int timeoutSec, int maxRetries) {
        this.baseUrl = baseUrl.endsWith("/") ? baseUrl : baseUrl + "/";
        this.apiKey = apiKey;
        this.limiter = new RateLimiter(ratePerSec);
        this.timeoutMs = timeoutSec * 1000;
        this.maxRetries = maxRetries;
    }

    @Override
    public String name() { return "NCBI"; }

    @Override
    public ApiResult lookup(String chr, long pos, String ea, String nea, String build) {
        // NCBI SPDI format: NC_00000X.Y:pos:ref:alt
        // Simpler: use the refsnp endpoint by position
        // GET /variation/v0/beta/genome/{build}/{chr}/{pos}/all
        String assembly = build.equalsIgnoreCase("hg38") ? "GCF_000001405.40" : "GCF_000001405.25";
        String ncbiChr = chr.startsWith("chr") ? chr : "chr" + chr;

        // Try NCBI SPDI overlap endpoint
        String url = baseUrl + "spdi/" + assembly + "/" + ncbiChr + "/" + pos + "/overlapping_variants";
        if (!apiKey.isEmpty()) url += "?api_key=" + apiKey;

        for (int attempt = 0; attempt <= maxRetries; attempt++) {
            try {
                limiter.acquire();
                String json = httpGet(url);
                if (json == null) return new ApiResult(null, "api_no_match", name());

                // Parse rsIDs from response and verify alleles
                return parseAndVerify(json, ea, nea);

            } catch (InterruptedException e) {
                Thread.currentThread().interrupt();
                return new ApiResult(null, "api_error", name());
            } catch (IOException e) {
                if (attempt < maxRetries && isRetryable(e)) {
                    try { Thread.sleep((long) Math.pow(2, attempt) * 1000); }
                    catch (InterruptedException ie) { Thread.currentThread().interrupt(); break; }
                    continue;
                }
                return new ApiResult(null, "api_error", name());
            }
        }
        return new ApiResult(null, "api_error", name());
    }

    private ApiResult parseAndVerify(String json, String ea, String nea) {
        String eaUp = ea.toUpperCase(), neaUp = nea.toUpperCase();

        // Look for rsID patterns in the response
        int idx = 0;
        String bestRsid = null;
        String bestReason = "api_no_match";

        while (true) {
            int rsPos = json.indexOf("\"rs", idx);
            if (rsPos < 0) break;

            // Extract rsID
            int rsStart = rsPos + 1;
            int rsEnd = rsStart;
            while (rsEnd < json.length() && json.charAt(rsEnd) != '"' && json.charAt(rsEnd) != ',') rsEnd++;
            String rsid = json.substring(rsStart, rsEnd);

            if (rsid.matches("rs\\d+")) {
                // Look for allele info near this rsID
                int searchEnd = Math.min(json.length(), rsEnd + 500);
                String context = json.substring(Math.max(0, rsPos - 200), searchEnd);

                // Try to extract ref/alt alleles from context
                String ref = extractAllele(context, "deleted_sequence");
                if (ref == null) ref = extractAllele(context, "ref");
                String alt = extractAllele(context, "inserted_sequence");
                if (alt == null) alt = extractAllele(context, "alt");

                if (ref != null && alt != null) {
                    ref = ref.toUpperCase();
                    alt = alt.toUpperCase();
                    // Forward: nea==ref, ea==alt
                    if (neaUp.equals(ref) && eaUp.equals(alt)) {
                        return new ApiResult(rsid, "api_matched_forward", name());
                    }
                    // Reverse: ea==ref, nea==alt
                    if (eaUp.equals(ref) && neaUp.equals(alt)) {
                        return new ApiResult(rsid, "api_matched_reverse", name());
                    }
                    bestRsid = rsid;
                    bestReason = "api_pos_only_mismatch";
                } else {
                    if (bestRsid == null) bestRsid = rsid;
                }
            }
            idx = rsEnd;
        }

        if (bestRsid != null && bestReason.equals("api_no_match")) {
            return new ApiResult(bestRsid, "api_pos_only_mismatch", name());
        }
        return new ApiResult(bestRsid, bestReason, name());
    }

    private String extractAllele(String context, String field) {
        String search = "\"" + field + "\":\"";
        int i = context.indexOf(search);
        if (i < 0) return null;
        i += search.length();
        int end = context.indexOf('"', i);
        if (end < 0) return null;
        return context.substring(i, end);
    }

    private String httpGet(String url) throws IOException {
        HttpURLConnection conn = (HttpURLConnection) new URL(url).openConnection();
        conn.setRequestMethod("GET");
        conn.setConnectTimeout(timeoutMs);
        conn.setReadTimeout(timeoutMs);
        conn.setRequestProperty("Accept", "application/json");
        int code = conn.getResponseCode();
        if (code == 404) return null;
        if (code == 429 || code >= 500) throw new IOException("HTTP " + code);
        if (code != 200) return null;
        try (InputStream is = conn.getInputStream()) {
            return new String(is.readAllBytes(), StandardCharsets.UTF_8);
        }
    }

    private boolean isRetryable(IOException e) {
        String msg = e.getMessage();
        return msg != null && (msg.contains("429") || msg.contains("500") || msg.contains("503") || msg.contains("timed out"));
    }
}
