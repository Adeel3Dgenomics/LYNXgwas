package rsid;

import java.io.*;
import java.net.*;
import java.nio.charset.StandardCharsets;

/**
 * gnomAD GraphQL provider for rsID lookup.
 * Queries by chr-pos-ref-alt variant ID, verifies alleles via ID matching.
 * // TODO: BUILD — verify gnomAD dataset version matches project build
 */
public class GnomadProvider implements RsidApiProvider {

    private final String graphqlUrl;
    private final String datasetHg19;
    private final String datasetHg38;
    private final RateLimiter limiter;
    private final int timeoutMs;
    private final int maxRetries;

    public GnomadProvider(String graphqlUrl, String datasetHg19, String datasetHg38,
                          double ratePerSec, int timeoutSec, int maxRetries) {
        this.graphqlUrl = graphqlUrl;
        this.datasetHg19 = datasetHg19;
        this.datasetHg38 = datasetHg38;
        this.limiter = new RateLimiter(ratePerSec);
        this.timeoutMs = timeoutSec * 1000;
        this.maxRetries = maxRetries;
    }

    @Override
    public String name() { return "gnomAD"; }

    @Override
    public ApiResult lookup(String chr, long pos, String ea, String nea, String build) {
        String dataset = build.equalsIgnoreCase("hg38") ? datasetHg38 : datasetHg19;

        // Try both orientations: chr-pos-nea-ea (forward) and chr-pos-ea-nea (reverse)
        ApiResult forward = queryVariant(chr, pos, nea.toUpperCase(), ea.toUpperCase(), dataset, "api_matched_forward");
        if (forward != null && forward.rsid != null) return forward;

        ApiResult reverse = queryVariant(chr, pos, ea.toUpperCase(), nea.toUpperCase(), dataset, "api_matched_reverse");
        if (reverse != null && reverse.rsid != null) return reverse;

        return new ApiResult(null, "api_no_match", name());
    }

    private ApiResult queryVariant(String chr, long pos, String ref, String alt,
                                   String dataset, String matchReason) {
        // gnomAD variant ID format: chr-pos-ref-alt
        String variantId = chr + "-" + pos + "-" + ref + "-" + alt;
        String query = "{\"query\":\"{ variant(variantId: \\\"" + variantId
            + "\\\", dataset: " + dataset + ") { rsids variant_id } }\"}";

        for (int attempt = 0; attempt <= maxRetries; attempt++) {
            try {
                limiter.acquire();
                String json = httpPost(graphqlUrl, query);
                if (json == null) return null;

                // Extract rsid from response
                String rsid = extractRsid(json);
                if (rsid != null) {
                    return new ApiResult(rsid, matchReason, name());
                }
                return null;

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

    private String extractRsid(String json) {
        // Look for "rsids":["rs12345"] or "rsids":["rs12345","rs67890"]
        int rsidsIdx = json.indexOf("\"rsids\"");
        if (rsidsIdx < 0) return null;
        int arrStart = json.indexOf('[', rsidsIdx);
        if (arrStart < 0) return null;
        int arrEnd = json.indexOf(']', arrStart);
        if (arrEnd < 0) return null;
        String arr = json.substring(arrStart + 1, arrEnd);
        // Extract first rsID
        int q1 = arr.indexOf('"');
        if (q1 < 0) return null;
        int q2 = arr.indexOf('"', q1 + 1);
        if (q2 < 0) return null;
        String rsid = arr.substring(q1 + 1, q2);
        if (rsid.matches("rs\\d+")) return rsid;
        return null;
    }

    private String httpPost(String url, String body) throws IOException {
        HttpURLConnection conn = (HttpURLConnection) new URL(url).openConnection();
        conn.setRequestMethod("POST");
        conn.setDoOutput(true);
        conn.setConnectTimeout(timeoutMs);
        conn.setReadTimeout(timeoutMs);
        conn.setRequestProperty("Content-Type", "application/json");
        conn.setRequestProperty("Accept", "application/json");
        try (OutputStream os = conn.getOutputStream()) {
            os.write(body.getBytes(StandardCharsets.UTF_8));
        }
        int code = conn.getResponseCode();
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
