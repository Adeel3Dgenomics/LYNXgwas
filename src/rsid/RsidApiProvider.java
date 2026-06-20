package rsid;

import java.util.List;

/**
 * Interface for external rsID lookup providers (NCBI, gnomAD).
 * Each provider queries by chr:pos:ref:alt and returns allele-verified rsIDs.
 */
public interface RsidApiProvider {

    String name();

    /**
     * Look up rsID for a single SNP. Returns null if no allele-verified match.
     * Must handle its own rate limiting, retries, and timeouts.
     */
    ApiResult lookup(String chr, long pos, String ea, String nea, String build);

    class ApiResult {
        public String rsid;
        public String matchReason;  // api_matched_forward, api_matched_reverse, api_pos_only_mismatch, api_no_match, api_error
        public String provider;

        public ApiResult(String rsid, String matchReason, String provider) {
            this.rsid = rsid;
            this.matchReason = matchReason;
            this.provider = provider;
        }
    }
}
