package rsid;

import java.util.*;
import java.util.concurrent.atomic.*;

/**
 * Orchestrates API-based rsID completion for unmatched SNPs.
 * Processes in batches, tries providers in order, caches results.
 */
public class RsidApiCompleter {

    public static class Config {
        public boolean enabled = false;
        public int maxUnmatched = 2000;
        public int batchSize = 100;
        public String[] providerOrder = {"ncbi", "gnomad"};
        public String ncbiBaseUrl = "https://api.ncbi.nlm.nih.gov/variation/v0/";
        public String ncbiApiKey = "";
        public double ncbiRate = 3;
        public String gnomadUrl = "https://gnomad.broadinstitute.org/api";
        public String gnomadDatasetHg19 = "gnomad_r2_1";
        public String gnomadDatasetHg38 = "gnomad_r4";
        public double gnomadRate = 5;
        public int timeoutSec = 20;
        public int maxRetries = 3;
    }

    public static class UnmatchedSnp {
        public String chr;
        public long pos;
        public String ea, nea, varid;
    }

    public static class CompletionResult {
        public int attempted = 0;
        public int apiMatched = 0;
        public int apiNoMatch = 0;
        public int apiError = 0;
        public Map<String, String> recovered = new LinkedHashMap<>(); // chr:pos → rsid
    }

    /**
     * Run API completion on a list of unmatched SNPs.
     * @return map of chr:pos → rsid for successfully recovered SNPs
     */
    public static CompletionResult complete(List<UnmatchedSnp> unmatched, String build,
                                            Config cfg, RsidApiCache cache,
                                            RsidProgress progress) {
        CompletionResult result = new CompletionResult();

        // Threshold gate
        if (unmatched.size() >= cfg.maxUnmatched) {
            System.out.printf("[RsidApi] API completion SKIPPED: %,d unmatched >= cap %d.%n" +
                "               Likely a data/build issue — review before enabling API.%n",
                unmatched.size(), cfg.maxUnmatched);
            return result;
        }

        System.out.printf("[RsidApi] API completion: ENABLED · cap %d · %d unmatched → eligible%n",
            cfg.maxUnmatched, unmatched.size());

        // Build providers
        List<RsidApiProvider> providers = new ArrayList<>();
        for (String name : cfg.providerOrder) {
            if ("ncbi".equalsIgnoreCase(name)) {
                providers.add(new NcbiDbSnpProvider(cfg.ncbiBaseUrl, cfg.ncbiApiKey,
                    cfg.ncbiRate, cfg.timeoutSec, cfg.maxRetries));
            } else if ("gnomad".equalsIgnoreCase(name)) {
                providers.add(new GnomadProvider(cfg.gnomadUrl, cfg.gnomadDatasetHg19,
                    cfg.gnomadDatasetHg38, cfg.gnomadRate, cfg.timeoutSec, cfg.maxRetries));
            }
        }

        if (providers.isEmpty()) {
            System.out.println("[RsidApi] No providers configured, skipping");
            return result;
        }

        String gnomadDataset = build.equalsIgnoreCase("hg38") ? cfg.gnomadDatasetHg38 : cfg.gnomadDatasetHg19;
        System.out.printf("[RsidApi] Provider order: %s · batch %d · build %s (gnomad %s)%n",
            String.join(" → ", cfg.providerOrder), cfg.batchSize, build, gnomadDataset);

        // Process in batches
        int totalBatches = (unmatched.size() + cfg.batchSize - 1) / cfg.batchSize;
        for (int b = 0; b < totalBatches; b++) {
            if (progress != null && progress.cancelled) break;

            int from = b * cfg.batchSize;
            int to = Math.min(from + cfg.batchSize, unmatched.size());
            List<UnmatchedSnp> batch = unmatched.subList(from, to);

            System.out.printf("[RsidApi] Batch %d/%d (%d SNPs)%n", b + 1, totalBatches, batch.size());

            for (UnmatchedSnp snp : batch) {
                result.attempted++;

                // Check cache first
                String cached = cache.get(build, snp.chr, snp.pos, snp.ea, snp.nea);
                if (cached != null) {
                    if (!cached.isEmpty()) {
                        result.apiMatched++;
                        result.recovered.put(snp.chr + ":" + snp.pos, cached);
                    } else {
                        result.apiNoMatch++;
                    }
                    continue;
                }

                // Try each provider in order
                boolean found = false;
                for (RsidApiProvider provider : providers) {
                    try {
                        RsidApiProvider.ApiResult ar = provider.lookup(snp.chr, snp.pos, snp.ea, snp.nea, build);
                        if (ar != null && ar.rsid != null &&
                                (ar.matchReason.contains("matched_forward") || ar.matchReason.contains("matched_reverse"))) {
                            result.apiMatched++;
                            result.recovered.put(snp.chr + ":" + snp.pos, ar.rsid);
                            cache.put(build, snp.chr, snp.pos, snp.ea, snp.nea, ar.rsid);
                            found = true;
                            System.out.printf("[RsidApi] %s:%d %s/%s → %s via %s (%s)%n",
                                snp.chr, snp.pos, snp.ea, snp.nea, ar.rsid, provider.name(), ar.matchReason);
                            break;
                        }
                    } catch (Exception e) {
                        System.err.printf("[RsidApi] %s error for %s:%d: %s%n",
                            provider.name(), snp.chr, snp.pos, e.getMessage());
                    }
                }
                if (!found) {
                    result.apiNoMatch++;
                    cache.put(build, snp.chr, snp.pos, snp.ea, snp.nea, ""); // cache negative result
                }
            }

            System.out.printf("[RsidApi] Batch %d/%d: %d matched, %d no match%n",
                b + 1, totalBatches, result.apiMatched, result.apiNoMatch);
        }

        cache.save();
        System.out.printf("[RsidApi] API recovered %d / %d · still unmatched: %d%n",
            result.apiMatched, result.attempted, result.apiNoMatch);
        return result;
    }
}
