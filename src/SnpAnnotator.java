import java.io.*;
import java.net.*;
import java.util.*;

/**
 * Optional: resolves chr:pos:ref:alt varids to rsids via NCBI E-utilities.
 * Only used when the GWAS file has no rsid column and no top_snp_file is given.
 * Results are cached in memory.
 */
public class SnpAnnotator {

    private final Map<String, String> cache = new HashMap<>();

    /**
     * For each locus, tries to look up rsids for the top SNP only (to minimise API calls).
     * Bulk annotation of all SNPs is intentionally skipped — the display ID falls back to varid.
     */
    public void annotateTopSnps(Map<Integer, Snp> topSnps, Config config) {
        if (!config.colRsid.isEmpty()) return; // rsid column is already present in GWAS data

        System.out.println("[SnpAnnotator] Attempting rsid lookup for top SNPs via NCBI...");
        int resolved = 0;
        for (Snp snp : topSnps.values()) {
            if (snp == null || !snp.id.contains(":")) continue; // already an rsid
            String rsid = lookup(snp.chr, snp.pos, snp.ea, snp.nea);
            if (rsid != null) {
                snp.id = rsid;
                resolved++;
            }
        }
        System.out.printf("[SnpAnnotator] Resolved %d / %d top SNPs to rsids%n",
            resolved, topSnps.size());
    }

    private String lookup(String chr, long pos, String ea, String nea) {
        String key = chr + ":" + pos + ":" + nea + ":" + ea;
        if (cache.containsKey(key)) return cache.get(key);

        // NCBI E-utilities: search by chr, pos, alleles
        // Uses the variant search endpoint
        String term = "\"" + chr + "[Chromosome]\" AND "
                    + pos + "[Base Position] AND "
                    + "\"Homo sapiens\"[Organism]";
        try {
            String url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
                       + "?db=snp&term=" + URLEncoder.encode(term, "UTF-8")
                       + "&retmax=5&retmode=json";
            String json = httpGet(url);
            String rsid = extractFirstRsid(json);
            cache.put(key, rsid);
            Thread.sleep(350); // respect NCBI rate limit (~3 req/s)
            return rsid;
        } catch (Exception e) {
            System.err.printf("[WARN] dbSNP lookup failed for %s: %s%n", key, e.getMessage());
            cache.put(key, null);
            return null;
        }
    }

    private static String httpGet(String url) throws Exception {
        HttpURLConnection conn = (HttpURLConnection) new URL(url).openConnection();
        conn.setRequestProperty("User-Agent", "LociCatalogue/1.0 (sammanmahmoud@gmail.com)");
        conn.setConnectTimeout(10_000);
        conn.setReadTimeout(15_000);
        try (BufferedReader br = new BufferedReader(new InputStreamReader(conn.getInputStream()))) {
            StringBuilder sb = new StringBuilder();
            String line;
            while ((line = br.readLine()) != null) sb.append(line);
            return sb.toString();
        }
    }

    private static String extractFirstRsid(String json) {
        // Simple extraction from NCBI JSON: "idlist":["123456789"]
        int idx = json.indexOf("\"idlist\"");
        if (idx < 0) return null;
        int bracketOpen = json.indexOf('[', idx);
        int bracketClose = json.indexOf(']', bracketOpen);
        if (bracketOpen < 0 || bracketClose < 0) return null;
        String ids = json.substring(bracketOpen + 1, bracketClose).trim();
        if (ids.isEmpty()) return null;
        // First id (may be quoted)
        String id = ids.replace("\"", "").split(",")[0].trim();
        return id.isEmpty() ? null : "rs" + id;
    }
}
