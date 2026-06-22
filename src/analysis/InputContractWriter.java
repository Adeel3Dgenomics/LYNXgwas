
import java.io.*;
import java.util.*;

/**
 * Writes the versioned input manifest that every analysis tool reads.
 * Contains: paths to harmonized data + LD matrix, sample sizes,
 * ancestry, genome build, trait type, and locus metadata.
 *
 * Schema version is bumped when the contract changes in a breaking way.
 */
public class InputContractWriter {

    public static final String SCHEMA_VERSION = "1.0";

    public static class DatasetMeta {
        public String projectId = "";
        public String locusId = "";
        public int locusIndex;
        public String chr = "";
        public long start, end;
        public String genomeBuild = "GRCh37";
        public String ancestry = "";
        public String traitType = "";        // "quantitative" or "binary"
        public String effectType = "";       // "beta", "OR", "logOR"
        public int sampleN;
        public int nCases;
        public int nControls;
    }

    public static void write(File runDir, File harmonizedDir, File ldDir,
                              DatasetMeta meta) throws IOException {
        runDir.mkdirs();
        File outFile = new File(runDir, "run_input.manifest.json");

        try (PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(outFile)))) {
            pw.println("{");
            pw.printf("  \"schema_version\": \"%s\",%n", SCHEMA_VERSION);
            pw.printf("  \"locus_id\": \"%s\",%n", esc(meta.locusId));
            pw.printf("  \"locus_index\": %d,%n", meta.locusIndex);
            pw.printf("  \"chr\": \"%s\",%n", esc(meta.chr));
            pw.printf("  \"start\": %d,%n", meta.start);
            pw.printf("  \"end\": %d,%n", meta.end);

            // Dataset metadata
            pw.printf("  \"project_id\": \"%s\",%n", esc(meta.projectId));
            pw.printf("  \"genome_build\": \"%s\",%n", esc(meta.genomeBuild));
            pw.printf("  \"ancestry\": \"%s\",%n", esc(meta.ancestry));
            pw.printf("  \"trait_type\": \"%s\",%n", esc(meta.traitType));
            pw.printf("  \"effect_type\": \"%s\",%n", esc(meta.effectType));
            pw.printf("  \"sample_n\": %d,%n", meta.sampleN);
            pw.printf("  \"n_cases\": %d,%n", meta.nCases);
            pw.printf("  \"n_controls\": %d,%n", meta.nControls);

            // Paths to shared artifacts
            pw.printf("  \"harmonized_gwas\": \"%s\",%n",
                esc(new File(harmonizedDir, "harmonized_gwas.tsv").getAbsolutePath()));
            pw.printf("  \"allele_alignment\": \"%s\",%n",
                esc(new File(harmonizedDir, "allele_alignment.tsv").getAbsolutePath()));
            pw.printf("  \"ld_r_matrix\": \"%s\",%n",
                esc(new File(ldDir, "ld_r.matrix").getAbsolutePath()));
            pw.printf("  \"ld_r2_matrix\": \"%s\",%n",
                esc(new File(ldDir, "ld_r2.matrix").getAbsolutePath()));
            pw.printf("  \"ld_dprime_matrix\": \"%s\",%n",
                esc(new File(ldDir, "ld_dprime.matrix").getAbsolutePath()));
            pw.printf("  \"ld_snp_order\": \"%s\",%n",
                esc(new File(ldDir, "ld_snp_order.txt").getAbsolutePath()));

            // Consistency diagnostic
            File diagFile = new File(ldDir, "consistency_summary.json");
            if (diagFile.exists()) {
                String diagJson = new String(java.nio.file.Files.readAllBytes(diagFile.toPath()), "UTF-8");
                pw.printf("  \"ld_diagnostic\": %s,%n", diagJson.trim());
            }

            pw.printf("  \"created_at\": %d%n", System.currentTimeMillis());
            pw.println("}");
        }
    }

    /**
     * Read and validate an input manifest. Returns null + prints error if invalid.
     */
    public static Map<String, String> readAndValidate(File manifestFile) throws IOException {
        if (!manifestFile.exists()) return null;
        String json = new String(java.nio.file.Files.readAllBytes(manifestFile.toPath()), "UTF-8");

        Map<String, String> fields = new LinkedHashMap<>();
        String[] required = {"schema_version", "locus_id", "harmonized_gwas",
                             "ld_r2_matrix", "ld_snp_order"};

        for (String key : required) {
            String val = extractField(json, key);
            if (val == null || val.isEmpty()) {
                System.err.printf("[InputContract] Missing required field: %s%n", key);
                return null;
            }
            fields.put(key, val);
        }

        // Validate schema version
        String ver = fields.get("schema_version");
        if (!SCHEMA_VERSION.equals(ver)) {
            System.err.printf("[InputContract] Schema version mismatch: expected %s, got %s%n",
                SCHEMA_VERSION, ver);
            return null;
        }

        // Validate files exist
        for (String pathKey : new String[]{"harmonized_gwas", "ld_r2_matrix", "ld_snp_order"}) {
            String path = fields.get(pathKey);
            if (path != null && !new File(path).exists()) {
                System.err.printf("[InputContract] Referenced file not found: %s = %s%n", pathKey, path);
                return null;
            }
        }

        // Add optional fields
        String[] optional = {"project_id", "chr", "genome_build", "ancestry",
                             "trait_type", "effect_type", "sample_n", "n_cases", "n_controls",
                             "ld_r_matrix", "ld_dprime_matrix"};
        for (String key : optional) {
            String val = extractField(json, key);
            if (val != null) fields.put(key, val);
        }

        return fields;
    }

    private static String extractField(String json, String key) {
        String marker = "\"" + key + "\":";
        int i = json.indexOf(marker);
        if (i < 0) return null;
        i += marker.length();
        while (i < json.length() && json.charAt(i) == ' ') i++;
        if (i >= json.length()) return null;
        if (json.charAt(i) == '"') {
            int start = ++i;
            while (i < json.length() && json.charAt(i) != '"') {
                if (json.charAt(i) == '\\') i++;
                i++;
            }
            return json.substring(start, i);
        }
        if (json.charAt(i) == '{') return null; // nested object, skip
        int end = i;
        while (end < json.length() && ",}\n\r".indexOf(json.charAt(end)) < 0) end++;
        return json.substring(i, end).trim();
    }

    private static String esc(String s) {
        if (s == null) return "";
        return s.replace("\\", "\\\\").replace("\"", "\\\"");
    }
}
