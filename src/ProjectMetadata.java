import java.io.*;
import java.nio.file.*;
import java.security.*;
import java.time.*;
import java.time.format.*;
import java.util.*;

/**
 * Manages project.json: fingerprint computation, staleness detection,
 * and metadata persistence for the multi-project pipeline cache.
 */
public class ProjectMetadata {

    public String id                    = "";
    public String name                  = "";
    public String description           = "";
    public int    lociCount             = 0;
    public int    totalSnps             = 0;
    public int    snpAnnotationSources  = 0;
    public int    locusAnnotationSources = 0;
    public int    snpsWithAnyAnnotation = 0;
    public String lastProcessed         = "";
    public String pipelineVersion       = "";
    public String coreInputFingerprint  = "";
    public String annotationFingerprint = "";

    // rsID recovery fields
    public boolean rsidColumnPresent     = false;
    public String  rsidRecoveryStatus    = "not_started";  // not_started | completed | failed
    public String  rsidRecoveryRate      = "";
    public String  rsidRecoveryDate      = "";
    public String  selectedSnpDatabase   = "";

    // ── Staleness detection ──────────────────────────────────────────────

    public enum StaleReason { NOT_STALE, NO_METADATA, VERSION_CHANGED, CORE_INPUT_CHANGED, ANNOTATION_CHANGED }

    /**
     * Checks whether a project needs reprocessing.
     * Returns NOT_STALE if up to date, or the first reason found (cheap checks first).
     */
    public static StaleReason checkStaleness(String projectDir, Config config) {
        File metaFile = new File(projectDir, "project.json");
        if (!metaFile.exists()) return StaleReason.NO_METADATA;

        ProjectMetadata stored = load(projectDir);
        if (stored == null) return StaleReason.NO_METADATA;

        if (!Config.PIPELINE_VERSION.equals(stored.pipelineVersion))
            return StaleReason.VERSION_CHANGED;

        String currentCore = computeCoreInputFingerprint(config, projectDir);
        if (!currentCore.equals(stored.coreInputFingerprint))
            return StaleReason.CORE_INPUT_CHANGED;

        String annotPath = new File(projectDir, "annotations.yaml").getAbsolutePath();
        String currentAnnot = computeAnnotationFingerprint(annotPath);
        if (!currentAnnot.equals(stored.annotationFingerprint))
            return StaleReason.ANNOTATION_CHANGED;

        return StaleReason.NOT_STALE;
    }

    public boolean isStale() {
        return !Config.PIPELINE_VERSION.equals(pipelineVersion);
    }

    // ── Core input fingerprint ───────────────────────────────────────────

    /**
     * Covers: GWAS file, loci file, GFF3 file, config.properties itself,
     * and reference panel files (size+mtime only for .bed/.bim/.fam).
     */
    public static String computeCoreInputFingerprint(Config config, String projectDir) {
        Map<String, String> cache = loadFingerprintCache(projectDir);
        boolean cacheChanged = false;

        try {
            MessageDigest digest = MessageDigest.getInstance("SHA-256");

            // Full SHA-256 for these files
            String[] fullHashFiles = {
                config.gwasFile,
                config.lociFile,
                config.gff3File,
                new File(projectDir, "config.properties").getAbsolutePath()
            };

            for (String path : fullHashFiles) {
                File f = new File(path);
                if (!f.exists()) {
                    digest.update(("missing:" + path).getBytes());
                    continue;
                }
                String newHash = cachedFileHash(f, cache);
                if (newHash != null) cacheChanged = true;
                // cachedFileHash always leaves the hash in the cache map
                String hash = cache.getOrDefault(cacheKey(f), "");
                digest.update(hash.getBytes());
            }

            // Reference panel: size+mtime fingerprint (files can be multi-GB)
            if (!config.refPanelPath.isEmpty()) {
                String[] exts = {".bed", ".bim", ".fam"};
                for (String ext : exts) {
                    File f = new File(config.refPanelPath + ext);
                    if (f.exists()) {
                        String sizeTime = f.length() + ":" + f.lastModified();
                        digest.update(sizeTime.getBytes());
                    } else {
                        digest.update(("missing:" + config.refPanelPath + ext).getBytes());
                    }
                }
            }

            if (cacheChanged) saveFingerprintCache(projectDir, cache);
            return "sha256:" + bytesToHex(digest.digest());

        } catch (NoSuchAlgorithmException e) {
            throw new RuntimeException("SHA-256 not available", e);
        }
    }

    // ── Annotation fingerprint ───────────────────────────────────────────

    /**
     * Covers: annotations.yaml content + every file referenced by a `file:` entry.
     */
    public static String computeAnnotationFingerprint(String annotationsYamlPath) {
        File yamlFile = new File(annotationsYamlPath);
        if (!yamlFile.exists()) return "sha256:none";

        try {
            MessageDigest digest = MessageDigest.getInstance("SHA-256");

            // Hash the YAML itself
            digest.update(sha256File(yamlFile).getBytes());

            // Extract and hash all file: references
            List<String> referencedFiles = extractFileReferences(yamlFile);
            File yamlDir = yamlFile.getParentFile();

            for (String ref : referencedFiles) {
                File refFile = new File(ref);
                if (!refFile.isAbsolute()) refFile = new File(yamlDir, ref);
                if (refFile.exists()) {
                    digest.update(sha256File(refFile).getBytes());
                } else {
                    digest.update(("missing:" + ref).getBytes());
                }
            }

            return "sha256:" + bytesToHex(digest.digest());

        } catch (NoSuchAlgorithmException e) {
            throw new RuntimeException("SHA-256 not available", e);
        }
    }

    /**
     * Minimal YAML parser: extracts values of `file:` keys from annotations.yaml.
     */
    static List<String> extractFileReferences(File yamlFile) {
        List<String> files = new ArrayList<>();
        try (BufferedReader br = new BufferedReader(new FileReader(yamlFile))) {
            String line;
            while ((line = br.readLine()) != null) {
                String trimmed = line.trim();
                if (trimmed.startsWith("file:")) {
                    String value = trimmed.substring(5).trim();
                    // Strip quotes
                    if (value.startsWith("\"") && value.endsWith("\""))
                        value = value.substring(1, value.length() - 1);
                    else if (value.startsWith("'") && value.endsWith("'"))
                        value = value.substring(1, value.length() - 1);
                    if (!value.isEmpty()) files.add(value);
                }
            }
        } catch (IOException e) {
            System.err.println("[ProjectMetadata] Could not read " + yamlFile + ": " + e.getMessage());
        }
        return files;
    }

    // ── Fingerprint cache (size+mtime → skip rehashing unchanged files) ──

    /**
     * Returns the new hash if the file changed since last cache entry, or null if cache is valid.
     */
    private static String cachedFileHash(File f, Map<String, String> cache) {
        String key = cacheKey(f);
        String cachedSize  = cache.get(key + ".size");
        String cachedMtime = cache.get(key + ".mtime");

        if (cachedSize != null && cachedMtime != null
                && cachedSize.equals(String.valueOf(f.length()))
                && cachedMtime.equals(String.valueOf(f.lastModified()))
                && cache.containsKey(key)) {
            return null; // cache hit — no change
        }

        // File changed or not cached yet — recompute
        String hash = sha256File(f);
        cache.put(key, hash);
        cache.put(key + ".size", String.valueOf(f.length()));
        cache.put(key + ".mtime", String.valueOf(f.lastModified()));
        return hash;
    }

    private static String cacheKey(File f) {
        try { return f.getCanonicalPath(); }
        catch (IOException e) { return f.getAbsolutePath(); }
    }

    private static Map<String, String> loadFingerprintCache(String projectDir) {
        Map<String, String> cache = new LinkedHashMap<>();
        File cacheFile = new File(projectDir, ".fingerprint_cache");
        if (!cacheFile.exists()) return cache;
        try (BufferedReader br = new BufferedReader(new FileReader(cacheFile))) {
            String line;
            while ((line = br.readLine()) != null) {
                int eq = line.indexOf('=');
                if (eq > 0) cache.put(line.substring(0, eq), line.substring(eq + 1));
            }
        } catch (IOException e) { /* fresh cache */ }
        return cache;
    }

    private static void saveFingerprintCache(String projectDir, Map<String, String> cache) {
        File cacheFile = new File(projectDir, ".fingerprint_cache");
        try (PrintWriter pw = new PrintWriter(new FileWriter(cacheFile))) {
            for (Map.Entry<String, String> e : cache.entrySet())
                pw.println(e.getKey() + "=" + e.getValue());
        } catch (IOException e) {
            System.err.println("[ProjectMetadata] Could not save fingerprint cache: " + e.getMessage());
        }
    }

    // ── SHA-256 helper ───────────────────────────────────────────────────

    static String sha256File(File f) {
        try {
            MessageDigest md = MessageDigest.getInstance("SHA-256");
            try (InputStream is = new BufferedInputStream(new FileInputStream(f), 65536)) {
                byte[] buf = new byte[65536];
                int n;
                while ((n = is.read(buf)) != -1) md.update(buf, 0, n);
            }
            return bytesToHex(md.digest());
        } catch (Exception e) {
            return "error:" + e.getMessage();
        }
    }

    private static String bytesToHex(byte[] bytes) {
        StringBuilder sb = new StringBuilder(bytes.length * 2);
        for (byte b : bytes) sb.append(String.format("%02x", b & 0xff));
        return sb.toString();
    }

    // ── Read/write project.json ──────────────────────────────────────────

    public static ProjectMetadata load(String projectDir) {
        File f = new File(projectDir, "project.json");
        if (!f.exists()) return null;
        try {
            String json = new String(Files.readAllBytes(f.toPath()), "UTF-8");
            ProjectMetadata pm = new ProjectMetadata();
            pm.id                      = extractJsonStr(json, "id");
            pm.name                    = extractJsonStr(json, "name");
            pm.description             = extractJsonStr(json, "description");
            pm.lociCount               = extractJsonInt(json, "loci_count");
            pm.totalSnps               = extractJsonInt(json, "total_snps");
            pm.snpAnnotationSources    = extractJsonInt(json, "snp_annotation_sources");
            pm.locusAnnotationSources  = extractJsonInt(json, "locus_annotation_sources");
            pm.snpsWithAnyAnnotation   = extractJsonInt(json, "snps_with_any_annotation");
            pm.lastProcessed           = extractJsonStr(json, "last_processed");
            pm.pipelineVersion         = extractJsonStr(json, "pipeline_version");
            pm.coreInputFingerprint    = extractJsonStr(json, "core_input_fingerprint");
            pm.annotationFingerprint   = extractJsonStr(json, "annotation_fingerprint");
            pm.rsidColumnPresent       = "true".equals(extractJsonStr(json, "rsid_column_present"));
            pm.rsidRecoveryStatus      = extractJsonStr(json, "rsid_recovery_status");
            if (pm.rsidRecoveryStatus.isEmpty()) pm.rsidRecoveryStatus = "not_started";
            pm.rsidRecoveryRate        = extractJsonStr(json, "rsid_recovery_rate");
            pm.rsidRecoveryDate        = extractJsonStr(json, "rsid_recovery_date");
            pm.selectedSnpDatabase     = extractJsonStr(json, "selected_snp_database");
            return pm;
        } catch (IOException e) {
            System.err.println("[ProjectMetadata] Could not read " + f + ": " + e.getMessage());
            return null;
        }
    }

    public void save(String projectDir) throws IOException {
        this.lastProcessed = Instant.now().atOffset(ZoneOffset.UTC)
            .format(DateTimeFormatter.ISO_OFFSET_DATE_TIME);
        this.pipelineVersion = Config.PIPELINE_VERSION;

        File f = new File(projectDir, "project.json");
        try (PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(f)))) {
            pw.println("{");
            pw.printf("  \"id\": \"%s\",%n", escJson(id));
            pw.printf("  \"name\": \"%s\",%n", escJson(name));
            pw.printf("  \"description\": \"%s\",%n", escJson(description));
            pw.printf("  \"loci_count\": %d,%n", lociCount);
            pw.printf("  \"total_snps\": %d,%n", totalSnps);
            pw.printf("  \"snp_annotation_sources\": %d,%n", snpAnnotationSources);
            pw.printf("  \"locus_annotation_sources\": %d,%n", locusAnnotationSources);
            pw.printf("  \"snps_with_any_annotation\": %d,%n", snpsWithAnyAnnotation);
            pw.printf("  \"last_processed\": \"%s\",%n", escJson(lastProcessed));
            pw.printf("  \"pipeline_version\": \"%s\",%n", escJson(pipelineVersion));
            pw.printf("  \"core_input_fingerprint\": \"%s\",%n", escJson(coreInputFingerprint));
            pw.printf("  \"annotation_fingerprint\": \"%s\",%n", escJson(annotationFingerprint));
            pw.printf("  \"rsid_column_present\": %s,%n", rsidColumnPresent);
            pw.printf("  \"rsid_recovery_status\": \"%s\",%n", escJson(rsidRecoveryStatus));
            pw.printf("  \"rsid_recovery_rate\": \"%s\",%n", escJson(rsidRecoveryRate));
            pw.printf("  \"rsid_recovery_date\": \"%s\",%n", escJson(rsidRecoveryDate));
            pw.printf("  \"selected_snp_database\": \"%s\"%n", escJson(selectedSnpDatabase));
            pw.println("}");
        }
    }

    // ── Annotation counting ──────────────────────────────────────────────

    /**
     * Counts annotation sources by level from annotations.yaml.
     * Sets snpAnnotationSources and locusAnnotationSources.
     */
    public void countAnnotationSources(String annotationsYamlPath) {
        File yamlFile = new File(annotationsYamlPath);
        if (!yamlFile.exists()) return;
        try (BufferedReader br = new BufferedReader(new FileReader(yamlFile))) {
            String line;
            while ((line = br.readLine()) != null) {
                String trimmed = line.trim();
                if (trimmed.startsWith("level:")) {
                    String level = trimmed.substring(6).trim()
                        .replace("\"", "").replace("'", "");
                    if ("snp".equalsIgnoreCase(level)) snpAnnotationSources++;
                    else if ("locus".equalsIgnoreCase(level)) locusAnnotationSources++;
                }
            }
        } catch (IOException e) {
            System.err.println("[ProjectMetadata] Could not count annotations: " + e.getMessage());
        }
    }

    // ── Minimal JSON helpers (no external dependencies) ──────────────────

    private static String extractJsonStr(String json, String key) {
        String search = "\"" + key + "\":";
        int i = json.indexOf(search);
        if (i < 0) return "";
        i += search.length();
        while (i < json.length() && json.charAt(i) == ' ') i++;
        if (i >= json.length()) return "";
        // Handle unquoted values (booleans, numbers, null)
        if (json.charAt(i) != '"') {
            int end = i;
            while (end < json.length() && ",}\n\r\t ".indexOf(json.charAt(end)) < 0) end++;
            return json.substring(i, end).trim();
        }
        i++;
        StringBuilder sb = new StringBuilder();
        while (i < json.length()) {
            char c = json.charAt(i++);
            if (c == '"') break;
            if (c == '\\' && i < json.length()) { sb.append(json.charAt(i++)); continue; }
            sb.append(c);
        }
        return sb.toString();
    }

    private static int extractJsonInt(String json, String key) {
        String search = "\"" + key + "\":";
        int i = json.indexOf(search);
        if (i < 0) return 0;
        i += search.length();
        while (i < json.length() && json.charAt(i) == ' ') i++;
        int end = i;
        while (end < json.length() && "0123456789-".indexOf(json.charAt(end)) >= 0) end++;
        if (end == i) return 0;
        try { return Integer.parseInt(json.substring(i, end)); }
        catch (NumberFormatException e) { return 0; }
    }

    private static String escJson(String s) {
        if (s == null) return "";
        return s.replace("\\", "\\\\").replace("\"", "\\\"")
                .replace("\n", "\\n").replace("\r", "\\r").replace("\t", "\\t");
    }
}
