package rsid;

import java.io.*;
import java.nio.file.*;
import java.util.*;

/**
 * Manages config/global.json: reference panels, SNP databases, NCBI API settings.
 * Shared across all projects — resources are referenced by id.
 */
public class GlobalConfig {

    public List<RefPanel> referencePanels = new ArrayList<>();
    public List<SnpDatabase> snpDatabases = new ArrayList<>();
    public NcbiApiConfig ncbiApi = new NcbiApiConfig();
    public SharedStorage sharedStorage = new SharedStorage();
    public HpcConfig hpc = new HpcConfig();

    public static class RefPanel {
        public String id = "";
        public String label = "";
        public String population = "";
        public String build = "";
        public String plinkPath = "";
        public List<String> validationErrors = new ArrayList<>();

        public void validate() {
            validationErrors.clear();
            if (plinkPath.isEmpty()) { validationErrors.add("PLINK path is empty"); return; }
            for (String ext : new String[]{".bed", ".bim", ".fam"}) {
                if (!new File(plinkPath + ext).exists())
                    validationErrors.add("Missing: " + plinkPath + ext);
            }
        }

        public boolean isValid() { return validationErrors.isEmpty(); }
    }

    public static class SnpDatabase {
        public String id = "";
        public String label = "";
        public String build = "";
        public String folder = "";
        public String filePattern = "{chr}.vcf.gz";
        public List<String> presentChromosomes = new ArrayList<>();
        public List<String> missingChromosomes = new ArrayList<>();
        public List<String> missingIndex = new ArrayList<>();

        public void validate() {
            presentChromosomes.clear();
            missingChromosomes.clear();
            missingIndex.clear();
            if (folder.isEmpty()) return;
            String[] chroms = {"1","2","3","4","5","6","7","8","9","10","11","12",
                "13","14","15","16","17","18","19","20","21","22","X","Y","MT"};
            for (String chr : chroms) {
                String vcf = filePattern.replace("{chr}", chr);
                File vcfFile = new File(folder, vcf);
                File tbiFile = new File(folder, vcf + ".tbi");
                if (vcfFile.exists()) {
                    presentChromosomes.add(chr);
                    if (!tbiFile.exists()) missingIndex.add(chr);
                } else {
                    missingChromosomes.add(chr);
                }
            }
        }

        public boolean isValid() {
            return !presentChromosomes.isEmpty() && missingIndex.isEmpty();
        }
    }

    public static class NcbiApiConfig {
        public boolean enabled = false;
        public String baseUrl = "https://api.ncbi.nlm.nih.gov/variation/v0/";
        public String apiKey = "";
        public int rateLimitPerSec = 3;
        public boolean restrictToLeadSnps = true;
    }

    /** Alternative to per-user download for large reference data (PLINK panels, GFF3
     *  annotation): a shared local/network folder (also covers a Google Drive/Dropbox desktop
     *  sync folder — those just look like an ordinary path once synced, no OAuth needed), or a
     *  git repo the user already has their own credentials configured for. This app never
     *  handles a token or password itself — "git" mode only ever shells out to the user's own
     *  installed git binary. */
    public static class SharedStorage {
        public String mode = "none"; // "none" | "folder" | "git"
        public String folderPath = "";
        public String gitUrl = "";
        public String gitLocalClone = "";
    }

    /** Optional, opt-in remote execution on an HPC cluster via the user's own SSH access. Disabled
     *  by default — local execution is completely unaffected either way. sshKeyPath is a path to
     *  an existing private key file; there is deliberately no password field anywhere in this
     *  config, key-based auth only. */
    public static class HpcConfig {
        public boolean enabled = false;
        public String sshHost = "";
        public String sshUser = "";
        public String sshKeyPath = "";
        public String remoteWorkDir = "";
        public Map<String, String> moduleNames = new LinkedHashMap<>();
        public int pollIntervalMinMinutes = 5;
        public int pollIntervalMaxMinutes = 20;
    }

    // ── Persistence ──────────────────────────────────────────────

    private static final String CONFIG_PATH = "config/global.json";

    public static GlobalConfig load() {
        File f = new File(CONFIG_PATH);
        if (!f.exists()) return createDefault();
        try {
            String json = new String(Files.readAllBytes(f.toPath()), "UTF-8");
            return parse(json);
        } catch (Exception e) {
            System.err.println("[GlobalConfig] Error loading: " + e.getMessage());
            return createDefault();
        }
    }

    public void save() throws IOException {
        new File("config").mkdirs();
        try (PrintWriter pw = new PrintWriter(new BufferedWriter(
                new FileWriter(CONFIG_PATH)))) {
            pw.println("{");
            pw.println("  \"reference_panels\": [");
            for (int i = 0; i < referencePanels.size(); i++) {
                RefPanel rp = referencePanels.get(i);
                pw.printf("    {%n      \"id\": \"%s\",%n      \"label\": \"%s\",%n      \"population\": \"%s\",%n      \"build\": \"%s\",%n      \"plink_path\": \"%s\"%n    }",
                    esc(rp.id), esc(rp.label), esc(rp.population), esc(rp.build), esc(rp.plinkPath));
                if (i < referencePanels.size() - 1) pw.println(",");
                else pw.println();
            }
            pw.println("  ],");
            pw.println("  \"snp_databases\": [");
            for (int i = 0; i < snpDatabases.size(); i++) {
                SnpDatabase sd = snpDatabases.get(i);
                pw.printf("    {%n      \"id\": \"%s\",%n      \"label\": \"%s\",%n      \"build\": \"%s\",%n      \"folder\": \"%s\",%n      \"file_pattern\": \"%s\"%n    }",
                    esc(sd.id), esc(sd.label), esc(sd.build), esc(sd.folder), esc(sd.filePattern));
                if (i < snpDatabases.size() - 1) pw.println(",");
                else pw.println();
            }
            pw.println("  ],");
            pw.println("  \"ncbi_api\": {");
            pw.printf("    \"enabled\": %s,%n", ncbiApi.enabled);
            pw.printf("    \"base_url\": \"%s\",%n", esc(ncbiApi.baseUrl));
            pw.printf("    \"api_key\": \"%s\",%n", esc(ncbiApi.apiKey));
            pw.printf("    \"rate_limit_per_sec\": %d,%n", ncbiApi.rateLimitPerSec);
            pw.printf("    \"restrict_to_lead_snps\": %s%n", ncbiApi.restrictToLeadSnps);
            pw.println("  }");
            pw.println("}");
        }
    }

    public void validateAll() {
        for (RefPanel rp : referencePanels) rp.validate();
        for (SnpDatabase sd : snpDatabases) sd.validate();
    }

    public RefPanel findPanel(String id) {
        for (RefPanel rp : referencePanels) if (rp.id.equals(id)) return rp;
        return null;
    }

    public SnpDatabase findDatabase(String id) {
        for (SnpDatabase sd : snpDatabases) if (sd.id.equals(id)) return sd;
        return null;
    }

    public String toJson() {
        StringBuilder j = new StringBuilder();
        j.append("{\"reference_panels\":[");
        for (int i = 0; i < referencePanels.size(); i++) {
            if (i > 0) j.append(',');
            RefPanel rp = referencePanels.get(i);
            rp.validate();
            j.append("{\"id\":\"").append(esc(rp.id))
             .append("\",\"label\":\"").append(esc(rp.label))
             .append("\",\"population\":\"").append(esc(rp.population))
             .append("\",\"build\":\"").append(esc(rp.build))
             .append("\",\"plink_path\":\"").append(esc(rp.plinkPath))
             .append("\",\"valid\":").append(rp.isValid())
             .append(",\"errors\":[");
            for (int e = 0; e < rp.validationErrors.size(); e++) {
                if (e > 0) j.append(',');
                j.append('"').append(esc(rp.validationErrors.get(e))).append('"');
            }
            j.append("]}");
        }
        j.append("],\"snp_databases\":[");
        for (int i = 0; i < snpDatabases.size(); i++) {
            if (i > 0) j.append(',');
            SnpDatabase sd = snpDatabases.get(i);
            sd.validate();
            j.append("{\"id\":\"").append(esc(sd.id))
             .append("\",\"label\":\"").append(esc(sd.label))
             .append("\",\"build\":\"").append(esc(sd.build))
             .append("\",\"folder\":\"").append(esc(sd.folder))
             .append("\",\"file_pattern\":\"").append(esc(sd.filePattern))
             .append("\",\"valid\":").append(sd.isValid())
             .append(",\"present_chromosomes\":").append(sd.presentChromosomes.size())
             .append(",\"missing_index\":[");
            for (int e = 0; e < sd.missingIndex.size(); e++) {
                if (e > 0) j.append(',');
                j.append('"').append(sd.missingIndex.get(e)).append('"');
            }
            j.append("]}");
        }
        j.append("],\"ncbi_api\":{");
        j.append("\"enabled\":").append(ncbiApi.enabled);
        j.append(",\"base_url\":\"").append(esc(ncbiApi.baseUrl)).append('"');
        j.append(",\"api_key\":\"").append(esc(ncbiApi.apiKey)).append('"');
        j.append(",\"rate_limit_per_sec\":").append(ncbiApi.rateLimitPerSec);
        j.append(",\"restrict_to_lead_snps\":").append(ncbiApi.restrictToLeadSnps);
        j.append("},\"shared_storage\":{");
        j.append("\"mode\":\"").append(esc(sharedStorage.mode)).append('"');
        j.append(",\"folder_path\":\"").append(esc(sharedStorage.folderPath)).append('"');
        j.append(",\"git_url\":\"").append(esc(sharedStorage.gitUrl)).append('"');
        j.append(",\"git_local_clone\":\"").append(esc(sharedStorage.gitLocalClone)).append('"');
        j.append("},\"hpc\":{");
        j.append("\"enabled\":").append(hpc.enabled);
        j.append(",\"ssh_host\":\"").append(esc(hpc.sshHost)).append('"');
        j.append(",\"ssh_user\":\"").append(esc(hpc.sshUser)).append('"');
        j.append(",\"ssh_key_path\":\"").append(esc(hpc.sshKeyPath)).append('"');
        j.append(",\"remote_work_dir\":\"").append(esc(hpc.remoteWorkDir)).append('"');
        j.append(",\"poll_interval_min_minutes\":").append(hpc.pollIntervalMinMinutes);
        j.append(",\"poll_interval_max_minutes\":").append(hpc.pollIntervalMaxMinutes);
        j.append(",\"module_names\":{");
        int mi = 0;
        for (Map.Entry<String, String> en : hpc.moduleNames.entrySet()) {
            if (mi++ > 0) j.append(',');
            j.append('"').append(esc(en.getKey())).append("\":\"").append(esc(en.getValue())).append('"');
        }
        j.append("}}}");
        return j.toString();
    }

    // ── Parsing ──────────────────────────────────────────────────

    private static GlobalConfig createDefault() {
        return new GlobalConfig();
    }

    /** Test-only entry point into the otherwise-private parse() — this project has no test
     *  framework to grant package-private access across the default/rsid package boundary. */
    public static GlobalConfig parseForTest(String json) {
        return parse(json);
    }

    private static GlobalConfig parse(String json) {
        GlobalConfig gc = new GlobalConfig();
        // Parse reference_panels array
        int rpStart = json.indexOf("\"reference_panels\"");
        if (rpStart >= 0) {
            int arrStart = json.indexOf('[', rpStart);
            int arrEnd = findMatchingBracket(json, arrStart);
            if (arrStart >= 0 && arrEnd > arrStart) {
                String arr = json.substring(arrStart + 1, arrEnd);
                gc.referencePanels = parseRefPanels(arr);
            }
        }
        // Parse snp_databases array
        int sdStart = json.indexOf("\"snp_databases\"");
        if (sdStart >= 0) {
            int arrStart = json.indexOf('[', sdStart);
            int arrEnd = findMatchingBracket(json, arrStart);
            if (arrStart >= 0 && arrEnd > arrStart) {
                String arr = json.substring(arrStart + 1, arrEnd);
                gc.snpDatabases = parseSnpDatabases(arr);
            }
        }
        // Parse ncbi_api
        int ncbiStart = json.indexOf("\"ncbi_api\"");
        if (ncbiStart >= 0) {
            int objStart = json.indexOf('{', ncbiStart);
            int objEnd = findMatchingBrace(json, objStart);
            if (objStart >= 0 && objEnd > objStart) {
                String obj = json.substring(objStart, objEnd + 1);
                gc.ncbiApi.enabled = "true".equals(jsonStr(obj, "enabled"));
                String bu = jsonStr(obj, "base_url");
                if (!bu.isEmpty()) gc.ncbiApi.baseUrl = bu;
                gc.ncbiApi.apiKey = jsonStr(obj, "api_key");
                String rl = jsonStr(obj, "rate_limit_per_sec");
                if (!rl.isEmpty()) gc.ncbiApi.rateLimitPerSec = Integer.parseInt(rl);
                gc.ncbiApi.restrictToLeadSnps = !"false".equals(jsonStr(obj, "restrict_to_lead_snps"));
            }
        }
        // Parse shared_storage
        int ssStart = json.indexOf("\"shared_storage\"");
        if (ssStart >= 0) {
            int objStart = json.indexOf('{', ssStart);
            int objEnd = findMatchingBrace(json, objStart);
            if (objStart >= 0 && objEnd > objStart) {
                String obj = json.substring(objStart, objEnd + 1);
                String mode = jsonStr(obj, "mode");
                if (!mode.isEmpty()) gc.sharedStorage.mode = mode;
                gc.sharedStorage.folderPath = jsonStr(obj, "folder_path");
                gc.sharedStorage.gitUrl = jsonStr(obj, "git_url");
                gc.sharedStorage.gitLocalClone = jsonStr(obj, "git_local_clone");
            }
        }
        // Parse hpc
        int hpcStart = json.indexOf("\"hpc\"");
        if (hpcStart >= 0) {
            int objStart = json.indexOf('{', hpcStart);
            int objEnd = findMatchingBrace(json, objStart);
            if (objStart >= 0 && objEnd > objStart) {
                String obj = json.substring(objStart, objEnd + 1);
                gc.hpc.enabled = "true".equals(jsonStr(obj, "enabled"));
                gc.hpc.sshHost = jsonStr(obj, "ssh_host");
                gc.hpc.sshUser = jsonStr(obj, "ssh_user");
                gc.hpc.sshKeyPath = jsonStr(obj, "ssh_key_path");
                gc.hpc.remoteWorkDir = jsonStr(obj, "remote_work_dir");
                String pMin = jsonStr(obj, "poll_interval_min_minutes");
                if (!pMin.isEmpty()) gc.hpc.pollIntervalMinMinutes = Integer.parseInt(pMin);
                String pMax = jsonStr(obj, "poll_interval_max_minutes");
                if (!pMax.isEmpty()) gc.hpc.pollIntervalMaxMinutes = Integer.parseInt(pMax);
                int mnStart = obj.indexOf("\"module_names\"");
                if (mnStart >= 0) {
                    int mnObjStart = obj.indexOf('{', mnStart);
                    int mnObjEnd = findMatchingBrace(obj, mnObjStart);
                    if (mnObjStart >= 0 && mnObjEnd > mnObjStart) {
                        gc.hpc.moduleNames = parseStringMap(obj.substring(mnObjStart + 1, mnObjEnd));
                    }
                }
            }
        }
        return gc;
    }

    private static Map<String, String> parseStringMap(String inner) {
        Map<String, String> map = new LinkedHashMap<>();
        int i = 0;
        while (i < inner.length()) {
            int keyStart = inner.indexOf('"', i);
            if (keyStart < 0) break;
            int keyEnd = skipString(inner, keyStart);
            String key = inner.substring(keyStart + 1, keyEnd);
            int colon = inner.indexOf(':', keyEnd);
            if (colon < 0) break;
            int valStart = inner.indexOf('"', colon);
            if (valStart < 0) break;
            int valEnd = skipString(inner, valStart);
            String val = inner.substring(valStart + 1, valEnd)
                .replace("\\\"", "\"").replace("\\\\", "\\");
            map.put(key, val);
            i = valEnd + 1;
        }
        return map;
    }

    private static List<RefPanel> parseRefPanels(String arr) {
        List<RefPanel> list = new ArrayList<>();
        int idx = 0;
        while (true) {
            int objStart = arr.indexOf('{', idx);
            if (objStart < 0) break;
            int objEnd = findMatchingBrace(arr, objStart);
            if (objEnd < 0) break;
            String obj = arr.substring(objStart, objEnd + 1);
            RefPanel rp = new RefPanel();
            rp.id = jsonStr(obj, "id");
            rp.label = jsonStr(obj, "label");
            rp.population = jsonStr(obj, "population");
            rp.build = jsonStr(obj, "build");
            rp.plinkPath = jsonStr(obj, "plink_path");
            list.add(rp);
            idx = objEnd + 1;
        }
        return list;
    }

    private static List<SnpDatabase> parseSnpDatabases(String arr) {
        List<SnpDatabase> list = new ArrayList<>();
        int idx = 0;
        while (true) {
            int objStart = arr.indexOf('{', idx);
            if (objStart < 0) break;
            int objEnd = findMatchingBrace(arr, objStart);
            if (objEnd < 0) break;
            String obj = arr.substring(objStart, objEnd + 1);
            SnpDatabase sd = new SnpDatabase();
            sd.id = jsonStr(obj, "id");
            sd.label = jsonStr(obj, "label");
            sd.build = jsonStr(obj, "build");
            sd.folder = jsonStr(obj, "folder");
            String fp = jsonStr(obj, "file_pattern");
            if (!fp.isEmpty()) sd.filePattern = fp;
            list.add(sd);
            idx = objEnd + 1;
        }
        return list;
    }

    private static int findMatchingBracket(String s, int openPos) {
        if (openPos < 0 || s.charAt(openPos) != '[') return -1;
        int depth = 1;
        for (int i = openPos + 1; i < s.length(); i++) {
            char c = s.charAt(i);
            if (c == '[') depth++;
            else if (c == ']') { depth--; if (depth == 0) return i; }
            else if (c == '"') i = skipString(s, i);
        }
        return -1;
    }

    private static int findMatchingBrace(String s, int openPos) {
        if (openPos < 0 || s.charAt(openPos) != '{') return -1;
        int depth = 1;
        for (int i = openPos + 1; i < s.length(); i++) {
            char c = s.charAt(i);
            if (c == '{') depth++;
            else if (c == '}') { depth--; if (depth == 0) return i; }
            else if (c == '"') i = skipString(s, i);
        }
        return -1;
    }

    private static int skipString(String s, int quotePos) {
        for (int i = quotePos + 1; i < s.length(); i++) {
            if (s.charAt(i) == '\\') i++;
            else if (s.charAt(i) == '"') return i;
        }
        return s.length() - 1;
    }

    private static String jsonStr(String json, String key) {
        String search = "\"" + key + "\":";
        int i = json.indexOf(search);
        if (i < 0) return "";
        i += search.length();
        while (i < json.length() && json.charAt(i) == ' ') i++;
        if (i >= json.length()) return "";
        if (json.charAt(i) == '"') {
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
        int end = i;
        while (end < json.length() && ",}] \n\r\t".indexOf(json.charAt(end)) < 0) end++;
        return json.substring(i, end).trim();
    }

    private static String esc(String s) {
        if (s == null) return "";
        // Normalize Windows backslashes to forward slashes for paths
        return s.replace("\\", "/").replace("\"", "\\\"")
                .replace("\n", "\\n").replace("\r", "\\r");
    }
}
