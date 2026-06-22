
import java.io.*;
import java.nio.file.*;
import java.util.*;

/**
 * Tracks inputs, content hashes, and build metadata for a single base step.
 * Written as JSON; read back to decide if the step needs re-running.
 */
public class StepManifest {

    public String stepName;
    public String inputHash;
    public String genomeBuild = "GRCh37";
    public String refPanelPopulation = "";
    public long   builtAt;
    public Map<String, String> inputs = new LinkedHashMap<>();
    public Map<String, String> outputs = new LinkedHashMap<>();

    public boolean isCurrent(String currentInputHash) {
        return inputHash != null && inputHash.equals(currentInputHash);
    }

    public void write(File dir) throws IOException {
        dir.mkdirs();
        File f = new File(dir, stepName + ".manifest.json");
        try (PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(f)))) {
            pw.println("{");
            pw.printf("  \"step\": \"%s\",%n", esc(stepName));
            pw.printf("  \"input_hash\": \"%s\",%n", esc(inputHash));
            pw.printf("  \"genome_build\": \"%s\",%n", esc(genomeBuild));
            pw.printf("  \"ref_panel_population\": \"%s\",%n", esc(refPanelPopulation));
            pw.printf("  \"built_at\": %d,%n", builtAt);
            pw.println("  \"inputs\": {");
            writeMap(pw, inputs);
            pw.println("  },");
            pw.println("  \"outputs\": {");
            writeMap(pw, outputs);
            pw.println("  }");
            pw.println("}");
        }
    }

    public static StepManifest read(File dir, String stepName) {
        File f = new File(dir, stepName + ".manifest.json");
        if (!f.exists()) return null;
        try {
            String json = new String(Files.readAllBytes(f.toPath()), "UTF-8");
            StepManifest m = new StepManifest();
            m.stepName = stepName;
            m.inputHash = extractField(json, "input_hash");
            m.genomeBuild = extractField(json, "genome_build");
            m.refPanelPopulation = extractField(json, "ref_panel_population");
            String builtStr = extractField(json, "built_at");
            if (builtStr != null) {
                try { m.builtAt = Long.parseLong(builtStr); } catch (NumberFormatException e) {}
            }
            return m;
        } catch (IOException e) {
            return null;
        }
    }

    private void writeMap(PrintWriter pw, Map<String, String> map) {
        int i = 0;
        for (Map.Entry<String, String> e : map.entrySet()) {
            pw.printf("    \"%s\": \"%s\"%s%n",
                esc(e.getKey()), esc(e.getValue()), (++i < map.size()) ? "," : "");
        }
    }

    private static String extractField(String json, String key) {
        String marker = "\"" + key + "\":";
        int i = json.indexOf(marker);
        if (i < 0) return null;
        i += marker.length();
        while (i < json.length() && (json.charAt(i) == ' ' || json.charAt(i) == '"')) i++;
        if (i >= json.length()) return null;
        // Could be a number (no quotes) or a string
        int start = i;
        if (i > 0 && json.charAt(i - 1) == '"') {
            int end = json.indexOf('"', i);
            return end > i ? json.substring(i, end) : null;
        }
        int end = i;
        while (end < json.length() && ",}\n\r".indexOf(json.charAt(end)) < 0) end++;
        return json.substring(start, end).trim();
    }

    private static String esc(String s) {
        if (s == null) return "";
        return s.replace("\\", "\\\\").replace("\"", "\\\"");
    }
}
