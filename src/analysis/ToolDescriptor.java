
import java.io.*;
import java.util.*;

/**
 * Declarative tool plugin descriptor, parsed from YAML.
 * Each descriptor defines: tool identity, execution command, required base artifacts,
 * user-configurable parameters (auto-rendered as a form), and output column mappings.
 *
 * Adding a new tool = drop a YAML descriptor + optional adapter script.
 * Zero core engine changes.
 */
public class ToolDescriptor {

    public String tool;           // e.g. "cojo_conditional"
    public String version;        // e.g. "1.0"
    public String label;          // human-readable name for UI
    public String description;    // short help text
    public String language;       // "binary", "R", "python"
    public String command;        // command template with {placeholders}
    public List<String> requires = new ArrayList<>();  // which base artifacts needed
    public List<Param> params = new ArrayList<>();
    public List<OutputColumn> outputMapping = new ArrayList<>();
    public String outputKey = "snp_id";  // key column in raw output
    public String inputAdapter;   // optional: class name for custom input prep
    public String outputAdapter;  // optional: class name for custom output parsing

    public static class Param {
        public String name;
        public String type;      // "float", "int", "string", "boolean", "select"
        public String defaultValue;
        public String label;
        public String description;
        public List<String> options;  // for "select" type

        public Map<String, Object> toMap() {
            Map<String, Object> m = new LinkedHashMap<>();
            m.put("name", name);
            m.put("type", type);
            m.put("default", defaultValue);
            m.put("label", label != null ? label : name);
            if (description != null) m.put("description", description);
            if (options != null) m.put("options", options);
            return m;
        }
    }

    public static class OutputColumn {
        public String raw;       // column name in tool's raw output
        public String as;        // unified column name
        public String scope;     // "per_snp" or "per_credible_set"
        public String type;      // "double", "int", "string", "boolean"
    }

    /**
     * Parse a YAML descriptor file (simple key-value + list parsing, no library needed).
     */
    public static ToolDescriptor parse(File yamlFile) throws IOException {
        ToolDescriptor td = new ToolDescriptor();
        String content = new String(java.nio.file.Files.readAllBytes(yamlFile.toPath()), "UTF-8");
        String[] lines = content.split("\n");

        String section = null;
        Param currentParam = null;
        OutputColumn currentCol = null;

        for (String rawLine : lines) {
            String line = rawLine.replace("\r", "");
            String trimmed = line.trim();
            if (trimmed.isEmpty() || trimmed.startsWith("#")) continue;

            int indent = countIndent(line);

            // Top-level fields
            if (indent == 0 && trimmed.contains(":")) {
                section = null;
                currentParam = null;
                currentCol = null;
                String[] kv = splitFirst(trimmed, ":");
                String key = kv[0].trim();
                String val = kv.length > 1 ? kv[1].trim() : "";
                val = unquote(val);

                switch (key) {
                    case "tool": td.tool = val; break;
                    case "version": td.version = val; break;
                    case "label": td.label = val; break;
                    case "description": td.description = val; break;
                    case "language": td.language = val; break;
                    case "command": td.command = val; break;
                    case "output_key": td.outputKey = val; break;
                    case "input_adapter": td.inputAdapter = val; break;
                    case "output_adapter": td.outputAdapter = val; break;
                    case "requires": section = "requires"; break;
                    case "params": section = "params"; break;
                    case "output_mapping": section = "output_mapping"; break;
                }
                // Inline list: requires: [harmonized_gwas, matched_ref, ld_r]
                if (key.equals("requires") && val.startsWith("[")) {
                    td.requires.clear();
                    String inner = val.substring(1, val.lastIndexOf(']'));
                    for (String item : inner.split(","))
                        td.requires.add(item.trim());
                    section = null;
                }
                continue;
            }

            // List items
            if ("requires".equals(section) && trimmed.startsWith("- ")) {
                td.requires.add(trimmed.substring(2).trim());
                continue;
            }

            if ("params".equals(section)) {
                if (trimmed.startsWith("- ")) {
                    currentParam = new Param();
                    td.params.add(currentParam);
                    // Inline object: - {name: x, type: float, default: 5e-8}
                    if (trimmed.contains("{")) {
                        parseInlineObject(trimmed.substring(2).trim(), currentParam);
                        continue;
                    }
                    String sub = trimmed.substring(2).trim();
                    if (sub.contains(":")) {
                        String[] kv = splitFirst(sub, ":");
                        applyParamField(currentParam, kv[0].trim(), kv.length > 1 ? unquote(kv[1].trim()) : "");
                    }
                } else if (currentParam != null && trimmed.contains(":")) {
                    String[] kv = splitFirst(trimmed, ":");
                    applyParamField(currentParam, kv[0].trim(), kv.length > 1 ? unquote(kv[1].trim()) : "");
                }
                continue;
            }

            if ("output_mapping".equals(section)) {
                if (trimmed.startsWith("- ")) {
                    currentCol = new OutputColumn();
                    td.outputMapping.add(currentCol);
                    if (trimmed.contains("{")) {
                        parseInlineOutputCol(trimmed.substring(2).trim(), currentCol);
                        continue;
                    }
                    String sub = trimmed.substring(2).trim();
                    if (sub.contains(":")) {
                        String[] kv = splitFirst(sub, ":");
                        applyOutputField(currentCol, kv[0].trim(), kv.length > 1 ? unquote(kv[1].trim()) : "");
                    }
                } else if (currentCol != null && trimmed.contains(":")) {
                    String[] kv = splitFirst(trimmed, ":");
                    applyOutputField(currentCol, kv[0].trim(), kv.length > 1 ? unquote(kv[1].trim()) : "");
                }
                continue;
            }
        }

        if (td.label == null || td.label.isEmpty()) td.label = td.tool;
        return td;
    }

    /**
     * Produce a JSON schema for auto-generating the UI form.
     */
    public String toFormJson() {
        StringBuilder sb = new StringBuilder();
        sb.append("{\"tool\":\"").append(esc(tool)).append('"');
        sb.append(",\"version\":\"").append(esc(version)).append('"');
        sb.append(",\"label\":\"").append(esc(label)).append('"');
        if (description != null)
            sb.append(",\"description\":\"").append(esc(description)).append('"');
        sb.append(",\"language\":\"").append(esc(language)).append('"');

        sb.append(",\"params\":[");
        for (int i = 0; i < params.size(); i++) {
            if (i > 0) sb.append(',');
            Param p = params.get(i);
            sb.append("{\"name\":\"").append(esc(p.name)).append('"');
            sb.append(",\"type\":\"").append(esc(p.type)).append('"');
            sb.append(",\"default\":\"").append(esc(p.defaultValue)).append('"');
            sb.append(",\"label\":\"").append(esc(p.label != null ? p.label : p.name)).append('"');
            if (p.description != null)
                sb.append(",\"description\":\"").append(esc(p.description)).append('"');
            if (p.options != null && !p.options.isEmpty()) {
                sb.append(",\"options\":[");
                for (int j = 0; j < p.options.size(); j++) {
                    if (j > 0) sb.append(',');
                    sb.append('"').append(esc(p.options.get(j))).append('"');
                }
                sb.append(']');
            }
            sb.append('}');
        }
        sb.append("]}");
        return sb.toString();
    }

    // ── Parsing helpers ──

    private static void parseInlineObject(String text, Param p) {
        text = text.replaceAll("^\\{|\\}$", "");
        for (String part : text.split(",")) {
            String[] kv = splitFirst(part.trim(), ":");
            if (kv.length == 2) applyParamField(p, kv[0].trim(), unquote(kv[1].trim()));
        }
    }

    private static void parseInlineOutputCol(String text, OutputColumn c) {
        text = text.replaceAll("^\\{|\\}$", "");
        for (String part : text.split(",")) {
            String[] kv = splitFirst(part.trim(), ":");
            if (kv.length == 2) applyOutputField(c, kv[0].trim(), unquote(kv[1].trim()));
        }
    }

    private static void applyParamField(Param p, String key, String val) {
        switch (key) {
            case "name": p.name = val; break;
            case "type": p.type = val; break;
            case "default": p.defaultValue = val; break;
            case "label": p.label = val; break;
            case "description": p.description = val; break;
        }
    }

    private static void applyOutputField(OutputColumn c, String key, String val) {
        switch (key) {
            case "raw": c.raw = val; break;
            case "as": c.as = val; break;
            case "scope": c.scope = val; break;
            case "type": c.type = val; break;
        }
    }

    private static int countIndent(String line) {
        int n = 0;
        for (char c : line.toCharArray()) {
            if (c == ' ') n++; else break;
        }
        return n;
    }

    private static String[] splitFirst(String s, String delim) {
        int idx = s.indexOf(delim);
        if (idx < 0) return new String[]{s};
        return new String[]{s.substring(0, idx), s.substring(idx + delim.length())};
    }

    private static String unquote(String s) {
        if (s.length() >= 2 &&
            ((s.startsWith("\"") && s.endsWith("\"")) ||
             (s.startsWith("'") && s.endsWith("'")))) {
            return s.substring(1, s.length() - 1);
        }
        return s;
    }

    private static String esc(String s) {
        if (s == null) return "";
        return s.replace("\\", "\\\\").replace("\"", "\\\"");
    }
}
