
import java.io.*;
import java.util.*;

/**
 * Validates tool output against the unified output contract.
 *
 * Every tool must produce:
 *   result.tsv — long-format, keyed by stable SNP id (chr:pos:ref:alt)
 *   result.manifest.json — declares each value column (name, type, scope, method)
 *
 * The validator rejects nonconforming output with a clear error message
 * before it reaches the viewer or exports.
 */
public class OutputContractValidator {

    public static final String OUTPUT_SCHEMA_VERSION = "1.0";

    public static class ValidationResult {
        public boolean valid;
        public List<String> errors = new ArrayList<>();
        public int rowCount;
        public List<ColumnDecl> columns = new ArrayList<>();
    }

    public static class ColumnDecl {
        public String name;
        public String type;    // "double", "int", "string", "boolean"
        public String scope;   // "per_snp" or "per_credible_set"
        public String method;
        public String methodVersion;

        @Override public String toString() {
            return name + "(" + type + ", " + scope + ")";
        }
    }

    /**
     * Validate result.tsv + result.manifest.json in the given run directory.
     */
    public static ValidationResult validate(File runDir) {
        ValidationResult vr = new ValidationResult();

        File manifestFile = new File(runDir, "result.manifest.json");
        File resultFile = new File(runDir, "result.tsv");

        // ── Check files exist ──
        if (!manifestFile.exists()) {
            vr.errors.add("result.manifest.json not found in " + runDir.getAbsolutePath());
            return vr;
        }
        if (!resultFile.exists()) {
            vr.errors.add("result.tsv not found in " + runDir.getAbsolutePath());
            return vr;
        }

        // ── Parse manifest ──
        String manifestJson;
        try {
            manifestJson = new String(java.nio.file.Files.readAllBytes(manifestFile.toPath()), "UTF-8");
        } catch (IOException e) {
            vr.errors.add("Cannot read result.manifest.json: " + e.getMessage());
            return vr;
        }

        String schemaVer = extractField(manifestJson, "schema_version");
        if (schemaVer == null) {
            vr.errors.add("Missing schema_version in result.manifest.json");
        } else if (!OUTPUT_SCHEMA_VERSION.equals(schemaVer)) {
            vr.errors.add("Schema version mismatch: expected " + OUTPUT_SCHEMA_VERSION + ", got " + schemaVer);
        }

        // Parse column declarations from the "columns" array
        vr.columns = parseColumns(manifestJson);
        if (vr.columns.isEmpty()) {
            vr.errors.add("No columns declared in result.manifest.json");
            return vr;
        }

        Set<String> declaredCols = new HashSet<>();
        for (ColumnDecl cd : vr.columns) {
            if (cd.name == null || cd.name.isEmpty()) {
                vr.errors.add("Column with empty name in manifest");
            }
            if (cd.type == null || !isValidType(cd.type)) {
                vr.errors.add("Column '" + cd.name + "' has invalid type: " + cd.type);
            }
            if (cd.scope == null || (!cd.scope.equals("per_snp") && !cd.scope.equals("per_credible_set"))) {
                vr.errors.add("Column '" + cd.name + "' has invalid scope: " + cd.scope + " (must be per_snp or per_credible_set)");
            }
            declaredCols.add(cd.name);
        }

        // ── Validate result.tsv ──
        try (BufferedReader br = new BufferedReader(new FileReader(resultFile))) {
            String headerLine = br.readLine();
            if (headerLine == null) {
                vr.errors.add("result.tsv is empty");
                return vr;
            }

            String[] headers = headerLine.split("\t", -1);

            // Must have snp_id as first column
            if (headers.length == 0 || !headers[0].trim().equals("snp_id")) {
                vr.errors.add("First column of result.tsv must be 'snp_id', got: '" +
                    (headers.length > 0 ? headers[0].trim() : "") + "'");
            }

            // Check all non-key columns are declared in manifest
            Set<String> tsvCols = new HashSet<>();
            for (int i = 1; i < headers.length; i++) {
                String col = headers[i].trim();
                tsvCols.add(col);
                if (!declaredCols.contains(col)) {
                    vr.errors.add("Column '" + col + "' in result.tsv not declared in manifest");
                }
            }

            // Check all declared columns are present in TSV
            for (String dc : declaredCols) {
                if (!tsvCols.contains(dc)) {
                    vr.errors.add("Declared column '" + dc + "' missing from result.tsv");
                }
            }

            // Validate rows
            int rowCount = 0;
            String line;
            while ((line = br.readLine()) != null) {
                if (line.trim().isEmpty()) continue;
                rowCount++;
                String[] fields = line.split("\t", -1);

                if (fields.length != headers.length) {
                    if (rowCount <= 3) {
                        vr.errors.add("Row " + rowCount + " has " + fields.length +
                            " columns, expected " + headers.length);
                    }
                    continue;
                }

                // Validate snp_id is non-empty (accept rsIDs, chr:pos:ref:alt, or BIM IDs)
                String snpId = fields[0].trim();
                if (snpId.isEmpty() && rowCount <= 3) {
                    vr.errors.add("Row " + rowCount + " has empty snp_id");
                }
            }
            vr.rowCount = rowCount;

        } catch (IOException e) {
            vr.errors.add("Cannot read result.tsv: " + e.getMessage());
        }

        vr.valid = vr.errors.isEmpty();
        return vr;
    }

    /**
     * Write a result manifest for a tool's output.
     */
    public static void writeManifest(File runDir, String method, String methodVersion,
                                      Map<String, String> params,
                                      List<ColumnDecl> columns) throws IOException {
        File outFile = new File(runDir, "result.manifest.json");
        try (PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(outFile)))) {
            pw.println("{");
            pw.printf("  \"schema_version\": \"%s\",%n", OUTPUT_SCHEMA_VERSION);
            pw.printf("  \"method\": \"%s\",%n", esc(method));
            pw.printf("  \"method_version\": \"%s\",%n", esc(methodVersion));

            // Parameters
            pw.println("  \"parameters\": {");
            int pi = 0;
            for (Map.Entry<String, String> e : params.entrySet()) {
                pw.printf("    \"%s\": \"%s\"%s%n",
                    esc(e.getKey()), esc(e.getValue()), (++pi < params.size()) ? "," : "");
            }
            pw.println("  },");

            // Columns
            pw.println("  \"columns\": [");
            for (int i = 0; i < columns.size(); i++) {
                ColumnDecl c = columns.get(i);
                pw.printf("    {\"name\":\"%s\",\"type\":\"%s\",\"scope\":\"%s\",\"method\":\"%s\",\"method_version\":\"%s\"}%s%n",
                    esc(c.name), esc(c.type), esc(c.scope),
                    esc(c.method != null ? c.method : method),
                    esc(c.methodVersion != null ? c.methodVersion : methodVersion),
                    (i < columns.size() - 1) ? "," : "");
            }
            pw.println("  ],");

            pw.printf("  \"created_at\": %d%n", System.currentTimeMillis());
            pw.println("}");
        }
    }

    private static List<ColumnDecl> parseColumns(String json) {
        List<ColumnDecl> cols = new ArrayList<>();
        int arrStart = json.indexOf("\"columns\"");
        if (arrStart < 0) return cols;
        int bracketStart = json.indexOf('[', arrStart);
        if (bracketStart < 0) return cols;

        // Find matching close bracket
        int depth = 0;
        int bracketEnd = -1;
        for (int i = bracketStart; i < json.length(); i++) {
            if (json.charAt(i) == '[') depth++;
            else if (json.charAt(i) == ']') { depth--; if (depth == 0) { bracketEnd = i; break; } }
        }
        if (bracketEnd < 0) return cols;

        String arrContent = json.substring(bracketStart + 1, bracketEnd);

        // Split by },{ to get each column object
        int objStart = arrContent.indexOf('{');
        while (objStart >= 0 && objStart < arrContent.length()) {
            int objEnd = arrContent.indexOf('}', objStart);
            if (objEnd < 0) break;
            String obj = arrContent.substring(objStart + 1, objEnd);

            ColumnDecl cd = new ColumnDecl();
            cd.name = extractInline(obj, "name");
            cd.type = extractInline(obj, "type");
            cd.scope = extractInline(obj, "scope");
            cd.method = extractInline(obj, "method");
            cd.methodVersion = extractInline(obj, "method_version");
            if (cd.name != null) cols.add(cd);

            objStart = arrContent.indexOf('{', objEnd);
        }
        return cols;
    }

    private static String extractInline(String obj, String key) {
        String marker = "\"" + key + "\":\"";
        int i = obj.indexOf(marker);
        if (i < 0) return null;
        i += marker.length();
        int end = obj.indexOf('"', i);
        return end > i ? obj.substring(i, end) : null;
    }

    private static boolean isValidType(String type) {
        return "double".equals(type) || "int".equals(type) ||
               "string".equals(type) || "boolean".equals(type);
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
        int end = i;
        while (end < json.length() && ",}\n\r".indexOf(json.charAt(end)) < 0) end++;
        return json.substring(i, end).trim();
    }

    private static String esc(String s) {
        if (s == null) return "";
        return s.replace("\\", "\\\\").replace("\"", "\\\"");
    }
}
