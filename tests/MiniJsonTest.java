import java.util.*;

/**
 * Regression test for MiniJson (DECISIONS_PHASE4.md section 2.1), covering: primitives, escapes
 * (including a unicode escape), nesting, and — the actual reason this parser exists — a real-shaped
 * OpenAI-compatible chat-completions tool-call response (nested objects/arrays with a JSON-encoded
 * string as one of the values, exactly how `arguments` arrives on a real tool call).
 */
public class MiniJsonTest {

    public static void main(String[] args) throws Exception {
        int failures = 0;

        failures += check("parses a plain string", "hello".equals(MiniJson.parse("\"hello\"")));
        failures += check("parses an integer as Double", MiniJson.parse("42").equals(42.0));
        failures += check("parses a negative decimal", MiniJson.parse("-3.5").equals(-3.5));
        failures += check("parses exponent notation", MiniJson.parse("1.5e2").equals(150.0));
        failures += check("parses true", Boolean.TRUE.equals(MiniJson.parse("true")));
        failures += check("parses false", Boolean.FALSE.equals(MiniJson.parse("false")));
        failures += check("parses null", MiniJson.parse("null") == null);

        failures += check("parses an empty object", MiniJson.asObject(MiniJson.parse("{}")).isEmpty());
        failures += check("parses an empty array", MiniJson.asArray(MiniJson.parse("[]")).isEmpty());

        // ── Escapes ──
        String escaped = MiniJson.parse("\"line1\\nline2\\t\\\"quoted\\\"\\u0041\"").toString();
        failures += check("handles \\n, \\t, \\\", and \\u unicode escape",
            escaped.equals("line1\nline2\t\"quoted\"A"));

        // ── Nesting ──
        Map<String, Object> obj = MiniJson.asObject(MiniJson.parse(
            "{\"a\": 1, \"b\": [1, 2, 3], \"c\": {\"d\": \"x\"}, \"e\": null}"));
        failures += check("object has 4 keys", obj.size() == 4);
        failures += check("nested array has 3 elements", MiniJson.asArray(obj.get("b")).size() == 3);
        failures += check("nested object field 'd' == 'x'",
            "x".equals(MiniJson.asObject(obj.get("c")).get("d")));
        failures += check("explicit null value preserved as null", obj.containsKey("e") && obj.get("e") == null);

        // ── Real-shaped OpenAI-compatible tool-call response ──
        String toolCallResponse = "{"
            + "\"choices\": [{"
            + "  \"message\": {"
            + "    \"role\": \"assistant\","
            + "    \"content\": null,"
            + "    \"tool_calls\": [{"
            + "      \"id\": \"call_abc123\","
            + "      \"type\": \"function\","
            + "      \"function\": {"
            + "        \"name\": \"create_or_update_project\","
            + "        \"arguments\": \"{\\\"project_id\\\":\\\"demo\\\",\\\"gwas_file\\\":\\\"input\\\\\\\\gwas.tsv\\\"}\""
            + "      }"
            + "    }]"
            + "  }"
            + "}]"
            + "}";
        Map<String, Object> resp = MiniJson.asObject(MiniJson.parse(toolCallResponse));
        List<Object> choices = MiniJson.asArray(resp.get("choices"));
        failures += check("one choice present", choices.size() == 1);
        Map<String, Object> message = MiniJson.asObject(MiniJson.asObject(choices.get(0)).get("message"));
        failures += check("message content is explicit null", message.containsKey("content") && message.get("content") == null);
        List<Object> toolCalls = MiniJson.asArray(message.get("tool_calls"));
        failures += check("exactly one tool call", toolCalls.size() == 1);
        Map<String, Object> call = MiniJson.asObject(toolCalls.get(0));
        failures += check("tool call id extracted", "call_abc123".equals(call.get("id")));
        Map<String, Object> fn = MiniJson.asObject(call.get("function"));
        failures += check("function name extracted", "create_or_update_project".equals(fn.get("name")));
        String argsJson = (String) fn.get("arguments");
        failures += check("arguments is a JSON-encoded string (double-escaped quotes decoded once)",
            argsJson.equals("{\"project_id\":\"demo\",\"gwas_file\":\"input\\\\gwas.tsv\"}"));
        // The arguments string is itself JSON and must parse cleanly a second time (this is how a
        // real tool-call payload works: the outer JSON's string value is inner JSON text).
        Map<String, Object> parsedArgs = MiniJson.asObject(MiniJson.parse(argsJson));
        failures += check("re-parsed arguments has project_id=demo", "demo".equals(parsedArgs.get("project_id")));
        failures += check("re-parsed arguments preserves a Windows-style backslash path",
            "input\\gwas.tsv".equals(parsedArgs.get("gwas_file")));

        // ── Malformed input must throw, not silently return a wrong value ──
        boolean threw = false;
        try { MiniJson.parse("{\"a\": }"); } catch (IllegalArgumentException e) { threw = true; }
        failures += check("malformed object throws IllegalArgumentException", threw);

        threw = false;
        try { MiniJson.parse("{\"a\": 1} trailing"); } catch (IllegalArgumentException e) { threw = true; }
        failures += check("trailing content after a valid value throws", threw);

        if (failures == 0) {
            System.out.println("PASS: all MiniJson tests passed");
        } else {
            System.out.println("FAIL: " + failures + " MiniJson test(s) failed");
            System.exit(1);
        }
    }

    private static int check(String label, boolean cond) {
        if (!cond) {
            System.out.println("FAIL: " + label);
            return 1;
        }
        System.out.println("PASS: " + label);
        return 0;
    }
}
