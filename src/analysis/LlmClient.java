import java.io.*;
import java.net.*;
import java.net.http.*;
import java.nio.charset.StandardCharsets;
import java.time.Duration;
import java.util.*;

/**
 * OpenAI-compatible Chat Completions client with tool-calling (DECISIONS_PHASE4.md section 1/2.2).
 * One implementation covers both requested cases — a locally-run open-source model and a
 * user-supplied cloud API key — because Ollama, LM Studio, and llama.cpp's own {@code server} binary
 * all expose this exact same request/response shape at their own {@code /v1} endpoint, and so do
 * cloud providers using it directly (OpenAI itself, Groq, Together, etc.). Only {@code base_url},
 * {@code api_key} (may be blank for a local runtime that doesn't check one), and {@code model}
 * differ between "local" and "cloud" — see {@link AgentConfig}.
 *
 * Anthropic's native Messages API uses a different shape and is intentionally not implemented here;
 * this class exists as the seam a future {@code AnthropicClient} would sit behind instead of being
 * special-cased into callers.
 */
public class LlmClient {

    public static class ToolCall {
        public final String id;
        public final String name;
        public final String argumentsJson; // raw JSON text, caller re-parses with MiniJson
        public ToolCall(String id, String name, String argumentsJson) {
            this.id = id; this.name = name; this.argumentsJson = argumentsJson;
        }
    }

    public static class ChatResult {
        public final String assistantText;      // may be null if the model only returned tool calls
        public final List<ToolCall> toolCalls;   // empty if the model returned a final text answer
        public ChatResult(String assistantText, List<ToolCall> toolCalls) {
            this.assistantText = assistantText; this.toolCalls = toolCalls;
        }
    }

    private final String baseUrl;
    private final String apiKey;
    private final String model;
    private final HttpClient http;

    public LlmClient(String baseUrl, String apiKey, String model) {
        // Tolerate a trailing slash so both "http://host:port/v1" and "http://host:port/v1/" work.
        this.baseUrl = baseUrl.endsWith("/") ? baseUrl.substring(0, baseUrl.length() - 1) : baseUrl;
        this.apiKey = apiKey;
        this.model = model;
        this.http = HttpClient.newBuilder().connectTimeout(Duration.ofSeconds(10)).build();
    }

    /**
     * @param messages  each a Map with at least "role" and "content" (content may be null for an
     *                  assistant message that only carried tool_calls, or a "tool" role message
     *                  whose content is the tool's JSON result string).
     * @param toolSchemas  the OpenAI "tools" array (function name/description/parameters JSON-schema),
     *                  as raw pre-built JSON objects from {@link AgentToolRegistry}.
     */
    public ChatResult chat(List<Map<String, Object>> messages, List<Map<String, Object>> toolSchemas)
            throws IOException, InterruptedException {
        String body = buildRequestJson(messages, toolSchemas);

        HttpRequest.Builder reqBuilder = HttpRequest.newBuilder()
            .uri(URI.create(baseUrl + "/chat/completions"))
            .timeout(Duration.ofSeconds(120))
            .header("Content-Type", "application/json");
        if (apiKey != null && !apiKey.isEmpty())
            reqBuilder.header("Authorization", "Bearer " + apiKey);
        HttpRequest req = reqBuilder.POST(HttpRequest.BodyPublishers.ofString(body, StandardCharsets.UTF_8)).build();

        HttpResponse<String> resp = http.send(req, HttpResponse.BodyHandlers.ofString(StandardCharsets.UTF_8));
        if (resp.statusCode() / 100 != 2) {
            throw new IOException("LLM endpoint returned HTTP " + resp.statusCode() + ": "
                + truncate(resp.body(), 500));
        }
        return parseResponse(resp.body());
    }

    // ── Request building (hand-rolled JSON, matching this project's existing kv/esc convention) ──

    private String buildRequestJson(List<Map<String, Object>> messages, List<Map<String, Object>> toolSchemas) {
        StringBuilder sb = new StringBuilder();
        sb.append("{");
        sb.append("\"model\":").append(jsonStr(model)).append(",");
        sb.append("\"messages\":[");
        for (int i = 0; i < messages.size(); i++) {
            if (i > 0) sb.append(",");
            sb.append(messageToJson(messages.get(i)));
        }
        sb.append("]");
        if (toolSchemas != null && !toolSchemas.isEmpty()) {
            sb.append(",\"tools\":[");
            for (int i = 0; i < toolSchemas.size(); i++) {
                if (i > 0) sb.append(",");
                sb.append(jsonEncodeAny(toolSchemas.get(i)));
            }
            sb.append("],\"tool_choice\":\"auto\"");
        }
        sb.append("}");
        return sb.toString();
    }

    @SuppressWarnings("unchecked")
    private String messageToJson(Map<String, Object> m) {
        StringBuilder sb = new StringBuilder("{");
        sb.append("\"role\":").append(jsonStr((String) m.get("role")));
        Object content = m.get("content");
        sb.append(",\"content\":").append(content == null ? "null" : jsonStr(String.valueOf(content)));
        if (m.get("tool_call_id") != null)
            sb.append(",\"tool_call_id\":").append(jsonStr(String.valueOf(m.get("tool_call_id"))));
        if (m.get("name") != null)
            sb.append(",\"name\":").append(jsonStr(String.valueOf(m.get("name"))));
        if (m.get("tool_calls") != null) {
            sb.append(",\"tool_calls\":[");
            List<Map<String, Object>> calls = (List<Map<String, Object>>) m.get("tool_calls");
            for (int i = 0; i < calls.size(); i++) {
                if (i > 0) sb.append(",");
                sb.append(jsonEncodeAny(calls.get(i)));
            }
            sb.append("]");
        }
        sb.append("}");
        return sb.toString();
    }

    // Tool schemas and echoed-back tool_calls are plain Map/List/String/Number/Boolean structures —
    // reuse MiniJson's general encoder rather than duplicating it here.
    private String jsonEncodeAny(Object v) {
        return MiniJson.encode(v);
    }

    private String jsonStr(String s) {
        if (s == null) return "null";
        StringBuilder sb = new StringBuilder("\"");
        for (int i = 0; i < s.length(); i++) {
            char c = s.charAt(i);
            switch (c) {
                case '"':  sb.append("\\\""); break;
                case '\\': sb.append("\\\\"); break;
                case '\n': sb.append("\\n"); break;
                case '\r': sb.append("\\r"); break;
                case '\t': sb.append("\\t"); break;
                default:
                    if (c < 0x20) sb.append(String.format("\\u%04x", (int) c));
                    else sb.append(c);
            }
        }
        return sb.append("\"").toString();
    }

    // ── Response parsing (via MiniJson) ──

    @SuppressWarnings("unchecked")
    private ChatResult parseResponse(String body) {
        Map<String, Object> root = MiniJson.asObject(MiniJson.parse(body));
        List<Object> choices = MiniJson.asArray(root.get("choices"));
        if (choices.isEmpty()) throw new IllegalArgumentException("LLM response had no choices: " + truncate(body, 300));
        Map<String, Object> message = MiniJson.asObject(MiniJson.asObject(choices.get(0)).get("message"));

        String text = message.get("content") == null ? null : String.valueOf(message.get("content"));
        List<ToolCall> calls = new ArrayList<>();
        Object rawCalls = message.get("tool_calls");
        if (rawCalls != null) {
            for (Object o : MiniJson.asArray(rawCalls)) {
                Map<String, Object> call = MiniJson.asObject(o);
                Map<String, Object> fn = MiniJson.asObject(call.get("function"));
                calls.add(new ToolCall(
                    String.valueOf(call.get("id")),
                    String.valueOf(fn.get("name")),
                    String.valueOf(fn.get("arguments"))
                ));
            }
        }
        return new ChatResult(text, calls);
    }

    private static String truncate(String s, int max) {
        return s.length() <= max ? s : s.substring(0, max) + "...";
    }
}
