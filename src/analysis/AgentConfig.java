import java.io.*;
import java.nio.charset.StandardCharsets;
import java.nio.file.*;
import java.util.*;

/**
 * Persisted agent LLM settings (DECISIONS_PHASE4.md section 1) — {@code mode} is purely descriptive
 * for the UI ("local" vs "cloud"); functionally both modes just populate {@code baseUrl}/{@code
 * apiKey}/{@code model} for {@link LlmClient}, since both speak the same OpenAI-compatible protocol.
 * Stored at {@code <workspace>/.agent_config.json} (gitignored). The API key is written to that local
 * file and returned as-is by the config-read endpoint (this is a single-user, single-machine
 * application per the project's existing security model — see DECISIONS.md's "Single-user,
 * single-machine design" limitation); it is never logged or included in any other response.
 */
public class AgentConfig {

    public String mode = "local";                       // "local" | "cloud" — UI hint only
    public String baseUrl = "http://localhost:11434/v1"; // Ollama's own default OpenAI-compatible endpoint
    public String apiKey = "";
    public String model = "";

    private static final String FILE_NAME = ".agent_config.json";

    public static AgentConfig load() {
        AgentConfig cfg = new AgentConfig();
        File f = new File(FILE_NAME);
        if (!f.exists()) return cfg;
        try {
            String json = new String(Files.readAllBytes(f.toPath()), StandardCharsets.UTF_8);
            Map<String, Object> obj = MiniJson.asObject(MiniJson.parse(json));
            if (obj.get("mode") != null) cfg.mode = String.valueOf(obj.get("mode"));
            if (obj.get("base_url") != null) cfg.baseUrl = String.valueOf(obj.get("base_url"));
            if (obj.get("api_key") != null) cfg.apiKey = String.valueOf(obj.get("api_key"));
            if (obj.get("model") != null) cfg.model = String.valueOf(obj.get("model"));
        } catch (Exception e) {
            System.err.println("[AgentConfig] Failed to read " + FILE_NAME + " (using defaults): " + e.getMessage());
        }
        return cfg;
    }

    public void save() throws IOException {
        String json = "{"
            + "\"mode\":" + jsonStr(mode) + ","
            + "\"base_url\":" + jsonStr(baseUrl) + ","
            + "\"api_key\":" + jsonStr(apiKey) + ","
            + "\"model\":" + jsonStr(model)
            + "}";
        Files.write(new File(FILE_NAME).toPath(), json.getBytes(StandardCharsets.UTF_8));
    }

    public String toJson() {
        return "{"
            + "\"mode\":" + jsonStr(mode) + ","
            + "\"base_url\":" + jsonStr(baseUrl) + ","
            + "\"api_key\":" + jsonStr(apiKey) + ","
            + "\"model\":" + jsonStr(model)
            + "}";
    }

    public void applyJson(String json) {
        Map<String, Object> obj = MiniJson.asObject(MiniJson.parse(json));
        if (obj.get("mode") != null) mode = String.valueOf(obj.get("mode"));
        if (obj.get("base_url") != null) baseUrl = String.valueOf(obj.get("base_url"));
        if (obj.get("api_key") != null) apiKey = String.valueOf(obj.get("api_key"));
        if (obj.get("model") != null) model = String.valueOf(obj.get("model"));
    }

    private static String jsonStr(String s) {
        if (s == null) s = "";
        StringBuilder sb = new StringBuilder("\"");
        for (int i = 0; i < s.length(); i++) {
            char c = s.charAt(i);
            if (c == '"') sb.append("\\\"");
            else if (c == '\\') sb.append("\\\\");
            else if (c == '\n') sb.append("\\n");
            else if (c < 0x20) sb.append(String.format("\\u%04x", (int) c));
            else sb.append(c);
        }
        return sb.append("\"").toString();
    }
}
