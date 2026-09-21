import java.io.IOException;
import java.util.*;

/**
 * The agent's tool-calling loop (DECISIONS_PHASE4.md sections 1/2.4). Stateless on the server side by
 * design: the caller (the browser tab) resends the full prior message history each turn and receives
 * the updated history back, rather than this class keeping per-session state in memory — simpler, and
 * consistent with this app being a single-user local tool with no login/session concept to hang state
 * off of.
 *
 * Every tool call the model makes is recorded as a {@link Step} in call order, in full, so the
 * frontend can show the user exactly what the agent did — no hidden actions (DECISIONS_PHASE4.md
 * section 1, "Transparency"). The loop is capped at {@link #MAX_TOOL_ROUNDS} round-trips to the LLM so
 * a model that keeps calling tools without ever producing a final answer cannot loop forever.
 */
public class AgentOrchestrator {

    private static final int MAX_TOOL_ROUNDS = 6;

    private static final String SYSTEM_PROMPT =
        "You are the built-in assistant for LYNXgwas, a local GWAS locus-analysis tool. You can " +
        "inspect and configure the user's projects and run its pipeline/analysis tools on their " +
        "behalf using the tools made available to you. Prefer calling a tool over guessing when the " +
        "answer depends on the user's actual data (project list, config, manifest, tool parameters). " +
        "When creating or editing a project, only pass the config fields the user actually specified " +
        "or that are clearly implied — do not invent file paths or column names. Long-running actions " +
        "(run_pipeline, run_analysis_tool) run in the background; tell the user you started them and " +
        "that they can ask you to check progress. Be concise.";

    public static class Step {
        public final String tool;
        public final String argumentsJson;
        public final String summary;
        public final boolean ok;
        Step(String tool, String argumentsJson, String summary, boolean ok) {
            this.tool = tool; this.argumentsJson = argumentsJson; this.summary = summary; this.ok = ok;
        }
    }

    public static class ChatOutcome {
        public final String reply;
        public final List<Step> steps;
        public final List<Map<String, Object>> messages; // full updated history — echoed back to the client
        ChatOutcome(String reply, List<Step> steps, List<Map<String, Object>> messages) {
            this.reply = reply; this.steps = steps; this.messages = messages;
        }
    }

    public static ChatOutcome chat(List<Map<String, Object>> history, String userMessage,
                                    AgentConfig cfg, int serverPort) throws IOException, InterruptedException {
        List<Map<String, Object>> messages = new ArrayList<>(history);
        if (messages.isEmpty()) messages.add(msg("system", SYSTEM_PROMPT));
        messages.add(msg("user", userMessage));

        AgentToolRegistry registry = new AgentToolRegistry(serverPort);
        LlmClient client = new LlmClient(cfg.baseUrl, cfg.apiKey, cfg.model);
        List<Step> steps = new ArrayList<>();
        List<Map<String, Object>> toolSchemas = registry.toolSchemas();

        for (int round = 0; round < MAX_TOOL_ROUNDS; round++) {
            LlmClient.ChatResult result = client.chat(messages, toolSchemas);

            if (result.toolCalls.isEmpty()) {
                String reply = result.assistantText == null ? "" : result.assistantText;
                messages.add(msg("assistant", reply));
                return new ChatOutcome(reply, steps, messages);
            }

            messages.add(assistantToolCallsMessage(result.toolCalls));
            for (LlmClient.ToolCall call : result.toolCalls) {
                Map<String, Object> args;
                try {
                    args = MiniJson.asObject(MiniJson.parse(call.argumentsJson));
                } catch (Exception e) {
                    args = new LinkedHashMap<>();
                }
                AgentToolRegistry.ToolExecutionResult execResult = registry.dispatch(call.name, args);
                steps.add(new Step(call.name, call.argumentsJson, execResult.summary, execResult.ok));
                messages.add(toolResultMessage(call.id, call.name, execResult.resultJson));
            }
        }

        String reply = "I stopped after " + MAX_TOOL_ROUNDS + " tool-call rounds without reaching a final " +
            "answer — this is a safety limit, not an error. You can ask me to continue.";
        messages.add(msg("assistant", reply));
        return new ChatOutcome(reply, steps, messages);
    }

    private static Map<String, Object> msg(String role, String content) {
        Map<String, Object> m = new LinkedHashMap<>();
        m.put("role", role);
        m.put("content", content);
        return m;
    }

    private static Map<String, Object> assistantToolCallsMessage(List<LlmClient.ToolCall> calls) {
        Map<String, Object> m = new LinkedHashMap<>();
        m.put("role", "assistant");
        m.put("content", null);
        List<Object> calls_ = new ArrayList<>();
        for (LlmClient.ToolCall c : calls) {
            Map<String, Object> fn = new LinkedHashMap<>();
            fn.put("name", c.name);
            fn.put("arguments", c.argumentsJson);
            Map<String, Object> call = new LinkedHashMap<>();
            call.put("id", c.id);
            call.put("type", "function");
            call.put("function", fn);
            calls_.add(call);
        }
        m.put("tool_calls", calls_);
        return m;
    }

    private static Map<String, Object> toolResultMessage(String toolCallId, String name, String resultJson) {
        Map<String, Object> m = new LinkedHashMap<>();
        m.put("role", "tool");
        m.put("tool_call_id", toolCallId);
        m.put("name", name);
        m.put("content", resultJson);
        return m;
    }
}
