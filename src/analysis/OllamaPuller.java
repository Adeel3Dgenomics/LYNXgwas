import java.net.URI;
import java.net.http.*;
import java.nio.charset.StandardCharsets;
import java.time.Duration;
import java.util.*;
import java.util.stream.Stream;

/**
 * Drives Ollama's own model-download API (DECISIONS_PHASE4.md follow-up: "make setting up a local
 * model easier"). Ollama exposes {@code POST /api/pull} — a real, stable, documented endpoint,
 * separate from the OpenAI-compatible surface {@link LlmClient} talks to — that streams
 * newline-delimited JSON progress objects while it downloads a named model. This lets the agent
 * settings panel offer a plain "Download" button instead of asking the user to open a terminal and
 * run {@code ollama pull <model>} themselves; it does not, and cannot, install Ollama itself (that is
 * a real system-level installer and is intentionally left to the user via the official link in the
 * UI — never silently automated).
 *
 * Only one pull runs at a time (a second call replaces the tracked status, matching this being a
 * single, small settings-panel affordance rather than a queueing system).
 */
public class OllamaPuller {

    public static class PullStatus {
        public volatile String status = "idle";
        public volatile long completed = 0;
        public volatile long total = 0;
        public volatile boolean done = true;
        public volatile String error = null;
    }

    private static final PullStatus CURRENT = new PullStatus();

    public static PullStatus current() { return CURRENT; }

    /** Starts a background pull; returns immediately. Progress is polled via {@link #current()}. */
    public static void startPull(String ollamaBaseUrl, String model) {
        synchronized (CURRENT) {
            CURRENT.status = "starting";
            CURRENT.completed = 0;
            CURRENT.total = 0;
            CURRENT.done = false;
            CURRENT.error = null;
        }
        new Thread(() -> runPull(ollamaBaseUrl, model), "ollama-pull-" + model).start();
    }

    private static void runPull(String ollamaBaseUrl, String model) {
        try {
            HttpClient client = HttpClient.newBuilder().connectTimeout(Duration.ofSeconds(5)).build();
            String body = "{\"name\":" + jsonStr(model) + ",\"stream\":true}";
            HttpRequest req = HttpRequest.newBuilder(URI.create(ollamaBaseUrl + "/api/pull"))
                .timeout(Duration.ofMinutes(30))
                .header("Content-Type", "application/json")
                .POST(HttpRequest.BodyPublishers.ofString(body, StandardCharsets.UTF_8))
                .build();
            HttpResponse<Stream<String>> resp = client.send(req, HttpResponse.BodyHandlers.ofLines());
            if (resp.statusCode() / 100 != 2) {
                CURRENT.error = "Ollama returned HTTP " + resp.statusCode();
                CURRENT.done = true;
                return;
            }
            resp.body().forEach(OllamaPuller::applyLine);
            synchronized (CURRENT) {
                if (CURRENT.error == null && !"success".equals(CURRENT.status)) {
                    // Stream ended without an explicit "success" or "error" line — treat as done anyway
                    // rather than leaving the UI stuck at done=false forever.
                }
                CURRENT.done = true;
            }
        } catch (Exception e) {
            CURRENT.error = e.getMessage() == null ? e.getClass().getSimpleName() : e.getMessage();
            CURRENT.done = true;
        }
    }

    private static void applyLine(String line) {
        if (line == null || line.isBlank()) return;
        try {
            Map<String, Object> obj = MiniJson.asObject(MiniJson.parse(line));
            if (obj.get("status") != null) CURRENT.status = String.valueOf(obj.get("status"));
            // Ollama's stream keeps sending lines after the model layer itself is fully downloaded
            // ("verifying sha256 digest", "writing manifest", "removing any unused layers", "success")
            // and some of those carry their own small, unrelated total/completed pair (e.g. the
            // manifest file's own byte count) — applying those would stomp the real download's
            // progress right at the end (observed directly: a real pull's final reported pct came out
            // as a nonsensical multi-hundred-million percent). Only "pulling ..." lines describe the
            // actual model download, so only those are allowed to update total/completed.
            String status = CURRENT.status;
            boolean isDownloadLine = status != null && status.startsWith("pulling");
            if (isDownloadLine && obj.get("total") != null) CURRENT.total = ((Number) obj.get("total")).longValue();
            if (isDownloadLine && obj.get("completed") != null) CURRENT.completed = ((Number) obj.get("completed")).longValue();
            if (obj.get("error") != null) CURRENT.error = String.valueOf(obj.get("error"));
        } catch (Exception ignore) {
            // A non-JSON or partial line should never crash the whole pull — skip it, the next line
            // (or the final status) still gets through.
        }
    }

    private static String jsonStr(String s) {
        StringBuilder sb = new StringBuilder("\"");
        for (int i = 0; i < s.length(); i++) {
            char c = s.charAt(i);
            if (c == '"') sb.append("\\\"");
            else if (c == '\\') sb.append("\\\\");
            else sb.append(c);
        }
        return sb.append("\"").toString();
    }
}
