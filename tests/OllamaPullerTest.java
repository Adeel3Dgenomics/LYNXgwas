import com.sun.net.httpserver.*;
import java.io.*;
import java.net.InetSocketAddress;
import java.nio.charset.StandardCharsets;

/**
 * Regression test for OllamaPuller, against a small fake HTTP server that streams a realistic
 * newline-delimited-JSON {@code /api/pull} response (manifest -> incrementally growing progress ->
 * success), plus a separate error case (Ollama reports a bad model name as a JSON {"error":...} line
 * within an HTTP 200 response, not as an HTTP error status). Bound to an ephemeral port so this can
 * never collide with a real, already-running Ollama instance.
 */
public class OllamaPullerTest {

    public static void main(String[] progArgs) throws Exception {
        int failures = 0;

        // ── Successful pull, with incremental progress ──
        HttpServer ok = HttpServer.create(new InetSocketAddress("localhost", 0), 0);
        int okPort = ok.getAddress().getPort();
        ok.createContext("/api/pull", ex -> {
            ex.getResponseHeaders().add("Content-Type", "application/x-ndjson");
            ex.sendResponseHeaders(200, 0); // chunked, unknown length
            try (OutputStream os = ex.getResponseBody()) {
                // Status strings match what a real Ollama instance actually sends (confirmed against
                // a real pull of llama3.2:1b), not a made-up shape — see the trailing bookkeeping
                // lines below, which is exactly what previously corrupted the final reported progress.
                writeLine(os, "{\"status\":\"pulling manifest\"}");
                sleep(50);
                writeLine(os, "{\"status\":\"pulling abc123digest\",\"total\":1000,\"completed\":250}");
                sleep(50);
                writeLine(os, "{\"status\":\"pulling abc123digest\",\"total\":1000,\"completed\":1000}");
                sleep(50);
                // Trailing bookkeeping lines, each with their own small, unrelated total/completed —
                // these must NOT overwrite the real download's progress (the bug this regression
                // guards against: a real pull's final reported pct came out as a nonsensical
                // multi-hundred-million percent because a later small total/completed pair stomped
                // the real ones).
                writeLine(os, "{\"status\":\"verifying sha256 digest\"}");
                writeLine(os, "{\"status\":\"writing manifest\",\"total\":485,\"completed\":485}");
                writeLine(os, "{\"status\":\"removing any unused layers\"}");
                writeLine(os, "{\"status\":\"success\"}");
            }
        });
        ok.start();
        try {
            OllamaPuller.startPull("http://localhost:" + okPort, "fake-model");

            // Poll until done (bounded wait so a real bug can't hang the test suite forever).
            OllamaPuller.PullStatus last = waitUntilDone(5000);
            failures += check("pull completed (done=true within timeout)", last.done);
            failures += check("no error on a successful pull", last.error == null);
            failures += check("final status is 'success'", "success".equals(last.status));
            failures += check("final completed/total reflect the real download, not a trailing " +
                "bookkeeping line's unrelated small total/completed (485)",
                last.completed == 1000 && last.total == 1000);
        } finally {
            ok.stop(0);
        }

        // ── Error case: HTTP 200 but the stream contains an {"error":...} line (bad model name) ──
        HttpServer err = HttpServer.create(new InetSocketAddress("localhost", 0), 0);
        int errPort = err.getAddress().getPort();
        err.createContext("/api/pull", ex -> {
            ex.getResponseHeaders().add("Content-Type", "application/x-ndjson");
            ex.sendResponseHeaders(200, 0);
            try (OutputStream os = ex.getResponseBody()) {
                writeLine(os, "{\"status\":\"pulling manifest\"}");
                writeLine(os, "{\"error\":\"model 'nonexistent-model' not found\"}");
            }
        });
        err.start();
        try {
            OllamaPuller.startPull("http://localhost:" + errPort, "nonexistent-model");
            OllamaPuller.PullStatus last = waitUntilDone(5000);
            failures += check("errored pull still reaches done=true", last.done);
            failures += check("error message captured", last.error != null && last.error.contains("not found"));
        } finally {
            err.stop(0);
        }

        // ── HTTP-level failure (server down / non-2xx) ──
        OllamaPuller.startPull("http://localhost:1", "whatever"); // nothing listens on port 1
        OllamaPuller.PullStatus last = waitUntilDone(5000);
        failures += check("connection failure still reaches done=true (does not hang)", last.done);
        failures += check("connection failure is reported as an error", last.error != null);

        if (failures == 0) {
            System.out.println("PASS: all OllamaPuller tests passed");
        } else {
            System.out.println("FAIL: " + failures + " OllamaPuller test(s) failed");
            System.exit(1);
        }
    }

    private static OllamaPuller.PullStatus waitUntilDone(long timeoutMs) throws InterruptedException {
        long deadline = System.currentTimeMillis() + timeoutMs;
        OllamaPuller.PullStatus s = OllamaPuller.current();
        while (!s.done && System.currentTimeMillis() < deadline) Thread.sleep(20);
        return s;
    }

    private static void writeLine(OutputStream os, String json) throws IOException {
        os.write((json + "\n").getBytes(StandardCharsets.UTF_8));
        os.flush();
    }

    private static void sleep(long ms) {
        try { Thread.sleep(ms); } catch (InterruptedException ignored) {}
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
