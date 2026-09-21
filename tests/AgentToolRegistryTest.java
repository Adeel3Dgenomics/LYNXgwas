import com.sun.net.httpserver.*;
import java.io.*;
import java.net.InetSocketAddress;
import java.nio.charset.StandardCharsets;
import java.util.*;
import java.util.concurrent.ConcurrentHashMap;

/**
 * Regression test for AgentToolRegistry (DECISIONS_PHASE4.md section 2.3), against a small,
 * purpose-built fake HTTP backend standing in for LocalServer's own project-scoped endpoints (bound
 * to an ephemeral port so this can never collide with a real, already-running LYNXgwas server) —
 * fast and deterministic, and lets this test assert exactly what AgentToolRegistry sent, which a real
 * server's real side effects would make much harder to inspect directly.
 *
 * The single most important thing this test verifies is the merge-before-write behavior described in
 * AgentToolRegistry's own class javadoc: updating one field of an *existing* project's config must
 * not silently reset every other field to Config's hard-coded defaults.
 */
public class AgentToolRegistryTest {

    private static final Map<String, String> capturedBodies = new ConcurrentHashMap<>();

    public static void main(String[] progArgs) throws Exception {
        int failures = 0;
        HttpServer server = HttpServer.create(new InetSocketAddress("localhost", 0), 0);
        int port = server.getAddress().getPort();
        registerFakeRoutes(server);
        server.start();
        try {
            AgentToolRegistry registry = new AgentToolRegistry(port);

            // ── list_projects ──
            AgentToolRegistry.ToolExecutionResult r = registry.dispatch("list_projects", new HashMap<>());
            failures += check("list_projects ok", r.ok);
            failures += check("list_projects summary mentions count", r.summary.contains("2"));

            // ── get_project_config: existing project ──
            Map<String, Object> callArgs = new LinkedHashMap<>();
            callArgs.put("project_id", "existing");
            r = registry.dispatch("get_project_config", callArgs);
            failures += check("get_project_config(existing) ok", r.ok);
            failures += check("get_project_config(existing) returns gwas.file", r.resultJson.contains("input/x.tsv"));

            // ── get_project_config: missing required arg ──
            r = registry.dispatch("get_project_config", new HashMap<>());
            failures += check("get_project_config() without project_id fails cleanly (no exception escapes)", !r.ok);

            // ── create_or_update_project: brand-new project (no prior config) ──
            callArgs = new LinkedHashMap<>();
            callArgs.put("project_id", "brandnew");
            callArgs.put("gwas.file", "input/new.tsv");
            callArgs.put("threads", "8");
            r = registry.dispatch("create_or_update_project", callArgs);
            failures += check("create_or_update_project(brandnew) ok", r.ok);
            failures += check("create_or_update_project(brandnew) summary says Created",
                r.summary.startsWith("Created project 'brandnew'"));
            String postedNew = capturedBodies.get("POST /api/project/brandnew/config");
            failures += check("brandnew POST body contains gwas.file",
                postedNew != null && postedNew.contains("input/new.tsv"));
            failures += check("brandnew POST body contains threads",
                postedNew != null && postedNew.contains("\"threads\":\"8\""));

            // ── create_or_update_project: EXISTING project, changing only one field ──
            // "existing" fake-backend config has gwas.file=input/x.tsv, col.chr=chrom, threads=4,
            // ld.enabled=false. We change only threads. The other three fields MUST survive in the
            // POSTed body — this is the merge-before-write behavior the whole class exists to get right.
            callArgs = new LinkedHashMap<>();
            callArgs.put("project_id", "existing");
            callArgs.put("threads", "16");
            r = registry.dispatch("create_or_update_project", callArgs);
            failures += check("create_or_update_project(existing, 1 field) ok", r.ok);
            failures += check("create_or_update_project(existing) summary says Updated",
                r.summary.startsWith("Updated project 'existing'"));
            String postedExisting = capturedBodies.get("POST /api/project/existing/config");
            failures += check("merge preserves untouched gwas.file",
                postedExisting != null && postedExisting.contains("input/x.tsv"));
            failures += check("merge preserves untouched col.chr",
                postedExisting != null && postedExisting.contains("\"col.chr\":\"chrom\""));
            failures += check("merge applies the new threads value",
                postedExisting != null && postedExisting.contains("\"threads\":\"16\""));
            failures += check("merge does not still contain the old threads value",
                postedExisting != null && !postedExisting.contains("\"threads\":\"4\""));

            // ── run_pipeline ──
            callArgs = new LinkedHashMap<>();
            callArgs.put("project_id", "existing");
            r = registry.dispatch("run_pipeline", callArgs);
            failures += check("run_pipeline ok", r.ok);
            String processBody = capturedBodies.get("POST /api/projects/process");
            failures += check("run_pipeline posted the right project id",
                processBody != null && processBody.contains("\"id\":\"existing\""));

            // ── get_pipeline_progress ──
            r = registry.dispatch("get_pipeline_progress", callArgs);
            failures += check("get_pipeline_progress ok", r.ok);

            // ── list_analysis_tools ──
            r = registry.dispatch("list_analysis_tools", new HashMap<>());
            failures += check("list_analysis_tools ok", r.ok);
            failures += check("list_analysis_tools returns susie", r.resultJson.contains("susie"));

            // ── run_analysis_tool ──
            callArgs = new LinkedHashMap<>();
            callArgs.put("project_id", "existing");
            callArgs.put("tool", "susie");
            callArgs.put("locus_id", "1");
            r = registry.dispatch("run_analysis_tool", callArgs);
            failures += check("run_analysis_tool ok", r.ok);
            String analysisBody = capturedBodies.get("POST /api/project/existing/analysis/run");
            failures += check("run_analysis_tool posted tool=susie",
                analysisBody != null && analysisBody.contains("\"tool\":\"susie\""));
            failures += check("run_analysis_tool did not post project_id in the body (only in the path)",
                analysisBody != null && !analysisBody.contains("project_id"));

            // ── run_analysis_tool: missing required arg ──
            r = registry.dispatch("run_analysis_tool", new HashMap<>());
            failures += check("run_analysis_tool() with no args fails cleanly", !r.ok);

            // ── get_project_manifest ──
            callArgs = new LinkedHashMap<>();
            callArgs.put("project_id", "existing");
            r = registry.dispatch("get_project_manifest", callArgs);
            failures += check("get_project_manifest ok", r.ok);

            // ── unknown tool name ──
            r = registry.dispatch("delete_everything", new HashMap<>());
            failures += check("unknown tool name is rejected, not silently ignored", !r.ok);

            // ── tool schemas are well-formed ──
            List<Map<String, Object>> schemas = registry.toolSchemas();
            failures += check("exposes exactly 8 tools", schemas.size() == 8);
            boolean anyDelete = false;
            for (Map<String, Object> s : schemas) {
                Map<?, ?> fn = (Map<?, ?>) s.get("function");
                if (String.valueOf(fn.get("name")).contains("delete")) anyDelete = true;
            }
            failures += check("no delete-project tool is exposed (deliberate scope limit)", !anyDelete);

        } finally {
            server.stop(0);
        }

        if (failures == 0) {
            System.out.println("PASS: all AgentToolRegistry tests passed");
        } else {
            System.out.println("FAIL: " + failures + " AgentToolRegistry test(s) failed");
            System.exit(1);
        }
    }

    private static void registerFakeRoutes(HttpServer server) {
        server.createContext("/api/projects", ex -> {
            respond(ex, 200, "[{\"id\":\"existing\"},{\"id\":\"brandnew\"}]");
        });
        server.createContext("/api/projects/process", ex -> {
            String body = readBody(ex);
            capturedBodies.put("POST /api/projects/process", body);
            respond(ex, 200, "{\"processing\":1}");
        });
        server.createContext("/api/project/", ex -> {
            String path = ex.getRequestURI().getPath();
            String method = ex.getRequestMethod();
            if (path.equals("/api/project/existing/config") && method.equals("GET")) {
                respond(ex, 200, "{\"gwas.file\":\"input/x.tsv\",\"col.chr\":\"chrom\",\"threads\":4,\"ld.enabled\":false}");
            } else if (path.equals("/api/project/brandnew/config") && method.equals("GET")) {
                respond(ex, 404, "{\"error\":\"No config.properties\"}");
            } else if (path.endsWith("/config") && method.equals("POST")) {
                String body = readBody(ex);
                capturedBodies.put("POST " + path, body);
                respond(ex, 200, "{\"ok\":true}");
            } else if (path.equals("/api/project/existing/progress") && method.equals("GET")) {
                respond(ex, 200, "{\"pct\":50,\"done\":false}");
            } else if (path.equals("/api/project/_agent/analysis/tools") && method.equals("GET")) {
                respond(ex, 200, "[{\"name\":\"susie\"},{\"name\":\"coloc\"}]");
            } else if (path.equals("/api/project/existing/analysis/run") && method.equals("POST")) {
                String body = readBody(ex);
                capturedBodies.put("POST " + path, body);
                respond(ex, 200, "{\"job_id\":\"abc123\"}");
            } else if (path.equals("/api/project/existing/manifest") && method.equals("GET")) {
                respond(ex, 200, "{\"loci\":[]}");
            } else {
                respond(ex, 404, "{\"error\":\"no fake route for " + path + "\"}");
            }
        });
    }

    private static String readBody(HttpExchange ex) throws IOException {
        ByteArrayOutputStream buf = new ByteArrayOutputStream();
        byte[] chunk = new byte[4096];
        int n;
        try (InputStream in = ex.getRequestBody()) {
            while ((n = in.read(chunk)) >= 0) buf.write(chunk, 0, n);
        }
        return buf.toString(StandardCharsets.UTF_8);
    }

    private static void respond(HttpExchange ex, int status, String json) throws IOException {
        byte[] bytes = json.getBytes(StandardCharsets.UTF_8);
        ex.getResponseHeaders().add("Content-Type", "application/json");
        ex.sendResponseHeaders(status, bytes.length);
        try (OutputStream os = ex.getResponseBody()) { os.write(bytes); }
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
