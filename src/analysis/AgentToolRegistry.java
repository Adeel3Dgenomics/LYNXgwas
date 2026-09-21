import java.io.*;
import java.net.*;
import java.net.http.*;
import java.nio.charset.StandardCharsets;
import java.time.Duration;
import java.util.*;

/**
 * The agent's fixed, fully-transparent tool surface (DECISIONS_PHASE4.md sections 1/2.3). Every tool
 * here is a thin loop-back HTTP call to this same server's own existing REST API — the agent can
 * therefore never do anything a human user couldn't already do by clicking through the UI; this is
 * automation of the existing surface, not a new one. Deliberately no delete-project tool.
 *
 * One important correctness detail this class owns rather than leaving to the LLM: {@code
 * POST /api/project/{id}/config} applies its JSON body onto a *fresh* {@code Config} object (see
 * {@code Config.applyJson}), so any field the caller omits reverts to that field's hard-coded
 * default — it is not a merge against whatever is already on disk. A naive "change one field" tool
 * call would therefore silently reset every other field of an existing project. {@link
 * #createOrUpdateProject} avoids this by reading the project's current config first (if it already
 * exists) and only overlaying the fields the model actually supplied before writing back the full,
 * merged document — exactly what a human editing the same field through the existing config-edit UI
 * form would submit.
 */
public class AgentToolRegistry {

    private final int serverPort;
    private final HttpClient http = HttpClient.newBuilder().connectTimeout(Duration.ofSeconds(5)).build();

    public AgentToolRegistry(int serverPort) {
        this.serverPort = serverPort;
    }

    public static class ToolExecutionResult {
        public final boolean ok;
        public final String resultJson; // fed back to the LLM as the tool message content
        public final String summary;    // one-line, shown to the user as a "step"
        public ToolExecutionResult(boolean ok, String resultJson, String summary) {
            this.ok = ok; this.resultJson = resultJson; this.summary = summary;
        }
    }

    // ── Tool schemas (OpenAI function-calling "tools" array shape) ──

    public List<Map<String, Object>> toolSchemas() {
        List<Map<String, Object>> tools = new ArrayList<>();
        tools.add(fn("list_projects",
            "List every project already known to LYNXgwas, with basic status (loci count, whether it needs reprocessing).",
            params()));
        tools.add(fn("get_project_config",
            "Read a project's current configuration (input files, column mapping, LD settings, dataset metadata).",
            params(req("project_id", "string", "The project's id (its folder name under projects/).")) ));
        tools.add(fn("create_or_update_project",
            "Create a new project, or change one or more settings of an existing one. Only include the fields you " +
            "want to set; anything you omit keeps its current value on an existing project, or the sensible default " +
            "on a new one. Field names match LYNXgwas's own config.properties keys exactly (e.g. 'gwas.file', " +
            "'col.chr', 'ref.panel.path', 'ld.enabled').",
            paramsMap(reqStr("project_id", "The project's id (folder name); created if it doesn't exist yet."),
                optStr("gwas.file", "Path to the raw GWAS summary-statistics file."),
                optStr("loci.file", "Path to a pre-defined loci list, if not using automatic clumping."),
                optStr("gff3.file", "Path to the GENCODE GFF3 gene-annotation file."),
                optStr("ref.panel.path", "Path (prefix) to the PLINK reference panel for LD computation."),
                optStr("ref.panel.population", "Reference panel population code, e.g. EUR, EAS, AFR."),
                optStr("col.chr", "GWAS file column name for chromosome."),
                optStr("col.pos", "GWAS file column name for base-pair position."),
                optStr("col.pvalue", "GWAS file column name for the p-value."),
                optStr("col.rsid", "GWAS file column name for rsID (blank if not present)."),
                optStr("col.varid", "GWAS file column name for a variant id, if any."),
                optStr("col.ea", "GWAS file column name for the effect allele."),
                optStr("col.nea", "GWAS file column name for the non-effect allele."),
                optStr("col.beta", "GWAS file column name for the effect size (beta), if present."),
                optStr("col.or", "GWAS file column name for the odds ratio, if present instead of beta."),
                optStr("col.se", "GWAS file column name for the standard error, if present."),
                optStr("col.n", "GWAS file column name for per-SNP sample size, if present."),
                optStr("ld.enabled", "\"true\" or \"false\" — whether to compute LD for this project."),
                optStr("threads", "Number of CPU threads to use, as a string, e.g. \"4\"."),
                optStr("sample.n", "Total GWAS sample size, as a string."),
                optStr("n.cases", "Number of cases, as a string, for a case/control study."),
                optStr("n.controls", "Number of controls, as a string, for a case/control study."),
                optStr("trait.type", "\"quantitative\" or \"binary\"."),
                optStr("effect.type", "\"beta\", \"OR\", or \"logOR\" — what col.beta/col.or actually contains."),
                optStr("genome.build", "\"GRCh37\" or \"GRCh38\"."),
                optStr("ancestry", "Free-text ancestry label for this dataset."),
                optStr("disease.name", "Free-text trait/disease name."))));
        tools.add(fn("run_pipeline",
            "Start (or restart) the base pipeline for a project: loci discovery, rsID recovery, LD computation. " +
            "Runs in the background; use get_pipeline_progress to check on it.",
            params(req("project_id", "string", "The project's id."))));
        tools.add(fn("get_pipeline_progress",
            "Check how far the base pipeline has gotten for a project that run_pipeline was called on.",
            params(req("project_id", "string", "The project's id."))));
        tools.add(fn("list_analysis_tools",
            "List the fine-mapping/analysis tools available to run on a locus (SuSiE, FINEMAP-style, GCTA-COJO, " +
            "colocalization, GWAMA, MAGMA, GCTA-GREML), with their parameters.",
            params()));
        tools.add(fn("run_analysis_tool",
            "Run one analysis tool (see list_analysis_tools for valid tool names and parameters) on one locus of " +
            "an already-processed project. Runs in the background.",
            paramsMixed(
                req("project_id", "string", "The project's id."),
                req("tool", "string", "The tool name, exactly as returned by list_analysis_tools."),
                req("locus_id", "string", "The locus id to run this tool on."))));
        tools.add(fn("get_project_manifest",
            "Read a project's manifest: its list of loci (position, top SNP, nearest gene) once the base pipeline " +
            "has finished.",
            params(req("project_id", "string", "The project's id."))));
        return tools;
    }

    // ── Dispatch ──

    public ToolExecutionResult dispatch(String toolName, Map<String, Object> args) {
        try {
            switch (toolName) {
                case "list_projects":              return listProjects();
                case "get_project_config":          return getProjectConfig(str(args, "project_id"));
                case "create_or_update_project":    return createOrUpdateProject(args);
                case "run_pipeline":                return runPipeline(str(args, "project_id"));
                case "get_pipeline_progress":       return getPipelineProgress(str(args, "project_id"));
                case "list_analysis_tools":          return listAnalysisTools();
                case "run_analysis_tool":            return runAnalysisTool(args);
                case "get_project_manifest":         return getProjectManifest(str(args, "project_id"));
                default:
                    return new ToolExecutionResult(false,
                        "{\"error\":\"Unknown tool: " + toolName + "\"}",
                        "Unknown tool requested: " + toolName);
            }
        } catch (Exception e) {
            String msg = e.getMessage() == null ? e.getClass().getSimpleName() : e.getMessage();
            return new ToolExecutionResult(false,
                "{\"error\":" + jsonStr(msg) + "}",
                "Tool '" + toolName + "' failed: " + msg);
        }
    }

    private ToolExecutionResult listProjects() throws IOException, InterruptedException {
        HttpResp r = get("/api/projects");
        int n = countTopLevelObjects(r.body);
        return new ToolExecutionResult(r.ok, r.body, "Listed projects (" + n + " found).");
    }

    private ToolExecutionResult getProjectConfig(String projectId) throws IOException, InterruptedException {
        requireStr(projectId, "project_id");
        HttpResp r = get("/api/project/" + enc(projectId) + "/config");
        return new ToolExecutionResult(r.ok, r.body,
            r.ok ? "Read config for project '" + projectId + "'."
                 : "No existing config for project '" + projectId + "' (it may not exist yet).");
    }

    @SuppressWarnings("unchecked")
    private ToolExecutionResult createOrUpdateProject(Map<String, Object> args) throws IOException, InterruptedException {
        String projectId = str(args, "project_id");
        requireStr(projectId, "project_id");

        // Merge onto the current on-disk config (if any) so fields the model didn't mention survive —
        // see the class-level javadoc for why this matters.
        Map<String, Object> merged = new LinkedHashMap<>();
        HttpResp current = get("/api/project/" + enc(projectId) + "/config");
        boolean isNew = !current.ok;
        if (current.ok) {
            Object parsed = MiniJson.parse(current.body);
            merged.putAll(MiniJson.asObject(parsed));
        }
        for (Map.Entry<String, Object> e : args.entrySet()) {
            if (e.getKey().equals("project_id")) continue;
            merged.put(e.getKey(), e.getValue());
        }

        String body = encodeFlatObject(merged);
        HttpResp r = post("/api/project/" + enc(projectId) + "/config", body);
        int changed = args.size() - 1;
        String summary = (isNew ? "Created project '" : "Updated project '") + projectId + "' ("
            + Math.max(changed, 0) + " field(s) set).";
        return new ToolExecutionResult(r.ok, r.body, summary);
    }

    private ToolExecutionResult runPipeline(String projectId) throws IOException, InterruptedException {
        requireStr(projectId, "project_id");
        HttpResp r = post("/api/projects/process", "{\"id\":" + jsonStr(projectId) + "}");
        return new ToolExecutionResult(r.ok, r.body,
            r.ok ? "Started pipeline processing for '" + projectId + "'."
                 : "Could not start pipeline for '" + projectId + "': " + r.body);
    }

    private ToolExecutionResult getPipelineProgress(String projectId) throws IOException, InterruptedException {
        requireStr(projectId, "project_id");
        HttpResp r = get("/api/project/" + enc(projectId) + "/progress");
        return new ToolExecutionResult(r.ok, r.body, "Checked pipeline progress for '" + projectId + "'.");
    }

    private ToolExecutionResult listAnalysisTools() throws IOException, InterruptedException {
        // This handler ignores the project id entirely (it lists globally-available tool descriptors),
        // so any placeholder segment works — see LocalServer's own comment at its route registration.
        HttpResp r = get("/api/project/_agent/analysis/tools");
        return new ToolExecutionResult(r.ok, r.body, "Listed available analysis tools.");
    }

    private ToolExecutionResult runAnalysisTool(Map<String, Object> args) throws IOException, InterruptedException {
        String projectId = str(args, "project_id");
        String tool = str(args, "tool");
        String locusId = str(args, "locus_id");
        requireStr(projectId, "project_id");
        requireStr(tool, "tool");
        requireStr(locusId, "locus_id");

        Map<String, Object> body = new LinkedHashMap<>(args);
        body.remove("project_id");
        HttpResp r = post("/api/project/" + enc(projectId) + "/analysis/run", encodeFlatObject(body));
        return new ToolExecutionResult(r.ok, r.body,
            r.ok ? "Started '" + tool + "' on locus '" + locusId + "' for project '" + projectId + "'."
                 : "Could not start '" + tool + "': " + r.body);
    }

    private ToolExecutionResult getProjectManifest(String projectId) throws IOException, InterruptedException {
        requireStr(projectId, "project_id");
        HttpResp r = get("/api/project/" + enc(projectId) + "/manifest");
        return new ToolExecutionResult(r.ok, r.body,
            r.ok ? "Read manifest for project '" + projectId + "'."
                 : "No manifest yet for '" + projectId + "' (run_pipeline may not have completed).");
    }

    // ── Loop-back HTTP plumbing ──

    private static class HttpResp {
        final boolean ok;
        final String body;
        HttpResp(boolean ok, String body) { this.ok = ok; this.body = body; }
    }

    private HttpResp get(String path) throws IOException, InterruptedException {
        HttpRequest req = HttpRequest.newBuilder(URI.create("http://localhost:" + serverPort + path))
            .timeout(Duration.ofSeconds(30)).GET().build();
        HttpResponse<String> resp = http.send(req, HttpResponse.BodyHandlers.ofString(StandardCharsets.UTF_8));
        return new HttpResp(resp.statusCode() / 100 == 2, resp.body());
    }

    private HttpResp post(String path, String jsonBody) throws IOException, InterruptedException {
        HttpRequest req = HttpRequest.newBuilder(URI.create("http://localhost:" + serverPort + path))
            .timeout(Duration.ofSeconds(30))
            .header("Content-Type", "application/json")
            .POST(HttpRequest.BodyPublishers.ofString(jsonBody, StandardCharsets.UTF_8)).build();
        HttpResponse<String> resp = http.send(req, HttpResponse.BodyHandlers.ofString(StandardCharsets.UTF_8));
        return new HttpResp(resp.statusCode() / 100 == 2, resp.body());
    }

    private static String enc(String s) {
        return URLEncoder.encode(s, StandardCharsets.UTF_8).replace("+", "%20");
    }

    // ── Small helpers ──

    private static String str(Map<String, Object> args, String key) {
        Object v = args.get(key);
        return v == null ? null : String.valueOf(v);
    }

    private static void requireStr(String v, String name) {
        if (v == null || v.isEmpty()) throw new IllegalArgumentException("Missing required argument: " + name);
    }

    private static int countTopLevelObjects(String jsonArray) {
        try {
            List<Object> list = MiniJson.asArray(MiniJson.parse(jsonArray));
            return list.size();
        } catch (Exception e) {
            return -1;
        }
    }

    /** Config.applyJson's own jsonStr() extractor reads a value's inner text whether it's quoted or
     *  a bare numeric/boolean token, so MiniJson's general encoder produces a document it reads
     *  correctly either way. */
    private static String encodeFlatObject(Map<String, Object> m) {
        return MiniJson.encode(m);
    }

    private static String jsonStr(String s) {
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

    // ── Tool-schema construction helpers ──

    private static Map<String, Object> fn(String name, String description, Map<String, Object> parameters) {
        Map<String, Object> function = new LinkedHashMap<>();
        function.put("name", name);
        function.put("description", description);
        function.put("parameters", parameters);
        Map<String, Object> wrapper = new LinkedHashMap<>();
        wrapper.put("type", "function");
        wrapper.put("function", function);
        return wrapper;
    }

    private static class P { String name, type, desc; boolean required; }

    private static P req(String name, String type, String desc) {
        P p = new P(); p.name = name; p.type = type; p.desc = desc; p.required = true; return p;
    }
    private static P reqStr(String name, String desc) { return req(name, "string", desc); }
    private static P optStr(String name, String desc) {
        P p = new P(); p.name = name; p.type = "string"; p.desc = desc; p.required = false; return p;
    }

    private static Map<String, Object> params(P... ps) { return paramsMixed(ps); }
    private static Map<String, Object> paramsMap(P... ps) { return paramsMixed(ps); }

    private static Map<String, Object> paramsMixed(P... ps) {
        Map<String, Object> properties = new LinkedHashMap<>();
        List<Object> required = new ArrayList<>();
        for (P p : ps) {
            Map<String, Object> prop = new LinkedHashMap<>();
            prop.put("type", p.type);
            prop.put("description", p.desc);
            properties.put(p.name, prop);
            if (p.required) required.add(p.name);
        }
        Map<String, Object> schema = new LinkedHashMap<>();
        schema.put("type", "object");
        schema.put("properties", properties);
        if (!required.isEmpty()) schema.put("required", required);
        return schema;
    }
}
