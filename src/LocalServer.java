import com.sun.net.httpserver.*;
import javax.swing.*;
import java.io.*;
import java.net.*;
import java.nio.file.*;
import java.util.*;
import java.util.concurrent.*;
import rsid.*;
import loci.*;
import export.*;

public class LocalServer {

    public static final int PORT = 8765;

    /** Global progress — used by legacy single-project mode. */
    public static final ProgressTracker progress = new ProgressTracker();

    private final HttpServer http;
    private final String     outputDir;
    private volatile String  lastFolder;

    // Legacy single-project state (set by Main after pipeline)
    private Config             config;
    private GffParser          gff;
    private List<Locus>        loci;
    private List<LocusOutput>  outputs;

    // Per-project progress for background processing
    private final Map<String, ProgressTracker> projectProgress = new ConcurrentHashMap<>();
    // Per-project pipeline state for live locus updates
    private final Map<String, ProjectState> projectStates = new ConcurrentHashMap<>();
    // Projects currently being processed (prevent double-processing)
    private final Set<String> processing = ConcurrentHashMap.newKeySet();

    static class ProjectState {
        Config config;
        GffParser gff;
        List<Locus> loci;
        List<LocusOutput> outputs;
        LociMutationService mutationService;
    }

    public void setPipelineState(Config config, GffParser gff, List<Locus> loci, List<LocusOutput> outputs) {
        this.config  = config;
        this.gff     = gff;
        this.loci    = loci;
        this.outputs = outputs;
    }

    public void setProjectState(String projectId, Config config, GffParser gff,
                                List<Locus> loci, List<LocusOutput> outputs) {
        ProjectState ps = new ProjectState();
        ps.config = config; ps.gff = gff; ps.loci = loci; ps.outputs = outputs;
        ps.mutationService = new LociMutationService(config, gff, loci, outputs);
        projectStates.put(projectId, ps);
    }

    public LocalServer(String outputDir) throws IOException {
        this.outputDir  = outputDir;
        this.lastFolder = outputDir + "/plots";
        new File(this.lastFolder).mkdirs();

        http = HttpServer.create(new InetSocketAddress("localhost", PORT), 32);

        // ── Project-scoped API (Step 6) ──────────────────────────────────
        http.createContext("/api/projects/process", this::processProjects);
        http.createContext("/api/projects",         this::listProjects);
        http.createContext("/api/delete-project",    this::projectDelete);
        http.createContext("/api/project/",         this::projectRouter);
        http.createContext("/api/peek-file-header", this::peekFileHeader);

        // ── rsID recovery API ────────────────────────────────────────────
        http.createContext("/api/global-config",     this::globalConfigEndpoint);
        http.createContext("/api/rsid-recover",      this::rsidRecover);
        http.createContext("/api/rsid-progress",     this::rsidProgressEndpoint);
        http.createContext("/api/rsid-detect",       this::rsidDetect);
        http.createContext("/api/loci-identify",     this::lociIdentify);
        http.createContext("/api/loci-progress",     this::lociProgressEndpoint);
        http.createContext("/api/export-excel",      this::exportExcel);
        http.createContext("/api/export-progress",   this::exportProgressEndpoint);
        http.createContext("/pick-folder-native",    this::pickFolder);

        // ── Legacy endpoints (kept for backward compatibility) ───────────
        http.createContext("/manifest",                    this::manifest);
        http.createContext("/locus/",                      this::locus);
        http.createContext("/pick-folder",                 this::pickFolder);
        http.createContext("/save-pdf",                    this::savePdf);
        http.createContext("/progress",                    this::progressEndpoint);
        http.createContext("/api/update-locus",            this::updateLocus);
        http.createContext("/api/create-locus",            this::createLocus);
        http.createContext("/api/split-locus",             this::splitLocus);
        http.createContext("/api/validate-split",          this::validateSplit);
        http.createContext("/pick-file",                   this::pickFile);
        http.createContext("/api/annotation-file",         this::annotationFile);
        http.createContext("/api/save-annotation-config",  this::saveAnnotationConfig);
        http.createContext("/",                             this::staticFiles);
        http.setExecutor(Executors.newFixedThreadPool(8));
    }

    public void start() { http.start(); }
    public void stop()  { http.stop(0); }

    // ══════════════════════════════════════════════════════════════════════
    //  PROJECT-SCOPED ENDPOINTS (Step 6)
    // ══════════════════════════════════════════════════════════════════════

    // GET /api/projects — list all projects with live-computed status
    private void listProjects(HttpExchange ex) throws IOException {
        cors(ex); if (preflight(ex)) return;
        List<File> dirs = Main.discoverProjects(new File("projects"));
        StringBuilder json = new StringBuilder("[");
        boolean first = true;
        for (File dir : dirs) {
            if (!first) json.append(',');
            first = false;
            json.append('{');
            String id = dir.getName();
            json.append("\"id\":\"").append(escJ(id)).append('"');

            ProjectMetadata pm = ProjectMetadata.load(dir.getAbsolutePath());
            if (pm != null) {
                json.append(",\"name\":\"").append(escJ(pm.name)).append('"');
                json.append(",\"description\":\"").append(escJ(pm.description)).append('"');
                json.append(",\"loci_count\":").append(pm.lociCount);
                json.append(",\"total_snps\":").append(pm.totalSnps);
                json.append(",\"snp_annotation_sources\":").append(pm.snpAnnotationSources);
                json.append(",\"locus_annotation_sources\":").append(pm.locusAnnotationSources);
                json.append(",\"snps_with_any_annotation\":").append(pm.snpsWithAnyAnnotation);
                json.append(",\"last_processed\":\"").append(escJ(pm.lastProcessed)).append('"');
                json.append(",\"pipeline_version\":\"").append(escJ(pm.pipelineVersion)).append('"');
                json.append(",\"rsid_recovery_status\":\"").append(escJ(pm.rsidRecoveryStatus)).append('"');
                json.append(",\"rsid_recovery_rate\":\"").append(escJ(pm.rsidRecoveryRate)).append('"');
            }

            // rsID status from project.json (always available, even if config fails)
            boolean rsidPresent = pm != null && pm.rsidColumnPresent;
            boolean hasLoci = false;

            // Live-computed pipeline status
            String status;
            if (processing.contains(id)) {
                status = "processing";
            } else {
                try {
                    Config cfg = Config.loadFromProject(dir.getAbsolutePath());
                    ProjectMetadata.StaleReason reason =
                        ProjectMetadata.checkStaleness(dir.getAbsolutePath(), cfg);
                    status = reason == ProjectMetadata.StaleReason.NOT_STALE
                        ? "up_to_date" : "needs_reprocessing";
                    // Also check config for rsID column (covers projects that had rsIDs from the start)
                    if (cfg.colRsid != null && !cfg.colRsid.isEmpty()) rsidPresent = true;
                    // Check if loci file exists
                    if (cfg.lociFile != null && !cfg.lociFile.isEmpty() && new File(cfg.lociFile).exists()) {
                        hasLoci = true;
                    }
                } catch (Exception e) {
                    status = "error";
                }
            }
            json.append(",\"status\":\"").append(status).append('"');
            json.append(",\"rsid_column_present\":").append(rsidPresent);
            json.append(",\"has_loci\":").append(hasLoci);
            json.append('}');
        }
        json.append(']');
        respond(ex, 200, "application/json", json.toString().getBytes("UTF-8"));
    }

    // POST /api/projects/process — { "id": "..." } or { "all_stale": true }
    private void processProjects(HttpExchange ex) throws IOException {
        cors(ex); if (preflight(ex)) return;
        if (!"POST".equalsIgnoreCase(ex.getRequestMethod())) {
            respond(ex, 405, "application/json", "{\"error\":\"POST required\"}".getBytes()); return;
        }
        String body = new String(readAll(ex.getRequestBody()), "UTF-8");
        String id = extractStr(body, "id");
        String allStale = extractStr(body, "all_stale");

        List<String> toProcess = new ArrayList<>();
        if (id != null && !id.isEmpty()) {
            File dir = new File("projects/" + id);
            if (!dir.isDirectory() || !new File(dir, "config.properties").exists()) {
                respond(ex, 404, "application/json",
                    ("{\"error\":\"Project not found: " + escJ(id) + "\"}").getBytes()); return;
            }
            toProcess.add(id);
        } else if ("true".equals(allStale)) {
            for (File dir : Main.discoverProjects(new File("projects"))) {
                String pid = dir.getName();
                try {
                    Config cfg = Config.loadFromProject(dir.getAbsolutePath());
                    if (ProjectMetadata.checkStaleness(dir.getAbsolutePath(), cfg)
                            != ProjectMetadata.StaleReason.NOT_STALE)
                        toProcess.add(pid);
                } catch (Exception ignored) {}
            }
        } else {
            respond(ex, 400, "application/json",
                "{\"error\":\"Provide 'id' or 'all_stale':true\"}".getBytes()); return;
        }

        // Launch processing in background
        for (String pid : toProcess) {
            if (processing.contains(pid)) continue;
            processing.add(pid);
            ProgressTracker pt = new ProgressTracker();
            projectProgress.put(pid, pt);
            new Thread(() -> processProjectBackground(pid, pt), "process-" + pid).start();
        }

        respond(ex, 200, "application/json",
            ("{\"processing\":" + toProcess.size() + "}").getBytes());
    }

    private void processProjectBackground(String id, ProgressTracker pt) {
        String projectDir = new File("projects", id).getAbsolutePath();
        try {
            Config cfg = Config.loadFromProject(projectDir);
            new File(cfg.outputDir + "/data").mkdirs();
            new File(cfg.outputDir + "/plots").mkdirs();
            new File(cfg.outputDir + "/tmp").mkdirs();
            if (cfg.ldEnabled) {
                new File(cfg.plinkSubsetsDir()).mkdirs();
                new File(cfg.ldResultsDir()).mkdirs();
            }

            Main.PipelineResult result = Main.runPipeline(cfg, projectDir, pt);
            pt.done = true;

            // Write project.json
            Set<String> uniqueSnps = new HashSet<>();
            for (LocusOutput o : result.outputs)
                for (Snp s : o.gwasSnps)
                    uniqueSnps.add(s.chr + ":" + s.pos);

            ProjectMetadata pm = ProjectMetadata.load(projectDir);
            ProjectMetadata meta = new ProjectMetadata();
            meta.id   = id;
            meta.name = (pm != null && !pm.name.isEmpty()) ? pm.name : id;
            meta.description = (pm != null) ? pm.description : "";
            meta.lociCount = result.outputs.size();
            meta.totalSnps = uniqueSnps.size();
            String annotPath = new File(projectDir, "annotations.yaml").getAbsolutePath();
            meta.countAnnotationSources(annotPath);
            meta.coreInputFingerprint = ProjectMetadata.computeCoreInputFingerprint(cfg, projectDir);
            meta.annotationFingerprint = ProjectMetadata.computeAnnotationFingerprint(annotPath);
            meta.save(projectDir);

            setProjectState(id, result.config, result.gff, result.loci, result.outputs);
            System.out.printf("[Server] Project '%s' processed: %d loci%n", id, meta.lociCount);
        } catch (Exception e) {
            pt.phase = "Error: " + e.getMessage();
            pt.done = true;
            System.err.printf("[Server] Project '%s' processing failed: %s%n", id, e.getMessage());
        } finally {
            processing.remove(id);
        }
    }

    // ── Project router: /api/project/{id}/... ────────────────────────────

    // POST /api/project-delete — { "id": "..." }
    private void projectDelete(HttpExchange ex) throws IOException {
        cors(ex); if (preflight(ex)) return;
        if (!"POST".equalsIgnoreCase(ex.getRequestMethod())) {
            respond(ex, 405, "application/json", "{\"error\":\"POST required\"}".getBytes()); return;
        }
        String body = new String(readAll(ex.getRequestBody()), "UTF-8");
        String id = extractStr(body, "id");
        if (id == null || id.isEmpty()) {
            respond(ex, 400, "application/json", "{\"error\":\"missing id\"}".getBytes()); return;
        }
        File dir = new File("projects", id);
        if (!dir.isDirectory()) {
            respond(ex, 404, "application/json", "{\"error\":\"Project not found\"}".getBytes()); return;
        }
        deleteDirectory(dir);
        projectStates.remove(id);
        respond(ex, 200, "application/json", "{\"ok\":true}".getBytes());
        System.out.printf("[Server] Deleted project '%s'%n", id);
    }

    private void projectRouter(HttpExchange ex) throws IOException {
        cors(ex); if (preflight(ex)) return;
        String path = ex.getRequestURI().getPath();
        // path = /api/project/{id}/...
        String rest = path.substring("/api/project/".length());
        int slash = rest.indexOf('/');
        String projectId = slash >= 0 ? rest.substring(0, slash) : rest;
        String action = slash >= 0 ? rest.substring(slash + 1) : "";

        String projectDir = new File("projects", projectId).getAbsolutePath();
        String dataDir = projectDir + "/data";

        // DELETE /api/project/{id} — delete entire project
        if (action.isEmpty() && "DELETE".equalsIgnoreCase(ex.getRequestMethod())) {
            File dir = new File("projects", projectId);
            if (!dir.isDirectory()) {
                respond(ex, 404, "application/json", "{\"error\":\"Project not found\"}".getBytes());
                return;
            }
            deleteDirectory(dir);
            projectStates.remove(projectId);
            respond(ex, 200, "application/json", "{\"ok\":true}".getBytes());
            System.out.printf("[Server] Deleted project '%s'%n", projectId);
            return;
        }

        if (action.equals("manifest")) {
            serveJson(ex, dataDir + "/manifest.json");
        } else if (action.startsWith("locus/")) {
            String n = action.substring("locus/".length()).replaceAll("\\D", "");
            serveJson(ex, dataDir + "/locus_" + n + ".json");
        } else if (action.equals("config")) {
            if ("POST".equalsIgnoreCase(ex.getRequestMethod()))
                projectConfigPost(ex, projectId, projectDir);
            else
                projectConfigGet(ex, projectDir);
        } else if (action.equals("annotations-yaml")) {
            File yamlFile = new File(projectDir, "annotations.yaml");
            if (yamlFile.exists()) {
                byte[] bytes = Files.readAllBytes(yamlFile.toPath());
                respond(ex, 200, "text/yaml", bytes);
            } else {
                respond(ex, 404, "text/plain", "annotations.yaml not found".getBytes());
            }
        } else if (action.equals("annotation-file")) {
            projectAnnotationFile(ex, projectDir);
        } else if (action.equals("save-annotation-config")) {
            projectSaveAnnotationConfig(ex, projectDir);
        } else if (action.equals("update-locus")) {
            projectUpdateLocus(ex, projectId, projectDir);
        } else if (action.equals("create-locus")) {
            projectCreateLocus(ex, projectId, projectDir);
        } else if (action.equals("split-locus")) {
            projectSplitLocus(ex, projectId, projectDir);
        } else if (action.equals("validate-split")) {
            projectValidateSplit(ex, projectId, projectDir);
        } else if (action.equals("merge-loci")) {
            projectMergeLoci(ex, projectId, projectDir);
        } else if (action.equals("undo-mutation")) {
            projectUndoMutation(ex, projectId, projectDir);
        } else if (action.equals("progress")) {
            projectProgressEndpoint(ex, projectId);
        } else {
            respond(ex, 404, "application/json",
                ("{\"error\":\"Unknown action: " + escJ(action) + "\"}").getBytes());
        }
    }

    // GET /api/project/{id}/config
    private void projectConfigGet(HttpExchange ex, String projectDir) throws IOException {
        File configFile = new File(projectDir, "config.properties");
        if (!configFile.exists()) {
            respond(ex, 404, "application/json", "{\"error\":\"No config.properties\"}".getBytes());
            return;
        }
        try {
            Config cfg = Config.loadFromProject(projectDir);
            respond(ex, 200, "application/json", cfg.toJson().getBytes("UTF-8"));
        } catch (Exception e) {
            respond(ex, 500, "application/json",
                ("{\"error\":\"" + escJ(e.getMessage()) + "\"}").getBytes());
        }
    }

    // POST /api/project/{id}/config
    private void projectConfigPost(HttpExchange ex, String projectId,
                                   String projectDir) throws IOException {
        String body = new String(readAll(ex.getRequestBody()), "UTF-8");
        try {
            new File(projectDir).mkdirs();
            Config cfg = new Config();
            cfg.applyJson(body);
            cfg.writeProperties(new File(projectDir, "config.properties").getAbsolutePath());
            respond(ex, 200, "application/json", "{\"ok\":true}".getBytes());
        } catch (Exception e) {
            respond(ex, 500, "application/json",
                ("{\"error\":\"" + escJ(e.getMessage()) + "\"}").getBytes());
        }
    }

    // GET /api/peek-file-header?path=...
    private void peekFileHeader(HttpExchange ex) throws IOException {
        cors(ex); if (preflight(ex)) return;
        String query = ex.getRequestURI().getQuery();
        String filePath = null;
        if (query != null) {
            for (String param : query.split("&")) {
                String[] kv = param.split("=", 2);
                if (kv.length == 2 && kv[0].equals("path"))
                    filePath = URLDecoder.decode(kv[1], "UTF-8");
            }
        }
        if (filePath == null) {
            respond(ex, 400, "application/json", "{\"error\":\"Missing path parameter\"}".getBytes());
            return;
        }
        String normalized = filePath.replace("\\", "/");
        if (normalized.contains("/../") || normalized.startsWith("../")) {
            respond(ex, 400, "application/json", "{\"error\":\"Path traversal not allowed\"}".getBytes());
            return;
        }
        File f = new File(filePath);
        if (!f.exists()) {
            respond(ex, 404, "application/json",
                ("{\"error\":\"File not found: " + escJ(filePath) + "\"}").getBytes());
            return;
        }
        try (BufferedReader br = new BufferedReader(new FileReader(f))) {
            String header = br.readLine();
            if (header == null) {
                respond(ex, 200, "application/json", "{\"columns\":[]}".getBytes());
                return;
            }
            String[] cols = header.split("\t");
            StringBuilder json = new StringBuilder("{\"columns\":[");
            for (int i = 0; i < cols.length; i++) {
                if (i > 0) json.append(',');
                json.append('"').append(cols[i].replace("\"", "\\\"")).append('"');
            }
            json.append("]}");
            respond(ex, 200, "application/json", json.toString().getBytes("UTF-8"));
        }
    }

    // GET /api/project/{id}/annotation-file?path=...
    private void projectAnnotationFile(HttpExchange ex, String projectDir) throws IOException {
        String query = ex.getRequestURI().getQuery();
        String filePath = null;
        if (query != null) {
            for (String param : query.split("&")) {
                String[] kv = param.split("=", 2);
                if (kv.length == 2 && kv[0].equals("path"))
                    filePath = URLDecoder.decode(kv[1], "UTF-8");
            }
        }
        if (filePath == null) {
            respond(ex, 400, "application/json", "{\"error\":\"Missing path parameter\"}".getBytes());
            return;
        }
        String normalized = filePath.replace("\\", "/");
        if (normalized.contains("/../") || normalized.startsWith("../")) {
            respond(ex, 400, "application/json", "{\"error\":\"Path traversal not allowed\"}".getBytes());
            return;
        }
        File f = new File(filePath);
        if (!f.exists()) f = new File(projectDir + "/" + filePath);
        if (!f.exists()) {
            File projectRoot = new File(projectDir).getParentFile();
            if (projectRoot != null) f = new File(projectRoot.getParent(), filePath);
        }
        if (!f.exists()) {
            respond(ex, 404, "application/json",
                ("{\"error\":\"File not found: " + filePath + "\"}").getBytes());
            return;
        }
        List<String> lines = Files.readAllLines(f.toPath());
        if (lines.isEmpty()) {
            respond(ex, 200, "application/json", "{\"columns\":[],\"rows\":[],\"row_count\":0}".getBytes());
            return;
        }
        String[] header = lines.get(0).split("\t");
        StringBuilder sb = new StringBuilder();
        sb.append("{\"columns\":[");
        for (int i = 0; i < header.length; i++) {
            if (i > 0) sb.append(',');
            sb.append('"').append(header[i].replace("\"","\\\"")).append('"');
        }
        sb.append("],\"rows\":[");
        for (int r = 1; r < lines.size(); r++) {
            if (r > 1) sb.append(',');
            String[] vals = lines.get(r).split("\t", -1);
            sb.append('{');
            for (int c = 0; c < header.length && c < vals.length; c++) {
                if (c > 0) sb.append(',');
                sb.append('"').append(header[c].replace("\"","\\\"")).append("\":\"");
                sb.append(vals[c].replace("\\","\\\\").replace("\"","\\\"")).append('"');
            }
            sb.append('}');
        }
        sb.append("],\"row_count\":").append(lines.size() - 1).append('}');
        respond(ex, 200, "application/json", sb.toString().getBytes("UTF-8"));
    }

    // POST /api/project/{id}/save-annotation-config
    private void projectSaveAnnotationConfig(HttpExchange ex, String projectDir) throws IOException {
        if (!"POST".equalsIgnoreCase(ex.getRequestMethod())) {
            respond(ex, 405, "application/json", "{\"error\":\"POST required\"}".getBytes()); return;
        }
        String body = new String(ex.getRequestBody().readAllBytes(), "UTF-8");
        String path = extractJsonString(body, "path");
        String content = extractJsonString(body, "content");
        if (path == null || content == null || path.contains("..")) {
            respond(ex, 400, "application/json", "{\"error\":\"Invalid request\"}".getBytes()); return;
        }
        File target = new File(projectDir + "/" + path);
        Files.writeString(target.toPath(), content);
        respond(ex, 200, "application/json", "{\"ok\":true}".getBytes());
    }

    // POST /api/project/{id}/update-locus
    private void projectUpdateLocus(HttpExchange ex, String projectId,
                                    String projectDir) throws IOException {
        if (!"POST".equalsIgnoreCase(ex.getRequestMethod())) {
            respond(ex, 405, "text/plain", "Method Not Allowed".getBytes()); return;
        }
        ProjectState ps = projectStates.get(projectId);
        if (ps == null || ps.config == null || ps.gff == null || ps.loci == null) {
            respond(ex, 503, "application/json",
                "{\"error\":\"Pipeline state not available. Process the project first.\"}".getBytes());
            return;
        }
        String body = new String(readAll(ex.getRequestBody()), "UTF-8");
        String sIdx   = extractStr(body, "locus_index");
        String sStart = extractStr(body, "new_start");
        String sEnd   = extractStr(body, "new_end");
        if (sIdx == null || sStart == null || sEnd == null) {
            respond(ex, 400, "application/json",
                "{\"error\":\"missing locus_index, new_start, or new_end\"}".getBytes()); return;
        }
        LocusUpdater.UpdateResult ur = LocusUpdater.update(
            Integer.parseInt(sIdx), Long.parseLong(sStart), Long.parseLong(sEnd),
            ps.config, ps.gff, ps.loci);
        if (ur.ok) {
            byte[] json = Files.readAllBytes(new File(ur.jsonPath).toPath());
            respond(ex, 200, "application/json", json);
        } else {
            respond(ex, 500, "application/json",
                ("{\"error\":\"" + escJ(ur.error) + "\"}").getBytes());
        }
    }

    // POST /api/project/{id}/split-locus
    private void projectSplitLocus(HttpExchange ex, String projectId,
                                   String projectDir) throws IOException {
        if (!"POST".equalsIgnoreCase(ex.getRequestMethod())) {
            respond(ex, 405, "text/plain", "Method Not Allowed".getBytes()); return;
        }
        ProjectState ps = projectStates.get(projectId);
        if (ps == null || ps.config == null || ps.gff == null || ps.loci == null || ps.outputs == null) {
            respond(ex, 503, "application/json",
                "{\"error\":\"Pipeline state not available\"}".getBytes()); return;
        }
        String body = new String(readAll(ex.getRequestBody()), "UTF-8");
        String sOrigIdx = extractStr(body, "locus_index");
        if (sOrigIdx == null) {
            respond(ex, 400, "application/json", "{\"error\":\"missing locus_index\"}".getBytes()); return;
        }
        List<LocusUpdater.SplitRegion> regions = parseSplitRegions(body);
        if (regions.size() < 2) {
            respond(ex, 400, "application/json", "{\"error\":\"need at least 2 regions\"}".getBytes()); return;
        }
        LocusUpdater.UpdateResult ur = LocusUpdater.split(
            Integer.parseInt(sOrigIdx), regions, ps.config, ps.gff, ps.loci, ps.outputs);
        if (ur.ok) {
            String dataDir = projectDir + "/data";
            byte[] manifestBytes = Files.readAllBytes(new File(dataDir + "/manifest.json").toPath());
            respond(ex, 200, "application/json", manifestBytes);
        } else {
            respond(ex, 500, "application/json",
                ("{\"error\":\"" + escJ(ur.error) + "\"}").getBytes());
        }
    }

    // POST /api/project/{id}/create-locus
    private void projectCreateLocus(HttpExchange ex, String projectId,
                                    String projectDir) throws IOException {
        if (!"POST".equalsIgnoreCase(ex.getRequestMethod())) {
            respond(ex, 405, "text/plain", "Method Not Allowed".getBytes()); return;
        }
        ProjectState ps = projectStates.get(projectId);
        if (ps == null || ps.config == null || ps.gff == null || ps.loci == null || ps.outputs == null) {
            respond(ex, 503, "application/json",
                "{\"error\":\"Pipeline state not available. Process the project first.\"}".getBytes());
            return;
        }
        String body  = new String(readAll(ex.getRequestBody()), "UTF-8");
        String chr   = extractStr(body, "chr");
        String sStart = extractStr(body, "start");
        String sEnd   = extractStr(body, "end");
        String name   = extractStr(body, "locus_name");
        if (chr == null || sStart == null || sEnd == null) {
            respond(ex, 400, "application/json",
                "{\"error\":\"missing chr, start, or end\"}".getBytes()); return;
        }
        LocusUpdater.UpdateResult ur = LocusUpdater.create(
            chr, Long.parseLong(sStart), Long.parseLong(sEnd), name,
            ps.config, ps.gff, ps.loci, ps.outputs);
        if (ur.ok) {
            byte[] manifestBytes = Files.readAllBytes(
                new File(projectDir + "/data/manifest.json").toPath());
            respond(ex, 200, "application/json", manifestBytes);
        } else {
            respond(ex, 500, "application/json",
                ("{\"error\":\"" + escJ(ur.error) + "\"}").getBytes());
        }
    }

    // POST /api/project/{id}/validate-split
    private void projectValidateSplit(HttpExchange ex, String projectId,
                                     String projectDir) throws IOException {
        if (!"POST".equalsIgnoreCase(ex.getRequestMethod())) {
            respond(ex, 405, "text/plain", "Method Not Allowed".getBytes()); return;
        }
        ProjectState ps = projectStates.get(projectId);
        if (ps == null || ps.config == null || ps.loci == null) {
            respond(ex, 503, "application/json",
                "{\"error\":\"Pipeline state not available\"}".getBytes()); return;
        }
        String body = new String(readAll(ex.getRequestBody()), "UTF-8");
        String sOrigIdx = extractStr(body, "locus_index");
        if (sOrigIdx == null) {
            respond(ex, 400, "application/json", "{\"error\":\"missing locus_index\"}".getBytes()); return;
        }
        List<LocusUpdater.SplitRegion> regions = parseSplitRegions(body);
        LocusUpdater.ValidationResult vr = LocusUpdater.validateSplit(
            Integer.parseInt(sOrigIdx), regions, ps.config, ps.loci);
        StringBuilder json = new StringBuilder();
        json.append("{\"valid\":").append(vr.valid);
        json.append(",\"split_ld_threshold\":").append(ps.config.splitLdThreshold);
        json.append(",\"split_min_distance_bp\":").append(ps.config.splitMinDistBp);
        json.append(",\"warnings\":[");
        for (int i = 0; i < vr.warnings.size(); i++) {
            if (i > 0) json.append(',');
            json.append('"').append(escJ(vr.warnings.get(i))).append('"');
        }
        json.append("],\"ld_violations\":[");
        for (int i = 0; i < vr.ldViolations.size(); i++) {
            if (i > 0) json.append(',');
            LocusUpdater.LdViolation v = vr.ldViolations.get(i);
            json.append("{\"snp_a\":\"").append(escJ(v.snpA))
                .append("\",\"snp_b\":\"").append(escJ(v.snpB))
                .append("\",\"pos_a\":").append(v.posA)
                .append(",\"pos_b\":").append(v.posB)
                .append(",\"region_a\":").append(v.regionA)
                .append(",\"region_b\":").append(v.regionB)
                .append(",\"r2\":").append(String.format("%.4f", v.r2))
                .append('}');
        }
        json.append("]}");
        respond(ex, 200, "application/json", json.toString().getBytes());
    }

    // POST /api/project/{id}/merge-loci
    private void projectMergeLoci(HttpExchange ex, String projectId,
                                  String projectDir) throws IOException {
        if (!"POST".equalsIgnoreCase(ex.getRequestMethod())) {
            respond(ex, 405, "text/plain", "Method Not Allowed".getBytes()); return;
        }
        ProjectState ps = projectStates.get(projectId);
        if (ps == null || ps.mutationService == null) {
            respond(ex, 503, "application/json",
                "{\"error\":\"Pipeline state not available\"}".getBytes()); return;
        }
        String body = new String(readAll(ex.getRequestBody()), "UTF-8");

        // Parse locus_indices array
        List<Integer> indices = new ArrayList<>();
        int arrStart = body.indexOf("\"locus_indices\"");
        if (arrStart >= 0) {
            int s = body.indexOf('[', arrStart);
            int e = body.indexOf(']', s);
            if (s >= 0 && e > s) {
                for (String part : body.substring(s+1, e).split(",")) {
                    part = part.trim();
                    if (!part.isEmpty()) indices.add(Integer.parseInt(part));
                }
            }
        }
        if (indices.size() < 2) {
            respond(ex, 400, "application/json",
                "{\"error\":\"need at least 2 locus_indices\"}".getBytes()); return;
        }
        String mergedName = extractStr(body, "merged_name");

        LociMutationService.MutationResult mr = ps.mutationService.merge(indices, mergedName);
        if (mr.ok) {
            respond(ex, 200, "application/json", mr.manifestJson.getBytes("UTF-8"));
        } else {
            respond(ex, 500, "application/json",
                ("{\"error\":\"" + escJ(mr.error) + "\"}").getBytes());
        }
    }

    // POST /api/project/{id}/undo-mutation
    private void projectUndoMutation(HttpExchange ex, String projectId,
                                     String projectDir) throws IOException {
        if (!"POST".equalsIgnoreCase(ex.getRequestMethod())) {
            respond(ex, 405, "text/plain", "Method Not Allowed".getBytes()); return;
        }
        ProjectState ps = projectStates.get(projectId);
        if (ps == null || ps.mutationService == null) {
            respond(ex, 503, "application/json",
                "{\"error\":\"Pipeline state not available\"}".getBytes()); return;
        }

        LociMutationService.MutationResult mr = ps.mutationService.undo();
        if (mr.ok) {
            String json = "{\"ok\":true,\"manifest\":" + mr.manifestJson + "}";
            respond(ex, 200, "application/json", json.getBytes("UTF-8"));
        } else {
            respond(ex, 500, "application/json",
                ("{\"error\":\"" + escJ(mr.error) + "\"}").getBytes());
        }
    }

    // GET /api/project/{id}/progress
    private void projectProgressEndpoint(HttpExchange ex, String projectId) throws IOException {
        ProgressTracker pt = projectProgress.get(projectId);
        if (pt == null) pt = new ProgressTracker();
        String json = String.format(
            "{\"phase\":\"%s\",\"locus\":%d,\"total\":%d,\"pct\":%d,\"done\":%s}",
            escJ(pt.phase), pt.locusIndex, pt.totalLoci, pt.pct(),
            pt.done ? "true" : "false");
        respond(ex, 200, "application/json", json.getBytes());
    }

    // ══════════════════════════════════════════════════════════════════════
    //  rsID RECOVERY ENDPOINTS
    // ══════════════════════════════════════════════════════════════════════

    private final Map<String, RsidProgress> rsidProgressMap = new ConcurrentHashMap<>();
    private final Map<String, LociProgress> lociProgressMap = new ConcurrentHashMap<>();
    private final Map<String, ExcelExporter.ExportProgress> exportProgressMap = new ConcurrentHashMap<>();

    // GET/POST /api/global-config
    private void globalConfigEndpoint(HttpExchange ex) throws IOException {
        cors(ex); if (preflight(ex)) return;
        if ("POST".equalsIgnoreCase(ex.getRequestMethod())) {
            String body = new String(readAll(ex.getRequestBody()), "UTF-8");
            try {
                // Parse and save - for now, accept the full JSON
                GlobalConfig gc = GlobalConfig.load();
                // Apply updates from body (simplified: just save the body as global config)
                new File("config").mkdirs();
                // Normalize Windows backslashes in paths to forward slashes
                body = body.replace("\\\\", "/").replace("\\", "/");
                Files.writeString(Path.of("config/global.json"), body);
                respond(ex, 200, "application/json", "{\"ok\":true}".getBytes());
            } catch (Exception e) {
                respond(ex, 500, "application/json",
                    ("{\"error\":\"" + escJ(e.getMessage()) + "\"}").getBytes());
            }
        } else {
            GlobalConfig gc = GlobalConfig.load();
            gc.validateAll();
            respond(ex, 200, "application/json", gc.toJson().getBytes("UTF-8"));
        }
    }

    // POST /api/rsid-recover — { "project_id": "...", "snp_database_id": "..." }
    private void rsidRecover(HttpExchange ex) throws IOException {
        cors(ex); if (preflight(ex)) return;
        if (!"POST".equalsIgnoreCase(ex.getRequestMethod())) {
            respond(ex, 405, "application/json", "{\"error\":\"POST required\"}".getBytes()); return;
        }
        String body = new String(readAll(ex.getRequestBody()), "UTF-8");
        String projectId = extractStr(body, "project_id");
        String snpDbId   = extractStr(body, "snp_database_id");

        if (projectId == null || projectId.isEmpty()) {
            respond(ex, 400, "application/json", "{\"error\":\"missing project_id\"}".getBytes()); return;
        }

        File projDir = new File("projects/" + projectId);
        if (!projDir.isDirectory()) {
            respond(ex, 404, "application/json", "{\"error\":\"Project not found\"}".getBytes()); return;
        }

        // Load project config
        Config projConfig;
        try { projConfig = Config.loadFromProject(projDir.getAbsolutePath()); }
        catch (Exception e) {
            respond(ex, 500, "application/json",
                ("{\"error\":\"Config load failed: " + escJ(e.getMessage()) + "\"}").getBytes()); return;
        }

        // Resolve SNP database
        GlobalConfig gc = GlobalConfig.load();
        gc.validateAll();
        GlobalConfig.SnpDatabase snpDb = null;
        if (snpDbId != null && !snpDbId.isEmpty()) {
            snpDb = gc.findDatabase(snpDbId);
        } else if (!gc.snpDatabases.isEmpty()) {
            snpDb = gc.snpDatabases.get(0);
        }
        if (snpDb == null) {
            respond(ex, 400, "application/json",
                "{\"error\":\"No SNP database configured. Add one in Resources & Settings.\"}".getBytes());
            return;
        }

        // Block if another heavy process is running for this project
        if (processing.contains(projectId)) {
            respond(ex, 409, "application/json",
                "{\"error\":\"Project is currently being processed. Wait for it to finish.\"}".getBytes());
            return;
        }

        // Parse API completion config from request
        RsidApiCompleter.Config apiCfg = null;
        String apiEnabled = extractStr(body, "api_enabled");
        if ("true".equals(apiEnabled)) {
            apiCfg = new RsidApiCompleter.Config();
            apiCfg.enabled = true;
            String maxU = extractStr(body, "api_max_unmatched");
            if (maxU != null && !maxU.isEmpty()) apiCfg.maxUnmatched = Integer.parseInt(maxU);
            String batch = extractStr(body, "api_batch_size");
            if (batch != null && !batch.isEmpty()) apiCfg.batchSize = Math.min(200, Integer.parseInt(batch));
            String ncbiKey = extractStr(body, "ncbi_api_key");
            if (ncbiKey != null) apiCfg.ncbiApiKey = ncbiKey;
        }

        // Start recovery in background
        RsidProgress rp = new RsidProgress();
        rsidProgressMap.put(projectId, rp);

        String snpFolder = snpDb.folder;
        String build = snpDb.build;
        RsidApiCompleter.Config finalApiCfg = apiCfg;
        new Thread(() -> RsidPipeline.run(
            projDir.getAbsolutePath(), projConfig.gwasFile,
            projConfig.lociFile, snpFolder, build,
            projConfig.colChr, projConfig.colPos, projConfig.colEa, projConfig.colNea,
            rp, finalApiCfg
        ), "rsid-" + projectId).start();

        respond(ex, 200, "application/json", "{\"started\":true}".getBytes());
    }

    // GET /api/rsid-progress?project=...
    private void rsidProgressEndpoint(HttpExchange ex) throws IOException {
        cors(ex); if (preflight(ex)) return;
        String query = ex.getRequestURI().getQuery();
        String projectId = null;
        if (query != null) {
            for (String param : query.split("&")) {
                String[] kv = param.split("=", 2);
                if (kv.length == 2 && kv[0].equals("project"))
                    projectId = URLDecoder.decode(kv[1], "UTF-8");
            }
        }
        RsidProgress rp = projectId != null ? rsidProgressMap.get(projectId) : null;
        if (rp == null) {
            respond(ex, 200, "application/json",
                "{\"pct\":0,\"current_step\":\"Not started\",\"step_index\":0,\"total_steps\":5,\"done\":false,\"error\":null,\"matched\":0,\"forward\":0,\"reverse\":0,\"unmatched\":0,\"total_snps\":0,\"current_locus\":0,\"total_loci\":0,\"recovery_rate\":0,\"output_file\":\"\"}".getBytes());
            return;
        }
        respond(ex, 200, "application/json", rp.toJson().getBytes("UTF-8"));
    }

    // GET /api/rsid-detect?gwas_file=...
    private void rsidDetect(HttpExchange ex) throws IOException {
        cors(ex); if (preflight(ex)) return;
        String query = ex.getRequestURI().getQuery();
        String gwasFile = null;
        if (query != null) {
            for (String param : query.split("&")) {
                String[] kv = param.split("=", 2);
                if (kv.length == 2 && kv[0].equals("gwas_file"))
                    gwasFile = URLDecoder.decode(kv[1], "UTF-8");
            }
        }
        if (gwasFile == null || gwasFile.isEmpty()) {
            respond(ex, 400, "application/json", "{\"error\":\"missing gwas_file\"}".getBytes()); return;
        }
        try {
            RsidDetector.DetectionResult dr = RsidDetector.detectFromFile(gwasFile);
            String json = String.format("{\"present\":%s,\"column\":%s,\"confidence\":%.3f}",
                dr.present,
                dr.column != null ? "\"" + escJ(dr.column) + "\"" : "null",
                dr.confidence);
            respond(ex, 200, "application/json", json.getBytes("UTF-8"));
        } catch (Exception e) {
            respond(ex, 500, "application/json",
                ("{\"error\":\"" + escJ(e.getMessage()) + "\"}").getBytes());
        }
    }

    // POST /api/loci-identify — { "project_id": "...", params... }
    private void lociIdentify(HttpExchange ex) throws IOException {
        cors(ex); if (preflight(ex)) return;
        if (!"POST".equalsIgnoreCase(ex.getRequestMethod())) {
            respond(ex, 405, "application/json", "{\"error\":\"POST required\"}".getBytes()); return;
        }
        String body = new String(readAll(ex.getRequestBody()), "UTF-8");
        String projectId = extractStr(body, "project_id");
        if (projectId == null || projectId.isEmpty()) {
            respond(ex, 400, "application/json", "{\"error\":\"missing project_id\"}".getBytes()); return;
        }
        File projDir = new File("projects/" + projectId);
        if (!projDir.isDirectory()) {
            respond(ex, 404, "application/json", "{\"error\":\"Project not found\"}".getBytes()); return;
        }

        Config projConfig;
        try { projConfig = Config.loadFromProject(projDir.getAbsolutePath()); }
        catch (Exception e) {
            respond(ex, 500, "application/json",
                ("{\"error\":\"" + escJ(e.getMessage()) + "\"}").getBytes()); return;
        }

        // Parse optional parameters
        LociIdentifier.Params params = new LociIdentifier.Params();
        String v;
        v = extractStr(body, "pre_filter_p");
        if (v != null && !v.isEmpty()) params.preFilterP = Double.parseDouble(v);
        v = extractStr(body, "lead_p_threshold");
        if (v != null && !v.isEmpty()) params.leadPThreshold = Double.parseDouble(v);
        v = extractStr(body, "merge_distance_kb");
        if (v != null && !v.isEmpty()) params.mergeDistanceBp = Integer.parseInt(v) * 1000;
        v = extractStr(body, "min_snps_per_locus");
        if (v != null && !v.isEmpty()) params.minSnpsPerLocus = Integer.parseInt(v);

        // Reference panel path for chr:pos matching and PLINK clumping
        params.refPanelPath = projConfig.refPanelPath;

        // Column mapping
        String colChr = projConfig.colChr;
        String colPos = projConfig.colPos;
        String colP   = projConfig.colPvalue;
        String colId  = projConfig.colRsid.isEmpty() ? projConfig.colVarid : projConfig.colRsid;

        LociProgress lp = new LociProgress();
        lociProgressMap.put(projectId, lp);
        String gwasFile = projConfig.gwasFile;
        String projectDir = projDir.getAbsolutePath();

        new Thread(() -> {
            try {
                List<LociIdentifier.IdentifiedLocus> loci = LociIdentifier.identify(
                    gwasFile, colChr, colPos, colP, colId, params, lp);

                String lociPath = projectDir + "/loci.txt";
                String detailPath = projectDir + "/loci_detail.tsv";
                LociIdentifier.writeLociFiles(loci, lociPath, detailPath);

                // Update project config to point to loci file
                File configFile = new File(projectDir, "config.properties");
                List<String> lines = Files.readAllLines(configFile.toPath());
                boolean wroteLoci = false;
                try (PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(configFile)))) {
                    for (String line : lines) {
                        if (line.trim().startsWith("loci.file=")) {
                            pw.println("loci.file=" + lociPath);
                            wroteLoci = true;
                        } else {
                            pw.println(line);
                        }
                    }
                    if (!wroteLoci) pw.println("loci.file=" + lociPath);
                }

                // Clear fingerprint cache so project is stale for reprocessing
                Files.deleteIfExists(new File(projectDir, ".fingerprint_cache").toPath());

                lp.done = true;
                lp.currentStep = "Complete";
                System.out.printf("[LociServer] Identified %d loci for project '%s'%n",
                    loci.size(), projectId);
            } catch (Exception e) {
                lp.error = e.getMessage();
                lp.done = true;
                lp.currentStep = "Error";
                System.err.printf("[LociServer] Failed for '%s': %s%n", projectId, e.getMessage());
                e.printStackTrace(System.err);
            }
        }, "loci-" + projectId).start();

        respond(ex, 200, "application/json", "{\"started\":true}".getBytes());
    }

    // GET /api/loci-progress?project=...
    private void lociProgressEndpoint(HttpExchange ex) throws IOException {
        cors(ex); if (preflight(ex)) return;
        String query = ex.getRequestURI().getQuery();
        String projectId = null;
        if (query != null) {
            for (String param : query.split("&")) {
                String[] kv = param.split("=", 2);
                if (kv.length == 2 && kv[0].equals("project"))
                    projectId = URLDecoder.decode(kv[1], "UTF-8");
            }
        }
        LociProgress lp = projectId != null ? lociProgressMap.get(projectId) : null;
        if (lp == null) {
            respond(ex, 200, "application/json",
                "{\"pct\":0,\"current_step\":\"Not started\",\"step_index\":0,\"total_steps\":4,\"done\":false,\"error\":null,\"candidate_snps\":0,\"seed_snps\":0,\"loci_found\":0,\"total_chromosomes\":0,\"current_chromosome\":0}".getBytes());
            return;
        }
        respond(ex, 200, "application/json", lp.toJson().getBytes("UTF-8"));
    }

    // POST /api/export-excel — { project_id, scope: "dataset"|"locus", locus_index? }
    private void exportExcel(HttpExchange ex) throws IOException {
        cors(ex); if (preflight(ex)) return;
        if (!"POST".equalsIgnoreCase(ex.getRequestMethod())) {
            respond(ex, 405, "application/json", "{\"error\":\"POST required\"}".getBytes()); return;
        }
        String body = new String(readAll(ex.getRequestBody()), "UTF-8");
        String projectId = extractStr(body, "project_id");
        String scope = extractStr(body, "scope");
        if (scope == null) scope = "dataset";

        if (projectId == null || projectId.isEmpty()) {
            respond(ex, 400, "application/json", "{\"error\":\"missing project_id\"}".getBytes()); return;
        }
        File projDir = new File("projects/" + projectId);
        if (!projDir.isDirectory()) {
            respond(ex, 404, "application/json", "{\"error\":\"Project not found\"}".getBytes()); return;
        }

        String projectDir = projDir.getAbsolutePath();
        String outputFolder = extractStr(body, "output_folder");
        File exportsDir;
        if (outputFolder != null && !outputFolder.isEmpty()) {
            exportsDir = new File(outputFolder);
        } else {
            exportsDir = new File(projectDir, "exports");
        }
        exportsDir.mkdirs();
        String date = java.time.LocalDate.now().toString();
        String outputPath;

        ExcelExporter.ExportProgress ep = new ExcelExporter.ExportProgress();
        exportProgressMap.put(projectId, ep);

        if ("locus".equals(scope)) {
            String sIdx = extractStr(body, "locus_index");
            int locusIdx = sIdx != null ? Integer.parseInt(sIdx) : 1;
            outputPath = new File(exportsDir, projectId + "_locus" + locusIdx + "_" + date + ".xlsx").getAbsolutePath();
            String finalPath = outputPath;
            new Thread(() -> ExcelExporter.exportLocus(projectDir, locusIdx, finalPath, ep),
                "excel-" + projectId).start();
        } else {
            outputPath = new File(exportsDir, projectId + "_" + date + ".xlsx").getAbsolutePath();
            String finalPath = outputPath;
            new Thread(() -> ExcelExporter.exportDataset(projectDir, finalPath, ep),
                "excel-" + projectId).start();
        }

        respond(ex, 200, "application/json",
            ("{\"started\":true,\"output\":\"" + escJ(outputPath) + "\"}").getBytes());
    }

    // GET /api/export-progress?project=...
    private void exportProgressEndpoint(HttpExchange ex) throws IOException {
        cors(ex); if (preflight(ex)) return;
        String query = ex.getRequestURI().getQuery();
        String projectId = null;
        if (query != null) for (String p : query.split("&")) {
            String[] kv = p.split("=", 2);
            if (kv.length == 2 && kv[0].equals("project"))
                projectId = URLDecoder.decode(kv[1], "UTF-8");
        }
        ExcelExporter.ExportProgress ep = projectId != null ? exportProgressMap.get(projectId) : null;
        if (ep == null) {
            respond(ex, 200, "application/json",
                "{\"pct\":0,\"step_index\":0,\"total_steps\":6,\"current_step\":\"Not started\",\"done\":false,\"error\":null,\"current_locus\":0,\"total_loci\":0}".getBytes());
            return;
        }
        respond(ex, 200, "application/json", ep.toJson().getBytes("UTF-8"));
    }

    // ══════════════════════════════════════════════════════════════════════
    //  LEGACY ENDPOINTS (unchanged behavior)
    // ══════════════════════════════════════════════════════════════════════

    private void manifest(HttpExchange ex) throws IOException {
        cors(ex); if (preflight(ex)) return;
        serveJson(ex, outputDir + "/data/manifest.json");
    }

    private void locus(HttpExchange ex) throws IOException {
        cors(ex); if (preflight(ex)) return;
        String n = ex.getRequestURI().getPath().replaceAll("\\D", "");
        serveJson(ex, outputDir + "/data/locus_" + n + ".json");
    }

    private void pickFolder(HttpExchange ex) throws IOException {
        cors(ex); if (preflight(ex)) return;
        String chosen = folderPicker(lastFolder);
        if (chosen != null) lastFolder = chosen;
        String json = chosen != null
            ? "{\"path\":\"" + escJ(chosen) + "\"}"
            : "{\"path\":null}";
        respond(ex, 200, "application/json", json.getBytes());
    }

    private void savePdf(HttpExchange ex) throws IOException {
        cors(ex); if (preflight(ex)) return;
        if (!"POST".equalsIgnoreCase(ex.getRequestMethod())) {
            respond(ex, 405, "text/plain", "Method Not Allowed".getBytes()); return;
        }
        String json     = new String(readAll(ex.getRequestBody()), "UTF-8");
        String filename = extractStr(json, "filename");
        String folder   = extractStr(json, "folder");
        String data     = extractStr(json, "data");
        if (filename == null || folder == null || data == null) {
            respond(ex, 400, "application/json", "{\"error\":\"missing fields\"}".getBytes()); return;
        }
        int comma = data.indexOf(',');
        byte[] pdf = Base64.getDecoder().decode(
            (comma >= 0 ? data.substring(comma + 1) : data).replaceAll("\\s", ""));
        new File(folder).mkdirs();
        Files.write(Paths.get(folder, filename), pdf);
        respond(ex, 200, "application/json", "{\"ok\":true}".getBytes());
    }

    private void progressEndpoint(HttpExchange ex) throws IOException {
        cors(ex); if (preflight(ex)) return;
        String json = String.format(
            "{\"phase\":\"%s\",\"locus\":%d,\"total\":%d,\"pct\":%d,\"done\":%s}",
            escJ(progress.phase), progress.locusIndex,
            progress.totalLoci, progress.pct(),
            progress.done ? "true" : "false");
        respond(ex, 200, "application/json", json.getBytes());
    }

    private void annotationFile(HttpExchange ex) throws IOException {
        cors(ex); if (preflight(ex)) return;
        projectAnnotationFile(ex, outputDir);
    }

    private void saveAnnotationConfig(HttpExchange ex) throws IOException {
        cors(ex); if (preflight(ex)) return;
        projectSaveAnnotationConfig(ex, outputDir);
    }

    private void updateLocus(HttpExchange ex) throws IOException {
        cors(ex); if (preflight(ex)) return;
        if (!"POST".equalsIgnoreCase(ex.getRequestMethod())) {
            respond(ex, 405, "text/plain", "Method Not Allowed".getBytes()); return;
        }
        if (config == null || gff == null || loci == null) {
            respond(ex, 503, "application/json",
                "{\"error\":\"Pipeline state not available\"}".getBytes()); return;
        }
        String body = new String(readAll(ex.getRequestBody()), "UTF-8");
        String sIdx   = extractStr(body, "locus_index");
        String sStart = extractStr(body, "new_start");
        String sEnd   = extractStr(body, "new_end");
        if (sIdx == null || sStart == null || sEnd == null) {
            respond(ex, 400, "application/json",
                "{\"error\":\"missing locus_index, new_start, or new_end\"}".getBytes()); return;
        }
        LocusUpdater.UpdateResult ur = LocusUpdater.update(
            Integer.parseInt(sIdx), Long.parseLong(sStart), Long.parseLong(sEnd),
            config, gff, loci);
        if (ur.ok) {
            respond(ex, 200, "application/json", Files.readAllBytes(new File(ur.jsonPath).toPath()));
        } else {
            respond(ex, 500, "application/json",
                ("{\"error\":\"" + escJ(ur.error) + "\"}").getBytes());
        }
    }

    private void createLocus(HttpExchange ex) throws IOException {
        cors(ex); if (preflight(ex)) return;
        if (!"POST".equalsIgnoreCase(ex.getRequestMethod())) {
            respond(ex, 405, "text/plain", "Method Not Allowed".getBytes()); return;
        }
        if (config == null || gff == null || loci == null || outputs == null) {
            respond(ex, 503, "application/json",
                "{\"error\":\"Pipeline state not available\"}".getBytes()); return;
        }
        String body  = new String(readAll(ex.getRequestBody()), "UTF-8");
        String chr   = extractStr(body, "chr");
        String sStart = extractStr(body, "start");
        String sEnd   = extractStr(body, "end");
        String name   = extractStr(body, "locus_name");
        if (chr == null || sStart == null || sEnd == null) {
            respond(ex, 400, "application/json",
                "{\"error\":\"missing chr, start, or end\"}".getBytes()); return;
        }
        LocusUpdater.UpdateResult ur = LocusUpdater.create(
            chr, Long.parseLong(sStart), Long.parseLong(sEnd), name,
            config, gff, loci, outputs);
        if (ur.ok) {
            respond(ex, 200, "application/json",
                Files.readAllBytes(new File(outputDir + "/data/manifest.json").toPath()));
        } else {
            respond(ex, 500, "application/json",
                ("{\"error\":\"" + escJ(ur.error) + "\"}").getBytes());
        }
    }

    private void splitLocus(HttpExchange ex) throws IOException {
        cors(ex); if (preflight(ex)) return;
        if (!"POST".equalsIgnoreCase(ex.getRequestMethod())) {
            respond(ex, 405, "text/plain", "Method Not Allowed".getBytes()); return;
        }
        if (config == null || gff == null || loci == null || outputs == null) {
            respond(ex, 503, "application/json",
                "{\"error\":\"Pipeline state not available\"}".getBytes()); return;
        }
        String body = new String(readAll(ex.getRequestBody()), "UTF-8");
        String sOrigIdx = extractStr(body, "locus_index");
        if (sOrigIdx == null) {
            respond(ex, 400, "application/json", "{\"error\":\"missing locus_index\"}".getBytes()); return;
        }
        List<LocusUpdater.SplitRegion> regions = parseSplitRegions(body);
        if (regions.size() < 2) {
            respond(ex, 400, "application/json", "{\"error\":\"need at least 2 regions\"}".getBytes()); return;
        }
        LocusUpdater.UpdateResult ur = LocusUpdater.split(
            Integer.parseInt(sOrigIdx), regions, config, gff, loci, outputs);
        if (ur.ok) {
            respond(ex, 200, "application/json",
                Files.readAllBytes(new File(outputDir + "/data/manifest.json").toPath()));
        } else {
            respond(ex, 500, "application/json",
                ("{\"error\":\"" + escJ(ur.error) + "\"}").getBytes());
        }
    }

    private void validateSplit(HttpExchange ex) throws IOException {
        cors(ex); if (preflight(ex)) return;
        if (!"POST".equalsIgnoreCase(ex.getRequestMethod())) {
            respond(ex, 405, "text/plain", "Method Not Allowed".getBytes()); return;
        }
        if (config == null || loci == null) {
            respond(ex, 503, "application/json",
                "{\"error\":\"Pipeline state not available\"}".getBytes()); return;
        }
        String body = new String(readAll(ex.getRequestBody()), "UTF-8");
        String sOrigIdx = extractStr(body, "locus_index");
        if (sOrigIdx == null) {
            respond(ex, 400, "application/json", "{\"error\":\"missing locus_index\"}".getBytes()); return;
        }
        List<LocusUpdater.SplitRegion> regions = parseSplitRegions(body);
        LocusUpdater.ValidationResult vr = LocusUpdater.validateSplit(
            Integer.parseInt(sOrigIdx), regions, config, loci);
        StringBuilder json = new StringBuilder();
        json.append("{\"valid\":").append(vr.valid);
        json.append(",\"split_ld_threshold\":").append(config.splitLdThreshold);
        json.append(",\"split_min_distance_bp\":").append(config.splitMinDistBp);
        json.append(",\"warnings\":[");
        for (int i = 0; i < vr.warnings.size(); i++) {
            if (i > 0) json.append(',');
            json.append('"').append(escJ(vr.warnings.get(i))).append('"');
        }
        json.append("],\"ld_violations\":[");
        for (int i = 0; i < vr.ldViolations.size(); i++) {
            if (i > 0) json.append(',');
            LocusUpdater.LdViolation v = vr.ldViolations.get(i);
            json.append("{\"snp_a\":\"").append(escJ(v.snpA))
                .append("\",\"snp_b\":\"").append(escJ(v.snpB))
                .append("\",\"pos_a\":").append(v.posA)
                .append(",\"pos_b\":").append(v.posB)
                .append(",\"region_a\":").append(v.regionA)
                .append(",\"region_b\":").append(v.regionB)
                .append(",\"r2\":").append(String.format("%.4f", v.r2))
                .append('}');
        }
        json.append("]}");
        respond(ex, 200, "application/json", json.toString().getBytes());
    }

    private void staticFiles(HttpExchange ex) throws IOException {
        cors(ex);
        String p = ex.getRequestURI().getPath();
        if (p.equals("/")) p = "/index.html";
        // Serve from project root first (index.html, viewer.html, assets/)
        File f = new File("." + p);
        if (!f.exists() || !f.isFile()) {
            // Fall back to outputDir for legacy single-project files
            f = new File(outputDir + p);
        }
        if (!f.exists() || !f.isFile()) {
            respond(ex, 404, "text/plain", ("Not found: " + p).getBytes()); return;
        }
        byte[] bytes = Files.readAllBytes(f.toPath());
        ex.getResponseHeaders().set("Content-Type", guessMime(p));
        ex.sendResponseHeaders(200, bytes.length);
        try (OutputStream os = ex.getResponseBody()) { os.write(bytes); }
    }

    private void pickFile(HttpExchange ex) throws IOException {
        cors(ex); if (preflight(ex)) return;
        String chosen = filePicker(".");
        String json = chosen != null
            ? "{\"path\":\"" + escJ(chosen) + "\"}"
            : "{\"path\":null}";
        respond(ex, 200, "application/json", json.getBytes());
    }

    // ── Helpers ────────────────────────────────────────────────────────────

    private void serveJson(HttpExchange ex, String path) throws IOException {
        File f = new File(path);
        if (!f.exists()) { respond(ex, 404, "application/json", "{\"error\":\"not found\"}".getBytes()); return; }
        respond(ex, 200, "application/json", Files.readAllBytes(f.toPath()));
    }

    private void respond(HttpExchange ex, int code, String mime, byte[] body) throws IOException {
        ex.getResponseHeaders().set("Content-Type", mime);
        ex.sendResponseHeaders(code, body.length);
        try (OutputStream os = ex.getResponseBody()) { os.write(body); }
    }

    private void cors(HttpExchange ex) {
        ex.getResponseHeaders().add("Access-Control-Allow-Origin",  "*");
        ex.getResponseHeaders().add("Access-Control-Allow-Methods", "GET, POST, DELETE, OPTIONS");
        ex.getResponseHeaders().add("Access-Control-Allow-Headers", "Content-Type");
    }

    private boolean preflight(HttpExchange ex) throws IOException {
        if ("OPTIONS".equalsIgnoreCase(ex.getRequestMethod())) {
            ex.sendResponseHeaders(204, -1); return true;
        }
        return false;
    }

    private static byte[] readAll(InputStream is) throws IOException {
        ByteArrayOutputStream buf = new ByteArrayOutputStream();
        byte[] tmp = new byte[8192]; int n;
        while ((n = is.read(tmp)) != -1) buf.write(tmp, 0, n);
        return buf.toByteArray();
    }

    private static String escJ(String s) {
        if (s == null) return "";
        return s.replace("\\", "\\\\").replace("\"", "\\\"");
    }

    private static String extractStr(String json, String key) {
        String marker = "\"" + key + "\":";
        int i = json.indexOf(marker);
        if (i < 0) return null;
        i += marker.length();
        while (i < json.length() && json.charAt(i) == ' ') i++;
        if (i >= json.length()) return null;
        if (json.charAt(i) == '"') {
            StringBuilder sb = new StringBuilder(); i++;
            while (i < json.length()) {
                char c = json.charAt(i++);
                if (c == '"') break;
                if (c == '\\' && i < json.length()) { sb.append(json.charAt(i++)); continue; }
                sb.append(c);
            }
            return sb.toString();
        }
        if (json.startsWith("null", i)) return null;
        int end = i;
        while (end < json.length() && ",}]".indexOf(json.charAt(end)) < 0) end++;
        return json.substring(i, end).trim();
    }

    private static String extractJsonString(String json, String key) {
        String search = "\"" + key + "\":\"";
        int idx = json.indexOf(search);
        if (idx < 0) return null;
        idx += search.length();
        StringBuilder sb = new StringBuilder();
        for (int i = idx; i < json.length(); i++) {
            char c = json.charAt(i);
            if (c == '"') break;
            if (c == '\\' && i + 1 < json.length()) { sb.append(json.charAt(++i)); continue; }
            sb.append(c);
        }
        return sb.toString();
    }

    private List<LocusUpdater.SplitRegion> parseSplitRegions(String body) {
        List<LocusUpdater.SplitRegion> regions = new ArrayList<>();
        int ri = body.indexOf("\"regions\"");
        if (ri < 0) return regions;
        int arrStart = body.indexOf('[', ri);
        int arrEnd   = body.indexOf(']', arrStart);
        if (arrStart < 0 || arrEnd < 0) return regions;
        String arrStr = body.substring(arrStart, arrEnd + 1);
        String[] parts = arrStr.replaceAll("^\\[\\{", "").replaceAll("\\}\\]$", "").split("\\},\\{");
        for (String part : parts) {
            String obj = "{" + part + "}";
            String name  = extractStr(obj, "name");
            String start = extractStr(obj, "start");
            String end   = extractStr(obj, "end");
            if (start == null || end == null) continue;
            regions.add(new LocusUpdater.SplitRegion(
                name != null ? name : "", Long.parseLong(start), Long.parseLong(end)));
        }
        return regions;
    }

    private static void deleteDirectory(File dir) {
        File[] files = dir.listFiles();
        if (files != null) {
            for (File f : files) {
                if (f.isDirectory()) deleteDirectory(f);
                else f.delete();
            }
        }
        dir.delete();
    }

    private static String guessMime(String path) {
        if (path.endsWith(".html")) return "text/html";
        if (path.endsWith(".js"))   return "application/javascript";
        if (path.endsWith(".json")) return "application/json";
        if (path.endsWith(".css"))  return "text/css";
        if (path.endsWith(".png"))  return "image/png";
        if (path.endsWith(".svg"))  return "image/svg+xml";
        if (path.endsWith(".yaml") || path.endsWith(".yml")) return "text/yaml";
        if (path.endsWith(".tsv"))  return "text/tab-separated-values";
        return "application/octet-stream";
    }

    private static String filePicker(String startDir) {
        final String[] result = {null};
        try {
            Thread edtThread = new Thread(() -> {
                try {
                    SwingUtilities.invokeAndWait(() -> {
                        try { UIManager.setLookAndFeel(UIManager.getSystemLookAndFeelClassName()); } catch (Exception ignore) {}
                        // Create an always-on-top frame so the dialog appears in front of the browser
                        JFrame frame = new JFrame();
                        frame.setAlwaysOnTop(true);
                        frame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
                        // Make it tiny and invisible but focusable
                        frame.setSize(0, 0);
                        frame.setLocationRelativeTo(null);
                        frame.setVisible(true);
                        frame.toFront();
                        frame.requestFocus();

                        JFileChooser fc = new JFileChooser(startDir);
                        fc.setFileSelectionMode(JFileChooser.FILES_ONLY);
                        fc.setDialogTitle("Select file");
                        fc.setAcceptAllFileFilterUsed(true);
                        fc.addChoosableFileFilter(new javax.swing.filechooser.FileNameExtensionFilter(
                            "Data files (*.tsv, *.csv, *.txt, *.gz, *.gff3, *.bed)", "tsv", "csv", "txt", "gz", "gff3", "bed", "bim", "fam"));
                        if (fc.showOpenDialog(frame) == JFileChooser.APPROVE_OPTION)
                            result[0] = fc.getSelectedFile().getAbsolutePath();
                        frame.dispose();
                    });
                } catch (Exception e) {
                    System.err.println("[LocalServer] File picker EDT error: " + e.getMessage());
                }
            });
            edtThread.setDaemon(true);
            edtThread.start();
            edtThread.join(120_000);
        } catch (Exception e) {
            System.err.println("[LocalServer] File picker error: " + e.getMessage());
        }
        return result[0];
    }

    private static String folderPicker(String startDir) {
        final String[] result = {null};
        try {
            SwingUtilities.invokeAndWait(() -> {
                try { UIManager.setLookAndFeel(UIManager.getSystemLookAndFeelClassName()); } catch (Exception ignore) {}
                JFrame frame = new JFrame();
                frame.setAlwaysOnTop(true);
                frame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
                frame.setSize(0, 0);
                frame.setLocationRelativeTo(null);
                frame.setVisible(true);
                frame.toFront();
                frame.requestFocus();

                JFileChooser fc = new JFileChooser(startDir);
                fc.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
                fc.setDialogTitle("Select folder");
                if (fc.showSaveDialog(frame) == JFileChooser.APPROVE_OPTION)
                    result[0] = fc.getSelectedFile().getAbsolutePath();
                frame.dispose();
            });
        } catch (Exception e) {
            System.err.println("[LocalServer] Folder picker error: " + e.getMessage());
        }
        return result[0];
    }
}
