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
import catalog.*;
import opentargets.*;

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

    // GWAS Catalog cross-phenotype lookup — one shared local snapshot for all projects, not
    // per-project config (mirrors config/global.json's reference-panel/SNP-database registry;
    // TODO: move this path into global.json once the snapshot location needs to be user-configurable).
    private static final String GWAS_CATALOG_TSV = "gwascatalog_data/gwas-catalog-download-associations-alt-full.tsv";
    private static final String GWAS_CATALOG_IDX_DIR = "gwascatalog_data";
    private static volatile GwasCatalogLocalIndex gwasCatalogIndex;
    private static volatile String gwasCatalogInitError;

    /** Lazily builds/loads the local GWAS Catalog index on first request. Returns null if unavailable. */
    private static GwasCatalogLocalIndex getGwasCatalogIndex() {
        if (gwasCatalogIndex != null || gwasCatalogInitError != null) return gwasCatalogIndex;
        synchronized (LocalServer.class) {
            if (gwasCatalogIndex != null || gwasCatalogInitError != null) return gwasCatalogIndex;
            if (!new File(GWAS_CATALOG_TSV).isFile()) {
                gwasCatalogInitError = "GWAS Catalog snapshot not found at " + GWAS_CATALOG_TSV +
                    " — download gwas-catalog-associations_ontology-annotated-full.zip from " +
                    "ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/ and unzip it there.";
                System.err.println("[LocalServer] " + gwasCatalogInitError);
                return null;
            }
            try {
                GwasCatalogLocalIndex idx = new GwasCatalogLocalIndex(GWAS_CATALOG_TSV, GWAS_CATALOG_IDX_DIR);
                if (!new File(GWAS_CATALOG_IDX_DIR, "rsid_index.tsv").isFile() ||
                    !new File(GWAS_CATALOG_IDX_DIR, "region_index.tsv").isFile()) {
                    idx.buildIndex();
                }
                idx.loadIndex();
                gwasCatalogIndex = idx;
            } catch (IOException e) {
                gwasCatalogInitError = "Failed to load GWAS Catalog index: " + e.getMessage();
                System.err.println("[LocalServer] " + gwasCatalogInitError);
            }
        }
        return gwasCatalogIndex;
    }

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

    private ProjectState ensureProjectState(String projectId, String projectDir) {
        ProjectState ps = projectStates.get(projectId);
        if (ps != null && ps.config != null && ps.loci != null) return ps;
        try {
            Config cfg = Config.loadFromProject(projectDir);
            GffParser gff = GffParser.parse(cfg);
            List<Locus> loci = new ArrayList<>();
            List<LocusOutput> outputs = new ArrayList<>();
            // Reconstruct state from manifest + locus JSON files
            File manifestFile = new File(projectDir, "data/manifest.json");
            if (manifestFile.exists()) {
                String mJson = new String(Files.readAllBytes(manifestFile.toPath()), "UTF-8");
                List<Integer> indices = new ArrayList<>();
                int lociArr = mJson.indexOf("\"loci\"");
                if (lociArr >= 0) {
                    int s = mJson.indexOf('[', lociArr);
                    int e = mJson.lastIndexOf(']');
                    if (s >= 0 && e > s) {
                        String arr = mJson.substring(s + 1, e);
                        int pos = 0;
                        while (pos < arr.length()) {
                            int ik = arr.indexOf("\"index\":", pos);
                            if (ik < 0) break;
                            int vs = ik + 8;
                            while (vs < arr.length() && arr.charAt(vs) == ' ') vs++;
                            int ve = vs;
                            while (ve < arr.length() && Character.isDigit(arr.charAt(ve))) ve++;
                            if (ve > vs) indices.add(Integer.parseInt(arr.substring(vs, ve)));
                            pos = ve;
                        }
                    }
                }
                for (int idx : indices) {
                    File lf = new File(projectDir, "data/locus_" + idx + ".json");
                    if (!lf.exists()) continue;
                    String json = new String(Files.readAllBytes(lf.toPath()), "UTF-8");
                    String chr = extractStr(json, "chr");
                    String startS = extractStr(json, "start");
                    String endS = extractStr(json, "end");
                    String id = extractStr(json, "id");
                    String name = extractStr(json, "locus_name");
                    if (chr == null || startS == null || endS == null) continue;
                    Locus l = id != null
                        ? new Locus(id, idx, chr, Long.parseLong(startS), Long.parseLong(endS), cfg.locusPadding)
                        : new Locus(idx, chr, Long.parseLong(startS), Long.parseLong(endS), cfg.locusPadding);
                    loci.add(l);
                    LocusOutput lo = new LocusOutput();
                    lo.id = l.id; lo.locusIndex = idx;
                    lo.locusName = (name != null && !name.isEmpty()) ? name : "Locus " + idx;
                    lo.chr = chr; lo.start = l.start; lo.end = l.end;
                    lo.paddedStart = l.paddedStart; lo.paddedEnd = l.paddedEnd;
                    lo.nearestGenes = extractStringArray(json, "nearest_genes");
                    lo.genes = gff.overlapping(chr, l.paddedStart, l.paddedEnd);

                    // Restore ref panel label for manifest accuracy after a lazy reload
                    // (LocusOutput.refPanel defaults to "" otherwise).
                    String refPanel = extractStr(json, "ref_panel");
                    if (refPanel != null) lo.refPanel = refPanel;

                    // Parse top_snp so manifest summaries (top_snp/top_snp_pval) stay
                    // correct after a live mutation (create/split/delete/merge/reorder)
                    // that re-exports the manifest from this reconstructed state.
                    int topSnpIdx = json.indexOf("\"top_snp\":");
                    if (topSnpIdx >= 0) {
                        String topId  = extractNestedStr(json, topSnpIdx, "id");
                        String topPos = extractNestedStr(json, topSnpIdx, "pos");
                        String topP   = extractNestedStr(json, topSnpIdx, "pvalue");
                        if (topId != null && topPos != null && topP != null) {
                            try {
                                lo.topSnp = new Snp(topId, chr, Long.parseLong(topPos),
                                    Double.parseDouble(topP), "", "");
                            } catch (NumberFormatException ignore) {}
                        }
                    }

                    outputs.add(lo);
                }
            }
            loci.sort(Comparator.comparingInt(Locus::chrInt).thenComparingLong(l -> l.start));
            outputs.sort(Comparator.comparingInt((LocusOutput o) -> Locus.chrToInt(o.chr))
                .thenComparingLong(o -> o.start));
            setProjectState(projectId, cfg, gff, loci, outputs);
            System.out.printf("[Server] Lazy-loaded project state for '%s': %d loci%n", projectId, loci.size());
            return projectStates.get(projectId);
        } catch (Exception e) {
            System.err.printf("[Server] Failed to lazy-load project '%s': %s%n", projectId, e.getMessage());
            return null;
        }
    }

    private LociMutationService ensureMutationService(ProjectState ps) {
        if (ps.mutationService == null && ps.config != null && ps.gff != null
                && ps.loci != null && ps.outputs != null) {
            ps.mutationService = new LociMutationService(ps.config, ps.gff, ps.loci, ps.outputs);
        }
        return ps.mutationService;
    }

    public LocalServer(String outputDir) throws IOException {
        this.outputDir  = outputDir;
        this.lastFolder = outputDir + "/plots";
        new File(this.lastFolder).mkdirs();

        http = HttpServer.create(new InetSocketAddress("localhost", PORT), 32);
        loadPersistedLocusMatrixJobs();

        // ── Project-scoped API (Step 6) ──────────────────────────────────
        http.createContext("/api/projects/process", this::processProjects);
        http.createContext("/api/projects",         this::listProjects);
        http.createContext("/api/delete-project",    this::projectDelete);
        http.createContext("/api/project/",         this::projectRouter);
        http.createContext("/api/peek-file-header", this::peekFileHeader);

        // ── rsID recovery API ────────────────────────────────────────────
        http.createContext("/api/global-config",     this::globalConfigEndpoint);
        http.createContext("/api/shared-storage/test", this::sharedStorageTest);
        http.createContext("/api/rsid-recover",      this::rsidRecover);
        http.createContext("/api/rsid-progress",     this::rsidProgressEndpoint);
        http.createContext("/api/rsid-detect",       this::rsidDetect);
        http.createContext("/api/missing-rsids",     this::missingRsids);
        http.createContext("/api/submit-rsid",       this::submitRsid);
        http.createContext("/api/loci-identify",     this::lociIdentify);
        http.createContext("/api/loci-progress",     this::lociProgressEndpoint);
        http.createContext("/api/export-excel",      this::exportExcel);
        http.createContext("/api/export-progress",   this::exportProgressEndpoint);
        http.createContext("/api/export-projects-info", this::exportProjectsInfo);
        http.createContext("/api/locus-matrix-run",      this::locusMatrixRun);
        http.createContext("/api/locus-matrix-progress", this::locusMatrixProgressEndpoint);
        http.createContext("/api/locus-matrix-jobs",     this::locusMatrixJobsList);
        http.createContext("/api/locus-matrix-result",   this::locusMatrixResult);
        http.createContext("/api/locus-matrix-export",   this::locusMatrixExport);
        http.createContext("/api/locus-matrix-delete",   this::locusMatrixDelete);
        http.createContext("/api/gene-constellation",    this::geneConstellation);
        http.createContext("/api/gene-constellation-export", this::geneConstellationExport);
        http.createContext("/api/serpent-plot",          this::serpentPlot);
        http.createContext("/api/serpent-locus",         this::serpentLocus);
        http.createContext("/api/susiex-run",        this::susiexRun);
        http.createContext("/api/susiex-progress",   this::susiexProgress);
        http.createContext("/api/susiex-result",     this::susiexResult);
        http.createContext("/api/mvsusie-run",       this::mvsusieRun);
        http.createContext("/api/mvsusie-progress",  this::mvsusieProgress);
        http.createContext("/api/mvsusie-result",    this::mvsusieResult);
        http.createContext("/api/mvsusie-latest",    this::mvsusieLatest);
        http.createContext("/api/mvsusie-for-dataset", this::mvsusieForDataset);
        http.createContext("/api/mvsusie-run-all",   this::mvsusieRunAll);
        http.createContext("/api/mvsusie-run-all-progress", this::mvsusieRunAllProgress);
        http.createContext("/api/mvsusie-run-all-cancel",   this::mvsusieRunAllCancel);
        http.createContext("/api/mvsusie-summary",   this::mvsusieSummary);
        http.createContext("/api/search/rebuild-index", this::searchRebuildIndex);
        http.createContext("/api/search",            this::globalSearch);
        http.createContext("/pick-folder-native",    this::pickFolder);

        // ── AI Agent panel (Phase 4) ──────────────────────────────────────
        http.createContext("/api/agent/chat",          this::agentChat);
        http.createContext("/api/agent/config",        this::agentConfigEndpoint);
        http.createContext("/api/agent/detect-local",  this::agentDetectLocal);
        http.createContext("/api/agent/pull-model",    this::agentPullModel);
        http.createContext("/api/agent/pull-progress", this::agentPullProgress);

        // ── GitHub project save/extract (Phase 6) ──────────────────────────
        http.createContext("/api/github/token-status", this::githubTokenStatus);
        http.createContext("/api/github/token",        this::githubToken);
        http.createContext("/api/github/import",       this::githubImport);

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
        } else if (action.equals("genome-skyline")) {
            serveJson(ex, dataDir + "/genome_skyline.json");
        } else if (action.matches("locus/\\d+/catalog")) {
            String n = action.replaceAll("\\D+", "");
            projectLocusCatalog(ex, dataDir, n);
        } else if (action.matches("locus/\\d+/novelty")) {
            String n = action.replaceAll("\\D+", "");
            projectLocusNovelty(ex, dataDir, projectDir, n);
        } else if (action.equals("catalog-report")) {
            projectCatalogReport(ex, dataDir, projectDir);
        } else if (action.matches("locus/\\d+/l2g")) {
            // NOT replaceAll("\\D+","") like the catalog/novelty routes above — "l2g" itself contains
            // a digit, so stripping all non-digits from e.g. "locus/2/l2g" corrupts "2" into "22".
            String n = action.replaceFirst("locus/(\\d+)/l2g", "$1");
            projectLocusL2G(ex, dataDir, projectDir, n);
        } else if (action.matches("locus/\\d+/regulatory")) {
            String n = action.replaceAll("\\D+", "");
            projectLocusRegulatory(ex, dataDir, projectId, n);
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
        } else if (action.equals("delete-locus")) {
            projectDeleteLocus(ex, projectId, projectDir);
        } else if (action.equals("merge-loci")) {
            projectMergeLoci(ex, projectId, projectDir);
        } else if (action.equals("undo-mutation")) {
            projectUndoMutation(ex, projectId, projectDir);
        } else if (action.equals("reorder-loci")) {
            projectReorderLoci(ex, projectId, projectDir);
        } else if (action.equals("analysis/base-status")) {
            projectAnalysisBaseStatus(ex, projectId, projectDir);
        } else if (action.equals("analysis/build-base")) {
            projectAnalysisBuildBase(ex, projectId, projectDir);
        } else if (action.equals("analysis/build-base-progress")) {
            projectAnalysisBuildBaseProgress(ex, projectId);
        } else if (action.startsWith("analysis/locus-log/")) {
            projectLocusBuildLog(ex, projectDir, action);
        } else if (action.equals("analysis/tools")) {
            projectAnalysisTools(ex);  // doesn't need project state
        } else if (action.equals("analysis/latest")) {
            projectAnalysisLatest(ex, projectId, projectDir);
        } else if (action.equals("analysis/run-result")) {
            projectAnalysisRunResult(ex, projectId, projectDir);
        } else if (action.equals("analysis/summary")) {
            projectAnalysisSummary(ex, projectId, projectDir);
        } else if (action.equals("analysis/run-all")) {
            projectAnalysisRunAll(ex, projectId, projectDir);
        } else if (action.equals("analysis/run-all-progress")) {
            projectAnalysisRunAllProgress(ex, projectId);
        } else if (action.equals("analysis/run-all-cancel")) {
            BatchRun br = batchRuns.get(projectId + "|" + queryParam(ex, "tool"));
            if (br != null) br.cancelled = true;
            respond(ex, 200, "application/json", "{\"ok\":true}".getBytes());
        } else if (action.equals("analysis/run")) {
            projectAnalysisRun(ex, projectId, projectDir);
        } else if (action.startsWith("analysis/job/")) {
            projectAnalysisJob(ex, projectId, projectDir, action);
        } else if (action.equals("analysis/history")) {
            projectAnalysisHistory(ex, projectId, projectDir);
        } else if (action.equals("progress")) {
            projectProgressEndpoint(ex, projectId);
        } else if (action.equals("qc")) {
            projectQc(ex, projectId, projectDir);
        } else if (action.equals("evidence")) {
            if ("POST".equalsIgnoreCase(ex.getRequestMethod())) {
                projectEvidenceUpload(ex, projectId, projectDir);
            } else {
                projectEvidenceList(ex, projectDir);
            }
        } else if (action.startsWith("evidence/")) {
            String evName = action.substring("evidence/".length());
            if ("DELETE".equalsIgnoreCase(ex.getRequestMethod())) {
                projectEvidenceDelete(ex, projectDir, evName);
            } else {
                projectEvidenceGet(ex, projectDir, evName);
            }
        } else if (action.equals("enrichment")) {
            projectEnrichment(ex, projectId, projectDir);
        } else if (action.equals("regulatory-enrichment")) {
            projectRegulatoryEnrichment(ex, projectId, projectDir);
        } else if (action.equals("github/push")) {
            projectGithubPush(ex, projectDir);
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

    // GET /api/project/{id}/qc — genomic inflation (lambda_GC), cached by GWAS-file content hash;
    // optionally accepts a POST with a user-supplied LDSC intercept to compute the attenuation ratio.
    private void projectQc(HttpExchange ex, String projectId, String projectDir) throws IOException {
        Config cfg;
        try {
            cfg = Config.loadFromProject(projectDir);
        } catch (Exception e) {
            respond(ex, 404, "application/json", ("{\"error\":\"" + escJ(e.getMessage()) + "\"}").getBytes());
            return;
        }

        Double ldscIntercept = null;
        if ("POST".equalsIgnoreCase(ex.getRequestMethod())) {
            String body = new String(readAll(ex.getRequestBody()), "UTF-8");
            String v = extractStr(body, "ldsc_intercept");
            if (v != null && !v.trim().isEmpty()) {
                try { ldscIntercept = Double.parseDouble(v.trim()); } catch (NumberFormatException ignored) {}
            }
        }

        File cacheFile = new File(projectDir, "data/qc.json");
        String gwasHash;
        try { gwasHash = ContentHasher.hashFile(new File(cfg.gwasFile)); }
        catch (IOException e) {
            respond(ex, 400, "application/json", ("{\"error\":\"GWAS file not found: " + escJ(cfg.gwasFile) + "\"}").getBytes());
            return;
        }

        GwasQc.Result r = null;
        if (cacheFile.exists()) {
            try {
                String cached = new String(Files.readAllBytes(cacheFile.toPath()), "UTF-8");
                // "n_snps_scanned" was added when distance-based pruning was introduced; a cache
                // written before that fix has the same gwas_hash (the GWAS file didn't change) but
                // holds an unpruned lambda_GC — without this check it would be served forever.
                boolean cacheIsCurrentFormat = cached.contains("\"n_snps_scanned\"");
                if (cacheIsCurrentFormat && gwasHash.equals(extractStr(cached, "gwas_hash"))) {
                    if (ldscIntercept == null) {
                        respond(ex, 200, "application/json", cached.getBytes("UTF-8"));
                        return;
                    }
                    r = new GwasQc.Result();
                    r.ok = true;
                    r.nSnps = (int) parseDoubleOr(extractStr(cached, "n_snps"), 0);
                    r.nSnpsScanned = (int) parseDoubleOr(extractStr(cached, "n_snps_scanned"), 0);
                    r.medianChi2 = parseDoubleOr(extractStr(cached, "median_chi2"), Double.NaN);
                    r.meanChi2 = parseDoubleOr(extractStr(cached, "mean_chi2"), Double.NaN);
                    r.lambdaGC = parseDoubleOr(extractStr(cached, "lambda_gc"), Double.NaN);
                    r.lambdaFlagged = "true".equals(extractStr(cached, "lambda_flagged"));
                }
            } catch (Exception ignored) {}
        }
        if (r == null) r = GwasQc.compute(cfg);
        StringBuilder json = new StringBuilder(r.toJson());
        json.setLength(json.length() - 1); // drop closing brace, append more fields
        json.append(",\"gwas_hash\":\"").append(escJ(gwasHash)).append('"');
        if (ldscIntercept != null && r.ok) {
            double ratio = GwasQc.attenuationRatio(ldscIntercept, r.meanChi2);
            json.append(",\"ldsc_intercept\":").append(ldscIntercept);
            json.append(",\"ldsc_flagged\":").append(ldscIntercept > GwasQc.LDSC_INTERCEPT_FLAG);
            json.append(",\"attenuation_ratio\":").append(Double.isNaN(ratio) ? "null" : ratio);
        }
        json.append('}');

        try {
            new File(projectDir, "data").mkdirs();
            Files.write(cacheFile.toPath(), json.toString().getBytes("UTF-8"));
        } catch (IOException ignored) {}

        respond(ex, r.ok ? 200 : 400, "application/json", json.toString().getBytes("UTF-8"));
    }

    // GET /api/project/{id}/locus/{n}/catalog — cross-phenotype lookup for this locus's top SNP
    // against a local GWAS Catalog snapshot (see GwasCatalogLocalIndex). No network call at request
    // time — the snapshot + rsID index are built once, offline, ahead of server start.
    private void projectLocusCatalog(HttpExchange ex, String dataDir, String locusIndex) throws IOException {
        File locusFile = new File(dataDir, "locus_" + locusIndex + ".json");
        if (!locusFile.exists()) {
            respond(ex, 404, "application/json", "{\"error\":\"Locus not found\"}".getBytes());
            return;
        }
        String locusJson = new String(Files.readAllBytes(locusFile.toPath()), "UTF-8");
        int topSnpIdx = locusJson.indexOf("\"top_snp\":");
        String rsid = topSnpIdx >= 0 ? extractNestedStr(locusJson, topSnpIdx, "id") : null;
        if (rsid == null || rsid.isEmpty() || !rsid.startsWith("rs")) {
            respond(ex, 200, "application/json",
                "{\"ok\":false,\"error\":\"Locus has no rsID-form top SNP to look up\"}".getBytes());
            return;
        }

        GwasCatalogLocalIndex idx = getGwasCatalogIndex();
        if (idx == null) {
            respond(ex, 503, "application/json",
                ("{\"ok\":false,\"error\":\"" + escJ(gwasCatalogInitError) + "\"}").getBytes());
            return;
        }

        List<GwasCatalogLocalIndex.Hit> hits;
        try {
            hits = idx.lookup(rsid);
        } catch (IOException e) {
            respond(ex, 500, "application/json",
                ("{\"ok\":false,\"error\":\"" + escJ(e.getMessage()) + "\"}").getBytes());
            return;
        }

        StringBuilder json = new StringBuilder();
        json.append("{\"ok\":true,\"rsid\":\"").append(escJ(rsid)).append("\",\"hits\":[");
        for (int i = 0; i < hits.size(); i++) {
            if (i > 0) json.append(',');
            json.append(hitToJson(hits.get(i)));
        }
        json.append("]}");
        respond(ex, 200, "application/json", json.toString().getBytes("UTF-8"));
    }

    // GET /api/project/{id}/locus/{n}/l2g — most likely causal GENE for this locus's top SNP, via
    // Open Targets' Locus-to-Gene model (live API — L2G is only published as Parquet in bulk, which
    // would need a real Parquet-decoding dependency this codebase doesn't have; the live query is a
    // single round trip per locus, cached to disk, so it never repeats the GWAS-Catalog-REST N+1 problem).
    // Returns an empty gene list (ok:true) rather than an error for any locus Open Targets hasn't seen,
    // including every genuinely novel locus by construction.
    private void projectLocusL2G(HttpExchange ex, String dataDir, String projectDir, String locusIndex) throws IOException {
        File locusFile = new File(dataDir, "locus_" + locusIndex + ".json");
        if (!locusFile.exists()) {
            respond(ex, 404, "application/json", "{\"error\":\"Locus not found\"}".getBytes());
            return;
        }
        String locusJson = new String(Files.readAllBytes(locusFile.toPath()), "UTF-8");
        int topSnpIdx = locusJson.indexOf("\"top_snp\":");
        String rsid = topSnpIdx >= 0 ? extractNestedStr(locusJson, topSnpIdx, "id") : null;
        if (rsid == null || rsid.isEmpty() || !rsid.startsWith("rs")) {
            respond(ex, 200, "application/json",
                "{\"ok\":false,\"error\":\"Locus has no rsID-form top SNP to look up\"}".getBytes());
            return;
        }

        OpenTargetsL2GClient client = new OpenTargetsL2GClient(projectDir);
        OpenTargetsL2GClient.Result r;
        try {
            r = client.lookup(rsid);
        } catch (IOException | InterruptedException e) {
            respond(ex, 502, "application/json",
                ("{\"ok\":false,\"error\":\"Open Targets lookup failed: " + escJ(e.getMessage()) + "\"}").getBytes());
            return;
        }

        StringBuilder json = new StringBuilder();
        json.append("{\"ok\":true,\"rsid\":\"").append(escJ(r.rsid)).append('"');
        json.append(",\"variant_id\":").append(r.variantId == null ? "null" : "\"" + escJ(r.variantId) + "\"");
        json.append(",\"credible_set_count\":").append(r.credibleSetCount);
        json.append(",\"genes\":[");
        for (int i = 0; i < r.genes.size(); i++) {
            if (i > 0) json.append(',');
            OpenTargetsL2GClient.GenePrediction gp = r.genes.get(i);
            json.append("{\"gene\":\"").append(escJ(gp.gene)).append('"');
            json.append(",\"score\":").append(gp.maxScore);
            json.append(",\"supporting_studies\":").append(gp.supportingStudies);
            json.append('}');
        }
        json.append(']');
        json.append(",\"enhancer_genes\":[");
        for (int i = 0; i < r.enhancerGenes.size(); i++) {
            if (i > 0) json.append(',');
            OpenTargetsL2GClient.EnhancerGenePrediction egp = r.enhancerGenes.get(i);
            json.append("{\"gene\":\"").append(escJ(egp.gene)).append('"');
            json.append(",\"score\":").append(egp.maxScore);
            json.append(",\"evidence\":[");
            for (int k = 0; k < egp.evidence.size(); k++) {
                if (k > 0) json.append(',');
                OpenTargetsL2GClient.EnhancerEvidence ev = egp.evidence.get(k);
                json.append("{\"biosample\":\"").append(escJ(ev.biosample)).append('"');
                json.append(",\"score\":").append(ev.score);
                json.append(",\"distance_to_tss\":").append(ev.distanceToTss);
                json.append(",\"pmid\":").append(ev.pmid == null ? "null" : "\"" + escJ(ev.pmid) + "\"");
                json.append('}');
            }
            json.append(']');
            json.append('}');
        }
        json.append("]}");
        respond(ex, 200, "application/json", json.toString().getBytes("UTF-8"));
    }

    // GET /api/project/{id}/locus/{n}/regulatory — Option C stacked regulatory tracks (Phase 3):
    // real histone ChIP-seq peaks (H3K27ac, H3K4me1, H3K4me3) from the disease-mapped reference
    // epigenome (RegulatoryPeakIndex), restricted to this locus's own already-computed padded
    // window — the SAME window the gene track already uses (padded_start/padded_end straight from
    // locus_N.json, not a newly invented padding value). Empty "marks" (still 200, ok, not an
    // error) for a project whose disease doesn't map to anything in DECISIONS_PHASE3.md's table.
    private void projectLocusRegulatory(HttpExchange ex, String dataDir, String projectId, String locusIndex) throws IOException {
        File locusFile = new File(dataDir, "locus_" + locusIndex + ".json");
        if (!locusFile.exists()) {
            respond(ex, 404, "application/json", "{\"error\":\"Locus not found\"}".getBytes());
            return;
        }
        String locusJson = new String(Files.readAllBytes(locusFile.toPath()), "UTF-8");
        String chr = extractStr(locusJson, "chr");
        String paddedStartStr = extractStr(locusJson, "padded_start");
        String paddedEndStr = extractStr(locusJson, "padded_end");

        String eid = RegulatoryPeakIndex.resolveEid(projectId);
        String tissue = eid == null ? null : RegulatoryPeakIndex.EID_TO_TISSUE.get(eid);

        StringBuilder json = new StringBuilder("{");
        json.append("\"eid\":").append(eid == null ? "null" : "\"" + escJ(eid) + "\"").append(",");
        json.append("\"tissue\":").append(tissue == null ? "null" : "\"" + escJ(tissue) + "\"").append(",");
        json.append("\"marks\":{");

        long start = -1, end = -1;
        if (paddedStartStr != null && paddedEndStr != null) {
            try {
                start = Long.parseLong(paddedStartStr);
                end = Long.parseLong(paddedEndStr);
            } catch (NumberFormatException nfe) { start = -1; end = -1; }
        }

        if (eid != null && chr != null && start >= 0 && end >= start) {
            RegulatoryPeakIndex idx = RegulatoryPeakIndex.instance();
            for (int mi = 0; mi < RegulatoryPeakIndex.MARKS.length; mi++) {
                String mark = RegulatoryPeakIndex.MARKS[mi];
                if (mi > 0) json.append(",");
                json.append("\"").append(mark).append("\":[");
                List<RegulatoryPeakIndex.Peak> peaks = idx.peaksOverlapping(eid, mark, chr, start, end);
                for (int pi = 0; pi < peaks.size(); pi++) {
                    if (pi > 0) json.append(",");
                    RegulatoryPeakIndex.Peak p = peaks.get(pi);
                    json.append("{\"start\":").append(p.start)
                        .append(",\"end\":").append(p.end)
                        .append(",\"signal\":").append(p.signalValue)
                        .append("}");
                }
                json.append("]");
            }
        }
        json.append("}}");
        respond(ex, 200, "application/json", json.toString().getBytes("UTF-8"));
    }

    // GET /api/project/{id}/regulatory-enrichment?threshold=5e-8 — real interval-overlap enrichment
    // (DECISIONS_PHASE3.md section 3.3/3.4): the project's own genome-wide-significant SNPs
    // (foreground) vs. its own below-significance SNPs (background), read from its own GWAS file
    // the same way GwasQc/MultiLocusScanner already do (BufferedReader + GwasParser.colIdx/splitTab,
    // MultiLocusScanner.normalizeChr), against the disease-mapped reference epigenome's real peaks.
    private void projectRegulatoryEnrichment(HttpExchange ex, String projectId, String projectDir) throws IOException {
        String thresholdParam = queryParam(ex, "threshold");
        double threshold = GeneConstellationBuilder.DEFAULT_THRESHOLD;
        if (thresholdParam != null && !thresholdParam.trim().isEmpty()) {
            try { threshold = Double.parseDouble(thresholdParam.trim()); } catch (NumberFormatException ignored) {}
        }

        Config cfg;
        try {
            cfg = Config.loadFromProject(projectDir);
        } catch (Exception e) {
            respond(ex, 404, "application/json",
                ("{\"error\":\"Project config not found: " + escJ(e.getMessage()) + "\"}").getBytes());
            return;
        }

        String eid = RegulatoryPeakIndex.resolveEid(projectId);
        String tissue = eid == null ? null : RegulatoryPeakIndex.EID_TO_TISSUE.get(eid);

        List<RegulatoryEnrichmentAnalyzer.SnpPos> foreground = new ArrayList<>();
        List<RegulatoryEnrichmentAnalyzer.SnpPos> background = new ArrayList<>();

        if (eid != null) {
            File gwasFile = new File(cfg.gwasFile);
            if (!gwasFile.isFile()) {
                respond(ex, 404, "application/json",
                    ("{\"error\":\"GWAS file not found: " + escJ(cfg.gwasFile) + "\"}").getBytes());
                return;
            }
            double sigThreshold = threshold;
            try (BufferedReader br = new BufferedReader(new FileReader(gwasFile), 1 << 20)) {
                String header = br.readLine();
                if (header != null) {
                    String[] cols = header.trim().split("\t");
                    int iChr = GwasParser.colIdx(cols, cfg.colChr);
                    int iPos = GwasParser.colIdx(cols, cfg.colPos);
                    int iP   = GwasParser.colIdx(cols, cfg.colPvalue);
                    if (iChr >= 0 && iPos >= 0 && iP >= 0) {
                        String line;
                        while ((line = br.readLine()) != null) {
                            if (line.isEmpty()) continue;
                            String[] f = GwasParser.splitTab(line);
                            if (f.length <= Math.max(iChr, Math.max(iPos, iP))) continue;
                            String chr = f[iChr].trim();
                            long pos; double p;
                            try {
                                pos = Long.parseLong(f[iPos].trim());
                                p   = Double.parseDouble(f[iP].trim());
                            } catch (NumberFormatException nfe) { continue; }
                            if (!(p > 0 && p < 1)) continue;
                            RegulatoryEnrichmentAnalyzer.SnpPos snp =
                                new RegulatoryEnrichmentAnalyzer.SnpPos(chr, pos);
                            if (p <= sigThreshold) foreground.add(snp); else background.add(snp);
                        }
                    }
                }
            }
        }

        RegulatoryEnrichmentAnalyzer.Result result = RegulatoryEnrichmentAnalyzer.run(
            eid, tissue, threshold, foreground, background, RegulatoryPeakIndex.instance());
        respond(ex, 200, "application/json", result.toJson().getBytes("UTF-8"));
    }

    private static String hitToJson(GwasCatalogLocalIndex.Hit h) {
        StringBuilder j = new StringBuilder();
        j.append("{\"trait\":\"").append(escJ(h.trait)).append('"');
        j.append(",\"mapped_trait\":\"").append(escJ(h.mappedTrait)).append('"');
        j.append(",\"mapped_gene\":\"").append(escJ(h.mappedGene)).append('"');
        j.append(",\"pvalue\":").append(Double.isFinite(h.pvalue) ? h.pvalue : 0);
        j.append(",\"pvalue_display\":\"").append(escJ(h.pvalueDisplay())).append('"');
        j.append(",\"study_accession\":\"").append(escJ(h.studyAccession)).append('"');
        j.append(",\"pubmed_id\":\"").append(escJ(h.pubmedId)).append('"');
        j.append(",\"initial_sample_size\":\"").append(escJ(h.initialSampleSize)).append('"');
        j.append(",\"link\":\"").append(escJ(h.link)).append('"');
        j.append(",\"chr\":\"").append(escJ(h.chr)).append('"');
        j.append(",\"pos\":").append(h.pos);
        j.append('}');
        return j.toString();
    }

    // GET /api/project/{id}/locus/{n}/novelty?trait=... — known-vs-novel verdict for this locus's
    // genomic region (not just its top SNP): does ANY previously reported association overlap this
    // locus's window at all, and if so, is it for the same trait (via a simple substring match against
    // the optional ?trait= keyword) or a different one. See GwasCatalogLocalIndex#classify.
    private void projectLocusNovelty(HttpExchange ex, String dataDir, String projectDir, String locusIndex) throws IOException {
        File locusFile = new File(dataDir, "locus_" + locusIndex + ".json");
        if (!locusFile.exists()) {
            respond(ex, 404, "application/json", "{\"error\":\"Locus not found\"}".getBytes());
            return;
        }
        String locusJson = new String(Files.readAllBytes(locusFile.toPath()), "UTF-8");
        String chr = extractStr(locusJson, "chr");
        String startStr = extractStr(locusJson, "start");
        String endStr = extractStr(locusJson, "end");
        if (chr == null || startStr == null || endStr == null) {
            respond(ex, 200, "application/json",
                "{\"ok\":false,\"error\":\"Locus JSON missing chr/start/end\"}".getBytes());
            return;
        }

        // ?trait= overrides the project's own configured disease name; neither given means "any trait counts".
        String trait = queryParam(ex, "trait");
        Config cfg = loadConfigQuiet(projectDir);
        if (trait == null || trait.trim().isEmpty())
            trait = (cfg != null && cfg.diseaseName != null && !cfg.diseaseName.trim().isEmpty()) ? cfg.diseaseName.trim() : null;

        // The GWAS Catalog index is GRCh38-only. A GRCh37 project's locus coordinates need converting
        // before comparing against it, or every check would silently compare mismatched coordinate
        // systems — see GenomeLiftover's class comment for how that was verified.
        long queryStart, queryEnd;
        String queryChr = chr;
        boolean lifted = false;
        if (cfg != null && GenomeLiftover.needsLiftover(cfg.genomeBuild)) {
            GenomeLiftover.Result lo = GenomeLiftover.toGRCh38(chr, Long.parseLong(startStr), Long.parseLong(endStr));
            if (!lo.ok) {
                respond(ex, 200, "application/json",
                    ("{\"ok\":false,\"error\":\"GRCh37->GRCh38 liftover failed: " + escJ(lo.error) + "\"}").getBytes());
                return;
            }
            queryChr = lo.chr;
            queryStart = lo.start;
            queryEnd = lo.end;
            lifted = true;
        } else {
            queryStart = Long.parseLong(startStr);
            queryEnd = Long.parseLong(endStr);
        }

        GwasCatalogLocalIndex idx = getGwasCatalogIndex();
        if (idx == null) {
            respond(ex, 503, "application/json",
                ("{\"ok\":false,\"error\":\"" + escJ(gwasCatalogInitError) + "\"}").getBytes());
            return;
        }

        GwasCatalogLocalIndex.RegionResult r;
        try {
            r = idx.classify(queryChr, queryStart, queryEnd, trait);
        } catch (Exception e) {
            respond(ex, 500, "application/json", ("{\"ok\":false,\"error\":\"" + escJ(e.getMessage()) + "\"}").getBytes());
            return;
        }

        StringBuilder json = new StringBuilder();
        json.append("{\"ok\":true,\"verdict\":\"").append(r.verdict.name().toLowerCase(java.util.Locale.ROOT)).append('"');
        json.append(",\"total_count\":").append(r.totalCount);
        json.append(",\"trait_used\":").append(trait == null ? "null" : "\"" + escJ(trait) + "\"");
        json.append(",\"lifted_to_grch38\":").append(lifted);
        if (lifted) json.append(",\"grch38_region\":\"").append(escJ(queryChr)).append(':').append(queryStart).append('-').append(queryEnd).append('"');
        json.append(",\"hits\":[");
        for (int i = 0; i < r.hits.size(); i++) {
            if (i > 0) json.append(',');
            json.append(hitToJson(r.hits.get(i)));
        }
        json.append("]}");
        respond(ex, 200, "application/json", json.toString().getBytes("UTF-8"));
    }

    /** Loads a project's Config, or null if that fails — used where a missing/unreadable config
     *  should just mean "no extra context available" rather than failing the whole request. */
    private static Config loadConfigQuiet(String projectDir) {
        try { return Config.loadFromProject(projectDir); }
        catch (Exception e) { return null; }
    }

    // GET /api/project/{id}/catalog-report?trait=... — genome-wide TSV export: one row per (locus,
    // overlapping Catalog hit), plus a verdict column per locus, for every locus in the project.
    private void projectCatalogReport(HttpExchange ex, String dataDir, String projectDir) throws IOException {
        File manifestFile = new File(dataDir, "manifest.json");
        if (!manifestFile.exists()) {
            respond(ex, 404, "application/json", "{\"error\":\"No manifest — process the project first\"}".getBytes());
            return;
        }
        GwasCatalogLocalIndex idx = getGwasCatalogIndex();
        if (idx == null) {
            respond(ex, 503, "application/json", ("{\"error\":\"" + escJ(gwasCatalogInitError) + "\"}").getBytes());
            return;
        }
        Config cfg = loadConfigQuiet(projectDir);
        String trait = queryParam(ex, "trait");
        if (trait == null || trait.trim().isEmpty())
            trait = (cfg != null && cfg.diseaseName != null && !cfg.diseaseName.trim().isEmpty()) ? cfg.diseaseName.trim() : null;
        boolean needsLift = cfg != null && GenomeLiftover.needsLiftover(cfg.genomeBuild);

        String manifestJson = new String(Files.readAllBytes(manifestFile.toPath()), "UTF-8");
        List<Integer> locusIndices = extractLocusIndices(manifestJson);

        StringBuilder tsv = new StringBuilder();
        tsv.append("locus_index\tchr\tstart\tend\tverdict_region_grch38\ttop_snp\tverdict\ttotal_known_hits\thit_trait\thit_mapped_trait\thit_gene\thit_pvalue\thit_study\thit_pmid\thit_url\n");
        for (int locusIndex : locusIndices) {
            File locusFile = new File(dataDir, "locus_" + locusIndex + ".json");
            if (!locusFile.exists()) continue;
            String locusJson = new String(Files.readAllBytes(locusFile.toPath()), "UTF-8");
            String chr = extractStr(locusJson, "chr");
            String startStr = extractStr(locusJson, "start");
            String endStr = extractStr(locusJson, "end");
            int topSnpIdx = locusJson.indexOf("\"top_snp\":");
            String topSnp = topSnpIdx >= 0 ? extractNestedStr(locusJson, topSnpIdx, "id") : "";
            if (chr == null || startStr == null || endStr == null) continue;

            long qStart, qEnd;
            String qChr = chr, liftedNote = "";
            if (needsLift) {
                GenomeLiftover.Result lo = GenomeLiftover.toGRCh38(chr, Long.parseLong(startStr), Long.parseLong(endStr));
                if (!lo.ok) {
                    tsv.append(locusIndex).append('\t').append(chr).append('\t').append(startStr).append('\t')
                       .append(endStr).append("\tLIFTOVER_FAILED\t").append(nz(topSnp)).append("\terror\t0\t\t\t\t\t\t\t\n");
                    continue;
                }
                qChr = lo.chr; qStart = lo.start; qEnd = lo.end;
                liftedNote = qChr + ":" + qStart + "-" + qEnd;
            } else {
                qStart = Long.parseLong(startStr); qEnd = Long.parseLong(endStr);
            }

            GwasCatalogLocalIndex.RegionResult r;
            try {
                r = idx.classify(qChr, qStart, qEnd, trait);
            } catch (Exception e) { continue; }

            if (r.hits.isEmpty()) {
                tsv.append(locusIndex).append('\t').append(chr).append('\t').append(startStr).append('\t')
                   .append(endStr).append('\t').append(nz(liftedNote)).append('\t').append(nz(topSnp)).append('\t')
                   .append(r.verdict.name().toLowerCase(Locale.ROOT))
                   .append('\t').append(r.totalCount).append("\t\t\t\t\t\t\t\n");
            } else {
                // r.hits is capped at GwasCatalogLocalIndex.MAX_OVERLAP_HITS (a handful of loci — APOE,
                // the MHC region — have tens of thousands of associations; total_known_hits is still the
                // exact uncapped count, so nothing here is silently lying, just not exhaustively listed.
                for (GwasCatalogLocalIndex.Hit h : r.hits) {
                    tsv.append(locusIndex).append('\t').append(chr).append('\t').append(startStr).append('\t')
                       .append(endStr).append('\t').append(nz(liftedNote)).append('\t').append(nz(topSnp)).append('\t')
                       .append(r.verdict.name().toLowerCase(Locale.ROOT)).append('\t')
                       .append(r.totalCount).append('\t')
                       .append(tsvEsc(h.trait)).append('\t').append(tsvEsc(h.mappedTrait)).append('\t')
                       .append(tsvEsc(h.mappedGene)).append('\t').append(Double.isFinite(h.pvalue) ? h.pvalue : "").append('\t')
                       .append(tsvEsc(h.studyAccession)).append('\t').append(tsvEsc(h.pubmedId)).append('\t')
                       .append(tsvEsc(h.link)).append('\n');
                }
            }
        }

        ex.getResponseHeaders().set("Content-Disposition", "attachment; filename=\"catalog_report.tsv\"");
        respond(ex, 200, "text/tab-separated-values", tsv.toString().getBytes("UTF-8"));
    }

    private static String nz(String s) { return s == null ? "" : s; }
    private static String tsvEsc(String s) { return s == null ? "" : s.replace("\t", " ").replace("\n", " "); }

    private static List<Integer> extractLocusIndices(String manifestJson) {
        List<Integer> out = new ArrayList<>();
        java.util.regex.Matcher m = java.util.regex.Pattern.compile("\"locusIndex\"\\s*:\\s*(\\d+)").matcher(manifestJson);
        while (m.find()) out.add(Integer.parseInt(m.group(1)));
        if (out.isEmpty()) {
            m = java.util.regex.Pattern.compile("\"index\"\\s*:\\s*(\\d+)").matcher(manifestJson);
            while (m.find()) out.add(Integer.parseInt(m.group(1)));
        }
        return out;
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
        ProjectState ps = ensureProjectState(projectId, projectDir);
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
        ProjectState ps = ensureProjectState(projectId, projectDir);
        if (ps == null || ps.config == null || ps.gff == null || ps.loci == null || ps.outputs == null) {
            respond(ex, 503, "application/json",
                "{\"error\":\"Pipeline state not available\"}".getBytes()); return;
        }
        LociMutationService ms = ensureMutationService(ps);
        if (ms == null) {
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
        LociMutationService.MutationResult mr = ms.split(
            Integer.parseInt(sOrigIdx), regions);
        if (mr.ok) {
            respond(ex, 200, "application/json", mr.manifestJson.getBytes("UTF-8"));
        } else {
            respond(ex, 500, "application/json",
                ("{\"error\":\"" + escJ(mr.error) + "\"}").getBytes());
        }
    }

    // POST /api/project/{id}/create-locus
    private void projectCreateLocus(HttpExchange ex, String projectId,
                                    String projectDir) throws IOException {
        if (!"POST".equalsIgnoreCase(ex.getRequestMethod())) {
            respond(ex, 405, "text/plain", "Method Not Allowed".getBytes()); return;
        }
        ProjectState ps = ensureProjectState(projectId, projectDir);
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
        ProjectState ps = ensureProjectState(projectId, projectDir);
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

    // POST /api/project/{id}/delete-locus
    private void projectDeleteLocus(HttpExchange ex, String projectId,
                                    String projectDir) throws IOException {
        if (!"POST".equalsIgnoreCase(ex.getRequestMethod())) {
            respond(ex, 405, "text/plain", "Method Not Allowed".getBytes()); return;
        }
        ProjectState ps = ensureProjectState(projectId, projectDir);
        if (ps == null) {
            respond(ex, 503, "application/json",
                "{\"error\":\"Pipeline state not available. Process the project first.\"}".getBytes()); return;
        }
        LociMutationService ms = ensureMutationService(ps);
        if (ms == null) {
            respond(ex, 503, "application/json",
                "{\"error\":\"Pipeline state not available\"}".getBytes()); return;
        }
        String body = new String(readAll(ex.getRequestBody()), "UTF-8");
        String idxStr = extractStr(body, "locus_index");
        if (idxStr == null || idxStr.isEmpty()) {
            respond(ex, 400, "application/json",
                "{\"error\":\"locus_index required\"}".getBytes()); return;
        }
        int locusIndex = Integer.parseInt(idxStr);

        LociMutationService.MutationResult mr = ms.delete(locusIndex);
        if (mr.ok) {
            respond(ex, 200, "application/json", mr.manifestJson.getBytes("UTF-8"));
        } else {
            respond(ex, 500, "application/json",
                ("{\"error\":\"" + escJ(mr.error) + "\"}").getBytes());
        }
    }

    // POST /api/project/{id}/merge-loci
    private void projectMergeLoci(HttpExchange ex, String projectId,
                                  String projectDir) throws IOException {
        if (!"POST".equalsIgnoreCase(ex.getRequestMethod())) {
            respond(ex, 405, "text/plain", "Method Not Allowed".getBytes()); return;
        }
        ProjectState ps = ensureProjectState(projectId, projectDir);
        if (ps == null) {
            respond(ex, 503, "application/json",
                "{\"error\":\"Pipeline state not available. Process the project first.\"}".getBytes()); return;
        }
        LociMutationService ms = ensureMutationService(ps);
        if (ms == null) {
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

        LociMutationService.MutationResult mr = ms.merge(indices, mergedName);
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
        ProjectState ps = ensureProjectState(projectId, projectDir);
        if (ps == null) {
            respond(ex, 503, "application/json",
                "{\"error\":\"Pipeline state not available. Process the project first.\"}".getBytes()); return;
        }
        LociMutationService ms = ensureMutationService(ps);
        if (ms == null) {
            respond(ex, 503, "application/json",
                "{\"error\":\"Pipeline state not available\"}".getBytes()); return;
        }

        LociMutationService.MutationResult mr = ms.undo();
        if (mr.ok) {
            String json = "{\"ok\":true,\"manifest\":" + mr.manifestJson + "}";
            respond(ex, 200, "application/json", json.getBytes("UTF-8"));
        } else {
            respond(ex, 500, "application/json",
                ("{\"error\":\"" + escJ(mr.error) + "\"}").getBytes());
        }
    }

    // POST /api/project/{id}/reorder-loci
    private void projectReorderLoci(HttpExchange ex, String projectId,
                                    String projectDir) throws IOException {
        if (!"POST".equalsIgnoreCase(ex.getRequestMethod())) {
            respond(ex, 405, "text/plain", "Method Not Allowed".getBytes()); return;
        }
        ProjectState ps = ensureProjectState(projectId, projectDir);
        if (ps == null) {
            respond(ex, 503, "application/json",
                "{\"error\":\"Pipeline state not available. Process the project first.\"}".getBytes()); return;
        }
        LociMutationService ms = ensureMutationService(ps);
        if (ms == null) {
            respond(ex, 503, "application/json",
                "{\"error\":\"Pipeline state not available\"}".getBytes()); return;
        }

        LociMutationService.MutationResult mr = ms.reorder();
        if (mr.ok) {
            respond(ex, 200, "application/json", mr.manifestJson.getBytes("UTF-8"));
        } else {
            respond(ex, 500, "application/json",
                ("{\"error\":\"" + escJ(mr.error) + "\"}").getBytes());
        }
    }

    // ── Analysis base pipeline endpoints ──────────────────────────────────

    private final Map<String, ProgressTracker> analysisProgress = new ConcurrentHashMap<>();

    // GET /api/project/{id}/analysis/base-status
    private void projectAnalysisBaseStatus(HttpExchange ex, String projectId,
                                            String projectDir) throws IOException {
        ProjectState ps = ensureProjectState(projectId, projectDir);
        if (ps == null || ps.loci == null) {
            respond(ex, 503, "application/json",
                "{\"error\":\"Pipeline state not available\"}".getBytes()); return;
        }

        StringBuilder json = new StringBuilder("{\"loci\":[");
        boolean first = true;
        for (Locus locus : ps.loci) {
            if (!first) json.append(',');
            first = false;
            Map<String, Boolean> status = BaseStepPipeline.checkStatus(ps.config, locus);
            json.append("{\"id\":\"").append(escJ(locus.id)).append('"');
            json.append(",\"index\":").append(locus.index);
            json.append(",\"base\":").append(status.get("base"));
            json.append(",\"matched\":").append(status.get("matched"));
            json.append(",\"harmonized\":").append(status.get("harmonized"));
            json.append(",\"ld\":").append(status.get("ld"));
            boolean allDone = status.values().stream().allMatch(v -> v);
            json.append(",\"ready\":").append(allDone);
            json.append(",\"mhc_overlap\":").append(BaseStepPipeline.overlapsMhc(locus, ps.config.genomeBuild));
            // Read consistency diagnostic summary if available
            File analysisRoot = BaseStepPipeline.analysisDir(ps.config, locus);
            File diagFile = new File(analysisRoot, "ld/consistency_summary.json");
            if (diagFile.exists()) {
                try {
                    String diagJson = new String(java.nio.file.Files.readAllBytes(diagFile.toPath()), "UTF-8");
                    String verdict = extractStr(diagJson, "verdict");
                    String flagged = extractStr(diagJson, "flagged_snps");
                    if (verdict != null) json.append(",\"diagnostic\":\"").append(escJ(verdict)).append('"');
                    if (flagged != null) json.append(",\"diagnostic_flagged\":").append(flagged);
                } catch (Exception e) {}
            }
            json.append(",\"has_log\":").append(new File(analysisRoot, "build.log").exists());
            json.append('}');
        }
        json.append("]}");
        respond(ex, 200, "application/json", json.toString().getBytes("UTF-8"));
    }

    // GET /api/project/{id}/analysis/locus-log/{locusId}
    private void projectLocusBuildLog(HttpExchange ex, String projectDir, String action)
            throws IOException {
        String locusId = action.substring("analysis/locus-log/".length()).trim();
        if (locusId.isEmpty()) {
            respond(ex, 400, "text/plain", "locus_id required".getBytes()); return;
        }
        File logFile = new File(projectDir, "loci_analysis/" + locusId + "/build.log");
        if (logFile.exists()) {
            byte[] bytes = java.nio.file.Files.readAllBytes(logFile.toPath());
            respond(ex, 200, "text/plain", bytes);
        } else {
            respond(ex, 404, "text/plain", "No build log found for this locus.".getBytes());
        }
    }

    // POST /api/project/{id}/analysis/build-base
    private void projectAnalysisBuildBase(HttpExchange ex, String projectId,
                                           String projectDir) throws IOException {
        if (!"POST".equalsIgnoreCase(ex.getRequestMethod())) {
            respond(ex, 405, "text/plain", "Method Not Allowed".getBytes()); return;
        }
        ProjectState ps = ensureProjectState(projectId, projectDir);
        if (ps == null || ps.config == null || ps.loci == null) {
            respond(ex, 503, "application/json",
                "{\"error\":\"Pipeline state not available. Process the project first.\"}".getBytes()); return;
        }
        if (ps.config.refPanelPath.isEmpty()) {
            respond(ex, 400, "application/json",
                "{\"error\":\"No reference panel configured\"}".getBytes()); return;
        }

        String body = new String(readAll(ex.getRequestBody()), "UTF-8");
        String locusId = extractStr(body, "locus_id");
        String threadsStr = extractStr(body, "threads");
        String ldWindowStr = extractStr(body, "ld_window");
        final int nThreads = (threadsStr != null && !threadsStr.isEmpty())
            ? Math.max(1, Math.min(Integer.parseInt(threadsStr), 16)) : 2;
        final int ldWindow = (ldWindowStr != null && !ldWindowStr.isEmpty())
            ? Math.max(10, Math.min(Integer.parseInt(ldWindowStr), 5000)) : LdMatrixComputer.DEFAULT_LD_WINDOW;

        ProgressTracker pt = new ProgressTracker();
        analysisProgress.put(projectId, pt);

        // Run in background thread
        final Config cfg = ps.config;
        final List<Locus> loci = ps.loci;
        new Thread(() -> {
            try {
                if (locusId != null && !locusId.isEmpty()) {
                    Locus target = null;
                    for (Locus l : loci) if (l.id.equals(locusId)) { target = l; break; }
                    if (target != null) {
                        pt.update("Building base artifacts", 0, 1);
                        BaseStepPipeline.runAll(cfg, target, ldWindow);
                        pt.update("Complete", 1, 1);
                    }
                } else {
                    BaseStepPipeline.runAllLoci(cfg, loci, pt, nThreads, ldWindow);
                }
                pt.done = true;
            } catch (Exception e) {
                pt.phase = "Error: " + e.getMessage();
                pt.done = true;
                System.err.printf("[Analysis] Build base failed: %s%n", e.getMessage());
            }
        }, "analysis-build-" + projectId).start();

        respond(ex, 202, "application/json",
            "{\"status\":\"started\",\"message\":\"Base artifact build started\"}".getBytes());
    }

    // GET /api/project/{id}/analysis/build-base-progress
    private void projectAnalysisBuildBaseProgress(HttpExchange ex, String projectId)
            throws IOException {
        ProgressTracker pt = analysisProgress.get(projectId);
        if (pt == null) {
            respond(ex, 200, "application/json",
                "{\"phase\":\"idle\",\"pct\":0,\"done\":true,\"completed_loci\":[]}".getBytes()); return;
        }
        StringBuilder json = new StringBuilder();
        json.append("{\"phase\":\"").append(escJ(pt.phase)).append('"');
        json.append(",\"locus\":").append(pt.locusIndex);
        json.append(",\"total\":").append(pt.totalLoci);
        json.append(",\"pct\":").append(pt.pct());
        json.append(",\"done\":").append(pt.done);
        json.append(",\"completed_loci\":[");
        boolean first = true;
        for (int idx : pt.completedLoci) {
            if (!first) json.append(',');
            json.append(idx);
            first = false;
        }
        json.append("]}");
        respond(ex, 200, "application/json", json.toString().getBytes());
    }

    // ── Analysis tool/run endpoints ────────────────────────────────

    private final Map<String, PluginEngine.RunResult> analysisJobs = new ConcurrentHashMap<>();
    private final Map<String, ProgressTracker> analysisJobProgress = new ConcurrentHashMap<>();

    private final Map<String, MultiLocusProgress> locusMatrixJobProgress = new ConcurrentHashMap<>();
    private final Map<String, MultiLocusResult>   locusMatrixJobs        = new ConcurrentHashMap<>();
    private final Map<String, LocusMatrixJobMeta> locusMatrixJobMeta     = new ConcurrentHashMap<>();
    private final List<String> locusMatrixJobOrder = new CopyOnWriteArrayList<>();

    private static class LocusMatrixJobMeta {
        String jobId;
        String name;
        List<String> projectIds;
        List<String> datasetNames;
        String refPanelId;
        String refPanelLabel;
        String createdAt;
        String kind = "constellation"; // "constellation" | "serpent" — which view the run was started for
    }
    private final Set<String> cancelledJobs = ConcurrentHashMap.newKeySet();
    /** Serpent Plot JSON per completed job (the category analysis runs permutations, so build once). */
    private final Map<String, String> serpentPlotCache = new ConcurrentHashMap<>();

    // GET /api/project/{id}/analysis/tools — discovered descriptors + param schemas
    private void projectAnalysisTools(HttpExchange ex) throws IOException {
        List<ToolDescriptor> tools = PluginEngine.discoverTools();
        StringBuilder json = new StringBuilder("{\"tools\":[");
        for (int i = 0; i < tools.size(); i++) {
            if (i > 0) json.append(',');
            json.append(tools.get(i).toFormJson());
        }
        json.append("]}");
        respond(ex, 200, "application/json", json.toString().getBytes("UTF-8"));
    }

    // POST /api/project/{id}/analysis/run — {tool, params, locus_id} → job id
    private void projectAnalysisRun(HttpExchange ex, String projectId,
                                     String projectDir) throws IOException {
        if (!"POST".equalsIgnoreCase(ex.getRequestMethod())) {
            respond(ex, 405, "text/plain", "Method Not Allowed".getBytes()); return;
        }
        ProjectState ps = ensureProjectState(projectId, projectDir);
        if (ps == null || ps.config == null || ps.loci == null) {
            respond(ex, 503, "application/json",
                "{\"error\":\"Pipeline state not available\"}".getBytes()); return;
        }

        String body = new String(readAll(ex.getRequestBody()), "UTF-8");
        String toolName = extractStr(body, "tool");
        String locusId = extractStr(body, "locus_id");
        if (toolName == null || locusId == null) {
            respond(ex, 400, "application/json",
                "{\"error\":\"tool and locus_id required\"}".getBytes()); return;
        }

        Locus locus = null;
        for (Locus l : ps.loci) if (l.id.equals(locusId)) { locus = l; break; }
        if (locus == null) {
            respond(ex, 404, "application/json",
                "{\"error\":\"Locus not found\"}".getBytes()); return;
        }

        // Parse user params from body
        Map<String, String> params = new LinkedHashMap<>();
        ToolDescriptor td = PluginEngine.findTool(toolName);
        if (td != null) {
            for (ToolDescriptor.Param p : td.params) {
                String val = extractStr(body, p.name);
                if (val != null) params.put(p.name, val);
                else if (p.defaultValue != null) params.put(p.name, p.defaultValue);
            }
        }

        // Prepare COJO input if this is a COJO tool
        final Locus targetLocus = locus;
        final Config cfg = ps.config;
        String jobId = UUID.randomUUID().toString();

        ProgressTracker pt = new ProgressTracker();
        analysisJobProgress.put(jobId, pt);

        // Run in background
        new Thread(() -> runAnalysisJob(projectId, cfg, targetLocus, toolName, params, jobId, pt),
            "analysis-run-" + jobId).start();

        respond(ex, 202, "application/json",
            ("{\"job_id\":\"" + jobId + "\",\"status\":\"started\"}").getBytes());
    }

    /**
     * Runs one tool on one locus and blocks until it finishes; the result is also recorded under
     * {@code jobId} for the job endpoints. Shared by the per-locus "Run" button and batch runs.
     */
    private PluginEngine.RunResult runAnalysisJob(String projectId, Config cfg, Locus targetLocus, String toolName,
                                                  Map<String, String> params, String jobId, ProgressTracker pt) {
        String locusId = targetLocus.id;
        try {
            pt.update("Preparing", 0, 4);

            // Tool-specific input preparation
            File analysisRoot = BaseStepPipeline.analysisDir(cfg, targetLocus);
            File harmonizedDir = new File(analysisRoot, "harmonized");
            File matchedDir = new File(analysisRoot, "matched");
            File ldDir = new File(analysisRoot, "ld");
            File runDir = new File(analysisRoot, "runs/" + jobId);
            runDir.mkdirs();

            int sampleN = cfg.sampleN;

            if (toolName.startsWith("cojo")) {
                double pCutoff = 5e-8;
                double collinear = 0.9;
                try { pCutoff = Double.parseDouble(params.getOrDefault("p_cutoff", "5e-8")); } catch (NumberFormatException e) {}
                try { collinear = Double.parseDouble(params.getOrDefault("collinear", "0.9")); } catch (NumberFormatException e) {}
                String gctaBin = params.getOrDefault("gcta_path", "bin/gcta64.exe");
                CojoAdapter.prepareRun(harmonizedDir, matchedDir, runDir, sampleN, pCutoff, collinear, gctaBin);
            } else if (toolName.equals("susie_finemapping")) {
                int maxCausal = 10;
                double coverage = 0.95, ldShrink = 0.1;
                int windowKb = 250;
                try { maxCausal = Integer.parseInt(params.getOrDefault("max_causal", "10")); } catch (NumberFormatException e) {}
                try { coverage = Double.parseDouble(params.getOrDefault("coverage", "0.95")); } catch (NumberFormatException e) {}
                try { ldShrink = Double.parseDouble(params.getOrDefault("ld_shrink", "0.1")); } catch (NumberFormatException e) {}
                try { windowKb = Integer.parseInt(params.getOrDefault("window_kb", "250")); } catch (NumberFormatException e) {}
                SusieAdapter.prepareRun(harmonizedDir, matchedDir, runDir, sampleN, maxCausal, coverage, ldShrink, windowKb);
            } else if (toolName.equals("finemap")) {
                int maxCausal = 5;
                try { maxCausal = Integer.parseInt(params.getOrDefault("max_causal", "5")); } catch (NumberFormatException e) {}
                FinemapAdapter.prepareRun(harmonizedDir, ldDir, matchedDir, runDir, sampleN, maxCausal);
            } else if (toolName.equals("coloc")) {
                String trait2File = params.get("trait2_file");
                String trait2Type = params.getOrDefault("trait2_type", "quant");
                int trait2N = 0, trait2NCases = 0;
                double p1 = 1e-4, p2 = 1e-4, p12 = 1e-5, pipThreshold = 0.1;
                boolean restrictToFinemapped = "true".equalsIgnoreCase(params.getOrDefault("restrict_to_finemapped", "false"));
                try { trait2N = Integer.parseInt(params.getOrDefault("trait2_n", "0")); } catch (NumberFormatException e) {}
                try { trait2NCases = Integer.parseInt(params.getOrDefault("trait2_n_cases", "0")); } catch (NumberFormatException e) {}
                try { p1 = Double.parseDouble(params.getOrDefault("p1", "1e-4")); } catch (NumberFormatException e) {}
                try { p2 = Double.parseDouble(params.getOrDefault("p2", "1e-4")); } catch (NumberFormatException e) {}
                try { p12 = Double.parseDouble(params.getOrDefault("p12", "1e-5")); } catch (NumberFormatException e) {}
                try { pipThreshold = Double.parseDouble(params.getOrDefault("finemap_pip_threshold", "0.1")); } catch (NumberFormatException e) {}
                ColocAdapter.prepareRun(harmonizedDir, runDir, targetLocus, cfg,
                    trait2File, trait2Type, trait2N, trait2NCases, p1, p2, p12,
                    analysisRoot, restrictToFinemapped, pipThreshold);
            } else if (toolName.equals("gwama_meta")) {
                GwamaAdapter.prepareRun(harmonizedDir, runDir);
            } else if (toolName.equals("magma_gene")) {
                PluginEngine.RunResult result = runMagmaGene(cfg, targetLocus, params, jobId, runDir, harmonizedDir, matchedDir, pt);
                analysisJobs.put(jobId, result);
                pt.done = true;
                pt.phase = result.ok ? "Complete" : "Error: " + result.error;
                return result;
            }

            PluginEngine.RunRequest req = new PluginEngine.RunRequest();
            req.tool = toolName;
            req.params = params;
            req.locusId = locusId;
            req.projectId = projectId;
            req.projectDir = cfg.outputDir;
            req.preCreatedRunDir = runDir.getAbsolutePath();

            pt.update("Running " + toolName, 1, 4);
            PluginEngine.RunResult result = PluginEngine.execute(req, cfg, targetLocus, pt);
            result.jobId = jobId;

            // R scripts write result.tsv + result.manifest.json directly

            analysisJobs.put(jobId, result);
            pt.done = true;
            pt.phase = result.ok ? "Complete" : "Error: " + result.error;
            return result;

        } catch (Exception e) {
            PluginEngine.RunResult errResult = new PluginEngine.RunResult();
            errResult.jobId = jobId;
            errResult.error = e.getMessage();
            analysisJobs.put(jobId, errResult);
            pt.done = true;
            pt.phase = "Error: " + e.getMessage();
            return errResult;
        }
    }

    /**
     * MAGMA gene-based test on one locus. MAGMA is a native binary driven by MagmaAdapter in Java
     * (annotate, then the gene test against the locus's matched reference panel), so it bypasses the
     * R-script plugin path but leaves the same run folder behind: result.tsv (one row per gene, best
     * first) and provenance.json, which is what saved-result listing and batch runs look for.
     */
    private PluginEngine.RunResult runMagmaGene(Config cfg, Locus locus, Map<String, String> params, String jobId,
                                                File runDir, File harmonizedDir, File matchedDir,
                                                ProgressTracker pt) throws Exception {
        PluginEngine.RunResult res = new PluginEngine.RunResult();
        res.jobId = jobId; res.tool = "magma_gene"; res.runDir = runDir.getAbsolutePath();
        long windowKb = 10;
        try { windowKb = Long.parseLong(params.getOrDefault("window_kb", "10")); } catch (NumberFormatException ignored) {}
        if (cfg.gff3File == null || !new File(cfg.gff3File).isFile())
            throw new IOException("No gene annotation (GFF3) configured for this project");
        String magmaBin = MagmaAdapter.findMagmaBinary(SusieAdapter.appRoot());
        pt.update("Loading gene annotation", 1, 4);
        GffParser gff = gffFuture(new File(cfg.gff3File).getAbsolutePath()).get();

        pt.update("MAGMA: assigning SNPs to genes", 2, 4);
        long start = Math.max(0, locus.start - windowKb * 1000), end = locus.end + windowKb * 1000;
        MagmaAdapter.AnnotateResult ann = MagmaAdapter.annotate(matchedDir, gff, runDir, locus.chr, start, end, windowKb, magmaBin);
        if (!ann.ok) throw new IOException(ann.error);

        pt.update("MAGMA: gene-based test", 3, 4);
        MagmaAdapter.GeneAnalysisResult ga = MagmaAdapter.geneAnalysis(harmonizedDir, matchedDir, ann, runDir, cfg.sampleN, magmaBin);
        if (!ga.ok) {
            if (ga.logTail != null && ga.logTail.contains("analysis failed for all genes"))
                throw new IOException("No gene at this locus has SNPs left after MAGMA's QC within " + windowKb
                    + " kb (" + ann.nGenes + " gene(s) annotated); nothing to test");
            throw new IOException(ga.error);
        }

        Map<String, String> nameOf = new HashMap<>();
        for (Gene g : gff.overlapping(locus.chr.replaceFirst("^chr", ""), start, end)) nameOf.put(g.geneId, g.geneName);
        List<MagmaAdapter.GeneResult> genes = new ArrayList<>(ga.genes);
        genes.sort(Comparator.comparingDouble(g -> g.p));
        int nTested = genes.size();
        try (PrintWriter pw = new PrintWriter(new OutputStreamWriter(new FileOutputStream(new File(runDir, "result.tsv")), "UTF-8"))) {
            pw.println("gene_id\tgene\tchr\tstart\tstop\tnsnps\tmagma_z\tmagma_p\tsignificant");
            for (MagmaAdapter.GeneResult g : genes) {
                String name = nameOf.getOrDefault(g.gene, g.gene);
                // Bonferroni across the genes tested at this locus
                boolean sig = g.p * nTested < 0.05;
                pw.printf(Locale.ROOT, "%s\t%s\t%s\t%d\t%d\t%d\t%.4f\t%.4g\t%s%n",
                    g.gene, name == null || name.isEmpty() ? g.gene : name, g.chr, g.start, g.stop, g.nsnps, g.zstat, g.p, sig ? "yes" : "no");
            }
        }
        try (PrintWriter pw = new PrintWriter(new OutputStreamWriter(new FileOutputStream(new File(runDir, "provenance.json")), "UTF-8"))) {
            pw.printf(Locale.ROOT, "{%n  \"job_id\": \"%s\",%n  \"tool\": \"magma_gene\",%n  \"tool_version\": \"1.10\",%n"
                + "  \"timestamp\": %d,%n  \"params\": {\"window_kb\": \"%d\", \"sample_n\": \"%d\"},%n"
                + "  \"snps_written\": %d,%n  \"snps_dropped_no_ref_id\": %d,%n  \"genes_tested\": %d%n}%n",
                escJ(jobId), System.currentTimeMillis(), windowKb, cfg.sampleN, ga.snpsWritten, ga.snpsDroppedNoRefId, nTested);
        }
        res.ok = true;
        res.resultRows = nTested;
        return res;
    }

    // ── Saved results: read finished runs back from disk ─────────────────────

    /** One finished run of one tool at one locus, summarised for display. */
    static class RunSummary {
        String tool, runId;
        long timestamp;
        int rows;
        String topSnp = "", topPos = "";
        double topPip = Double.NaN;
        int nSets;                                   // SuSiE/FINEMAP credible sets, COJO independent signals
        String reliability = "";                     // COJO: ok | suspect | likely_artifact
        String note = "";                            // COJO: why it is not reliable
        final List<String> setLines = new ArrayList<>();   // "CS1: 9 SNPs, lead rs123 (PIP 0.26)" / "rs123 (pJ 1e-12)"

        String toJson() {
            StringBuilder j = new StringBuilder("{");
            j.append("\"tool\":\"").append(escJ(tool)).append('"');
            j.append(",\"run_id\":\"").append(escJ(runId)).append('"');
            j.append(",\"timestamp\":").append(timestamp);
            j.append(",\"rows\":").append(rows);
            j.append(",\"top_snp\":\"").append(escJ(topSnp)).append('"');
            j.append(",\"top_pos\":\"").append(escJ(topPos)).append('"');
            j.append(",\"top_pip\":").append(Double.isNaN(topPip) ? "null" : String.format(Locale.ROOT, "%.4f", topPip));
            j.append(",\"n_sets\":").append(nSets);
            j.append(",\"reliability\":\"").append(escJ(reliability)).append('"');
            j.append(",\"note\":\"").append(escJ(note)).append('"');
            j.append(",\"sets\":[");
            for (int i = 0; i < setLines.size(); i++) { if (i > 0) j.append(','); j.append('"').append(escJ(setLines.get(i))).append('"'); }
            return j.append("]}").toString();
        }
    }

    /** Newest successful run per tool at this locus (a run counts when provenance.json and result.tsv exist). */
    static Map<String, RunSummary> latestRuns(Config cfg, Locus locus) {
        Map<String, RunSummary> latest = new TreeMap<>();
        File[] dirs = new File(BaseStepPipeline.analysisDir(cfg, locus), "runs").listFiles(File::isDirectory);
        if (dirs == null) return latest;
        for (File d : dirs) {
            File prov = new File(d, "provenance.json"), res = new File(d, "result.tsv");
            if (!prov.isFile() || !res.isFile()) continue;
            try {
                Map<String, Object> p = MiniJson.asObject(MiniJson.parse(new String(Files.readAllBytes(prov.toPath()), "UTF-8")));
                String tool = MiniJson.getStr(p, "tool", "");
                long ts = p.get("timestamp") instanceof Number ? ((Number) p.get("timestamp")).longValue() : d.lastModified();
                RunSummary prev = latest.get(tool);
                if (prev != null && prev.timestamp >= ts) continue;
                RunSummary rs = summarise(tool, d.getName(), ts, res);
                if (rs != null) latest.put(tool, rs);
            } catch (Exception ignored) {}
        }
        return latest;
    }

    static RunSummary summarise(String tool, String runId, long ts, File resultTsv) throws IOException {
        RunSummary rs = new RunSummary();
        rs.tool = tool; rs.runId = runId; rs.timestamp = ts;
        List<String> lines = Files.readAllLines(resultTsv.toPath());
        if (lines.isEmpty()) return rs;
        List<String> h = Arrays.asList(lines.get(0).split("\t", -1));
        if (tool.equals("magma_gene")) {
            // One row per gene, best first: lead = top gene; "sets" = genes significant after
            // Bonferroni across the genes tested at this locus, listed (up to 5) with their p-values.
            int iGene = h.indexOf("gene"), iP = h.indexOf("magma_p"), iN = h.indexOf("nsnps"), iSig = h.indexOf("significant");
            for (int r = 1; r < lines.size(); r++) {
                String[] f = lines.get(r).split("\t", -1);
                if (f.length < h.size()) continue;
                rs.rows++;
                if (rs.topSnp.isEmpty()) rs.topSnp = f[iGene];
                boolean sig = "yes".equals(f[iSig]);
                if (sig) rs.nSets++;
                if (rs.setLines.size() < 5)
                    rs.setLines.add(f[iGene] + ": p " + f[iP] + ", " + f[iN] + " SNPs" + (sig ? " (significant)" : ""));
            }
            if (rs.rows > 0) rs.note = rs.rows + " genes tested; Bonferroni threshold " + String.format(Locale.ROOT, "%.2g", 0.05 / rs.rows);
            return rs;
        }
        int iSnp = h.indexOf("snp_id"), iPos = h.indexOf("pos");
        int iPip = h.indexOf(tool.startsWith("susie") ? "susie_pip" : "finemap_pip");
        int iCs = h.indexOf(tool.startsWith("susie") ? "susie_cs" : "finemap_cs");
        int iSel = h.indexOf("cojo_selected"), iPj = h.indexOf("cojo_pJ");
        Map<String, List<String[]>> sets = new TreeMap<>();
        for (int r = 1; r < lines.size(); r++) {
            String[] f = lines.get(r).split("\t", -1);
            if (f.length < h.size()) continue;
            rs.rows++;
            if (iPip >= 0) {
                try {
                    double pip = Double.parseDouble(f[iPip]);
                    if (Double.isNaN(rs.topPip) || pip > rs.topPip) { rs.topPip = pip; rs.topSnp = f[iSnp]; rs.topPos = iPos >= 0 ? f[iPos] : ""; }
                } catch (NumberFormatException ignored) {}
            }
            if (iCs >= 0 && !f[iCs].isEmpty() && !f[iCs].equals("NA")) sets.computeIfAbsent(f[iCs], k -> new ArrayList<>()).add(f);
            if (iSel >= 0 && f[iSel].equalsIgnoreCase("selected")) {
                rs.nSets++;
                rs.setLines.add(f[iSnp] + (iPj >= 0 && !f[iPj].equals("NA") ? " (joint p " + f[iPj] + ")" : ""));
                if (rs.topSnp.isEmpty()) { rs.topSnp = f[iSnp]; rs.topPos = iPos >= 0 ? f[iPos] : ""; }
            }
        }
        File log = new File(resultTsv.getParentFile(), "run.log");
        if (tool.startsWith("cojo") && log.isFile()) {
            for (String line : Files.readAllLines(log.toPath())) {
                if (line.startsWith("Reliability: ")) rs.reliability = line.substring(13).trim();
                else if (line.trim().startsWith("Reasons: ")) rs.note = line.trim().substring(9).trim();
            }
        }
        for (Map.Entry<String, List<String[]>> e : sets.entrySet()) {
            String[] lead = null; double best = -1;
            for (String[] f : e.getValue()) {
                try { double v = Double.parseDouble(f[iPip]); if (v > best) { best = v; lead = f; } } catch (NumberFormatException ignored) {}
            }
            rs.nSets++;
            rs.setLines.add("Set " + e.getKey() + ": " + e.getValue().size() + " SNPs, lead " + (lead == null ? "?" : lead[iSnp])
                + String.format(Locale.ROOT, " (PIP %.2f)", Math.max(best, 0)));
        }
        if (tool.startsWith("susie")) {
            // Every allowed effect forming a credible set usually means SuSiE is explaining an LD
            // mismatch between the reference panel and the GWAS, not that many causal variants exist.
            int maxL = 10;
            try {
                String prov = new String(Files.readAllBytes(new File(resultTsv.getParentFile(), "provenance.json").toPath()), "UTF-8");
                java.util.regex.Matcher m = java.util.regex.Pattern.compile("\"max_causal\"\\s*:\\s*\"?(\\d+)").matcher(prov);
                if (m.find()) maxL = Integer.parseInt(m.group(1));
            } catch (Exception ignored) {}
            // The run log records the measured GWAS-vs-reference LD inconsistency (estimate_s_rss)
            double sRss = Double.NaN;
            File runLog = new File(resultTsv.getParentFile(), "run.log");
            if (runLog.isFile())
                for (String line : Files.readAllLines(runLog.toPath()))
                    if (line.contains("inconsistency s = "))
                        try { sRss = Double.parseDouble(line.substring(line.indexOf("s = ") + 4).trim()); } catch (NumberFormatException ignored) {}
            boolean saturated = rs.nSets >= maxL, inconsistent = sRss > 0.3;
            rs.reliability = saturated || inconsistent ? "check" : "ok";
            if (saturated)
                rs.note = "All " + maxL + " SuSiE effects formed credible sets: this usually means the reference-panel LD "
                    + "does not match the GWAS in this region (common in large, dense regions). Interpret with caution.";
            else if (inconsistent)
                rs.note = String.format(Locale.ROOT, "The GWAS z-scores are inconsistent with the reference LD here (s = %.2f > 0.3): "
                    + "credible sets may include spurious signals. Interpret with caution.", sRss);
        }
        return rs;
    }

    private Locus findLocus(ProjectState ps, String locusId) {
        if (ps == null || ps.loci == null || locusId == null) return null;
        for (Locus l : ps.loci) if (l.id.equals(locusId)) return l;
        return null;
    }

    // GET /api/project/{id}/analysis/latest?locus_id=... — newest finished run per tool, summarised
    private void projectAnalysisLatest(HttpExchange ex, String projectId, String projectDir) throws IOException {
        ProjectState ps = ensureProjectState(projectId, projectDir);
        Locus l = findLocus(ps, queryParam(ex, "locus_id"));
        if (l == null) { respond(ex, 404, "application/json", "{\"error\":\"Locus not found\"}".getBytes()); return; }
        StringBuilder j = new StringBuilder("{\"runs\":[");
        boolean first = true;
        for (RunSummary rs : latestRuns(ps.config, l).values()) {
            if (!first) j.append(',');
            first = false;
            j.append(rs.toJson());
        }
        respond(ex, 200, "application/json", j.append("]}").toString().getBytes("UTF-8"));
    }

    // GET /api/project/{id}/analysis/run-result?locus_id=...&run_id=... — a saved run's result.tsv
    private void projectAnalysisRunResult(HttpExchange ex, String projectId, String projectDir) throws IOException {
        ProjectState ps = ensureProjectState(projectId, projectDir);
        Locus l = findLocus(ps, queryParam(ex, "locus_id"));
        String runId = queryParam(ex, "run_id");
        if (l == null || runId == null || !runId.matches("[A-Za-z0-9-]{1,64}")) {
            respond(ex, 404, "application/json", "{\"error\":\"Run not found\"}".getBytes()); return;
        }
        File f = new File(new File(new File(BaseStepPipeline.analysisDir(ps.config, l), "runs"), runId), "result.tsv");
        if (!f.isFile()) { respond(ex, 404, "application/json", "{\"error\":\"Run not found\"}".getBytes()); return; }
        respond(ex, 200, "text/tab-separated-values", Files.readAllBytes(f.toPath()));
    }

    // GET /api/project/{id}/analysis/summary[?format=tsv] — every locus x tool, latest run each
    private void projectAnalysisSummary(HttpExchange ex, String projectId, String projectDir) throws IOException {
        ProjectState ps = ensureProjectState(projectId, projectDir);
        if (ps == null || ps.loci == null) { respond(ex, 503, "application/json", "{\"error\":\"Pipeline state not available\"}".getBytes()); return; }
        boolean tsv = "tsv".equals(queryParam(ex, "format"));
        StringBuilder out = new StringBuilder();
        if (tsv) out.append("locus\tlocus_id\tchr\tstart\tend\ttool\ttop_snp\ttop_pos\ttop_pip\tn_sets\tsets\treliability\tnote\trun_id\n");
        else out.append("{\"loci\":[");
        boolean first = true;
        for (Locus l : ps.loci) {
            Map<String, RunSummary> runs = latestRuns(ps.config, l);
            if (tsv) {
                for (RunSummary rs : runs.values())
                    out.append(l.index).append('\t').append(l.id).append('\t').append(l.chr).append('\t').append(l.start).append('\t').append(l.end)
                       .append('\t').append(rs.tool).append('\t').append(rs.topSnp).append('\t').append(rs.topPos)
                       .append('\t').append(Double.isNaN(rs.topPip) ? "" : String.format(Locale.ROOT, "%.4f", rs.topPip))
                       .append('\t').append(rs.nSets).append('\t').append(String.join("; ", rs.setLines))
                       .append('\t').append(rs.reliability).append('\t').append(rs.note.replace('\t', ' ')).append('\t').append(rs.runId).append('\n');
                continue;
            }
            if (!first) out.append(',');
            first = false;
            out.append("{\"index\":").append(l.index).append(",\"id\":\"").append(escJ(l.id)).append("\",\"chr\":\"").append(escJ(l.chr))
               .append("\",\"start\":").append(l.start).append(",\"end\":").append(l.end)
               .append(",\"mhc\":").append(BaseStepPipeline.overlapsMhc(l, ps.config.genomeBuild)).append(",\"runs\":[");
            boolean f2 = true;
            for (RunSummary rs : runs.values()) { if (!f2) out.append(','); f2 = false; out.append(rs.toJson()); }
            out.append("]}");
        }
        if (tsv) {
            ex.getResponseHeaders().set("Content-Disposition", "attachment; filename=\"" + projectId + "_finemapping_summary.tsv\"");
            respond(ex, 200, "text/tab-separated-values", out.toString().getBytes("UTF-8"));
        } else {
            respond(ex, 200, "application/json", out.append("]}").toString().getBytes("UTF-8"));
        }
    }

    // ── Joint fine-mapping across datasets (mvSuSiE) on an aligned locus ────
    //
    // POST /api/mvsusie-run {"job": <cross-dataset run id>, "locus": <aligned locus index>,
    //                        "mode": "within" | "cross", "disease": <optional, within mode>}
    // Safeguards: only datasets genome-wide significant at the locus with built inputs; datasets whose
    // SuSiE result is flagged "check" (LD mismatch) are left out; "within" (default) keeps one disease;
    // "cross" keeps only datasets whose SuSiE credible sets overlap the consensus signal, since pooling
    // different causal variants misleads the fit. Every exclusion is reported with its reason; the R side
    // reports the between-dataset correlation (sample overlap) and, per credible set, the datasets it is
    // active in.
    //
    // Every run is saved under output/mvsusie/<id>/ (manifest, result.json, meta.json) and indexed in
    // memory, so a locus shows its latest joint result when reopened, each dataset's regional plot can
    // draw it, and a whole cross-dataset run can be fine-mapped jointly in one batch.

    static class MvsusieJob {
        String ownerId;
        volatile String status = "running", phase = "Selecting datasets", error;
        volatile String selectionJson = "{}", resultJson, metaJson;
        final long startedAt = System.currentTimeMillis();
    }
    private final Map<String, MvsusieJob> mvsusieJobs = new ConcurrentHashMap<>();
    static final long MVSUSIE_HALF_WINDOW = 250_000L;
    static final File MVSUSIE_DIR = new File("output/mvsusie");

    /** A saved run's summary (from its meta.json), for listing without reading the result. */
    static class MvsusieSaved {
        String id, job, mode, disease = "", status, error, ownerId, leadSnp = "", gene = "", chr = "";
        int locus, nSets, nDatasets;
        long created;
        double leadPip = Double.NaN;
        final Map<String, Integer> datasetLoci = new LinkedHashMap<>();   // dataset -> its own locus index

        String toJson() {
            StringBuilder j = new StringBuilder("{");
            j.append("\"id\":\"").append(escJ(id)).append("\",\"job\":\"").append(escJ(job)).append("\",\"locus\":").append(locus)
             .append(",\"mode\":\"").append(escJ(mode)).append("\",\"disease\":\"").append(escJ(disease))
             .append("\",\"status\":\"").append(escJ(status)).append("\",\"created\":").append(created)
             .append(",\"n_sets\":").append(nSets).append(",\"n_datasets\":").append(nDatasets)
             .append(",\"lead_snp\":\"").append(escJ(leadSnp)).append('"')
             .append(",\"lead_pip\":").append(Double.isNaN(leadPip) ? "null" : String.format(Locale.ROOT, "%.4f", leadPip))
             .append(",\"gene\":\"").append(escJ(gene)).append('"');
            if (error != null) j.append(",\"error\":\"").append(escJ(error)).append('"');
            j.append(",\"datasets\":{");
            int i = 0;
            for (Map.Entry<String, Integer> e : datasetLoci.entrySet())
                j.append(i++ > 0 ? "," : "").append('"').append(escJ(e.getKey())).append("\":").append(e.getValue());
            return j.append("}}").toString();
        }
    }
    private final Map<String, MvsusieSaved> mvsusieSaved = new ConcurrentHashMap<>();
    private volatile boolean mvsusieIndexed = false;

    // Access hooks (desktop: single user, everything visible; the server edition narrows these)
    private String mvsusieUser(HttpExchange ex) { return null; }
    private boolean mvsusieCanSeeRun(HttpExchange ex, String job) { return job != null; }
    private boolean mvsusieReadable(HttpExchange ex, String pid) { return true; }
    private boolean mvsusieVisible(HttpExchange ex, MvsusieSaved s) {
        return s.ownerId == null || s.ownerId.equals(mvsusieUser(ex));
    }

    private static MvsusieSaved parseSavedMeta(Map<String, Object> m) {
        MvsusieSaved s = new MvsusieSaved();
        s.id = MiniJson.getStr(m, "id", "");
        s.job = MiniJson.getStr(m, "job", "");
        s.locus = (int) num(m.get("locus"), -1);
        s.mode = MiniJson.getStr(m, "mode", "within");
        s.disease = MiniJson.getStr(m, "disease", "");
        s.status = MiniJson.getStr(m, "status", "error");
        s.error = m.get("error") == null ? null : String.valueOf(m.get("error"));
        s.ownerId = m.get("owner") == null ? null : String.valueOf(m.get("owner"));
        s.created = (long) num(m.get("created"), 0);
        s.nSets = (int) num(m.get("n_sets"), 0);
        s.leadSnp = MiniJson.getStr(m, "lead_snp", "");
        s.leadPip = num(m.get("lead_pip"), Double.NaN);
        s.gene = MiniJson.getStr(m, "gene", "");
        s.chr = MiniJson.getStr(m, "chr", "");
        Object ds = m.get("datasets");
        if (ds instanceof List)
            for (Object o : (List<?>) ds)
                if (o instanceof Map) {
                    Map<String, Object> d = MiniJson.asObject(o);
                    s.datasetLoci.put(MiniJson.getStr(d, "id", ""), (int) num(d.get("locus_index"), -1));
                }
        s.nDatasets = s.datasetLoci.size();
        return s;
    }

    private static double num(Object o, double dflt) {
        if (o instanceof Number) return ((Number) o).doubleValue();
        try { return o == null ? dflt : Double.parseDouble(String.valueOf(o)); } catch (NumberFormatException e) { return dflt; }
    }

    /** Loads every saved run's meta.json once (lazily, on first use). */
    private void ensureMvsusieIndex() {
        if (mvsusieIndexed) return;
        synchronized (mvsusieSaved) {
            if (mvsusieIndexed) return;
            File[] dirs = MVSUSIE_DIR.listFiles(File::isDirectory);
            if (dirs != null)
                for (File d : dirs) {
                    File meta = new File(d, "meta.json");
                    if (!meta.isFile()) continue;
                    try {
                        MvsusieSaved s = parseSavedMeta(MiniJson.asObject(MiniJson.parse(
                            new String(Files.readAllBytes(meta.toPath()), "UTF-8"))));
                        if (s.id.equals(d.getName())) mvsusieSaved.put(s.id, s);
                    } catch (Exception e) {
                        System.err.println("[mvSuSiE] skipping unreadable " + meta + ": " + e.getMessage());
                    }
                }
            mvsusieIndexed = true;
        }
    }

    /** Saved runs at one aligned locus of one cross-dataset run, newest first. */
    private List<MvsusieSaved> savedAtLocus(String job, int locus) {
        ensureMvsusieIndex();
        List<MvsusieSaved> out = new ArrayList<>();
        for (MvsusieSaved s : mvsusieSaved.values()) if (s.job.equals(job) && s.locus == locus) out.add(s);
        out.sort((a, b) -> Long.compare(b.created, a.created));
        return out;
    }

    private String diseaseOf(MultiLocusResult result, String pid) {
        // Same grouping the Serpent Plot shows, so "within one disease" matches what the user sees
        for (MultiLocusResult.DatasetInfo di : result.datasets)
            if (di.id.equals(pid)) return GeneConstellationBuilder.diseaseGroupOf(pid, di.diseaseName);
        return GeneConstellationBuilder.diseaseGroupOf(pid);
    }

    /**
     * The dataset's latest SuSiE credible sets at this locus as chr:pos sets, the primary set first:
     * the set containing the dataset's strongest GWAS association ({@code leadKey}), else the set
     * holding the highest-PIP SNP. Empty when there is no run or no set.
     */
    private static List<Set<String>> susieCredibleSets(Config cfg, Locus l, RunSummary rs, String leadKey) {
        Map<String, Set<String>> sets = new LinkedHashMap<>();
        String primary = null; double best = -1;
        if (rs == null) return new ArrayList<>();
        File f = new File(new File(new File(BaseStepPipeline.analysisDir(cfg, l), "runs"), rs.runId), "result.tsv");
        try {
            List<String> lines = Files.readAllLines(f.toPath());
            List<String> h = Arrays.asList(lines.get(0).split("\t", -1));
            int ic = h.indexOf("chr"), ip = h.indexOf("pos"), ics = h.indexOf("susie_cs"), ipip = h.indexOf("susie_pip");
            for (int i = 1; i < lines.size(); i++) {
                String[] x = lines.get(i).split("\t", -1);
                if (ics < 0 || ics >= x.length || x[ics].isEmpty() || x[ics].equals("NA")) continue;
                sets.computeIfAbsent(x[ics], k -> new HashSet<>()).add(x[ic] + ":" + x[ip]);
                try { double v = Double.parseDouble(x[ipip]); if (v > best) { best = v; primary = x[ics]; } } catch (Exception ignored) {}
            }
        } catch (Exception ignored) {}
        for (Map.Entry<String, Set<String>> e : sets.entrySet())
            if (e.getValue().contains(leadKey)) { primary = e.getKey(); break; }   // strongest association wins
        List<Set<String>> out = new ArrayList<>();
        if (primary != null) out.add(sets.remove(primary));
        out.addAll(sets.values());
        return out;
    }

    /** Single-dataset SuSiE PIPs (position -> PIP >= 0.001) from the dataset's latest run at this locus. */
    private static Map<String, Object> susiePips(Config cfg, Locus l, RunSummary rs, long from, long to) {
        Map<String, Object> out = new LinkedHashMap<>();
        if (rs == null) return out;
        File f = new File(new File(new File(BaseStepPipeline.analysisDir(cfg, l), "runs"), rs.runId), "result.tsv");
        try {
            List<String> lines = Files.readAllLines(f.toPath());
            List<String> h = Arrays.asList(lines.get(0).split("\t", -1));
            int ip = h.indexOf("pos"), ipip = h.indexOf("susie_pip");
            if (ip < 0 || ipip < 0) return out;
            for (int i = 1; i < lines.size(); i++) {
                String[] x = lines.get(i).split("\t", -1);
                try {
                    long pos = Long.parseLong(x[ip]);
                    double v = Double.parseDouble(x[ipip]);
                    if (v >= 0.001 && pos >= from && pos <= to) out.put(String.valueOf(pos), Math.round(v * 10000) / 10000.0);
                } catch (Exception ignored) {}
            }
        } catch (Exception ignored) {}
        return out;
    }

    /** Datasets genome-wide significant at an aligned locus, by disease group (the Serpent Plot's grouping). */
    private Map<String, List<String>> significantByDisease(MultiLocusResult result, MultiLocusResult.LocusRow row) {
        Map<String, List<String>> out = new TreeMap<>();
        for (Map.Entry<String, MultiLocusResult.DatasetLocusStat> e : row.cells.entrySet())
            if (e.getValue() != null && e.getValue().nP5e8 > 0)
                out.computeIfAbsent(diseaseOf(result, e.getKey()), k -> new ArrayList<>()).add(e.getKey());
        return out;
    }

    private MultiLocusResult.LocusRow locusRow(MultiLocusResult result, String locusStr) {
        try {
            int li = Integer.parseInt(locusStr);
            for (MultiLocusResult.LocusRow r : result.loci) if (r.index == li) return r;
        } catch (Exception ignored) {}
        return null;
    }

    private void mvsusieRun(HttpExchange ex) throws IOException {
        cors(ex); if (preflight(ex)) return;
        if (!"POST".equalsIgnoreCase(ex.getRequestMethod())) {
            respond(ex, 405, "application/json", "{\"error\":\"POST required\"}".getBytes()); return;
        }
        String body = new String(readAll(ex.getRequestBody()), "UTF-8");
        String jobRef = extractStr(body, "job");
        String mode = "cross".equals(extractStr(body, "mode")) ? "cross" : "within";
        String wantDisease = extractStr(body, "disease");
        if (wantDisease != null && wantDisease.length() > 80) wantDisease = null;
        MultiLocusResult result = jobRef == null || !mvsusieCanSeeRun(ex, jobRef) ? null : locusMatrixJobs.get(jobRef);
        if (result == null) { respond(ex, 404, "application/json", "{\"error\":\"Cross-dataset run not found\"}".getBytes()); return; }
        MultiLocusResult.LocusRow row = locusRow(result, extractStr(body, "locus"));
        if (row == null) { respond(ex, 404, "application/json", "{\"error\":\"Locus not found\"}".getBytes()); return; }
        if (mvsusieBatches.containsKey(jobRef) && !mvsusieBatches.get(jobRef).finished) {
            respond(ex, 409, "application/json", "{\"error\":\"A joint fine-mapping batch is running on this run; wait for it or cancel it\"}".getBytes()); return;
        }
        Set<String> readable = new HashSet<>();
        for (String pid : row.cells.keySet()) if (mvsusieReadable(ex, pid)) readable.add(pid);

        String id = UUID.randomUUID().toString();
        MvsusieJob job = new MvsusieJob();
        job.ownerId = mvsusieUser(ex);
        mvsusieJobs.put(id, job);
        final MultiLocusResult.LocusRow r = row;
        final String disease = wantDisease;
        new Thread(() -> runMvsusieJob(job, id, jobRef, result, r, mode, disease, readable), "mvsusie-" + id).start();
        respond(ex, 202, "application/json", ("{\"id\":\"" + id + "\"}").getBytes());
    }

    private void runMvsusieJob(MvsusieJob job, String id, String jobRef, MultiLocusResult result, MultiLocusResult.LocusRow row,
                               String mode, String wantDisease, Set<String> readable) {
        class Cand { String pid, disease; Config cfg; Locus locus; File harmonized; double p; long pos; List<Set<String>> cs; RunSummary su; }
        List<Cand> chosen = new ArrayList<>();
        String disease = "";
        long from = 0, to = 0;
        File work = new File(MVSUSIE_DIR, id).getAbsoluteFile();
        try {
            // 1. Candidates: significant here, inputs built, not flagged
            List<Cand> eligible = new ArrayList<>();
            List<String[]> excluded = new ArrayList<>();
            for (Map.Entry<String, MultiLocusResult.DatasetLocusStat> e : row.cells.entrySet()) {
                MultiLocusResult.DatasetLocusStat st = e.getValue();
                if (st == null || st.nP5e8 <= 0) continue;
                String pid = e.getKey();
                if (!readable.contains(pid)) continue;
                ProjectState ps = ensureProjectState(pid, new File("projects", pid).getAbsolutePath());
                if (ps == null || ps.loci == null) { excluded.add(new String[]{pid, "project not available"}); continue; }
                Locus own = null;
                for (Locus l : ps.loci)
                    if (l.chr.equals(row.chr) && l.start <= row.end && l.end >= row.start
                            && !BaseStepPipeline.overlapsMhc(l, ps.config.genomeBuild)
                            && (own == null || (l.start <= st.bestPos && st.bestPos <= l.end))) own = l;
                if (own == null) { excluded.add(new String[]{pid, "no fine-mappable locus of its own here (MHC or not identified)"}); continue; }
                File harm = new File(BaseStepPipeline.analysisDir(ps.config, own), "harmonized/harmonized_gwas.tsv");
                if (!harm.isFile()) { excluded.add(new String[]{pid, "inputs not built for this locus (run Analyze first)"}); continue; }
                RunSummary su = latestRuns(ps.config, own).get("susie_finemapping");
                if (su != null && "check".equals(su.reliability)) {
                    excluded.add(new String[]{pid, "SuSiE flagged this locus (reference LD does not match the GWAS)"}); continue;
                }
                Cand c = new Cand();
                c.pid = pid; c.disease = diseaseOf(result, pid); c.cfg = ps.config; c.locus = own; c.harmonized = harm; c.su = su;
                c.p = Double.isNaN(st.bestP) ? 1 : st.bestP; c.pos = st.bestPos; c.cs = susieCredibleSets(ps.config, own, su, row.chr + ":" + st.bestPos);
                eligible.add(c);
            }

            // 2. Within one disease, or colocalising datasets across diseases
            if (mode.equals("within")) {
                Map<String, Integer> counts = new TreeMap<>();
                for (Cand c : eligible) counts.merge(c.disease, 1, Integer::sum);
                disease = wantDisease != null && counts.containsKey(wantDisease)
                    ? wantDisease
                    : wantDisease != null && !wantDisease.isEmpty() ? wantDisease
                    : counts.entrySet().stream().max(Map.Entry.comparingByValue()).map(Map.Entry::getKey).orElse("");
                for (Cand c : eligible) {
                    if (c.disease.equals(disease)) chosen.add(c);
                    else excluded.add(new String[]{c.pid, "different disease (" + c.disease + "); within-disease mode"});
                }
            } else {
                List<Cand> withCs = new ArrayList<>();
                for (Cand c : eligible) {
                    if (c.cs.isEmpty()) excluded.add(new String[]{c.pid, "no SuSiE credible set to test colocalisation (run SuSiE first)"});
                    else withCs.add(c);
                }
                // Anchor = the consensus signal: the credible set that overlaps credible sets in the most other
                // datasets (ties: the dataset with the stronger association). A dataset joins only if one of
                // its own sets shares a SNP with that anchor; no chaining through secondary signals, which
                // would pull in datasets whose signal here is a different causal variant.
                withCs.sort(Comparator.comparingDouble(c -> c.p));
                Set<String> anchor = null; String anchorFrom = null; int bestSupport = -1;
                for (Cand c : withCs)
                    for (Set<String> set : c.cs) {
                        int support = 0;
                        for (Cand o : withCs) {
                            if (o == c) continue;
                            for (Set<String> os : o.cs) if (!Collections.disjoint(os, set)) { support++; break; }
                        }
                        if (support > bestSupport) { bestSupport = support; anchor = set; anchorFrom = c.pid; }
                    }
                if (anchor != null) {
                    for (Cand c : withCs) {
                        boolean overlaps = false;
                        for (Set<String> set : c.cs) if (!Collections.disjoint(set, anchor)) { overlaps = true; break; }
                        if (overlaps) chosen.add(c);
                        else excluded.add(new String[]{c.pid, "no credible set overlaps the shared signal (the "
                            + anchor.size() + "-SNP set of " + anchorFrom + " supported by " + bestSupport
                            + " other datasets): likely a different causal variant"});
                    }
                }
            }

            StringBuilder sel = new StringBuilder("{\"mode\":\"").append(mode).append("\",\"disease\":\"").append(escJ(disease)).append("\",\"chosen\":[");
            for (int i = 0; i < chosen.size(); i++) sel.append(i > 0 ? "," : "").append('"').append(escJ(chosen.get(i).pid)).append('"');
            sel.append("],\"excluded\":[");
            for (int i = 0; i < excluded.size(); i++)
                sel.append(i > 0 ? "," : "").append("{\"dataset\":\"").append(escJ(excluded.get(i)[0])).append("\",\"reason\":\"").append(escJ(excluded.get(i)[1])).append("\"}");
            sel.append("]}");
            job.selectionJson = sel.toString();
            if (chosen.size() < 2) throw new IOException("Need at least 2 datasets for joint fine-mapping here; " + chosen.size() + " qualified");

            // 3. Manifest for the R script
            Cand lead = Collections.min(chosen, Comparator.comparingDouble(c -> c.p));
            from = Math.max(1, lead.pos - MVSUSIE_HALF_WINDOW);
            to = lead.pos + MVSUSIE_HALF_WINDOW;
            work.mkdirs();
            String plink = PlinkSubsetter.findPlink(lead.cfg);
            if (plink == null) throw new IOException("PLINK not found");
            Map<String, Object> m = new LinkedHashMap<>();
            m.put("chr", row.chr);
            m.put("from", from);
            m.put("to", to);
            m.put("ref_panel", lead.cfg.refPanelPath);
            m.put("plink", new File(plink).getAbsolutePath());
            m.put("ukbb_dir", SusieAdapter.ukbbLdDir());
            m.put("ukbb_script", new File(SusieAdapter.appRoot(), "scripts/ukbb_ld.py").getAbsolutePath());
            m.put("extra_lib", mvsusieRLib().replace("\\", "/"));
            m.put("L", 10); m.put("coverage", 0.95); m.put("ld_shrink", 0.1);
            m.put("work_dir", work.getAbsolutePath().replace("\\", "/"));
            m.put("out", new File(work, "result.json").getAbsolutePath().replace("\\", "/"));
            List<Object> ds = new ArrayList<>();
            for (Cand c : chosen) {
                Map<String, Object> d = new LinkedHashMap<>();
                d.put("id", c.pid); d.put("disease", c.disease); d.put("N", c.cfg.sampleN); d.put("p", c.p);
                Map<String, Object> gw = new LinkedHashMap<>();
                gw.put("file", c.cfg.gwasFile.replace("\\", "/")); gw.put("chr", c.cfg.colChr); gw.put("pos", c.cfg.colPos);
                gw.put("ea", c.cfg.colEa); gw.put("nea", c.cfg.colNea); gw.put("se", c.cfg.colSe);
                gw.put("beta", c.cfg.colBeta); gw.put("or", c.cfg.colOr);
                d.put("gwas", gw);
                d.put("null_cache", mvsusieNullCache(c.pid, c.cfg).getAbsolutePath().replace("\\", "/"));
                d.put("files", new ArrayList<Object>(Collections.singletonList(c.harmonized.getAbsolutePath().replace("\\", "/"))));
                ds.add(d);
            }
            m.put("datasets", ds);
            File manifest = new File(work, "manifest.json");
            Files.write(manifest.toPath(), MiniJson.encode(m).getBytes("UTF-8"));

            // 4. Run
            job.phase = "Running mvSuSiE on " + chosen.size() + " datasets";
            ProcessBuilder pb = new ProcessBuilder(ToolLocator.rscript(),
                new File(SusieAdapter.appRoot(), "tools/r/mvsusie_run.R").getAbsolutePath(), manifest.getAbsolutePath());
            pb.redirectErrorStream(true);
            pb.redirectOutput(new File(work, "run.log"));
            Process proc = pb.start();
            if (!proc.waitFor(20, TimeUnit.MINUTES)) { proc.destroyForcibly(); throw new IOException("mvSuSiE did not finish within 20 minutes"); }
            File out = new File(work, "result.json");
            if (!out.isFile()) throw new IOException("mvSuSiE produced no result; see " + new File(work, "run.log"));
            String res = new String(Files.readAllBytes(out.toPath()), "UTF-8");
            Map<String, Object> parsed = MiniJson.asObject(MiniJson.parse(res));
            if (Boolean.FALSE.equals(parsed.get("ok"))) throw new IOException(String.valueOf(parsed.get("error")));
            job.resultJson = res;
            job.phase = "Complete";
            job.status = "done";
        } catch (Exception e) {
            job.error = mvsusiePublicError(e.getMessage());
            job.status = "error";
            job.phase = "Error";
            System.err.printf("[mvSuSiE] %s failed: %s%n", id, e.getMessage());
        }
        try {
            saveMvsusieRun(job, id, jobRef, row, mode, disease, chosen.stream().map(c -> new Object[]{c.pid, c.disease, c.locus, c.cfg, c.su, c.cs})
                .collect(java.util.stream.Collectors.toList()), from, to);
        } catch (Exception e) {
            System.err.printf("[mvSuSiE] %s could not be saved: %s%n", id, e.getMessage());
        }
    }

    /** R library holding mvsusieR (needs a newer susieR than the main library). */
    private String mvsusieRLib() {
        String lib = System.getenv("LYNXGWAS_R_LIB");
        if (lib != null && !lib.isEmpty()) return lib;
        File appLib = new File(SusieAdapter.appRoot(), "rlib"), userLib = new File(System.getProperty("user.home"), "R/mvsusie-lib");
        return appLib.isDirectory() ? appLib.getAbsolutePath() : userLib.isDirectory() ? userLib.getAbsolutePath() : "";
    }

    /** Where a dataset's genome-wide null z sample is cached (rebuilt by R when the GWAS file changes). */
    private File mvsusieNullCache(String pid, Config cfg) {
        return new File(cfg.outputDir, "cache/mvsusie_null_z.rds");
    }

    private String mvsusiePublicError(String msg) { return msg == null ? "mvSuSiE failed" : msg; }

    /** Removes every saved joint fine-mapping run made on a cross-dataset run (when that run is deleted). */
    private void purgeMvsusieRunsOf(String jobId) {
        ensureMvsusieIndex();
        for (MvsusieSaved sv : new ArrayList<>(mvsusieSaved.values())) {
            if (!sv.job.equals(jobId)) continue;
            mvsusieSaved.remove(sv.id);
            File dir = new File(MVSUSIE_DIR, sv.id);
            File[] fs = dir.listFiles();
            if (fs != null) for (File f : fs) f.delete();
            dir.delete();
        }
        mvsusieBatches.remove(jobId);
    }

    /**
     * Writes meta.json next to the run (selection, datasets with their own locus, single-dataset SuSiE PIPs
     * and credible-set sizes for the comparison panel) and indexes it. Error runs are saved too, so a
     * batch does not retry loci that cannot be fine-mapped jointly.
     */
    @SuppressWarnings("unchecked")
    private void saveMvsusieRun(MvsusieJob job, String id, String jobRef, MultiLocusResult.LocusRow row, String mode,
                                String disease, List<Object[]> chosen, long from, long to) throws IOException {
        File work = new File(MVSUSIE_DIR, id);
        work.mkdirs();
        Map<String, Object> meta = new LinkedHashMap<>();
        meta.put("id", id); meta.put("job", jobRef); meta.put("locus", row.index); meta.put("chr", row.chr);
        meta.put("start", row.start); meta.put("end", row.end); meta.put("gene", row.nearestGene == null ? "" : row.nearestGene);
        meta.put("mode", mode); meta.put("disease", disease);
        meta.put("created", System.currentTimeMillis());
        meta.put("status", job.status); meta.put("error", job.error); meta.put("owner", job.ownerId);
        try { meta.put("selection", MiniJson.parse(job.selectionJson)); } catch (Exception e) { meta.put("selection", new LinkedHashMap<>()); }

        Map<String, Object> result = null;
        if (job.resultJson != null) result = MiniJson.asObject(MiniJson.parse(job.resultJson));
        List<Object> sets = result != null && result.get("credible_sets") instanceof List ? (List<Object>) result.get("credible_sets") : new ArrayList<>();
        meta.put("n_sets", sets.size());
        if (!sets.isEmpty()) {
            Map<String, Object> s0 = MiniJson.asObject(sets.get(0));
            meta.put("lead_snp", s0.get("lead_snp")); meta.put("lead_pip", s0.get("lead_pip"));
        }
        // Datasets the R step actually pooled (it can drop one that is too sparse in the window)
        if (result != null && result.get("datasets") instanceof List) {
            Set<String> used = new HashSet<>();
            for (Object o : (List<Object>) result.get("datasets")) used.add(String.valueOf(o));
            chosen = chosen.stream().filter(c -> used.contains((String) c[0])).collect(java.util.stream.Collectors.toList());
        }
        List<Object> ds = new ArrayList<>();
        Map<String, Object> singlePips = new LinkedHashMap<>();
        for (Object[] c : chosen) {
            Map<String, Object> d = new LinkedHashMap<>();
            Locus l = (Locus) c[2];
            d.put("id", c[0]); d.put("disease", c[1]); d.put("locus_index", l.index); d.put("locus_id", l.id);
            ds.add(d);
            if (job.resultJson != null) singlePips.put((String) c[0], susiePips((Config) c[3], l, (RunSummary) c[4], from, to));
        }
        meta.put("datasets", ds);
        meta.put("single_pips", singlePips);
        // Per joint credible set: each dataset's smallest single-dataset set holding the joint lead SNP
        List<Object> setCmp = new ArrayList<>();
        for (Object so : sets) {
            Map<String, Object> s = MiniJson.asObject(so);
            String key = row.chr + ":" + (long) num(s.get("lead_pos"), -1);
            Map<String, Object> sizes = new LinkedHashMap<>();
            for (Object[] c : chosen) {
                int best = Integer.MAX_VALUE;
                for (Set<String> cs : (List<Set<String>>) c[5]) if (cs.contains(key)) best = Math.min(best, cs.size());
                sizes.put((String) c[0], best == Integer.MAX_VALUE ? null : best);
            }
            setCmp.add(sizes);
        }
        meta.put("single_set_size_at_lead", setCmp);
        String metaJson = MiniJson.encode(meta);
        Files.write(new File(work, "meta.json").toPath(), metaJson.getBytes("UTF-8"));
        job.metaJson = metaJson;
        ensureMvsusieIndex();
        mvsusieSaved.put(id, parseSavedMeta(meta));
    }

    // GET /api/mvsusie-progress?id=...
    private void mvsusieProgress(HttpExchange ex) throws IOException {
        cors(ex); if (preflight(ex)) return;
        MvsusieJob job = mvsusieJobs.get(String.valueOf(queryParam(ex, "id")));
        if (job != null && job.ownerId != null && !job.ownerId.equals(mvsusieUser(ex))) job = null;
        if (job == null) { respond(ex, 404, "application/json", "{\"error\":\"Job not found\"}".getBytes()); return; }
        respond(ex, 200, "application/json", ("{\"status\":\"" + job.status + "\",\"phase\":\"" + escJ(job.phase) + "\""
            + ",\"elapsed_s\":" + (System.currentTimeMillis() - job.startedAt) / 1000
            + (job.error != null ? ",\"error\":\"" + escJ(job.error) + "\"" : "") + ",\"selection\":" + job.selectionJson + "}").getBytes("UTF-8"));
    }

    // GET /api/mvsusie-result?id=... -> {selection, result (null for a failed run), error?, meta}
    private void mvsusieResult(HttpExchange ex) throws IOException {
        cors(ex); if (preflight(ex)) return;
        String id = String.valueOf(queryParam(ex, "id"));
        ensureMvsusieIndex();
        MvsusieSaved saved = mvsusieSaved.get(id);
        if (saved == null || !mvsusieVisible(ex, saved) || !mvsusieCanSeeRun(ex, saved.job)) {
            respond(ex, 404, "application/json", "{\"error\":\"Result not available\"}".getBytes()); return;
        }
        File work = new File(MVSUSIE_DIR, saved.id);
        String meta = new String(Files.readAllBytes(new File(work, "meta.json").toPath()), "UTF-8");
        File res = new File(work, "result.json");
        String result = "ok".equals(saved.status) || "done".equals(saved.status)
            ? (res.isFile() ? new String(Files.readAllBytes(res.toPath()), "UTF-8") : "null") : "null";
        Map<String, Object> mm = MiniJson.asObject(MiniJson.parse(meta));
        respond(ex, 200, "application/json", ("{\"selection\":" + MiniJson.encode(mm.get("selection"))
            + ",\"result\":" + result + (saved.error != null ? ",\"error\":\"" + escJ(saved.error) + "\"" : "")
            + ",\"meta\":" + meta + "}").getBytes("UTF-8"));
    }

    // GET /api/mvsusie-latest?job=...&locus=... -> {"runs":[saved run summaries, newest first]}
    private void mvsusieLatest(HttpExchange ex) throws IOException {
        cors(ex); if (preflight(ex)) return;
        String job = queryParam(ex, "job");
        int locus;
        try { locus = Integer.parseInt(String.valueOf(queryParam(ex, "locus"))); }
        catch (NumberFormatException e) { respond(ex, 400, "application/json", "{\"error\":\"locus required\"}".getBytes()); return; }
        if (!mvsusieCanSeeRun(ex, job)) { respond(ex, 404, "application/json", "{\"error\":\"Cross-dataset run not found\"}".getBytes()); return; }
        StringBuilder j = new StringBuilder("{\"runs\":[");
        int i = 0;
        for (MvsusieSaved s : savedAtLocus(job, locus))
            if (mvsusieVisible(ex, s)) j.append(i++ > 0 ? "," : "").append(s.toJson());
        respond(ex, 200, "application/json", j.append("]}").toString().getBytes("UTF-8"));
    }

    // GET /api/mvsusie-for-dataset?project=...&locus_index=... -> successful joint runs that pooled this
    // dataset's locus (latest per cross-dataset run, mode and disease), for the regional-plot strip
    private void mvsusieForDataset(HttpExchange ex) throws IOException {
        cors(ex); if (preflight(ex)) return;
        String pid = queryParam(ex, "project");
        int li;
        try { li = Integer.parseInt(String.valueOf(queryParam(ex, "locus_index"))); }
        catch (NumberFormatException e) { respond(ex, 400, "application/json", "{\"error\":\"locus_index required\"}".getBytes()); return; }
        if (pid == null || !mvsusieReadable(ex, pid)) { respond(ex, 404, "application/json", "{\"error\":\"Project not found\"}".getBytes()); return; }
        ensureMvsusieIndex();
        Map<String, MvsusieSaved> latest = new LinkedHashMap<>();
        for (MvsusieSaved s : mvsusieSaved.values()) {
            if (!"done".equals(s.status) || !Integer.valueOf(li).equals(s.datasetLoci.get(pid))) continue;
            if (!mvsusieVisible(ex, s) || !mvsusieCanSeeRun(ex, s.job)) continue;
            String k = s.job + "|" + s.mode + "|" + s.disease;
            MvsusieSaved cur = latest.get(k);
            if (cur == null || s.created > cur.created) latest.put(k, s);
        }
        StringBuilder j = new StringBuilder("{\"runs\":[");
        int i = 0;
        for (MvsusieSaved s : latest.values()) j.append(i++ > 0 ? "," : "").append(s.toJson());
        respond(ex, 200, "application/json", j.append("]}").toString().getBytes("UTF-8"));
    }

    // ── Joint fine-mapping over every aligned locus of a cross-dataset run ──
    //
    // POST /api/mvsusie-run-all {"job":..., "mode":"within"|"cross", "only_missing":true, "parallel":3}
    //   within: every (locus, disease) with >= 2 genome-wide significant datasets of that disease
    //   cross:  every locus genome-wide significant in >= 2 diseases
    // GET  /api/mvsusie-run-all-progress?job=...   POST /api/mvsusie-run-all-cancel {"job":...}
    // GET  /api/mvsusie-summary?job=...&format=tsv  (one row per joint credible set, latest run per unit)

    static class MvsusieBatch {
        final String mode; String ownerId;
        volatile int total, done, ok, failed, skipped;
        volatile String current = "";
        volatile boolean finished, cancelled;
        final long startedAt = System.currentTimeMillis();
        MvsusieBatch(String mode) { this.mode = mode; }
        String toJson() {
            return "{\"mode\":\"" + mode + "\",\"total\":" + total + ",\"done\":" + done + ",\"ok\":" + ok + ",\"failed\":" + failed
                + ",\"skipped\":" + skipped + ",\"current\":\"" + escJ(current) + "\",\"finished\":" + finished
                + ",\"cancelled\":" + cancelled + ",\"elapsed_s\":" + (System.currentTimeMillis() - startedAt) / 1000 + "}";
        }
    }
    private final Map<String, MvsusieBatch> mvsusieBatches = new ConcurrentHashMap<>();

    private void mvsusieRunAll(HttpExchange ex) throws IOException {
        cors(ex); if (preflight(ex)) return;
        if (!"POST".equalsIgnoreCase(ex.getRequestMethod())) {
            respond(ex, 405, "application/json", "{\"error\":\"POST required\"}".getBytes()); return;
        }
        String body = new String(readAll(ex.getRequestBody()), "UTF-8");
        String jobRef = extractStr(body, "job");
        String mode = "cross".equals(extractStr(body, "mode")) ? "cross" : "within";
        boolean onlyMissing = !"false".equals(extractStr(body, "only_missing"));
        int parallel = (int) Math.max(1, Math.min(4, parseDoubleOr(extractStr(body, "parallel"), 3)));
        MultiLocusResult result = jobRef == null || !mvsusieCanSeeRun(ex, jobRef) ? null : locusMatrixJobs.get(jobRef);
        if (result == null) { respond(ex, 404, "application/json", "{\"error\":\"Cross-dataset run not found\"}".getBytes()); return; }
        if (!mvsusieBatchAllowed(ex, jobRef)) return;
        MvsusieBatch prev = mvsusieBatches.get(jobRef);
        if (prev != null && !prev.finished) {
            respond(ex, 409, "application/json", "{\"error\":\"A joint fine-mapping batch is already running on this run\"}".getBytes()); return;
        }
        Set<String> readable = new HashSet<>();
        for (MultiLocusResult.DatasetInfo di : result.datasets) if (mvsusieReadable(ex, di.id)) readable.add(di.id);

        // Units of work
        List<Object[]> units = new ArrayList<>();   // {row, disease}
        for (MultiLocusResult.LocusRow row : result.loci) {
            Map<String, List<String>> byDis = significantByDisease(result, row);
            if (mode.equals("within")) {
                for (Map.Entry<String, List<String>> e : byDis.entrySet())
                    if (e.getValue().size() >= 2) units.add(new Object[]{row, e.getKey()});
            } else if (byDis.size() >= 2) units.add(new Object[]{row, ""});
        }
        MvsusieBatch b = new MvsusieBatch(mode);
        b.ownerId = mvsusieUser(ex);
        b.total = units.size();
        mvsusieBatches.put(jobRef, b);
        final String owner = b.ownerId;
        Runnable all = () -> {
            java.util.concurrent.ExecutorService pool = Executors.newFixedThreadPool(parallel);
            List<java.util.concurrent.Future<?>> fs = new ArrayList<>();
            for (Object[] u : units) {
                MultiLocusResult.LocusRow row = (MultiLocusResult.LocusRow) u[0];
                String dis = (String) u[1];
                fs.add(pool.submit(() -> {
                    if (b.cancelled) return;
                    boolean have = false;
                    if (onlyMissing)
                        for (MvsusieSaved s : savedAtLocus(jobRef, row.index))
                            if (s.mode.equals(mode) && s.disease.equals(dis) && Objects.equals(s.ownerId, owner)) { have = true; break; }
                    if (have) { synchronized (b) { b.skipped++; b.done++; } return; }
                    b.current = (row.nearestGene == null || row.nearestGene.isEmpty() ? row.chr + ":" + row.start : row.nearestGene)
                        + (dis.isEmpty() ? "" : " (" + dis + ")");
                    String id = UUID.randomUUID().toString();
                    MvsusieJob job = new MvsusieJob();
                    job.ownerId = owner;
                    runMvsusieJob(job, id, jobRef, result, row, mode, dis.isEmpty() ? null : dis, readable);
                    synchronized (b) { if ("done".equals(job.status)) b.ok++; else b.failed++; b.done++; }
                }));
            }
            for (java.util.concurrent.Future<?> f : fs) try { f.get(); } catch (Exception ignored) {}
            pool.shutdown();
            b.current = "";
            b.finished = true;
        };
        String refused = mvsusieStartBatch(ex, owner, jobRef, all);
        if (refused != null) { mvsusieBatches.remove(jobRef); respond(ex, 429, "application/json", ("{\"error\":\"" + escJ(refused) + "\"}").getBytes("UTF-8")); return; }
        respond(ex, 202, "application/json", b.toJson().getBytes("UTF-8"));
    }

    /** Desktop: any batch may run. */
    private boolean mvsusieBatchAllowed(HttpExchange ex, String jobRef) throws IOException { return true; }

    /** Starts the batch; returns a refusal message, or null when it started. */
    private String mvsusieStartBatch(HttpExchange ex, String owner, String jobRef, Runnable r) {
        new Thread(r, "mvsusie-batch-" + jobRef).start();
        return null;
    }

    private void mvsusieRunAllProgress(HttpExchange ex) throws IOException {
        cors(ex); if (preflight(ex)) return;
        String job = queryParam(ex, "job");
        MvsusieBatch b = job == null ? null : mvsusieBatches.get(job);
        if (b != null && b.ownerId != null && !b.ownerId.equals(mvsusieUser(ex))) b = null;
        respond(ex, 200, "application/json", (b == null ? "{\"running\":false}" : b.toJson()).getBytes("UTF-8"));
    }

    private void mvsusieRunAllCancel(HttpExchange ex) throws IOException {
        cors(ex); if (preflight(ex)) return;
        if (!"POST".equalsIgnoreCase(ex.getRequestMethod())) {
            respond(ex, 405, "application/json", "{\"error\":\"POST required\"}".getBytes()); return;
        }
        String job = extractStr(new String(readAll(ex.getRequestBody()), "UTF-8"), "job");
        MvsusieBatch b = job == null ? null : mvsusieBatches.get(job);
        if (b != null && b.ownerId != null && !b.ownerId.equals(mvsusieUser(ex))) b = null;
        if (b == null) { respond(ex, 404, "application/json", "{\"error\":\"No batch running\"}".getBytes()); return; }
        b.cancelled = true;
        respond(ex, 200, "application/json", "{\"ok\":true}".getBytes());
    }

    // GET /api/mvsusie-summary?job=...[&format=tsv]
    @SuppressWarnings("unchecked")
    private void mvsusieSummary(HttpExchange ex) throws IOException {
        cors(ex); if (preflight(ex)) return;
        String job = queryParam(ex, "job");
        boolean tsv = "tsv".equals(queryParam(ex, "format"));
        if (!mvsusieCanSeeRun(ex, job)) { respond(ex, 404, "application/json", "{\"error\":\"Cross-dataset run not found\"}".getBytes()); return; }
        ensureMvsusieIndex();
        // Latest visible run per (locus, mode, disease)
        Map<String, MvsusieSaved> latest = new TreeMap<>();
        for (MvsusieSaved s : mvsusieSaved.values()) {
            if (!s.job.equals(job) || !mvsusieVisible(ex, s)) continue;
            String k = String.format(Locale.ROOT, "%08d|%s|%s", s.locus, s.mode, s.disease);
            MvsusieSaved cur = latest.get(k);
            if (cur == null || s.created > cur.created) latest.put(k, s);
        }
        if (!tsv) {
            int ok = 0, sets = 0, single = 0;
            for (MvsusieSaved s : latest.values()) if ("done".equals(s.status)) { ok++; sets += s.nSets; if (s.nSets == 1) single++; }
            respond(ex, 200, "application/json", String.format(Locale.ROOT,
                "{\"units\":%d,\"ok\":%d,\"failed\":%d,\"credible_sets\":%d,\"single_signal_units\":%d}",
                latest.size(), ok, latest.size() - ok, sets, single).getBytes("UTF-8"));
            return;
        }
        StringBuilder out = new StringBuilder("locus\tchr\tstart\tend\tnearest_gene\tmode\tdisease\tstatus\tn_datasets\tdatasets\tn_sets"
            + "\tset\tset_size\tpip_min\tpip_max\tlead_snp\tlead_pos\tlead_pip\tpurity\tn_active\tactive_in"
            + "\tsmallest_single_dataset_set_at_lead\tld_source\tnote\trun_id\n");
        for (MvsusieSaved s : latest.values()) {
            File work = new File(MVSUSIE_DIR, s.id);
            Map<String, Object> meta, res = null;
            try {
                meta = MiniJson.asObject(MiniJson.parse(new String(Files.readAllBytes(new File(work, "meta.json").toPath()), "UTF-8")));
                File rf = new File(work, "result.json");
                if ("done".equals(s.status) && rf.isFile()) res = MiniJson.asObject(MiniJson.parse(new String(Files.readAllBytes(rf.toPath()), "UTF-8")));
            } catch (Exception e) { continue; }
            String head = s.locus + "\t" + s.chr + "\t" + (long) num(meta.get("start"), 0) + "\t" + (long) num(meta.get("end"), 0) + "\t"
                + tsvEsc(s.gene) + "\t" + s.mode + "\t" + tsvEsc(s.disease) + "\t" + ("done".equals(s.status) ? "ok" : "failed") + "\t"
                + s.nDatasets + "\t" + tsvEsc(String.join(",", s.datasetLoci.keySet())) + "\t" + s.nSets;
            String ld = res == null ? "" : tsvEsc(MiniJson.getStr(res, "ld_source", ""));
            List<Object> sets = res != null && res.get("credible_sets") instanceof List ? (List<Object>) res.get("credible_sets") : new ArrayList<>();
            if (sets.isEmpty()) {
                out.append(head).append("\t\t\t\t\t\t\t\t\t\t\t\t").append(ld).append('\t')
                   .append(tsvEsc(s.error != null ? s.error : "no credible set reached 95% coverage")).append('\t').append(s.id).append('\n');
                continue;
            }
            List<Object> cmp = meta.get("single_set_size_at_lead") instanceof List ? (List<Object>) meta.get("single_set_size_at_lead") : new ArrayList<>();
            for (int i = 0; i < sets.size(); i++) {
                Map<String, Object> c = MiniJson.asObject(sets.get(i));
                List<Object> act = c.get("active_in") instanceof List ? (List<Object>) c.get("active_in") : new ArrayList<>();
                String smallest = "";
                if (i < cmp.size() && cmp.get(i) instanceof Map) {
                    int best = Integer.MAX_VALUE; String who = "";
                    for (Map.Entry<String, Object> e : MiniJson.asObject(cmp.get(i)).entrySet())
                        if (e.getValue() instanceof Number && ((Number) e.getValue()).intValue() < best) { best = ((Number) e.getValue()).intValue(); who = e.getKey(); }
                    if (best != Integer.MAX_VALUE) smallest = who + ":" + best;
                }
                List<String> actS = new ArrayList<>();
                for (Object a : act) actS.add(String.valueOf(a));
                out.append(head).append('\t').append(i + 1).append('\t').append((long) num(c.get("size"), 0))
                   .append('\t').append(MiniJson.getStr(c, "pip_min", "")).append('\t').append(MiniJson.getStr(c, "pip_max", ""))
                   .append('\t').append(tsvEsc(MiniJson.getStr(c, "lead_snp", ""))).append('\t').append((long) num(c.get("lead_pos"), 0))
                   .append('\t').append(MiniJson.getStr(c, "lead_pip", "")).append('\t').append(MiniJson.getStr(c, "min_abs_corr", ""))
                   .append('\t').append(actS.size()).append('\t').append(tsvEsc(String.join(",", actS)))
                   .append('\t').append(tsvEsc(smallest)).append('\t').append(ld).append("\t\t").append(s.id).append('\n');
            }
        }
        ex.getResponseHeaders().set("Content-Disposition", "attachment; filename=\"joint_finemapping_" + (job.length() > 8 ? job.substring(0, 8) : job) + ".tsv\"");
        respond(ex, 200, "text/tab-separated-values", out.toString().getBytes("UTF-8"));
    }

    // ── Batch runs: one tool over every ready locus of a project ─────────────

    static class BatchRun {
        final String tool;
        volatile int total, done, ok, failed, skipped;
        volatile String current = "";
        volatile boolean finished, cancelled;
        final long startedAt = System.currentTimeMillis();
        volatile long finishedAt;
        final List<String> errors = Collections.synchronizedList(new ArrayList<>());
        BatchRun(String tool) { this.tool = tool; }

        String toJson() {
            StringBuilder j = new StringBuilder("{");
            j.append("\"tool\":\"").append(escJ(tool)).append('"');
            j.append(",\"total\":").append(total).append(",\"done\":").append(done);
            j.append(",\"ok\":").append(ok).append(",\"failed\":").append(failed).append(",\"skipped\":").append(skipped);
            j.append(",\"current\":\"").append(escJ(current)).append('"');
            j.append(",\"finished\":").append(finished).append(",\"cancelled\":").append(cancelled);
            j.append(",\"elapsed_s\":").append(((finished ? finishedAt : System.currentTimeMillis()) - startedAt) / 1000);
            j.append(",\"errors\":[");
            synchronized (errors) {
                for (int i = 0; i < errors.size(); i++) { if (i > 0) j.append(','); j.append('"').append(escJ(errors.get(i))).append('"'); }
            }
            return j.append("]}").toString();
        }
    }
    private final Map<String, BatchRun> batchRuns = new ConcurrentHashMap<>();

    /** True when this locus already has a successful run of {@code tool} (provenance + result written). */
    static boolean hasSuccessfulRun(Config cfg, Locus locus, String tool) {
        File runs = new File(BaseStepPipeline.analysisDir(cfg, locus), "runs");
        File[] dirs = runs.listFiles(File::isDirectory);
        if (dirs == null) return false;
        for (File d : dirs) {
            File prov = new File(d, "provenance.json");
            if (!prov.isFile() || !new File(d, "result.tsv").isFile()) continue;
            try {
                String j = new String(Files.readAllBytes(prov.toPath()), "UTF-8");
                if (j.contains("\"tool\": \"" + tool + "\"") || j.contains("\"tool\":\"" + tool + "\"")) return true;
            } catch (IOException ignored) {}
        }
        return false;
    }

    // POST /api/project/{id}/analysis/run-all
    //   {"tool": "...", "only_missing": "true", "include_mhc": "false", "parallel": "2", <tool params>}
    // Runs the tool on every locus whose base artifacts are ready, skipping the extended MHC
    // (dense LD, not fine-mappable with a reference panel) and, by default, loci already done.
    private void projectAnalysisRunAll(HttpExchange ex, String projectId, String projectDir) throws IOException {
        if (!"POST".equalsIgnoreCase(ex.getRequestMethod())) {
            respond(ex, 405, "text/plain", "Method Not Allowed".getBytes()); return;
        }
        ProjectState ps = ensureProjectState(projectId, projectDir);
        if (ps == null || ps.config == null || ps.loci == null) {
            respond(ex, 503, "application/json", "{\"error\":\"Pipeline state not available\"}".getBytes()); return;
        }
        String body = new String(readAll(ex.getRequestBody()), "UTF-8");
        String toolName = extractStr(body, "tool");
        ToolDescriptor td = toolName == null ? null : PluginEngine.findTool(toolName);
        if (td == null) {
            respond(ex, 400, "application/json", "{\"error\":\"Unknown tool\"}".getBytes()); return;
        }
        if (toolName.equals("coloc") || toolName.equals("gwama_meta")) {
            respond(ex, 400, "application/json", ("{\"error\":\"" + escJ(toolName)
                + " compares datasets; run it per locus or as a cross-dataset run\"}").getBytes()); return;
        }
        String key = projectId + "|" + toolName;
        BatchRun existing = batchRuns.get(key);
        if (existing != null && !existing.finished) {
            respond(ex, 409, "application/json", ("{\"error\":\"A batch run of " + escJ(toolName) + " is already running\"}").getBytes()); return;
        }
        boolean onlyMissing = !"false".equalsIgnoreCase(extractStr(body, "only_missing"));
        boolean includeMhc = "true".equalsIgnoreCase(extractStr(body, "include_mhc"));
        int parallel = 2;
        try { parallel = Math.max(1, Math.min(8, Integer.parseInt(extractStr(body, "parallel")))); } catch (Exception ignored) {}

        Map<String, String> params = new LinkedHashMap<>();
        for (ToolDescriptor.Param p : td.params) {
            String val = extractStr(body, p.name);
            if (val != null) params.put(p.name, val);
            else if (p.defaultValue != null) params.put(p.name, p.defaultValue);
        }

        final Config cfg = ps.config;
        List<Locus> todo = new ArrayList<>();
        BatchRun br = new BatchRun(toolName);
        for (Locus l : ps.loci) {
            Map<String, Boolean> st = BaseStepPipeline.checkStatus(cfg, l);
            boolean ready = !st.isEmpty() && !st.containsValue(false);
            if (!ready || (!includeMhc && BaseStepPipeline.overlapsMhc(l, cfg.genomeBuild))
                    || (onlyMissing && hasSuccessfulRun(cfg, l, toolName))) { br.skipped++; continue; }
            todo.add(l);
        }
        br.total = todo.size();
        batchRuns.put(key, br);
        final int threads = parallel;
        new Thread(() -> {
            ExecutorService pool = Executors.newFixedThreadPool(threads);
            try {
                List<Future<?>> fs = new ArrayList<>();
                for (Locus l : todo) fs.add(pool.submit(() -> {
                    if (br.cancelled) return;
                    br.current = "Locus " + l.index + " (chr" + l.chr + ")";
                    String jobId = UUID.randomUUID().toString();
                    ProgressTracker pt = new ProgressTracker();
                    analysisJobProgress.put(jobId, pt);
                    PluginEngine.RunResult r = runAnalysisJob(projectId, cfg, l, toolName, new LinkedHashMap<>(params), jobId, pt);
                    synchronized (br) {
                        br.done++;
                        if (r.ok) br.ok++;
                        else {
                            br.failed++;
                            if (br.errors.size() < 100) br.errors.add("Locus " + l.index + ": " + (r.error == null ? "failed" : r.error));
                        }
                    }
                }));
                for (Future<?> f : fs) { try { f.get(); } catch (Exception ignored) {} }
            } finally {
                pool.shutdown();
                br.current = "";
                br.finished = true;
                br.finishedAt = System.currentTimeMillis();
                System.out.printf("[Batch] %s on '%s': %d ok, %d failed, %d skipped in %ds%n",
                    toolName, projectId, br.ok, br.failed, br.skipped, (br.finishedAt - br.startedAt) / 1000);
            }
        }, "batch-" + toolName + "-" + projectId).start();

        respond(ex, 202, "application/json", br.toJson().getBytes("UTF-8"));
    }

    // GET /api/project/{id}/analysis/run-all-progress?tool=...
    private void projectAnalysisRunAllProgress(HttpExchange ex, String projectId) throws IOException {
        BatchRun br = batchRuns.get(projectId + "|" + queryParam(ex, "tool"));
        if (br == null) {
            respond(ex, 404, "application/json", "{\"error\":\"No batch run for this tool\"}".getBytes()); return;
        }
        respond(ex, 200, "application/json", br.toJson().getBytes("UTF-8"));
    }

    // GET /api/project/{id}/analysis/job/{jobId}/progress or /result or /cancel
    private void projectAnalysisJob(HttpExchange ex, String projectId,
                                     String projectDir, String action) throws IOException {
        // action = "analysis/job/{jobId}/progress" or "/result" or "/cancel" or "/log"
        String rest = action.substring("analysis/job/".length());
        int slash = rest.indexOf('/');
        String jobId = slash >= 0 ? rest.substring(0, slash) : rest;
        String subAction = slash >= 0 ? rest.substring(slash + 1) : "progress";

        if ("progress".equals(subAction)) {
            ProgressTracker pt = analysisJobProgress.get(jobId);
            PluginEngine.RunResult result = analysisJobs.get(jobId);
            if (pt == null) {
                respond(ex, 404, "application/json", "{\"error\":\"Job not found\"}".getBytes());
                return;
            }
            StringBuilder json = new StringBuilder("{");
            json.append("\"job_id\":\"").append(escJ(jobId)).append('"');
            json.append(",\"phase\":\"").append(escJ(pt.phase)).append('"');
            json.append(",\"pct\":").append(pt.pct());
            json.append(",\"done\":").append(pt.done);
            if (result != null) {
                json.append(",\"ok\":").append(result.ok);
                if (result.error != null) json.append(",\"error\":\"").append(escJ(result.error)).append('"');
                if (result.ok) json.append(",\"result_rows\":").append(result.resultRows);
            }
            json.append('}');
            respond(ex, 200, "application/json", json.toString().getBytes("UTF-8"));

        } else if ("result".equals(subAction)) {
            PluginEngine.RunResult result = analysisJobs.get(jobId);
            if (result == null || !result.ok) {
                respond(ex, 404, "application/json",
                    "{\"error\":\"Result not available\"}".getBytes()); return;
            }
            File resultFile = new File(result.runDir, "result.tsv");
            if (resultFile.exists()) {
                byte[] bytes = java.nio.file.Files.readAllBytes(resultFile.toPath());
                respond(ex, 200, "text/tab-separated-values", bytes);
            } else {
                respond(ex, 404, "application/json",
                    "{\"error\":\"result.tsv not found\"}".getBytes());
            }

        } else if ("cancel".equals(subAction)) {
            cancelledJobs.add(jobId);
            respond(ex, 200, "application/json", "{\"ok\":true}".getBytes());

        } else if ("log".equals(subAction)) {
            PluginEngine.RunResult result = analysisJobs.get(jobId);
            if (result != null && result.runDir != null) {
                File logFile = new File(result.runDir, "run.log");
                if (logFile.exists()) {
                    byte[] bytes = java.nio.file.Files.readAllBytes(logFile.toPath());
                    respond(ex, 200, "text/plain", bytes);
                    return;
                }
            }
            respond(ex, 404, "text/plain", "Log not found".getBytes());
        }
    }

    // GET /api/project/{id}/analysis/history?locus_id=...
    private void projectAnalysisHistory(HttpExchange ex, String projectId,
                                         String projectDir) throws IOException {
        String query = ex.getRequestURI().getQuery();
        String locusId = null;
        if (query != null) {
            for (String p : query.split("&")) {
                String[] kv = p.split("=", 2);
                if (kv.length == 2 && kv[0].equals("locus_id")) locusId = kv[1];
            }
        }

        if (locusId == null) {
            respond(ex, 400, "application/json",
                "{\"error\":\"locus_id query param required\"}".getBytes()); return;
        }

        File runsDir = new File(projectDir, "loci_analysis/" + locusId + "/runs");
        StringBuilder json = new StringBuilder("{\"runs\":[");
        boolean first = true;

        if (runsDir.isDirectory()) {
            File[] runDirs = runsDir.listFiles(File::isDirectory);
            if (runDirs != null) {
                Arrays.sort(runDirs, Comparator.comparingLong(File::lastModified).reversed());
                for (File rd : runDirs) {
                    File provFile = new File(rd, "provenance.json");
                    if (!provFile.exists()) continue;
                    if (!first) json.append(',');
                    first = false;
                    String prov = new String(java.nio.file.Files.readAllBytes(provFile.toPath()), "UTF-8");
                    // Add run status
                    boolean hasResult = new File(rd, "result.tsv").exists();
                    json.append("{\"run_dir\":\"").append(escJ(rd.getName())).append('"');
                    json.append(",\"has_result\":").append(hasResult);
                    json.append(",\"provenance\":").append(prov.trim());
                    json.append('}');
                }
            }
        }
        json.append("]}");
        respond(ex, 200, "application/json", json.toString().getBytes("UTF-8"));
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

    // GET /api/agent/config — current agent LLM settings (mode/base_url/api_key/model).
    // POST /api/agent/config — replace them (DECISIONS_PHASE4.md section 1).
    private void agentConfigEndpoint(HttpExchange ex) throws IOException {
        cors(ex); if (preflight(ex)) return;
        if ("POST".equalsIgnoreCase(ex.getRequestMethod())) {
            String body = new String(readAll(ex.getRequestBody()), "UTF-8");
            try {
                AgentConfig cfg = new AgentConfig();
                cfg.applyJson(body);
                cfg.save();
                respond(ex, 200, "application/json", "{\"ok\":true}".getBytes());
            } catch (Exception e) {
                respond(ex, 500, "application/json",
                    ("{\"error\":\"" + escJ(e.getMessage()) + "\"}").getBytes());
            }
        } else {
            respond(ex, 200, "application/json", AgentConfig.load().toJson().getBytes("UTF-8"));
        }
    }

    // POST /api/agent/pull-model — { "model": "llama3.2" } — proxies to Ollama's own POST /api/pull
    // (a real, stable, documented endpoint) so the settings panel can offer a plain "Download" button
    // instead of asking the user to open a terminal. Does not, and cannot, install Ollama itself — see
    // OllamaPuller's javadoc. Only one pull tracked at a time (a settings-panel affordance, not a queue).
    private void agentPullModel(HttpExchange ex) throws IOException {
        cors(ex); if (preflight(ex)) return;
        if (!"POST".equalsIgnoreCase(ex.getRequestMethod())) {
            respond(ex, 405, "application/json", "{\"error\":\"POST required\"}".getBytes()); return;
        }
        String body = new String(readAll(ex.getRequestBody()), "UTF-8");
        String model = extractStr(body, "model");
        if (model == null || model.isEmpty()) {
            respond(ex, 400, "application/json", "{\"error\":\"model is required\"}".getBytes()); return;
        }
        OllamaPuller.startPull("http://localhost:11434", model);
        respond(ex, 200, "application/json", "{\"started\":true}".getBytes());
    }

    // GET /api/agent/pull-progress — current status of the most recently started model download.
    private void agentPullProgress(HttpExchange ex) throws IOException {
        cors(ex); if (preflight(ex)) return;
        OllamaPuller.PullStatus s = OllamaPuller.current();
        long pct = s.total > 0 ? (100 * s.completed / s.total) : 0;
        String json = "{"
            + "\"status\":" + jsonStrStatic(s.status) + ","
            + "\"completed\":" + s.completed + ","
            + "\"total\":" + s.total + ","
            + "\"pct\":" + pct + ","
            + "\"done\":" + s.done + ","
            + "\"error\":" + (s.error == null ? "null" : jsonStrStatic(s.error))
            + "}";
        respond(ex, 200, "application/json", json.getBytes("UTF-8"));
    }

    // ── GitHub project save/extract (Phase 6) ───────────────────────────────

    // GET /api/github/token-status — reports only whether a credential currently resolves for
    // github.com; never the value itself (see GithubTokenStore's own javadoc for why).
    private void githubTokenStatus(HttpExchange ex) throws IOException {
        cors(ex); if (preflight(ex)) return;
        boolean configured = GithubTokenStore.isConfigured();
        respond(ex, 200, "application/json", ("{\"configured\":" + configured + "}").getBytes("UTF-8"));
    }

    // POST /api/github/token — { "token": "..." } — saves it via the user's own configured git
    // credential helper. DELETE /api/github/token — forgets it. Neither ever echoes the token back.
    private void githubToken(HttpExchange ex) throws IOException {
        cors(ex); if (preflight(ex)) return;
        if ("DELETE".equalsIgnoreCase(ex.getRequestMethod())) {
            GithubTokenStore.Result r = GithubTokenStore.forget();
            respond(ex, r.ok ? 200 : 500, "application/json",
                ("{\"ok\":" + r.ok + ",\"message\":" + jsonStrStatic(r.message) + "}").getBytes("UTF-8"));
            return;
        }
        if (!"POST".equalsIgnoreCase(ex.getRequestMethod())) {
            respond(ex, 405, "application/json", "{\"error\":\"POST or DELETE required\"}".getBytes()); return;
        }
        String body = new String(readAll(ex.getRequestBody()), "UTF-8");
        String token = extractStr(body, "token");
        GithubTokenStore.Result r = GithubTokenStore.save(token);
        respond(ex, r.ok ? 200 : 400, "application/json",
            ("{\"ok\":" + r.ok + ",\"message\":" + jsonStrStatic(r.message) + "}").getBytes("UTF-8"));
    }

    // POST /api/project/{id}/github/push — { "repo_url", "branch"?, "commit_message"?, "token"? } —
    // pushes this project's own directory. "token" here, if given, is a one-time credential
    // (DECISIONS_PHASE6.md section 2) never saved anywhere; omit it to use whatever credential is
    // already configured via /api/github/token.
    private void projectGithubPush(HttpExchange ex, String projectDir) throws IOException {
        if (!"POST".equalsIgnoreCase(ex.getRequestMethod())) {
            respond(ex, 405, "application/json", "{\"error\":\"POST required\"}".getBytes()); return;
        }
        String body = new String(readAll(ex.getRequestBody()), "UTF-8");
        String repoUrl = extractStr(body, "repo_url");
        String branch = extractStr(body, "branch");
        String commitMessage = extractStr(body, "commit_message");
        String token = extractStr(body, "token");
        GithubProjectSync.Result r = GithubProjectSync.push(new File(projectDir), repoUrl, branch, commitMessage, token);
        respond(ex, r.ok ? 200 : 400, "application/json",
            ("{\"ok\":" + r.ok + ",\"message\":" + jsonStrStatic(r.message) + "}").getBytes("UTF-8"));
    }

    // POST /api/github/import — { "repo_url", "branch"?, "project_id", "token"? } — extracts a
    // project from a GitHub repo into projects/{project_id}. Refuses to overwrite an existing,
    // unrelated directory (see GithubProjectSync.pullInto).
    private void githubImport(HttpExchange ex) throws IOException {
        cors(ex); if (preflight(ex)) return;
        if (!"POST".equalsIgnoreCase(ex.getRequestMethod())) {
            respond(ex, 405, "application/json", "{\"error\":\"POST required\"}".getBytes()); return;
        }
        String body = new String(readAll(ex.getRequestBody()), "UTF-8");
        String repoUrl = extractStr(body, "repo_url");
        String branch = extractStr(body, "branch");
        String projectId = extractStr(body, "project_id");
        String token = extractStr(body, "token");
        if (projectId == null || projectId.isEmpty()) {
            respond(ex, 400, "application/json", "{\"error\":\"project_id is required\"}".getBytes()); return;
        }
        File targetDir = new File("projects", projectId);
        GithubProjectSync.Result r = GithubProjectSync.pullInto(targetDir, repoUrl, branch, token);
        respond(ex, r.ok ? 200 : 400, "application/json",
            ("{\"ok\":" + r.ok + ",\"message\":" + jsonStrStatic(r.message) + "}").getBytes("UTF-8"));
    }

    // GET /api/agent/detect-local — probes common local LLM runtime ports (Ollama 11434, LM Studio
    // 1234) with a short timeout so the agent settings UI can offer a live model dropdown instead of
    // asking the user to type a model name blind. Absence of either is normal, not an error.
    private void agentDetectLocal(HttpExchange ex) throws IOException {
        cors(ex); if (preflight(ex)) return;
        StringBuilder json = new StringBuilder("{");
        json.append("\"ollama\":").append(detectOllama()).append(",");
        json.append("\"lmstudio\":").append(detectLmStudio());
        json.append("}");
        respond(ex, 200, "application/json", json.toString().getBytes("UTF-8"));
    }

    @SuppressWarnings("unchecked")
    private String detectOllama() {
        try {
            String body = shortGet("http://localhost:11434/api/tags");
            Map<String, Object> obj = MiniJson.asObject(MiniJson.parse(body));
            List<Object> models = MiniJson.asArray(obj.getOrDefault("models", new ArrayList<>()));
            List<String> names = new ArrayList<>();
            for (Object m : models) names.add(MiniJson.getStr(MiniJson.asObject(m), "name", ""));
            return "{\"available\":true,\"base_url\":\"http://localhost:11434/v1\",\"models\":"
                + MiniJson.encode(names) + "}";
        } catch (Exception e) {
            return "{\"available\":false}";
        }
    }

    @SuppressWarnings("unchecked")
    private String detectLmStudio() {
        try {
            String body = shortGet("http://localhost:1234/v1/models");
            Map<String, Object> obj = MiniJson.asObject(MiniJson.parse(body));
            List<Object> data = MiniJson.asArray(obj.getOrDefault("data", new ArrayList<>()));
            List<String> ids = new ArrayList<>();
            for (Object m : data) ids.add(MiniJson.getStr(MiniJson.asObject(m), "id", ""));
            return "{\"available\":true,\"base_url\":\"http://localhost:1234/v1\",\"models\":"
                + MiniJson.encode(ids) + "}";
        } catch (Exception e) {
            return "{\"available\":false}";
        }
    }

    private String shortGet(String url) throws IOException, InterruptedException {
        java.net.http.HttpClient client = java.net.http.HttpClient.newBuilder()
            .connectTimeout(java.time.Duration.ofMillis(700)).build();
        java.net.http.HttpRequest req = java.net.http.HttpRequest.newBuilder(URI.create(url))
            .timeout(java.time.Duration.ofMillis(1200)).GET().build();
        java.net.http.HttpResponse<String> resp = client.send(req,
            java.net.http.HttpResponse.BodyHandlers.ofString(java.nio.charset.StandardCharsets.UTF_8));
        if (resp.statusCode() / 100 != 2) throw new IOException("HTTP " + resp.statusCode());
        return resp.body();
    }

    // POST /api/agent/chat — { "message": "...", "history": [...] } -> { "reply", "steps", "history" }
    // Runs the tool-calling loop (AgentOrchestrator) against whichever LLM endpoint is configured via
    // /api/agent/config. Stateless server-side by design — see AgentOrchestrator's own javadoc.
    @SuppressWarnings("unchecked")
    private void agentChat(HttpExchange ex) throws IOException {
        cors(ex); if (preflight(ex)) return;
        if (!"POST".equalsIgnoreCase(ex.getRequestMethod())) {
            respond(ex, 405, "application/json", "{\"error\":\"POST required\"}".getBytes()); return;
        }
        String body = new String(readAll(ex.getRequestBody()), "UTF-8");
        try {
            Map<String, Object> req = MiniJson.asObject(MiniJson.parse(body));
            String message = MiniJson.getStr(req, "message", "");
            if (message.isEmpty()) {
                respond(ex, 400, "application/json", "{\"error\":\"message is required\"}".getBytes());
                return;
            }
            List<Map<String, Object>> history = new ArrayList<>();
            Object rawHistory = req.get("history");
            if (rawHistory != null) {
                for (Object m : MiniJson.asArray(rawHistory)) history.add(MiniJson.asObject(m));
            }

            AgentConfig cfg = AgentConfig.load();
            if (cfg.model == null || cfg.model.isEmpty()) {
                respond(ex, 400, "application/json",
                    "{\"error\":\"No model configured yet — open Agent Settings and pick a local or cloud model first.\"}".getBytes());
                return;
            }

            AgentOrchestrator.ChatOutcome outcome = AgentOrchestrator.chat(history, message, cfg, PORT);

            List<Object> stepsJson = new ArrayList<>();
            for (AgentOrchestrator.Step s : outcome.steps) {
                Map<String, Object> sj = new LinkedHashMap<>();
                sj.put("tool", s.tool);
                sj.put("arguments", s.argumentsJson);
                sj.put("summary", s.summary);
                sj.put("ok", s.ok);
                stepsJson.add(sj);
            }
            Map<String, Object> respObj = new LinkedHashMap<>();
            respObj.put("reply", outcome.reply);
            respObj.put("steps", stepsJson);
            respObj.put("history", outcome.messages);
            respond(ex, 200, "application/json", MiniJson.encode(respObj).getBytes("UTF-8"));
        } catch (Exception e) {
            String msg = e.getMessage() == null ? e.getClass().getSimpleName() : e.getMessage();
            respond(ex, 500, "application/json", ("{\"error\":" + jsonStrStatic(msg) + "}").getBytes());
        }
    }

    private static String jsonStrStatic(String s) {
        return MiniJson.encode(s);
    }

    // POST /api/shared-storage/test — validates the currently-saved shared-storage config (folder
    // reachability, or a git clone/pull), for the Resources & Settings "Test" button.
    private void sharedStorageTest(HttpExchange ex) throws IOException {
        cors(ex); if (preflight(ex)) return;
        if (!"POST".equalsIgnoreCase(ex.getRequestMethod())) {
            respond(ex, 405, "application/json", "{\"error\":\"POST required\"}".getBytes()); return;
        }
        GlobalConfig gc = GlobalConfig.load();
        rsid.SharedStorageResolver.Result r = rsid.SharedStorageResolver.test(gc.sharedStorage);
        String json = "{\"ok\":" + r.ok + ",\"message\":\"" + escJ(r.message) + "\"}";
        respond(ex, 200, "application/json", json.getBytes("UTF-8"));
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

    // GET /api/missing-rsids — aggregates unmatched SNPs (chr/pos/ea/nea) across every
    // project's rsid_recovery_results.csv, so the user can manually supply the rsID.
    private void missingRsids(HttpExchange ex) throws IOException {
        cors(ex); if (preflight(ex)) return;
        List<File> dirs = Main.discoverProjects(new File("projects"));
        StringBuilder json = new StringBuilder("{\"missing\":[");
        boolean first = true;
        for (File dir : dirs) {
            String projectId = dir.getName();
            File csvFile = new File(dir, "rsid_recovery_results.csv");
            if (!csvFile.exists()) continue;
            ProjectMetadata pm = ProjectMetadata.load(dir.getAbsolutePath());
            String projectName = pm != null && !pm.name.isEmpty() ? pm.name : projectId;
            try (BufferedReader br = new BufferedReader(new FileReader(csvFile))) {
                br.readLine(); // header: chr,pos,gwas_ea,gwas_nea,assigned_rsid,match_reason,pos_only_rsid,n_candidates_at_pos
                String line;
                while ((line = br.readLine()) != null) {
                    if (line.isEmpty()) continue;
                    String[] f = line.split(",", -1);
                    if (f.length < 6 || !f[4].isEmpty()) continue; // already has an rsid
                    if (!first) json.append(',');
                    first = false;
                    json.append('{')
                        .append("\"project_id\":\"").append(escJ(projectId)).append("\",")
                        .append("\"project_name\":\"").append(escJ(projectName)).append("\",")
                        .append("\"chr\":\"").append(escJ(f[0])).append("\",")
                        .append("\"pos\":").append(f[1]).append(',')
                        .append("\"ea\":\"").append(escJ(f[2])).append("\",")
                        .append("\"nea\":\"").append(escJ(f[3])).append("\",")
                        .append("\"match_reason\":\"").append(escJ(f[5])).append("\",")
                        .append("\"pos_only_rsid\":\"").append(f.length > 6 ? escJ(f[6]) : "").append('"')
                        .append('}');
                }
            } catch (IOException e) {
                System.err.println("[Server] Failed reading " + csvFile + ": " + e.getMessage());
            }
        }
        json.append("]}");
        respond(ex, 200, "application/json", json.toString().getBytes("UTF-8"));
    }

    // POST /api/submit-rsid — { "project_id", "chr", "pos", "rsid" }
    // Manually supplies an rsID for one SNP: patches the locus JSON(s) and marks the
    // row in rsid_recovery_results.csv as resolved so it drops out of /api/missing-rsids.
    private void submitRsid(HttpExchange ex) throws IOException {
        cors(ex); if (preflight(ex)) return;
        if (!"POST".equalsIgnoreCase(ex.getRequestMethod())) {
            respond(ex, 405, "text/plain", "Method Not Allowed".getBytes()); return;
        }
        String body = new String(readAll(ex.getRequestBody()), "UTF-8");
        String projectId = extractStr(body, "project_id");
        String chr = extractStr(body, "chr");
        String posStr = extractStr(body, "pos");
        String rsid = extractStr(body, "rsid");

        if (projectId == null || chr == null || posStr == null || rsid == null || rsid.trim().isEmpty()) {
            respond(ex, 400, "application/json",
                "{\"error\":\"missing project_id, chr, pos, or rsid\"}".getBytes()); return;
        }
        rsid = rsid.trim();

        File projDir = new File("projects", projectId);
        if (!projDir.isDirectory()) {
            respond(ex, 404, "application/json", "{\"error\":\"Project not found\"}".getBytes()); return;
        }

        try {
            Map<String, String> posToRsid = new HashMap<>();
            posToRsid.put(chr + ":" + posStr, rsid);
            int patched = RsidPipeline.patchLocusJsonsWithRsids(projDir.getAbsolutePath(), posToRsid);

            File csvFile = new File(projDir, "rsid_recovery_results.csv");
            if (csvFile.exists()) {
                List<String> lines = Files.readAllLines(csvFile.toPath());
                List<String> updated = new ArrayList<>();
                for (String line : lines) {
                    String[] f = line.split(",", -1);
                    if (f.length >= 6 && f[0].equals(chr) && f[1].equals(posStr) && f[4].isEmpty()) {
                        f[4] = rsid;
                        f[5] = "user_provided";
                        updated.add(String.join(",", f));
                    } else {
                        updated.add(line);
                    }
                }
                Files.write(csvFile.toPath(), updated);
            }

            String json = "{\"ok\":true,\"patched_files\":" + patched + "}";
            respond(ex, 200, "application/json", json.getBytes("UTF-8"));
        } catch (Exception e) {
            respond(ex, 500, "application/json",
                ("{\"error\":\"" + escJ(e.getMessage()) + "\"}").getBytes());
        }
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
        v = extractStr(body, "clump_p2");
        if (v != null && !v.isEmpty()) params.clumpP2 = Double.parseDouble(v);
        v = extractStr(body, "clump_r2");
        if (v != null && !v.isEmpty()) params.clumpR2 = Double.parseDouble(v);
        v = extractStr(body, "clump_kb");
        if (v != null && !v.isEmpty()) params.clumpKb = Integer.parseInt(v);

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

        export.GwasSchema schema = null;
        try {
            Config cfg = Config.loadFromProject(projectDir);
            schema = new export.GwasSchema();
            schema.gwasFile = cfg.gwasFile;
            schema.colChr = cfg.colChr; schema.colPos = cfg.colPos;
            schema.colEa = cfg.colEa; schema.colNea = cfg.colNea; schema.colPvalue = cfg.colPvalue;
            schema.colRsid = cfg.colRsid; schema.colVarid = cfg.colVarid;
            schema.colBeta = cfg.colBeta; schema.colOr = cfg.colOr; schema.colSe = cfg.colSe;
            schema.colN = cfg.colN; schema.colMaf = cfg.colMaf; schema.colInfo = cfg.colInfo;
        } catch (Exception e) {
            System.err.println("[LocalServer] Could not load config for extra GWAS columns: " + e.getMessage());
        }
        final export.GwasSchema finalSchema = schema;

        if ("locus".equals(scope)) {
            String sIdx = extractStr(body, "locus_index");
            int locusIdx = sIdx != null ? Integer.parseInt(sIdx) : 1;
            outputPath = new File(exportsDir, projectId + "_locus" + locusIdx + "_" + date + ".xlsx").getAbsolutePath();
            String finalPath = outputPath;
            new Thread(() -> ExcelExporter.exportLocus(projectDir, locusIdx, finalPath, ep, finalSchema),
                "excel-" + projectId).start();
        } else {
            outputPath = new File(exportsDir, projectId + "_" + date + ".xlsx").getAbsolutePath();
            String finalPath = outputPath;
            new Thread(() -> ExcelExporter.exportDataset(projectDir, finalPath, ep, finalSchema),
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

    // GET /api/export-projects-info — one-sheet cross-project summary table
    // (loci, SNPs, cases/controls/total, ancestry, ref panel, rsID status, etc.)
    private void exportProjectsInfo(HttpExchange ex) throws IOException {
        cors(ex); if (preflight(ex)) return;
        try {
            List<File> dirs = Main.discoverProjects(new File("projects"));

            ByteArrayOutputStream baos = new ByteArrayOutputStream();
            try (XlsxWriter xlsx = new XlsxWriter(baos)) {
                XlsxWriter.Sheet s = xlsx.addSheet("Projects Info").freezeHeader().autoFilter();
                s.addRow("Project ID", "Name", "Status", "Loci", "Total SNPs",
                    "Cases", "Controls", "Total N", "Ancestry", "Genome Build",
                    "Ref Panel Population", "Ref Panel Path", "RSID Present",
                    "RSID Recovery Status", "RSID Recovery Rate", "Last Processed",
                    "GWAS File");

                for (File dir : dirs) {
                    String id = dir.getName();
                    String projectDir = dir.getAbsolutePath();
                    ProjectMetadata pm = ProjectMetadata.load(projectDir);

                    Config cfg = null;
                    String status;
                    try {
                        cfg = Config.loadFromProject(projectDir);
                        ProjectMetadata.StaleReason reason =
                            ProjectMetadata.checkStaleness(projectDir, cfg);
                        status = reason == ProjectMetadata.StaleReason.NOT_STALE
                            ? "up_to_date" : "needs_reprocessing";
                    } catch (Exception e) {
                        status = "error";
                    }

                    boolean rsidPresent = (pm != null && pm.rsidColumnPresent)
                        || (cfg != null && cfg.colRsid != null && !cfg.colRsid.isEmpty());

                    s.addRow(
                        id,
                        pm != null && !pm.name.isEmpty() ? pm.name : id,
                        status,
                        pm != null ? pm.lociCount : 0,
                        pm != null ? pm.totalSnps : 0,
                        cfg != null ? cfg.nCases : 0,
                        cfg != null ? cfg.nControls : 0,
                        cfg != null ? cfg.sampleN : 0,
                        cfg != null ? cfg.ancestry : "",
                        cfg != null ? cfg.genomeBuild : "",
                        cfg != null ? cfg.refPanelPopulation : "",
                        cfg != null ? cfg.refPanelPath : "",
                        rsidPresent,
                        pm != null ? pm.rsidRecoveryStatus : "not_started",
                        pm != null ? pm.rsidRecoveryRate : "",
                        pm != null ? pm.lastProcessed : "",
                        cfg != null ? cfg.gwasFile : ""
                    );
                }
                xlsx.finish();
            }

            byte[] bytes = baos.toByteArray();
            ex.getResponseHeaders().set("Content-Type",
                "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet");
            ex.getResponseHeaders().set("Content-Disposition",
                "attachment; filename=\"lynxgwas_projects_info.xlsx\"");
            ex.sendResponseHeaders(200, bytes.length);
            try (OutputStream os = ex.getResponseBody()) { os.write(bytes); }
        } catch (Exception e) {
            respond(ex, 500, "application/json",
                ("{\"error\":\"" + escJ(e.getMessage()) + "\"}").getBytes());
        }
    }

    // ══════════════════════════════════════════════════════════════════════
    //  LOCUS MATRIX (cross-dataset locus significance grid)
    // ══════════════════════════════════════════════════════════════════════

    // POST /api/locus-matrix-run — {"project_ids":[...], "ref_panel_id":"..."}
    private void locusMatrixRun(HttpExchange ex) throws IOException {
        cors(ex); if (preflight(ex)) return;
        if (!"POST".equalsIgnoreCase(ex.getRequestMethod())) {
            respond(ex, 405, "application/json", "{\"error\":\"POST required\"}".getBytes()); return;
        }
        String body = new String(readAll(ex.getRequestBody()), "UTF-8");
        List<String> projectIds = extractStringArray(body, "project_ids");
        String refPanelId = extractStr(body, "ref_panel_id");
        String requestedName = extractStr(body, "name");
        String requestedKind = extractStr(body, "kind");

        if (projectIds.isEmpty()) {
            respond(ex, 400, "application/json", "{\"error\":\"project_ids required\"}".getBytes()); return;
        }
        if (refPanelId == null || refPanelId.isEmpty()) {
            respond(ex, 400, "application/json", "{\"error\":\"ref_panel_id required\"}".getBytes()); return;
        }

        GlobalConfig gc = GlobalConfig.load();
        GlobalConfig.RefPanel panel = gc.findPanel(refPanelId);
        if (panel == null || !panel.isValid()) {
            respond(ex, 400, "application/json",
                ("{\"error\":\"Reference panel not found or invalid: " + escJ(refPanelId) + "\"}").getBytes());
            return;
        }

        List<Config> configs = new ArrayList<>();
        for (String pid : projectIds) {
            File dir = new File("projects/" + pid);
            if (!dir.isDirectory()) {
                respond(ex, 400, "application/json",
                    ("{\"error\":\"Project not found: " + escJ(pid) + "\"}").getBytes());
                return;
            }
            try {
                configs.add(Config.loadFromProject(dir.getAbsolutePath()));
            } catch (Exception e) {
                respond(ex, 400, "application/json",
                    ("{\"error\":\"Failed to load project '" + escJ(pid) + "': " + escJ(e.getMessage()) + "\"}").getBytes());
                return;
            }
        }

        List<String> datasetNames = new ArrayList<>();
        for (String pid : projectIds) {
            ProjectMetadata pm = ProjectMetadata.load(new File("projects/" + pid).getAbsolutePath());
            datasetNames.add(pm != null && !pm.name.isEmpty() ? pm.name : pid);
        }

        String jobId = UUID.randomUUID().toString();
        MultiLocusProgress progress = new MultiLocusProgress();
        locusMatrixJobProgress.put(jobId, progress);

        LocusMatrixJobMeta meta = new LocusMatrixJobMeta();
        meta.jobId = jobId;
        meta.name = (requestedName != null && !requestedName.trim().isEmpty())
            ? requestedName.trim() : String.join(" + ", datasetNames);
        meta.projectIds = projectIds;
        meta.datasetNames = datasetNames;
        meta.refPanelId = panel.id;
        meta.refPanelLabel = panel.label;
        meta.createdAt = java.time.Instant.now().toString();
        if ("serpent".equals(requestedKind)) meta.kind = "serpent";
        locusMatrixJobMeta.put(jobId, meta);
        locusMatrixJobOrder.add(jobId);
        saveLocusMatrixJobMeta(meta, "running", null);

        final GlobalConfig.RefPanel finalPanel = panel;
        final List<String> finalProjectIds = projectIds;
        final String finalJobName = meta.name;
        new Thread(() -> {
            try {
                File jobDir = new File("output/multi_locus/" + jobId);
                jobDir.mkdirs();
                File mergedFile = new File(jobDir, "merged.tsv");

                MultiLocusMerger.merge(configs, 5e-3, mergedFile, progress);

                progress.phase = "identify";
                LociIdentifier.Params params = new LociIdentifier.Params();
                params.refPanelPath = finalPanel.plinkPath;
                List<LociIdentifier.IdentifiedLocus> loci = LociIdentifier.identify(
                    mergedFile.getAbsolutePath(), "chrom", "pos", "p", "", params, progress.identifyProgress);
                progress.lociFound = loci.size();

                GffParser gff = GffParser.parse(configs.get(0));

                MultiLocusResult result = MultiLocusScanner.scan(configs, finalProjectIds, loci, gff, progress);
                result.name = finalJobName;
                result.refPanelId = finalPanel.id;
                result.refPanelLabel = finalPanel.label;
                result.createdAt = java.time.Instant.now().toString();

                locusMatrixJobs.put(jobId, result);
                saveLocusMatrixJobResult(jobId, result);
                saveLocusMatrixJobMeta(meta, "done", null);
                progress.phase = "done";
                progress.done = true;
            } catch (Exception e) {
                progress.phase = "error";
                progress.error = e.getMessage();
                progress.done = true;
                saveLocusMatrixJobMeta(meta, "error", e.getMessage());
                System.err.printf("[LocusMatrix] Job '%s' failed: %s%n", jobId, e.getMessage());
                e.printStackTrace(System.err);
            }
        }, "locus-matrix-" + jobId).start();

        respond(ex, 202, "application/json", ("{\"job_id\":\"" + jobId + "\",\"status\":\"started\"}").getBytes());
    }

    // GET /api/locus-matrix-progress?job=<jobId>
    private void locusMatrixProgressEndpoint(HttpExchange ex) throws IOException {
        cors(ex); if (preflight(ex)) return;
        String jobId = queryParam(ex, "job");
        MultiLocusProgress p = jobId != null ? locusMatrixJobProgress.get(jobId) : null;
        if (p == null) {
            respond(ex, 200, "application/json",
                "{\"pct\":0,\"phase\":\"merge\",\"dataset_index\":0,\"dataset_total\":0,\"current_dataset\":\"\",\"loci_found\":0,\"done\":false,\"error\":null}".getBytes());
            return;
        }
        respond(ex, 200, "application/json", p.toJson().getBytes("UTF-8"));
    }

    // GET /api/locus-matrix-jobs — list all Locus Matrix runs this session, most recent first
    private void locusMatrixJobsList(HttpExchange ex) throws IOException {
        cors(ex); if (preflight(ex)) return;
        StringBuilder j = new StringBuilder("[");
        for (int i = locusMatrixJobOrder.size() - 1; i >= 0; i--) {
            String jobId = locusMatrixJobOrder.get(i);
            LocusMatrixJobMeta meta = locusMatrixJobMeta.get(jobId);
            if (meta == null) continue;
            MultiLocusProgress progress = locusMatrixJobProgress.get(jobId);
            MultiLocusResult result = locusMatrixJobs.get(jobId);

            String status = "running";
            if (progress != null && progress.done) status = progress.error != null ? "error" : "done";

            if (j.length() > 1) j.append(",");
            j.append("{");
            j.append("\"job_id\":\"").append(escJ(jobId)).append("\",");
            j.append("\"name\":\"").append(escJ(meta.name)).append("\",");
            j.append("\"datasets\":\"").append(escJ(String.join(", ", meta.datasetNames))).append("\",");
            j.append("\"ref_panel_label\":\"").append(escJ(meta.refPanelLabel)).append("\",");
            j.append("\"created_at\":\"").append(escJ(meta.createdAt)).append("\",");
            j.append("\"kind\":\"").append(escJ(meta.kind)).append("\",");
            j.append("\"status\":\"").append(status).append("\",");
            j.append("\"loci_found\":").append(result != null ? result.loci.size() : 0).append(",");
            j.append("\"error\":").append(progress != null && progress.error != null
                ? "\"" + escJ(progress.error) + "\"" : "null");
            j.append("}");
        }
        j.append("]");
        respond(ex, 200, "application/json", j.toString().getBytes("UTF-8"));
    }

    // POST /api/locus-matrix-delete — {"job_id":"..."}
    private void locusMatrixDelete(HttpExchange ex) throws IOException {
        cors(ex); if (preflight(ex)) return;
        if (!"POST".equalsIgnoreCase(ex.getRequestMethod())) {
            respond(ex, 405, "application/json", "{\"error\":\"POST required\"}".getBytes()); return;
        }
        String body = new String(readAll(ex.getRequestBody()), "UTF-8");
        String jobId = extractStr(body, "job_id");
        if (jobId == null || jobId.isEmpty()) {
            respond(ex, 400, "application/json", "{\"error\":\"job_id required\"}".getBytes()); return;
        }

        locusMatrixJobOrder.remove(jobId);
        locusMatrixJobMeta.remove(jobId);
        locusMatrixJobProgress.remove(jobId);
        locusMatrixJobs.remove(jobId);
        serpentPlotCache.remove(jobId);
        if (jobId.matches("[A-Za-z0-9-]+")) purgeMvsusieRunsOf(jobId);

        try {
            File jobDir = new File(MULTI_LOCUS_DIR, jobId);
            if (jobDir.isDirectory() && jobId.matches("[A-Za-z0-9-]+")) {
                File[] files = jobDir.listFiles();
                if (files != null) for (File f : files) f.delete();
                jobDir.delete();
            }
        } catch (Exception ignored) {}

        respond(ex, 200, "application/json", "{\"ok\":true}".getBytes());
    }

    // GET /api/locus-matrix-result?job=<jobId>
    private void locusMatrixResult(HttpExchange ex) throws IOException {
        cors(ex); if (preflight(ex)) return;
        String jobId = queryParam(ex, "job");
        MultiLocusResult result = jobId != null ? locusMatrixJobs.get(jobId) : null;
        if (result == null) {
            respond(ex, 404, "application/json", "{\"error\":\"Job not found or not complete\"}".getBytes());
            return;
        }
        respond(ex, 200, "application/json", result.toJson().getBytes("UTF-8"));
    }

    // GET /api/gene-constellation?job=<jobId>&threshold=<p>  — derived gene-level view of an
    // already-completed Locus Matrix job (reuses the same locusMatrixJobs lookup, no new pipeline run).
    private void geneConstellation(HttpExchange ex) throws IOException {
        cors(ex); if (preflight(ex)) return;
        String jobId = queryParam(ex, "job");
        MultiLocusResult result = jobId != null ? locusMatrixJobs.get(jobId) : null;
        if (result == null) {
            respond(ex, 404, "application/json", "{\"error\":\"Job not found or not complete\"}".getBytes());
            return;
        }
        double threshold = GeneConstellationBuilder.DEFAULT_THRESHOLD;
        String thresholdParam = queryParam(ex, "threshold");
        if (thresholdParam != null && !thresholdParam.isEmpty()) {
            try { threshold = Double.parseDouble(thresholdParam); } catch (NumberFormatException ignored) {}
        }
        int minEdgeCount = GeneConstellationBuilder.DEFAULT_MIN_EDGE_COUNT;
        String minEdgeParam = queryParam(ex, "min_edge_count");
        if (minEdgeParam != null && !minEdgeParam.isEmpty()) {
            try { minEdgeCount = Math.max(1, Integer.parseInt(minEdgeParam)); } catch (NumberFormatException ignored) {}
        }
        int maxEdges = GeneConstellationBuilder.DEFAULT_MAX_EDGES;
        String maxEdgesParam = queryParam(ex, "max_edges");
        if (maxEdgesParam != null && !maxEdgesParam.isEmpty()) {
            try { maxEdges = Math.max(0, Integer.parseInt(maxEdgesParam)); } catch (NumberFormatException ignored) {}
        }
        try {
            GeneConstellationResult gcr = GeneConstellationBuilder.build(result, threshold, minEdgeCount, maxEdges);
            respond(ex, 200, "application/json", gcr.toJson().getBytes("UTF-8"));
        } catch (Exception e) {
            respond(ex, 500, "application/json", ("{\"error\":\"" + escJ(e.getMessage()) + "\"}").getBytes());
        }
    }

    // ── Persisted Locus Matrix runs (Gene Constellation / Locus Serpent Plot) ─────────────────────
    // Each run lives in output/multi_locus/<jobId>/: merged.tsv (pipeline input), job.json (name,
    // kind, datasets, ref panel, created, status) and, once finished, result.json
    // (MultiLocusResult.toJson()). Runs are reloaded at startup so they survive a restart; a run that
    // was still going when LYNXgwas stopped comes back as "interrupted" instead of vanishing.
    private static final String MULTI_LOCUS_DIR = "output/multi_locus";

    private void saveLocusMatrixJobMeta(LocusMatrixJobMeta meta, String status, String error) {
        try {
            Map<String, Object> m = new LinkedHashMap<>();
            m.put("job_id", meta.jobId);
            m.put("name", meta.name);
            m.put("kind", meta.kind);
            m.put("project_ids", new ArrayList<Object>(meta.projectIds));
            m.put("dataset_names", new ArrayList<Object>(meta.datasetNames));
            m.put("ref_panel_id", meta.refPanelId);
            m.put("ref_panel_label", meta.refPanelLabel);
            m.put("created_at", meta.createdAt);
            m.put("status", status);
            m.put("error", error);
            writeAtomically(new File(new File(MULTI_LOCUS_DIR, meta.jobId), "job.json"), MiniJson.encode(m));
        } catch (Exception e) {
            System.err.printf("[LocusMatrix] Could not save job.json for '%s': %s%n", meta.jobId, e.getMessage());
        }
    }

    private void saveLocusMatrixJobResult(String jobId, MultiLocusResult result) {
        try {
            writeAtomically(new File(new File(MULTI_LOCUS_DIR, jobId), "result.json"), result.toJson());
        } catch (Exception e) {
            System.err.printf("[LocusMatrix] Could not save result.json for '%s': %s%n", jobId, e.getMessage());
        }
    }

    /** Write to a temp file then rename, so a crash mid-write never leaves a truncated file. */
    private static void writeAtomically(File target, String content) throws IOException {
        target.getParentFile().mkdirs();
        File tmp = new File(target.getParentFile(), target.getName() + ".tmp");
        Files.write(tmp.toPath(), content.getBytes(java.nio.charset.StandardCharsets.UTF_8));
        try {
            Files.move(tmp.toPath(), target.toPath(), java.nio.file.StandardCopyOption.REPLACE_EXISTING,
                       java.nio.file.StandardCopyOption.ATOMIC_MOVE);
        } catch (java.nio.file.AtomicMoveNotSupportedException e) {
            Files.move(tmp.toPath(), target.toPath(), java.nio.file.StandardCopyOption.REPLACE_EXISTING);
        }
    }

    private void loadPersistedLocusMatrixJobs() {
        File root = new File(MULTI_LOCUS_DIR);
        File[] dirs = root.listFiles(File::isDirectory);
        if (dirs == null) return;
        List<LocusMatrixJobMeta> loaded = new ArrayList<>();
        int ok = 0, failed = 0;
        for (File dir : dirs) {
            File metaFile = new File(dir, "job.json");
            if (!metaFile.isFile()) continue; // runs from before persistence existed: nothing to restore
            try {
                Map<String, Object> m = MiniJson.asObject(MiniJson.parse(
                    new String(Files.readAllBytes(metaFile.toPath()), java.nio.charset.StandardCharsets.UTF_8)));
                LocusMatrixJobMeta meta = new LocusMatrixJobMeta();
                meta.jobId = MiniJson.getStr(m, "job_id", dir.getName());
                meta.name = MiniJson.getStr(m, "name", meta.jobId);
                meta.kind = MiniJson.getStr(m, "kind", "constellation");
                meta.projectIds = new ArrayList<>();
                meta.datasetNames = new ArrayList<>();
                if (m.get("project_ids") instanceof List) for (Object x : MiniJson.asArray(m.get("project_ids"))) meta.projectIds.add(String.valueOf(x));
                if (m.get("dataset_names") instanceof List) for (Object x : MiniJson.asArray(m.get("dataset_names"))) meta.datasetNames.add(String.valueOf(x));
                meta.refPanelId = MiniJson.getStr(m, "ref_panel_id", "");
                meta.refPanelLabel = MiniJson.getStr(m, "ref_panel_label", "");
                meta.createdAt = MiniJson.getStr(m, "created_at", "");
                String status = MiniJson.getStr(m, "status", "");

                MultiLocusProgress progress = new MultiLocusProgress();
                progress.done = true;
                File resultFile = new File(dir, "result.json");
                if ("done".equals(status) && resultFile.isFile()) {
                    MultiLocusResult result = MultiLocusResult.fromJson(
                        new String(Files.readAllBytes(resultFile.toPath()), java.nio.charset.StandardCharsets.UTF_8));
                    locusMatrixJobs.put(meta.jobId, result);
                    progress.phase = "done";
                    progress.lociFound = result.loci.size();
                    ok++;
                } else {
                    progress.phase = "error";
                    progress.error = "done".equals(status) ? "Saved result is missing (result.json)."
                        : "error".equals(status) ? MiniJson.getStr(m, "error", "Run failed.")
                        : "Interrupted: LYNXgwas stopped before this run finished. Delete it and run again.";
                    failed++;
                }
                locusMatrixJobMeta.put(meta.jobId, meta);
                locusMatrixJobProgress.put(meta.jobId, progress);
                loaded.add(meta);
            } catch (Exception e) {
                failed++;
                System.err.printf("[LocusMatrix] Could not restore run '%s': %s%n", dir.getName(), e.getMessage());
            }
        }
        loaded.sort(Comparator.comparing(x -> x.createdAt == null ? "" : x.createdAt));
        for (LocusMatrixJobMeta meta : loaded) locusMatrixJobOrder.add(meta.jobId);
        if (!loaded.isEmpty())
            System.out.printf("Restored %d saved cross-dataset run(s) (%d complete, %d failed/interrupted)%n", loaded.size(), ok, failed);
    }

    // GET /api/serpent-plot?job=<jobId> — Locus Serpent Plot view + category analysis, derived from an
    // already-completed Locus Matrix job (same job store as Gene Constellation; no pipeline re-run).
    private void serpentPlot(HttpExchange ex) throws IOException {
        cors(ex); if (preflight(ex)) return;
        String jobId = queryParam(ex, "job");
        MultiLocusResult result = jobId != null ? locusMatrixJobs.get(jobId) : null;
        if (result == null) {
            respond(ex, 404, "application/json", "{\"error\":\"Job not found or not complete\"}".getBytes());
            return;
        }
        warmGeneAnnotation(jobId); // so the first click on a serpent does not wait for the GFF parse
        try {
            String json = serpentPlotCache.computeIfAbsent(jobId,
                k -> SerpentPlotBuilder.toJson(SerpentPlotBuilder.build(result)));
            respond(ex, 200, "application/json", json.getBytes("UTF-8"));
        } catch (Exception e) {
            respond(ex, 500, "application/json", ("{\"error\":\"" + escJ(e.getMessage()) + "\"}").getBytes());
        }
    }

    // GET /api/serpent-locus?job=<jobId>&locus=<index> — table shown when the user clicks a serpent:
    // each dataset's lead SNP (position, alleles, OR + 95% CI, p) with its nearest gene.
    private void serpentLocus(HttpExchange ex) throws IOException {
        cors(ex); if (preflight(ex)) return;
        String jobId = queryParam(ex, "job");
        MultiLocusResult result = jobId != null ? locusMatrixJobs.get(jobId) : null;
        if (result == null) {
            respond(ex, 404, "application/json", "{\"error\":\"Job not found or not complete\"}".getBytes());
            return;
        }
        int locus;
        try { locus = Integer.parseInt(queryParam(ex, "locus")); }
        catch (Exception e) { respond(ex, 400, "application/json", "{\"error\":\"locus index required\"}".getBytes()); return; }
        try {
            GffParser gff = geneAnnotationFor(jobId);
            SerpentPlotBuilder.GeneLookup lookup = gff == null ? null : (chr, pos) -> nearestGene(gff, chr, pos);
            String json = SerpentPlotBuilder.locusDetailJson(result, locus, lookup);
            if (json == null) { respond(ex, 404, "application/json", "{\"error\":\"Locus not found\"}".getBytes()); return; }
            respond(ex, 200, "application/json", json.getBytes("UTF-8"));
        } catch (Exception e) {
            respond(ex, 500, "application/json", ("{\"error\":\"" + escJ(e.getMessage()) + "\"}").getBytes());
        }
    }

    /** Parsed gene annotation per GFF3 path, loaded once (the GENCODE file is ~1.3 GB). */
    private static final Map<String, java.util.concurrent.Future<GffParser>> gffByPath = new ConcurrentHashMap<>();

    /** The GFF3 of the first project in the run (all projects in a run share one annotation in practice). */
    private String gffPathForJob(String jobId) {
        LocusMatrixJobMeta meta = locusMatrixJobMeta.get(jobId);
        if (meta == null || meta.projectIds == null) return null;
        for (String pid : meta.projectIds) {
            try {
                Config c = Config.loadFromProject(new File("projects/" + pid).getAbsolutePath());
                if (c.gff3File != null && new File(c.gff3File).isFile()) return c.gff3File;
            } catch (Exception ignored) {}
        }
        return null;
    }

    private java.util.concurrent.Future<GffParser> gffFuture(String path) {
        return gffByPath.computeIfAbsent(path, k -> {
            java.util.concurrent.FutureTask<GffParser> task = new java.util.concurrent.FutureTask<>(() -> {
                Config c = new Config();
                c.gff3File = k;
                return GffParser.parse(c);
            });
            Thread t = new Thread(task, "gff-load"); t.setDaemon(true); t.start();
            return task;
        });
    }

    private void warmGeneAnnotation(String jobId) {
        try { String path = gffPathForJob(jobId); if (path != null) gffFuture(path); } catch (Exception ignored) {}
    }

    private GffParser geneAnnotationFor(String jobId) {
        String path = gffPathForJob(jobId);
        if (path == null) return null;
        try { return gffFuture(path).get(); }
        catch (Exception e) {
            gffByPath.remove(path); // let a later request retry
            System.err.printf("[Serpent] Could not load gene annotation %s: %s%n", path, e.getMessage());
            return null;
        }
    }

    /** Nearest protein-coding gene within 1 Mb of pos (falling back to any annotated gene), as {name, distance}. */
    private static Object[] nearestGene(GffParser gff, String chr, long pos) {
        String c = chr == null ? "" : (chr.startsWith("chr") ? chr.substring(3) : chr);
        List<Gene> near = gff.overlapping(c, Math.max(0, pos - 1_000_000), pos + 1_000_000);
        Gene best = null, bestAny = null;
        for (Gene g : near) {
            if (bestAny == null || g.distanceTo(pos) < bestAny.distanceTo(pos)) bestAny = g;
            Transcript t = g.canonical();
            boolean coding = t != null && "protein_coding".equals(t.transcriptType);
            if (coding && (best == null || g.distanceTo(pos) < best.distanceTo(pos))) best = g;
        }
        Gene g = best != null ? best : bestAny;
        return g == null ? null : new Object[]{ g.geneName, g.distanceTo(pos) };
    }

    // GET /api/locus-matrix-export?job=<jobId>
    private void locusMatrixExport(HttpExchange ex) throws IOException {
        cors(ex); if (preflight(ex)) return;
        String jobId = queryParam(ex, "job");
        MultiLocusResult result = jobId != null ? locusMatrixJobs.get(jobId) : null;
        if (result == null) {
            respond(ex, 404, "application/json", "{\"error\":\"Job not found or not complete\"}".getBytes());
            return;
        }
        try {
            ByteArrayOutputStream baos = new ByteArrayOutputStream();
            MultiLocusExcelWriter.write(result, baos);
            byte[] bytes = baos.toByteArray();
            ex.getResponseHeaders().set("Content-Type",
                "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet");
            ex.getResponseHeaders().set("Content-Disposition",
                "attachment; filename=\"locus_matrix_" + jobId + ".xlsx\"");
            ex.sendResponseHeaders(200, bytes.length);
            try (OutputStream os = ex.getResponseBody()) { os.write(bytes); }
        } catch (Exception e) {
            respond(ex, 500, "application/json", ("{\"error\":\"" + escJ(e.getMessage()) + "\"}").getBytes());
        }
    }

    // GET /api/gene-constellation-export?job=<jobId>&threshold=<p> — gene-level Excel export,
    // mirrors locusMatrixExport's structure exactly (same job lookup, same response headers).
    private void geneConstellationExport(HttpExchange ex) throws IOException {
        cors(ex); if (preflight(ex)) return;
        String jobId = queryParam(ex, "job");
        MultiLocusResult result = jobId != null ? locusMatrixJobs.get(jobId) : null;
        if (result == null) {
            respond(ex, 404, "application/json", "{\"error\":\"Job not found or not complete\"}".getBytes());
            return;
        }
        double threshold = GeneConstellationBuilder.DEFAULT_THRESHOLD;
        String thresholdParam = queryParam(ex, "threshold");
        if (thresholdParam != null && !thresholdParam.isEmpty()) {
            try { threshold = Double.parseDouble(thresholdParam); } catch (NumberFormatException ignored) {}
        }
        try {
            GeneConstellationResult gcr = GeneConstellationBuilder.build(result, threshold);
            ByteArrayOutputStream baos = new ByteArrayOutputStream();
            GeneConstellationExcelWriter.write(gcr, baos);
            byte[] bytes = baos.toByteArray();
            ex.getResponseHeaders().set("Content-Type",
                "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet");
            ex.getResponseHeaders().set("Content-Disposition",
                "attachment; filename=\"gene_constellation_" + jobId + ".xlsx\"");
            ex.sendResponseHeaders(200, bytes.length);
            try (OutputStream os = ex.getResponseBody()) { os.write(bytes); }
        } catch (Exception e) {
            respond(ex, 500, "application/json", ("{\"error\":\"" + escJ(e.getMessage()) + "\"}").getBytes());
        }
    }

    // ══════════════════════════════════════════════════════════════════════
    //  CROSS-PROJECT GENE/SNP/POSITION SEARCH
    // ══════════════════════════════════════════════════════════════════════

    // GET /api/search?type=gene|snp|position&q=<term>&chr=<c>&pos=<p>&build=GRCh37|GRCh38&threshold=<p>
    private void globalSearch(HttpExchange ex) throws IOException {
        cors(ex); if (preflight(ex)) return;
        String type = queryParam(ex, "type");
        if (type == null || type.isEmpty()) {
            respond(ex, 400, "application/json", "{\"error\":\"type is required (gene|snp|position)\"}".getBytes());
            return;
        }
        double threshold = GlobalSearchIndex.DEFAULT_THRESHOLD;
        String thresholdParam = queryParam(ex, "threshold");
        if (thresholdParam != null && !thresholdParam.isEmpty()) {
            try { threshold = Double.parseDouble(thresholdParam); } catch (NumberFormatException ignored) {}
        }

        GlobalSearchIndex idx = GlobalSearchIndex.getOrBuild(new File("projects"), threshold);
        GlobalSearchIndex.Result result;
        String q = queryParam(ex, "q");

        switch (type) {
            case "gene":
                result = idx.searchGene(q);
                break;
            case "snp":
                result = idx.searchSnp(q);
                break;
            case "position": {
                String chr = queryParam(ex, "chr");
                String posParam = queryParam(ex, "pos");
                String build = queryParam(ex, "build");
                long pos;
                try {
                    pos = Long.parseLong(posParam);
                } catch (Exception e) {
                    respond(ex, 400, "application/json", "{\"error\":\"pos must be an integer\"}".getBytes());
                    return;
                }
                if (chr == null || chr.isEmpty()) {
                    respond(ex, 400, "application/json", "{\"error\":\"chr is required for position search\"}".getBytes());
                    return;
                }
                result = idx.searchPosition(chr, pos, build);
                break;
            }
            default:
                respond(ex, 400, "application/json", ("{\"error\":\"unknown type: " + escJ(type) + "\"}").getBytes());
                return;
        }

        StringBuilder j = new StringBuilder();
        j.append('{');
        j.append("\"type\":\"").append(escJ(type)).append("\",");
        j.append("\"q\":\"").append(escJ(q == null ? "" : q)).append("\",");
        j.append("\"threshold\":").append(threshold).append(',');
        j.append("\"count\":").append(result.hits.size()).append(',');
        j.append("\"hits\":[");
        for (int i = 0; i < result.hits.size(); i++) {
            if (i > 0) j.append(',');
            j.append(result.hits.get(i).toJson());
        }
        j.append("],");
        j.append("\"not_compared\":[");
        for (int i = 0; i < result.notCompared.size(); i++) {
            if (i > 0) j.append(',');
            j.append(result.notCompared.get(i).toJson());
        }
        j.append("],");
        j.append("\"index_summary\":").append(idx.summaryJson());
        j.append('}');
        respond(ex, 200, "application/json", j.toString().getBytes("UTF-8"));
    }

    // POST /api/search/rebuild-index — forces a fresh build, returns the same summary stats logged
    // to stdout during the build (project/loci/SNP counts, build time).
    private void searchRebuildIndex(HttpExchange ex) throws IOException {
        cors(ex); if (preflight(ex)) return;
        if (!"POST".equalsIgnoreCase(ex.getRequestMethod())) {
            respond(ex, 405, "application/json", "{\"error\":\"POST required\"}".getBytes());
            return;
        }
        double threshold = GlobalSearchIndex.DEFAULT_THRESHOLD;
        String thresholdParam = queryParam(ex, "threshold");
        if (thresholdParam != null && !thresholdParam.isEmpty()) {
            try { threshold = Double.parseDouble(thresholdParam); } catch (NumberFormatException ignored) {}
        }
        GlobalSearchIndex idx = GlobalSearchIndex.forceRebuild(new File("projects"), threshold);
        respond(ex, 200, "application/json", idx.summaryJson().getBytes("UTF-8"));
    }

    // ══════════════════════════════════════════════════════════════════════
    //  EVIDENCE (Feature A) — per-project gene-level evidence tables
    //  Stored at projects/<id>/evidence/<name>.tsv + <name>.meta.json.
    //  No multipart parsing in this codebase (checked: every existing upload-shaped endpoint
    //  reads a raw request body), so uploads are POST with the raw file body plus
    //  ?name=&kind=&gene_col= query params, per the project's established simple-HTTP-handling
    //  convention.
    // ══════════════════════════════════════════════════════════════════════

    private static String sanitizeEvidenceName(String name) {
        return name == null ? "" : name.trim().replaceAll("[^A-Za-z0-9_.-]", "_");
    }

    /** Minimal CSV/TSV line splitter: tab-delimited lines split on tab as-is; comma-delimited
     *  lines get basic double-quote handling (quoted commas don't split, "" is an escaped quote). */
    private static String[] splitDelimited(String line, String delim) {
        if ("\t".equals(delim)) return line.split("\t", -1);
        List<String> out = new ArrayList<>();
        StringBuilder cur = new StringBuilder();
        boolean inQuotes = false;
        for (int i = 0; i < line.length(); i++) {
            char c = line.charAt(i);
            if (inQuotes) {
                if (c == '"') {
                    if (i + 1 < line.length() && line.charAt(i + 1) == '"') { cur.append('"'); i++; }
                    else inQuotes = false;
                } else cur.append(c);
            } else {
                if (c == '"') inQuotes = true;
                else if (c == ',') { out.add(cur.toString()); cur.setLength(0); }
                else cur.append(c);
            }
        }
        out.add(cur.toString());
        return out.toArray(new String[0]);
    }

    // GET /api/project/{id}/evidence — list attached evidence tables (each meta.json is already
    // valid JSON in the exact shape we want to return, so this just concatenates them into an array).
    private void projectEvidenceList(HttpExchange ex, String projectDir) throws IOException {
        File evDir = new File(projectDir, "evidence");
        File[] metas = evDir.isDirectory() ? evDir.listFiles((d, n) -> n.endsWith(".meta.json")) : null;
        StringBuilder j = new StringBuilder("[");
        if (metas != null) {
            Arrays.sort(metas, Comparator.comparing(File::getName));
            for (int i = 0; i < metas.length; i++) {
                if (i > 0) j.append(",");
                j.append(new String(Files.readAllBytes(metas[i].toPath()), "UTF-8"));
            }
        }
        j.append("]");
        respond(ex, 200, "application/json", j.toString().getBytes("UTF-8"));
    }

    // POST /api/project/{id}/evidence?name=<name>&kind=<ppi|expression|other>&gene_col=<col>
    // Body = raw CSV or TSV file content (header row required; delimiter auto-detected: a tab
    // in the header line means TSV, else comma).
    private void projectEvidenceUpload(HttpExchange ex, String projectId, String projectDir) throws IOException {
        if (!"POST".equalsIgnoreCase(ex.getRequestMethod())) {
            respond(ex, 405, "application/json", "{\"error\":\"POST required\"}".getBytes()); return;
        }
        String name = sanitizeEvidenceName(queryParam(ex, "name"));
        if (name.isEmpty()) {
            respond(ex, 400, "application/json", "{\"error\":\"name is required\"}".getBytes()); return;
        }
        String kind = queryParam(ex, "kind");
        if (kind == null || !(kind.equals("ppi") || kind.equals("expression") || kind.equals("other"))) kind = "other";
        String geneColParam = queryParam(ex, "gene_col");

        String text = new String(readAll(ex.getRequestBody()), "UTF-8");
        if (text.trim().isEmpty()) {
            respond(ex, 400, "application/json", "{\"error\":\"Empty upload\"}".getBytes()); return;
        }
        String[] lines = text.split("\r\n|\n|\r");
        int hIdx = 0;
        while (hIdx < lines.length && lines[hIdx].trim().isEmpty()) hIdx++;
        if (hIdx >= lines.length) {
            respond(ex, 400, "application/json", "{\"error\":\"No header row found\"}".getBytes()); return;
        }
        String headerLine = lines[hIdx];
        String delim = headerLine.contains("\t") ? "\t" : ",";
        String[] headers = splitDelimited(headerLine, delim);
        for (int i = 0; i < headers.length; i++) headers[i] = headers[i].trim();

        int geneColIdx = -1;
        if (geneColParam != null && !geneColParam.trim().isEmpty()) {
            for (int i = 0; i < headers.length; i++) {
                if (headers[i].equalsIgnoreCase(geneColParam.trim())) { geneColIdx = i; break; }
            }
            if (geneColIdx < 0) {
                respond(ex, 400, "application/json",
                    ("{\"error\":\"gene_col '" + escJ(geneColParam) + "' not found in header\"}").getBytes()); return;
            }
        } else {
            String[] candidates = {"gene", "gene_symbol", "symbol", "gene_name"};
            search:
            for (String cand : candidates) {
                for (int i = 0; i < headers.length; i++) {
                    if (headers[i].equalsIgnoreCase(cand)) { geneColIdx = i; break search; }
                }
            }
            if (geneColIdx < 0) {
                respond(ex, 400, "application/json",
                    "{\"error\":\"Could not auto-detect a gene-symbol column (expected one of gene/gene_symbol/symbol/gene_name); specify gene_col explicitly\"}".getBytes());
                return;
            }
        }
        String geneColName = headers[geneColIdx];

        List<String[]> rows = new ArrayList<>();
        for (int i = hIdx + 1; i < lines.length; i++) {
            if (lines[i].trim().isEmpty()) continue;
            rows.add(splitDelimited(lines[i], delim));
        }
        if (rows.isEmpty()) {
            respond(ex, 400, "application/json", "{\"error\":\"No data rows found\"}".getBytes()); return;
        }

        List<String> otherCols = new ArrayList<>();
        for (String h : headers) if (!h.equals(geneColName)) otherCols.add(h);
        Map<String, Boolean> numericFlag = new LinkedHashMap<>();
        for (String col : otherCols) numericFlag.put(col, true);
        for (String[] r : rows) {
            for (int i = 0; i < headers.length; i++) {
                String h = headers[i];
                if (h.equals(geneColName)) continue;
                String v = i < r.length ? r[i].trim() : "";
                if (v.isEmpty()) continue;
                if (Boolean.TRUE.equals(numericFlag.get(h))) {
                    try { Double.parseDouble(v); } catch (NumberFormatException e) { numericFlag.put(h, false); }
                }
            }
        }

        File evDir = new File(projectDir, "evidence");
        evDir.mkdirs();
        File tsvFile = new File(evDir, name + ".tsv");
        File metaFile = new File(evDir, name + ".meta.json");
        try (PrintWriter pw = new PrintWriter(new OutputStreamWriter(new FileOutputStream(tsvFile), "UTF-8"))) {
            pw.println(String.join("\t", headers));
            for (String[] r : rows) {
                StringBuilder line = new StringBuilder();
                for (int i = 0; i < headers.length; i++) {
                    if (i > 0) line.append('\t');
                    String v = i < r.length ? r[i].trim() : "";
                    line.append(v.replace("\t", " "));
                }
                pw.println(line);
            }
        }

        StringBuilder meta = new StringBuilder();
        meta.append("{");
        meta.append("\"name\":\"").append(escJ(name)).append("\",");
        meta.append("\"kind\":\"").append(escJ(kind)).append("\",");
        meta.append("\"gene_col\":\"").append(escJ(geneColName)).append("\",");
        meta.append("\"row_count\":").append(rows.size()).append(",");
        meta.append("\"columns\":[");
        for (int i = 0; i < otherCols.size(); i++) {
            if (i > 0) meta.append(",");
            meta.append("{\"name\":\"").append(escJ(otherCols.get(i))).append("\",\"type\":\"")
                .append(Boolean.TRUE.equals(numericFlag.get(otherCols.get(i))) ? "numeric" : "categorical")
                .append("\"}");
        }
        meta.append("]}");
        try (PrintWriter pw = new PrintWriter(new OutputStreamWriter(new FileOutputStream(metaFile), "UTF-8"))) {
            pw.print(meta);
        }

        respond(ex, 200, "application/json", meta.toString().getBytes("UTF-8"));
        System.out.printf("[Server] Uploaded evidence '%s' (%s) to project '%s': %d rows%n",
            name, kind, projectId, rows.size());
    }

    // GET /api/project/{id}/evidence/{name} — full parsed content as JSON: gene -> {column: value}
    private void projectEvidenceGet(HttpExchange ex, String projectDir, String rawName) throws IOException {
        String name = sanitizeEvidenceName(rawName);
        File evDir = new File(projectDir, "evidence");
        File tsvFile = new File(evDir, name + ".tsv");
        File metaFile = new File(evDir, name + ".meta.json");
        if (!tsvFile.exists() || !metaFile.exists()) {
            respond(ex, 404, "application/json", "{\"error\":\"Evidence table not found\"}".getBytes()); return;
        }
        String metaJson = new String(Files.readAllBytes(metaFile.toPath()), "UTF-8");
        String geneCol = extractStr(metaJson, "gene_col");

        StringBuilder j = new StringBuilder("{");
        boolean first = true;
        try (BufferedReader br = new BufferedReader(new InputStreamReader(new FileInputStream(tsvFile), "UTF-8"))) {
            String header = br.readLine();
            if (header != null) {
                String[] headers = header.split("\t", -1);
                int geneIdx = -1;
                for (int i = 0; i < headers.length; i++) if (headers[i].equals(geneCol)) { geneIdx = i; break; }
                String line;
                while (geneIdx >= 0 && (line = br.readLine()) != null) {
                    if (line.trim().isEmpty()) continue;
                    String[] fields = line.split("\t", -1);
                    if (geneIdx >= fields.length) continue;
                    String gene = fields[geneIdx].trim();
                    if (gene.isEmpty()) continue;
                    if (!first) j.append(",");
                    first = false;
                    j.append("\"").append(escJ(gene)).append("\":{");
                    boolean f2 = true;
                    for (int i = 0; i < headers.length; i++) {
                        if (i == geneIdx) continue;
                        if (!f2) j.append(",");
                        f2 = false;
                        String v = i < fields.length ? fields[i].trim() : "";
                        j.append("\"").append(escJ(headers[i])).append("\":\"").append(escJ(v)).append("\"");
                    }
                    j.append("}");
                }
            }
        }
        j.append("}");
        respond(ex, 200, "application/json", j.toString().getBytes("UTF-8"));
    }

    // DELETE /api/project/{id}/evidence/{name}
    private void projectEvidenceDelete(HttpExchange ex, String projectDir, String rawName) throws IOException {
        String name = sanitizeEvidenceName(rawName);
        File evDir = new File(projectDir, "evidence");
        File tsvFile = new File(evDir, name + ".tsv");
        File metaFile = new File(evDir, name + ".meta.json");
        boolean existed = tsvFile.exists() || metaFile.exists();
        tsvFile.delete();
        metaFile.delete();
        if (!existed) {
            respond(ex, 404, "application/json", "{\"error\":\"Evidence table not found\"}".getBytes()); return;
        }
        respond(ex, 200, "application/json", "{\"ok\":true}".getBytes());
    }

    // GET /api/project/{id}/enrichment?evidence=<name> — locus-based enrichment test (Fisher's
    // exact for categorical columns, Mann-Whitney U for numeric columns) of this project's
    // identified-loci genes against a background gene set, run live (not cached/precomputed).
    private void projectEnrichment(HttpExchange ex, String projectId, String projectDir) throws IOException {
        String evidenceName = queryParam(ex, "evidence");
        if (evidenceName == null || evidenceName.trim().isEmpty()) {
            respond(ex, 400, "application/json", "{\"error\":\"evidence query param is required\"}".getBytes()); return;
        }
        String name = sanitizeEvidenceName(evidenceName);
        File evDir = new File(projectDir, "evidence");
        File tsvFile = new File(evDir, name + ".tsv");
        File metaFile = new File(evDir, name + ".meta.json");
        if (!tsvFile.exists() || !metaFile.exists()) {
            respond(ex, 404, "application/json", "{\"error\":\"Evidence table not found\"}".getBytes()); return;
        }

        ProjectState ps = ensureProjectState(projectId, projectDir);
        if (ps == null || ps.gff == null || ps.outputs == null) {
            respond(ex, 400, "application/json", "{\"error\":\"Project not available\"}".getBytes()); return;
        }

        try {
            String metaJson = new String(Files.readAllBytes(metaFile.toPath()), "UTF-8");
            String geneCol = extractStr(metaJson, "gene_col");
            Map<String, String> columnTypes = new LinkedHashMap<>();
            for (String c : extractObjectArray(metaJson, "columns")) {
                String cname = extractStr(c, "name");
                String ctype = extractStr(c, "type");
                if (cname != null) columnTypes.put(cname, ctype);
            }

            Map<String, Map<String, String>> evidence = new LinkedHashMap<>();
            try (BufferedReader br = new BufferedReader(new InputStreamReader(new FileInputStream(tsvFile), "UTF-8"))) {
                String header = br.readLine();
                if (header != null) {
                    String[] headers = header.split("\t", -1);
                    int geneIdx = -1;
                    for (int i = 0; i < headers.length; i++) if (headers[i].equals(geneCol)) { geneIdx = i; break; }
                    String line;
                    while (geneIdx >= 0 && (line = br.readLine()) != null) {
                        if (line.trim().isEmpty()) continue;
                        String[] fields = line.split("\t", -1);
                        if (geneIdx >= fields.length) continue;
                        String gene = fields[geneIdx].trim();
                        if (gene.isEmpty()) continue;
                        Map<String, String> row = new LinkedHashMap<>();
                        for (int i = 0; i < headers.length; i++) {
                            if (i == geneIdx) continue;
                            row.put(headers[i], i < fields.length ? fields[i].trim() : "");
                        }
                        evidence.put(gene, row);
                    }
                }
            }

            // In-loci genes: genes overlapping this project's own identified loci (already-parsed
            // annotation, reused as-is — see LocusOutput.genes / ensureProjectState).
            Set<String> inLociGenes = new LinkedHashSet<>();
            Set<String> chrs = new LinkedHashSet<>();
            for (LocusOutput lo : ps.outputs) {
                chrs.add(lo.chr);
                for (Gene g : lo.genes) if (g.geneName != null) inLociGenes.add(g.geneName);
            }

            // Background: every gene in the parsed GFF3 on the same chromosome(s) as those loci.
            Set<String> backgroundGenes = new LinkedHashSet<>();
            for (String chr : chrs) {
                for (Gene g : ps.gff.overlapping(chr, 0, Long.MAX_VALUE)) {
                    if (g.geneName != null) backgroundGenes.add(g.geneName);
                }
            }

            EnrichmentAnalyzer.Result result =
                EnrichmentAnalyzer.run(inLociGenes, backgroundGenes, evidence, columnTypes);
            respond(ex, 200, "application/json", result.toJson().getBytes("UTF-8"));
        } catch (Exception e) {
            respond(ex, 500, "application/json", ("{\"error\":\"" + escJ(e.getMessage()) + "\"}").getBytes());
        }
    }

    private static String queryParam(HttpExchange ex, String key) throws IOException {
        String query = ex.getRequestURI().getQuery();
        if (query == null) return null;
        for (String param : query.split("&")) {
            String[] kv = param.split("=", 2);
            if (kv.length == 2 && kv[0].equals(key)) return URLDecoder.decode(kv[1], "UTF-8");
        }
        return null;
    }

    // ══════════════════════════════════════════════════════════════════════
    //  SUSIEX (cross-ancestry, cross-project fine-mapping)
    // ══════════════════════════════════════════════════════════════════════

    private static class SusiexJob {
        String status = "running"; // running | done | error
        String error;
        SusiexAdapter.RunResult result;
        String runDir;
    }
    private final Map<String, SusiexJob> susiexJobs = new ConcurrentHashMap<>();

    // POST /api/susiex-run — {"members":[{"project_id":"x","locus_index":1}, ...], "params":{...}}
    private void susiexRun(HttpExchange ex) throws IOException {
        cors(ex); if (preflight(ex)) return;
        if (!"POST".equalsIgnoreCase(ex.getRequestMethod())) {
            respond(ex, 405, "application/json", "{\"error\":\"POST required\"}".getBytes()); return;
        }
        String body = new String(readAll(ex.getRequestBody()), "UTF-8");
        List<String> memberObjs = extractObjectArray(body, "members");
        if (memberObjs.size() < 2) {
            respond(ex, 400, "application/json",
                "{\"error\":\"At least two project/locus members are required\"}".getBytes()); return;
        }

        String paramsObj = extractRawObject(body, "params");
        double pvalThresh = parseDoubleOr(extractStr(paramsObj, "pval_thresh"), 1e-5);
        double mafThresh  = parseDoubleOr(extractStr(paramsObj, "maf_thresh"), 0.005);
        int maxCausal      = (int) parseDoubleOr(extractStr(paramsObj, "max_causal"), 5);
        String susiexPath = extractStr(paramsObj, "susiex_path");
        String plinkPath  = extractStr(paramsObj, "plink_path");
        if (susiexPath == null || susiexPath.trim().isEmpty()) {
            respond(ex, 400, "application/json",
                "{\"error\":\"susiex_path (path to SuSiEx.py) is required\"}".getBytes()); return;
        }

        // Resolve each member: project + locus + its base-pipeline artifacts
        List<SusiexAdapter.Member> members = new ArrayList<>();
        List<Config> memberConfigs = new ArrayList<>();
        List<Locus> memberLoci = new ArrayList<>();
        String sharedChr = null;
        long regionStart = Long.MAX_VALUE, regionEnd = Long.MIN_VALUE;
        for (String mo : memberObjs) {
            String pid = extractStr(mo, "project_id");
            String idxStr = extractStr(mo, "locus_index");
            if (pid == null || idxStr == null) {
                respond(ex, 400, "application/json",
                    "{\"error\":\"Each member needs project_id and locus_index\"}".getBytes()); return;
            }
            File dir = new File("projects/" + pid);
            ProjectState ps = ensureProjectState(pid, dir.getAbsolutePath());
            if (ps == null || ps.config == null || ps.loci == null) {
                respond(ex, 400, "application/json",
                    ("{\"error\":\"Project not available: " + escJ(pid) + "\"}").getBytes()); return;
            }
            int locusIndex;
            try { locusIndex = Integer.parseInt(idxStr); } catch (NumberFormatException e) {
                respond(ex, 400, "application/json", "{\"error\":\"Invalid locus_index\"}".getBytes()); return;
            }
            Locus locus = null;
            for (Locus l : ps.loci) if (l.index == locusIndex) { locus = l; break; }
            if (locus == null) {
                respond(ex, 400, "application/json",
                    ("{\"error\":\"Locus " + locusIndex + " not found in project " + escJ(pid) + "\"}").getBytes()); return;
            }
            if (ps.config.sampleN <= 0) {
                respond(ex, 400, "application/json",
                    ("{\"error\":\"Project '" + escJ(pid) + "' has no sample size (N) configured\"}").getBytes()); return;
            }
            memberConfigs.add(ps.config);
            memberLoci.add(locus);
            if (sharedChr == null) sharedChr = locus.chr;
            regionStart = Math.min(regionStart, locus.start);
            regionEnd = Math.max(regionEnd, locus.end);

            SusiexAdapter.Member m = new SusiexAdapter.Member();
            ProjectMetadata pm = ProjectMetadata.load(dir.getAbsolutePath());
            m.label = (pm != null && !pm.name.isEmpty()) ? pm.name : pid;
            m.sampleN = ps.config.sampleN;
            File analysisRoot = BaseStepPipeline.analysisDir(ps.config, locus);
            m.harmonizedDir = new File(analysisRoot, "harmonized");
            m.matchedDir = new File(analysisRoot, "matched");
            members.add(m);
        }

        String jobId = UUID.randomUUID().toString();
        SusiexJob job = new SusiexJob();
        job.runDir = "output/susiex/" + jobId;
        susiexJobs.put(jobId, job);

        final String finalChr = sharedChr;
        final long finalStart = regionStart, finalEnd = regionEnd;
        final double finalPval = pvalThresh, finalMaf = mafThresh;
        final int finalMaxCausal = maxCausal;
        final String finalSusiexPath = susiexPath.trim();
        final String finalPlinkPath = (plinkPath != null && !plinkPath.trim().isEmpty()) ? plinkPath.trim() : "plink";

        new Thread(() -> {
            try {
                // Ensure each member's base artifacts (harmonized GWAS + matched ref panel) exist.
                for (int i = 0; i < memberConfigs.size(); i++) {
                    BaseStepPipeline.runAll(memberConfigs.get(i), memberLoci.get(i), LdMatrixComputer.DEFAULT_LD_WINDOW);
                }
                File runDir = new File(job.runDir);
                SusiexAdapter.RunResult r = SusiexAdapter.prepareAndRun(members, runDir,
                    finalChr, finalStart, finalEnd, finalSusiexPath, finalPlinkPath,
                    finalMaxCausal, finalPval, finalMaf);
                job.result = r;
                job.status = r.ok ? "done" : "error";
                job.error = r.error;
            } catch (Exception e) {
                job.status = "error";
                job.error = e.getMessage();
                System.err.printf("[SuSiEx] Job '%s' failed: %s%n", jobId, e.getMessage());
                e.printStackTrace(System.err);
            }
        }, "susiex-" + jobId).start();

        respond(ex, 202, "application/json", ("{\"job_id\":\"" + jobId + "\",\"status\":\"started\"}").getBytes());
    }

    // GET /api/susiex-progress?job=<jobId>
    private void susiexProgress(HttpExchange ex) throws IOException {
        cors(ex); if (preflight(ex)) return;
        String jobId = queryParam(ex, "job");
        SusiexJob job = jobId != null ? susiexJobs.get(jobId) : null;
        if (job == null) {
            respond(ex, 404, "application/json", "{\"error\":\"Job not found\"}".getBytes()); return;
        }
        StringBuilder json = new StringBuilder("{");
        json.append("\"status\":\"").append(escJ(job.status)).append('"');
        if (job.error != null) json.append(",\"error\":\"").append(escJ(job.error)).append('"');
        json.append('}');
        respond(ex, 200, "application/json", json.toString().getBytes("UTF-8"));
    }

    // GET /api/susiex-result?job=<jobId>
    private void susiexResult(HttpExchange ex) throws IOException {
        cors(ex); if (preflight(ex)) return;
        String jobId = queryParam(ex, "job");
        SusiexJob job = jobId != null ? susiexJobs.get(jobId) : null;
        if (job == null || job.result == null) {
            respond(ex, 404, "application/json", "{\"error\":\"Job not found or not complete\"}".getBytes()); return;
        }
        SusiexAdapter.RunResult r = job.result;
        StringBuilder json = new StringBuilder("{");
        json.append("\"ok\":").append(r.ok);
        if (r.error != null) json.append(",\"error\":\"").append(escJ(r.error)).append('"');
        json.append(",\"n_credible_sets\":").append(r.nCredibleSets);
        json.append(",\"rows\":[");
        for (int i = 0; i < r.rows.size(); i++) {
            if (i > 0) json.append(',');
            SusiexAdapter.ResultRow row = r.rows.get(i);
            json.append("{\"snp_id\":\"").append(escJ(row.snpId)).append('"');
            json.append(",\"chr\":\"").append(escJ(row.chr != null ? row.chr : "")).append('"');
            json.append(",\"pos\":").append(row.pos);
            json.append(",\"pip\":").append(Double.isNaN(row.pip) ? "null" : row.pip);
            json.append(",\"cs_id\":\"").append(escJ(row.csId)).append('"');
            json.append('}');
        }
        json.append("]}");
        respond(ex, 200, "application/json", json.toString().getBytes("UTF-8"));
    }

    private static double parseDoubleOr(String s, double def) {
        if (s == null || s.trim().isEmpty()) return def;
        try { return Double.parseDouble(s.trim()); } catch (NumberFormatException e) { return def; }
    }

    /** Returns the raw {...} substring for a top-level object-valued key, or "{}" if absent. */
    private static String extractRawObject(String json, String key) {
        String marker = "\"" + key + "\":";
        int i = json.indexOf(marker);
        if (i < 0) return "{}";
        i += marker.length();
        while (i < json.length() && json.charAt(i) != '{') i++;
        if (i >= json.length()) return "{}";
        int depth = 0, start = i;
        for (int j = i; j < json.length(); j++) {
            char c = json.charAt(j);
            if (c == '{') depth++;
            else if (c == '}') { depth--; if (depth == 0) return json.substring(start, j + 1); }
        }
        return "{}";
    }

    /** Splits a top-level array-of-objects value into raw {...} substrings. */
    private static List<String> extractObjectArray(String json, String key) {
        List<String> result = new ArrayList<>();
        String marker = "\"" + key + "\":";
        int i = json.indexOf(marker);
        if (i < 0) return result;
        i += marker.length();
        while (i < json.length() && json.charAt(i) != '[') i++;
        if (i >= json.length()) return result;
        int arrEnd = i, depth = 0;
        for (int j = i; j < json.length(); j++) {
            char c = json.charAt(j);
            if (c == '[') depth++;
            else if (c == ']') { depth--; if (depth == 0) { arrEnd = j; break; } }
        }
        String inner = json.substring(i + 1, arrEnd);
        int objDepth = 0, objStart = -1;
        for (int j = 0; j < inner.length(); j++) {
            char c = inner.charAt(j);
            if (c == '{') { if (objDepth == 0) objStart = j; objDepth++; }
            else if (c == '}') { objDepth--; if (objDepth == 0 && objStart >= 0) result.add(inner.substring(objStart, j + 1)); }
        }
        return result;
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
        StringBuilder b = new StringBuilder(s.length() + 16);
        for (int i = 0; i < s.length(); i++) {
            char c = s.charAt(i);
            switch (c) {
                case '\\': b.append("\\\\"); break;
                case '"':  b.append("\\\""); break;
                case '\n': b.append("\\n"); break;
                case '\r': b.append("\\r"); break;
                case '\t': b.append("\\t"); break;
                default:
                    // Other control characters are not allowed raw inside a JSON string
                    if (c < 0x20) b.append(String.format("\\u%04x", (int) c)); else b.append(c);
            }
        }
        return b.toString();
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

    /** Like extractStr, but only searches from a given offset (e.g. inside a nested object). */
    private static String extractNestedStr(String json, int fromIndex, String key) {
        String marker = "\"" + key + "\":";
        int i = json.indexOf(marker, fromIndex);
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

    private static List<String> extractStringArray(String json, String key) {
        List<String> result = new ArrayList<>();
        String marker = "\"" + key + "\":";
        int i = json.indexOf(marker);
        if (i < 0) return result;
        i += marker.length();
        while (i < json.length() && json.charAt(i) != '[') i++;
        if (i >= json.length()) return result;
        int end = json.indexOf(']', i);
        if (end < 0) return result;
        String arrContent = json.substring(i + 1, end);
        for (String part : arrContent.split(",")) {
            part = part.trim();
            if (part.startsWith("\"") && part.endsWith("\""))
                result.add(part.substring(1, part.length() - 1));
        }
        return result;
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
