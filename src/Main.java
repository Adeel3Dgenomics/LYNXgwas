import java.io.*;
import java.nio.file.*;
import java.util.*;

public class Main {

    static class PipelineResult {
        Config             config;
        GffParser          gff;
        List<Locus>        loci;
        List<LocusOutput>  outputs;
    }

    public static void main(String[] args) throws Exception {
        System.out.println("=== LYNXgwas – Locus Analysis and Genomic Explorer ===");
        System.out.printf("Pipeline version: %s%n%n", Config.PIPELINE_VERSION);

        // Parse orchestration flags
        String forceProjectId = null;
        boolean forceAll = false;
        for (String arg : args) {
            if (arg.startsWith("--project=")) forceProjectId = arg.substring("--project=".length());
            else if (arg.equals("--all"))     forceAll = true;
        }

        // Migrate legacy single-project layout if needed (Step 11)
        migrateLegacyIfNeeded();

        // Discover projects
        List<File> projectDirs = discoverProjects(new File("projects"));

        if (!projectDirs.isEmpty()) {
            runMultiProject(projectDirs, forceProjectId, forceAll);
        } else if (new File("config.properties").exists()) {
            System.out.println("[INFO] No projects found after migration check. Running single-project mode.");
            System.out.println("[TIP]  Create projects/{id}/config.properties for multi-project mode.\n");
            runLegacySingleProject(args);
        } else {
            System.err.println("[ERROR] No projects found in projects/ and no config.properties at root.");
            System.err.println("Create projects/{id}/config.properties or place config.properties in the root.");
            System.exit(1);
        }
    }

    // ── Multi-project orchestration ──────────────────────────────────────

    private static void runMultiProject(List<File> projectDirs, String forceProjectId, boolean forceAll) {
        // Filter to specific project if --project= given
        if (forceProjectId != null) {
            String id = forceProjectId;
            projectDirs.removeIf(d -> !d.getName().equals(id));
            if (projectDirs.isEmpty()) {
                System.err.printf("[ERROR] Project '%s' not found in projects/%n", id);
                System.exit(1);
            }
        }

        // Determine staleness for each project
        List<String> staleIds   = new ArrayList<>();
        List<String> currentIds = new ArrayList<>();
        Map<String, Config> configs = new LinkedHashMap<>();
        Map<String, ProjectMetadata.StaleReason> reasons = new LinkedHashMap<>();

        for (File dir : projectDirs) {
            String id = dir.getName();
            try {
                Config config = Config.loadFromProject(dir.getAbsolutePath());
                configs.put(id, config);

                if (forceAll || forceProjectId != null) {
                    staleIds.add(id);
                    reasons.put(id, ProjectMetadata.StaleReason.NO_METADATA);
                } else {
                    ProjectMetadata.StaleReason reason =
                        ProjectMetadata.checkStaleness(dir.getAbsolutePath(), config);
                    reasons.put(id, reason);
                    if (reason == ProjectMetadata.StaleReason.NOT_STALE) currentIds.add(id);
                    else staleIds.add(id);
                }
            } catch (Exception e) {
                System.err.printf("[ERROR] Could not load config for project '%s': %s%n", id, e.getMessage());
                reasons.put(id, ProjectMetadata.StaleReason.NO_METADATA);
            }
        }

        // Pre-processing summary
        System.out.printf("Discovered %d project(s): %d stale, %d up to date%n",
            configs.size(), staleIds.size(), currentIds.size());
        for (String id : staleIds)
            System.out.printf("  %-20s  → needs reprocessing (%s)%n", id, reasonLabel(reasons.get(id)));
        for (String id : currentIds)
            System.out.printf("  %-20s  → up to date%n", id);
        System.out.println();

        // Start server (serves home page + API regardless of whether any projects are stale)
        ProgressTracker progress = LocalServer.progress;
        LocalServer server = null;
        try {
            server = new LocalServer("output");
            server.start();
            System.out.printf("Server at http://localhost:%d/%n%n", LocalServer.PORT);
        } catch (Exception e) {
            System.err.println("[WARN] Could not start server: " + e.getMessage());
        }

        // Only auto-process if explicitly requested via --project= or --all
        if (forceProjectId == null && !forceAll) {
            if (!staleIds.isEmpty()) {
                System.out.printf("%d project(s) need reprocessing. Use 'Reprocess' on the home page or run with --all.%n", staleIds.size());
            } else {
                System.out.println("All projects are up to date.");
            }
            if (server != null) {
                System.out.printf("Home page: http://localhost:%d/%nPress Ctrl+C to stop.%n", LocalServer.PORT);
            }
            return;
        }

        if (staleIds.isEmpty()) {
            System.out.println("All projects are up to date. Nothing to process.");
            if (server != null) {
                System.out.printf("Home page: http://localhost:%d/%nPress Ctrl+C to stop.%n", LocalServer.PORT);
            }
            return;
        }

        // Process stale projects sequentially (only when explicitly requested)
        int processed = 0, succeeded = 0, failed = 0;
        List<String[]> summaryRows = new ArrayList<>();
        PipelineResult lastResult = null;

        for (String id : staleIds) {
            Config config = configs.get(id);
            if (config == null) {
                summaryRows.add(new String[]{id, "SKIP", "config load failed"});
                failed++;
                continue;
            }

            String projectDir = new File("projects", id).getAbsolutePath();

            // Skip projects that have no loci file — they need "Get Loci" first
            if (config.lociFile.isEmpty() || !new File(config.lociFile).exists()) {
                summaryRows.add(new String[]{id, "SKIP", "no loci file — use Get Loci"});
                System.out.printf("[SKIP] Project '%s': no loci file. Use 'Get Loci' on the home page.%n", id);
                continue;
            }

            System.out.printf("────────────────────────────────────────────────────%n");
            System.out.printf("Processing project: %s%n", id);
            System.out.printf("────────────────────────────────────────────────────%n");
            progress.update("Project: " + id, 0, 1);
            progress.done = false;

            long t0 = System.currentTimeMillis();
            try {
                PipelineResult result = runPipeline(config, projectDir, progress);
                long elapsed = System.currentTimeMillis() - t0;

                // Write project.json with fresh metadata
                writeProjectMetadata(id, projectDir, config, result);

                summaryRows.add(new String[]{id, "OK",
                    String.format("%.1fs, %d loci", elapsed / 1000.0, result.outputs.size())});
                succeeded++;
                lastResult = result;

                printLociTable(result.outputs);

            } catch (Exception e) {
                long elapsed = System.currentTimeMillis() - t0;
                System.err.printf("[ERROR] Project '%s' failed after %.1fs: %s%n",
                    id, elapsed / 1000.0, e.getMessage());
                e.printStackTrace(System.err);
                summaryRows.add(new String[]{id, "FAIL", e.getMessage()});
                failed++;
                // Do NOT write project.json for failed runs — stale state is preserved
            }
            processed++;
        }

        progress.done = true;

        // Final summary
        System.out.printf("%n════════════════════════════════════════════════════%n");
        System.out.printf("Processing complete: %d processed, %d succeeded, %d failed%n",
            processed, succeeded, failed);
        System.out.printf("════════════════════════════════════════════════════%n");
        System.out.printf("%-20s %-6s %s%n", "Project", "Status", "Details");
        System.out.println("-".repeat(60));
        for (String[] row : summaryRows)
            System.out.printf("%-20s %-6s %s%n", row[0], row[1], row[2]);

        // Give the server the last successfully processed project's state
        if (server != null && lastResult != null) {
            server.setPipelineState(lastResult.config, lastResult.gff,
                lastResult.loci, lastResult.outputs);
            System.out.printf("%nViewer server running. Press Ctrl+C to stop.%n");
        } else if (server != null) {
            System.out.printf("%nServer running (no projects processed). Press Ctrl+C to stop.%n");
        }
    }

    // ── Pipeline execution (one project) ─────────────────────────────────

    static PipelineResult runPipeline(Config config, String projectDir,
                                      ProgressTracker progress) throws Exception {
        System.out.printf("GWAS file  : %s%n", config.gwasFile);
        System.out.printf("Loci file  : %s%n", config.lociFile);
        System.out.printf("GFF3 file  : %s%n", config.gff3File);
        System.out.printf("LD enabled : %s%n", config.ldEnabled);
        System.out.printf("Ref panel  : %s (%s)%n",
            config.refPanelPath.isEmpty() ? "(none)" : config.refPanelPath,
            config.refPanelPopulation);
        System.out.printf("Output     : %s%n%n", config.outputDir);

        // Create output directories
        new File(config.outputDir + "/data").mkdirs();
        new File(config.outputDir + "/plots").mkdirs();
        new File(config.outputDir + "/tmp").mkdirs();

        if (config.ldEnabled) {
            new File(config.plinkSubsetsDir()).mkdirs();
            new File(config.ldResultsDir()).mkdirs();
        }

        // Set up phased progress: parsing, [plink subsets], [LD], export
        int numPhases = config.ldEnabled ? 4 : 2;
        progress.setPhases(numPhases);

        // ── Phase 1: Parse inputs ────────────────────────────────────────
        progress.nextPhase("Parsing loci & GFF3", 4);
        List<Locus> loci = LociParser.parse(config);
        progress.advance(1);

        System.out.println("[Phase 1] Parsing GFF3 annotation...");
        GffParser gff = GffParser.parse(config);
        progress.advance(2);

        System.out.println("[Phase 1] Streaming GWAS summary statistics...");
        GwasParser.parse(config, loci);
        progress.advance(3);

        Map<Integer, Snp> topSnps = new LinkedHashMap<>();
        for (Locus locus : loci) {
            Snp top = locus.snps.stream()
                .min(Comparator.comparingDouble(s -> s.pvalue))
                .orElse(null);
            topSnps.put(locus.index, top);
        }

        System.out.println("[Phase 1] Building genome-wide skyline...");
        GenomeSkyline.generate(config);
        progress.advance(4);

        if (config.colRsid.isEmpty() && config.topSnpFile.isEmpty())
            new SnpAnnotator().annotateTopSnps(topSnps, config);

        // ── Phase 2: PLINK subset extraction ─────────────────────────────
        Map<Integer, PlinkSubsetter.SubsetResult> subsets = new LinkedHashMap<>();
        String plinkBin = null;
        if (config.ldEnabled) {
            plinkBin = PlinkSubsetter.findPlink(config);
            if (plinkBin == null) {
                System.err.println("[WARN] PLINK executable not found — disabling LD");
                config.ldEnabled = false;
            }
        }
        if (config.ldEnabled) {
            progress.nextPhase("Extracting PLINK subsets", loci.size());
            System.out.printf("[Phase 2] Extracting PLINK subsets for %d loci...%n", loci.size());
            subsets = PlinkSubsetter.subsetAll(loci, topSnps, plinkBin, config, progress);
        }

        // ── Phase 3: LD computation ──────────────────────────────────────
        Map<Integer, LdCalculator.LdResult> ldResults = new LinkedHashMap<>();
        if (config.ldEnabled) {
            progress.nextPhase("Computing LD", loci.size());
            System.out.printf("[Phase 3] Computing LD (%d parallel jobs)...%n",
                config.ldParallelJobs);
            ldResults = LdCalculator.computeAll(
                loci, topSnps, subsets, plinkBin, config, progress);
        }

        // ── Phase 4: Build output + export ───────────────────────────────
        progress.nextPhase("Exporting JSON", loci.size());
        System.out.println("[Phase 4] Building locus output...");
        List<LocusOutput> outputs = new ArrayList<>();

        for (int i = 0; i < loci.size(); i++) {
            Locus locus = loci.get(i);
            LocusOutput lo = new LocusOutput();
            lo.id          = locus.id;
            lo.locusIndex  = locus.index;
            lo.locusName   = "Locus " + locus.index;
            lo.chr         = locus.chr;
            lo.start       = locus.start;
            lo.end         = locus.end;
            lo.paddedStart = locus.paddedStart;
            lo.paddedEnd   = locus.paddedEnd;
            lo.refPanel    = config.ldEnabled ? config.refPanelPopulation : "";
            lo.topSnp      = topSnps.get(locus.index);

            lo.genes        = gff.overlapping(locus.chr, locus.paddedStart, locus.paddedEnd);
            lo.nearestGenes = gff.nearestGeneNames(locus.mid(), new ArrayList<>(lo.genes));

            lo.gwasSnps = new ArrayList<>(locus.snps);
            lo.gwasSnps.sort(Comparator.comparingLong(s -> s.pos));

            LdCalculator.LdResult ld = ldResults.get(locus.index);
            if (ld != null && !ld.ldFailed) {
                for (Snp snp : lo.gwasSnps) {
                    Double r2 = ld.r2ByPos.get(snp.chr + ":" + snp.pos);
                    if (r2 != null) snp.r2 = r2;
                }
                lo.ldTriangle = ld.triangle;
            } else if (lo.topSnp != null && !config.ldEnabled) {
                for (Snp snp : lo.gwasSnps)
                    if (snp.id.equals(lo.topSnp.id)) { snp.r2 = 1.0; break; }
            }

            lo.locusContext = new LocusOutput.LocusContext();
            if (i > 0) {
                Locus prev = loci.get(i - 1);
                long dist = locus.chr.equals(prev.chr)
                    ? Math.abs(locus.start - prev.end) : Long.MAX_VALUE;
                lo.locusContext.prevLocus = new LocusOutput.LocusRef(
                    prev.index, prev.chr, prev.start, prev.end, prev.mid(), dist);
            }
            if (i < loci.size() - 1) {
                Locus next = loci.get(i + 1);
                long dist = locus.chr.equals(next.chr)
                    ? Math.abs(next.start - locus.end) : Long.MAX_VALUE;
                lo.locusContext.nextLocus = new LocusOutput.LocusRef(
                    next.index, next.chr, next.start, next.end, next.mid(), dist);
            }

            outputs.add(lo);
            progress.advance(i + 1);
        }

        JsonExporter.export(loci, outputs, config);

        PipelineResult result = new PipelineResult();
        result.config  = config;
        result.gff     = gff;
        result.loci    = loci;
        result.outputs = outputs;
        return result;
    }

    // ── Project metadata write ───────────────────────────────────────────

    private static void writeProjectMetadata(String id, String projectDir,
                                             Config config, PipelineResult result) {
        try {
            // Preserve user-set fields from previous project.json
            ProjectMetadata existing = ProjectMetadata.load(projectDir);

            ProjectMetadata pm = new ProjectMetadata();
            pm.id          = id;
            pm.name        = (existing != null && !existing.name.isEmpty()) ? existing.name : id;
            pm.description = (existing != null) ? existing.description : "";
            pm.lociCount   = result.outputs.size();

            // Count unique SNPs across all loci
            Set<String> uniqueSnps = new HashSet<>();
            for (LocusOutput o : result.outputs)
                for (Snp s : o.gwasSnps)
                    uniqueSnps.add(s.chr + ":" + s.pos);
            pm.totalSnps = uniqueSnps.size();

            String annotPath = new File(projectDir, "annotations.yaml").getAbsolutePath();
            pm.countAnnotationSources(annotPath);
            pm.coreInputFingerprint  = ProjectMetadata.computeCoreInputFingerprint(config, projectDir);
            pm.annotationFingerprint = ProjectMetadata.computeAnnotationFingerprint(annotPath);

            // Preserve rsID recovery fields from previous project.json
            if (existing != null) {
                pm.rsidColumnPresent   = existing.rsidColumnPresent;
                pm.rsidRecoveryStatus  = existing.rsidRecoveryStatus;
                pm.rsidRecoveryRate    = existing.rsidRecoveryRate;
                pm.rsidRecoveryDate    = existing.rsidRecoveryDate;
                pm.selectedSnpDatabase = existing.selectedSnpDatabase;
            }

            // Auto-detect rsID presence from the config's rsid column mapping
            if (!config.colRsid.isEmpty()) {
                pm.rsidColumnPresent = true;
                if ("not_started".equals(pm.rsidRecoveryStatus)) {
                    pm.rsidRecoveryStatus = "completed";
                }
            }

            pm.save(projectDir);
            System.out.printf("[Metadata] Wrote project.json for '%s': %d loci, %d unique SNPs, rsid=%s%n",
                id, pm.lociCount, pm.totalSnps, pm.rsidColumnPresent ? "present" : "missing");
        } catch (Exception e) {
            System.err.printf("[WARN] Could not write project.json for '%s': %s%n", id, e.getMessage());
        }
    }

    // ── Project discovery ────────────────────────────────────────────────

    static List<File> discoverProjects(File projectsRoot) {
        List<File> projects = new ArrayList<>();
        if (!projectsRoot.isDirectory()) return projects;

        File[] children = projectsRoot.listFiles();
        if (children == null) return projects;

        for (File child : children) {
            if (!child.isDirectory()) continue;
            if (new File(child, "config.properties").exists())
                projects.add(child);
        }

        projects.sort(Comparator.comparing(File::getName));
        return projects;
    }

    // ── Legacy migration (Step 11) ─────────────────────────────────────

    private static void migrateLegacyIfNeeded() {
        File projectsDir = new File("projects");
        // Only migrate if projects/ has no project subdirectories yet
        // (the template file may already exist from Step 7)
        if (!discoverProjects(projectsDir).isEmpty()) return;

        File legacyManifest = new File("output/data/manifest.json");
        File legacyConfig   = new File("config.properties");
        if (!legacyManifest.exists() || !legacyConfig.exists()) return;

        System.out.println("═══════════════════════════════════════════════════");
        System.out.println("  Migrating legacy single-project layout...");
        System.out.println("═══════════════════════════════════════════════════");

        File defaultDir = new File("projects/default");
        try {
            // 1. Create projects/default/
            defaultDir.mkdirs();

            // 2. Copy config.properties
            Files.copy(legacyConfig.toPath(),
                new File(defaultDir, "config.properties").toPath(),
                StandardCopyOption.REPLACE_EXISTING);
            System.out.println("  Copied config.properties → projects/default/config.properties");

            // 3. Copy annotations.yaml if it exists
            File legacyAnnot = new File("output/annotations.yaml");
            if (!legacyAnnot.exists()) legacyAnnot = new File("annotations.yaml");
            if (legacyAnnot.exists()) {
                Files.copy(legacyAnnot.toPath(),
                    new File(defaultDir, "annotations.yaml").toPath(),
                    StandardCopyOption.REPLACE_EXISTING);
                System.out.println("  Copied " + legacyAnnot.getPath() + " → projects/default/annotations.yaml");
            }

            // 4. Copy output/data/ → projects/default/data/
            File legacyData = new File("output/data");
            File newData = new File(defaultDir, "data");
            newData.mkdirs();
            File[] dataFiles = legacyData.listFiles();
            if (dataFiles != null) {
                int count = 0;
                for (File f : dataFiles) {
                    if (f.isFile()) {
                        Files.copy(f.toPath(), new File(newData, f.getName()).toPath(),
                            StandardCopyOption.REPLACE_EXISTING);
                        count++;
                    }
                }
                System.out.printf("  Copied %d files from output/data/ → projects/default/data/%n", count);
            }

            // 5. Generate project.json from existing manifest (no pipeline rerun)
            Config config = Config.loadFromProject(defaultDir.getAbsolutePath());
            String manifestJson = new String(Files.readAllBytes(legacyManifest.toPath()), "UTF-8");
            int lociCount = extractManifestLociCount(manifestJson);

            ProjectMetadata pm = new ProjectMetadata();
            pm.id   = "default";
            pm.name = "Default Project";
            pm.description = "Migrated from legacy single-project layout";
            pm.lociCount = lociCount;

            // Count SNPs from locus files
            Set<String> uniqueSnps = new HashSet<>();
            for (int i = 1; i <= lociCount; i++) {
                File locusFile = new File(newData, "locus_" + i + ".json");
                if (locusFile.exists()) {
                    String lj = new String(Files.readAllBytes(locusFile.toPath()), "UTF-8");
                    countSnpsFromLocusJson(lj, uniqueSnps);
                }
            }
            pm.totalSnps = uniqueSnps.size();

            String annotPath = new File(defaultDir, "annotations.yaml").getAbsolutePath();
            pm.countAnnotationSources(annotPath);
            pm.coreInputFingerprint  = ProjectMetadata.computeCoreInputFingerprint(config, defaultDir.getAbsolutePath());
            pm.annotationFingerprint = ProjectMetadata.computeAnnotationFingerprint(annotPath);
            pm.save(defaultDir.getAbsolutePath());

            System.out.printf("  Generated project.json: %d loci, %d SNPs%n", pm.lociCount, pm.totalSnps);
            System.out.println("  Migration complete. Project available as 'default'.");
            System.out.println("═══════════════════════════════════════════════════\n");

        } catch (Exception e) {
            System.err.println("[ERROR] Migration failed: " + e.getMessage());
            e.printStackTrace(System.err);
            // Clean up partial migration so it can be retried
            if (defaultDir.exists() && new File(defaultDir, "config.properties").exists()
                    && !new File(defaultDir, "project.json").exists()) {
                System.err.println("[WARN] Cleaning up partial migration in projects/default/");
                deleteDir(defaultDir);
            }
        }
    }

    private static int extractManifestLociCount(String json) {
        String key = "\"total_loci\":";
        int i = json.indexOf(key);
        if (i < 0) return 0;
        i += key.length();
        int end = i;
        while (end < json.length() && Character.isDigit(json.charAt(end))) end++;
        try { return Integer.parseInt(json.substring(i, end)); }
        catch (NumberFormatException e) { return 0; }
    }

    private static void countSnpsFromLocusJson(String json, Set<String> uniqueSnps) {
        // Extract chr:pos pairs from gwas_snps array entries
        int idx = 0;
        while (true) {
            int chrPos = json.indexOf("\"chr\":\"", idx);
            if (chrPos < 0) break;
            chrPos += 7;
            int chrEnd = json.indexOf('"', chrPos);
            if (chrEnd < 0) break;
            String chr = json.substring(chrPos, chrEnd);

            int posPos = json.indexOf("\"pos\":", chrEnd);
            if (posPos < 0 || posPos - chrEnd > 50) { idx = chrEnd + 1; continue; }
            posPos += 6;
            int posEnd = posPos;
            while (posEnd < json.length() && (Character.isDigit(json.charAt(posEnd)) || json.charAt(posEnd) == '-'))
                posEnd++;
            if (posEnd > posPos) {
                uniqueSnps.add(chr + ":" + json.substring(posPos, posEnd));
            }
            idx = posEnd;
        }
    }

    private static void deleteDir(File dir) {
        File[] files = dir.listFiles();
        if (files != null) for (File f : files) {
            if (f.isDirectory()) deleteDir(f);
            else f.delete();
        }
        dir.delete();
    }

    // ── Legacy single-project mode ───────────────────────────────────────

    private static void runLegacySingleProject(String[] args) throws Exception {
        if (args.length == 0 && new File("config.properties").exists())
            args = new String[]{"--config", "config.properties"};

        Config config = Config.load(args);
        config.validate();

        // Copy viewer HTML into output directory
        new File(config.outputDir + "/data").mkdirs();
        new File(config.outputDir + "/plots").mkdirs();
        new File(config.outputDir + "/tmp").mkdirs();
        Path viewerSrc = Paths.get("index.html");
        if (Files.exists(viewerSrc)) {
            Files.copy(viewerSrc, Paths.get(config.outputDir, "index.html"),
                StandardCopyOption.REPLACE_EXISTING);
        }

        // Start server
        LocalServer server = null;
        try {
            server = new LocalServer(config.outputDir);
            server.start();
            System.out.printf("Viewer available at http://localhost:%d/index.html%n%n",
                LocalServer.PORT);
        } catch (Exception e) {
            System.err.println("[WARN] Could not start local server: " + e.getMessage());
        }

        long t0 = System.currentTimeMillis();
        ProgressTracker progress = LocalServer.progress;
        PipelineResult result = runPipeline(config, config.outputDir, progress);
        progress.done = true;

        long elapsed = System.currentTimeMillis() - t0;
        System.out.printf("%n=== Done in %.1f s ===%n", elapsed / 1000.0);
        printLociTable(result.outputs);

        if (server != null) {
            server.setPipelineState(result.config, result.gff, result.loci, result.outputs);
            System.out.printf("%nPress Ctrl+C to stop the viewer server.%n");
        } else {
            System.out.printf("Open: %s/index.html%n", config.outputDir);
        }
    }

    // ── Helpers ──────────────────────────────────────────────────────────

    private static void printLociTable(List<LocusOutput> outputs) {
        long noSnps  = outputs.stream().filter(o -> o.gwasSnps.isEmpty()).count();
        long noGenes = outputs.stream().filter(o -> o.genes.isEmpty()).count();
        long hasLd   = outputs.stream().filter(o -> o.ldTriangle != null).count();
        System.out.printf("%d loci processed", outputs.size());
        if (noSnps  > 0) System.out.printf(", %d with no SNPs", noSnps);
        if (noGenes > 0) System.out.printf(", %d with no genes", noGenes);
        if (hasLd   > 0) System.out.printf(", %d with LD triangle", hasLd);
        System.out.println("\n");

        System.out.printf("%-6s %-24s %-16s %-12s %s%n",
            "Locus", "Region", "Top SNP", "P-value", "Nearest genes");
        System.out.println("-".repeat(82));
        for (LocusOutput lo : outputs) {
            String region = "chr" + lo.chr + ":" + lo.start + "-" + lo.end;
            String topId  = lo.topSnp != null ? lo.topSnp.id : "—";
            String pval   = lo.topSnp != null ? String.format("%.2e", lo.topSnp.pvalue) : "—";
            String genes  = String.join(", ", lo.nearestGenes);
            System.out.printf("%-6d %-24s %-16s %-12s %s%n",
                lo.locusIndex, region, topId, pval, genes);
        }
    }

    private static String reasonLabel(ProjectMetadata.StaleReason reason) {
        if (reason == null) return "unknown";
        switch (reason) {
            case NO_METADATA:         return "never processed";
            case VERSION_CHANGED:     return "pipeline version changed";
            case CORE_INPUT_CHANGED:  return "input files changed";
            case ANNOTATION_CHANGED:  return "annotations changed";
            default:                  return reason.name();
        }
    }
}
