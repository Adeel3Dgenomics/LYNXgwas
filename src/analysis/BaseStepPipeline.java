
import java.io.*;
import java.util.*;
import java.util.concurrent.*;

/**
 * Orchestrates the four base steps for a single locus:
 *   1. Extract  → base/
 *   2. Match    → matched/
 *   3. Harmonize → harmonized/
 *   4. LD       → ld/
 *
 * Each step is cached by content hash and only re-runs if inputs changed.
 * Call {@link #runAll} for the full pipeline or individual step methods
 * if partial re-execution is needed.
 */
public class BaseStepPipeline {

    public static class PipelineResult {
        public boolean ok;
        public String error;
        public String locusId;
        public int gwasSnps;
        public int matchedSnps;
        public int harmonizedSnps;
        public int ldSnps;
        public boolean allCached;
        public String diagnosticVerdict;
        public int diagnosticFlagged;

        public Map<String, Object> toMap() {
            Map<String, Object> m = new LinkedHashMap<>();
            m.put("ok", ok);
            if (error != null) m.put("error", error);
            m.put("locus_id", locusId);
            m.put("gwas_snps", gwasSnps);
            m.put("matched_snps", matchedSnps);
            m.put("harmonized_snps", harmonizedSnps);
            m.put("ld_snps", ldSnps);
            m.put("all_cached", allCached);
            if (diagnosticVerdict != null) {
                m.put("diagnostic_verdict", diagnosticVerdict);
                m.put("diagnostic_flagged", diagnosticFlagged);
            }
            return m;
        }
    }

    /**
     * Get the analysis root directory for a locus.
     */
    public static File analysisDir(Config config, Locus locus) {
        return new File(config.outputDir, "loci_analysis/" + locus.id);
    }

    public static File analysisDir(String projectDir, String locusId) {
        return new File(projectDir, "loci_analysis/" + locusId);
    }

    /**
     * Run all four base steps for a single locus (default LD window).
     */
    public static PipelineResult runAll(Config config, Locus locus) throws IOException {
        return runAll(config, locus, LdMatrixComputer.DEFAULT_LD_WINDOW);
    }

    /**
     * Run all four base steps for a single locus with configurable LD window.
     */
    public static PipelineResult runAll(Config config, Locus locus, int ldWindow) throws IOException {
        PipelineResult pr = new PipelineResult();
        pr.locusId = locus.id;

        File root = analysisDir(config, locus);
        File baseDir       = new File(root, "base");
        File matchedDir    = new File(root, "matched");
        File harmonizedDir = new File(root, "harmonized");
        File ldDir         = new File(root, "ld");

        boolean allCached = true;

        // Step 1: Extract
        LocusGwasExtractor.Result ext = LocusGwasExtractor.run(config, locus, baseDir);
        if (!ext.ok) {
            pr.error = "Extract failed: " + ext.error;
            return pr;
        }
        pr.gwasSnps = ext.gwasSnpCount;

        // Step 2: Match (requires ref panel)
        if (config.refPanelPath.isEmpty()) {
            pr.error = "No reference panel configured";
            return pr;
        }

        SnpMatcher.Result match = SnpMatcher.run(baseDir, matchedDir, config);
        if (!match.ok) {
            pr.error = "Match failed: " + match.error;
            return pr;
        }
        pr.matchedSnps = match.matchedCount;

        // Step 3: Harmonize
        AlleleHarmonizer.Result harm = AlleleHarmonizer.run(matchedDir, harmonizedDir);
        if (!harm.ok) {
            pr.error = "Harmonize failed: " + harm.error;
            return pr;
        }
        pr.harmonizedSnps = harm.kept + harm.flipped + harm.complemented;

        // Step 4: LD
        LdMatrixComputer.Result ld = LdMatrixComputer.run(matchedDir, harmonizedDir, ldDir, config, ldWindow);
        if (!ld.ok) {
            pr.error = "LD computation failed: " + ld.error;
            return pr;
        }
        pr.ldSnps = ld.snpCount;

        // Step 5: LD-GWAS consistency diagnostic
        try {
            LdGwasDiagnostic.DiagnosticResult diag = LdGwasDiagnostic.run(harmonizedDir, ldDir);
            if (diag.ok) {
                pr.diagnosticVerdict = diag.verdict;
                pr.diagnosticFlagged = diag.flaggedSnps;
            }
        } catch (Exception e) {
            System.err.printf("[BaseStepPipeline] Diagnostic failed (non-fatal): %s%n", e.getMessage());
        }

        pr.ok = true;
        pr.allCached = allCached;
        return pr;
    }

    /**
     * Run all base steps for every locus in the project (single-threaded, default window).
     */
    public static List<PipelineResult> runAllLoci(Config config, List<Locus> loci,
                                                   ProgressTracker progress) throws IOException {
        return runAllLoci(config, loci, progress, 1, LdMatrixComputer.DEFAULT_LD_WINDOW);
    }

    /**
     * Run all base steps for every locus in the project with parallel threads.
     */
    public static List<PipelineResult> runAllLoci(Config config, List<Locus> loci,
                                                   ProgressTracker progress,
                                                   int threads) throws IOException {
        return runAllLoci(config, loci, progress, threads, LdMatrixComputer.DEFAULT_LD_WINDOW);
    }

    public static List<PipelineResult> runAllLoci(Config config, List<Locus> loci,
                                                   ProgressTracker progress,
                                                   int threads, int ldWindow) throws IOException {
        if (threads <= 1) {
            List<PipelineResult> results = new ArrayList<>();
            for (int i = 0; i < loci.size(); i++) {
                Locus locus = loci.get(i);
                if (progress != null)
                    progress.update("Building base artifacts", i, loci.size());
                results.add(runOneSafe(config, locus, ldWindow));
            }
            if (progress != null) progress.update("Base artifacts complete", loci.size(), loci.size());
            return results;
        }

        ExecutorService pool = Executors.newFixedThreadPool(threads);
        List<Future<PipelineResult>> futures = new ArrayList<>();
        java.util.concurrent.atomic.AtomicInteger done = new java.util.concurrent.atomic.AtomicInteger(0);
        int total = loci.size();

        for (Locus locus : loci) {
            futures.add(pool.submit(() -> {
                PipelineResult pr = runOneSafe(config, locus, ldWindow);
                int completed = done.incrementAndGet();
                if (progress != null)
                    progress.update("Building base artifacts (" + threads + " threads)",
                        completed, total);
                return pr;
            }));
        }

        List<PipelineResult> results = new ArrayList<>();
        for (Future<PipelineResult> f : futures) {
            try {
                results.add(f.get(10, TimeUnit.MINUTES));
            } catch (Exception e) {
                PipelineResult pr = new PipelineResult();
                pr.error = e.getMessage();
                results.add(pr);
            }
        }
        pool.shutdown();
        if (progress != null) progress.update("Base artifacts complete", total, total);
        return results;
    }

    private static PipelineResult runOneSafe(Config config, Locus locus, int ldWindow) {
        try {
            PipelineResult pr = runAll(config, locus, ldWindow);
            if (!pr.ok) {
                System.err.printf("[BaseStepPipeline] Locus %s failed: %s%n",
                    locus.id, pr.error);
            }
            return pr;
        } catch (Exception e) {
            PipelineResult pr = new PipelineResult();
            pr.locusId = locus.id;
            pr.error = e.getMessage();
            System.err.printf("[BaseStepPipeline] Locus %s exception: %s%n",
                locus.id, e.getMessage());
            return pr;
        }
    }

    /**
     * Invalidate all analysis artifacts for a locus (called after mutation).
     */
    public static void invalidate(Config config, Locus locus) {
        File root = analysisDir(config, locus);
        if (root.exists()) {
            deleteRecursive(root);
            System.out.printf("[BaseStepPipeline] Invalidated analysis for locus %s%n", locus.id);
        }
    }

    public static void invalidate(String projectDir, String locusId) {
        File root = analysisDir(projectDir, locusId);
        if (root.exists()) {
            deleteRecursive(root);
            System.out.printf("[BaseStepPipeline] Invalidated analysis for locus %s%n", locusId);
        }
    }

    /**
     * Check if base artifacts are current for a locus.
     */
    public static Map<String, Boolean> checkStatus(Config config, Locus locus) {
        Map<String, Boolean> status = new LinkedHashMap<>();
        File root = analysisDir(config, locus);

        status.put("base", checkStep(new File(root, "base"), "base"));
        status.put("matched", checkStep(new File(root, "matched"), "matched"));
        status.put("harmonized", checkStep(new File(root, "harmonized"), "harmonized"));
        status.put("ld", checkStep(new File(root, "ld"), "ld"));

        return status;
    }

    private static boolean checkStep(File stepDir, String stepName) {
        if (!stepDir.exists()) return false;
        StepManifest m = StepManifest.read(stepDir, stepName);
        return m != null && m.builtAt > 0;
    }

    private static void deleteRecursive(File file) {
        if (file.isDirectory()) {
            File[] children = file.listFiles();
            if (children != null) {
                for (File child : children) deleteRecursive(child);
            }
        }
        file.delete();
    }
}
