
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

    // MHC/extended-HLA region (chr6) has extreme, unusual LD structure that breaks the statistical
    // assumptions behind fine-mapping/colocalization (SuSiE, FINEMAP, COJO, coloc, SuSiEx) — results
    // there aren't trustworthy without special handling, and its LD density can also OOM naive
    // LD-matrix computation. We don't block analysis (a user may deliberately want to look here),
    // just surface a clear warning wherever this locus is analyzed.
    private static final long MHC_GRCH37_START = 28_477_897L, MHC_GRCH37_END = 33_448_354L;
    private static final long MHC_GRCH38_START = 28_510_120L, MHC_GRCH38_END = 33_480_577L;

    public static boolean overlapsMhc(Locus locus, String genomeBuild) {
        String chr = locus.chr.replaceFirst("^chr", "");
        if (!chr.equals("6")) return false;
        boolean grch38 = "GRCh38".equalsIgnoreCase(genomeBuild);
        long mhcStart = grch38 ? MHC_GRCH38_START : MHC_GRCH37_START;
        long mhcEnd   = grch38 ? MHC_GRCH38_END   : MHC_GRCH37_END;
        return locus.start <= mhcEnd && locus.end >= mhcStart;
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
        root.mkdirs();
        File baseDir       = new File(root, "base");
        File matchedDir    = new File(root, "matched");
        File harmonizedDir = new File(root, "harmonized");
        File ldDir         = new File(root, "ld");

        boolean allCached = true;
        List<String> log = new ArrayList<>();
        log.add(String.format("=== Base pipeline: Locus %d (id=%s) ===", locus.index, locus.id));
        log.add(String.format("Region: chr%s:%d-%d (padded %d-%d)", locus.chr, locus.start, locus.end, locus.paddedStart, locus.paddedEnd));
        log.add(String.format("Ref panel: %s (%s)", config.refPanelPath, config.refPanelPopulation));
        log.add(String.format("LD window: %d SNPs each side", ldWindow));
        if (overlapsMhc(locus, config.genomeBuild)) {
            log.add("WARNING: This locus overlaps the MHC/extended-HLA region. Its unusual LD structure");
            log.add("  breaks the statistical assumptions behind fine-mapping/colocalization (SuSiE, FINEMAP,");
            log.add("  COJO, coloc, SuSiEx) -- treat PIPs/credible sets/coloc results here with caution.");
        }
        log.add("");

        // Step 1: Extract
        log.add("[Step 1] GWAS extraction...");
        LocusGwasExtractor.Result ext = LocusGwasExtractor.run(config, locus, baseDir);
        if (!ext.ok) {
            pr.error = "Extract failed: " + ext.error;
            log.add("  FAILED: " + ext.error);
            writeLog(root, log);
            return pr;
        }
        pr.gwasSnps = ext.gwasSnpCount;
        log.add(String.format("  OK: %d GWAS SNPs extracted", ext.gwasSnpCount));

        // Step 2: Match (requires ref panel)
        if (config.refPanelPath.isEmpty()) {
            pr.error = "No reference panel configured";
            log.add("  FAILED: No reference panel configured");
            writeLog(root, log);
            return pr;
        }

        log.add("[Step 2] SNP matching...");
        SnpMatcher.Result match = SnpMatcher.run(baseDir, matchedDir, config);
        if (!match.ok) {
            pr.error = "Match failed: " + match.error;
            log.add("  FAILED: " + match.error);
            writeLog(root, log);
            return pr;
        }
        pr.matchedSnps = match.matchedCount;
        log.add(String.format("  OK: %d SNPs matched to ref panel", match.matchedCount));

        // Step 2.5: Compute ref panel allele frequencies
        File refFreqFile = new File(matchedDir, "ref_freq.tsv");
        if (!refFreqFile.exists()) {
            log.add("[Step 2.5] Computing ref panel allele frequencies...");
            String plinkBin = PlinkSubsetter.findPlink(config);
            if (plinkBin != null) {
                try {
                    String bfile = new File(matchedDir, "matched_ref").getAbsolutePath();
                    String freqPrefix = new File(matchedDir, "ref_freq_tmp").getAbsolutePath();
                    ProcessBuilder pb = new ProcessBuilder(plinkBin,
                        "--bfile", bfile, "--freq", "--out", freqPrefix, "--silent");
                    pb.redirectErrorStream(true);
                    Process proc = pb.start();
                    proc.getInputStream().transferTo(java.io.OutputStream.nullOutputStream());
                    proc.waitFor(120, java.util.concurrent.TimeUnit.SECONDS);

                    File plinkFrq = new File(freqPrefix + ".frq");
                    if (plinkFrq.exists()) {
                        // Convert PLINK .frq to clean TSV: snp_id \t a1 \t freq
                        int freqCount = 0;
                        try (BufferedReader fbr = new BufferedReader(new FileReader(plinkFrq));
                             PrintWriter fpw = new PrintWriter(new BufferedWriter(new FileWriter(refFreqFile)))) {
                            fpw.println("snp_id\ta1\tfreq");
                            fbr.readLine(); // skip header
                            String fl;
                            while ((fl = fbr.readLine()) != null) {
                                String[] ff = fl.trim().split("\\s+");
                                if (ff.length >= 5) {
                                    fpw.printf("%s\t%s\t%s%n", ff[1], ff[2], ff[4]);
                                    freqCount++;
                                }
                            }
                        }
                        plinkFrq.delete();
                        new File(freqPrefix + ".log").delete();
                        new File(freqPrefix + ".nosex").delete();
                        log.add(String.format("  OK: %d SNP frequencies computed from ref panel", freqCount));
                    } else {
                        log.add("  Skipped: PLINK --freq produced no output");
                    }
                } catch (Exception e) {
                    log.add("  Skipped: " + e.getMessage());
                }
            } else {
                log.add("  Skipped: PLINK not found");
            }
        } else {
            log.add("[Step 2.5] Ref panel frequencies already computed.");
        }

        // Step 3: Harmonize
        log.add("[Step 3] Allele harmonization...");
        AlleleHarmonizer.Result harm = AlleleHarmonizer.run(matchedDir, harmonizedDir);
        if (!harm.ok) {
            pr.error = "Harmonize failed: " + harm.error;
            log.add("  FAILED: " + harm.error);
            writeLog(root, log);
            return pr;
        }
        pr.harmonizedSnps = harm.kept + harm.flipped + harm.complemented;
        log.add(String.format("  OK: %d kept, %d flipped, %d complemented, %d dropped (%d palindromic/unresolvable)",
            harm.kept, harm.flipped, harm.complemented, harm.dropped, harm.droppedPalindromic));

        // Step 4: LD
        log.add("[Step 4] LD matrix computation...");
        LdMatrixComputer.Result ld = LdMatrixComputer.run(matchedDir, harmonizedDir, ldDir, config, ldWindow);
        if (!ld.ok) {
            pr.error = "LD computation failed: " + ld.error;
            log.add("  FAILED: " + ld.error);
            writeLog(root, log);
            return pr;
        }
        pr.ldSnps = ld.snpCount;
        log.add(String.format("  OK: %d SNPs in LD matrix", ld.snpCount));

        // Step 5: LD-GWAS consistency diagnostic
        log.add("[Step 5] LD-GWAS consistency diagnostic...");
        try {
            LdGwasDiagnostic.DiagnosticResult diag = LdGwasDiagnostic.run(harmonizedDir, ldDir);
            if (diag.ok) {
                pr.diagnosticVerdict = diag.verdict;
                pr.diagnosticFlagged = diag.flaggedSnps;
                log.add(String.format("  Verdict: %s (%d/%d SNPs flagged)",
                    diag.verdict, diag.flaggedSnps, diag.totalSnps));
                if ("high_warn".equals(diag.verdict)) {
                    log.add("  NOTE: High warning — >5% of SNPs show LD-GWAS inconsistency.");
                    log.add("  This typically indicates LD reference panel population mismatch.");
                    log.add("  Fine-mapping results (COJO/SuSiE) may be less reliable for this locus.");
                }
            } else if (diag.error != null) {
                log.add("  Skipped: " + diag.error);
            }
        } catch (Exception e) {
            log.add("  Error (non-fatal): " + e.getMessage());
            System.err.printf("[BaseStepPipeline] Diagnostic failed (non-fatal): %s%n", e.getMessage());
        }

        log.add("");
        log.add("Pipeline complete.");
        pr.ok = true;
        pr.allCached = allCached;
        writeLog(root, log);
        return pr;
    }

    private static void writeLog(File root, List<String> lines) {
        try (PrintWriter pw = new PrintWriter(new BufferedWriter(
                new FileWriter(new File(root, "build.log"))))) {
            for (String line : lines) pw.println(line);
        } catch (IOException e) {
            System.err.printf("[BaseStepPipeline] Failed to write build.log: %s%n", e.getMessage());
        }
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
                if (progress != null) progress.completedLoci.add(locus.index);
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
                if (progress != null) {
                    progress.update("Building base artifacts (" + threads + " threads)",
                        completed, total);
                    progress.completedLoci.add(locus.index);
                }
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
