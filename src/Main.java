import java.io.*;
import java.nio.file.*;
import java.util.*;

public class Main {

    public static void main(String[] args) throws Exception {
        System.out.println("=== LYNXgwas – Locus Analysis and Genomic Explorer ===");

        if (args.length == 0 && new File("config.properties").exists())
            args = new String[]{"--config", "config.properties"};

        Config config = Config.load(args);
        config.validate();

        System.out.printf("GWAS file  : %s%n", config.gwasFile);
        System.out.printf("Loci file  : %s%n", config.lociFile);
        System.out.printf("GFF3 file  : %s%n", config.gff3File);
        System.out.printf("LD enabled : %s%n", config.ldEnabled);
        System.out.printf("Ref panel  : %s (%s)%n",
            config.refPanelPath.isEmpty() ? "(none)" : config.refPanelPath,
            config.refPanelPopulation);
        System.out.printf("Output     : %s%n%n", config.outputDir);

        // Create all output directories
        new File(config.outputDir + "/data").mkdirs();
        new File(config.outputDir + "/plots").mkdirs();
        new File(config.outputDir + "/tmp").mkdirs();

        // Copy viewer HTML into output directory
        Path viewerSrc = Paths.get("index.html");
        if (Files.exists(viewerSrc)) {
            Files.copy(viewerSrc, Paths.get(config.outputDir, "index.html"),
                StandardCopyOption.REPLACE_EXISTING);
        }
        if (config.ldEnabled) {
            new File(config.plinkSubsetsDir()).mkdirs();
            new File(config.ldResultsDir()).mkdirs();
        }

        // Start HTTP server immediately so the loading screen can poll /progress
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

        // ── Phase 1: Parse inputs ──────────────────────────────────────────
        progress.update("Parsing loci", 0, 1);
        List<Locus> loci = LociParser.parse(config);

        progress.update("Parsing GFF3", 0, 1);
        System.out.println("[Phase 1] Parsing GFF3 annotation...");
        GffParser gff = GffParser.parse(config);

        progress.update("Streaming GWAS", 0, loci.size());
        System.out.println("[Phase 1] Streaming GWAS summary statistics...");
        GwasParser.parse(config, loci);

        // Identify top SNPs
        Map<Integer, Snp> topSnps = new LinkedHashMap<>();
        for (Locus locus : loci) {
            Snp top = locus.snps.stream()
                .min(Comparator.comparingDouble(s -> s.pvalue))
                .orElse(null);
            topSnps.put(locus.index, top);
        }

        // Genome-wide skyline (precomputed binned Manhattan overview)
        progress.update("Building genome skyline", 0, 1);
        System.out.println("[Phase 1] Building genome-wide skyline...");
        GenomeSkyline.generate(config);

        // Optional rsid annotation
        if (config.colRsid.isEmpty() && config.topSnpFile.isEmpty())
            new SnpAnnotator().annotateTopSnps(topSnps, config);

        // ── Phase 2: PLINK subset extraction (sequential) ─────────────────
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
            System.out.printf("[Phase 2] Extracting PLINK subsets for %d loci...%n", loci.size());
            subsets = PlinkSubsetter.subsetAll(loci, topSnps, plinkBin, config, progress);
        }

        // ── Phase 3: LD computation (parallel) ────────────────────────────
        Map<Integer, LdCalculator.LdResult> ldResults = new LinkedHashMap<>();
        if (config.ldEnabled) {
            System.out.printf("[Phase 3] Computing LD (%d parallel jobs)...%n",
                config.ldParallelJobs);
            ldResults = LdCalculator.computeAll(
                loci, topSnps, subsets, plinkBin, config, progress);
        }

        // ── Phase 4: Build output + export ────────────────────────────────
        progress.update("Exporting JSON", 0, loci.size());
        System.out.println("[Phase 4] Building locus output...");
        List<LocusOutput> outputs = new ArrayList<>();

        for (int i = 0; i < loci.size(); i++) {
            Locus locus = loci.get(i);
            LocusOutput lo = new LocusOutput();
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

            // Annotate SNPs with r² from LD results
            LdCalculator.LdResult ld = ldResults.get(locus.index);
            if (ld != null && !ld.ldFailed) {
                for (Snp snp : lo.gwasSnps) {
                    Double r2 = ld.r2ByPos.get(snp.chr + ":" + snp.pos);
                    if (r2 != null) snp.r2 = r2;
                }
                lo.ldTriangle = ld.triangle;
            } else if (lo.topSnp != null && !config.ldEnabled) {
                // Tag top SNP with r2=1.0 so it still displays purple
                for (Snp snp : lo.gwasSnps)
                    if (snp.id.equals(lo.topSnp.id)) { snp.r2 = 1.0; break; }
            }

            // Locus context (boundary-to-boundary distances)
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
            progress.update("Exporting JSON", i + 1, loci.size());
        }

        JsonExporter.export(loci, outputs, config);
        progress.done = true;

        long elapsed = System.currentTimeMillis() - t0;
        System.out.printf("%n=== Done in %.1f s ===%n", elapsed / 1000.0);
        printSummary(outputs);

        if (server != null) {
            server.setPipelineState(config, gff, loci, outputs);
            System.out.printf("%nPress Ctrl+C to stop the viewer server.%n");
        } else {
            System.out.printf("Open: %s/index.html%n", config.outputDir);
        }
    }

    private static void printSummary(List<LocusOutput> outputs) {
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
}
