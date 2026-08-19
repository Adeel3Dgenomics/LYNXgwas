import java.io.*;
import java.nio.file.Files;
import java.util.*;

/**
 * Recomputes a single locus with new boundaries: re-streams GWAS, re-runs GFF overlap,
 * PLINK subset, LD, and re-exports the JSON.
 */
public class LocusUpdater {

    public static class UpdateResult {
        public boolean ok;
        public String error;
        public String jsonPath;
    }

    public static UpdateResult update(int locusIndex, long newStart, long newEnd,
                                      Config config, GffParser gff,
                                      List<Locus> allLoci) {
        UpdateResult result = new UpdateResult();
        try {
            // Find the original locus
            Locus original = null;
            int locusPos = -1;
            for (int i = 0; i < allLoci.size(); i++) {
                if (allLoci.get(i).index == locusIndex) {
                    original = allLoci.get(i);
                    locusPos = i;
                    break;
                }
            }
            if (original == null) {
                result.error = "Locus " + locusIndex + " not found";
                return result;
            }

            // Create updated locus with new boundaries, preserving stable ID
            Locus updated = new Locus(original.id, original.index, original.chr, newStart, newEnd, config.locusPadding);

            // Re-stream GWAS for this single locus
            streamGwasForLocus(updated, config);

            // Find top SNP
            Snp topSnp = updated.snps.stream()
                .min(Comparator.comparingDouble(s -> s.pvalue))
                .orElse(null);

            // Downsample if needed
            if (updated.snps.size() > config.maxSnpsPerLocus) {
                List<Snp> significant = new ArrayList<>();
                List<Snp> rest = new ArrayList<>();
                for (Snp s : updated.snps) {
                    if (s.pvalue < 1e-4) significant.add(s);
                    else rest.add(s);
                }
                int remaining = config.maxSnpsPerLocus - significant.size();
                if (remaining > 0 && !rest.isEmpty()) {
                    Collections.shuffle(rest, new Random(42));
                    significant.addAll(rest.subList(0, Math.min(remaining, rest.size())));
                }
                updated.snps.clear();
                updated.snps.addAll(significant);
            }

            // Optional rsid annotation
            if (topSnp != null && config.colRsid.isEmpty() && config.topSnpFile.isEmpty()) {
                Map<Integer, Snp> topMap = new LinkedHashMap<>();
                topMap.put(locusIndex, topSnp);
                new SnpAnnotator().annotateTopSnps(topMap, config);
            }

            // PLINK subset + LD
            PlinkSubsetter.SubsetResult subset = null;
            LdCalculator.LdResult ldResult = null;
            String plinkBin = null;

            if (config.ldEnabled) {
                plinkBin = PlinkSubsetter.findPlink(config);
                if (plinkBin != null) {
                    // Subset
                    new File(config.plinkSubsetsDir()).mkdirs();
                    List<Locus> singleList = Collections.singletonList(updated);
                    Map<Integer, Snp> topMap = new LinkedHashMap<>();
                    topMap.put(locusIndex, topSnp);
                    Map<Integer, PlinkSubsetter.SubsetResult> subsets =
                        PlinkSubsetter.subsetAll(singleList, topMap, plinkBin, config, null);
                    subset = subsets.get(locusIndex);

                    // LD
                    if (subset != null && subset.ok) {
                        new File(config.ldResultsDir()).mkdirs();
                        Map<Integer, LdCalculator.LdResult> ldResults =
                            LdCalculator.computeAll(singleList, topMap, subsets, plinkBin, config, null);
                        ldResult = ldResults.get(locusIndex);
                    }
                }
            }

            // Build output
            LocusOutput lo = new LocusOutput();
            lo.id          = updated.id;
            lo.locusIndex  = updated.index;
            lo.locusName   = "Locus " + updated.index;
            lo.chr         = updated.chr;
            lo.start       = updated.start;
            lo.end         = updated.end;
            lo.paddedStart = updated.paddedStart;
            lo.paddedEnd   = updated.paddedEnd;
            lo.refPanel    = config.ldEnabled ? config.refPanelPopulation : "";
            lo.topSnp      = topSnp;

            lo.genes        = gff.overlapping(updated.chr, updated.paddedStart, updated.paddedEnd);
            lo.nearestGenes = gff.nearestGeneNames(updated.mid(), new ArrayList<>(lo.genes));

            lo.gwasSnps = new ArrayList<>(updated.snps);
            lo.gwasSnps.sort(Comparator.comparingLong(s -> s.pos));

            if (ldResult != null && !ldResult.ldFailed) {
                for (Snp snp : lo.gwasSnps) {
                    Double r2 = ldResult.r2ByPos.get(snp.chr + ":" + snp.pos);
                    if (r2 != null) snp.r2 = r2;
                }
                lo.ldTriangle = ldResult.triangle;
            } else if (lo.topSnp != null && !config.ldEnabled) {
                for (Snp snp : lo.gwasSnps)
                    if (snp.id.equals(lo.topSnp.id)) { snp.r2 = 1.0; break; }
            }

            // Locus context
            lo.locusContext = new LocusOutput.LocusContext();
            if (locusPos > 0) {
                Locus prev = allLoci.get(locusPos - 1);
                long dist = updated.chr.equals(prev.chr)
                    ? Math.abs(updated.start - prev.end) : Long.MAX_VALUE;
                lo.locusContext.prevLocus = new LocusOutput.LocusRef(
                    prev.index, prev.chr, prev.start, prev.end, prev.mid(), dist);
            }
            if (locusPos < allLoci.size() - 1) {
                Locus next = allLoci.get(locusPos + 1);
                long dist = updated.chr.equals(next.chr)
                    ? Math.abs(next.start - updated.end) : Long.MAX_VALUE;
                lo.locusContext.nextLocus = new LocusOutput.LocusRef(
                    next.index, next.chr, next.start, next.end, next.mid(), dist);
            }

            // Export JSON
            String dataDir = config.outputDir + "/data";
            String json = JsonExporter.locusToJson(lo);
            String jsonPath = dataDir + "/locus_" + locusIndex + ".json";
            try (PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(jsonPath)))) {
                pw.print(json);
            }
            try (PrintWriter pw = new PrintWriter(new BufferedWriter(
                    new FileWriter(dataDir + "/locus_" + locusIndex + ".js")))) {
                pw.print("(function(){window.LOCUS_DATA=window.LOCUS_DATA||{};");
                pw.print("window.LOCUS_DATA[" + locusIndex + "]=" + json + ";");
                pw.print("})();");
            }

            // Update the loci list so subsequent updates see the new boundaries
            allLoci.set(locusPos, updated);

            result.ok = true;
            result.jsonPath = jsonPath;
            System.out.printf("[LocusUpdater] Locus %d updated: chr%s:%d-%d → %d SNPs%n",
                locusIndex, updated.chr, newStart, newEnd, lo.gwasSnps.size());

        } catch (Exception e) {
            result.error = e.getMessage();
            e.printStackTrace();
        }
        return result;
    }

    public static UpdateResult create(String chr, long newStart, long newEnd, String locusName,
                                       Config config, GffParser gff,
                                       List<Locus> allLoci, List<LocusOutput> allOutputs) {
        UpdateResult result = new UpdateResult();
        try {
            int newIndex = allLoci.stream().mapToInt(l -> l.index).max().orElse(0) + 1;

            Locus locus = new Locus(newIndex, chr, newStart, newEnd, config.locusPadding);

            // Stream GWAS
            streamGwasForLocus(locus, config);

            // Top SNP
            Snp topSnp = locus.snps.stream()
                .min(Comparator.comparingDouble(s -> s.pvalue))
                .orElse(null);

            // Downsample
            if (locus.snps.size() > config.maxSnpsPerLocus) {
                List<Snp> significant = new ArrayList<>();
                List<Snp> rest = new ArrayList<>();
                for (Snp s : locus.snps) {
                    if (s.pvalue < 1e-4) significant.add(s);
                    else rest.add(s);
                }
                int remaining = config.maxSnpsPerLocus - significant.size();
                if (remaining > 0 && !rest.isEmpty()) {
                    Collections.shuffle(rest, new Random(42));
                    significant.addAll(rest.subList(0, Math.min(remaining, rest.size())));
                }
                locus.snps.clear();
                locus.snps.addAll(significant);
            }

            // rsid annotation
            if (topSnp != null && config.colRsid.isEmpty() && config.topSnpFile.isEmpty()) {
                Map<Integer, Snp> topMap = new LinkedHashMap<>();
                topMap.put(newIndex, topSnp);
                new SnpAnnotator().annotateTopSnps(topMap, config);
            }

            // PLINK + LD
            PlinkSubsetter.SubsetResult subset = null;
            LdCalculator.LdResult ldResult = null;
            if (config.ldEnabled) {
                String plinkBin = PlinkSubsetter.findPlink(config);
                if (plinkBin != null) {
                    new File(config.plinkSubsetsDir()).mkdirs();
                    List<Locus> singleList = Collections.singletonList(locus);
                    Map<Integer, Snp> topMap = new LinkedHashMap<>();
                    topMap.put(newIndex, topSnp);
                    Map<Integer, PlinkSubsetter.SubsetResult> subsets =
                        PlinkSubsetter.subsetAll(singleList, topMap, plinkBin, config, null);
                    subset = subsets.get(newIndex);
                    // subsetAll may have replaced topSnp with a fallback
                    topSnp = topMap.get(newIndex);

                    if (subset != null && subset.ok) {
                        new File(config.ldResultsDir()).mkdirs();
                        Map<Integer, LdCalculator.LdResult> ldResults =
                            LdCalculator.computeAll(singleList, topMap, subsets, plinkBin, config, null);
                        ldResult = ldResults.get(newIndex);
                    }
                }
            }

            // Build output
            LocusOutput lo = new LocusOutput();
            lo.id          = locus.id;
            lo.locusIndex  = newIndex;
            lo.locusName   = locusName != null && !locusName.isEmpty() ? locusName : "Locus " + newIndex;
            lo.chr         = chr;
            lo.start       = newStart;
            lo.end         = newEnd;
            lo.paddedStart = locus.paddedStart;
            lo.paddedEnd   = locus.paddedEnd;
            lo.refPanel    = config.ldEnabled ? config.refPanelPopulation : "";
            lo.topSnp      = topSnp;
            lo.genes       = gff.overlapping(chr, locus.paddedStart, locus.paddedEnd);
            lo.nearestGenes = gff.nearestGeneNames(locus.mid(), new ArrayList<>(lo.genes));
            lo.gwasSnps    = new ArrayList<>(locus.snps);
            lo.gwasSnps.sort(Comparator.comparingLong(s -> s.pos));

            if (ldResult != null && !ldResult.ldFailed) {
                for (Snp snp : lo.gwasSnps) {
                    Double r2 = ldResult.r2ByPos.get(snp.chr + ":" + snp.pos);
                    if (r2 != null) snp.r2 = r2;
                }
                lo.ldTriangle = ldResult.triangle;
            } else if (lo.topSnp != null && !config.ldEnabled) {
                for (Snp snp : lo.gwasSnps)
                    if (snp.id.equals(lo.topSnp.id)) { snp.r2 = 1.0; break; }
            }

            lo.locusContext = new LocusOutput.LocusContext();

            // Export locus JSON
            String dataDir = config.outputDir + "/data";
            String json = JsonExporter.locusToJson(lo);
            String jsonPath = dataDir + "/locus_" + newIndex + ".json";
            try (PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(jsonPath)))) {
                pw.print(json);
            }
            try (PrintWriter pw = new PrintWriter(new BufferedWriter(
                    new FileWriter(dataDir + "/locus_" + newIndex + ".js")))) {
                pw.print("(function(){window.LOCUS_DATA=window.LOCUS_DATA||{};");
                pw.print("window.LOCUS_DATA[" + newIndex + "]=" + json + ";");
                pw.print("})();");
            }

            // Add to lists, sort by chromosome then start, rebuild navigation
            allLoci.add(locus);
            allOutputs.add(lo);
            allLoci.sort(Comparator.comparingInt(Locus::chrInt).thenComparingLong(l -> l.start));
            allOutputs.sort(Comparator.comparingInt((LocusOutput o) -> Locus.chrToInt(o.chr))
                .thenComparingLong(o -> o.start));
            rebuildLocusContext(allOutputs, config);
            ProjectMetadata.syncLociCount(config.outputDir, allOutputs.size());

            result.ok = true;
            result.jsonPath = jsonPath;
            System.out.printf("[LocusUpdater] Created locus %d: chr%s:%d-%d → %d SNPs%n",
                newIndex, chr, newStart, newEnd, lo.gwasSnps.size());

        } catch (Exception e) {
            result.error = e.getMessage();
            e.printStackTrace();
        }
        return result;
    }

    public static class SplitRegion {
        public String name;
        public long start, end;
        public SplitRegion(String name, long start, long end) {
            this.name = name; this.start = start; this.end = end;
        }
    }

    public static UpdateResult split(int origIndex, List<SplitRegion> regions,
                                     Config config, GffParser gff,
                                     List<Locus> allLoci, List<LocusOutput> allOutputs) {
        UpdateResult result = new UpdateResult();
        try {
            // Find and remove original
            Locus original = null;
            int origPos = -1;
            for (int i = 0; i < allLoci.size(); i++) {
                if (allLoci.get(i).index == origIndex) {
                    original = allLoci.get(i);
                    origPos = i;
                    break;
                }
            }
            if (original == null) {
                result.error = "Locus " + origIndex + " not found";
                return result;
            }

            // Create each sub-locus
            List<Integer> newIndices = new ArrayList<>();
            for (SplitRegion r : regions) {
                UpdateResult sub = create(original.chr, r.start, r.end, r.name,
                    config, gff, allLoci, allOutputs);
                if (!sub.ok) {
                    result.error = "Failed creating sub-locus '" + r.name + "': " + sub.error;
                    return result;
                }
                newIndices.add(allLoci.get(allLoci.size() - 1).index);
            }

            // Remove original locus from lists
            allLoci.removeIf(l -> l.index == origIndex);
            allOutputs.removeIf(o -> o.locusIndex == origIndex);

            // Sort by chromosome then start position
            allLoci.sort(Comparator.comparingInt(Locus::chrInt).thenComparingLong(l -> l.start));
            allOutputs.sort(Comparator.comparingInt((LocusOutput o) -> Locus.chrToInt(o.chr))
                .thenComparingLong(o -> o.start));

            // Rebuild prev/next context for all loci and re-export their JSONs
            rebuildLocusContext(allOutputs, config);
            ProjectMetadata.syncLociCount(config.outputDir, allOutputs.size());

            result.ok = true;
            System.out.printf("[LocusUpdater] Split locus %d into %d sub-loci: %s%n",
                origIndex, regions.size(), newIndices);

        } catch (Exception e) {
            result.error = e.getMessage();
            e.printStackTrace();
        }
        return result;
    }

    public static class ValidationResult {
        public boolean valid = true;
        public List<String> warnings = new ArrayList<>();
        public List<LdViolation> ldViolations = new ArrayList<>();
    }

    public static class LdViolation {
        public String snpA, snpB;
        public long posA, posB;
        public int regionA, regionB;
        public double r2;
    }

    public static ValidationResult validateSplit(int origIndex, List<SplitRegion> regions,
                                                  Config config, List<Locus> allLoci) {
        ValidationResult vr = new ValidationResult();

        // Check minimum distance between adjacent regions
        List<SplitRegion> sorted = new ArrayList<>(regions);
        sorted.sort(Comparator.comparingLong(r -> r.start));
        for (int i = 0; i < sorted.size() - 1; i++) {
            long gap = sorted.get(i + 1).start - sorted.get(i).end;
            if (gap < config.splitMinDistBp) {
                vr.warnings.add("Regions '" + sorted.get(i).name + "' and '" +
                    sorted.get(i + 1).name + "' are only " + gap +
                    " bp apart (minimum: " + config.splitMinDistBp + " bp).");
                vr.valid = false;
            }
        }

        // Check overlapping regions
        for (int i = 0; i < sorted.size() - 1; i++) {
            if (sorted.get(i).end > sorted.get(i + 1).start) {
                vr.warnings.add("Regions '" + sorted.get(i).name + "' and '" +
                    sorted.get(i + 1).name + "' overlap.");
                vr.valid = false;
            }
        }

        // Cross-region LD check using PLINK
        if (!config.ldEnabled) return vr;
        String plinkBin = PlinkSubsetter.findPlink(config);
        if (plinkBin == null) return vr;

        // Find original locus to get the PLINK subset prefix
        Locus original = null;
        for (Locus l : allLoci) {
            if (l.index == origIndex) { original = l; break; }
        }
        if (original == null) return vr;

        String subPrefix = config.plinkSubsetsDir() + "/locus_" + origIndex;
        File bimFile = new File(subPrefix + ".bim");
        if (!bimFile.exists()) return vr;

        try {
            // Read BIM to get SNP positions
            Map<Long, String> posToBimId = new LinkedHashMap<>();
            try (BufferedReader br = new BufferedReader(new FileReader(bimFile))) {
                String line;
                while ((line = br.readLine()) != null) {
                    String[] f = line.trim().split("\t", -1);
                    if (f.length < 4) continue;
                    long pos = Long.parseLong(f[3].trim());
                    posToBimId.put(pos, f[1]);
                }
            }

            // For each region, collect SNP BIM IDs within that region
            List<List<String>> regionSnps = new ArrayList<>();
            for (SplitRegion r : regions) {
                List<String> ids = new ArrayList<>();
                for (Map.Entry<Long, String> e : posToBimId.entrySet()) {
                    if (e.getKey() >= r.start && e.getKey() <= r.end) ids.add(e.getValue());
                }
                regionSnps.add(ids);
            }

            // For each pair of regions, pick up to 10 random SNPs from each and compute LD
            Random rng = new Random(42);
            for (int i = 0; i < regions.size(); i++) {
                for (int j = i + 1; j < regions.size(); j++) {
                    List<String> snpsA = regionSnps.get(i);
                    List<String> snpsB = regionSnps.get(j);
                    if (snpsA.isEmpty() || snpsB.isEmpty()) continue;

                    // Sample up to 10 from each
                    List<String> sampleA = new ArrayList<>(snpsA);
                    List<String> sampleB = new ArrayList<>(snpsB);
                    Collections.shuffle(sampleA, rng);
                    Collections.shuffle(sampleB, rng);
                    sampleA = sampleA.subList(0, Math.min(10, sampleA.size()));
                    sampleB = sampleB.subList(0, Math.min(10, sampleB.size()));

                    // Write SNP list for extraction
                    String ldDir = config.ldResultsDir();
                    new File(ldDir).mkdirs();
                    String extractFile = ldDir + "/split_validate_snps.txt";
                    try (PrintWriter pw = new PrintWriter(extractFile)) {
                        for (String s : sampleA) pw.println(s);
                        for (String s : sampleB) pw.println(s);
                    }

                    // Run PLINK --r2
                    String ldPrefix = ldDir + "/split_validate";
                    long windowKb = (original.paddedEnd - original.paddedStart) / 1000 + 2;
                    Process proc = new ProcessBuilder(Arrays.asList(
                        plinkBin,
                        "--bfile", subPrefix,
                        "--r2",
                        "--extract", extractFile,
                        "--ld-window", "999999",
                        "--ld-window-kb", String.valueOf(windowKb),
                        "--ld-window-r2", "0.0",
                        "--out", ldPrefix,
                        "--silent"
                    )).redirectErrorStream(true).start();
                    try (BufferedReader br = new BufferedReader(
                            new InputStreamReader(proc.getInputStream()))) {
                        while (br.readLine() != null) {}
                    }
                    proc.waitFor();

                    // Parse results — look for cross-region pairs with r² > threshold
                    Set<String> setA = new HashSet<>(sampleA);
                    Set<String> setB = new HashSet<>(sampleB);

                    File ldFile = new File(ldPrefix + ".ld");
                    if (ldFile.exists()) {
                        try (BufferedReader br = new BufferedReader(new FileReader(ldFile))) {
                            br.readLine(); // header
                            String line;
                            while ((line = br.readLine()) != null) {
                                String[] f = line.trim().split("\\s+");
                                if (f.length < 7) continue;
                                String snpA = f[2], snpB = f[5];
                                double r2 = Double.parseDouble(f[6]);
                                // Check if cross-region
                                boolean crossRegion = (setA.contains(snpA) && setB.contains(snpB))
                                    || (setA.contains(snpB) && setB.contains(snpA));
                                if (crossRegion && r2 > config.splitLdThreshold) {
                                    LdViolation v = new LdViolation();
                                    v.snpA = snpA; v.snpB = snpB;
                                    v.posA = Long.parseLong(f[1]); v.posB = Long.parseLong(f[4]);
                                    v.regionA = i; v.regionB = j;
                                    v.r2 = r2;
                                    vr.ldViolations.add(v);
                                    vr.valid = false;
                                }
                            }
                        }
                        ldFile.delete();
                    }
                    // Cleanup
                    new File(extractFile).delete();
                    new File(ldPrefix + ".log").delete();
                    new File(ldPrefix + ".nosex").delete();
                }
            }

            if (!vr.ldViolations.isEmpty()) {
                vr.warnings.add(vr.ldViolations.size() + " SNP pair(s) exceed r² threshold of " +
                    config.splitLdThreshold + " across regions.");
            }

        } catch (Exception e) {
            vr.warnings.add("LD validation error: " + e.getMessage());
        }

        return vr;
    }

    private static void rebuildLocusContext(List<LocusOutput> allOutputs, Config config) {
        for (int i = 0; i < allOutputs.size(); i++) {
            LocusOutput lo = allOutputs.get(i);
            lo.locusContext = new LocusOutput.LocusContext();
            if (i > 0) {
                LocusOutput prev = allOutputs.get(i - 1);
                long dist = lo.chr.equals(prev.chr)
                    ? Math.abs(lo.start - prev.end) : Long.MAX_VALUE;
                lo.locusContext.prevLocus = new LocusOutput.LocusRef(
                    prev.locusIndex, prev.chr, prev.start, prev.end,
                    (prev.start + prev.end) / 2, dist);
            }
            if (i < allOutputs.size() - 1) {
                LocusOutput next = allOutputs.get(i + 1);
                long dist = lo.chr.equals(next.chr)
                    ? Math.abs(next.start - lo.end) : Long.MAX_VALUE;
                lo.locusContext.nextLocus = new LocusOutput.LocusRef(
                    next.locusIndex, next.chr, next.start, next.end,
                    (next.start + next.end) / 2, dist);
            }
        }
        // Re-export locus JSONs with updated context + manifest.
        // For skeleton LocusOutputs (lazy-loaded, missing gwasSnps/topSnp/ldTriangle —
        // see LocalServer.ensureProjectState), only patch the locus_context field in the
        // existing JSON to avoid overwriting LD matrix, SNP data, etc. with empty values.
        try {
            String dataDir = config.outputDir + "/data";
            for (LocusOutput lo : allOutputs) {
                boolean hasFullData = lo.topSnp != null || !lo.gwasSnps.isEmpty();
                if (hasFullData) {
                    String json = JsonExporter.locusToJson(lo);
                    try (PrintWriter pw = new PrintWriter(new BufferedWriter(
                            new FileWriter(dataDir + "/locus_" + lo.locusIndex + ".json")))) {
                        pw.print(json);
                    }
                    try (PrintWriter pw = new PrintWriter(new BufferedWriter(
                            new FileWriter(dataDir + "/locus_" + lo.locusIndex + ".js")))) {
                        pw.print("(function(){window.LOCUS_DATA=window.LOCUS_DATA||{};");
                        pw.print("window.LOCUS_DATA[" + lo.locusIndex + "]=" + json + ";");
                        pw.print("})();");
                    }
                } else {
                    patchLocusContextOnDisk(dataDir, lo);
                }
            }
            JsonExporter.exportManifest(allOutputs, config);
        } catch (IOException e) {
            System.err.println("[LocusUpdater] Failed to re-export locus context: " + e.getMessage());
        }
    }

    private static void patchLocusContextOnDisk(String dataDir, LocusOutput lo) throws IOException {
        File jsonFile = new File(dataDir + "/locus_" + lo.locusIndex + ".json");
        if (!jsonFile.exists()) return;

        String existing = new String(Files.readAllBytes(jsonFile.toPath()), "UTF-8");

        // Build the replacement locus_context JSON fragment
        StringBuilder ctx = new StringBuilder("\"locus_context\":");
        if (lo.locusContext != null) {
            ctx.append("{\"prev_locus\":");
            appendLocusRef(ctx, lo.locusContext.prevLocus);
            ctx.append(",\"next_locus\":");
            appendLocusRef(ctx, lo.locusContext.nextLocus);
            ctx.append('}');
        } else {
            ctx.append("null");
        }

        int ctxStart = existing.lastIndexOf("\"locus_context\"");
        if (ctxStart < 0) return;

        String patched = existing.substring(0, ctxStart) + ctx.toString() + "}";
        try (PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(jsonFile)))) {
            pw.print(patched);
        }
        File jsFile = new File(dataDir + "/locus_" + lo.locusIndex + ".js");
        try (PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(jsFile)))) {
            pw.print("(function(){window.LOCUS_DATA=window.LOCUS_DATA||{};");
            pw.print("window.LOCUS_DATA[" + lo.locusIndex + "]=" + patched + ";");
            pw.print("})();");
        }
    }

    private static void appendLocusRef(StringBuilder sb, LocusOutput.LocusRef ref) {
        if (ref == null) { sb.append("null"); return; }
        sb.append("{\"index\":").append(ref.index);
        sb.append(",\"chr\":\"").append(ref.chr).append('"');
        sb.append(",\"start\":").append(ref.start);
        sb.append(",\"end\":").append(ref.end);
        sb.append(",\"mid\":").append(ref.mid);
        sb.append(",\"distance_bp\":").append(ref.distanceBp);
        sb.append('}');
    }

    private static void streamGwasForLocus(Locus locus, Config config) throws IOException {
        try (BufferedReader br = new BufferedReader(new FileReader(config.gwasFile), 1024 * 1024)) {
            String header = br.readLine();
            if (header == null) return;

            String[] cols = header.trim().split("\t");
            String avail = String.join(", ", cols);

            int iChr  = GwasParser.colIdx(cols, config.colChr);
            int iPos  = GwasParser.colIdx(cols, config.colPos);
            int iPval = GwasParser.colIdx(cols, config.colPvalue);
            int iEa   = GwasParser.colIdx(cols, config.colEa);
            int iNea  = GwasParser.colIdx(cols, config.colNea);
            int iVarid = config.colVarid.isEmpty() ? -1 : GwasParser.colIdx(cols, config.colVarid);
            int iRsid  = config.colRsid.isEmpty()  ? -1 : GwasParser.colIdx(cols, config.colRsid);

            int iBeta = config.colBeta.isEmpty() ? -1 : GwasParser.colIdx(cols, config.colBeta);
            int iOr   = config.colOr.isEmpty()   ? -1 : GwasParser.colIdx(cols, config.colOr);
            int iSe   = config.colSe.isEmpty()   ? -1 : GwasParser.colIdx(cols, config.colSe);
            int iN    = config.colN.isEmpty()    ? -1 : GwasParser.colIdx(cols, config.colN);
            int iMaf  = config.colMaf.isEmpty()  ? -1 : GwasParser.colIdx(cols, config.colMaf);
            int iInfo = config.colInfo.isEmpty() ? -1 : GwasParser.colIdx(cols, config.colInfo);

            String line;
            while ((line = br.readLine()) != null) {
                if (line.isEmpty()) continue;
                String[] f = GwasParser.splitTab(line);
                if (f.length <= Math.max(iChr, Math.max(iPos, iPval))) continue;

                String chr = f[iChr].trim();
                if (!chr.equals(locus.chr)) continue;

                long pos;
                double pval;
                try {
                    pos  = Long.parseLong(f[iPos].trim());
                    pval = Double.parseDouble(f[iPval].trim());
                } catch (NumberFormatException e) { continue; }

                if (pos < locus.paddedStart || pos > locus.paddedEnd) continue;

                String varid = (iVarid >= 0 && iVarid < f.length) ? f[iVarid].trim() : chr + ":" + pos;
                String rsid  = (iRsid  >= 0 && iRsid  < f.length) ? f[iRsid].trim()  : "";
                String id    = (!rsid.isEmpty() && !rsid.equals(".")) ? rsid : varid;
                String ea    = (iEa  < f.length) ? f[iEa].trim()  : ".";
                String nea   = (iNea < f.length) ? f[iNea].trim() : ".";

                Snp snp = new Snp(id, chr, pos, pval, ea, nea);
                if (iBeta >= 0) snp.beta      = parseOptionalDouble(f, iBeta);
                if (iOr   >= 0) snp.oddsRatio = parseOptionalDouble(f, iOr);
                if (iSe   >= 0) snp.se        = parseOptionalDouble(f, iSe);
                if (iN    >= 0) snp.sampleN   = parseOptionalDouble(f, iN);
                if (iMaf  >= 0) snp.maf       = parseOptionalDouble(f, iMaf);
                if (iInfo >= 0) snp.infoScore  = parseOptionalDouble(f, iInfo);

                locus.snps.add(snp);
            }
        }
    }

    private static double parseOptionalDouble(String[] fields, int idx) {
        if (idx >= fields.length) return Double.NaN;
        String v = fields[idx].trim();
        if (v.isEmpty() || v.equals(".") || v.equalsIgnoreCase("NA") || v.equalsIgnoreCase("nan"))
            return Double.NaN;
        try { return Double.parseDouble(v); }
        catch (NumberFormatException e) { return Double.NaN; }
    }
}
