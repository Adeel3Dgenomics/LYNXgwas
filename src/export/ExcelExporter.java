package export;

import java.io.*;
import java.nio.file.*;
import java.time.*;
import java.time.format.*;
import java.util.*;

/**
 * Builds Excel workbooks by iterating registered column providers.
 * No hardcoded columns — everything flows through ExportRegistry.
 */
public class ExcelExporter {

    public static class ExportProgress {
        public volatile int stepIndex = 0, totalSteps = 6;
        public volatile String currentStep = "Initializing";
        public volatile boolean done = false, cancelled = false;
        public volatile String error = null;
        public volatile int currentLocus = 0, totalLoci = 0;

        public int pct() {
            if (done) return 100;
            double pw = 100.0 / totalSteps;
            double base = Math.max(0, stepIndex - 1) * pw;
            double within = totalLoci > 0 ? currentLocus * pw / totalLoci : 0;
            return Math.min(99, (int)(base + within));
        }

        public String toJson() {
            return String.format("{\"pct\":%d,\"step_index\":%d,\"total_steps\":%d," +
                "\"current_step\":\"%s\",\"done\":%s,\"error\":%s," +
                "\"current_locus\":%d,\"total_loci\":%d}",
                pct(), stepIndex, totalSteps, esc(currentStep), done,
                error != null ? "\"" + esc(error) + "\"" : "null",
                currentLocus, totalLoci);
        }
        private String esc(String s) { return s == null ? "" : s.replace("\\","\\\\").replace("\"","\\\""); }
    }

    /**
     * Export the whole dataset (all loci).
     */
    public static void exportDataset(String projectDir, String outputPath,
                                     ExportProgress progress) {
        try {
            progress.stepIndex = 1;
            progress.currentStep = "Loading project & manifest";

            String dataDir = projectDir + "/data";
            File manifestFile = new File(dataDir, "manifest.json");
            if (!manifestFile.exists()) throw new IOException("No manifest.json — process the project first");

            String manifestJson = new String(Files.readAllBytes(manifestFile.toPath()), "UTF-8");
            List<Integer> locusIndices = parseLocusIndices(manifestJson);
            progress.totalLoci = locusIndices.size();

            // Collect providers
            List<SnpColumnProvider> snpProvs = new ArrayList<>();
            for (SnpColumnProvider p : ExportRegistry.snpProviders())
                if (p.isAvailable(projectDir)) snpProvs.add(p);
            List<LocusColumnProvider> locProvs = new ArrayList<>();
            for (LocusColumnProvider p : ExportRegistry.locusProviders())
                if (p.isAvailable(projectDir)) locProvs.add(p);

            // Step 2-4: Collect data from locus JSONs
            progress.stepIndex = 2;
            progress.currentStep = "Collecting locus data";

            List<LocusContext> lociCtx = new ArrayList<>();
            List<SnpContext> allSnps = new ArrayList<>();

            for (int li = 0; li < locusIndices.size(); li++) {
                if (progress.cancelled) return;
                progress.currentLocus = li + 1;
                int locIdx = locusIndices.get(li);

                File locusFile = new File(dataDir, "locus_" + locIdx + ".json");
                if (!locusFile.exists()) continue;
                String json = new String(Files.readAllBytes(locusFile.toPath()), "UTF-8");

                LocusContext lc = parseLocusContext(json, locIdx);
                lociCtx.add(lc);

                List<SnpContext> snps = parseSnpContexts(json, lc);
                allSnps.addAll(snps);
            }

            progress.stepIndex = 5;
            progress.currentStep = "Computing summary";

            // Compute IMD (inter-marker distance) for each locus
            for (int i = 0; i < lociCtx.size(); i++) {
                LocusContext cur = lociCtx.get(i);
                if (i == 0) {
                    cur.imd = "start";
                } else {
                    LocusContext prev = lociCtx.get(i - 1);
                    if (!cur.chr.equals(prev.chr)) {
                        cur.imd = "start";
                    } else {
                        long dist = cur.start - prev.end;
                        cur.imd = String.valueOf(Math.max(0, dist));
                    }
                }
            }

            // Summary stats
            int totalSnps = allSnps.size();
            int nP5e3 = 0, nP5e5 = 0, nP5e8 = 0;
            for (SnpContext s : allSnps) {
                if (s.pvalue <= 5e-3) nP5e3++;
                if (s.pvalue <= 5e-5) nP5e5++;
                if (s.pvalue <= 5e-8) nP5e8++;
            }

            // Step 6: Write workbook
            progress.stepIndex = 6;
            progress.currentStep = "Writing workbook";

            new File(outputPath).getParentFile().mkdirs();
            try (XlsxWriter xlsx = new XlsxWriter(new FileOutputStream(outputPath))) {
                // Sheet 1: Dataset Summary
                writeSummarySheet(xlsx, projectDir, lociCtx.size(), totalSnps, nP5e3, nP5e5, nP5e8);

                // Sheet 2: Loci Summary
                writeLociSheet(xlsx, lociCtx, locProvs);

                // Sheet 3: SNPs
                writeSnpsSheet(xlsx, allSnps, snpProvs);

                xlsx.finish();
            }

            progress.done = true;
            progress.currentStep = "Complete";
            System.out.printf("[ExcelExporter] Wrote %s (%d loci, %,d SNPs)%n",
                outputPath, lociCtx.size(), totalSnps);

        } catch (Exception e) {
            progress.error = e.getMessage();
            progress.done = true;
            System.err.println("[ExcelExporter] Failed: " + e.getMessage());
        }
    }

    /**
     * Export a single locus.
     */
    public static void exportLocus(String projectDir, int locusIndex,
                                   String outputPath, ExportProgress progress) {
        try {
            progress.totalSteps = 3;
            progress.stepIndex = 1;
            progress.currentStep = "Loading locus data";

            String dataDir = projectDir + "/data";
            File locusFile = new File(dataDir, "locus_" + locusIndex + ".json");
            if (!locusFile.exists()) throw new IOException("Locus file not found: " + locusFile);
            String json = new String(Files.readAllBytes(locusFile.toPath()), "UTF-8");

            List<SnpColumnProvider> snpProvs = new ArrayList<>();
            for (SnpColumnProvider p : ExportRegistry.snpProviders())
                if (p.isAvailable(projectDir)) snpProvs.add(p);
            List<LocusColumnProvider> locProvs = new ArrayList<>();
            for (LocusColumnProvider p : ExportRegistry.locusProviders())
                if (p.isAvailable(projectDir)) locProvs.add(p);

            LocusContext lc = parseLocusContext(json, locusIndex);
            List<SnpContext> snps = parseSnpContexts(json, lc);

            progress.stepIndex = 2;
            progress.currentStep = "Writing workbook";

            new File(outputPath).getParentFile().mkdirs();
            try (XlsxWriter xlsx = new XlsxWriter(new FileOutputStream(outputPath))) {
                writeLociSheet(xlsx, Collections.singletonList(lc), locProvs);
                writeSnpsSheet(xlsx, snps, snpProvs);
                xlsx.finish();
            }

            progress.done = true;
            progress.currentStep = "Complete";
            System.out.printf("[ExcelExporter] Locus %d → %s (%d SNPs)%n",
                locusIndex, outputPath, snps.size());

        } catch (Exception e) {
            progress.error = e.getMessage();
            progress.done = true;
        }
    }

    // ── Sheet writers ────────────────────────────────────────────

    private static void writeSummarySheet(XlsxWriter xlsx, String projectDir,
                                          int nLoci, int nSnps, int nP3, int nP5, int nP8) {
        XlsxWriter.Sheet s = xlsx.addSheet("Dataset Summary");
        String projId = new File(projectDir).getName();
        String date = Instant.now().atOffset(ZoneOffset.UTC).format(DateTimeFormatter.ISO_LOCAL_DATE);

        s.addRow("Project", projId);
        s.addRow("Export Date", date);
        s.addRow("Number of Loci", nLoci);
        s.addRow("Total SNPs in Loci", nSnps);
        s.addRow("SNPs with p < 5e-3", nP3);
        s.addRow("SNPs with p < 5e-5", nP5);
        s.addRow("SNPs with p < 5e-8", nP8);
    }

    private static void writeLociSheet(XlsxWriter xlsx, List<LocusContext> loci,
                                       List<LocusColumnProvider> providers) {
        XlsxWriter.Sheet s = xlsx.addSheet("Loci Summary").freezeHeader().autoFilter();
        // Header row
        List<ColumnSpec> allCols = new ArrayList<>();
        for (LocusColumnProvider p : providers)
            allCols.addAll(p.columns());
        Object[] header = new Object[allCols.size()];
        for (int i = 0; i < allCols.size(); i++) header[i] = allCols.get(i).header;
        s.addRow(header);
        // Data rows
        for (LocusContext lc : loci) {
            Object[] row = new Object[allCols.size()];
            int ci = 0;
            for (LocusColumnProvider p : providers)
                for (ColumnSpec col : p.columns())
                    row[ci++] = p.value(lc, col);
            s.addRow(row);
        }
    }

    private static void writeSnpsSheet(XlsxWriter xlsx, List<SnpContext> snps,
                                       List<SnpColumnProvider> providers) {
        XlsxWriter.Sheet s = xlsx.addSheet("SNPs").freezeHeader().autoFilter();
        List<ColumnSpec> allCols = new ArrayList<>();
        for (SnpColumnProvider p : providers)
            allCols.addAll(p.columns());
        Object[] header = new Object[allCols.size()];
        for (int i = 0; i < allCols.size(); i++) header[i] = allCols.get(i).header;
        s.addRow(header);
        for (SnpContext snp : snps) {
            Object[] row = new Object[allCols.size()];
            int ci = 0;
            for (SnpColumnProvider p : providers)
                for (ColumnSpec col : p.columns())
                    row[ci++] = p.value(snp, col);
            s.addRow(row);
        }
    }

    // ── JSON parsing helpers ─────────────────────────────────────

    static List<Integer> parseLocusIndices(String manifest) {
        List<Integer> indices = new ArrayList<>();
        int idx = 0;
        while (true) {
            int pos = manifest.indexOf("\"index\":", idx);
            if (pos < 0) break;
            pos += 8;
            int end = pos;
            while (end < manifest.length() && Character.isDigit(manifest.charAt(end))) end++;
            if (end > pos) indices.add(Integer.parseInt(manifest.substring(pos, end)));
            idx = end;
        }
        return indices;
    }

    static LocusContext parseLocusContext(String json, int locIdx) {
        LocusContext lc = new LocusContext();
        lc.index = locIdx;
        lc.name = jStr(json, "locus_name");
        lc.chr = jStr(json, "chr");
        lc.start = jLong(json, "start");
        lc.end = jLong(json, "end");
        lc.sizeBp = lc.end - lc.start + 1;
        lc.refPanel = jStr(json, "ref_panel");
        lc.ldComputed = json.contains("\"ld_triangle\"") && !json.contains("\"ld_triangle\":null");

        // Lead SNP
        int tsIdx = json.indexOf("\"top_snp\"");
        if (tsIdx >= 0) {
            int tsEnd = json.indexOf("}", tsIdx);
            if (tsEnd > 0) {
                String ts = json.substring(tsIdx, tsEnd + 1);
                lc.leadSnpId = jStr(ts, "id");
                lc.leadSnpPos = jLong(ts, "pos");
                lc.leadP = jDouble(ts, "pvalue");
                lc.leadBeta = jDoubleOpt(ts, "beta");
                lc.leadOr = jDoubleOpt(ts, "or");
                // Compute OR from beta if not directly present
                if (Double.isNaN(lc.leadOr) && !Double.isNaN(lc.leadBeta)) {
                    lc.leadOr = Math.exp(lc.leadBeta);
                }
            }
        }

        // Nearest genes (locus-level, from manifest)
        int ngIdx = json.indexOf("\"nearest_genes\"");
        if (ngIdx >= 0) {
            int arrStart = json.indexOf('[', ngIdx);
            int arrEnd = json.indexOf(']', arrStart);
            if (arrStart >= 0 && arrEnd > arrStart) {
                String arr = json.substring(arrStart + 1, arrEnd);
                int q1 = arr.indexOf('"');
                if (q1 >= 0) {
                    int q2 = arr.indexOf('"', q1 + 1);
                    if (q2 > q1) lc.nearestGene = arr.substring(q1 + 1, q2);
                }
            }
        }

        // Compute locus nearest gene distance from lead SNP to gene positions
        List<GenePos> genes = parseGenes(json);
        if (!genes.isEmpty() && lc.leadSnpPos > 0) {
            String[] ng = nearestGene(lc.leadSnpPos, genes);
            if (!ng[0].isEmpty()) {
                lc.nearestGene = ng[0];
                try { lc.nearestGeneDist = Long.parseLong(ng[1]); }
                catch (NumberFormatException e) {}
            }
        }

        // Count SNPs by p-value tier
        int snpScan = 0;
        while (true) {
            int pIdx = json.indexOf("\"pvalue\":", snpScan);
            if (pIdx < 0) break;
            pIdx += 9;
            int pEnd = pIdx;
            while (pEnd < json.length() && "0123456789.eE+-".indexOf(json.charAt(pEnd)) >= 0) pEnd++;
            if (pEnd > pIdx) {
                try {
                    double p = Double.parseDouble(json.substring(pIdx, pEnd));
                    lc.nSnps++;
                    if (p <= 5e-3) lc.nP5e3++;
                    if (p <= 5e-5) lc.nP5e5++;
                    if (p <= 5e-8) lc.nP5e8++;
                } catch (NumberFormatException ignored) {}
            }
            snpScan = pEnd;
        }
        // Subtract 1 for the top_snp pvalue counted above
        if (json.contains("\"top_snp\"") && lc.nSnps > 0) lc.nSnps--;

        return lc;
    }

    // Simple gene record for nearest-gene computation
    static class GenePos { String name; long start; long end; }

    static List<GenePos> parseGenes(String json) {
        List<GenePos> genes = new ArrayList<>();
        int genesIdx = json.indexOf("\"genes\"");
        if (genesIdx < 0) return genes;
        int arrStart = json.indexOf('[', genesIdx);
        if (arrStart < 0) return genes;

        int idx = arrStart;
        while (true) {
            int nameIdx = json.indexOf("\"gene_name\":\"", idx);
            if (nameIdx < 0) break;
            int nameStart = nameIdx + 13;
            int nameEnd = json.indexOf('"', nameStart);
            if (nameEnd < 0) break;
            String name = json.substring(nameStart, nameEnd);

            // Find the gene-level start/end (not transcript)
            int startIdx = json.indexOf("\"start\":", nameEnd);
            int endIdx = json.indexOf("\"end\":", nameEnd);
            if (startIdx < 0 || endIdx < 0) { idx = nameEnd; continue; }

            long start = jLong(json.substring(startIdx, Math.min(startIdx + 30, json.length())), "start");
            long end = jLong(json.substring(endIdx, Math.min(endIdx + 30, json.length())), "end");

            if (start > 0 && end > 0) {
                GenePos gp = new GenePos();
                gp.name = name; gp.start = start; gp.end = end;
                genes.add(gp);
            }
            idx = nameEnd + 1;
        }
        return genes;
    }

    static String[] nearestGene(long pos, List<GenePos> genes) {
        String bestName = "";
        long bestDist = Long.MAX_VALUE;
        for (GenePos g : genes) {
            long dist;
            if (pos >= g.start && pos <= g.end) dist = 0;
            else if (pos < g.start) dist = g.start - pos;
            else dist = pos - g.end;
            if (dist < bestDist) { bestDist = dist; bestName = g.name; }
        }
        return new String[]{bestName, String.valueOf(bestDist)};
    }

    static List<SnpContext> parseSnpContexts(String json, LocusContext lc) {
        List<SnpContext> snps = new ArrayList<>();
        // Parse genes for per-SNP nearest gene computation
        List<GenePos> genes = parseGenes(json);

        int gwasIdx = json.indexOf("\"gwas_snps\"");
        if (gwasIdx < 0) return snps;
        int arrStart = json.indexOf('[', gwasIdx);
        int arrEnd = findMatchingBracket(json, arrStart);
        if (arrStart < 0 || arrEnd < 0) return snps;

        String arr = json.substring(arrStart + 1, arrEnd);
        int idx = 0;
        while (true) {
            int objStart = arr.indexOf('{', idx);
            if (objStart < 0) break;
            int objEnd = findMatchingBrace(arr, objStart);
            if (objEnd < 0) break;
            String obj = arr.substring(objStart, objEnd + 1);

            SnpContext sc = new SnpContext();
            sc.id = jStr(obj, "id");
            sc.chr = jStr(obj, "chr");
            sc.pos = jLong(obj, "pos");
            sc.pvalue = jDouble(obj, "pvalue");
            sc.negLog10P = jDouble(obj, "neg_log10_p");
            sc.r2WithLead = jDouble(obj, "r2_with_index");
            sc.ea = jStr(obj, "ea");
            sc.nea = jStr(obj, "nea");
            sc.beta = jDoubleOpt(obj, "beta");
            sc.oddsRatio = jDoubleOpt(obj, "or");
            sc.se = jDoubleOpt(obj, "se");
            sc.sampleN = jDoubleOpt(obj, "n");
            sc.maf = jDoubleOpt(obj, "maf");
            sc.info = jDoubleOpt(obj, "info");

            // Compute OR from beta if not directly available
            if (Double.isNaN(sc.oddsRatio) && !Double.isNaN(sc.beta)) {
                sc.oddsRatio = Math.exp(sc.beta);
            }

            sc.locusIndex = lc.index;
            sc.locusName = lc.name;
            sc.locusChr = lc.chr;
            sc.locusStart = lc.start;
            sc.locusEnd = lc.end;
            sc.isLead = lc.leadSnpId != null && lc.leadSnpId.equals(sc.id);
            sc.distanceToLead = Math.abs(sc.pos - lc.leadSnpPos);

            // Per-SNP nearest gene
            if (!genes.isEmpty()) {
                String[] ng = nearestGene(sc.pos, genes);
                sc.nearestGene = ng[0];
                try { sc.nearestGeneDist = Long.parseLong(ng[1]); }
                catch (NumberFormatException e) { sc.nearestGeneDist = Long.MAX_VALUE; }
            }

            snps.add(sc);
            idx = objEnd + 1;
        }

        // Backfill locus lead beta/OR from the matching SNP in gwas_snps
        for (SnpContext sc : snps) {
            if (sc.isLead) {
                if (Double.isNaN(lc.leadBeta) && !Double.isNaN(sc.beta)) lc.leadBeta = sc.beta;
                if (Double.isNaN(lc.leadOr) && !Double.isNaN(sc.oddsRatio)) lc.leadOr = sc.oddsRatio;
                break;
            }
        }

        return snps;
    }

    // ── Minimal JSON helpers ─────────────────────────────────────

    private static String jStr(String json, String key) {
        String search = "\"" + key + "\":\"";
        int i = json.indexOf(search);
        if (i < 0) return "";
        i += search.length();
        int end = json.indexOf('"', i);
        return end > i ? json.substring(i, end) : "";
    }

    private static long jLong(String json, String key) {
        String search = "\"" + key + "\":";
        int i = json.indexOf(search);
        if (i < 0) return 0;
        i += search.length();
        while (i < json.length() && json.charAt(i) == ' ') i++;
        int end = i;
        while (end < json.length() && "0123456789-".indexOf(json.charAt(end)) >= 0) end++;
        if (end == i) return 0;
        try { return Long.parseLong(json.substring(i, end)); }
        catch (NumberFormatException e) { return 0; }
    }

    private static double jDouble(String json, String key) {
        String search = "\"" + key + "\":";
        int i = json.indexOf(search);
        if (i < 0) return Double.NaN;
        i += search.length();
        while (i < json.length() && json.charAt(i) == ' ') i++;
        if (i < json.length() && json.charAt(i) == 'n') return Double.NaN; // null
        int end = i;
        while (end < json.length() && "0123456789.eE+-".indexOf(json.charAt(end)) >= 0) end++;
        if (end == i) return Double.NaN;
        try { return Double.parseDouble(json.substring(i, end)); }
        catch (NumberFormatException e) { return Double.NaN; }
    }

    private static double jDoubleOpt(String json, String key) {
        String search = "\"" + key + "\":";
        int i = json.indexOf(search);
        if (i < 0) return Double.NaN;
        i += search.length();
        while (i < json.length() && json.charAt(i) == ' ') i++;
        if (i < json.length() && json.charAt(i) == 'n') return Double.NaN;
        int end = i;
        while (end < json.length() && "0123456789.eE+-".indexOf(json.charAt(end)) >= 0) end++;
        if (end == i) return Double.NaN;
        try { return Double.parseDouble(json.substring(i, end)); }
        catch (NumberFormatException e) { return Double.NaN; }
    }

    private static int findMatchingBracket(String s, int pos) {
        if (pos < 0 || s.charAt(pos) != '[') return -1;
        int d = 1;
        for (int i = pos + 1; i < s.length(); i++) {
            char c = s.charAt(i);
            if (c == '[') d++; else if (c == ']') { d--; if (d == 0) return i; }
            else if (c == '"') { i++; while (i < s.length() && s.charAt(i) != '"') { if (s.charAt(i) == '\\') i++; i++; } }
        }
        return -1;
    }

    private static int findMatchingBrace(String s, int pos) {
        if (pos < 0 || s.charAt(pos) != '{') return -1;
        int d = 1;
        for (int i = pos + 1; i < s.length(); i++) {
            char c = s.charAt(i);
            if (c == '{') d++; else if (c == '}') { d--; if (d == 0) return i; }
            else if (c == '"') { i++; while (i < s.length() && s.charAt(i) != '"') { if (s.charAt(i) == '\\') i++; i++; } }
        }
        return -1;
    }
}
