import java.io.*;
import java.util.*;
import java.util.concurrent.*;

/**
 * Phase 3 — LD computation via PLINK 1.9, using subset bfiles from PlinkSubsetter.
 *
 * Per locus (in parallel):
 *   Step A  Select 100-SNP window from .bim  → extract file
 *   Step B  LD with index (top SNP vs all)   → r² map by chr:pos
 *   Step C  Pairwise LD square matrix         → ld_triangle data
 */
public class LdCalculator {

    public static class LdResult {
        /** "chr:pos" → r² with the top SNP. Used for Manhattan coloring. */
        public Map<String, Double> r2ByPos  = new HashMap<>();
        /** Pairwise LD triangle (null if unavailable or fewer than 2 SNPs). */
        public LocusOutput.LdTriangle triangle;
        public boolean ldFailed = false;
    }

    public static Map<Integer, LdResult> computeAll(
            List<Locus>                             loci,
            Map<Integer, Snp>                       topSnps,
            Map<Integer, PlinkSubsetter.SubsetResult> subsets,
            String                                  plinkBin,
            Config                                  config,
            ProgressTracker                         progress) {

        Map<Integer, LdResult> results = new ConcurrentHashMap<>();
        ExecutorService pool = Executors.newFixedThreadPool(
            Math.max(1, config.ldParallelJobs));
        List<Future<?>> futures = new ArrayList<>();

        for (Locus locus : loci) {
            Snp topSnp = topSnps.get(locus.index);
            PlinkSubsetter.SubsetResult subset = subsets.get(locus.index);

            if (topSnp == null || subset == null || !subset.ok || subset.topSnpBimId == null) {
                LdResult r = new LdResult(); r.ldFailed = true;
                results.put(locus.index, r); continue;
            }

            futures.add(pool.submit(() -> {
                try {
                    if (progress != null)
                        progress.update("Computing LD", locus.index, loci.size());
                    LdResult r = compute(locus, topSnp, subset, plinkBin, config);
                    results.put(locus.index, r);
                } catch (Exception e) {
                    System.err.printf("[WARN] LD failed locus %d: %s%n",
                        locus.index, e.getMessage());
                    LdResult r = new LdResult(); r.ldFailed = true;
                    results.put(locus.index, r);
                }
                return null;
            }));
        }

        for (Future<?> f : futures) {
            try { f.get(5, TimeUnit.MINUTES); }
            catch (TimeoutException e) { f.cancel(true); System.err.println("[WARN] LD job timed out"); }
            catch (Exception ignored) {}
        }
        pool.shutdown();
        return results;
    }

    // ─────────────────────────────────────────────────────────────────────────

    private static LdResult compute(Locus locus, Snp topSnp,
                                    PlinkSubsetter.SubsetResult subset,
                                    String plinkBin, Config config)
            throws IOException, InterruptedException {

        LdResult result = new LdResult();
        String subPrefix = config.plinkSubsetsDir() + "/locus_" + locus.index;
        String ldDir     = config.ldResultsDir();
        new File(ldDir).mkdirs();

        // ── Step A: select 100-SNP window around top SNP ─────────────────────
        List<String[]> bimRows = readBimSorted(subPrefix + ".bim");
        int topIdx = -1;
        for (int i = 0; i < bimRows.size(); i++) {
            if (bimRows.get(i)[1].equals(subset.topSnpBimId)) { topIdx = i; break; }
        }

        int lo = topIdx < 0 ? 0 : Math.max(0, topIdx - config.ldTriangleBoundary);
        int hi = topIdx < 0
            ? Math.min(bimRows.size() - 1, config.ldTriangleBoundary * 2)
            : Math.min(bimRows.size() - 1, topIdx + config.ldTriangleBoundary);
        List<String[]> window = bimRows.subList(lo, hi + 1);

        // Write extract file for Step C
        String extractFile = subPrefix + "_ld_snps.txt";
        try (PrintWriter pw = new PrintWriter(extractFile)) {
            for (String[] row : window) pw.println(row[1]);
        }

        // ── Step B: LD with index SNP (all locus SNPs) ───────────────────────
        String ldIdxPrefix = ldDir + "/locus_" + locus.index + "_index";
        long   windowKb    = (locus.paddedEnd - locus.paddedStart) / 1000 + 2;
        runPlink(Arrays.asList(
            plinkBin,
            "--bfile",        subPrefix,
            "--r2",
            "--ld-snp",       subset.topSnpBimId,
            "--ld-window",    "999999",
            "--ld-window-kb", String.valueOf(windowKb),
            "--ld-window-r2", "0.0",
            "--out",          ldIdxPrefix,
            "--silent"
        ));

        File idxFile = new File(ldIdxPrefix + ".ld");
        if (idxFile.exists()) {
            // CHR_A  BP_A  SNP_A  CHR_B  BP_B  SNP_B  R2
            try (BufferedReader br = new BufferedReader(new FileReader(idxFile))) {
                br.readLine(); // header
                String line;
                while ((line = br.readLine()) != null) {
                    String[] f = line.trim().split("\\s+");
                    if (f.length < 7) continue;
                    try {
                        double r2 = Double.parseDouble(f[6]);
                        result.r2ByPos.put(f[3] + ":" + f[4], r2);
                    } catch (NumberFormatException ignored) {}
                }
            }
            idxFile.delete();
        }
        cleanupPlink(ldIdxPrefix);

        // ── Step C: pairwise LD square matrix (window only) ─────────────────
        String ldPairPrefix = ldDir + "/locus_" + locus.index + "_pairwise";
        String plinkOut = runPlinkCapture(Arrays.asList(
            plinkBin,
            "--bfile",   subPrefix,
            "--r2",      "square",
            "--extract", extractFile,
            "--out",     ldPairPrefix,
            "--silent"
        ));

        File pairFile = new File(ldPairPrefix + ".ld");
        System.out.printf("[LD] Locus %d pairwise: exists=%b size=%d%n",
            locus.index, pairFile.exists(), pairFile.exists() ? pairFile.length() : -1L);
        if (!plinkOut.trim().isEmpty()) {
            String preview = plinkOut.length() > 400 ? plinkOut.substring(0, 400) : plinkOut;
            System.out.printf("[LD] Locus %d PLINK output: %s%n", locus.index, preview);
        }
        if (pairFile.exists() && pairFile.length() > 0) {
            try (BufferedReader dbr = new BufferedReader(new FileReader(pairFile))) {
                for (int di = 0; di < 3; di++) {
                    String dl = dbr.readLine();
                    if (dl == null) break;
                    String dlPrev = dl.length() > 120 ? dl.substring(0, 120) + "..." : dl;
                    System.out.printf("[LD] Locus %d pairwise line %d: [%s]%n",
                        locus.index, di, dlPrev);
                }
            }
        }

        // Determine PLINK's actual SNP order: read BIM in natural file order, keep only
        // SNPs present in the extract file.  PLINK always outputs in BIM-file order, so
        // this is the authoritative ordering for matrix rows/columns.
        Set<String> extractSet = new LinkedHashSet<>();
        try (BufferedReader er = new BufferedReader(new FileReader(extractFile))) {
            String el;
            while ((el = er.readLine()) != null) {
                el = el.trim();
                if (!el.isEmpty()) extractSet.add(el);
            }
        }
        List<String[]> plinkOrder = new ArrayList<>();
        try (BufferedReader br = new BufferedReader(new FileReader(subPrefix + ".bim"))) {
            String bl;
            while ((bl = br.readLine()) != null) {
                String[] parts = bl.trim().split("\t", -1);
                if (parts.length >= 4 && extractSet.contains(parts[1])) plinkOrder.add(parts);
            }
        }

        // Ordering diagnostic
        System.out.printf("[LD] Locus %d snp ordering: window=%d  plinkBim=%d%n",
            locus.index, window.size(), plinkOrder.size());
        int diagN = Math.min(3, Math.min(window.size(), plinkOrder.size()));
        for (int di = 0; di < diagN; di++) {
            String wId = window.get(di)[1], pId = plinkOrder.get(di)[1];
            System.out.printf("[LD] Locus %d  snp[%d]: window=%s  plink=%s  match=%b%n",
                locus.index, di, wId, pId, wId.equals(pId));
        }
        if (window.size() != plinkOrder.size())
            System.out.printf("[LD] Locus %d SIZE MISMATCH — some SNPs excluded by PLINK QC%n",
                locus.index);

        // Build triangle: parse full ref-panel matrix, then filter to GWAS-only SNPs
        if (pairFile.exists() && plinkOrder.size() >= 2) {
            double[][] fullMatrix = parsePairwiseMatrix(pairFile, plinkOrder.size());

            // Part 0: build set of GWAS SNP positions for this locus
            Set<Long> gwasPositions = new HashSet<>();
            for (Snp gwas : locus.snps) gwasPositions.add(gwas.pos);

            // Identify which plinkOrder indices have a matching GWAS SNP
            List<Integer> keepIdx = new ArrayList<>();
            for (int wi = 0; wi < plinkOrder.size(); wi++) {
                long pos = Long.parseLong(plinkOrder.get(wi)[3].trim());
                if (gwasPositions.contains(pos)) keepIdx.add(wi);
            }

            // Always keep the top SNP even if it somehow lacks a GWAS entry
            if (subset.topSnpBimId != null) {
                for (int wi = 0; wi < plinkOrder.size(); wi++) {
                    if (plinkOrder.get(wi)[1].equals(subset.topSnpBimId)
                            && !keepIdx.contains(wi)) {
                        keepIdx.add(wi);
                        break;
                    }
                }
                Collections.sort(keepIdx);
            }

            int refTotal = plinkOrder.size();
            int refOnly  = refTotal - keepIdx.size();
            int gwasOnly = 0;
            Set<Long> triPositions = new HashSet<>();
            for (int wi : keepIdx)
                triPositions.add(Long.parseLong(plinkOrder.get(wi)[3].trim()));
            for (Snp gwas : locus.snps)
                if (!triPositions.contains(gwas.pos)) gwasOnly++;
            System.out.printf("[LD] Locus %d filtering: refPanel=%d, inBoth=%d, refOnly=%d, gwasOnly=%d%n",
                locus.index, refTotal, keepIdx.size(), refOnly, gwasOnly);

            // Extract subsetted dense matrix (keepIdx × keepIdx)
            int M = keepIdx.size();
            double[][] matrix = new double[M][M];
            for (int i = 0; i < M; i++) {
                for (int j = 0; j < M; j++) {
                    matrix[i][j] = fullMatrix[keepIdx.get(i)][keepIdx.get(j)];
                }
            }

            if (M >= 2) {
                LocusOutput.LdTriangle tri = new LocusOutput.LdTriangle();
                tri.snpCount = M;
                tri.matrix   = matrix;
                tri.snps     = new ArrayList<>();
                for (int ni = 0; ni < M; ni++) {
                    String[] row = plinkOrder.get(keepIdx.get(ni));
                    LocusOutput.LdTriangle.LdSnp s = new LocusOutput.LdTriangle.LdSnp();
                    s.bimId    = row[1];
                    s.pos      = Long.parseLong(row[3].trim());
                    s.isTopSnp = row[1].equals(subset.topSnpBimId);
                    s.rank     = ni;
                    s.id       = s.bimId;
                    for (Snp gwas : locus.snps) {
                        if (gwas.pos == s.pos) { s.id = gwas.id; break; }
                    }
                    if (s.isTopSnp) tri.topSnpRank = ni;
                    tri.snps.add(s);
                }
                result.triangle = tri;
            }
            pairFile.delete();
        }
        cleanupPlink(ldPairPrefix);

        return result;
    }

    // ─────────────────────────────────────────────────────────────────────────

    private static List<String[]> readBimSorted(String path) throws IOException {
        List<String[]> rows = new ArrayList<>();
        File f = new File(path);
        if (!f.exists()) return rows;
        try (BufferedReader br = new BufferedReader(new FileReader(f))) {
            String line;
            while ((line = br.readLine()) != null) {
                if (line.isEmpty()) continue;
                String[] parts = line.split("\t", -1);
                if (parts.length >= 4) rows.add(parts);
            }
        }
        rows.sort(Comparator.comparingLong(r -> {
            try { return Long.parseLong(r[3].trim()); } catch (NumberFormatException e) { return 0L; }
        }));
        return rows;
    }

    private static double[][] parsePairwiseMatrix(File f, int n) throws IOException {
        double[][] m = new double[n][n];
        for (double[] row : m) Arrays.fill(row, Double.NaN);

        List<String> lines = new ArrayList<>();
        try (BufferedReader br = new BufferedReader(new FileReader(f))) {
            String line;
            while ((line = br.readLine()) != null) {
                if (!line.trim().isEmpty()) lines.add(line);
            }
        }
        if (lines.isEmpty()) return m;

        // Auto-detect delimiter: tab takes priority, else split on whitespace
        String firstLine = lines.get(0);
        String sep = firstLine.contains("\t") ? "\t" : "\\s+";

        // Skip header row if first token is non-numeric
        int startLine = 0;
        String[] probe = firstLine.trim().split(sep, 2);
        if (probe.length > 0) {
            try { Double.parseDouble(probe[0].trim()); }
            catch (NumberFormatException e) { startLine = 1; }
        }

        for (int i = startLine; i < lines.size() && (i - startLine) < n; i++) {
            String[] parts = lines.get(i).trim().split(sep, -1);
            int row = i - startLine;
            for (int j = 0; j < parts.length && j < n; j++) {
                String val = parts[j].trim();
                if (val.isEmpty()
                        || val.equalsIgnoreCase("nan")
                        || val.equalsIgnoreCase("na")
                        || val.equals(".")) {
                    m[row][j] = Double.NaN;
                } else {
                    try { m[row][j] = Double.parseDouble(val); }
                    catch (NumberFormatException e) { m[row][j] = Double.NaN; }
                }
            }
        }

        // Guarantee a full symmetric N×N matrix: mirror upper→lower, diagonal = 1.0
        for (int i = 0; i < n; i++) {
            m[i][i] = 1.0;
            for (int j = i + 1; j < n; j++) {
                if (!Double.isNaN(m[i][j])) {
                    m[j][i] = m[i][j];
                } else if (!Double.isNaN(m[j][i])) {
                    m[i][j] = m[j][i];
                }
            }
        }
        return m;
    }

    private static void runPlink(List<String> cmd) throws IOException, InterruptedException {
        Process p = new ProcessBuilder(cmd).redirectErrorStream(true).start();
        try (BufferedReader br = new BufferedReader(new InputStreamReader(p.getInputStream()))) {
            while (br.readLine() != null) {}
        }
        p.waitFor();
    }

    private static String runPlinkCapture(List<String> cmd) throws IOException, InterruptedException {
        Process p = new ProcessBuilder(cmd).redirectErrorStream(true).start();
        StringBuilder sb = new StringBuilder();
        try (BufferedReader br = new BufferedReader(new InputStreamReader(p.getInputStream()))) {
            String line;
            while ((line = br.readLine()) != null) sb.append(line).append('\n');
        }
        p.waitFor();
        return sb.toString();
    }

    private static void cleanupPlink(String prefix) {
        for (String ext : new String[]{".log", ".nosex", ".noperm"})
            new File(prefix + ext).delete();
    }
}
