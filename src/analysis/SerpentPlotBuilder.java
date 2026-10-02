import java.util.*;

/**
 * Derives the Locus Serpent Plot view from an already-computed {@link MultiLocusResult} (the same
 * Locus Matrix job output Gene Constellation reads), plus the plot's category analysis.
 *
 * Per dataset x locus it reports a significance tier and the lead-SNP offset. The tiers match the
 * plot's color/thickness steps:
 *   0: p >= 5e-3 or missing   1: 5e-5 <= p < 5e-3   2: 5e-8 <= p < 5e-5   3: p < 5e-8
 * The offset is the dataset's best-SNP position minus the position of the locus's overall best SNP
 * (the strongest dataset's lead variant), so a serpent that bends means that dataset peaks elsewhere.
 *
 * Datasets with fewer than MIN_GWS_LOCI genome-wide significant loci carry almost no shared-loci
 * information (they sit at distance ~1 from everything and would dominate any cut of the tree), so
 * they are drawn but left out of the clustering and the category tests, and reported as such.
 * The remaining datasets are clustered on their tier vectors (cosine distance, average linkage / UPGMA). The
 * category analysis then asks which dataset attribute (disease, ancestry, sample-size tier, or any
 * user-supplied category.<name>) best explains that shared-loci structure:
 *   - PERMANOVA (Anderson 2001) on the same distance matrix: R^2, pseudo-F, permutation p-value;
 *   - adjusted Rand index between the dendrogram's K clusters and the category's labels;
 *   - per category level, loci that are consistently genome-wide significant within the level
 *     (>= half of its datasets, at least 2) and specifically so: more often than in the other
 *     datasets with a two-sided Fisher's exact p < 0.05 (uncorrected). A locus significant almost
 *     everywhere (e.g. an extended-MHC LD artifact) is therefore not reported as group-specific.
 *
 * Pure, read-only derived view: no pipeline re-run.
 */
public class SerpentPlotBuilder {

    public static final int DEFAULT_PERMUTATIONS = 999;
    public static final long PERMUTATION_SEED = 20260928L;
    /** Consistent-locus rule: at least this fraction of a level's datasets at p < 5e-8. */
    public static final double CONSISTENT_FRACTION = 0.5;
    public static final int MAX_CONSISTENT_LOCI_PER_LEVEL = 25;
    /** Consistent-locus rule: group-vs-rest two-sided Fisher's exact p below this (uncorrected). */
    public static final double CONSISTENT_FISHER_P = 0.05;
    /** Datasets with fewer genome-wide significant loci than this are left out of clustering/tests. */
    public static final int MIN_GWS_LOCI = 3;

    private SerpentPlotBuilder() {}

    static int tierOf(double p) {
        if (Double.isNaN(p)) return 0;
        if (p < 5e-8) return 3;
        if (p < 5e-5) return 2;
        if (p < 5e-3) return 1;
        return 0;
    }

    // ── result model ─────────────────────────────────────────────────────────────────────────

    public static class Dataset {
        public String id = "", name = "", disease = "", ancestry = "";
        public int sampleN = 0;
        public Map<String, String> categories = new LinkedHashMap<>(); // category name -> level (all categories)
    }

    public static class Cell {
        public int tier = 0;
        public double p = Double.NaN;
        public long pos = 0;
        public long offsetBp = 0;
        public String snp = "";
    }

    public static class Locus {
        public int index;
        public String chr = "", gene = "";
        public long start, end, refPos;
        public int nDatasetsGws = 0;
        public Cell[] cells;
    }

    public static class Merge { public int left, right, size; public double height; }

    public static class ConsistentLocus {
        public int locusIndex; public String gene = "", chr = "";
        public long pos;
        public int nIn, nInSig, nOut, nOutSig;
        public double fisherP;
    }

    public static class Level {
        public String level = "";
        public int nDatasets = 0;
        public double cohesion = Double.NaN;   // 1 - mean within-level distance
        public double separation = Double.NaN; // mean distance to datasets outside the level
        public int nConsistentLoci = 0;
        public List<ConsistentLocus> consistentLoci = new ArrayList<>();
    }

    public static class CategoryAnalysis {
        public String category = "";
        public boolean testable = false;
        public String note = "";
        public int nDatasets = 0, nLevels = 0;
        public double r2 = Double.NaN, pseudoF = Double.NaN, permP = Double.NaN, ari = Double.NaN;
        public int nPermutations = 0;
        public List<Level> levels = new ArrayList<>();
    }

    public static class Result {
        public String name = "", refPanelLabel = "", createdAt = "";
        public List<Dataset> datasets = new ArrayList<>();
        public List<Locus> loci = new ArrayList<>();
        public List<String> categoryNames = new ArrayList<>();
        public int[] leafOrder = new int[0];   // dataset indices: clustered datasets in dendrogram order, then excluded ones
        public int[] clusterOf = new int[0];   // per dataset; -1 = left out (too few significant loci)
        public int[] members = new int[0];     // dendrogram leaf i -> dataset index
        public List<Integer> excluded = new ArrayList<>();
        public int nClusters = 0;
        public List<Merge> merges = new ArrayList<>();
        public List<CategoryAnalysis> analyses = new ArrayList<>();
    }

    // ── build ────────────────────────────────────────────────────────────────────────────────

    public static Result build(MultiLocusResult mlr) { return build(mlr, DEFAULT_PERMUTATIONS); }

    public static Result build(MultiLocusResult mlr, int nPermutations) {
        Result r = new Result();
        r.name = mlr.name; r.refPanelLabel = mlr.refPanelLabel; r.createdAt = mlr.createdAt;
        int n = mlr.datasets.size();

        for (MultiLocusResult.DatasetInfo d : mlr.datasets) {
            Dataset ds = new Dataset();
            ds.id = d.id; ds.name = d.name;
            ds.disease = GeneConstellationBuilder.diseaseGroupOf(d.id, d.diseaseName);
            ds.ancestry = normalizeAncestry(d.ancestry);
            ds.sampleN = d.sampleN;
            r.datasets.add(ds);
        }

        // categories: Disease, Ancestry, sample-size tertile (power control), then user categories
        r.categoryNames.add("Disease");
        r.categoryNames.add("Ancestry");
        for (Dataset ds : r.datasets) { ds.categories.put("Disease", ds.disease); ds.categories.put("Ancestry", ds.ancestry); }
        String[] nTier = sampleSizeTertiles(r.datasets);
        if (nTier != null) {
            r.categoryNames.add("Sample size");
            for (int i = 0; i < n; i++) r.datasets.get(i).categories.put("Sample size", nTier[i]);
        }
        Set<String> custom = new TreeSet<>();
        for (MultiLocusResult.DatasetInfo d : mlr.datasets) custom.addAll(d.categories.keySet());
        for (String c : custom) {
            r.categoryNames.add(c);
            for (int i = 0; i < n; i++) {
                String v = mlr.datasets.get(i).categories.get(c);
                r.datasets.get(i).categories.put(c, v == null ? "" : v.trim());
            }
        }

        // loci + cells
        for (MultiLocusResult.LocusRow row : mlr.loci) {
            Locus L = new Locus();
            L.index = row.index; L.chr = row.chr; L.start = row.start; L.end = row.end;
            L.gene = row.nearestGene == null ? "" : row.nearestGene;
            MultiLocusResult.DatasetLocusStat best = row.cells.get(row.overallBestDatasetId);
            L.refPos = best != null && best.bestPos > 0 ? best.bestPos : (row.start + row.end) / 2;
            L.cells = new Cell[n];
            for (int i = 0; i < n; i++) {
                Cell c = new Cell();
                MultiLocusResult.DatasetLocusStat s = row.cells.get(mlr.datasets.get(i).id);
                if (s != null) {
                    c.p = s.bestP; c.tier = tierOf(s.bestP); c.pos = s.bestPos; c.snp = s.bestSnpId == null ? "" : s.bestSnpId;
                    c.offsetBp = s.bestPos > 0 ? s.bestPos - L.refPos : 0;
                }
                if (c.tier == 3) L.nDatasetsGws++;
                L.cells[i] = c;
            }
            r.loci.add(L);
        }

        if (n == 0) return r;

        // clustering on informative datasets only
        double[][] dist = tierDistances(r.loci, n);
        boolean[] use = new boolean[n];
        List<Integer> mem = new ArrayList<>();
        for (int i = 0; i < n; i++) {
            int gws = 0; for (Locus L : r.loci) if (L.cells[i].tier == 3) gws++;
            use[i] = gws >= MIN_GWS_LOCI;
            if (use[i]) mem.add(i); else r.excluded.add(i);
        }
        if (mem.size() < 3) { // too few informative datasets to cluster: keep everyone rather than nothing
            mem.clear(); r.excluded.clear();
            for (int i = 0; i < n; i++) { use[i] = true; mem.add(i); }
        }
        int m = mem.size();
        r.members = new int[m];
        for (int a = 0; a < m; a++) r.members[a] = mem.get(a);
        double[][] sub = new double[m][m];
        for (int a = 0; a < m; a++) for (int b = 0; b < m; b++) sub[a][b] = dist[r.members[a]][r.members[b]];
        r.merges = upgma(sub);
        r.nClusters = chooseClusterCount(r.merges, sub, m);
        int[] subOrder = leafOrder(r.merges, m), subCluster = cutTree(r.merges, m, r.nClusters);
        r.leafOrder = new int[n];
        r.clusterOf = new int[n];
        Arrays.fill(r.clusterOf, -1);
        int w = 0;
        for (int leaf : subOrder) r.leafOrder[w++] = r.members[leaf];
        for (int i : r.excluded) r.leafOrder[w++] = i;
        for (int a = 0; a < m; a++) r.clusterOf[r.members[a]] = subCluster[a];

        // category analyses (same informative datasets)
        for (String cat : r.categoryNames) {
            String[] labels = new String[n];
            for (int i = 0; i < n; i++) labels[i] = r.datasets.get(i).categories.getOrDefault(cat, "");
            r.analyses.add(analyzeCategory(cat, labels, use, dist, r.clusterOf, r.loci, nPermutations));
        }
        return r;
    }

    /** "EUR (Spanish)" -> "EUR", "Multi (EUR+AFR...)" -> "MULTI", "" -> "". */
    static String normalizeAncestry(String a) {
        if (a == null) return "";
        String t = a.trim();
        if (t.isEmpty()) return "";
        int cut = t.length();
        for (int k = 0; k < t.length(); k++) { char ch = t.charAt(k); if (ch == ' ' || ch == '(' || ch == ',' || ch == ';') { cut = k; break; } }
        return t.substring(0, cut).toUpperCase(Locale.ROOT);
    }

    /** Low/Mid/High tertiles of sample size, or null when fewer than 6 datasets report N. */
    static String[] sampleSizeTertiles(List<Dataset> ds) {
        List<Integer> ns = new ArrayList<>();
        for (Dataset d : ds) if (d.sampleN > 0) ns.add(d.sampleN);
        if (ns.size() < 6) return null;
        Collections.sort(ns);
        int t1 = ns.get(ns.size() / 3), t2 = ns.get((2 * ns.size()) / 3);
        String[] out = new String[ds.size()];
        for (int i = 0; i < ds.size(); i++) {
            int v = ds.get(i).sampleN;
            out[i] = v <= 0 ? "" : v < t1 ? "Low N" : v < t2 ? "Mid N" : "High N";
        }
        return out;
    }

    // ── distances + clustering ───────────────────────────────────────────────────────────────

    /** Cosine distance between datasets' tier vectors; an all-zero vector is at distance 1 from every other dataset. */
    static double[][] tierDistances(List<Locus> loci, int n) {
        double[][] v = new double[n][loci.size()];
        for (int k = 0; k < loci.size(); k++) for (int i = 0; i < n; i++) v[i][k] = loci.get(k).cells[i].tier;
        double[] norm = new double[n];
        for (int i = 0; i < n; i++) { double s = 0; for (double x : v[i]) s += x * x; norm[i] = Math.sqrt(s); }
        double[][] d = new double[n][n];
        for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) {
            double dd;
            if (norm[i] == 0 || norm[j] == 0) dd = 1.0;
            else { double dot = 0; for (int k = 0; k < v[i].length; k++) dot += v[i][k] * v[j][k]; dd = 1.0 - dot / (norm[i] * norm[j]); }
            dd = Math.max(0, Math.min(1, dd));
            d[i][j] = d[j][i] = dd;
        }
        return d;
    }

    /** Average-linkage (UPGMA) agglomerative clustering via the Lance-Williams update. Leaves are
     *  0..n-1; the k-th merge creates node n+k. Ties break on the lowest cluster ids, so it is deterministic. */
    static List<Merge> upgma(double[][] d0) {
        int n = d0.length;
        List<Merge> merges = new ArrayList<>();
        if (n < 2) return merges;
        int total = 2 * n - 1;
        double[][] d = new double[total][total];
        for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) d[i][j] = d0[i][j];
        int[] size = new int[total];
        boolean[] active = new boolean[total];
        for (int i = 0; i < n; i++) { size[i] = 1; active[i] = true; }
        int next = n;
        while (next < total) {
            int ba = -1, bb = -1; double best = Double.POSITIVE_INFINITY;
            for (int a = 0; a < next; a++) if (active[a]) for (int b = a + 1; b < next; b++) if (active[b]) {
                if (d[a][b] < best - 1e-15) { best = d[a][b]; ba = a; bb = b; }
            }
            Merge m = new Merge(); m.left = ba; m.right = bb; m.height = best; m.size = size[ba] + size[bb];
            merges.add(m);
            int c = next++;
            size[c] = m.size; active[c] = true; active[ba] = false; active[bb] = false;
            for (int k = 0; k < c; k++) if (active[k] && k != c) {
                double v = (size[ba] * d[ba][k] + size[bb] * d[bb][k]) / (double) (size[ba] + size[bb]);
                d[c][k] = d[k][c] = v;
            }
        }
        return merges;
    }

    static int[] leafOrder(List<Merge> merges, int n) {
        if (n == 0) return new int[0];
        if (merges.isEmpty()) { int[] o = new int[n]; for (int i = 0; i < n; i++) o[i] = i; return o; }
        List<Integer> out = new ArrayList<>();
        Deque<Integer> stack = new ArrayDeque<>();
        stack.push(n + merges.size() - 1);
        while (!stack.isEmpty()) {
            int node = stack.pop();
            if (node < n) { out.add(node); continue; }
            Merge m = merges.get(node - n);
            stack.push(m.right); stack.push(m.left);
        }
        int[] o = new int[out.size()];
        for (int i = 0; i < o.length; i++) o[i] = out.get(i);
        return o;
    }

    /**
     * Number of clusters to cut the tree into: the K in [2, max(2, n/2)] with the highest mean
     * silhouette width (Rousseeuw 1987) on the same distance matrix; ties go to the smaller K.
     * Merging two genuinely different groups lowers the silhouette, while an outlier dataset that
     * shares almost no loci with anything (e.g. an X-chromosome-only GWAS, or one with no genome-wide
     * hits) can split off as a singleton at no cost (a singleton's silhouette is 0), so outliers do
     * not force real groups together, and a corpus with three real groups is not split into four.
     */
    static int chooseClusterCount(List<Merge> merges, double[][] d, int n) {
        if (n < 3) return n;
        int maxK = Math.max(2, n / 2);
        int bestK = 2; double best = Double.NEGATIVE_INFINITY;
        for (int k = 2; k <= maxK && k < n; k++) {
            double sw = meanSilhouette(d, cutTree(merges, n, k), k);
            if (sw > best + 1e-9) { best = sw; bestK = k; }
        }
        return bestK;
    }

    /** Mean silhouette width; singleton clusters contribute 0, as in the standard definition. */
    static double meanSilhouette(double[][] d, int[] cl, int k) {
        int n = cl.length;
        int[] size = new int[k];
        for (int c : cl) size[c]++;
        double total = 0;
        for (int i = 0; i < n; i++) {
            if (size[cl[i]] < 2) continue;
            double[] sum = new double[k];
            for (int j = 0; j < n; j++) if (j != i) sum[cl[j]] += d[i][j];
            double a = sum[cl[i]] / (size[cl[i]] - 1);
            double b = Double.POSITIVE_INFINITY;
            for (int c = 0; c < k; c++) if (c != cl[i] && size[c] > 0) b = Math.min(b, sum[c] / size[c]);
            double m = Math.max(a, b);
            total += m > 0 ? (b - a) / m : 0;
        }
        return total / n;
    }

    /** Cluster ids 0..k-1 after applying the first n-k merges; ids are numbered in leaf order. */
    static int[] cutTree(List<Merge> merges, int n, int k) {
        int[] parent = new int[2 * n];
        for (int i = 0; i < parent.length; i++) parent[i] = i;
        int apply = Math.max(0, Math.min(merges.size(), n - k));
        for (int t = 0; t < apply; t++) {
            Merge m = merges.get(t);
            int node = n + t;
            parent[find(parent, m.left)] = node;
            parent[find(parent, m.right)] = node;
        }
        int[] order = leafOrder(merges, n);
        Map<Integer, Integer> ids = new HashMap<>();
        int[] out = new int[n];
        for (int leaf : order) {
            int root = find(parent, leaf);
            Integer id = ids.get(root);
            if (id == null) { id = ids.size(); ids.put(root, id); }
            out[leaf] = id;
        }
        return out;
    }

    private static int find(int[] p, int x) { while (p[x] != x) { p[x] = p[p[x]]; x = p[x]; } return x; }

    // ── statistics ───────────────────────────────────────────────────────────────────────────

    public static class Permanova { public double ssTotal, ssWithin, ssBetween, r2, pseudoF, p; public int nPerm; }

    /** One-way PERMANOVA (Anderson 2001) on a distance matrix; groups[i] in 0..a-1. p = (#{F* >= F} + 1)/(nPerm + 1). */
    static Permanova permanova(double[][] d, int[] groups, int nPerm, long seed) {
        int N = groups.length;
        int a = 0; for (int g : groups) a = Math.max(a, g + 1);
        Permanova res = new Permanova();
        double ssT = 0;
        for (int i = 0; i < N; i++) for (int j = i + 1; j < N; j++) ssT += d[i][j] * d[i][j];
        ssT /= N;
        res.ssTotal = ssT;
        res.ssWithin = ssWithin(d, groups, a);
        res.ssBetween = ssT - res.ssWithin;
        res.r2 = ssT > 0 ? res.ssBetween / ssT : Double.NaN;
        res.pseudoF = fStat(ssT, res.ssWithin, a, N);
        res.nPerm = nPerm;
        if (nPerm <= 0 || Double.isNaN(res.pseudoF)) { res.p = Double.NaN; return res; }
        Random rnd = new Random(seed);
        int[] perm = groups.clone();
        int ge = 0;
        for (int t = 0; t < nPerm; t++) {
            for (int i = N - 1; i > 0; i--) { int j = rnd.nextInt(i + 1); int tmp = perm[i]; perm[i] = perm[j]; perm[j] = tmp; }
            double f = fStat(ssT, ssWithin(d, perm, a), a, N);
            if (f >= res.pseudoF - 1e-12) ge++;
        }
        res.p = (ge + 1.0) / (nPerm + 1.0);
        return res;
    }

    private static double ssWithin(double[][] d, int[] groups, int a) {
        double[] ss = new double[a]; int[] cnt = new int[a];
        for (int g : groups) cnt[g]++;
        for (int i = 0; i < groups.length; i++) for (int j = i + 1; j < groups.length; j++)
            if (groups[i] == groups[j]) ss[groups[i]] += d[i][j] * d[i][j];
        double tot = 0;
        for (int g = 0; g < a; g++) if (cnt[g] > 0) tot += ss[g] / cnt[g];
        return tot;
    }

    private static double fStat(double ssT, double ssW, int a, int N) {
        if (a < 2 || N - a < 1) return Double.NaN;
        double ssB = ssT - ssW;
        if (ssW <= 1e-15) return ssB > 1e-15 ? Double.POSITIVE_INFINITY : Double.NaN;
        return (ssB / (a - 1)) / (ssW / (N - a));
    }

    /** Adjusted Rand index (Hubert & Arabie 1985) between two labelings of the same items. */
    static double adjustedRandIndex(int[] x, int[] y) {
        int n = x.length;
        if (n < 2) return Double.NaN;
        Map<Long, Integer> nij = new HashMap<>();
        Map<Integer, Integer> ai = new HashMap<>(), bj = new HashMap<>();
        for (int i = 0; i < n; i++) {
            nij.merge(((long) x[i] << 32) | (y[i] & 0xffffffffL), 1, Integer::sum);
            ai.merge(x[i], 1, Integer::sum); bj.merge(y[i], 1, Integer::sum);
        }
        double sumNij = 0, sumA = 0, sumB = 0;
        for (int v : nij.values()) sumNij += c2(v);
        for (int v : ai.values()) sumA += c2(v);
        for (int v : bj.values()) sumB += c2(v);
        double expected = sumA * sumB / c2(n);
        double max = 0.5 * (sumA + sumB);
        if (Math.abs(max - expected) < 1e-15) return 1.0; // both partitions trivial and identical in structure
        return (sumNij - expected) / (max - expected);
    }

    private static double c2(int v) { return v * (v - 1) / 2.0; }

    // ── category analysis ────────────────────────────────────────────────────────────────────

    static CategoryAnalysis analyzeCategory(String cat, String[] labels, boolean[] use, double[][] dist, int[] clusterOf,
                                            List<Locus> loci, int nPerm) {
        CategoryAnalysis ca = new CategoryAnalysis();
        ca.category = cat;
        List<Integer> idx = new ArrayList<>();
        int unlabeled = 0;
        for (int i = 0; i < labels.length; i++) {
            if (!use[i]) continue;
            if (labels[i] != null && !labels[i].isEmpty()) idx.add(i); else unlabeled++;
        }
        Map<String, Integer> levelId = new LinkedHashMap<>();
        for (int i : idx) levelId.putIfAbsent(labels[i], levelId.size());
        ca.nDatasets = idx.size();
        ca.nLevels = levelId.size();

        int[] g = new int[idx.size()];
        int[] cl = new int[idx.size()];
        double[][] sub = new double[idx.size()][idx.size()];
        int[] levelSize = new int[levelId.size()];
        for (int a = 0; a < idx.size(); a++) {
            g[a] = levelId.get(labels[idx.get(a)]);
            cl[a] = clusterOf[idx.get(a)];
            levelSize[g[a]]++;
            for (int b = 0; b < idx.size(); b++) sub[a][b] = dist[idx.get(a)][idx.get(b)];
        }
        boolean anyReplicated = false; for (int s : levelSize) if (s >= 2) anyReplicated = true;

        if (ca.nLevels < 2) ca.note = "Only one " + cat.toLowerCase(Locale.ROOT) + " level among the selected datasets, so it cannot separate them.";
        else if (!anyReplicated) ca.note = "Every " + cat.toLowerCase(Locale.ROOT) + " level has a single dataset, so within-level consistency cannot be measured.";
        else if (ca.nDatasets - ca.nLevels < 1) ca.note = "Too few datasets for the number of levels.";
        else {
            ca.testable = true;
            int replicated = 0; for (int s2 : levelSize) if (s2 >= 2) replicated++;
            if (replicated < 2) ca.note = "Only one " + cat.toLowerCase(Locale.ROOT)
                + " group has two or more datasets, so this test has little power; treat it as uninformative, not as evidence of no effect.";
            Permanova pm = permanova(sub, g, nPerm, PERMUTATION_SEED + cat.hashCode());
            ca.r2 = pm.r2; ca.pseudoF = pm.pseudoF; ca.permP = pm.p; ca.nPermutations = nPerm;
            ca.ari = adjustedRandIndex(g, cl);
            if (unlabeled > 0) ca.note = (ca.note.isEmpty() ? "" : ca.note + " ")
                + unlabeled + " dataset(s) without a " + cat.toLowerCase(Locale.ROOT) + " value were left out.";
        }

        // per-level cohesion / separation / consistently significant loci
        for (Map.Entry<String, Integer> e : levelId.entrySet()) {
            Level L = new Level();
            L.level = e.getKey();
            int lv = e.getValue();
            L.nDatasets = levelSize[lv];
            double win = 0, out = 0; int nw = 0, no = 0;
            for (int a = 0; a < g.length; a++) for (int b = a + 1; b < g.length; b++) {
                boolean ia = g[a] == lv, ib = g[b] == lv;
                if (ia && ib) { win += sub[a][b]; nw++; } else if (ia || ib) { out += sub[a][b]; no++; }
            }
            if (nw > 0) L.cohesion = 1.0 - win / nw;
            if (no > 0) L.separation = out / no;
            if (L.nDatasets >= 2) {
                List<ConsistentLocus> found = new ArrayList<>();
                for (Locus loc : loci) {
                    int inSig = 0, outSig = 0, nIn = 0, nOut = 0;
                    for (int a = 0; a < idx.size(); a++) {
                        boolean sig = loc.cells[idx.get(a)].tier == 3;
                        if (g[a] == lv) { nIn++; if (sig) inSig++; } else { nOut++; if (sig) outSig++; }
                    }
                    double fin = nIn > 0 ? inSig / (double) nIn : 0, fout = nOut > 0 ? outSig / (double) nOut : 0;
                    if (inSig >= 2 && fin >= CONSISTENT_FRACTION && fin > fout) {
                        double fp = EnrichmentAnalyzer.fisherExactTwoSided(inSig, nIn - inSig, outSig, nOut - outSig);
                        if (fp >= CONSISTENT_FISHER_P) continue;
                        ConsistentLocus c = new ConsistentLocus();
                        c.locusIndex = loc.index; c.gene = loc.gene; c.chr = loc.chr; c.pos = loc.refPos;
                        c.nIn = nIn; c.nInSig = inSig; c.nOut = nOut; c.nOutSig = outSig;
                        c.fisherP = fp;
                        found.add(c);
                    }
                }
                found.sort(Comparator.comparingDouble((ConsistentLocus c) -> c.fisherP)
                        .thenComparing(c -> -(double) c.nInSig / Math.max(1, c.nIn)));
                L.nConsistentLoci = found.size();
                L.consistentLoci = new ArrayList<>(found.subList(0, Math.min(MAX_CONSISTENT_LOCI_PER_LEVEL, found.size())));
            }
            ca.levels.add(L);
        }
        ca.levels.sort(Comparator.comparingInt((Level l) -> -l.nConsistentLoci).thenComparing(l -> l.level));
        return ca;
    }

    // ── locus detail (click on a serpent) ────────────────────────────────────────────────────

    /** Nearest gene to a position: returns {name, distanceBp} or null when nothing is annotated nearby. */
    public interface GeneLookup { Object[] nearest(String chr, long pos); }

    /**
     * Per-dataset lead SNPs for one locus, as shown when the user clicks a serpent: lead SNP id,
     * position, alleles, OR with 95% CI, p-value, tier, offset from the strongest lead SNP, and the
     * nearest gene to that lead SNP. Rows are sorted by p-value. The OR is the reported OR when the
     * dataset has one, otherwise exp(beta); the CI is exp(ln OR +/- 1.96 SE) when an SE exists.
     * The header's gene list is the nearest genes of all lead SNPs, deduplicated, in the order of the
     * strongest dataset first, with how many datasets' lead SNPs point to each.
     */
    public static String locusDetailJson(MultiLocusResult mlr, int locusIndex, GeneLookup genes) {
        MultiLocusResult.LocusRow row = null;
        for (MultiLocusResult.LocusRow r : mlr.loci) if (r.index == locusIndex) { row = r; break; }
        if (row == null) return null;
        MultiLocusResult.DatasetLocusStat bestCell = row.cells.get(row.overallBestDatasetId);
        long refPos = bestCell != null && bestCell.bestPos > 0 ? bestCell.bestPos : (row.start + row.end) / 2;

        class R { MultiLocusResult.DatasetInfo d; MultiLocusResult.DatasetLocusStat s; double or, lo, hi; String gene = ""; long dist = -1; }
        List<R> rows = new ArrayList<>();
        for (MultiLocusResult.DatasetInfo d : mlr.datasets) {
            MultiLocusResult.DatasetLocusStat s = row.cells.get(d.id);
            if (s == null || s.bestPos <= 0 || Double.isNaN(s.bestP)) continue;
            R x = new R(); x.d = d; x.s = s;
            double ln;
            if (!Double.isNaN(s.or) && s.or > 0) ln = Math.log(s.or);
            else ln = s.beta;
            x.or = Double.isNaN(ln) ? Double.NaN : Math.exp(ln);
            x.lo = x.hi = Double.NaN;
            if (!Double.isNaN(ln) && !Double.isNaN(s.se) && s.se > 0) { x.lo = Math.exp(ln - 1.959964 * s.se); x.hi = Math.exp(ln + 1.959964 * s.se); }
            String chr = s.bestChr == null || s.bestChr.isEmpty() ? row.chr : s.bestChr;
            Object[] g = genes == null ? null : genes.nearest(chr, s.bestPos);
            if (g != null) { x.gene = String.valueOf(g[0]); x.dist = ((Number) g[1]).longValue(); }
            rows.add(x);
        }
        rows.sort(Comparator.comparingDouble(x -> x.s.bestP));

        Map<String, Integer> geneCount = new LinkedHashMap<>();
        Map<String, Long> geneMinDist = new HashMap<>();
        for (R x : rows) if (!x.gene.isEmpty()) {
            geneCount.merge(x.gene, 1, Integer::sum);
            geneMinDist.merge(x.gene, x.dist, Math::min);
        }

        StringBuilder j = new StringBuilder(4096);
        j.append("{\"locus\":{\"index\":").append(row.index).append(",");
        kv(j, "chr", row.chr).append(",\"start\":").append(row.start).append(",\"end\":").append(row.end)
            .append(",\"ref_pos\":").append(refPos).append(",");
        kv(j, "locus_nearest_genes", row.nearestGene == null ? "" : row.nearestGene).append(",");
        kv(j, "best_dataset_id", row.overallBestDatasetId).append(",");
        int gws = 0; for (R x : rows) if (x.s.bestP < 5e-8) gws++;
        j.append("\"n_datasets\":").append(rows.size()).append(",\"n_gws\":").append(gws).append("},\"genes\":[");
        int k = 0;
        for (Map.Entry<String, Integer> e : geneCount.entrySet()) {
            if (k++ > 0) j.append(",");
            j.append("{"); kv(j, "gene", e.getKey()).append(",\"n_lead_snps\":").append(e.getValue())
             .append(",\"min_distance\":").append(geneMinDist.get(e.getKey())).append("}");
        }
        j.append("],\"rows\":[");
        for (int i = 0; i < rows.size(); i++) {
            R x = rows.get(i);
            if (i > 0) j.append(",");
            j.append("{"); kv(j, "dataset_id", x.d.id).append(","); kv(j, "dataset", x.d.name).append(",");
            kv(j, "disease", GeneConstellationBuilder.diseaseGroupOf(x.d.id, x.d.diseaseName)).append(",");
            kv(j, "ancestry", normalizeAncestry(x.d.ancestry)).append(",");
            kv(j, "snp", x.s.bestSnpId).append(",");
            kv(j, "chr", x.s.bestChr == null || x.s.bestChr.isEmpty() ? row.chr : x.s.bestChr).append(",");
            j.append("\"pos\":").append(x.s.bestPos).append(",\"offset\":").append(x.s.bestPos - refPos).append(",");
            kv(j, "ea", x.s.ea).append(","); kv(j, "nea", x.s.nea).append(",");
            j.append("\"p\":").append(num(x.s.bestP)).append(",\"tier\":").append(tierOf(x.s.bestP))
             .append(",\"or\":").append(num(x.or)).append(",\"or_lo\":").append(num(x.lo)).append(",\"or_hi\":").append(num(x.hi))
             .append(",\"or_from_beta\":").append(Double.isNaN(x.s.or) || x.s.or <= 0).append(",");
            kv(j, "nearest_gene", x.gene).append(",\"gene_distance\":").append(x.dist).append("}");
        }
        j.append("]}");
        return j.toString();
    }

    // ── JSON ─────────────────────────────────────────────────────────────────────────────────

    public static String toJson(Result r) {
        StringBuilder j = new StringBuilder(1 << 16);
        j.append("{");
        kv(j, "name", r.name).append(",");
        kv(j, "ref_panel_label", r.refPanelLabel).append(",");
        kv(j, "created_at", r.createdAt).append(",");
        j.append("\"categories\":[");
        for (int i = 0; i < r.categoryNames.size(); i++) { if (i > 0) j.append(","); str(j, r.categoryNames.get(i)); }
        j.append("],\"datasets\":[");
        for (int i = 0; i < r.datasets.size(); i++) {
            Dataset d = r.datasets.get(i);
            if (i > 0) j.append(",");
            j.append("{"); kv(j, "id", d.id).append(","); kv(j, "name", d.name).append(",");
            kv(j, "disease", d.disease).append(","); kv(j, "ancestry", d.ancestry).append(",");
            j.append("\"sample_n\":").append(d.sampleN).append(",\"categories\":{");
            int k = 0;
            for (Map.Entry<String, String> e : d.categories.entrySet()) { if (k++ > 0) j.append(","); kv(j, e.getKey(), e.getValue()); }
            j.append("}}");
        }
        j.append("],\"loci\":[");
        for (int i = 0; i < r.loci.size(); i++) {
            Locus L = r.loci.get(i);
            if (i > 0) j.append(",");
            j.append("{\"index\":").append(L.index).append(",");
            kv(j, "chr", L.chr).append(",");
            j.append("\"start\":").append(L.start).append(",\"end\":").append(L.end).append(",\"ref_pos\":").append(L.refPos).append(",");
            kv(j, "gene", L.gene).append(",");
            j.append("\"n_gws\":").append(L.nDatasetsGws).append(",");
            // compact per-dataset arrays, in dataset order
            j.append("\"tier\":["); for (int k = 0; k < L.cells.length; k++) { if (k > 0) j.append(","); j.append(L.cells[k].tier); }
            j.append("],\"p\":["); for (int k = 0; k < L.cells.length; k++) { if (k > 0) j.append(","); j.append(num(L.cells[k].p)); }
            j.append("],\"off\":["); for (int k = 0; k < L.cells.length; k++) { if (k > 0) j.append(","); j.append(L.cells[k].offsetBp); }
            j.append("],\"snp\":["); for (int k = 0; k < L.cells.length; k++) { if (k > 0) j.append(","); str(j, L.cells[k].snp); }
            j.append("]}");
        }
        j.append("],\"clustering\":{\"n_clusters\":").append(r.nClusters).append(",\"min_gws_loci\":").append(MIN_GWS_LOCI)
         .append(",\"leaf_order\":");
        ints(j, r.leafOrder); j.append(",\"cluster_of\":"); ints(j, r.clusterOf);
        j.append(",\"members\":"); ints(j, r.members);
        int[] ex = new int[r.excluded.size()]; for (int i = 0; i < ex.length; i++) ex[i] = r.excluded.get(i);
        j.append(",\"excluded\":"); ints(j, ex);
        j.append(",\"merges\":[");
        for (int i = 0; i < r.merges.size(); i++) {
            Merge m = r.merges.get(i);
            if (i > 0) j.append(",");
            j.append("[").append(m.left).append(",").append(m.right).append(",").append(num(m.height)).append(",").append(m.size).append("]");
        }
        j.append("]},\"analyses\":[");
        for (int i = 0; i < r.analyses.size(); i++) {
            CategoryAnalysis a = r.analyses.get(i);
            if (i > 0) j.append(",");
            j.append("{"); kv(j, "category", a.category).append(",");
            j.append("\"testable\":").append(a.testable).append(","); kv(j, "note", a.note).append(",");
            j.append("\"n_datasets\":").append(a.nDatasets).append(",\"n_levels\":").append(a.nLevels).append(",");
            j.append("\"r2\":").append(num(a.r2)).append(",\"pseudo_f\":").append(num(a.pseudoF)).append(",");
            j.append("\"perm_p\":").append(num(a.permP)).append(",\"n_permutations\":").append(a.nPermutations).append(",");
            j.append("\"ari\":").append(num(a.ari)).append(",\"levels\":[");
            for (int k = 0; k < a.levels.size(); k++) {
                Level L = a.levels.get(k);
                if (k > 0) j.append(",");
                j.append("{"); kv(j, "level", L.level).append(",");
                j.append("\"n_datasets\":").append(L.nDatasets).append(",\"cohesion\":").append(num(L.cohesion))
                 .append(",\"separation\":").append(num(L.separation)).append(",\"n_consistent_loci\":").append(L.nConsistentLoci)
                 .append(",\"consistent_loci\":[");
                for (int q = 0; q < L.consistentLoci.size(); q++) {
                    ConsistentLocus c = L.consistentLoci.get(q);
                    if (q > 0) j.append(",");
                    j.append("{\"locus_index\":").append(c.locusIndex).append(","); kv(j, "gene", c.gene).append(",");
                    kv(j, "chr", c.chr).append(",\"pos\":").append(c.pos).append(",\"n_in\":").append(c.nIn)
                     .append(",\"n_in_sig\":").append(c.nInSig).append(",\"n_out\":").append(c.nOut)
                     .append(",\"n_out_sig\":").append(c.nOutSig).append(",\"fisher_p\":").append(num(c.fisherP)).append("}");
                }
                j.append("]}");
            }
            j.append("]}");
        }
        j.append("]}");
        return j.toString();
    }

    private static StringBuilder kv(StringBuilder j, String k, String v) { str(j, k); j.append(":"); str(j, v); return j; }

    private static void str(StringBuilder j, String s) {
        j.append('"');
        if (s != null) for (int i = 0; i < s.length(); i++) {
            char c = s.charAt(i);
            switch (c) {
                case '"': j.append("\\\""); break;
                case '\\': j.append("\\\\"); break;
                case '\n': j.append("\\n"); break;
                case '\r': j.append("\\r"); break;
                case '\t': j.append("\\t"); break;
                default: if (c < 0x20) j.append(String.format("\\u%04x", (int) c)); else j.append(c);
            }
        }
        j.append('"');
    }

    private static void ints(StringBuilder j, int[] a) {
        j.append("["); for (int i = 0; i < a.length; i++) { if (i > 0) j.append(","); j.append(a[i]); } j.append("]");
    }

    private static String num(double d) {
        if (Double.isNaN(d)) return "null";
        if (Double.isInfinite(d)) return d > 0 ? "1e308" : "-1e308";
        return String.valueOf(d);
    }
}
