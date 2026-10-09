import java.io.*;
import java.nio.charset.StandardCharsets;
import java.nio.file.*;
import java.util.*;
import java.util.concurrent.ConcurrentHashMap;
import java.util.zip.GZIPInputStream;

/**
 * scE2G enhancer-gene links, bundled with LYNXgwas (regulatory_links/sce2g/, built by
 * scripts/regulatory/build_sce2g_pack.py from the public predictions of Sheth, Qiu et al., Engreitz lab).
 *
 * Per cell type: the element-gene pairs scE2G_multiome calls (E2G.Score.qnorm >= 0.177), in GRCh37 (lifted) and the
 * original GRCh38. Loaded per cell type on first use into compact arrays.
 *
 * Two uses: a locus view (the elements in a window, merged across cell types, the element nearest a SNP, its target
 * genes) and project statistics (where the lead SNPs sit relative to the elements, and per-cell-type enrichment against
 * random SNPs matched on the distance to the nearest gene TSS).
 *
 * Positions follow the scE2G files: an element is BED-like (start, end]; a SNP at pos is inside when start < pos <= end.
 */
public class Sce2gLinks {

    public static final File DIR = new File("regulatory_links/sce2g");
    static final String[] CLASSES = {"promoter", "genic", "intergenic"};

    public static class CellType {
        public String id, label, group, lineage, file;
        public int pairs;
    }

    /** One cell type, one chromosome: link rows, plus the distinct elements sorted by start for each build. */
    static class Chrom {
        int n;
        int[] s37, e37, s38, e38, tss37, tss38, gene;
        byte[] cls;
        float[] e2g, abc, kendall;
        boolean[] self;
        int[][] el37, el38;               // [start, end, class, firstRow] per distinct element, sorted by start
    }

    private static final List<String> GENES = new ArrayList<>();
    private static final Map<String, Integer> GENE_IDX = new HashMap<>();
    private static volatile List<CellType> cellTypes;
    private static final Map<String, Map<String, Chrom>> LOADED = new ConcurrentHashMap<>();
    private static final Map<String, Object> LOCKS = new ConcurrentHashMap<>();   // one per cell type: files load in parallel
    private static String citation = "", threshold = "";

    // ── Manifest and selection ──────────────────────────────────────────

    public static synchronized List<CellType> cellTypes() {
        if (cellTypes != null) return cellTypes;
        List<CellType> out = new ArrayList<>();
        File mf = new File(DIR, "manifest.json");
        if (mf.isFile()) {
            try {
                Map<String, Object> m = MiniJson.asObject(MiniJson.parse(new String(Files.readAllBytes(mf.toPath()), StandardCharsets.UTF_8)));
                citation = MiniJson.getStr(m, "citation", "");
                threshold = MiniJson.getStr(m, "threshold", "");
                for (Object o : MiniJson.asArray(m.get("cell_types"))) {
                    Map<String, Object> e = MiniJson.asObject(o);
                    CellType c = new CellType();
                    c.id = MiniJson.getStr(e, "id", "");
                    c.label = MiniJson.getStr(e, "label", c.id);
                    c.group = MiniJson.getStr(e, "group", "");
                    c.lineage = MiniJson.getStr(e, "lineage", "");
                    c.file = MiniJson.getStr(e, "file", "");
                    Object p = e.get("pairs");
                    c.pairs = p instanceof Number ? ((Number) p).intValue() : 0;
                    if (new File(DIR, c.file).isFile()) out.add(c);
                }
            } catch (Exception e) {
                System.err.println("[scE2G] could not read the manifest: " + e.getMessage());
            }
        }
        cellTypes = out;
        return out;
    }

    public static boolean available() { return !cellTypes().isEmpty(); }

    public static String citation() { cellTypes(); return citation; }

    public static String thresholdText() { cellTypes(); return threshold; }

    static final Set<String> IMMUNE_DISEASES = new HashSet<>(Arrays.asList(
        "ra", "ibd", "sle", "lupus", "ms", "t1d", "celiac", "psoriasis", "crohn", "uc", "colitis", "asthma", "aitd",
        "jia", "pbc", "psc", "as", "ss", "sjogren", "vitiligo", "rheumatoid arthritis", "inflammatory bowel disease",
        "systemic lupus erythematosus", "multiple sclerosis", "type 1 diabetes", "crohn's disease", "ulcerative colitis"));
    static final Set<String> ISLET_DISEASES = new HashSet<>(Arrays.asList("t2d", "type 2 diabetes", "diabetes", "t2dm"));

    /** Named selections offered in the viewer and summary page. */
    public static List<String[]> presets() {
        List<String[]> p = new ArrayList<>();
        p.add(new String[]{"immune", "Immune cells (PBMC + bone marrow)"});
        p.add(new String[]{"pbmc", "PBMC"});
        p.add(new String[]{"bm", "Bone marrow"});
        p.add(new String[]{"islets", "Pancreatic islets"});
        p.add(new String[]{"all", "All cell types"});
        return p;
    }

    /** The default selection for a disease group: islets for type 2 diabetes, immune cells for immune diseases, else all. */
    public static String defaultPreset(String diseaseGroup) {
        String d = diseaseGroup == null ? "" : diseaseGroup.toLowerCase(Locale.ROOT).trim();
        if (ISLET_DISEASES.contains(d)) return "islets";
        if (IMMUNE_DISEASES.contains(d)) return "immune";
        return "all";
    }

    /** Cell types for a selection: a preset name, or a comma-separated list of cell-type ids. */
    public static List<CellType> select(String sel) {
        List<CellType> all = cellTypes(), out = new ArrayList<>();
        String s = sel == null || sel.trim().isEmpty() ? "all" : sel.trim();
        for (CellType c : all) {
            boolean keep;
            switch (s) {
                case "all": keep = true; break;
                case "immune": keep = c.group.equals("PBMC") || c.group.equals("Bone marrow"); break;
                case "pbmc": keep = c.group.equals("PBMC"); break;
                case "bm": keep = c.group.equals("Bone marrow"); break;
                case "islets": keep = c.group.equals("Pancreatic islets"); break;
                default: keep = Arrays.asList(s.split(",")).contains(c.id);
            }
            if (keep) out.add(c);
        }
        return out;
    }

    // ── Loading ─────────────────────────────────────────────────────────

    static String chrKey(String chr) {
        String c = chr.trim();
        return c.startsWith("chr") ? c : "chr" + c;
    }

    static Map<String, Chrom> load(CellType ct) {
        Map<String, Chrom> got = LOADED.get(ct.id);
        if (got != null) return got;
        synchronized (LOCKS.computeIfAbsent(ct.id, k -> new Object())) {
            got = LOADED.get(ct.id);
            if (got != null) return got;
            long t0 = System.currentTimeMillis();
            Map<String, List<String[]>> byChr = new LinkedHashMap<>();
            String[] head = null;
            try (BufferedReader br = new BufferedReader(new InputStreamReader(
                    new GZIPInputStream(new FileInputStream(new File(DIR, ct.file)), 1 << 16), StandardCharsets.UTF_8), 1 << 16)) {
                String line;
                while ((line = br.readLine()) != null) {
                    if (line.startsWith("#")) continue;
                    String[] f = line.split("\t", -1);
                    if (head == null) { head = f; continue; }
                    byChr.computeIfAbsent(f[0], k -> new ArrayList<>()).add(f);
                }
            } catch (IOException e) {
                System.err.println("[scE2G] could not read " + ct.file + ": " + e.getMessage());
                return Collections.emptyMap();
            }
            List<String> h = Arrays.asList(head);
            int iS = h.indexOf("start"), iE = h.indexOf("end"), iC = h.indexOf("class"), iG = h.indexOf("gene"),
                iT = h.indexOf("tss"), iSelf = h.indexOf("self_promoter"), iQ = h.indexOf("e2g"), iA = h.indexOf("abc"),
                iK = h.indexOf("kendall"), iS8 = h.indexOf("start38"), iE8 = h.indexOf("end38"), iT8 = h.indexOf("tss38"),
                // GRCh38 stored as offsets from GRCh37 (scripts/regulatory/delta38.py) instead of absolute columns
                iD = h.indexOf("d38"), iDT = h.indexOf("dtss38"), iDL = h.indexOf("dlen38");
            boolean delta = iS8 < 0 && iD >= 0;
            Map<String, Chrom> out = new HashMap<>();
            // genes are numbered within the file, then mapped to the shared table in one step (no lock per row)
            Map<String, Integer> localIdx = new HashMap<>();
            List<String> localGenes = new ArrayList<>();
            for (Map.Entry<String, List<String[]>> en : byChr.entrySet()) {
                List<String[]> rows = en.getValue();
                Chrom c = new Chrom();
                c.n = rows.size();
                c.s37 = new int[c.n]; c.e37 = new int[c.n]; c.s38 = new int[c.n]; c.e38 = new int[c.n];
                c.tss37 = new int[c.n]; c.tss38 = new int[c.n]; c.gene = new int[c.n]; c.cls = new byte[c.n];
                c.e2g = new float[c.n]; c.abc = new float[c.n]; c.kendall = new float[c.n]; c.self = new boolean[c.n];
                for (int r = 0; r < c.n; r++) {
                    String[] f = rows.get(r);
                    c.s37[r] = Integer.parseInt(f[iS]); c.e37[r] = Integer.parseInt(f[iE]);
                    c.tss37[r] = Integer.parseInt(f[iT]);
                    if (delta) {
                        int d = Integer.parseInt(f[iD]);
                        c.s38[r] = c.s37[r] + d;
                        c.e38[r] = c.s38[r] + (c.e37[r] - c.s37[r]) + Integer.parseInt(f[iDL]);
                        c.tss38[r] = c.tss37[r] + d + Integer.parseInt(f[iDT]);
                    } else {
                        c.s38[r] = Integer.parseInt(f[iS8]); c.e38[r] = Integer.parseInt(f[iE8]);
                        c.tss38[r] = Integer.parseInt(f[iT8]);
                    }
                    Integer gi = localIdx.get(f[iG]);
                    if (gi == null) { gi = localGenes.size(); localGenes.add(f[iG]); localIdx.put(f[iG], gi); }
                    c.gene[r] = gi;
                    c.cls[r] = (byte) Math.max(0, Arrays.asList(CLASSES).indexOf(f[iC]));
                    c.e2g[r] = num(f[iQ]); c.abc[r] = num(f[iA]); c.kendall[r] = num(f[iK]);
                    c.self[r] = "1".equals(f[iSelf]);
                }
                c.el37 = elements(c, false);
                c.el38 = elements(c, true);
                out.put(en.getKey(), c);
            }
            int[] map = geneIndices(localGenes);
            for (Chrom c : out.values()) for (int r = 0; r < c.n; r++) c.gene[r] = map[c.gene[r]];
            LOADED.put(ct.id, out);
            System.out.printf("[scE2G] loaded %s (%d ms)%n", ct.label, System.currentTimeMillis() - t0);
            return out;
        }
    }

    private static synchronized int[] geneIndices(List<String> genes) {
        int[] out = new int[genes.size()];
        for (int k = 0; k < out.length; k++) {
            Integer i = GENE_IDX.get(genes.get(k));
            if (i == null) { i = GENES.size(); GENES.add(genes.get(k)); GENE_IDX.put(genes.get(k), i); }
            out[k] = i;
        }
        return out;
    }

    /** Gene name of a shared index (synchronized: another cell type may be adding genes meanwhile). */
    private static synchronized String geneName(int i) { return GENES.get(i); }

    /** Loads the given cell types, 4 files at a time (each file once); returns when all are in memory.
     *  Capped because reading a file briefly takes ~100 MB of text; more at once only adds memory, not much speed. */
    private static final java.util.concurrent.ExecutorService POOL = java.util.concurrent.Executors.newFixedThreadPool(
        Math.max(1, Math.min(4, Runtime.getRuntime().availableProcessors())), r -> {
            Thread t = new Thread(r, "sce2g-load"); t.setDaemon(true); return t; });

    static void loadAll(List<CellType> cts) {
        List<java.util.concurrent.Future<?>> fs = new ArrayList<>();
        for (CellType ct : cts) if (!LOADED.containsKey(ct.id)) fs.add(POOL.submit(() -> load(ct)));
        for (java.util.concurrent.Future<?> f : fs) {
            try { f.get(); } catch (Exception e) { System.err.println("[scE2G] load failed: " + e.getMessage()); }
        }
    }

    private static float num(String s) {
        try { return s.isEmpty() ? Float.NaN : Float.parseFloat(s); } catch (NumberFormatException e) { return Float.NaN; }
    }

    /** Distinct elements of a chromosome, sorted by start: [start, end, class, first row]. */
    private static int[][] elements(Chrom c, boolean b38) {
        Map<Long, int[]> seen = new LinkedHashMap<>();
        for (int r = 0; r < c.n; r++) {
            int s = b38 ? c.s38[r] : c.s37[r], e = b38 ? c.e38[r] : c.e37[r];
            seen.putIfAbsent(((long) s << 32) | (e & 0xffffffffL), new int[]{s, e, c.cls[r], r});
        }
        int[][] el = seen.values().toArray(new int[0][]);
        Arrays.sort(el, Comparator.comparingInt(a -> a[0]));
        return el;
    }

    // ── Distances ───────────────────────────────────────────────────────

    /** Signed distance from pos to the nearest element of this chromosome (0 = inside); the element via out[0]. */
    static long nearest(int[][] el, long pos, int[] out) {
        if (el == null || el.length == 0) { if (out != null) out[0] = -1; return Long.MAX_VALUE; }
        int n = el.length, lo = 0, hi = n;
        while (lo < hi) {                              // first element with start >= pos (n if none)
            int mid = (lo + hi) >>> 1;
            if (el[mid][0] < pos) lo = mid + 1; else hi = mid;
        }
        long best = Long.MAX_VALUE;
        int bi = -1;
        if (lo < n) { best = el[lo][0] - pos + 1; bi = lo; }
        // elements starting before pos: the one containing it, or the closest end (a few back, for long elements)
        for (int k = lo - 1; k >= 0 && k >= lo - 50; k--) {
            long d = pos <= el[k][1] ? 0 : pos - el[k][1];
            if (d < best) { best = d; bi = k; }
            if (best == 0) break;
        }
        if (out != null) out[0] = bi;
        return best;
    }

    // ── Locus view ──────────────────────────────────────────────────────

    /**
     * The elements in [from, to] for the selected cell types, merged across cell types into regions; for the region
     * nearest each query SNP, its target genes. snps: {pos, pip} (pip < 0 for the lead SNP).
     */
    public static Map<String, Object> locus(List<CellType> cts, boolean b38, String chr, long from, long to,
                                            long leadPos, List<double[]> csSnps) {
        loadAll(cts);
        String ck = chrKey(chr);
        List<int[]> calls = new ArrayList<>();          // [start, end, class, ctIndex, row]
        for (int ci = 0; ci < cts.size(); ci++) {
            Chrom c = load(cts.get(ci)).get(ck);
            if (c == null) continue;
            int[][] el = b38 ? c.el38 : c.el37;
            for (int[] e : el) if (e[1] >= from && e[0] <= to) calls.add(new int[]{e[0], e[1], e[2], ci, e[3]});
        }
        calls.sort(Comparator.comparingInt(a -> a[0]));
        // merge overlapping calls of different cell types into regions
        List<List<int[]>> regions = new ArrayList<>();
        int curEnd = Integer.MIN_VALUE;
        for (int[] cl : calls) {
            if (regions.isEmpty() || cl[0] >= curEnd) { regions.add(new ArrayList<>()); curEnd = cl[1]; }
            else curEnd = Math.max(curEnd, cl[1]);
            regions.get(regions.size() - 1).add(cl);
        }
        List<Object> regionOut = new ArrayList<>();
        long[] rs = new long[regions.size()], re = new long[regions.size()];
        for (int k = 0; k < regions.size(); k++) {
            List<int[]> g = regions.get(k);
            int s = Integer.MAX_VALUE, e = Integer.MIN_VALUE;
            int[] clsCount = new int[3];
            Set<Integer> ctSet = new TreeSet<>();
            float best = 0;
            for (int[] cl : g) {
                s = Math.min(s, cl[0]); e = Math.max(e, cl[1]);
                clsCount[cl[2]]++;
                ctSet.add(cl[3]);
                best = Math.max(best, maxScore(cts.get(cl[3]), ck, cl[0], cl[1], b38));
            }
            rs[k] = s; re[k] = e;
            int cls = clsCount[0] > 0 ? 0 : (clsCount[1] >= clsCount[2] ? 1 : 2);   // promoter wins, else majority
            Map<String, Object> r = new LinkedHashMap<>();
            r.put("s", s); r.put("e", e); r.put("cls", CLASSES[cls]); r.put("n", ctSet.size());
            r.put("score", round(best));
            List<Object> ctl = new ArrayList<>();
            for (int ci : ctSet) ctl.add(cts.get(ci).label);
            r.put("cts", ctl);
            regionOut.add(r);
        }
        Map<String, Object> out = new LinkedHashMap<>();
        out.put("regions", regionOut);
        // the region nearest the lead SNP, and its target genes
        int lead = nearestRegion(rs, re, leadPos);
        if (lead >= 0) {
            Map<String, Object> nr = new LinkedHashMap<>(MiniJson.asObject(regionOut.get(lead)));
            nr.put("index", lead);
            nr.put("dist", dist(rs[lead], re[lead], leadPos));
            nr.put("links", targetGenes(cts, ck, rs[lead], re[lead], b38));
            out.put("nearest", nr);
        }
        // credible-set SNPs inside a region: which genes they would regulate, weighted by PIP
        Map<String, double[]> geneScore = new TreeMap<>();          // gene -> [sum PIP x E2G, best E2G, n SNPs]
        List<Object> csOut = new ArrayList<>();
        for (double[] sn : csSnps) {
            long pos = (long) sn[0];
            int k = nearestRegion(rs, re, pos);
            if (k < 0 || dist(rs[k], re[k], pos) != 0) continue;
            List<Object> links = targetGenes(cts, ck, rs[k], re[k], b38);
            Map<String, Object> c = new LinkedHashMap<>();
            c.put("pos", pos); c.put("pip", round(sn[1])); c.put("region", k); c.put("links", links);
            csOut.add(c);
            for (Object lo : links) {
                Map<String, Object> l = MiniJson.asObject(lo);
                double e2g = ((Number) l.get("e2g")).doubleValue();
                double[] gs = geneScore.computeIfAbsent((String) l.get("gene"), x -> new double[3]);
                gs[0] += sn[1] * e2g; gs[1] = Math.max(gs[1], e2g); gs[2]++;
            }
        }
        out.put("cs_in_elements", csOut);
        List<Object> gsOut = new ArrayList<>();
        geneScore.entrySet().stream().sorted((a, b) -> Double.compare(b.getValue()[0], a.getValue()[0])).forEach(en -> {
            Map<String, Object> g = new LinkedHashMap<>();
            g.put("gene", en.getKey()); g.put("score", round(en.getValue()[0])); g.put("best_e2g", round(en.getValue()[1]));
            g.put("snps", (int) en.getValue()[2]);
            gsOut.add(g);
        });
        out.put("cs_gene_scores", gsOut);
        // the nearest element to the lead SNP in each cell type
        List<Object> strip = new ArrayList<>();
        int[] idx = new int[1];
        for (CellType ct : cts) {
            Chrom c = load(ct).get(ck);
            Map<String, Object> s = new LinkedHashMap<>();
            s.put("ct", ct.label); s.put("group", ct.group); s.put("lineage", ct.lineage);
            long d = c == null ? Long.MAX_VALUE : nearest(b38 ? c.el38 : c.el37, leadPos, idx);
            if (c != null && idx[0] >= 0 && d != Long.MAX_VALUE) {
                int[] el = (b38 ? c.el38 : c.el37)[idx[0]];
                s.put("dist", d); s.put("cls", CLASSES[el[2]]);
                int br = bestRow(c, el, b38);
                s.put("gene", geneName(c.gene[br])); s.put("e2g", round(c.e2g[br]));
            }
            strip.add(s);
        }
        out.put("cell_types", strip);
        return out;
    }

    /** The best-scoring link row of the element [s, e]: its rows are contiguous from firstRow. */
    static int bestRow(Chrom c, int[] el, boolean b38) {
        int best = el[3];
        for (int r = el[3]; r < c.n && (b38 ? c.s38[r] : c.s37[r]) == el[0] && (b38 ? c.e38[r] : c.e37[r]) == el[1]; r++)
            if (c.e2g[r] > c.e2g[best]) best = r;
        return best;
    }

    private static int nearestRegion(long[] rs, long[] re, long pos) {
        int best = -1;
        long bd = Long.MAX_VALUE;
        for (int k = 0; k < rs.length; k++) {
            long d = dist(rs[k], re[k], pos);
            if (d < bd) { bd = d; best = k; }
        }
        return best;
    }

    static long dist(long s, long e, long pos) {
        if (s < pos && pos <= e) return 0;
        return pos <= s ? s - pos + 1 : pos - e;
    }

    private static float maxScore(CellType ct, String ck, int s, int e, boolean b38) {
        Chrom c = load(ct).get(ck);
        float best = 0;
        if (c == null) return best;
        for (int r = 0; r < c.n; r++)
            if ((b38 ? c.s38[r] : c.s37[r]) == s && (b38 ? c.e38[r] : c.e37[r]) == e) best = Math.max(best, c.e2g[r]);
        return best;
    }

    /** Target genes of a region, best link per gene across the selected cell types, strongest first. */
    static List<Object> targetGenes(List<CellType> cts, String ck, long s, long e, boolean b38) {
        Map<String, Object[]> best = new HashMap<>();                 // gene -> [e2g, abc, kendall, tss, n, topCt, self]
        for (CellType ct : cts) {
            Chrom c = load(ct).get(ck);
            if (c == null) continue;
            Set<Integer> genesHere = new HashSet<>();
            for (int r = 0; r < c.n; r++) {
                int rs = b38 ? c.s38[r] : c.s37[r], re = b38 ? c.e38[r] : c.e37[r];
                if (re <= s || rs >= e) continue;
                String g = geneName(c.gene[r]);
                Object[] b = best.get(g);
                if (b == null) { b = new Object[]{0f, 0f, 0f, 0, 0, "", false}; best.put(g, b); }
                if (genesHere.add(c.gene[r])) b[4] = (int) b[4] + 1;
                if (c.e2g[r] > (float) b[0]) {
                    b[0] = c.e2g[r]; b[1] = c.abc[r]; b[2] = c.kendall[r]; b[3] = b38 ? c.tss38[r] : c.tss37[r];
                    b[5] = ct.label; b[6] = c.self[r];
                }
            }
        }
        List<Map.Entry<String, Object[]>> list = new ArrayList<>(best.entrySet());
        list.sort((a, b) -> Float.compare((float) b.getValue()[0], (float) a.getValue()[0]));
        List<Object> out = new ArrayList<>();
        for (Map.Entry<String, Object[]> en : list) {
            Object[] b = en.getValue();
            Map<String, Object> l = new LinkedHashMap<>();
            l.put("gene", en.getKey()); l.put("e2g", round((float) b[0])); l.put("abc", round((float) b[1]));
            l.put("kendall", round((float) b[2])); l.put("tss", b[3]); l.put("n", b[4]); l.put("top_ct", b[5]);
            l.put("self_promoter", b[6]);
            out.add(l);
        }
        return out;
    }

    static double round(double v) {
        if (Double.isNaN(v) || Double.isInfinite(v)) return 0;
        return Math.round(v * 10000.0) / 10000.0;
    }

    // ── Project statistics ──────────────────────────────────────────────

    /** A lead SNP (or any query SNP) for the statistics. */
    public static class Snp {
        public String locus, id, chr, gene;
        public long pos;
        public Snp(String locus, String id, String chr, long pos, String gene) {
            this.locus = locus; this.id = id; this.chr = chr; this.pos = pos; this.gene = gene;
        }
    }

    static final long[] BANDS = {0, 1_000, 10_000, 50_000, 250_000};

    /**
     * Where the SNPs sit relative to the elements of each selected cell type, and enrichment against nNull sets of random
     * SNPs matched on the signed distance to the nearest gene TSS (tssByChr: sorted TSS positions per chromosome, the
     * project's own gene annotation). One-sided empirical P = (1 + #null >= observed) / (1 + nNull); q by Benjamini-Hochberg
     * across cell types, per definition (inside, within 1 kb, within 10 kb).
     */
    public static Map<String, Object> stats(List<CellType> cts, boolean b38, List<Snp> snps, Map<String, long[]> tssByChr,
                                            int nNull, long seed) {
        loadAll(cts);
        int nS = snps.size(), nC = cts.size();
        List<Map<String, Chrom>> data = new ArrayList<>();
        for (CellType ct : cts) data.add(load(ct));
        // observed distances, per SNP x cell type
        long[][] obs = new long[nS][nC];
        int[] idx = new int[1];
        List<Object> perSnp = new ArrayList<>();
        int[] classCount = new int[3];
        long[] minDist = new long[nS];
        for (int i = 0; i < nS; i++) {
            Snp s = snps.get(i);
            String ck = chrKey(s.chr);
            long best = Long.MAX_VALUE;
            int bestCls = -1, nInside = 0;
            List<Object> inside = new ArrayList<>();
            String bestGene = "";
            float bestE2g = 0;
            for (int c = 0; c < nC; c++) {
                Chrom ch = data.get(c).get(ck);
                long d = ch == null ? Long.MAX_VALUE : nearest(b38 ? ch.el38 : ch.el37, s.pos, idx);
                obs[i][c] = d;
                if (d == 0) { nInside++; inside.add(cts.get(c).label); }
                if (ch != null && idx[0] >= 0 && d < best) {
                    best = d;
                    int[] el = (b38 ? ch.el38 : ch.el37)[idx[0]];
                    bestCls = el[2];
                    int r = bestRow(ch, el, b38);
                    bestGene = geneName(ch.gene[r]); bestE2g = ch.e2g[r];
                }
            }
            minDist[i] = best;
            if (bestCls >= 0) classCount[bestCls]++;
            Map<String, Object> m = new LinkedHashMap<>();
            m.put("locus", s.locus); m.put("snp", s.id); m.put("chr", s.chr); m.put("pos", s.pos); m.put("gene", s.gene);
            m.put("dist", best == Long.MAX_VALUE ? null : best); m.put("cls", bestCls < 0 ? "" : CLASSES[bestCls]);
            m.put("n_ct_inside", nInside); m.put("cts_inside", inside);
            m.put("nearest_gene_link", bestGene); m.put("nearest_e2g", round(bestE2g));
            perSnp.add(m);
        }
        // null: each SNP keeps its signed offset to its nearest TSS, placed at a random TSS of the genome
        List<String> tssChr = new ArrayList<>();
        List<long[]> tssPos = new ArrayList<>();
        long total = 0;
        for (Map.Entry<String, long[]> en : tssByChr.entrySet()) {
            if (en.getValue().length == 0 || !data.isEmpty() && !anyChrom(data, chrKey(en.getKey()))) continue;
            tssChr.add(en.getKey()); tssPos.add(en.getValue()); total += en.getValue().length;
        }
        long[] cum = new long[tssPos.size()];
        long acc = 0;
        for (int k = 0; k < cum.length; k++) { acc += tssPos.get(k).length; cum[k] = acc; }
        long[] offset = new long[nS];
        for (int i = 0; i < nS; i++) {
            long[] t = tssByChr.get(snps.get(i).chr);
            offset[i] = nearestSigned(t, snps.get(i).pos);
        }
        Random rnd = new Random(seed);
        int D = 3;                                                       // inside, <= 1 kb, <= 10 kb
        long[] lim = {0, 1_000, 10_000};
        int[][] obsCount = new int[nC + 1][D];                           // last row: any cell type
        for (int i = 0; i < nS; i++) {
            boolean[] any = new boolean[D];
            for (int c = 0; c < nC; c++)
                for (int d = 0; d < D; d++) if (obs[i][c] <= lim[d]) { obsCount[c][d]++; any[d] = true; }
            for (int d = 0; d < D; d++) if (any[d]) obsCount[nC][d]++;
        }
        int[][][] nullCount = new int[nNull][nC + 1][D];
        long[] nullMedian = new long[nNull];
        List<Long> nullDistPool = new ArrayList<>();
        if (total > 0) {
            long[] setDist = new long[nS];
            for (int r = 0; r < nNull; r++) {
                for (int i = 0; i < nS; i++) {
                    long pos = -1;
                    String chr = null;
                    for (int tries = 0; tries < 20 && pos <= 0; tries++) {
                        long pick = (long) (rnd.nextDouble() * total);
                        int k = Arrays.binarySearch(cum, pick + 1);
                        if (k < 0) k = -k - 1;
                        long[] t = tssPos.get(k);
                        long base = t[(int) (pick - (k == 0 ? 0 : cum[k - 1]))];
                        chr = tssChr.get(k);
                        pos = base + offset[i];
                    }
                    String ck = chrKey(chr);
                    boolean[] any = new boolean[D];
                    long mn = Long.MAX_VALUE;
                    for (int c = 0; c < nC; c++) {
                        Chrom ch = data.get(c).get(ck);
                        long d = ch == null ? Long.MAX_VALUE : nearest(b38 ? ch.el38 : ch.el37, pos, null);
                        mn = Math.min(mn, d);
                        for (int dd = 0; dd < D; dd++) if (d <= lim[dd]) { nullCount[r][c][dd]++; any[dd] = true; }
                    }
                    for (int dd = 0; dd < D; dd++) if (any[dd]) nullCount[r][nC][dd]++;
                    setDist[i] = mn;
                    if (r < 50) nullDistPool.add(mn);
                }
                long[] sd = setDist.clone();
                Arrays.sort(sd);
                nullMedian[r] = sd[nS / 2];
            }
        }
        String[] defs = {"inside", "kb1", "kb10"};
        List<Object> ctOut = new ArrayList<>();
        double[][] p = new double[nC + 1][D];
        double[][] fold = new double[nC + 1][D];
        double[][] exp = new double[nC + 1][D];
        for (int c = 0; c <= nC; c++) for (int d = 0; d < D; d++) {
            int ge = 0;
            double sum = 0;
            for (int r = 0; r < nNull; r++) { sum += nullCount[r][c][d]; if (nullCount[r][c][d] >= obsCount[c][d]) ge++; }
            exp[c][d] = nNull == 0 ? 0 : sum / nNull;
            fold[c][d] = exp[c][d] > 0 ? obsCount[c][d] / exp[c][d] : 0;
            p[c][d] = nNull == 0 ? 1 : (1.0 + ge) / (1.0 + nNull);
        }
        double[][] q = new double[nC][D];
        for (int d = 0; d < D; d++) {
            double[] pv = new double[nC];
            for (int c = 0; c < nC; c++) pv[c] = p[c][d];
            double[] qq = bh(pv);
            for (int c = 0; c < nC; c++) q[c][d] = qq[c];
        }
        for (int c = 0; c <= nC; c++) {
            Map<String, Object> m = new LinkedHashMap<>();
            if (c < nC) {
                m.put("id", cts.get(c).id); m.put("label", cts.get(c).label); m.put("group", cts.get(c).group);
                m.put("lineage", cts.get(c).lineage);
            } else { m.put("id", "any"); m.put("label", "Any selected cell type"); }
            for (int d = 0; d < D; d++) {
                Map<String, Object> v = new LinkedHashMap<>();
                v.put("obs", obsCount[c][d]); v.put("exp", round(exp[c][d])); v.put("fold", round(fold[c][d]));
                v.put("p", p[c][d]); if (c < nC) v.put("q", round(q[c][d]));
                m.put(defs[d], v);
            }
            if (c < nC) ctOut.add(m);
            else ctOut.add(0, m);
        }
        // where the SNPs sit: distance bands, class of the nearest element, distance distribution vs null
        int[] band = new int[BANDS.length + 1];
        for (long d : minDist) {
            int b = BANDS.length;
            for (int k = 0; k < BANDS.length; k++) if (d <= BANDS[k]) { b = k; break; }
            band[b]++;
        }
        long[] md = minDist.clone();
        Arrays.sort(md);
        Collections.sort(nullDistPool);
        Map<String, Object> out = new LinkedHashMap<>();
        out.put("resource", "scE2G");
        out.put("citation", citation());
        out.put("threshold", thresholdText());
        out.put("build", b38 ? "GRCh38" : "GRCh37");
        out.put("n_snps", nS);
        out.put("n_cell_types", nC);
        out.put("n_null", nNull);
        out.put("null_method", "random SNPs matched on the signed distance to the nearest gene TSS (each lead SNP's offset placed at a random TSS of the genome)");
        out.put("cell_types", ctOut);
        Map<String, Object> bands = new LinkedHashMap<>();
        String[] bl = {"inside", "1 bp - 1 kb", "1 - 10 kb", "10 - 50 kb", "50 - 250 kb", "> 250 kb"};
        for (int k = 0; k < bl.length; k++) bands.put(bl[k], band[k]);
        out.put("bands", bands);
        Map<String, Object> cls = new LinkedHashMap<>();
        for (int k = 0; k < 3; k++) cls.put(CLASSES[k], classCount[k]);
        out.put("nearest_class", cls);
        out.put("median_dist", nS == 0 ? null : md[nS / 2]);
        long[] nm = nullMedian.clone();
        Arrays.sort(nm);
        out.put("median_dist_null", nNull == 0 ? null : nm[nNull / 2]);
        out.put("cdf_obs", quantiles(md));
        long[] pool = new long[nullDistPool.size()];
        for (int k = 0; k < pool.length; k++) pool[k] = nullDistPool.get(k);
        out.put("cdf_null", quantiles(pool));
        out.put("snps", perSnp);
        return out;
    }

    private static boolean anyChrom(List<Map<String, Chrom>> data, String ck) {
        for (Map<String, Chrom> m : data) if (m.containsKey(ck)) return true;
        return false;
    }

    /** pos minus the nearest TSS (signed); 0 when there is no TSS on the chromosome. */
    static long nearestSigned(long[] tss, long pos) {
        if (tss == null || tss.length == 0) return 0;
        int k = Arrays.binarySearch(tss, pos);
        if (k >= 0) return 0;
        k = -k - 1;
        long best = Long.MAX_VALUE;
        for (int j = Math.max(0, k - 1); j <= Math.min(tss.length - 1, k); j++)
            if (Math.abs(pos - tss[j]) < Math.abs(best)) best = pos - tss[j];
        return best;
    }

    /** log-spaced distance quantiles for the cumulative curve: share of SNPs within each distance. */
    private static List<Object> quantiles(long[] sorted) {
        List<Object> out = new ArrayList<>();
        long[] at = {0, 100, 300, 1_000, 3_000, 10_000, 30_000, 100_000, 250_000};
        for (long a : at) {
            int k = 0;
            while (k < sorted.length && sorted[k] <= a) k++;
            Map<String, Object> m = new LinkedHashMap<>();
            m.put("d", a); m.put("share", sorted.length == 0 ? 0 : round((double) k / sorted.length));
            out.add(m);
        }
        return out;
    }

    static double[] bh(double[] p) {
        int n = p.length;
        Integer[] o = new Integer[n];
        for (int i = 0; i < n; i++) o[i] = i;
        Arrays.sort(o, (a, b) -> Double.compare(p[b], p[a]));
        double[] q = new double[n];
        double prev = 1;
        for (int r = 0; r < n; r++) {
            int i = o[r];
            double v = Math.min(prev, p[i] * n / (n - r));
            q[i] = Math.min(1, v);
            prev = q[i];
        }
        return q;
    }
}
