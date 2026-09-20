import java.io.*;
import java.util.*;
import java.util.concurrent.ConcurrentHashMap;
import java.util.zip.GZIPInputStream;

/**
 * Process-lifetime, lazily-built-per-(EID,mark), in-memory interval index over the NIH Roadmap
 * Epigenomics consolidated narrowPeak histone ChIP-seq peak files at
 * {@code regulatory_data/<EID>/<EID>-<mark>.narrowPeak.gz} (see DECISIONS_PHASE3.md section 1/2).
 *
 * Mirrors the interval-lookup convention already used by {@link GwasParser}/{@link MultiLocusScanner}
 * for GWAS-SNP-to-locus matching: peaks are grouped per chromosome (normalized via
 * {@code MultiLocusScanner.normalizeChr} so "1"/"chr1" key the same bucket), sorted by start
 * position, then queried by binary search for the first peak whose start exceeds the query
 * window's end, followed by a linear scan back over the candidates for genomic overlap with the
 * window's start — the same two-step pattern those classes use for a single point, generalized
 * here to a start/end window.
 *
 * Storage-conscious by design (DECISIONS_PHASE3.md section 2): the gzipped narrowPeak files stay
 * on disk exactly as downloaded; each (EID, mark) pair is decompressed and parsed into memory only
 * once, on first request, and cached for the rest of the process's lifetime — never pre-expanded
 * into a larger on-disk form, and never duplicated per project (one shared copy per tissue,
 * resolved by disease at query time).
 */
public class RegulatoryPeakIndex {

    /** One narrowPeak record: genomic span plus its real (not synthetic) ChIP-seq signal value. */
    public static class Peak {
        public final String chr;
        public final long start;
        public final long end;
        public final double signalValue;

        public Peak(String chr, long start, long end, double signalValue) {
            this.chr = chr;
            this.start = start;
            this.end = end;
            this.signalValue = signalValue;
        }
    }

    // ── Disease -> reference epigenome (EID) mapping, from DECISIONS_PHASE3.md section 1. Derive a
    // project's disease key the same best-effort way GeneConstellationBuilder.diseaseGroupOf already
    // does (substring before the project id's first '-', else "ungrouped") before looking it up here.
    // A disease key with no entry (including "ungrouped") means "no regulatory data" — empty, not
    // an error — since this table only covers the 30-dataset/6-disease validation corpus by design.
    public static final Map<String, String> DISEASE_TO_EID;
    public static final Map<String, String> EID_TO_TISSUE;
    /** The three histone marks downloaded for every mapped tissue, in the order tracks are shown. */
    public static final String[] MARKS = {"H3K27ac", "H3K4me1", "H3K4me3"};

    static {
        Map<String, String> d2e = new LinkedHashMap<>();
        d2e.put("scz", "E073");
        d2e.put("alz", "E073");
        d2e.put("t2d", "E098");
        d2e.put("cad", "E065");
        d2e.put("ra",  "E034");
        d2e.put("ibd", "E106");
        DISEASE_TO_EID = Collections.unmodifiableMap(d2e);

        Map<String, String> e2t = new LinkedHashMap<>();
        e2t.put("E073", "Brain, Dorsolateral Prefrontal Cortex");
        e2t.put("E098", "Pancreas");
        e2t.put("E065", "Aorta");
        e2t.put("E034", "Primary T cells, peripheral blood");
        e2t.put("E106", "Sigmoid colon");
        EID_TO_TISSUE = Collections.unmodifiableMap(e2t);
    }

    /** Resolves a project id's EID via its disease group, or null if unmapped (empty tracks/enrichment). */
    public static String resolveEid(String projectId) {
        String disease = GeneConstellationBuilder.diseaseGroupOf(projectId);
        return DISEASE_TO_EID.get(disease);
    }

    private static final String DEFAULT_DATA_ROOT = "regulatory_data";
    private static final RegulatoryPeakIndex INSTANCE = new RegulatoryPeakIndex(DEFAULT_DATA_ROOT);

    /** The process-lifetime shared instance the server uses; tests build their own with a fixture root. */
    public static RegulatoryPeakIndex instance() { return INSTANCE; }

    private final String dataRoot;
    // (eid + "|" + mark) -> (normalized chr -> peaks sorted by start). An empty map is cached for a
    // missing/corrupt file too, so a bad input is never re-parsed on every subsequent request.
    private final Map<String, Map<String, List<Peak>>> cache = new ConcurrentHashMap<>();

    public RegulatoryPeakIndex(String dataRoot) {
        this.dataRoot = dataRoot;
    }

    /**
     * Peaks of {@code mark} for reference epigenome {@code eid} overlapping [start, end] on
     * {@code chr} (inclusive on both ends, so a peak exactly touching a query window edge counts).
     * Returns an empty list, never throws, for an unknown EID/mark, a missing/corrupt file, or a
     * chromosome with no peaks — building/parsing happens lazily on first use and is cached after.
     */
    public List<Peak> peaksOverlapping(String eid, String mark, String chr, long start, long end) {
        if (eid == null || mark == null || chr == null) return Collections.emptyList();
        Map<String, List<Peak>> byChr = cache.computeIfAbsent(key(eid, mark), k -> build(eid, mark));
        List<Peak> chrPeaks = byChr.get(MultiLocusScanner.normalizeChr(chr));
        if (chrPeaks == null || chrPeaks.isEmpty()) return Collections.emptyList();

        // Binary search: first index whose start exceeds the query window's end.
        int bsLo = 0, bsHi = chrPeaks.size();
        while (bsLo < bsHi) {
            int mid = (bsLo + bsHi) >>> 1;
            if (chrPeaks.get(mid).start <= end) bsLo = mid + 1; else bsHi = mid;
        }
        // Every peak at indices [0..bsLo-1] starts at or before the window's end; keep those whose
        // own end also reaches into the window (i.e. genuinely overlaps, not just starts before it).
        List<Peak> hits = new ArrayList<>();
        for (int i = 0; i < bsLo; i++) {
            Peak p = chrPeaks.get(i);
            if (p.end >= start) hits.add(p);
        }
        return hits;
    }

    private static String key(String eid, String mark) { return eid + "|" + mark; }

    private Map<String, List<Peak>> build(String eid, String mark) {
        File f = new File(dataRoot, eid + File.separator + eid + "-" + mark + ".narrowPeak.gz");
        Map<String, List<Peak>> byChr = new HashMap<>();
        if (!f.isFile()) {
            System.err.printf("[RegulatoryPeakIndex] Missing peak file (treating as no data): %s%n", f.getPath());
            return byChr;
        }
        int kept = 0, skipped = 0;
        try (BufferedReader br = new BufferedReader(new InputStreamReader(
                new GZIPInputStream(new FileInputStream(f)), "UTF-8"), 1 << 20)) {
            String line;
            while ((line = br.readLine()) != null) {
                if (line.isEmpty()) continue;
                // narrowPeak (BED6+4): chrom, start, end, name, score, strand, signalValue, pValue, qValue, peak
                String[] cols = line.split("\t");
                if (cols.length < 7) { skipped++; continue; }
                try {
                    String chr = MultiLocusScanner.normalizeChr(cols[0].trim());
                    long start = Long.parseLong(cols[1].trim());
                    long end = Long.parseLong(cols[2].trim());
                    double signal = Double.parseDouble(cols[6].trim());
                    byChr.computeIfAbsent(chr, k -> new ArrayList<>()).add(new Peak(chr, start, end, signal));
                    kept++;
                } catch (NumberFormatException nfe) {
                    skipped++;
                }
            }
        } catch (IOException e) {
            System.err.printf("[RegulatoryPeakIndex] Failed to read/parse %s (treating as no data): %s%n",
                f.getPath(), e.getMessage());
            return new HashMap<>();
        }
        for (List<Peak> list : byChr.values()) {
            list.sort(Comparator.comparingLong(p -> p.start));
        }
        System.out.printf("[RegulatoryPeakIndex] Built %s/%s: %d peaks kept, %d skipped, %d chromosomes%n",
            eid, mark, kept, skipped, byChr.size());
        return byChr;
    }
}
