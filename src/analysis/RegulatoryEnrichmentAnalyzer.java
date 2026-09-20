import java.util.*;

/**
 * Real interval-overlap enrichment test (DECISIONS_PHASE3.md section 3.3): for a project's own
 * genome-wide-significant SNPs (foreground) versus its own below-significance SNPs (background,
 * same GWAS file), does each histone mark's set of real Roadmap Epigenomics peaks
 * ({@link RegulatoryPeakIndex}) overlap the foreground more than the background?
 *
 * This is a genuinely different foreground/background definition and classifier from
 * {@link EnrichmentAnalyzer}'s existing gene-based test (a SNP's genomic position falling inside a
 * peak interval, vs. a gene's evidence-table category/value), so this is new code — it does not
 * touch EnrichmentAnalyzer's own per-gene logic. It does directly reuse
 * {@link EnrichmentAnalyzer#fisherExactTwoSided(int, int, int, int)} for the 2x2 significance test,
 * since that is already a real, tested, exact hypergeometric implementation and this project's
 * convention (StatsUtil/AnovaUtil) is not to reimplement shared math a second time.
 */
public class RegulatoryEnrichmentAnalyzer {

    /** A single SNP's genomic position (chromosome + 1-based position), independent of its p-value. */
    public static class SnpPos {
        public final String chr;
        public final long pos;
        public SnpPos(String chr, long pos) { this.chr = chr; this.pos = pos; }
    }

    public static class MarkResult {
        public String mark;
        public double oddsRatio = Double.NaN;
        public double pValue = Double.NaN;
        public int nForegroundInPeak;
        public int nForegroundTotal;
        public int nBackgroundInPeak;
        public int nBackgroundTotal;

        public String toJson() {
            StringBuilder j = new StringBuilder("{");
            kv(j, "mark", mark); j.append(",");
            j.append("\"odds_ratio\":").append(num(oddsRatio)).append(",");
            j.append("\"p_value\":").append(num(pValue)).append(",");
            j.append("\"n_foreground_in_peak\":").append(nForegroundInPeak).append(",");
            j.append("\"n_foreground_total\":").append(nForegroundTotal).append(",");
            j.append("\"n_background_in_peak\":").append(nBackgroundInPeak).append(",");
            j.append("\"n_background_total\":").append(nBackgroundTotal);
            j.append("}");
            return j.toString();
        }
    }

    public static class Result {
        public String eid;          // null when the project's disease has no mapped reference epigenome
        public String tissue;       // null alongside eid
        public double threshold;
        public List<MarkResult> marks = new ArrayList<>();

        public String toJson() {
            StringBuilder j = new StringBuilder("{");
            j.append("\"eid\":").append(eid == null ? "null" : "\"" + esc(eid) + "\"").append(",");
            j.append("\"tissue\":").append(tissue == null ? "null" : "\"" + esc(tissue) + "\"").append(",");
            j.append("\"threshold\":").append(threshold).append(",");
            j.append("\"marks\":[");
            for (int i = 0; i < marks.size(); i++) {
                if (i > 0) j.append(",");
                j.append(marks.get(i).toJson());
            }
            j.append("]}");
            return j.toString();
        }
    }

    /**
     * Runs the test for one project's already-mapped reference epigenome. Returns a Result with an
     * empty {@code marks} list (not an error) when {@code eid} is null, matching the "unmapped
     * disease = no regulatory data" contract used everywhere else in this feature.
     */
    public static Result run(String eid, String tissue, double threshold,
                              List<SnpPos> foreground, List<SnpPos> background,
                              RegulatoryPeakIndex index) {
        Result out = new Result();
        out.eid = eid;
        out.tissue = tissue;
        out.threshold = threshold;
        if (eid == null) return out;

        for (String mark : RegulatoryPeakIndex.MARKS) {
            MarkResult mr = new MarkResult();
            mr.mark = mark;
            mr.nForegroundTotal = foreground.size();
            mr.nBackgroundTotal = background.size();
            mr.nForegroundInPeak = countInPeak(foreground, eid, mark, index);
            mr.nBackgroundInPeak = countInPeak(background, eid, mark, index);

            int a = mr.nForegroundInPeak;
            int b = mr.nForegroundTotal - a;
            int c = mr.nBackgroundInPeak;
            int d = mr.nBackgroundTotal - c;
            mr.pValue = EnrichmentAnalyzer.fisherExactTwoSided(a, b, c, d);
            mr.oddsRatio = oddsRatio(a, b, c, d);
            out.marks.add(mr);
        }
        return out;
    }

    private static int countInPeak(List<SnpPos> snps, String eid, String mark, RegulatoryPeakIndex index) {
        int n = 0;
        for (SnpPos s : snps) {
            if (!index.peaksOverlapping(eid, mark, s.chr, s.pos, s.pos).isEmpty()) n++;
        }
        return n;
    }

    /**
     * (a*d)/(b*c), with the standard Haldane-Anscombe correction (+0.5 to every cell) when any cell
     * is zero — the usual, defensible way to avoid an undefined (0/0) or infinite (b=0 or c=0) odds
     * ratio for a sparse table, rather than special-casing each zero combination separately.
     */
    static double oddsRatio(int a, int b, int c, int d) {
        if (a == 0 || b == 0 || c == 0 || d == 0) {
            return ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5));
        }
        return ((double) a * d) / ((double) b * c);
    }

    private static void kv(StringBuilder j, String key, String value) {
        j.append("\"").append(key).append("\":\"").append(esc(value)).append("\"");
    }
    private static String num(double d) {
        return Double.isNaN(d) || Double.isInfinite(d) ? "null" : String.valueOf(d);
    }
    private static String esc(String s) {
        if (s == null) return "";
        return s.replace("\\", "\\\\").replace("\"", "\\\"");
    }
}
