import java.util.*;

/**
 * Locus-based enrichment testing (hand-rolled, no external math library — matches
 * StatsUtil/AnovaUtil's existing dependency-free convention): compares evidence-table values
 * for genes found inside a project's identified loci ("in-loci") against a background gene set.
 *
 * Two test types, chosen per evidence column by its already-inferred type (see
 * LocalServer's evidence upload: a column is "numeric" if every non-empty value in it parses as
 * a double, "categorical" otherwise):
 *
 *   - Categorical column -&gt; Fisher's exact test on a 2x2 contingency table
 *         [[a, b], [c, d]] = [[in-loci & category, in-loci & not-category],
 *                             [background & category, background & not-category]]
 *     computed via the exact hypergeometric tail-sum definition of the two-sided p-value (sum
 *     the probabilities of every 2x2 table with the same margins whose probability is no
 *     greater than the observed table's), using log-binomial-coefficients (backed by
 *     StatsUtil.logGamma) to avoid overflow for large gene sets.
 *
 *     "Category of interest" for columns with more than two categories: the category most
 *     frequent among the in-loci genes' own values (ties broken alphabetically) — documented
 *     here since the spec calls for a single category-vs-rest 2x2 table per column, and the
 *     in-loci set's own dominant category is the most defensible parameter-free choice of which
 *     category the run is "testing for".
 *
 *   - Numeric column -&gt; Mann-Whitney U test (rank-sum), with ties handled via average ranks
 *     and a normal approximation (with the standard tie-correction term in the variance) for the
 *     two-sided p-value, via StatsUtil.zToP.
 *
 * Background gene set: by design here, "every gene from the parsed GFF3 on the same
 * chromosome(s) as at least one of the project's identified loci" (see LocalServer's
 * projectEnrichment). This is a whole-chromosome background, not a strict foreground/complement
 * split — in-loci genes are (by construction) a subset of the background set, mirroring a
 * common gene-set-enrichment convention (foreground vs. genome/chromosome background) rather
 * than excluding the foreground from its own comparison population. This is a documented
 * simplification, not an oversight.
 */
public class EnrichmentAnalyzer {

    public static class ColumnResult {
        public String column;
        public String testType;              // "fisher_exact" | "mann_whitney_u"
        public double statistic = Double.NaN;
        public double pValue = Double.NaN;
        public int nInLoci = 0;
        public int nBackground = 0;
        public String categoryOfInterest;    // only set for fisher_exact
    }

    public static class Result {
        public List<ColumnResult> columns = new ArrayList<>();

        public String toJson() {
            StringBuilder j = new StringBuilder("{\"columns\":[");
            for (int i = 0; i < columns.size(); i++) {
                if (i > 0) j.append(",");
                ColumnResult c = columns.get(i);
                j.append("{");
                kv(j, "column", c.column); j.append(",");
                kv(j, "test_type", c.testType); j.append(",");
                j.append("\"statistic\":").append(num(c.statistic)).append(",");
                j.append("\"p_value\":").append(num(c.pValue)).append(",");
                j.append("\"n_in_loci\":").append(c.nInLoci).append(",");
                j.append("\"n_background\":").append(c.nBackground);
                if (c.categoryOfInterest != null) {
                    j.append(","); kv(j, "category_of_interest", c.categoryOfInterest);
                }
                j.append("}");
            }
            j.append("]}");
            return j.toString();
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

    /**
     * @param inLociGenes     gene symbols found in this project's identified loci (duplicates
     *                        and background-overlap are fine; deduped internally via a Set)
     * @param backgroundGenes gene symbols making up the background set
     * @param evidence        gene symbol -&gt; {column name -&gt; raw string value}
     * @param columnTypes     column name -&gt; "numeric" | "categorical"
     */
    public static Result run(Collection<String> inLociGenes, Collection<String> backgroundGenes,
                              Map<String, Map<String, String>> evidence,
                              Map<String, String> columnTypes) {
        Result out = new Result();
        Set<String> inLoci = new LinkedHashSet<>(inLociGenes);
        Set<String> background = new LinkedHashSet<>(backgroundGenes);

        for (Map.Entry<String, String> ce : columnTypes.entrySet()) {
            String col = ce.getKey();
            boolean numeric = "numeric".equals(ce.getValue());
            ColumnResult cr = new ColumnResult();
            cr.column = col;
            cr.testType = numeric ? "mann_whitney_u" : "fisher_exact";

            if (numeric) {
                double[] a = valuesFor(inLoci, evidence, col);
                double[] b = valuesFor(background, evidence, col);
                cr.nInLoci = a.length;
                cr.nBackground = b.length;
                if (a.length == 0 || b.length == 0) { out.columns.add(cr); continue; }
                MannWhitneyResult mw = mannWhitneyU(a, b);
                cr.statistic = mw.u;
                cr.pValue = mw.pValue;
            } else {
                Map<String, Integer> inLociCounts = new TreeMap<>();
                Map<String, Integer> bgCounts = new HashMap<>();
                int nInLoci = 0, nBg = 0;
                for (String g : inLoci) {
                    Map<String, String> row = evidence.get(g);
                    if (row == null) continue;
                    String v = row.get(col);
                    if (v == null || v.isEmpty()) continue;
                    nInLoci++;
                    inLociCounts.merge(v, 1, Integer::sum);
                }
                for (String g : background) {
                    Map<String, String> row = evidence.get(g);
                    if (row == null) continue;
                    String v = row.get(col);
                    if (v == null || v.isEmpty()) continue;
                    nBg++;
                    bgCounts.merge(v, 1, Integer::sum);
                }
                cr.nInLoci = nInLoci;
                cr.nBackground = nBg;
                if (nInLoci == 0 || nBg == 0 || inLociCounts.isEmpty()) { out.columns.add(cr); continue; }

                String best = null; int bestCount = -1;
                for (Map.Entry<String, Integer> e2 : inLociCounts.entrySet()) {
                    if (e2.getValue() > bestCount) { bestCount = e2.getValue(); best = e2.getKey(); }
                }
                cr.categoryOfInterest = best;
                int a = inLociCounts.getOrDefault(best, 0);
                int bCell = nInLoci - a;
                int c = bgCounts.getOrDefault(best, 0);
                int dCell = nBg - c;
                cr.statistic = a;
                cr.pValue = fisherExactTwoSided(a, bCell, c, dCell);
            }
            out.columns.add(cr);
        }
        return out;
    }

    private static double[] valuesFor(Set<String> genes, Map<String, Map<String, String>> evidence, String col) {
        List<Double> vals = new ArrayList<>();
        for (String g : genes) {
            Map<String, String> row = evidence.get(g);
            if (row == null) continue;
            String v = row.get(col);
            if (v == null || v.isEmpty()) continue;
            try { vals.add(Double.parseDouble(v)); } catch (NumberFormatException ignored) {}
        }
        double[] a = new double[vals.size()];
        for (int i = 0; i < a.length; i++) a[i] = vals.get(i);
        return a;
    }

    // ── Mann-Whitney U (rank-sum), average ranks for ties, normal approximation p-value ──

    public static class MannWhitneyResult {
        public double u;       // U for group a (in-loci)
        public double pValue;
    }

    /**
     * Rank-sum (Mann-Whitney) U test. Ranks the pooled sample, breaking ties with the average
     * rank of the tied block, then derives U1 = R1 - n1(n1+1)/2 for group a. The p-value uses the
     * standard normal approximation with the tie-correction term in the variance:
     *   Var(U) = n1*n2/12 * [(N+1) - sum(t_j^3 - t_j) / (N(N-1))]
     * where t_j is the size of the j-th tied block and N = n1+n2 (Numerical Recipes /
     * standard nonparametric-statistics formula for the tie-corrected normal approximation).
     */
    static MannWhitneyResult mannWhitneyU(double[] a, double[] b) {
        int n1 = a.length, n2 = b.length;
        double[] all = new double[n1 + n2];
        System.arraycopy(a, 0, all, 0, n1);
        System.arraycopy(b, 0, all, n1, n2);
        Integer[] order = new Integer[all.length];
        for (int i = 0; i < order.length; i++) order[i] = i;
        Arrays.sort(order, (x, y) -> Double.compare(all[x], all[y]));

        double[] ranks = new double[all.length];
        double tieSum = 0;
        int i = 0;
        while (i < order.length) {
            int j = i;
            while (j + 1 < order.length && all[order[j + 1]] == all[order[i]]) j++;
            double avgRank = (i + j) / 2.0 + 1; // 1-based ranks
            for (int k = i; k <= j; k++) ranks[order[k]] = avgRank;
            double t = j - i + 1;
            tieSum += (t * t * t - t);
            i = j + 1;
        }

        double r1 = 0;
        for (int k = 0; k < n1; k++) r1 += ranks[k];
        double u1 = r1 - n1 * (n1 + 1) / 2.0;

        MannWhitneyResult r = new MannWhitneyResult();
        r.u = u1;
        double meanU = n1 * (double) n2 / 2.0;
        int N = n1 + n2;
        double varU = (N > 1)
            ? (n1 * (double) n2 / 12.0) * ((N + 1) - tieSum / (double) (N * (N - 1)))
            : 0;
        if (varU <= 0) {
            r.pValue = 1.0; // no variance to distinguish groups (e.g. every value tied)
            return r;
        }
        double z = (u1 - meanU) / Math.sqrt(varU);
        r.pValue = StatsUtil.zToP(z);
        return r;
    }

    // ── Fisher's exact test (two-sided), via exact hypergeometric tail-sum ──

    /** log(nCk) via StatsUtil.logGamma, avoiding overflow for large gene sets. */
    static double logChoose(int n, int k) {
        if (k < 0 || k > n) return Double.NEGATIVE_INFINITY;
        return StatsUtil.logGamma(n + 1) - StatsUtil.logGamma(k + 1) - StatsUtil.logGamma(n - k + 1);
    }

    /**
     * Two-sided Fisher's exact test p-value for a 2x2 table [[a,b],[c,d]], via direct summation
     * of the hypergeometric probabilities of every table sharing the same row/column margins
     * whose probability is no greater than the observed table's — the standard exact definition
     * of Fisher's two-sided p-value (not an approximation).
     */
    static double fisherExactTwoSided(int a, int b, int c, int d) {
        int row1 = a + b, row2 = c + d;
        int col1 = a + c, col2 = b + d;
        int n = row1 + row2;
        if (row1 == 0 || row2 == 0 || col1 == 0 || col2 == 0) return 1.0;

        double logDenom = logChoose(n, col1);
        double logObserved = logChoose(row1, a) + logChoose(row2, c) - logDenom;

        int aMin = Math.max(0, col1 - row2);
        int aMax = Math.min(row1, col1);
        double pSum = 0;
        final double EPS = 1e-9;
        for (int ai = aMin; ai <= aMax; ai++) {
            int ci = col1 - ai;
            double logP = logChoose(row1, ai) + logChoose(row2, ci) - logDenom;
            if (logP <= logObserved + EPS) pSum += Math.exp(logP);
        }
        return Math.min(1.0, pSum);
    }
}
