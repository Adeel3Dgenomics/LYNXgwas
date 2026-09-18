import java.util.List;

/**
 * One-way ANOVA (hand-rolled, no external math library — matches StatsUtil's existing
 * dependency-free convention). Used to test whether a gene's significance varies meaningfully
 * (a) across the datasets within one disease, and (b) across the six disease groups.
 */
public class AnovaUtil {

    public static class Result {
        public int k;               // number of groups
        public int n;                // total observations
        public double grandMean;
        public double ssBetween, ssWithin, ssTotal;
        public int dfBetween, dfWithin;
        public double msBetween, msWithin;
        public double fStat;
        public double pValue;
    }

    /**
     * @param groups each element is the set of observations (e.g. -log10(p) per dataset) for one
     *                group (e.g. one disease, or one dataset within a disease). Groups with zero
     *                observations are ignored; at least 2 non-empty groups are required.
     */
    public static Result oneWay(List<double[]> groups) {
        int k = 0;
        int n = 0;
        double sum = 0;
        for (double[] g : groups) {
            if (g.length == 0) continue;
            k++;
            n += g.length;
            for (double v : g) sum += v;
        }
        if (k < 2 || n <= k) {
            throw new IllegalArgumentException("ANOVA requires at least 2 non-empty groups and n > k observations");
        }
        double grandMean = sum / n;

        double ssBetween = 0, ssWithin = 0;
        for (double[] g : groups) {
            if (g.length == 0) continue;
            double gsum = 0;
            for (double v : g) gsum += v;
            double gmean = gsum / g.length;
            ssBetween += g.length * (gmean - grandMean) * (gmean - grandMean);
            for (double v : g) ssWithin += (v - gmean) * (v - gmean);
        }

        Result r = new Result();
        r.k = k;
        r.n = n;
        r.grandMean = grandMean;
        r.ssBetween = ssBetween;
        r.ssWithin = ssWithin;
        r.ssTotal = ssBetween + ssWithin;
        r.dfBetween = k - 1;
        r.dfWithin = n - k;
        r.msBetween = ssBetween / r.dfBetween;
        r.msWithin = r.dfWithin > 0 ? ssWithin / r.dfWithin : 0;

        if (r.msWithin == 0) {
            // All within-group variance is exactly zero (degenerate/synthetic data): any
            // between-group difference is infinitely significant, no difference is undefined-but-
            // reported as non-significant rather than propagating NaN into callers/exports.
            r.fStat = r.msBetween > 0 ? Double.POSITIVE_INFINITY : 0.0;
            r.pValue = r.msBetween > 0 ? 0.0 : 1.0;
        } else {
            r.fStat = r.msBetween / r.msWithin;
            r.pValue = StatsUtil.fDistPValue(r.fStat, r.dfBetween, r.dfWithin);
        }
        return r;
    }
}
