
/**
 * Small self-contained statistics helpers (no external math library dependency,
 * matching this project's existing hand-rolled-over-dependency convention).
 */
public class StatsUtil {

    /**
     * Inverse standard normal CDF (quantile function), Acklam's rational approximation
     * (relative error &lt; 1.15e-9 everywhere). Used to back-derive beta/SE from OR + p-value
     * when a GWAS file reports only the former.
     */
    public static double qnorm(double p) {
        if (p <= 0) return Double.NEGATIVE_INFINITY;
        if (p >= 1) return Double.POSITIVE_INFINITY;

        double[] a = {-3.969683028665376e+01, 2.209460984245205e+02, -2.759285104469687e+02,
                       1.383577518672690e+02, -3.066479806614716e+01, 2.506628277459239e+00};
        double[] b = {-5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02,
                       6.680131188771972e+01, -1.328068155288572e+01};
        double[] c = {-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00,
                      -2.549732539343734e+00, 4.374664141464968e+00, 2.938163982698783e+00};
        double[] d = {7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00,
                      3.754408661907416e+00};
        double pLow = 0.02425, pHigh = 1 - pLow;

        if (p < pLow) {
            double q = Math.sqrt(-2 * Math.log(p));
            return (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) /
                   ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1);
        } else if (p <= pHigh) {
            double q = p - 0.5, r = q * q;
            return (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q /
                   (((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1);
        } else {
            double q = Math.sqrt(-2 * Math.log(1 - p));
            return -((((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) /
                    ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1));
        }
    }

    /**
     * Complementary error function via the Numerical Recipes Chebyshev approximation
     * (fractional error &lt; 1.2e-7 everywhere). Unlike computing 1-CDF by subtraction, this
     * evaluates the tail directly (through exp(-z^2+...)) so it stays accurate for extreme
     * z rather than catastrophically cancelling to 0 once the true tail probability is smaller
     * than double's ~1e-16 subtractive resolution.
     */
    public static double erfc(double x) {
        double z = Math.abs(x);
        double t = 1.0 / (1.0 + 0.5 * z);
        double ans = t * Math.exp(-z*z - 1.26551223 +
            t*(1.00002368 +
            t*(0.37409196 +
            t*(0.09678418 +
            t*(-0.18628806 +
            t*(0.27886807 +
            t*(-1.13520398 +
            t*(1.48851587 +
            t*(-0.82215223 +
            t*0.17087277)))))))));
        return x >= 0.0 ? ans : 2.0 - ans;
    }

    private static final double SQRT2 = Math.sqrt(2.0);

    /** Accurate two-sided p-value from a z-score (z = beta/se), for extreme z where 1-CDF underflows. */
    public static double zToP(double z) {
        return erfc(Math.abs(z) / SQRT2);
    }

    /**
     * Derives beta = ln(OR) and se = |beta| / z, z = qnorm(1 - p/2), for a GWAS row that reports
     * OR and p-value but not beta/SE directly. Returns null if any input is invalid (OR &lt;= 0,
     * p outside (0,1), non-finite) rather than propagating a bad value.
     */
    public static double[] deriveBetaSeFromOrP(double or, double pval) {
        if (!(or > 0) || !Double.isFinite(or)) return null;
        if (!(pval > 0 && pval < 1) || !Double.isFinite(pval)) return null;
        double beta = Math.log(or);
        double z = qnorm(1.0 - pval / 2.0);
        if (!Double.isFinite(z) || z == 0) return null;
        double se = Math.abs(beta) / z;
        if (!Double.isFinite(se) || se <= 0) return null;
        return new double[]{beta, se};
    }
}
