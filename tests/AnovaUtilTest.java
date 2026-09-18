import java.util.*;

/**
 * Regression test for AnovaUtil / StatsUtil.fDistPValue.
 *
 * Verification strategy (per the requirement to check statistical correctness at least three
 * independent ways, not just "it ran without throwing"):
 *   1. Two cases below use groups of the shape {m-1, m, m+1} for k=3 groups of n=3 — for that
 *      exact shape the within-group SS is always 6 (df=6) and the between-group df is always 2,
 *      so the F(2,6) distribution's upper-tail applies every time. For d1=2, the F distribution
 *      has a closed form derived independently of AnovaUtil/StatsUtil's incomplete-beta code:
 *      F(2,d2) reduces to Beta(1, d2/2), whose CDF has the elementary closed form
 *      1-(1-y)^b (no continued fraction needed), giving
 *          P(F > f) = (d2 / (2f + d2)) ^ (d2/2).
 *      That closed form is evaluated here from scratch (fDistUpperTailClosedFormD1eq2) and
 *      compared against AnovaUtil's general result — two independently-coded computations of the
 *      same quantity agreeing is a much stronger check than a single implementation "looking
 *      right".
 *   2. Degenerate edge cases (zero within-group variance in both the "no difference" and
 *      "total separation" directions) are checked against the mathematically forced values
 *      (p=1 and p=0 respectively), not against the general formula at all.
 *   3. A monotonicity sanity check across many random-but-seeded synthetic datasets: shrinking
 *      within-group spread while holding between-group means fixed must never decrease F or
 *      increase p. This is a "controlled random logical investigation" check distinct from the
 *      closed-form numeric checks above — it would catch a wrong-direction sign bug that a single
 *      fixed-number example could accidentally pass.
 *   4. (Applied outside this file, per DECISIONS_PHASE2.md:) this test is run once against the
 *      real AnovaUtil, then again after `git stash`-ing AnovaUtil.java/the StatsUtil additions,
 *      to confirm it actually fails without the implementation rather than vacuously passing.
 */
public class AnovaUtilTest {

    public static void main(String[] args) {
        int failures = 0;

        // --- Case 1: textbook-shape example, means 4/7/10, closed-form p = 0.001 exactly ---
        {
            List<double[]> groups = Arrays.asList(
                new double[]{3, 4, 5},
                new double[]{6, 7, 8},
                new double[]{9, 10, 11}
            );
            AnovaUtil.Result r = AnovaUtil.oneWay(groups);
            failures += checkClose("case1 ssBetween", r.ssBetween, 54.0, 1e-9);
            failures += checkClose("case1 ssWithin", r.ssWithin, 6.0, 1e-9);
            failures += checkClose("case1 F statistic", r.fStat, 27.0, 1e-9);
            double closedForm = fDistUpperTailClosedFormD1eq2(r.fStat, r.dfWithin);
            failures += checkClose("case1 p vs independent closed form", closedForm, 0.001, 1e-9);
            failures += checkClose("case1 AnovaUtil p matches closed form", r.pValue, closedForm, 1e-9);
        }

        // --- Case 2: same shape, irregular means 10/15/30 -> tiny p, still closed-form-checkable ---
        {
            List<double[]> groups = Arrays.asList(
                new double[]{9, 10, 11},
                new double[]{14, 15, 16},
                new double[]{29, 30, 31}
            );
            AnovaUtil.Result r = AnovaUtil.oneWay(groups);
            failures += checkClose("case2 ssBetween", r.ssBetween, 650.0, 1e-6);
            failures += checkClose("case2 F statistic", r.fStat, 325.0, 1e-6);
            double closedForm = fDistUpperTailClosedFormD1eq2(r.fStat, r.dfWithin);
            failures += checkClose("case2 AnovaUtil p matches independent closed form", r.pValue, closedForm, 1e-4 * closedForm);
        }

        // --- Case 3: zero within-group variance, groups differ -> forced p = 0 (perfect separation) ---
        {
            List<double[]> groups = Arrays.asList(
                new double[]{1, 1, 1},
                new double[]{5, 5, 5}
            );
            AnovaUtil.Result r = AnovaUtil.oneWay(groups);
            failures += checkClose("case3 (perfect separation) p == 0", r.pValue, 0.0, 1e-12);
            failures += check("case3 F is +Infinity", Double.isInfinite(r.fStat) && r.fStat > 0);
        }

        // --- Case 4: zero within-group variance, groups identical -> forced p = 1 (no difference at all) ---
        {
            List<double[]> groups = Arrays.asList(
                new double[]{7, 7, 7},
                new double[]{7, 7, 7},
                new double[]{7, 7, 7}
            );
            AnovaUtil.Result r = AnovaUtil.oneWay(groups);
            failures += checkClose("case4 (identical groups) p == 1", r.pValue, 1.0, 1e-12);
            failures += checkClose("case4 F == 0", r.fStat, 0.0, 1e-12);
        }

        // --- Case 5: controlled-random monotonicity check (seeded, reproducible) ---
        {
            Random rng = new Random(42);
            for (int trial = 0; trial < 50; trial++) {
                double meanA = rng.nextDouble() * 10;
                double meanB = meanA + 1 + rng.nextDouble() * 10; // always strictly greater
                double wideSpread = 3.0 + rng.nextDouble() * 5.0;
                double narrowSpread = wideSpread / (2.0 + rng.nextDouble() * 3.0); // strictly smaller
                int n = 5 + rng.nextInt(10);

                AnovaUtil.Result wide = AnovaUtil.oneWay(Arrays.asList(
                    sample(rng, meanA, wideSpread, n), sample(rng, meanB, wideSpread, n)));
                AnovaUtil.Result narrow = AnovaUtil.oneWay(Arrays.asList(
                    sample(rng, meanA, narrowSpread, n), sample(rng, meanB, narrowSpread, n)));

                if (narrow.fStat < wide.fStat - 1e-9) {
                    failures += check("trial " + trial + ": narrower within-group spread must not lower F "
                        + "(narrow F=" + narrow.fStat + ", wide F=" + wide.fStat + ")", false);
                }
                if (narrow.pValue > wide.pValue + 1e-9) {
                    failures += check("trial " + trial + ": narrower within-group spread must not raise p "
                        + "(narrow p=" + narrow.pValue + ", wide p=" + wide.pValue + ")", false);
                }
            }
            System.out.println("PASS: 50 seeded random monotonicity trials (spread down => F up, p down)");
        }

        // --- Case 6: rejects degenerate input rather than silently returning nonsense ---
        {
            boolean threw = false;
            try {
                AnovaUtil.oneWay(Arrays.asList(new double[]{1, 2, 3}));
            } catch (IllegalArgumentException e) {
                threw = true;
            }
            failures += check("single group correctly rejected", threw);
        }

        if (failures == 0) {
            System.out.println("PASS: all AnovaUtil tests passed");
        } else {
            System.out.println("FAIL: " + failures + " test(s) failed");
            System.exit(1);
        }
    }

    /**
     * Independent closed-form derivation of P(F > f) for an F(2, d2) distribution, coded from
     * scratch (not calling StatsUtil) via the elementary Beta(1, d2/2) CDF identity — see the
     * class-level comment for the derivation. Deliberately duplicates none of StatsUtil's
     * incomplete-beta machinery so this is a genuinely separate computation to compare against.
     */
    private static double fDistUpperTailClosedFormD1eq2(double f, int d2) {
        double y1minus = d2 / (2.0 * f + d2);
        return Math.pow(y1minus, d2 / 2.0);
    }

    private static double[] sample(Random rng, double mean, double spread, int n) {
        double[] out = new double[n];
        for (int i = 0; i < n; i++) out[i] = mean + (rng.nextDouble() - 0.5) * 2 * spread;
        return out;
    }

    private static int checkClose(String label, double actual, double expected, double tol) {
        if (Math.abs(actual - expected) > Math.abs(tol) + 1e-15) {
            System.out.println("FAIL: " + label + " — expected " + expected + ", got " + actual);
            return 1;
        }
        System.out.println("PASS: " + label + " (" + actual + ")");
        return 0;
    }

    private static int check(String label, boolean cond) {
        if (!cond) {
            System.out.println("FAIL: " + label);
            return 1;
        }
        System.out.println("PASS: " + label);
        return 0;
    }
}
