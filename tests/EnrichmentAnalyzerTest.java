import java.util.*;

/**
 * Regression test for EnrichmentAnalyzer (Fisher's exact test + Mann-Whitney U test).
 *
 * Verification strategy (independent cross-checks, not just "it ran without throwing"):
 *   1. Fisher's exact test against the classic "lady tasting tea" 2x2 example (Fisher, 1935):
 *      an 8-cup experiment, 4 milk-first / 4 tea-first, table [[3,1],[1,3]]. The two-sided exact
 *      p-value for this table is a widely-cited textbook value, independently re-derivable by
 *      hand from the hypergeometric distribution: p = (C(4,0)C(4,4) + C(4,1)C(4,3) + C(4,3)C(4,1)
 *      + C(4,4)C(4,0)) / C(8,4) = (1+16+16+1)/70 = 34/70 = 17/35 = 0.485714285714... (see the
 *      per-case comment below for the by-hand table enumeration).
 *   2. Mann-Whitney U against a small no-tie example (n=3 vs n=3, fully separated: {1,2,3} vs
 *      {4,5,6}) with ranks and U computed by hand in the comment, and a tie-handling example
 *      ({1,2,2} vs {2,3,4}) with average-rank computation shown by hand.
 *   3. Because the two-sided normal-approximation p-value itself isn't elementary to verify
 *      fully by hand, StatsUtil.zToP is cross-checked here against an independently-coded normal
 *      CDF approximation (Abramowitz &amp; Stegun 26.2.17) — a second, differently-derived
 *      implementation agreeing with StatsUtil's erfc-based one is a much stronger check than
 *      trusting a single implementation, matching AnovaUtilTest's existing "second independent
 *      formula" strategy in this repo.
 *   4. A full-pipeline integration case (EnrichmentAnalyzer.run) with a hand-computed categorical
 *      (Fisher) result against a background set that includes the in-loci genes (this project's
 *      documented background convention), plus edge cases: empty in-loci set, empty background,
 *      and a zero-variance numeric column, none of which may crash or produce NaN/garbage.
 *   5. (Applied outside this file, per the task's fail-then-pass requirement:) this test is run
 *      once against the real EnrichmentAnalyzer, then again after deliberately breaking a
 *      computation, to confirm it actually fails rather than vacuously passing.
 */
public class EnrichmentAnalyzerTest {

    public static void main(String[] args) {
        int failures = 0;

        // ── Fisher's exact test: classic tea-tasting 2x2, table [[3,1],[1,3]] ──
        // By-hand hypergeometric enumeration for margins row=(4,4), col=(4,4), n=8, C(8,4)=70:
        //   a=0: C(4,0)*C(4,4)=1*1=1     -> P=1/70=0.0142857  (<=P(3), included)
        //   a=1: C(4,1)*C(4,3)=4*4=16    -> P=16/70=0.2285714 (<=P(3), included)
        //   a=2: C(4,2)*C(4,2)=6*6=36    -> P=36/70=0.5142857 (> P(3), excluded)
        //   a=3: C(4,3)*C(4,1)=4*4=16    -> P=16/70=0.2285714 (observed table itself)
        //   a=4: C(4,4)*C(4,0)=1*1=1     -> P=1/70=0.0142857  (<=P(3), included)
        // two-sided p = (1+16+16+1)/70 = 34/70 = 17/35 = 0.4857142857142857
        {
            double p = EnrichmentAnalyzer.fisherExactTwoSided(3, 1, 1, 3);
            failures += checkClose("fisher tea-tasting [[3,1],[1,3]] p == 17/35", p, 17.0 / 35.0, 1e-9);
        }

        // ── Fisher's exact test: symmetric extreme case, [[4,0],[0,4]] -> maximally significant ──
        // Same margins (4,4)/(4,4) as the tea-tasting case above, C(8,4)=70, per-table
        // probabilities P(a=0..4) = 1,16,36,16,1 (over 70). Observed a=4 has P=1/70. The
        // two-sided definition sums every table with P <= observed P: that's a=4 (1/70, the
        // observed table itself) AND the symmetric opposite extreme a=0 (also 1/70) — a=1,2,3 all
        // have strictly larger probabilities (16/70, 36/70, 16/70) and are excluded.
        // p = (1+1)/70 = 2/70 = 1/35 = 0.028571428571428571
        {
            double p = EnrichmentAnalyzer.fisherExactTwoSided(4, 0, 0, 4);
            failures += checkClose("fisher perfect-separation [[4,0],[0,4]] p == 1/35", p, 1.0 / 35.0, 1e-9);
        }

        // ── Fisher's exact test: degenerate margins (a whole row/col is zero) must not crash ──
        {
            double p = EnrichmentAnalyzer.fisherExactTwoSided(0, 0, 5, 5);
            failures += checkClose("fisher degenerate zero-row margin returns p=1", p, 1.0, 1e-12);
        }

        // ── Mann-Whitney: no ties, {1,2,3} vs {4,5,6} ──
        // Pooled ranks: 1->1, 2->2, 3->3, 4->4, 5->5, 6->6 (no ties).
        // R_a = 1+2+3 = 6.  U_a = R_a - n1(n1+1)/2 = 6 - 3 = 0 (a is uniformly the smaller group).
        // meanU = n1*n2/2 = 4.5.  varU (no ties) = n1*n2*(N+1)/12 = 9*7/12 = 5.25.
        // z = (0 - 4.5) / sqrt(5.25) = -4.5 / 2.291287847... = -1.963961012...
        {
            double[] a = {1, 2, 3}, b = {4, 5, 6};
            EnrichmentAnalyzer.MannWhitneyResult r = EnrichmentAnalyzer.mannWhitneyU(a, b);
            failures += checkClose("MWU no-ties U_a == 0", r.u, 0.0, 1e-9);
            double z = (0.0 - 4.5) / Math.sqrt(5.25);
            double independentP = twoSidedNormalP(z);
            failures += checkClose("MWU no-ties p matches independent normal-CDF cross-check", r.pValue, independentP, 1e-6);
        }

        // ── Mann-Whitney: tie handling, {1,2,2} vs {2,3,4} ──
        // Pooled values sorted: 1,2,2,2,3,4 (three 2's: one from a, two from b).
        // Ranks: 1->1; the three 2's occupy positions 2,3,4 -> average rank (2+3+4)/3=3 each; 3->5; 4->6.
        // Group a = {1,2,2} -> ranks {1,3,3} -> R_a = 7.  U_a = 7 - 3*4/2 = 7-6 = 1.
        // meanU = 3*3/2 = 4.5.  tie term: one tied block of size 3 -> (3^3-3)=24.
        // varU = (3*3/12) * (7 - 24/(6*5)) = 0.75 * (7-0.8) = 0.75*6.2 = 4.65
        // z = (1-4.5)/sqrt(4.65) = -3.5/2.156385... = -1.623242...
        {
            double[] a = {1, 2, 2}, b = {2, 3, 4};
            EnrichmentAnalyzer.MannWhitneyResult r = EnrichmentAnalyzer.mannWhitneyU(a, b);
            failures += checkClose("MWU tie-handling U_a == 1", r.u, 1.0, 1e-9);
            double z = (1.0 - 4.5) / Math.sqrt(4.65);
            double independentP = twoSidedNormalP(z);
            failures += checkClose("MWU tie-handling p matches independent normal-CDF cross-check", r.pValue, independentP, 1e-6);
        }

        // ── Mann-Whitney: zero-variance (every value identical) must not produce NaN ──
        {
            double[] a = {5, 5, 5}, b = {5, 5, 5, 5};
            EnrichmentAnalyzer.MannWhitneyResult r = EnrichmentAnalyzer.mannWhitneyU(a, b);
            failures += check("MWU all-tied p is not NaN", !Double.isNaN(r.pValue));
            failures += checkClose("MWU all-tied p == 1 (no distinguishing information)", r.pValue, 1.0, 1e-12);
        }

        // ── Full-pipeline integration: EnrichmentAnalyzer.run with categorical + numeric columns ──
        // In-loci genes: G1..G4.  Background (per this analyzer's documented convention) = G1..G8,
        // i.e. includes the in-loci genes themselves.
        //   risk column: G1=yes,G2=yes,G3=no,G4=yes,G5=no,G6=no,G7=no,G8=yes
        //   in-loci counts: yes=3, no=1 -> category of interest = "yes" (most frequent in-loci)
        //   in-loci table row: a=3 (yes), b=1 (no)
        //   background counts: yes=4 (G1,G2,G4,G8), no=4 (G3,G5,G6,G7) -> c=4, d=4
        //   2x2 = [[3,1],[4,4]], row=(4,8), col=(7,5), n=12, C(12,7)=792
        //   By-hand hypergeometric enumeration (a=in-loci&yes, c=7-a):
        //     a=0: C(4,0)C(8,7)=1*8=8         a=1: C(4,1)C(8,6)=4*28=112
        //     a=2: C(4,2)C(8,5)=6*56=336      a=3: C(4,3)C(8,4)=4*70=280 (observed)
        //     a=4: C(4,4)C(8,3)=1*56=56
        //   sum=792 (checks out against C(12,7)).  Tables with P<=P(3)=280/792: a=0,1,3,4
        //   (a=2's 336/792 is strictly greater, excluded).  p = (8+112+280+56)/792 = 456/792 = 19/33
        {
            Set<String> inLoci = new LinkedHashSet<>(Arrays.asList("G1", "G2", "G3", "G4"));
            Set<String> background = new LinkedHashSet<>(Arrays.asList(
                "G1", "G2", "G3", "G4", "G5", "G6", "G7", "G8"));
            Map<String, Map<String, String>> evidence = new LinkedHashMap<>();
            String[][] rows = {
                {"G1", "yes", "5"}, {"G2", "yes", "6"}, {"G3", "no", "1"}, {"G4", "yes", "7"},
                {"G5", "no", "2"}, {"G6", "no", "1.5"}, {"G7", "no", "0.5"}, {"G8", "yes", "8"}
            };
            for (String[] row : rows) {
                Map<String, String> m = new LinkedHashMap<>();
                m.put("risk", row[1]);
                m.put("expr", row[2]);
                evidence.put(row[0], m);
            }
            Map<String, String> columnTypes = new LinkedHashMap<>();
            columnTypes.put("risk", "categorical");
            columnTypes.put("expr", "numeric");

            EnrichmentAnalyzer.Result result = EnrichmentAnalyzer.run(inLoci, background, evidence, columnTypes);
            failures += check("integration: 2 columns returned", result.columns.size() == 2);
            EnrichmentAnalyzer.ColumnResult risk = findColumn(result, "risk");
            EnrichmentAnalyzer.ColumnResult expr = findColumn(result, "expr");
            failures += check("integration: risk column found", risk != null);
            failures += check("integration: expr column found", expr != null);
            if (risk != null) {
                failures += check("integration: risk test type is fisher_exact", "fisher_exact".equals(risk.testType));
                failures += check("integration: risk category_of_interest == 'yes'", "yes".equals(risk.categoryOfInterest));
                failures += check("integration: risk n_in_loci == 4", risk.nInLoci == 4);
                failures += check("integration: risk n_background == 8", risk.nBackground == 8);
                failures += checkClose("integration: risk p matches by-hand Fisher computation (19/33)",
                    risk.pValue, 19.0 / 33.0, 1e-9);
            }
            if (expr != null) {
                failures += check("integration: expr test type is mann_whitney_u", "mann_whitney_u".equals(expr.testType));
                failures += check("integration: expr p is a valid probability, not NaN", !Double.isNaN(expr.pValue) && expr.pValue >= 0 && expr.pValue <= 1);
                failures += check("integration: expr n_in_loci == 4", expr.nInLoci == 4);
                failures += check("integration: expr n_background == 8", expr.nBackground == 8);
            }
        }

        // ── Edge case: empty in-loci set must not crash or produce garbage ──
        {
            Set<String> inLoci = Collections.emptySet();
            Set<String> background = new LinkedHashSet<>(Arrays.asList("G1", "G2"));
            Map<String, Map<String, String>> evidence = new LinkedHashMap<>();
            Map<String, String> m1 = new LinkedHashMap<>(); m1.put("expr", "1"); evidence.put("G1", m1);
            Map<String, String> m2 = new LinkedHashMap<>(); m2.put("expr", "2"); evidence.put("G2", m2);
            Map<String, String> columnTypes = new LinkedHashMap<>();
            columnTypes.put("expr", "numeric");
            EnrichmentAnalyzer.Result result = EnrichmentAnalyzer.run(inLoci, background, evidence, columnTypes);
            EnrichmentAnalyzer.ColumnResult expr = findColumn(result, "expr");
            failures += check("edge case: empty in-loci does not crash", expr != null);
            failures += check("edge case: empty in-loci gives n_in_loci == 0", expr != null && expr.nInLoci == 0);
            failures += check("edge case: empty in-loci gives NaN p (not garbage 0/Infinity)", expr != null && Double.isNaN(expr.pValue));
        }

        // ── Edge case: empty background must not crash ──
        {
            Set<String> inLoci = new LinkedHashSet<>(Arrays.asList("G1"));
            Set<String> background = Collections.emptySet();
            Map<String, Map<String, String>> evidence = new LinkedHashMap<>();
            Map<String, String> m1 = new LinkedHashMap<>(); m1.put("expr", "1"); evidence.put("G1", m1);
            Map<String, String> columnTypes = new LinkedHashMap<>();
            columnTypes.put("expr", "numeric");
            EnrichmentAnalyzer.Result result = EnrichmentAnalyzer.run(inLoci, background, evidence, columnTypes);
            EnrichmentAnalyzer.ColumnResult expr = findColumn(result, "expr");
            failures += check("edge case: empty background does not crash", expr != null);
            failures += check("edge case: empty background gives n_background == 0", expr != null && expr.nBackground == 0);
        }

        // ── Edge case: evidence column with no variance at all (constant value) ──
        {
            Set<String> inLoci = new LinkedHashSet<>(Arrays.asList("G1", "G2"));
            Set<String> background = new LinkedHashSet<>(Arrays.asList("G1", "G2", "G3", "G4"));
            Map<String, Map<String, String>> evidence = new LinkedHashMap<>();
            for (String g : Arrays.asList("G1", "G2", "G3", "G4")) {
                Map<String, String> m = new LinkedHashMap<>(); m.put("expr", "3.5"); evidence.put(g, m);
            }
            Map<String, String> columnTypes = new LinkedHashMap<>();
            columnTypes.put("expr", "numeric");
            EnrichmentAnalyzer.Result result = EnrichmentAnalyzer.run(inLoci, background, evidence, columnTypes);
            EnrichmentAnalyzer.ColumnResult expr = findColumn(result, "expr");
            failures += check("edge case: zero-variance column p is not NaN", expr != null && !Double.isNaN(expr.pValue));
            failures += checkClose("edge case: zero-variance column p == 1", expr != null ? expr.pValue : -1, 1.0, 1e-9);
        }

        if (failures == 0) {
            System.out.println("PASS: all EnrichmentAnalyzer tests passed");
        } else {
            System.out.println("FAIL: " + failures + " test(s) failed");
            System.exit(1);
        }
    }

    private static EnrichmentAnalyzer.ColumnResult findColumn(EnrichmentAnalyzer.Result r, String name) {
        for (EnrichmentAnalyzer.ColumnResult c : r.columns) if (c.column.equals(name)) return c;
        return null;
    }

    /**
     * Independent two-sided normal-distribution p-value from a z-score, via the Abramowitz &amp;
     * Stegun 26.2.17 polynomial approximation of Phi(z) (max error ~7.5e-8) — coded from scratch
     * here, sharing no code with StatsUtil.erfc/zToP, specifically so it can serve as an
     * independent cross-check of StatsUtil's erfc-based implementation.
     */
    private static double twoSidedNormalP(double z) {
        double az = Math.abs(z);
        double t = 1.0 / (1.0 + 0.2316419 * az);
        double poly = t * (0.319381530
            + t * (-0.356563782
            + t * (1.781477937
            + t * (-1.821255978
            + t * 1.330274429))));
        double phiUpper = (1.0 / Math.sqrt(2 * Math.PI)) * Math.exp(-az * az / 2.0) * poly;
        // phiUpper approximates P(Z > az) = 1 - Phi(az); two-sided p = 2 * that upper tail.
        return 2.0 * phiUpper;
    }

    private static int checkClose(String label, double actual, double expected, double tol) {
        if (Double.isNaN(actual) || Math.abs(actual - expected) > Math.abs(tol) + 1e-15) {
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
