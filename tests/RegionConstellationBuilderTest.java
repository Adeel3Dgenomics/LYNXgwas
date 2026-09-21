import java.util.*;

/**
 * Regression test for RegionConstellationBuilder (DECISIONS_PHASE5.md section 3.2/3.6).
 *
 * Fixture: 4 datasets across 2 diseases, using the EXPLICIT disease-name field (not the id-prefix
 * convention — dataset ids here are "study1".."study4", deliberately without hyphens, so this test
 * also proves the Phase 5.1 explicit-disease wiring actually works, not just the id-prefix fallback
 * already covered by GeneConstellationBuilderTest):
 *   study1, study2 -> diseaseName "Disease A"
 *   study3, study4 -> diseaseName "Disease B"
 *
 * Two pooled regions (LocusRows):
 *   Region 1 (chr1:1000-2000, nearestGene="GENE_X"): all 4 datasets testable.
 *     study1: p=1e-9  (-log10p=9),  beta=0.5
 *     study2: p=1e-10 (-log10p=10), beta=0.6
 *     study3: p=1e-9  (-log10p=9),  beta=-0.4
 *     study4: p=1e-8  (-log10p=8),  beta=-0.5
 *     -> between-disease p-groups:      {DiseaseA:[9,10], DiseaseB:[9,8]}
 *     -> between-disease effect-groups: {DiseaseA:[0.5,0.6], DiseaseB:[-0.4,-0.5]}
 *   Region 2 (chr1:5000-6000, nearestGene="" — no assignable nearest gene): only study1 testable.
 *     study1: p=1e-9 (-log10p=9), beta=0.3
 *     -> only one disease present at all -> between-disease ANOVA must be null (not attempted on 1 group)
 *     -> that one disease has only 1 contributing dataset -> within-disease ANOVA for it must also be null
 *
 * The exact F/p-value correctness of AnovaUtil.oneWay itself is already covered by AnovaUtilTest; this
 * test verifies RegionConstellationBuilder assembles the *correct groups* by cross-checking its output
 * against a direct AnovaUtil.oneWay call on the same hand-derived expected groups, and verifies the
 * region-identity/no-gene-still-gets-a-node/null-degenerate-case plumbing around that.
 */
public class RegionConstellationBuilderTest {

    public static void main(String[] progArgs) throws Exception {
        int failures = 0;

        MultiLocusResult mlr = buildFixture();
        List<GeneConstellationResult.RegionEntry> regions = RegionConstellationBuilder.build(mlr, 5e-8);

        failures += check("exactly 2 regions produced (one per LocusRow, none dropped)", regions.size(), 2);

        GeneConstellationResult.RegionEntry r1 = findByRegionId(regions, "1:1000-2000");
        GeneConstellationResult.RegionEntry r2 = findByRegionId(regions, "1:5000-6000");
        failures += check("region 1 found", r1 != null, true);
        failures += check("region 2 found", r2 != null, true);

        // ── Region 1: basic aggregates ──
        failures += check("region1 nearestGene", r1.nearestGene, "GENE_X");
        failures += check("region1 n_datasets_total == 4", r1.nDatasetsTotal, 4);
        failures += check("region1 n_datasets_significant == 4 (all <= 5e-8)", r1.nDatasetsSignificant, 4);
        failures += check("region1 agg_neglog10_p == 10 (min p = 1e-10)", r1.aggNegLog10P, 10.0);
        // direction: study1(+),study2(+),study3(-),study4(-) -> 2 of 4 risk
        failures += check("region1 direction_fraction_risk == 0.5", r1.directionFractionRisk, 0.5);

        // ── Region 1: between-disease ANOVA on -log10(p), cross-checked against a direct AnovaUtil call ──
        AnovaUtil.Result expectedP = AnovaUtil.oneWay(Arrays.asList(
            new double[]{9.0, 10.0}, new double[]{9.0, 8.0}));
        failures += check("region1 betweenDiseaseAnova present", r1.betweenDiseaseAnova != null, true);
        failures += check("region1 betweenDiseaseAnova.fStat matches direct AnovaUtil call",
            r1.betweenDiseaseAnova.fStat, expectedP.fStat);
        failures += check("region1 betweenDiseaseAnova.pValue matches direct AnovaUtil call",
            r1.betweenDiseaseAnova.pValue, expectedP.pValue);

        // ── Region 1: between-disease ANOVA on ln(effect) (beta here, so no ln() applied) ──
        AnovaUtil.Result expectedEffect = AnovaUtil.oneWay(Arrays.asList(
            new double[]{0.5, 0.6}, new double[]{-0.4, -0.5}));
        failures += check("region1 betweenDiseaseAnovaEffect present", r1.betweenDiseaseAnovaEffect != null, true);
        failures += check("region1 betweenDiseaseAnovaEffect.fStat matches direct AnovaUtil call",
            r1.betweenDiseaseAnovaEffect.fStat, expectedEffect.fStat);

        // ── Within-disease ANOVA is never attempted at region granularity (see the class javadoc):
        // a single region gives each dataset exactly one observation, so "group by dataset within a
        // disease" always has zero residual degrees of freedom by construction, not just here. ──
        failures += check("region1 withinDiseaseAnova is left empty (never attempted at region level)",
            r1.withinDiseaseAnova.isEmpty(), true);
        failures += check("region1 withinDiseaseAnovaEffect is left empty (never attempted at region level)",
            r1.withinDiseaseAnovaEffect.isEmpty(), true);

        // ── Region 2: no nearest gene still produces a node (unlike the gene-level view) ──
        failures += check("region2 nearestGene is empty (no gene assigned)", r2.nearestGene, "");
        failures += check("region2 n_datasets_total == 1", r2.nDatasetsTotal, 1);

        // ── Region 2: degenerate cases must stay null/empty, not throw or fabricate a result ──
        failures += check("region2 betweenDiseaseAnova is null (only 1 disease present)",
            r2.betweenDiseaseAnova == null, true);
        failures += check("region2 withinDiseaseAnova is empty (never attempted at region level)",
            r2.withinDiseaseAnova.isEmpty(), true);

        if (failures == 0) {
            System.out.println("PASS: all RegionConstellationBuilder tests passed");
        } else {
            System.out.println("FAIL: " + failures + " RegionConstellationBuilder test(s) failed");
            System.exit(1);
        }
    }

    private static MultiLocusResult buildFixture() {
        MultiLocusResult mlr = new MultiLocusResult();
        mlr.datasets = new ArrayList<>();
        mlr.datasets.add(dataset("study1", "Disease A"));
        mlr.datasets.add(dataset("study2", "Disease A"));
        mlr.datasets.add(dataset("study3", "Disease B"));
        mlr.datasets.add(dataset("study4", "Disease B"));

        MultiLocusResult.LocusRow row1 = new MultiLocusResult.LocusRow();
        row1.index = 0; row1.chr = "1"; row1.start = 1000; row1.end = 2000; row1.nearestGene = "GENE_X";
        row1.cells = new LinkedHashMap<>();
        row1.cells.put("study1", stat(1e-9, 0.5));
        row1.cells.put("study2", stat(1e-10, 0.6));
        row1.cells.put("study3", stat(1e-9, -0.4));
        row1.cells.put("study4", stat(1e-8, -0.5));
        row1.overallBestP = 1e-10;

        MultiLocusResult.LocusRow row2 = new MultiLocusResult.LocusRow();
        row2.index = 1; row2.chr = "1"; row2.start = 5000; row2.end = 6000; row2.nearestGene = "";
        row2.cells = new LinkedHashMap<>();
        row2.cells.put("study1", stat(1e-9, 0.3));
        row2.overallBestP = 1e-9;

        mlr.loci = Arrays.asList(row1, row2);
        return mlr;
    }

    private static MultiLocusResult.DatasetInfo dataset(String id, String diseaseName) {
        MultiLocusResult.DatasetInfo d = new MultiLocusResult.DatasetInfo(id, id, "Beta");
        d.diseaseName = diseaseName;
        return d;
    }

    private static MultiLocusResult.DatasetLocusStat stat(double p, double beta) {
        MultiLocusResult.DatasetLocusStat s = new MultiLocusResult.DatasetLocusStat();
        s.bestP = p;
        s.beta = beta;
        return s;
    }

    private static GeneConstellationResult.RegionEntry findByRegionId(
            List<GeneConstellationResult.RegionEntry> regions, String id) {
        for (GeneConstellationResult.RegionEntry r : regions) if (r.regionId.equals(id)) return r;
        return null;
    }

    private static int check(String label, Object actual, Object expected) {
        boolean ok = Objects.equals(actual, expected)
            || (actual instanceof Double && expected instanceof Double
                && Math.abs((Double) actual - (Double) expected) < 1e-9);
        if (!ok) {
            System.out.println("FAIL: " + label + " (expected " + expected + ", got " + actual + ")");
            return 1;
        }
        System.out.println("PASS: " + label);
        return 0;
    }

    private static int check(String label, boolean cond) {
        if (!cond) { System.out.println("FAIL: " + label); return 1; }
        System.out.println("PASS: " + label);
        return 0;
    }
}
