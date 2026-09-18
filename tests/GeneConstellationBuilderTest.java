import java.util.*;

/**
 * Regression test for GeneConstellationBuilder's gene-pair "co-significance" edges: the count of
 * datasets in which two genes are *both* independently significant.
 *
 * Fixture: 3 datasets (d1,d2,d3), 3 genes (GENE_A, GENE_B, GENE_C), each with a hand-chosen best-p
 * per dataset (1e-9 = significant at the 5e-8 default threshold, 0.5/0.9 = not). Per-dataset
 * significant-gene sets, worked out by hand:
 *   d1: {GENE_A, GENE_B}          (GENE_C's p=0.9 in d1 -> not significant)
 *   d2: {GENE_A, GENE_B, GENE_C}  (all three significant)
 *   d3: {GENE_B, GENE_C}          (GENE_A's p=0.5 in d3 -> not significant)
 * Expected co-significance counts (a pair counts once per dataset where BOTH are in that
 * dataset's significant set):
 *   GENE_A-GENE_B: d1 yes, d2 yes, d3 no  -> 2
 *   GENE_A-GENE_C: d1 no,  d2 yes, d3 no  -> 1
 *   GENE_B-GENE_C: d1 no,  d2 yes, d3 yes -> 2
 */
public class GeneConstellationBuilderTest {

    public static void main(String[] args) {
        int failures = 0;

        MultiLocusResult mlr = buildFixture();

        // --- Default min-edge-count (2): only the two count=2 pairs should appear ---
        GeneConstellationResult r2 = GeneConstellationBuilder.build(mlr, 5e-8, 2, 300);
        failures += check("default threshold: exactly 2 edges reported (A-C's count=1 filtered out)",
            r2.edges.size(), 2);
        failures += checkHasEdge(r2, "GENE_A", "GENE_B", 2);
        failures += checkHasEdge(r2, "GENE_B", "GENE_C", 2);
        failures += check("GENE_A-GENE_C is NOT present at min_edge_count=2",
            !hasEdge(r2, "GENE_A", "GENE_C"), true);

        // --- min-edge-count=1: all three pairs, with exact hand-computed counts ---
        GeneConstellationResult r1 = GeneConstellationBuilder.build(mlr, 5e-8, 1, 300);
        failures += check("min_edge_count=1: all 3 pairs reported", r1.edges.size(), 3);
        failures += checkHasEdge(r1, "GENE_A", "GENE_B", 2);
        failures += checkHasEdge(r1, "GENE_A", "GENE_C", 1);
        failures += checkHasEdge(r1, "GENE_B", "GENE_C", 2);

        // --- edges are sorted highest-count-first ---
        failures += check("edges sorted descending by count", r1.edges.get(0).count >= r1.edges.get(1).count
            && r1.edges.get(1).count >= r1.edges.get(2).count, true);

        // --- maxEdges cap is respected ---
        GeneConstellationResult rCapped = GeneConstellationBuilder.build(mlr, 5e-8, 1, 1);
        failures += check("max_edges=1 caps the list to 1, keeping the highest count",
            rCapped.edges.size() == 1 && rCapped.edges.get(0).count == 2, true);

        // --- a gene alone in every dataset (no co-significance partner) produces zero edges ---
        MultiLocusResult solo = new MultiLocusResult();
        solo.datasets = Collections.singletonList(new MultiLocusResult.DatasetInfo("d1", "D1", "Beta"));
        MultiLocusResult.LocusRow row = new MultiLocusResult.LocusRow();
        row.index = 0; row.chr = "1"; row.start = 100; row.end = 200; row.nearestGene = "LONE_GENE";
        MultiLocusResult.DatasetLocusStat s = new MultiLocusResult.DatasetLocusStat();
        s.bestP = 1e-9;
        row.cells = new LinkedHashMap<>();
        row.cells.put("d1", s);
        row.overallBestP = 1e-9;
        solo.loci = Collections.singletonList(row);
        GeneConstellationResult soloResult = GeneConstellationBuilder.build(solo, 5e-8, 1, 300);
        failures += check("a single gene with no co-significant partner produces zero edges",
            soloResult.edges.size(), 0);

        if (failures == 0) {
            System.out.println("PASS: all GeneConstellationBuilder edge tests passed");
        } else {
            System.out.println("FAIL: " + failures + " test(s) failed");
            System.exit(1);
        }
    }

    private static MultiLocusResult buildFixture() {
        MultiLocusResult mlr = new MultiLocusResult();
        mlr.name = "edge-test";
        mlr.datasets = Arrays.asList(
            new MultiLocusResult.DatasetInfo("d1", "Dataset 1", "Beta"),
            new MultiLocusResult.DatasetInfo("d2", "Dataset 2", "Beta"),
            new MultiLocusResult.DatasetInfo("d3", "Dataset 3", "Beta"));
        mlr.loci = Arrays.asList(
            makeLocus(0, "GENE_A", 1_000, new double[]{1e-9, 1e-9, 0.5}),
            makeLocus(1, "GENE_B", 3_000, new double[]{1e-9, 1e-9, 1e-9}),
            makeLocus(2, "GENE_C", 5_000, new double[]{0.9, 1e-9, 1e-9}));
        return mlr;
    }

    /** dsPs[0..2] = best-p in d1,d2,d3 respectively. */
    private static MultiLocusResult.LocusRow makeLocus(int index, String gene, long start, double[] dsPs) {
        MultiLocusResult.LocusRow row = new MultiLocusResult.LocusRow();
        row.index = index; row.chr = "1"; row.start = start; row.end = start + 500;
        row.nearestGene = gene;
        row.cells = new LinkedHashMap<>();
        String[] ids = {"d1", "d2", "d3"};
        double best = Double.POSITIVE_INFINITY;
        for (int i = 0; i < ids.length; i++) {
            MultiLocusResult.DatasetLocusStat s = new MultiLocusResult.DatasetLocusStat();
            s.bestP = dsPs[i];
            row.cells.put(ids[i], s);
            best = Math.min(best, dsPs[i]);
        }
        row.overallBestP = best;
        return row;
    }

    private static boolean hasEdge(GeneConstellationResult r, String a, String b) {
        for (GeneConstellationResult.GeneEdge e : r.edges) {
            if ((e.geneA.equals(a) && e.geneB.equals(b)) || (e.geneA.equals(b) && e.geneB.equals(a))) return true;
        }
        return false;
    }

    private static int checkHasEdge(GeneConstellationResult r, String a, String b, int expectedCount) {
        for (GeneConstellationResult.GeneEdge e : r.edges) {
            if ((e.geneA.equals(a) && e.geneB.equals(b)) || (e.geneA.equals(b) && e.geneB.equals(a))) {
                return check(a + "-" + b + " count", e.count, expectedCount);
            }
        }
        System.out.println("FAIL: edge " + a + "-" + b + " not found at all");
        return 1;
    }

    private static int check(String label, Object actual, Object expected) {
        if (!Objects.equals(actual, expected)) {
            System.out.println("FAIL: " + label + " — expected " + expected + ", got " + actual);
            return 1;
        }
        System.out.println("PASS: " + label + " (" + actual + ")");
        return 0;
    }
}
