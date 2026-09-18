import java.util.*;

/**
 * Derives a gene-level "constellation" view from an already-computed {@link MultiLocusResult}
 * (the Locus Matrix job output): one row per nearest-gene instead of one row per locus, with
 * cross-dataset replication/significance/direction summaries plus one-way ANOVA (via the
 * already-tested {@code AnovaUtil.oneWay}) both across the six disease groups and within each
 * disease group's own datasets.
 *
 * This is a pure, read-only derived view — it does not re-run any expensive pipeline, it just
 * re-groups the rows of a completed MultiLocusResult by gene.
 */
public class GeneConstellationBuilder {

    /** Default genome-wide significance threshold used to count a dataset "significant" for a gene. */
    public static final double DEFAULT_THRESHOLD = 5e-8;

    private GeneConstellationBuilder() {}

    public static GeneConstellationResult build(MultiLocusResult mlr, double genomewideThreshold) {
        GeneConstellationResult out = new GeneConstellationResult();
        out.name = mlr.name;
        out.refPanelLabel = mlr.refPanelLabel;
        out.createdAt = mlr.createdAt;
        out.threshold = genomewideThreshold;
        out.datasets = mlr.datasets;

        // dataset id -> disease group (best-effort: substring before first '-', else "ungrouped")
        Map<String, String> diseaseOf = new LinkedHashMap<>();
        for (MultiLocusResult.DatasetInfo d : mlr.datasets) {
            diseaseOf.put(d.id, diseaseGroupOf(d.id));
        }
        Map<String, String> nameOf = new LinkedHashMap<>();
        for (MultiLocusResult.DatasetInfo d : mlr.datasets) nameOf.put(d.id, d.name);

        LinkedHashMap<String, GeneAgg> genes = new LinkedHashMap<>();
        int unassigned = 0;

        for (MultiLocusResult.LocusRow row : mlr.loci) {
            String gene = row.nearestGene == null ? "" : row.nearestGene.trim();
            if (gene.isEmpty()) { unassigned++; continue; }

            GeneAgg agg = genes.computeIfAbsent(gene, g -> new GeneAgg(g));
            if (Double.isNaN(agg.repP) || (!Double.isNaN(row.overallBestP) && row.overallBestP < agg.repP)) {
                agg.repP = row.overallBestP;
                agg.chr = row.chr;
                agg.pos = row.start;
            }

            for (Map.Entry<String, MultiLocusResult.DatasetLocusStat> e : row.cells.entrySet()) {
                String dsId = e.getKey();
                MultiLocusResult.DatasetLocusStat s = e.getValue();
                if (s == null || Double.isNaN(s.bestP)) continue; // not testable in this dataset for this locus

                GeneAgg.Obs obs = new GeneAgg.Obs();
                obs.datasetId = dsId;
                obs.datasetName = nameOf.getOrDefault(dsId, dsId);
                obs.disease = diseaseOf.getOrDefault(dsId, "ungrouped");
                obs.bestP = s.bestP;
                obs.beta = s.beta;
                obs.or = s.or;
                agg.observations.add(obs);
            }
        }

        for (GeneAgg agg : genes.values()) {
            GeneConstellationResult.GeneEntry ge = new GeneConstellationResult.GeneEntry();
            ge.gene = agg.gene;
            ge.chr = agg.chr;
            ge.pos = agg.pos;

            double minP = Double.NaN;
            int nSig = 0, nTotal = 0;
            int riskCount = 0, dirTotal = 0;
            for (GeneAgg.Obs o : agg.observations) {
                nTotal++;
                if (o.bestP <= genomewideThreshold) nSig++;
                if (Double.isNaN(minP) || o.bestP < minP) minP = o.bestP;

                Boolean risk = null;
                if (!Double.isNaN(o.beta)) risk = o.beta > 0;
                else if (!Double.isNaN(o.or)) risk = o.or > 1;
                if (risk != null) { dirTotal++; if (risk) riskCount++; }

                GeneConstellationResult.PerDataset pd = new GeneConstellationResult.PerDataset();
                pd.id = o.datasetId; pd.name = o.datasetName; pd.disease = o.disease;
                pd.bestP = o.bestP; pd.beta = o.beta; pd.or = o.or;
                ge.perDataset.add(pd);
            }
            ge.nDatasetsSignificant = nSig;
            ge.nDatasetsTotal = nTotal;
            ge.aggNegLog10P = (Double.isNaN(minP) || minP <= 0) ? Double.NaN : -Math.log10(minP);
            ge.directionFractionRisk = dirTotal > 0 ? (double) riskCount / dirTotal : Double.NaN;

            // ── Between-disease ANOVA: group this gene's -log10(p) observations by disease ──
            Map<String, List<Double>> byDisease = new LinkedHashMap<>();
            for (GeneAgg.Obs o : agg.observations) {
                double v = neglog10(o.bestP);
                if (Double.isNaN(v)) continue;
                byDisease.computeIfAbsent(o.disease, k -> new ArrayList<>()).add(v);
            }
            if (byDisease.size() >= 2) {
                List<double[]> groups = new ArrayList<>();
                for (List<Double> vals : byDisease.values()) groups.add(toArr(vals));
                try {
                    AnovaUtil.Result r = AnovaUtil.oneWay(groups);
                    GeneConstellationResult.AnovaSummary a = new GeneConstellationResult.AnovaSummary();
                    a.fStat = r.fStat; a.pValue = r.pValue; a.dfBetween = r.dfBetween; a.dfWithin = r.dfWithin;
                    ge.betweenDiseaseAnova = a;
                } catch (IllegalArgumentException skip) {
                    // too few non-empty groups / degenerate n<=k — leave null, don't fail the whole build
                }
            }

            // ── Within-disease ANOVA: for each disease, group by dataset id ──
            Map<String, Map<String, List<Double>>> byDiseaseThenDataset = new LinkedHashMap<>();
            for (GeneAgg.Obs o : agg.observations) {
                double v = neglog10(o.bestP);
                if (Double.isNaN(v)) continue;
                byDiseaseThenDataset
                    .computeIfAbsent(o.disease, k -> new LinkedHashMap<>())
                    .computeIfAbsent(o.datasetId, k -> new ArrayList<>())
                    .add(v);
            }
            for (Map.Entry<String, Map<String, List<Double>>> de : byDiseaseThenDataset.entrySet()) {
                String disease = de.getKey();
                Map<String, List<Double>> byDataset = de.getValue();
                if (byDataset.size() < 2) {
                    ge.withinDiseaseAnova.put(disease, null); // only one dataset in this disease group for this gene
                    continue;
                }
                List<double[]> groups = new ArrayList<>();
                for (List<Double> vals : byDataset.values()) groups.add(toArr(vals));
                try {
                    AnovaUtil.Result r = AnovaUtil.oneWay(groups);
                    GeneConstellationResult.AnovaSummary a = new GeneConstellationResult.AnovaSummary();
                    a.fStat = r.fStat; a.pValue = r.pValue; a.dfBetween = r.dfBetween; a.dfWithin = r.dfWithin;
                    ge.withinDiseaseAnova.put(disease, a);
                } catch (IllegalArgumentException skip) {
                    ge.withinDiseaseAnova.put(disease, null);
                }
            }

            out.genes.add(ge);
        }

        out.unassignedLociCount = unassigned;
        return out;
    }

    /** Best-effort disease-group key: substring of the project id before its first '-'.
     *  Falls back cleanly to "ungrouped" for ids that don't follow the convention. */
    static String diseaseGroupOf(String projectId) {
        if (projectId == null || projectId.isEmpty()) return "ungrouped";
        int dash = projectId.indexOf('-');
        if (dash <= 0) return "ungrouped";
        return projectId.substring(0, dash);
    }

    private static double neglog10(double p) {
        if (Double.isNaN(p) || p <= 0) return Double.NaN;
        return -Math.log10(p);
    }

    private static double[] toArr(List<Double> vals) {
        double[] a = new double[vals.size()];
        for (int i = 0; i < a.length; i++) a[i] = vals.get(i);
        return a;
    }

    /** Internal accumulator while grouping locus rows by gene. */
    private static class GeneAgg {
        final String gene;
        String chr = "";
        long pos = 0;
        double repP = Double.NaN; // smallest overallBestP seen so far, drives chr/pos pick
        List<Obs> observations = new ArrayList<>();

        GeneAgg(String gene) { this.gene = gene; }

        static class Obs {
            String datasetId, datasetName, disease;
            double bestP = Double.NaN, beta = Double.NaN, or = Double.NaN;
        }
    }
}
