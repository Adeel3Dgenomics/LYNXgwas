import java.util.*;

/**
 * Derives a region-level "constellation" view from an already-computed {@link MultiLocusResult}
 * (DECISIONS_PHASE5.md section 3.2) — one entry per already-pooled cross-dataset locus, i.e. one
 * entry per {@link MultiLocusResult.LocusRow}, *not* collapsed by nearest gene the way
 * {@link GeneConstellationBuilder} collapses them. Region pooling itself is not new work: it already
 * happens upstream (`MultiLocusScanner`/`MultiLocusMerger` build `LocusRow`s by genomic position
 * across every dataset) — this class only avoids throwing that granularity away, so a region with no
 * assignable nearest gene still gets a node instead of being dropped from `unassignedLociCount`.
 *
 * Reuses {@link GeneConstellationBuilder#effectValue} for the ln(effect) scale and
 * {@code AnovaUtil.oneWay} for the actual test — the same statistical machinery as gene-level
 * ANOVA, applied at a different, finer granularity, not a second implementation of either.
 *
 * <p><b>Only between-disease ANOVA is computed at this granularity — never within-disease.</b> This
 * was discovered while testing, not decided up front: gene-level within-disease ANOVA works because a
 * gene can span several pooled regions, giving one dataset several observations to group by; a single
 * region is, by definition, one {@code LocusRow}, so every dataset contributes at most one observation
 * to it. "Group by dataset within a disease" therefore always has exactly one observation per group —
 * zero residual degrees of freedom, every time, not just when data happens to be sparse. Computing it
 * anyway would mean calling {@code AnovaUtil.oneWay}, catching the resulting
 * {@code IllegalArgumentException}, and storing {@code null} on literally every region — dead work
 * with a guaranteed-null result, not a real analysis. {@code RegionEntry.withinDiseaseAnova}/
 * {@code withinDiseaseAnovaEffect} are left at their empty-map default here for this reason.</p>
 */
public class RegionConstellationBuilder {

    private RegionConstellationBuilder() {}

    public static List<GeneConstellationResult.RegionEntry> build(
            MultiLocusResult mlr, double genomewideThreshold) {

        Map<String, String> diseaseOf = new LinkedHashMap<>();
        Map<String, String> nameOf = new LinkedHashMap<>();
        for (MultiLocusResult.DatasetInfo d : mlr.datasets) {
            diseaseOf.put(d.id, GeneConstellationBuilder.diseaseGroupOf(d.id, d.diseaseName));
            nameOf.put(d.id, d.name);
        }

        List<GeneConstellationResult.RegionEntry> out = new ArrayList<>();

        for (MultiLocusResult.LocusRow row : mlr.loci) {
            GeneConstellationResult.RegionEntry re = new GeneConstellationResult.RegionEntry();
            re.regionId = row.chr + ":" + row.start + "-" + row.end;
            re.chr = row.chr;
            re.start = row.start;
            re.end = row.end;
            re.nearestGene = row.nearestGene == null ? "" : row.nearestGene.trim();

            double minP = Double.NaN;
            int nSig = 0, nTotal = 0;
            int riskCount = 0, dirTotal = 0;
            Map<String, List<Double>> byDiseaseP = new LinkedHashMap<>();
            Map<String, List<Double>> byDiseaseEffect = new LinkedHashMap<>();

            for (Map.Entry<String, MultiLocusResult.DatasetLocusStat> e : row.cells.entrySet()) {
                String dsId = e.getKey();
                MultiLocusResult.DatasetLocusStat s = e.getValue();
                if (s == null || Double.isNaN(s.bestP)) continue; // not testable in this dataset for this region

                String disease = diseaseOf.getOrDefault(dsId, "ungrouped");

                nTotal++;
                if (s.bestP <= genomewideThreshold) nSig++;
                if (Double.isNaN(minP) || s.bestP < minP) minP = s.bestP;

                Boolean risk = null;
                if (!Double.isNaN(s.beta)) risk = s.beta > 0;
                else if (!Double.isNaN(s.or)) risk = s.or > 1;
                if (risk != null) { dirTotal++; if (risk) riskCount++; }

                GeneConstellationResult.PerDataset pd = new GeneConstellationResult.PerDataset();
                pd.id = dsId; pd.name = nameOf.getOrDefault(dsId, dsId); pd.disease = disease;
                pd.bestP = s.bestP; pd.beta = s.beta; pd.or = s.or;
                re.perDataset.add(pd);

                // Only the between-disease grouping is collected here — see the class-level note on
                // withinDiseaseAnova below for why a within-disease ANOVA is never attempted at this
                // granularity at all, not just sometimes null.
                double pVal = neglog10(s.bestP);
                if (!Double.isNaN(pVal)) byDiseaseP.computeIfAbsent(disease, k -> new ArrayList<>()).add(pVal);
                double eVal = GeneConstellationBuilder.effectValue(s.beta, s.or);
                if (!Double.isNaN(eVal)) byDiseaseEffect.computeIfAbsent(disease, k -> new ArrayList<>()).add(eVal);
            }

            re.nDatasetsSignificant = nSig;
            re.nDatasetsTotal = nTotal;
            re.aggNegLog10P = (Double.isNaN(minP) || minP <= 0) ? Double.NaN : -Math.log10(minP);
            re.directionFractionRisk = dirTotal > 0 ? (double) riskCount / dirTotal : Double.NaN;

            re.betweenDiseaseAnova = betweenAnova(byDiseaseP);
            re.betweenDiseaseAnovaEffect = betweenAnova(byDiseaseEffect);
            // withinDiseaseAnova/withinDiseaseAnovaEffect are deliberately left as their empty-map
            // default (unlike GeneEntry's, which can be genuinely non-null) — see the class-level
            // javadoc: a single region gives each dataset exactly one observation, so a "within-disease,
            // grouped by dataset" comparison always has 0 residual degrees of freedom by construction,
            // not merely for lack of data. Attempting it and always getting IllegalArgumentException
            // back would be dead code exercised on every region, every time, for no result.

            out.add(re);
        }

        return out;
    }

    private static GeneConstellationResult.AnovaSummary betweenAnova(Map<String, List<Double>> byDisease) {
        if (byDisease.size() < 2) return null;
        List<double[]> groups = new ArrayList<>();
        for (List<Double> vals : byDisease.values()) groups.add(toArr(vals));
        try {
            AnovaUtil.Result r = AnovaUtil.oneWay(groups);
            GeneConstellationResult.AnovaSummary a = new GeneConstellationResult.AnovaSummary();
            a.fStat = r.fStat; a.pValue = r.pValue; a.dfBetween = r.dfBetween; a.dfWithin = r.dfWithin;
            return a;
        } catch (IllegalArgumentException skip) {
            return null;
        }
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
}
