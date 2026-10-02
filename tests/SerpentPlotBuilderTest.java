import java.util.*;

/**
 * Tests for SerpentPlotBuilder (Locus Serpent Plot view + category analysis).
 *
 * Statistics are checked against independent computations, not just "it ran":
 *   1. Tier boundaries at exactly 5e-3 / 5e-5 / 5e-8.
 *   2. PERMANOVA on 1-D points with Euclidean distance must equal the classical centroid-based
 *      sums of squares (Anderson 2001 shows the two coincide for Euclidean distances), computed
 *      here from the coordinates, never from the distance matrix.
 *   3. Permutation p-value for a perfectly separated 4+4 design: only 2 of the C(8,4)=70 labelings
 *      reach the observed F, so the exact p is 2/70; the Monte-Carlo p must land near it.
 *   4. Adjusted Rand index: identical and relabeled partitions give 1, and a hand-computed 2x3
 *      contingency example gives 0.8/3.3.
 *   5. UPGMA on a 4-point matrix: merge heights equal the hand-computed averages.
 *   6. End-to-end on a synthetic corpus where disease drives the loci and ancestry does not:
 *      disease must explain far more of the structure than ancestry, the clusters must match
 *      disease exactly, the consistent loci must be exactly each disease's own loci, and the JSON
 *      must parse.
 */
public class SerpentPlotBuilderTest {

    static int failures = 0;

    public static void main(String[] args) {
        // 1. tiers
        check("tier NaN", SerpentPlotBuilder.tierOf(Double.NaN) == 0);
        check("tier 5e-3 is grey", SerpentPlotBuilder.tierOf(5e-3) == 0);
        check("tier 4.9e-3 is green", SerpentPlotBuilder.tierOf(4.9e-3) == 1);
        check("tier 5e-5 is green", SerpentPlotBuilder.tierOf(5e-5) == 1);
        check("tier 4.9e-5 is yellow", SerpentPlotBuilder.tierOf(4.9e-5) == 2);
        check("tier 5e-8 is yellow", SerpentPlotBuilder.tierOf(5e-8) == 2);
        check("tier 4.9e-8 is red", SerpentPlotBuilder.tierOf(4.9e-8) == 3);

        // 2. PERMANOVA == centroid sums of squares for Euclidean 1-D points
        {
            double[] x = {1, 2, 3, 10, 11, 12.5, 5, 6};
            int[] g = {0, 0, 0, 1, 1, 1, 2, 2};
            double[][] d = new double[x.length][x.length];
            for (int i = 0; i < x.length; i++) for (int j = 0; j < x.length; j++) d[i][j] = Math.abs(x[i] - x[j]);
            double mean = Arrays.stream(x).average().getAsDouble();
            double ssT = 0; for (double v : x) ssT += (v - mean) * (v - mean);
            double[] gs = new double[3]; int[] gc = new int[3];
            for (int i = 0; i < x.length; i++) { gs[g[i]] += x[i]; gc[g[i]]++; }
            double ssW = 0; for (int i = 0; i < x.length; i++) { double m = gs[g[i]] / gc[g[i]]; ssW += (x[i] - m) * (x[i] - m); }
            double F = ((ssT - ssW) / 2) / (ssW / (x.length - 3));
            SerpentPlotBuilder.Permanova pm = SerpentPlotBuilder.permanova(d, g, 0, 1);
            close("permanova ssTotal", pm.ssTotal, ssT, 1e-9);
            close("permanova ssWithin", pm.ssWithin, ssW, 1e-9);
            close("permanova R2", pm.r2, (ssT - ssW) / ssT, 1e-12);
            close("permanova pseudo-F", pm.pseudoF, F, 1e-9);
        }

        // 3. permutation p for a perfectly separated 4+4 design ~ 2/70
        {
            double[] x = {0, 0.1, 0.2, 0.3, 5, 5.1, 5.2, 5.3};
            int[] g = {0, 0, 0, 0, 1, 1, 1, 1};
            double[][] d = new double[8][8];
            for (int i = 0; i < 8; i++) for (int j = 0; j < 8; j++) d[i][j] = Math.abs(x[i] - x[j]);
            SerpentPlotBuilder.Permanova pm = SerpentPlotBuilder.permanova(d, g, 4999, 7);
            close("permutation p near exact 2/70", pm.p, 2.0 / 70.0, 0.01);
            SerpentPlotBuilder.Permanova again = SerpentPlotBuilder.permanova(d, g, 4999, 7);
            close("permutation p reproducible with same seed", again.p, pm.p, 0);
        }

        // 4. adjusted Rand index
        close("ARI identical", SerpentPlotBuilder.adjustedRandIndex(new int[]{0,0,1,1,2,2}, new int[]{0,0,1,1,2,2}), 1.0, 1e-12);
        close("ARI relabeled", SerpentPlotBuilder.adjustedRandIndex(new int[]{0,0,1,1,2,2}, new int[]{2,2,0,0,1,1}), 1.0, 1e-12);
        close("ARI hand example", SerpentPlotBuilder.adjustedRandIndex(new int[]{0,0,0,1,1,1}, new int[]{0,0,1,1,2,2}), 0.8 / 3.3, 1e-12);

        // 5. UPGMA heights
        {
            double[][] d = {
                {0, 0.1, 0.9, 0.9},
                {0.1, 0, 0.7, 0.9},
                {0.9, 0.7, 0, 0.2},
                {0.9, 0.9, 0.2, 0}};
            List<SerpentPlotBuilder.Merge> m = SerpentPlotBuilder.upgma(d);
            check("UPGMA 3 merges", m.size() == 3);
            close("UPGMA first merge height", m.get(0).height, 0.1, 1e-12);
            close("UPGMA second merge height", m.get(1).height, 0.2, 1e-12);
            close("UPGMA root height = mean of 4 cross distances", m.get(2).height, (0.9 + 0.9 + 0.7 + 0.9) / 4, 1e-12);
            int[] cl = SerpentPlotBuilder.cutTree(m, 4, 2);
            check("cut k=2 keeps pairs", cl[0] == cl[1] && cl[2] == cl[3] && cl[0] != cl[2]);
            int[] order = SerpentPlotBuilder.leafOrder(m, 4);
            check("leaf order keeps pairs adjacent", Math.abs(pos(order, 0) - pos(order, 1)) == 1 && Math.abs(pos(order, 2) - pos(order, 3)) == 1);
        }

        // 5b. silhouette: two tight, well-separated pairs -> mean silhouette by hand
        {
            double[][] d = {
                {0, 0.1, 0.9, 0.9},
                {0.1, 0, 0.9, 0.9},
                {0.9, 0.9, 0, 0.2},
                {0.9, 0.9, 0.2, 0}};
            // points 0,1: a=0.1, b=0.9 -> 0.8/0.9 ; points 2,3: a=0.2, b=0.9 -> 0.7/0.9
            double want = (2 * (0.8 / 0.9) + 2 * (0.7 / 0.9)) / 4;
            close("mean silhouette hand example", SerpentPlotBuilder.meanSilhouette(d, new int[]{0, 0, 1, 1}, 2), want, 1e-12);
            close("singleton contributes 0", SerpentPlotBuilder.meanSilhouette(d, new int[]{0, 0, 1, 2}, 3), (2 * (0.8 / 0.9)) / 4, 1e-12);
        }

        // 6. end-to-end synthetic corpus
        {
            MultiLocusResult mlr = syntheticCorpus();
            SerpentPlotBuilder.Result r = SerpentPlotBuilder.build(mlr, 999);
            SerpentPlotBuilder.CategoryAnalysis dis = find(r, "Disease"), anc = find(r, "Ancestry");
            check("disease testable", dis != null && dis.testable);
            check("ancestry testable", anc != null && anc.testable);
            check("disease R2 > 0.8 (" + dis.r2 + ")", dis.r2 > 0.8);
            check("ancestry R2 < 0.3 (" + anc.r2 + ")", anc.r2 < 0.3);
            check("disease permutation p <= 0.01 (" + dis.permP + ")", dis.permP <= 0.01);
            check("ancestry permutation p > 0.05 (" + anc.permP + ")", anc.permP > 0.05);
            check("disease ARI == 1 with 3 clusters (" + dis.ari + ", k=" + r.nClusters + ")", r.nClusters == 3 && Math.abs(dis.ari - 1.0) < 1e-12);
            // consistent loci of disease "aaa" are exactly loci 0..4
            SerpentPlotBuilder.Level la = null;
            for (SerpentPlotBuilder.Level l : dis.levels) if (l.level.equals("aaa")) la = l;
            check("level aaa present", la != null);
            if (la != null) {
                Set<Integer> got = new TreeSet<>();
                for (SerpentPlotBuilder.ConsistentLocus c : la.consistentLoci) got.add(c.locusIndex);
                check("aaa consistent loci = {0..4} (" + got + ")", got.equals(new TreeSet<>(Arrays.asList(0, 1, 2, 3, 4))));
                SerpentPlotBuilder.ConsistentLocus c0 = la.consistentLoci.get(0);
                close("consistent-locus Fisher p matches EnrichmentAnalyzer",
                    c0.fisherP, EnrichmentAnalyzer.fisherExactTwoSided(4, 0, 0, 8), 1e-12);
            }
            // lead-SNP offset is relative to the overall best SNP
            SerpentPlotBuilder.Locus L0 = r.loci.get(0);
            check("offset of best dataset is 0", L0.cells[0].offsetBp == 0);
            check("offset of dataset 1 is +5000 (" + L0.cells[1].offsetBp + ")", L0.cells[1].offsetBp == 5000);
            // sample-size category present (>= 6 datasets with N)
            check("sample-size control category present", r.categoryNames.contains("Sample size"));
            // custom category carried through
            check("custom category present", r.categoryNames.contains("tissue"));
            String json = SerpentPlotBuilder.toJson(r);
            Object parsed = MiniJson.parse(json);
            Map<String, Object> obj = MiniJson.asObject(parsed);
            check("JSON parses with loci/analyses", obj != null && obj.containsKey("loci") && obj.containsKey("analyses") && obj.containsKey("clustering"));
        }

        // 7. outlier datasets must not use up the clusters, and a locus significant in every dataset
        //    must not be reported as group-specific
        {
            MultiLocusResult mlr = syntheticCorpus();
            for (int o = 0; o < 2; o++) { // two datasets with no genome-wide hits at all
                MultiLocusResult.DatasetInfo di = new MultiLocusResult.DatasetInfo("zzz-out" + o, "outlier " + o, "OR");
                di.diseaseName = "zzz"; mlr.datasets.add(di);
                for (MultiLocusResult.LocusRow row : mlr.loci) {
                    MultiLocusResult.DatasetLocusStat s = new MultiLocusResult.DatasetLocusStat();
                    s.bestP = (row.index + o) % 7 == 0 ? 1e-3 : 0.5; s.bestPos = 1_100_000; row.cells.put(di.id, s);
                }
            }
            MultiLocusResult.LocusRow every = new MultiLocusResult.LocusRow();
            every.index = 99; every.chr = "6"; every.start = 30_000_000; every.end = 31_000_000; every.nearestGene = "EVERYWHERE";
            for (MultiLocusResult.DatasetInfo di : mlr.datasets) {
                MultiLocusResult.DatasetLocusStat s = new MultiLocusResult.DatasetLocusStat();
                s.bestP = di.id.startsWith("zzz") ? 0.5 : 1e-20; s.bestPos = 30_500_000; every.cells.put(di.id, s);
            }
            every.overallBestDatasetId = mlr.datasets.get(0).id;
            mlr.loci.add(every);
            SerpentPlotBuilder.Result r = SerpentPlotBuilder.build(mlr, 199);
            boolean intact = true;
            for (int d = 0; d < 3; d++) for (int k = 1; k < 4; k++) if (r.clusterOf[d * 4 + k] != r.clusterOf[d * 4]) intact = false;
            boolean distinct = r.clusterOf[0] != r.clusterOf[4] && r.clusterOf[4] != r.clusterOf[8] && r.clusterOf[0] != r.clusterOf[8];
            check("outliers do not merge the three disease clusters (k=" + r.nClusters + ")", intact && distinct);
            check("datasets with < MIN_GWS_LOCI significant loci are left out of clustering",
                r.clusterOf[12] == -1 && r.clusterOf[13] == -1 && r.excluded.size() == 2);
            check("three clusters among the informative datasets (k=" + r.nClusters + ")", r.nClusters == 3);
            check("excluded datasets come last in the row order", r.leafOrder[12] >= 12 && r.leafOrder[13] >= 12);
            boolean everywhereListed = false;
            for (SerpentPlotBuilder.Level l : find(r, "Disease").levels)
                for (SerpentPlotBuilder.ConsistentLocus c : l.consistentLoci) if (c.locusIndex == 99) everywhereListed = true;
            check("locus significant in every dataset is not reported as group-specific", !everywhereListed);
        }

        // 8. locus detail (click on a serpent): per-dataset lead SNPs, OR/CI, dedup nearest genes
        {
            MultiLocusResult mlr = syntheticCorpus();
            MultiLocusResult.DatasetLocusStat s0 = mlr.loci.get(0).cells.get("aaa-ds0");
            s0.or = 1.5; s0.se = 0.1;                    // reported OR
            MultiLocusResult.DatasetLocusStat s1 = mlr.loci.get(0).cells.get("aaa-ds1");
            s1.beta = Math.log(2.0); s1.se = 0.2;        // OR must come from exp(beta)
            // fake gene annotation: GENE_A before 1,103,000, GENE_B after
            SerpentPlotBuilder.GeneLookup lookup = (chr, pos) -> pos < 1_103_000 ? new Object[]{"GENE_A", 0L} : new Object[]{"GENE_B", pos - 1_103_000};
            String js = SerpentPlotBuilder.locusDetailJson(mlr, 0, lookup);
            Map<String, Object> o = MiniJson.asObject(MiniJson.parse(js));
            List<Object> rows = MiniJson.asArray(o.get("rows"));
            check("detail has one row per dataset (" + rows.size() + ")", rows.size() == 12);
            Map<String, Object> r0 = MiniJson.asObject(rows.get(0));
            check("rows sorted by p (first is aaa-ds0)", "aaa-ds0".equals(r0.get("dataset_id")));
            close("reported OR kept", ((Number) r0.get("or")).doubleValue(), 1.5, 1e-12);
            close("CI lower = exp(ln1.5 - 1.96*0.1)", ((Number) r0.get("or_lo")).doubleValue(), Math.exp(Math.log(1.5) - 1.959964 * 0.1), 1e-9);
            Map<String, Object> r1 = null;
            for (Object x : rows) { Map<String, Object> m = MiniJson.asObject(x); if ("aaa-ds1".equals(m.get("dataset_id"))) r1 = m; }
            close("OR from beta = exp(beta)", ((Number) r1.get("or")).doubleValue(), 2.0, 1e-12);
            check("or_from_beta flagged", Boolean.TRUE.equals(r1.get("or_from_beta")));
            List<Object> genes = MiniJson.asArray(o.get("genes"));
            Set<String> names = new LinkedHashSet<>();
            for (Object g : genes) names.add(String.valueOf(MiniJson.asObject(g).get("gene")));
            check("genes deduplicated (" + names + ")", genes.size() == names.size() && names.equals(new LinkedHashSet<>(Arrays.asList("GENE_A", "GENE_B"))));
            check("unknown locus returns null", SerpentPlotBuilder.locusDetailJson(mlr, 999, lookup) == null);
        }

        if (failures == 0) System.out.println("PASS: all SerpentPlotBuilder tests passed");
        else { System.out.println("FAIL: " + failures + " SerpentPlotBuilder check(s) failed"); System.exit(1); }
    }

    /** 3 diseases x 4 datasets; each disease has 5 loci significant in all its datasets and nowhere
     *  else; ancestry alternates EUR/EAS across diseases (unrelated to the loci). */
    static MultiLocusResult syntheticCorpus() {
        MultiLocusResult mlr = new MultiLocusResult();
        mlr.name = "synthetic";
        String[] dis = {"aaa", "bbb", "ccc"};
        int[] ns = {1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000, 11000, 12000};
        int k = 0;
        for (String d : dis) for (int r = 0; r < 4; r++) {
            MultiLocusResult.DatasetInfo di = new MultiLocusResult.DatasetInfo(d + "-ds" + r, d + " " + r, "OR");
            di.diseaseName = d;
            di.ancestry = (r % 2 == 0) ? "EUR" : "EAS (Japan)";
            di.sampleN = ns[k++];
            di.categories.put("tissue", r < 2 ? "blood" : "brain");
            mlr.datasets.add(di);
        }
        Random rnd = new Random(3);
        for (int li = 0; li < 15; li++) {
            MultiLocusResult.LocusRow row = new MultiLocusResult.LocusRow();
            row.index = li; row.chr = String.valueOf(1 + li); row.start = 1_000_000; row.end = 1_200_000;
            row.nearestGene = "G" + li;
            int owner = li / 5;
            for (int di = 0; di < mlr.datasets.size(); di++) {
                MultiLocusResult.DatasetLocusStat s = new MultiLocusResult.DatasetLocusStat();
                boolean mine = di / 4 == owner;
                s.bestP = mine ? 1e-12 * (1 + di % 4) : 0.2 + 0.5 * rnd.nextDouble();
                s.bestPos = 1_100_000 + (di % 4) * 5000;
                s.bestSnpId = "rs" + li + "_" + di;
                row.cells.put(mlr.datasets.get(di).id, s);
            }
            row.overallBestDatasetId = mlr.datasets.get(owner * 4).id;
            mlr.loci.add(row);
        }
        return mlr;
    }

    static SerpentPlotBuilder.CategoryAnalysis find(SerpentPlotBuilder.Result r, String c) {
        for (SerpentPlotBuilder.CategoryAnalysis a : r.analyses) if (a.category.equals(c)) return a;
        return null;
    }

    static int pos(int[] a, int v) { for (int i = 0; i < a.length; i++) if (a[i] == v) return i; return -1; }

    static void check(String name, boolean ok) {
        if (ok) System.out.println("PASS: " + name);
        else { System.out.println("FAIL: " + name); failures++; }
    }

    static void close(String name, double got, double want, double tol) {
        boolean ok = !Double.isNaN(got) && Math.abs(got - want) <= tol;
        if (ok) System.out.println("PASS: " + name + " (" + got + ")");
        else { System.out.println("FAIL: " + name + " got " + got + " want " + want); failures++; }
    }
}
