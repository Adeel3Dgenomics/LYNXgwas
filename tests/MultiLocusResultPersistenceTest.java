import java.util.*;

/**
 * Saved cross-dataset runs (output/multi_locus/<job>/result.json) must reload exactly:
 *   1. toJson -> fromJson -> toJson is byte-identical (every field survives, including NaN/missing
 *      values, categories, ancestry, sample size and extra columns);
 *   2. awkward characters in free-text fields (quotes, backslashes, tabs, newlines) produce valid JSON
 *      and come back unchanged;
 *   3. a reloaded result gives the same Locus Serpent Plot as the original (same tiers, clusters, R^2).
 */
public class MultiLocusResultPersistenceTest {

    static int failures = 0;

    public static void main(String[] args) {
        MultiLocusResult r = fixture();
        String json1 = r.toJson();
        MultiLocusResult back;
        try {
            back = MultiLocusResult.fromJson(json1);
        } catch (Exception e) {
            System.out.println("FAIL: fromJson threw " + e);
            System.exit(1);
            return;
        }
        String json2 = back.toJson();
        check("round trip is byte-identical", json1.equals(json2));

        MultiLocusResult.DatasetInfo d0 = back.datasets.get(0);
        check("ancestry survives", "EUR (Spanish)".equals(d0.ancestry));
        check("sample size survives", d0.sampleN == 12345);
        check("categories survive", "blood".equals(d0.categories.get("tissue")));
        MultiLocusResult.DatasetLocusStat c = back.loci.get(0).cells.get("aaa-1");
        check("NaN stays NaN", Double.isNaN(c.beta) && Double.isNaN(c.maf));
        check("p-value exact", c.bestP == 3.25e-12);
        check("position exact", c.bestPos == 123_456_789L);
        check("awkward extra value unchanged", "a \"quoted\"\tvalue\\with\nnewline".equals(c.extra.get("note")));
        check("empty cell map entry kept", back.loci.get(1).cells.containsKey("bbb-1"));

        SerpentPlotBuilder.Result s1 = SerpentPlotBuilder.build(r, 99), s2 = SerpentPlotBuilder.build(back, 99);
        check("serpent view identical after reload", SerpentPlotBuilder.toJson(s1).equals(SerpentPlotBuilder.toJson(s2)));

        if (failures == 0) System.out.println("PASS: all MultiLocusResult persistence tests passed");
        else { System.out.println("FAIL: " + failures + " persistence check(s) failed"); System.exit(1); }
    }

    static MultiLocusResult fixture() {
        MultiLocusResult r = new MultiLocusResult();
        r.name = "run \"one\""; r.refPanelId = "eur_1000g"; r.refPanelLabel = "1000 Genomes EUR"; r.createdAt = "2026-09-28T20:00:00Z";
        String[] ids = {"aaa-1", "aaa-2", "bbb-1", "bbb-2"};
        for (int i = 0; i < ids.length; i++) {
            MultiLocusResult.DatasetInfo d = new MultiLocusResult.DatasetInfo(ids[i], "Dataset " + i, i % 2 == 0 ? "OR" : "Beta");
            d.diseaseName = i < 2 ? "Alpha" : "";
            d.ancestry = i == 0 ? "EUR (Spanish)" : "EAS";
            d.sampleN = i == 0 ? 12345 : 1000 * (i + 1);
            if (i == 0) d.categories.put("tissue", "blood");
            d.extraColumns.add("note");
            r.datasets.add(d);
        }
        for (int li = 0; li < 3; li++) {
            MultiLocusResult.LocusRow row = new MultiLocusResult.LocusRow();
            row.index = li; row.chr = li == 2 ? "X" : String.valueOf(li + 1);
            row.start = 1_000_000L * (li + 1); row.end = row.start + 50_000; row.sizeBp = 50_000;
            row.nearestGene = li == 0 ? "GENE1, AC012345.1" : ""; row.geneDist = li == 0 ? 0 : -1;
            row.overallBestDatasetId = "aaa-1"; row.overallBestSnpId = "rs" + li; row.overallBestP = 3.25e-12;
            for (String id : ids) {
                MultiLocusResult.DatasetLocusStat s = new MultiLocusResult.DatasetLocusStat();
                s.bestSnpId = "rs" + li + id; s.bestChr = row.chr; s.bestPos = li == 0 && id.equals("aaa-1") ? 123_456_789L : row.start + 10;
                s.bestP = id.equals("aaa-1") ? 3.25e-12 : 0.01 * (li + 1);
                s.ea = "A"; s.nea = "G"; s.or = id.startsWith("aaa") ? 1.2 : Double.NaN; s.se = 0.05;
                s.n = 5000; s.nP5e5 = 2; s.nP5e8 = id.equals("aaa-1") ? 1 : 0;
                if (li == 0 && id.equals("aaa-1")) s.extra.put("note", "a \"quoted\"\tvalue\\with\nnewline");
                row.cells.put(id, s);
            }
            r.loci.add(row);
        }
        return r;
    }

    static void check(String name, boolean ok) {
        if (ok) System.out.println("PASS: " + name);
        else { System.out.println("FAIL: " + name); failures++; }
    }
}
