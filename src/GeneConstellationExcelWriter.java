import java.io.*;
import java.util.*;
import export.XlsxWriter;

/** Builds the Gene Constellation export workbook from an already-computed
 *  GeneConstellationResult. Mirrors MultiLocusExcelWriter's exact style (static write() using
 *  export.XlsxWriter). Two sheets: one row per gene ("Genes"), and one row per (gene, dataset)
 *  pair drawn from GeneEntry.perDataset ("Per-dataset detail"). */
public class GeneConstellationExcelWriter {

    public static void write(GeneConstellationResult result, OutputStream out) throws IOException {
        try (XlsxWriter xlsx = new XlsxWriter(out)) {
            writeGenesSheet(xlsx, result);
            writePerDatasetSheet(xlsx, result);
            xlsx.finish();
        }
    }

    private static void writeGenesSheet(XlsxWriter xlsx, GeneConstellationResult result) {
        XlsxWriter.Sheet s = xlsx.addSheet("Genes").freezeHeader().autoFilter();

        // Collect the set of diseases that have at least one within-disease ANOVA entry
        // anywhere in the result, so every gene row gets the same column set.
        LinkedHashSet<String> diseases = new LinkedHashSet<>();
        for (GeneConstellationResult.GeneEntry g : result.genes) {
            diseases.addAll(g.withinDiseaseAnova.keySet());
        }

        List<Object> header = new ArrayList<>();
        header.add("Gene");
        header.add("Chr");
        header.add("Pos");
        header.add("N Datasets Significant");
        header.add("N Datasets Total");
        header.add("Agg -log10(P)");
        header.add("Direction Fraction Risk");
        header.add("Between-Disease ANOVA F");
        header.add("Between-Disease ANOVA P");
        header.add("Between-Disease df1");
        header.add("Between-Disease df2");
        for (String disease : diseases) {
            header.add(disease + " Within-Disease ANOVA F");
            header.add(disease + " Within-Disease ANOVA P");
        }
        s.addRow(header.toArray());

        for (GeneConstellationResult.GeneEntry g : result.genes) {
            List<Object> row = new ArrayList<>();
            row.add(dash(g.gene));
            row.add(dash(g.chr));
            row.add(g.pos);
            row.add(g.nDatasetsSignificant);
            row.add(g.nDatasetsTotal);
            row.add(num(g.aggNegLog10P));
            row.add(num(g.directionFractionRisk));
            if (g.betweenDiseaseAnova != null) {
                row.add(num(g.betweenDiseaseAnova.fStat));
                row.add(num(g.betweenDiseaseAnova.pValue));
                row.add(g.betweenDiseaseAnova.dfBetween);
                row.add(g.betweenDiseaseAnova.dfWithin);
            } else {
                // Not computed for this gene (fewer than 2 non-empty disease groups) —
                // leave genuinely blank rather than writing 0 or a fabricated value.
                row.add(""); row.add(""); row.add(""); row.add("");
            }
            for (String disease : diseases) {
                GeneConstellationResult.AnovaSummary a = g.withinDiseaseAnova.get(disease);
                if (a != null) {
                    row.add(num(a.fStat));
                    row.add(num(a.pValue));
                } else {
                    // Either not present for this gene, or present-but-null (only one dataset
                    // in that disease group for this gene) — both cases mean "not computed".
                    row.add(""); row.add("");
                }
            }
            s.addRow(row.toArray());
        }
    }

    private static void writePerDatasetSheet(XlsxWriter xlsx, GeneConstellationResult result) {
        XlsxWriter.Sheet s = xlsx.addSheet("Per-dataset detail").freezeHeader().autoFilter();

        List<Object> header = new ArrayList<>();
        header.add("Gene");
        header.add("Chr");
        header.add("Dataset");
        header.add("Disease");
        header.add("Best P");
        header.add("Beta");
        header.add("OR");
        s.addRow(header.toArray());

        for (GeneConstellationResult.GeneEntry g : result.genes) {
            for (GeneConstellationResult.PerDataset pd : g.perDataset) {
                List<Object> row = new ArrayList<>();
                row.add(dash(g.gene));
                row.add(dash(g.chr));
                row.add(dash(pd.name != null && !pd.name.isEmpty() ? pd.name : pd.id));
                row.add(dash(pd.disease));
                row.add(num(pd.bestP));
                row.add(num(pd.beta));
                row.add(num(pd.or));
                s.addRow(row.toArray());
            }
        }
    }

    private static Object dash(String s) {
        return (s == null || s.isEmpty()) ? "-" : s;
    }

    private static Object num(double d) {
        return Double.isNaN(d) ? "-" : d;
    }
}
