import java.util.*;

/** JSON-serializable result of {@link GeneConstellationBuilder#build}: one entry per gene,
 *  derived from an already-computed MultiLocusResult. Hand-rolled JSON, same style as
 *  MultiLocusResult.toJson() (kv/num/esc helpers). */
public class GeneConstellationResult {
    public String name          = "";
    public String refPanelLabel = "";
    public String createdAt     = "";
    public double threshold     = GeneConstellationBuilder.DEFAULT_THRESHOLD;
    public int unassignedLociCount = 0;
    public List<MultiLocusResult.DatasetInfo> datasets = new ArrayList<>();
    public List<GeneEntry> genes = new ArrayList<>();

    public static class PerDataset {
        public String id, name, disease;
        public double bestP = Double.NaN, beta = Double.NaN, or = Double.NaN;
    }

    public static class AnovaSummary {
        public double fStat, pValue;
        public int dfBetween, dfWithin;
    }

    public static class GeneEntry {
        public String gene = "";
        public String chr  = "";
        public long   pos  = 0;
        public int    nDatasetsSignificant = 0;
        public int    nDatasetsTotal       = 0;
        public double aggNegLog10P         = Double.NaN;
        public double directionFractionRisk = Double.NaN; // fraction of cells with beta>0/OR>1; NaN = no direction data
        public AnovaSummary betweenDiseaseAnova = null;
        public Map<String, AnovaSummary> withinDiseaseAnova = new LinkedHashMap<>(); // disease -> summary or null
        public List<PerDataset> perDataset = new ArrayList<>();
    }

    public String toJson() {
        StringBuilder j = new StringBuilder();
        j.append("{");
        kv(j, "name", name); j.append(",");
        kv(j, "ref_panel_label", refPanelLabel); j.append(",");
        kv(j, "created_at", createdAt); j.append(",");
        j.append("\"threshold\":").append(num(threshold)).append(",");
        j.append("\"unassigned_loci_count\":").append(unassignedLociCount).append(",");

        j.append("\"datasets\":[");
        for (int i = 0; i < datasets.size(); i++) {
            if (i > 0) j.append(",");
            MultiLocusResult.DatasetInfo d = datasets.get(i);
            j.append("{");
            kv(j, "id", d.id); j.append(",");
            kv(j, "name", d.name); j.append(",");
            kv(j, "effect_type", d.effectType);
            j.append("}");
        }
        j.append("],");

        j.append("\"genes\":[");
        for (int i = 0; i < genes.size(); i++) {
            if (i > 0) j.append(",");
            appendGene(j, genes.get(i));
        }
        j.append("]");
        j.append("}");
        return j.toString();
    }

    private void appendGene(StringBuilder j, GeneEntry g) {
        j.append("{");
        kv(j, "gene", g.gene); j.append(",");
        kv(j, "chr", g.chr); j.append(",");
        j.append("\"pos\":").append(g.pos).append(",");
        j.append("\"n_datasets_significant\":").append(g.nDatasetsSignificant).append(",");
        j.append("\"n_datasets_total\":").append(g.nDatasetsTotal).append(",");
        j.append("\"agg_neglog10_p\":").append(num(g.aggNegLog10P)).append(",");
        j.append("\"direction_fraction_risk\":").append(num(g.directionFractionRisk)).append(",");

        j.append("\"between_disease_anova\":");
        appendAnova(j, g.betweenDiseaseAnova);
        j.append(",");

        j.append("\"within_disease_anova\":{");
        int wi = 0;
        for (Map.Entry<String, AnovaSummary> e : g.withinDiseaseAnova.entrySet()) {
            if (wi++ > 0) j.append(",");
            j.append("\"").append(esc(e.getKey())).append("\":");
            appendAnova(j, e.getValue());
        }
        j.append("},");

        j.append("\"per_dataset\":[");
        for (int i = 0; i < g.perDataset.size(); i++) {
            if (i > 0) j.append(",");
            PerDataset p = g.perDataset.get(i);
            j.append("{");
            kv(j, "id", p.id); j.append(",");
            kv(j, "name", p.name); j.append(",");
            kv(j, "disease", p.disease); j.append(",");
            j.append("\"best_p\":").append(num(p.bestP)).append(",");
            j.append("\"beta\":").append(num(p.beta)).append(",");
            j.append("\"or\":").append(num(p.or));
            j.append("}");
        }
        j.append("]");
        j.append("}");
    }

    private void appendAnova(StringBuilder j, AnovaSummary a) {
        if (a == null) { j.append("null"); return; }
        j.append("{");
        j.append("\"f_stat\":").append(num(a.fStat)).append(",");
        j.append("\"p_value\":").append(num(a.pValue)).append(",");
        j.append("\"df_between\":").append(a.dfBetween).append(",");
        j.append("\"df_within\":").append(a.dfWithin);
        j.append("}");
    }

    private static void kv(StringBuilder j, String key, String value) {
        j.append("\"").append(key).append("\":\"").append(esc(value)).append("\"");
    }

    private static String num(double d) {
        return Double.isNaN(d) || Double.isInfinite(d) ? "null" : String.valueOf(d);
    }

    private static String esc(String s) {
        if (s == null) return "";
        return s.replace("\\", "\\\\").replace("\"", "\\\"");
    }
}
