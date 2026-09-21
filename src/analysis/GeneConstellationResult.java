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
    public List<GeneEdge> edges = new ArrayList<>();
    /** Region-level view (DECISIONS_PHASE5.md section 3.2): one entry per already-pooled
     *  cross-dataset locus, not collapsed by nearest gene — see {@link RegionConstellationBuilder}. */
    public List<RegionEntry> regions = new ArrayList<>();

    /** A "significant together" edge between two genes: the number of datasets in which both
     *  independently reach the significance threshold (not a correlation or shared-variant claim). */
    public static class GeneEdge {
        public String geneA = "", geneB = "";
        public int count = 0;
    }

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
        /** Same ANOVA machinery, applied to ln(effect) instead of -log10(p) — see DECISIONS_PHASE5.md
         *  section 2 for why ln(effect) rather than raw OR, and the stated cross-study-allele-coding
         *  caveat this carries. Null when fewer than 2 non-NaN-effect groups exist, same convention as
         *  betweenDiseaseAnova/withinDiseaseAnova. */
        public AnovaSummary betweenDiseaseAnovaEffect = null;
        public Map<String, AnovaSummary> withinDiseaseAnovaEffect = new LinkedHashMap<>();
        public List<PerDataset> perDataset = new ArrayList<>();
    }

    /** One already-pooled cross-dataset locus (DECISIONS_PHASE5.md section 3.2) — a genomic *region*,
     *  not collapsed by nearest gene, so a region with no assignable nearest gene still gets a node
     *  instead of being dropped (unlike the gene-level view, where such loci are only counted in
     *  {@code unassignedLociCount}). */
    public static class RegionEntry {
        public String regionId = "";   // "chrN:start-end"
        public String chr = "";
        public long   start = 0, end = 0;
        public String nearestGene = ""; // "" if none within this app's gene-assignment window
        public int    nDatasetsSignificant = 0;
        public int    nDatasetsTotal       = 0;
        public double aggNegLog10P          = Double.NaN;
        public double directionFractionRisk = Double.NaN;
        public AnovaSummary betweenDiseaseAnova       = null;
        public AnovaSummary betweenDiseaseAnovaEffect = null;
        public Map<String, AnovaSummary> withinDiseaseAnova       = new LinkedHashMap<>();
        public Map<String, AnovaSummary> withinDiseaseAnovaEffect = new LinkedHashMap<>();
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
        j.append("],");

        j.append("\"edges\":[");
        for (int i = 0; i < edges.size(); i++) {
            if (i > 0) j.append(",");
            GeneEdge e = edges.get(i);
            j.append("{");
            kv(j, "gene_a", e.geneA); j.append(",");
            kv(j, "gene_b", e.geneB); j.append(",");
            j.append("\"count\":").append(e.count);
            j.append("}");
        }
        j.append("],");

        j.append("\"regions\":[");
        for (int i = 0; i < regions.size(); i++) {
            if (i > 0) j.append(",");
            appendRegion(j, regions.get(i));
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

        j.append("\"between_disease_anova_effect\":");
        appendAnova(j, g.betweenDiseaseAnovaEffect);
        j.append(",");

        j.append("\"within_disease_anova_effect\":{");
        wi = 0;
        for (Map.Entry<String, AnovaSummary> e : g.withinDiseaseAnovaEffect.entrySet()) {
            if (wi++ > 0) j.append(",");
            j.append("\"").append(esc(e.getKey())).append("\":");
            appendAnova(j, e.getValue());
        }
        j.append("},");

        j.append("\"per_dataset\":[");
        appendPerDatasetArray(j, g.perDataset);
        j.append("]");
        j.append("}");
    }

    private void appendRegion(StringBuilder j, RegionEntry r) {
        j.append("{");
        kv(j, "region_id", r.regionId); j.append(",");
        kv(j, "chr", r.chr); j.append(",");
        j.append("\"start\":").append(r.start).append(",");
        j.append("\"end\":").append(r.end).append(",");
        kv(j, "nearest_gene", r.nearestGene); j.append(",");
        j.append("\"n_datasets_significant\":").append(r.nDatasetsSignificant).append(",");
        j.append("\"n_datasets_total\":").append(r.nDatasetsTotal).append(",");
        j.append("\"agg_neglog10_p\":").append(num(r.aggNegLog10P)).append(",");
        j.append("\"direction_fraction_risk\":").append(num(r.directionFractionRisk)).append(",");

        j.append("\"between_disease_anova\":");
        appendAnova(j, r.betweenDiseaseAnova);
        j.append(",");
        j.append("\"between_disease_anova_effect\":");
        appendAnova(j, r.betweenDiseaseAnovaEffect);
        j.append(",");

        j.append("\"within_disease_anova\":{");
        int wi = 0;
        for (Map.Entry<String, AnovaSummary> e : r.withinDiseaseAnova.entrySet()) {
            if (wi++ > 0) j.append(",");
            j.append("\"").append(esc(e.getKey())).append("\":");
            appendAnova(j, e.getValue());
        }
        j.append("},");
        j.append("\"within_disease_anova_effect\":{");
        wi = 0;
        for (Map.Entry<String, AnovaSummary> e : r.withinDiseaseAnovaEffect.entrySet()) {
            if (wi++ > 0) j.append(",");
            j.append("\"").append(esc(e.getKey())).append("\":");
            appendAnova(j, e.getValue());
        }
        j.append("},");

        j.append("\"per_dataset\":[");
        appendPerDatasetArray(j, r.perDataset);
        j.append("]");
        j.append("}");
    }

    private void appendPerDatasetArray(StringBuilder j, List<PerDataset> list) {
        for (int i = 0; i < list.size(); i++) {
            if (i > 0) j.append(",");
            PerDataset p = list.get(i);
            j.append("{");
            kv(j, "id", p.id); j.append(",");
            kv(j, "name", p.name); j.append(",");
            kv(j, "disease", p.disease); j.append(",");
            j.append("\"best_p\":").append(num(p.bestP)).append(",");
            j.append("\"beta\":").append(num(p.beta)).append(",");
            j.append("\"or\":").append(num(p.or));
            j.append("}");
        }
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
