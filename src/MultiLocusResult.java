import java.util.*;

/** Result of a Locus Matrix job: the dataset x locus significance grid. */
public class MultiLocusResult {
    public String name          = "";
    public String refPanelId    = "";
    public String refPanelLabel = "";
    public String createdAt     = "";
    public List<DatasetInfo> datasets = new ArrayList<>();
    public List<LocusRow> loci        = new ArrayList<>();

    public static class DatasetInfo {
        public String id;
        public String name;
        public String effectType; // "OR" or "Beta"
        public String diseaseName = ""; // explicit Config.diseaseName, "" if the user never set it
        public String ancestry    = ""; // Config.ancestry (falls back to the ref-panel population)
        public int    sampleN     = 0;  // Config.sampleN, 0 if unknown
        public Map<String, String> categories = new LinkedHashMap<>(); // Config.categories (category.<name>=...)
        public List<String> extraColumns = new ArrayList<>(); // that dataset's own unmapped header columns

        public DatasetInfo(String id, String name, String effectType) {
            this.id = id; this.name = name; this.effectType = effectType;
        }
    }

    public static class DatasetLocusStat {
        public String bestSnpId = "";
        public String bestChr   = "";
        public long   bestPos   = 0;
        public double bestP     = Double.NaN;
        public String ea        = "";
        public String nea       = "";
        public double beta      = Double.NaN;
        public double or        = Double.NaN;
        public double se        = Double.NaN;
        public double maf       = Double.NaN;
        public double n         = Double.NaN;
        public double info      = Double.NaN;
        public int    nP5e5     = 0;
        public int    nP5e8     = 0;
        public Map<String, String> extra = new LinkedHashMap<>(); // dataset's own unmapped columns, for the best SNP
    }

    public static class LocusRow {
        public int    index;
        public String chr;
        public long   start, end;
        public long   sizeBp;
        public String nearestGene = "";
        public long   geneDist    = -1;
        public String overallBestDatasetId = "";
        public String overallBestSnpId     = "";
        public double overallBestP         = Double.NaN;
        public Map<String, DatasetLocusStat> cells = new LinkedHashMap<>();
    }

    public String toJson() {
        StringBuilder j = new StringBuilder();
        j.append("{");
        kv(j, "name", name); j.append(",");
        kv(j, "ref_panel_id", refPanelId); j.append(",");
        kv(j, "ref_panel_label", refPanelLabel); j.append(",");
        kv(j, "created_at", createdAt); j.append(",");

        j.append("\"datasets\":[");
        for (int i = 0; i < datasets.size(); i++) {
            if (i > 0) j.append(",");
            DatasetInfo d = datasets.get(i);
            j.append("{");
            kv(j, "id", d.id); j.append(",");
            kv(j, "name", d.name); j.append(",");
            kv(j, "effect_type", d.effectType); j.append(",");
            kv(j, "disease_name", d.diseaseName); j.append(",");
            kv(j, "ancestry", d.ancestry); j.append(",");
            j.append("\"sample_n\":").append(d.sampleN).append(",");
            j.append("\"categories\":{");
            int ci = 0;
            for (Map.Entry<String, String> c : d.categories.entrySet()) {
                if (ci++ > 0) j.append(",");
                kv(j, c.getKey(), c.getValue());
            }
            j.append("},");
            j.append("\"extra_columns\":[");
            for (int k = 0; k < d.extraColumns.size(); k++) {
                if (k > 0) j.append(",");
                j.append("\"").append(esc(d.extraColumns.get(k))).append("\"");
            }
            j.append("]");
            j.append("}");
        }
        j.append("],");

        j.append("\"loci\":[");
        for (int i = 0; i < loci.size(); i++) {
            if (i > 0) j.append(",");
            appendLocus(j, loci.get(i));
        }
        j.append("]");
        j.append("}");
        return j.toString();
    }

    private void appendLocus(StringBuilder j, LocusRow l) {
        j.append("{");
        j.append("\"index\":").append(l.index).append(",");
        kv(j, "chr", l.chr); j.append(",");
        j.append("\"start\":").append(l.start).append(",");
        j.append("\"end\":").append(l.end).append(",");
        j.append("\"size_bp\":").append(l.sizeBp).append(",");
        kv(j, "nearest_gene", l.nearestGene); j.append(",");
        j.append("\"gene_dist\":").append(l.geneDist).append(",");
        kv(j, "overall_best_dataset_id", l.overallBestDatasetId); j.append(",");
        kv(j, "overall_best_snp_id", l.overallBestSnpId); j.append(",");
        j.append("\"overall_best_p\":").append(num(l.overallBestP)).append(",");
        j.append("\"cells\":{");
        int i = 0;
        for (Map.Entry<String, DatasetLocusStat> e : l.cells.entrySet()) {
            if (i++ > 0) j.append(",");
            j.append("\"").append(esc(e.getKey())).append("\":");
            appendStat(j, e.getValue());
        }
        j.append("}");
        j.append("}");
    }

    private void appendStat(StringBuilder j, DatasetLocusStat s) {
        j.append("{");
        kv(j, "best_snp_id", s.bestSnpId); j.append(",");
        kv(j, "best_chr", s.bestChr); j.append(",");
        j.append("\"best_pos\":").append(s.bestPos).append(",");
        j.append("\"best_p\":").append(num(s.bestP)).append(",");
        kv(j, "ea", s.ea); j.append(",");
        kv(j, "nea", s.nea); j.append(",");
        j.append("\"beta\":").append(num(s.beta)).append(",");
        j.append("\"or\":").append(num(s.or)).append(",");
        j.append("\"se\":").append(num(s.se)).append(",");
        j.append("\"maf\":").append(num(s.maf)).append(",");
        j.append("\"n\":").append(num(s.n)).append(",");
        j.append("\"info\":").append(num(s.info)).append(",");
        j.append("\"n_p5e5\":").append(s.nP5e5).append(",");
        j.append("\"n_p5e8\":").append(s.nP5e8).append(",");
        j.append("\"extra\":{");
        int ei = 0;
        for (Map.Entry<String, String> e : s.extra.entrySet()) {
            if (ei++ > 0) j.append(",");
            kv(j, e.getKey(), e.getValue());
        }
        j.append("}");
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
        StringBuilder b = new StringBuilder(s.length() + 8);
        for (int i = 0; i < s.length(); i++) {
            char c = s.charAt(i);
            switch (c) {
                case '\\': b.append("\\\\"); break;
                case '"':  b.append("\\\""); break;
                case '\n': b.append("\\n"); break;
                case '\r': b.append("\\r"); break;
                case '\t': b.append("\\t"); break;
                default:
                    if (c < 0x20) b.append(String.format("\\u%04x", (int) c)); else b.append(c);
            }
        }
        return b.toString();
    }

    // ── reading back (persisted runs: output/multi_locus/<job>/result.json) ──────────────────

    /** Inverse of {@link #toJson()}: rebuilds a result saved to disk, so runs survive a restart. */
    public static MultiLocusResult fromJson(String json) {
        Map<String, Object> o = MiniJson.asObject(MiniJson.parse(json));
        MultiLocusResult r = new MultiLocusResult();
        r.name          = str(o, "name");
        r.refPanelId    = str(o, "ref_panel_id");
        r.refPanelLabel = str(o, "ref_panel_label");
        r.createdAt     = str(o, "created_at");
        for (Object dv : list(o, "datasets")) {
            Map<String, Object> d = MiniJson.asObject(dv);
            DatasetInfo di = new DatasetInfo(str(d, "id"), str(d, "name"), str(d, "effect_type"));
            di.diseaseName = str(d, "disease_name");
            di.ancestry    = str(d, "ancestry");
            di.sampleN     = (int) lng(d, "sample_n");
            Object cats = d.get("categories");
            if (cats instanceof Map) for (Map.Entry<String, Object> e : MiniJson.asObject(cats).entrySet())
                di.categories.put(e.getKey(), e.getValue() == null ? "" : String.valueOf(e.getValue()));
            for (Object x : list(d, "extra_columns")) di.extraColumns.add(String.valueOf(x));
            r.datasets.add(di);
        }
        for (Object lv : list(o, "loci")) {
            Map<String, Object> l = MiniJson.asObject(lv);
            LocusRow row = new LocusRow();
            row.index  = (int) lng(l, "index");
            row.chr    = str(l, "chr");
            row.start  = lng(l, "start");
            row.end    = lng(l, "end");
            row.sizeBp = lng(l, "size_bp");
            row.nearestGene = str(l, "nearest_gene");
            row.geneDist    = lng(l, "gene_dist");
            row.overallBestDatasetId = str(l, "overall_best_dataset_id");
            row.overallBestSnpId     = str(l, "overall_best_snp_id");
            row.overallBestP         = dbl(l, "overall_best_p");
            Object cells = l.get("cells");
            if (cells instanceof Map) for (Map.Entry<String, Object> e : MiniJson.asObject(cells).entrySet()) {
                Map<String, Object> c = MiniJson.asObject(e.getValue());
                DatasetLocusStat st = new DatasetLocusStat();
                st.bestSnpId = str(c, "best_snp_id");
                st.bestChr   = str(c, "best_chr");
                st.bestPos   = lng(c, "best_pos");
                st.bestP     = dbl(c, "best_p");
                st.ea = str(c, "ea"); st.nea = str(c, "nea");
                st.beta = dbl(c, "beta"); st.or = dbl(c, "or"); st.se = dbl(c, "se");
                st.maf = dbl(c, "maf"); st.n = dbl(c, "n"); st.info = dbl(c, "info");
                st.nP5e5 = (int) lng(c, "n_p5e5"); st.nP5e8 = (int) lng(c, "n_p5e8");
                Object ex = c.get("extra");
                if (ex instanceof Map) for (Map.Entry<String, Object> x : MiniJson.asObject(ex).entrySet())
                    st.extra.put(x.getKey(), x.getValue() == null ? "" : String.valueOf(x.getValue()));
                row.cells.put(e.getKey(), st);
            }
            r.loci.add(row);
        }
        return r;
    }

    private static String str(Map<String, Object> o, String k) {
        Object v = o.get(k); return v == null ? "" : String.valueOf(v);
    }
    private static double dbl(Map<String, Object> o, String k) {
        Object v = o.get(k); return v instanceof Number ? ((Number) v).doubleValue() : Double.NaN;
    }
    private static long lng(Map<String, Object> o, String k) {
        Object v = o.get(k); return v instanceof Number ? ((Number) v).longValue() : 0L;
    }
    private static List<Object> list(Map<String, Object> o, String k) {
        Object v = o.get(k); return v instanceof List ? MiniJson.asArray(v) : Collections.emptyList();
    }
}
