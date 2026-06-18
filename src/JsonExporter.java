import java.io.*;
import java.util.*;

/** Writes per-locus JSON + .js wrapper files and manifest.json/manifest.js */
public class JsonExporter {

    public static void export(List<Locus> loci,
                              List<LocusOutput> outputs,
                              Config config) throws IOException {
        String dataDir = config.outputDir + "/data";
        new File(dataDir).mkdirs();

        for (LocusOutput lo : outputs) {
            String json = locusToJson(lo);
            try (PrintWriter pw = new PrintWriter(new BufferedWriter(
                    new FileWriter(dataDir + "/locus_" + lo.locusIndex + ".json")))) {
                pw.print(json);
            }
            try (PrintWriter pw = new PrintWriter(new BufferedWriter(
                    new FileWriter(dataDir + "/locus_" + lo.locusIndex + ".js")))) {
                pw.print("(function(){window.LOCUS_DATA=window.LOCUS_DATA||{};");
                pw.print("window.LOCUS_DATA[" + lo.locusIndex + "]=" + json + ";");
                pw.print("})();");
            }
        }

        String manifestJson = manifestToJson(outputs);
        try (PrintWriter pw = new PrintWriter(new BufferedWriter(
                new FileWriter(dataDir + "/manifest.json")))) {
            pw.print(manifestJson);
        }
        try (PrintWriter pw = new PrintWriter(new BufferedWriter(
                new FileWriter(dataDir + "/manifest.js")))) {
            pw.print("(function(){window.MANIFEST=" + manifestJson + ";})();");
        }

        System.out.printf("[JsonExporter] Wrote %d locus files + manifest to %s/%n",
            outputs.size(), dataDir);
    }

    // ── Per-locus JSON ────────────────────────────────────────────────────────

    static String locusToJson(LocusOutput lo) {
        Jb j = new Jb();
        j.obj(() -> {
            j.kv("locus_index",   lo.locusIndex);
            j.kv("locus_name",    lo.locusName);
            j.kv("chr",           lo.chr);
            j.kv("start",         lo.start);
            j.kv("end",           lo.end);
            j.kv("padded_start",  lo.paddedStart);
            j.kv("padded_end",    lo.paddedEnd);
            j.kv("ref_panel",     lo.refPanel);

            j.key("top_snp");
            if (lo.topSnp != null) {
                j.obj(() -> {
                    j.kv("id",             lo.topSnp.id);
                    j.kv("bim_id",         lo.topSnp.bimId != null ? lo.topSnp.bimId : "");
                    j.kv("chr",            lo.topSnp.chr);
                    j.kv("pos",            lo.topSnp.pos);
                    j.kv("pvalue",         lo.topSnp.pvalue);
                    j.kv("neg_log10_p",    lo.topSnp.negLog10P);
                    j.kv("triangle_index", lo.ldTriangle != null ? lo.ldTriangle.topSnpRank : -1);
                });
            } else { j.raw("null"); }

            j.key("nearest_genes");
            j.arr(() -> lo.nearestGenes.forEach(g -> j.str(g)));

            j.key("gwas_snps");
            j.arr(() -> lo.gwasSnps.forEach(s -> j.obj(() -> {
                j.kv("id",           s.id);
                j.kv("chr",          s.chr);
                j.kv("pos",          s.pos);
                j.kv("pvalue",       s.pvalue);
                j.kv("neg_log10_p",  s.negLog10P);
                j.kv("r2_with_index", s.r2);   // -1 = no LD data
                j.kv("ea",           s.ea);
                j.kv("nea",          s.nea);
            })));

            j.key("genes");
            j.arr(() -> lo.genes.forEach(g -> j.obj(() -> {
                j.kv("gene_name", g.geneName);
                j.kv("gene_id",   g.geneId);
                j.kv("strand",    g.strand);
                j.kv("start",     g.start);
                j.kv("end",       g.end);
                j.key("transcripts");
                j.arr(() -> g.transcripts.forEach(tx -> j.obj(() -> {
                    j.kv("transcript_id",   tx.transcriptId);
                    j.kv("transcript_type", tx.transcriptType);
                    j.kv("strand",          tx.strand);
                    j.kv("start",           tx.start);
                    j.kv("end",             tx.end);
                    j.key("exons"); j.arr(() -> tx.exons.forEach(e -> exonObj(j, e)));
                    j.key("cds");   j.arr(() -> tx.cds.forEach(e   -> exonObj(j, e)));
                    j.key("utrs");  j.arr(() -> tx.utrs.forEach(e  -> exonObj(j, e)));
                })));
            })));

            // Diagnostic: count real r² values before export
            if (lo.ldTriangle != null && lo.ldTriangle.matrix != null) {
                int nonNull = 0, total = 0;
                int sz = lo.ldTriangle.snps.size();
                for (int i = 0; i < sz; i++) {
                    for (int jj = i; jj < sz; jj++) {
                        total++;
                        if (!Double.isNaN(lo.ldTriangle.matrix[i][jj])) nonNull++;
                    }
                }
                System.out.printf("[JSON] Locus %d matrix: %d/%d upper-triangle cells have real r² values%n",
                    lo.locusIndex, nonNull, total);
            }

            j.key("ld_triangle");
            writeLdTriangle(j, lo.ldTriangle);

            j.key("locus_context");
            if (lo.locusContext != null) {
                j.obj(() -> {
                    j.key("prev_locus"); writeLocusRef(j, lo.locusContext.prevLocus);
                    j.key("next_locus"); writeLocusRef(j, lo.locusContext.nextLocus);
                });
            } else { j.raw("null"); }
        });
        return j.toString();
    }

    private static void exonObj(Jb j, Exon e) {
        j.obj(() -> {
            j.kv("start", e.start);
            j.kv("end",   e.end);
            if (e.type != null) j.kv("type", e.type);
        });
    }

    private static void writeLdTriangle(Jb j, LocusOutput.LdTriangle tri) {
        if (tri == null || tri.snps == null || tri.matrix == null) { j.raw("null"); return; }
        j.obj(() -> {
            j.kv("cell_size_px", tri.cellSizePx);
            j.kv("snp_count",    tri.snpCount);
            j.key("snps");
            j.arr(() -> tri.snps.forEach(s -> j.obj(() -> {
                j.kv("id",         s.id);
                j.kv("rank",       s.rank);
                j.kv("bim_id",     s.bimId);
                j.kv("pos",        s.pos);
                j.kv("is_top_snp", s.isTopSnp);
            })));
            j.key("matrix");
            j.arr(() -> {
                for (int i = 0; i < tri.snpCount; i++) {
                    final int ii = i;
                    j.arr(() -> {
                        for (int jj = 0; jj < tri.snpCount; jj++) {
                            double v = tri.matrix[ii][jj];
                            j.raw(Double.isNaN(v) ? "null" : String.format("%.4f", v));
                        }
                    });
                }
            });
        });
    }

    private static void writeLocusRef(Jb j, LocusOutput.LocusRef ref) {
        if (ref == null) { j.raw("null"); return; }
        j.obj(() -> {
            j.kv("index",       ref.index);
            j.kv("chr",         ref.chr);
            j.kv("start",       ref.start);
            j.kv("end",         ref.end);
            j.kv("mid",         ref.mid);
            j.kv("distance_bp", ref.distanceBp);
        });
    }

    // ── Manifest JSON ─────────────────────────────────────────────────────────

    static void exportManifest(List<LocusOutput> outputs, Config config) throws IOException {
        String dataDir = config.outputDir + "/data";
        String manifestJson = manifestToJson(outputs);
        try (PrintWriter pw = new PrintWriter(new BufferedWriter(
                new FileWriter(dataDir + "/manifest.json")))) {
            pw.print(manifestJson);
        }
        try (PrintWriter pw = new PrintWriter(new BufferedWriter(
                new FileWriter(dataDir + "/manifest.js")))) {
            pw.print("(function(){window.MANIFEST=" + manifestJson + ";})();");
        }
    }

    private static String manifestToJson(List<LocusOutput> outputs) {
        List<LocusOutput> sorted = new ArrayList<>(outputs);
        sorted.sort(Comparator.comparingInt((LocusOutput o) -> Locus.chrToInt(o.chr))
            .thenComparingLong(o -> o.start));
        Jb j = new Jb();
        j.obj(() -> {
            j.kv("total_loci", sorted.size());
            j.key("loci");
            j.arr(() -> sorted.forEach(lo -> j.obj(() -> {
                j.kv("index",        lo.locusIndex);
                j.kv("chr",          lo.chr);
                j.kv("start",        lo.start);
                j.kv("end",          lo.end);
                j.kv("locus_name",   lo.locusName);
                j.kv("top_snp",      lo.topSnp != null ? lo.topSnp.id : "");
                j.kv("top_snp_pval", lo.topSnp != null ? lo.topSnp.pvalue : Double.NaN);
                j.key("nearest_genes");
                j.arr(() -> lo.nearestGenes.forEach(g -> j.str(g)));
            })));
        });
        return j.toString();
    }

    // ── Minimal JSON builder ──────────────────────────────────────────────────

    static class Jb {
        private final StringBuilder sb = new StringBuilder(4096);
        private boolean needComma = false;

        Jb obj(Runnable body) {
            comma(); sb.append('{'); needComma = false; body.run(); sb.append('}'); needComma = true; return this;
        }
        Jb arr(Runnable body) {
            comma(); sb.append('['); needComma = false; body.run(); sb.append(']'); needComma = true; return this;
        }
        void key(String k) { comma(); sb.append('"').append(escapeJson(k)).append("\":"); needComma = false; }
        void str(String v) { comma(); sb.append('"').append(escapeJson(v)).append('"'); needComma = true; }
        void kv(String k, String v)  { key(k); str(v != null ? v : ""); }
        void kv(String k, long v)    { key(k); raw(String.valueOf(v)); }
        void kv(String k, int v)     { key(k); raw(String.valueOf(v)); }
        void kv(String k, boolean v) { key(k); raw(String.valueOf(v)); }
        void kv(String k, double v)  {
            key(k);
            if (Double.isNaN(v) || Double.isInfinite(v)) raw("null");
            else if (v == (long)v && Math.abs(v) < 1e15) raw(String.valueOf((long)v));
            else raw(String.format("%.6g", v));
        }
        void raw(String v) { comma(); sb.append(v); needComma = true; }
        void comma() { if (needComma) sb.append(','); }

        @Override public String toString() { return sb.toString(); }

        private static String escapeJson(String s) {
            if (s == null) return "";
            StringBuilder out = new StringBuilder(s.length());
            for (int i = 0; i < s.length(); i++) {
                char c = s.charAt(i);
                switch (c) {
                    case '"':  out.append("\\\""); break;
                    case '\\': out.append("\\\\"); break;
                    case '\n': out.append("\\n");  break;
                    case '\r': out.append("\\r");  break;
                    case '\t': out.append("\\t");  break;
                    default:   out.append(c);
                }
            }
            return out.toString();
        }
    }
}
