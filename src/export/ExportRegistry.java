package export;

import java.util.*;

import static export.ColumnSpec.ColType.*;

/**
 * Central registry of column providers. Adding a new annotation =
 * implement SnpColumnProvider or LocusColumnProvider and register it here.
 * The ExcelExporter iterates all registered providers — no hardcoded columns.
 */
public class ExportRegistry {

    private static final List<SnpColumnProvider> snpProviders = new ArrayList<>();
    private static final List<LocusColumnProvider> locusProviders = new ArrayList<>();

    static {
        // Register built-in providers in display order
        snpProviders.add(new LocusIdentitySnpProvider());
        snpProviders.add(new GwasPassthroughProvider());
        snpProviders.add(new RsidProvider());
        snpProviders.add(new LeadFlagsProvider());
        snpProviders.add(new LdWithLeadProvider());
        snpProviders.add(new NearestGeneSnpProvider());

        locusProviders.add(new LocusCoreProvider());
        locusProviders.add(new LocusLeadProvider());
        locusProviders.add(new LocusNearestGeneProvider());
        locusProviders.add(new LocusCountsProvider());
    }

    public static void registerSnpProvider(SnpColumnProvider p) { snpProviders.add(p); }
    public static void registerLocusProvider(LocusColumnProvider p) { locusProviders.add(p); }

    public static List<SnpColumnProvider> snpProviders() { return Collections.unmodifiableList(snpProviders); }
    public static List<LocusColumnProvider> locusProviders() { return Collections.unmodifiableList(locusProviders); }

    // ── Built-in SNP providers ───────────────────────────────────

    static class LocusIdentitySnpProvider implements SnpColumnProvider {
        final List<ColumnSpec> cols = Arrays.asList(
            new ColumnSpec("locus_index", "Locus", "Locus", INT, 0),
            new ColumnSpec("locus_name",  "Locus Name", "Locus", TEXT, 1),
            new ColumnSpec("locus_chr",   "Locus Chr", "Locus", TEXT, 2),
            new ColumnSpec("locus_start", "Locus Start", "Locus", INT, 3),
            new ColumnSpec("locus_end",   "Locus End", "Locus", INT, 4));
        public List<ColumnSpec> columns() { return cols; }
        public Object value(SnpContext c, ColumnSpec col) {
            switch (col.key) {
                case "locus_index": return c.locusIndex;
                case "locus_name":  return c.locusName;
                case "locus_chr":   return c.locusChr;
                case "locus_start": return c.locusStart;
                case "locus_end":   return c.locusEnd;
                default: return "";
            }
        }
        public boolean isAvailable(String p) { return true; }
    }

    static class GwasPassthroughProvider implements SnpColumnProvider {
        final List<ColumnSpec> cols = Arrays.asList(
            new ColumnSpec("chr",    "Chr",     "GWAS", TEXT,       10),
            new ColumnSpec("pos",    "Position","GWAS", INT,        11),
            new ColumnSpec("ea",     "EA",      "GWAS", TEXT,       12),
            new ColumnSpec("nea",    "NEA",     "GWAS", TEXT,       13),
            new ColumnSpec("pvalue", "P-value", "GWAS", SCIENTIFIC, 14),
            new ColumnSpec("neglog", "-log10(P)","GWAS",DOUBLE,     15),
            new ColumnSpec("beta",   "Beta",    "GWAS", DOUBLE,     16),
            new ColumnSpec("or",     "OR",      "GWAS", DOUBLE,     17),
            new ColumnSpec("se",     "SE",      "GWAS", DOUBLE,     18),
            new ColumnSpec("n",      "N",       "GWAS", INT,        19),
            new ColumnSpec("maf",    "MAF",     "GWAS", DOUBLE,     20),
            new ColumnSpec("info",   "Info",    "GWAS", DOUBLE,     21));
        public List<ColumnSpec> columns() { return cols; }
        public Object value(SnpContext c, ColumnSpec col) {
            switch (col.key) {
                case "chr":    return c.chr;
                case "pos":    return c.pos;
                case "ea":     return c.ea;
                case "nea":    return c.nea;
                case "pvalue": return c.pvalue;
                case "neglog": return c.negLog10P;
                case "beta":   return Double.isNaN(c.beta) ? "-" : c.beta;
                case "or":     return Double.isNaN(c.oddsRatio) ? "-" : c.oddsRatio;
                case "se":     return Double.isNaN(c.se) ? "-" : c.se;
                case "n":      return Double.isNaN(c.sampleN) ? "-" : (int) c.sampleN;
                case "maf":    return Double.isNaN(c.maf) ? "-" : c.maf;
                case "info":   return Double.isNaN(c.info) ? "-" : c.info;
                default: return "";
            }
        }
        public boolean isAvailable(String p) { return true; }
    }

    static class RsidProvider implements SnpColumnProvider {
        final List<ColumnSpec> cols = Collections.singletonList(
            new ColumnSpec("rsid", "rsID", "Identity", TEXT, 30));
        public List<ColumnSpec> columns() { return cols; }
        public Object value(SnpContext c, ColumnSpec col) {
            return (c.id != null && !c.id.isEmpty()) ? c.id : "-";
        }
        public boolean isAvailable(String p) { return true; }
    }

    static class LeadFlagsProvider implements SnpColumnProvider {
        final List<ColumnSpec> cols = Arrays.asList(
            new ColumnSpec("is_lead",  "Is Lead", "Computed", BOOL, 40),
            new ColumnSpec("dist_lead","Dist to Lead (bp)", "Computed", INT, 41));
        public List<ColumnSpec> columns() { return cols; }
        public Object value(SnpContext c, ColumnSpec col) {
            switch (col.key) {
                case "is_lead":   return c.isLead;
                case "dist_lead": return c.distanceToLead;
                default: return "";
            }
        }
        public boolean isAvailable(String p) { return true; }
    }

    static class LdWithLeadProvider implements SnpColumnProvider {
        final List<ColumnSpec> cols = Collections.singletonList(
            new ColumnSpec("r2_lead", "r² with Lead", "LD", DOUBLE, 50));
        public List<ColumnSpec> columns() { return cols; }
        public Object value(SnpContext c, ColumnSpec col) {
            return c.r2WithLead >= 0 ? c.r2WithLead : "-";
        }
        public boolean isAvailable(String p) { return true; }
    }

    static class NearestGeneSnpProvider implements SnpColumnProvider {
        final List<ColumnSpec> cols = Arrays.asList(
            new ColumnSpec("nearest_gene", "Nearest Gene", "Gene", TEXT, 60),
            new ColumnSpec("gene_dist",    "Gene Dist (bp)", "Gene", INT, 61));
        public List<ColumnSpec> columns() { return cols; }
        public Object value(SnpContext c, ColumnSpec col) {
            switch (col.key) {
                case "nearest_gene": return c.nearestGene.isEmpty() ? "-" : c.nearestGene;
                case "gene_dist":    return c.nearestGeneDist == Long.MAX_VALUE ? "-" : c.nearestGeneDist;
                default: return "";
            }
        }
        public boolean isAvailable(String p) { return true; }
    }

    // ── Built-in Locus providers ─────────────────────────────────

    static class LocusCoreProvider implements LocusColumnProvider {
        final List<ColumnSpec> cols = Arrays.asList(
            new ColumnSpec("index", "Locus", "Core", INT, 0),
            new ColumnSpec("name",  "Name",  "Core", TEXT, 1),
            new ColumnSpec("chr",   "Chr",   "Core", TEXT, 2),
            new ColumnSpec("start", "Start", "Core", INT, 3),
            new ColumnSpec("end",   "End",   "Core", INT, 4),
            new ColumnSpec("size",  "Size (bp)", "Core", INT, 5));
        public List<ColumnSpec> columns() { return cols; }
        public Object value(LocusContext c, ColumnSpec col) {
            switch (col.key) {
                case "index": return c.index;
                case "name":  return c.name;
                case "chr":   return c.chr;
                case "start": return c.start;
                case "end":   return c.end;
                case "size":  return c.sizeBp;
                default: return "";
            }
        }
        public boolean isAvailable(String p) { return true; }
    }

    static class LocusLeadProvider implements LocusColumnProvider {
        final List<ColumnSpec> cols = Arrays.asList(
            new ColumnSpec("lead_snp", "Lead SNP",   "Lead", TEXT,       10),
            new ColumnSpec("lead_pos", "Lead Pos",   "Lead", INT,        11),
            new ColumnSpec("lead_p",   "Lead P",     "Lead", SCIENTIFIC, 12),
            new ColumnSpec("lead_beta","Lead Beta",  "Lead", DOUBLE,     13),
            new ColumnSpec("lead_or",  "Lead OR",    "Lead", DOUBLE,     14));
        public List<ColumnSpec> columns() { return cols; }
        public Object value(LocusContext c, ColumnSpec col) {
            switch (col.key) {
                case "lead_snp":  return c.leadSnpId;
                case "lead_pos":  return c.leadSnpPos;
                case "lead_p":    return Double.isNaN(c.leadP) ? "-" : c.leadP;
                case "lead_beta": return Double.isNaN(c.leadBeta) ? "-" : c.leadBeta;
                case "lead_or":   return Double.isNaN(c.leadOr) ? "-" : c.leadOr;
                default: return "";
            }
        }
        public boolean isAvailable(String p) { return true; }
    }

    static class LocusNearestGeneProvider implements LocusColumnProvider {
        final List<ColumnSpec> cols = Arrays.asList(
            new ColumnSpec("nearest_gene", "Nearest Gene", "Gene", TEXT, 20),
            new ColumnSpec("gene_dist",    "Gene Dist (bp)","Gene",INT,  21));
        public List<ColumnSpec> columns() { return cols; }
        public Object value(LocusContext c, ColumnSpec col) {
            switch (col.key) {
                case "nearest_gene": return c.nearestGene.isEmpty() ? "-" : c.nearestGene;
                case "gene_dist":    return c.nearestGeneDist == Long.MAX_VALUE ? "-" : c.nearestGeneDist;
                default: return "";
            }
        }
        public boolean isAvailable(String p) { return true; }
    }

    static class LocusCountsProvider implements LocusColumnProvider {
        final List<ColumnSpec> cols = Arrays.asList(
            new ColumnSpec("n_snps",   "N SNPs",       "Counts", INT, 30),
            new ColumnSpec("n_p5e3",   "N p<5e-3",     "Counts", INT, 31),
            new ColumnSpec("n_p5e5",   "N p<5e-5",     "Counts", INT, 32),
            new ColumnSpec("n_p5e8",   "N p<5e-8",     "Counts", INT, 33),
            new ColumnSpec("ld_done",  "LD Computed",   "LD",     BOOL, 34),
            new ColumnSpec("ref_panel","Ref Panel",     "LD",     TEXT, 35));
        public List<ColumnSpec> columns() { return cols; }
        public Object value(LocusContext c, ColumnSpec col) {
            switch (col.key) {
                case "n_snps":    return c.nSnps;
                case "n_p5e3":    return c.nP5e3;
                case "n_p5e5":    return c.nP5e5;
                case "n_p5e8":    return c.nP5e8;
                case "ld_done":   return c.ldComputed;
                case "ref_panel": return c.refPanel;
                default: return "";
            }
        }
        public boolean isAvailable(String p) { return true; }
    }
}
