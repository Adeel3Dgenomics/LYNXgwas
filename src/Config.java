import java.io.*;
import java.util.*;

public class Config {
    // Bump this after any pipeline behavior change to force all projects to reprocess.
    public static final String PIPELINE_VERSION = "1.0.0";

    public String gwasFile           = "";
    public String lociFile           = "input/loci.txt";
    public String gff3File           = "resources/gencode.v37.annotation.gff3";
    // Reference panel — accept both old key (ref.panel) and new (ref.panel.path)
    public String refPanelPath       = "";
    public String refPanelPopulation = "EAS";
    public String colChr             = "chrom";
    public String colPos             = "pos";
    public String colPvalue          = "p";
    public String colRsid            = "";
    public String colVarid           = "varid";
    public String colEa              = "ea";
    public String colNea             = "nea";
    // Optional GWAS columns — empty string means not mapped
    public String colBeta            = "";
    public String colOr              = "";
    public String colSe              = "";
    public String colN               = "";
    public String colMaf             = "";
    public String colInfo            = "";
    public String topSnpFile         = "";
    public long   locusPadding       = 200000;
    public boolean ldEnabled         = false;
    public int    ldTriangleBoundary = 100;   // SNPs each side of top for triangle
    public double ldR2Threshold      = 0.0;
    public int    ldParallelJobs     = 4;
    public String outputDir          = "output";
    public int    threads            = 4;
    public int    maxSnpsPerLocus    = 5000;
    public double splitLdThreshold  = 0.2;
    public long   splitMinDistBp    = 250000;

    /**
     * Load config from a project directory's config.properties.
     * Output is directed into the project directory itself.
     */
    public static Config loadFromProject(String projectDir) throws IOException {
        Config c = new Config();
        String path = new File(projectDir, "config.properties").getAbsolutePath();
        c.loadProperties(path);
        c.outputDir = projectDir;
        c.validate();
        return c;
    }

    public static Config load(String[] args) throws IOException {
        Config c = new Config();
        for (int i = 0; i < args.length - 1; i++) {
            if (args[i].equals("--config")) { c.loadProperties(args[i + 1]); break; }
        }
        for (int i = 0; i < args.length; i++) {
            switch (args[i]) {
                case "--gwas":       if (i+1 < args.length) c.gwasFile       = args[++i]; break;
                case "--loci":       if (i+1 < args.length) c.lociFile       = args[++i]; break;
                case "--gff3":       if (i+1 < args.length) c.gff3File       = args[++i]; break;
                case "--ref-panel":  if (i+1 < args.length) c.refPanelPath   = args[++i]; break;
                case "--output":     if (i+1 < args.length) c.outputDir      = args[++i]; break;
                case "--padding":    if (i+1 < args.length) c.locusPadding   = Long.parseLong(args[++i]); break;
                case "--threads":    if (i+1 < args.length) c.threads        = Integer.parseInt(args[++i]); break;
                case "--no-ld":      c.ldEnabled = false; break;
                case "--ld":         c.ldEnabled = true;  break;
            }
        }
        return c;
    }

    private void loadProperties(String path) throws IOException {
        Properties p = new Properties();
        // Pre-process: escape lone backslashes (Windows paths) so Properties.load doesn't drop them
        StringBuilder sb = new StringBuilder();
        try (BufferedReader br = new BufferedReader(new FileReader(path))) {
            String line;
            while ((line = br.readLine()) != null) {
                if (!line.startsWith("#") && !line.startsWith("!") && line.contains("=")) {
                    int eq = line.indexOf('=');
                    String val = line.substring(eq + 1);
                    // Replace lone \ (not already doubled) with /  so paths work on Windows
                    val = val.replace("\\", "/");
                    line = line.substring(0, eq + 1) + val;
                }
                sb.append(line).append('\n');
            }
        }
        p.load(new java.io.StringReader(sb.toString()));

        gwasFile           = p.getProperty("gwas.file",            gwasFile);
        lociFile           = p.getProperty("loci.file",            lociFile);
        gff3File           = p.getProperty("gff3.file",            gff3File);
        // Accept both old key and new key; new key wins
        refPanelPath       = p.getProperty("ref.panel",            refPanelPath);
        refPanelPath       = p.getProperty("ref.panel.path",       refPanelPath);
        refPanelPopulation = p.getProperty("ref.panel.population", refPanelPopulation);
        colChr             = p.getProperty("col.chr",              colChr);
        colPos             = p.getProperty("col.pos",              colPos);
        colPvalue          = p.getProperty("col.pvalue",           colPvalue);
        colRsid            = p.getProperty("col.rsid",             colRsid);
        colVarid           = p.getProperty("col.varid",            colVarid);
        colEa              = p.getProperty("col.ea",               colEa);
        colNea             = p.getProperty("col.nea",              colNea);
        colBeta            = p.getProperty("col.beta",             colBeta);
        colOr              = p.getProperty("col.or",               colOr);
        colSe              = p.getProperty("col.se",               colSe);
        colN               = p.getProperty("col.n",                colN);
        colMaf             = p.getProperty("col.maf",              colMaf);
        colInfo            = p.getProperty("col.info",             colInfo);
        topSnpFile         = p.getProperty("top.snp.file",         topSnpFile);
        locusPadding       = Long.parseLong(p.getProperty("locus.padding", String.valueOf(locusPadding)));
        ldEnabled          = Boolean.parseBoolean(p.getProperty("ld.enabled", String.valueOf(ldEnabled)));
        ldTriangleBoundary = Integer.parseInt(p.getProperty("ld.triangle.snp.boundary",
                             p.getProperty("ld.triangle.boundary", String.valueOf(ldTriangleBoundary))));
        ldR2Threshold      = Double.parseDouble(p.getProperty("ld.r2.threshold", String.valueOf(ldR2Threshold)));
        ldParallelJobs     = Integer.parseInt(p.getProperty("ld.parallel.jobs", String.valueOf(ldParallelJobs)));
        outputDir          = p.getProperty("output.dir",           outputDir);
        threads            = Integer.parseInt(p.getProperty("threads", String.valueOf(threads)));
        maxSnpsPerLocus    = Integer.parseInt(p.getProperty("max.snps.per.locus", String.valueOf(maxSnpsPerLocus)));
        splitLdThreshold   = Double.parseDouble(p.getProperty("split.ld.threshold", String.valueOf(splitLdThreshold)));
        splitMinDistBp     = Long.parseLong(p.getProperty("split.min.distance.bp", String.valueOf(splitMinDistBp)));

        // Auto-enable LD if ref.panel.path is set and ld.enabled not explicitly false
        if (!refPanelPath.isEmpty() && !p.containsKey("ld.enabled")) ldEnabled = true;
    }

    public void validate() {
        if (gwasFile.isEmpty()) throw new IllegalStateException("gwas.file must be specified");
        if (lociFile.isEmpty()) throw new IllegalStateException("loci.file must be specified");
        if (gff3File.isEmpty()) throw new IllegalStateException("gff3.file must be specified");
        if (ldEnabled && refPanelPath.isEmpty()) {
            System.err.println("[WARN] ld.enabled=true but ref.panel.path is not set — disabling LD");
            ldEnabled = false;
        }
    }

    public String plinkSubsetsDir() { return outputDir + "/plink_subsets"; }
    public String ldResultsDir()    { return outputDir + "/ld_results"; }

    /** Serialize config as JSON (powers GET /api/project/{id}/config). */
    public String toJson() {
        StringBuilder j = new StringBuilder();
        j.append('{');
        kv(j, "gwas.file", gwasFile, true);
        kv(j, "loci.file", lociFile, false);
        kv(j, "gff3.file", gff3File, false);
        kv(j, "ref.panel.path", refPanelPath, false);
        kv(j, "ref.panel.population", refPanelPopulation, false);
        kv(j, "col.chr", colChr, false);
        kv(j, "col.pos", colPos, false);
        kv(j, "col.pvalue", colPvalue, false);
        kv(j, "col.rsid", colRsid, false);
        kv(j, "col.varid", colVarid, false);
        kv(j, "col.ea", colEa, false);
        kv(j, "col.nea", colNea, false);
        kv(j, "col.beta", colBeta, false);
        kv(j, "col.or", colOr, false);
        kv(j, "col.se", colSe, false);
        kv(j, "col.n", colN, false);
        kv(j, "col.maf", colMaf, false);
        kv(j, "col.info", colInfo, false);
        kv(j, "top.snp.file", topSnpFile, false);
        j.append(String.format(",\"locus.padding\":%d", locusPadding));
        j.append(String.format(",\"ld.enabled\":%s", ldEnabled));
        j.append(String.format(",\"ld.triangle.boundary\":%d", ldTriangleBoundary));
        j.append(String.format(",\"ld.r2.threshold\":%.4f", ldR2Threshold));
        j.append(String.format(",\"ld.parallel.jobs\":%d", ldParallelJobs));
        j.append(String.format(",\"threads\":%d", threads));
        j.append(String.format(",\"max.snps.per.locus\":%d", maxSnpsPerLocus));
        j.append(String.format(",\"split.ld.threshold\":%.4f", splitLdThreshold));
        j.append(String.format(",\"split.min.distance.bp\":%d", splitMinDistBp));
        j.append('}');
        return j.toString();
    }

    private static void kv(StringBuilder j, String k, String v, boolean first) {
        if (!first) j.append(',');
        j.append('"').append(k).append("\":\"");
        j.append(v.replace("\\", "\\\\").replace("\"", "\\\""));
        j.append('"');
    }

    /**
     * Write config to a properties file in the canonical format.
     * This is the ONE writer shared by both hand-editing and the wizard POST endpoint.
     */
    public void writeProperties(String path) throws IOException {
        try (PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(path)))) {
            pw.println("# LYNXgwas project configuration");
            pw.println("# Generated by LYNXgwas " + PIPELINE_VERSION);
            pw.println();
            pw.println("# Input files");
            pw.println("gwas.file=" + gwasFile);
            pw.println("loci.file=" + lociFile);
            pw.println("gff3.file=" + gff3File);
            pw.println();
            pw.println("# Reference panel");
            pw.println("ref.panel.path=" + refPanelPath);
            pw.println("ref.panel.population=" + refPanelPopulation);
            pw.println();
            pw.println("# GWAS column mapping");
            pw.println("col.chr=" + colChr);
            pw.println("col.pos=" + colPos);
            pw.println("col.pvalue=" + colPvalue);
            pw.println("col.rsid=" + colRsid);
            pw.println("col.varid=" + colVarid);
            pw.println("col.ea=" + colEa);
            pw.println("col.nea=" + colNea);
            if (!colBeta.isEmpty()) pw.println("col.beta=" + colBeta);
            if (!colOr.isEmpty())   pw.println("col.or=" + colOr);
            if (!colSe.isEmpty())   pw.println("col.se=" + colSe);
            if (!colN.isEmpty())    pw.println("col.n=" + colN);
            if (!colMaf.isEmpty())  pw.println("col.maf=" + colMaf);
            if (!colInfo.isEmpty()) pw.println("col.info=" + colInfo);
            pw.println();
            pw.println("# Top SNP override file");
            pw.println("top.snp.file=" + topSnpFile);
            pw.println();
            pw.println("# Locus parameters");
            pw.println("locus.padding=" + locusPadding);
            pw.println();
            pw.println("# LD computation");
            pw.println("ld.enabled=" + ldEnabled);
            pw.println("ld.triangle.boundary=" + ldTriangleBoundary);
            pw.println("ld.r2.threshold=" + ldR2Threshold);
            pw.println("ld.parallel.jobs=" + ldParallelJobs);
            pw.println();
            pw.println("# Locus splitting");
            pw.println("split.ld.threshold=" + splitLdThreshold);
            pw.println("split.min.distance.bp=" + splitMinDistBp);
            pw.println();
            pw.println("# Performance");
            pw.println("threads=" + threads);
            pw.println("max.snps.per.locus=" + maxSnpsPerLocus);
        }
    }

    /**
     * Apply JSON values onto this Config. Keys match the properties file format.
     */
    public void applyJson(String json) {
        gwasFile           = jsonStr(json, "gwas.file",            gwasFile);
        lociFile           = jsonStr(json, "loci.file",            lociFile);
        gff3File           = jsonStr(json, "gff3.file",            gff3File);
        refPanelPath       = jsonStr(json, "ref.panel.path",       refPanelPath);
        refPanelPopulation = jsonStr(json, "ref.panel.population", refPanelPopulation);
        colChr             = jsonStr(json, "col.chr",              colChr);
        colPos             = jsonStr(json, "col.pos",              colPos);
        colPvalue          = jsonStr(json, "col.pvalue",           colPvalue);
        colRsid            = jsonStr(json, "col.rsid",             colRsid);
        colVarid           = jsonStr(json, "col.varid",            colVarid);
        colEa              = jsonStr(json, "col.ea",               colEa);
        colNea             = jsonStr(json, "col.nea",              colNea);
        colBeta            = jsonStr(json, "col.beta",             colBeta);
        colOr              = jsonStr(json, "col.or",               colOr);
        colSe              = jsonStr(json, "col.se",               colSe);
        colN               = jsonStr(json, "col.n",                colN);
        colMaf             = jsonStr(json, "col.maf",              colMaf);
        colInfo            = jsonStr(json, "col.info",             colInfo);
        topSnpFile         = jsonStr(json, "top.snp.file",         topSnpFile);
        String v;
        v = jsonStr(json, "locus.padding", "");
        if (!v.isEmpty()) locusPadding = Long.parseLong(v);
        v = jsonStr(json, "ld.enabled", "");
        if (!v.isEmpty()) ldEnabled = Boolean.parseBoolean(v);
        v = jsonStr(json, "ld.triangle.boundary", "");
        if (!v.isEmpty()) ldTriangleBoundary = Integer.parseInt(v);
        v = jsonStr(json, "ld.r2.threshold", "");
        if (!v.isEmpty()) ldR2Threshold = Double.parseDouble(v);
        v = jsonStr(json, "ld.parallel.jobs", "");
        if (!v.isEmpty()) ldParallelJobs = Integer.parseInt(v);
        v = jsonStr(json, "threads", "");
        if (!v.isEmpty()) threads = Integer.parseInt(v);
        v = jsonStr(json, "max.snps.per.locus", "");
        if (!v.isEmpty()) maxSnpsPerLocus = Integer.parseInt(v);
        v = jsonStr(json, "split.ld.threshold", "");
        if (!v.isEmpty()) splitLdThreshold = Double.parseDouble(v);
        v = jsonStr(json, "split.min.distance.bp", "");
        if (!v.isEmpty()) splitMinDistBp = Long.parseLong(v);

        if (!refPanelPath.isEmpty() && !ldEnabled) ldEnabled = true;
    }

    private static String jsonStr(String json, String key, String fallback) {
        String search = "\"" + key + "\":";
        int i = json.indexOf(search);
        if (i < 0) return fallback;
        i += search.length();
        while (i < json.length() && json.charAt(i) == ' ') i++;
        if (i >= json.length()) return fallback;
        if (json.charAt(i) == '"') {
            i++;
            StringBuilder sb = new StringBuilder();
            while (i < json.length()) {
                char c = json.charAt(i++);
                if (c == '"') break;
                if (c == '\\' && i < json.length()) { sb.append(json.charAt(i++)); continue; }
                sb.append(c);
            }
            return sb.toString();
        }
        if (json.startsWith("null", i)) return fallback;
        int end = i;
        while (end < json.length() && ",}]".indexOf(json.charAt(end)) < 0) end++;
        return json.substring(i, end).trim();
    }
}
