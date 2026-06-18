import java.io.*;
import java.util.*;

public class Config {
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
}
