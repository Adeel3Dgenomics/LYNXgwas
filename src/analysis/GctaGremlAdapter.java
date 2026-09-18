
import java.io.*;
import java.nio.file.Files;
import java.util.*;

/**
 * GCTA-GREML heritability adapter — separate from CojoAdapter, which only wraps GCTA's
 * --cojo-slct conditional analysis. This adapter builds a genetic relationship matrix (GRM)
 * from a genotype panel and runs REML variance-component estimation against a phenotype file.
 *
 * ── Data-availability caveat (read before wiring this into anything) ──────────────────────
 * LYNXgwas's data model is built around GWAS *summary statistics* (per-SNP beta/se/p) plus a
 * PLINK reference/genotype panel subset+matched to a locus (matched_ref.{bed,bim,fam}). It does
 * not currently ingest per-individual PHENOTYPE values for the genotyped samples: grepping the
 * whole src/ tree for "pheno"/"phenotype" (case-insensitive) turns up nothing outside this
 * file — no phenotype file path in Config.java, no column-mapping UI, no per-sample data
 * structure anywhere.
 *
 * GREML fundamentally requires individual-level phenotypes: --reml partitions phenotypic
 * variance using the GRM *plus* an actual measured trait value per individual (a GCTA --pheno
 * file: FID, IID, value). A GWAS summary-statistic p-value or beta is not a substitute and
 * cannot be turned into one — there is no honest way to derive per-individual trait values from
 * per-SNP association statistics. Consequently:
 *   - makeGrm(...) works today: it only needs the genotype panel that already exists.
 *   - reml(...) requires a caller-supplied phenotype file. If none is given (or the given path
 *     doesn't exist), it throws a clear, specific IOException explaining exactly why, rather
 *     than fabricating a phenotype column or returning a placeholder/mocked h² — silently
 *     degrading to a fake number would be a correctness violation.
 * The moment a project has genuine per-individual phenotype data (of whatever form this app
 * eventually decides to ingest it in), reml() is ready to run against it unchanged.
 */
public class GctaGremlAdapter {

    public static class GrmResult {
        public boolean ok;
        public String error;
        public String grmPrefix;
        public String logTail;
    }

    /** One row of GCTA's .hsq output (Source, Variance, SE — SE is absent on some rows, e.g. logL/df/n). */
    public static class HsqResult {
        public String source;
        public double variance = Double.NaN;
        public double se = Double.NaN;
    }

    public static class RemlResult {
        public boolean ok;
        public String error;
        public List<HsqResult> components = new ArrayList<>();
        public double h2 = Double.NaN;
        public double h2Se = Double.NaN;
        public File hsqFile;
        public String manifestPath;
        public String logTail;
    }

    /**
     * Builds the genetic relationship matrix: gcta64 --bfile <ref> --make-grm --out <prefix>.
     * Needs only the genotype panel this app already builds for a locus (matched_ref.bim/.bed/.fam)
     * — no phenotype data required for this step.
     */
    public static GrmResult makeGrm(File matchedDir, File runDir, String gctaBin) throws IOException {
        File gctaFile = GctaBinaryResolver.verify(gctaBin);
        runDir.mkdirs();

        File bimFile = new File(matchedDir, "matched_ref.bim");
        if (!bimFile.exists())
            throw new IOException("matched_ref.bim not found in " + matchedDir.getAbsolutePath()
                + ". Run the base pipeline first to build the matched reference/genotype panel.");

        String bfilePrefix = new File(matchedDir, "matched_ref").getAbsolutePath();
        File grmPrefix = new File(runDir, "gcta_grm");
        List<String> cmd = Arrays.asList(gctaFile.getAbsolutePath(),
            "--bfile", bfilePrefix, "--make-grm", "--out", grmPrefix.getAbsolutePath());

        File logFile = new File(runDir, "gcta_make_grm.log");
        GrmResult result = new GrmResult();
        int exit = runProcess(cmd, runDir, logFile);
        result.logTail = tailOf(logFile, 60);

        File grmBin = new File(grmPrefix.getAbsolutePath() + ".grm.bin");
        if (exit != 0 || !grmBin.exists()) {
            result.ok = false;
            result.error = "GCTA --make-grm exited with code " + exit + ". Check " + logFile.getAbsolutePath()
                + (result.logTail != null ? "\n---\n" + result.logTail : "");
            return result;
        }
        result.grmPrefix = grmPrefix.getAbsolutePath();
        result.ok = true;
        return result;
    }

    /**
     * Runs REML: gcta64 --reml --grm <prefix> --pheno <phenoFile> --out <prefix>, then parses
     * the .hsq output for h² (the V(G)/Vp row) and its SE.
     *
     * @param phenoFile a GCTA-format phenotype file (FID IID VALUE, tab/space-delimited). Must
     *                  be supplied explicitly by the caller — see the class-level comment for
     *                  why this app cannot construct one itself from what it currently ingests.
     */
    public static RemlResult reml(String grmPrefix, File phenoFile, File runDir, String gctaBin) throws IOException {
        File gctaFile = GctaBinaryResolver.verify(gctaBin);
        if (grmPrefix == null || !new File(grmPrefix + ".grm.bin").exists())
            throw new IOException("GRM not found at " + grmPrefix + ".grm.bin — run makeGrm() first.");

        // Honest handling of the missing-input-type problem described in the class comment:
        // never fabricate a phenotype file or a mocked h² value.
        if (phenoFile == null) {
            throw new IOException("GREML requires a phenotype file (FID IID VALUE, GCTA --pheno format) "
                + "with one measured trait value per individual. LYNXgwas does not currently ingest or "
                + "derive per-individual phenotypes from GWAS summary statistics, so none can be built "
                + "automatically here. This adapter is ready to run the moment a genuine phenotype file "
                + "exists for this project/dataset — pass its path explicitly to reml().");
        }
        if (!phenoFile.exists()) {
            throw new IOException("Phenotype file not found: " + phenoFile.getAbsolutePath()
                + ". GREML requires individual-level phenotype data, which LYNXgwas does not currently "
                + "ingest for GWAS-summary-statistics projects — see the class-level comment on "
                + "GctaGremlAdapter for detail.");
        }
        runDir.mkdirs();

        File outPrefix = new File(runDir, "gcta_reml");
        List<String> cmd = Arrays.asList(gctaFile.getAbsolutePath(),
            "--reml", "--grm", grmPrefix, "--pheno", phenoFile.getAbsolutePath(),
            "--out", outPrefix.getAbsolutePath());

        File logFile = new File(runDir, "gcta_reml.log");
        RemlResult result = new RemlResult();
        int exit = runProcess(cmd, runDir, logFile);
        result.logTail = tailOf(logFile, 60);

        File hsqFile = new File(outPrefix.getAbsolutePath() + ".hsq");
        if (exit != 0 || !hsqFile.exists()) {
            result.ok = false;
            result.error = "GCTA --reml exited with code " + exit + ". Check " + logFile.getAbsolutePath()
                + (result.logTail != null ? "\n---\n" + result.logTail : "");
            return result;
        }
        result.hsqFile = hsqFile;
        result.components = parseHsq(hsqFile);
        for (HsqResult c : result.components) {
            if ("V(G)/Vp".equals(c.source)) {
                result.h2 = c.variance;
                result.h2Se = c.se;
            }
        }
        if (Double.isNaN(result.h2)) {
            result.ok = false;
            result.error = "GCTA .hsq output did not contain a V(G)/Vp row — cannot extract h². "
                + "See " + hsqFile.getAbsolutePath();
            return result;
        }
        result.ok = true;

        File manifestFile = new File(runDir, "gcta_reml.manifest.json");
        writeRemlManifest(manifestFile, result.h2, result.h2Se);
        result.manifestPath = manifestFile.getAbsolutePath();

        return result;
    }

    /**
     * Parses GCTA's tab-separated .hsq output: header line "Source Variance SE" followed by
     * rows such as V(G), V(e), Vp, V(G)/Vp (all three columns), and logL/LRT/df/Pval/n (which
     * GCTA writes with only Source+Variance, no SE) — SE is left NaN for those.
     */
    static List<HsqResult> parseHsq(File f) throws IOException {
        List<HsqResult> out = new ArrayList<>();
        try (BufferedReader br = new BufferedReader(new FileReader(f))) {
            br.readLine(); // header: "Source\tVariance\tSE"
            String line;
            while ((line = br.readLine()) != null) {
                if (line.trim().isEmpty()) continue;
                String[] cols = line.split("\t", -1);
                if (cols.length < 2) continue;
                HsqResult r = new HsqResult();
                r.source = cols[0].trim();
                try {
                    r.variance = Double.parseDouble(cols[1].trim());
                } catch (NumberFormatException e) {
                    continue;
                }
                if (cols.length > 2 && !cols[2].trim().isEmpty()) {
                    try { r.se = Double.parseDouble(cols[2].trim()); } catch (NumberFormatException ignored) {}
                }
                out.add(r);
            }
        }
        return out;
    }

    private static void writeRemlManifest(File f, double h2, double h2Se) throws IOException {
        try (PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(f)))) {
            pw.print("{\"schema_version\":\"1.0\",\"method\":\"gcta_greml\",\"method_version\":\"1.0\",");
            pw.print("\"parameters\":{},");
            pw.print("\"columns\":[");
            pw.print("{\"name\":\"h2\",\"type\":\"double\",\"scope\":\"per_dataset\",\"method\":\"gcta_greml\",\"method_version\":\"1.0\"},");
            pw.print("{\"name\":\"h2_se\",\"type\":\"double\",\"scope\":\"per_dataset\",\"method\":\"gcta_greml\",\"method_version\":\"1.0\"}");
            pw.print("],");
            pw.printf("\"h2\":%.6f,\"h2_se\":%.6f,\"created_at\":%d}", h2, h2Se, System.currentTimeMillis());
        }
    }

    private static int runProcess(List<String> cmd, File workDir, File logFile) throws IOException {
        System.out.printf("[GctaGremlAdapter] Running: %s%n", String.join(" ", cmd));
        try {
            ProcessBuilder pb = new ProcessBuilder(cmd).directory(workDir).redirectErrorStream(true);
            Process proc = pb.start();
            try (BufferedReader br = new BufferedReader(new InputStreamReader(proc.getInputStream()));
                 PrintWriter log = new PrintWriter(new BufferedWriter(new FileWriter(logFile)))) {
                String line;
                while ((line = br.readLine()) != null) log.println(line);
            }
            return proc.waitFor();
        } catch (InterruptedException e) {
            Thread.currentThread().interrupt();
            throw new IOException("GCTA process interrupted", e);
        }
    }

    private static String tailOf(File logFile, int maxLines) {
        if (!logFile.exists()) return null;
        try {
            List<String> lines = Files.readAllLines(logFile.toPath());
            int from = Math.max(0, lines.size() - maxLines);
            return String.join("\n", lines.subList(from, lines.size()));
        } catch (IOException e) {
            return null;
        }
    }
}
