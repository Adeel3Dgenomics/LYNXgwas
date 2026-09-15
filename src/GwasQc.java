import java.io.*;
import java.util.*;

/**
 * Genome-wide GWAS inflation/confounding triage: genomic inflation factor (lambda_GC), computed
 * directly from the project's own summary statistics (no external tool needed), plus an optional
 * user-supplied LDSC intercept to derive the attenuation ratio that separates genuine polygenicity
 * from confounding/population stratification.
 *
 * Caveat surfaced to the user: lambda_GC here is computed from every SNP in the file, not an
 * LD-pruned independent subset, so it runs slightly hot versus a properly pruned estimate — a quick
 * triage signal, not a substitute for LDSC's regression-based intercept.
 */
public class GwasQc {

    private static final double MEDIAN_CHI2_DF1 = 0.454936423;
    public static final double LAMBDA_GC_FLAG = 1.1;
    public static final double LDSC_INTERCEPT_FLAG = 1.05;

    public static class Result {
        public boolean ok;
        public String error;
        public int nSnps;
        public double medianChi2;
        public double meanChi2;
        public double lambdaGC;
        public boolean lambdaFlagged;

        public String toJson() {
            StringBuilder sb = new StringBuilder("{");
            sb.append("\"ok\":").append(ok);
            if (error != null) sb.append(",\"error\":\"").append(esc(error)).append('"');
            sb.append(",\"n_snps\":").append(nSnps);
            sb.append(",\"median_chi2\":").append(medianChi2);
            sb.append(",\"mean_chi2\":").append(meanChi2);
            sb.append(",\"lambda_gc\":").append(lambdaGC);
            sb.append(",\"lambda_flagged\":").append(lambdaFlagged);
            sb.append('}');
            return sb.toString();
        }

        private static String esc(String s) { return s == null ? "" : s.replace("\\", "\\\\").replace("\"", "\\\""); }
    }

    public static Result compute(Config config) {
        Result r = new Result();
        File f = new File(config.gwasFile);
        if (!f.isFile()) { r.error = "GWAS file not found: " + config.gwasFile; return r; }

        List<Double> chi2s = new ArrayList<>();
        double sum = 0;
        try (BufferedReader br = new BufferedReader(new FileReader(f), 1 << 20)) {
            String headerLine = br.readLine();
            if (headerLine == null) { r.error = "Empty GWAS file"; return r; }
            String[] header = headerLine.trim().split("\t");
            int iPval = colIdx(header, config.colPvalue);
            if (iPval < 0) { r.error = "P-value column not found: " + config.colPvalue; return r; }

            String line;
            while ((line = br.readLine()) != null) {
                if (line.isEmpty()) continue;
                String[] row = line.split("\t", -1);
                if (iPval >= row.length) continue;
                double p;
                try { p = Double.parseDouble(row[iPval].trim()); } catch (NumberFormatException e) { continue; }
                if (!(p > 0 && p < 1)) continue;
                double z = StatsUtil.qnorm(1.0 - p / 2.0);
                if (!Double.isFinite(z)) continue;
                double chi2 = z * z;
                chi2s.add(chi2);
                sum += chi2;
            }
        } catch (IOException e) {
            r.error = "Failed to read GWAS file: " + e.getMessage();
            return r;
        }

        if (chi2s.isEmpty()) { r.error = "No usable p-values found"; return r; }
        Collections.sort(chi2s);
        int n = chi2s.size();
        double median = (n % 2 == 1) ? chi2s.get(n / 2) : (chi2s.get(n/2 - 1) + chi2s.get(n/2)) / 2.0;

        r.nSnps = n;
        r.medianChi2 = median;
        r.meanChi2 = sum / n;
        r.lambdaGC = median / MEDIAN_CHI2_DF1;
        r.lambdaFlagged = r.lambdaGC > LAMBDA_GC_FLAG;
        r.ok = true;
        return r;
    }

    /** (intercept - 1) / (mean_chi2 - 1) — near 0 means genuine polygenicity, well above 0 means confounding. */
    public static double attenuationRatio(double ldscIntercept, double meanChi2) {
        double denom = meanChi2 - 1.0;
        if (denom == 0) return Double.NaN;
        return (ldscIntercept - 1.0) / denom;
    }

    private static int colIdx(String[] header, String name) {
        for (int i = 0; i < header.length; i++) if (header[i].trim().equalsIgnoreCase(name)) return i;
        return -1;
    }
}
