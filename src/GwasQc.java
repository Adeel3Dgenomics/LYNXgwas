import java.io.*;
import java.util.*;

/**
 * Genome-wide GWAS inflation/confounding triage: genomic inflation factor (lambda_GC), computed
 * directly from the project's own summary statistics (no external tool needed), plus an optional
 * user-supplied LDSC intercept to derive the attenuation ratio that separates genuine polygenicity
 * from confounding/population stratification.
 *
 * Before computing lambda_GC, SNPs are thinned to an approximately-independent set via simple
 * distance-based pruning (at most one SNP kept per PRUNE_DISTANCE_BP window per chromosome) — this
 * project deliberately requires no reference panel or external tool for this quick check, so it
 * cannot run genotype-based LD pruning (PLINK --indep-pairwise); distance-based thinning is a
 * lighter proxy that still substantially reduces the effect of dense LD blocks contributing many
 * highly-correlated chi-square values to the median. It is still not equivalent to true LD-based
 * pruning or LDSC's regression intercept — this remains a fast triage signal, not a replacement.
 */
public class GwasQc {

    private static final double MEDIAN_CHI2_DF1 = 0.454936423;
    public static final double LAMBDA_GC_FLAG = 1.1;
    public static final double LDSC_INTERCEPT_FLAG = 1.05;
    /** Keep at most one SNP per this many bp per chromosome before computing lambda_GC. */
    public static final long PRUNE_DISTANCE_BP = 250_000;

    public static class Result {
        public boolean ok;
        public String error;
        public int nSnps;
        public int nSnpsScanned;
        public double medianChi2;
        public double meanChi2;
        public double lambdaGC;
        public boolean lambdaFlagged;

        public String toJson() {
            StringBuilder sb = new StringBuilder("{");
            sb.append("\"ok\":").append(ok);
            if (error != null) sb.append(",\"error\":\"").append(esc(error)).append('"');
            sb.append(",\"n_snps\":").append(nSnps);
            sb.append(",\"n_snps_scanned\":").append(nSnpsScanned);
            sb.append(",\"median_chi2\":").append(medianChi2);
            sb.append(",\"mean_chi2\":").append(meanChi2);
            sb.append(",\"lambda_gc\":").append(lambdaGC);
            sb.append(",\"lambda_flagged\":").append(lambdaFlagged);
            sb.append('}');
            return sb.toString();
        }

        private static String esc(String s) { return s == null ? "" : s.replace("\\", "\\\\").replace("\"", "\\\""); }
    }

    private static class SnpRec {
        String chr; long pos; double chi2;
        SnpRec(String chr, long pos, double chi2) { this.chr = chr; this.pos = pos; this.chi2 = chi2; }
    }

    public static Result compute(Config config) {
        Result r = new Result();
        File f = new File(config.gwasFile);
        if (!f.isFile()) { r.error = "GWAS file not found: " + config.gwasFile; return r; }

        List<SnpRec> recs = new ArrayList<>();
        try (BufferedReader br = new BufferedReader(new FileReader(f), 1 << 20)) {
            String headerLine = br.readLine();
            if (headerLine == null) { r.error = "Empty GWAS file"; return r; }
            String[] header = headerLine.trim().split("\t");
            int iPval = colIdx(header, config.colPvalue);
            int iChr  = colIdx(header, config.colChr);
            int iPos  = colIdx(header, config.colPos);
            if (iPval < 0) { r.error = "P-value column not found: " + config.colPvalue; return r; }
            // Chr/pos are needed only for distance-based pruning; if either is unmapped, fall back
            // to computing over every SNP (matches this method's pre-pruning behavior) rather than
            // failing the whole QC check.
            boolean canPrune = iChr >= 0 && iPos >= 0;

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

                String chr = ".";
                long pos = 0;
                if (canPrune) {
                    if (iChr >= row.length || iPos >= row.length) continue;
                    try { pos = Long.parseLong(row[iPos].trim()); } catch (NumberFormatException e) { continue; }
                    chr = MultiLocusScanner.normalizeChr(row[iChr].trim());
                }
                recs.add(new SnpRec(chr, pos, chi2));
            }
        } catch (IOException e) {
            r.error = "Failed to read GWAS file: " + e.getMessage();
            return r;
        }

        if (recs.isEmpty()) { r.error = "No usable p-values found"; return r; }
        r.nSnpsScanned = recs.size();

        List<Double> chi2s = prune(recs);

        Collections.sort(chi2s);
        int n = chi2s.size();
        double median = (n % 2 == 1) ? chi2s.get(n / 2) : (chi2s.get(n/2 - 1) + chi2s.get(n/2)) / 2.0;
        double sum = 0;
        for (double c : chi2s) sum += c;

        r.nSnps = n;
        r.medianChi2 = median;
        r.meanChi2 = sum / n;
        r.lambdaGC = median / MEDIAN_CHI2_DF1;
        r.lambdaFlagged = r.lambdaGC > LAMBDA_GC_FLAG;
        r.ok = true;
        return r;
    }

    /** Distance-based approximate independence pruning: sort each chromosome by position and keep
     *  a SNP only if it is at least PRUNE_DISTANCE_BP past the last kept SNP on that chromosome. */
    private static List<Double> prune(List<SnpRec> recs) {
        Map<String, List<SnpRec>> byChr = new HashMap<>();
        for (SnpRec s : recs) byChr.computeIfAbsent(s.chr, k -> new ArrayList<>()).add(s);

        List<Double> kept = new ArrayList<>();
        for (List<SnpRec> chrRecs : byChr.values()) {
            if (chrRecs.get(0).chr.equals(".")) { // chr/pos unavailable — no pruning possible
                for (SnpRec s : chrRecs) kept.add(s.chi2);
                continue;
            }
            chrRecs.sort(Comparator.comparingLong(s -> s.pos));
            long lastKeptPos = Long.MIN_VALUE / 2;
            for (SnpRec s : chrRecs) {
                if (s.pos - lastKeptPos >= PRUNE_DISTANCE_BP) {
                    kept.add(s.chi2);
                    lastKeptPos = s.pos;
                }
            }
        }
        return kept;
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
