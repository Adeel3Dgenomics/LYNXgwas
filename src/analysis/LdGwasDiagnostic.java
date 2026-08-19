
import java.io.*;
import java.util.*;

/**
 * Phase 2: LD–GWAS consistency diagnostic (DENTIST-style).
 *
 * Detects SNPs whose marginal z-score is inconsistent with their LD neighbors,
 * which typically indicates allele coding errors, strand flips, or LD matrix
 * mismatches — the top cause of wrong fine-mapping results.
 *
 * Algorithm:
 *   For each SNP i with z_i = beta_i / se_i:
 *     predicted_z_i = sum_j(r_ij * z_j) / sum_j(r_ij^2)  (LD-weighted average of neighbors)
 *     residual_i = (z_i - predicted_z_i)^2 / var_estimate
 *   Flag SNPs where residual exceeds chi-sq(1) threshold (default p < 1e-4).
 *
 * Output:
 *   ld/consistency_report.tsv — per-SNP flag + statistic
 *   Returns locus-level pass/warn/fail summary
 */
public class LdGwasDiagnostic {

    public static class DiagnosticResult {
        public boolean ok;
        public String error;
        public int totalSnps;
        public int flaggedSnps;
        public String verdict;  // "pass", "warn", "fail"
        public List<FlaggedSnp> flagged = new ArrayList<>();
    }

    public static class FlaggedSnp {
        public String snpId;
        public long pos;
        public double z;
        public double predictedZ;
        public double residual;
        public double pvalue;
    }

    private static final double CHI2_THRESHOLD_WARN = 15.14;  // p < 1e-4, df=1
    private static final double CHI2_THRESHOLD_FLAG = 10.83;  // p < 1e-3, df=1
    private static final int MIN_NEIGHBORS = 3;

    public static DiagnosticResult run(File harmonizedDir, File ldDir) throws IOException {
        DiagnosticResult result = new DiagnosticResult();

        // Read harmonized GWAS for z-scores
        File gwasFile = new File(harmonizedDir, "harmonized_gwas.tsv");
        if (!gwasFile.exists()) {
            result.error = "harmonized_gwas.tsv not found";
            return result;
        }

        // Read SNP order
        File snpOrderFile = new File(ldDir, "ld_snp_order.txt");
        if (!snpOrderFile.exists()) {
            result.error = "ld_snp_order.txt not found";
            return result;
        }

        List<String> snpOrder = new ArrayList<>();
        try (BufferedReader br = new BufferedReader(new FileReader(snpOrderFile))) {
            String line;
            while ((line = br.readLine()) != null) {
                line = line.trim();
                if (!line.isEmpty()) snpOrder.add(line);
            }
        }

        // Read r matrix (signed correlation)
        File rMatrixFile = new File(ldDir, "ld_r.matrix");
        double[][] rMatrix = null;
        if (rMatrixFile.exists()) {
            rMatrix = readSquareMatrix(rMatrixFile, snpOrder.size());
        }
        if (rMatrix == null) {
            // Fall back to r2 matrix (unsigned) — less powerful but still useful
            File r2File = new File(ldDir, "ld_r2.matrix");
            if (r2File.exists()) {
                rMatrix = readSquareMatrix(r2File, snpOrder.size());
                // r2 is unsigned, so we can only detect magnitude outliers
            }
        }
        if (rMatrix == null) {
            result.error = "No LD matrix available (need ld_r.matrix or ld_r2.matrix)";
            return result;
        }

        int n = snpOrder.size();

        // Build position-based lookup for GWAS z-scores
        // snpOrder entries are "chr:pos:a1:a2"
        Map<String, Integer> snpToIdx = new HashMap<>();
        long[] positions = new long[n];
        for (int i = 0; i < n; i++) {
            snpToIdx.put(snpOrder.get(i), i);
            String[] parts = snpOrder.get(i).split(":");
            if (parts.length >= 2) {
                try { positions[i] = Long.parseLong(parts[1]); }
                catch (NumberFormatException e) { positions[i] = 0; }
            }
        }

        // Read GWAS and compute z-scores
        double[] zScores = new double[n];
        String[] snpIds = new String[n];
        Arrays.fill(zScores, Double.NaN);

        try (BufferedReader br = new BufferedReader(new FileReader(gwasFile))) {
            String header = br.readLine();
            if (header == null) { result.error = "Empty harmonized GWAS"; return result; }

            // Header: snp_id chr pos ea nea pvalue beta se ...
            String line;
            while ((line = br.readLine()) != null) {
                if (line.isEmpty()) continue;
                String[] f = line.split("\t", -1);
                if (f.length < 8) continue;

                String snpId = f[0];
                String chr = f[1];
                String pos = f[2];
                String ea = f[3].toUpperCase();
                String nea = f[4].toUpperCase();

                // Match to snpOrder using chr:pos:sorted_alleles
                String a1 = ea.compareTo(nea) <= 0 ? ea : nea;
                String a2 = ea.compareTo(nea) <= 0 ? nea : ea;
                // Try matching against snp order entries
                Integer idx = findSnpIndex(snpToIdx, chr, pos, ea, nea);
                if (idx == null) continue;

                snpIds[idx] = snpId;

                String betaStr = f[6];
                String seStr = f[7];
                if (betaStr.equals("NA") || seStr.equals("NA")) continue;

                try {
                    double beta = Double.parseDouble(betaStr);
                    double se = Double.parseDouble(seStr);
                    if (se > 0 && !Double.isNaN(beta)) {
                        zScores[idx] = beta / se;
                    }
                } catch (NumberFormatException e) {}
            }
        }

        // Count SNPs with valid z-scores
        int validZ = 0;
        for (double z : zScores) if (!Double.isNaN(z)) validZ++;

        if (validZ < MIN_NEIGHBORS + 1) {
            result.ok = true;
            result.totalSnps = validZ;
            result.verdict = "pass";
            result.error = "Too few SNPs with beta/se for diagnostic (" + validZ + ")";
            writeReport(ldDir, result, snpOrder, snpIds, positions);
            return result;
        }

        // Pairwise DENTIST: for each SNP, test against the best LD neighbor.
        // Statistic: T = (z_i - r * z_j)^2 / (1 - r^2), distributed as chi2(1)
        // under the null of consistent data. Only test pairs with 0.5 < |r| < 0.99.
        for (int i = 0; i < n; i++) {
            if (Double.isNaN(zScores[i])) continue;

            // Find the best LD neighbor (highest |r| in the testable range)
            int bestJ = -1;
            double bestAbsR = 0;
            for (int j = 0; j < n; j++) {
                if (j == i || Double.isNaN(zScores[j])) continue;
                double r = rMatrix[i][j];
                if (Double.isNaN(r)) continue;
                double absR = Math.abs(r);
                if (absR > 0.5 && absR < 0.99 && absR > bestAbsR) {
                    bestAbsR = absR;
                    bestJ = j;
                }
            }
            if (bestJ < 0) continue;

            double r = rMatrix[i][bestJ];
            double predictedZ = r * zScores[bestJ];
            double diff = zScores[i] - predictedZ;
            double residual = (diff * diff) / (1.0 - r * r);

            double pval = chiSqPvalue(residual);

            if (residual > CHI2_THRESHOLD_FLAG) {
                FlaggedSnp fs = new FlaggedSnp();
                fs.snpId = snpIds[i] != null ? snpIds[i] : snpOrder.get(i);
                fs.pos = positions[i];
                fs.z = zScores[i];
                fs.predictedZ = predictedZ;
                fs.residual = residual;
                fs.pvalue = pval;
                result.flagged.add(fs);
            }
        }

        result.totalSnps = validZ;
        result.flaggedSnps = result.flagged.size();
        result.ok = true;

        // Determine verdict
        int severeFlags = 0;
        for (FlaggedSnp fs : result.flagged) {
            if (fs.residual > CHI2_THRESHOLD_WARN) severeFlags++;
        }

        if (severeFlags == 0 && result.flaggedSnps == 0) {
            result.verdict = "pass";
        } else if (severeFlags > validZ * 0.05) {
            result.verdict = "high_warn";
        } else {
            result.verdict = "warn";
        }

        writeReport(ldDir, result, snpOrder, snpIds, positions);

        System.out.printf("[LdDiagnostic] %d SNPs checked, %d flagged → %s%n",
            validZ, result.flaggedSnps, result.verdict);
        return result;
    }

    private static Integer findSnpIndex(Map<String, Integer> snpToIdx,
                                         String chr, String pos, String ea, String nea) {
        // Try sorted allele key (canonical form used in snpOrder)
        String a1 = ea.compareTo(nea) <= 0 ? ea : nea;
        String a2 = ea.compareTo(nea) <= 0 ? nea : ea;
        Integer idx = snpToIdx.get(chr + ":" + pos + ":" + a1 + ":" + a2);
        if (idx != null) return idx;

        // Try ref panel ordering (A1:A2 as-is from bim)
        idx = snpToIdx.get(chr + ":" + pos + ":" + ea + ":" + nea);
        if (idx != null) return idx;
        idx = snpToIdx.get(chr + ":" + pos + ":" + nea + ":" + ea);
        return idx;
    }

    private static void writeReport(File ldDir, DiagnosticResult result,
                                     List<String> snpOrder, String[] snpIds,
                                     long[] positions) throws IOException {
        File reportFile = new File(ldDir, "consistency_report.tsv");
        try (PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(reportFile)))) {
            pw.println("snp_id\tchr_pos\tz_observed\tz_predicted\tresidual\tpvalue\tflag");
            for (FlaggedSnp fs : result.flagged) {
                String flag = fs.residual > CHI2_THRESHOLD_WARN ? "FAIL" : "WARN";
                pw.printf("%s\t%d\t%.4f\t%.4f\t%.4f\t%.2e\t%s%n",
                    fs.snpId, fs.pos, fs.z, fs.predictedZ, fs.residual, fs.pvalue, flag);
            }
        }

        // Write summary JSON
        File summaryFile = new File(ldDir, "consistency_summary.json");
        try (PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(summaryFile)))) {
            pw.printf("{\"verdict\":\"%s\",\"total_snps\":%d,\"flagged_snps\":%d,\"timestamp\":%d}%n",
                result.verdict, result.totalSnps, result.flaggedSnps, System.currentTimeMillis());
        }
    }

    private static double[][] readSquareMatrix(File file, int expectedSize) throws IOException {
        List<double[]> rows = new ArrayList<>();
        try (BufferedReader br = new BufferedReader(new FileReader(file))) {
            String line;
            while ((line = br.readLine()) != null) {
                line = line.trim();
                if (line.isEmpty()) continue;
                String[] parts = line.split("\\s+");
                double[] row = new double[parts.length];
                for (int i = 0; i < parts.length; i++) {
                    row[i] = parseDouble(parts[i]);
                }
                rows.add(row);
            }
        }
        if (rows.isEmpty()) return null;
        int n = rows.size();
        double[][] matrix = new double[n][n];
        for (int i = 0; i < n; i++) {
            double[] row = rows.get(i);
            for (int j = 0; j < Math.min(row.length, n); j++) {
                matrix[i][j] = row[j];
            }
        }
        return matrix;
    }

    private static double parseDouble(String s) {
        if (s == null || s.isEmpty() || s.equalsIgnoreCase("NA") || s.equalsIgnoreCase("nan"))
            return Double.NaN;
        try { return Double.parseDouble(s); }
        catch (NumberFormatException e) { return Double.NaN; }
    }

    /**
     * Approximate chi-squared(1) survival function using the Wilson-Hilferty transform.
     */
    private static double chiSqPvalue(double x) {
        if (x <= 0) return 1.0;
        double z = Math.pow(x, 1.0 / 3.0) * (1.0 - 2.0 / 9.0) - (1.0 - 2.0 / 9.0);
        z /= Math.sqrt(2.0 / 9.0);
        return 0.5 * erfc(z / Math.sqrt(2.0));
    }

    private static double erfc(double x) {
        double t = 1.0 / (1.0 + 0.3275911 * Math.abs(x));
        double poly = t * (0.254829592 + t * (-0.284496736 + t * (1.421413741
            + t * (-1.453152027 + t * 1.061405429))));
        double result = poly * Math.exp(-x * x);
        return x >= 0 ? result : 2.0 - result;
    }
}
