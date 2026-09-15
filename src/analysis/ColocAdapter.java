
import java.io.*;
import java.util.*;

/**
 * Colocalisation (coloc.abf) adapter.
 * Joins this locus's own harmonized GWAS (trait 1) against an externally supplied
 * second trait's summary statistics file (trait 2) by chr:pos, harmonizes alleles,
 * and generates a self-contained R script that runs coloc::coloc.abf.
 *
 * Unlike SusieAdapter/FinemapAdapter this adapter does NOT write result.tsv directly —
 * it writes a raw output file (coloc_raw.tsv) with columns matching coloc.yaml's
 * output_mapping (snp_id, SNP.PP.H4, SNP.PP.H3, SNP.PP.H0) and lets PluginEngine's
 * generic mapOutput() do the column mapping, same as any declarative tool.
 */
public class ColocAdapter {

    public static void prepareRun(File harmonizedDir, File runDir, Locus locus, Config config,
                                   String trait2File, String trait2Type, int trait2N, int trait2NCases,
                                   double p1, double p2, double p12,
                                   File analysisRoot, boolean restrictToFinemapped, double pipThreshold) throws IOException {
        if (trait2File == null || trait2File.trim().isEmpty())
            throw new IOException("Trait 2 summary statistics file is required for colocalisation.");
        File t2 = new File(trait2File.trim());
        if (!t2.isFile())
            throw new IOException("Trait 2 file not found: " + t2.getAbsolutePath());
        if (config.sampleN <= 0)
            throw new IOException("Sample size (N) is required for coloc (trait 1). Set it in the project configuration.");
        runDir.mkdirs();

        // ── Load trait 1 (own harmonized GWAS for this locus) ──
        // Column order fixed by LocusGwasExtractor: snp_id chr pos ea nea pvalue beta se or n maf info rsid varid
        Map<String, double[]> t1 = new LinkedHashMap<>();   // chr:pos -> {beta, se, maf}
        Map<String, String[]> t1Alleles = new LinkedHashMap<>(); // chr:pos -> {ea, nea}
        Map<String, String> t1Snp = new LinkedHashMap<>();  // chr:pos -> snp_id
        File gwas1 = new File(harmonizedDir, "harmonized_gwas.tsv");
        if (!gwas1.exists()) throw new IOException("harmonized_gwas.tsv not found. Run the base pipeline first.");
        try (BufferedReader br = new BufferedReader(new FileReader(gwas1))) {
            br.readLine();
            String line;
            while ((line = br.readLine()) != null) {
                if (line.isEmpty()) continue;
                String[] f = line.split("\t", -1);
                if (f.length < 8) continue;
                String key = f[1] + ":" + f[2];
                try {
                    double beta = Double.parseDouble(f[6]);
                    double se = Double.parseDouble(f[7]);
                    if (!Double.isFinite(beta) || !Double.isFinite(se) || se <= 0) continue;
                    double maf = (f.length > 10) ? parseOrNaN(f[10]) : Double.NaN;
                    t1.put(key, new double[]{beta, se, maf});
                    t1Alleles.put(key, new String[]{f[3].toUpperCase(), f[4].toUpperCase()});
                    t1Snp.put(key, f[0]);
                } catch (NumberFormatException ignored) {}
            }
        }
        if (t1.isEmpty()) throw new IOException("No usable trait-1 SNPs (beta/se) in harmonized_gwas.tsv");

        Set<String> allowedSnpIds = null;
        if (restrictToFinemapped) {
            allowedSnpIds = findFinemappedSnpIds(analysisRoot, pipThreshold);
            if (allowedSnpIds == null) {
                System.out.println("[ColocAdapter] restrict_to_finemapped requested but no prior SuSiE/FINEMAP/COJO "
                    + "run found for this locus — falling back to the whole-locus SNP set.");
            } else if (allowedSnpIds.isEmpty()) {
                throw new IOException("restrict_to_finemapped found a prior fine-mapping run, but no SNP passed "
                    + "its threshold — nothing to test. Lower the PIP threshold or run without restriction.");
            }
        }

        // ── Auto-detect and load trait 2 file ──
        try (BufferedReader br = new BufferedReader(new FileReader(t2))) {
            String header = br.readLine();
            if (header == null) throw new IOException("Trait 2 file is empty: " + t2.getAbsolutePath());
            String delim = header.contains("\t") ? "\t" : (header.contains(",") ? "," : "\\s+");
            String[] cols = header.split(delim, -1);
            int chrIdx = findCol(cols, "chr", "chrom", "chromosome");
            int posIdx = findCol(cols, "pos", "bp", "position", "base_pair_location");
            int eaIdx  = findCol(cols, "ea", "a1", "effect_allele", "allele1");
            int neaIdx = findCol(cols, "nea", "a2", "other_allele", "allele2", "non_effect_allele");
            int betaIdx = findCol(cols, "beta", "b", "effect");
            int seIdx  = findCol(cols, "se", "stderr", "standard_error");
            int mafIdx = findCol(cols, "maf", "eaf", "freq", "af", "effect_allele_frequency");
            if (chrIdx < 0 || posIdx < 0 || eaIdx < 0 || neaIdx < 0 || betaIdx < 0 || seIdx < 0) {
                throw new IOException("Could not auto-detect chr/pos/ea/nea/beta/se columns in trait 2 file header: " + header);
            }

            File inputTsv = new File(runDir, "coloc_input.tsv");
            int matched = 0, fallbackMaf = 0;
            try (PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(inputTsv)))) {
                pw.println("snp_id\tchr\tpos\tbeta1\tse1\tmaf1\tbeta2\tse2\tmaf2");
                String line;
                while ((line = br.readLine()) != null) {
                    if (line.trim().isEmpty()) continue;
                    String[] f = line.split(delim, -1);
                    int maxIdx = max(chrIdx, posIdx, eaIdx, neaIdx, betaIdx, seIdx, mafIdx);
                    if (f.length <= maxIdx) continue;
                    String chr = f[chrIdx].replaceFirst("^chr", "").trim();
                    String pos = f[posIdx].trim();
                    String key = chr + ":" + pos;
                    double[] row1 = t1.get(key);
                    if (row1 == null) continue;
                    if (allowedSnpIds != null && !allowedSnpIds.contains(t1Snp.get(key))) continue;

                    double beta2, se2;
                    try {
                        beta2 = Double.parseDouble(f[betaIdx]);
                        se2 = Double.parseDouble(f[seIdx]);
                        if (!Double.isFinite(beta2) || !Double.isFinite(se2) || se2 <= 0) continue;
                    } catch (NumberFormatException e) { continue; }

                    // Harmonize trait2 alleles against trait1's effect/other allele
                    String ea2 = f[eaIdx].toUpperCase().trim();
                    String nea2 = f[neaIdx].toUpperCase().trim();
                    String[] alleles1 = t1Alleles.get(key);
                    if (alleles1[0].equals(nea2) && alleles1[1].equals(ea2)) {
                        beta2 = -beta2; // swapped relative to trait 1 -> flip sign
                    } else if (!(alleles1[0].equals(ea2) && alleles1[1].equals(nea2))) {
                        continue; // allele mismatch (strand issue / different variant) — skip rather than guess
                    }

                    double maf2 = mafIdx >= 0 ? parseOrNaN(f[mafIdx]) : Double.NaN;
                    if (Double.isNaN(maf2)) { maf2 = row1[2]; fallbackMaf++; }

                    pw.printf("%s\t%s\t%s\t%.8g\t%.8g\t%s\t%.8g\t%.8g\t%s%n",
                        t1Snp.get(key), chr, pos, row1[0], row1[1],
                        Double.isNaN(row1[2]) ? "NA" : String.format("%.4f", row1[2]),
                        beta2, se2, Double.isNaN(maf2) ? "NA" : String.format("%.4f", maf2));
                    matched++;
                }
            }
            if (matched == 0)
                throw new IOException("No overlapping, allele-consistent SNPs between trait 1 and trait 2 in this locus.");
            System.out.printf("[ColocAdapter] Matched %d SNPs (%d trait-2 MAF fallbacks to ref panel/trait-1)%n",
                matched, fallbackMaf);

            writeRScript(runDir, inputTsv, config, trait2Type, trait2N, trait2NCases, p1, p2, p12);
        }
    }

    private static void writeRScript(File runDir, File inputTsv, Config config,
                                      String trait2Type, int trait2N, int trait2NCases,
                                      double p1, double p2, double p12) throws IOException {
        boolean t1Binary = "binary".equalsIgnoreCase(config.traitType);
        boolean t2Binary = "cc".equalsIgnoreCase(trait2Type);
        double t1S = (t1Binary && config.sampleN > 0) ? (double) config.nCases / config.sampleN : Double.NaN;
        double t2S = (t2Binary && trait2N > 0) ? (double) trait2NCases / trait2N : Double.NaN;

        String inputPath = inputTsv.getAbsolutePath().replace("\\", "/");
        String rawOut = new File(runDir, "coloc_raw.tsv").getAbsolutePath().replace("\\", "/");
        String summaryOut = new File(runDir, "coloc_summary.json").getAbsolutePath().replace("\\", "/");

        File rScript = new File(runDir, "coloc_run.R");
        try (PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(rScript)))) {
            pw.println("#!/usr/bin/env Rscript");
            pw.println("suppressPackageStartupMessages({ library(coloc); library(data.table) })");
            pw.println();
            pw.printf("df <- fread('%s', sep='\\t', data.table=FALSE)%n", inputPath);
            pw.println("df$beta1 <- as.numeric(df$beta1); df$se1 <- as.numeric(df$se1)");
            pw.println("df$beta2 <- as.numeric(df$beta2); df$se2 <- as.numeric(df$se2)");
            pw.println("df$maf1 <- suppressWarnings(as.numeric(df$maf1))");
            pw.println("df$maf2 <- suppressWarnings(as.numeric(df$maf2))");
            pw.printf("n1 <- %d%n", config.sampleN);
            pw.printf("n2 <- %d%n", trait2N);
            pw.println();

            pw.println("d1 <- list(snp=df$snp_id, beta=df$beta1, varbeta=df$se1^2, position=df$pos, N=n1,");
            pw.printf("           type='%s')%n", t1Binary ? "cc" : "quant");
            if (t1Binary && !Double.isNaN(t1S)) pw.printf("d1$s <- %.6f%n", t1S);
            pw.println("if (all(!is.na(df$maf1))) d1$MAF <- df$maf1");
            pw.println();
            pw.println("d2 <- list(snp=df$snp_id, beta=df$beta2, varbeta=df$se2^2, position=df$pos, N=n2,");
            pw.printf("           type='%s')%n", t2Binary ? "cc" : "quant");
            if (t2Binary && !Double.isNaN(t2S)) pw.printf("d2$s <- %.6f%n", t2S);
            pw.println("if (all(!is.na(df$maf2))) d2$MAF <- df$maf2");
            pw.println();

            pw.printf("res <- coloc.abf(dataset1=d1, dataset2=d2, p1=%s, p2=%s, p12=%s)%n",
                fmt(p1), fmt(p2), fmt(p12));
            pw.println();

            pw.println("per_snp <- res$results");
            pw.println("out <- data.frame(snp_id=per_snp$snp,");
            pw.println("  'SNP.PP.H4'=round(per_snp$SNP.PP.H4, 6),");
            pw.println("  'SNP.PP.H3'=round(per_snp$SNP.PP.H3, 6),");
            pw.println("  'SNP.PP.H0'=round(per_snp$SNP.PP.H0, 6),");
            pw.println("  check.names=FALSE)");
            pw.printf("write.table(out, '%s', sep='\\t', quote=FALSE, row.names=FALSE)%n", rawOut);
            pw.println();

            pw.println("pp <- res$summary");
            // Standard PP.H4 interpretation scale (one shared causal variant between the two traits):
            // >=0.75 strong, >=0.50 moderate (common literature threshold), 0.10-0.50 weak/exploratory, <0.10 none.
            pw.println("h4 <- pp['PP.H4.abf']");
            pw.println("interp <- if (h4 >= 0.75) 'strong' else if (h4 >= 0.50) 'moderate' else if (h4 >= 0.10) 'weak' else 'none'");
            pw.printf("writeLines(sprintf('{\"nsnps\":%%d,\"PP.H0\":%%.6f,\"PP.H1\":%%.6f,\"PP.H2\":%%.6f,\"PP.H3\":%%.6f,\"PP.H4\":%%.6f,\"PP.H4.interpretation\":\"%%s\"}',%n");
            pw.println("  pp['nsnps'], pp['PP.H0.abf'], pp['PP.H1.abf'], pp['PP.H2.abf'], pp['PP.H3.abf'], pp['PP.H4.abf'], interp),");
            pw.printf("  '%s')%n", summaryOut);
            pw.println();
            pw.println("cat(sprintf('coloc: PP.H4 (shared causal variant) = %.4f [%s evidence]\\n', h4, interp))");
            pw.println("cat('Scale: >=0.75 strong, >=0.50 moderate, 0.10-0.50 weak, <0.10 none\\n')");
        }
    }

    private static String fmt(double d) {
        return String.valueOf(d);
    }

    /**
     * Finds SNP ids passing a prior SuSiE/FINEMAP/COJO run's threshold for this locus, by scanning each
     * subdirectory of analysisRoot's "runs" directory for a provenance.json, picking the most recently
     * modified run whose tool is one of the three, and reading its result.tsv. Returns null if no prior
     * run of any of the three tools exists (caller falls back to no restriction), or an empty (non-null)
     * set if a run exists but nothing passed the threshold.
     */
    private static Set<String> findFinemappedSnpIds(File analysisRoot, double pipThreshold) {
        if (analysisRoot == null) return null;
        File runsDir = new File(analysisRoot, "runs");
        File[] runDirs = runsDir.listFiles(File::isDirectory);
        if (runDirs == null) return null;

        File bestRun = null;
        String bestTool = null;
        long bestTime = -1;
        for (File rd : runDirs) {
            File prov = new File(rd, "provenance.json");
            if (!prov.exists()) continue;
            String tool;
            try {
                String json = new String(java.nio.file.Files.readAllBytes(prov.toPath()), "UTF-8");
                int i = json.indexOf("\"tool\":");
                if (i < 0) continue;
                int start = json.indexOf('"', i + 7) + 1;
                int end = json.indexOf('"', start);
                tool = json.substring(start, end);
            } catch (IOException e) { continue; }
            if (!tool.equals("susie_finemapping") && !tool.equals("finemap") && !tool.equals("cojo_conditional")) continue;
            File result = new File(rd, "result.tsv");
            if (!result.exists() || rd.lastModified() <= bestTime) continue;
            bestRun = rd; bestTool = tool; bestTime = rd.lastModified();
        }
        if (bestRun == null) return null;

        Set<String> ids = new LinkedHashSet<>();
        File result = new File(bestRun, "result.tsv");
        try (BufferedReader br = new BufferedReader(new FileReader(result))) {
            String header = br.readLine();
            if (header == null) return ids;
            String[] cols = header.split("\t", -1);
            int idIdx = 0; // all three tools key their result.tsv on the first column
            int valIdx = -1;
            boolean isCojo = bestTool.equals("cojo_conditional");
            String pipCol = bestTool.equals("finemap") ? "finemap_pip" : "susie_pip";
            for (int i = 0; i < cols.length; i++) {
                if (cols[i].equals("snp_id")) idIdx = i;
                if (cols[i].equals(isCojo ? "cojo_selected" : pipCol)) valIdx = i;
            }
            if (valIdx < 0) return ids;
            String line;
            while ((line = br.readLine()) != null) {
                if (line.trim().isEmpty()) continue;
                String[] f = line.split("\t", -1);
                if (f.length <= Math.max(idIdx, valIdx)) continue;
                boolean qualifies = isCojo
                    ? f[valIdx].trim().equalsIgnoreCase("selected")
                    : parseDoubleOrZero(f[valIdx]) >= pipThreshold;
                if (qualifies) ids.add(f[idIdx].trim());
            }
        } catch (IOException e) {
            return new LinkedHashSet<>();
        }
        System.out.printf("[ColocAdapter] Restricting to %d SNP(s) from prior %s run (%s)%n",
            ids.size(), bestTool, bestRun.getName());
        return ids;
    }

    private static double parseDoubleOrZero(String s) {
        try { return Double.parseDouble(s.trim()); } catch (Exception e) { return 0; }
    }

    private static double parseOrNaN(String s) {
        if (s == null) return Double.NaN;
        s = s.trim();
        if (s.isEmpty() || s.equalsIgnoreCase("NA") || s.equals(".")) return Double.NaN;
        try {
            double v = Double.parseDouble(s);
            return (v > 0 && v < 1) ? v : Double.NaN;
        } catch (NumberFormatException e) { return Double.NaN; }
    }

    private static int findCol(String[] cols, String... aliases) {
        for (int i = 0; i < cols.length; i++) {
            String c = cols[i].trim().toLowerCase();
            for (String a : aliases) if (c.equals(a)) return i;
        }
        return -1;
    }

    private static int max(int... vals) {
        int m = Integer.MIN_VALUE;
        for (int v : vals) if (v > m) m = v;
        return m;
    }
}
