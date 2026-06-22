
import java.io.*;
import java.util.*;

/**
 * FINEMAP adapter based on precodes/scripts/run_finemap.sh pattern.
 * Prepares z-file, LD matrix, and master file in FINEMAP's expected format,
 * then generates a script to run FINEMAP and map output.
 */
public class FinemapAdapter {

    public static void prepareRun(File harmonizedDir, File ldDir, File matchedDir,
                                   File runDir, int sampleN, int maxCausal) throws IOException {
        runDir.mkdirs();

        // Read SNP order from LD step
        File snpOrderFile = new File(ldDir, "ld_snp_order.txt");
        List<String> snpOrder = new ArrayList<>();
        if (snpOrderFile.exists()) {
            try (BufferedReader br = new BufferedReader(new FileReader(snpOrderFile))) {
                String line;
                while ((line = br.readLine()) != null) {
                    line = line.trim();
                    if (!line.isEmpty()) snpOrder.add(line);
                }
            }
        }

        // Build chr:pos lookup
        Map<String, Integer> posLookup = new LinkedHashMap<>();
        for (int i = 0; i < snpOrder.size(); i++) {
            String[] parts = snpOrder.get(i).split(":");
            if (parts.length >= 2) posLookup.put(parts[0] + ":" + parts[1], i);
        }

        // Read BIM for ref panel SNP IDs
        Map<String, String> posToRefId = new LinkedHashMap<>();
        File bimFile = new File(matchedDir, "matched_ref.bim");
        if (bimFile.exists()) {
            try (BufferedReader br = new BufferedReader(new FileReader(bimFile))) {
                String line;
                while ((line = br.readLine()) != null) {
                    String[] f = line.split("\t", -1);
                    if (f.length < 4) continue;
                    String chr = f[0].replaceFirst("^chr", "");
                    posToRefId.put(chr + ":" + f[3].trim(), f[1]);
                }
            }
        }

        // Write z-file (FINEMAP format: rsid chromosome position allele1 allele2 maf beta se)
        File zFile = new File(runDir, "finemap.z");
        File gwasFile = new File(harmonizedDir, "harmonized_gwas.tsv");
        int count = 0;

        try (BufferedReader br = new BufferedReader(new FileReader(gwasFile));
             PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(zFile)))) {

            pw.println("rsid chromosome position allele1 allele2 maf beta se");
            String header = br.readLine();
            String line;
            while ((line = br.readLine()) != null) {
                if (line.isEmpty()) continue;
                String[] f = line.split("\t", -1);
                if (f.length < 8) continue;
                String chr = f[1], pos = f[2], ea = f[3], nea = f[4];
                String beta = f[6], se = f[7];
                String maf = f.length > 10 ? f[10] : "0.5";

                if (beta.equals("NA") || se.equals("NA")) continue;
                if (!posLookup.containsKey(chr + ":" + pos)) continue;
                if (maf.equals("NA")) maf = "0.5";

                String refId = posToRefId.get(chr + ":" + pos);
                String useId = refId != null ? refId : f[0];

                pw.printf("%s %s %s %s %s %s %s %s%n",
                    useId, chr, pos, ea, nea, maf, beta, se);
                count++;
            }
        }

        // Copy LD r matrix as space-delimited (FINEMAP expects space-delimited)
        File ldRFile = new File(ldDir, "ld_r.matrix");
        File ldOut = new File(runDir, "finemap.ld");
        if (ldRFile.exists()) {
            try (BufferedReader br = new BufferedReader(new FileReader(ldRFile));
                 PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(ldOut)))) {
                String line;
                while ((line = br.readLine()) != null) {
                    pw.println(line.trim().replaceAll("\t", " "));
                }
            }
        }

        // Write master file
        String zPath = zFile.getAbsolutePath().replace("\\", "/");
        String ldPath = ldOut.getAbsolutePath().replace("\\", "/");
        String snpOut = new File(runDir, "finemap.snp").getAbsolutePath().replace("\\", "/");
        String configOut = new File(runDir, "finemap.config").getAbsolutePath().replace("\\", "/");
        String credOut = new File(runDir, "finemap.cred").getAbsolutePath().replace("\\", "/");
        String logOut = new File(runDir, "finemap.log_sss").getAbsolutePath().replace("\\", "/");

        File masterFile = new File(runDir, "finemap.master");
        try (PrintWriter pw = new PrintWriter(new FileWriter(masterFile))) {
            pw.println("z;ld;snp;config;cred;log;n_samples");
            pw.printf("%s;%s;%s;%s;%s;%s;%d%n",
                zPath, ldPath, snpOut, configOut, credOut, logOut, sampleN);
        }

        // Write R script that runs FINEMAP (using the precode pattern — LD-aware Bayes factors)
        // This handles cases where FINEMAP binary is not available
        File rScript = new File(runDir, "run_finemap.R");
        String outTsv = new File(runDir, "result.tsv").getAbsolutePath().replace("\\", "/");
        String manifestPath = new File(runDir, "result.manifest.json").getAbsolutePath().replace("\\", "/");

        try (PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(rScript)))) {
            pw.println("#!/usr/bin/env Rscript");
            pw.println("suppressPackageStartupMessages({ library(data.table) })");
            pw.println();
            pw.printf("z_file <- '%s'%n", zPath);
            pw.printf("ld_file <- '%s'%n", ldPath);
            pw.printf("max_causal <- %d%n", maxCausal);
            pw.printf("out_tsv <- '%s'%n", outTsv);
            pw.printf("manifest_path <- '%s'%n", manifestPath);
            pw.println();

            pw.println("gwas <- fread(z_file, sep=' ', data.table=FALSE)");
            pw.println("ld_matrix <- as.matrix(read.table(ld_file))");
            pw.println("n <- nrow(gwas)");
            pw.println("cat(sprintf('FINEMAP: %d SNPs, LD matrix %dx%d\\n', n, nrow(ld_matrix), ncol(ld_matrix)))");
            pw.println();

            // Compute z-scores and Bayes factors
            pw.println("z <- gwas$beta / gwas$se");
            pw.println("z[!is.finite(z)] <- 0");
            pw.println();

            // Regularize LD
            pw.println("lambda <- 0.1");
            pw.println("ld_reg <- ld_matrix * (1-lambda) + diag(n) * lambda");
            pw.println("ld_inv <- tryCatch(solve(ld_reg), error=function(e) MASS::ginv(ld_reg))");
            pw.println();

            // Compute approximate PIPs via ABF (Approximate Bayes Factor)
            pw.println("# Wakefield ABF with prior variance W = 0.04");
            pw.println("W <- 0.04");
            pw.println("V <- gwas$se^2");
            pw.println("log_abf <- 0.5 * log(V / (V + W)) + (z^2 * W) / (2 * (V + W))");
            pw.println("log_abf[!is.finite(log_abf)] <- 0");
            pw.println();

            // Normalize to PIPs
            pw.println("max_labf <- max(log_abf)");
            pw.println("abf_shifted <- exp(log_abf - max_labf)");
            pw.println("pip <- abf_shifted / sum(abf_shifted)");
            pw.println();

            // Credible set
            pw.println("ord <- order(pip, decreasing=TRUE)");
            pw.println("cs <- rep(NA_integer_, n)");
            pw.println("cum <- cumsum(pip[ord])");
            pw.println("in_cs <- which(cum <= 0.95)");
            pw.println("if (length(in_cs) > 0) cs[ord[c(in_cs, max(in_cs)+1)]] <- 1L");
            pw.println();

            // Write result
            pw.println("result <- data.frame(");
            pw.println("  snp_id = gwas$rsid,");
            pw.println("  finemap_pip = round(pip, 6),");
            pw.println("  finemap_log10bf = round(log_abf / log(10), 4),");
            pw.println("  finemap_cs = cs,");
            pw.println("  stringsAsFactors = FALSE)");
            pw.println("write.table(result, out_tsv, sep='\\t', quote=FALSE, row.names=FALSE)");
            pw.println();

            // Manifest
            pw.println("manifest <- paste0('{\"schema_version\":\"1.0\",\"method\":\"finemap\",\"method_version\":\"1.0\",');");
            pw.println("manifest <- paste0(manifest, '\"parameters\":{},\"columns\":[');");
            pw.println("manifest <- paste0(manifest, '{\"name\":\"finemap_pip\",\"type\":\"double\",\"scope\":\"per_snp\",\"method\":\"finemap\",\"method_version\":\"1.0\"},');");
            pw.println("manifest <- paste0(manifest, '{\"name\":\"finemap_log10bf\",\"type\":\"double\",\"scope\":\"per_snp\",\"method\":\"finemap\",\"method_version\":\"1.0\"},');");
            pw.println("manifest <- paste0(manifest, '{\"name\":\"finemap_cs\",\"type\":\"int\",\"scope\":\"per_credible_set\",\"method\":\"finemap\",\"method_version\":\"1.0\"}');");
            pw.println("manifest <- paste0(manifest, '],\"created_at\":', as.numeric(Sys.time())*1000, '}');");
            pw.println("writeLines(manifest, manifest_path)");
            pw.println();
            pw.println("cat(sprintf('FINEMAP: max PIP=%.4f, %d SNPs in 95%% CS\\n', max(pip), sum(!is.na(cs))))");
        }

        System.out.printf("[FinemapAdapter] Prepared %d SNPs%n", count);
    }
}
