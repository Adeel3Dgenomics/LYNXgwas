
import java.io.*;
import java.util.*;

/**
 * SuSiE adapter based on precodes/scripts/run_susier.R pattern.
 * Uses bigsnpr + susieR. Computes LD internally from PLINK genotypes.
 * Generates a self-contained R script that the plugin engine executes.
 */
public class SusieAdapter {

    public static void prepareRun(File harmonizedDir, File matchedDir, File runDir,
                                   int sampleN, int maxCausal, double coverage,
                                   double ldShrink, int windowKb) throws IOException {
        if (sampleN <= 0)
            throw new IOException("Sample size (N) is required for SuSiE. Set it in the project configuration.");
        runDir.mkdirs();

        // Write matched GWAS in the format the R script expects (SNP BP A1 A2 BETA SE P)
        File gwasFile = new File(harmonizedDir, "harmonized_gwas.tsv");
        File inputTsv = new File(runDir, "susie_input.tsv");
        Map<String, String> posToRefId = readBimIds(new File(matchedDir, "matched_ref.bim"));

        int count = 0;
        try (BufferedReader br = new BufferedReader(new FileReader(gwasFile));
             PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(inputTsv)))) {

            pw.println("SNP\tBP\tA1\tA2\tBETA\tSE\tP");
            String header = br.readLine();
            String line;
            while ((line = br.readLine()) != null) {
                if (line.isEmpty()) continue;
                String[] f = line.split("\t", -1);
                if (f.length < 8) continue;
                String chr = f[1], pos = f[2], ea = f[3], nea = f[4];
                String beta = f[6], se = f[7];
                if (beta.equals("NA") || se.equals("NA")) continue;
                try {
                    double b = Double.parseDouble(beta), s = Double.parseDouble(se);
                    if (Double.isNaN(b) || Double.isNaN(s) || s <= 0) continue;
                } catch (NumberFormatException e) { continue; }

                String refId = posToRefId.get(chr + ":" + pos);
                String useId = refId != null ? refId : f[0];
                pw.printf("%s\t%s\t%s\t%s\t%s\t%s\t%s%n", useId, pos, ea, nea, beta, se, f[5]);
                count++;
            }
        }

        String plinkPrefix = new File(matchedDir, "matched_ref").getAbsolutePath().replace("\\", "/");
        String inputPath = inputTsv.getAbsolutePath().replace("\\", "/");
        String outTsv = new File(runDir, "result.tsv").getAbsolutePath().replace("\\", "/");
        String diagJson = new File(runDir, "susie_diag.json").getAbsolutePath().replace("\\", "/");
        String diagLog = new File(runDir, "susie_diag.log").getAbsolutePath().replace("\\", "/");
        String plotOverview = new File(runDir, "susie_overview.png").getAbsolutePath().replace("\\", "/");
        String manifestPath = new File(runDir, "result.manifest.json").getAbsolutePath().replace("\\", "/");

        File rScript = new File(runDir, "susie_run.R");
        try (PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(rScript)))) {
            pw.println("#!/usr/bin/env Rscript");
            pw.println("suppressPackageStartupMessages({");
            pw.println("  library(bigsnpr)");
            pw.println("  library(susieR)");
            pw.println("  library(data.table)");
            pw.println("  library(Matrix)");
            pw.println("})");
            pw.println();
            pw.printf("input_tsv <- '%s'%n", inputPath);
            pw.printf("plink_prefix <- '%s'%n", plinkPrefix);
            pw.printf("n_samples <- %d%n", sampleN);
            pw.printf("L <- %d%n", maxCausal);
            pw.printf("coverage <- %.2f%n", coverage);
            pw.printf("ld_shrink <- %.3f%n", ldShrink);
            pw.printf("window_kb <- %d%n", windowKb);
            pw.printf("out_tsv <- '%s'%n", outTsv);
            pw.printf("diag_json <- '%s'%n", diagJson);
            pw.printf("manifest_path <- '%s'%n", manifestPath);
            pw.println();

            // Load GWAS
            pw.println("df <- fread(input_tsv, sep='\\t', data.table=FALSE)");
            pw.println("df$BETA <- as.numeric(df$BETA); df$SE <- as.numeric(df$SE)");
            pw.println("df$P <- as.numeric(df$P); df$BP <- as.numeric(df$BP)");
            pw.println("bad <- !is.finite(df$BETA) | !is.finite(df$SE) | df$SE<=0 | !is.finite(df$P) | df$P<=0");
            pw.println("df <- df[!bad, , drop=FALSE]");
            pw.println();

            // Window around index SNP
            pw.println("index_bp <- df$BP[which.min(df$P)]");
            pw.println("bp_lo <- index_bp - window_kb*1000; bp_hi <- index_bp + window_kb*1000");
            pw.println("df <- df[df$BP >= bp_lo & df$BP <= bp_hi, , drop=FALSE]");
            pw.println("df <- df[order(df$BP), ]");
            pw.println("cat(sprintf('SuSiE: %d SNPs in window\\n', nrow(df)))");
            pw.println();

            // PLINK intersection
            pw.println("bim <- fread(paste0(plink_prefix, '.bim'), header=FALSE,");
            pw.println("  col.names=c('CHR','SNP','CM','BP','A1','A2'), data.table=FALSE)");
            pw.println("bim_win <- bim[bim$BP >= bp_lo & bim$BP <= bp_hi, , drop=FALSE]");
            pw.println("common <- intersect(df$SNP, bim_win$SNP)");
            pw.println("if (length(common) < 10) stop('Too few SNPs matched between GWAS and PLINK')");
            pw.println("df <- df[df$SNP %in% common, , drop=FALSE]");
            pw.println("df <- df[match(bim_win$SNP[bim_win$SNP %in% common], df$SNP), , drop=FALSE]");
            pw.println("cat(sprintf('SuSiE: %d SNPs matched with PLINK\\n', nrow(df)))");
            pw.println();

            // Allele harmonization
            pw.println("bim_m <- bim_win[match(df$SNP, bim_win$SNP), , drop=FALSE]");
            pw.println("flip <- (df$A1 == bim_m$A2 & df$A2 == bim_m$A1)");
            pw.println("df$BETA[flip] <- -df$BETA[flip]");
            pw.println("cat(sprintf('SuSiE: %d allele flips applied\\n', sum(flip)))");
            pw.println();

            // bigsnpr LD computation
            pw.println("tmp_dir <- tempdir()");
            pw.println("backing <- file.path(tmp_dir, paste0('bigsnp_', basename(plink_prefix)))");
            pw.println("rds_path <- paste0(backing, '.rds')");
            pw.println("if (!file.exists(rds_path)) {");
            pw.println("  rds_created <- snp_readBed(paste0(plink_prefix, '.bed'), backingfile=backing)");
            pw.println("  snp_obj <- snp_attach(rds_created)");
            pw.println("} else { snp_obj <- snp_attach(rds_path) }");
            pw.println("G <- snp_obj$genotypes; bim_full <- snp_obj$map");
            pw.println("plink_idx <- match(df$SNP, as.character(bim_full$marker.ID))");
            pw.println("if (any(is.na(plink_idx))) stop('SNPs not found in PLINK bim')");
            pw.println();
            pw.println("cat('Computing LD matrix with bigsnpr...\\n')");
            pw.println("R_raw <- as.matrix(snp_cor(Gna=G, ind.col=plink_idx, ncores=max(1L, parallel::detectCores()-1L)))");
            pw.println("R <- (R_raw + t(R_raw)) / 2; diag(R) <- 1.0");
            pw.println();

            // Eigenvalue check + nearPD
            pw.println("eig_min <- min(eigen(R, symmetric=TRUE, only.values=TRUE)$values)");
            pw.println("if (eig_min < 0) { R <- as.matrix(nearPD(R, corr=TRUE, keepDiag=TRUE)$mat); diag(R) <- 1.0 }");
            pw.println("if (ld_shrink > 0) R <- (1-ld_shrink)*R + ld_shrink*diag(nrow(R))");
            pw.println();

            // Run SuSiE
            pw.println("z <- df$BETA / df$SE");
            pw.println("cat(sprintf('Running susie_rss: %d SNPs, N=%d, L=%d\\n', length(z), n_samples, L))");
            pw.println("fit <- susie_rss(z=z, R=R, n=n_samples, L=L, coverage=coverage, verbose=FALSE)");
            pw.println("pips <- susie_get_pip(fit)");
            pw.println();

            // Credible sets
            pw.println("cs_args <- names(formals(susie_get_cs))");
            pw.println("cs_obj <- if ('Rr' %in% cs_args) susie_get_cs(fit, coverage=coverage, Rr=R)");
            pw.println("  else if ('Xcorr' %in% cs_args) susie_get_cs(fit, coverage=coverage, Xcorr=R)");
            pw.println("  else susie_get_cs(fit, coverage=coverage)");
            pw.println("cs_list <- cs_obj$cs; n_cs <- length(cs_list)");
            pw.println("cs_member <- rep(NA_integer_, length(pips))");
            pw.println("cs_cover <- rep(NA_real_, length(pips))");
            pw.println("for (ci in seq_along(cs_list)) {");
            pw.println("  idx <- cs_list[[ci]]; cs_member[idx] <- ci; cs_cover[idx] <- sum(pips[idx])");
            pw.println("}");
            pw.println();

            // Write result
            pw.println("result <- data.frame(snp_id=df$SNP, chr=bim_m$CHR, pos=df$BP,");
            pw.println("  susie_pip=round(pips,6),");
            pw.println("  susie_cs=cs_member, susie_cs_coverage=round(cs_cover,4), stringsAsFactors=FALSE)");
            pw.println("write.table(result, out_tsv, sep='\\t', quote=FALSE, row.names=FALSE)");
            pw.println();

            // Write manifest
            pw.println("manifest <- paste0('{\"schema_version\":\"1.0\",\"method\":\"susie_finemapping\",\"method_version\":\"1.0\",');");
            pw.println("manifest <- paste0(manifest, '\"parameters\":{},\"columns\":[');");
            pw.println("manifest <- paste0(manifest, '{\"name\":\"chr\",\"type\":\"string\",\"scope\":\"per_snp\",\"method\":\"susie_finemapping\",\"method_version\":\"1.0\"},');");
            pw.println("manifest <- paste0(manifest, '{\"name\":\"pos\",\"type\":\"int\",\"scope\":\"per_snp\",\"method\":\"susie_finemapping\",\"method_version\":\"1.0\"},');");
            pw.println("manifest <- paste0(manifest, '{\"name\":\"susie_pip\",\"type\":\"double\",\"scope\":\"per_snp\",\"method\":\"susie_finemapping\",\"method_version\":\"1.0\"},');");
            pw.println("manifest <- paste0(manifest, '{\"name\":\"susie_cs\",\"type\":\"int\",\"scope\":\"per_credible_set\",\"method\":\"susie_finemapping\",\"method_version\":\"1.0\"},');");
            pw.println("manifest <- paste0(manifest, '{\"name\":\"susie_cs_coverage\",\"type\":\"double\",\"scope\":\"per_credible_set\",\"method\":\"susie_finemapping\",\"method_version\":\"1.0\"}');");
            pw.println("manifest <- paste0(manifest, '],\"created_at\":', as.numeric(Sys.time())*1000, '}');");
            pw.println("writeLines(manifest, manifest_path)");
            pw.println();
            pw.println("cat(sprintf('SuSiE: %d credible sets, max PIP=%.4f\\n', n_cs, max(pips)))");
        }

        System.out.printf("[SusieAdapter] Prepared %d SNPs, R script at %s%n", count, rScript.getName());
    }

    private static Map<String, String> readBimIds(File bimFile) throws IOException {
        Map<String, String> map = new LinkedHashMap<>();
        if (!bimFile.exists()) return map;
        try (BufferedReader br = new BufferedReader(new FileReader(bimFile))) {
            String line;
            while ((line = br.readLine()) != null) {
                String[] f = line.split("\t", -1);
                if (f.length < 4) continue;
                String chr = f[0].replaceFirst("^chr", "");
                map.put(chr + ":" + f[3].trim(), f[1]);
            }
        }
        return map;
    }
}
