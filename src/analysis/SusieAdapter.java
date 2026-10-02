
import java.io.*;
import java.util.*;

/**
 * SuSiE adapter based on precodes/scripts/run_susier.R pattern.
 * Uses bigsnpr + susieR. Computes LD internally from PLINK genotypes.
 * Generates a self-contained R script that the plugin engine executes.
 */
public class SusieAdapter {

    /** App folder: LYNXGWAS_HOME, else the folder containing the compiled classes' bin/ directory. */
    static File appRoot() {
        String env = System.getenv("LYNXGWAS_HOME");
        if (env != null && !env.isEmpty()) return new File(env);
        try {
            File classes = new File(SusieAdapter.class.getProtectionDomain().getCodeSource().getLocation().toURI());
            return classes.getParentFile() != null ? classes.getParentFile() : new File(".");
        } catch (Exception e) {
            return new File(".");
        }
    }

    /** UK Biobank LD windows (resources/ukbb_ld under the app folder, or LYNXGWAS_UKBB_LD); "" if absent. */
    static String ukbbLdDir() {
        String env = System.getenv("LYNXGWAS_UKBB_LD");
        File d = env != null && !env.isEmpty() ? new File(env) : new File(appRoot(), "resources/ukbb_ld");
        return d.isDirectory() ? d.getAbsolutePath() : "";
    }

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

            pw.println("SNP\tBP\tA1\tA2\tBETA\tSE\tP\tN");
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
                String nSnp = f.length > 9 && !f[9].isEmpty() ? f[9] : "NA";   // per-SNP sample size (meta-analyses vary)
                pw.printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s%n", useId, pos, ea, nea, beta, se, f[5], nSnp);
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
        try (PrintWriter pw = new PrintWriter(new BufferedWriter(new OutputStreamWriter(new FileOutputStream(rScript), java.nio.charset.StandardCharsets.UTF_8)))) {
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
            // Large-sample LD reference (UK Biobank, scripts/ukbb_ld.py); empty dir = not installed
            pw.printf("ukbb_dir <- '%s'%n", ukbbLdDir().replace("\\", "/"));
            pw.printf("ukbb_script <- '%s'%n", new File(appRoot(), "scripts/ukbb_ld.py").getAbsolutePath().replace("\\", "/"));
            pw.println();

            // Load GWAS
            pw.println("df <- fread(input_tsv, sep='\\t', data.table=FALSE)");
            pw.println("df$BETA <- as.numeric(df$BETA); df$SE <- as.numeric(df$SE)");
            pw.println("df$P <- as.numeric(df$P); df$BP <- as.numeric(df$BP)");
            pw.println("bad <- !is.finite(df$BETA) | !is.finite(df$SE) | df$SE<=0 | !is.finite(df$P) | df$P<=0");
            pw.println("df <- df[!bad, , drop=FALSE]");
            // Meta-analyses report a different N per SNP (not every cohort has every SNP). SuSiE-RSS
            // assumes one N, and low-N SNPs then look inconsistent with the LD (extra spurious signals).
            // Keep SNPs with N >= 80% of the locus maximum when per-SNP N is available.
            pw.println("if ('N' %in% names(df)) {");
            pw.println("  df$N <- suppressWarnings(as.numeric(df$N))");
            pw.println("  if (sum(is.finite(df$N)) > 0.5 * nrow(df)) {");
            pw.println("    n_max <- max(df$N, na.rm=TRUE); low_n <- !is.finite(df$N) | df$N < 0.8 * n_max");
            pw.println("    if (any(low_n)) cat(sprintf('SuSiE: removed %d of %d SNPs with per-SNP N < 80%% of the maximum (%.0f): varying N breaks the single-N model\\n', sum(low_n), nrow(df), n_max))");
            pw.println("    df <- df[!low_n, , drop=FALSE]");
            pw.println("  }");
            pw.println("}");
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
            pw.println("if (length(common) < 10) stop(sprintf('Only %d SNPs in this locus window match the reference panel (at least 10 needed): the summary statistics are too sparse here to fine-map', length(common)))");
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
            // Prefer UK Biobank LD (337k Europeans) when a window covers the locus: a 503-person panel's
            // LD errors look like extra signals in GWAS with hundreds of thousands of samples.
            pw.println("ukb <- NULL");
            pw.println("py_bin <- Sys.getenv('LYNXGWAS_PYTHON'); if (!nzchar(py_bin)) py_bin <- Sys.which('python'); if (!nzchar(py_bin)) py_bin <- Sys.which('python3')");
            pw.println("if (nzchar(ukbb_dir) && dir.exists(ukbb_dir) && file.exists(ukbb_script) && nzchar(py_bin)) {");
            pw.println("  snp_in <- tempfile(fileext='.tsv'); out_pre <- tempfile()");
            pw.println("  write.table(data.frame(chr=bim_m$CHR, pos=bim_m$BP, a1=bim_m$A1, a2=bim_m$A2), snp_in, sep='\\t', quote=FALSE, row.names=FALSE)");
            pw.println("  msg <- suppressWarnings(system2(py_bin, c(shQuote(ukbb_script), shQuote(ukbb_dir), shQuote(snp_in), shQuote(out_pre)), stdout=TRUE, stderr=TRUE))");
            pw.println("  if (is.null(attr(msg, 'status')) && file.exists(paste0(out_pre, '.ld'))) {");
            pw.println("    ukb_idx <- scan(paste0(out_pre, '.idx'), quiet=TRUE)");
            pw.println("    if (length(ukb_idx) >= 10) {");
            pw.println("      ukb <- as.matrix(read.table(paste0(out_pre, '.ld')))");
            pw.println("      cat(sprintf('SuSiE: LD from UK Biobank (%d of %d SNPs matched; %s)\\n', length(ukb_idx), nrow(df), paste(msg, collapse=' ')))");
            pw.println("      df <- df[ukb_idx, , drop=FALSE]; bim_m <- bim_m[ukb_idx, , drop=FALSE]");
            pw.println("      R <- (ukb + t(ukb)) / 2; R[!is.finite(R)] <- 0; diag(R) <- 1.0; dimnames(R) <- NULL");
            pw.println("    }");
            pw.println("  }");
            pw.println("}");
            pw.println("if (is.null(ukb)) {");
            pw.println("cat('Computing LD matrix with bigsnpr (1000 Genomes reference panel)...\\n')");
            // Full LD matrix: snp_cor's default only correlates SNPs within 500 positions of each other and
            // leaves 0 elsewhere, which is not a valid correlation matrix for wide windows (hundreds of
            // negative eigenvalues -> nearPD for tens of minutes, and wrong LD for distant pairs).
            // One core per process: batch runs already run several loci in parallel.
            pw.println("R_raw <- as.matrix(snp_cor(Gna=G, ind.col=plink_idx, size=length(plink_idx), ncores=1L))");
            pw.println("R <- (R_raw + t(R_raw)) / 2; diag(R) <- 1.0");
            // SNPs that don't vary in the reference panel have undefined LD (NaN): drop them from the
            // matrix and the association data together instead of failing the whole locus.
            // Identify them by allele frequency (MAF 0 in the panel) or an undefined self-correlation:
            // counting NaNs per row is unreliable when many SNPs in the window are monomorphic.
            pw.println("maf_ld <- snp_MAF(G, ind.col=plink_idx)");
            pw.println("bad_ld <- sort(union(which(!is.finite(maf_ld) | maf_ld <= 0), which(!is.finite(diag(R_raw)))))");
            pw.println("if (length(bad_ld) > 0) {");
            pw.println("  cat(sprintf('SuSiE: dropped %d SNPs with undefined LD (monomorphic in the reference panel)\\n', length(bad_ld)))");
            pw.println("  keep_ld <- setdiff(seq_len(nrow(R)), bad_ld); R <- R[keep_ld, keep_ld, drop=FALSE]");
            pw.println("  df <- df[keep_ld, , drop=FALSE]; bim_m <- bim_m[keep_ld, , drop=FALSE]");
            pw.println("}");
            pw.println("n_undef <- sum(!is.finite(R))");
            pw.println("if (n_undef > 0) { cat(sprintf('SuSiE: %d remaining undefined LD entries set to 0\\n', n_undef)); R[!is.finite(R)] <- 0; diag(R) <- 1.0 }");
            pw.println("}  # end 1000 Genomes LD");
            pw.println("if (nrow(R) < 2) stop('Fewer than 2 SNPs with defined LD in this window')");
            pw.println();

            // Eigenvalue check + nearPD
            pw.println("eig_min <- min(eigen(R, symmetric=TRUE, only.values=TRUE)$values)");
            // A full reference-panel LD matrix is positive semi-definite up to rounding; repair only real violations
            pw.println("if (eig_min < -1e-6) { cat(sprintf('SuSiE: LD not positive definite (min eigenvalue %.3g); projecting to nearest correlation matrix\\n', eig_min)); R <- as.matrix(nearPD(R, corr=TRUE, keepDiag=TRUE, maxit=20)$mat); diag(R) <- 1.0 }");
            pw.println("R_ld <- R  # LD as estimated; shrinkage below is for fitting only, purity uses the real LD");
            // LD consistency check (reference panel vs GWAS). kriging_rss predicts each SNP's z from the
            // others through the LD; SNPs that contradict it (allele flips, imputation errors, panel
            // differences) are removed a few at a time, then estimate_s_rss measures the remaining
            // inconsistency, which becomes the LD regularisation when it exceeds the default.
            pw.println("z_qc <- df$BETA / df$SE");
            pw.println("ld_qc_removed <- 0L");
            pw.println("for (qc_iter in 1:3) {");
            pw.println("  kr <- tryCatch(kriging_rss(z_qc, R_ld, n_samples)$conditional_dist, error=function(e) NULL)");
            pw.println("  if (is.null(kr)) break");
            pw.println("  bad <- which(kr$logLR > 2 & abs(kr$z) > 2)");
            pw.println("  if (length(bad) == 0) break");
            pw.println("  bad <- head(bad[order(-kr$logLR[bad])], max(1L, ceiling(0.05 * length(z_qc))))");
            pw.println("  keep_qc <- setdiff(seq_along(z_qc), bad)");
            pw.println("  R_ld <- R_ld[keep_qc, keep_qc, drop=FALSE]; df <- df[keep_qc, , drop=FALSE]; bim_m <- bim_m[keep_qc, , drop=FALSE]; z_qc <- z_qc[keep_qc]");
            pw.println("  ld_qc_removed <- ld_qc_removed + length(bad)");
            pw.println("}");
            pw.println("if (ld_qc_removed > 0) cat(sprintf('SuSiE LD check: removed %d SNPs whose z-scores contradict the reference LD\\n', ld_qc_removed))");
            pw.println("s_rss <- tryCatch(estimate_s_rss(z_qc, R_ld, n_samples), error=function(e) NA_real_)");
            pw.println("if (is.finite(s_rss)) cat(sprintf('SuSiE LD check: GWAS vs reference LD inconsistency s = %.3f\\n', s_rss))");
            pw.println("if (is.finite(s_rss) && s_rss > ld_shrink) { ld_shrink <- min(s_rss, 0.5); cat(sprintf('SuSiE LD check: LD regularisation raised to %.3f\\n', ld_shrink)) }");
            pw.println("R <- R_ld");
            pw.println("if (ld_shrink > 0) R <- (1-ld_shrink)*R + ld_shrink*diag(nrow(R))");
            pw.println();

            // Run SuSiE
            pw.println("z <- df$BETA / df$SE");
            pw.println("cat(sprintf('Running susie_rss: %d SNPs, N=%d, L=%d\\n', length(z), n_samples, L))");
            // estimate_residual_variance=TRUE assumes R is the true in-sample LD matrix; R here comes
            // from a reference panel (bigsnpr), not in-sample genotypes, so it must be FALSE or
            // susie_rss can fail with "Estimating residual variance failed: the estimated value is negative".
            pw.println("fit <- susie_rss(z=z, R=R, n=n_samples, L=L, coverage=coverage, estimate_residual_variance=FALSE, verbose=FALSE)");
            pw.println("pips <- susie_get_pip(fit)");
            pw.println();

            // Credible sets
            pw.println("cs_args <- names(formals(susie_get_cs))");
            // Purity (min |r| within a credible set) must be judged on the real LD: the shrunk matrix
            // scales every correlation by (1-ld_shrink) and drops sets sitting just above 0.5.
            pw.println("cs_obj <- if ('Rr' %in% cs_args) { susie_get_cs(fit, coverage=coverage, Rr=R_ld) } else if ('Xcorr' %in% cs_args) {");
            pw.println("  susie_get_cs(fit, coverage=coverage, Xcorr=R_ld) } else { susie_get_cs(fit, coverage=coverage) }");
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
