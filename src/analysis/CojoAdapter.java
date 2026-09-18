
import java.io.*;
import java.util.*;

/**
 * COJO (GCTA) adapter — v2.1: delegates selection to GCTA --cojo-slct.
 *
 * Fixes over v2.0:
 *   - Allele/beta/freq orientation explicitly verified per SNP (Fix 1)
 *   - SNPs without ref panel ID dropped, never fallback to GWAS ID (Fix 2)
 *   - Per-SNP N accepted without requiring global N (Fix 3)
 *   - LD-GWAS consistency diagnostic wired in before COJO (Fix 4)
 *   - GCTA freq-discrepancy warnings parsed and reported (Fix 5)
 *   - Empty .jma.cojo handled as 0 signals cleanly (Fix 6)
 */
public class CojoAdapter {

    public static int prepareRun(File harmonizedDir, File matchedDir, File runDir,
                                  int sampleN, double pCutoff, double collinear,
                                  String gctaBin) throws IOException {

        // ── Verify GCTA binary (shared with GctaGremlAdapter — see GctaBinaryResolver) ──
        File gctaFile = GctaBinaryResolver.verify(gctaBin);

        runDir.mkdirs();

        File gwasFile = new File(harmonizedDir, "harmonized_gwas.tsv");
        if (!gwasFile.exists()) return 0;

        // ── Read ref panel BIM: chr:pos -> (snpId, A1, A2) ──
        Map<String, String> posToRefId = new LinkedHashMap<>();
        Map<String, String> refBimA1 = new LinkedHashMap<>();  // snpId -> A1 (col5)
        Map<String, String> refBimA2 = new LinkedHashMap<>();  // snpId -> A2 (col6)
        File bimFile = new File(matchedDir, "matched_ref.bim");
        if (bimFile.exists()) {
            try (BufferedReader br = new BufferedReader(new FileReader(bimFile))) {
                String line;
                while ((line = br.readLine()) != null) {
                    String[] f = line.split("\t", -1);
                    if (f.length < 6) continue;
                    String chr = f[0].replaceFirst("^chr", "");
                    String snpId = f[1];
                    posToRefId.put(chr + ":" + f[3].trim(), snpId);
                    refBimA1.put(snpId, f[4].toUpperCase());
                    refBimA2.put(snpId, f[5].toUpperCase());
                }
            }
        }

        // ── Load ref panel A1 frequencies (base pipeline Step 2.5) ──
        Map<String, Double> refFreq = new LinkedHashMap<>();
        Map<String, String> refFreqA1 = new LinkedHashMap<>();
        File refFreqFile = new File(matchedDir, "ref_freq.tsv");
        if (refFreqFile.exists()) {
            try (BufferedReader br = new BufferedReader(new FileReader(refFreqFile))) {
                br.readLine();
                String line;
                while ((line = br.readLine()) != null) {
                    String[] f = line.split("\t", -1);
                    if (f.length < 3) continue;
                    try {
                        double freq = Double.parseDouble(f[2]);
                        if (freq > 0 && freq < 1) {
                            refFreq.put(f[0], freq);
                            refFreqA1.put(f[0], f[1].toUpperCase());
                        }
                    } catch (NumberFormatException ignored) {}
                }
            }
            System.out.printf("[CojoAdapter] Loaded %d ref panel frequencies%n", refFreq.size());
        } else {
            throw new IOException("ref_freq.tsv not found in matched directory. "
                + "Run the base pipeline first to compute ref panel frequencies.");
        }

        // ── Write .ma file with per-SNP audit ──
        File maFile = new File(runDir, "cojo_input.ma");
        File auditFile = new File(runDir, "ma_audit.tsv");
        int written = 0, droppedNoRefId = 0, droppedNoFreq = 0;
        int droppedOrientation = 0, droppedStats = 0;
        int freqFlipped = 0;
        boolean hasPerSnpN = false;
        boolean anySnpHasN = false;

        // First pass: check if per-SNP N exists
        try (BufferedReader br = new BufferedReader(new FileReader(gwasFile))) {
            String header = br.readLine();
            if (header != null) {
                String[] hcols = header.split("\t", -1);
                int iN = -1;
                for (int i = 0; i < hcols.length; i++) {
                    if (hcols[i].trim().equalsIgnoreCase("n")) { iN = i; break; }
                }
                if (iN >= 0) {
                    String line = br.readLine();
                    if (line != null) {
                        String[] f = line.split("\t", -1);
                        if (iN < f.length) {
                            String nVal = f[iN].trim();
                            if (!nVal.isEmpty() && !nVal.equals("NA")) {
                                try {
                                    if ((int) Double.parseDouble(nVal) > 0) anySnpHasN = true;
                                } catch (NumberFormatException ignored) {}
                            }
                        }
                    }
                }
            }
        }

        // Fix 3: require either global N or per-SNP N
        if (sampleN <= 0 && !anySnpHasN)
            throw new IOException("Sample size (N) is required for COJO. "
                + "Set it in the project configuration or ensure the GWAS has a per-SNP N column.");

        try (BufferedReader br = new BufferedReader(new FileReader(gwasFile));
             PrintWriter ma = new PrintWriter(new BufferedWriter(new FileWriter(maFile)));
             PrintWriter audit = new PrintWriter(new BufferedWriter(new FileWriter(auditFile)))) {

            ma.println("SNP\tA1\tA2\tfreq\tb\tse\tp\tN");
            audit.println("snp_id\tref_id\tstatus\tdetail");

            String header = br.readLine();
            if (header == null) return 0;

            String[] hcols = header.split("\t", -1);
            int iN = -1;
            for (int i = 0; i < hcols.length; i++) {
                if (hcols[i].trim().equalsIgnoreCase("n")) { iN = i; break; }
            }

            String line;
            while ((line = br.readLine()) != null) {
                if (line.isEmpty()) continue;
                String[] f = line.split("\t", -1);
                if (f.length < 8) continue;

                String snpId = f[0];
                String chr   = f[1];
                String pos   = f[2];
                String ea    = f[3].toUpperCase();
                String nea   = f[4].toUpperCase();
                String pval  = f[5];
                String beta  = f[6];
                String se    = f[7];

                if (beta.equals("NA") || se.equals("NA")) { droppedStats++; continue; }
                double betaVal, seVal;
                try {
                    betaVal = Double.parseDouble(beta);
                    seVal = Double.parseDouble(se);
                } catch (NumberFormatException e) { droppedStats++; continue; }
                if (Double.isNaN(betaVal) || Double.isNaN(seVal) || seVal <= 0) {
                    droppedStats++; continue;
                }

                // Fix 2: require ref panel ID — never fall back to GWAS ID
                String refId = posToRefId.get(chr + ":" + pos);
                if (refId == null) {
                    droppedNoRefId++;
                    audit.printf("%s\t\tdropped_no_ref_id\tchr%s:%s not in ref panel BIM%n",
                        snpId, chr, pos);
                    continue;
                }

                // Get ref panel freq
                Double refF = refFreq.get(refId);
                if (refF == null) {
                    droppedNoFreq++;
                    audit.printf("%s\t%s\tdropped_no_freq\tno frequency in ref_freq.tsv%n",
                        snpId, refId);
                    continue;
                }

                // Fix 1: verify allele orientation
                // Convention: A1 = harmonized EA (beta refers to this allele)
                //             freq = frequency of A1 (EA) in the ref panel
                String rA1 = refBimA1.get(refId);
                String rA2 = refBimA2.get(refId);
                String refFa1 = refFreqA1.get(refId);

                double freqEA;
                String orientStatus;

                if (ea.equals(rA1) || ea.equals(refFa1)) {
                    // EA matches ref A1 — freq is already EA's frequency
                    freqEA = refF;
                    orientStatus = "kept";
                } else if (ea.equals(rA2) || (rA2 != null && ea.equals(rA2))) {
                    // EA matches ref A2 — flip freq to get EA's frequency
                    freqEA = 1.0 - refF;
                    orientStatus = "freq_flipped";
                    freqFlipped++;
                } else {
                    // EA matches neither ref allele — drop
                    droppedOrientation++;
                    audit.printf("%s\t%s\tdropped_orientation\tEA=%s NEA=%s vs ref A1=%s A2=%s%n",
                        snpId, refId, ea, nea, rA1 != null ? rA1 : "?", rA2 != null ? rA2 : "?");
                    continue;
                }

                // Per-SNP N (Fix 3)
                int snpN = sampleN > 0 ? sampleN : 0;
                if (iN >= 0 && iN < f.length && !f[iN].trim().isEmpty()
                        && !f[iN].trim().equals("NA")) {
                    try {
                        int parsed = (int) Double.parseDouble(f[iN].trim());
                        if (parsed > 0) {
                            snpN = parsed;
                            hasPerSnpN = true;
                        }
                    } catch (NumberFormatException ignored) {}
                }
                if (snpN <= 0) {
                    audit.printf("%s\t%s\tdropped_no_n\tno valid sample size%n", snpId, refId);
                    continue;
                }

                // Write .ma row: A1=EA, freq=freq(EA), b=beta(EA), unchanged
                ma.printf("%s\t%s\t%s\t%.6f\t%.6g\t%.6g\t%s\t%d%n",
                    refId, ea, nea, freqEA, betaVal, seVal, pval, snpN);
                audit.printf("%s\t%s\t%s\tfreq=%.4f%n", snpId, refId, orientStatus, freqEA);
                written++;
            }
        }

        String nSource = hasPerSnpN ? "per_snp" : "global";
        System.out.printf("[CojoAdapter] .ma: %d written, dropped: %d no-ref-id, %d no-freq, "
            + "%d orientation, %d bad-stats, %d freq-flipped. N-source: %s%n",
            written, droppedNoRefId, droppedNoFreq, droppedOrientation, droppedStats,
            freqFlipped, nSource);

        // ── Run LD-GWAS consistency diagnostic before COJO (Fix 4) ──
        File ldDir = new File(harmonizedDir.getParentFile(), "ld");
        String ldConsistencyVerdict = "not_run";
        int ldConsistencyFlagged = 0;
        if (ldDir.exists()) {
            File summaryFile = new File(ldDir, "consistency_summary.json");
            if (summaryFile.exists()) {
                try {
                    String sumJson = new String(java.nio.file.Files.readAllBytes(summaryFile.toPath()), "UTF-8");
                    // Simple parse
                    int vi = sumJson.indexOf("\"verdict\"");
                    if (vi >= 0) {
                        int vs = sumJson.indexOf('"', vi + 10);
                        int ve = sumJson.indexOf('"', vs + 1);
                        if (vs >= 0 && ve > vs) ldConsistencyVerdict = sumJson.substring(vs + 1, ve);
                    }
                    int fi = sumJson.indexOf("\"flagged_snps\"");
                    if (fi >= 0) {
                        int fs = fi + 15;
                        while (fs < sumJson.length() && !Character.isDigit(sumJson.charAt(fs))) fs++;
                        int fe = fs;
                        while (fe < sumJson.length() && Character.isDigit(sumJson.charAt(fe))) fe++;
                        if (fe > fs) ldConsistencyFlagged = Integer.parseInt(sumJson.substring(fs, fe));
                    }
                } catch (Exception ignored) {}
            }
        }

        // ── Generate R script ──
        String bfilePrefix = new File(matchedDir, "matched_ref").getAbsolutePath().replace("\\", "/");
        String maPath = maFile.getAbsolutePath().replace("\\", "/");
        String outTsv = new File(runDir, "result.tsv").getAbsolutePath().replace("\\", "/");
        String diagJson = new File(runDir, "cojo_diag.json").getAbsolutePath().replace("\\", "/");
        String manifestPath = new File(runDir, "result.manifest.json").getAbsolutePath().replace("\\", "/");
        String gctaAbsolute = gctaFile.getAbsolutePath().replace("\\", "/");
        String cojoPrefix = new File(runDir, "cojo").getAbsolutePath().replace("\\", "/");

        File rScript = new File(runDir, "run_cojo.R");
        try (PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(rScript)))) {
            pw.println("#!/usr/bin/env Rscript");
            pw.printf("gcta_bin <- '%s'%n", gctaAbsolute);
            pw.printf("ma_file <- '%s'%n", maPath);
            pw.printf("bfile_prefix <- '%s'%n", bfilePrefix);
            pw.printf("p_cutoff <- %s%n", String.valueOf(pCutoff));
            pw.printf("collinear_thr <- %s%n", String.valueOf(collinear));
            pw.printf("out_prefix <- '%s'%n", cojoPrefix);
            pw.printf("out_tsv <- '%s'%n", outTsv);
            pw.printf("diag_json <- '%s'%n", diagJson);
            pw.printf("manifest_path <- '%s'%n", manifestPath);
            pw.printf("sample_size <- %d%n", sampleN > 0 ? sampleN : 0);
            pw.printf("freq_source <- 'ref_panel'%n");
            pw.printf("n_source <- '%s'%n", nSource);
            pw.printf("snps_written <- %d%n", written);
            pw.printf("snps_dropped_no_ref_id <- %d%n", droppedNoRefId);
            pw.printf("snps_dropped_no_freq <- %d%n", droppedNoFreq);
            pw.printf("snps_dropped_orientation <- %d%n", droppedOrientation);
            pw.printf("snps_freq_flipped <- %d%n", freqFlipped);
            pw.printf("ld_consistency_verdict <- '%s'%n", ldConsistencyVerdict);
            pw.printf("ld_consistency_flagged <- %d%n", ldConsistencyFlagged);
            pw.printf("artifact_inflation_factor <- 5%n");
            pw.println();

            // Verify + version
            pw.println("if (!file.exists(gcta_bin)) stop(paste0('GCTA binary not found: ', gcta_bin))");
            pw.println("gcta_ver <- tryCatch({ paste(system2(gcta_bin, '--version', stdout=TRUE, stderr=TRUE), collapse=' ') }, error=function(e) 'unknown')");
            pw.println("cat(sprintf('GCTA: %s\\nVersion: %s\\n', gcta_bin, gcta_ver))");
            pw.println();

            // LD-consistency warning
            pw.println("if (ld_consistency_verdict == 'high_warn') {");
            pw.println("  cat(sprintf('WARNING: LD-GWAS consistency = %s (%d SNPs flagged)\\n', ld_consistency_verdict, ld_consistency_flagged))");
            pw.println("  cat('  Joint estimates assume LD panel matches GWAS population — cross-ancestry mismatch may produce artifacts.\\n')");
            pw.println("}");
            pw.println();

            // ═══ Step 1: GCTA --cojo-slct ═══
            pw.println("cat(sprintf('SNPs in .ma: %d\\n', snps_written))");
            pw.println("cojo_args <- c('--bfile', bfile_prefix, '--cojo-file', ma_file, '--cojo-slct',");
            pw.println("  '--cojo-p', format(p_cutoff, scientific=TRUE), '--cojo-collinear', as.character(collinear_thr), '--out', out_prefix)");
            pw.println("cat(sprintf('Command: %s %s\\n', gcta_bin, paste(cojo_args, collapse=' ')))");
            pw.println("cojo_out <- system2(gcta_bin, args=cojo_args, stdout=TRUE, stderr=TRUE)");
            pw.println("cojo_exit <- attr(cojo_out, 'status'); if (is.null(cojo_exit)) cojo_exit <- 0L");
            pw.println("writeLines(cojo_out, paste0(out_prefix, '_stdout.log'))");
            pw.println("cat(paste(cojo_out, collapse='\\n'), '\\n')");
            pw.println();

            // Build result from .ma
            pw.println("ma_data <- read.table(ma_file, header=TRUE, stringsAsFactors=FALSE)");
            pw.println("bim <- read.table(paste0(bfile_prefix, '.bim'), stringsAsFactors=FALSE)");
            pw.println("result <- data.frame(snp_id=ma_data$SNP, stringsAsFactors=FALSE)");
            pw.println("bim_idx <- match(result$snp_id, bim$V2)");
            pw.println("result$chr <- ifelse(is.na(bim_idx), NA, bim$V1[bim_idx])");
            pw.println("result$pos <- ifelse(is.na(bim_idx), NA, bim$V4[bim_idx])");
            pw.println("result$cojo_p_marginal <- ma_data$p");
            pw.println("result$cojo_b_marginal <- ma_data$b");
            pw.println("result$cojo_freq <- ma_data$freq");
            pw.println("result$cojo_n <- ma_data$N");
            pw.println("result$cojo_pJ <- NA_real_; result$cojo_bJ <- NA_real_; result$cojo_bJ_se <- NA_real_");
            pw.println("result$cojo_pC <- NA_real_; result$cojo_bC <- NA_real_; result$cojo_bC_se <- NA_real_");
            pw.println("result$cojo_selected <- ''");
            pw.println();

            // ═══ Parse .jma.cojo (joint estimates for selected SNPs) ═══
            pw.println("jma_file <- paste0(out_prefix, '.jma.cojo')");
            pw.println("n_selected <- 0L; selected_snps <- character(0); freq_discrepancies <- data.frame()");
            pw.println("if (file.exists(jma_file) && file.info(jma_file)$size > 0) {");
            pw.println("  jma <- tryCatch(read.table(jma_file, header=TRUE, stringsAsFactors=FALSE), error=function(e) data.frame())");
            pw.println("  if (nrow(jma) > 0) {");
            pw.println("    n_selected <- nrow(jma); selected_snps <- jma$SNP");
            pw.println("    cat(sprintf('GCTA selected %d independent SNP(s)\\n', n_selected))");
            pw.println("    for (j in seq_len(nrow(jma))) {");
            pw.println("      idx <- which(result$snp_id == jma$SNP[j])");
            pw.println("      if (length(idx) > 0) {");
            pw.println("        result$cojo_pJ[idx] <- jma$pJ[j]; result$cojo_bJ[idx] <- jma$bJ[j]; result$cojo_bJ_se[idx] <- jma$bJ_se[j]");
            pw.println("        result$cojo_selected[idx] <- 'selected'");
            pw.println("    } }");
            pw.println("    if ('freq_geno' %in% colnames(jma) && 'freq' %in% colnames(jma)) {");
            pw.println("      jma$freq_diff <- abs(jma$freq - jma$freq_geno)");
            pw.println("      bad_freq <- jma[jma$freq_diff > 0.15, , drop=FALSE]");
            pw.println("      if (nrow(bad_freq) > 0) freq_discrepancies <- bad_freq[, c('SNP','freq','freq_geno','freq_diff')]");
            pw.println("    }");
            pw.println("  } else { cat('Header-only .jma.cojo — 0 signals\\n') }");
            pw.println("} else { cat(if (cojo_exit==0) 'No .jma.cojo — 0 signals\\n' else sprintf('GCTA exited %d\\n', cojo_exit)) }");
            pw.println();

            // ═══ Step 2: --cojo-cond on selected set → pC for ALL SNPs ═══
            pw.println("cond_exit <- NA_integer_");
            pw.println("if (n_selected > 0) {");
            pw.println("  cat('Running --cojo-cond with selected SNPs for conditional p-values...\\n')");
            pw.println("  cond_snp_file <- paste0(out_prefix, '_selected.snps')");
            pw.println("  writeLines(selected_snps, cond_snp_file)");
            pw.println("  cond_prefix <- paste0(out_prefix, '_final_cond')");
            pw.println("  cond_args <- c('--bfile', bfile_prefix, '--cojo-file', ma_file,");
            pw.println("    '--cojo-cond', cond_snp_file, '--out', cond_prefix)");
            pw.println("  cond_out <- system2(gcta_bin, args=cond_args, stdout=TRUE, stderr=TRUE)");
            pw.println("  cond_exit <- attr(cond_out, 'status'); if (is.null(cond_exit)) cond_exit <- 0L");
            pw.println("  cma_file <- paste0(cond_prefix, '.cma.cojo')");
            pw.println("  if (file.exists(cma_file)) {");
            pw.println("    cma <- tryCatch(read.table(cma_file, header=TRUE, stringsAsFactors=FALSE), error=function(e) data.frame())");
            pw.println("    if (nrow(cma) > 0) {");
            pw.println("      p_col <- if ('pC' %in% colnames(cma)) 'pC' else 'p'");
            pw.println("      b_col <- if ('bC' %in% colnames(cma)) 'bC' else 'b'");
            pw.println("      se_col <- if ('bC_se' %in% colnames(cma)) 'bC_se' else 'se'");
            pw.println("      for (j in seq_len(nrow(cma))) {");
            pw.println("        idx <- which(result$snp_id == cma$SNP[j])");
            pw.println("        if (length(idx) > 0) {");
            pw.println("          result$cojo_pC[idx] <- cma[[p_col]][j]");
            pw.println("          result$cojo_bC[idx] <- cma[[b_col]][j]");
            pw.println("          result$cojo_bC_se[idx] <- cma[[se_col]][j]");
            pw.println("        }");
            pw.println("      }");
            pw.println("      cat(sprintf('  Conditional p-values computed for %d SNPs\\n', nrow(cma)))");
            pw.println("    }");
            pw.println("  } else { cat('  WARNING: --cojo-cond produced no output\\n') }");
            pw.println("  # Selected SNPs: their conditional p = their joint p (they are the model)");
            pw.println("  for (s in selected_snps) {");
            pw.println("    idx <- which(result$snp_id == s)");
            pw.println("    if (length(idx) > 0) {");
            pw.println("      result$cojo_pC[idx] <- result$cojo_pJ[idx]");
            pw.println("      result$cojo_bC[idx] <- result$cojo_bJ[idx]");
            pw.println("      result$cojo_bC_se[idx] <- result$cojo_bJ_se[idx]");
            pw.println("    }");
            pw.println("  }");
            pw.println("}");
            pw.println();

            // ═══ Locus-fully-explained check ═══
            pw.println("non_sel_pC <- result$cojo_pC[result$cojo_selected != 'selected']");
            pw.println("non_sel_pC <- non_sel_pC[!is.na(non_sel_pC)]");
            pw.println("max_residual_pC <- if (length(non_sel_pC) > 0) min(non_sel_pC) else 1.0");
            pw.println("locus_fully_explained <- (n_selected == 0 && all(is.na(result$cojo_pC))) || (max_residual_pC >= p_cutoff)");
            pw.println("residual_snps <- character(0)");
            pw.println("if (!locus_fully_explained && length(non_sel_pC) > 0) {");
            pw.println("  res_idx <- which(result$cojo_selected != 'selected' & !is.na(result$cojo_pC) & result$cojo_pC < p_cutoff)");
            pw.println("  residual_snps <- result$snp_id[res_idx]");
            pw.println("  cat(sprintf('WARNING: Locus NOT fully explained — %d SNPs still significant after conditioning (min pC = %.2e)\\n', length(residual_snps), max_residual_pC))");
            pw.println("} else if (n_selected > 0) {");
            pw.println("  cat(sprintf('Locus fully explained: all conditional p >= %.0e (min pC = %.2e)\\n', p_cutoff, max_residual_pC))");
            pw.println("}");
            pw.println();

            // ═══ Artifact detection ═══
            pw.println("artifact_flags <- list()");
            pw.println("if (n_selected > 0) {");
            pw.println("  sel_idx <- which(result$cojo_selected == 'selected')");
            pw.println("  for (i in sel_idx) {");
            pw.println("    snp <- result$snp_id[i]; flags <- character(0)");
            pw.println("    bM <- result$cojo_b_marginal[i]; bJ <- result$cojo_bJ[i]");
            pw.println("    pM <- result$cojo_p_marginal[i]; pJ <- result$cojo_pJ[i]");
            // Effect inflation
            pw.println("    if (!is.na(bM) && !is.na(bJ) && abs(bM) > 0 && abs(bJ) > artifact_inflation_factor * abs(bM))");
            pw.println("      flags <- c(flags, sprintf('effect_inflation: |bJ|=%.4f > %dx|bM|=%.4f', abs(bJ), artifact_inflation_factor, abs(bM)))");
            // Sign flip with non-significant marginal
            pw.println("    if (!is.na(bM) && !is.na(bJ) && !is.na(pM) && sign(bM) != sign(bJ) && pM > 0.01)");
            pw.println("      flags <- c(flags, sprintf('sign_flip: bM=%.4f bJ=%.4f, marginal p=%.2e', bM, bJ, pM))");
            // Marginal-vs-joint inconsistency
            pw.println("    if (!is.na(pM) && !is.na(pJ) && pM > 1e-3 && pJ < 1e-20)");
            pw.println("      flags <- c(flags, sprintf('marginal_joint_inconsistency: pM=%.2e pJ=%.2e', pM, pJ))");
            pw.println("    if (length(flags) > 0) artifact_flags[[snp]] <- flags");
            pw.println("  }");
            pw.println("}");
            pw.println();

            // Reliability status
            pw.println("n_artifact_snps <- length(artifact_flags)");
            pw.println("cojo_reliability <- 'ok'");
            pw.println("reliability_reasons <- character(0)");
            pw.println("if (n_artifact_snps > 0) {");
            pw.println("  has_inflation <- any(sapply(artifact_flags, function(f) any(grepl('effect_inflation', f))))");
            pw.println("  has_inconsistency <- any(sapply(artifact_flags, function(f) any(grepl('marginal_joint_inconsistency', f))))");
            pw.println("  if (has_inflation || has_inconsistency) {");
            pw.println("    cojo_reliability <- 'likely_artifact'");
            pw.println("    reliability_reasons <- c(reliability_reasons, sprintf('%d SNP(s) with artifact flags', n_artifact_snps))");
            pw.println("  } else { cojo_reliability <- 'suspect'; reliability_reasons <- c(reliability_reasons, 'sign flips detected') }");
            pw.println("}");
            pw.println("if (ld_consistency_verdict == 'high_warn' && n_artifact_snps > 0) {");
            pw.println("  cojo_reliability <- 'likely_artifact'");
            pw.println("  reliability_reasons <- c(reliability_reasons, 'LD-GWAS consistency = high_warn + artifact flags — reference panel ancestry may not match GWAS')");
            pw.println("}");
            pw.println("if (nrow(freq_discrepancies) > 0) reliability_reasons <- c(reliability_reasons, sprintf('%d freq discrepancies', nrow(freq_discrepancies)))");
            pw.println("cat(sprintf('Reliability: %s\\n', cojo_reliability))");
            pw.println("if (length(reliability_reasons) > 0) cat(paste('  Reasons:', paste(reliability_reasons, collapse='; '), '\\n'))");
            pw.println("if (n_artifact_snps > 0) {");
            pw.println("  cat('Artifact details:\\n')");
            pw.println("  for (snp in names(artifact_flags)) cat(sprintf('  %s: %s\\n', snp, paste(artifact_flags[[snp]], collapse='; ')))");
            pw.println("}");
            pw.println();

            // ═══ Automatic safe fallback when artifacts detected ═══
            pw.println("safe_fallback_used <- FALSE");
            pw.println("safe_n_selected <- n_selected");
            pw.println("safe_selected_snps <- selected_snps");
            pw.println("if (cojo_reliability == 'likely_artifact' && collinear_thr > 0.5) {");
            pw.println("  cat('\\n=== ARTIFACT DETECTED — running safe fallback with collinear=0.5 ===\\n')");
            pw.println("  safe_prefix <- paste0(out_prefix, '_safe')");
            pw.println("  safe_args <- c('--bfile', bfile_prefix, '--cojo-file', ma_file, '--cojo-slct',");
            pw.println("    '--cojo-p', format(p_cutoff, scientific=TRUE), '--cojo-collinear', '0.5', '--out', safe_prefix)");
            pw.println("  safe_out <- system2(gcta_bin, args=safe_args, stdout=TRUE, stderr=TRUE)");
            pw.println("  safe_jma <- paste0(safe_prefix, '.jma.cojo')");
            pw.println("  if (file.exists(safe_jma) && file.info(safe_jma)$size > 0) {");
            pw.println("    safe_jma_data <- tryCatch(read.table(safe_jma, header=TRUE, stringsAsFactors=FALSE), error=function(e) data.frame())");
            pw.println("    if (nrow(safe_jma_data) > 0) {");
            pw.println("      safe_fallback_used <- TRUE");
            pw.println("      safe_n_selected <- nrow(safe_jma_data)");
            pw.println("      safe_selected_snps <- safe_jma_data$SNP");
            pw.println("      # Overwrite result with safe estimates");
            pw.println("      result$cojo_pJ <- NA_real_; result$cojo_bJ <- NA_real_; result$cojo_bJ_se <- NA_real_");
            pw.println("      result$cojo_selected <- ''");
            pw.println("      for (j in seq_len(nrow(safe_jma_data))) {");
            pw.println("        idx <- which(result$snp_id == safe_jma_data$SNP[j])");
            pw.println("        if (length(idx) > 0) {");
            pw.println("          result$cojo_pJ[idx] <- safe_jma_data$pJ[j]; result$cojo_bJ[idx] <- safe_jma_data$bJ[j]");
            pw.println("          result$cojo_bJ_se[idx] <- safe_jma_data$bJ_se[j]; result$cojo_selected[idx] <- 'selected'");
            pw.println("        }");
            pw.println("      }");
            pw.println("      # Re-run conditional pass with safe selected set");
            pw.println("      result$cojo_pC <- NA_real_; result$cojo_bC <- NA_real_; result$cojo_bC_se <- NA_real_");
            pw.println("      safe_cond_snp <- paste0(safe_prefix, '_sel.snps')");
            pw.println("      writeLines(safe_selected_snps, safe_cond_snp)");
            pw.println("      safe_cond_pfx <- paste0(safe_prefix, '_cond')");
            pw.println("      system2(gcta_bin, args=c('--bfile', bfile_prefix, '--cojo-file', ma_file,");
            pw.println("        '--cojo-cond', safe_cond_snp, '--out', safe_cond_pfx), stdout=TRUE, stderr=TRUE)");
            pw.println("      safe_cma <- paste0(safe_cond_pfx, '.cma.cojo')");
            pw.println("      if (file.exists(safe_cma)) {");
            pw.println("        sc <- tryCatch(read.table(safe_cma, header=TRUE, stringsAsFactors=FALSE), error=function(e) data.frame())");
            pw.println("        if (nrow(sc)>0) {");
            pw.println("          pc <- if ('pC' %in% colnames(sc)) 'pC' else 'p'; bc <- if ('bC' %in% colnames(sc)) 'bC' else 'b'; sec <- if ('bC_se' %in% colnames(sc)) 'bC_se' else 'se'");
            pw.println("          for (j in seq_len(nrow(sc))) { idx <- which(result$snp_id==sc$SNP[j]); if(length(idx)>0) { result$cojo_pC[idx]<-sc[[pc]][j]; result$cojo_bC[idx]<-sc[[bc]][j]; result$cojo_bC_se[idx]<-sc[[sec]][j] } }");
            pw.println("        }");
            pw.println("      }");
            pw.println("      for (s in safe_selected_snps) { idx<-which(result$snp_id==s); if(length(idx)>0){ result$cojo_pC[idx]<-result$cojo_pJ[idx]; result$cojo_bC[idx]<-result$cojo_bJ[idx]; result$cojo_bC_se[idx]<-result$cojo_bJ_se[idx] } }");
            pw.println("      # Recompute explained check");
            pw.println("      non_sel_pC <- result$cojo_pC[result$cojo_selected != 'selected']");
            pw.println("      non_sel_pC <- non_sel_pC[!is.na(non_sel_pC)]");
            pw.println("      max_residual_pC <- if (length(non_sel_pC) > 0) min(non_sel_pC) else 1.0");
            pw.println("      locus_fully_explained <- max_residual_pC >= p_cutoff");
            pw.println("      residual_snps <- character(0)");
            pw.println("      if (!locus_fully_explained) residual_snps <- result$snp_id[which(result$cojo_selected != 'selected' & !is.na(result$cojo_pC) & result$cojo_pC < p_cutoff)]");
            pw.println("      n_selected <- safe_n_selected; selected_snps <- safe_selected_snps");
            pw.println("      cat(sprintf('Safe fallback: %d signal(s), max|bJ/b|=%.1f\\n', safe_n_selected,");
            pw.println("        max(abs(safe_jma_data$bJ / safe_jma_data$b), na.rm=TRUE)))");
            pw.println("      reliability_reasons <- c(reliability_reasons, sprintf('safe fallback used (collinear=0.5): %d signal(s) — original had artifacts', safe_n_selected))");
            pw.println("    }");
            pw.println("  }");
            pw.println("}");
            pw.println();

            // Parse GCTA log warnings
            pw.println("gcta_warnings <- character(0); n_freq_badsnps <- 0L");
            pw.println("gcta_log_file <- paste0(out_prefix, '.log')");
            pw.println("if (file.exists(gcta_log_file)) {");
            pw.println("  log_lines <- readLines(gcta_log_file)");
            pw.println("  gcta_warnings <- grep('Warning|collinear|Error|dropped|removed', log_lines, value=TRUE, ignore.case=TRUE)");
            pw.println("  freq_line <- grep('large difference of allele frequency', log_lines, value=TRUE)");
            pw.println("  if (length(freq_line) > 0) { m <- regmatches(freq_line[1], regexpr('[0-9]+', freq_line[1])); if (length(m)>0) n_freq_badsnps <- as.integer(m) }");
            pw.println("}");
            pw.println("ldr_max_r2 <- NA_real_");
            pw.println("ldr_file <- paste0(out_prefix, '.ldr.cojo')");
            pw.println("if (file.exists(ldr_file)) { ldr <- tryCatch(read.table(ldr_file, header=TRUE, stringsAsFactors=FALSE), error=function(e) NULL)");
            pw.println("  if (!is.null(ldr) && nrow(ldr)>0 && 'r' %in% colnames(ldr)) ldr_max_r2 <- max(abs(ldr$r), na.rm=TRUE)^2 }");
            pw.println();

            // Write result.tsv
            pw.println("write.table(result, out_tsv, sep='\\t', quote=FALSE, row.names=FALSE)");
            pw.println();

            // Write manifest (marginal + joint + conditional clearly labeled)
            pw.println("v <- '2.2'");
            pw.println("m <- 'cojo_conditional'");
            pw.println("cd <- function(n,t) sprintf('{\"name\":\"%s\",\"type\":\"%s\",\"scope\":\"per_snp\",\"method\":\"%s\",\"method_version\":\"%s\"}', n, t, m, v)");
            pw.println("col_defs <- paste(c(cd('chr','string'), cd('pos','int'),");
            pw.println("  cd('cojo_p_marginal','double'), cd('cojo_b_marginal','double'),");
            pw.println("  cd('cojo_freq','double'), cd('cojo_n','int'),");
            pw.println("  cd('cojo_pJ','double'), cd('cojo_bJ','double'), cd('cojo_bJ_se','double'),");
            pw.println("  cd('cojo_pC','double'), cd('cojo_bC','double'), cd('cojo_bC_se','double'),");
            pw.println("  cd('cojo_selected','string')), collapse=',')");
            pw.println("manifest <- paste0('{\"schema_version\":\"1.0\",\"method\":\"', m, '\",\"method_version\":\"', v, '\",',");
            pw.println("  '\"parameters\":{},\"columns\":[', col_defs, '],\"created_at\":', as.numeric(Sys.time())*1000, '}')");
            pw.println("writeLines(manifest, manifest_path)");
            pw.println();

            // Write diagnostics JSON
            pw.println("esc <- function(s) gsub('\"', '\\\\\"', s)");
            pw.println("freq_disc_json <- '[]'");
            pw.println("if (nrow(freq_discrepancies) > 0) freq_disc_json <- paste0('[', paste(sprintf(");
            pw.println("  '{\"snp\":\"%s\",\"ma_freq\":%.4f,\"geno_freq\":%.4f,\"diff\":%.4f}',");
            pw.println("  freq_discrepancies$SNP, freq_discrepancies$freq, freq_discrepancies$freq_geno, freq_discrepancies$freq_diff), collapse=','), ']')");
            pw.println("artifact_json <- '{}'");
            pw.println("if (length(artifact_flags) > 0) {");
            pw.println("  af_entries <- sapply(names(artifact_flags), function(s) sprintf('\"%s\":[%s]', s, paste0('\"', esc(artifact_flags[[s]]), '\"', collapse=',')))");
            pw.println("  artifact_json <- paste0('{', paste(af_entries, collapse=','), '}')");
            pw.println("}");
            pw.println("residual_json <- if (length(residual_snps)>0) paste0('[', paste0('\"', residual_snps, '\"', collapse=','), ']') else '[]'");
            pw.println("warn_json <- if (length(gcta_warnings)>0) paste0('\"', esc(gcta_warnings), '\"', collapse=',') else ''");
            pw.println("reason_json <- if (length(reliability_reasons)>0) paste0('\"', esc(reliability_reasons), '\"', collapse=',') else ''");
            pw.println("diag <- sprintf(paste0(");
            pw.println("  '{\"method\":\"cojo_slct\",\"method_version\":\"2.2\",',");
            pw.println("  '\"gcta_binary\":\"%s\",\"gcta_version\":\"%s\",',");
            pw.println("  '\"n_selected\":%d,\"p_cutoff\":\"%s\",\"collinear_threshold\":%s,',");
            pw.println("  '\"sample_size\":%d,\"n_source\":\"%s\",\"freq_source\":\"%s\",',");
            pw.println("  '\"snps_in_ma\":%d,\"snps_dropped_no_ref_id\":%d,\"snps_dropped_no_freq\":%d,',");
            pw.println("  '\"snps_dropped_orientation\":%d,\"snps_freq_flipped\":%d,\"gcta_freq_badsnps\":%d,',");
            pw.println("  '\"ld_consistency_verdict\":\"%s\",\"ld_consistency_flagged\":%d,',");
            pw.println("  '\"locus_fully_explained\":%s,\"max_residual_pC\":%s,\"residual_signal_snps\":%s,',");
            pw.println("  '\"cojo_reliability\":\"%s\",\"reliability_reasons\":[%s],',");
            pw.println("  '\"artifact_flags\":%s,',");
            pw.println("  '\"gcta_exit_code\":%d,\"gcta_warnings\":[%s],',");
            pw.println("  '\"safe_fallback_used\":%s,',");
            pw.println("  '\"freq_discrepancies\":%s,\"ldr_max_r2\":%s,',");
            pw.println("  '\"selected_snps\":[%s]}'),");
            pw.println("  esc(gcta_bin), esc(gcta_ver),");
            pw.println("  n_selected, format(p_cutoff, scientific=TRUE), collinear_thr,");
            pw.println("  sample_size, n_source, freq_source,");
            pw.println("  snps_written, snps_dropped_no_ref_id, snps_dropped_no_freq,");
            pw.println("  snps_dropped_orientation, snps_freq_flipped, n_freq_badsnps,");
            pw.println("  ld_consistency_verdict, ld_consistency_flagged,");
            pw.println("  tolower(as.character(locus_fully_explained)),");
            pw.println("  sprintf('%.2e', max_residual_pC), residual_json,");
            pw.println("  cojo_reliability, reason_json,");
            pw.println("  artifact_json, cojo_exit, warn_json,");
            pw.println("  tolower(as.character(safe_fallback_used)),");
            pw.println("  freq_disc_json,");
            pw.println("  if (!is.na(ldr_max_r2)) sprintf('%.4f', ldr_max_r2) else 'null',");
            pw.println("  if (n_selected>0) paste0('\"', selected_snps, '\"', collapse=',') else '')");
            pw.println("writeLines(diag, diag_json)");
            pw.println();

            pw.println("cat(sprintf('COJO complete: %d signal(s), reliability=%s, explained=%s\\n', n_selected, cojo_reliability, locus_fully_explained))");
        }

        System.out.printf("[CojoAdapter] Prepared %d SNPs. R script at %s%n", written, rScript.getName());
        return written;
    }
}
