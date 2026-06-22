
import java.io.*;
import java.util.*;

/**
 * COJO (GCTA) adapter: prepares the matched GWAS in GCTA .ma format,
 * generates the R script that runs iterative COJO (from precodes pattern),
 * and maps output to the unified contract.
 */
public class CojoAdapter {

    /**
     * Write GCTA .ma input AND the R runner script.
     */
    public static int prepareRun(File harmonizedDir, File matchedDir, File runDir,
                                  int sampleN, double pCutoff, String gctaBin) throws IOException {
        runDir.mkdirs();

        // ── Write .ma file from harmonized GWAS ──
        File gwasFile = new File(harmonizedDir, "harmonized_gwas.tsv");
        if (!gwasFile.exists()) return 0;

        File maFile = new File(runDir, "cojo_input.ma");
        // Also write a matched TSV in the column format the R script expects
        File matchedTsv = new File(runDir, "matched_gwas_cojo.tsv");
        int count = 0;

        // Read matched .bim for SNP IDs that GCTA can reference
        Map<String, String> posToRefId = new LinkedHashMap<>();
        File bimFile = new File(matchedDir, "matched_ref.bim");
        if (bimFile.exists()) {
            try (BufferedReader br = new BufferedReader(new FileReader(bimFile))) {
                String line;
                while ((line = br.readLine()) != null) {
                    String[] f = line.split("\t", -1);
                    if (f.length < 6) continue;
                    String chr = f[0].replaceFirst("^chr", "");
                    posToRefId.put(chr + ":" + f[3].trim(), f[1]);
                }
            }
        }

        try (BufferedReader br = new BufferedReader(new FileReader(gwasFile));
             PrintWriter ma = new PrintWriter(new BufferedWriter(new FileWriter(maFile)));
             PrintWriter tsv = new PrintWriter(new BufferedWriter(new FileWriter(matchedTsv)))) {

            ma.println("SNP\tA1\tA2\tfreq\tb\tse\tp\tN");
            tsv.println("SNP\tA1\tA2\tBETA\tSE\tP\tBP\tMAF");

            String header = br.readLine();
            if (header == null) return 0;

            String line;
            while ((line = br.readLine()) != null) {
                if (line.isEmpty()) continue;
                String[] f = line.split("\t", -1);
                if (f.length < 8) continue;

                String snpId = f[0];
                String chr   = f[1];
                String pos   = f[2];
                String ea    = f[3];
                String nea   = f[4];
                String pval  = f[5];
                String beta  = f[6];
                String se    = f[7];
                String maf   = f.length > 10 ? f[10] : "NA";

                if (beta.equals("NA") || se.equals("NA")) continue;
                double betaVal, seVal;
                try {
                    betaVal = Double.parseDouble(beta);
                    seVal = Double.parseDouble(se);
                } catch (NumberFormatException e) { continue; }
                if (Double.isNaN(betaVal) || Double.isNaN(seVal) || seVal <= 0) continue;

                double freq = 0.2;
                if (!maf.equals("NA")) {
                    try { freq = Double.parseDouble(maf); } catch (NumberFormatException e) {}
                    if (freq <= 0 || freq >= 1) freq = 0.2;
                }

                // Use ref panel SNP ID if available (GCTA matches against bfile)
                String refId = posToRefId.get(chr + ":" + pos);
                String useId = refId != null ? refId : snpId;

                ma.printf("%s\t%s\t%s\t%.4f\t%.6g\t%.6g\t%s\t%d%n",
                    useId, ea, nea, freq, betaVal, seVal, pval, sampleN);
                tsv.printf("%s\t%s\t%s\t%.6g\t%.6g\t%s\t%s\t%.4f%n",
                    useId, ea, nea, betaVal, seVal, pval, pos, freq);
                count++;
            }
        }

        // ── Write R runner script (based on precodes/scripts/run_cojo_iterative_gcta.R) ──
        String bfilePrefix = new File(matchedDir, "matched_ref").getAbsolutePath().replace("\\", "/");
        String maPath = maFile.getAbsolutePath().replace("\\", "/");
        String outTsv = new File(runDir, "result.tsv").getAbsolutePath().replace("\\", "/");
        String diagJson = new File(runDir, "cojo_diag.json").getAbsolutePath().replace("\\", "/");

        // Resolve GCTA binary to absolute path
        File gctaFile = new File(gctaBin);
        String gctaAbsolute = gctaFile.isAbsolute() ? gctaFile.getAbsolutePath()
            : new File(gctaBin).getAbsoluteFile().getAbsolutePath();

        File rScript = new File(runDir, "run_cojo.R");
        try (PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(rScript)))) {
            pw.println("#!/usr/bin/env Rscript");
            pw.printf("gcta_bin <- '%s'%n", gctaAbsolute.replace("\\", "/"));
            pw.printf("ma_file <- '%s'%n", maPath);
            pw.printf("bfile_prefix <- '%s'%n", bfilePrefix);
            pw.printf("p_cutoff <- %s%n", String.valueOf(pCutoff));
            pw.printf("sample_size <- %d%n", sampleN);
            pw.printf("out_tsv <- '%s'%n", outTsv);
            pw.printf("diag_json <- '%s'%n", diagJson);
            pw.println();
            pw.println("tmp_dir <- file.path(dirname(out_tsv), 'gcta_tmp')");
            pw.println("dir.create(tmp_dir, recursive=TRUE, showWarnings=FALSE)");
            pw.println("out_prefix <- file.path(tmp_dir, 'slct')");
            pw.println();
            pw.println("# Run GCTA --cojo-slct");
            pw.println("cmd_args <- c('--bfile', bfile_prefix,");
            pw.println("  '--cojo-file', ma_file,");
            pw.println("  '--cojo-p', format(p_cutoff, scientific=TRUE),");
            pw.println("  '--cojo-slct', '--cojo-collinear', '0.9',");
            pw.println("  '--out', out_prefix)");
            pw.println();
            pw.println("status <- system2(gcta_bin, args=cmd_args, stdout=TRUE, stderr=TRUE)");
            pw.println("status_text <- paste(status, collapse='\\n')");
            pw.println();
            pw.println("jma_file <- paste0(out_prefix, '.jma.cojo')");
            pw.println("if (!file.exists(jma_file)) {");
            pw.println("  # No SNPs selected — write empty result");
            pw.println("  writeLines('snp_id\\tcojo_pJ\\tcojo_bJ\\tcojo_bJ_se\\tcojo_selected', out_tsv)");
            pw.println("  writeLines('{\"method\":\"cojo\",\"n_selected\":0,\"status\":\"no_snps_selected\"}', diag_json)");
            pw.println("  cat('COJO: No SNPs passed p-value threshold\\n')");
            pw.println("  quit(save='no', status=0)");
            pw.println("}");
            pw.println();
            pw.println("jma <- read.table(jma_file, header=TRUE, stringsAsFactors=FALSE)");
            pw.println("selected_snps <- jma$SNP");
            pw.println();
            pw.println("# Also run --cojo-cond to get conditional p for ALL snps");
            pw.println("if (length(selected_snps) > 0) {");
            pw.println("  cond_file <- file.path(tmp_dir, 'cond.snps')");
            pw.println("  writeLines(selected_snps, cond_file)");
            pw.println("  cond_prefix <- file.path(tmp_dir, 'cond')");
            pw.println("  system2(gcta_bin, args=c('--bfile', bfile_prefix, '--cojo-file', ma_file,");
            pw.println("    '--cojo-cond', cond_file, '--out', cond_prefix), stdout=TRUE, stderr=TRUE)");
            pw.println("  cma_file <- paste0(cond_prefix, '.cma.cojo')");
            pw.println("  if (file.exists(cma_file)) {");
            pw.println("    cma <- read.table(cma_file, header=TRUE, stringsAsFactors=FALSE)");
            pw.println("  } else { cma <- NULL }");
            pw.println("} else { cma <- NULL }");
            pw.println();
            pw.println("# Build unified result");
            pw.println("ma_data <- read.table(ma_file, header=TRUE, stringsAsFactors=FALSE)");
            pw.println("result <- data.frame(snp_id=ma_data$SNP, stringsAsFactors=FALSE)");
            pw.println("result$cojo_pJ <- NA_real_");
            pw.println("result$cojo_bJ <- NA_real_");
            pw.println("result$cojo_bJ_se <- NA_real_");
            pw.println("result$cojo_selected <- FALSE");
            pw.println();
            pw.println("# Fill from jma (selected SNPs)");
            pw.println("for (i in seq_len(nrow(jma))) {");
            pw.println("  idx <- which(result$snp_id == jma$SNP[i])");
            pw.println("  if (length(idx) > 0) {");
            pw.println("    result$cojo_pJ[idx] <- jma$pJ[i]");
            pw.println("    result$cojo_bJ[idx] <- jma$bJ[i]");
            pw.println("    result$cojo_bJ_se[idx] <- jma$bJ_se[i]");
            pw.println("    result$cojo_selected[idx] <- TRUE");
            pw.println("  }");
            pw.println("}");
            pw.println();
            pw.println("# Fill conditional p-values from cma");
            pw.println("if (!is.null(cma)) {");
            pw.println("  p_col <- if ('pC' %in% colnames(cma)) 'pC' else 'p'");
            pw.println("  for (i in seq_len(nrow(cma))) {");
            pw.println("    idx <- which(result$snp_id == cma$SNP[i])");
            pw.println("    if (length(idx) > 0 && is.na(result$cojo_pJ[idx])) {");
            pw.println("      result$cojo_pJ[idx] <- cma[[p_col]][i]");
            pw.println("    }");
            pw.println("  }");
            pw.println("}");
            pw.println();
            pw.println("write.table(result, out_tsv, sep='\\t', quote=FALSE, row.names=FALSE)");
            pw.println();
            pw.println("# Write result manifest");
            pw.printf("manifest_path <- '%s'%n",
                new File(runDir, "result.manifest.json").getAbsolutePath().replace("\\", "/"));
            pw.println("manifest <- paste0('{\"schema_version\":\"1.0\",\"method\":\"cojo_conditional\",\"method_version\":\"1.0\",');");
            pw.println("manifest <- paste0(manifest, '\"parameters\":{},\"columns\":[');");
            pw.println("manifest <- paste0(manifest, '{\"name\":\"cojo_pJ\",\"type\":\"double\",\"scope\":\"per_snp\",\"method\":\"cojo_conditional\",\"method_version\":\"1.0\"},');");
            pw.println("manifest <- paste0(manifest, '{\"name\":\"cojo_bJ\",\"type\":\"double\",\"scope\":\"per_snp\",\"method\":\"cojo_conditional\",\"method_version\":\"1.0\"},');");
            pw.println("manifest <- paste0(manifest, '{\"name\":\"cojo_bJ_se\",\"type\":\"double\",\"scope\":\"per_snp\",\"method\":\"cojo_conditional\",\"method_version\":\"1.0\"},');");
            pw.println("manifest <- paste0(manifest, '{\"name\":\"cojo_selected\",\"type\":\"boolean\",\"scope\":\"per_snp\",\"method\":\"cojo_conditional\",\"method_version\":\"1.0\"}');");
            pw.println("manifest <- paste0(manifest, '],\"created_at\":', as.numeric(Sys.time())*1000, '}');");
            pw.println("writeLines(manifest, manifest_path)");
            pw.println();
            pw.println("# Diagnostics");
            pw.printf("writeLines(sprintf('{\"method\":\"cojo\",\"n_selected\":%%d,\"p_cutoff\":%%s,\"sample_size\":%%d}',");
            pw.println(" length(selected_snps), format(p_cutoff, scientific=TRUE), sample_size), diag_json)");
            pw.println();
            pw.println("cat(sprintf('COJO: %%d SNPs selected at p < %%s\\n', length(selected_snps), format(p_cutoff, scientific=TRUE)))");
        }

        System.out.printf("[CojoAdapter] Prepared %d SNPs, R script at %s%n", count, rScript.getName());
        return count;
    }
}
