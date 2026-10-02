# LYNXgwas: joint fine-mapping of one aligned locus across several datasets with mvSuSiE.
#
#   Rscript mvsusie_run.R <manifest.json>
#
# The manifest (written by LocalServer) lists the datasets LYNXgwas selected for this locus (same
# disease, or colocalising datasets across diseases; datasets flagged for LD mismatch already
# excluded), each with its harmonised locus files and sample size, plus the window, reference
# panel, PLINK binary and optional UK Biobank LD. This script:
#   1. reads each dataset's z-scores in the window (beta/se), drops SNPs with per-SNP N < 80% of the
#      dataset's maximum (meta-analyses), and orients every dataset to the reference panel's A1
#      (palindromic SNPs dropped);
#   2. keeps datasets strongest-first while the SNPs shared by all stay >= max(300, 10% of the densest);
#   3. builds LD for the shared SNPs: UK Biobank when available (337k Europeans), else the panel;
#   4. removes SNPs whose z contradicts the LD in any dataset (susieR::kriging_rss) and sets the LD
#      regularisation from the largest remaining inconsistency (estimate_s_rss, >= 0.1, <= 0.5);
#   5. estimates the between-dataset correlation of z at null SNPs (sample overlap) and runs
#      mvSuSiE-RSS with it as residual covariance;
#   6. writes result.json: credible sets (size, PIP range, lead, purity, datasets each set is active
#      in, lfsr < 0.05), the correlation matrix, the LD source, and every exclusion with its reason.
args <- commandArgs(trailingOnly = TRUE)
suppressPackageStartupMessages(library(jsonlite))
`%||%` <- function(a, b) if (is.null(a)) b else a
man <- fromJSON(args[1], simplifyVector = FALSE)
if (nzchar(man$extra_lib %||% "")) .libPaths(c(man$extra_lib, .libPaths()))
need <- c("susieR", "mvsusieR", "bigsnpr", "data.table")
missing_pkgs <- need[!vapply(need, requireNamespace, TRUE, quietly = TRUE)]
out_json <- man$out
fail <- function(msg) { writeLines(toJSON(list(ok = FALSE, error = msg), auto_unbox = TRUE), out_json); cat("ERROR:", msg, "\n"); quit(status = 1) }
if (length(missing_pkgs)) fail(paste("R packages not installed:", paste(missing_pkgs, collapse = ", "),
                                    "(mvsusieR needs susieR >= 0.15.54; install both into the folder given by LYNXGWAS_R_LIB)"))
suppressPackageStartupMessages({ library(susieR); library(mvsusieR); library(bigsnpr); library(data.table) })

t_start <- Sys.time()
chr <- as.character(man$chr); lo <- as.numeric(man$from); hi <- as.numeric(man$to)
L <- man$L %||% 10; coverage <- man$coverage %||% 0.95
work <- man$work_dir
excluded <- list()
exclude <- function(id, reason) excluded[[length(excluded) + 1]] <<- list(dataset = id, reason = reason)
COMP <- c(A = "T", T = "A", C = "G", G = "C")

# ── 1. Reference panel window ────────────────────────────────────────────
win <- file.path(work, "win")
st <- system2(man$plink, c("--bfile", shQuote(man$ref_panel), "--chr", chr, "--from-bp", lo, "--to-bp", hi,
                           "--make-bed", "--memory", "4096", "--threads", "1", "--out", shQuote(win), "--silent"),
              stdout = TRUE, stderr = TRUE)
if (!file.exists(paste0(win, ".bim"))) fail("PLINK could not extract the reference panel window")
bim <- fread(paste0(win, ".bim"), header = FALSE, col.names = c("CHR", "SNP", "CM", "BP", "A1", "A2"), data.table = FALSE)
bim$key <- paste0(bim$CHR, ":", bim$BP)
dup <- bim$key[duplicated(bim$key)]
bim <- bim[!bim$key %in% dup, ]                                       # multi-allelic positions: ambiguous
bim <- bim[!(paste0(bim$A1, bim$A2) %in% c("AT", "TA", "CG", "GC")), ]  # palindromic: strand unknown
rownames(bim) <- bim$key

# ── 2. Each dataset's z-scores, oriented to the panel's A1 ───────────────
zs <- list(); ns <- c()
for (d in man$datasets) {
  frames <- lapply(unlist(d$files), function(f) if (file.exists(f)) fread(f, data.table = FALSE) else NULL)
  g <- do.call(rbind, frames[!vapply(frames, is.null, TRUE)])
  if (is.null(g) || !nrow(g)) { exclude(d$id, "no harmonised data in this window"); next }
  g <- g[as.character(g$chr) == chr & g$pos >= lo & g$pos <= hi, , drop = FALSE]
  g$beta <- suppressWarnings(as.numeric(g$beta)); g$se <- suppressWarnings(as.numeric(g$se))
  g <- g[is.finite(g$beta) & is.finite(g$se) & g$se > 0, , drop = FALSE]
  if ("n" %in% names(g)) {
    nn <- suppressWarnings(as.numeric(g$n))
    if (sum(is.finite(nn)) > 0.5 * nrow(g)) g <- g[is.finite(nn) & nn >= 0.8 * max(nn, na.rm = TRUE), , drop = FALSE]
  }
  g$key <- paste0(g$chr, ":", g$pos)
  g <- g[!duplicated(g$key) & g$key %in% bim$key, , drop = FALSE]
  b <- bim[g$key, ]
  ea <- toupper(g$ea); nea <- toupper(g$nea)
  same <- (ea == b$A1 & nea == b$A2) | (COMP[ea] == b$A1 & COMP[nea] == b$A2)
  swap <- (ea == b$A2 & nea == b$A1) | (COMP[ea] == b$A2 & COMP[nea] == b$A1)
  same[is.na(same)] <- FALSE; swap[is.na(swap)] <- FALSE
  z <- g$beta / g$se
  z <- ifelse(same, z, ifelse(swap, -z, NA))
  keep <- is.finite(z)
  if (sum(keep) < 50) { exclude(d$id, sprintf("only %d SNPs align with the reference panel", sum(keep))); next }
  zs[[d$id]] <- setNames(z[keep], g$key[keep]); ns[d$id] <- d$N
}

# Strongest association first (the datasets carrying the signal matter most), then add each dataset
# while the SNPs shared by all stay >= max(300, 10% of the densest dataset)
pv <- vapply(names(zs), function(id) { x <- Filter(function(d) d$id == id, man$datasets)[[1]]$p; if (is.null(x)) 1 else x }, 1)
ord <- names(zs)[order(pv, -vapply(zs, length, 1L))]
min_shared <- max(300, 0.1 * max(vapply(zs, length, 1L)))
kept <- character(0); shared <- NULL
for (id in ord) {
  nxt <- if (is.null(shared)) names(zs[[id]]) else intersect(shared, names(zs[[id]]))
  if (length(nxt) >= min_shared) { kept <- c(kept, id); shared <- nxt }
  else exclude(id, sprintf("too sparse here: would cut the shared SNP set to %d", length(nxt)))
}
if (length(kept) < 2) fail("fewer than 2 datasets have enough shared SNPs at this locus")
shared <- bim$key[bim$key %in% shared]                                  # panel order (by position)
Z <- sapply(kept, function(id) zs[[id]][shared]); rownames(Z) <- shared
N <- ns[kept]

# ── 3. LD: UK Biobank when available, else the reference panel ───────────
R <- NULL; ld_source <- "1000 Genomes reference panel"
py <- Sys.getenv("LYNXGWAS_PYTHON"); if (!nzchar(py)) py <- Sys.which("python"); if (!nzchar(py)) py <- Sys.which("python3")
if (nzchar(man$ukbb_dir %||% "") && dir.exists(man$ukbb_dir) && file.exists(man$ukbb_script %||% "") && nzchar(py)) {
  snp_in <- file.path(work, "ukb_snps.tsv"); out_pre <- file.path(work, "ukb")
  b <- bim[shared, ]
  write.table(data.frame(chr = b$CHR, pos = b$BP, a1 = b$A1, a2 = b$A2), snp_in, sep = "\t", quote = FALSE, row.names = FALSE)
  msg <- suppressWarnings(system2(py, c(shQuote(man$ukbb_script), shQuote(man$ukbb_dir), shQuote(snp_in), shQuote(out_pre)), stdout = TRUE, stderr = TRUE))
  if (is.null(attr(msg, "status")) && file.exists(paste0(out_pre, ".ld"))) {
    idx <- scan(paste0(out_pre, ".idx"), quiet = TRUE)
    if (length(idx) >= 50) {
      R <- as.matrix(read.table(paste0(out_pre, ".ld"))); dimnames(R) <- NULL
      Z <- Z[idx, , drop = FALSE]; shared <- shared[idx]
      ld_source <- sprintf("UK Biobank (%s)", paste(msg, collapse = " "))
    }
  }
}
if (is.null(R)) {
  writeLines(bim[shared, "SNP"], file.path(work, "extract.txt"))
  sub <- file.path(work, "ref")
  system2(man$plink, c("--bfile", shQuote(win), "--extract", shQuote(file.path(work, "extract.txt")), "--make-bed",
                       "--memory", "2048", "--threads", "1", "--out", shQuote(sub), "--silent"), stdout = TRUE, stderr = TRUE)
  unlink(file.path(work, c("ref_bigsnp.bk", "ref_bigsnp.rds")))
  rds <- snp_readBed(paste0(sub, ".bed"), backingfile = file.path(work, "ref_bigsnp"))
  obj <- snp_attach(rds)
  ix <- match(bim[shared, "SNP"], obj$map$marker.ID)
  maf <- snp_MAF(obj$genotypes, ind.col = ix)
  ok <- is.finite(maf) & maf > 0
  Z <- Z[ok, , drop = FALSE]; shared <- shared[ok]; ix <- ix[ok]
  R <- as.matrix(snp_cor(obj$genotypes, ind.col = ix, size = length(ix), ncores = 1L))
}
R <- (R + t(R)) / 2; R[!is.finite(R)] <- 0; diag(R) <- 1

# ── 4. LD consistency check across datasets ──────────────────────────────
removed <- 0L
for (it in 1:3) {
  bad <- integer(0)
  for (k in seq_len(ncol(Z))) {
    kr <- tryCatch(kriging_rss(Z[, k], R, N[k])$conditional_dist, error = function(e) NULL)
    if (!is.null(kr)) bad <- union(bad, which(kr$logLR > 2 & abs(kr$z) > 2))
  }
  if (!length(bad)) break
  bad <- head(bad, max(1L, ceiling(0.05 * nrow(Z))))
  keep <- setdiff(seq_len(nrow(Z)), bad)
  Z <- Z[keep, , drop = FALSE]; R <- R[keep, keep, drop = FALSE]; shared <- shared[keep]; removed <- removed + length(bad)
}
s_by <- vapply(seq_len(ncol(Z)), function(k) tryCatch(estimate_s_rss(Z[, k], R, N[k]), error = function(e) NA_real_), 0)
names(s_by) <- colnames(Z)
lambda <- min(0.5, max(c(man$ld_shrink %||% 0.1, s_by), na.rm = TRUE))
Rs <- (1 - lambda) * R + lambda * diag(nrow(R))

# ── 5. Between-dataset correlation (sample overlap) and mvSuSiE ─────────
K <- ncol(Z)
# Estimated genome-wide (mvSuSiE paper): z at a thinned set of SNPs on the other chromosomes with
# |z| < 2 in both datasets. Inside the locus window the shared true signal leaks into "null" SNPs and
# inflates the estimate, which makes mvSuSiE split one signal into dataset-specific sets. Each
# dataset's thinned z table is cached next to its project (rebuilt when the GWAS file changes).
null_z <- function(d) {
  g <- d$gwas; if (is.null(g) || !nzchar(g$file %||% "") || !file.exists(g$file)) return(NULL)
  cache <- d$null_cache
  if (!is.null(cache) && file.exists(cache) && file.mtime(cache) > file.mtime(g$file)) return(readRDS(cache))
  cols <- c(g$chr, g$pos, g$ea, g$nea, g$se, if (nzchar(g$beta %||% "")) g$beta else g$or)
  x <- tryCatch(fread(g$file, select = cols, data.table = FALSE, showProgress = FALSE), error = function(e) NULL)
  if (is.null(x) || !nrow(x)) return(NULL)
  pos <- suppressWarnings(as.numeric(x[[g$pos]]))
  keep <- is.finite(pos) & pos %% 53 == 7                                   # ~2% of SNPs, same ones in every dataset
  x <- x[keep, , drop = FALSE]; pos <- pos[keep]
  eff <- suppressWarnings(as.numeric(x[[cols[6]]])); if (!nzchar(g$beta %||% "")) eff <- log(eff)
  z <- eff / suppressWarnings(as.numeric(x[[g$se]]))
  ea <- toupper(x[[g$ea]]); nea <- toupper(x[[g$nea]])
  ok <- is.finite(z) & nchar(ea) == 1 & nchar(nea) == 1 & !(paste0(ea, nea) %in% c("AT", "TA", "CG", "GC"))
  # orient to the alphabetically first allele, strand-collapsed (A/T -> A, C/G -> C)
  canon <- c(A = "A", T = "A", C = "C", G = "C")
  flip <- canon[ea] > canon[nea]
  ch <- sub("^chr", "", as.character(x[[g$chr]]), ignore.case = TRUE)
  out <- data.frame(chr = ch[ok], key = paste0(ch, ":", pos)[ok], z = ifelse(flip, -z, z)[ok], stringsAsFactors = FALSE)
  out <- out[!duplicated(out$key) & !out$key %in% out$key[duplicated(out$key)], ]
  # write-then-rename: parallel runs may build the same dataset's cache at once
  if (!is.null(cache)) tryCatch({
    dir.create(dirname(cache), showWarnings = FALSE, recursive = TRUE)
    tmp <- paste0(cache, ".", Sys.getpid(), ".tmp"); saveRDS(out, tmp)
    if (!file.rename(tmp, cache)) unlink(tmp)
  }, error = function(e) NULL)
  out
}
ids <- colnames(Z); nz <- list()
for (d in man$datasets) if (d$id %in% ids) {
  t <- null_z(d)
  if (!is.null(t)) { t <- t[t$chr != chr & abs(t$z) < 2, ]; nz[[d$id]] <- setNames(t$z, t$key) }
}
V <- diag(K); n_pair <- matrix(0L, K, K); v_source <- "genome-wide null SNPs"
for (a in seq_len(K - 1)) for (b in (a + 1):K) {
  A <- nz[[ids[a]]]; B <- nz[[ids[b]]]
  if (!is.null(A) && !is.null(B)) {
    k <- intersect(names(A), names(B)); n_pair[a, b] <- n_pair[b, a] <- length(k)
    if (length(k) >= 1000) { V[a, b] <- V[b, a] <- cor(A[k], B[k]); next }
  }
  # fallback: SNPs in this window that are null in both (less reliable, see above)
  k <- intersect(names(zs[[ids[a]]]), names(zs[[ids[b]]]))
  za <- zs[[ids[a]]][k]; zb <- zs[[ids[b]]][k]; nul <- abs(za) < 2 & abs(zb) < 2
  n_pair[a, b] <- n_pair[b, a] <- sum(nul); v_source <- "genome-wide null SNPs (some pairs: locus window only)"
  if (sum(nul) >= 100) V[a, b] <- V[b, a] <- cor(za[nul], zb[nul])
}
V[!is.finite(V)] <- 0; V <- (V + t(V)) / 2; diag(V) <- 1
V_raw <- V; null <- if (K > 1) min(n_pair[upper.tri(n_pair)]) else 0
ev <- eigen(V, symmetric = TRUE); V <- cov2cor(ev$vectors %*% diag(pmax(ev$values, 0.05)) %*% t(ev$vectors))
dimnames(V) <- list(colnames(Z), colnames(Z))

fit <- mvsusie_rss(Z = Z, R = Rs, N = round(median(N)), prior_variance = create_mixture_prior(R = K),
                   residual_variance = V, L = L, coverage = coverage, estimate_prior_variance = TRUE, verbose = FALSE)

# ── 6. Result ─────────────────────────────────────────────────────────────
snp_id <- bim[shared, "SNP"]; snp_pos <- bim[shared, "BP"]
pip <- as.numeric(fit$pip)
cs <- if (!is.null(fit$sets$cs)) fit$sets$cs else list()
purity <- fit$sets$purity
lf <- fit$single_effect_lfsr
sets <- lapply(seq_along(cs), function(i) {
  m <- cs[[i]]; l <- as.integer(sub("L", "", names(cs)[i]))
  lead <- m[which.max(pip[m])]
  active <- if (!is.null(lf)) colnames(Z)[lf[l, ] < 0.05] else character(0)
  list(set = i, size = length(m), pip_min = round(min(pip[m]), 4), pip_max = round(max(pip[m]), 4),
       lead_snp = snp_id[lead], lead_pos = snp_pos[lead], lead_pip = round(pip[lead], 4),
       min_abs_corr = if (!is.null(purity)) round(purity[i, 1], 3) else NA,
       active_in = as.list(active),
       snps = lapply(m[order(-pip[m])], function(j) list(id = snp_id[j], pos = snp_pos[j], pip = round(pip[j], 4))))
})
res <- list(ok = TRUE, chr = chr, from = lo, to = hi, n_snps = nrow(Z), datasets = as.list(colnames(Z)),
            N = as.list(N), ld_source = ld_source, ld_regularisation = round(lambda, 3),
            ld_inconsistency = as.list(round(s_by, 3)), snps_removed_ld_check = removed,
            null_correlation = lapply(seq_len(K), function(i) unname(round(V_raw[i, ], 3))), n_null_snps = null, null_correlation_source = v_source,
            credible_sets = sets, max_pip = round(max(pip), 4), excluded = excluded,
            snp_pips = { k <- which(pip >= 0.001); list(id = as.list(snp_id[k]), pos = as.list(snp_pos[k]), pip = as.list(round(pip[k], 4))) },
            seconds = round(as.numeric(difftime(Sys.time(), t_start, units = "secs")), 1))
writeLines(toJSON(res, auto_unbox = TRUE, null = "null", na = "null", digits = NA), out_json)
cat(sprintf("mvSuSiE: %d datasets, %d SNPs, %d credible sets, LD: %s\n", K, nrow(Z), length(sets), ld_source))
