#!/usr/bin/env Rscript
# post_freeze_exploratory
#
# Refits the three resistance/recurrence differential-expression models with
# IDENTICAL parameters to the frozen analyses, and exports the per-gene
# quantities that the frozen output tables do not carry: a model-based
# standard error for the log2 fold change, the total degrees of freedom the
# test actually used, the gene-wise dispersion, and the QL variance scaling.
#
# Nothing here changes a frozen result. The refits are verified against the
# committed frozen tables gene-by-gene (log2fc and p-value) and the run aborts
# if any dataset fails to reproduce.
#
# Data sources and versions (all already in this repository, no downloads):
#   GSE118713  data/processed/gse118713_gene_tpm_filtered.tsv.gz (config
#              gse118713_phase2b.filtering.filtered_gene_tpm_tsv); limma
#              3.66.0, lmFit/contrasts.fit/eBayes on log2(TPM + 1), design
#              ~0 + group, contrast TAMR - MCF7. Verified against the
#              UNREDACTED frozen table (14,838 genes), not the published
#              redacted one (14,836): KDM1A and RCOR1 were withheld from the
#              redacted table by the preregistered blinding, which was
#              retired on 2026-08-10 (CLAUDE.md; PREANALYSIS.md amendments
#              log). Both genes were always fitted; only their reporting was
#              withheld, so using the unredacted table changes no estimate.
#              KDM1A is one of the 13 candidates, so this matters here.
#   GSE111151  results/tables/gse111151/counts_matrix.tsv.gz; edgeR 4.8.2,
#              filterByExpr + TMM + estimateDisp + glmQLFit/glmQLFTest,
#              design ~cell_line + resistance_status.
#   GSE240112  results/tables/gse240112_pseudobulk/tumor_cell_counts.tsv.gz;
#              edgeR 4.8.2, same workflow, design ~group (RT vs PT).
#
# WHAT THE STANDARD ERROR IS, AND IS NOT
# --------------------------------------
# edgeR's glmQLFTest is a deviance-based quasi-likelihood F test, not a Wald
# test, so no standard error can be recovered exactly by inverting its
# p-value. What is exported instead is the model-based Wald standard error
# implied by the same fit:
#
#     Var(beta_hat) = var.post * (X' W X)^-1 ,   W = mu / (1 + phi * mu)
#
# on the natural-log scale, divided by log(2) for log2. This is exact
# CONDITIONAL on the fitted means, the gene-wise dispersion and the posterior
# QL variance -- it is not an exact inversion of glmQLFTest. The script
# therefore reports, per dataset, how closely t_wald^2 tracks the reported
# quasi-F statistic, so a mismatched covariance reconstruction cannot pass
# silently. Downstream consumers must treat the Wald p-value as a diagnostic,
# not as a replacement for the frozen p-value.
#
# Output: results/post_poster/de_refit/*.tsv.gz and refit_validation.tsv

suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
  library(yaml)
  library(data.table)
})

args <- commandArgs(trailingOnly = FALSE)
repo_root <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", args, value = TRUE))), ".."))
if (length(repo_root) == 0 || !dir.exists(file.path(repo_root, "config"))) repo_root <- getwd()

config <- yaml::read_yaml(file.path(repo_root, "config", "config.yaml"))
OUT_DIR <- file.path(repo_root, "results", "post_poster", "de_refit")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

validation <- list()

log_rows <- function(label, stage, n_in, n_out) {
  message(sprintf("[%s] %-38s rows in = %7d  rows out = %7d  lost = %6d",
                  label, stage, n_in, n_out, n_in - n_out))
}

# ---------------------------------------------------------------------------
# Model-based Wald standard errors from a fitted glmQLFit object.
#
# For a negative-binomial GLM with log link, no observation weights and a
# per-gene dispersion phi, the IRLS working weight for observation i is
# mu_i / (1 + phi * mu_i). These are the same weights edgeR itself uses when
# fitting, so (X' W X)^-1 is the fit's own unscaled coefficient covariance.
# ---------------------------------------------------------------------------
wald_se_from_qlfit <- function(fit, coef_name) {
  design <- fit$design
  j <- match(coef_name, colnames(design))
  stopifnot(!is.na(j))
  mu <- fit$fitted.values
  phi <- fit$dispersion
  if (length(phi) == 1L) phi <- rep(phi, nrow(mu))
  # edgeR >= 4.4 renamed the QL variance scaling var.post -> s2.post; accept
  # either rather than silently producing an unscaled SE.
  s2_post <- if (!is.null(fit$s2.post)) fit$s2.post else fit$var.post
  if (is.null(s2_post)) stop("fit has neither s2.post nor var.post -- not a glmQLFit object")

  n_genes <- nrow(mu)
  se_nat <- numeric(n_genes)
  for (g in seq_len(n_genes)) {
    w <- mu[g, ] / (1 + phi[g] * mu[g, ])
    xtwx <- crossprod(design, design * w)
    inv <- tryCatch(solve(xtwx), error = function(e) NULL)
    se_nat[g] <- if (is.null(inv)) NA_real_ else sqrt(s2_post[g] * inv[j, j])
  }
  se_nat / log(2)
}

# The residual df glmQLFTest actually used: df.residual.zeros when the fit
# carries it, otherwise df.residual.adj (edgeR >= 4.4 QL pipeline). Mirrors
# the branch inside glmQLFTest rather than assuming one of them.
ql_df_residual <- function(fit) {
  if (!is.null(fit$df.residual.zeros)) fit$df.residual.zeros else fit$df.residual.adj
}

run_edger_refit <- function(label, counts_path, metadata_path, gene_col,
                            design_fn, coef_pattern, frozen_path, frozen_gene_col) {
  message(sprintf("=== %s ===", label))
  counts_dt <- read.delim(gzfile(counts_path), header = TRUE, sep = "\t",
                          stringsAsFactors = FALSE, check.names = FALSE)
  id_cols <- intersect(c("gene_id", "gene_name", "gene"), colnames(counts_dt))
  gene_ids <- counts_dt[[gene_col]]
  counts <- as.matrix(counts_dt[, setdiff(colnames(counts_dt), id_cols)])
  rownames(counts) <- gene_ids
  gene_name_map <- if ("gene_name" %in% colnames(counts_dt)) {
    setNames(counts_dt$gene_name, counts_dt[[gene_col]])
  } else {
    setNames(counts_dt[[gene_col]], counts_dt[[gene_col]])
  }

  metadata <- read.delim(metadata_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  stopifnot(identical(colnames(counts), metadata$sample_id))
  design <- design_fn(metadata)

  dge <- DGEList(counts = counts, samples = metadata)
  n_before <- nrow(dge)
  keep <- filterByExpr(dge, design = design)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  log_rows(label, "filterByExpr", n_before, nrow(dge))
  dge <- calcNormFactors(dge, method = "TMM")
  dge <- estimateDisp(dge, design)
  fit <- glmQLFit(dge, design)
  coef_name <- grep(coef_pattern, colnames(design), value = TRUE)
  stopifnot(length(coef_name) == 1)
  qlf <- glmQLFTest(fit, coef = coef_name)
  tab <- topTags(qlf, n = Inf, sort.by = "none")$table

  # Take df.total straight from the test object -- that is, by construction,
  # the df the reported p-value was computed against; do not re-derive it.
  df_total <- qlf$df.total
  if (length(df_total) == 1L) df_total <- rep(df_total, nrow(dge))
  df_res <- ql_df_residual(fit)
  if (length(df_res) == 1L) df_res <- rep(df_res, nrow(dge))
  s2_post <- if (!is.null(fit$s2.post)) fit$s2.post else fit$var.post
  se_log2 <- wald_se_from_qlfit(fit, coef_name)
  disp <- if (length(fit$dispersion) == 1L) rep(fit$dispersion, nrow(dge)) else fit$dispersion

  out <- data.table(
    gene = rownames(dge),
    gene_symbol = unname(gene_name_map[rownames(dge)]),
    log2fc = tab$logFC,
    se_log2fc_wald = se_log2,
    df_total = df_total,
    df_prior = if (length(fit$df.prior) == 1L) rep(fit$df.prior, nrow(dge)) else fit$df.prior,
    df_residual = df_res,
    dispersion = disp,
    var_post = s2_post,
    avg_log_cpm = aveLogCPM(dge),
    quasi_f = tab$F,
    p_value_frozen_method = tab$PValue,
    fdr_frozen_method = tab$FDR,
    engine = "edgeR_glmQLF",
    post_freeze_exploratory = TRUE
  )
  out[, t_wald := log2fc / se_log2fc_wald]
  out[, p_value_wald := 2 * pt(-abs(t_wald), df = df_total)]

  # --- verification against the committed frozen table -----------------
  frozen <- fread(cmd = paste("zcat", shQuote(frozen_path)))
  setnames(frozen, frozen_gene_col, "gene")
  m <- merge(out[, .(gene, log2fc, p_value_frozen_method)],
             frozen[, .(gene, log2fc_frozen = log2fc, p_frozen = p_value)],
             by = "gene")
  log_rows(label, "join to frozen table", nrow(out), nrow(m))
  stopifnot(nrow(m) == nrow(out), nrow(m) == nrow(frozen))
  max_dlfc <- max(abs(m$log2fc - m$log2fc_frozen))
  max_dp <- max(abs(m$p_value_frozen_method - m$p_frozen))
  message(sprintf("  reproduction: max |dlog2fc| = %.3e, max |dp| = %.3e", max_dlfc, max_dp))
  if (max_dlfc > 1e-8 || max_dp > 1e-10) {
    stop(sprintf("%s refit does NOT reproduce the frozen table -- refusing to write", label))
  }

  # --- diagnostic: does t_wald^2 track the reported quasi-F? ------------
  ok <- is.finite(out$t_wald) & is.finite(out$quasi_f) & out$quasi_f > 0
  ratio <- (out$t_wald[ok]^2) / out$quasi_f[ok]
  message(sprintf("  t_wald^2 / quasi_F: median = %.4f, IQR = [%.4f, %.4f], n = %d",
                  median(ratio), quantile(ratio, 0.25), quantile(ratio, 0.75), sum(ok)))

  validation[[label]] <<- data.table(
    dataset = label, engine = "edgeR_glmQLF", n_genes = nrow(out),
    # Row accounting persisted, not merely logged: how many genes the counts
    # matrix carried, how many filterByExpr removed, and how many joined the
    # frozen table.
    n_genes_before_filter = n_before,
    n_genes_removed_by_filter = n_before - nrow(dge),
    n_genes_after_filter = nrow(dge),
    n_genes_joined_to_frozen = nrow(m),
    n_samples = ncol(dge), design = paste(colnames(design), collapse = " + "),
    coef = coef_name, df_residual = df_res[1],
    df_prior_median = median(out$df_prior), df_total_median = median(df_total),
    max_abs_log2fc_diff_vs_frozen = max_dlfc, max_abs_p_diff_vs_frozen = max_dp,
    t_wald_sq_over_quasi_f_median = median(ratio),
    t_wald_sq_over_quasi_f_q25 = quantile(ratio, 0.25),
    t_wald_sq_over_quasi_f_q75 = quantile(ratio, 0.75)
  )

  path <- file.path(OUT_DIR, sprintf("%s_model_stats.tsv.gz", label))
  fwrite(out[order(gene)], path, sep = "\t")
  message(sprintf("  wrote %s (%d genes)", path, nrow(out)))
  invisible(out)
}

# ---------------------------------------------------------------------------
# GSE111151 -- edgeR, ~cell_line + resistance_status
# ---------------------------------------------------------------------------
cfg111 <- config$gse111151
run_edger_refit(
  label = "gse111151",
  counts_path = file.path(repo_root, cfg111$output$counts_tsv),
  metadata_path = file.path(repo_root, cfg111$output$metadata_tsv),
  gene_col = "gene_id",
  design_fn = function(md) {
    md$cell_line <- factor(md$cell_line)
    md$resistance_status <- factor(md$resistance_status, levels = c("parental", "resistant"))
    model.matrix(~cell_line + resistance_status, data = md)
  },
  coef_pattern = "^resistance_statusresistant$",
  frozen_path = file.path(repo_root, cfg111$output$de$genomewide_tsv),
  frozen_gene_col = "gene_id"
)

# ---------------------------------------------------------------------------
# GSE240112 -- edgeR, ~group (RT vs PT), tumour-cell pseudobulk
# ---------------------------------------------------------------------------
cfg240 <- config$gse240112
run_edger_refit(
  label = "gse240112",
  counts_path = file.path(repo_root, cfg240$output$tumor_cell$counts_tsv),
  metadata_path = file.path(repo_root, cfg240$output$tumor_cell$metadata_tsv),
  gene_col = "gene",
  design_fn = function(md) {
    md$group <- factor(md$group, levels = c("PT", "RT"))
    model.matrix(~group, data = md)
  },
  coef_pattern = "^groupRT$",
  frozen_path = file.path(repo_root, cfg240$output$de$genomewide_tsv),
  frozen_gene_col = "gene"
)

# ---------------------------------------------------------------------------
# GSE118713 -- limma on log2(TPM + 1), contrast TAMR - MCF7
#
# The frozen table already carries the exact moderated standard error
# (sqrt(s2.post) * stdev.unscaled). What it does not carry is df.prior, and
# therefore not the total df the moderated t was referred to. The refit
# recovers those and is verified to reproduce the frozen se and p-value.
# ---------------------------------------------------------------------------
message("=== gse118713 ===")
source(file.path(repo_root, "scripts", "analysis", "gse118713_limma_lib.R"))
cfg118 <- config$gse118713
cfg118b <- config$gse118713_phase2b
expr <- read.delim(gzfile(file.path(repo_root, cfg118b$filtering$filtered_gene_tpm_tsv)),
                   header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
meta <- read.delim(file.path(repo_root, cfg118$output$sample_metadata_tsv),
                   header = TRUE, sep = "\t", stringsAsFactors = FALSE)
log2_mat <- build_log2_matrix(expr, meta)
design118 <- build_design(meta)
fit118 <- fit_limma(log2_mat, design118)
CONTRAST <- "TAMR_vs_MCF7"

df_prior118 <- fit118$df.prior
df_res118 <- fit118$df.residual
if (length(df_res118) == 1L) df_res118 <- rep(df_res118, nrow(log2_mat))
if (length(df_prior118) == 1L) df_prior118 <- rep(df_prior118, nrow(log2_mat))
tt118 <- topTable(fit118, coef = CONTRAST, number = Inf, sort.by = "none")
se118 <- sqrt(fit118$s2.post) * fit118$stdev.unscaled[, CONTRAST]

out118 <- data.table(
  gene = expr$gene_id,
  gene_symbol = expr$gene_symbol,
  log2fc = tt118$logFC,
  se_log2fc_wald = se118,
  df_total = df_res118 + df_prior118,
  df_prior = df_prior118,
  df_residual = df_res118,
  dispersion = NA_real_,          # not a count model; s2.post plays this role
  var_post = fit118$s2.post,
  avg_log_cpm = tt118$AveExpr,    # mean log2(TPM+1), not CPM -- see report
  quasi_f = tt118$t^2,
  p_value_frozen_method = tt118$P.Value,
  fdr_frozen_method = tt118$adj.P.Val,
  engine = "limma_eBayes_log2TPM1",
  post_freeze_exploratory = TRUE
)
out118[, t_wald := log2fc / se_log2fc_wald]
out118[, p_value_wald := 2 * pt(-abs(t_wald), df = df_total)]

frozen118_path <- file.path(repo_root, cfg118b$unredaction$differential_expression_unredacted_tsv_gz)
frozen118 <- fread(cmd = paste("zcat", shQuote(frozen118_path)))
frozen118 <- frozen118[contrast == CONTRAST]
m118 <- merge(out118[, .(gene, log2fc, se_log2fc_wald, p_value_frozen_method)],
              frozen118[, .(gene = gene_id, log2fc_frozen = log2fc, se_frozen = se, p_frozen = p_value)],
              by = "gene")
log_rows("gse118713", "join to frozen table", nrow(out118), nrow(m118))
stopifnot(nrow(m118) == nrow(out118))
max_dlfc118 <- max(abs(m118$log2fc - m118$log2fc_frozen))
max_dse118 <- max(abs(m118$se_log2fc_wald - m118$se_frozen))
max_dp118 <- max(abs(m118$p_value_frozen_method - m118$p_frozen))
message(sprintf("  reproduction: max |dlog2fc| = %.3e, max |dse| = %.3e, max |dp| = %.3e",
                max_dlfc118, max_dse118, max_dp118))
if (max_dlfc118 > 1e-8 || max_dse118 > 1e-8 || max_dp118 > 1e-10) {
  stop("gse118713 refit does NOT reproduce the frozen table -- refusing to write")
}

# GSE118713 starts from an already-filtered TPM matrix, so its filter loss is
# recorded upstream rather than here; the counts are carried across so the
# validation table is comparable across all three datasets.
n_before_118 <- config$gse118713_phase2b$expected_n_genes
validation[["gse118713"]] <- data.table(
  dataset = "gse118713", engine = "limma_eBayes_log2TPM1", n_genes = nrow(out118),
  n_genes_before_filter = n_before_118,
  n_genes_removed_by_filter = n_before_118 - nrow(out118),
  n_genes_after_filter = nrow(out118),
  n_genes_joined_to_frozen = nrow(m118),
  n_samples = ncol(log2_mat), design = paste(colnames(design118), collapse = " + "),
  coef = CONTRAST, df_residual = df_res118[1],
  df_prior_median = median(df_prior118), df_total_median = median(out118$df_total),
  max_abs_log2fc_diff_vs_frozen = max_dlfc118, max_abs_p_diff_vs_frozen = max_dp118,
  t_wald_sq_over_quasi_f_median = 1, t_wald_sq_over_quasi_f_q25 = 1,
  t_wald_sq_over_quasi_f_q75 = 1
)
path118 <- file.path(OUT_DIR, "gse118713_model_stats.tsv.gz")
fwrite(out118[order(gene)], path118, sep = "\t")
message(sprintf("  wrote %s (%d genes)", path118, nrow(out118)))

val <- rbindlist(validation, use.names = TRUE)
val[, post_freeze_exploratory := TRUE]
fwrite(val, file.path(OUT_DIR, "refit_validation.tsv"), sep = "\t")
message("wrote refit_validation.tsv")
message("Done.")
