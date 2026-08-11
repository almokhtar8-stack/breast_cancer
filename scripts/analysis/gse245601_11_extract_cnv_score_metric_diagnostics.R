#!/usr/bin/env Rscript
# GSE245601 -- per-cell CNV signal-extent diagnostics (docs/CNV_METHOD_AUDIT.md
# Point 1 follow-up: is the whole-genome mean-squared CNV score an
# appropriate malignancy metric when abnormalities affect only part of the
# genome?). Read-only reload of the already-frozen run.final.infercnv_obj
# per sample -- InferCNV is NOT rerun; this only recomputes diagnostic
# summary statistics from its already-saved output. No malignancy label is
# changed, no threshold is changed, no candidate gene is read anywhere.
#
# For each of the 12 selected samples:
#   1. epi_cells = the frozen gse245601_malignant_cell_labels.tsv rows for
#      that sample (already verified elsewhere to equal InferCNV's own
#      observation-cell set exactly).
#   2. expr = run.final.infercnv_obj@expr.data[, epi_cells] -- the FULL
#      InferCNV gene set for this sample (unlike the earlier heatmap-
#      comparison extraction, genes are NOT intersected with CopyKAT here:
#      this analysis is about InferCNV's own signal on its own terms).
#   3. centered = expr - 1; abs_dev = abs(centered).
#   4. cnv_score_recomputed = colMeans(centered^2) -- cross-checked against
#      the frozen cnv_score column; the script stops if any cell's
#      recomputed score does not match the frozen one (a mismatch would
#      mean this script is reading the wrong object/cells, and nothing
#      downstream should be trusted).
#   5. Per cell: fraction of genes with abs_dev exceeding each of
#      config deviation_levels (0.05, 0.10, 0.15, pre-declared, not tuned
#      to any outcome); max/p95/p99 abs_dev across genes.
#   6. Per (cell, chromosome): mean abs_dev across that chromosome's genes
#      (from InferCNV's own @gene_order). Per cell: count of chromosomes
#      whose own mean abs_dev exceeds each deviation level.
#
# Requires the sc245601 micromamba env (infercnv, yaml):
#   micromamba run -n sc245601 Rscript scripts/analysis/gse245601_11_extract_cnv_score_metric_diagnostics.R
#
# Input:  data/processed/gse245601/infercnv/<sample>/run.final.infercnv_obj
#         results/tables/gse245601_malignant_cell_labels.tsv
# Output: results/tables/gse245601_infercnv_score_metric_diagnostics.tsv (committed)

suppressPackageStartupMessages({
  library(infercnv)
  library(yaml)
})

CNV_SCORE_TOLERANCE <- 1e-9
N_CHR_AUTOSOMES <- 22

args <- commandArgs(trailingOnly = FALSE)
repo_root <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", args, value = TRUE))), "..", ".."))
if (length(repo_root) == 0 || !dir.exists(file.path(repo_root, "config"))) repo_root <- getwd()

config <- yaml::read_yaml(file.path(repo_root, "config", "config.yaml"))
cfg <- config$gse245601_infercnv_score_metric_diagnostics

selected_samples <- cfg$selected_samples
infercnv_dir <- file.path(repo_root, cfg$inputs$infercnv_dir)
labels_path <- file.path(repo_root, cfg$inputs$infercnv_labels_tsv)
deviation_levels <- as.numeric(cfg$deviation_levels)
out_path <- file.path(repo_root, cfg$output$diagnostics_tsv)

dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)

labels <- read.delim(labels_path, stringsAsFactors = FALSE)

all_rows <- list()

for (sample_id in selected_samples) {
  message("=== ", sample_id, " ===")
  sample_labels <- labels[labels$sample_id == sample_id, ]
  epi_cells <- sample_labels$cell_id
  n_epi <- length(epi_cells)
  if (n_epi == 0) stop(sprintf("%s: 0 epithelial cells in the frozen label table", sample_id))
  if (anyDuplicated(epi_cells) > 0) stop(sprintf("%s: duplicate cell_id in the frozen label table", sample_id))

  obj <- readRDS(file.path(infercnv_dir, sample_id, "run.final.infercnv_obj"))
  missing_cells <- setdiff(epi_cells, colnames(obj@expr.data))
  if (length(missing_cells) > 0) {
    stop(sprintf("%s: %d epithelial cell(s) absent from run.final.infercnv_obj", sample_id, length(missing_cells)))
  }

  expr <- obj@expr.data[, epi_cells, drop = FALSE]
  centered <- expr - 1
  abs_dev <- abs(centered)
  n_genes <- nrow(expr)

  chr_raw <- as.character(obj@gene_order$chr)
  chr_num <- suppressWarnings(as.integer(sub("^chr", "", chr_raw)))
  if (any(is.na(chr_num))) stop(sprintf("%s: non-numeric chromosome value in gene_order", sample_id))
  stopifnot(length(chr_num) == n_genes)
  chrs_present <- sort(unique(chr_num))
  n_chr_total <- length(chrs_present)

  cnv_score_recomputed <- colMeans(centered^2)
  frozen_score <- setNames(sample_labels$cnv_score, sample_labels$cell_id)[epi_cells]
  score_diff <- abs(cnv_score_recomputed - frozen_score)
  if (any(score_diff > CNV_SCORE_TOLERANCE)) {
    stop(sprintf(
      "%s: recomputed cnv_score does not match the frozen cnv_score for %d cell(s) (max diff=%.3g) -- refusing to proceed on a mismatched premise",
      sample_id, sum(score_diff > CNV_SCORE_TOLERANCE), max(score_diff)
    ))
  }

  # per-(chromosome, cell) mean abs deviation via rowsum (fast, vectorized)
  chr_sums <- rowsum(abs_dev, group = chr_num) # rows = sorted chromosomes, cols = cells
  chr_gene_counts <- as.numeric(table(chr_num)[as.character(chrs_present)])
  chr_means <- chr_sums / chr_gene_counts # row-wise divide (length matches nrow)

  frac_cols <- list()
  n_chr_affected_cols <- list()
  for (lvl in deviation_levels) {
    frac_cols[[paste0("fraction_genes_dev_gt_", lvl)]] <- colMeans(abs_dev > lvl)
    n_chr_affected_cols[[paste0("n_chromosomes_affected_", lvl)]] <- colSums(chr_means > lvl)
  }

  max_abs_dev <- apply(abs_dev, 2, max)
  p95_abs_dev <- apply(abs_dev, 2, quantile, probs = 0.95, names = FALSE)
  p99_abs_dev <- apply(abs_dev, 2, quantile, probs = 0.99, names = FALSE)

  chr_mean_df <- as.data.frame(t(chr_means))
  colnames(chr_mean_df) <- paste0("chr", chrs_present, "_mean_abs_dev")
  # pad to a fixed chr1..chr22 column set (NA where a chromosome had 0 genes
  # passing filtering for this particular sample) so every sample's row can
  # be row-bound into one table without a column-set mismatch.
  for (k in seq_len(N_CHR_AUTOSOMES)) {
    nm <- paste0("chr", k, "_mean_abs_dev")
    if (!nm %in% names(chr_mean_df)) chr_mean_df[[nm]] <- NA_real_
  }
  chr_mean_df <- chr_mean_df[, paste0("chr", seq_len(N_CHR_AUTOSOMES), "_mean_abs_dev")]

  sample_df <- data.frame(
    cell_id = epi_cells,
    sample_id = sample_id,
    patient = sample_labels$patient,
    condition = sample_labels$condition,
    threshold_group = sample_labels$threshold_group,
    primary_malignancy_label = sample_labels$primary_malignancy_label,
    cnv_correlation_to_seed = sample_labels$cnv_correlation_to_seed,
    cnv_score = frozen_score,
    cnv_score_recomputed = cnv_score_recomputed,
    n_genes_total = n_genes,
    n_chromosomes_total = n_chr_total,
    max_abs_deviation = max_abs_dev,
    p95_abs_deviation = p95_abs_dev,
    p99_abs_deviation = p99_abs_dev,
    stringsAsFactors = FALSE
  )
  for (nm in names(frac_cols)) sample_df[[nm]] <- frac_cols[[nm]]
  for (nm in names(n_chr_affected_cols)) sample_df[[nm]] <- n_chr_affected_cols[[nm]]
  sample_df <- cbind(sample_df, chr_mean_df)

  message(sprintf(
    "  n_cells=%d, n_genes=%d, n_chromosomes_present=%d, cnv_score reconstruction verified (max diff=%.3g)",
    n_epi, n_genes, n_chr_total, max(score_diff)
  ))

  all_rows[[sample_id]] <- sample_df
  rm(obj, expr, centered, abs_dev, chr_sums, chr_means)
  gc(verbose = FALSE)
}

out <- do.call(rbind, all_rows)
rownames(out) <- NULL
write.table(out, out_path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
message(sprintf("Wrote %s (%d rows, %d cols)", out_path, nrow(out), ncol(out)))
