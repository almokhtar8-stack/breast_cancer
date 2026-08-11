#!/usr/bin/env Rscript
# GSE245601 -- Step 8: per-sample InferCNV + malignant-cell classification.
#
# Implements the frozen procedure in docs/gse245601_PREANALYSIS.md section 9
# ("ng_2021_and_thresholding"), reconstructed from the paper's Methods and
# the author code (run_infercnv.R, identify_reference_clusters.R) at pinned
# commit ceabf3f331c88f464e6a57b0ad1f9c500bedde85:
#
#   1. Reference (diploid) cells per sample = immune (t_nk + b_plasma +
#      myeloid) + endothelial cells, from that sample only.
#   2. Observation cells per sample = epithelial cells, from that sample only.
#   3. infercnv::run(), "vignettes" (non-HMM) preset: cutoff=0.1,
#      min_cells_per_gene=3, cluster_by_groups=TRUE, scale_data=FALSE,
#      HMM=FALSE, analysis_mode="samples", denoise=TRUE, no_prelim_plot=TRUE.
#   4. Per epithelial cell, CNV score = mean((expr - 1)^2) across genes
#      (infercnv's default output is centered at 1.0 = no relative change).
#   5. Within each (sample, Seurat cluster) group of epithelial cells, take
#      the top 5% (minimum 2) of cells by CNV score as a malignant "seed";
#      correlate (Kendall tau) every cell's CNV profile to the seed's mean
#      profile.
#   6. Adaptive per-group thresholds: th_value = mean(CNV score) - 2*SD,
#      clamped to [0.01, 0.05]; th_corr = mean(correlation) - 1.5*SD,
#      clamped to [0.2, 0.4]. Malignant if CNV score > th_value AND
#      correlation > th_corr; else non-malignant epithelial.
#
# Treatment condition (Control/Tamoxifen) is NEVER an input to any step of
# this procedure -- reference-cell selection, infercnv itself, CNV-score
# computation, and thresholds are all computed identically regardless of
# condition. No candidate gene (13 candidates + PAICS) is read, plotted, or
# used anywhere in this script.
#
# A sample with fewer than 20 combined reference cells after QC is flagged
# `insufficient_reference` and excluded from CNV/malignancy calling for
# that sample only (frozen fallback rule, PREANALYSIS section 9).
# A (sample, cluster) group with fewer than 10 epithelial cells falls back
# to whole-sample (all epithelial clusters pooled) thresholding for that
# sample -- frozen here, amendment 3, before any malignancy label was
# computed.
#
# Input:  data/processed/gse245601/seurat_clustered/annotated.rds
#         data/processed/gse245601/infercnv/gene_order.tsv
# Output: data/processed/gse245601/infercnv/<sample_id>/ (gitignored, per-sample infercnv working dirs)
#         results/tables/gse245601_malignant_cell_labels.tsv (committed; per-cell)
#         results/tables/gse245601_malignant_summary_per_sample.tsv (committed)

suppressPackageStartupMessages({
  library(Seurat)
  library(infercnv)
  library(Matrix)
  library(pcaPP)  # cor.fk(): exact Kendall tau via a fast O(n log n) algorithm,
                  # numerically identical to stats::cor(method="kendall") (verified
                  # to machine precision) -- required because base R's cor() uses
                  # an O(n^2) algorithm for Kendall that is intractable at the
                  # ~8,000-gene vector length used here. This is a performance-only
                  # substitution, not a change to the statistic being computed.
})

set.seed(42)
MIN_REFERENCE_CELLS <- 20
MIN_CLUSTER_CELLS_FOR_OWN_THRESHOLD <- 10
SEED_TOP_FRACTION <- 0.05

args <- commandArgs(trailingOnly = FALSE)
repo_root <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", args, value = TRUE))), "..", ".."))
if (length(repo_root) == 0 || !dir.exists(file.path(repo_root, "config"))) repo_root <- getwd()

merged <- readRDS(file.path(repo_root, "data", "processed", "gse245601", "seurat_clustered", "annotated.rds"))
gene_order_path <- file.path(repo_root, "data", "processed", "gse245601", "infercnv", "gene_order.tsv")
gene_order <- read.delim(gene_order_path, header = FALSE, col.names = c("gene", "chr", "start", "end"), stringsAsFactors = FALSE)

immune_types <- c("t_nk", "b_plasma", "myeloid")
endothelial_type <- "endothelial"
epithelial_type <- "epithelial"

samples <- sort(unique(merged$sample_id))
message("Samples to process: ", length(samples))

dir_infercnv_root <- file.path(repo_root, "data", "processed", "gse245601", "infercnv")

cell_label_rows <- list()
sample_summary_rows <- list()

kendall_cor_to_seed <- function(mat_centered, seed_mean_profile) {
  apply(mat_centered, 2, function(cell_vec) {
    suppressWarnings(pcaPP::cor.fk(cell_vec, seed_mean_profile))
  })
}

for (sample_id in samples) {
  message("=== ", sample_id, " ===")
  sample_cells <- colnames(merged)[merged$sample_id == sample_id]
  sample_obj <- subset(merged, cells = sample_cells)

  ref_cells <- colnames(sample_obj)[sample_obj$broad_cell_type %in% c(immune_types, endothelial_type)]
  obs_cells <- colnames(sample_obj)[sample_obj$broad_cell_type == epithelial_type]

  n_ref <- length(ref_cells)
  n_obs <- length(obs_cells)
  message(sprintf("  reference cells=%d, epithelial (observation) cells=%d", n_ref, n_obs))

  if (n_ref < MIN_REFERENCE_CELLS || n_obs < 2) {
    sample_summary_rows[[sample_id]] <- data.frame(
      sample_id = sample_id, n_reference_cells = n_ref, n_epithelial_cells = n_obs,
      status = "insufficient_reference", n_malignant = NA, n_nonmalignant_epithelial = NA,
      stringsAsFactors = FALSE
    )
    next
  }

  cnv_input_cells <- c(ref_cells, obs_cells)
  raw_counts <- GetAssayData(sample_obj, layer = "counts")[, cnv_input_cells]
  common_genes <- intersect(rownames(raw_counts), gene_order$gene)
  raw_counts <- raw_counts[common_genes, ]
  ordered_genes <- gene_order[match(common_genes, gene_order$gene), ]

  annotations_df <- data.frame(
    cell = cnv_input_cells,
    group = ifelse(cnv_input_cells %in% ref_cells, "reference", "observation"),
    row.names = cnv_input_cells
  )

  sample_dir <- file.path(dir_infercnv_root, sample_id)
  dir.create(sample_dir, recursive = TRUE, showWarnings = FALSE)

  gene_order_tmp <- file.path(sample_dir, "gene_order_used.tsv")
  write.table(ordered_genes, gene_order_tmp, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  annotations_tmp <- file.path(sample_dir, "annotations.tsv")
  write.table(annotations_df[, "group", drop = FALSE], annotations_tmp, sep = "\t", quote = FALSE, col.names = FALSE)

  infercnv_obj <- CreateInfercnvObject(
    raw_counts_matrix = raw_counts,
    annotations_file = annotations_tmp,
    gene_order_file = gene_order_tmp,
    ref_group_names = "reference"
  )

  run_result <- tryCatch(
    {
      infercnv::run(
        infercnv_obj, cutoff = 0.1, min_cells_per_gene = 3, out_dir = sample_dir,
        cluster_by_groups = TRUE, scale_data = FALSE, HMM = FALSE, analysis_mode = "samples",
        denoise = TRUE, no_prelim_plot = TRUE, png_res = 60, num_threads = 8,
        plot_steps = FALSE, resume_mode = FALSE, save_rds = FALSE
      )
    },
    error = function(e) {
      message("  InferCNV run FAILED: ", conditionMessage(e))
      NULL
    }
  )

  if (is.null(run_result)) {
    sample_summary_rows[[sample_id]] <- data.frame(
      sample_id = sample_id, n_reference_cells = n_ref, n_epithelial_cells = n_obs,
      status = "infercnv_run_failed", n_malignant = NA, n_nonmalignant_epithelial = NA,
      stringsAsFactors = FALSE
    )
    next
  }

  expr_mat <- run_result@expr.data[, obs_cells, drop = FALSE]
  centered <- expr_mat - 1
  cnv_score <- colMeans(centered^2)

  epi_clusters <- as.character(sample_obj$seurat_clusters[obs_cells])
  names(epi_clusters) <- obs_cells

  cell_status <- setNames(rep(NA_character_, length(obs_cells)), obs_cells)
  cell_cnv_corr <- setNames(rep(NA_real_, length(obs_cells)), obs_cells)
  cell_threshold_group <- setNames(rep(NA_character_, length(obs_cells)), obs_cells)

  cluster_tab <- table(epi_clusters)
  small_clusters <- names(cluster_tab)[cluster_tab < MIN_CLUSTER_CELLS_FOR_OWN_THRESHOLD]
  groups_to_process <- if (length(small_clusters) > 0) {
    grp <- epi_clusters
    grp[epi_clusters %in% small_clusters] <- "whole_sample_pooled"
    grp
  } else {
    epi_clusters
  }

  for (grp in unique(groups_to_process)) {
    grp_cells <- names(groups_to_process)[groups_to_process == grp]
    grp_scores <- cnv_score[grp_cells]
    n_seed <- max(round(length(grp_cells) * SEED_TOP_FRACTION), 2)
    n_seed <- min(n_seed, length(grp_cells))
    seed_cells <- names(sort(grp_scores, decreasing = TRUE))[1:n_seed]
    seed_mean_profile <- rowMeans(centered[, seed_cells, drop = FALSE])

    grp_corr <- kendall_cor_to_seed(centered[, grp_cells, drop = FALSE], seed_mean_profile)

    th_value <- mean(grp_scores) - 2 * sd(grp_scores)
    th_value <- max(min(th_value, 0.05), 0.01)
    th_corr <- mean(grp_corr, na.rm = TRUE) - 1.5 * sd(grp_corr, na.rm = TRUE)
    th_corr <- max(min(th_corr, 0.4), 0.2)

    is_malignant <- (grp_scores > th_value) & (grp_corr > th_corr)
    is_malignant[is.na(is_malignant)] <- FALSE

    cell_status[grp_cells] <- ifelse(is_malignant, "malignant", "non-malignant epithelial")
    cell_cnv_corr[grp_cells] <- grp_corr
    cell_threshold_group[grp_cells] <- grp
  }

  cell_label_rows[[sample_id]] <- data.frame(
    cell_id = obs_cells,
    sample_id = sample_id,
    patient = sample_obj$patient[obs_cells],
    condition = sample_obj$condition[obs_cells],
    seurat_cluster = epi_clusters[obs_cells],
    threshold_group = cell_threshold_group[obs_cells],
    cnv_score = cnv_score[obs_cells],
    cnv_correlation_to_seed = cell_cnv_corr[obs_cells],
    primary_malignancy_label = cell_status[obs_cells],
    stringsAsFactors = FALSE
  )

  n_malignant <- sum(cell_status == "malignant", na.rm = TRUE)
  n_nonmalignant <- sum(cell_status == "non-malignant epithelial", na.rm = TRUE)
  sample_summary_rows[[sample_id]] <- data.frame(
    sample_id = sample_id, n_reference_cells = n_ref, n_epithelial_cells = n_obs,
    status = "ok", n_malignant = n_malignant, n_nonmalignant_epithelial = n_nonmalignant,
    stringsAsFactors = FALSE
  )
  message(sprintf("  malignant=%d, non-malignant epithelial=%d", n_malignant, n_nonmalignant))
}

cell_labels <- do.call(rbind, cell_label_rows)
rownames(cell_labels) <- NULL
out1 <- file.path(repo_root, "results", "tables", "gse245601_malignant_cell_labels.tsv")
write.table(cell_labels, out1, sep = "\t", quote = FALSE, row.names = FALSE)
message(sprintf("Wrote %s (%d cells)", out1, nrow(cell_labels)))

sample_summary <- do.call(rbind, sample_summary_rows)
rownames(sample_summary) <- NULL
out2 <- file.path(repo_root, "results", "tables", "gse245601_malignant_summary_per_sample.tsv")
write.table(sample_summary, out2, sep = "\t", quote = FALSE, row.names = FALSE)
message(sprintf("Wrote %s (%d samples)", out2, nrow(sample_summary)))
