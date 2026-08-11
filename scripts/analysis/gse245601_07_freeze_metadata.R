#!/usr/bin/env Rscript
# GSE245601 -- Step 10: freeze immutable per-cell metadata and compute
# pseudobulk pair-eligibility flags. Does NOT run pseudobulk itself.
#
# Combines: cluster/broad_cell_type/epithelial flag (Step 6-7), CNV
# metrics + primary (InferCNV) malignancy label (Step 8), sensitivity
# (CopyKAT) malignancy label (Step 9) into one per-cell table covering
# every QC-passed cell in the 20 primary tumor samples (not just
# epithelial cells, so the full compositional context is preserved).
#
# Eligibility rule (frozen in docs/gse245601_PREANALYSIS.md section 11,
# BEFORE any malignant-cell count was known): a tumor pair (Control +
# Tamoxifen) is provisionally eligible for a future pseudobulk phase only
# if BOTH arms contain >=50 primary (InferCNV) malignant cells.
#
# No candidate gene is read, plotted, or used anywhere in this script. No
# pseudobulk aggregation, no differential expression, no candidate ranking.
#
# Output: results/tables/gse245601_cell_metadata_frozen.tsv.gz (committed)
#         results/tables/gse245601_malignant_counts_per_sample.tsv (committed)
#         results/tables/gse245601_pair_eligibility.tsv (committed)

suppressPackageStartupMessages({
  library(Seurat)
})

MALIGNANT_CELL_MIN_COUNT <- 50

args <- commandArgs(trailingOnly = FALSE)
repo_root <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", args, value = TRUE))), "..", ".."))
if (length(repo_root) == 0 || !dir.exists(file.path(repo_root, "config"))) repo_root <- getwd()

merged <- readRDS(file.path(repo_root, "data", "processed", "gse245601", "seurat_clustered", "annotated.rds"))
malignant_labels <- read.delim(file.path(repo_root, "results", "tables", "gse245601_malignant_cell_labels.tsv"), stringsAsFactors = FALSE)
sensitivity_path <- file.path(repo_root, "results", "tables", "gse245601_copykat_sensitivity_labels.tsv")
sensitivity_labels <- if (file.exists(sensitivity_path)) read.delim(sensitivity_path, stringsAsFactors = FALSE) else NULL

meta <- merged@meta.data
meta$cell_id <- rownames(meta)

full_meta <- data.frame(
  cell_id = meta$cell_id, GSM = meta$GSM, patient = meta$patient, condition = meta$condition,
  sample_id = meta$sample_id, nCount_RNA = meta$nCount_RNA, nFeature_RNA = meta$nFeature_RNA,
  percent.mt = meta$percent.mt, seurat_cluster = meta$seurat_clusters, broad_cell_type = meta$broad_cell_type,
  is_epithelial = meta$broad_cell_type == "epithelial",
  stringsAsFactors = FALSE
)

full_meta <- merge(full_meta, malignant_labels[, c("cell_id", "cnv_score", "cnv_correlation_to_seed", "threshold_group", "primary_malignancy_label")], by = "cell_id", all.x = TRUE)
if (!is.null(sensitivity_labels)) {
  full_meta <- merge(full_meta, sensitivity_labels[, c("cell_id", "sensitivity_malignancy_label")], by = "cell_id", all.x = TRUE)
} else {
  full_meta$sensitivity_malignancy_label <- NA
}

out1 <- file.path(repo_root, "results", "tables", "gse245601_cell_metadata_frozen.tsv.gz")
gz <- gzfile(out1, "w")
write.table(full_meta, gz, sep = "\t", quote = FALSE, row.names = FALSE)
close(gz)
message(sprintf("Wrote %s (%d cells, %d columns)", out1, nrow(full_meta), ncol(full_meta)))

# --- Per-sample malignant counts + eligibility -----------------------------
samples <- sort(unique(full_meta$sample_id))
count_rows <- list()
for (s in samples) {
  sub <- full_meta[full_meta$sample_id == s, ]
  count_rows[[s]] <- data.frame(
    sample_id = s,
    patient = unique(sub$patient),
    condition = unique(sub$condition),
    n_qc_passed_cells = nrow(sub),
    n_epithelial_cells = sum(sub$is_epithelial, na.rm = TRUE),
    n_primary_malignant = sum(sub$primary_malignancy_label == "malignant", na.rm = TRUE),
    n_primary_nonmalignant_epithelial = sum(sub$primary_malignancy_label == "non-malignant epithelial", na.rm = TRUE),
    malignant_fraction_of_epithelial = {
      denom <- sum(sub$is_epithelial, na.rm = TRUE)
      if (denom > 0) sum(sub$primary_malignancy_label == "malignant", na.rm = TRUE) / denom else NA
    },
    stringsAsFactors = FALSE
  )
}
counts <- do.call(rbind, count_rows)
rownames(counts) <- NULL
out2 <- file.path(repo_root, "results", "tables", "gse245601_malignant_counts_per_sample.tsv")
write.table(counts, out2, sep = "\t", quote = FALSE, row.names = FALSE)
message(sprintf("Wrote %s (%d samples)", out2, nrow(counts)))

# --- Pair eligibility (frozen rule; pseudobulk itself NOT run here) --------
patients <- sort(unique(counts$patient))
elig_rows <- list()
for (p in patients) {
  psub <- counts[counts$patient == p, ]
  ctrl <- psub$n_primary_malignant[psub$condition == "Control"]
  tam <- psub$n_primary_malignant[psub$condition == "Tamoxifen"]
  ctrl <- if (length(ctrl) == 1) ctrl else NA
  tam <- if (length(tam) == 1) tam else NA
  eligible <- !is.na(ctrl) && !is.na(tam) && ctrl >= MALIGNANT_CELL_MIN_COUNT && tam >= MALIGNANT_CELL_MIN_COUNT
  elig_rows[[p]] <- data.frame(
    patient = p, n_malignant_control = ctrl, n_malignant_tamoxifen = tam,
    min_required_per_arm = MALIGNANT_CELL_MIN_COUNT, eligible_for_pseudobulk = eligible,
    stringsAsFactors = FALSE
  )
}
eligibility <- do.call(rbind, elig_rows)
rownames(eligibility) <- NULL
out3 <- file.path(repo_root, "results", "tables", "gse245601_pair_eligibility.tsv")
write.table(eligibility, out3, sep = "\t", quote = FALSE, row.names = FALSE)
message(sprintf("Wrote %s: %d/%d tumor pairs eligible", out3, sum(eligibility$eligible_for_pseudobulk), nrow(eligibility)))
