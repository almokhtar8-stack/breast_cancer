#!/usr/bin/env Rscript
# GSE245601 -- malignant vs non-malignant epithelial context check
# (gse245601_PREANALYSIS.md section 13, Phase 8). Validation/context only,
# NOT the treatment-response analysis: this script never touches
# Control-vs-Tamoxifen, and candidate-gene results are never used to
# modify the frozen InferCNV labels here or anywhere else.
#
# Two outputs:
#   1. A per-epithelial-cell descriptive table (CNV score already frozen
#      in gse245601_malignant_cell_labels.tsv; existing UMAP/PCA
#      coordinates; existing epithelial-program module score
#      (score_epithelial1, computed during the original candidate-blind
#      broad cell-type annotation); existing QC metrics
#      (nFeature_RNA, percent.mt); and cell-cycle/proliferation scores via
#      Seurat's standard, built-in CellCycleScoring() using the
#      field-standard cc.genes.updated.2019 gene sets -- not a candidate
#      gene, not project-specific). Purely descriptive; no p-value is
#      computed per cell (that would be pseudoreplication).
#   2. A malignant-vs-non-malignant-epithelial pseudobulk (raw counts,
#      pooled across BOTH treatment arms per patient, patient x
#      malignancy-status as the unit), for patients with >=50 pooled
#      malignant cells (config gse245601_pseudobulk.malignant_vs_nonmalignant_min_pooled_cells
#      -- the SAME threshold value as the frozen Track B per-arm rule,
#      reapplied to a different, separate, non-preregistered grouping;
#      this is not the frozen Track B eligibility rule). Used downstream
#      (Python) only for a secondary, exploratory 13-candidate + PAICS
#      context table, kept explicitly separate from the Control-vs-
#      Tamoxifen question.
#
# Input:  data/processed/gse245601/seurat_clustered/annotated.rds
#         results/tables/gse245601_malignant_cell_labels.tsv
# Output: results/tables/gse245601_pseudobulk/malignant_vs_nonmalignant_cell_level_summary.tsv (committed)
#         results/tables/gse245601_pseudobulk/malignant_vs_nonmalignant_counts.tsv.gz (committed)
#         results/tables/gse245601_pseudobulk/malignant_vs_nonmalignant_metadata.tsv (committed)

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(yaml)
  library(data.table)
})

args <- commandArgs(trailingOnly = FALSE)
repo_root <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", args, value = TRUE))), "..", ".."))
if (length(repo_root) == 0 || !dir.exists(file.path(repo_root, "config"))) repo_root <- getwd()

config <- yaml::read_yaml(file.path(repo_root, "config", "config.yaml"))
cfg <- config$gse245601_pseudobulk
min_pooled <- as.integer(cfg$malignant_vs_nonmalignant_min_pooled_cells)

out_cell_level <- file.path(repo_root, cfg$output$malignant_vs_nonmalignant$cell_level_summary_tsv)
out_counts <- file.path(repo_root, cfg$output$malignant_vs_nonmalignant$counts_tsv)
out_meta <- file.path(repo_root, cfg$output$malignant_vs_nonmalignant$metadata_tsv)
for (p in c(out_cell_level, out_counts, out_meta)) dir.create(dirname(p), recursive = TRUE, showWarnings = FALSE)

merged <- readRDS(file.path(repo_root, cfg$inputs$seurat_object))
malignant_labels <- read.delim(file.path(repo_root, cfg$inputs$infercnv_labels_tsv), stringsAsFactors = FALSE)

epi <- subset(merged, subset = broad_cell_type == "epithelial")
message(sprintf("epithelial cells: %d", ncol(epi)))

label_map <- setNames(malignant_labels$primary_malignancy_label, malignant_labels$cell_id)
epi_labels <- label_map[colnames(epi)]
if (any(is.na(epi_labels))) stop("some epithelial cells have no frozen malignancy label -- refusing to proceed")
epi$malignancy_status <- epi_labels

# --- 1. cell-level descriptive table ---
epi <- CellCycleScoring(
  epi,
  s.features = cc.genes.updated.2019$s.genes,
  g2m.features = cc.genes.updated.2019$g2m.genes,
  set.ident = FALSE
)

umap_coords <- Embeddings(epi, reduction = "umap")
pca_coords <- Embeddings(epi, reduction = "pca")[, 1:2]
cnv_score_map <- setNames(malignant_labels$cnv_score, malignant_labels$cell_id)

cell_level <- data.frame(
  cell_id = colnames(epi),
  sample_id = epi$sample_id,
  patient = epi$patient,
  condition = epi$condition,
  malignancy_status = epi$malignancy_status,
  cnv_score = cnv_score_map[colnames(epi)],
  score_epithelial_program = epi$score_epithelial1,
  nFeature_RNA = epi$nFeature_RNA,
  nCount_RNA = epi$nCount_RNA,
  percent_mt = epi$percent.mt,
  S_Score = epi$S.Score,
  G2M_Score = epi$G2M.Score,
  cell_cycle_phase = epi$Phase,
  umap_1 = umap_coords[, 1],
  umap_2 = umap_coords[, 2],
  pca_1 = pca_coords[, 1],
  pca_2 = pca_coords[, 2],
  stringsAsFactors = FALSE
)
write.table(cell_level, out_cell_level, sep = "\t", quote = FALSE, row.names = FALSE)
message(sprintf("wrote %s (%d cells)", out_cell_level, nrow(cell_level)))

# --- 2. malignant vs non-malignant pseudobulk, pooled across both arms, per patient ---
pooled_malignant_n <- table(epi$patient[epi$malignancy_status == "malignant"])
eligible_patients <- names(pooled_malignant_n)[pooled_malignant_n >= min_pooled]
message(sprintf("patients with >=%d pooled malignant cells: %s", min_pooled, paste(sort(eligible_patients), collapse = ", ")))

counts <- GetAssayData(merged, layer = "counts")
stopifnot(!anyDuplicated(rownames(counts))) # gene symbols must be unique before any pseudobulk row-sum aggregation
cell_groups <- list()
for (p in sort(eligible_patients)) {
  mal_cells <- colnames(epi)[epi$patient == p & epi$malignancy_status == "malignant"]
  nonmal_cells <- colnames(epi)[epi$patient == p & epi$malignancy_status == "non-malignant epithelial"]
  cell_groups[[sprintf("%s_malignant", p)]] <- mal_cells
  cell_groups[[sprintf("%s_nonmalignant", p)]] <- nonmal_cells
}
stopifnot(length(cell_groups) == length(eligible_patients) * 2)

pb_mat <- sapply(cell_groups, function(cell_ids) Matrix::rowSums(counts[, cell_ids, drop = FALSE]))
rownames(pb_mat) <- rownames(counts)
pb_meta <- data.frame(
  sample_id = names(cell_groups),
  patient = sub("_(malignant|nonmalignant)$", "", names(cell_groups)),
  malignancy_status = sub("^.*_", "", names(cell_groups)),
  n_contributing_cells = sapply(cell_groups, length),
  total_library_size = colSums(pb_mat),
  n_detected_genes = colSums(pb_mat > 0),
  stringsAsFactors = FALSE
)

out_dt <- data.table(gene = rownames(pb_mat), as.data.frame(pb_mat))
fwrite(out_dt, out_counts, sep = "\t")
write.table(pb_meta, out_meta, sep = "\t", quote = FALSE, row.names = FALSE)
message(sprintf("wrote %s (%d genes x %d samples) and %s", out_counts, nrow(pb_mat), ncol(pb_mat), out_meta))
message("Done. No candidate-gene expression was used to derive or modify any malignancy label in this script.")
