#!/usr/bin/env Rscript
# GSE245601 candidate deep-dive -- per-cell extraction only (read-only
# follow-up, same discipline as gse245601_17_extract_usp34_expression.R).
# Does NOT rerun QC/normalization/clustering/InferCNV/CopyKAT/pseudobulk,
# does NOT change any frozen label, does NOT recompute UMAP/PCA. Loads
# the already-frozen Seurat object only to read the four candidate genes'
# already-computed raw counts ("counts" layer) and log-normalized
# expression ("data" layer) for the already-frozen epithelial cell set,
# plus the existing seurat_clusters assignment (the only clustering that
# exists in this project -- computed on the full multi-lineage object,
# not a separate epithelial-only subclustering; epithelial cells occupy a
# subset of these cluster IDs). Malignancy status, UMAP coordinates,
# library size, and patient/condition metadata are reused as-is from the
# already-frozen malignant_vs_nonmalignant_cell_level_summary.tsv --
# never recomputed here.
#
# Input:  data/processed/gse245601/seurat_clustered/annotated.rds
#         results/tables/gse245601_pseudobulk/malignant_vs_nonmalignant_cell_level_summary.tsv
# Output: results/tables/gse245601_candidate_deepdive/candidate_per_cell_expression.tsv (committed)

suppressPackageStartupMessages({
  library(Seurat)
  library(yaml)
})

args <- commandArgs(trailingOnly = FALSE)
repo_root <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", args, value = TRUE))), "..", ".."))
if (length(repo_root) == 0 || !dir.exists(file.path(repo_root, "config"))) repo_root <- getwd()

config <- yaml::read_yaml(file.path(repo_root, "config", "config.yaml"))
cfg <- config$gse245601_candidate_deepdive
genes <- unlist(cfg$genes)

out_path <- file.path(repo_root, cfg$output$per_cell_tsv)
dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)

cell_table <- read.delim(file.path(repo_root, cfg$inputs$cell_level_summary_tsv), stringsAsFactors = FALSE)
epi_cell_ids <- cell_table$cell_id
message(sprintf("frozen epithelial cell table: %d cells", length(epi_cell_ids)))

merged <- readRDS(file.path(repo_root, cfg$inputs$seurat_object))

missing_genes <- setdiff(genes, rownames(merged))
if (length(missing_genes) > 0) stop(sprintf("gene(s) not found in the Seurat object feature space: %s", paste(missing_genes, collapse = ", ")))
for (g in genes) {
  n_found <- sum(rownames(merged) == g)
  message(sprintf("%s found in feature space: %d occurrence(s)", g, n_found))
  if (n_found != 1) stop(sprintf("%s is not present exactly once -- refusing to proceed", g))
}

missing_cells <- setdiff(epi_cell_ids, colnames(merged))
if (length(missing_cells) > 0) stop(sprintf("%d epithelial cell(s) from the frozen table are absent from the Seurat object", length(missing_cells)))

epi <- subset(merged, cells = epi_cell_ids)
epi <- epi[, epi_cell_ids]  # enforce exact frozen order

counts_layer <- GetAssayData(epi, layer = "counts")
data_layer <- GetAssayData(epi, layer = "data")

out <- data.frame(cell_id = epi_cell_ids, seurat_clusters = as.character(epi@meta.data[epi_cell_ids, "seurat_clusters"]), stringsAsFactors = FALSE)
for (g in genes) {
  out[[paste0(g, "_raw_count")]] <- as.numeric(counts_layer[g, epi_cell_ids])
  out[[paste0(g, "_log_norm")]] <- as.numeric(data_layer[g, epi_cell_ids])
}

n_mismatch <- sum(out$cell_id != epi_cell_ids)
if (n_mismatch > 0) stop("cell_id order mismatch after extraction -- refusing to write a misaligned table")

write.table(out, out_path, sep = "\t", quote = FALSE, row.names = FALSE)
message(sprintf("wrote %s (%d cells x %d genes)", out_path, nrow(out), length(genes)))
