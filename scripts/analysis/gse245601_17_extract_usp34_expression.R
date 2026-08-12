#!/usr/bin/env Rscript
# GSE245601 -- USP34 expression extraction for visualization only
# (read-only follow-up). Does NOT rerun QC/normalization/clustering/
# InferCNV/CopyKAT/pseudobulk, does NOT change any label or eligibility
# rule. Loads the already-frozen Seurat object only to read USP34's
# already-computed log-normalized expression ("data" layer) for the
# already-frozen epithelial cell set -- clustering, annotation, and UMAP
# coordinates are reused as-is (from the already-frozen
# gse245601_pseudobulk_..._cell_level_summary.tsv table joined in Python,
# not recomputed here).
#
# Input:  data/processed/gse245601/seurat_clustered/annotated.rds
#         results/tables/gse245601_pseudobulk/malignant_vs_nonmalignant_cell_level_summary.tsv
#           (read-only, only to get the exact frozen epithelial cell_id set)
# Output: results/tables/gse245601_usp34_expression/usp34_per_cell_expression.tsv (committed)

suppressPackageStartupMessages({
  library(Seurat)
  library(yaml)
})

args <- commandArgs(trailingOnly = FALSE)
repo_root <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", args, value = TRUE))), "..", ".."))
if (length(repo_root) == 0 || !dir.exists(file.path(repo_root, "config"))) repo_root <- getwd()

config <- yaml::read_yaml(file.path(repo_root, "config", "config.yaml"))
cfg <- config$gse245601_usp34_expression
gene <- cfg$gene

out_path <- file.path(repo_root, cfg$output$expression_tsv)
dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)

epi_cells_table <- read.delim(file.path(repo_root, cfg$inputs$cell_level_summary_tsv), stringsAsFactors = FALSE)
epi_cell_ids <- epi_cells_table$cell_id

merged <- readRDS(file.path(repo_root, cfg$inputs$seurat_object))

n_found <- sum(rownames(merged) == gene)
message(sprintf("%s found in feature space: %d occurrence(s)", gene, n_found))
if (n_found != 1) stop(sprintf("%s is not present exactly once in the Seurat object's feature space -- refusing to proceed", gene))

missing_cells <- setdiff(epi_cell_ids, colnames(merged))
if (length(missing_cells) > 0) stop(sprintf("%d epithelial cell(s) from the frozen table are absent from the Seurat object", length(missing_cells)))

data_layer <- GetAssayData(merged, layer = "data")
expr <- data_layer[gene, epi_cell_ids]

out <- data.frame(cell_id = epi_cell_ids, usp34_log_norm_expression = as.numeric(expr), stringsAsFactors = FALSE)
write.table(out, out_path, sep = "\t", quote = FALSE, row.names = FALSE)
message(sprintf("wrote %s (%d cells, %d expressing >0)", out_path, nrow(out), sum(out$usp34_log_norm_expression > 0)))
