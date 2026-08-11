#!/usr/bin/env Rscript
# GSE245601 -- Step 5: normalization, PCA, clustering, UMAP.
#
# Frozen rules (docs/gse245601_PREANALYSIS.md section 6): LogNormalize
# (scale.factor=10000), 2000 variable features, ScaleData regressing
# percent.mt, RunPCA(npcs=50), neighbor graph + Louvain clustering on the
# first 30 PCs at resolution 0.8, UMAP (n_neighbors=30, min_dist=0.3,
# metric="cosine"). Raw UMI counts are preserved in the "RNA" assay's
# counts layer, untouched, for future pseudobulk use; normalized/scaled
# values are used only for clustering/visualization here.
#
# This script performs no candidate-gene-driven parameter tuning and does
# not inspect any of the 13 candidate genes or PAICS.
#
# Input:  data/processed/gse245601/seurat_qc/<patient>_<condition>.rds (20 files)
# Output: data/processed/gse245601/seurat_clustered/merged.rds (gitignored)
#         results/tables/gse245601_cluster_summary.tsv (committed)

suppressPackageStartupMessages({
  library(Seurat)
  library(yaml)
})

set.seed(42)

args <- commandArgs(trailingOnly = FALSE)
repo_root <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", args, value = TRUE))), "..", ".."))
if (length(repo_root) == 0 || !dir.exists(file.path(repo_root, "config"))) repo_root <- getwd()

config <- yaml::read_yaml(file.path(repo_root, "config", "config.yaml"))
cfg <- config$gse245601

dir_seurat_qc <- file.path(repo_root, "data", "processed", "gse245601", "seurat_qc")
dir_clustered <- file.path(repo_root, "data", "processed", "gse245601", "seurat_clustered")
dir.create(dir_clustered, recursive = TRUE, showWarnings = FALSE)

rds_files <- list.files(dir_seurat_qc, pattern = "\\.rds$", full.names = TRUE)
stopifnot(length(rds_files) == 20)

message("Loading ", length(rds_files), " QC'd sample objects ...")
obj_list <- lapply(rds_files, readRDS)

merged <- if (length(obj_list) > 1) {
  merge(obj_list[[1]], y = obj_list[-1])
} else {
  obj_list[[1]]
}
merged <- JoinLayers(merged)

message("Merged object: ", ncol(merged), " cells, ", nrow(merged), " genes")

merged <- NormalizeData(merged, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
merged <- FindVariableFeatures(merged, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
merged <- ScaleData(merged, vars.to.regress = "percent.mt", verbose = FALSE)
merged <- RunPCA(merged, npcs = 50, seed.use = 42, verbose = FALSE)
merged <- FindNeighbors(merged, dims = 1:30, verbose = FALSE)
merged <- FindClusters(merged, resolution = 0.8, verbose = FALSE)
merged <- RunUMAP(merged, dims = 1:30, n.neighbors = 30, min.dist = 0.3, metric = "cosine", seed.use = 42, verbose = FALSE)

saveRDS(merged, file.path(dir_clustered, "merged.rds"))

cluster_summary <- as.data.frame(table(cluster = merged$seurat_clusters, patient = merged$patient, condition = merged$condition))
colnames(cluster_summary) <- c("cluster", "patient", "condition", "n_cells")
cluster_summary <- cluster_summary[cluster_summary$n_cells > 0, ]

output_tsv <- file.path(repo_root, "results", "tables", "gse245601_cluster_summary.tsv")
write.table(cluster_summary, output_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
message(sprintf("Wrote %s (%d rows); %d total clusters, %d total cells", output_tsv, nrow(cluster_summary), length(unique(merged$seurat_clusters)), ncol(merged)))
