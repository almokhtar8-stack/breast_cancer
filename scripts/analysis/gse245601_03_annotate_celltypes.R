#!/usr/bin/env Rscript
# GSE245601 -- Step 6: broad cell-type annotation (candidate-gene-blind).
#
# Strategy: per-cluster canonical-marker module scoring, preferring the
# authors' documented epithelial marker panel [A] (basal KRT5/6A/7/14/17,
# luminal KRT8/18/19 + FOXA1) over introducing an unrelated automated
# classifier. Non-epithelial lineage markers are standard, field-canonical
# marker sets [C] (frozen in docs/gse245601_PREANALYSIS.md section 7) since
# the authors' exact non-epithelial marker panels were not recoverable
# verbatim from the paper text or the excerpted code (which defers most
# non-epithelial typing to a SingleR/BlueprintEncodeData step we do not
# reproduce here -- see PREANALYSIS amendment). NONE of the 13 candidate
# genes or PAICS appears in any marker set used here.
#
# For each cluster, AddModuleScore is computed per lineage; the cluster is
# assigned to the lineage with the highest mean module score, UNLESS the
# top score does not clearly dominate (margin over the second-best score
# below MARGIN_THRESHOLD), in which case the cluster is labeled
# "uncertain" rather than forced into a lineage.
#
# Input:  data/processed/gse245601/seurat_clustered/merged.rds
# Output: data/processed/gse245601/seurat_clustered/annotated.rds (gitignored)
#         results/tables/gse245601_celltype_marker_evidence.tsv (committed)
#         results/tables/gse245601_celltype_counts_per_sample.tsv (committed)

suppressPackageStartupMessages({
  library(Seurat)
  library(yaml)
  library(dplyr)
})

set.seed(42)
MARGIN_THRESHOLD <- 0.05  # minimum module-score margin to call a lineage confidently

args <- commandArgs(trailingOnly = FALSE)
repo_root <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", args, value = TRUE))), "..", ".."))
if (length(repo_root) == 0 || !dir.exists(file.path(repo_root, "config"))) repo_root <- getwd()

dir_clustered <- file.path(repo_root, "data", "processed", "gse245601", "seurat_clustered")
merged <- readRDS(file.path(dir_clustered, "merged.rds"))

# --- Frozen marker panels (docs/gse245601_PREANALYSIS.md section 7) --------
marker_sets <- list(
  epithelial = c("KRT5", "KRT6A", "KRT7", "KRT14", "KRT17", "KRT8", "KRT18", "KRT19", "FOXA1"),
  t_nk = c("CD3D", "CD3E", "CD2", "NKG7", "GNLY"),
  b_plasma = c("MS4A1", "CD79A", "CD19", "JCHAIN", "MZB1"),
  myeloid = c("CD68", "LYZ", "ITGAM", "FCGR3A", "CSF1R"),
  fibroblast = c("COL1A1", "COL1A2", "PDGFRB", "DCN"),
  endothelial = c("PECAM1", "VWF", "CLDN5")
)

candidate_genes <- c("USP34", "CTDNEP1", "EIF4ENIF1", "HMGB1", "KDM1A", "PET117", "TADA2B", "VEZF1", "ICK", "SUPT4H1", "TLK2", "TSR3", "USP17L29", "PAICS")
stopifnot(!any(unlist(marker_sets) %in% candidate_genes))

for (lineage in names(marker_sets)) {
  genes_present <- intersect(marker_sets[[lineage]], rownames(merged))
  if (length(genes_present) == 0) next
  merged <- AddModuleScore(merged, features = list(genes_present), name = paste0("score_", lineage), seed = 42)
}
score_cols <- grep("^score_.*1$", colnames(merged@meta.data), value = TRUE)
names(score_cols) <- sub("^score_", "", sub("1$", "", score_cols))

cluster_scores <- merged@meta.data %>%
  group_by(seurat_clusters) %>%
  summarise(across(all_of(unname(score_cols)), mean), n_cells = n(), .groups = "drop")

lineage_call <- character(nrow(cluster_scores))
top_score <- numeric(nrow(cluster_scores))
second_score <- numeric(nrow(cluster_scores))
for (i in seq_len(nrow(cluster_scores))) {
  vals <- as.numeric(cluster_scores[i, unname(score_cols)])
  ord <- order(vals, decreasing = TRUE)
  top_score[i] <- vals[ord[1]]
  second_score[i] <- ifelse(length(vals) > 1, vals[ord[2]], -Inf)
  margin <- top_score[i] - second_score[i]
  lineage_call[i] <- if (margin >= MARGIN_THRESHOLD) names(score_cols)[ord[1]] else "uncertain"
}
cluster_scores$broad_cell_type <- lineage_call
cluster_scores$top_score <- top_score
cluster_scores$second_score <- second_score
cluster_scores$margin <- top_score - second_score

cluster_to_type <- setNames(cluster_scores$broad_cell_type, as.character(cluster_scores$seurat_clusters))
merged$broad_cell_type <- unname(cluster_to_type[as.character(merged$seurat_clusters)])

marker_evidence_rows <- list()
for (lineage in names(marker_sets)) {
  genes_present <- intersect(marker_sets[[lineage]], rownames(merged))
  genes_absent <- setdiff(marker_sets[[lineage]], rownames(merged))
  clusters_assigned <- cluster_scores$seurat_clusters[cluster_scores$broad_cell_type == lineage]
  marker_evidence_rows[[lineage]] <- data.frame(
    lineage = lineage,
    markers_used = paste(genes_present, collapse = ","),
    markers_absent_from_reference = paste(genes_absent, collapse = ","),
    n_clusters_assigned = length(clusters_assigned),
    clusters_assigned = paste(clusters_assigned, collapse = ","),
    n_cells_assigned = sum(cluster_scores$n_cells[cluster_scores$broad_cell_type == lineage]),
    stringsAsFactors = FALSE
  )
}
marker_evidence <- do.call(rbind, marker_evidence_rows)
uncertain_row <- data.frame(
  lineage = "uncertain", markers_used = NA, markers_absent_from_reference = NA,
  n_clusters_assigned = sum(cluster_scores$broad_cell_type == "uncertain"),
  clusters_assigned = paste(cluster_scores$seurat_clusters[cluster_scores$broad_cell_type == "uncertain"], collapse = ","),
  n_cells_assigned = sum(cluster_scores$n_cells[cluster_scores$broad_cell_type == "uncertain"]),
  stringsAsFactors = FALSE
)
marker_evidence <- rbind(marker_evidence, uncertain_row)

saveRDS(merged, file.path(dir_clustered, "annotated.rds"))

out1 <- file.path(repo_root, "results", "tables", "gse245601_celltype_marker_evidence.tsv")
write.table(marker_evidence, out1, sep = "\t", quote = FALSE, row.names = FALSE)
message(sprintf("Wrote %s (%d lineages)", out1, nrow(marker_evidence)))

counts_per_sample <- as.data.frame(table(sample_id = merged$sample_id, broad_cell_type = merged$broad_cell_type))
colnames(counts_per_sample) <- c("sample_id", "broad_cell_type", "n_cells")
counts_per_sample <- counts_per_sample[counts_per_sample$n_cells > 0, ]
out2 <- file.path(repo_root, "results", "tables", "gse245601_celltype_counts_per_sample.tsv")
write.table(counts_per_sample, out2, sep = "\t", quote = FALSE, row.names = FALSE)
message(sprintf("Wrote %s (%d rows)", out2, nrow(counts_per_sample)))

message("Cluster -> broad_cell_type summary:")
print(cluster_scores[, c("seurat_clusters", "broad_cell_type", "n_cells", "top_score", "margin")])
