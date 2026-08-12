#!/usr/bin/env Rscript
# GSE240112 -- tumor-cell pseudobulk construction (Phase 9;
# docs/GSE240112_PREANALYSIS.md section F). Primary population = the
# author-defined tumor/malignant cells (TTs_cancer_060223.h5seurat, in
# which cell.annot == "Breast cancer cells" for all cells -- Phase 7
# Case A). Raw UMI counts summed per sample (orig.ident: PT1-3/RT1-3).
# Individual cells are never treated as biological replicates -- one
# pseudobulk profile per sample, six total.
#
# Input:  data/processed/gse240112/tt_cancer_extracted.rds
# Output: results/tables/gse240112_pseudobulk/tumor_cell_counts.tsv.gz (committed)
#         results/tables/gse240112_pseudobulk/tumor_cell_metadata.tsv (committed)

suppressPackageStartupMessages({
  library(Matrix)
  library(yaml)
  library(data.table)
})

args <- commandArgs(trailingOnly = FALSE)
repo_root <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", args, value = TRUE))), "..", ".."))
if (length(repo_root) == 0 || !dir.exists(file.path(repo_root, "config"))) repo_root <- getwd()

config <- yaml::read_yaml(file.path(repo_root, "config", "config.yaml"))
cfg <- config$gse240112

dir.create(dirname(file.path(repo_root, cfg$output$tumor_cell$counts_tsv)), recursive = TRUE, showWarnings = FALSE)

tt <- readRDS(file.path(repo_root, cfg$output$tt_cancer_rds))
counts <- tt$counts
meta <- tt$meta
stopifnot(identical(colnames(counts), meta$cell_id))
stopifnot(all(meta$cell.annot == "Breast cancer cells")) # Phase 7 Case A population check -- every aggregated cell is author-labeled tumor

sample_groups <- split(meta$cell_id, meta$orig.ident)
stopifnot(length(sample_groups) == 6)
stopifnot(setequal(names(sample_groups), c("PT1", "PT2", "PT3", "RT1", "RT2", "RT3")))
stopifnot(sum(lengths(sample_groups)) == nrow(meta)) # every tumor cell assigned to exactly one pseudobulk sample, none lost or duplicated

pb_mat <- sapply(sample_groups, function(cell_ids) Matrix::rowSums(counts[, cell_ids, drop = FALSE]))
rownames(pb_mat) <- rownames(counts)
pb_mat <- pb_mat[, c("PT1", "PT2", "PT3", "RT1", "RT2", "RT3")] # fixed column order

group_of <- substr(colnames(pb_mat), 1, 2)
pb_meta <- data.frame(
  sample_id = colnames(pb_mat),
  group = group_of,
  n_contributing_cells = lengths(sample_groups)[colnames(pb_mat)],
  total_library_size = colSums(pb_mat),
  n_detected_genes = colSums(pb_mat > 0),
  stringsAsFactors = FALSE
)

out_counts <- data.table(gene = rownames(pb_mat), as.data.frame(pb_mat))
fwrite(out_counts, file.path(repo_root, cfg$output$tumor_cell$counts_tsv), sep = "\t")
write.table(pb_meta, file.path(repo_root, cfg$output$tumor_cell$metadata_tsv), sep = "\t", quote = FALSE, row.names = FALSE)

message(sprintf("wrote tumor-cell pseudobulk: %d genes x %d samples", nrow(pb_mat), ncol(pb_mat)))
message(sprintf("total tumor cells aggregated: %d (matches source: %s)", sum(pb_meta$n_contributing_cells), sum(pb_meta$n_contributing_cells) == nrow(meta)))
print(pb_meta)
message("Done.")
