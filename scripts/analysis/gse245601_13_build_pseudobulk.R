#!/usr/bin/env Rscript
# GSE245601 -- Track A / Track B pseudobulk count construction
# (gse245601_PREANALYSIS.md section 13). Requires
# gse245601_12_pseudobulk_design_audit.R to have already passed (this
# script does not re-run those checks, but re-derives the same patient/
# condition/eligibility values from the same frozen sources so the two
# scripts cannot silently drift apart). Raw UMI counts only (the Seurat
# object's "counts" layer) -- never normalized/scaled expression. Patient
# x treatment is the pseudobulk unit; individual cells are never treated
# as independent replicates.
#
# Track A: every epithelial cell (broad_cell_type == "epithelial") from
#   all 10 paired tumors, summed per (patient, condition) -> 20 pseudobulk
#   samples.
# Track B: frozen InferCNV-malignant cells only
#   (primary_malignancy_label == "malignant" in the frozen per-cell label
#   table), restricted to the 3 tumor pairs already meeting the frozen
#   >=50-malignant-cells-per-arm rule (Tumor_02, Tumor_03, Tumor_07) ->
#   6 pseudobulk samples. No other tumor is included in Track B, and the
#   eligibility set is read from config (itself cross-checked against the
#   frozen gse245601_pair_eligibility.tsv in the audit script), not
#   re-decided here.
#
# Input:  data/processed/gse245601/seurat_clustered/annotated.rds
#         results/tables/gse245601_malignant_cell_labels.tsv
# Output: results/tables/gse245601_pseudobulk/track_a_epithelial_counts.tsv.gz (committed)
#         results/tables/gse245601_pseudobulk/track_a_epithelial_metadata.tsv (committed)
#         results/tables/gse245601_pseudobulk/track_b_malignant_counts.tsv.gz (committed)
#         results/tables/gse245601_pseudobulk/track_b_malignant_metadata.tsv (committed)

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
track_b_patients <- cfg$track_b_eligible_patients

for (p in c(cfg$output$track_a$counts_tsv, cfg$output$track_b$counts_tsv)) {
  dir.create(dirname(file.path(repo_root, p)), recursive = TRUE, showWarnings = FALSE)
}

merged <- readRDS(file.path(repo_root, cfg$inputs$seurat_object))
malignant_labels <- read.delim(file.path(repo_root, cfg$inputs$infercnv_labels_tsv), stringsAsFactors = FALSE)
counts <- GetAssayData(merged, layer = "counts")
stopifnot(is(counts, "dgCMatrix"))
stopifnot(!anyDuplicated(rownames(counts))) # gene symbols must be unique before any pseudobulk row-sum aggregation

build_pseudobulk <- function(cell_groups) {
  # cell_groups: named list, name = pseudobulk sample_id ("Tumor_XX_Condition"), value = cell_id vector
  pb_mat <- sapply(cell_groups, function(cell_ids) {
    stopifnot(length(cell_ids) > 0, all(cell_ids %in% colnames(counts)))
    Matrix::rowSums(counts[, cell_ids, drop = FALSE])
  })
  rownames(pb_mat) <- rownames(counts)
  meta <- data.frame(
    sample_id = names(cell_groups),
    patient = sub("_(Control|Tamoxifen)$", "", names(cell_groups)),
    condition = sub("^Tumor_[0-9]+_", "", names(cell_groups)),
    n_contributing_cells = sapply(cell_groups, length),
    total_library_size = colSums(pb_mat),
    n_detected_genes = colSums(pb_mat > 0),
    stringsAsFactors = FALSE
  )
  list(counts = pb_mat, metadata = meta)
}

write_pseudobulk <- function(pb, counts_path, metadata_path) {
  out <- data.table(gene = rownames(pb$counts), as.data.frame(pb$counts))
  fwrite(out, counts_path, sep = "\t")
  write.table(pb$metadata, metadata_path, sep = "\t", quote = FALSE, row.names = FALSE)
  message(sprintf("  wrote %s (%d genes x %d samples) and %s", counts_path, nrow(pb$counts), ncol(pb$counts), metadata_path))
}

# --- Track A: all epithelial cells, all 10 patients ---
message("=== Track A (all epithelial) ===")
epi_meta <- merged@meta.data[merged$broad_cell_type == "epithelial", c("cell_id", "sample_id")]
track_a_groups <- split(epi_meta$cell_id, epi_meta$sample_id)
stopifnot(length(track_a_groups) == 20)
stopifnot(sum(lengths(track_a_groups)) == nrow(epi_meta)) # every epithelial cell assigned to exactly one group
track_a <- build_pseudobulk(track_a_groups)
message(sprintf("  total epithelial cells aggregated: %d (matches Seurat epithelial count: %s)",
                 sum(track_a$metadata$n_contributing_cells), sum(track_a$metadata$n_contributing_cells) == nrow(epi_meta)))
write_pseudobulk(
  track_a,
  file.path(repo_root, cfg$output$track_a$counts_tsv),
  file.path(repo_root, cfg$output$track_a$metadata_tsv)
)

# --- Track B: frozen malignant cells, eligible patients only ---
message("=== Track B (strict malignant, eligible pairs only) ===")
malignant_cells <- malignant_labels[malignant_labels$primary_malignancy_label == "malignant", c("cell_id", "sample_id", "patient")]
malignant_cells <- malignant_cells[malignant_cells$patient %in% track_b_patients, ]
stopifnot(identical(sort(unique(malignant_cells$patient)), sort(track_b_patients)))
track_b_groups <- split(malignant_cells$cell_id, malignant_cells$sample_id)
stopifnot(length(track_b_groups) == length(track_b_patients) * 2)
stopifnot(sum(lengths(track_b_groups)) == nrow(malignant_cells))
# every Track B cell must also be one of the Track A epithelial cells for its sample (malignant cells are a subset of epithelial cells)
stopifnot(all(malignant_cells$cell_id %in% epi_meta$cell_id))
track_b <- build_pseudobulk(track_b_groups)
message(sprintf("  total malignant cells aggregated: %d", sum(track_b$metadata$n_contributing_cells)))
write_pseudobulk(
  track_b,
  file.path(repo_root, cfg$output$track_b$counts_tsv),
  file.path(repo_root, cfg$output$track_b$metadata_tsv)
)

message("Done. Track A and Track B pseudobulk counts and metadata written.")
