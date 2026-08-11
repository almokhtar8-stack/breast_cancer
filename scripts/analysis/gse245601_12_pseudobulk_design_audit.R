#!/usr/bin/env Rscript
# GSE245601 -- pseudobulk design/eligibility audit (gse245601_PREANALYSIS.md
# section 13 follow-up). Read-only: loads the frozen Seurat object's
# METADATA and gene NAMES only (no count/expression value is read or
# inspected anywhere in this script) plus the frozen InferCNV label and
# eligibility tables. Verifies every design invariant the pseudobulk phase
# depends on BEFORE any pseudobulk aggregation or candidate-gene
# expression is touched, and writes the result to an audit table. Stops
# (does not silently continue) if any check fails.
#
# Checks:
#   1. sample_id set == exactly Tumor_01..10 x {Control, Tamoxifen} (20).
#   2. Each patient has exactly one Control and one Tamoxifen sample_id.
#   3. Per-sample epithelial cell counts (from the Seurat object's own
#      broad_cell_type column) reproduce the frozen
#      gse245601_malignant_summary_per_sample.tsv n_epithelial_cells
#      column exactly.
#   4. Per-sample strict-malignant cell counts (from the frozen
#      gse245601_malignant_cell_labels.tsv) reproduce the frozen
#      gse245601_malignant_summary_per_sample.tsv n_malignant column
#      exactly.
#   5. Track B eligibility (>=50 malignant cells in BOTH arms), recomputed
#      independently from those frozen malignant counts, reproduces the
#      frozen gse245601_pair_eligibility.tsv eligible_for_pseudobulk
#      column exactly, and the eligible set equals exactly
#      {Tumor_02, Tumor_03, Tumor_07} (config track_b_eligible_patients).
#   6. Each of the 13 candidate genes + PAICS is present exactly once in
#      the Seurat object's gene (row) names -- a mapping/uniqueness check
#      only; no expression value for any gene is read.
#
# Input:  data/processed/gse245601/seurat_clustered/annotated.rds (metadata + gene names only)
#         results/tables/gse245601_malignant_cell_labels.tsv
#         results/tables/gse245601_malignant_summary_per_sample.tsv
#         results/tables/gse245601_pair_eligibility.tsv
# Output: results/tables/gse245601_pseudobulk/design_eligibility_audit.tsv (committed)

suppressPackageStartupMessages({
  library(Seurat)
  library(yaml)
})

args <- commandArgs(trailingOnly = FALSE)
repo_root <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", args, value = TRUE))), "..", ".."))
if (length(repo_root) == 0 || !dir.exists(file.path(repo_root, "config"))) repo_root <- getwd()

config <- yaml::read_yaml(file.path(repo_root, "config", "config.yaml"))
cfg <- config$gse245601_pseudobulk

candidates_13 <- cfg$candidates$thirteen
paics <- cfg$candidates$paics
min_per_arm <- as.integer(cfg$min_malignant_cells_per_arm)
expected_track_b <- sort(cfg$track_b_eligible_patients)
audit_out <- file.path(repo_root, cfg$output$design_audit_tsv)
dir.create(dirname(audit_out), recursive = TRUE, showWarnings = FALSE)

merged <- readRDS(file.path(repo_root, cfg$inputs$seurat_object))
malignant_labels <- read.delim(file.path(repo_root, cfg$inputs$infercnv_labels_tsv), stringsAsFactors = FALSE)
malignant_summary <- read.delim(file.path(repo_root, cfg$inputs$malignant_summary_tsv), stringsAsFactors = FALSE)
pair_eligibility <- read.delim(file.path(repo_root, cfg$inputs$pair_eligibility_tsv), stringsAsFactors = FALSE)

checks <- list()
add_check <- function(name, expected, actual, pass) {
  checks[[length(checks) + 1]] <<- data.frame(
    check = name,
    expected = as.character(expected),
    actual = as.character(actual),
    pass = pass,
    stringsAsFactors = FALSE
  )
}

# --- 1. sample_id set ---
expected_samples <- sort(as.vector(outer(sprintf("Tumor_%02d", 1:10), c("Control", "Tamoxifen"), paste, sep = "_")))
actual_samples <- sort(unique(merged$sample_id))
add_check(
  "sample_id_set_is_exactly_20_tumor_pairs",
  paste(expected_samples, collapse = ","),
  paste(actual_samples, collapse = ","),
  identical(expected_samples, actual_samples)
)

# --- 2. exactly one Control + one Tamoxifen per patient ---
tab <- table(merged$patient, merged$condition)
one_each <- all(tab[, "Control"] > 0) && all(tab[, "Tamoxifen"] > 0) &&
  length(unique(merged$sample_id[merged$condition == "Control"])) == 10 &&
  length(unique(merged$sample_id[merged$condition == "Tamoxifen"])) == 10
add_check("exactly_one_control_and_one_tamoxifen_sample_per_patient", TRUE, one_each, isTRUE(one_each))

# --- 3. epithelial cell counts reproduce frozen summary ---
epi_counts <- table(merged$sample_id[merged$broad_cell_type == "epithelial"])
epi_counts <- epi_counts[malignant_summary$sample_id]
epi_match <- identical(as.integer(epi_counts), as.integer(malignant_summary$n_epithelial_cells))
add_check(
  "epithelial_cell_counts_reproduce_frozen_summary_exactly",
  paste(malignant_summary$n_epithelial_cells, collapse = ","),
  paste(as.integer(epi_counts), collapse = ","),
  epi_match
)

# --- 4. strict malignant counts reproduce frozen summary ---
malignant_counts_recomputed <- sapply(malignant_summary$sample_id, function(s) {
  sum(malignant_labels$sample_id == s & malignant_labels$primary_malignancy_label == "malignant")
})
malignant_match <- identical(as.integer(malignant_counts_recomputed), as.integer(malignant_summary$n_malignant))
add_check(
  "strict_malignant_counts_reproduce_frozen_summary_exactly",
  paste(malignant_summary$n_malignant, collapse = ","),
  paste(as.integer(malignant_counts_recomputed), collapse = ","),
  malignant_match
)

# --- 5. Track B eligibility ---
malignant_summary$patient <- sub("_(Control|Tamoxifen)$", "", malignant_summary$sample_id)
malignant_summary$condition <- sub("^Tumor_[0-9]+_", "", malignant_summary$sample_id)
patients <- sort(unique(malignant_summary$patient))
control_n <- setNames(malignant_summary$n_malignant[malignant_summary$condition == "Control"],
                       malignant_summary$patient[malignant_summary$condition == "Control"])
tamoxifen_n <- setNames(malignant_summary$n_malignant[malignant_summary$condition == "Tamoxifen"],
                         malignant_summary$patient[malignant_summary$condition == "Tamoxifen"])
recomputed_eligible <- (control_n[patients] >= min_per_arm) & (tamoxifen_n[patients] >= min_per_arm)
names(recomputed_eligible) <- patients

frozen_eligible <- setNames(pair_eligibility$eligible_for_pseudobulk, pair_eligibility$patient)
eligibility_match <- identical(as.logical(recomputed_eligible[patients]), as.logical(frozen_eligible[patients]))
add_check(
  "recomputed_track_b_eligibility_matches_frozen_pair_eligibility_table",
  paste(sprintf("%s=%s", patients, frozen_eligible[patients]), collapse = ";"),
  paste(sprintf("%s=%s", patients, recomputed_eligible[patients]), collapse = ";"),
  eligibility_match
)

actual_track_b <- sort(patients[recomputed_eligible])
add_check(
  "track_b_eligible_set_is_exactly_configured_patients",
  paste(expected_track_b, collapse = ","),
  paste(actual_track_b, collapse = ","),
  identical(expected_track_b, actual_track_b)
)

# --- 6. candidate + PAICS unique mapping in the raw gene feature space ---
# NOTE: presence/uniqueness of the gene NAME only -- no expression value
# for any gene is read here.
gene_names <- rownames(merged)
for (g in c(candidates_13, paics)) {
  n_found <- sum(gene_names == g)
  add_check(sprintf("candidate_gene_unique_in_feature_space__%s", g), 1, n_found, n_found == 1)
}

audit <- do.call(rbind, checks)
write.table(audit, audit_out, sep = "\t", quote = FALSE, row.names = FALSE)
message(sprintf("Wrote %s (%d checks, %d passed, %d failed)", audit_out, nrow(audit), sum(audit$pass), sum(!audit$pass)))

if (any(!audit$pass)) {
  print(audit[!audit$pass, ])
  stop("gse245601_12_pseudobulk_design_audit.R: one or more design/eligibility checks failed -- see printed rows above. Refusing to proceed to pseudobulk aggregation on a mismatched premise.")
}
message("All design/eligibility checks passed. No expression value or candidate-gene treatment effect was inspected in this script.")
