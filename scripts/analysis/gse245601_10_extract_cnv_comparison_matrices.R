#!/usr/bin/env Rscript
# GSE245601 -- visualization prep for the matched InferCNV/CopyKAT heatmap
# comparison (docs/CNV_METHOD_AUDIT.md Point 1 follow-up). Read-only
# against already-frozen InferCNV (run.final.infercnv_obj) and CopyKAT
# (<sample>_copykat_raw_results_gene_by_cell.txt) outputs for the 6
# priority tumors (config gse245601_cnv_heatmap_comparison.selected_samples)
# -- neither method is rerun here, no malignancy label is changed, no
# candidate gene is read anywhere in this script.
#
# For each selected sample:
#   1. Epithelial cell_id set = the frozen
#      results/tables/gse245601_malignant_cell_labels.tsv rows for that
#      sample (already verified, in the prior audit step, to be identical
#      to InferCNV's own observation-cell set and a subset of CopyKAT's
#      per-cell prediction).
#   2. InferCNV matrix = run.final.infercnv_obj@expr.data, columns
#      restricted to that cell_id set, rows = gene symbols.
#   3. CopyKAT matrix = <sample>_copykat_raw_results_gene_by_cell.txt,
#      columns restricted to that cell_id set via data.table::fread's
#      select= (the file is 60MB-860MB; select= avoids ever materializing
#      the full width in memory).
#   4. Gene symbols duplicated WITHIN either matrix (InferCNV rownames or
#      CopyKAT hgnc_symbol), or missing/blank, are excluded from BOTH
#      before intersecting -- a deterministic exclude-not-merge rule,
#      never a silent pick-first. The excluded-row counts (not distinct
#      symbol counts -- a symbol appearing on 3 rows contributes 3 to the
#      count) are written to gene_exclusion_report.tsv.
#   5. The shared, de-duplicated gene set is sorted by (numeric chromosome
#      1-22, start position) using InferCNV's own @gene_order as the
#      single position authority -- InferCNV's on-disk gene order is
#      chromosome-blocked but not chromosome-numeric (e.g. chr1, chr12,
#      chr2, ...), so this is a display-only re-sort applied identically
#      to both matrices, not a change to either method's output. CopyKAT's
#      own chromosome_name for the same genes is cross-checked against
#      InferCNV's (not silently trusted to agree); start_position is not
#      cross-checked -- InferCNV's start is used as the single ordering
#      authority without comparing it to CopyKAT's own start_position
#      column.
#
# Requires the sc245601 micromamba env (infercnv, data.table, yaml):
#   micromamba run -n sc245601 Rscript scripts/analysis/gse245601_10_extract_cnv_comparison_matrices.R
#
# Input:  data/processed/gse245601/infercnv/<sample>/run.final.infercnv_obj
#         data/processed/gse245601/copykat/<sample>/<sample>_copykat_raw_results_gene_by_cell.txt
#         results/tables/gse245601_malignant_cell_labels.tsv
# Output: data/processed/gse245601/cnv_heatmap_comparison/<sample>/infercnv_matrix.tsv.gz (gitignored)
#         data/processed/gse245601/cnv_heatmap_comparison/<sample>/copykat_matrix.tsv.gz (gitignored)
#         data/processed/gse245601/cnv_heatmap_comparison/<sample>/gene_order.tsv (gitignored)
#         results/tables/gse245601_cnv_heatmap_comparison/gene_exclusion_report.tsv (committed)

suppressPackageStartupMessages({
  library(infercnv)
  library(data.table)
  library(yaml)
})

CHR_MISMATCH_FRACTION_ERROR_THRESHOLD <- 0.05

args <- commandArgs(trailingOnly = FALSE)
repo_root <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", args, value = TRUE))), "..", ".."))
if (length(repo_root) == 0 || !dir.exists(file.path(repo_root, "config"))) repo_root <- getwd()

config <- yaml::read_yaml(file.path(repo_root, "config", "config.yaml"))
cfg <- config$gse245601_cnv_heatmap_comparison

selected_samples <- cfg$selected_samples
infercnv_dir <- file.path(repo_root, cfg$inputs$infercnv_dir)
copykat_dir <- file.path(repo_root, cfg$inputs$copykat_dir)
labels_path <- file.path(repo_root, cfg$inputs$infercnv_labels_tsv)
working_dir <- file.path(repo_root, cfg$working_dir)
gene_exclusion_out <- file.path(repo_root, cfg$output$gene_exclusion_tsv)

dir.create(working_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(gene_exclusion_out), recursive = TRUE, showWarnings = FALSE)

labels <- read.delim(labels_path, stringsAsFactors = FALSE)

exclusion_rows <- list()

for (sample_id in selected_samples) {
  message("=== ", sample_id, " ===")
  epi_cells <- labels$cell_id[labels$sample_id == sample_id]
  n_epi <- length(epi_cells)
  message(sprintf("  epithelial cells (frozen label table): %d", n_epi))
  if (n_epi == 0) stop(sprintf("%s: 0 epithelial cells in the frozen label table", sample_id))
  if (anyDuplicated(epi_cells) > 0) {
    stop(sprintf("%s: duplicate cell_id value(s) in the frozen label table for this sample", sample_id))
  }

  # --- InferCNV ---
  infercnv_obj_path <- file.path(infercnv_dir, sample_id, "run.final.infercnv_obj")
  obj <- readRDS(infercnv_obj_path)
  missing_from_rds <- setdiff(epi_cells, colnames(obj@expr.data))
  if (length(missing_from_rds) > 0) {
    stop(sprintf(
      "%s: %d epithelial cell(s) in the frozen label table are absent from run.final.infercnv_obj",
      sample_id, length(missing_from_rds)
    ))
  }
  ic_mat <- obj@expr.data[, epi_cells, drop = FALSE]
  ic_genes_all <- rownames(ic_mat)
  n_ic_genes_all <- length(ic_genes_all)
  ic_invalid <- is.na(ic_genes_all) | !nzchar(ic_genes_all)
  ic_dup <- duplicated(ic_genes_all) | duplicated(ic_genes_all, fromLast = TRUE) | ic_invalid
  n_ic_dup <- sum(ic_dup)
  ic_gene_order <- obj@gene_order # same row order as expr.data (chr, start, stop)

  # --- CopyKAT ---
  ck_path <- file.path(copykat_dir, sample_id, paste0(sample_id, "_copykat_raw_results_gene_by_cell.txt"))
  ck_header <- names(fread(ck_path, nrows = 0))
  id_cols <- c("hgnc_symbol", "chromosome_name", "start_position")
  missing_id_cols <- setdiff(id_cols, ck_header)
  if (length(missing_id_cols) > 0) {
    stop(sprintf("%s: CopyKAT gene_by_cell file missing expected columns: %s", sample_id, paste(missing_id_cols, collapse = ", ")))
  }
  missing_from_ck <- setdiff(epi_cells, ck_header)
  if (length(missing_from_ck) > 0) {
    stop(sprintf(
      "%s: %d epithelial cell(s) in the frozen label table are absent from the CopyKAT gene_by_cell file",
      sample_id, length(missing_from_ck)
    ))
  }
  ck_dt <- fread(ck_path, select = c(id_cols, epi_cells))
  ck_dt <- ck_dt[, c(id_cols, epi_cells), with = FALSE] # enforce column order regardless of fread's internal order
  ck_genes_all <- ck_dt$hgnc_symbol
  n_ck_genes_all <- length(ck_genes_all)
  ck_invalid <- is.na(ck_genes_all) | !nzchar(ck_genes_all)
  ck_dup <- duplicated(ck_genes_all) | duplicated(ck_genes_all, fromLast = TRUE) | ck_invalid
  n_ck_dup <- sum(ck_dup)

  # --- deterministic shared, de-duplicated, chromosome-sorted gene set ---
  ic_genes_unique <- ic_genes_all[!ic_dup]
  ck_genes_unique <- ck_genes_all[!ck_dup]
  shared_genes <- intersect(ic_genes_unique, ck_genes_unique)
  n_shared <- length(shared_genes)
  message(sprintf(
    "  genes: infercnv=%d (duplicate/invalid-symbol rows excluded=%d), copykat=%d (duplicate/invalid-symbol rows excluded=%d), shared after exclusion=%d",
    n_ic_genes_all, n_ic_dup, n_ck_genes_all, n_ck_dup, n_shared
  ))
  if (n_shared == 0) stop(sprintf("%s: 0 shared genes after duplicate exclusion", sample_id))

  chr_num <- suppressWarnings(as.integer(sub("^chr", "", ic_gene_order[shared_genes, "chr"])))
  start_pos <- ic_gene_order[shared_genes, "start"]
  if (any(is.na(chr_num))) stop(sprintf("%s: non-numeric chromosome value among shared genes (InferCNV side)", sample_id))
  if (any(is.na(start_pos))) stop(sprintf("%s: missing start position among shared genes (InferCNV side)", sample_id))
  ord <- order(chr_num, start_pos)
  shared_genes_sorted <- shared_genes[ord]

  # Cross-check InferCNV's chromosome assignment against CopyKAT's own, for
  # the same gene symbols -- do not silently trust the two annotation
  # sources agree.
  ck_chr_lookup <- setNames(as.integer(as.character(ck_dt$chromosome_name)), ck_dt$hgnc_symbol)
  ck_chr_for_shared <- unname(ck_chr_lookup[shared_genes_sorted])
  n_chr_mismatch <- sum(chr_num[ord] != ck_chr_for_shared, na.rm = TRUE) + sum(is.na(ck_chr_for_shared))
  mismatch_frac <- n_chr_mismatch / n_shared
  message(sprintf("  chromosome cross-check: %d/%d shared genes mismatch InferCNV vs CopyKAT (%.4f%%)", n_chr_mismatch, n_shared, 100 * mismatch_frac))
  if (mismatch_frac > CHR_MISMATCH_FRACTION_ERROR_THRESHOLD) {
    stop(sprintf(
      "%s: %.2f%% of shared genes have a chromosome disagreement between InferCNV and CopyKAT annotation -- exceeds the %.0f%% tolerance, refusing to proceed silently",
      sample_id, 100 * mismatch_frac, 100 * CHR_MISMATCH_FRACTION_ERROR_THRESHOLD
    ))
  }

  gene_order_df <- data.frame(
    gene = shared_genes_sorted,
    chr = chr_num[ord],
    start = start_pos[ord],
    stringsAsFactors = FALSE
  )

  # --- write compact per-sample outputs (same gene rows, same cell columns, same order, both matrices) ---
  sample_dir <- file.path(working_dir, sample_id)
  dir.create(sample_dir, recursive = TRUE, showWarnings = FALSE)

  ic_out <- as.data.frame(ic_mat[shared_genes_sorted, , drop = FALSE])
  ic_out <- cbind(gene = shared_genes_sorted, ic_out)
  fwrite(ic_out, file.path(sample_dir, "infercnv_matrix.tsv.gz"), sep = "\t")

  ck_mat_full <- as.matrix(ck_dt[, ..epi_cells])
  rownames(ck_mat_full) <- ck_dt$hgnc_symbol
  ck_out <- as.data.frame(ck_mat_full[shared_genes_sorted, epi_cells, drop = FALSE])
  ck_out <- cbind(gene = shared_genes_sorted, ck_out)
  fwrite(ck_out, file.path(sample_dir, "copykat_matrix.tsv.gz"), sep = "\t")

  fwrite(gene_order_df, file.path(sample_dir, "gene_order.tsv"), sep = "\t")

  exclusion_rows[[sample_id]] <- data.frame(
    sample_id = sample_id,
    n_epithelial_cells = n_epi,
    n_infercnv_genes = n_ic_genes_all,
    n_infercnv_duplicate_or_invalid_symbol_rows_excluded = n_ic_dup,
    n_copykat_genes = n_ck_genes_all,
    n_copykat_duplicate_or_invalid_symbol_rows_excluded = n_ck_dup,
    n_shared_genes_final = n_shared,
    n_chr_mismatch_infercnv_vs_copykat = n_chr_mismatch,
    stringsAsFactors = FALSE
  )

  rm(obj, ic_mat, ck_dt, ck_mat_full)
  gc(verbose = FALSE)
}

exclusion_report <- do.call(rbind, exclusion_rows)
rownames(exclusion_report) <- NULL
write.table(exclusion_report, gene_exclusion_out, sep = "\t", quote = FALSE, row.names = FALSE)
message(sprintf("Wrote %s (%d samples)", gene_exclusion_out, nrow(exclusion_report)))
