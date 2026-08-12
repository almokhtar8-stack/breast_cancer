#!/usr/bin/env Rscript
# GSE240112 -- Phase 15 OPTIONAL/SECONDARY normal-tissue context. Purely
# descriptive: mean log-normalized expression of the testable candidates
# in NT (normal breast) epithelial cells, for comparison against the
# already-computed PT/RT epithelial pseudobulk means (Phase 14). No
# statistical test is run here (n=2 NT samples) and no "tumor
# specificity"/safety claim is made from it -- see
# docs/GSE240112_PREANALYSIS.md and Phase 15's explicit instruction not to
# let this distract from the primary RT-vs-PT question.
#
# USP34 is absent from the author-processed NTs.h5seurat object's own
# 25,543-gene feature set (a smaller/differently-filtered set than the
# 27,161-gene TTs_cancer object) -- confirmed directly, not assumed -- so
# USP34 cannot be included in this NT comparison; this is reported
# explicitly, not silently skipped.
#
# Epithelial cells in NT are called with a simple, bounded, non-tuned
# rule (not a malignancy classifier, not optimized against any candidate
# gene): a cell is "epithelial" if it has nonzero counts for at least 2
# of {EPCAM, KRT8, KRT18, KRT19} and zero counts for PTPRC (immune
# marker). This rule is not iterated on or adjusted after inspecting any
# candidate-gene result.
#
# Input:  data/processed/gse240112/nt_extracted.rds
# Output: results/tables/gse240112/nt_epithelial_candidate_means.tsv (committed)

suppressPackageStartupMessages({
  library(Matrix)
  library(yaml)
})

args <- commandArgs(trailingOnly = FALSE)
repo_root <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", args, value = TRUE))), "..", ".."))
if (length(repo_root) == 0 || !dir.exists(file.path(repo_root, "config"))) repo_root <- getwd()

config <- yaml::read_yaml(file.path(repo_root, "config", "config.yaml"))
cfg <- config$gse240112
candidates <- cfg$candidates$thirteen
paics_gene <- cfg$candidates$paics

nt <- readRDS(file.path(repo_root, cfg$output$nt_rds))
counts <- nt$counts
meta <- nt$meta
stopifnot(identical(colnames(counts), meta$cell_id))

epi_markers <- c("EPCAM", "KRT8", "KRT18", "KRT19")
epi_markers_present <- intersect(epi_markers, rownames(counts))
n_epi_markers_detected <- Matrix::colSums(counts[epi_markers_present, , drop = FALSE] > 0)
ptprc_detected <- if ("PTPRC" %in% rownames(counts)) counts["PTPRC", ] > 0 else rep(FALSE, ncol(counts))
is_epithelial <- (n_epi_markers_detected >= 2) & !ptprc_detected

message(sprintf("NT: %d/%d cells called epithelial by the bounded marker rule", sum(is_epithelial), ncol(counts)))
meta$is_epithelial <- is_epithelial

epi_counts <- counts[, is_epithelial, drop = FALSE]
epi_sample <- meta$orig.ident[is_epithelial]

lib_size <- Matrix::colSums(epi_counts)
all_genes <- c(candidates, paics_gene)
rows <- list()
for (g in all_genes) {
  if (!(g %in% rownames(epi_counts))) {
    rows[[g]] <- data.frame(gene = g, present_in_nt_feature_space = FALSE, n_epithelial_cells = sum(is_epithelial),
                             mean_log_norm_expression = NA_real_, pct_expressing = NA_real_, stringsAsFactors = FALSE)
    next
  }
  raw <- epi_counts[g, ]
  log_norm <- log1p(raw / lib_size * 1e4)
  rows[[g]] <- data.frame(gene = g, present_in_nt_feature_space = TRUE, n_epithelial_cells = sum(is_epithelial),
                           mean_log_norm_expression = mean(log_norm), pct_expressing = 100 * mean(raw > 0), stringsAsFactors = FALSE)
}
out <- do.call(rbind, rows)

out_path <- file.path(repo_root, "results", "tables", "gse240112", "nt_epithelial_candidate_means.tsv")
dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
write.table(out, out_path, sep = "\t", quote = FALSE, row.names = FALSE)
message(sprintf("wrote %s", out_path))
print(out)
message("Done.")
