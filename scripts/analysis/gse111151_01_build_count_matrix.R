#!/usr/bin/env Rscript
# GSE111151 -- build a genes x samples raw-count matrix from the 11
# per-sample supplementary files (docs/GSE111151_DATA_AUDIT.md section 4).
# Each file has identical Ensembl-gene-ID row order (verified in the data
# audit); this script re-verifies that identity itself before merging,
# rather than assuming it, and uses raw counts only (column 3) -- not the
# publisher's precomputed CPM columns (docs/GSE111151_PREANALYSIS.md
# section E).
#
# Input:  data/raw/gse111151/extracted/*.txt.gz
# Output: results/tables/gse111151/counts_matrix.tsv.gz (committed)
#         results/tables/gse111151/sample_metadata.tsv (committed)

suppressPackageStartupMessages({
  library(yaml)
  library(data.table)
})

args <- commandArgs(trailingOnly = FALSE)
repo_root <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", args, value = TRUE))), "..", ".."))
if (length(repo_root) == 0 || !dir.exists(file.path(repo_root, "config"))) repo_root <- getwd()

config <- yaml::read_yaml(file.path(repo_root, "config", "config.yaml"))
cfg <- config$gse111151
samples <- cfg$samples

raw_dir <- file.path(repo_root, cfg$raw$per_sample_dir)

message(sprintf("Reading %d per-sample files...", length(samples)))
ref_ids <- NULL
ref_names <- NULL
counts_list <- list()
lib_sizes <- numeric(length(samples))
sample_ids <- character(length(samples))

for (i in seq_along(samples)) {
  s <- samples[[i]]
  path <- file.path(raw_dir, s$filename)
  dt <- read.delim(gzfile(path), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  dt$counts <- as.integer(dt$counts)
  stopifnot(nrow(dt) == 60619)
  stopifnot(!anyDuplicated(dt$EnsEMBL_GenID))

  if (is.null(ref_ids)) {
    ref_ids <- dt$EnsEMBL_GenID
    ref_names <- dt$gene_name
  } else {
    stopifnot(identical(dt$EnsEMBL_GenID, ref_ids)) # verify identical gene ordering across every file, not assumed
    stopifnot(identical(dt$gene_name, ref_names))
  }

  counts_list[[s$sample_id]] <- dt$counts
  lib_sizes[i] <- sum(dt$counts)
  sample_ids[i] <- s$sample_id
  message(sprintf("  %s (%s): %d genes, library size %d", s$sample_id, s$gsm, nrow(dt), sum(dt$counts)))
}

counts_mat <- as.data.frame(counts_list, check.names = FALSE)
stopifnot(nrow(counts_mat) == 60619, ncol(counts_mat) == length(samples))
stopifnot(!anyNA(counts_mat))

out_counts <- data.table(gene_id = ref_ids, gene_name = ref_names, counts_mat)
dir.create(file.path(repo_root, "results", "tables", "gse111151"), recursive = TRUE, showWarnings = FALSE)
fwrite(out_counts, file.path(repo_root, cfg$output$counts_tsv), sep = "\t")
message(sprintf("wrote %s (%d genes x %d samples)", cfg$output$counts_tsv, nrow(out_counts), length(samples)))

metadata <- data.frame(
  sample_id = sample_ids,
  gsm = sapply(samples, function(s) s$gsm),
  cell_line = sapply(samples, function(s) s$parental_line),
  resistance_status = sapply(samples, function(s) s$status),
  paired_parental_sample_id = sapply(samples, function(s) if (is.null(s$paired_parental_sample_id)) NA_character_ else s$paired_parental_sample_id),
  library_size = lib_sizes,
  n_detected_genes = sapply(counts_list, function(x) sum(x > 0)),
  stringsAsFactors = FALSE
)
write.table(metadata, file.path(repo_root, cfg$output$metadata_tsv), sep = "\t", quote = FALSE, row.names = FALSE)
message(sprintf("wrote %s", cfg$output$metadata_tsv))
print(metadata)
message("Done. Every gene row present in every file, identical order verified, 0 NAs, 0 duplicate Ensembl IDs.")
