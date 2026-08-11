#!/usr/bin/env Rscript
# GSE245601 -- paired-design edgeR differential expression for Track A
# (all epithelial, n=10 patient pairs) and Track B (strict malignant,
# n=3 patient pairs, explicitly exploratory), per
# gse245601_PREANALYSIS.md section 13.
#
# Design: ~ patient + treatment (patient as a blocking factor so each
# patient is its own paired control; treatment coefficient = Tamoxifen
# vs Control, Control as the reference level). No batch covariate is
# added beyond patient -- nothing else is supported by the metadata.
#
# Method: standard modern edgeR workflow -- filterByExpr() (using the
# design matrix) to remove lowly-expressed genes, calcNormFactors() (TMM),
# estimateDisp(), glmQLFit()/glmQLFTest() (quasi-likelihood F-test) for
# the treatment coefficient. Genome-wide results are written in full (not
# filtered to significant genes) for context; no threshold is loosened for
# Track B to manufacture significance -- its small n is reported as-is.
#
# Input:  results/tables/gse245601_pseudobulk/track_a_epithelial_counts.tsv.gz (+ metadata)
#         results/tables/gse245601_pseudobulk/track_b_malignant_counts.tsv.gz (+ metadata)
# Output: results/tables/gse245601_pseudobulk/track_a_genomewide_de.tsv.gz (committed)
#         results/tables/gse245601_pseudobulk/track_b_genomewide_de.tsv.gz (committed)
#         results/tables/gse245601_pseudobulk/track_a_edger_filtering_summary.tsv (committed)
#         results/tables/gse245601_pseudobulk/track_b_edger_filtering_summary.tsv (committed)

suppressPackageStartupMessages({
  library(edgeR)
  library(yaml)
  library(data.table)
})

args <- commandArgs(trailingOnly = FALSE)
repo_root <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", args, value = TRUE))), "..", ".."))
if (length(repo_root) == 0 || !dir.exists(file.path(repo_root, "config"))) repo_root <- getwd()

config <- yaml::read_yaml(file.path(repo_root, "config", "config.yaml"))
cfg <- config$gse245601_pseudobulk

run_edger_track <- function(counts_path, metadata_path, de_out_path, summary_out_path, track_label, n_patient_pairs) {
  message(sprintf("=== %s (n=%d patient pairs) ===", track_label, n_patient_pairs))
  counts <- read.delim(gzfile(counts_path), stringsAsFactors = FALSE, check.names = FALSE)
  gene_names <- counts$gene
  counts <- as.matrix(counts[, setdiff(colnames(counts), "gene")])
  rownames(counts) <- gene_names
  metadata <- read.delim(metadata_path, stringsAsFactors = FALSE)
  stopifnot(identical(colnames(counts), metadata$sample_id))

  metadata$patient <- factor(metadata$patient)
  metadata$treatment <- factor(metadata$condition, levels = c("Control", "Tamoxifen"))
  design <- model.matrix(~ patient + treatment, data = metadata)

  dge <- DGEList(counts = counts, samples = metadata)
  n_genes_before <- nrow(dge)
  keep <- filterByExpr(dge, design = design)
  n_genes_after <- sum(keep)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- calcNormFactors(dge, method = "TMM")

  dge <- estimateDisp(dge, design)
  fit <- glmQLFit(dge, design)
  # treatment coefficient is the last column of the design matrix (Control is the reference level)
  treatment_coef <- grep("^treatmentTamoxifen$", colnames(design), value = TRUE)
  stopifnot(length(treatment_coef) == 1)
  qlf <- glmQLFTest(fit, coef = treatment_coef)

  results <- topTags(qlf, n = Inf, sort.by = "none")$table
  results$gene <- rownames(results)
  avg_log_cpm <- aveLogCPM(dge)
  results$avg_log_cpm <- avg_log_cpm[match(results$gene, rownames(dge))]
  out <- results[, c("gene", "logFC", "avg_log_cpm", "PValue", "FDR")]
  colnames(out) <- c("gene", "log2fc", "avg_log_cpm", "p_value", "fdr")
  out <- out[order(out$gene), ]

  dir.create(dirname(de_out_path), recursive = TRUE, showWarnings = FALSE)
  fwrite(out, de_out_path, sep = "\t")
  message(sprintf("  wrote %s (%d genes tested)", de_out_path, nrow(out)))

  filtering_summary <- data.frame(
    track = track_label,
    n_patient_pairs = n_patient_pairs,
    n_pseudobulk_samples = ncol(dge),
    n_genes_before_filterByExpr = n_genes_before,
    n_genes_after_filterByExpr = n_genes_after,
    normalization_method = "TMM",
    de_method = "edgeR glmQLFit/glmQLFTest",
    design_formula = "~ patient + treatment",
    contrast = "Tamoxifen vs Control",
    residual_df = fit$df.residual[1],
    stringsAsFactors = FALSE
  )
  write.table(filtering_summary, summary_out_path, sep = "\t", quote = FALSE, row.names = FALSE)
  message(sprintf("  wrote %s", summary_out_path))
  message(sprintf("  filtering: %d -> %d genes; residual df = %d", n_genes_before, n_genes_after, fit$df.residual[1]))

  invisible(out)
}

run_edger_track(
  file.path(repo_root, cfg$output$track_a$counts_tsv),
  file.path(repo_root, cfg$output$track_a$metadata_tsv),
  file.path(repo_root, cfg$output$de$track_a_genomewide_tsv),
  file.path(repo_root, cfg$output$de$track_a_edger_summary_tsv),
  "track_a_epithelial",
  10
)

run_edger_track(
  file.path(repo_root, cfg$output$track_b$counts_tsv),
  file.path(repo_root, cfg$output$track_b$metadata_tsv),
  file.path(repo_root, cfg$output$de$track_b_genomewide_tsv),
  file.path(repo_root, cfg$output$de$track_b_edger_summary_tsv),
  "track_b_malignant",
  length(cfg$track_b_eligible_patients)
)

message("Done.")
