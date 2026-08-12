#!/usr/bin/env Rscript
# GSE240112 -- unpaired edgeR differential expression, tumor-cell
# pseudobulk, recurrent (RT) vs primary (PT), per
# docs/GSE240112_PREANALYSIS.md sections C/F/H. Design: ~ group (PT
# reference level), contrast = groupRT (RT vs PT). No blocking factor --
# PT and RT samples come from different, unrelated patients (data audit
# section 2), so no patient/pair structure exists to model.
#
# Method: standard modern edgeR workflow -- filterByExpr() (using the
# design matrix) to remove lowly-expressed genes, calcNormFactors() (TMM),
# estimateDisp(), glmQLFit()/glmQLFTest() (quasi-likelihood F-test). No
# threshold is loosened or tightened after inspecting any result.
#
# Input:  results/tables/gse240112_pseudobulk/tumor_cell_counts.tsv.gz (+ metadata)
# Output: results/tables/gse240112_pseudobulk/tumor_cell_genomewide_de.tsv.gz (committed)
#         results/tables/gse240112_pseudobulk/tumor_cell_edger_filtering_summary.tsv (committed)

suppressPackageStartupMessages({
  library(edgeR)
  library(yaml)
  library(data.table)
})

args <- commandArgs(trailingOnly = FALSE)
repo_root <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", args, value = TRUE))), "..", ".."))
if (length(repo_root) == 0 || !dir.exists(file.path(repo_root, "config"))) repo_root <- getwd()

config <- yaml::read_yaml(file.path(repo_root, "config", "config.yaml"))
cfg <- config$gse240112

run_edger <- function(counts_path, metadata_path, de_out_path, summary_out_path, label) {
  message(sprintf("=== %s ===", label))
  counts <- read.delim(gzfile(counts_path), stringsAsFactors = FALSE, check.names = FALSE)
  gene_names <- counts$gene
  counts <- as.matrix(counts[, setdiff(colnames(counts), "gene")])
  rownames(counts) <- gene_names
  metadata <- read.delim(metadata_path, stringsAsFactors = FALSE)
  stopifnot(identical(colnames(counts), metadata$sample_id))

  metadata$group <- factor(metadata$group, levels = c("PT", "RT")) # PT = reference level
  design <- model.matrix(~group, data = metadata)

  dge <- DGEList(counts = counts, samples = metadata)
  n_genes_before <- nrow(dge)
  keep <- filterByExpr(dge, design = design)
  n_genes_after <- sum(keep)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- calcNormFactors(dge, method = "TMM")

  dge <- estimateDisp(dge, design)
  fit <- glmQLFit(dge, design)
  group_coef <- grep("^groupRT$", colnames(design), value = TRUE)
  stopifnot(length(group_coef) == 1)
  qlf <- glmQLFTest(fit, coef = group_coef)

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
    label = label,
    n_samples = ncol(dge),
    n_genes_before_filterByExpr = n_genes_before,
    n_genes_after_filterByExpr = n_genes_after,
    normalization_method = "TMM",
    de_method = "edgeR glmQLFit/glmQLFTest",
    design_formula = "~ group",
    contrast = "RT vs PT (RT - PT), unpaired",
    residual_df = fit$df.residual[1],
    stringsAsFactors = FALSE
  )
  write.table(filtering_summary, summary_out_path, sep = "\t", quote = FALSE, row.names = FALSE)
  message(sprintf("  wrote %s", summary_out_path))
  message(sprintf("  filtering: %d -> %d genes; residual df = %d", n_genes_before, n_genes_after, fit$df.residual[1]))

  invisible(out)
}

run_edger(
  file.path(repo_root, cfg$output$tumor_cell$counts_tsv),
  file.path(repo_root, cfg$output$tumor_cell$metadata_tsv),
  file.path(repo_root, cfg$output$de$genomewide_tsv),
  file.path(repo_root, cfg$output$de$edger_summary_tsv),
  "tumor_cell_RT_vs_PT"
)

message("Done.")
