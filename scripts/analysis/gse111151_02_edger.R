#!/usr/bin/env Rscript
# GSE111151 -- cell-line-blocked edgeR differential expression, tamoxifen-
# resistant vs parental, per docs/GSE111151_PREANALYSIS.md sections C/E/F/G.
# Design: ~ cell_line + resistance_status (cell_line as a 4-level blocking
# factor; each resistant subline compared against its own cell line's
# single parental sample). Not a naive unpaired two-group test -- that
# would confound cell-line identity with the resistance effect.
#
# Input:  results/tables/gse111151/counts_matrix.tsv.gz (+ sample_metadata.tsv)
# Output: results/tables/gse111151/genomewide_de.tsv.gz (committed)
#         results/tables/gse111151/edger_filtering_summary.tsv (committed)

suppressPackageStartupMessages({
  library(edgeR)
  library(yaml)
  library(data.table)
})

args <- commandArgs(trailingOnly = FALSE)
repo_root <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", args, value = TRUE))), "..", ".."))
if (length(repo_root) == 0 || !dir.exists(file.path(repo_root, "config"))) repo_root <- getwd()

config <- yaml::read_yaml(file.path(repo_root, "config", "config.yaml"))
cfg <- config$gse111151

counts_dt <- read.delim(gzfile(file.path(repo_root, cfg$output$counts_tsv)), header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
gene_ids <- counts_dt$gene_id
gene_names <- counts_dt$gene_name
counts <- as.matrix(counts_dt[, setdiff(colnames(counts_dt), c("gene_id", "gene_name"))])
rownames(counts) <- gene_ids
stopifnot(!anyDuplicated(rownames(counts)))

metadata <- read.delim(file.path(repo_root, cfg$output$metadata_tsv), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
stopifnot(identical(colnames(counts), metadata$sample_id))

metadata$cell_line <- factor(metadata$cell_line)
metadata$resistance_status <- factor(metadata$resistance_status, levels = c("parental", "resistant")) # parental = reference level
design <- model.matrix(~cell_line + resistance_status, data = metadata)
message("Design matrix:")
print(design)
stopifnot(qr(design)$rank == ncol(design)) # design must be full rank (no perfect confound between cell_line and resistance_status)

dge <- DGEList(counts = counts, samples = metadata)
n_genes_before <- nrow(dge)
keep <- filterByExpr(dge, design = design)
n_genes_after <- sum(keep)
dge <- dge[keep, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge, method = "TMM")

dge <- estimateDisp(dge, design)
fit <- glmQLFit(dge, design)
resistance_coef <- grep("^resistance_statusresistant$", colnames(design), value = TRUE)
stopifnot(length(resistance_coef) == 1)
qlf <- glmQLFTest(fit, coef = resistance_coef)

results <- topTags(qlf, n = Inf, sort.by = "none")$table
results$gene_id <- rownames(results)
avg_log_cpm <- aveLogCPM(dge)
results$avg_log_cpm <- avg_log_cpm[match(results$gene_id, rownames(dge))]
gene_name_map <- setNames(gene_names, gene_ids)
results$gene_name <- gene_name_map[results$gene_id]
out <- results[, c("gene_id", "gene_name", "logFC", "avg_log_cpm", "PValue", "FDR")]
colnames(out) <- c("gene_id", "gene_name", "log2fc", "avg_log_cpm", "p_value", "fdr")
out <- out[order(out$gene_id), ]

dir.create(dirname(file.path(repo_root, cfg$output$de$genomewide_tsv)), recursive = TRUE, showWarnings = FALSE)
fwrite(out, file.path(repo_root, cfg$output$de$genomewide_tsv), sep = "\t")
message(sprintf("wrote %s (%d genes tested)", cfg$output$de$genomewide_tsv, nrow(out)))

# Per-sample TMM normalization factors, exported so that Python-side
# sample-level descriptive values (used throughout Phase 7/8/11) use the
# same normalization as the statistical model, not a naive library-size-
# only CPM -- these differ substantially for some samples here (e.g. MCF-7
# parental TMM factor 1.30 vs MCF-7_Tam1 0.97), and using naive CPM was
# found to produce a spurious apparent direction for several candidates.
norm_factors <- data.frame(
  sample_id = colnames(dge),
  library_size = dge$samples$lib.size,
  norm_factor = dge$samples$norm.factors,
  effective_library_size = dge$samples$lib.size * dge$samples$norm.factors,
  stringsAsFactors = FALSE
)
norm_factors_path <- file.path(repo_root, cfg$output$de$tmm_norm_factors_tsv)
write.table(norm_factors, norm_factors_path, sep = "\t", quote = FALSE, row.names = FALSE)
message(sprintf("wrote %s", norm_factors_path))

filtering_summary <- data.frame(
  label = "gse111151_resistant_vs_parental",
  n_samples = ncol(dge),
  n_cell_lines = nlevels(metadata$cell_line),
  n_genes_before_filterByExpr = n_genes_before,
  n_genes_after_filterByExpr = n_genes_after,
  normalization_method = "TMM",
  de_method = "edgeR glmQLFit/glmQLFTest",
  design_formula = "~ cell_line + resistance_status",
  contrast = "resistant vs parental, blocked by cell_line",
  residual_df = fit$df.residual[1],
  stringsAsFactors = FALSE
)
write.table(filtering_summary, file.path(repo_root, cfg$output$de$edger_summary_tsv), sep = "\t", quote = FALSE, row.names = FALSE)
message(sprintf("wrote %s", cfg$output$de$edger_summary_tsv))
message(sprintf("filtering: %d -> %d genes; residual df = %d", n_genes_before, n_genes_after, fit$df.residual[1]))
message("Done.")
