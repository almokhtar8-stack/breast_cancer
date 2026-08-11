#!/usr/bin/env Rscript
# GSE245601 -- Step 9: independent malignancy-classification sensitivity
# check (CopyKAT), pre-registered in docs/gse245601_PREANALYSIS.md section
# 10 BEFORE any malignancy label -- primary or sensitivity -- was
# cross-referenced against any candidate gene.
#
# CopyKAT (Gao et al. 2021, Nat Biotechnol) is run PER SAMPLE on the same
# mixed cell population (reference immune/endothelial + epithelial) that
# InferCNV used as its input, using CopyKAT's OWN internal
# confident-normal-cell inference (we do not supply our marker-based
# immune/endothelial labels as `norm.cell.names`) -- this is what makes it
# a methodologically independent check rather than a re-run of the same
# reference assumption. Default parameters throughout (ngene.chr=5,
# win.size=25, KS.cut=0.1 (package default), distance="euclidean", genome="hg20").
#
# Concordance against the primary InferCNV-based labels is reported per
# sample and overall. Discordant cells are reported, not adjudicated.
# No candidate gene is read, plotted, or used anywhere in this script.
#
# Input:  data/processed/gse245601/seurat_clustered/annotated.rds
#         results/tables/gse245601_malignant_cell_labels.tsv (primary labels)
# Output: data/processed/gse245601/copykat/<sample_id>/ (gitignored)
#         results/tables/gse245601_copykat_sensitivity_labels.tsv (committed)
#         results/tables/gse245601_malignancy_concordance.tsv (committed)

suppressPackageStartupMessages({
  library(Seurat)
  library(copykat)
})

set.seed(42)
MIN_CELLS_FOR_COPYKAT <- 30

args <- commandArgs(trailingOnly = FALSE)
repo_root <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", args, value = TRUE))), "..", ".."))
if (length(repo_root) == 0 || !dir.exists(file.path(repo_root, "config"))) repo_root <- getwd()

merged <- readRDS(file.path(repo_root, "data", "processed", "gse245601", "seurat_clustered", "annotated.rds"))
primary_labels <- read.delim(file.path(repo_root, "results", "tables", "gse245601_malignant_cell_labels.tsv"), stringsAsFactors = FALSE)

immune_types <- c("t_nk", "b_plasma", "myeloid")
endothelial_type <- "endothelial"
epithelial_type <- "epithelial"

samples_with_primary_labels <- unique(primary_labels$sample_id)
message("Samples with primary InferCNV labels: ", length(samples_with_primary_labels))

dir_copykat_root <- file.path(repo_root, "data", "processed", "gse245601", "copykat")
dir.create(dir_copykat_root, recursive = TRUE, showWarnings = FALSE)

sensitivity_rows <- list()
concordance_rows <- list()

orig_wd <- getwd()

for (sample_id in samples_with_primary_labels) {
  message("=== ", sample_id, " (CopyKAT) ===")
  sample_cells_all <- colnames(merged)[merged$sample_id == sample_id]
  sample_obj <- subset(merged, cells = sample_cells_all)

  ref_cells <- colnames(sample_obj)[sample_obj$broad_cell_type %in% c(immune_types, endothelial_type)]
  epi_cells <- colnames(sample_obj)[sample_obj$broad_cell_type == epithelial_type]
  copykat_input_cells <- c(ref_cells, epi_cells)

  if (length(copykat_input_cells) < MIN_CELLS_FOR_COPYKAT) {
    sensitivity_rows[[sample_id]] <- NULL
    concordance_rows[[sample_id]] <- data.frame(sample_id = sample_id, status = "too_few_cells", n_concordant = NA, n_discordant = NA, n_compared = NA, concordance_rate = NA)
    next
  }

  raw_counts <- as.matrix(GetAssayData(sample_obj, layer = "counts")[, copykat_input_cells])

  sample_dir <- file.path(dir_copykat_root, sample_id)
  dir.create(sample_dir, recursive = TRUE, showWarnings = FALSE)
  setwd(sample_dir)

  ck_result <- tryCatch(
    copykat(
      rawmat = raw_counts, id.type = "S", cell.line = "no", ngene.chr = 5,
      win.size = 25, KS.cut = 0.1, sam.name = sample_id, distance = "euclidean",
      output.seg = FALSE, plot.genes = FALSE, genome = "hg20", n.cores = 8
    ),
    error = function(e) {
      message("  CopyKAT FAILED: ", conditionMessage(e))
      NULL
    }
  )
  setwd(orig_wd)

  if (is.null(ck_result) || is.null(ck_result$prediction)) {
    concordance_rows[[sample_id]] <- data.frame(sample_id = sample_id, status = "copykat_failed", n_concordant = NA, n_discordant = NA, n_compared = NA, concordance_rate = NA)
    next
  }

  pred <- ck_result$prediction
  colnames(pred)[colnames(pred) == "cell.names"] <- "cell_id"
  pred$sensitivity_malignancy_label <- ifelse(pred$copykat.pred == "aneuploid", "malignant",
                                        ifelse(pred$copykat.pred == "diploid", "non-malignant epithelial", "not_defined"))

  epi_pred <- pred[pred$cell_id %in% epi_cells, c("cell_id", "sensitivity_malignancy_label")]
  epi_pred$sample_id <- sample_id
  sensitivity_rows[[sample_id]] <- epi_pred

  merged_compare <- merge(
    primary_labels[primary_labels$sample_id == sample_id, c("cell_id", "primary_malignancy_label")],
    epi_pred[, c("cell_id", "sensitivity_malignancy_label")],
    by = "cell_id"
  )
  comparable <- merged_compare[merged_compare$sensitivity_malignancy_label != "not_defined", ]
  n_concordant <- sum(comparable$primary_malignancy_label == comparable$sensitivity_malignancy_label)
  n_compared <- nrow(comparable)
  concordance_rows[[sample_id]] <- data.frame(
    sample_id = sample_id, status = "ok",
    n_concordant = n_concordant, n_discordant = n_compared - n_concordant, n_compared = n_compared,
    concordance_rate = ifelse(n_compared > 0, n_concordant / n_compared, NA)
  )
  message(sprintf("  concordance: %d/%d (%.1f%%)", n_concordant, n_compared, 100 * n_concordant / max(n_compared, 1)))
}

sensitivity_labels <- do.call(rbind, sensitivity_rows)
out1 <- file.path(repo_root, "results", "tables", "gse245601_copykat_sensitivity_labels.tsv")
write.table(sensitivity_labels, out1, sep = "\t", quote = FALSE, row.names = FALSE)
message(sprintf("Wrote %s (%d cells)", out1, nrow(sensitivity_labels)))

concordance <- do.call(rbind, concordance_rows)
rownames(concordance) <- NULL
out2 <- file.path(repo_root, "results", "tables", "gse245601_malignancy_concordance.tsv")
write.table(concordance, out2, sep = "\t", quote = FALSE, row.names = FALSE)
message(sprintf("Wrote %s (%d samples)", out2, nrow(concordance)))

overall_concordant <- sum(concordance$n_concordant, na.rm = TRUE)
overall_compared <- sum(concordance$n_compared, na.rm = TRUE)
message(sprintf("OVERALL concordance: %d/%d (%.1f%%)", overall_concordant, overall_compared, 100 * overall_concordant / max(overall_compared, 1)))
