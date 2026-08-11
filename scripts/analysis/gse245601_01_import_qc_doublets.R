#!/usr/bin/env Rscript
# GSE245601 -- Step 3+4: Seurat import, QC, and doublet detection.
#
# Data source: GSE245601 (Kim, Whitman et al., Clin Cancer Res
# 2023;29(23):4894-4907, PMID 37747807), Cell Ranger filtered
# feature-barcode H5 matrices, GRCh38 10x reference (33,538 features,
# confirmed identical across all 20 primary tumor samples).
#
# QC and doublet-detection rules are frozen in docs/gse245601_PREANALYSIS.md
# BEFORE this script is run: nCount_RNA >= 5000, nFeature_RNA >= 2000,
# percent.mt <= 25 (paper ABlock3 thresholds); doublets removed via
# DoubletFinder (paper-stated method). Treatment condition (Control vs
# Tamoxifen) is never used as a QC criterion -- identical thresholds are
# applied to every cell regardless of condition. Raw UMI counts are
# preserved unmodified in a separate assay layer for future pseudobulk use.
#
# This script does not read, plot, or otherwise use any of the 13
# candidate genes or PAICS. It performs no differential expression, no
# pseudobulk, and no candidate-level analysis of any kind.
#
# Outputs:
#   data/processed/gse245601/seurat_qc/<patient>_<condition>.rds  (gitignored, not committed)
#   results/tables/gse245601_qc_summary.tsv  (committed)
#   results/tables/gse245601_session_info.txt (committed)

suppressPackageStartupMessages({
  library(Seurat)
  library(DoubletFinder)
  library(yaml)
  library(dplyr)
})

set.seed(42)

repo_root <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE))), "..", ".."))
if (length(repo_root) == 0 || !dir.exists(file.path(repo_root, "config"))) {
  repo_root <- getwd()
}
config <- yaml::read_yaml(file.path(repo_root, "config", "config.yaml"))
cfg <- config$gse245601

manifest <- read.delim(file.path(repo_root, cfg$output$sample_manifest_tsv), stringsAsFactors = FALSE)
primary <- manifest[manifest$primary_analysis == "True" | manifest$primary_analysis == TRUE, ]
stopifnot(nrow(primary) == 20)

# --- Frozen QC thresholds (docs/gse245601_PREANALYSIS.md section 4) --------
QC_MIN_NCOUNT <- 5000
QC_MIN_NFEATURE <- 2000
QC_MAX_PERCENT_MT <- 25
MITO_PATTERN <- "^MT-"

dir_seurat_qc <- file.path(repo_root, "data", "processed", "gse245601", "seurat_qc")
dir.create(dir_seurat_qc, recursive = TRUE, showWarnings = FALSE)

qc_rows <- list()

for (i in seq_len(nrow(primary))) {
  row <- primary[i, ]
  sample_id <- paste(row$patient, row$condition, sep = "_")
  h5_path <- file.path(repo_root, cfg$h5_dir, row$h5_filename)
  message(sprintf("[%d/%d] %s (%s)", i, nrow(primary), sample_id, row$GSM))

  counts <- Read10X_h5(h5_path)
  obj <- CreateSeuratObject(counts = counts, project = sample_id, min.cells = 0, min.features = 0)

  n_cells_raw <- ncol(obj)

  obj$GSM <- row$GSM
  obj$patient <- row$patient
  obj$condition <- row$condition
  obj$paired_patient <- row$patient
  obj$sample_id <- sample_id
  obj$original_barcode <- colnames(obj)
  obj <- RenameCells(obj, new.names = paste0(sample_id, "_", colnames(obj)))
  obj$cell_id <- colnames(obj)

  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = MITO_PATTERN)

  # --- QC filter (frozen thresholds; identical for Control and Tamoxifen) --
  n_fail_ncount <- sum(obj$nCount_RNA < QC_MIN_NCOUNT)
  n_fail_nfeature <- sum(obj$nFeature_RNA < QC_MIN_NFEATURE)
  n_fail_mt <- sum(obj$percent.mt > QC_MAX_PERCENT_MT)

  keep <- obj$nCount_RNA >= QC_MIN_NCOUNT & obj$nFeature_RNA >= QC_MIN_NFEATURE & obj$percent.mt <= QC_MAX_PERCENT_MT
  n_cells_postqc <- sum(keep)
  obj_qc <- subset(obj, cells = colnames(obj)[keep])

  # --- Doublet detection (DoubletFinder), post-QC ---------------------------
  n_doublets <- NA_integer_
  n_cells_final <- n_cells_postqc
  doublet_error <- NA_character_

  if (n_cells_postqc >= 50) {
    result <- tryCatch(
      {
        o <- NormalizeData(obj_qc, verbose = FALSE)
        o <- FindVariableFeatures(o, nfeatures = 2000, verbose = FALSE)
        o <- ScaleData(o, vars.to.regress = "percent.mt", verbose = FALSE)
        o <- RunPCA(o, npcs = 30, seed.use = 42, verbose = FALSE)
        o <- FindNeighbors(o, dims = 1:30, verbose = FALSE)
        o <- FindClusters(o, resolution = 0.8, verbose = FALSE)

        sweep_res <- paramSweep(o, PCs = 1:30, sct = FALSE)
        sweep_stats <- summarizeSweep(sweep_res, GT = FALSE)
        pdf(file = NULL)  # find.pK() plots as a side effect; suppress in batch mode
        bcmvn <- find.pK(sweep_stats)
        dev.off()
        pK_opt <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))

        annotations <- o$seurat_clusters
        homotypic_prop <- modelHomotypic(annotations)
        # Standard 10x multiplet-rate scaling: ~0.8% doublets per 1000 cells recovered.
        expected_rate <- min(0.076, 0.008 * (ncol(o) / 1000))
        nExp_poi <- round(expected_rate * ncol(o))
        nExp_poi_adj <- round(nExp_poi * (1 - homotypic_prop))

        pN_val <- 0.25
        pK_col_name <- paste0("pANN_", pN_val, "_", pK_opt, "_", nExp_poi)
        o <- doubletFinder(o, PCs = 1:30, pN = pN_val, pK = pK_opt, nExp = nExp_poi_adj, reuse.pANN = NULL, sct = FALSE)
        classification_col <- grep("^DF.classifications", colnames(o@meta.data), value = TRUE)
        classification_col <- classification_col[length(classification_col)]
        singlets <- colnames(o)[o@meta.data[[classification_col]] == "Singlet"]
        n_doublets_local <- ncol(o) - length(singlets)
        list(obj = subset(obj_qc, cells = singlets), n_doublets = n_doublets_local)
      },
      error = function(e) list(obj = NULL, n_doublets = NA_integer_, error = conditionMessage(e))
    )
    if (!is.null(result$obj)) {
      obj_final <- result$obj
      n_doublets <- result$n_doublets
      n_cells_final <- ncol(obj_final)
    } else {
      obj_final <- obj_qc
      doublet_error <- result$error
      message(sprintf("  DoubletFinder failed for %s: %s -- proceeding without doublet removal", sample_id, doublet_error))
    }
  } else {
    obj_final <- obj_qc
    message(sprintf("  %s: only %d cells post-QC (<50), skipping doublet detection", sample_id, n_cells_postqc))
  }

  saveRDS(obj_final, file.path(dir_seurat_qc, paste0(sample_id, ".rds")))

  qc_rows[[sample_id]] <- data.frame(
    sample_id = sample_id, GSM = row$GSM, patient = row$patient, condition = row$condition,
    n_cells_raw = n_cells_raw,
    n_fail_ncount_rna = n_fail_ncount, n_fail_nfeature_rna = n_fail_nfeature, n_fail_percent_mt = n_fail_mt,
    n_cells_post_qc = n_cells_postqc,
    n_doublets_removed = n_doublets,
    n_cells_final = n_cells_final,
    doublet_detection_error = doublet_error,
    stringsAsFactors = FALSE
  )
}

qc_summary <- do.call(rbind, qc_rows)
rownames(qc_summary) <- NULL

output_tsv <- file.path(repo_root, "results", "tables", "gse245601_qc_summary.tsv")
write.table(qc_summary, output_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
message(sprintf("Wrote %s (%d rows)", output_tsv, nrow(qc_summary)))

session_txt <- file.path(repo_root, "results", "tables", "gse245601_session_info.txt")
writeLines(capture.output(sessionInfo()), session_txt)
message(sprintf("Wrote %s", session_txt))
