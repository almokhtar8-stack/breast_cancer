#!/usr/bin/env Rscript
# GSE240112 -- extraction of the author-processed h5Seurat objects
# (docs/GSE240112_DATA_AUDIT.md section 4; docs/GSE240112_PREANALYSIS.md
# section G). SeuratDisk::LoadH5Seurat() fails hard on Seurat >=5.0.0
# (GetAssayData(slot=...) is defunct; SeuratDisk:::AssembleAssay() still
# calls it) -- this exactly matches the author GitHub README's own
# Seurat-v4-vs-v5 compatibility warning. This script bypasses
# LoadH5Seurat()/as.Seurat() entirely and reads the h5Seurat file's raw
# on-disk HDF5 structure directly (hdf5r has no Seurat-API dependency),
# reconstructing the raw-counts sparse matrix, cell metadata, and UMAP
# embedding as plain R objects. No cell, gene, or metadata value is
# dropped, subset, or recomputed by this script -- it is a pure format
# conversion, verified against the file's own declared dimensions.
#
# Input:  data/raw/gse240112/NTs.h5seurat
#         data/raw/gse240112/TTs_cancer_060223.h5seurat
# Output: data/processed/gse240112/tt_cancer_extracted.rds (not committed; data/processed is gitignored)
#         data/processed/gse240112/nt_extracted.rds (not committed)
#         data/processed/gse240112/tt_cancer_metadata.tsv (not committed)
#         data/processed/gse240112/tt_cancer_candidate_raw_counts.tsv (not committed)
#         data/processed/gse240112/tt_cancer_candidate_lognorm.tsv (not committed)
#         results/tables/gse240112/extraction_log.tsv (committed -- rows-in/rows-out provenance)

suppressPackageStartupMessages({
  library(hdf5r)
  library(Matrix)
  library(yaml)
})

args <- commandArgs(trailingOnly = FALSE)
repo_root <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", args, value = TRUE))), "..", ".."))
if (length(repo_root) == 0 || !dir.exists(file.path(repo_root, "config"))) repo_root <- getwd()

config <- yaml::read_yaml(file.path(repo_root, "config", "config.yaml"))
cfg <- config$gse240112

dir.create(file.path(repo_root, cfg$output$processed_dir), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(file.path(repo_root, cfg$output$tt_cancer_extraction_log_tsv)), recursive = TRUE, showWarnings = FALSE)

# --- generic h5Seurat-internal readers -------------------------------------

read_meta_field <- function(h5, field) {
  # An R factor is stored as an HDF5 group with "levels" (character) and
  # "values" (0-based integer codes); NA is encoded as a value of -1.
  # A plain vector (numeric/character) is stored as a bare dataset.
  obj <- h5[[paste0("meta.data/", field)]]
  if (is(obj, "H5Group")) {
    levels <- obj[["levels"]]$read()
    values <- obj[["values"]]$read()
    out <- ifelse(values < 0, NA_character_, levels[values + 1L])
    return(out)
  }
  obj$read()
}

read_counts_matrix <- function(h5, assay = "RNA") {
  grp <- h5[[paste0("assays/", assay, "/counts")]]
  dims <- h5attributes(grp)$dims # c(n_genes, n_cells), CSC sparse (genes are rows)
  data <- grp[["data"]]$read()
  indices <- grp[["indices"]]$read()
  indptr <- grp[["indptr"]]$read()
  mat <- new("dgCMatrix", i = as.integer(indices), p = as.integer(indptr), x = as.numeric(data), Dim = as.integer(dims))
  features <- h5[[paste0("assays/", assay, "/features")]]$read()
  cell_names <- h5[["cell.names"]]$read()
  rownames(mat) <- features
  colnames(mat) <- cell_names
  mat
}

read_umap <- function(h5) {
  emb <- h5[["reductions/umap/cell.embeddings"]]$read()
  # hdf5r reads row-major C storage as [n_cells, 2] already oriented correctly for h5Seurat's dense-matrix convention
  if (nrow(emb) != length(h5[["cell.names"]]$read())) emb <- t(emb)
  colnames(emb) <- c("umap_1", "umap_2")
  rownames(emb) <- h5[["cell.names"]]$read()
  emb
}

extract_object <- function(h5_path, meta_fields) {
  h5 <- H5File$new(h5_path, mode = "r")
  on.exit(h5$close_all())

  counts <- read_counts_matrix(h5, assay = "RNA")
  cell_names <- h5[["cell.names"]]$read()
  stopifnot(identical(colnames(counts), cell_names))
  stopifnot(!anyDuplicated(rownames(counts)))

  meta <- data.frame(cell_id = cell_names, stringsAsFactors = FALSE)
  for (f in meta_fields) {
    if (h5$link_exists(paste0("meta.data/", f))) {
      meta[[f]] <- read_meta_field(h5, f)
    }
  }

  umap <- tryCatch(as.data.frame(read_umap(h5)), error = function(e) NULL)
  if (!is.null(umap)) {
    umap$cell_id <- rownames(umap)
    meta <- merge(meta, umap, by = "cell_id", sort = FALSE)
  }
  stopifnot(nrow(meta) == ncol(counts))

  list(counts = counts, meta = meta, n_genes = nrow(counts), n_cells = ncol(counts))
}

extraction_log <- list()

# --- TT (cancer cells only): primary tumor-cell object ---------------------
message("=== extracting TT cancer-cell object ===")
tt_path <- file.path(repo_root, cfg$raw$tt_cancer_h5seurat)
tt_fields <- c(
  "orig.ident", "ident", "cell.annot", "seurat_clusters",
  "integrated_snn_res.0.5", "integrated_snn_res.0.7", "integrated_snn_res.0.9",
  "nCount_RNA", "nFeature_RNA", "nCount_SCT", "nFeature_SCT", "percent.mt",
  "scDblFinder.class", "scDblFinder.score",
  "AURKA_module1", "CASP3_module1", "ERBB2_module1", "ESR1_module1",
  "PLAU_module1", "STAT1_module1", "TR_module1", "VEGF_module1"
)
tt <- extract_object(tt_path, tt_fields)
stopifnot(tt$n_cells == 9942, tt$n_genes == 27161) # matches h5 file's own declared dims (data audit section 4)
tt$meta$group <- substr(tt$meta$orig.ident, 1, 2)
stopifnot(all(tt$meta$group %in% c("PT", "RT")))
saveRDS(tt, file.path(repo_root, cfg$output$tt_cancer_rds))
write.table(tt$meta, file.path(repo_root, cfg$output$tt_cancer_metadata_tsv), sep = "\t", quote = FALSE, row.names = FALSE)
message(sprintf("  TT: %d genes x %d cells extracted (0 lost; matches source h5 dims)", tt$n_genes, tt$n_cells))
extraction_log[["tt_cancer"]] <- data.frame(
  object = "tt_cancer", source_file = cfg$raw$tt_cancer_h5seurat,
  n_genes_in_source = 27161, n_cells_in_source = 9942,
  n_genes_extracted = tt$n_genes, n_cells_extracted = tt$n_cells, n_cells_lost = 9942 - tt$n_cells,
  stringsAsFactors = FALSE
)

# --- Candidate gene raw + log-normalized per-cell expression (Phase 8/12) --
candidates <- c(cfg$candidates$thirteen, cfg$candidates$paics)
present <- candidates[candidates %in% rownames(tt$counts)]
absent <- setdiff(candidates, present)
message(sprintf("  candidate genes present in TT feature space: %d/%d (absent: %s)", length(present), length(candidates), paste(absent, collapse = ", ")))

raw_sub <- as.matrix(tt$counts[present, , drop = FALSE])
raw_df <- data.frame(cell_id = colnames(raw_sub), t(raw_sub), check.names = FALSE, stringsAsFactors = FALSE)
write.table(raw_df, file.path(repo_root, cfg$output$tt_cancer_candidate_raw_tsv), sep = "\t", quote = FALSE, row.names = FALSE)

# log-normalized (per-cell library-size normalization, natural-log(1 + counts-per-10k)) -- Seurat's standard "LogNormalize", computed directly from the raw counts we already extracted, not read from the SCT assay (whose "data" slot uses a different, non-count-depth-comparable Pearson-residual-based correction).
lib_size <- Matrix::colSums(tt$counts)
norm_sub <- sweep(raw_sub, 2, lib_size, "/") * 1e4
lognorm_sub <- log1p(norm_sub)
lognorm_df <- data.frame(cell_id = colnames(lognorm_sub), t(lognorm_sub), check.names = FALSE, stringsAsFactors = FALSE)
write.table(lognorm_df, file.path(repo_root, cfg$output$tt_cancer_candidate_lognorm_tsv), sep = "\t", quote = FALSE, row.names = FALSE)
message(sprintf("  wrote candidate raw + log-normalized per-cell expression for %d cells", ncol(raw_sub)))

# --- NT (normal tissue): optional Phase 15 context --------------------------
message("=== extracting NT object ===")
nt_path <- file.path(repo_root, cfg$raw$nt_h5seurat)
nt_fields <- c("orig.ident", "ident", "SCT_snn_res.0.5", "SCT_snn_res.0.7", "nCount_RNA", "nFeature_RNA", "percent.mt", "scDblFinder.class", "scDblFinder.score")
nt <- extract_object(nt_path, nt_fields)
stopifnot(nt$n_cells == 7529, nt$n_genes == 25543)
saveRDS(nt, file.path(repo_root, cfg$output$nt_rds))
message(sprintf("  NT: %d genes x %d cells extracted (0 lost; matches source h5 dims)", nt$n_genes, nt$n_cells))
extraction_log[["nt"]] <- data.frame(
  object = "nt", source_file = cfg$raw$nt_h5seurat,
  n_genes_in_source = 25543, n_cells_in_source = 7529,
  n_genes_extracted = nt$n_genes, n_cells_extracted = nt$n_cells, n_cells_lost = 7529 - nt$n_cells,
  stringsAsFactors = FALSE
)

log_out <- do.call(rbind, extraction_log)
write.table(log_out, file.path(repo_root, cfg$output$tt_cancer_extraction_log_tsv), sep = "\t", quote = FALSE, row.names = FALSE)
message("Done. Extraction log written to results/tables/gse240112/extraction_log.tsv")
