#!/usr/bin/env Rscript
# GSE240112 -- Phase 6/14 secondary analysis: broad cell-type compartment
# calling and all-epithelial-cell pseudobulk, from the raw GEO Cell Ranger
# matrices (docs/GSE240112_DATA_AUDIT.md section 4 point 2). This exists
# only because no author-provided all-cell-type PT/RT object exists (the
# only PT/RT-covering processed object, TTs_cancer_060223.h5seurat, is
# restricted to cancer cells by construction). This is explicitly a
# SECONDARY/sensitivity analysis (docs/GSE240112_PREANALYSIS.md section G)
# -- it identifies broad EPITHELIAL vs non-epithelial compartments via
# standard canonical markers only, and never redefines or replaces the
# Phase 7 Case A tumor-cell labels used for the primary analysis.
#
# GSM-to-sample mapping is by GEO's declared replicate order (the only
# unambiguous source; see data audit section 3) -- raw internal file
# prefixes (PT1/PT2/PT5, RT3/RT4/RT6) differ from the clean PT1-3/RT1-3
# labels used throughout this analysis.
#
# QC: standard fixed thresholds (min.cells=3, min.features=200 at load,
# consistent with the paper's own stated CreateSeuratObject call; then
# nFeature_RNA>200 & percent.mt<20 & doublet removal via scDblFinder) --
# not the paper's per-sample adaptive 3xMAD procedure, which is
# unnecessary complexity for a compartment-calling-only secondary
# analysis. No sample is merged/integrated across the batch beyond simple
# concatenation -- Phase 6 explicitly discourages integrating away real
# PT/RT biology, and none is scientifically justified here since the goal
# is only a broad epithelial/immune/stromal/endothelial split, not
# fine subclustering.
#
# Input:  data/raw/gse240112/cellranger/*.tsv.gz, *.mtx.gz
# Output: results/tables/gse240112/epithelial_compartment_composition.tsv (committed)
#         results/tables/gse240112_pseudobulk/epithelial_counts.tsv.gz (committed)
#         results/tables/gse240112_pseudobulk/epithelial_metadata.tsv (committed)
#         results/tables/gse240112_pseudobulk/epithelial_genomewide_de.tsv.gz (committed)
#         results/tables/gse240112_pseudobulk/epithelial_edger_filtering_summary.tsv (committed)
#         results/figures/gse240112_pseudobulk/broad_celltype_umap.png (committed)

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(scDblFinder)
  library(SingleCellExperiment)
  library(edgeR)
  library(yaml)
  library(data.table)
})

args <- commandArgs(trailingOnly = FALSE)
repo_root <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", args, value = TRUE))), "..", ".."))
if (length(repo_root) == 0 || !dir.exists(file.path(repo_root, "config"))) repo_root <- getwd()

config <- yaml::read_yaml(file.path(repo_root, "config", "config.yaml"))
cfg <- config$gse240112
cr_dir <- file.path(repo_root, cfg$raw$cellranger_dir)

set.seed(1)

# GEO declared replicate order -> raw internal file prefix (data audit section 3)
sample_map <- list(
  PT1 = list(gsm = "GSM7681687", prefix = "PT1"),
  PT2 = list(gsm = "GSM7681688", prefix = "PT2"),
  PT3 = list(gsm = "GSM7681689", prefix = "PT5"),
  RT1 = list(gsm = "GSM7681690", prefix = "RT3"),
  RT2 = list(gsm = "GSM7681691", prefix = "RT4"),
  RT3 = list(gsm = "GSM7681692", prefix = "RT6")
)

message("=== loading raw Cell Ranger matrices ===")
objs <- list()
for (sample_id in names(sample_map)) {
  m <- sample_map[[sample_id]]
  prefix <- file.path(cr_dir, paste0(m$gsm, "_", m$prefix, "_"))
  # Read10X() expects a directory of standard-named files; ours are
  # GSM-prefixed, so the barcodes/features/matrix triplet is read directly.
  barcodes <- read.delim(paste0(prefix, "barcodes.tsv.gz"), header = FALSE, stringsAsFactors = FALSE)$V1
  features <- read.delim(paste0(prefix, "features.tsv.gz"), header = FALSE, stringsAsFactors = FALSE)
  mat <- Matrix::readMM(gzfile(paste0(prefix, "matrix.mtx.gz")))
  rownames(mat) <- make.unique(features$V2)
  colnames(mat) <- paste0(sample_id, "_", barcodes)
  n_cells_raw <- ncol(mat)

  obj <- CreateSeuratObject(counts = mat, project = sample_id, min.cells = 3, min.features = 200)
  n_cells_min_filter <- ncol(obj)
  obj$sample_id <- sample_id
  obj$group <- substr(sample_id, 1, 2)
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")

  qc_keep <- obj$nFeature_RNA > 200 & obj$percent.mt < 20
  obj <- obj[, qc_keep]
  n_cells_qc <- ncol(obj)

  sce <- as.SingleCellExperiment(obj)
  sce <- scDblFinder(sce, verbose = FALSE)
  obj$scDblFinder.class <- sce$scDblFinder.class
  obj <- obj[, obj$scDblFinder.class == "singlet"]
  n_cells_final <- ncol(obj)

  message(sprintf("  %s (%s_%s): raw=%d -> min.cells/min.features=%d -> QC=%d -> singlets=%d",
                   sample_id, m$gsm, m$prefix, n_cells_raw, n_cells_min_filter, n_cells_qc, n_cells_final))
  objs[[sample_id]] <- obj
}

message("=== merging (no cross-sample integration) ===")
merged <- merge(objs[[1]], y = objs[-1], add.cell.ids = NULL)
merged <- JoinLayers(merged)
message(sprintf("merged: %d genes x %d cells", nrow(merged), ncol(merged)))

message("=== normalize / cluster / UMAP (standard, no integration) ===")
merged <- NormalizeData(merged, verbose = FALSE)
merged <- FindVariableFeatures(merged, verbose = FALSE)
merged <- ScaleData(merged, verbose = FALSE)
merged <- RunPCA(merged, npcs = 30, verbose = FALSE)
merged <- FindNeighbors(merged, dims = 1:30, verbose = FALSE)
merged <- FindClusters(merged, resolution = 0.5, verbose = FALSE)
merged <- RunUMAP(merged, dims = 1:30, verbose = FALSE)

checkpoint_path <- file.path(repo_root, cfg$output$processed_dir, "epithelial_pipeline_checkpoint.rds")
dir.create(dirname(checkpoint_path), recursive = TRUE, showWarnings = FALSE)
saveRDS(merged, checkpoint_path)
message(sprintf("wrote checkpoint %s (post-clustering, pre-compartment-calling)", checkpoint_path))

message("=== broad compartment marker scoring ===")
# Direct, transparent per-cluster mean log-normalized expression per
# marker set (not Seurat::AddModuleScore, whose background-correction and
# automatic column-naming produced an unverifiable, and in a first run
# incorrect -- every cluster called "epithelial" -- classification).
marker_sets <- list(
  epithelial = c("EPCAM", "KRT8", "KRT18", "KRT19"),
  immune = c("PTPRC"),
  stromal = c("PDGFRB", "COL1A1", "DCN"),
  endothelial = c("PECAM1", "VWF")
)
log_norm <- GetAssayData(merged, assay = "RNA", layer = "data")
cell_scores <- sapply(names(marker_sets), function(nm) {
  genes_present <- intersect(marker_sets[[nm]], rownames(log_norm))
  stopifnot(length(genes_present) > 0)
  Matrix::colMeans(log_norm[genes_present, , drop = FALSE])
})
rownames(cell_scores) <- colnames(merged)
stopifnot(identical(colnames(cell_scores), names(marker_sets)))

cluster_mean_scores <- aggregate(cell_scores, by = list(cluster = merged$seurat_clusters), FUN = mean)
message("per-cluster mean marker-set scores:")
print(cluster_mean_scores)
score_matrix <- as.matrix(cluster_mean_scores[, names(marker_sets)])
cluster_call <- names(marker_sets)[apply(score_matrix, 1, which.max)]
names(cluster_call) <- cluster_mean_scores$cluster
merged$broad_cell_type <- unname(cluster_call[as.character(merged$seurat_clusters)])
stopifnot(!anyNA(merged$broad_cell_type))

message("cluster -> broad_cell_type calls:")
print(table(merged$seurat_clusters, merged$broad_cell_type))

composition <- as.data.frame(table(merged$sample_id, merged$broad_cell_type))
colnames(composition) <- c("sample_id", "broad_cell_type", "n_cells")
composition$group <- substr(as.character(composition$sample_id), 1, 2)
comp_path <- file.path(repo_root, cfg$output$epithelial$compartment_composition_tsv)
dir.create(dirname(comp_path), recursive = TRUE, showWarnings = FALSE)
write.table(composition, comp_path, sep = "\t", quote = FALSE, row.names = FALSE)
message(sprintf("wrote %s", comp_path))

umap_df <- as.data.frame(Embeddings(merged, "umap"))
umap_df$broad_cell_type <- merged$broad_cell_type
umap_df$sample_id <- merged$sample_id
fig_dir <- file.path(repo_root, cfg$output$qc$figures_dir)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
png(file.path(fig_dir, "broad_celltype_umap.png"), width = 1400, height = 1100, res = 150)
plot(umap_df$umap_1, umap_df$umap_2, col = as.factor(umap_df$broad_cell_type), pch = 16, cex = 0.3,
     xlab = "UMAP_1", ylab = "UMAP_2", main = "GSE240112 PT/RT raw Cell Ranger: broad compartment (marker-based)")
legend("topright", legend = levels(as.factor(umap_df$broad_cell_type)), col = 1:length(levels(as.factor(umap_df$broad_cell_type))), pch = 16, cex = 0.8)
dev.off()
message("wrote broad_celltype_umap.png")

message("=== epithelial pseudobulk ===")
epi_cells <- colnames(merged)[merged$broad_cell_type == "epithelial"]
epi_counts <- GetAssayData(merged, assay = "RNA", layer = "counts")[, epi_cells, drop = FALSE]
epi_meta <- merged@meta.data[epi_cells, c("sample_id", "group")]

sample_groups <- split(rownames(epi_meta), epi_meta$sample_id)
stopifnot(setequal(names(sample_groups), names(sample_map)))
pb_mat <- sapply(sample_groups, function(cell_ids) Matrix::rowSums(epi_counts[, cell_ids, drop = FALSE]))
rownames(pb_mat) <- rownames(epi_counts)
pb_mat <- pb_mat[, names(sample_map)]

pb_meta <- data.frame(
  sample_id = colnames(pb_mat),
  group = substr(colnames(pb_mat), 1, 2),
  n_contributing_cells = lengths(sample_groups)[colnames(pb_mat)],
  total_library_size = colSums(pb_mat),
  n_detected_genes = colSums(pb_mat > 0),
  stringsAsFactors = FALSE
)

epi_counts_path <- file.path(repo_root, cfg$output$epithelial$counts_tsv)
epi_meta_path <- file.path(repo_root, cfg$output$epithelial$metadata_tsv)
dir.create(dirname(epi_counts_path), recursive = TRUE, showWarnings = FALSE)
fwrite(data.table(gene = rownames(pb_mat), as.data.frame(pb_mat)), epi_counts_path, sep = "\t")
write.table(pb_meta, epi_meta_path, sep = "\t", quote = FALSE, row.names = FALSE)
message(sprintf("wrote epithelial pseudobulk: %d genes x %d samples", nrow(pb_mat), ncol(pb_mat)))
print(pb_meta)

message("=== epithelial edgeR DE (RT vs PT, unpaired, same design as primary) ===")
pb_meta$group <- factor(pb_meta$group, levels = c("PT", "RT"))
design <- model.matrix(~group, data = pb_meta)
dge <- DGEList(counts = pb_mat, samples = pb_meta)
n_before <- nrow(dge)
keep <- filterByExpr(dge, design = design)
n_after <- sum(keep)
dge <- dge[keep, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge, method = "TMM")
dge <- estimateDisp(dge, design)
fit <- glmQLFit(dge, design)
qlf <- glmQLFTest(fit, coef = "groupRT")
results <- topTags(qlf, n = Inf, sort.by = "none")$table
results$gene <- rownames(results)
avg_log_cpm <- aveLogCPM(dge)
results$avg_log_cpm <- avg_log_cpm[match(results$gene, rownames(dge))]
out <- results[, c("gene", "logFC", "avg_log_cpm", "PValue", "FDR")]
colnames(out) <- c("gene", "log2fc", "avg_log_cpm", "p_value", "fdr")
out <- out[order(out$gene), ]

de_path <- file.path(repo_root, cfg$output$epithelial$genomewide_de_tsv)
summary_path <- file.path(repo_root, cfg$output$epithelial$edger_summary_tsv)
fwrite(out, de_path, sep = "\t")
filtering_summary <- data.frame(
  label = "epithelial_RT_vs_PT", n_samples = ncol(dge),
  n_genes_before_filterByExpr = n_before, n_genes_after_filterByExpr = n_after,
  normalization_method = "TMM", de_method = "edgeR glmQLFit/glmQLFTest",
  design_formula = "~ group", contrast = "RT vs PT (RT - PT), unpaired",
  residual_df = fit$df.residual[1], stringsAsFactors = FALSE
)
write.table(filtering_summary, summary_path, sep = "\t", quote = FALSE, row.names = FALSE)
message(sprintf("wrote %s (%d genes tested); filtering %d -> %d", de_path, nrow(out), n_before, n_after))
message("Done.")
