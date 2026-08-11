#!/usr/bin/env Rscript
# GSE245601 -- Step 11: required QC/annotation figures. NO candidate gene
# (13 candidates or PAICS) is plotted anywhere in this script -- every
# figure here is a structural/QC/annotation figure only.
#
# Output: results/figures/gse245601_preprocessing/*.pdf (committed)

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
})

args <- commandArgs(trailingOnly = FALSE)
repo_root <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", args, value = TRUE))), "..", ".."))
if (length(repo_root) == 0 || !dir.exists(file.path(repo_root, "config"))) repo_root <- getwd()

dir_fig <- file.path(repo_root, "results", "figures", "gse245601_preprocessing")
dir.create(dir_fig, recursive = TRUE, showWarnings = FALSE)

merged <- readRDS(file.path(repo_root, "data", "processed", "gse245601", "seurat_clustered", "annotated.rds"))
qc_summary <- read.delim(file.path(repo_root, "results", "tables", "gse245601_qc_summary.tsv"), stringsAsFactors = FALSE)
malignant_labels <- read.delim(file.path(repo_root, "results", "tables", "gse245601_malignant_cell_labels.tsv"), stringsAsFactors = FALSE)
counts_per_sample <- read.delim(file.path(repo_root, "results", "tables", "gse245601_malignant_counts_per_sample.tsv"), stringsAsFactors = FALSE)
concordance_path <- file.path(repo_root, "results", "tables", "gse245601_malignancy_concordance.tsv")

theme_set(theme_minimal(base_size = 10))

# 1. QC distributions (post-QC, since raw per-cell values before filtering
# are only available in the transient raw object; post-QC distributions by
# sample are shown here alongside the QC summary table's before/after counts).
p1 <- VlnPlot(merged, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), group.by = "sample_id", pt.size = 0, ncol = 1) &
  theme(axis.text.x = element_text(angle = 90, size = 6))
ggsave(file.path(dir_fig, "01_qc_distributions_post_filter.pdf"), p1, width = 12, height = 10)

qc_long <- data.frame(
  sample_id = rep(qc_summary$sample_id, 3),
  stage = rep(c("raw", "post_qc", "final_post_doublet"), each = nrow(qc_summary)),
  n_cells = c(qc_summary$n_cells_raw, qc_summary$n_cells_post_qc, qc_summary$n_cells_final)
)
qc_long$stage <- factor(qc_long$stage, levels = c("raw", "post_qc", "final_post_doublet"))
p1b <- ggplot(qc_long, aes(x = sample_id, y = n_cells, fill = stage)) +
  geom_col(position = "dodge") + coord_flip() + labs(title = "Cells before/after each QC stage", x = NULL, y = "n cells")
ggsave(file.path(dir_fig, "01b_qc_cells_before_after.pdf"), p1b, width = 8, height = 8)

# 2. UMAP by broad cell type
p2 <- DimPlot(merged, group.by = "broad_cell_type", label = TRUE, raster = FALSE) + labs(title = "UMAP by broad cell type")
ggsave(file.path(dir_fig, "02_umap_broad_cell_type.pdf"), p2, width = 8, height = 6)

# 3. UMAP by patient
p3 <- DimPlot(merged, group.by = "patient", raster = FALSE) + labs(title = "UMAP by patient")
ggsave(file.path(dir_fig, "03_umap_by_patient.pdf"), p3, width = 8, height = 6)

# 4. UMAP by treatment
p4 <- DimPlot(merged, group.by = "condition", raster = FALSE, cols = c(Control = "#94A3B8", Tamoxifen = "#6E44FF")) + labs(title = "UMAP by treatment condition")
ggsave(file.path(dir_fig, "04_umap_by_treatment.pdf"), p4, width = 7, height = 6)

# 5. Epithelial-cell UMAP (highlighted)
merged$is_epithelial_label <- ifelse(merged$broad_cell_type == "epithelial", "epithelial", "other")
p5 <- DimPlot(merged, group.by = "is_epithelial_label", cols = c(epithelial = "#F2650B", other = "#D1D5DB"), raster = FALSE) +
  labs(title = "Epithelial vs. non-epithelial cells")
ggsave(file.path(dir_fig, "05_umap_epithelial.pdf"), p5, width = 7, height = 6)

# 6. InferCNV / CNV diagnostic figure: CNV score distribution by primary label
# (local copy + local column, so the shared malignant_labels$primary_malignancy_label
# character column used later by Figure 7 is never mutated into a factor)
malignant_labels_fig6 <- malignant_labels
malignant_labels_fig6$primary_malignancy_label <- factor(malignant_labels_fig6$primary_malignancy_label, levels = c("non-malignant epithelial", "malignant"))
p6 <- ggplot(malignant_labels_fig6, aes(x = cnv_score, y = cnv_correlation_to_seed, color = primary_malignancy_label)) +
  geom_point(size = 0.6, alpha = 0.6) +
  scale_color_manual(values = c("non-malignant epithelial" = "#94A3B8", "malignant" = "#6E44FF")) +
  labs(title = "InferCNV diagnostic: CNV score vs. correlation-to-seed", x = "CNV score (mean squared deviation)", y = "Kendall correlation to malignant seed") +
  facet_wrap(~sample_id, scales = "free")
ggsave(file.path(dir_fig, "06_infercnv_diagnostic_scatter.pdf"), p6, width = 16, height = 12)

# 7. UMAP malignant vs non-malignant epithelial (epithelial cells only, in full UMAP context)
# Built directly from the UMAP embedding + explicit geom_point/scale_color_manual
# (bypassing DimPlot's internal grouping/color machinery) for reliable color mapping.
malignancy_umap_label <- setNames(rep("non-epithelial", ncol(merged)), colnames(merged))
malignancy_umap_label[malignant_labels$cell_id] <- as.character(malignant_labels$primary_malignancy_label)
umap_df <- as.data.frame(Embeddings(merged, "umap"))
umap_df$malignancy_umap_label <- factor(malignancy_umap_label[rownames(umap_df)], levels = c("non-epithelial", "non-malignant epithelial", "malignant"))
stopifnot(sum(is.na(umap_df$malignancy_umap_label)) == 0)
stopifnot(sum(umap_df$malignancy_umap_label == "malignant") == sum(malignant_labels$primary_malignancy_label == "malignant"))
umap_df <- umap_df[order(umap_df$malignancy_umap_label), ]  # plot non-epithelial first so malignant points are drawn on top
p7 <- ggplot(umap_df, aes(x = umap_1, y = umap_2, color = malignancy_umap_label)) +
  geom_point(size = 0.4) +
  scale_color_manual(values = c("non-epithelial" = "#E5E7EB", "non-malignant epithelial" = "#0F9B8E", "malignant" = "#6E44FF")) +
  labs(title = "Malignant vs. non-malignant epithelial cells", color = NULL) +
  guides(color = guide_legend(override.aes = list(size = 3)))
ggsave(file.path(dir_fig, "07_umap_malignant_vs_nonmalignant.pdf"), p7, width = 8, height = 6)

# 8. Malignant-cell counts by patient x condition
counts_per_sample$condition <- factor(counts_per_sample$condition, levels = c("Control", "Tamoxifen"))
p8 <- ggplot(counts_per_sample, aes(x = patient, y = n_primary_malignant, fill = condition)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = c(Control = "#94A3B8", Tamoxifen = "#6E44FF")) +
  geom_hline(yintercept = 50, linetype = "dashed", color = "#D946EF") +
  labs(title = "Primary (InferCNV) malignant-cell counts by patient x condition", subtitle = "Dashed line = frozen eligibility threshold (>=50 cells)", y = "n malignant cells")
ggsave(file.path(dir_fig, "08_malignant_counts_by_patient_condition.pdf"), p8, width = 9, height = 5)

# 9. Primary-vs-sensitivity malignancy concordance
if (file.exists(concordance_path)) {
  concordance <- read.delim(concordance_path, stringsAsFactors = FALSE)
  conc_ok <- concordance[concordance$status == "ok", ]
  if (nrow(conc_ok) > 0) {
    p9 <- ggplot(conc_ok, aes(x = sample_id, y = concordance_rate)) +
      geom_col(fill = "#0F9B8E") + coord_flip() + ylim(0, 1) +
      labs(title = "Primary (InferCNV) vs. sensitivity (CopyKAT) malignancy concordance", y = "concordance rate", x = NULL)
    ggsave(file.path(dir_fig, "09_malignancy_concordance.pdf"), p9, width = 8, height = 6)
  }
}

message("Wrote QC/annotation figures to ", dir_fig)
