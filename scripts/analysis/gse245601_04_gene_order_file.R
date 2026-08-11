#!/usr/bin/env Rscript
# GSE245601 -- InferCNV prerequisite: build the gene-order (genomic
# position) file required by infercnv::CreateInfercnvObject().
#
# Built from TxDb.Hsapiens.UCSC.hg38.knownGene + org.Hs.eg.db (a standard,
# documented way of constructing an infercnv gene_order_file when one is
# not separately downloaded) -- restricted to genes present in the merged
# Seurat object's expression matrix, sorted by chromosome and start
# position, canonical chromosomes only (chr1-22, X, Y). This file contains
# gene symbol + genomic coordinates only; no expression values.
#
# Output: data/processed/gse245601/infercnv/gene_order.tsv (gitignored;
#         regenerable from this script)

suppressPackageStartupMessages({
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(org.Hs.eg.db)
  library(GenomicFeatures)
  library(Seurat)
})

args <- commandArgs(trailingOnly = FALSE)
repo_root <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", args, value = TRUE))), "..", ".."))
if (length(repo_root) == 0 || !dir.exists(file.path(repo_root, "config"))) repo_root <- getwd()

merged <- readRDS(file.path(repo_root, "data", "processed", "gse245601", "seurat_clustered", "annotated.rds"))
expr_genes <- rownames(merged)
message("Genes in expression matrix: ", length(expr_genes))

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
genes_gr <- genes(txdb)

entrez_to_symbol <- AnnotationDbi::select(org.Hs.eg.db, keys = names(genes_gr), columns = "SYMBOL", keytype = "ENTREZID")
entrez_to_symbol <- entrez_to_symbol[!duplicated(entrez_to_symbol$ENTREZID), ]
rownames(entrez_to_symbol) <- entrez_to_symbol$ENTREZID

genes_gr$symbol <- entrez_to_symbol[names(genes_gr), "SYMBOL"]
genes_gr <- genes_gr[!is.na(genes_gr$symbol)]

canonical_chr <- paste0("chr", c(1:22, "X", "Y"))
genes_gr <- genes_gr[as.character(GenomeInfoDb::seqnames(genes_gr)) %in% canonical_chr]

# One row per gene symbol present in our expression matrix. If a symbol
# maps to >1 Entrez/genomic locus (rare), keep the first (alphabetically
# first ENTREZID) for a single, deterministic gene-order file -- flagged
# here rather than silently averaged/merged.
genes_gr <- genes_gr[genes_gr$symbol %in% expr_genes]
dup_symbols <- unique(genes_gr$symbol[duplicated(genes_gr$symbol)])
if (length(dup_symbols) > 0) {
  message(sprintf("Note: %d gene symbols map to >1 genomic locus in this TxDb; keeping first locus for each.", length(dup_symbols)))
  genes_gr <- genes_gr[!duplicated(genes_gr$symbol)]
}

gene_order <- data.frame(
  gene = genes_gr$symbol,
  chr = as.character(GenomeInfoDb::seqnames(genes_gr)),
  start = BiocGenerics::start(genes_gr),
  end = BiocGenerics::end(genes_gr),
  stringsAsFactors = FALSE
)

chr_order <- factor(gene_order$chr, levels = canonical_chr)
gene_order <- gene_order[order(chr_order, gene_order$start), ]
gene_order <- gene_order[!duplicated(gene_order$gene), ]

dir_infercnv <- file.path(repo_root, "data", "processed", "gse245601", "infercnv")
dir.create(dir_infercnv, recursive = TRUE, showWarnings = FALSE)
out_path <- file.path(dir_infercnv, "gene_order.tsv")
write.table(gene_order, out_path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

message(sprintf("Wrote %s: %d genes (of %d in expression matrix; %.1f%% coverage)",
                 out_path, nrow(gene_order), length(expr_genes), 100 * nrow(gene_order) / length(expr_genes)))

# Candidate genes + PAICS must be present in the gene-order file too, else
# they would be silently dropped from the CNV-input matrix (they are NOT
# used for CNV inference or malignancy calling -- inferCNV requires a
# fixed, cohort-wide gene-ordered matrix, so ALL genes including the 13
# candidates naturally pass through this step; this check only confirms
# they were not accidentally lost by the coordinate lookup).
candidate_genes <- c("USP34", "CTDNEP1", "EIF4ENIF1", "HMGB1", "KDM1A", "PET117", "TADA2B", "VEZF1", "ICK", "SUPT4H1", "TLK2", "TSR3", "USP17L29", "PAICS")
missing_candidates <- setdiff(candidate_genes, gene_order$gene)
if (length(missing_candidates) > 0) {
  message("Note: the following candidate genes have no chr/start/end in this TxDb build and will be absent from the CNV gene-order matrix (informational only -- CNV inference itself never selects on or reports these genes): ", paste(missing_candidates, collapse = ", "))
}
