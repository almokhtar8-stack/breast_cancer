# Project Status

Last updated: 2026-08-11. See [`../README.md`](../README.md) for the full
narrative, [`CODE_MAP.md`](CODE_MAP.md) for where each piece of code
lives, and [`gse245601_PREANALYSIS.md`](gse245601_PREANALYSIS.md) for the
single-cell pre-analysis plan and blinding rule referenced below.

## COMPLETE

- CRISPR reanalysis of Hany et al. 2023 (19,103 genes fitted, 28 Gate-1
  hits at FDR<0.1)
- GSE118713 bulk RNA analysis (QC, limma DE across MCF7/TAMR/FASR)
- CRISPR x bulk-RNA integration and evidence-class hierarchy for the 13
  sensitising candidates
- PAICS published-benchmark check (not a Gate-1 hit in this reanalysis)
- NEBULA poster figures (`src/nebula_plots_final.py`)
- GSE245601 download, checksum-verified manifest, candidate feature-space
  check
- GSE245601 preprocessing: Seurat import, QC, doublet removal,
  normalization, clustering, UMAP
- GSE245601 broad cell-type annotation and epithelial identification
  (candidate-gene-blind)
- GSE245601 InferCNV primary malignant-cell classification
- GSE245601 CopyKAT independent sensitivity check
- GSE245601 per-cell metadata frozen (44,140 cells); pseudobulk
  pair-eligibility flags computed per tumor (pseudobulk aggregation itself
  has not been run)

## CURRENT

- Deciding the statistical design for the GSE245601 candidate-level
  analysis, given that only 3/10 paired tumors (Tumor_02, Tumor_03,
  Tumor_07) currently meet the frozen ≥50-malignant-cell-per-arm
  eligibility rule, and InferCNV/CopyKAT concordance is highly variable
  between samples (mean ~56%, range ~0-100%).

## NEXT

- GSE245601 pseudobulk aggregation and candidate-level differential
  expression (blocked on the statistical-design decision above)
- Single-cell candidate ranking

## LATER

- GSE111151 independent confirmation
- Mechanism / pathway analysis
- Druggability assessment
- Normal-tissue expression analysis
- Final candidate prioritization

No candidate-gene expression value has been inspected in GSE245601 at any
point — the blinding rule in
[`gse245601_PREANALYSIS.md`](gse245601_PREANALYSIS.md) remains
intact through the end of the CURRENT preprocessing phase.
