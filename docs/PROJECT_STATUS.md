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

## COMPLETE (continued)

- GSE245601 InferCNV/CopyKAT method audit and CNV-score-metric diagnostic
  (`docs/CNV_METHOD_AUDIT.md`,
  `docs/GSE245601_INFERCNV_THRESHOLD_DIAGNOSTIC.md`,
  `docs/GSE245601_INFERCNV_SCORE_METRIC_DIAGNOSTIC.md`): established that
  InferCNV/CopyKAT disagreement is real and sample-dependent, and that the
  InferCNV downstream classifier's whole-genome CNV score is
  extent-weighted with an effectively fixed 0.01 floor -- a known,
  documented limitation, not a coding error. CNV-method optimization is
  stopped; the frozen InferCNV labels remain the primary malignant-cell
  definition (`gse245601_PREANALYSIS.md` section 13).
- Two-track pseudobulk design frozen (section 13): Track A
  (epithelial-compartment, all 10 paired tumors) and Track B (strict
  malignant, Tumor_02/03/07 only, n=3, exploratory) -- decided before any
  pseudobulk aggregation or candidate expression was inspected.
- GSE245601 pseudobulk aggregation, QC, edgeR differential expression
  (Track A and Track B), 13-candidate + PAICS extraction, patient-level
  visualization, malignant-vs-non-malignant epithelial context check, and
  integration with the frozen CRISPR/bulk-RNA candidate evidence.

## CURRENT

- Reviewing the integrated single-cell + CRISPR + bulk-RNA candidate
  evidence for next steps (see the GSE245601 pseudobulk candidate
  integration report for current findings).

## NEXT

- GSE111151 independent confirmation, pathway-level analysis,
  druggability, normal-tissue expression context, and poster redesign
  (explicitly out of scope for the pseudobulk phase just completed).
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
