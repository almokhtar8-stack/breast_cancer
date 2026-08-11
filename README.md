# Breast Cancer Tamoxifen Response Project

Computational reanalysis identifying and prioritizing genes whose loss
modulates tamoxifen response in ER-positive breast cancer, by combining a
public genome-wide CRISPR screen with public bulk and single-cell
transcriptomics.

**Quick links:** [Code map](docs/CODE_MAP.md) ·
[Project status](docs/PROJECT_STATUS.md) ·
[Root pre-analysis plan](PREANALYSIS.md) ·
[GSE245601 pre-analysis plan](docs/gse245601_PREANALYSIS.md) ·
[Result tables](results/tables/) ·
[Figures](results/figures/) ·
[Tests](tests/)

## Overview

This project computationally reanalyses published functional and
transcriptomic datasets to identify and prioritize genes whose inhibition
may modulate tamoxifen response in ER+ breast cancer. It is a
**computational reanalysis of public data — no wet-lab work was performed
by this project.** A public genome-wide CRISPR knockout screen (Hany et
al. 2023) supplies functional perturbation evidence; public bulk RNA-seq
of tamoxifen-resistant cell lines (GSE118713) supplies resistance-associated
expression evidence; public human single-cell RNA-seq of tumors treated ex
vivo with tamoxifen (GSE245601) supplies additional human-tumor response
context. No claim of therapeutic efficacy or safety is made or implied by
any result in this repository.

## Research Workflow

```
Genome-wide CRISPR screen (Hany et al. 2023)              [COMPLETE]
        |
Bulk RNA resistance analysis (GSE118713)                   [COMPLETE]
        |
CRISPR x RNA candidate prioritization                       [COMPLETE]
        |
Human single-cell response context (GSE245601)
   preprocessing + malignant-cell ID                        [COMPLETE]
   candidate-level expression / pseudobulk DE / ranking      [NOT STARTED]
        |
Independent validation (GSE111151)                          [NOT STARTED]
        |
Pathway / mechanism                                         [NOT STARTED]
        |
Druggability                                                 [NOT STARTED]
        |
Normal-tissue expression                                     [NOT STARTED]
        |
Final candidate prioritization                               [NOT STARTED]
```

**Complete:**
- CRISPR reanalysis of the Hany et al. 2023 screen
- GSE118713 bulk RNA analysis
- CRISPR x bulk-RNA integration and candidate evidence classification
- PAICS published-benchmark check
- NEBULA poster figures
- GSE245601 single-cell preprocessing (QC, doublet removal, normalization, clustering)
- Broad cell-type annotation and epithelial identification
- InferCNV primary malignant-cell classification
- CopyKAT independent sensitivity analysis
- Malignant-cell labels frozen

**Current decision point:** statistical design for the GSE245601
candidate-level analysis (only 3/10 paired tumors currently meet the
frozen pseudobulk-eligibility threshold — see
[Limitations](#limitations--current-methodological-question)).

**Not yet done:** GSE245601 candidate expression analysis, pseudobulk
differential expression, single-cell candidate ranking, GSE111151
independent confirmation, mechanism/pathway analysis, druggability
assessment, normal-tissue expression analysis, final candidate
prioritization.

## Datasets

| Dataset | Role | Biological system | Status |
|---|---|---|---|
| Hany et al. 2023 CRISPR screen (*Sci Adv* 9:eadd3685) | Functional labels (gene-by-treatment interaction effect, E2+4-OHT vs E2) | MCF7-V, a drug-tolerant parental clone | Complete |
| GSE118713 | Resistance-associated bulk expression features | MCF7 parental / TAMR / FASR cell lines | Complete |
| GSE245601 | Human ex vivo tamoxifen-response single-cell context | 10 paired primary ER+/HER2− tumors, 10 µM tamoxifen vs control media, 12 h ex vivo, scRNA-seq | Preprocessing + malignant-cell ID complete; candidate analysis not started |
| GSE111151 | Independent post-hoc confirmation only (excluded from model training) | 4 parental lines + 7 TamR derivatives | Not yet used |

GSE245601's treatment contrast is **10 µM tamoxifen vs control media for
12 hours, ex vivo, on primary tumors** — this is deliberately not
described with the Hany-screen's "E2 vs E2+4-OHT" wording, which belongs
only to that unrelated CRISPR screen (see
[`docs/gse245601_PREANALYSIS.md` §0](docs/gse245601_PREANALYSIS.md#0-scientific-correction-to-the-prior-audit-frozen-here)).

## Current Findings

**CRISPR x bulk RNA (complete):**
- 19,103 genes fitted in the CRISPR reanalysis; 28 Gate-1 hits at the
  frozen FDR<0.1 threshold
- 13 of the 28 hits are sensitising knockout candidates (knockout depletes
  under tamoxifen, i.e. associated with increased tamoxifen sensitivity)
- USP34 currently has the strongest combined CRISPR + TAMR-vs-MCF7
  bulk-expression support among the 13 candidates
- PAICS remains a published benchmark gene, checked separately — it is
  **not** a Gate-1 CRISPR hit in this reanalysis (CRISPR FDR = 0.85)

**GSE245601 single-cell (preprocessing complete, candidates not yet examined):**
- 44,140 cells retained after QC and doublet removal across 20 primary
  tumor samples (10 patients, control + tamoxifen)
- 29,175 of those cells classified as epithelial
- InferCNV-based malignant-vs-non-malignant classification completed for
  all epithelial cells, with CopyKAT run per sample as an independent
  sensitivity check
- Only **Tumor_02, Tumor_03, and Tumor_07** currently satisfy the frozen
  ≥50-malignant-cells-per-arm pseudobulk-eligibility rule
  (`results/tables/gse245601_pair_eligibility.tsv`)
- Candidate-gene expression has **not** been examined in GSE245601 — no
  pseudobulk aggregation or differential expression has been run

No result in this repository should be read as "validated therapeutic
target," "restores sensitivity," or "safe target" — none of those claims
is established by the analyses completed so far.

## Repository Guide — Where is the code?

| Task | Main code | Main outputs | Notes |
|---|---|---|---|
| CRISPR screen reanalysis | [`src/labels.py`](src/labels.py), [`src/gate1_checks.py`](src/gate1_checks.py) | [`results/tables/gate1_decision.tsv`](results/tables/gate1_decision.tsv) | Gate-1 FDR<0.1 hit decision, CEG2 essentiality check |
| Bulk RNA — GSE118713 | [`src/gse118713_prep.py`](src/gse118713_prep.py), [`src/gse118713_qc.py`](src/gse118713_qc.py), [`scripts/analysis/gse118713_limma.R`](scripts/analysis/gse118713_limma.R) | [`results/tables/gse118713_de_summary.tsv`](results/tables/gse118713_de_summary.tsv) | limma DE across TAMR/FASR/MCF7 |
| CRISPR x bulk integration | [`src/crispr_gse118713_integration.py`](src/crispr_gse118713_integration.py) | [`results/tables/crispr_gse118713_master_table.tsv`](results/tables/crispr_gse118713_master_table.tsv) | Joins Gate-1 hits to GSE118713 expression |
| Candidate prioritization | [`src/candidate_evidence_summary.py`](src/candidate_evidence_summary.py) | [`results/tables/candidate_evidence_summary.tsv`](results/tables/candidate_evidence_summary.tsv), [`candidate_sensitisation_candidates.tsv`](results/tables/candidate_sensitisation_candidates.tsv) | Evidence-class hierarchy for the 13 sensitising candidates |
| PAICS benchmark | [`src/candidate_evidence_summary.py`](src/candidate_evidence_summary.py) | [`results/tables/candidate_paics_benchmark.tsv`](results/tables/candidate_paics_benchmark.tsv) | Published-benchmark check, separate from Gate-1 |
| NEBULA poster figures | [`src/nebula_plots_final.py`](src/nebula_plots_final.py) | [`results/figures/nebula_final/`](results/figures/nebula_final/) | Current, self-contained figure module |
| GSE245601 download + manifest | [`scripts/download/download_gse245601.sh`](scripts/download/download_gse245601.sh), [`src/gse245601_manifest.py`](src/gse245601_manifest.py) | [`results/tables/gse245601_sample_manifest.tsv`](results/tables/gse245601_sample_manifest.tsv) | Checksum-verified 26-sample manifest |
| GSE245601 candidate feature-space check | [`src/gse245601_feature_check.py`](src/gse245601_feature_check.py) | [`results/tables/gse245601_candidate_feature_availability.tsv`](results/tables/gse245601_candidate_feature_availability.tsv) | Gene-symbol presence only, no expression values |
| GSE245601 QC + doublets | [`scripts/analysis/gse245601_01_import_qc_doublets.R`](scripts/analysis/gse245601_01_import_qc_doublets.R) | [`results/tables/gse245601_qc_summary.tsv`](results/tables/gse245601_qc_summary.tsv) | Frozen thresholds, treatment-blind filter mask |
| GSE245601 normalization/clustering | [`scripts/analysis/gse245601_02_normalize_cluster.R`](scripts/analysis/gse245601_02_normalize_cluster.R) | [`results/tables/gse245601_cluster_summary.tsv`](results/tables/gse245601_cluster_summary.tsv) | LogNormalize, PCA, Louvain, UMAP |
| GSE245601 cell-type annotation | [`scripts/analysis/gse245601_03_annotate_celltypes.R`](scripts/analysis/gse245601_03_annotate_celltypes.R) | [`results/tables/gse245601_celltype_counts_per_sample.tsv`](results/tables/gse245601_celltype_counts_per_sample.tsv) | Candidate-gene-blind marker scoring |
| GSE245601 malignant-cell classification | [`scripts/analysis/gse245601_05_infercnv_malignant.R`](scripts/analysis/gse245601_05_infercnv_malignant.R) | [`results/tables/gse245601_malignant_cell_labels.tsv`](results/tables/gse245601_malignant_cell_labels.tsv) | InferCNV-based, treatment-blind |
| GSE245601 sensitivity check | [`scripts/analysis/gse245601_06_copykat_sensitivity.R`](scripts/analysis/gse245601_06_copykat_sensitivity.R) | [`results/tables/gse245601_malignancy_concordance.tsv`](results/tables/gse245601_malignancy_concordance.tsv) | Independent CopyKAT call, same cells |
| GSE245601 metadata freeze | [`scripts/analysis/gse245601_07_freeze_metadata.R`](scripts/analysis/gse245601_07_freeze_metadata.R) | [`results/tables/gse245601_cell_metadata_frozen.tsv.gz`](results/tables/gse245601_cell_metadata_frozen.tsv.gz), [`gse245601_pair_eligibility.tsv`](results/tables/gse245601_pair_eligibility.tsv) | Per-cell metadata + pseudobulk-eligibility flags; no pseudobulk run |
| Tests | [`tests/`](tests/) | — | One test module per analysis module — see [Code map](docs/CODE_MAP.md) |

## Key Figures

1. **[`results/figures/nebula_final/fig1_crispr_landscape.png`](results/figures/nebula_final/fig1_crispr_landscape.png)** — CRISPR hit landscape: all 28 Gate-1 hits, sensitising vs tolerance-associated, with the PAICS benchmark inset.
2. **[`results/figures/gse245601_preprocessing/02_umap_broad_cell_type.pdf`](results/figures/gse245601_preprocessing/02_umap_broad_cell_type.pdf)** — GSE245601 broad cell-type UMAP (epithelial, fibroblast, myeloid, endothelial, B/plasma, T/NK) across all 20 primary tumor samples.
3. **[`results/figures/gse245601_preprocessing/07_umap_malignant_vs_nonmalignant.pdf`](results/figures/gse245601_preprocessing/07_umap_malignant_vs_nonmalignant.pdf)** — Malignant vs. non-malignant epithelial cells, InferCNV-derived primary classification.
4. **[`results/figures/gse245601_preprocessing/08_malignant_counts_by_patient_condition.pdf`](results/figures/gse245601_preprocessing/08_malignant_counts_by_patient_condition.pdf)** — Malignant-cell counts per patient x condition, with the frozen ≥50-cell eligibility line — shows directly why only 3/10 tumor pairs currently qualify for pseudobulk analysis.

InferCNV supplied the **primary** malignant-cell classification; CopyKAT
was run as an **independent sensitivity check** on the same cells, not as
a second vote — see [Limitations](#limitations--current-methodological-question)
for how well the two agree. The classification should not be read as a
ground-truth cancer-cell label.

## Single-Cell Workflow

```
Cell Ranger filtered H5 matrices
  -> Seurat import
  -> QC + doublet removal
  -> normalization (LogNormalize, PCA, UMAP)
  -> broad cell-type annotation (candidate-gene-blind marker scoring)
  -> epithelial subset
  -> InferCNV malignant-cell classification (primary)
  -> CopyKAT sensitivity analysis (independent check)
  -> frozen malignant-cell set + pseudobulk-eligibility flags
```

**Pseudobulk aggregation and candidate-gene expression analysis have not
yet been performed** on GSE245601. The pipeline stops at a frozen,
per-cell metadata table; no gene-expression value for any of the 13
sensitising candidates (or PAICS) has been inspected in this dataset.

## Reproducibility

- Analyses run on KAUST Ibex (SLURM); raw sequencing data are **not**
  committed to this repository. Public processed matrices are downloaded
  from source by scripts under [`scripts/download/`](scripts/download/)
  (e.g. [`download_gse245601.sh`](scripts/download/download_gse245601.sh)).
- All file paths come from [`config/config.yaml`](config/config.yaml); no
  path is hardcoded in analysis code.
- Every module has a corresponding pytest module under
  [`tests/`](tests/) exercising its logic, not merely that it runs.
- Thresholds are declared in advance in [`PREANALYSIS.md`](PREANALYSIS.md)
  (CRISPR/bulk phase) and
  [`docs/gse245601_PREANALYSIS.md`](docs/gse245601_PREANALYSIS.md)
  (single-cell phase), with dated, append-only amendment logs — never
  edited in place after analysis begins.
- The conda/micromamba environment is pinned in
  [`environment.yml`](environment.yml).

## Limitations / Current Methodological Question

- Malignant-cell yield varies strongly by tumor: some samples (e.g.
  Tumor_01, Tumor_04) have only a handful of InferCNV-classified malignant
  cells, while others (Tumor_03) have hundreds.
- Only **3 of 10** paired tumors (Tumor_02, Tumor_03, Tumor_07) currently
  meet the frozen ≥50-malignant-cells-per-arm pseudobulk-eligibility
  threshold.
- InferCNV/CopyKAT concordance is highly variable between samples (from
  near-0% to near-100% agreement per sample; ~56% on average) — the two
  methods do not always agree on which epithelial cells are malignant.
- As a result, candidate-level single-cell inference is **not yet
  finalized**, and the statistical design for how to handle
  low-malignant-cell tumors (exclude, pool, or model uncertainty
  explicitly) is the current open decision point for this phase.

## HOW TO SHOW THIS REPO TO A RESEARCHER

- **Open first:** this README, then [`docs/PROJECT_STATUS.md`](docs/PROJECT_STATUS.md) for the one-page status.
- **CRISPR code:** [`src/labels.py`](src/labels.py) and [`src/gate1_checks.py`](src/gate1_checks.py) (see the Repository Guide table above).
- **Bulk RNA code:** [`src/gse118713_prep.py`](src/gse118713_prep.py), [`src/gse118713_qc.py`](src/gse118713_qc.py), [`scripts/analysis/gse118713_limma.R`](scripts/analysis/gse118713_limma.R).
- **Single-cell code:** [`scripts/analysis/gse245601_01_import_qc_doublets.R`](scripts/analysis/gse245601_01_import_qc_doublets.R) through [`_08_qc_figures.R`](scripts/analysis/gse245601_08_qc_figures.R), in numeric order.
- **Malignant-cell code:** [`scripts/analysis/gse245601_05_infercnv_malignant.R`](scripts/analysis/gse245601_05_infercnv_malignant.R) (primary) and [`scripts/analysis/gse245601_06_copykat_sensitivity.R`](scripts/analysis/gse245601_06_copykat_sensitivity.R) (independent check).
- **Latest figures:** [`results/figures/nebula_final/`](results/figures/nebula_final/) (CRISPR/bulk poster figures) and [`results/figures/gse245601_preprocessing/`](results/figures/gse245601_preprocessing/) (single-cell QC/annotation figures).
