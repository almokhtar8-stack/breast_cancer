# Code Map

Practical lookup sheet: for a given analysis task, where is the code, what
does it read, what does it write, and what tests cover it. Built from the
files that actually exist in this repository as of 2026-08-11 — see
[`../README.md`](../README.md) for the narrative overview and
[`PROJECT_STATUS.md`](PROJECT_STATUS.md) for current status.

All paths below are relative to the repository root.

---

## CRISPR (Hany et al. 2023 screen reanalysis)

**Purpose:** Build per-gene tamoxifen-response labels from the published
screen and apply the pre-registered Gate-1 significance/essentiality
decision.

**Main script(s):**
- [`src/labels.py`](../src/labels.py) — guide- and gene-level label
  construction (MAGeCK-VISPR-normalised counts, E2 vs E2+4-OHT interaction
  term).
- [`src/gate1_checks.py`](../src/gate1_checks.py) — Gate-1 FDR<0.1
  decision, CEG2 essentiality-contamination check, direction-sanity
  extraction.
- [`scripts/download/download_ceg2.py`](../scripts/download/download_ceg2.py)
  — deterministic download of the Hart et al. 2017 CEG2 reference list.

**Inputs:** Hany et al. 2023 Data S1 (`config.yaml: data.raw.hany_data_s1`).

**Outputs:**
- `data/processed/labels.parquet` (frozen, force-added to git)
- [`results/tables/gate1_decision.tsv`](../results/tables/gate1_decision.tsv)
- [`results/tables/gate1_essentiality_summary.tsv`](../results/tables/gate1_essentiality_summary.tsv)
- [`results/tables/gate1_direction_sanity.tsv`](../results/tables/gate1_direction_sanity.tsv)

**Tests:** [`tests/test_labels.py`](../tests/test_labels.py),
[`tests/test_gate1_checks.py`](../tests/test_gate1_checks.py),
[`tests/test_download_ceg2.py`](../tests/test_download_ceg2.py)

---

## Bulk RNA — GSE118713

**Purpose:** Gene-level TPM preparation, QC, and limma differential
expression across MCF7 / TAMR / FASR to produce resistance-associated
expression evidence.

**Main script(s):**
- [`src/gse118713_prep.py`](../src/gse118713_prep.py) — transcript-to-gene
  TPM aggregation.
- [`src/gse118713_expression_filter.py`](../src/gse118713_expression_filter.py)
  — frozen-matrix validation + expression filtering.
- [`src/gse118713_qc.py`](../src/gse118713_qc.py) — sample QC, correlations, PCA.
- [`scripts/analysis/gse118713_limma.R`](../scripts/analysis/gse118713_limma.R)
  (+ [`gse118713_limma_lib.R`](../scripts/analysis/gse118713_limma_lib.R))
  — limma differential expression (subprocess call from
  [`src/gse118713_phase2b.py`](../src/gse118713_phase2b.py)).
- [`src/gse118713_de_summary.py`](../src/gse118713_de_summary.py) —
  per-contrast DE summary statistics.
- [`src/gse118713_tamr_specificity.py`](../src/gse118713_tamr_specificity.py)
  — joins the three preregistered contrasts.
- [`src/gse118713_unredact.py`](../src/gse118713_unredact.py) — controlled
  release of the RCOR1/KDM1A blind, per the 2026-08-10 `PREANALYSIS.md`
  amendment.
- [`src/gse118713_phase2b.py`](../src/gse118713_phase2b.py) — orchestrates
  the above, in order.

**Inputs:** GSE118713 transcript-level TPM archive
(`config.yaml: data.raw.gse118713_transcript_tpm`).

**Outputs:**
- `data/processed/gse118713_gene_tpm.parquet`, `gse118713_gene_tpm_filtered.tsv.gz` (frozen, force-added)
- [`results/tables/gse118713_de_summary.tsv`](../results/tables/gse118713_de_summary.tsv)
- [`results/tables/gse118713_differential_expression.tsv.gz`](../results/tables/gse118713_differential_expression.tsv.gz)
- [`results/figures/gse118713_pca.pdf`](../results/figures/gse118713_pca.pdf), [`gse118713_sample_correlation.pdf`](../results/figures/gse118713_sample_correlation.pdf)

**Tests:** [`tests/test_gse118713_prep.py`](../tests/test_gse118713_prep.py),
[`tests/test_gse118713_expression_filter.py`](../tests/test_gse118713_expression_filter.py),
[`tests/test_gse118713_qc.py`](../tests/test_gse118713_qc.py),
[`tests/test_gse118713_limma.py`](../tests/test_gse118713_limma.py),
[`tests/test_gse118713_de_summary.py`](../tests/test_gse118713_de_summary.py),
[`tests/test_gse118713_tamr_specificity.py`](../tests/test_gse118713_tamr_specificity.py),
[`tests/test_gse118713_unredact.py`](../tests/test_gse118713_unredact.py),
[`tests/test_gse118713_phase2b.py`](../tests/test_gse118713_phase2b.py)

---

## GSE111151 (prepared, not yet used downstream)

**Purpose:** Independent post-hoc confirmation dataset (4 parental lines,
7 TamR derivatives) — prepared but explicitly excluded from model
training and not yet integrated into any candidate analysis.

**Main script(s):** [`src/gse111151_prep.py`](../src/gse111151_prep.py)

**Outputs:** `data/processed/gse111151_log2cpm.parquet` (frozen, force-added),
[`results/tables/gse111151_sample_metadata.tsv`](../results/tables/gse111151_sample_metadata.tsv)

**Tests:** [`tests/test_gse111151_prep.py`](../tests/test_gse111151_prep.py)

---

## Candidate prioritization

**Purpose:** Join the 28 Gate-1 CRISPR hits to GSE118713 expression and
classify each of the 13 sensitising-knockout candidates into an
evidence-class hierarchy (direction-aware, no new statistical test).

**Main script(s):**
- [`src/crispr_gse118713_integration.py`](../src/crispr_gse118713_integration.py)
  — master integration table.
- [`src/candidate_evidence_summary.py`](../src/candidate_evidence_summary.py)
  — evidence-class hierarchy + PAICS benchmark check.

**Inputs:** `data/processed/labels.parquet`,
[`results/tables/gate1_decision.tsv`](../results/tables/gate1_decision.tsv),
GSE118713 DE table.

**Outputs:**
- [`results/tables/crispr_gse118713_master_table.tsv`](../results/tables/crispr_gse118713_master_table.tsv)
- [`results/tables/candidate_evidence_summary.tsv`](../results/tables/candidate_evidence_summary.tsv)
- [`results/tables/candidate_sensitisation_candidates.tsv`](../results/tables/candidate_sensitisation_candidates.tsv)
- [`results/tables/candidate_paics_benchmark.tsv`](../results/tables/candidate_paics_benchmark.tsv) (PAICS benchmark, kept separate from Gate-1)
- [`results/tables/candidate_tolerance_hits.tsv`](../results/tables/candidate_tolerance_hits.tsv), [`candidate_shortlist.tsv`](../results/tables/candidate_shortlist.tsv), [`candidate_secondary_context.tsv`](../results/tables/candidate_secondary_context.tsv)

**Tests:** [`tests/test_crispr_gse118713_integration.py`](../tests/test_crispr_gse118713_integration.py),
[`tests/test_candidate_evidence_summary.py`](../tests/test_candidate_evidence_summary.py)

---

## NEBULA figures

**Purpose:** Poster-ready figures for the CRISPR/bulk candidate
prioritization results.

**Main script(s):**
- [`src/nebula_plots_final.py`](../src/nebula_plots_final.py) — **current**,
  self-contained module (no dependency on the retired `nebula_plots_v2`
  through `v8` iteration modules), produces the visually-approved figure
  set (Figures 1-3 = approved "v6" design, Figures 4-5 = approved "v8"
  design).
- [`src/nebula_plots.py`](../src/nebula_plots.py) — earlier module that
  produces Figures 1-5 in the intermediate "v6" appearance; still valid
  and referenced by `config.yaml`, kept alongside `nebula_plots_final.py`.

**Inputs:** `results/tables/candidate_evidence_summary.tsv`,
`crispr_gse118713_master_table.tsv`, GSE118713 filtered TPM matrix — no
recomputation, plotting only.

**Outputs:**
- [`results/figures/nebula_final/`](../results/figures/nebula_final/) — current poster figure set + `figure_manifest.tsv` + contact sheet (tracked in git)
- `results/figures/nebula/` — earlier figure set, superseded by `nebula_final/`; regenerated locally by `src/nebula_plots.py`, not tracked in git

**Tests:** [`tests/test_nebula_plots_final.py`](../tests/test_nebula_plots_final.py),
[`tests/test_nebula_plots.py`](../tests/test_nebula_plots.py)

---

## Single-cell — GSE245601

Human ex vivo tamoxifen-response validation. Governed by
[`docs/gse245601_PREANALYSIS.md`](gse245601_PREANALYSIS.md). Pipeline
stops after malignant-cell freezing — **no pseudobulk or candidate
expression analysis exists yet** (see [`PROJECT_STATUS.md`](PROJECT_STATUS.md)).

### Download + manifest

**Main script(s):**
[`scripts/download/download_gse245601.sh`](../scripts/download/download_gse245601.sh),
[`src/gse245601_manifest.py`](../src/gse245601_manifest.py)

**Outputs:** [`results/tables/gse245601_sample_manifest.tsv`](../results/tables/gse245601_sample_manifest.tsv)
(checksum-verified 26-sample manifest, 20 primary + 6 excluded)

**Tests:** [`tests/test_gse245601_manifest.py`](../tests/test_gse245601_manifest.py)

### Candidate feature-space check

**Main script(s):** [`src/gse245601_feature_check.py`](../src/gse245601_feature_check.py)
— reads only `matrix/features/*` of a Cell Ranger H5 (gene metadata), never
expression values.

**Outputs:** [`results/tables/gse245601_candidate_feature_availability.tsv`](../results/tables/gse245601_candidate_feature_availability.tsv)

**Tests:** [`tests/test_gse245601_feature_check.py`](../tests/test_gse245601_feature_check.py)

### QC + doublets

**Main script(s):** [`scripts/analysis/gse245601_01_import_qc_doublets.R`](../scripts/analysis/gse245601_01_import_qc_doublets.R)
— Seurat import, frozen QC thresholds (treatment-blind filter mask), DoubletFinder.

**Outputs:** [`results/tables/gse245601_qc_summary.tsv`](../results/tables/gse245601_qc_summary.tsv)

### Seurat normalization / clustering

**Main script(s):** [`scripts/analysis/gse245601_02_normalize_cluster.R`](../scripts/analysis/gse245601_02_normalize_cluster.R)
— LogNormalize, 2000 variable features, PCA(50), Louvain clustering, UMAP.

**Outputs:** [`results/tables/gse245601_cluster_summary.tsv`](../results/tables/gse245601_cluster_summary.tsv)

### Cell-type annotation

**Main script(s):** [`scripts/analysis/gse245601_03_annotate_celltypes.R`](../scripts/analysis/gse245601_03_annotate_celltypes.R)
— candidate-gene-blind marker-set module scoring, per-cluster argmax
assignment.

**Outputs:** [`results/tables/gse245601_celltype_marker_evidence.tsv`](../results/tables/gse245601_celltype_marker_evidence.tsv),
[`gse245601_celltype_counts_per_sample.tsv`](../results/tables/gse245601_celltype_counts_per_sample.tsv)

### InferCNV (primary malignant-cell classification)

**Main script(s):**
[`scripts/analysis/gse245601_04_gene_order_file.R`](../scripts/analysis/gse245601_04_gene_order_file.R)
(prerequisite gene-order file), then
[`scripts/analysis/gse245601_05_infercnv_malignant.R`](../scripts/analysis/gse245601_05_infercnv_malignant.R)
— per-sample InferCNV run + adaptive CNV-score/correlation thresholding
("ng_2021_and_thresholding" procedure, reconstructed and frozen in
[`docs/gse245601_PREANALYSIS.md` §9](gse245601_PREANALYSIS.md)).

**Outputs:** [`results/tables/gse245601_malignant_cell_labels.tsv`](../results/tables/gse245601_malignant_cell_labels.tsv),
[`gse245601_malignant_summary_per_sample.tsv`](../results/tables/gse245601_malignant_summary_per_sample.tsv)

### CopyKAT (independent sensitivity check)

**Main script(s):** [`scripts/analysis/gse245601_06_copykat_sensitivity.R`](../scripts/analysis/gse245601_06_copykat_sensitivity.R)
— run per sample on the same mixed cell population, own internal
normal-cell inference (genuinely independent of the InferCNV call).

**Outputs:** [`results/tables/gse245601_copykat_sensitivity_labels.tsv`](../results/tables/gse245601_copykat_sensitivity_labels.tsv),
[`gse245601_malignancy_concordance.tsv`](../results/tables/gse245601_malignancy_concordance.tsv)

### Metadata freeze

**Main script(s):** [`scripts/analysis/gse245601_07_freeze_metadata.R`](../scripts/analysis/gse245601_07_freeze_metadata.R)
— combines annotation + InferCNV + CopyKAT into one per-cell table;
computes (does not run) pseudobulk pair-eligibility.

**Outputs:** [`results/tables/gse245601_cell_metadata_frozen.tsv.gz`](../results/tables/gse245601_cell_metadata_frozen.tsv.gz)
(44,140 cells), [`gse245601_malignant_counts_per_sample.tsv`](../results/tables/gse245601_malignant_counts_per_sample.tsv),
[`gse245601_pair_eligibility.tsv`](../results/tables/gse245601_pair_eligibility.tsv)

### Figures

**Main script(s):** [`scripts/analysis/gse245601_08_qc_figures.R`](../scripts/analysis/gse245601_08_qc_figures.R)
— QC/annotation/malignancy figures only; no candidate gene plotted anywhere.

**Outputs:** [`results/figures/gse245601_preprocessing/`](../results/figures/gse245601_preprocessing/)

### Tests (all GSE245601 preprocessing scripts)

[`tests/test_gse245601_preprocessing.py`](../tests/test_gse245601_preprocessing.py)
— structural candidate-blindness checks on every script above, plus
output-table structural tests (skipped, not failed, if a table hasn't
been generated yet in a given environment).
