# Project Workflow

Narrative, phase-by-phase walkthrough of the full 16-phase analysis
pipeline, with code, output, and report pointers for each phase. See
[`../README.md`](../README.md) for the one-page summary and
[`RESULTS_GUIDE.md`](RESULTS_GUIDE.md) for what each output table/figure
means and whether it is current or superseded.

All paths below are relative to the repository root. Phases 1-11 are
**frozen** — later phases read from them but never alter their values.
Phases 12-16 characterize the frozen shortlist's top two candidates
(USP34, VEZF1) without reopening candidate discovery.

---

## Phase 1 — CRISPR screen reanalysis (Hany et al. 2023)

**Question:** which genes' knockout changes tamoxifen response in the Hany
et al. genome-wide screen?

**Code:** [`src/labels.py`](../src/labels.py) (guide/gene-level label
construction), [`src/gate1_checks.py`](../src/gate1_checks.py) (Gate-1
FDR<0.1 decision, CEG2 essentiality check).

**Output:** [`results/tables/gate1_decision.tsv`](../results/tables/gate1_decision.tsv)
— 19,103 genes fitted; 28 Gate-1 hits at FDR<0.1; 13 are sensitising
knockouts (the direction of interest for a resensitisation strategy).

## Phase 2 — GSE118713 bulk RNA resistance analysis

**Question:** which of the CRISPR hits also show resistance-associated
expression changes in an independent bulk-RNA resistance model?

**Code:** [`src/gse118713_prep.py`](../src/gse118713_prep.py),
[`src/gse118713_qc.py`](../src/gse118713_qc.py),
[`scripts/analysis/gse118713_limma.R`](../scripts/analysis/gse118713_limma.R).

**Output:** [`results/tables/gse118713_de_summary.tsv`](../results/tables/gse118713_de_summary.tsv).

## Phase 3 — CRISPR x bulk-RNA integration

**Code:** [`src/crispr_gse118713_integration.py`](../src/crispr_gse118713_integration.py),
[`src/candidate_evidence_summary.py`](../src/candidate_evidence_summary.py).

**Output:** [`results/tables/candidate_evidence_summary.tsv`](../results/tables/candidate_evidence_summary.tsv)
— evidence-class hierarchy for the 13 sensitising candidates.

## Phase 4 — NEBULA poster figures

**Code:** [`src/nebula_plots_final.py`](../src/nebula_plots_final.py).
**Output:** [`results/figures/nebula_final/`](../results/figures/nebula_final/).

## Phase 5 — GSE245601 single-cell preprocessing + malignant-cell ID

**Code:** [`scripts/analysis/gse245601_01_import_qc_doublets.R`](../scripts/analysis/gse245601_01_import_qc_doublets.R)
through `_07_freeze_metadata.R`, in numeric order. InferCNV supplies the
**primary** malignant-cell classification; CopyKAT is an **independent
sensitivity check** on the same cells (see
[`CNV_METHOD_AUDIT.md`](CNV_METHOD_AUDIT.md)).

**Output:** [`results/tables/gse245601_malignant_cell_labels.tsv`](../results/tables/gse245601_malignant_cell_labels.tsv),
[`gse245601_pair_eligibility.tsv`](../results/tables/gse245601_pair_eligibility.tsv)
(only Tumor_02/03/07 meet the frozen >=50-malignant-cell pseudobulk rule).

## Phase 6 — GSE245601 candidate-level expression / pseudobulk

**Code:** [`scripts/analysis/gse245601_13_build_pseudobulk.R`](../scripts/analysis/gse245601_13_build_pseudobulk.R),
`_15_pseudobulk_edger.R`, [`src/gse245601_candidate_deepdive_data.py`](../src/gse245601_candidate_deepdive_data.py)
and the other `gse245601_candidate_deepdive_*.py` modules.

**Output:** [`results/tables/gse245601_candidate_deepdive/`](../results/tables/gse245601_candidate_deepdive/),
report: [`GSE245601_CANDIDATE_DEEPDIVE.md`](GSE245601_CANDIDATE_DEEPDIVE.md).

## Phase 7 — GSE240112 primary-vs-recurrent scRNA-seq (4th evidence layer)

**Code:** [`scripts/analysis/gse240112_01_extract_h5seurat.R`](../scripts/analysis/gse240112_01_extract_h5seurat.R)
through `_05_normal_tissue_context.R`, [`src/gse240112_candidate_integration.py`](../src/gse240112_candidate_integration.py).

**Output:** [`results/tables/gse240112/`](../results/tables/gse240112/),
report: [`GSE240112_ANALYSIS_REPORT.md`](GSE240112_ANALYSIS_REPORT.md).

## Phase 8 — GSE111151 independent resistance-model validation (5th evidence layer)

**Code:** [`scripts/analysis/gse111151_01_build_count_matrix.R`](../scripts/analysis/gse111151_01_build_count_matrix.R),
`_02_edger.R`, [`src/gse111151_candidate_integration.py`](../src/gse111151_candidate_integration.py).

**Output:** [`results/tables/gse111151/`](../results/tables/gse111151/),
report: [`GSE111151_ANALYSIS_REPORT.md`](GSE111151_ANALYSIS_REPORT.md).

## Phase 9 — Unbiased genome-wide cross-dataset integration

**Question:** across all five evidence layers (CRISPR, GSE118713, GSE245601,
GSE240112, GSE111151), which genes have the most consistent, coherent
evidence — checked genome-wide, not just within the 13-candidate shortlist?

**Code:** [`src/cross_dataset_genomewide_*.py`](../src/).

**Output:** [`results/tables/cross_dataset_genomewide/`](../results/tables/cross_dataset_genomewide/)
(`all_genes_cross_dataset_evidence_with_ranking.tsv` is the master table),
report: [`CROSS_DATASET_GENOMEWIDE_ANALYSIS_REPORT.md`](CROSS_DATASET_GENOMEWIDE_ANALYSIS_REPORT.md).

## Phase 10 — Candidate adjudication

**Code:** [`src/candidate_adjudication_*.py`](../src/).
**Output:** [`results/tables/candidate_adjudication/`](../results/tables/candidate_adjudication/)
— narrows to the 7 `MULTIMODAL_STRONG` genes.

## Phase 11 — Evidence freeze (FROZEN)

**Output:** [`results/tables/evidence_freeze/THERAPEUTIC_SHORTLIST_FREEZE.tsv`](../results/tables/evidence_freeze/THERAPEUTIC_SHORTLIST_FREEZE.tsv),
report: [`THERAPEUTIC_SHORTLIST_FREEZE.md`](THERAPEUTIC_SHORTLIST_FREEZE.md).
**Frozen four-gene shortlist: USP34 > VEZF1 > EML5 > CITED2.** No later
phase alters these values.

## Phase 12 — Systems / pathway network mapping

**Code:** [`src/systems_network_*.py`](../src/).
**Output:** [`results/networks/systems_network/`](../results/networks/systems_network/),
[`results/tables/systems_network/`](../results/tables/systems_network/),
report: [`results/reports/systems_network/four_candidate_network_audit.md`](../results/reports/systems_network/four_candidate_network_audit.md).

## Phase 13 — Literature mechanism review

**Output:** [`results/reports/literature_mechanism/four_candidate_mechanism_review.md`](../results/reports/literature_mechanism/four_candidate_mechanism_review.md),
[`verified_references.tsv`](../results/reports/literature_mechanism/verified_references.tsv).

## Phase 14 — Independent validation (TCGA-BRCA + DepMap 26Q1)

**Code:** [`src/independent_validation_*.py`](../src/).
**Output:** [`results/tables/independent_validation/`](../results/tables/independent_validation/),
report: [`results/reports/independent_validation/four_candidate_TCGA_DepMap_review.md`](../results/reports/independent_validation/four_candidate_TCGA_DepMap_review.md).
Explicitly **not** a tamoxifen-resistance cohort — see caveats in
[`../README.md`](../README.md#e-datasets).

## Phase 15 — Lead-target deep dive + druggability/safety

**Question:** for USP34 and VEZF1 only, what is the structural/chemical
druggability, and what safety liabilities exist in normal tissue?

**Code:** [`src/lead_target_deep_dive_*.py`](../src/),
[`src/druggability_safety_*.py`](../src/).
**Output:** [`results/tables/lead_target_deep_dive/`](../results/tables/lead_target_deep_dive/),
[`results/tables/druggability_safety/`](../results/tables/druggability_safety/),
reports: [`results/reports/lead_target_deep_dive/USP34_VEZF1_translational_deep_dive.md`](../results/reports/lead_target_deep_dive/USP34_VEZF1_translational_deep_dive.md),
[`results/reports/druggability_safety/four_candidate_druggability_safety_review.md`](../results/reports/druggability_safety/four_candidate_druggability_safety_review.md).

## Phase 16 — Final translational design + GDSC pharmacogenomics

**Code:** [`src/final_translational_*.py`](../src/),
[`src/final_pharmacogenomics_*.py`](../src/).
**Output:** [`results/tables/final_translational/`](../results/tables/final_translational/),
[`results/tables/final_pharmacogenomics/`](../results/tables/final_pharmacogenomics/),
reports: [`results/reports/final_translational/final_USP34_VEZF1_translational_plan.md`](../results/reports/final_translational/final_USP34_VEZF1_translational_plan.md),
[`results/reports/final_pharmacogenomics/USP34_VEZF1_GDSC_review.md`](../results/reports/final_pharmacogenomics/USP34_VEZF1_GDSC_review.md).
This phase designed (EXP-1 through EXP-5) but did **not run** any wet-lab
experiment, and performed a real fpocket structural pocket analysis
(concluding `DOCKING_NOT_YET_JUSTIFIED`) and a real GDSC Release 8.5
drug-response correlation lookup for USP34/VEZF1 only.

---

## What comes after Phase 16

Per this project's explicit scope discipline, the next steps (not yet
started) are: poster assembly -> final written project story -> a fresh
independent blind audit -> fixing any issues found -> a final Git
freeze/tag -> **only then** consider making the repository public.
