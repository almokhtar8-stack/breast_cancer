# Data Audit (poster-story-v1)

This phase reuses the exact data inventory already documented in
`results/reports/poster_exploration_v2/DATA_FOR_VISUALIZATION_AUDIT.md`
(genome-wide CRISPR, GSE118713/GSE111151/GSE240112/GSE245601 per-sample
transcriptomics, DepMap 26Q1, structures) -- not repeated here in full.
This document records only what is NEW or being used differently this
phase.

## New real layer: malignant vs. non-malignant (GSE245601)

| Need | Source | Unit |
|---|---|---|
| Per-gene delta + significance, all 13 sensitising hits | `results/tables/gse245601_candidate_integration/malignant_vs_nonmalignant_candidates.tsv` | gene (already-frozen test result) |
| Same delta, already in the master cross-dataset table | `results/tables/cross_dataset_genomewide/all_genes_cross_dataset_evidence_with_ranking.tsv`, columns `gse245601_malignant_log2fc/_p/_fdr` | gene, genome-wide |
| Real per-patient pseudobulk raw counts | `results/tables/gse245601_pseudobulk/malignant_vs_nonmalignant_counts.tsv.gz` + `malignant_vs_nonmalignant_metadata.tsv` (library sizes) | 5 patients (Tumor_02/03/07/09/10) x malignant/non-malignant |
| Cell-level malignancy calls (copyKAT-derived) | `results/tables/gse245601_copykat_sensitivity_labels.tsv` (29,175 cells) | cell |

Confirmed: `gse245601_malignant_log2fc` in the frozen master table and the
per-patient pseudobulk raw counts are two views of the SAME already-frozen
comparison (the master table is the summary statistic; the raw counts let
real per-patient points be plotted). Real per-patient log2(CPM+1) values
are computed here from the frozen raw counts + frozen library sizes using
the IDENTICAL disclosed formula already used and verified in
`poster_exploration_v2_data.py` for the acute-response track (verified
again this phase against the frozen master-table delta as a consistency
check -- see `tests/test_poster_story_v1.py`).

**copyKAT's role, stated precisely:** copyKAT is the tool that produced
the per-cell malignancy calls (`sensitivity_malignancy_label`) which
define the "malignant" vs. "non-malignant" grouping used above. It is not
visualized as its own QC figure in this phase (that would be
methods-supplement material) -- its provenance is stated once, in the
disease/clinical-context figure's caption.

## Network coverage (re-confirmed unchanged from v2)

`four_candidate_direct_neighbors.tsv` / `candidate_pathway_membership.tsv`
still cover only the original 4 candidates (USP34: 10 neighbor rows,
CITED2: 18, VEZF1: 1, EML5: 1); KDM1A and TLK2 have zero rows in either
table. Unchanged conclusion from v2: not usable as a main-sequence panel
for a 4-focus-gene story; kept as an explicit backup item only.

## TCGA (re-confirmed unchanged)

USP34/VEZF1 paired tumour-vs-normal associations remain weak/non-significant
(FDR 0.21 / 0.90, from `results/reports/post_audit/SCIENCE_FREEZE_REPORT.md`
Section 7 lineage); TCGA was never run for KDM1A/TLK2. Omitted entirely
from this phase's small, high-signal set (see `STORY_PLAN.md`).
