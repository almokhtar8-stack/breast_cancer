# GSE245601 candidate deep-dive: object and metadata audit

**Date written:** 2026-08-13. This is a read-only exploratory
decomposition of the acute 12h-tamoxifen GSE245601 pseudobulk result for
the four FROZEN therapeutic candidates (USP34, VEZF1, EML5, CITED2). It
does **not** modify `results/tables/evidence_freeze/`,
`results/figures/evidence_freeze/`, or `docs/THERAPEUTIC_SHORTLIST_FREEZE.md`,
does not rerun QC/clustering/InferCNV/CopyKAT/the pseudobulk-edgeR
pipeline, and does not alter any frozen label.

## Source object

`data/processed/gse245601/seurat_clustered/annotated.rds` -- the
already-frozen, fully-annotated Seurat object (5,5 `SeuratObject`
version). 44,140 total cells across all lineages; RNA assay with
`counts`, `data` (log-normalized), and `scale.data` layers; `pca` and
`umap` reductions (reused as-is, never recomputed).

`merged@meta.data` columns used: `patient`, `condition`, `sample_id`,
`cell_id`, `seurat_clusters` (the project's only clustering -- computed
once, on the full multi-lineage object at `RNA_snn_res.0.8`; there is no
separate epithelial-only subclustering, so "epithelial clusters" in this
deep-dive means the subset of `seurat_clusters` IDs that contain
epithelial cells), `broad_cell_type` (`epithelial` for 29,175 cells, used
only to confirm the frozen epithelial cell set matches).

## Frozen per-cell metadata (reused, not recomputed)

`results/tables/gse245601_pseudobulk/malignant_vs_nonmalignant_cell_level_summary.tsv`
-- 29,175 rows, one per epithelial cell (matches
`broad_cell_type == "epithelial"` in the Seurat object exactly). Columns
reused as-is: `cell_id`, `sample_id`, `patient`, `condition`,
`malignancy_status` (`malignant` / `non-malignant epithelial`, from the
already-frozen InferCNV calling), `nCount_RNA` (per-cell library size),
`umap_1`/`umap_2` (frozen UMAP coordinates), `cnv_score`,
`score_epithelial_program`, `cell_cycle_phase`.

10 patients: `Tumor_01` ... `Tumor_10`. 2 conditions: `Control`,
`Tamoxifen` -- confirmed to be the acute 12-hour ex vivo tamoxifen
exposure (`gse245601_PREANALYSIS.md`), never an established/chronic
resistant state. `malignancy_status` breakdown: 1,739 malignant, 27,436
non-malignant epithelial.

## New extraction (this phase only)

`scripts/analysis/gse245601_18_extract_candidate_deepdive.R` loads the
frozen Seurat object, subsets to the exact frozen 29,175-cell epithelial
set (order-enforced to match the frozen metadata table), confirms each of
the four genes is present exactly once in the feature space, and writes
per-cell raw counts (`counts` layer) and log-normalized expression
(`data` layer) for each gene, plus each cell's `seurat_clusters`
assignment -- to
`results/tables/gse245601_candidate_deepdive/candidate_per_cell_expression.tsv`.
Joined with the frozen metadata table (by `cell_id`, one-to-one, no rows
lost or duplicated) by `src.gse245601_candidate_deepdive_data.load_per_cell_table`,
which is the single source every table/figure in this deep-dive reads
through.

## Pseudobulk sources (reused, not recomputed)

Track A (all epithelial): `results/tables/gse245601_pseudobulk/track_a_epithelial_counts.tsv.gz`
(33,538 genes x 20 samples) + `track_a_epithelial_metadata.tsv`.
Track B (strict malignant, >=50-cell eligibility already applied at
pseudobulk-construction time to Tumor_02/03/07 only):
`results/tables/gse245601_pseudobulk/track_b_malignant_counts.tsv.gz`
(33,538 genes x 6 samples) + `track_b_malignant_metadata.tsv`. Genome-wide
DE results: `track_{a,b}_genomewide_de.tsv.gz`.

## Cluster support rule (declared before inspecting any gene-specific result)

A cluster is "sufficiently represented" iff at least 3 tumors have >=10
cells in BOTH Control and Tamoxifen within that cluster. Of the 16
epithelial-populated `seurat_clusters` IDs, 5 meet this bar: clusters 3,
8, 17, 18, 25 (3-8 tumors supported each). The other 11 are reported in
`cluster_candidate_response.tsv` for transparency but excluded from any
interpretive claim.

## Phase 17 (response-state exploration): explicitly omitted

A rigorous, non-circular acute-response transcriptional-state score
(excluding the 4 candidate genes) is most defensibly computed with
Seurat's `AddModuleScore`, which corrects for each signature gene's own
expression-level control bin -- reproducing that correctly would require
a new R computation step beyond simple frozen-data extraction. A
simplified mean-expression proxy was considered and rejected as an
avoidable extra methodological choice not worth the added risk, given
that Phases 4-20 already directly and comprehensively answer the primary
question (what structure underlies the acute pseudobulk result) without
it. This phase is marked OPTIONAL in the task specification and is
omitted here with this documented justification, per that phase's own
explicit allowance.

## Statistical unit discipline (binding for every phase below)

The tumor/patient is the biological replicate. scRNA-seq is destructive:
Control and Tamoxifen cells are different cells from matched pieces of
the same tumor, never the same cell measured twice. Individual cells are
never connected by a line across conditions, never called a trajectory or
time series, and never used as independent n for a treatment-effect
p-value. Cell-level distributions are shown descriptively only.
