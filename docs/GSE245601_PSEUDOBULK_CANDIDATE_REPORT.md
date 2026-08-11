# GSE245601 pseudobulk single-cell validation: pipeline and results

**Status:** frozen InferCNV labels unchanged; frozen 0.01 CNV-score cutoff
unchanged; CopyKAT untouched as sensitivity-only. No candidate-gene
expression was inspected until after `docs/gse245601_PREANALYSIS.md`
section 13's design freeze and the Phase 2 design/eligibility audit
(below) passed all 20 checks. Design decisions per section 13 and the
task instructions.

## Pipeline

| Phase | Script/module | Output |
|---|---|---|
| 2. Design/eligibility audit | `gse245601_12_pseudobulk_design_audit.R` | `results/tables/gse245601_pseudobulk/design_eligibility_audit.tsv` (20/20 checks passed) |
| 3. Pseudobulk construction | `gse245601_13_build_pseudobulk.R` | Track A (20 samples, 29,175 cells) + Track B (6 samples, 1,401 cells) counts/metadata |
| 4. QC | `src/gse245601_pseudobulk_qc.py` | library size/detected genes, sample correlation, PCA (both tracks) |
| 5. edgeR DE | `gse245601_15_pseudobulk_edger.R` | genome-wide `~ patient + treatment` results, both tracks |
| 6. Candidate extraction | `src/gse245601_candidate_extraction.py` | 13-candidate + PAICS tables, both tracks |
| 7. Visualization | `src/gse245601_candidate_visualization.py` | paired grids, response heatmaps, Track A/B summary |
| 8. Malignant vs non-malignant | `gse245601_16_malignant_vs_nonmalignant.R` + `src/gse245601_malignant_vs_nonmalignant.py` | cell-level descriptive table + 5-patient secondary candidate table |
| 9. Integration | `src/gse245601_candidate_integration.py` | one row per candidate, CRISPR+bulk+SC, no composite score |

## Phase 2: design/eligibility audit (20/20 passed)

Verified directly from the frozen Seurat object's metadata and the frozen
label/eligibility tables, **before any pseudobulk aggregation or
candidate-gene expression was touched**: exactly 20 samples
(Tumor_01-10 x Control/Tamoxifen), one Control + one Tamoxifen per
patient, epithelial and strict-malignant cell counts reproduce the frozen
`gse245601_malignant_summary_per_sample.tsv` exactly, recomputed Track B
eligibility reproduces the frozen `gse245601_pair_eligibility.tsv` exactly
and equals exactly {Tumor_02, Tumor_03, Tumor_07}, and all 13 candidates +
PAICS map to exactly one gene-name row in the raw 33,538-gene feature
space.

## Phase 3-4: pseudobulk construction and QC

Track A: 29,175 epithelial cells -> 20 pseudobulk samples (library size
1.6M-83.8M, 15,117-23,061 detected genes). Track B: 1,401 frozen-malignant
cells (317 Tumor_02 + 858 Tumor_03 + 226 Tumor_07, matching the frozen
per-sample counts exactly) -> 6 pseudobulk samples.

**QC result: samples cluster overwhelmingly by patient in both tracks**
(PCA and hierarchical clustering on the correlation matrix both show every
patient's Control/Tamoxifen pair as nearest neighbors, well before any
other patient) -- exactly the structure that justifies the paired
`~ patient + treatment` design. No catastrophic sample-quality problem
was found; Tumor_04 (70-151 contributing cells) has visibly lower absolute
correlation with other samples, consistent with its low cell count, but
still pairs correctly with its own Control/Tamoxifen arm. No batch
correction was applied.

## Phase 5: edgeR differential expression

`~ patient + treatment` (Control reference level), `filterByExpr` ->
`calcNormFactors(method="TMM")` -> `estimateDisp` ->
`glmQLFit`/`glmQLFTest` on the treatment coefficient.

| Track | patient pairs | pseudobulk samples | genes before filter | genes after filter | residual df |
|---|---:|---:|---:|---:|---:|
| A (epithelial) | 10 | 20 | 33,538 | 17,987 | 9 |
| B (malignant, exploratory) | 3 | 6 | 33,538 | 13,863 | 2 |

## Phase 6: candidate extraction

BH correction applied **only to the candidates that survived
`filterByExpr`** in each track (12 of 13 -- USP17L29 was filtered out as
too lowly expressed in both tracks, correctly excluded from the BH family
rather than treated as non-significant). PAICS is never in that family.
**No candidate reaches candidate-set BH FDR < 0.05 in either track.**
Nominal (uncorrected) p<0.05 in Track B only: VEZF1 (p=0.0215), all 3
tested patients decreased.

## Phase 8: malignant vs non-malignant epithelial context

Descriptive: CNV score distribution differs by construction (malignant
cells were classified using this score -- not independent evidence).
Existing UMAP shows partial but real spatial separation (visible
malignant-enriched clusters, with expected intermixing given no batch
integration was applied). Cell-cycle scores (Seurat's standard
`CellCycleScoring`, field-standard gene sets) and the existing epithelial-
program module score were compared descriptively; dispersion metrics
(detected genes, %MT) are reported without a "messier = cancer" claim.

Secondary candidate table (5 patients with >=50 pooled malignant cells,
paired Wilcoxon signed-rank across patients, malignant vs non-malignant,
both treatment arms pooled): no candidate survives BH-13 correction
(minimum achievable p at n=5 is 0.0625; several hit that floor).

## Phase 9: integration

One row per candidate, `crispr_*` / `bulk_*` / `sc_track_a_*` /
`sc_track_b_*` / `malignant_context_*` columns, all raw effect
sizes/p-values/FDRs preserved, no composite score. Two descriptive
agreement flags (never evidence on their own): Track A/Track B direction
agrees for 9/12 tested candidates (disagrees: HMGB1, TADA2B, TSR3).

## Codex independent review

**Verdict: PASS WITH NOTES.** No fatal implementation bug found in
pseudobulk aggregation, patient-level replication, the edgeR design/
contrast, filtering/normalization, or the candidate BH-scope logic.
Corrections applied:

1. **Fixed (real bug):** the two integration "agreement" flags silently
   treated an untested gene as `False` (via `.fillna(False)`), which
   reads as "checked and disagrees" rather than "could not be checked."
   Changed to propagate `pd.NA` when either layer is untested --
   USP17L29 now correctly shows `<NA>`, not `False`, for both flags.
2. **Fixed (reproducibility):** the malignant-vs-nonmalignant paired
   Wilcoxon test now pins `zero_method="wilcox"`,
   `alternative="two-sided"`, `method="auto"` explicitly (previously
   relying on scipy's defaults, which are not guaranteed stable across
   versions).
3. **Fixed (defensive):** added an explicit `stopifnot(!anyDuplicated(rownames(counts)))`
   before pseudobulk row-summing in both the Track A/B and the
   malignant-vs-nonmalignant R scripts (Seurat already guarantees this;
   now verified, not assumed).
4. **Wording only (applied throughout this document and the code
   docstrings, not a code change):** the candidate BH correction must be
   described as "BH over the candidates that passed expression
   filtering" (12 of 13), never a blanket "BH-13"; Track A must never be
   called malignant-specific; Track B's raw p-values are nominal only
   (none survive correction); the malignant-vs-nonmalignant comparison is
   descriptive context, not evidence of malignant-state causation, and
   pools both treatment arms per patient (a design choice, stated
   explicitly, not the frozen Track B rule); the two integration
   agreement flags are descriptive navigation aids only, never validation
   or replication evidence.

All four are reflected in the code/docs as of this commit.
