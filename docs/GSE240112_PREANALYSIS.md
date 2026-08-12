# Pre-analysis plan: GSE240112 tumor-cell primary-vs-recurrent candidate analysis

**Date written:** 2026-08-12
**Status:** locked before any candidate-gene expression value is inspected
in GSE240112. Written after `docs/GSE240112_DATA_AUDIT.md` (data-availability
facts only, no candidate values) and before any pseudobulk aggregation,
differential expression, or candidate-gene inspection. Amendments after
that point must be appended, dated, below -- never silently edited in
place.

## A. Primary contrast

Recurrent tumor cells (RT) vs primary tumor cells (PT), in the tumor/
malignant-cell compartment. Direction convention (fixed here, restated in
every downstream table): **positive log2FC = higher in RT than PT**.
Consistent with the project's cross-dataset framing, results are reported
as "higher/lower in recurrent tumors" or "associated with recurrent/
tamoxifen-treated disease," never as "tamoxifen caused X to change" --
this is a cross-sectional, unpaired comparison and cannot support a causal
claim.

## B. Unit of replication

Patient/sample (one pseudobulk profile per PT/RT sample; n=3 vs n=3).
Individual cells are never treated as biological replicates in any
inferential (p-value/FDR) step. Cell-level plots (Phase 12/13's per-cell
UMAP and dot-plot panels) are descriptive only.

## C. Statistical design

PT (OriGene Technologies Inc.) and RT (Ontario Tumor Bank) tissues come
from different source institutions and no pairing statement exists
anywhere in the GEO metadata or paper Methods (`GSE240112_DATA_AUDIT.md`
section 2). The design is therefore **unpaired**: model formula `~ group`
(group in {PT, RT}), contrast `RT - PT`. No patient/sample is treated as
matched to another; sample-level plots (Phase 11) never connect PT and RT
points with a line implying pairing.

## D. Candidate family (frozen, verified against the repository's existing
frozen table before writing this document)

The 13 candidates are the `sensitising_knockout`-direction rows of the
existing, already-committed `results/tables/candidate_evidence_summary.tsv`
(negative CRISPR `effect_size` = sensitising direction), cross-checked
against the user-supplied verbatim list and found to match exactly:

USP34, CTDNEP1, EIF4ENIF1, HMGB1, KDM1A, PET117, TADA2B, VEZF1, ICK,
SUPT4H1, TLK2, TSR3, USP17L29

`PAICS` is a separate published benchmark gene, reported alongside but
never included in this family's multiple-testing correction.

This set is frozen as of this document and will not be changed after
inspecting any GSE240112 expression or DE value.

## E. Multiple testing

Benjamini-Hochberg FDR is applied to the candidate family only (not
genome-wide), restricted to whichever of the 13 candidates are
statistically testable in the tumor-cell pseudobulk DE (pass the
pre-specified expression filter in section H). Any candidate that fails
the filter is retained in every output table with FDR = NA and an
explicit reason string ("filtered_out_low_expression" or
"gene_not_detected" etc.) -- never silently dropped from a table.
`PAICS` is reported with its own nominal p-value only, excluded from the
13-gene BH family.

Given the small sample size (n=3 vs n=3), a genome-wide FDR value is also
reported for context (standard edgeR output) but is not used for any
candidate-level claim.

## F. Primary expression analysis

Sample-level pseudobulk of raw UMI counts, summed within each PT/RT
sample over the author-defined tumor/malignant-cell population (the
`TTs_cancer_060223.h5seurat` object, in which `cell.annot ==
"Breast cancer cells"` for all 9,942 cells -- see Phase 7 tumor-cell
definition below). One profile per sample (6 profiles: PT1, PT2, PT3,
RT1, RT2, RT3). Differential expression via edgeR (matching this
project's established convention from GSE118713/CRISPR analyses),
TMM normalization, `filterByExpr`-style low-count filtering applied
genome-wide before fitting, quasi-likelihood F-test (`glmQLFit`/
`glmQLFTest`) for the `RT - PT` contrast.

## G. Tumor-cell definition (Phase 7 decision, frozen here)

**Case A applies**: author-provided malignant/tumor-cell labels exist
publicly. The `TTs_cancer_060223.h5seurat` object is, by construction,
already restricted to the authors' own cancer-cell calls (`cell.annot`
reads `"Breast cancer cells"` for every one of its 9,942 cells; per the
paper Methods this was generated via `Seurat::FindAllMarkers` against
curated reference marker panels -- PanglaoDB, CellMarker, ProteinAtlas,
literature -- explicitly not a CNV-based method). This is used directly
as the PRIMARY tumor-cell population, with no independent CNV-based
re-classification, no threshold tuning, and no comparison of alternative
malignancy-calling methods (per the user's explicit time-boxing
instruction and the prior GSE245601 inferCNV/CopyKAT experience, which is
not to be repeated).

RT3 has only 178 cancer cells in this object versus 1,029-2,721 for the
other five samples -- an order of magnitude fewer. This is flagged now,
before any DE result is inspected, as a sample to inspect carefully at
the pseudobulk-QC stage (Phase 9). It will not be excluded pre-emptively;
if pseudobulk QC (library size, PCA/MDS placement, correlation with other
RT samples) shows it behaves as a normal (if small) biological replicate,
it is kept. It will only be dropped if it shows a catastrophic technical
signature, and if that happens the reason is reported explicitly before
any candidate result is shown, per the user's Rule 2 (no threshold/QC
changes made after seeing candidate effects).

**Sensitivity population (Phase 14)**: all epithelial cells (not
malignant-restricted), reconstructed from the raw Cell Ranger PT/RT
matrices (`data/raw/gse240112/cellranger/`) via standard canonical
epithelial markers (EPCAM, KRT8, KRT18, KRT19) versus immune (PTPRC),
stromal (PDGFRB, COL1A1), and endothelial (PECAM1) markers -- broad
compartment calling only, not a malignancy classifier, and not used to
retroactively alter the Case A tumor-cell labels above.

## H. Significance

FDR < 0.05 (candidate-family BH, section E) is the only threshold treated
as statistically supported. Nominal (uncorrected) p-values are shown
alongside FDR in every candidate table but are always explicitly labeled
"nominal" and never described as significant on their own. A
candidate-set-testable gene is defined as one that survives edgeR's
standard `filterByExpr` low-count filter applied to the tumor-cell
pseudobulk count matrix; genes that do not pass this filter are reported
as NA/untestable (section E), not silently excluded.

## I. Freeze condition

No threshold in sections E-H will be changed after any candidate's
GSE240112 pseudobulk expression, log2FC, p-value, or FDR is inspected.
Any genuinely necessary correction discovered later (e.g. at Codex review,
Phase 23) will be appended as a dated amendment below, not edited in
place, matching the convention already used in the root `PREANALYSIS.md`
and `docs/gse245601_PREANALYSIS.md`.

---

*(No amendments yet.)*
