# Pre-analysis plan: GSE111151 independent tamoxifen-resistance validation

**Date written:** 2026-08-12. Locked before any candidate-gene expression
value is inspected in GSE111151. Written after `docs/GSE111151_DATA_AUDIT.md`
(data-availability facts only, no candidate values). Amendments after this
point must be appended, dated, below -- never silently edited in place.

## A. Primary contrast

Tamoxifen-resistant sublines vs. their isogenic parental cell line,
**within cell-line blocks**. Direction convention (fixed here, restated
in every downstream table): **positive log2FC = higher in
tamoxifen-resistant cells**. This is a genuine tamoxifen-resistance
contrast (confirmed from primary sources, not assumed -- see data audit
section 2), so "tamoxifen-resistance-associated" language is used; this
remains a cell-line-model, cross-sectional comparison, and does not by
itself establish that any single gene's expression change causes
resistance.

## B. Sample inclusion

All 11 samples (4 parental, 7 resistant) are included. No sample is
excluded: library sizes range 30.9M-75.0M raw counts with no catastrophic
outlier (data audit section 5, item 12). If pseudobulk/sample-level QC
(Phase 4) reveals a genuine catastrophic technical problem, it will be
reported and any exclusion decision made *before* any candidate value is
inspected, exactly as for GSE240112.

## C. Unit of biological replication and paired/unpaired design

**Cell-line-blocked design, not a naive unpaired two-group comparison.**
Each of the 4 cell lines (MCF-7, T-47D, ZR-75-1, BT-474) has exactly one
parental sample and 1-2 independently-derived resistant sublines. A naive
unpaired two-group test (all 7 resistant vs. all 4 parental, ignoring
cell-line identity) would confound cell-line-intrinsic expression
differences (these are molecularly distinct lines -- BT-474 is
HER2-amplified, the others are not) with the resistance effect, and is
explicitly rejected here. Model: `~ cell_line + resistance_status`
(cell_line as a 4-level blocking factor, analogous in spirit to this
project's `~ patient + treatment` GSE245601 design, but here the "patient"
is a cell-line background rather than an individual, and the design is
unbalanced -- 1 parental vs. 1 or 2 resistant per block, not a strict 1:1
pair). Each resistant subline is compared to its own cell line's single
parental sample; Tam1 and Tam2 (where both exist) are treated as two
independent resistant-status observations within the same block, not
averaged or dropped.

## D. Candidate family (frozen, read from the repository's existing
frozen config, not manually retyped)

Read from `config/config.yaml` → `gse245601_pseudobulk.candidates.thirteen`
(cross-checked identical to `gse240112.candidates.thirteen`, both already
frozen and committed): USP34, CTDNEP1, EIF4ENIF1, HMGB1, KDM1A, PET117,
TADA2B, VEZF1, ICK, SUPT4H1, TLK2, TSR3, USP17L29. `PAICS` is the
separate published benchmark, reported alongside but never included in
this family's multiple-testing correction. This set is not redefined,
reordered, or reduced after inspecting any GSE111151 value.

All 14 genes (13 candidates + PAICS) are confirmed present by exact
Ensembl-gene-ID match in the GSE111151 expression matrix (verified in the
data audit, before any count value was inspected) -- this is a
presence-in-reference check only, not a detection/expression check, which
is deferred to Phase 6.

## E. Normalization approach

Raw integer counts (the `counts` column in each per-sample file) are used
as input to an independently-run edgeR pipeline (TMM normalization within
edgeR), matching this project's established convention
(GSE118713/GSE245601/GSE240112 all use edgeR on raw counts). The
publisher-supplied `CPM_batch (log2)` column is **not** used as the
primary input: per the data audit (section 4), that batch correction
targets a "study origin" effect from merging this cell-line panel with an
external patient cohort in the original paper -- a correction that is not
relevant to, and would not be appropriately scoped for, an
independently-specified within-GSE111151 model.

## F. Filtering rule

Standard `edgeR::filterByExpr()` applied to the full 60,619-gene raw
count matrix using the `~ cell_line + resistance_status` design matrix,
before any candidate gene is inspected.

## G. Statistical model

`edgeR`: `DGEList` → `filterByExpr` → `calcNormFactors` (TMM) →
`estimateDisp` → `glmQLFit` → `glmQLFTest` on the `resistance_status`
(resistant vs. parental) coefficient, with `cell_line` as an additional
blocking term in the design matrix. Genome-wide log2FC, p-value, and
genome-wide BH FDR are reported for context; the candidate-level claim
never relies on genome-wide FDR alone (section H).

## H. Candidate multiple-testing rule

Benjamini-Hochberg FDR restricted to whichever of the 13 candidates pass
`filterByExpr` in this model ("candidate-set BH"), computed on exactly
that testable subset -- mirroring the established convention in
`src/gse245601_candidate_extraction.py` and
`src/gse240112_candidate_extraction.py`. Any candidate that fails the
filter (or is absent from the count matrix, which section D already
rules out) is retained in every output table with FDR = NA and an
explicit reason string, never silently dropped. `PAICS` is reported with
its own nominal p-value/genome-wide FDR only, excluded from the 13-gene
BH family.

## I. PAICS handling

Reported exactly like a candidate row in every output table but flagged
`is_paics_benchmark = True` and never included in the candidate-set BH
family, matching the convention already used throughout this project.

## J. Planned QC

Library size, detected-gene count, sample-sample correlation, and PCA on
log2(CPM+1) of the raw counts, labeled by `cell_line` and
`resistance_status`. Strong clustering by cell line is expected and is
not itself a QC failure or an "outlier" to remove -- it is exactly why
the model in section C blocks on cell line rather than ignoring it. No
sample is excluded based on how it affects any candidate gene's apparent
result.

## K. Planned figures

PCA/MDS by group, sample correlation heatmap, candidate effect-size plot
(13 genes), candidate expression heatmap across all 11 samples, USP34/
VEZF1/SUPT4H1 (and any other candidate that stands out) sample-level
plots grouped by cell line, optional volcano plot with the 13 candidates
highlighted.

## L. Interpretation limits (frozen before results are seen)

- A significant result in this panel does not prove causality for
  tamoxifen resistance in patients; these are long-term-cultured cell
  line derivatives, and cross-resistance to other stresses of prolonged
  culture cannot be excluded.
- The panel spans more than one molecular subtype (BT-474 is
  HER2-amplified); a candidate's effect that is consistent across all 4
  backgrounds is stronger evidence than one driven by a single cell line,
  and this must be checked explicitly at the sample level (Phase 7)
  before any "independently supported" claim is made.
- With only 4 blocks and 7 resistant vs. 4 parental observations, this
  panel has limited power; a nonsignificant candidate-set FDR does not
  refute prior CRISPR/GSE118713/GSE245601/GSE240112 evidence (Absolute
  Rule 8, carried over from the GSE240112 run), and a nominal p<0.05
  alone is not treated as "independently supported" (Phase 9's explicit
  instruction).
- This dataset is not used to redefine the candidate list, and no
  analysis choice (filtering, normalization, model form, threshold) is
  changed after any candidate value has been inspected.

---

*(No amendments yet.)*
