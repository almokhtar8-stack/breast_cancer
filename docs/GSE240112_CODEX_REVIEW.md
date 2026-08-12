# GSE240112 Phase 23: Codex independent review

**Date:** 2026-08-12. Reviewer: Codex (GPT-5.2-Codex via `mcp__codex__codex`),
read-only access to the full repository, instructed to be skeptical and
audit for genuine errors rather than confirm the existing conclusions.

**Initial verdict: FAIL.** Corrections were applied (below); the
technical/statistical implementation itself was found sound.

## What Codex found

1. **USP34 individual-sample-consistency claim was wrong.** The original
   report claimed "2 of 3 RT samples clearly exceed all PT samples."
   Verified directly: PT range is 6.760-7.088 (log2 CPM+1); RT1=7.060 and
   RT2=6.858 both fall *inside* that range; only RT3 (7.656) exceeds it.
   RT3 is also the pseudobulk sample built from the fewest cells (178 vs.
   1,029-2,721) and the smallest library (1.14M vs. 43-65M UMIs). The
   positive RT-vs-PT effect is therefore substantially driven by the
   smallest, noisiest sample -- the opposite of what "not driven by a
   single outlier" claimed.
2. **The tumor-cell and all-epithelial tracks were described as "two
   independent human single-cell populations."** They are not: the
   all-epithelial track reprocesses the same six sequencing libraries
   with a broader (overlapping) cell-selection pipeline. It is a
   sensitivity check on the cell-population definition, not independent
   replication, and should not be counted as a second independent
   dataset when assessing convergent evidence.
3. **"No dataset shows USP34 moving in the opposite direction" was
   factually false.** GSE245601 Track A (log2FC=-0.033) and Track B
   (log2FC=-0.181) are both negative -- opposite in sign to GSE240112's
   +0.400 -- even though neither is itself statistically significant.
   The defensible statement is that no dataset shows a *statistically
   supported* opposite-direction effect, not that no dataset shows the
   opposite sign at all.
4. **"MODESTLY STRENGTHENS" was not well supported once (1)-(3) are
   corrected.** Given the RT3-driven fragility, the non-independence of
   the sensitivity track, and GSE245601's (nonsignificant) opposite sign,
   the original upgrade from "no new significant evidence" to "modestly
   strengthens" relied on claims that did not hold up.
5. **Patient-unrelatedness was overstated as proven.** "Different,
   unrelated patients" and "no single patient contributes both" go
   beyond what the evidence (different tissue-source institutions, no
   pairing statement in the retrieved metadata, Table S1 not retrieved)
   actually establishes. The unpaired design decision itself is still
   correct -- it should rest on "no evidence of pairing was found," not
   on a stronger, unproven claim of patient unrelatedness.
6. **The evidence-classification scheme (`src/gse240112_evidence_classification.py`)
   has a conceptual gap**: for genes without a significant prior
   GSE118713 bulk direction (most of them, including VEZF1), the
   "directionally supportive" tier has no pre-specified direction to
   check consistency against -- it is measuring within-GSE240112 evidence
   strength only, not cross-dataset concordance, which should be stated
   explicitly rather than left implicit.

## What passed review (no changes needed)

- `sample_map` in `scripts/analysis/gse240112_04_cellranger_epithelial.R`
  exactly matches the audit's GSM/raw-prefix mapping.
- The direct-HDF5 CSC sparse-matrix reconstruction in
  `scripts/analysis/gse240112_01_extract_h5seurat.R` (zero-based `i`,
  column-pointer `p`, values `x`, declared `dims`) is correct.
- Tumor-cell pseudobulk aggregation (`scripts/analysis/gse240112_02_build_pseudobulk.R`)
  uses disjoint raw-UMI sums, exactly one profile per sample, with
  explicit no-loss/no-duplication checks.
- edgeR design and direction (`scripts/analysis/gse240112_03_pseudobulk_edger.R`):
  PT reference level, `~group`, TMM, `filterByExpr`, QL fit/test, positive
  logFC = higher in RT -- all correct.
- Candidate-set BH correction was independently reproduced to numerical
  precision across the 12 testable candidates; PAICS correctly excluded.
- USP17L29's absence was independently re-verified in all 6 raw Cell
  Ranger matrices: `ENSG00000231637`/`USP17L29` is present in the
  reference at feature row 8023 but has exactly zero nonzero counts in
  every cell of every sample -- the "genuinely undetected, not a
  symbol-mapping error" claim holds.
- Treating `TTs_cancer_060223.h5seurat` (`cell.annot` uniform across all
  cells) as the author-defined tumor-cell population is a reasonable
  reading of a marker-based (not CNV-based) author call, correctly
  disclosed as such rather than independently validated.
- The 4-layer integration table has exactly 13 unique genes, correct
  direct joins, and no invented composite score.
- Spot-checked figures (PCA, effect-size, USP34 per-sample distribution,
  broad-compartment UMAP) correctly reflect their underlying tables.
- No inappropriate causal ("tamoxifen caused...") wording was found.
- `python3 -m pytest tests/ -k gse240112 -q` reproduced: 48 passed.

## Corrections applied

All in `docs/GSE240112_ANALYSIS_REPORT.md` unless noted:

- Section 5 (USP34): replaced the incorrect "2/3 RT samples exceed PT
  range, not outlier-driven" claim with the corrected, RT3-driven
  description; reframed the all-epithelial result as a non-independent
  sensitivity check.
- Section 6/7: added an explicit caveat that the classification scheme
  is descriptive/post-hoc and, for most genes, measures only
  within-GSE240112 evidence strength (item 6 above).
- Section 8: left factually as-is (it already reported the GSE245601
  negative values correctly; only section 9's summary sentence was
  wrong).
- Section 9: **USP34 final status changed from MODESTLY STRENGTHENS to
  UNCHANGED-NEUTRAL**, with a rewritten justification that accurately
  states GSE245601's opposite sign, the RT3-driven fragility, and the
  non-independence of the sensitivity track.
- Section 11: added limitations entries for the RT3-driven USP34
  fragility and the non-independence of the epithelial sensitivity
  track; softened the patient-pairing limitation to match the corrected
  epistemic status.
- `docs/GSE240112_DATA_AUDIT.md` section 2: softened "different,
  unrelated patients" / "no single patient contributes both" to "no
  evidence of pairing was found," matching what the retrieved sources
  actually establish.

No source data, extraction, pseudobulk, or DE computation was rerun --
none of the corrections involved a computational error; they were all
interpretive/wording errors in the report layer. No test needed to be
added or changed as a result (the existing 48 GSE240112 tests already
passed and cover the correct underlying computations; the bug was in
report prose, not in code).
