# GSE111151 Phase 14: Codex independent review

**Date:** 2026-08-12. Reviewer: Codex (GPT-5.2-Codex via `mcp__codex__codex`),
read-only access to the full repository, instructed to be skeptical and
audit for genuine errors rather than confirm the existing conclusions.

**Verdict: PASS WITH NOTES.** The core statistical/computational
implementation was independently re-verified and found sound; three
report-wording issues were found and corrected.

## What Codex found

1. **Factual contradiction about RNA significance.** The report claimed
   VEZF1 was "the only candidate with any FDR<0.05 RNA result in this
   entire project" -- false, since USP34 has a candidate-set-significant
   GSE118713 bulk result (FDR=0.0073), stated elsewhere in the same
   report. Corrected to clarify which FDR family/layer each claim refers
   to (section 8).
2. **Overly strong prioritization language.** "Should be dropped from
   further prioritization" was stronger than the underpowered,
   heterogeneous RNA validation evidence justifies, and was in tension
   with the report's own statement (carried from the GSE240112 run's
   Absolute Rule 8) that nonsignificance does not refute prior CRISPR
   evidence. Corrected to "not prioritized by current RNA evidence"
   (section 8).
3. **Broken forward references** to this file (`docs/GSE111151_CODEX_REVIEW.md`)
   before it existed. Resolved by writing this file.

## What passed review (no changes needed)

- `~ cell_line + resistance_status` is the appropriate primary design (a
  naive 7-vs-4 comparison would be confounded by cell-line identity);
  the full-rank check (`qr(design)$rank == ncol(design)`,
  `scripts/analysis/gse111151_02_edger.R`) is real, not vacuous --
  confirmed 5 design columns, full rank, 6 residual degrees of freedom.
  Parental is the reference level; the tested coefficient
  (`resistance_statusresistant`) is uniquely matched.
- The pre-existing GSE111151 config block (committed 2026-08-05,
  `a7f3bbf1`) was correctly preserved: all 11 GSM/file/cell-line/status
  mappings match the data audit table; this session's additions were
  purely additive, no existing field was changed.
- The Ensembl-gene-ID row-order identity check
  (`scripts/analysis/gse111151_01_build_count_matrix.R`) uses vector-wise
  `identical()` against the first file and would catch a reordered,
  missing, or substituted gene ID; gene-name order is checked too.
- edgeR's `filterByExpr` -> TMM -> `estimateDisp` -> `glmQLFit`/
  `glmQLFTest` workflow is standard and appropriate.
- **The TMM normalization fix (caught during this run, before Codex
  review) was independently verified to be genuinely and consistently
  propagated**: the exported effective library size
  (`lib.size * norm.factor`) is used in every downstream sample-level
  calculation (candidate extraction, candidate heatmap, USP34/VEZF1/
  SUPT4H1 sample-level plots) -- confirmed by independently recomputing
  every saved candidate/PAICS sample-level value from the raw counts and
  TMM factors, matching to numerical precision (max absolute error
  8.9e-16).
- USP17L29 confirmed to have exactly zero counts in all 11 raw samples.
- Candidate lookup is by frozen Ensembl ID (not gene symbol, which has
  308 duplicates in this reference), confirmed in code.
- Candidate-set BH was independently reproduced over exactly the 12
  testable candidates, PAICS excluded; zero discrepancy from the
  committed table.
- Reported sample-level cell-line contrasts were independently
  recomputed and matched exactly: USP34 (+0.296, +0.211, -0.036, +0.200
  -> 3/4 consistent), VEZF1 (-0.605, -0.705, +0.075, +0.264 -> 2/4
  consistent), SUPT4H1 (-0.661, +0.329, +0.088, +0.129 -> 3/4 by sign but
  weak/mixed, p=0.916).
- The tightened classification rule (`p < 0.3 AND n_consistent >= 3`,
  `src/gse111151_evidence_classification.py`) is correctly implemented
  as an AND, with the separate `p < 0.05` route remaining a valid OR
  branch; the real output contains exactly USP34 and ICK in the
  `directionally_supportive_but_weak` tier, matching the report.
- "USP34 modestly strengthened" and "VEZF1 did not reproduce" were found
  to be defensible, appropriately-qualified interpretations.
- The 5-layer integration table performs correct gene-keyed joins,
  creates no composite score, and only *reads* the frozen upstream
  tables (CRISPR/GSE118713/GSE245601/GSE240112) -- confirmed it writes
  solely the new GSE111151 integration output.
- No causal ("tamoxifen caused...") language was found.
- `python3 -m pytest tests/ -k gse111151 -q` reproduced: 47 passed. Spot
  -checked tests contain real numerical normalization, Ensembl-lookup,
  BH-family, and classification-rule-regression checks, not only smoke
  tests.

## Corrections applied

All in `docs/GSE111151_ANALYSIS_REPORT.md` section 8 unless noted:
softened the "only candidate with FDR<0.05" claim to name both USP34
(GSE118713 layer) and VEZF1 (GSE240112 layer) with their respective
scopes; replaced "should be dropped" with "not prioritized by current RNA
evidence" plus an explicit non-finality caveat; created this file to
resolve the forward references.

No source data, extraction, or statistical computation was rerun -- all
three corrections were report-wording issues, not computational errors.
