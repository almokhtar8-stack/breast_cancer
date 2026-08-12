# Cross-dataset genome-wide integration: Codex independent review

**Date:** 2026-08-12. Reviewer: Codex (GPT-5.2-Codex via `mcp__codex__codex`),
read-only access to the full repository, explicitly instructed to try to
break the "no dataset gets extra weight" and "no gene name influences the
ranking" claims rather than confirm them.

**Initial verdict: FAIL.** Three genuine implementation bugs were found
that materially affected the coverage-tier and leave-one-out stability
results (though not the reported global Top 20's composition). All three
were fixed, the entire downstream pipeline was rerun, and the fixes were
independently re-verified. See "Corrections applied" below.

## What Codex found (material issues)

1. **GSE240112 epithelial-only genes received a phantom coverage vote.**
   `build_wide_matrix` (`src/cross_dataset_evidence_tables.py`) set
   `gse240112_testable = tumor_tested OR epithelial_tested`, but
   `compute_dataset_percentiles` (`src/cross_dataset_ranking.py`) uses
   **only** the tumor-cell track for the dataset's percentile (correctly,
   per the frozen Phase 7 rule). A gene tested only in the epithelial
   sensitivity track therefore got `gse240112_testable=True` with **no**
   actual percentile -- silently inflating its coverage-tier credit even
   though GSE240112 contributed nothing to its ranking. Codex measured
   this affecting 1,841 genes, with 628 entering the primary ranking who
   should not have, 302 of those meeting the >=3-dataset eligibility
   threshold only because of the phantom vote, and 67 wrongly placed in
   coverage Tier A. The bug propagated into the human-only and RNA-only
   rankings too (both also gate on `gse240112_testable`).
2. **Leave-one-dataset-out reruns did not replicate the main ranking
   hierarchy.** `_rank_subset`'s `"hierarchy"` scheme (`src/cross_dataset_stability.py`)
   sorted by `n_fdr05 -> n_top10 -> median -> mean -> gene`, omitting the
   coverage/testable-count key and the top-20%-count key that the main
   `build_global_ranking` hierarchy uses. Codex's independent
   reconstruction using the true main hierarchy changed most numerical
   leave-one-out ranks and flipped 3 stability labels in the Top 20
   (MSMO1, GABBR2, SLC4A10).
3. **`CONTEXT_DEPENDENT` category silently dropped its documented
   "or top-10%" branch.** The docstring said ">=2 datasets at FDR<0.05 OR
   top-10%," but `assign_evidence_category` (`src/cross_dataset_consensus_views.py`)
   checked only the FDR<0.05 count. Zero genes reached `CONTEXT_DEPENDENT`
   in the original output as a direct result; Codex identified at least
   9 genes (e.g. CTSZ, FOXD3, PAH, TFF2) that should have qualified via
   the missing branch.

## Smaller issues and cautions (also addressed)

- `compute_within_dataset_percentile`'s tie-break relied on the caller's
  input row order rather than an explicit key -- harmless in practice
  (the gene universe is built in alphabetical order upstream, so ties
  already broke alphabetically), but not self-evidently robust. Fixed by
  adding an explicit `gene` parameter used as a final ascending tie-break
  key.
- The word "unbiased" overstated what the analysis demonstrates (it is
  candidate-list-independent, not free of dataset-power or coverage
  effects on the ranking). Softened to "genome-wide,
  candidate-list-independent" in module docstrings.
- GSE245601's Track A/Track B mean-of-two-percentiles vs. single-track
  genes creates a real (documented, not hidden) missing-data asymmetry --
  noted as a limitation, not a bug requiring a code change.

## Corrections applied

1. `src/cross_dataset_evidence_tables.py`: `gse240112_testable` now
   reflects the tumor-cell track alone.
2. `src/cross_dataset_stability.py`: `_rank_subset`'s `"hierarchy"`
   scheme now uses the exact same 7-key sort as `build_global_ranking`
   (testable-count, FDR<0.05 count, top-10% count, top-20% count, median
   percentile, mean percentile, gene) -- verified by a new test asserting
   that applying it to all 5 datasets exactly reproduces the main global
   ranking's ranks.
3. `src/cross_dataset_consensus_views.py`: `CONTEXT_DEPENDENT` now checks
   `(n_datasets_fdr05 >= 2) OR (n_datasets_top10pct >= 2)`, matching the
   docstring; regression test added.
4. `src/cross_dataset_ranking.py`: `compute_within_dataset_percentile`
   gained an explicit `gene` tie-break parameter, now used everywhere it
   is called from `compute_dataset_percentiles`.
5. Wording: "unbiased" softened to "genome-wide,
   candidate-list-independent" in `src/cross_dataset_top20_lists.py` and
   `src/cross_dataset_surprise_discovery.py`.

**Effect of the corrections on results:** eligible genes in the primary
global ranking dropped from 15,557 to 15,255 (the 302 genes Codex
predicted would lose eligibility). `CONTEXT_DEPENDENT` now has 9 genes
(was 0). **The global Top 20's composition and order are unchanged** --
independently re-verified after the fix. The anonymization audit still
shows an exact 15,255/15,255 rank match (zero mismatches) on the
corrected pipeline. Ranking-stability numbers for the Top 20 changed
numerically (as expected, since the leave-one-out hierarchy fix changes
what "best/worst/median rank" means), though the qualitative picture
(most Top-20 genes MODERATELY_STABLE, a few DATASET_DEPENDENT, PLOD2 the
sole ROBUST gene) is similar to before.

## What passed review unchanged

- `build_full_gene_universe` is a genuine union with no candidate-list
  gate.
- Ambiguous duplicate symbols are excluded wholesale, never arbitrarily
  selected.
- The long-form evidence table has exactly 188,155 rows (37,631 genes x
  5 datasets, never 7).
- GSE118713 uses only the `TAMR_vs_MCF7` contrast.
- CRISPR functional direction and RNA up/down direction are never
  conflated (separate vocabularies, separate columns).
- The global ranking is a transparent lexicographic sort -- no hidden
  weighted score anywhere in `src/cross_dataset_*.py` (verified by Codex
  independently, and by a Phase 25 grep before the review).
- Resistance consensus uses only GSE118713, GSE240112 (tumor-cell), and
  GSE111151; never GSE245601; never requires CRISPR significance.
- CRISPR-independent discovery contains zero genes with CRISPR FDR<0.10.
- RNA-only excludes CRISPR entirely; human-only uses only GSE245601 and
  GSE240112.
- The anonymized-gene-ID mapping is genuinely non-alphabetical, and the
  anonymized ranking is rebuilt from scratch (percentiles recomputed on
  anonymized IDs, not just a relabeled copy of the named ranking).
- Figures spot-checked against their source tables were accurate.
- No causal or overclaiming language was found in the methods/audit docs
  beyond the "unbiased" wording already addressed.
- `python3 -m pytest tests/ -k cross_dataset -q`: 105 passed before the
  fixes (109 after, with the 4 new regression tests added for this
  review); full project suite 681 passed after the fixes, no regressions
  to any frozen prior-phase result.

## Explicit answers (from the review)

**"If every gene name were hidden and every dataset label were replaced
with Dataset A-E, would the algorithm still return the same ranking?"**
Yes for gene names (demonstrated empirically, zero mismatches under a
non-alphabetical shuffle). Dataset labels are structurally tied to each
dataset's own track-collapsing rule (e.g. GSE245601's two-track mean vs.
GSE240112's tumor-cell-only rule), so they are not interchangeable
symbols in the way gene names are -- this is an intentional design
choice (each dataset's biology is preserved, not genericized), not a
bias.

**"Does any dataset effectively receive more weight than another
because of duplicated analyses or hidden ranking logic?"** No dataset
receives two formal percentile/count votes. The GSE240112 phantom
-coverage bug (corrected above) *did* give it undue influence on
coverage-tier placement for genes it had no actual percentile for; that
is now fixed. GSE245601's single-track-vs-two-track missingness pattern
remains a documented, non-hidden statistical asymmetry (fewer genes
average two percentiles than use one), not a weighting bug.
