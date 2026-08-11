# InferCNV malignant-classification threshold diagnostics

**Status: audit-only, diagnostic-only. No InferCNV or CopyKAT rerun. No
label changed. No candidate gene read. No threshold optimized, selected,
or informed by CopyKAT agreement anywhere in this document.**

Scope: the 12 samples from the CNV method-comparison audit --
Tumor_01/02/03/04/08/10 x Control/Tamoxifen. Grouped per the request:

- **Disagreement cases** (InferCNV shows chromosome-scale abnormality on
  its own heatmap but calls almost no cells malignant): Tumor_01, Tumor_04,
  Tumor_08.
- **Good agreement controls**: Tumor_02, Tumor_03.
- **Signal-quality case, tracked separately**: Tumor_10.

All numbers below are reproducible from `results/tables/gse245601_infercnv_threshold_diagnostics.tsv`,
`results/tables/gse245601_infercnv_threshold_sensitivity_grid.tsv`, and
`results/tables/gse245601_infercnv_local_score_sensitivity.tsv`, all
written by `src/gse245601_infercnv_threshold_diagnostics.py`
(`run_threshold_diagnostics()`).

---

## 1. The rule, reconstructed exactly

Source: `scripts/analysis/gse245601_05_infercnv_malignant.R`, which itself
implements `docs/gse245601_PREANALYSIS.md` section 9 ("Malignant-cell
classification") and its 2026-08-11 amendment 3. Provenance tags below
follow that document's own convention (section 3): **[A]** = stated in the
published paper (PMID 37747807); **[B]** = found in the authors' code, not
(or not exactly) in the paper prose; **[A+B]** = both; **[C]** = neither
source fully specified it -- a decision this project made and froze.

| Step | Rule | Provenance |
|---|---|---|
| CNV score | `cnv_score = colMeans((expr_mat - 1)^2)`, per epithelial cell, across all genes in InferCNV's output | **[A+B]** -- paper: "CNV values"; code: `run_infercnv.R` |
| Grouping | Per (sample, Seurat cluster); clusters with <10 epithelial cells are pooled into one `whole_sample_pooled` group per sample | **[C]** -- amendment 3, frozen before InferCNV was run, because a tiny cluster gives an unstable mean/SD |
| Seed selection | Top `max(round(n_group * 0.05), 2)` cells by CNV score within the group | **[A+B]** -- paper: "top 5% of cells with high CNV values"; the minimum-2 floor is in the code |
| Correlation | Kendall tau (`pcaPP::cor.fk`, a fast, numerically-identical substitute for `stats::cor(method="kendall")`) between each cell's centered CNV profile and the seed group's mean profile | **[A+B]** -- paper: "correlation values"; code: exact formula |
| CNV-score threshold | `th_value = clamp(mean(score) - 2*SD(score), 0.01, 0.05)` | **[A+B]** for the mean-2SD form and the clamp bounds (both in the code); **[C]** for the decision to use this clamped variant (`"ng_2021_and_thresholding"`) rather than the simpler, unclamped `"ng_2021"` variant -- PREANALYSIS.md section 9 states this choice was not literally recoverable as a single default and was frozen by this project |
| Correlation threshold | `th_corr = clamp(mean(corr, na.rm=TRUE) - 1.5*SD(corr, na.rm=TRUE), 0.2, 0.4)` | **[A+B]**, same as above |
| Final call | `malignant` iff `score > th_value AND corr > th_corr` (strict `>`); NA on either side forces `non-malignant epithelial` | **[A+B]** |

**This diagnostic module's reconstruction is arithmetic, not a rerun.**
`gse245601_malignant_cell_labels.tsv` already stores each cell's
`cnv_score`, `cnv_correlation_to_seed`, `threshold_group`, and the
recorded `primary_malignancy_label` -- InferCNV itself is never touched
again. `recompute_group_thresholds()` reapplies the formula above to those
frozen per-cell values, and `verify_reconstruction()` checks that doing so
reproduces the frozen label for **every single one of the 11,211 cells**
in the 12 selected samples, exactly. It does (0 mismatches; also
cross-checked: the recomputed malignant count, summed per sample, exactly
equals the independently-written `gse245601_malignant_summary_per_sample.tsv`
`n_malignant` column for all 12 samples). One documented asymmetry: R's
correlation threshold uses `na.rm=TRUE` but its CNV-score threshold does
not; this module's pandas recomputation skips NaN for both, so the two are
only guaranteed equivalent when neither column has missing values --
`load_cnv_score_table()` asserts this (raises if not), rather than
assuming it.

---

## 2. Good vs. bad tumors: the numbers

One row per tumor arm, cells rolled up across all (sample, threshold_group)
groups (a sample may have several groups, each with its own adaptive
threshold -- see the full per-group table in the TSV for that detail).

| Tumor | Arm | n cells | malignant | fails CNV-score only | fails correlation only | fails both | malignant % |
|---|---|---:|---:|---:|---:|---:|---:|
| Tumor_01 | Control | 3334 | 1 | 3075 | 0 | 258 | 0.03% |
| Tumor_01 | Tamoxifen | 2774 | 2 | 2569 | 1 | 202 | 0.07% |
| Tumor_04 | Control | 151 | 0 | 90 | 0 | 61 | 0.00% |
| Tumor_04 | Tamoxifen | 70 | 1 | 48 | 0 | 21 | 1.43% |
| Tumor_08 | Control | 906 | 17 | 823 | 0 | 66 | 1.88% |
| Tumor_08 | Tamoxifen | 984 | 25 | 880 | 0 | 79 | 2.54% |
| Tumor_02 | Control | 355 | 134 | 112 | 0 | 109 | 37.75% |
| Tumor_02 | Tamoxifen | 325 | 183 | 54 | 3 | 85 | 56.31% |
| Tumor_03 | Control | 513 | 276 | 140 | 3 | 94 | 53.80% |
| Tumor_03 | Tamoxifen | 1105 | 582 | 283 | 2 | 238 | 52.67% |
| Tumor_10 | Control | 416 | 61 | 175 | 3 | 177 | 14.66% |
| Tumor_10 | Tamoxifen | 278 | 2 | 59 | 0 | 217 | 0.72% |

**CNV-score threshold: exactly 0.01 (the clamp floor) in every one of the
63 (sample, threshold_group) groups across these 12 samples** -- and,
checked separately, in all 115 groups across the entire 20-sample cohort.
The "adaptive" `mean(score) - 2*SD(score)` formula never once landed
inside its allowed `[0.01, 0.05]` range without being clamped up to the
floor; it behaves as a fixed 0.01 cutoff throughout this dataset, not as a
per-group-adaptive one. (Correlation threshold, by contrast, does vary
across its 0.2-0.4 range from group to group -- see the TSV.)

Median/range of CNV score and correlation, and seed diagnostics, are in
the full table; the headline pattern: Tumor_01/04/08's cells have a CNV
score distribution sitting almost entirely below 0.01 (sample medians
0.003-0.006; see the TSV's `median_cnv_score` column), while Tumor_02/03's
sample medians (0.008-0.011) straddle the floor, leaving roughly half the
population above it.

---

## 3. Diagnostic plots

Two figures per sample under `results/figures/gse245601_infercnv_threshold_diagnostics/`:

- **`{sample}_score_vs_correlation.png`** -- one panel per threshold_group
  actually used in that sample (a sample is not force-fit to one pair of
  thresholds; the real per-group adaptivity, such as it is, is shown
  faithfully). X = CNV score, Y = Kendall correlation to seed. The two
  **actual, recomputed** threshold lines for that group are drawn (dashed).
  Point color = the frozen InferCNV call (red = malignant, blue =
  non-malignant); point shape = CopyKAT's call (circle = aneuploid,
  triangle = diploid, x = not.defined) -- shown for visual context only,
  never used to place or move a threshold line.
- **`{sample}_sensitivity_grid.png`** -- Point 5's coarse sweep (see below),
  with the sample's actual per-group threshold(s) marked as red `+`.

One additional figure, **`local_score_sensitivity_all_samples.png`**,
built specifically to answer the "near or far from cutoff" question with
adequate resolution (see Point 4).

---

## 4. Which threshold, and how far?

**The CNV-score threshold is overwhelmingly the dominant, and for most
cells the sole, blocking criterion in Tumor_01/04/08.** In all three,
`fails_correlation_only` is 0 or 1 cell; nearly every non-malignant call
is either "fails CNV-score only" (the large majority) or "fails both"
(a non-trivial minority, especially in Tumor_04 where fails_both/fails_cnv_only
is roughly 2:3). Correlation-to-seed is essentially never the reason a
cell is excluded in these three tumors -- most cells that fail the
CNV-score bar would have passed the correlation bar on their own.

**How far below the cutoff, and is it a knife-edge right at 0.01?** Two
separate checks, because the coarse 0.005-spaced sweep in Point 5 cannot
resolve local behavior:

- Among cells that fail the CNV-score criterion, the median shortfall
  (`0.01 - score`) is **0.0070 (Tumor_01), 0.0070 (Tumor_04), 0.0039
  (Tumor_08)** -- 39-70% of the threshold's own value. This is not a
  population of cells sitting a hair's breadth below the line.
- The dense local sweep (`local_score_sensitivity_all_samples.png`, 17
  points from 0.006 to 0.014, i.e. roughly double the floor on either
  side) shows **three distinct shapes**, not one:
  - **Tumor_01 and Tumor_04**: already down to ~1-7% passing at the *left*
    edge of this band (0.006, i.e. 40% below the actual floor) and
    essentially flat at ~0% by 0.0095-0.01. The population isn't clustered
    at a cliff near 0.01 -- it has already thinned out well before reaching
    it. Moving the threshold slightly would not meaningfully change the
    result for these two; only a much larger relaxation would.
  - **Tumor_08**: a genuinely steep, gradual decline through this whole
    band (~50% at 0.006 down to ~2% at 0.01) -- this tumor's yield *is*
    locally sensitive to where the CNV-score line is drawn, more so than
    Tumor_01/04.
  - **Tumor_02/03 (good controls)**: gently sloping, staying at 40-70%
    across the entire band -- nowhere near a cliff.

**Seed-cell construction**: reconstructed independently (top-5%-by-score
per group, minimum 2 cells). This sharpens the picture rather than
softening it: in **Tumor_01, the seed's own mean CNV score is below the
0.01 floor in all 12 of its groups** (`seed_margin_above_cnv_threshold`
ranges -0.0070 to -0.0020); in **Tumor_04, the same is true in all 10 of
its groups** (-0.0068 to -0.0009). That is, even the most CNV-elevated 5%
of cells within every single cluster of these two tumors -- the cells the
rule itself picked as the best available "malignant-like" reference --
still do not clear the fixed floor. Tumor_08 is intermediate (4 of 6
groups negative, 2 positive, up to +0.0004). By contrast, in Tumor_02 and
Tumor_03 roughly half of groups have a *substantially* positive seed
margin (up to +0.021 and +0.010 respectively) -- some clusters have a
seed population clearly and comfortably above the floor, which Tumor_01
and Tumor_04 never do. The correlation-to-seed step is still not the
bottleneck (see above) -- but the seed-construction diagnostic shows the
CNV-score shortfall in Tumor_01/04 is not a property of "most" cells while
a distinct high-scoring subpopulation exists; no such subpopulation was
found in any cluster of either tumor.

**No evidence of a coding or implementation error was found.**
`verify_reconstruction` reproduces the frozen labels exactly for all 11,211
cells; the clamp and formula are being applied exactly as coded. The
finding is a **threshold-metric mismatch, not a bug**: a whole-genome mean
of squared deviations (`colMeans((expr-1)^2)` across ~8,000+ genes) is
diluted by however much of the genome is *not* affected in a given tumor,
and the hard 0.01 floor apparently sits above that diluted score for most
cells in three of these six tumors. This document does not determine
*why* Tumor_01/04/08's diluted scores are lower than Tumor_02/03's (signal
extent, InferCNV's internal scaling/centering/denoising, or something else
in reference-cell composition are all plausible and untested here) --
that is a mechanistic question, not a threshold-diagnostic one, and is
explicitly out of scope for this document.

---

## 5. Sensitivity sweep (diagnostic only)

`build_sensitivity_grid()` sweeps a **fixed** threshold pair (not
per-group-adaptive) over CNV score in `[0.0, 0.05]` (11 points) x
correlation in `[0.2, 0.4]` (9 points, the rule's own clamp bounds),
counting malignant cells pooled across each sample. This is never used to
choose, optimize, or recommend a threshold -- it exists only to show the
shape of the yield curve.

Findings, consistent with Point 4:

- Malignant fraction depends almost entirely on the CNV-score axis; moving
  the correlation threshold across its full 0.2-0.4 range changes the
  result by only a few percentage points at any fixed CNV-score value
  (the sensitivity-grid heatmaps are essentially horizontal bands).
- At `cnv_score_threshold=0.0` (correlation threshold varied across its
  full 0.2-0.4 range), malignant fraction is 41-60% for Tumor_04_Control
  and 81-98% for Tumor_08_Control; both fall sharply by 0.005 (16% and
  69-77% respectively) and are ~0-2% by 0.01 -- concentrated, substantial
  sensitivity in the 0.0-0.01 band for both, though starting from very
  different levels.
- For Tumor_02_Control, the same sweep still shows ~35-40% malignant at
  the actual 0.01 floor -- it does not collapse.
- The finer, local sweep in Point 4 refines this further: the *rate* of
  decline right around 0.01 differs by tumor (steep for Tumor_08, already
  flat-near-zero for Tumor_01/04), so "how sensitive" is not a single
  cohort-wide answer.

---

## Files

- `results/tables/gse245601_infercnv_threshold_diagnostics.tsv` -- one row
  per (sample, threshold_group): recomputed thresholds, score/correlation
  distribution, 4-way failure breakdown, seed diagnostics, gap-to-cutoff.
- `results/tables/gse245601_infercnv_threshold_sensitivity_grid.tsv` --
  coarse 11x9 grid, malignant fraction per sample.
- `results/tables/gse245601_infercnv_local_score_sensitivity.tsv` -- dense
  17-point local sweep around the actual 0.01 floor, per sample.
- `results/figures/gse245601_infercnv_threshold_diagnostics/` -- 2 figures
  per sample (24 total) + 1 combined local-sensitivity figure.
- `src/gse245601_infercnv_threshold_diagnostics.py` +
  `tests/test_gse245601_infercnv_threshold_diagnostics.py` (21 tests,
  including the real-data exact-reconstruction check).

## Codex independent review

Reviewed inline (sandbox filesystem access unavailable in this
environment, as in prior audit phases): the recomputation formula, the
failure-category logic, and the interpretive claims (a)-(e). Verdict:
core reconstruction and arithmetic confirmed correct and exactly matching
the frozen R formula on this data (the one documented `na.rm` asymmetry is
inert here because the data contain no NaNs, now asserted rather than
assumed). Two substantive corrections were requested and applied:

1. The original "not a knife-edge" claim was not supported by the coarse
   0.005-spaced grid alone -- **the dense local sweep (Point 4) was added
   specifically in response to this**, and shows the answer is
   tumor-specific rather than a single cohort-wide statement.
2. "Intrinsic property of the metric" was reframed as a hypothesis about
   *mechanism* (dilution by unaffected genome) that this document does not
   establish, versus the directly-supported, narrower claim that this is
   a threshold-metric mismatch and not a coding error.

Both corrections are reflected in the text above.
