---
title: Minimum detectable effects for the three resistance datasets
status: post_freeze_exploratory
analysis_date: 2026-08-18
branch: post-poster-strengthening
scope: GSE111151, GSE118713, GSE240112 (GSE245601 excluded by design)
---

# Characterising the GSE111151 null

**post_freeze_exploratory.** No frozen result is changed. No new ranking is
produced.

## The question

`docs/FINAL_PUBLIC_REPO_AUDIT.md` line 122 and `poster/README.md` line 99 both
record that "GSE111151 candidate differential expression was largely null".
Neither states an FDR; measured here from the frozen table, the two poster leads
present in that dataset are USP34 at FDR 0.632 and VEZF1 at 0.608. *"We found
nothing"* and *"we could not have found anything at this sample size"* are
different scientific statements. This report distinguishes them, per gene, and
does the same for the other two resistance datasets so the three are comparable.

## What was computed, and what was deliberately not

For each gene, the **minimum detectable effect (MDE)**: the smallest
|log2 fold change| the study would detect with 80% power at a stated alpha,
given **that gene's own fitted dispersion, the real library sizes and the real
design**, obtained by refitting each frozen model with identical parameters
(`src/post_poster_de_refit.R`, verified to reproduce the frozen tables to
machine precision — see `META_ANALYSIS_REPORT.md`).

**This is not observed power.** Power evaluated at the observed effect estimate
is a monotone function of the p-value and adds nothing to it. No observed
p-value is converted into a power value anywhere in this analysis, and
`tests/test_post_poster_power.py::test_mde_does_not_depend_on_the_observed_effect`
enforces that the MDE is a function only of SE, df and alpha.

Reproduce with:

```
Rscript src/post_poster_de_refit.R      # sc245601 R env
python -m src.post_poster_power
pytest tests/test_post_poster_power.py
```

### Significance levels

| Name | Value | Status |
|---|---|---|
| `nominal_0.05` | 0.05 | prespecified per-gene type-I error — the headline |
| `bonferroni_13` | 0.05/13 = 0.003846 | prespecified, correcting across the 13-candidate family |
| `observed_bh_threshold` | dataset-specific | **descriptive only** |

The third is the largest nominal p-value that actually attained BH FDR ≤ 0.05 in
each dataset (GSE111151 7.5e-05, GSE118713 0.0151, GSE240112 0.0144). It is
outcome-dependent — a function of the realised p-value distribution, including
the gene being assessed — so it is **not** called "power at FDR 0.05" and is not
used for the headline. There is no fixed nominal alpha equivalent to BH FDR
control; genuine transcriptome-wide FDR power would require simulating the whole
pipeline under an assumed alternative.

### Only the information-based standard error may define an MDE

The meta-analysis carries a second quantity, `f_calibrated` = |log2FC| / √F, as
a pooling sensitivity. It is **deliberately not available here**: it contains
the observed effect by construction, so an MDE built on it would be partly
outcome-derived — precisely the circularity this analysis exists to avoid.
`compute_power_table()` raises rather than accepting it, and a test asserts the
refusal. Every MDE below uses the model-based Wald standard error.

### The conditional-approximation caveat

The MDE holds each gene's standard error fixed at its fitted value. A
workflow-level calculation would resimulate filtering, TMM, dispersion
estimation, empirical-Bayes moderation and BH under an assumed alternative, and
the standard error itself moves with the mean and the effect. These numbers are
tight enough to separate "underpowered" from "well powered" and not tight enough
to quote to two decimals.

## Results

MDE at 80% power, alpha = 0.05, against the observed effect. Bold = observed
effect exceeds the MDE.

| Gene | GSE111151 MDE / obs | GSE118713 MDE / obs | GSE240112 MDE / obs |
|---|---|---|---|
| CTDNEP1 | 0.64 / 0.04 | 0.67 / 0.36 | 1.01 / 0.12 |
| EIF4ENIF1 | 0.47 / 0.13 | 0.58 / 0.53 | 0.77 / 0.06 |
| HMGB1 | 0.62 / 0.13 | 0.46 / 0.08 | 1.46 / 0.16 |
| ICK | 0.57 / 0.17 | 0.50 / 0.11 | 1.09 / 0.36 |
| KDM1A | 0.63 / 0.07 | 0.46 / 0.15 | 0.81 / 0.19 |
| PET117 | 0.80 / 0.02 | 0.42 / 0.06 | 1.00 / 0.25 |
| SUPT4H1 | 0.63 / 0.02 | 0.51 / 0.04 | 1.25 / 1.00 |
| TADA2B | 0.59 / 0.06 | 0.52 / 0.01 | 1.15 / 0.43 |
| TLK2 | 0.42 / 0.07 | 0.74 / 0.15 | 0.81 / 0.06 |
| TSR3 | 0.70 / 0.09 | 0.78 / 0.47 | 1.14 / 0.43 |
| USP34 | 0.46 / 0.16 | 0.42 / **0.59** | 0.77 / 0.40 |
| VEZF1 | 0.64 / 0.24 | 0.79 / 0.43 | 1.00 / **1.15** |
| USP17L29 | not tested | not tested | not tested |

Median MDE across the 12 testable candidates
(`dataset_sensitivity_summary.tsv`):

| alpha | GSE111151 | GSE118713 | GSE240112 |
|---|---|---|---|
| 0.05 | 0.627 | 0.514 | 1.007 |
| 0.05/13 | 0.938 | 0.779 | 1.555 |
| observed BH threshold | 1.480 | 0.635 | 1.266 |

## Answer to the GSE111151 question

**The GSE111151 candidate null is uninformative, for all 12 testable
candidates.** At alpha = 0.05 and 80% power, GSE111151 could detect a 1.5-fold
change (median MDE 0.627 log2). The largest candidate effect it actually
observed was 0.24 (VEZF1) — well under half its own detection floor, and every
one of the 12 sits below its own MDE. Nothing about "FDR 0.63 for USP34" tells
us that USP34 is unchanged in resistance; it tells us the study could not have
distinguished the observed change from zero at this sample size, whatever the
truth is.

Under the family-wise Bonferroni alpha the floor rises to 0.94 log2 — a 1.9-fold
change — and all 12 remain below it.

## Which dataset had the power to say anything at all

**None of them, at the effect sizes actually observed.** Across 36
(gene, dataset) cells at alpha = 0.05, only **two** have an observed effect
exceeding their own MDE:

| Gene | Dataset | observed | MDE | external benchmark | observed FDR |
|---|---|---|---|---|---|
| USP34 | GSE118713 | 0.590 | 0.416 | 0.400 | 0.0073 |
| VEZF1 | GSE240112 | 1.149 | 0.997 | 0.427 | 0.0195 |

Both of those are the datasets' *own* significant hits, so this is not new
information about them — it is confirmation that the two significant candidate
results are the only two the studies were sensitive enough to produce.

Ranking the three by sensitivity: **GSE118713 (0.514) > GSE111151 (0.627)
≫ GSE240112 (1.007)**. But GSE118713's apparent sensitivity is inflated: its
three replicates per arm are of a single MCF7/TAMR lineage, so its residual
variance is within-lineage technical variation, not biological variation between
independent resistance models. Its true sensitivity to a *biological* resistance
effect is worse than 0.514 by an unidentifiable factor. On that reading
**GSE111151 is the most informative of the three about resistance, and it still
could not detect any candidate effect it observed.**

GSE240112's floor of ~1.01 log2 (a 2.0-fold change) is the worst of the three,
which is consistent with n = 3 vs 3 unpaired pseudobulk. Its `group` is also
perfectly confounded with biobank, so even an effect above that floor would not
be attributable to recurrence.

## How to read the interpretation column

`candidate_minimum_detectable_effects.tsv` carries one interpretation per
(gene, dataset, alpha):

| Value | Meaning |
|---|---|
| `null_uninformative_observed_effect_below_mde80` | the observed effect is smaller than the smallest the study could detect; the null carries no information about this gene |
| `sensitive_to_observed_magnitude` | the study could detect an effect the size of the one observed |
| `sensitive_to_observed_magnitude_but_not_to_external_benchmark` | as above, but the MDE still exceeds the largest effect this gene shows in the other datasets |
| `not_tested` | the gene is absent from this dataset; no MDE and no verdict |
| `not_estimable` | tested, but the MDE could not be computed |

**No value asserts a genuine negative**, and this is enforced by test. Failing to
reject while nominally able to detect a hypothetical effect is not evidence of
absence: 80% power still misses one true effect in five, and the observed
estimate is itself noisy. Establishing absence would need a prespecified
smallest effect of interest and an equivalence test — this project never
prespecified one, so that analysis is not available and is not faked.

The `external_reference_abs_log2fc` column is the largest |log2FC| that gene
shows in the *other* datasets. It is an explicitly post-hoc external benchmark,
not a truth.

## Genes not tested

**USP17L29 is absent from all three fitted tables**, but not for the same reason
in each, and the `not_tested_reason` column keeps them apart:

| Dataset | Reason |
|---|---|
| GSE111151 | present before filtering; removed by `filterByExpr` |
| GSE118713 | present before filtering; removed by the expression filter |
| GSE240112 | not in the gene annotation at all; never measured |

It receives `not_tested` rows everywhere and no MDE, rather than a null verdict.
Note that `resistance_fdr05_count` records it as 0, which reads as "tested and
null" but is in two datasets "measured and filtered out for low expression" and
in the third "never measured".

## What this changes about how the evidence should be described

The phrase "GSE111151 candidate differential expression was largely null" should
be replaced by: *GSE111151 was underpowered for every candidate effect it
observed; its nulls do not constitute evidence against the candidates.* That is
a weaker statement about the data and a **more favourable** one for the
candidates — but it is not support for them either. It removes a piece of
apparent negative evidence rather than adding positive evidence.

Nothing here contradicts a frozen conclusion. The frozen shortlist was not built
on these nulls being informative.
