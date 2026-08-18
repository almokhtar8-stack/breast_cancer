---
title: Post-poster strengthening analyses — summary
status: post_freeze_exploratory
analysis_date: 2026-08-18
branch: post-poster-strengthening (unmerged; merging is a human decision)
---

# Summary

Three strengthening analyses were requested. Two were performed. One was
stopped at its own feasibility gate, as the brief required.

Everything here is **post-freeze exploratory**. No frozen value, threshold,
ranking, shortlist or conclusion was changed. Every new table carries a
`post_freeze_exploratory` column, and that is enforced by test.

---

## 1. Meta-analysis replacing vote counting

**What changed about how the evidence should be described.**

The frozen table's `resistance_fdr05_count` should not be read as a measure of
resistance evidence. Pooling the three datasets shows that **no candidate
reaches BH FDR 0.05 for a pooled resistance effect in any arm, under either
standard-error definition** — the smallest value anywhere is FDR 0.100, and at
k = 2-3 studies even that is a descriptive plug-in value rather than a
calibrated test. The correct
description of the resistance transcriptomics for all 13 candidates is "weak",
not "supported in N of 3 datasets".

Three specific things the count was hiding:

- **VEZF1** shares the top vote-count rank on a single significant dataset while
  its two contributing datasets point in **opposite directions** (I² = 92.5,
  Cochran Q p = 0.00088). It falls four places when direction and effect size
  are used. This is the count's clearest failure.
- **SUPT4H1** has a vote count of zero yet the largest pooled effect of any
  candidate (+0.451, consistent direction). The count treats it worst; the
  pooled estimate treats it best.
- **HMGB1** and **EIF4ENIF1** fall 9 and 8 places because their pooled effects
  point the *wrong way*. A count of zero is direction-blind.

Two reporting defects in the count itself: **USP17L29** is recorded as count 0,
which reads as "tested and null" but means "removed by the expression filter" in
two datasets and "never measured" in the third; and
**TSR3, CTDNEP1 and PET117** have no row in the frozen evidence table at all, so
no count exists for them.

**Does anything contradict the frozen conclusions? No.** USP34 — the frozen
rank-1 gene — is rank 1 here too, and is the only candidate with even a
nominally significant pooled effect (and even that only under a sensitivity SE
— see below). The frozen shortlist was built primarily on
CRISPR effect, which this analysis does not touch. The candidate ordering *below*
rank 1 does change; that is reported as a discrepancy and has not been acted on.

One result is explicitly **not** a conclusion: USP34's pooled p is 0.078 under
the primary SE and 0.049 under the `f_calibrated` sensitivity. It straddles 0.05,
and at k = 2 studies the plug-in random-effects p-value is anti-conservative
anyway. USP34 must not be described as significant.

Full report: `META_ANALYSIS_REPORT.md`. Tables:
`results/post_poster/meta_analysis/`.

---

## 2. Characterising the GSE111151 null

**What changed about how the evidence should be described.**

"GSE111151 candidate differential expression was largely null" should be
replaced by: **GSE111151 was underpowered for every candidate effect it
observed, so its nulls are uninformative.** At 80% power and alpha 0.05 it could
detect a 1.5-fold change (median minimum detectable effect 0.627 log2); the
largest candidate effect it saw was 0.24. All 12 testable candidates sit below
their own detection floor.

Extending the same calculation to the other two datasets answers "which of the
three had the power to say anything at all": **none of them, at the effect sizes
actually observed.** Across 36 (gene, dataset) cells, only two have an observed
effect exceeding their own minimum detectable effect, and both are the datasets'
own already-significant hits. GSE240112 is the least sensitive (median MDE 1.007
log2, a 2.0-fold change) and is additionally confounded with biobank.
GSE118713 looks the most sensitive, but its residual variance is within-lineage
technical variation, so its true sensitivity to a biological resistance effect is
worse than it appears.

**Does anything contradict the frozen conclusions? No.** This *removes* a piece
of apparent negative evidence rather than adding positive evidence. The frozen
shortlist was not built on those nulls being informative. No verdict in this
analysis asserts a genuine negative, and that is test-enforced.

Full report: `POWER_ANALYSIS_REPORT.md`. Tables: `results/post_poster/power/`.

---

## 3. Reconstruction gap against the published annotations

**Not performed. The feasibility gate failed, and the analysis was stopped
before any proxy was attempted.**

The authors' per-cell `Epi. Tumor` / `Epi. Nontumor` labels are not publicly
available. Checked, and recorded machine-readably: GEO carries 26 Cell Ranger
`.h5` matrices and nothing else; all 76 worksheets of the paper's seven
supplementary tables were opened and scanned, containing **zero** cell barcodes
(they are gene-level DE and signature tables, largest 910 rows); the authors'
repository at the pinned commit `ceabf3f` has only `.ipynb`, `.md`, `.png` and
`.r` files in its entire history, with no releases, no tags and no Git LFS. The
processed data sits behind dbGaP controlled access (`phs003186.v1.p1`), which is
not a public comparator.

Reconstructing the authors' labels in order to compare them against our
reconstruction of the authors' labels would measure agreement between two runs
of our own pipeline. That circularity is why no proxy was produced.

**Does anything contradict the frozen conclusions? No.** The gap remains
qualitative, as `docs/CNV_METHOD_AUDIT.md` left it. What can be said without any
comparison is that our own malignant fractions are strongly sample-dependent —
from 3 of 6,108 cells (Tumor_01) to 858 of 1,618 (Tumor_03) — and that the three
Track B eligible tumours (Tumor_02, Tumor_03, Tumor_07) are precisely the
high-malignant-fraction ones.

Full report: `ANNOTATION_CONCORDANCE_REPORT.md`. Tables:
`results/post_poster/annotation_concordance/`.

---

## Constraint verification

Recorded at the start of this work and again before the final commit.

| Check | Start | End |
|---|---|---|
| `git rev-parse science-freeze-2026-08-15^{commit}` | `9a1b7777d6c69c2be44f16f25bc950769dc2ffda` | `9a1b7777d6c69c2be44f16f25bc950769dc2ffda` |
| SHA-256 of `THERAPEUTIC_SHORTLIST_FREEZE.tsv` | `b6990d7e20f1191df7ebbca18be7542cc312e6e8489f05af6e4b4d31e66fb9dc` | `b6990d7e20f1191df7ebbca18be7542cc312e6e8489f05af6e4b4d31e66fb9dc` |
| `git diff science-freeze-2026-08-15 -- results/tables/evidence_freeze/` | empty | empty |
| `poster/`, root `README.md`, `PREANALYSIS.md`, `docs/*PREANALYSIS*.md` changed on this branch | no | no |
| working tree clean apart from the new work | yes | yes |

One caveat worth recording: running the existing poster figure tests re-renders
their PDF and SVG outputs, which are not byte-reproducible (documented by commit
`477c992`). Eight such files under `results/figures/poster_*` showed as modified
during the test run and were reverted before committing — no figure content
changed, and no file outside the permitted paths is part of this branch.

Two notes on the brief's stated checks:

- `git rev-parse science-freeze-2026-08-15` returns `cd9750d8…`, the **annotated
  tag object**. The commit needs `^{commit}` and is `9a1b7777…` as expected.
- `git diff science-freeze-2026-08-15 -- poster/` is **not** empty and cannot be:
  `poster/` did not exist at the freeze and was created by the seven
  post-freeze audit commits already on `main`. The meaningful check for this
  branch is against its base commit `477c992`, which is empty.

## Environment note

`environment.yml` was reported to be missing `h5py`. **Not confirmed** — `h5py`
is present at line 10, it imports (3.16.0), and `pytest --collect-only` gathers
1,770 tests with no errors. No fix commit was made.

A genuine, *pre-existing* packaging gap was found instead: the frozen R analyses
and this branch's refit require `r-yaml` and `r-data.table`, which the `bc`
environment lacks; they run in the `sc245601` micromamba environment. Recorded,
not fixed — changing `environment.yml` is outside this branch's scope.

## Test results

New tests: **77 passed** (`tests/test_post_poster_*.py`). Run together with the
release-integrity and freeze-shortlist tests: **104 passed, 0 failed**.

Wider suite (95 test files, run file by file): **89 files produced a result —
1,466 passed, 1 skipped**. Six files produced no result because they hang
indefinitely in this environment and survive their own `timeout`, blocking any
sequential run:

```
test_final_pharmacogenomics.py   test_poster_exploration_v2.py
test_independent_validation.py   test_poster_exploration_v3.py
test_post_audit_sensitivity.py   test_poster_figures.py
```

This reproduces the limitation already recorded for the previous session, where
the suite likewise could not be completed in one process. **The suite was not
fully run, and this branch is not claiming a green suite.**

One failure appeared and is not a regression:
`test_nebula_plots_final.py::test_final_figures_written_with_expected_filenames_and_dimensions`
failed while several figure-rendering tests were running concurrently against
the same output directory. Run on its own it passes (19 passed). It is a
concurrency artifact of how I ran the suite.

What the unrun files do and do not leave open: this branch **adds only new
files**. No existing tracked source, test, config or data file was modified —
`git status` before committing showed only the new paths, with the working tree
otherwise clean. The unrun tests therefore exercise code this branch did not
touch. That is an argument from the diff, not from a green suite, and is offered
as such.

## Independent review

Codex reviewed at checkpoints 1 and 2, and both materially changed the work.
Checkpoint 1 reshaped the plan (primary arm, alpha choice, verdict wording,
heterogeneity estimators). Checkpoint 2 rejected my central inferential
choice — I had made the statistic-derived `f_calibrated` the primary standard
error, which is outcome-dependent and cannot legitimately weight an
inverse-variance pooling; `wald` is now primary throughout, and that also
resolved a report-versus-table mismatch and a circularity in the power
calculation. Acting on Codex's drop-accounting point also exposed a factual
error in my own report about why USP17L29 is missing. All of it, including two
bugs my own tests caught, is logged in `CODEX_REVIEWS.md`.

## Status

The branch is left **unmerged**.
