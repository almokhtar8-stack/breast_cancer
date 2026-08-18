---
title: Candidate volcano figure — proposed Figure 2 replacement
status: post_freeze_exploratory
analysis_date: 2026-08-18
branch: figure2-volcano (unmerged; whether this replaces Figure 2 is a human decision)
---

# Figure 2 volcano candidate

**post_freeze_exploratory.** No statistic is recomputed; every plotted value
is read from a frozen genome-wide DE table. The existing Figure 2, its
renderer, and `poster/figure_manifest.tsv` are untouched. This is a candidate
replacement, not a swap.

Outputs: `results/figures/poster_candidate_volcano_v1/candidate_volcano_v1.{png,pdf,svg}`
plus `volcano_manifest.tsv`, `candidate_values_plotted.tsv` and
`verification_against_frozen.tsv` alongside. Code:
`src/poster_candidate_volcano_v1.py`, wrapper
`scripts/poster/02b_candidate_volcano.py`, tests
`tests/test_poster_candidate_volcano_v1.py`.

## Which GSE118713 contrast, and why

**`TAMR_vs_MCF7`**, determined from the file rather than assumed: the table
carries three contrasts (`TAMR_vs_MCF7`, `FASR_vs_MCF7`, `TAMR_vs_FASR`), and
the module selects the unique label whose `_vs_` endpoints are exactly
{TAMR, MCF7}, raising if that match is not unique. This is the contrast the
current Figure 2 and `candidate_evidence_summary.tsv` report (its
`tamr_vs_mcf7_*` columns), and `evidence_long.tsv` records it as the primary
GSE118713 contrast. Row accounting: 44,514 rows in, 14,838 kept, 29,676
excluded as the other two contrasts.

## Two source-file corrections, made before plotting and confirmed by the requester

The verification gate exposed that two of the four files named in the brief do
**not** reproduce the frozen reference values. Both substitutions were put to
the requester and confirmed before any figure was rendered:

1. **Panel A reads `gse118713_differential_expression_unredacted.tsv.gz`**
   (14,838 genes/contrast), not the redacted file the brief names (14,836).
   The redacted file still carries the KDM1A/RCOR1 blinding retired on
   2026-08-10, so KDM1A has no row in it at all — the gate's KDM1A FDR 0.494
   is only checkable, and only plottable, from the unredacted table. The
   USP34, VEZF1 and TLK2 rows are identical in both files.
2. **Panel C reads `gse240112_pseudobulk/tumor_cell_genomewide_de.tsv.gz`**
   (18,428 genes), not the epithelial track the brief names. The epithelial
   file fails the gate on **all four** candidates (KDM1A 0.332 vs 0.59, TLK2
   0.759 vs 0.88, USP34 0.196 vs 0.23, VEZF1 0.0074 vs 0.0195); the
   tumour-cell track reproduces all four, and `evidence_long.tsv` states that
   the tumour-cell track is the primary analysis while all-epithelial is
   sensitivity-only. The FDR<0.05 candidate count is 2 either way, but the
   annotated VEZF1 value would have read 0.007 instead of the frozen 0.0195.

Both paths are taken from `config/config.yaml`
(`cross_dataset_genomewide.inputs`), which already names exactly these four
files as the genome-wide evidence inputs. Track A is used for GSE245601 as
instructed; Track B is not read anywhere.

## Verification against the frozen values

All 16 candidate FDRs and both quoted log2FCs match the frozen reference
numbers (`verification_against_frozen.tsv`; enforced again by test, and the
module raises rather than plots on any mismatch). One tolerance note: the
brief's reference values are quoted at 2–4 decimal places, and at a strict
1e-3 seven of the 2-dp quotes would "fail" purely by rounding (e.g. KDM1A
GSE111151: table value 0.905931 vs quote "0.91", difference 0.0041). The
comparison therefore uses the larger of 1e-3 and half a unit in the last
quoted decimal place — i.e. "0.91" accepts exactly the values that round to
0.91. Every value quoted at ≥3 decimals matches within the strict 1e-3.

Exactly **2** candidate points are significant at FDR 0.05 across all four
panels — USP34 in GSE118713 (0.0073) and VEZF1 in GSE240112 (0.0195) — and
both the module and the tests assert that count is exactly 2.

## What this figure shows that the heatmap does not

- **Significance is on an axis.** The heatmap's z-score colour scale spans
  its full blue-to-red range regardless of significance, so KDM1A at FDR 0.49
  can be — and was — misread as strongly upregulated. Here a non-significant
  candidate is a hollow ring below a dashed line; there is no colour intensity
  to over-read.
- **The 2-of-16 count is legible without reading a number**: the only
  saturated, filled candidate points in the whole figure are the two
  significant ones, each annotated with its FDR — the only two numbers a
  reader needs.
- **Genome-wide context.** Each candidate is seen against every gene the
  dataset tested, so "significant in this dataset" is visibly calibrated
  against what significance looks like there (GSE111151 and GSE245601 have
  almost no signal genome-wide; GSE118713 and GSE240112 have thousands of
  significant genes among which the candidates are unremarkable except the
  two).

## What it loses relative to the heatmap

- **Sample-level detail.** The heatmap shows every biological observation
  (29 rows: individual sublines, tumours, and patients), which is how a reader
  sees the n behind each contrast and the within-group spread. The volcano
  collapses each dataset to one point per gene.
- **Absolute expression level.** The heatmap's per-sample values convey
  where each gene sits in each sample; the volcano carries only the contrast
  (log2FC) and its FDR — a gene barely expressed and a gene highly expressed
  look the same at equal fold change.

If the replacement is adopted, those two losses argue for keeping the
sample-level view available somewhere (supplement or a side panel), not for
keeping the z-score colour encoding.

## Axis and layout decisions

- Y is **−log10(FDR)** (BH-adjusted), labelled as such, because the
  significance claim is made on FDR; the dashed reference line sits at
  FDR 0.05, with a vertical line at log2FC 0. Y-limits are shared across all
  four panels so panels are directly comparable — GSE111151 and GSE245601
  are honestly flat.
- **No point is trimmed from view in any panel.** Panel C's x-range (±16.5)
  is wide because pseudobulk fold changes at low expression are extreme; the
  caption notes this rather than clipping. Candidates are never trimmed.
- Panel titles put the biological description first and the accession second,
  smaller and lighter, matching Figure 2 — including its corrected GSE245601
  headline ("Acute 12 h tamoxifen — per-tumour pseudobulk", not
  "single-cell"). Candidate colours are imported from the Figure 2 renderer
  (`FOCUS_COLORS`), not chosen here.
- PDF bytes are reproducible (`SOURCE_DATE_EPOCH=1600000000`, as
  `scripts/poster/build_all.py` pins it; covered by test); SVG ids are pinned
  with `svg.hashsalt`.

## Test results

New tests: **15 passed** (`tests/test_poster_candidate_volcano_v1.py`),
covering the verification gate (including that it fails loudly on a corrupted
value), the exact significant count of 2, the FDR 0.05 threshold asserted
rather than assumed, per-panel row accounting with no candidate dropped, all
three output formats, the `post_freeze_exploratory` metadata label, manifest
contents, and PDF byte-reproducibility.

Full suite: **1,437 passed, 1 skipped, 0 failed** in a single process
(7 m 13 s), with six test files excluded by name because they hang
indefinitely in this environment (documented on the
post-poster-strengthening branch, whose SUMMARY.md records the same
exclusion):
`test_final_pharmacogenomics.py`, `test_independent_validation.py`,
`test_post_audit_sensitivity.py`, `test_poster_exploration_v2.py`,
`test_poster_exploration_v3.py`, `test_poster_figures.py`.

Running the suite re-renders 19 PDF/SVG/adjustText-PNG files under
`results/figures/poster_*`, which are not byte-reproducible (documented by
commit `477c992`); they were reverted before committing, and no figure
content changed. The new volcano tests render only into pytest temporary
directories, so the committed volcano outputs are not churned by the suite.

## Freeze verification

| Check | Start | End |
|---|---|---|
| `git rev-parse science-freeze-2026-08-15^{commit}` | `9a1b7777d6c69c2be44f16f25bc950769dc2ffda` | `9a1b7777d6c69c2be44f16f25bc950769dc2ffda` |
| `git diff science-freeze-2026-08-15 -- results/tables/evidence_freeze/` | empty | empty |
| `poster/`, root `README.md`, `PREANALYSIS.md`, `docs/` changed on this branch (vs base `477c992`) | no | no |

Files changed on this branch, in full: the volcano module, wrapper, tests,
the six files under `results/figures/poster_candidate_volcano_v1/`, and this
note. Nothing else.

(The brief's combined `git diff science-freeze-2026-08-15 -- poster/` check
cannot return empty on any branch: `poster/` did not exist at the freeze tag
and was created by post-freeze commits already on `main`. The meaningful
check — `poster/` unchanged relative to this branch's base commit — is the
one recorded.)
