---
title: Figure 1 palette change — colour only
status: post_freeze_exploratory
analysis_date: 2026-08-19
branch: figure1-palette (unmerged)
---

# Figure 1: new candidate colours, nothing else

The original lollipop design is unchanged. Only the four candidate gene
colours differ. A previous session had rebuilt this figure as a scatter with
false discovery rate on the vertical axis; that redesign is not used here and
nothing from it is carried over.

**Output:** `results/figures/poster_crispr_v2/CRISPR_discovery_v2.{png,pdf,svg}`
**Renderer:** `src/poster_crispr_discovery_v2.py` · **Wrapper:**
`scripts/poster/01b_crispr_discovery_v2.py`

## The colour change

| Gene | From | To | Changed |
|---|---|---|---|
| KDM1A | `#D55E00` | `#D55E00` | no |
| TLK2 | `#CC79A7` | `#009E73` | **yes** |
| USP34 | `#0072B2` | `#6A3D9A` | **yes** |
| VEZF1 | `#E69F00` | `#56B4E9` | **yes** |

The previous set contained two oranges (`#D55E00`, `#E69F00`) and, at small
mark sizes, a blue that paired visually with the pink-purple (`#0072B2` with
`#CC79A7`), so the four genes read as two similar pairs. The new set gives four
clearly separated hues. Three are Okabe-Ito colour-vision-safe values; the
purple `#6A3D9A` is a deliberate addition, since Okabe-Ito contains no true
purple.

## How it was made without touching anything else

`src/post_audit_sensitivity_visualization.py` (which defines the shared
`FOCUS_COLORS` dict) and `src/poster_crispr_discovery_v1.py` are **untouched**.
Editing the shared dict in place would have silently recoloured every other
committed figure that imports it.

Instead: a new `src/poster_palette.py` holds `GENE_COLOURS` as the single
future source of truth, and `src/poster_crispr_discovery_v2.py` is v1 with one
substantive difference — it imports `GENE_COLOURS` from the palette module
rather than `FOCUS_COLORS` from the shared visualisation module. Everything
else is behaviourally identical: same sort order, same bar and marker geometry,
same title and subtitle strings, same typography, same greys for the nine
non-candidate genes, same axis treatment, same figure size, same margins.

Two additions beyond the colour import, both required by the task and neither
visual: the verification gate described below, and `SOURCE_DATE_EPOCH` pinning
plus a `post_freeze_exploratory` metadata label on each saved file.

## Verification against the frozen source

Asserted **before** rendering, in `verify_against_frozen()`; a mismatch raises
`VerificationError` and refuses to plot rather than substituting.

- the thirteen genes, **by name and in plotted order** (ascending effect size)
- each gene's **effect size**, to a tolerance of 1e-6
- the genome-wide fitted-gene count (19,103)

**Result: passed** — 13 genes, effect sizes and order all match.

## Diff against the committed figure 01

`scripts/poster/01b_crispr_v2_diff_report.py` compares v2 against
`poster/final_figures/01_crispr_discovery.png` two ways, and writes
`diff_report.json`.

**Structural** — both figures rebuilt in-process and their artists compared:

| Property | Result |
|---|---|
| Thirteen genes, same order | identical |
| Axis limits (x and y) | identical |
| Tick positions (x and y) | identical |
| Title and both subtitle strings | identical |
| Font sizes (title, subtitles, tick labels, axis label) | identical |
| Figure size | identical (12.5 × 7.6 in) |
| Margins / subplot parameters | identical |
| Marker sizes and line widths | identical |
| Spine visibility | identical |
| **Differences outside colour** | **none** |

The only properties that differ are `ytick_colors` and `hline_colors` — the two
that carry the gene colours, which is exactly the intended change.

**Pixel** — the two PNGs compared pixel by pixel:

- same image dimensions
- **103,019 pixels differ (1.333%)**
- **0 of those are unexplained.** Every differing pixel is accounted for by one
  of the three colour changes: near an old colour in the committed figure and
  near the corresponding new colour in v2, allowing for anti-aliased blends
  toward white and toward the axis grey.

`"PASS": true`.

## Colour-vision check

Deuteranopia and protanopia simulated with the Machado, Oliveira & Fernandes
(2009) severity-1.0 matrices applied in linear sRGB — the direct computation
already used elsewhere in this repository. No colour-vision library exists in
this environment and none was installed. Applied both to the four hex values
and pixel-wise to the rendered PNG
(`CRISPR_discovery_v2_deuteranopia.png`, `_protanopia.png`;
full table in `cvd_palette_simulation.tsv`).

**All four genes remain distinguishable under both simulations.**

Pairwise CIE76 ΔE (ΔL\* in brackets):

| Pair | normal | deuteranopia | protanopia |
|---|---|---|---|
| **USP34 / VEZF1** (purple / light blue — the at-risk pair) | 61.9 (34) | **32.3 (32)** | **40.9 (38)** |
| KDM1A / TLK2 | 102.4 | 53.3 | 37.2 |
| KDM1A / USP34 | 108.6 | 107.8 | 101.3 |
| KDM1A / VEZF1 | 113.6 | 101.1 | 88.7 |
| TLK2 / USP34 | 103.9 | 56.1 | 70.5 |
| TLK2 / VEZF1 | 59.5 | 48.9 | 51.7 |

Purple and light blue converge in hue under both simulations, as expected, and
separate on **lightness** by 32–38 L\* units — the mechanism the task
anticipated. The smallest separation anywhere is ΔE 32.3, comfortably above the
threshold at which two marks become confusable. In the simulated renders USP34
reads as a dark blue and VEZF1 as a pale blue; KDM1A as ochre and TLK2 as a
grey-olive.

**The purple was therefore not darkened, and `#6A3D9A` stands unchanged.**

For comparison, the palette this replaces had a worst-case CVD separation of
**ΔE 18.3**; the new set raises that to **32.3**, a 1.8× improvement in the
tightest pair.

## Reproducibility

`SOURCE_DATE_EPOCH` is pinned, and two renders produce **byte-identical PNG and
PDF**. SVG is not byte-reproducible — matplotlib emits per-run element ids —
which matches the behaviour already documented for the committed poster
figures.

## Freeze

Verified at start and end of this work:

| Check | Result |
|---|---|
| `git rev-parse science-freeze-2026-08-15^{commit}` | `9a1b7777d6c69c2be44f16f25bc950769dc2ffda` |
| `git diff science-freeze-2026-08-15 -- results/tables/evidence_freeze/` | empty |
| `poster/` changed on this branch | no |
| `README.md`, `PREANALYSIS.md`, `docs/` changed | no |
| `src/post_audit_sensitivity_visualization.py`, `src/poster_crispr_discovery_v1.py` changed | no |

One note on the freeze command as stated: `git diff science-freeze-2026-08-15 -- results/tables/evidence_freeze/ poster/`
cannot return empty on any branch, because `poster/` did not exist at the
freeze tag and was created by later commits already on `main`. The meaningful
check — `poster/` unchanged relative to this branch's base commit — is the one
recorded above, and it is empty.
