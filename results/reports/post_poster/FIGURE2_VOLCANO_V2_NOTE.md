---
title: Candidate volcano v2 — layout fix, project palette, two variants
status: post_freeze_exploratory
analysis_date: 2026-08-18
branch: figure2-volcano (unmerged; whether either variant replaces Figure 2 is a human decision)
supersedes: FIGURE2_VOLCANO_NOTE.md (v1, retained alongside for provenance)
---

# Figure 2 volcano candidate, v2

**post_freeze_exploratory.** v2 changes layout and colour only. No statistic
is recomputed, no source file changes, and no threshold or reference value
moves. v2 imports v1's data loaders and verification gate unchanged, so all
16 candidate values are still verified against the frozen reference numbers
before anything is plotted, and the significant count is still asserted to
be exactly 2 (USP34/GSE118713 FDR 0.0073, VEZF1/GSE240112 FDR 0.0195). The
two source-file corrections from v1 stand: unredacted GSE118713 table
(so KDM1A has a row) and the GSE240112 tumour-cell track (the one that
reproduces the frozen FDRs). v1 is left in place, untouched.

Outputs: `results/figures/poster_candidate_volcano_v2/`
`candidate_volcano_v2_variant_a_genomewide.{png,pdf,svg}`,
`candidate_volcano_v2_variant_b_candidates_only.{png,pdf,svg}`,
`volcano_v2_manifest.tsv` (both variants, SHA-256 per file),
`candidate_values_plotted.tsv`, `verification_against_frozen.tsv`,
`cosmetic_offsets.tsv`, `cvd_palette_simulation.tsv`, and
`cvd_simulation/` (deuteranopia and protanopia renderings of both PNGs).
Code: `src/poster_candidate_volcano_v2.py`, wrapper
`scripts/poster/02b_candidate_volcano_v2.py`, tests
`tests/test_poster_candidate_volcano_v2.py`.

## What changed from v1, and why

| v1 problem | v2 fix |
|---|---|
| Per-panel autoscaling: x ranges ±7 / ±10 / ±15 / ±3, so the same fold change sat in a different place in every panel and panel 4 was stretched while panel 3 was compressed | One shared x range and one shared y range per variant, derived from the data and asserted to contain every plotted point; identical panel geometry |
| Three of four candidates crowd near the origin in every panel | Variant A: a shared **zoom row** beneath the four main panels covering the candidate region, with a rectangle on each main panel marking it. Variant B: no backdrop, so the tight range itself is the zoom |
| Overlapping candidate rings not individually countable (GSE111151 KDM1A/TLK2 differ by 0.001 log2FC; GSE245601 TLK2/USP34/VEZF1 sit within 0.02) | Explicit draw order and thicker outlines; in A's zoom row a fixed, on-figure-disclosed horizontal spread; in B concentric rings at the true position with zero displacement (below) |
| Three-line caption unreadable at poster size | One-line caption on the figure; everything else moved here |
| Legend inside the header block | Legend strip below the panels |
| Panel titles one or two lines, so header heights differed | Two-line titles for all four; accession beneath, smaller and lighter; GSE245601 stays "acute 12 h tamoxifen, per-tumour pseudobulk" |
| Gene colours imported from the Figure 2 heatmap renderer | The project's data-tier gene palette defined **once** in `GENE_COLOURS` (`src/poster_candidate_volcano_v2.py`), for later figure passes to import; no colour is read from v1 or the heatmap renderer any more (test-enforced) |
| A second grey for "other genes at FDR<0.05" | Only the repository's neutral ladder is used; background genes are one grey and their significance is read from position above the line |

Main-panel labels in variant A: at ±16.5 log2FC any label sits on the
clouds, so the main panels carry **no** text labels; the zoom row directly
beneath labels all four candidates and annotates the two significant FDRs.
Variant B labels all four in the main panels.

## Colours

| Gene | Hex |
|---|---|
| USP34 | `#6A3D9A` |
| KDM1A | `#D55E00` |
| TLK2 | `#009E73` |
| VEZF1 | `#56B4E9` |

Neutral ladder: `#262626` figure text, `#555555` subtitles/axis names,
`#8c8c8c` tick labels, `#b0b0b0` threshold and reference lines, `#c9c9c9`
background points, `#e6e6e6` gridlines (declared; not drawn). No other
colour appears in the module — a test greps every hex literal against this
allow-list. White is used only as the page and as the thin edge of filled
markers.

### Colour-blindness simulation

No CVD library is installed (`colorspacious`, `daltonize`, `scikit-image`,
`colour` all absent) and none was added; the simulation is computed
directly from the Machado, Oliveira & Fernandes (2009) severity-1.0
matrices applied in linear sRGB, then CIE76 ΔE and ΔL* in CIELAB. Applied
both to the four hex values and, pixel-wise, to both rendered PNGs
(`cvd_simulation/`).

Pairwise ΔE76 (ΔL* in brackets):

| Pair | normal | deuteranopia | protanopia |
|---|---|---|---|
| USP34 / VEZF1 (purple / light blue — the flagged pair) | 61.9 (34) | **32.3 (32)** | **40.9 (38)** |
| KDM1A / TLK2 (orange / green) | 102.4 (4) | 53.3 (2) | 37.2 (13) |
| KDM1A / USP34 | 108.6 | 107.8 | 101.3 |
| KDM1A / VEZF1 | 113.6 | 101.1 | 88.7 |
| TLK2 / USP34 | 103.9 | 56.1 | 70.5 |
| TLK2 / VEZF1 | 59.5 | 48.9 | 51.7 |

**All four remain distinguishable under both simulations.** Purple and
light blue converge in hue (both read as blue) but keep a lightness gap of
32–38 L\* units, exactly the lightness separation the brief anticipated;
the smallest ΔE anywhere is 32. Visually in the simulated PNGs, USP34 is a
dark blue, VEZF1 a pale blue, KDM1A ochre and TLK2 a grey-tan. **The purple
was therefore not darkened and `#6A3D9A` stands unchanged.** The pair with
the least *lightness* separation is orange/green under deuteranopia (ΔL\*
2), which separates on chroma (ΔE 53); at poster size those two rings are
also always far apart in position.

## The inset decision

**A shared zoom row beneath the four main panels, not a per-panel inset.**
Reasons: (i) in panels 1 and 3 the genome-wide clouds occupy the upper
corners where an inset would have to sit, so an inset would occlude the
very backdrop variant A exists to show; (ii) a per-panel inset would be
~40% of a 3-inch panel — at final print size its labels would not be
readable; (iii) a row of four equal zoom panels shares one set of limits, so
the zoomed views are as comparable to each other as the main panels are.
Each main panel carries a rectangle marking the zoom region.

## Axis ranges

| | x (log2FC) | y (−log10 FDR) |
|---|---|---|
| Variant A, main panels | −16.53 to 16.53 | 0 to 8.73 |
| Variant A, zoom row | −1.5 to 1.5 | 0 to 2.5 |
| Variant B | −1.5 to 1.5 | 0 to 2.5 |

Variant A's limits are the data extremes (GSE240112 |log2FC| 15.59;
GSE118713 −log10 FDR 8.23) × 1.06, so every gene in every panel is inside;
that makes panels 1, 2 and 4 look empty at the edges, which is informative
and was kept. Panel 3's wide range reflects pseudobulk fold-change extremes
at low expression; nothing is trimmed from view in any panel, and a
candidate falling outside any range raises rather than clips (test-enforced).
The zoom region and variant B use the same limits, chosen to contain every
candidate (max |log2FC| 1.149, max −log10 FDR 2.136) with the FDR 0.05 line
at 1.30 inside; the tick labels ±1.5 are suppressed there so neighbouring
panels' edge labels do not touch. Variant B's range is thus 11× tighter in x
and 3.5× tighter in y than variant A's main panels.

## Coincident candidates: two different treatments, one per variant

Two clusters of coincident candidates exist (GSE111151 KDM1A/TLK2 differ by
0.001 log2FC; GSE245601 TLK2/USP34/VEZF1 sit within 0.02). Coincidence rule
(`COINCIDENCE_TOL = 0.06` in zoom/variant-B data units, both axes).

**Variant A (zoom row only): a fixed horizontal spread, disclosed on the
figure.** Cluster members are drawn `COINCIDENCE_STEP = 0.14` log2FC apart
in alphabetical order, centred on the cluster's mean x; **y is never moved,
so no point can change side of the threshold line** (test-enforced). Five
points are displaced; 0.14 is 0.85% of variant A's x range and invisible in
the main panels, which draw measured positions. Because displacing a data
point alters what the figure shows, the disclosure sits **in the legend
strip of the figure itself**, not only here: *"Zoom row: candidates within
0.06 log2FC of each other are drawn 0.14 apart horizontally so each can be
counted (5 points; y unchanged; measured values in the manifest)."*
(test-enforced to be present on A). Displacements are in
`cosmetic_offsets.tsv`; measured values in `candidate_values_plotted.tsv`.

| Dataset | Cluster | Gene | measured x | displayed x (A zoom row) |
|---|---|---|---|---|
| GSE111151 | KDM1A + TLK2 | KDM1A | 0.071 | 0.002 |
| | | TLK2 | 0.073 | 0.142 |
| GSE245601 | TLK2 + USP34 + VEZF1 | TLK2 | −0.036 | −0.180 |
| | | USP34 | −0.033 | −0.040 |
| | | VEZF1 | −0.051 | 0.100 |

**Variant B: zero displacement.** On B's ±1.5 axis the same 0.14 step
would be 9.3% of the range — three GSE245601 candidates that differ by 0.02
would be drawn nearly a tenth of the axis apart, and a reader taking fold
change off the axis could conclude they point in opposite directions. So
in B nothing moves: coincident candidates are drawn as **concentric rings
at their true shared position**, same centre and increasing radius
(`CONCENTRIC_SIZES`, alphabetical order, largest drawn beneath so none is
hidden). Every ring stays countable and no x is falsified. Tests assert
that B applies zero displacement, that **plotted x equals the source
log2FC exactly (`==`) for all sixteen candidate points**, that cluster
members are drawn at their own measured centres with strictly increasing
radius, and that no disclosure text is present on B because there is
nothing to disclose. Candidate rings are drawn unclipped so a ring centred
just above y = 0 keeps its lower arc.

## Text moved off the figure

- Hollow rings are candidates not reaching FDR 0.05; filled discs are the
  only saturated marks in the figure and are exactly the two significant
  points.
- Panel 4 (GSE245601) measures acute 12 h ex vivo tamoxifen response, not
  resistance.
- Panel 3 (GSE240112) has a wide fold-change range from pseudobulk extremes
  at low expression.
- No point is trimmed from view in any panel of either variant.
- (The offset disclosure for variant A is deliberately NOT moved off the
  figure — see above.)

## Reproducibility

Labels sit at fixed, hand-placed positions (`LABEL_POS`) with leader lines;
no collision solver. PDF bytes are pinned by `SOURCE_DATE_EPOCH=1600000000`
(as `scripts/poster/build_all.py` does), SVG ids by `svg.hashsalt`, and a
test asserts that **both PDF and PNG bytes are identical across two renders**
of each variant.

## Recommendation

**Take variant B to the supervisor as the primary proposal, with variant A
in hand as the supporting view.**

- *Variant B, candidates only.* Everything the figure exists to communicate
  is legible at poster distance: sixteen large marks, two of them filled,
  every label readable, and the FDR 0.05 line separating them. It fixes the
  misreading that motivated the whole exercise (a non-significant gene read
  as strongly regulated) with nothing competing for the eye. Its cost is
  that it looks like a small-n plot rather than a genome-wide analysis, and
  it does not show how ordinary the candidates are relative to the
  transcriptome.
- *Variant A, genome-wide backdrop.* Honest about scale — the reader sees
  that GSE118713 and GSE240112 have thousands of significant genes among
  which the candidates are unremarkable except two, and that GSE111151 and
  GSE245601 have almost none. But the clouds pull the eye upward while the
  candidates sit low; the counting has to happen in the zoom row, so the
  figure is taller and busier, and at poster size the zoom-row rings are
  small.

If only one can be on the poster, B carries the message; A belongs in the
supplement or as a backup slide, and its zoom row is what B is.

## Tests

`tests/test_poster_candidate_volcano_v2.py`: **22 passed** — gate passes on
all 16 values and fails loudly on a corrupted value; exactly 2 significant
at an asserted 0.05; the four gene colours match the dict exactly; v2 does
not import colours from v1 or the heatmap renderer; no hex outside the
palette + neutral ladder appears in the module; CVD simulation keeps all
pairs above ΔE 20 and the purple/light-blue pair above ΔL\* 25; shared axis
limits identical across the four panels of each variant (main and zoom rows
each); no candidate outside the visible range in either variant; coincidence
rule spreads horizontally only, never changes threshold side, and matches
the two documented clusters; variant B applies zero displacement with
plotted x == source log2FC for all 16 points, draws concentric rings of
strictly increasing radius at each member's own centre, and carries no
disclosure; variant A spreads 5 points and carries the disclosure text on
the figure; both variants build in all three formats; PDF
and PNG bytes reproducible; manifest covers both variants with distinct
SHA-256s.

Full suite (six known-hanging files excluded by name —
`test_final_pharmacogenomics.py`, `test_independent_validation.py`,
`test_post_audit_sensitivity.py`, `test_poster_exploration_v2.py`,
`test_poster_exploration_v3.py`, `test_poster_figures.py`): **1,459 passed, 1 skipped, 0 failed** in a single process (7 m 52 s) — v1's 1,437 plus the 22 new tests. Running the suite re-renders non-byte-reproducible PDF/SVG figures under `results/figures/poster_*` (documented by commit `477c992`); those were reverted before committing. The v2 tests render only into pytest temporary directories.

## Freeze verification

| Check | Start | End |
|---|---|---|
| `git rev-parse science-freeze-2026-08-15^{commit}` | `9a1b7777d6c69c2be44f16f25bc950769dc2ffda` | `9a1b7777d6c69c2be44f16f25bc950769dc2ffda` |
| `git diff science-freeze-2026-08-15 -- results/tables/evidence_freeze/` | empty | empty |
| `poster/`, root `README.md`, `PREANALYSIS.md`, `docs/` changed on this branch (vs base `477c992`) | no | no |

(As recorded for v1: `git diff science-freeze-2026-08-15 -- poster/` cannot
be empty on any branch because `poster/` post-dates the freeze tag; the
meaningful check is `poster/` unchanged relative to the branch base.)
