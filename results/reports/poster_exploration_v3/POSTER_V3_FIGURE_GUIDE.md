# Poster-Grade Figure Bank v3 -- Figure Guide

Six main figures + two alternates, rebuilt from the same frozen science
already used in v2, with a stricter poster-grade visual design system.
v2 was judged scientifically sound but visually too close to an internal
analysis/report deck (too many small subplots, too much grey, dense
footnotes, weak visual hierarchy). v3 does not change a single scientific
number -- it changes how those numbers are shown.

Every figure exports `.png` (300dpi) + `.pdf` + `.svg`, except
`F5_depmap_context`, `F6_structural_comparison`, and `ALT_B_pocket_closeups`
(PNG-only, since their content includes a rasterized PyMOL render or a
composited raster image -- a vector re-export would carry no additional
information). See `src/poster_exploration_v3_data.py` for the exact frozen
source of every plotted number.

## Style system actually applied

- White background throughout; no grey panel fills anywhere.
- Large type: hero titles ~19-22pt bold, panel titles ~15-17pt bold, axis
  labels ~14-16pt, tick labels ~12-13pt.
- At most one short takeaway line per figure -- no paragraph footnotes.
- 1-2 panels per figure, never a small-multiple grid larger than 4.
- Minimal gridlines; light, recessive spines; direct labeling over legend
  boxes wherever practical.
- Fixed candidate palette everywhere: KDM1A = orange (`#D55E00`), TLK2 =
  light blue (`#56B4E9`), USP34 = strong blue (`#0072B2`), VEZF1 = gold
  (`#E69F00`) -- identical hex values already used and colorblind-safety-
  validated in v2/poster_final, not redefined here.
- Every number on every figure is read from a loaded DataFrame/Series at
  render time -- see `tests/test_poster_exploration_v3.py` for the checks
  that pin this (including a direct grep guard against a hand-typed
  "19,103" string).

---

## The 6 main figures

### F1. `F1_crispr_discovery` -- CRISPR discovery
**Selected because:** it is the cleanest possible answer to "where do the
4 focus genes sit among the significant hits" -- a single large lollipop,
sorted by effect size, focus genes bold and colored, the other 9 muted
gray. Built from v2's A2 (13-hit lollipop), redesigned at poster scale: 2x
larger fonts, thicker stems, no small subplots. A small "19,103 genes
screened" annotation gives genome-wide context in one line rather than a
second panel, per the brief's "may incorporate a small amount of A1
context if useful."
**One takeaway:** 13 significant sensitising hits; the 4 focus genes span
the full range of functional strength -- KDM1A/TLK2 far stronger than
USP34/VEZF1 by this measure alone.

### F2. `F2_transcriptomic_corroboration` -- Transcriptomic corroboration
**Selected because:** it replaces v2's 4x4 = 16-panel micro-multiple dump
(B1-B4 combined) with ONE elegant 2-panel figure. All 4 genes appear on a
SHARED, comparable y-axis in each panel by using a disclosed centering
transform (delta relative to each gene's own parental line for GSE111151;
delta relative to each gene's own primary-tumour-group mean for
GSE240112) -- see `build_gse111151_delta_from_parental` /
`build_gse240112_delta_from_primary_mean` in `poster_exploration_v3_data.py`
for the exact, disclosed subtraction (no new statistical test, no p-value
computed). Every point is a real sample/tumour observation, not an
aggregate.
**One takeaway:** GSE111151 is a resistance-model cell-line panel;
GSE240112 is unpaired human recurrence data -- VEZF1 shows the clearest
real shift in the recurrence context.

### F3. `F3_pathway_convergence` -- Pathway convergence
**Selected because:** it was already the strongest visual hero in v2 (C2)
and needed only scaling up, not redesigning: larger fonts, bolder line
weights, direct end-of-line labels instead of a legend box, a one-line
takeaway instead of a caveat paragraph.
**One takeaway:** Estrogen response falls in every transcriptomic context
tested; EMT rises in resistance and recurrence but falls in the acute 12h
context -- the clearest pathway-level divergence in the whole dataset.

### F4. `F4_postaudit_interpretation` -- Why these genes / post-audit interpretation
**Selected because:** it is the figure that makes "this is an
interpretation framework, not a naive top-CRISPR ranking" visually
undeniable. Rebuilt from v2's E3 with the same underlying Rule 0 (original
frozen gate) vs. Rule 1 (CRISPR-only rank) comparison, but decluttered
dramatically: the 9 non-focus, gate-excluded genes get ONE shared
annotation instead of 9 repeated "excluded by RNA-eligibility gate"
labels, and only the two genes that actually passed the original gate
(USP34, VEZF1) get a connecting line into the right-hand column.
**One takeaway:** KDM1A and TLK2 rank 1st and 2nd by CRISPR strength alone
yet were never eligible under the original RNA-corroboration gate --
USP34 and VEZF1 were selected on different grounds entirely.

### F5. `F5_depmap_context` -- Human / DepMap context
**Selected because:** a heatmap of real per-line values with real
cell-line names is inherently stronger than any percentage bar chart --
this figure keeps that real-data heatmap (from v2's F1) and adds one very
small secondary element (a per-gene median-value strip along the top) so
a viewer gets both the aggregate signal and the underlying heterogeneity
in one glance, without a second full panel.
**One takeaway:** TLK2 shows the strongest and most consistent baseline
dependency across all 11 real ER+/luminal lines -- explicitly flagged as
not automatically an advantage for tamoxifen-specific sensitisation.

### F6. `F6_structural_comparison` -- Structural / pharmacological comparison
**Selected because:** it is very likely the single strongest visual hero
in the entire bank. Real, experimentally solved structures for KDM1A
(6NQU), TLK2 (5O0Y), and USP34 (7W3U), fetched and rendered in a matched
camera/lighting style, placed alongside an honest "no experimental
structure" panel for VEZF1 (a plain neutral box, never a fabricated
homology-model or AlphaFold render standing in for real data). Each
gene's pharmacological maturity is stated in one short line, not a
paragraph.
**One takeaway:** four candidates occupy four genuinely different points
on the structural/pharmacological maturity spectrum -- from KDM1A's
clinical-stage inhibitors to VEZF1's complete absence of a structure.

---

## Alternates (2)

### Alt A. `ALT_A_genomewide_landscape` -- alternate to F1
A genome-wide ranked-landscape view (all 19,103 genes) with a magnified
inset showing the 13 significant hits, instead of a lollipop. Kept as a
genuine alternate because there is real uncertainty about whether a
poster audience responds better to "here is the whole screen, zoomed in"
(this figure) or "here are the 13 hits, ranked" (F1) as the opening image.

### Alt B. `ALT_B_pocket_closeups` -- alternate to F6
The same three real structures, zoomed to the binding/catalytic site
instead of showing the whole protein. Kept as a genuine alternate because
whole-protein context (F6) and binding-site detail (this figure) trade off
differently at poster viewing distance -- F6 reads better from a few
meters away; this figure rewards close-up inspection.

---

## What was intentionally retired from v2

- All v2 small-multiple grids of 4+ tiny panels (B1, B3, B4, B5a/b, C1b)
  -- replaced by fewer, larger, denser-but-legible panels (F2 in
  particular replaces an effective 16-panel dump with 2 panels).
- v2's dense multi-line footnotes (NES-comparability caveats,
  denominator disclosures, etc.) -- the underlying facts are still true of
  v3's data (nothing scientific changed) but are no longer printed as
  paragraph text inside the poster figures; they remain documented in
  `results/reports/poster_exploration_v2/DATA_FOR_VISUALIZATION_AUDIT.md`
  and this guide.
- v2's UpSet plot (E1) and quantitative evidence map (E2) -- both real and
  informative, but judged report-grade rather than poster-grade at their
  original density; F4 replaces E3's core message more elegantly, and
  E2's message (sensitisation strength is not baseline dependency) is
  already implicit across F1/F5 without a dedicated panel.
- v2's network figures (D1, D2) -- not rebuilt for v3. They remain
  honestly weak for poster use (no KDM1A/TLK2 network data exists at all,
  a real coverage gap, not a redesign problem) and the brief explicitly
  said not to force a weak network figure into the main set.
- v2's H2 (experimental schematic) and A3/A4 (redundant CRISPR variants)
  -- not rebuilt; F1/Alt A already cover CRISPR discovery, and no
  translational-schematic figure was judged necessary for this smaller,
  denser-signal set.

## Likely visual heroes

**F6** (structural comparison) and **F3** (pathway convergence) are the
two strongest candidates for the single most memorable poster image --
both are bold, high-contrast, immediately legible in 5-10 seconds, and
carry real data. **F1** and **F5** are close behind as strong, clean
supporting heroes.
