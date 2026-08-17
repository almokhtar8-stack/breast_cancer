# Poster Layout Recommendation

This is a layout recommendation only -- **no poster has been assembled**.
It describes how the six figures in `results/figures/poster_final/` would
read as a physical/digital poster, for review before assembly is
authorized. See `POSTER_FINAL_FIGURE_GUIDE.md` for what each figure shows
and why it was selected.

## Recommended reading order and zone layout

A standard landscape scientific-poster grid (3 columns x 2 rows of
content, title band across the top) maps directly onto the six figures in
their numbered order -- no re-sequencing needed:

```
+-----------------------------------------------------------------+
|  TITLE: From a CRISPR Screen to Therapeutic Vulnerabilities:    |
|  Functional Genomics of Tamoxifen Sensitisation in ER+ Breast   |
|  Cancer                                    [affiliation / logo] |
+-----------------------------------------------------------------+
|   F1. Genome-wide   |   F2. Pathway        |   F3. Candidate    |
|   CRISPR discovery  |   systems view        |   evidence         |
|                      |                       |   integration      |
|   (landscape, wide)  |   (landscape, wide)   |   (2x2 grid,tall)  |
+-----------------------------------------------------------------+
|   F4. Human /        |   F5. USP34           |   F6. Final        |
|   DepMap validation   |   structure &         |   translational    |
|   (2-panel, wide)     |   tractability        |   framework        |
|                       |   (3-panel hero)      |   (2-row, wide)     |
+-----------------------------------------------------------------+
|  Caption strip / QR code to full repository / contact info       |
+-----------------------------------------------------------------+
```

**Reading order is left-to-right, top-to-bottom (F1->F2->F3->F4->F5->F6)**,
matching the project's own narrative arc: discovery, systems context, why
these four candidates, orthogonal validation, structural deep-dive on the
most novel lead, translational framework. This order requires no numbering
callouts beyond the existing "F1"-"F6" panel labels already baked into each
figure.

## Sizing guidance (relative, not absolute -- depends on final poster
canvas dimensions)

- **F1 and F2** are naturally wide (landscape) figures with dense
  right-side labeling (F1) or a wide dataset axis (F2, 4 columns). Give
  them the top-row's two wider slots.
- **F3** is a 2x2 grid and reads best as a tall, roughly-square panel --
  give it the top-row's narrower third slot, or promote it to full-width
  if poster space allows (it is the single most information-dense figure
  in the set and rewards standing-distance reading time).
- **F4** (two side-by-side panels) and **F6** (scatter + 4 role cards) are
  both naturally wide -- bottom-row's two wider slots.
- **F5** is the single most visually striking figure (real ray-traced
  structures) and is a strong candidate for a slightly larger "hero" slot
  if the poster template allows one panel to break the grid -- e.g. sized
  ~1.3x the others, since it is designed to be readable and eye-catching
  from several meters away.

## Typography and whitespace

- Every figure already carries its own panel letters (F1-F6, and A/B/C/D
  within multi-panel figures) at a poster-legible size (11-14pt bold at
  300dpi) -- no additional in-poster panel numbering is needed.
- Each figure's own italic caption strip (bottom of the PNG) is intended to
  be read at close range (poster-visitor distance, not across-the-hall
  distance) -- keep it visible, do not crop it off when placing the figure
  into a poster template.
- Leave a consistent gutter (recommend >=1.5% of poster width) between
  panels; none of the six figures have transparent backgrounds (all are
  flattened to white), so panels will show a visible white rectangle
  against any non-white poster background color -- either keep the poster
  background white/very light, or add a subtle drop shadow or thin border
  per panel so panel edges do not appear to float.

## What is deliberately NOT recommended

- Do not compress F3 or F6 below roughly 40% of poster width -- both rely
  on 4-column dot-matrices / role cards that become illegible at poster
  viewing distance below that size.
- Do not crop F5's caption strip -- the "no validated selective USP34
  inhibitor exists; docking was not pursued" sentence is load-bearing for
  not overclaiming druggability and must stay attached to the structural
  images.
- Do not add a 7th panel from the supporting/optional list
  (`POSTER_FINAL_FIGURE_GUIDE.md`) into the main grid without re-running
  this layout recommendation -- the six-figure grid above is sized
  assuming exactly six main panels.

## Open decision for the poster owner

F5 can be swapped for the alternative structural comparison
(`12b_USP34_structure_comparison`, in the retired-candidates contact sheet)
if a reviewer prefers a strict apo-vs-bound side-by-side over the current
apo / bound / close-up three-panel composition -- both are real PyMOL
renders of the same two frozen structures; this is a presentation
preference, not a scientific difference, and is left open rather than
decided here.
