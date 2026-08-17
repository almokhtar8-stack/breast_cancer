# Hero Heatmap v6 -- Short Data Note

**Why v6 exists:** v5 fixed hierarchy and orientation but still had far
too much unused white space around a comparatively small heatmap body.
v6 changes GEOMETRY ONLY -- it imports `DATASET_BUILDERS`,
`DIVERGING_CMAP`, `STATE_COLORS`, `CONTEXT_TITLE`, and `STATE_LABELS`
directly from `poster_hero_heatmap_v4`/`v5`, so the 29 biological rows,
the 4 genes, gene order, row order, and every z-score value are
byte-identical to v4/v5. Nothing about the data or the transform changed.

**What changed (geometry only):**
- Gene columns widened (`col_w` 1.00 -> 1.28) and rows made taller
  (`row_h` 0.30 -> 0.40) -- cells are visibly larger, substantial
  rectangles rather than thin strips.
- Left-side label/bracket/annotation zone compressed (bracket and label
  positions pulled in, subgroup/block gaps shrunk) so the heatmap now
  occupies ~74% of the plotted content width, up from ~63-65% in v5.
- Matplotlib's default ~12% axes margin was reclaimed via
  `fig.subplots_adjust(left=0.005, right=0.995, top=0.995, bottom=0.005)`
  so the canvas is sized to the actual content, not a default empty
  border.
- A real bug was caught and fixed during rendering: at the original title
  fontsize (21pt) the title text was measured at 13.7 inches wide --
  wider than the whole canvas -- which forced `bbox_inches='tight'` to
  expand the saved image far to the right of the actual heatmap, leaving
  a large dead zone. Fixed by measuring actual rendered text width and
  capping the title at 18.5pt (fits the reclaimed axes width with
  margin) rather than guessing a size.
- Legend condensed to one horizontal band (colorbar + two small text
  lines on the left, one compact 4-item state key on the right); its
  total height is well under 10% of the figure height.
- Landscape canvas, final image ~3988x2690 px (~1.48:1, within the
  requested 4:3-3:2 range).

**Confirmation the data matrix is identical to v5/v4:** every row's
z-score dict is produced by calling `poster_hero_heatmap_v4.DATASET_BUILDERS`
directly (see `tests/test_poster_hero_heatmap_v6.py::test_zscore_matrix_identical_to_v5`).
No value, row order, gene order, or pairing relationship was changed.
GSE240112 remains unpaired (no bracket, group divider + labels only);
GSE111151 and GSE245601 keep their genuine parent/patient-matched
brackets, redrawn more compactly. No clipping is applied to the color
scale (max |z| = 1.90 < 2, as documented in the v5 note).

**Confirmation no scientific result changed:** no CRISPR, DE, FDR,
DepMap, TCGA, pathway, network, or structural conclusion is shown in or
touched by this figure.
