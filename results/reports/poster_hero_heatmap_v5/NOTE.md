# Hero Heatmap v5 -- Short Data Note

**Why v5 exists:** v4 established the correct concept (every row a real
biological observation, gene-wise within-dataset z-score) but read as a
technical genomics heatmap rather than a polished poster centerpiece. v5
changes PRESENTATION ONLY -- it imports `DATASET_BUILDERS` (the row
construction functions) and `DIVERGING_CMAP`/`STATE_COLORS` directly from
`src/poster_hero_heatmap_v4.py`, so the 29 biological rows, the 4 genes,
and every z-score value are byte-identical to v4. Nothing was
recomputed, reordered for aesthetics, or re-derived.

**What changed (presentation only):**
- Each dataset block is now headlined by its biological context, not its
  accession: "Cell-line resistance model" / "Independent
  tamoxifen-resistant sublines" / "Primary vs recurrent tumours" /
  "Single-cell: before vs 12 h after tamoxifen". The GEO accession is
  kept as small, muted, italic secondary text directly beneath the
  title -- present for traceability, no longer the headline.
- Block headers are attached directly above their own rows (with a thin
  rule beneath), instead of floating in a vertically-centered side
  column, so each header reads as part of its block.
- Landscape canvas (~4170x3010 px, ~4:3) so the heatmap body dominates
  the figure; the 4-gene heatmap occupies ~63-65% of the figure width.
- One compact legend (a single diverging colorbar, "Low / 0 / High") plus
  one small 4-item state-annotation key, replacing v4's larger two-part
  legend footer.
- GSE118713 and GSE240112 (the two datasets with no real sample-to-sample
  pairing) use a plain thin dashed divider between their two conditions;
  GSE111151 and GSE245601 (genuine parent->derivative / patient-matched
  pairs) keep a bracket, redrawn tighter so it never crosses row-label
  text.

**Display limits / clipping:** no clipping was applied. The maximum
|z-score| across all 116 gene x sample observations is 1.90, already
inside +/-2, so a +/-2 (or +/-2.5) clip would not move a single cell's
color -- clipping is therefore not justified by the actual values and
was omitted rather than silently applied. The shared diverging color
scale spans exactly [-1.90, +1.90].

**Confirmation no scientific result changed:** the row set, the gene set,
the pairing structure, and every z-score value are identical to v4 (see
`tests/test_poster_hero_heatmap_v5.py::test_zscore_matrix_identical_to_v4`).
No CRISPR, DE, FDR, DepMap, TCGA, pathway, network, or structural
conclusion is shown in or touched by this figure.
