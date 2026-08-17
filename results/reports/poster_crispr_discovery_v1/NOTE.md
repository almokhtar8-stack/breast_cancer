# CRISPR Discovery Figure v1 -- Short Data Note

**Purpose:** the poster's entry/discovery figure -- answers "which genes
sensitize ER+ breast-cancer cells to tamoxifen, and where do the four
focus candidates sit within that discovery result?" before the
transcriptomic heatmap.

**Exact frozen data sources** (read unmodified, no new computation):
- `data/processed/labels.parquet` -- the genome-wide CRISPR fit (19,103
  fitted genes), loaded via `post_audit_sensitivity_data.load_genomewide_crispr()`.
- The same module's `load_significant_sensitising_hits()`, which applies
  the pre-specified, already-frozen significance gate (FDR < 0.10,
  PREANALYSIS.md Section 4) and the sensitising direction (effect_size <
  0), then sorts by effect size. This is the exact function already used
  and tested elsewhere in this project's post-audit sensitivity analysis
  -- nothing here re-derives, re-thresholds, or re-ranks it.

**Every number on the figure is computed at render time from these
tables, not hand-typed:** the "13 genes" and "19,103 genome-wide fitted
genes" counts are `len(hits)` and `len(genomewide)`; every stem length
and marker position is the loaded `effect_size` column; gene order is the
loaded `rank_by_effect` order.

**What the figure shows:** a horizontal ranked lollipop of all 13
significant sensitising hits, most negative (strongest sensitising
knockout) at the top. The four focus candidates (KDM1A, TLK2, USP34,
VEZF1) are drawn in their established poster colors and bold labels; the
other 9 significant hits are muted gray. This makes the real discovery
context legible at a glance: KDM1A (rank 1) and TLK2 (rank 4) are among
the strongest functional hits in the entire significant set; VEZF1 (rank
8) and USP34 (rank 12) are comparatively weaker by CRISPR effect size
alone but are genuine members of the same pre-specified significant
sensitising set, not later additions.

**Why a single ranked lollipop over alternatives:** a volcano or a full
19,103-gene landscape would either bury the 13 real hits in a dense,
mostly-irrelevant cloud or require heavy annotation to point out the four
focus genes -- both work against "legible at poster distance." The
ranked lollipop restricted to the 13-hit significant set puts every
labeled gene at a defined, readable row, keeps the four focus genes
visually prominent without a busy multi-panel layout, and directly
answers the five questions this figure needs to answer. No alternate
figure was produced -- the main plot already satisfies the figure's
selection standard (CRISPR discovery result, 13 hits, relative strength
of the four candidates) on its own, so a second panel would only add
scope without adding information.

**Confirmation no scientific result changed:** no CRISPR effect size,
FDR, significance gate, or ranking was recomputed, re-thresholded, or
re-ordered. The figure reuses `post_audit_sensitivity_data`'s already-frozen
loaders unmodified.
