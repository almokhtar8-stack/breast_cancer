# results/figures

**Looking for the poster? Go to [`poster/final_figures/`](../../poster/final_figures/).**
That directory holds the six canonical figures with a provenance manifest and
SHA-256 hashes.

This directory is the **working figure archive**: one subdirectory per analysis
phase, including superseded development iterations. Nothing here is deleted, so
that every published figure can be traced back through its own version history.

## Which directories back the canonical poster figures

| Poster figure | Source directory here |
|---|---|
| 01 CRISPR discovery | `poster_crispr_discovery_v1/` |
| 02 Candidate expression | `poster_hero_heatmap_v6/` |
| 03 Molecular networks | `poster_network_mechanism_v4/` |
| 04 Pathway remodeling | `poster_pathway_v2/` |
| 05 Baseline dependency | `poster_depmap_v2/` |
| 06 Structural tractability | `poster_druggability_v1/` (+ `renders/` PyMOL images) |

## Superseded development iterations (retained for provenance)

- `poster_hero_heatmap/`, `_v2/` … `_v5/` — earlier expression-heatmap designs.
  `v3` used group means rather than real samples and was rejected; `v4` introduced
  the true sample-level matrix; `v5`/`v6` are presentation-only refinements.
- `poster_network_mechanism_v1/` … `_v3/` — `v1` reused an older, asymmetric
  network built for a different candidate shortlist; `v2` is the standardized
  STRING rebuild; `v3` a compact single-canvas layout; `v4` is canonical.
- `poster_pathway_v1/` — earlier three-panel pathway figure, superseded by `v2`.
- `poster_depmap_v1/` — the full 11 × 4 Chronos heatmap. Scientifically complete
  and still useful for per-cell-line detail; `v2` is the poster version.
- `poster_exploration_v2/`, `_v3/`, `poster_story_v1/`, `poster_final/`,
  `nebula/` — earlier exploratory figure banks.

Other subdirectories correspond to upstream analysis phases (QC, CNV method
audits, per-dataset diagnostics) and are documented in
[`docs/RESULTS_GUIDE.md`](../../docs/RESULTS_GUIDE.md).

Some large superseded PDFs are gitignored (see `.gitignore`); regenerate them
with the modules listed in [`docs/analysis_map.md`](../../docs/analysis_map.md).
