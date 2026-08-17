# Hero Heatmap v3 -- Short Data Note

**Exact frozen data inputs** (reused unchanged via
`poster_story_v1_data.build_hero_heatmap_pairs()` -- no new computation,
same values already used in `poster_hero_heatmap` v1 and v2):
- GSE118713: `data/processed/gse118713_gene_tpm.parquet` (MCF7 vs TAMR, TPM)
- GSE111151: `data/processed/gse111151_log2cpm.parquet` (parental vs resistant, log2CPM)
- GSE240112: `results/tables/gse240112/candidate_sample_level_log2cpm.tsv` (primary vs recurrent, log2CPM)
- GSE245601: frozen raw counts + library sizes, log2(CPM+1) (control vs tamoxifen 12h)

**Reference + comparison rows (rows 1-2 of each block):** colored with the
SAME sequential palette (pale mist -> deep taupe) on ONE `Normalize` fit
jointly to the 8 values in that block (4 genes x {reference, comparison}).
Scaling the two rows together, instead of separately as in v2, means a
color difference between row 1 and row 2 for a given gene is now a real,
visually honest signal of that gene's own before/after shift -- not an
artifact of two independently-stretched scales.

**Delta strip (row 3 of each block):** log2 fold-change of the comparison
condition vs. that block's own reference mean -- the exact same `log2fc`
column already frozen in `build_hero_heatmap_pairs()`, unchanged from v1/v2.
Colored with a DIVERGING palette (slate blue = down, terracotta = up) on
ONE shared `TwoSlopeNorm` scale across all four dataset blocks.

**Why expression colors are only within-dataset comparable:** GSE118713 is
TPM, GSE111151/GSE240112 are log2CPM, and GSE245601 is computed
log2(CPM+1) -- these are different platforms and units, not on a shared
scale. Each block's sequential norm is fit only from that block's own 8
cells, and the legend caption says so explicitly ("scaled WITHIN each
dataset block"); no comparison of shading is implied or valid across the
four blocks.

**Why delta colors are cross-context comparable:** log2 fold-change is a
ratio, not a raw expression value, so it is the standard comparable unit
across independent studies (the same convention already used throughout
this project's frozen cross-dataset evidence tables). All four blocks'
delta cells sit on one shared diverging norm, so color intensity is
directly comparable block to block -- e.g. GSE245601's uniformly stronger
blue and GSE240112's standout VEZF1 red are real, unexaggerated signal.

**Confirmation no scientific result changed:** no new statistical test, no
recomputed p-value/FDR, no altered log2FC or expression mean. This figure
only changes how the same three frozen quantities (reference mean,
comparison mean, log2FC) are laid out and colored. TLK2 remains negative
in all four contexts, GSE245601 remains negative across all four genes,
USP34 remains positive in GSE118713/GSE240112 but not universally, VEZF1's
strongest positive shift remains in GSE240112, and GSE111151 remains
mostly weak/null -- exactly as in the frozen source data.
