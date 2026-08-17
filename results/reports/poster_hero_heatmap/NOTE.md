# Hero Heatmap -- Short Data Note

**Data sources** (all frozen, unmodified; reused via
`poster_story_v1_data.build_hero_heatmap_pairs()`, no new computation):
- GSE118713: `data/processed/gse118713_gene_tpm.parquet` (MCF7 vs TAMR, TPM)
- GSE111151: `data/processed/gse111151_log2cpm.parquet` (parental vs resistant, log2CPM)
- GSE240112: `results/tables/gse240112/candidate_sample_level_log2cpm.tsv` (primary vs recurrent, log2CPM)
- GSE245601: frozen raw counts + library sizes, log2(CPM+1) (control vs tamoxifen 12h)

**What the color encodes:** log2 fold-change of the comparison condition
(TAMR / resistant / recurrent / tamoxifen) relative to that same
dataset's own reference-condition mean (MCF7 / parental / primary /
control). The reference row itself is flat neutral gray -- not on the
color scale, since it IS the reference point, not a signal.

**Why this is valid across datasets:** GSE118713 is TPM, GSE111151/GSE240112
are log2CPM, and GSE245601 is computed log2(CPM+1) -- these raw units are
not on a shared scale, so plotting raw values with one colorbar would be
misleading. log2 fold-change is a ratio, not an absolute scale, and is the
standard comparable unit across independent studies (the same convention
already used throughout this project's frozen cross-dataset evidence
tables). No p-value or new statistical test was computed for this figure.
