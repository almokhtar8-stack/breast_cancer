# Hero Heatmap v2 -- Short Data Note

**Exact frozen data inputs** (reused unchanged via
`poster_story_v1_data.build_hero_heatmap_pairs()` -- no new computation,
same values already used in `poster_hero_heatmap` v1):
- GSE118713: `data/processed/gse118713_gene_tpm.parquet` (MCF7 vs TAMR, TPM)
- GSE111151: `data/processed/gse111151_log2cpm.parquet` (parental vs resistant, log2CPM)
- GSE240112: `results/tables/gse240112/candidate_sample_level_log2cpm.tsv` (primary vs recurrent, log2CPM)
- GSE245601: frozen raw counts + library sizes, log2(CPM+1) (control vs tamoxifen 12h)

**Top/reference row encodes:** that dataset's own reference-condition mean
expression (MCF7 / parental / primary / control), colored with a
SEQUENTIAL palette (pale mist -> deep taupe) scaled by min-max **within
that dataset block's 4 genes only**. This is real, data-bearing color --
not a placeholder -- but it is explicitly NOT comparable across the four
dataset blocks (different platforms/units), which is why it uses a
separate, differently-scaled legend from the comparison row.

**Bottom/comparison row encodes:** log2 fold-change of the comparison
condition (TAMR / resistant / recurrent / tamoxifen) vs. that same
dataset's own reference mean, colored with a DIVERGING palette (slate blue
= down, terracotta = up) on ONE shared scale across all four blocks. log2
fold-change is a ratio, not a raw scale, so this row -- and only this row
-- is the cross-context-comparable signal.

**Why this fixes the v1 design:** v1's reference row was flat neutral gray
in every cell, so it looked blank and wasted half the figure's visual
real estate. v2's reference row now carries real information (relative
baseline expression pattern within each context) while remaining honest
about what it is and is not comparable to -- fixing the "blank-looking"
problem without pretending raw values are comparable across studies. No
p-value, FDR, or biological result was recomputed or altered; TLK2 still
reads consistently down and the acute (GSE245601) context still reads
negative across all four genes, exactly as in the frozen source data.
