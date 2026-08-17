# Hero Heatmap v4 -- Short Data Note

**Why v4 exists:** v1-v3 collapsed every dataset to one reference-mean row
and one comparison-mean row per gene, which reads as a summary table, not
a genomics heatmap. v4 shows every real biological observation as its own
row, drawn from the same frozen per-sample loaders already used and
tested in `src/poster_exploration_v2_data.py` (`load_gse118713_focus_gene_samples`,
`load_gse111151_focus_gene_samples`, `load_gse240112_focus_gene_tumours`,
`load_gse245601_paired_focus_genes`). No new discovery, no new statistic,
no altered log2FC/FDR/p-value -- only a reshaping and a disclosed
visualization-only transform.

**Biological rows shown, per dataset (29 total):**
- GSE118713: 6 rows -- MCF7 replicates 1-3, TAMR replicates 1-3 (TPM).
  FASR (a fulvestrant-resistant derivative in the same series) is
  deliberately excluded -- this poster is scoped to tamoxifen
  sensitisation, and FASR does not bear on that story.
- GSE111151: 11 rows -- the real blocked design of 4 independent parental
  cell-line backgrounds (MCF-7, T-47D, ZR-75-1, BT-474) and their 7
  independently-derived TamR sublines (log2CPM). Parent -> derivative
  relationships are drawn exactly as encoded in the frozen
  `paired_parental_sample_id` metadata column; a line with two
  derivatives gets one bracket to both, never a fabricated 1:1 pairing.
- GSE240112: 6 rows -- Primary-1..3, Recurrent-1..3 tumour pseudobulk
  (log2CPM). These are UNPAIRED (different patients/biobanks): no bracket
  is drawn between primary and recurrent rows, only a group divider and
  label.
- GSE245601: 6 rows -- the 3 patients (Tumor_02/03/07) passing the
  project's own pre-declared pseudobulk eligibility filter, each with a
  genuinely patient-matched Control and 12h ex vivo Tamoxifen row
  (log2(CPM+1)). Bracketed within patient because the pairing is real.
  Labelled "acute tamoxifen (12h)", never called resistance.

**Visualization transform (the only transform applied):** for each
dataset and each gene independently, z = (value - dataset/gene mean) /
dataset/gene SD, computed across exactly the rows shown for that dataset
(6, 11, 6, and 6 observations respectively). This never touches a
p-value, FDR, effect size, or the frozen log2FC values used in v1-v3 --
it only re-expresses each gene's own within-dataset distribution as
unitless standard deviations. Because z-scores are dimensionless, one
shared diverging colorbar can honestly span all four dataset blocks
(unlike raw TPM/CPM/log2CPM, which are on different platforms/units and
are never plotted on one shared scale here). All N >= 6 with nonzero
variance in every block, so no degenerate (SD=0) case was encountered;
the code still guards for and would flag that case if it ever occurred.

**Row order / no clustering:** rows stay grouped dataset -> biological
condition -> sample/model, in the order encoded in the frozen metadata --
never hierarchically clustered or reordered to flatter a pattern. Gene
columns keep the fixed KDM1A / TLK2 / USP34 / VEZF1 order.

**Confirmation no scientific result changed:** no CRISPR, DE, FDR, DepMap,
TCGA, pathway, network, or structural conclusion was touched. This figure
only reshapes already-frozen per-sample expression values and applies the
disclosed z-score transform above for display purposes.
