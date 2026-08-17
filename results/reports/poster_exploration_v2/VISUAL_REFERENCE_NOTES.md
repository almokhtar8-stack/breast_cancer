# Visual Reference Notes (poster-exploration-v2)

Brief web research was performed before building this figure bank (two
`WebSearch` queries; general web search, not full-text access to specific
paywalled figures, so no specific published figure is copied or even
directly viewed here -- principles only, general and already
well-established genomics-visualization conventions, corroborated by the
sources below).

Sources consulted:
- [A novel pipeline for prioritizing cancer type-specific therapeutic vulnerabilities using DepMap (PMC10850805)](https://pmc.ncbi.nlm.nih.gov/articles/PMC10850805/)
- [A dependency map enhanced with next-generation 3D cancer models, Nature](https://www.nature.com/articles/s41586-026-10843-7)
- [Building a translational cancer dependency map for TCGA, Nature Cancer](https://www.nature.com/articles/s43018-024-00789-y)
- [VolcaNoseR: a web app for creating, exploring, labeling and sharing volcano plots, Scientific Reports](https://www.nature.com/articles/s41598-020-76603-3)
- [Visualization of functional enrichment results, clusterProfiler/enrichplot book](https://yulab-smu.top/biomedical-knowledge-mining-book/enrichplot.html)
- [GSEA running-sum enrichment-score methodology, Bioconductor discussion](https://support.bioconductor.org/p/83814/)

## Principles extracted (not copied from any single figure)

**Density and restraint.** Modern cancer-genomics figures (DepMap-style
dependency papers, CRISPR-screen papers) favor dense, information-rich
single panels over decorative summary boxes. A genome-wide screen is
almost always shown as one scatter/volcano with tens of thousands of
points rendered small and translucent, with only a handful of named hits
labeled directly on the plot -- never as a table pretending to be a
figure.

**Direct labeling over legends.** Gene/hit labels are placed next to their
point with a short leader line or none at all, in the same color as the
point. Large boxed legends are avoided when there are 4 or fewer
categories -- color + direct label carries identity.

**One or two accent colors against neutral gray.** The overwhelming
majority of a genome-wide or gene-set plot is neutral gray (background
genes, non-significant tests, other cell lines); accent color is reserved
for the specific genes/hits the figure is about. This is the strongest,
most consistently reproduced convention across the sources above and is
the direct reason the previous poster-final figure set (with saturated
fills and colored card backgrounds on every panel) read as a dashboard
rather than a genomics figure.

**Real biological units, not summary shapes.** DepMap-style dependency
figures plot actual per-cell-line values (heatmaps with real model names,
or strip/violin plots of real Chronos scores) rather than a single "% of
lines dependent" bar. CRISPR screen papers show the real rank/effect
distribution, not just a table of the top hits. This is the central
critique this exploration phase is responding to.

**GSEA enrichment curves are the classic running-sum (weighted
Kolmogorov-Smirnov) statistic**: walk down the ranked gene list, step up
at each gene inside the set, step down (weighted by set size) at each gene
outside it; the curve's peak deviation from zero is the enrichment score.
This is a standard, fully reconstructable statistic from a ranked gene
list and a gene set -- confirmed by the Bioconductor/enrichplot references
above -- and is used here (Section C3) exactly as intended: a
visualization of an already-computed statistic, not a new test.

**UpSet plots over Venn diagrams for >3 sets.** Confirmed as the standard
modern choice for showing multi-set overlap membership (e.g. "which of
these 13 genes have RNA support AND high dependency AND a structure")
without the combinatorial unreadability of a multi-circle Venn diagram.

**Structural biology panels favor consistent camera philosophy across
multiple targets** -- same background, same ray-tracing/lighting
approach, same relative scale logic -- specifically so a reader can
visually compare tractability across targets in one glance, rather than
each structure being rendered as an isolated illustration. This directly
informs Section G's approach (KDM1A/TLK2/USP34 rendered with matched
style).

## What this means concretely for this figure bank

- Backgrounds: white. Non-focal data: light neutral gray, never colored.
- Candidate accent colors used ONLY where gene identity is the point of
  the panel (kept from the previous phase's Okabe-Ito-derived palette,
  but never as a filled card/box background).
- No panel uses a tinted rectangle, colored dashboard card, or scorecard
  matrix as its primary content.
- Real per-sample/per-patient/per-cell-line values are shown wherever the
  data audit confirms they exist (see `DATA_FOR_VISUALIZATION_AUDIT.md`),
  in preference to any pre-aggregated percentage or summary statistic.
