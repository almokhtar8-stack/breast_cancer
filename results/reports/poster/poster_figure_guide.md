# Poster Figure Guide

**This is a full rebuild of the poster figure set (see "What changed and
why" at the bottom).** The prior version of this figure set was
schematic/report-style (icon matrices, a dashboard table, a matplotlib
wireframe). This version replaces every main figure except the final
schematic with a real, data-driven scientific visual: an actual
genome-wide CRISPR volcano plot, actual per-sample expression points,
actual DepMap Chronos distributions, an actual pharmacogenomic scatter,
and a real PyMOL structural render.

All plotted numbers are loaded programmatically from frozen source tables
or already-validated per-sample/per-cell-line data matrices -- see
`src/poster_figures_data.py` and the pinned-value regression tests in
`tests/test_poster_figures.py`. No new analysis, no new candidate
discovery, no re-ranking.

---

## 1. `results/figures/poster/01_crispr_discovery.png`

- **Type:** genome-wide volcano plot (effect size vs -log10 FDR),
  19,103 genes.
- **Scientific message:** the four frozen candidates emerged from a real,
  genome-wide functional CRISPR screen, not a curated list -- all four are
  Gate-1 hits (FDR<0.1) with a sensitising-knockout (negative effect)
  direction under 4-OHT.
- **Source:** `data/processed/labels.parquet` (the genome-wide per-gene
  Hany screen fit, same file the Gate-1 decision itself was built from).
- **Why this is better than the previous version:** the prior poster set
  had no genome-wide discovery figure at all (the workflow diagram named
  the screen but showed no data). This shows the actual screen result for
  all 19,103 genes, with the four candidates highlighted in place among
  every other tested gene.

## 2. `results/figures/poster/02_expression_evidence.png`

- **Type:** 3-panel real-data figure -- (A) per-sample strip+box plot,
  (B) per-tumor strip+box plot, (C) cross-dataset forest/dot plot.
- **Scientific message:** USP34 is elevated in both TAMR and FASR
  resistant derivatives vs parental MCF7 (GSE118713); VEZF1 is elevated in
  recurrent vs primary tumors (GSE240112); across all four datasets and
  four candidates, the direction of effect is mostly consistent but not
  uniformly significant.
- **Source:** `data/processed/gse118713_gene_tpm.parquet` (9 real samples),
  `results/tables/gse118713_differential_expression.tsv.gz`;
  `results/tables/gse240112/candidate_sample_level_log2cpm.tsv` (6 real
  tumor pseudobulk samples), `results/tables/gse240112/candidate_table.tsv`;
  `results/tables/cross_dataset_genomewide/all_genes_cross_dataset_evidence_with_ranking.tsv`.
- **Why this is better than the previous version:** the prior version 2
  ("four-candidate evidence matrix") was a symbol/icon grid -- triangles
  and diamonds standing in for data. This version shows the actual sample
  points (every replicate, every tumor) that those symbols were
  summarizing.

## 3. `results/figures/poster/03_depmap_distributions.png`

- **Type:** box + jittered strip plot of real per-cell-line values.
- **Scientific message:** VEZF1 shows real baseline dependency in
  ER+/luminal breast cancer cell lines (several lines with strongly
  negative Chronos scores); USP34, EML5, and CITED2 do not.
- **Source:** real Chronos gene-effect values from DepMap 26Q1
  `CRISPRGeneEffect.csv`, loaded via the already-validated
  `independent_validation_depmap_data.load_gene_effect()` for the n=11
  ER+/luminal screened lines; percentages annotated from
  `results/tables/independent_validation/DepMap_candidate_dependency.tsv`.
- **Why this is better than the previous version:** the prior version 3
  was a bar chart of a single percentage per gene. This shows the actual
  11 per-line data points behind that percentage, including the two
  individual VEZF1 lines with strongly negative scores that drive it.

## 4. `results/figures/poster/04_pharmacogenomics.png`

- **Type:** scatter plot with OLS regression line + 95% CI.
- **Scientific message:** USP34 expression is negatively associated with
  AZD7762 (a CHK1/CHK2 inhibitor) response across GDSC breast cancer cell
  lines -- a real, FDR-significant pharmacogenomic correlation.
- **Source:** real per-cell-line USP34 expression (DepMap 26Q1) joined to
  GDSC1 AZD7762 response (DRUG_ID 1402), via the already-validated
  `final_pharmacogenomics_gdsc_data.build_breast_expression_joined()` --
  the single strongest, FDR-significant USP34 GDSC hit from that phase.
- **Included by judgment call:** this was the strongest and cleanest of
  the candidate GDSC associations (n=44 real cell lines, visible negative
  trend, FDR=0.008) and was judged compelling enough to include as a core
  figure rather than a demoted/omitted supplement.
- **Why this is new:** the prior poster set had no pharmacogenomic figure.

## 5. `results/figures/poster/05_structure.png`

- **Type:** real structural-biology figure -- four ray-traced PyMOL
  cartoon panels (not a matplotlib wireframe) + a small clean summary line.
- **Scientific message:** USP34 has a real experimental structure with a
  confirmed catalytic Cys1903/His2164 dyad and demonstrated covalent
  reactivity (7W3U's ubiquitin-propargylamide probe), but no selective
  inhibitor exists and docking has not been pursued.
- **Source:** the same two already-frozen, already-verified PDB structures
  (7W3R apo, 7W3U covalent-probe-bound) used in the final_translational
  phase, re-rendered here with PyMOL (`src/poster_structure_render.py`)
  purely for visualization -- no new structural analysis, pocket
  detection, or docking was performed. Chain-D-to-chain-A probe pairing
  and the apo/bound structural alignment (RMSD ~0.61 Å) were verified
  directly from the coordinates before rendering.
- **Why this is a complete redesign:** the previous version 5 used a raw
  matplotlib Cα-backbone wireframe ("spaghetti"), explicitly called out as
  unacceptable. This version uses PyMOL cartoon/ribbon rendering with
  properly colored catalytic residues and the covalent probe, matching
  standard structural-biology poster conventions.

## 6. `results/figures/poster/06_experimental_strategy.png`

- **Type:** experimental schematic (the one schematic figure in the set).
- **Scientific message:** the proposed USP34 validation experiment (4 arms
  x 3 readout groups + 2 normal-cell comparators) is designed to
  distinguish direct dependency, tamoxifen sensitisation, and dual action;
  VEZF1 appears only as a smaller secondary follow-up note.
- **Source:** `results/tables/final_translational/final_experimental_design.tsv`,
  `final_normal_cell_comparators.tsv`.
- **Why this is better than the previous version:** functionally similar
  to the prior version, but VEZF1's box is now deliberately smaller and
  visually secondary (an outlined note, not a filled box of equal weight)
  per the instruction that the lead USP34 experiment should visually
  dominate.

---

## What changed and why (summary)

| # | Prior design | New design | Reason |
|---|---|---|---|
| 1 | Workflow flow diagram (no data) | Genome-wide CRISPR volcano plot | Show the actual screen, not just name it |
| 2 | Icon/symbol evidence matrix | Real per-sample strip+box + forest plot | Replace symbols with actual data points |
| 3 | Single-percentage bar chart | Real per-line Chronos distribution | Show the 11 cell lines behind each percentage |
| 4 | (did not exist) | Real GDSC pharmacogenomic scatter | New, judged compelling and worth including |
| 5 | Matplotlib Cα wireframe ("spaghetti") | Real PyMOL cartoon/ribbon render | Explicitly required redesign |
| 6 | Two-arm-equal schematic | Same schematic, VEZF1 de-emphasized | Keep USP34 experiment visually dominant |

The previous "four-candidate evidence matrix," "independent validation
DepMap" bar chart, and "USP34 vs VEZF1 head-to-head dashboard" figures
have been retired from the core poster set -- their content is now folded
into the real-data figures above (2 and 3) or considered a supporting
supplement rather than a main poster panel, per the instruction that a
head-to-head table should not be a main panel.
