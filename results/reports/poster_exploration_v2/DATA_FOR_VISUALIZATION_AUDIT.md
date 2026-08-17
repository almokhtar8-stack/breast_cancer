# Data-for-Visualization Audit (poster-exploration-v2)

Inventory of every REAL, plot-level (not just summary-statistic) data
source available for the new exploratory figure bank, without recomputing
any frozen biological result. Every figure built in this phase cites one
of the rows below as its source. Where a concept in the prompt turned out
to be impossible from frozen data, that is stated explicitly rather than
faked.

Two classes of transformation are used below and are NOT considered "new
science": (a) reading an already-frozen table/matrix and reshaping it for
plotting (selecting rows/columns, pivoting, melting); (b) a disclosed,
standard, deterministic normalization applied for visualization only
(e.g. CPM = raw_count / library_size * 1e6, log2(CPM+1)) that does not
touch any p-value, FDR, effect size, or significance call anywhere in the
project. Both classes were already used and documented in
`src/poster_figures_bank_data.py` and `src/post_audit_sensitivity_data.py`
in earlier phases; this phase follows the same convention.

## CRISPR (genome-wide)

| Need | Source | Unit | N |
|---|---|---|---|
| All fitted genes | `data/processed/labels.parquet` | gene | 19,103 |
| 28 Gate-1 hits (both directions) | same, filtered `fdr<0.1` | gene | 28 |
| 13 significant sensitising hits + ranks | `post_audit_sensitivity_data.load_significant_sensitising_hits()` | gene | 13 |
| Blind control (not recovered) | same table, gene=RCOR1 | gene | 1 |

Fully available. No gap.

## GSE118713 (MCF7 / TAMR / FASR, resistance-model cell line)

| Need | Source | Unit | N |
|---|---|---|---|
| Per-sample TPM, any gene | `data/processed/gse118713_gene_tpm.parquet` | sample | 9 (3 MCF7 + 3 TAMR + 3 FASR, replicate aliquots of ONE derivation event) |
| Validated DE (log2FC/FDR) | `results/tables/gse118713_differential_expression.tsv.gz` | gene x contrast | genome-wide |
| Sample QC / PCA | `results/tables/gse118713_pca_coordinates.tsv` | sample | 9 |

Fully available for all 4 focus genes (KDM1A, TLK2, USP34, VEZF1), not just
USP34 -- the frozen TPM matrix is genome-wide. **Caveat that must stay
visible in any figure:** this is ONE parental line -> one TAMR-derivation
event + one FASR-derivation event, each with 3 replicate RNA-seq samples
(technical/aliquot replicates of a single derivation, not 3 independent
derivations) -- never implied as population replication.

## GSE111151 (4 parental lines x 7 independent TamR sublines, resistance-model)

| Need | Source | Unit | N |
|---|---|---|---|
| Per-sample log2CPM, any gene | `data/processed/gse111151_log2cpm.parquet` | sample | 11 (4 parental + 7 resistant) |
| Explicit pairing | `results/tables/gse111151_sample_metadata.tsv` (`paired_parental_sample_id` column) | sample | 11 |
| Validated per-gene DE | `results/tables/systems_network/gse111151_ranked_genes.tsv` (ranking stat, log2FC, FDR) | gene | genome-wide |

Fully available, genome-wide, all 4 focus genes, with REAL parental->
derivative pairing already encoded in the metadata (MCF-7->MCF-7_Tam1;
T-47D->{Tam1,Tam2}; ZR-75-1->{Tam1,Tam2}; BT-474->{Tam1,Tam2}). This
supports a genuine slopegraph/trajectory plot using the actual blocked
design -- no invented pairing needed. BT-474 is HER2-amplified, a
different molecular subtype from the other three -- kept visible, not
hidden.

## GSE240112 (3 primary vs 3 recurrent human tumours, UNPAIRED, recurrence-associated)

| Need | Source | Unit | N |
|---|---|---|---|
| Per-tumour log2CPM, all 13 sensitising genes (incl. KDM1A, TLK2) | `results/tables/gse240112/candidate_sample_level_log2cpm.tsv` | tumour (pseudobulk) | 6 (3 PT + 3 RT) |
| Genome-wide per-tumour raw counts | `results/tables/gse240112_pseudobulk/tumor_cell_counts.tsv.gz` + `tumor_cell_metadata.tsv` | tumour | 6 |
| Validated genome-wide DE | `results/tables/gse240112_pseudobulk/tumor_cell_genomewide_de.tsv.gz` | gene | genome-wide |

Fully available for all 4 focus genes at true per-tumour resolution
(`candidate_sample_level_log2cpm.tsv` already includes KDM1A and TLK2, not
only the original 4 -- a frozen table from the post-audit evidence-matrix
build). **Caveat that must stay visible:** PT and RT samples are from
different source institutions/biobanks, unpaired, n=3 vs n=3 -- no
connecting line may be drawn between any PT and RT observation.

## GSE245601 (single-cell, acute 12h ex vivo tamoxifen, patient-matched Control/Tamoxifen)

| Need | Source | Unit | N |
|---|---|---|---|
| Paired per-patient pseudobulk, original 4 candidates | `results/tables/gse245601_candidate_deepdive/patient_malignant_pseudobulk.tsv` | patient (pseudobulk) | 10 patients present, but only 3 pass the pre-declared eligibility filter |
| Eligibility filter | `results/tables/gse245601_pair_eligibility.tsv` (`eligible_for_pseudobulk`) | patient | 3 eligible: Tumor_02, Tumor_03, Tumor_07 |
| Genome-wide raw counts, eligible patients only | `results/tables/gse245601_pseudobulk/track_b_malignant_counts.tsv.gz` + `track_b_malignant_metadata.tsv` (library sizes) | patient x condition | 6 (3 patients x Control/Tamoxifen) |
| Validated genome-wide DE | `results/tables/gse245601_pseudobulk/track_b_genomewide_de.tsv.gz` | gene | genome-wide |

KDM1A and TLK2 do NOT appear in the frozen `patient_malignant_pseudobulk.tsv`
(built only for the original 4). However, the genome-wide raw-count matrix
IS available for the same 3 eligible patients, so KDM1A/TLK2 log2(CPM+1)
values are computed here using the IDENTICAL formula already used to
produce the frozen `normalized_expression` column for the other 4 genes
(log2((raw_umi/library_size*1e6)+1) -- verified by reproducing USP34's
already-frozen Tumor_02 value from the raw counts before trusting this
for KDM1A/TLK2). This is a disclosed, formula-identical extension of an
existing frozen computation to 2 additional genes for plotting only -- no
DE/FDR/significance value is touched. **Caveat that must stay visible:**
this is real patient-matched pairing (same tumour, Control vs Tamoxifen,
12h ex vivo) -- connecting lines ARE appropriate here, unlike GSE240112.
n=3 patients only.

## Pathway biology (GSEA)

| Need | Source | Unit | N |
|---|---|---|---|
| NES/FDR per pathway per dataset | `results/tables/systems_network/gsea_{dataset}.tsv` | pathway x dataset | 4 datasets, hallmark + go_bp collections |
| Full pre-ranked gene list (for real running-ES curves) | `results/tables/systems_network/{dataset}_ranked_genes.tsv` (`ranking_stat` column) | gene | genome-wide per dataset |
| Full gene-set membership (for real running-ES curves) | `data/reference/genesets/hallmark.gmt` (already-frozen, locally cached MSigDB Hallmark collection) | gene set | 50 hallmark sets |

**This makes real GSEA running-enrichment-score curves buildable** (Section
C3): the standard weighted Kolmogorov-Smirnov-style running-sum statistic
can be reconstructed deterministically from the already-frozen rank list +
gene-set membership, reproducing the same curve shape the frozen NES/FDR
already summarize -- this recomputes a VISUALIZATION of an existing
statistic, not a new enrichment test, and the frozen NES/FDR values
displayed alongside it are read directly from `gsea_{dataset}.tsv`, never
recomputed.

## Network / systems biology

| Need | Source | Coverage |
|---|---|---|
| Candidate-pathway membership | `results/tables/systems_network/candidate_pathway_membership.tsv` | **Original 4 candidates only** (USP34, VEZF1, EML5, CITED2) |
| Candidate-candidate connections | `results/tables/systems_network/candidate_candidate_connections.tsv`, `four_candidate_convergence.tsv` | **Original 4 candidates only** |
| Direct STRING neighbors | `results/tables/systems_network/four_candidate_direct_neighbors.tsv` | **Original 4 candidates**: CITED2 18 rows, USP34 10 rows, VEZF1 1 row, EML5 1 row (VEZF1/EML5 too sparse for their own network figure) |

**Important gap, honestly documented:** the systems-network phase of this
project was run ONLY on the original 4-candidate set (USP34/VEZF1/EML5/
CITED2), before KDM1A/TLK2 were reconsidered post-audit. There is NO
frozen network table for KDM1A or TLK2 anywhere in this repository. Any
Section D network figure can therefore only show USP34/VEZF1 (the two
still-current focus genes with network data) plus EML5/CITED2 for context
-- it cannot be a "4 focus gene" figure. This is disclosed in the figure
guide, not hidden by silently substituting different genes.

## Human tumour / DepMap

| Need | Source | Unit | N |
|---|---|---|---|
| TCGA-BRCA paired tumour-vs-normal | `results/tables/independent_validation/TCGA_candidate_expression.tsv` | patient (paired) | **Original 4 candidates only** -- KDM1A/TLK2 never assessed |
| DepMap 26Q1 per-line Chronos, all 4 focus genes | DepMap `CRISPRGeneEffect.csv` raw matrix, generic symbol match (same method as `poster_final_data.load_f4_depmap_effect`) | ER+/luminal cell line | 11 |
| Cell-line identity | DepMap `Model.csv` (`StrippedCellLineName` / `CellLineName`) | cell line | 11 |

Fully available for the DepMap side (real per-line values + real model
names, not just a summary %). TCGA remains an original-4-only asset, same
gap already disclosed in the science-freeze and poster-final phases.

## Structural biology

| Gene | Structure | Local file | Ligand/probe |
|---|---|---|---|
| KDM1A | 6NQU (LSD1 + GSK2879552) | fetched this phase, `kdm1a_tlk2_structures_dir` (new config key, RCSB download, see `scripts/download/download_kdm1a_tlk2_structures.py`) | KWM (inhibitor, genuinely bound) |
| TLK2 | 5O0Y (kinase domain) | fetched this phase, same dir | AGS (ATP-gamma-S, NOT an inhibitor) |
| USP34 | 7W3R (apo) / 7W3U (probe-bound) | already frozen, `usp34_structures_dir` | AYE (covalent ubiquitin-propargylamide probe, NOT a small-molecule) |
| VEZF1 | none exists | -- | -- (only a published homology model + one weak IC50=20uM screening hit, already documented in Table 06b -- not rendered as a fabricated "structure") |

All three experimentally-solved structures are now locally available for
rendering in a consistent visual style. VEZF1 is intentionally NOT given a
structural panel with fabricated content -- Section G represents it as
"no experimental structure" rather than substituting a homology model or
AlphaFold prediction as if it were equivalent evidence.

## What was skipped, and why

- **scRNA-seq UMAP/cluster-level plots for GSE245601 malignant-cell
  composition**: real cluster tables exist (`gse245601_cluster_summary.tsv`,
  `malignant_enrichment.tsv`) but were judged out of scope for THIS figure
  bank (candidate-expression focus, not cell-typing methodology) --
  candidate for a future phase, not built here.
- **KDM1A/TLK2 network figures**: impossible from frozen data (see Network
  section above) -- not fabricated by re-running network analysis with a
  different gene list.
- **A true GSEA enrichment curve for GSE245601 (acute)**: the acute
  dataset's own `gsea_gse245601.tsv` / `gse245601_ranked_genes.tsv` exist
  and ARE used, but Section C prioritizes Estrogen Response Early/Late and
  EMT curves in the three non-acute contexts plus one acute contrast for
  comparison, not all 4x50 combinations (kept to a readable number of
  curves).
