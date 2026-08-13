# Systems-network phase: genome-wide input audit

Phase 0-1 of the pathway/network systems-biology phase. Purely additive over
all previously frozen phases (evidence-freeze, GSE245601 candidate deep-dive).
Nothing under `results/tables/evidence_freeze/` or
`docs/THERAPEUTIC_SHORTLIST_FREEZE.md` is read-write here; both are read-only
reference points confirming the frozen four candidates
(USP34, VEZF1, EML5, CITED2) and are not touched.

## Repository state at start of this phase

- Branch: `main`, working tree clean.
- Evidence-freeze commit present: `4d085bf` ("Freeze therapeutic candidate
  evidence and five-layer summaries").
- GSE245601 deep-dive commits present: `33d848b`, `cd20ad4`.
- HEAD at start of this phase: `cd20ad4`.

## Genome-wide source tables located

Full detail in `results/tables/systems_network/input_inventory.tsv`. Summary:

| Dataset | File | Genes tested |
|---|---|---|
| GSE118713 (TAMR vs MCF7) | `results/tables/gse118713_differential_expression_unredacted.tsv.gz` | 14,838 |
| GSE240112 (tumor_cell, recurrent vs primary) | `results/tables/gse240112_pseudobulk/tumor_cell_genomewide_de.tsv.gz` | 18,429 |
| GSE111151 (resistant vs parental) | `results/tables/gse111151/genomewide_de.tsv.gz` | 27,418 |
| GSE245601 Track A (12h acute, all epithelial) | `results/tables/gse245601_pseudobulk/track_a_genomewide_de.tsv.gz` | 17,988 |
| GSE245601 Track B (12h acute, malignant, exploratory secondary) | `results/tables/gse245601_pseudobulk/track_b_genomewide_de.tsv.gz` | 13,864 |
| CRISPR (E2+4-OHT vs E2) | `data/processed/labels.parquet` | 19,103 |

These are the same terminal, already-frozen genome-wide DE/screen tables used
by the prior `cross_dataset_genomewide` integration phase (see
`results/tables/cross_dataset_genomewide/dataset_inventory.tsv`, itself
read-only reference here). No upstream QC/alignment/clustering/edgeR/limma/
CRISPR-model step is rerun in this phase.

## Why GSE240112 uses the `tumor_cell` track, not `epithelial`

Consistent with the frozen dataset_inventory: `tumor_cell` is the primary,
pre-declared contrast (author marker-based cancer-cell annotation);
`all_epithelial` is a sensitivity-only secondary track built from the same
libraries and is not treated as a second independent dataset vote anywhere
in this phase either.

## Test-statistic availability

- **GSE118713**: limma export includes `moderated_t` — used directly as the
  ranking statistic (Phase 2).
- **GSE240112, GSE111151, GSE245601 (both tracks)**: the frozen edgeR
  genomewide tables export only `log2fc`, `avg_log_cpm`/`avg_expr`,
  `p_value`, `fdr` — no LR/QLF test statistic column was retained at export
  time. Re-running edgeR to recover that internal statistic would violate
  the "do NOT rerun upstream analyses" instruction for a non-essential gain
  (a fully transparent documented fallback is explicitly permitted instead).
  Ranking for these four uses the documented fallback
  `sign(log2fc) * -log10(p_value)` (Phase 2), never `log2fc` alone.
- **CRISPR**: `labels.parquet` has `effect_size` and `se` (no test statistic
  column exported either) — the Wald statistic `effect_size / se` is used as
  the ranking/enrichment statistic (Phases 2 and 11).

## Gene-symbol harmonization

**Correction (post Codex review):** an earlier draft of this section
claimed `results/tables/cross_dataset_genomewide/gene_mapping_audit.tsv`
was reused for symbol resolution. It is not — `src/systems_network_ranking.py`
never reads that file. The actual implementation is described below.

The per-dataset files above use dataset-native identifiers (Ensembl gene_id
for GSE118713/GSE111151, gene symbol for GSE240112/GSE245601/CRISPR). Each
already carries an upstream-provided `gene_symbol`/`gene_name` column
(populated by that dataset's own earlier, already-frozen phase, e.g.
`gse118713_limma.R`'s Ensembl-to-symbol annotation) -- this phase reads
that column directly rather than re-deriving Ensembl-to-symbol mapping.
Duplicate symbols within a single dataset (e.g. GSE118713's pseudoautosomal
ENSGR/ENSG duplicates, 81 rows; GSE111151, 310 rows) are resolved
independently here, not by reusing the cross-dataset-genomewide phase's
specific resolution choices: `src.systems_network_ranking._dedup_by_max_abs_stat`
keeps, per symbol, the row with the largest `|ranking_stat|` (ties broken by
smallest p-value, then row order) -- see the module's own docstring; each ranked vector is built
directly from its own dataset's genome-wide file (Phase 2).

## GSEA implementation

`fgsea`/`msigdbr`/`STRINGdb`/`reactome.db` are not installed in the project's
R environment (`sc245601` micromamba env) and installing them via BiocManager
in this environment is impractically slow/unreliable for this phase. `gseapy`
(Python, prerank/GSEA implementation of the same weighted Kolmogorov-Smirnov
algorithm underlying `fgsea`/the original Broad GSEA) is used instead and is
documented here explicitly as the implementation choice. No modeling
assumptions differ from a standard preranked GSEA.

## Gene sets

MSigDB Hallmark, Reactome, and GO Biological Process gene sets are not
present locally and are downloaded once via `scripts/download_genesets.sh`
into `data/reference/genesets/` (see `docs/SYSTEMS_NETWORK_GENESETS.md` for
exact source/version/date). All downstream GSEA modules read these local
files only — no network calls at analysis runtime.
