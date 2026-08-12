# Cross-dataset genome-wide integration: data audit

**Date written:** 2026-08-12. Written from direct inspection of the
already-frozen, already-committed genome-wide output files of all five
independent datasets, before any gene-level ranking is computed. No
prior candidate list (13 sensitising, 28 CRISPR hits, PAICS, or any
individually-discussed gene) is used to select or filter which genes
enter this audit.

## The five independent evidence sources

Per the task's Central Principle, each of the five sources below
contributes exactly **one** independent-dataset slot to every ranking
built in this analysis, regardless of how many tracks/contrasts/
sub-analyses exist within it.

### 1. CRISPR functional screen

| Item | Value |
|---|---|
| Source file | `data/processed/labels.parquet` (frozen; sha256-pinned in `results/tables/gate1_decision.tsv`) |
| Rows | 19,103 genes, **all with complete effect_size/p_value/fdr** (0 NA) |
| Gene identifier | Symbol only (unique -- 19,103 rows, 19,103 unique symbols, verified) |
| Effect column | `effect_size` |
| Significance | `p_value`, `fdr` |
| Filtering | None beyond the screen's own guide-fitting model (`n_guides` per gene retained as context) |
| Contrast | Hany-style MCF7-V, E2+4-OHT vs E2 (a genetic-perturbation screen, not an RNA contrast) |
| Direction semantics | **Negative** `effect_size` = knockout relatively depleted under 4-OHT = **sensitising KO** (loss increases tamoxifen sensitivity). **Positive** = knockout relatively favored under 4-OHT = **tolerance-associated KO**. **This sign is a functional-fitness direction and is never treated as equivalent to an RNA log2FC sign** -- a CRISPR-positive gene is not "up" and a CRISPR-negative gene is not "down" in any expression sense. |

### 2. GSE118713 (bulk RNA, acquired tamoxifen/endocrine resistance)

| Item | Value |
|---|---|
| Source file | `results/tables/gse118713_differential_expression_unredacted.tsv.gz` (the **unredacted** table is used -- the KDM1A/RCOR1 blinding referenced in the older `gse118713_differential_expression.tsv.gz` file was retired 2026-08-10 per `PREANALYSIS.md`'s amendments log and `CLAUDE.md`; using the still-redacted file here would silently exclude two now-unblinded genes from a genome-wide search, which is exactly the kind of silent omission this integration must avoid) |
| Contrasts present | `TAMR_vs_MCF7` (primary), `FASR_vs_MCF7`, `TAMR_vs_FASR` (both secondary) -- **only `TAMR_vs_MCF7` is used as this dataset's one independent contribution**, per the task's explicit instruction that TAMR-vs-FASR is secondary context, not a second independent dataset |
| Rows (TAMR_vs_MCF7) | 14,838 |
| Gene identifiers | Ensembl `gene_id` (unique -- 0 duplicates) **and** `gene_symbol` (**71 duplicated symbols** -- see gene-mapping audit) |
| Effect column | `log2fc` |
| Significance | `p_value`, `fdr` |
| Filtering | Genes below an expression threshold were not fit (limma); genes present in this table are exactly the testable set, no separate "present but untestable" distinction is recoverable from this file alone |
| Contrast direction | **Positive** `log2fc` = higher in TAMR (tamoxifen-resistant) vs. MCF7 parental. **Negative** = lower in TAMR. |
| Biological meaning | Bulk RNA-seq of one acquired-resistance cell-line model (established/chronic resistant state, not acute treatment). |

### 3. GSE245601 (acute human single-cell tamoxifen response)

| Item | Value |
|---|---|
| Source files | `results/tables/gse245601_pseudobulk/track_a_genomewide_de.tsv.gz` (all epithelial cells), `track_b_genomewide_de.tsv.gz` (strict malignant cells) -- **both tracks are ONE dataset**, collapsed to a single dataset-level contribution in every ranking (Phase 6/7) |
| Rows | Track A: 17,988; Track B: 13,864 (each track's own `filterByExpr` output -- these ARE the testable sets for that track) |
| Gene identifier | Symbol only (0 duplicates in either track) |
| Effect column | `log2fc` |
| Significance | `p_value`, `fdr` |
| Filtering | `edgeR::filterByExpr`, applied separately per track; a gene absent from a track's file was filtered out (untestable in that track) -- not distinguishable here from "not measured at all" without the pre-filter feature list, which is available (`data/processed/gse245601/...` per-sample matrices, not re-consulted here since testability is what matters for this integration, not raw detection) |
| Contrast direction | **Positive** `log2fc` = higher after 12h ex vivo tamoxifen vs. control media. |
| Biological meaning | **Acute (12h) treatment response, not established/chronic resistance.** Track A/B differ only in which epithelial-cell population is aggregated (all epithelial vs. strict malignant), from the *same* sequencing libraries -- not independent evidence of each other. |

### 4. GSE240112 (human primary-vs-recurrent tumor context)

| Item | Value |
|---|---|
| Source files | `results/tables/gse240112_pseudobulk/tumor_cell_genomewide_de.tsv.gz` (primary analysis), `epithelial_genomewide_de.tsv.gz` (sensitivity analysis, same 6 sequencing libraries reprocessed with a broader cell-selection pipeline -- see `docs/GSE240112_CELLTYPE_DEFINITIONS.md`) |
| Rows | Tumor-cell: 18,429; epithelial: 20,267 |
| Gene identifier | Symbol only (0 duplicates in either track) |
| Effect column | `log2fc` |
| Significance | `p_value`, `fdr` |
| Filtering | `edgeR::filterByExpr` per track |
| Contrast direction | **Positive** `log2fc` = higher in recurrent/tamoxifen-treated (RT) tumor cells vs. primary (PT). |
| Biological meaning | Human tumor tissue, primary-vs-recurrent/treatment-history context -- a *different* biological question from GSE245601's acute ex vivo response and from GSE118713/GSE111151's cell-line acquired-resistance models. **The tumor-cell track is this dataset's one independent contribution to every global ranking; the epithelial track is used only as a same-library robustness/sensitivity flag, never as a second dataset vote**, per the task's explicit instruction. |

### 5. GSE111151 (independent tamoxifen-resistant cell-line panel)

| Item | Value |
|---|---|
| Source file | `results/tables/gse111151/genomewide_de.tsv.gz` |
| Rows | 27,418 |
| Gene identifiers | Ensembl `gene_id` (unique -- 0 duplicates) **and** `gene_name` (**126 duplicated symbols**) |
| Effect column | `log2fc` |
| Significance | `p_value`, `fdr` |
| Filtering | `edgeR::filterByExpr`, cell-line-blocked model (`~ cell_line + resistance_status`) |
| Contrast direction | **Positive** `log2fc` = higher in tamoxifen-resistant sublines vs. their own parental cell line. |
| Biological meaning | A *second*, independent panel of acquired-resistance cell lines (4 different cell-line backgrounds, unlike GSE118713's single MCF7-derived panel) -- genuinely independent of GSE118713 (different cell lines, different laboratory, different study), so both may legitimately each contribute one dataset slot. |

## Dataset inventory table

See `results/tables/cross_dataset_genomewide/dataset_inventory.tsv`.
