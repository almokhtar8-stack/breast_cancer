# Cross-dataset genome-wide integration: final report

**Date written:** 2026-08-12. All numbers below are read directly from
the committed output tables (post-Codex-review corrections), not
hand-typed. See `docs/CROSS_DATASET_GENOMEWIDE_DATA_AUDIT.md`,
`docs/CROSS_DATASET_GENOMEWIDE_METHODS.md`, and
`docs/CROSS_DATASET_GENOMEWIDE_CODEX_REVIEW.md` for full methodology and
the review that found and fixed 3 material bugs before this report was
finalized.

## 1. Full universe

37,631 unique harmonized genes (union across all 5 datasets, never
seeded by any prior candidate list). Coverage breakdown:

| Datasets testable | Genes | Tier |
|---|---|---|
| 5/5 | 11,093 | A |
| 4/5 | 1,934 | B |
| 3/5 | 2,228 | C |
| 2/5 | 5,725 | D |
| 1/5 | 15,866 | E |
| 0/5 | 785 | (present in the union via one dataset's broader "presence" concept but with no usable inferential percentile from any dataset) |

15,255 genes (Tier A/B/C) are eligible for the primary global ranking.

## 2. Global Top 20

(`results/tables/cross_dataset_genomewide/top20_global.tsv`)

| Rank | Gene | Coverage | Datasets FDR<0.05 | Datasets top-10% | Resistance consensus | Evidence category |
|---|---|---|---|---|---|---|
| 1 | GREB1 | A | 3 | 4 | all_down | RNA_RESISTANCE_CONSENSUS |
| 2 | LARP6 | A | 3 | 4 | all_up | RNA_RESISTANCE_CONSENSUS |
| 3 | PLOD2 | A | 3 | 4 | all_up | RNA_RESISTANCE_CONSENSUS |
| 4 | MED13L | A | 3 | 4 | all_down | RNA_RESISTANCE_CONSENSUS |
| 5 | FGD3 | A | 3 | 3 | all_down | RNA_RESISTANCE_CONSENSUS |
| 6 | EPAS1 | A | 3 | 3 | all_up | RNA_RESISTANCE_CONSENSUS |
| 7 | MYBL1 | A | 3 | 3 | all_down | RNA_RESISTANCE_CONSENSUS |
| 8 | SYTL5 | A | 3 | 3 | majority_down | RNA_RESISTANCE_CONSENSUS |
| 9 | HRK | A | 3 | 3 | all_up | RNA_RESISTANCE_CONSENSUS |
| 10 | CUX1 | A | 3 | 3 | all_down | MULTIMODAL_STRONG |
| 11 | PTGER4 | A | 3 | 2 | all_up | RNA_RESISTANCE_CONSENSUS |
| 12 | C6orf141 | A | 3 | 2 | all_down | RNA_RESISTANCE_CONSENSUS |
| 13 | LAMB3 | A | 3 | 2 | all_up | RNA_RESISTANCE_CONSENSUS |
| 14 | RUNX1 | A | 3 | 2 | majority_up | RNA_RESISTANCE_CONSENSUS |
| 15 | DMRTA1 | A | 3 | 2 | all_up | RNA_RESISTANCE_CONSENSUS |
| 16 | SOX2 | A | 3 | 2 | all_up | MULTIMODAL_STRONG |
| 17 | MSMO1 | A | 3 | 2 | all_up | RNA_RESISTANCE_CONSENSUS |
| 18 | NEDD9 | A | 3 | 2 | majority_up | RNA_RESISTANCE_CONSENSUS |
| 19 | WBP2 | A | 3 | 2 | majority_down | RNA_RESISTANCE_CONSENSUS |
| 20 | HLA-DQB1 | A | 3 | 1 | all_up | RNA_RESISTANCE_CONSENSUS |

Main caveats: none flagged for any Top-20 gene (all Tier A, no
single-track-only or mixed-direction issues).

## 3. Resistance-state Top 20

(GSE118713 + GSE240112 tumor-cell + GSE111151 only; no CRISPR
requirement; GSE245601 excluded)

EPAS1, GABBR2, SLC4A10, DMRTA1, LAMB3, HLA-DQB1 (all 3/3 resistance
datasets FDR<0.05, clean direction consensus), then GREB1, SPRY4, GFRA1,
SUSD3, DLX2, SULF1, DSCAM, EYA2, GJB2, ACOT4, LARP6, MALRD1, PDGFC, IL20
(2/3 FDR<0.05, all_up or all_down). Notably, several of these (GREB1,
SUSD3, SULF1, GFRA1) are established ER-pathway genes in the breast
cancer literature.

## 4. RNA-only Top 20

(All 4 RNA datasets, CRISPR excluded entirely)

GREB1, LARP6, EPAS1, MYBL1, SLC4A10, GABBR2, SYTL5, PLOD2, FGD3, DMRTA1,
PTGER4, MED13L, HRK, C6orf141, LAMB3, RUNX1, HLA-DQB1, RRAGD, NEDD9,
MSMO1 -- nearly identical to the global Top 20 (RNA evidence dominates
the ranking since 4 of 5 datasets are RNA).

## 5. Functional CRISPR Top 20 (sensitising direction)

Ranks 1-13: KDM1A, TLK2, TADA2B, USP17L29, VEZF1, HMGB1, EIF4ENIF1,
SUPT4H1, ICK, USP34, PET117, TSR3, CTDNEP1 -- **exactly** the frozen 13
sensitising candidates from earlier phases of this project, in the same
relative order, confirming the functional-CRISPR-only ranking reproduces
the prior hand-curated list exactly (a consistency check, not a new
finding). Ranks 14-20 (DOT1L, MRGBP, F8A3, DCAF7, CITED2, TAF6L, SYNCRIP)
extend beyond that frozen set -- new sensitising-direction candidates
that did not previously meet the Gate-1 FDR<0.1 threshold used to define
the "28 CRISPR-significant hits" but rank immediately below them here.

## 6. CRISPR-nonsignificant RNA Top 20

**This is the key list the old CRISPR-first strategy would have missed
entirely.** Top entries (crispr_fdr, resistance_fdr05_count,
resistance_direction_consensus): EPAS1 (0.92, 3, all_up), SLC4A10 (0.97,
3, all_down), GABBR2 (0.93, 3, all_up), DMRTA1 (0.98, 3, all_up), LAMB3
(0.77, 3, all_up), HLA-DQB1 (0.92, 3, all_up), NCF2 (0.81, 2,
majority_up), GREB1 (0.82, 2, all_down), SPRY4 (0.79, 2, all_up), GFRA1
(0.79, 2, all_down), MALL (0.79, 2, majority_up), **SUSD3** (0.79, 2,
all_down), DLX2 (0.89, 2, all_up), HMGA2 (0.89, 2, majority_up), CXCL8
(0.95, 2, all_up), LOX (0.72, 2, all_up), HSPB8 (0.94, 2, all_down),
**SULF1** (0.80, 2, all_down), SSTR5-AS1 (untested, 2, all_up), DSCAM
(0.94, 2, all_down). SUSD3 and SULF1 in particular are well-known ER
-pathway/endocrine-context genes in the literature despite having no
CRISPR screen signal whatsoever.

## 7. Human-tumor Top 20

(GSE245601 + GSE240112 only) GREB1, SYTL5, MYBL1, LARP6, FGD3,
AL121578.3, PTGER4, HRK, MSMO1, AC105411.1, AL160408.1, RRAGD, MED13L,
RUNX1, TRPS1, NEU1, C6orf141, PALLD, LPIN1, IRS1.

## 8. Robust top genes

Only **PLOD2** is `ROBUST` (stays in the global Top 20 under all 5
leave-one-dataset-out variants; best_rank=3, worst_rank=42). 14 of the
remaining Top 20 are `MODERATELY_STABLE` (stay in Top 20 under >=3/5
leave-one-out variants); EPAS1, HRK, DMRTA1, NEDD9, HLA-DQB1 are
`DATASET_DEPENDENT`. Removing GSE245601 causes the largest disruption to
the Top 20 overall (`results/figures/cross_dataset_genomewide/final_review/08_leave_one_dataset_out.png`),
reflecting that dataset's larger testable-gene count and correspondingly
finer percentile resolution -- worth further methodological attention in
any follow-up, but does not by itself indicate a coding error (verified
in the Codex review).

## 9. Multimodal convergent genes

7 genes are `MULTIMODAL_STRONG` (strong CRISPR signal, FDR<0.10, AND
independent resistance-RNA support in >=2 datasets): **USP34**, **VEZF1**
(both previously studied in this project), plus **CUX1, DPP9, LZTR1,
SOX2, TFAP2C** -- all newly surfaced by this unbiased search. SOX2
(stemness) and TFAP2C (AP-2γ, a luminal-breast-cancer master
transcription factor) are particularly biologically plausible additions.

## 10. Functional-only genes

19 genes: strong CRISPR signal (FDR<0.10) but no resistance-RNA or acute
-RNA support at FDR<0.05 in any dataset (`results/tables/cross_dataset_genomewide/crispr_functional_all_genes.tsv`
filtered to this category) -- these remain candidates whose importance,
if real, is not reflected in RNA abundance and would need orthogonal
(e.g. protein-level or pathway) follow-up to interpret.

## 11. RNA-consensus / CRISPR-weak genes

5,726 genes are `RNA_RESISTANCE_CONSENSUS` (repeated resistance-RNA
support, CRISPR nonsignificant) -- by far the largest category, including
17 of the 20 global Top-20 genes. This is the single biggest finding of
this analysis: **the vast majority of the strongest cross-dataset RNA
signal in this project has no CRISPR screen counterpart**, meaning a
CRISPR-first strategy structurally cannot find these genes.

## 12. Context-dependent / discordant genes

9 genes are `CONTEXT_DEPENDENT` (substantial evidence, >=2 datasets at
FDR<0.05 or top-10%, but a "mixed" resistance-direction consensus):
C2orf54, CHRNG, CTSZ, FOXD3, LPPR2, PAH, PIPOX, SEPT6, TFF2.

## 13. Surprises

**18 of the global Top 20 (90%) were never part of the project's prior
28 CRISPR-significant hits.** Only SOX2 and CUX1 overlap with the old
list. GREB1 (#1), LARP6 (#2), PLOD2 (#3, the only ROBUST gene), MYBL1,
EPAS1, RUNX1 and 12 others are entirely new to this project's evidence
base, discovered only because the search stopped requiring CRISPR
significance as a gate.

## 14. Previously discussed genes

| Gene | Global rank (of 15,255) | Evidence category | Stability |
|---|---|---|---|
| USP34 | 994 | MULTIMODAL_STRONG | DATASET_DEPENDENT |
| VEZF1 | 1,138 | MULTIMODAL_STRONG | DATASET_DEPENDENT |
| SUPT4H1 | 2,518 | FUNCTIONAL_ONLY | DATASET_DEPENDENT |
| ICK | 2,584 | FUNCTIONAL_ONLY | DATASET_DEPENDENT |
| PAICS | 7,634 | LOW_EVIDENCE | DATASET_DEPENDENT |

USP34 and VEZF1 both hold up respectably (top ~7% of eligible genes,
genuine `MULTIMODAL_STRONG` classification) but neither reaches the
global Top 20 -- their prior prominence in this project reflected being
part of the pre-specified 13-candidate CRISPR-sensitising family, not
that they were the single strongest cross-dataset signals overall.

## 15. Did the old 13/28 filter miss anything important?

**Yes, unambiguously.** 18 of the global Top 20, and the entire
5,726-gene `RNA_RESISTANCE_CONSENSUS` category (including well-known
ER-pathway genes GREB1, SUSD3, SULF1, GFRA1), would never have been
examined under the old CRISPR-first strategy. The 28-gene CRISPR filter
correctly identified genes with strong functional screen evidence, but
that is only one of several ways a gene can matter for tamoxifen
resistance, and it structurally excludes any gene whose functional role
is not captured by this specific CRISPR screen's design (e.g. genes
whose knockout effect is buffered, genes not covered by the guide
library, or genes whose relevance is transcriptional/biomarker-like
rather than functional-fitness-like).

## 16. Best overall candidates (no single winner forced)

- **Strongest multimodal candidates:** USP34, VEZF1, SOX2, TFAP2C, CUX1,
  LZTR1, DPP9 (all `MULTIMODAL_STRONG`).
- **Strongest resistance-expression candidates:** GREB1, LARP6, PLOD2,
  EPAS1, MYBL1 (global Top 5, all with resistance-RNA consensus).
- **Strongest functional-sensitisation candidates:** KDM1A, TLK2, TADA2B
  (largest-magnitude, most significant sensitising CRISPR effects).
- **Strongest human-tumor candidates:** GREB1, SYTL5, MYBL1, LARP6, FGD3
  (top of the human-only view).

## 17. Figures

`results/figures/cross_dataset_genomewide/final_review/`: 12 figures
(`01_all_gene_evidence_heatmap_top100.png` through
`12_surprise_candidates.png`).

## 18. Codex verdict

Initial verdict: **FAIL** (3 material bugs: a phantom GSE240112 coverage
vote for epithelial-only genes, a leave-one-out hierarchy that didn't
match the main ranking, and a missing branch in the CONTEXT_DEPENDENT
category). All 3 were fixed, the full pipeline was rerun, and a follow
-up Codex pass independently re-verified each fix and gave a final
verdict of **PASS WITH NOTES**. Full detail:
`docs/CROSS_DATASET_GENOMEWIDE_CODEX_REVIEW.md`.

## 19. Tests

109 cross-dataset-specific tests passed
(`python3 -m pytest tests/ -k cross_dataset -q`); 681 passed for the full
project test suite, confirming no regression to any existing frozen
module.

## 20. Git commit

See final commit hash reported at the end of this run.
