# Candidate adjudication: MULTIMODAL_STRONG category audit

**Date written:** 2026-08-12. This is a verification exercise, not a
re-run: the rule below was re-derived independently in
`src/candidate_adjudication_category_audit.py` (written from the rule's
own docstring, not by importing `assign_evidence_category`), then checked
against the frozen `results/tables/cross_dataset_genomewide/evidence_categories.tsv`
output. Source: `all_genes_cross_dataset_evidence_with_ranking.tsv`,
`resistance_consensus_all_genes.tsv` (both frozen, read-only).

## 1. Exact rule

A gene is `MULTIMODAL_STRONG` iff, evaluated on the 15,255-gene primary
ranking table (coverage tier A/B/C only -- the category tree's first
branch, `LOW_COVERAGE`, removes tier D/E genes before this rule is ever
applied):

```
crispr_fdr < 0.10
AND (
      resistance_fdr05_count >= 2
      OR
      (resistance_fdr05_count >= 1  AND  n_datasets_fdr05_overall >= 2)
    )
```

where `resistance_fdr05_count` is the number of the three resistance
-state datasets (GSE118713, GSE240112 tumor-cell track, GSE111151) with
FDR<0.05, and `n_datasets_fdr05_overall` is the number of *all five*
datasets (CRISPR + the four RNA datasets, GSE245601 collapsed to its
representative track) with FDR<0.05.

## 2. Did all seven satisfy the same rule?

**No -- two distinct sub-paths, not one uniform pattern.** Per-gene FDRs:

| gene | crispr_fdr | gse118713_fdr | gse240112_tumor_fdr | gse111151_fdr | gse245601_epi_fdr |
|---|---|---|---|---|---|
| CUX1 | 0.0450 | **0.0008** | **0.0326** | 0.353 | 0.802 |
| SOX2 | 0.0122 | **0.0093** | **0.0001** | 0.632 | 0.694 |
| DPP9 | 0.0417 | **0.0021** | 0.802 | 0.954 | 0.971 |
| LZTR1 | 0.0124 | **0.0151** | 0.608 | 0.717 | 0.913 |
| TFAP2C | 0.0482 | **0.0010** | 0.461 | 0.376 | 0.162 |
| USP34 | 0.0417 | **0.0073** | 0.228 | 0.632 | 0.901 |
| VEZF1 | 0.0373 | 0.2375 | **0.0195** | 0.608 | 0.890 |

- **CUX1 and SOX2** qualify through the stronger path: **two independent
  resistance datasets** each at FDR<0.05 (`resistance_fdr05_count >= 2`).
- **DPP9, LZTR1, TFAP2C, USP34, VEZF1** qualify only through the weaker
  path: exactly **one** resistance dataset at FDR<0.05, with CRISPR's own
  FDR<0.05 (all five are <0.05, not just <0.10) supplying the second of
  the two datasets required by `n_datasets_fdr05_overall >= 2`. Without
  CRISPR being independently significant at the stricter 0.05 level (not
  merely the 0.10 category threshold), these five would not qualify.
- **GSE111151 and GSE245601 contribute zero FDR<0.05 hits to any of the
  seven genes.** All resistance support among the seven comes from
  GSE118713 (6/7 genes) and/or GSE240112 tumor-cell (3/7 genes:
  CUX1, SOX2, VEZF1).

## 3. Are any borderline?

Yes. DPP9, LZTR1, TFAP2C, USP34, and VEZF1 rest on a single resistance
-dataset FDR plus CRISPR's own significance -- remove that one resistance
dataset's significance and they would drop to `FUNCTIONAL_ONLY` (CRISPR
strong, no resistance-dataset FDR<0.05). CUX1 and SOX2 are more robust:
they would still satisfy `resistance_fdr05_count >= 1` even after losing
either one of their two significant resistance datasets, as long as
CRISPR remains significant. TFAP2C's CRISPR FDR (0.0482) is the closest
of the seven to the 0.05 boundary used implicitly by the weaker path.

## 4. How many genes nearly satisfied the category?

Applying the same rule, exactly **two** tier-A genes are one condition
away from `MULTIMODAL_STRONG`:

| gene | crispr_fdr | n_datasets_fdr05 | resistance_fdr05_count | missing condition |
|---|---|---|---|---|
| SPRED2 | 0.0783 | 1 | 1 | has CRISPR FDR<0.10, one resistance dataset FDR<0.05, but only 1 dataset overall at FDR<0.05 (CRISPR itself falls in [0.05,0.10) so does not count toward the FDR<0.05 tally) |
| EML5 | 0.1487 | 2 | 2 | has the full resistance-FDR pattern (2 datasets FDR<0.05) but CRISPR FDR (0.149) misses the 0.10 gate |

This narrow ("exactly one condition short") definition is intentionally
strict. A broader near-miss search (allowing comparable percentile
strength rather than a strict FDR-boundary miss) is done separately in
Phase 4 (`multimodal_near_misses.tsv`) and returns more candidates.

## 5. Genes just outside deserving inspection?

SPRED2 and EML5 both deserve a look in Phase 4/9: SPRED2 has real CRISPR
signal (nominally close to the 0.10 gate) plus one significant resistance
dataset -- essentially the same evidence pattern as DPP9/LZTR1/TFAP2C/
USP34/VEZF1, just missing the "2 datasets FDR<0.05 overall" technicality
because its own CRISPR FDR sits in [0.05, 0.10) rather than below 0.05.
EML5 is the mirror case: full resistance support, CRISPR just above the
0.10 gate. Neither is manufactured into the category -- they are reported
as near-misses, not promoted.

## Conclusion

The independently-reconstructed rule reproduces the frozen
`evidence_categories.tsv` MULTIMODAL_STRONG set **exactly**: 7 genes, no
more, no fewer (`src/candidate_adjudication_category_audit.py` raises on
any mismatch and none occurred). The category is real but not uniform --
it contains two evidentially stronger genes (CUX1, SOX2, doubly
resistance-supported) and five that rest on a single resistance dataset's
significance plus CRISPR's own strength. This distinction carries forward
into the Phase 9 head-to-head.
