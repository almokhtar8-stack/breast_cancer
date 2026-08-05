# Tamoxifen response modifiers in ER-positive breast cancer

Which gene knockouts change how ER-positive breast cancer cells respond to
tamoxifen, and which of those targets are least likely to add to the joint
pain, muscle pain and fatigue that cause patients to stop endocrine therapy.

Computational only. Public data. No wet lab.

## Data

| Source | Role |
|---|---|
| Hany 2023, Sci Adv 9:eadd3685 | Labels: MAGeCK-VISPR median-ratio-normalised sgRNA abundance; E2+4-OHT versus E2 |
| GSE118713 | Transcript-level TPM aggregated to gene level; MCF7 / TAMR / FASR, n=3 per group; primary replicated resistance dataset for feature construction |
| GSE111151 | Raw counts plus published uncorrected log2 CPM; 4 parental lines and 7 TamR derivatives, n=1 per condition; independent post-model confirmation only and excluded from model training |

## Claim

The screen was run in a drug-tolerant parental clone, not an acquired
resistant derivative. The claim is *modulates tamoxifen response*, not
restores sensitivity.

## Version-controlled outputs

`data/processed/labels.parquet` is force-added and tracked in git, unlike
the rest of `data/processed/*` (gitignored). It is the frozen per-gene
label table that every downstream analysis reads from, so its history
needs to be reviewable the same way code is.

`data/processed/gse118713_gene_tpm.parquet` and
`data/processed/gse111151_log2cpm.parquet` are likewise force-added. They
are the frozen, checksum-verified gene-level expression matrices (see
`results/tables/gse118713_preparation_qc.tsv` and
`results/tables/gse111151_preparation_qc.tsv` for provenance) that every
downstream Phase 2 analysis reads from, so they need the same reviewable
history as the labels.

`data/processed/gse118713_gene_tpm_filtered.tsv.gz` is also force-added.
It is frozen because it is the exact statistical input used for the
committed Phase 2B limma analysis (see
`results/tables/gse118713_expression_filtering.tsv` for the filtering
record and gene counts that produced it).
