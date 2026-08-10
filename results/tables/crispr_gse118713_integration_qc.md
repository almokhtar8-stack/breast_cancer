# CRISPR-to-GSE118713 master table QC summary

Descriptive, aggregate counts only. FDR values here are evidence for
candidate interpretation, not a predictive feature or a ranking rule.

- CRISPR hits (Gate-1 FDR<0.1): 28
- Uniquely mapped and passed the GSE118713 expression filter: 27
- Uniquely mapped but removed by the GSE118713 expression filter: 1
- Ambiguous mapping (excluded from RNA fields): 0
- Unmatched (no GSE118713 gene ID resolves to this symbol): 0
- Mapped and filter-passing but DE unavailable (RCOR1/KDM1A blinded at source, not recomputed): 1
- TAMR_vs_MCF7 FDR<0.05: 6
- TAMR_vs_FASR FDR<0.05: 15
- Significant in both comparisons: 4
- Directionality where DE is available: 14 up / 12 down in TAMR_vs_MCF7; 17 up / 9 down in TAMR_vs_FASR. Direction is descriptive only and not a causal claim.
