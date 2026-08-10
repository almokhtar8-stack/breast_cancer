# CRISPR-to-GSE118713 master table QC summary

Descriptive, aggregate counts only. FDR values here are evidence for
candidate interpretation, not a predictive feature or a ranking rule.

- CRISPR hits (Gate-1 FDR<0.1): 28
- Uniquely mapped and passed the GSE118713 expression filter: 27
- Uniquely mapped but removed by the GSE118713 expression filter: 1
- Ambiguous mapping (excluded from RNA fields): 0
- Unmatched (no GSE118713 gene ID resolves to this symbol): 0
- Mapped and filter-passing but DE unavailable (historical blind gap, expected to be zero when reading the post-unblinding specificity table): 0
- DE values available: 27
- TAMR_vs_MCF7 FDR<0.05: 6
- TAMR_vs_FASR FDR<0.05: 16
- Significant in both comparisons: 4
- Directionality where DE is available: 15 up / 12 down in TAMR_vs_MCF7; 18 up / 9 down in TAMR_vs_FASR. Direction is descriptive only and not a causal claim.
