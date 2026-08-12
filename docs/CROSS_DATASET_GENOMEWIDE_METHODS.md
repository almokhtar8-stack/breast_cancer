# Cross-dataset genome-wide integration: methods

**Date written:** 2026-08-12. Consolidates the fixed, deterministic rules
implemented across `src/cross_dataset_*.py` (each module's own docstring
is the authoritative detail; this is an index/summary written for
readability). No rule below was adjusted after inspecting which specific
genes it would surface -- every rule is a fixed, general procedure
applied uniformly to the full 37,631-gene union universe.

## Gene universe and harmonization (`src/cross_dataset_gene_mapping.py`)

Union of every gene symbol tested/testable in any of the five datasets'
own frozen output files -- never seeded by, or filtered through, any
prior candidate list (the 13 sensitising candidates, the 28 CRISPR
-significant hits, PAICS, or any individually-discussed gene).
Pseudoautosomal-region `ENSGR`-prefixed duplicate Ensembl IDs with
numerically identical statistics (parsed-value equality, not raw-byte equality) to their `ENSG` counterpart are resolved
(one row kept); any remaining duplicate symbol (a real, distinct second
gene) is excluded from that dataset's single-row table and documented in
`gene_mapping_audit.tsv`, never collapsed by guessing.

## One vote per independent dataset (`src/cross_dataset_evidence_tables.py`, `src/cross_dataset_ranking.py`)

- **CRISPR**: the full 19,103-gene fitted screen, one contribution.
- **GSE118713**: only the `TAMR_vs_MCF7` contrast; `TAMR_vs_FASR`/
  `FASR_vs_MCF7` are secondary context, never a second vote.
- **GSE245601**: Track A (all epithelial) and Track B (strict malignant)
  are the same libraries. The dataset-level evidence percentile is the
  arithmetic mean of both tracks' own within-track percentiles (or the
  single available track's percentile if only one is testable, flagged
  `gse245601_one_track_only`).
- **GSE240112**: the tumor-cell track is the dataset's sole inferential
  contribution; the all-epithelial track (same libraries) is retained
  only as a same-library robustness/direction-agreement flag
  (`gse240112_track_direction_agreement`, `gse240112_outlier_fragility`),
  never blended into the percentile.
- **GSE111151**: the cell-line-blocked `resistant vs parental` model, one
  contribution.

## Within-dataset evidence percentile (`compute_within_dataset_percentile`)

Deterministic rank -> percentile, computed independently within each
dataset's own testable-gene set: sort by FDR ascending, then nominal p
ascending, then |effect| descending; rank 1 (strongest) -> percentile
1.0; rank N (weakest) -> percentile 0.0. Untestable genes get `NaN`,
never `0` (which would falsely imply "tested and weakest"). No percentile
is compared in raw-effect-size units across datasets -- CRISPR effect
size, bulk log2FC, and single-cell pseudobulk log2FC are never treated as
numerically interchangeable.

## Coverage tiers (`assign_coverage_tier`)

Tier A = testable in 5/5 independent datasets, B = 4/5, C = 3/5, D = 2/5,
E = 1/5. The primary global ranking requires >=3/5 (Tier A/B/C);
Tier D/E genes are never discarded -- they remain in the full wide table
and, if they show an extreme (>=99th percentile) signal in their one or
two testable datasets, are surfaced in `high_signal_low_coverage.tsv`
rather than hidden.

## Global ranking (`build_global_ranking`)

A transparent, reconstructable sort hierarchy, not a weighted formula:
(1) coverage tier, (2) number of the 5 datasets with FDR<0.05 (using each
dataset's one representative track), (3) number of datasets in the
top 10% within their own dataset, (4) number in the top 20%, (5) median
evidence percentile across testable datasets, (6) mean evidence
percentile, (7) gene symbol ascending (deterministic tie-break only).
Every summary statistic that could be mistaken for a hidden score is
explicitly named (`equal_dataset_mean_percentile`,
`median_evidence_percentile`) -- no unlabeled composite "final score" is
computed anywhere in this codebase (verified by grep across
`src/cross_dataset_*.py` before the Codex review, Phase 25).

## Alternative views (`src/cross_dataset_consensus_views.py`)

**Resistance-state consensus** uses only GSE118713 + GSE240112 +
GSE111151 (established/chronic acquired-resistance or recurrence
context); GSE245601 (acute 12h ex vivo) is explicitly excluded; CRISPR
significance is never required. **Functional CRISPR ranking** covers the
full fitted screen (not just prior FDR<0.10 hits), sensitising and
tolerance-associated knockouts kept as separate rankings (opposite
biological meanings, never merged onto one scale). **Human-only**
(GSE245601 + GSE240112) and **RNA-only** (all four RNA datasets, CRISPR
excluded) use the same equal-percentile-mean logic as the global ranking,
restricted to their respective dataset subsets. **CRISPR-independent
discovery** surfaces genes with FDR>=0.10 (or untestable) in CRISPR but
repeated resistance-RNA support, explicitly labeled
`resistance_biomarker_or_pathway_candidate`, never dismissed.

## Evidence categories (`assign_evidence_category`)

One deterministic, precedence-ordered decision tree (LOW_COVERAGE ->
MULTIMODAL_STRONG -> RNA_RESISTANCE_CONSENSUS -> FUNCTIONAL_ONLY ->
ACUTE_RESPONSE -> HUMAN_TUMOR_SUPPORTED -> CONTEXT_DEPENDENT ->
LOW_EVIDENCE), evaluated top-to-bottom, first match wins, every gene
gets exactly one category. CRISPR significance threshold (FDR<0.10)
matches this project's pre-existing Gate-1 threshold
(`config.yaml gate1.fdr_threshold`), not a newly invented cutoff.

## Stability (`src/cross_dataset_stability.py`)

Four alternative equal-treatment sort schemes (median-percentile-first,
mean-percentile-first, FDR-count-first, top-10%-count-first) plus five
leave-one-dataset-out reruns of the main ranking hierarchy on the
remaining four datasets. `ROBUST` = stays in the global Top 20 under all
5 leave-one-out variants; `MODERATELY_STABLE` = under >=3/5;
`DATASET_DEPENDENT` = otherwise.

## Anonymization audit (`src/cross_dataset_anonymization_audit.py`)

Every gene reassigned a `GeneNNNNN` ID under a fixed-seed but
**non-alphabetical** shuffle (chosen specifically so the tie-break sort
key differs from the original run -- an alphabetical anonymization would
trivially reproduce the original order regardless of any hidden bias).
The full percentile + ranking pipeline is rerun from scratch on the
anonymized table. Result on the real, corrected data (after the Phase 27
Codex review's fixes -- see `docs/CROSS_DATASET_GENOMEWIDE_CODEX_REVIEW.md`):
**15,255/15,255 genes' global ranks matched exactly** between the named
and anonymized runs, zero mismatches.

## What this analysis explicitly does not do

No arbitrary weighted score (Phase 25); no gene-name-based branching
anywhere in the ranking logic (verified empirically by the anonymization
audit); no CRISPR sign compared directly to any RNA sign (CRISPR
direction uses `sensitising_KO`/`tolerance_associated_KO` vocabulary,
RNA directions use `up`/`down` vocabulary, kept in separate columns and
never merged); GSE245601's acute 12h response never counted as resistance
-state evidence; no dataset's multiple tracks/contrasts ever cast more
than one vote in any ranking.
