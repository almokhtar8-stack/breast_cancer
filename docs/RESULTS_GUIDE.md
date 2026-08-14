# Results Guide

What each major result table/figure/report means, and — critically —
whether it is **CURRENT** (the live, reported value) or **ARCHIVED**
(kept only for traceability, never the reported value). Read this before
citing any specific number from this repository.

For "where do I start," see [`../results/final/README.md`](../results/final/README.md)
instead — this document is a reference lookup, not a starting point.

---

## The single most important current-vs-archived distinction: DepMap release

**DepMap 26Q1 is the CURRENT, ACTIVE, REPORTED release. DepMap 24Q4 is
ARCHIVED and is NEVER the reported value for anything in this repository.**

- `results/tables/independent_validation/DepMap_candidate_dependency.tsv`
  — **CURRENT (26Q1)**. Every row's `depmap_release` column reads `26Q1`.
  This is the table cited everywhere else in the project.
- `results/tables/independent_validation/archive_24Q4/DepMap_candidate_dependency_24Q4.tsv`
  — **ARCHIVED**. Kept only so the 26Q1 values can be shown to reproduce
  the original 24Q4 classification tiers exactly (they do, for all four
  frozen candidates). **This archived file has no `depmap_release` column
  at all** — its path (`archive_24Q4/`) and its own `_24Q4` filename suffix
  are the only markers of its vintage, so never rely on column content
  alone to judge a DepMap table's currency; check the path first. If a copy
  of this file is ever extracted or renamed without its path, treat any
  unlabeled DepMap dependency table you can't trace back to a path as
  unknown-vintage, not as current.
- Same rule applies everywhere DepMap expression/dependency/codependency is
  used: `DepMap_candidate_expression.tsv`, `DepMap_candidate_codependency.tsv`,
  the GDSC pharmacogenomics phase's cell-line join (`final_pharmacogenomics`,
  26Q1 only) — all 26Q1, never 24Q4, unless the path explicitly says
  `archive_24Q4`.

If a number you're reading is from a **current** (non-archived) DepMap
table, it will carry a `depmap_release` column reading `26Q1` — check it
before citing. If the table has no such column, check its path for
`archive_24Q4` before assuming it's current.

---

## Frozen vs. non-frozen tables

**Frozen** (never altered after the phase that produced them; later phases
only ever *read* these):

- `results/tables/evidence_freeze/THERAPEUTIC_SHORTLIST_FREEZE.tsv` — the
  4-gene shortlist ranking (USP34 > VEZF1 > EML5 > CITED2).
- `results/tables/gate1_decision.tsv` — the CRISPR Gate-1 hit calls.
- `results/tables/cross_dataset_genomewide/` — the genome-wide integration.
- `results/tables/systems_network/` — the network/pathway mapping.
- `results/networks/systems_network/` — the Cytoscape network files.

**Non-frozen but stable** (extended, never contradicted, by later phases):

- `results/tables/independent_validation/` — TCGA/DepMap validation adds
  context but never changes the frozen ranking above.
- `results/tables/lead_target_deep_dive/`, `druggability_safety/` —
  structural/safety characterization of USP34/VEZF1 only.

**Most recent / actively updated this run:**

- `results/tables/final_translational/` — updated to incorporate PMID
  28499884 EMT/stemness counter-evidence and the split EXP-2A/EXP-2B
  normal-cell comparators.
- `results/tables/final_pharmacogenomics/` — new this run; the GDSC
  Release 8.5 drug-response lookup for USP34/VEZF1.

## Per-phase canonical report (read this, not the raw tables, for narrative)

| Phase | Canonical report |
|---|---|
| Evidence freeze | [`THERAPEUTIC_SHORTLIST_FREEZE.md`](THERAPEUTIC_SHORTLIST_FREEZE.md) |
| Cross-dataset genome-wide | [`CROSS_DATASET_GENOMEWIDE_ANALYSIS_REPORT.md`](CROSS_DATASET_GENOMEWIDE_ANALYSIS_REPORT.md) |
| Systems/network | [`../results/reports/systems_network/four_candidate_network_audit.md`](../results/reports/systems_network/four_candidate_network_audit.md) |
| Literature mechanism | [`../results/reports/literature_mechanism/four_candidate_mechanism_review.md`](../results/reports/literature_mechanism/four_candidate_mechanism_review.md) |
| Independent validation (TCGA/DepMap) | [`../results/reports/independent_validation/four_candidate_TCGA_DepMap_review.md`](../results/reports/independent_validation/four_candidate_TCGA_DepMap_review.md) |
| Lead-target deep dive | [`../results/reports/lead_target_deep_dive/USP34_VEZF1_translational_deep_dive.md`](../results/reports/lead_target_deep_dive/USP34_VEZF1_translational_deep_dive.md) |
| Druggability/safety | [`../results/reports/druggability_safety/four_candidate_druggability_safety_review.md`](../results/reports/druggability_safety/four_candidate_druggability_safety_review.md) |
| Final translational plan | [`../results/reports/final_translational/final_USP34_VEZF1_translational_plan.md`](../results/reports/final_translational/final_USP34_VEZF1_translational_plan.md) |
| GDSC pharmacogenomics | [`../results/reports/final_pharmacogenomics/USP34_VEZF1_GDSC_review.md`](../results/reports/final_pharmacogenomics/USP34_VEZF1_GDSC_review.md) |

## Other current-vs-superseded distinctions to watch for

- **GSE245601 malignant-cell classification:** InferCNV is the **primary**
  call; CopyKAT is an **independent sensitivity check** on the same cells,
  never a second vote to average against InferCNV. Concordance is reported
  transparently (~56% average agreement, highly variable per sample) — see
  [`CNV_METHOD_AUDIT.md`](CNV_METHOD_AUDIT.md).
- **GDSC1 vs GDSC2:** two separate screening campaigns, never pooled for
  FDR correction. A compound can carry more than one `DRUG_ID` even within
  one release/campaign (e.g. AZD7762 has two distinct GDSC1 `DRUG_ID`s) —
  always check `DRUG_ID`, not just `DRUG_NAME`, before treating two rows as
  independent.
- **Anonymized ranking table:** `results/tables/cross_dataset_genomewide/anonymized_ranking.tsv`
  and `anonymized_gene_mapping.tsv` exist for a specific ranking-stability
  check and are **not** the primary gene-identified ranking table — use
  `all_genes_cross_dataset_evidence_with_ranking.tsv` for actual gene names.
- **Follow-up rankings vs. the frozen ranking:** `four_candidate_followup_rankings.tsv`
  is explicitly labeled, in its own `note` column, as a follow-up-order
  signal only — it does not and must not be read as a revision of the
  frozen `THERAPEUTIC_SHORTLIST_FREEZE.tsv` ranking. See
  [`../README.md`](../README.md#f-the-frozen-candidates--three-distinct-rankings-not-one).
