# Final Results — Canonical Index

This page is the canonical entry point to this project's final findings. It
links to (never duplicates) the actual tables, figures, and reports that
carry the real numbers — see [`../../docs/RESULTS_GUIDE.md`](../../docs/RESULTS_GUIDE.md)
for what's current vs. archived across the whole project.

## Frozen four-gene shortlist

**USP34 > VEZF1 > EML5 > CITED2**, frozen at the evidence-freeze phase and
never altered afterward.

- Table: [`../tables/evidence_freeze/THERAPEUTIC_SHORTLIST_FREEZE.tsv`](../tables/evidence_freeze/THERAPEUTIC_SHORTLIST_FREEZE.tsv)
- Report: [`../../docs/THERAPEUTIC_SHORTLIST_FREEZE.md`](../../docs/THERAPEUTIC_SHORTLIST_FREEZE.md)

## Final translational conclusions: USP34 = lead, VEZF1 = second/backup

- Table: [`../tables/final_translational/final_translational_conclusions.tsv`](../tables/final_translational/final_translational_conclusions.tsv)
- Full report: [`../reports/final_translational/final_USP34_VEZF1_translational_plan.md`](../reports/final_translational/final_USP34_VEZF1_translational_plan.md)
- Proposed (not yet run) wet-lab plan: EXP-1 (USP34), EXP-2A/EXP-2B (USP34
  normal-cell comparators), EXP-3 (VEZF1), EXP-4 (VEZF1 comparator), EXP-5
  (VEZF1-TEAD1 hypothesis test) — see the report above, Section M of the
  [main README](../../README.md#m-proposed-next-experiment) for a summary.
- Structural/pocket analysis conclusion: **`DOCKING_NOT_YET_JUSTIFIED`** —
  table: [`../tables/final_translational/USP34_docking_decision.tsv`](../tables/final_translational/USP34_docking_decision.tsv).

## GDSC Release 8.5 pharmacogenomics (USP34/VEZF1 only, correlational)

- Full report: [`../reports/final_pharmacogenomics/USP34_VEZF1_GDSC_review.md`](../reports/final_pharmacogenomics/USP34_VEZF1_GDSC_review.md)
- USP34: 9 of 1,278 drug-tests reach FDR<0.05 (strongest: AZD7762/CHK1-CHK2,
  FDR=0.008) — table: [`../tables/final_pharmacogenomics/USP34_GDSC_drug_associations.tsv`](../tables/final_pharmacogenomics/USP34_GDSC_drug_associations.tsv).
- VEZF1: 0 FDR-significant hits — genuine negative result — table:
  [`../tables/final_pharmacogenomics/VEZF1_GDSC_drug_associations.tsv`](../tables/final_pharmacogenomics/VEZF1_GDSC_drug_associations.tsv).
- Final interpretation (does GDSC change USP34=LEAD or VEZF1=BACKUP? **No,
  for both**): [`../tables/final_pharmacogenomics/GDSC_final_interpretation.tsv`](../tables/final_pharmacogenomics/GDSC_final_interpretation.tsv).
- Figures: [`../figures/final_pharmacogenomics/`](../figures/final_pharmacogenomics/).

## Independent validation (TCGA-BRCA + DepMap 26Q1)

- Report: [`../reports/independent_validation/four_candidate_TCGA_DepMap_review.md`](../reports/independent_validation/four_candidate_TCGA_DepMap_review.md)
- Table: [`../tables/independent_validation/four_candidate_independent_validation.tsv`](../tables/independent_validation/four_candidate_independent_validation.tsv)

## Structural / druggability / safety characterization

- Lead-target deep dive: [`../reports/lead_target_deep_dive/USP34_VEZF1_translational_deep_dive.md`](../reports/lead_target_deep_dive/USP34_VEZF1_translational_deep_dive.md)
- Druggability/safety review: [`../reports/druggability_safety/four_candidate_druggability_safety_review.md`](../reports/druggability_safety/four_candidate_druggability_safety_review.md)

## Mechanism and pathway context

- Literature mechanism review: [`../reports/literature_mechanism/four_candidate_mechanism_review.md`](../reports/literature_mechanism/four_candidate_mechanism_review.md)
- Systems/network mapping: [`../reports/systems_network/four_candidate_network_audit.md`](../reports/systems_network/four_candidate_network_audit.md)

## Foundational evidence (frozen, upstream of the shortlist)

- Genome-wide cross-dataset integration: [`../../docs/CROSS_DATASET_GENOMEWIDE_ANALYSIS_REPORT.md`](../../docs/CROSS_DATASET_GENOMEWIDE_ANALYSIS_REPORT.md)
- CRISPR Gate-1 decision: [`../tables/gate1_decision.tsv`](../tables/gate1_decision.tsv)

---

For the full 16-phase workflow narrative, see
[`../../docs/PROJECT_WORKFLOW.md`](../../docs/PROJECT_WORKFLOW.md). For
what's current vs. archived (especially DepMap 24Q4 vs. 26Q1), see
[`../../docs/RESULTS_GUIDE.md`](../../docs/RESULTS_GUIDE.md).
