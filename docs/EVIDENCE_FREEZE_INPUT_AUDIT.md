# Evidence freeze: input audit

**Date written:** 2026-08-12. Lists the exact, already-committed candidate
-adjudication output files (commit `fdd1a44`) used as the sole inputs to
this evidence-freeze phase. All read-only here; nothing in this phase
reruns the CRISPR/scRNA/bulk upstream analyses.

## Source files located and used

| Purpose | File |
|---|---|
| Exact seven-gene evidence table | `results/tables/candidate_adjudication/multimodal7_exact_evidence.tsv` |
| Final candidate decision table (~30 genes) | `results/tables/candidate_adjudication/final_candidate_decision_table.tsv` |
| Final candidate summary (compact) | `results/tables/candidate_adjudication/final_candidate_summary.tsv` |
| **Provisional multimodal therapeutic shortlist** | `results/tables/candidate_adjudication/shortlist_A_multimodal_therapeutic.tsv` |
| Resistance-biomarker shortlist | `results/tables/candidate_adjudication/shortlist_B_resistance_biomarker.tsv` |
| Functional-sensitisation shortlist | `results/tables/candidate_adjudication/shortlist_C_functional_sensitisation.tsv` |
| Human-tumor shortlist | `results/tables/candidate_adjudication/shortlist_D_human_tumor.tsv` |
| Ranking stability (leave-one-out) | `results/tables/cross_dataset_genomewide/ranking_stability.tsv`, `results/tables/candidate_adjudication/candidate_rank_stability.tsv` |
| Near-miss multimodal genes | `results/tables/candidate_adjudication/multimodal_near_misses.tsv` |
| Sample-level tables/plots | `results/figures/candidate_adjudication/sample_level/*.png` (GSE118713, GSE111151, GSE240112, GSE245601 Track A/B) |
| Gene cards | `results/figures/candidate_adjudication/gene_cards/{USP34,VEZF1,CUX1,DPP9,LZTR1,SOX2,TFAP2C}.png` |
| Codex review record (candidate adjudication) | inline in the prior session; verdict PASS WITH NOTES after 2 material fixes (see conversation record; no separate doc file was written for that review) |
| Test results (candidate adjudication) | `tests/test_candidate_adjudication_*.py` (54 tests, part of 735 project-wide) |
| Category-audit doc | `docs/CANDIDATE_ADJUDICATION_CATEGORY_AUDIT.md` |
| Seven-gene head-to-head doc | `docs/CANDIDATE_ADJUDICATION_MULTIMODAL7.md` |
| GSE245601 influence doc | `docs/CANDIDATE_ADJUDICATION_GSE245601_INFLUENCE.md` |

## Upstream (frozen, read-only) sources these ultimately trace back to

`results/tables/cross_dataset_genomewide/all_genes_cross_dataset_evidence_with_ranking.tsv`
(the wide evidence matrix), `resistance_consensus_all_genes.tsv`,
`evidence_categories.tsv`, `crispr_functional_all_genes.tsv` -- all
committed in `2d7ab1d`, confirmed present in this branch's history
(Phase 0).

## What this phase does NOT touch

No file under `results/tables/cross_dataset_genomewide/`,
`results/tables/gse111151/`, `results/tables/gse240112_pseudobulk/`,
`results/tables/gse245601_pseudobulk/`, or the candidate-adjudication
directory itself is modified. All evidence-freeze outputs go into new
`results/tables/evidence_freeze/`, `results/figures/evidence_freeze/`,
and `docs/THERAPEUTIC_SHORTLIST_FREEZE.md`.
