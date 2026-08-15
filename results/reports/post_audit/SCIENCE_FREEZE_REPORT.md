# Science Freeze Report

**Status: authoritative final scientific-state document as of 2026-08-15.**
This report consolidates the original frozen analysis, the external
scientific audit, the post-audit sensitivity analysis, and the final
corrections made during the science-freeze pass. It supersedes no data —
it is a summary and pointer document.

---

## 1. Original historical conclusion

Under a predefined multimodal eligibility gate (CRISPR sensitisation +
at least one RNA/human-tumour corroboration dataset), the project froze a
four-gene therapeutic shortlist:

**USP34 > VEZF1 > EML5 > CITED2**

(`results/tables/evidence_freeze/THERAPEUTIC_SHORTLIST_FREEZE.tsv`,
`docs/THERAPEUTIC_SHORTLIST_FREEZE.md`). This result is **historically
preserved and byte-identical to the original commit** (Section 12 below)
— it remains valid as the answer to the specific, predefined question
that gate was designed to answer.

## 2. External-audit challenge

An external reviewer challenged the RNA-corroboration eligibility
requirement: a gene does not need to become transcriptionally
differentially expressed during acquired resistance for its knockout to
sensitise cells to tamoxifen, so the requirement may have excluded
stronger *functional* CRISPR hits — specifically named KDM1A and TLK2.
Secondary claims concerned GSE118713/GSE111151 design, GSE240112 pairing
language, GDSC's "9 significant rows," environment-specification gaps,
and a 2026 USP34 human-genetics paper. Every claim was independently
verified (`results/tables/post_audit_sensitivity/01_external_audit_claim_verification.tsv`);
all were confirmed correct except one genuine documentation error
(GSE240112 "matched" language — Section 6 below).

## 3. Post-audit sensitivity conclusion

A dedicated, separate analysis
(`results/reports/post_audit/POST_AUDIT_SENSITIVITY_REPORT.md`) found
**no single gene wins on every evidence axis**. Among the 13 genes
significant at the pre-registered Gate-1 threshold (FDR<0.1, sensitising
direction) in the Hany CRISPR screen:

- **KDM1A** ranks 1st of 13 by both CRISPR effect and FDR — the
  strongest functional sensitiser in the entire screen — but is null in
  every RNA dataset this project has.
- **TLK2** ranks 2nd of 13 by FDR — and has, by a wide margin, the
  strongest baseline ER+/luminal DepMap dependency **among the 4 focus
  genes** (81.8%, i.e. 9 of 11 dependency-evaluable lines, DepMap Public
  26Q1 -- provenance: 96 breast models in metadata -> 22 ER+/luminal ->
  53 breast / 11 ER+/luminal dependency-evaluable; SUPT4H1 is higher,
  90.9%, across the full 13-gene universe but is not one of the 4 focus
  genes). High dependency is not automatically an advantage for a
  tamoxifen-sensitisation strategy -- it may instead indicate a
  narrower, less tamoxifen-specific therapeutic window; this caveat
  applies equally to VEZF1's 27.3%.
- **USP34** ranks 12th of 13 by effect / 10th of 13 by FDR — the
  *weakest* CRISPR evidence of the four focus genes — but has exactly
  one significant non-acute (recurrence/resistance-context) RNA dataset
  (GSE118713, an established acquired-resistance cell-line model) and a
  real, crystallographically-proven covalently-reactive catalytic
  cysteine, demonstrated with a macromolecular activity-based ubiquitin
  probe (not a small-molecule inhibitor).
- **VEZF1** ranks 8th of 13 by effect / 6th of 13 by FDR, with exactly
  one significant non-acute RNA dataset (GSE240112 -- a recurrence
  *association* in unpaired human tumours, not an experimentally
  established resistance model; the only significant human-tumour
  signal of the four).
- GSE111151, the only dataset in this project with multiple independent
  resistance backgrounds, corroborates **zero** of the 12 testable
  significant sensitising hits (including all 4 focus genes) at
  FDR<0.05 — a real, disclosed null result.
- **KDM1A vs. USP34, pairwise (not a new ranking rule):** both are 0.0%
  strongly dependent in the same 11-line DepMap subset. KDM1A has
  substantially stronger CRISPR evidence and existing clinical-stage
  pharmacology; USP34's distinguishing advantages are novelty, its
  demonstrated catalytic-cysteine reactivity, and its single-model
  GSE118713 corroboration. Neither is stated to be superior overall.

## 4. Final candidate-role interpretation

| Role | Gene | Confidence |
|---|---|---|
| Strongest functional sensitiser | KDM1A | high — direct recomputation from `labels.parquet` |
| Most pharmacologically mature | KDM1A | high — multiple clinical-stage LSD1 inhibitors exist |
| Strongest baseline ER+/luminal dependency among the 4 focus genes | TLK2 | high — DepMap 26Q1, 11-line dependency-evaluable denominator; not automatically an advantage (see Section 3 above) |
| Strongest human recurrence-associated RNA signal | VEZF1 | moderate — single unpaired, small-n, source-confounded, recurrence-associated (not resistance-model) dataset |
| Target with the most direct evidence of covalent PROTEIN reactivity (activity probe, not a drug) | USP34 | moderate — hedged, not a proven comparative-tractability ranking, and distinct in kind from KDM1A's covalent small-molecule pharmacology; see Section 5 |
| Best-supported multimodal (CRISPR+RNA) target | Not a clean tie between USP34/VEZF1, not decidable | low — each has exactly 1 significant non-acute RNA dataset, but of different kinds/rigor (established resistance model vs. small unpaired recurrence association) — not interchangeable |
| **Overall universal winner** | **NO** | — |

## 5. Structural correction

The first draft of the post-audit report stated USP34's advantage as
"has a structure whereas KDM1A/TLK2 do not." **This was wrong and has
been corrected.** A facet-by-facet re-audit
(`results/tables/post_audit_sensitivity/06b_structural_tractability_audit.tsv`,
report Section 13A) independently verified (web search + RCSB PDB
fetch, 2026-08-15):

- **KDM1A**: many experimental structures, several inhibitor-bound
  co-crystals (e.g. PDB 6NQU with GSK2879552), multiple clinical-stage
  selective inhibitors (e.g. iadademstat/ORY-1001).
- **TLK2**: a real experimental kinase-domain structure exists (PDB
  5O0Y, 2.86Å) — bound to an ATP analog (AGS), **not** an inhibitor. No
  selective inhibitor exists; TLK1-paralog selectivity is an explicit
  unsolved medicinal-chemistry problem.
- **USP34**: real experimental catalytic-domain structures (7W3R apo,
  7W3U covalent-probe-bound), proven covalent reactivity at Cys1903 via
  a ubiquitin-propargylamide activity-based **probe** (not a drug). No
  selective small-molecule inhibitor exists.
- **VEZF1**: **no experimental structure exists** (a published
  inhibitor-screening study used a Zif268 homology model for lack of
  one). One real but weak (IC50=20µM), structure-unguided preliminary
  screening hit has been published — a correction to this analysis's
  earlier, too-strong claim of "no inhibitor program identified."

**Corrected precise USP34 claim, used throughout the README and reports
from this point forward:** *"a novel catalytic DUB target with
experimentally resolved covalent chemical addressability at Cys1903
(proven by a covalent activity-based probe, not an inhibitor) but no
validated selective small-molecule inhibitor identified."*

## 6. GSE240112 correction

Verified unpaired (3 primary + 3 recurrent ER+ tumours, different
patients, different source biobanks — OriGene for PT, Ontario Tumor Bank
for RT — no pairing statement in GEO metadata or the source paper). Two
occurrences of incorrect "matched primary/recurrent" language were found
and corrected this science-freeze pass: `docs/DATA_PROVENANCE.md` (found
and fixed in the prior post-audit phase) and **README.md's dataset table
(a second occurrence, missed by the prior pass, found and fixed in this
phase)**. A repo-wide regression test
(`tests/test_post_audit_sensitivity.py::TestDocumentationCorrection`)
now guards against a third recurrence. The frozen GSE240112 differential
-expression analysis itself was always correctly run as unpaired — only
documentation summary lines were ever wrong.

## 7. GDSC downgrade

`README.md` and the post-audit tables now state explicitly: 9
FDR-significant rows represent only **8 unique compounds** (AZD7762
counted twice, once per response metric), **all 9 are GDSC1 with zero
independent GDSC2 replication**, and 2 of the 7 LN_IC50-only hits
**reverse sign** in the AUC metric. Only AZD7762 (FDR=0.008, consistent
direction in both metrics) is a genuinely robust association. VEZF1 has
**zero** significant GDSC associations. GDSC is now explicitly labeled
exploratory/supporting evidence and explicitly **not** a reason USP34 was
selected, in both the README dataset table and narrative text.

## 8. USP34 safety update

Wigoda et al. 2026 (*Clinical Genetics*, DOI 10.1111/cge.70194 —
independently re-found via web search, matching this project's existing
PMID 42315110 citation) reports 6 individuals (5 confirmed de novo) with
heterozygous USP34 loss-of-function and a neurodevelopmental phenotype
(developmental delay, speech impairment, autism-spectrum features,
craniofacial dysmorphism), mediated through USP34's role stabilizing
Axin/Wnt-beta-catenin signalling in development. This is a **congenital,
lifelong, germline ~50%-dosage** phenotype. It is integrated into the
README and reports with the explicit, correct distinction: this is a
genuine, serious translational liability signal — it does **not** by
itself establish that adult, time-limited, partial pharmacological USP34
inhibition is unsafe. Neither the "ignore it" nor the "it proves danger"
reading is used anywhere in this repository's current text.

## 9. Reproducibility status

`environment.yml` now declares `gseapy`, `h5py`, `lifelines`, `networkx`,
`requests` (each independently re-verified this phase as genuinely
imported by real `src/`/`scripts/` code — not added speculatively) and
pins `pandas>=2,<3` with a documented rationale. Full test suite result:
see Section 11.

## 10. Codex final-review outcome

Independent adversarial review (read-only; verified live against RCSB PDB
for the TLK2 structural claim; overall verdict MEDIUM, no CRITICAL/HIGH).
Confirmed correct: the complete 13-gene sensitising universe; KDM1A/TLK2
fair representation; GSE245601 acute-only framing everywhere; DepMap
direction/denominators; README GDSC downgrade language; USP34
human-genetics framing (neither omitted nor overstated); no hidden
master/composite score. Valid MEDIUM/LOW issues found and **fixed this
pass**:

1. A second, surviving "recurrent tumor pairs" occurrence for GSE240112
   in `README.md`'s Limitations section (the prior pass's fix and
   regression test only covered `docs/DATA_PROVENANCE.md`) — fixed, and
   the regression test rewritten to sweep all markdown files at
   paragraph level (the old test only checked single physical lines and
   missed a wrapped continuation line).
2. GSE240112 was repeatedly labeled a "chronic-resistance" dataset in the
   post-audit report and README, equating a small, unpaired, source
   -confounded human recurrence *association* with an experimentally
   established resistance model (GSE118713/GSE111151) — reworded
   throughout to distinguish "resistance model" from "recurrence
   association" and never call GSE240112 "chronic resistance" evidence.
3. "Most structurally addressable novel target = USP34" and "best
   -supported multimodal target = tied USP34/VEZF1" were stated as clean
   superlatives/ties without acknowledging that (a) proven covalent
   reactivity to a macromolecular probe does not by itself establish
   comparative small-molecule tractability against TLK2's canonical
   kinase pocket, and (b) GSE118713 and GSE240112 are different kinds and
   rigor of evidence, not interchangeable "votes" — both hedged
   throughout the report, README, and this document.
4. A residual internal-consistency error: the audit-verification table
   (table 01) still said "zero of the 13" for GSE111151 corroboration
   instead of "zero of the 12 testable" (one gene, USP17L29, has no
   GSE111151 count-matrix entry) — the main report already had this
   right; the table's own text did not. Fixed.
5. A stale "39 passed" reproducibility count in the post-audit report,
   left over from before this science-freeze pass added 2 more
   regression tests — updated to 41.

No CRITICAL or HIGH issue was found or needed fixing.

## 11. Exact test result

`tests/test_post_audit_sensitivity.py`: **41 passed**. Full project
suite, run after all fixes in this report: **1,295 passed, 1 skipped, 0
failed** (`pytest -q`).

## 12. Historical frozen-file integrity

Verified via `git diff --stat` against every frozen result directory
(`evidence_freeze`, `cross_dataset_genomewide`, `independent_validation`,
`final_translational`, `final_pharmacogenomics`, `systems_network`,
`druggability_safety`, `lead_target_deep_dive`, `literature_mechanism`,
`candidate_adjudication`, the GSE-specific result directories,
`data/processed`) and the original gate code/pre-analysis plans
(`PREANALYSIS.md`, `GSE240112_PREANALYSIS.md`, `GSE111151_PREANALYSIS.md`,
`src/evidence_freeze_shortlist_freeze.py`, `src/labels.py`): **empty diff
— byte-identical to the repository's last commit (`9a4397f`) for every
one of these paths.** SHA-256 of the frozen shortlist file:
`b6990d7e20f1191df7ebbca18be7542cc312e6e8489f05af6e4b4d31e66fb9dc`
(matches the value independently computed by Codex against `HEAD` in the
prior post-audit review round). The post-audit analysis is an additional
sensitivity layer only — it has not silently replaced, edited, or
recomputed any part of the original historical analysis.

## 13. Known limitations that remain

- KDM1A and TLK2 were never carried into this project's structural/TCGA/
  pathway/GDSC/tissue-liability phases (those phases predate the
  post-audit analysis and were built specifically around USP34/VEZF1/
  EML5/CITED2). Their evidence base is CRISPR + DepMap + independently
  verified external literature only — a real, disclosed evidence-depth
  gap, not evidence they are weaker.
- USP34's and VEZF1's non-acute RNA corroboration is each a **single
  dataset** (USP34: GSE118713, an established acquired-resistance
  cell-line model; VEZF1: GSE240112, a recurrence *association* in
  unpaired human tumours, not an experimentally established resistance
  model -- not directly comparable in kind), not multiple independent
  confirmations; leave-one-dataset-out analysis shows each candidate's
  entire RNA support disappears if that one dataset is removed.
- A documentation tension between two GSE111151 pre-analysis plans (an
  earlier "no ordinary FDR, post-selection-confirmation-only" rule vs. a
  later, separately-dated, locked plan specifying a conventional blocked
  edgeR/FDR model, which is what the implemented eligibility gate
  actually consumes) was found and reported but **not resolved** — it
  requires the original analysis owner's intent, not a documentation
  typo fix, and is disclosed rather than silently adjudicated.
- A rebuild of a genuinely fresh conda environment from `environment.yml`
  was judged impractical within this project's session-based working
  model and was not performed; environment completeness was instead
  verified by an exhaustive static AST-level import audit (exact, not an
  estimate) cross-checked against the declared dependency list.
- No formal external-registry cross-check (e.g. ClinicalTrials.gov) was
  performed for the KDM1A/TLK2 pharmacology claims — they rely on
  literature/web-search synthesis, disclosed as such.
- Rule 4 of the selection-rule sensitivity analysis (human-evidence
  -first) was constructed, after reading the external audit, specifically
  to test that audit's own named hypothesis — this is disclosed
  explicitly in the report and is not presented as an independently
  -discovered neutral rule.

**No candidate discussed anywhere in this repository — USP34, VEZF1,
KDM1A, TLK2, EML5, or CITED2 — has been experimentally validated as a
therapy by this computational project. Every finding here is a
computational hypothesis awaiting wet-lab validation.**
