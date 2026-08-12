# Therapeutic shortlist freeze

**Date/time:** 2026-08-12, 21:40 UTC (revised after the Phase 20 Codex
review corrected a ranking-hierarchy implementation bug -- see below).
**Commit/HEAD used:** `fdd1a4446b8fab20f6637302165720490c14a1e0`
(candidate-adjudication phase; this freeze phase's own commit is layered
on top of it, see `docs/EVIDENCE_FREEZE_INPUT_AUDIT.md` for the exact
input files). Built by `src/evidence_freeze_shortlist_freeze.py`,
independently re-deriving eligibility and rank from
`results/tables/evidence_freeze/final_candidate_evidence.tsv` rather than
copying the candidate-adjudication phase's `shortlist_A_multimodal_therapeutic.tsv`
-- both methods converge on the identical 4-gene *membership*
(VEZF1/USP34/EML5/CITED2), a form of independent confirmation, though the
candidate-adjudication phase's List A used a different (also disclosed,
also banded) hierarchy and so is not expected to reproduce this list's
exact order.

## The frozen shortlist

| Rank | Gene | CRISPR effect/FDR | Full RNA pattern (↑\|↑\|↑ \|\| ↓, \* = FDR<0.05) |
|---|---|---|---|
| 1 | USP34 | −1.391 / 0.0417 | ↑\* \| ↑ \| ↑ \|\| ↓ |
| 2 | VEZF1 | −1.602 / 0.0373 | ↑ \| ↑\* \| ↓ \|\| ↓ |
| 3 | EML5 | −1.058 / 0.1488 | ↑\* \| ↑\* \| ↑ \|\| ↑ |
| 4 | CITED2 | −1.495 / 0.1100 | ↓ \| ↑\* \| ↑ \|\| ↑ |

USP34 and VEZF1 tie on the first two ranking criteria (both CRISPR-strength
band STRONG, both resistance-evidence band STRONG at 1/3 FDR<0.05); the
third criterion, resistance direction consistency, separates them --
USP34's resistance direction is fully concordant (`all_up`, all three
datasets point the same way) while VEZF1's is only `majority_up` (its
GSE111151 value points down, disagreeing with GSE118713/GSE240112). This
is why USP34 ranks first despite VEZF1 having a marginally lower CRISPR
FDR and (unlike USP34) significant human-tumor support -- direction
consistency is checked before human evidence in the declared hierarchy.

Full table: `results/tables/evidence_freeze/THERAPEUTIC_SHORTLIST_FREEZE.tsv`.

## Freeze rules (exact)

**Eligibility (hard gate, both required):**
1. CRISPR direction is `sensitising_KO` (negative effect size -- knockout
   relatively depleted under 4-OHT).
2. CRISPR evidence is real: FDR<0.25 (rejects a merely-negative-by-chance
   sign with no statistical support -- this excludes several
   resistance-leader genes that carry a `sensitising_KO` *label* with
   CRISPR FDR>0.7, e.g. GREB1, SLC4A10, DMRTA1, DLX2, GJB2, ACOT4, IL20).
3. At least one form of RNA/human corroboration: ≥1 resistance dataset
   (GSE118713/GSE240112/GSE111151) FDR<0.05, OR significant human-tumor
   support (GSE240112 or either GSE245601 track FDR<0.05). This is what
   distinguishes this list from the separate FUNCTIONAL_SENSITISATION
   list (class C, RNA support explicitly not required there) -- without
   it, classes A and C would be identical.

**Ranking** (banded comparison at every step, never a hidden numeric
score -- a tie within a band is expected and is exactly where the next
criterion decides the order): CRISPR-strength band → resistance-RNA
-evidence band → resistance direction consistency → human-tumor evidence
band → sample/replicate robustness band (GSE111151 cell-line consistency)
→ cross-dataset leave-one-out stability → gene symbol (deterministic
final tie-break only).

**12 candidates in the adjudicated pool met the sensitising-direction +
real-CRISPR-evidence gate; 4 also passed the multimodal-corroboration
gate.** The other 8 -- KDM1A, TLK2, TADA2B, USP17L29, HMGB1, EIF4ENIF1,
SUPT4H1, ICK -- are excluded from this list for lack of any resistance
-RNA or human-tumor support, and appear instead in class C,
`results/tables/evidence_freeze/frozen_candidate_classes.tsv`.

## Evidence sources

CRISPR functional screen (Hany-style MCF7-V, E2+4-OHT vs E2); GSE118713
(bulk RNA, TAMR vs parental MCF7); GSE240112 (human scRNA pseudobulk,
recurrent vs primary tumor, tumor-cell track); GSE111151 (bulk RNA,
tamoxifen-resistant vs parental, cell-line-blocked); GSE245601 (human
scRNA pseudobulk, acute 12h tamoxifen vs control, epithelial track as the
representative summary, malignant track also retained in the full
evidence table).

## Why GSE245601 is shown but excluded from the resistance consensus

GSE245601 measures an acute (12-hour) ex vivo tamoxifen exposure --
fundamentally a different biological question from the other three
datasets. GSE118713 and GSE111151 measure an established, chronic
resistant state reached after long-term drug exposure or clonal
selection; GSE240112 measures a recurrence-context *association*
(recurrent vs. primary tumor tissue) rather than a directly observed
selection history -- a real, meaningful distinction from acute exposure,
but not itself proof of a specific chronic-resistance mechanism. Treating
a 12-hour transcriptional response as equivalent evidence of "resistance"
would conflate acute pharmacodynamic response with the resistant/
recurrent phenotype these projects are ultimately about. GSE245601 is
displayed as the mandatory fourth arrow in every table and in the primary
RNA-comparison figures in this freeze (visually separated by a `||`
divider and, in
figures, by spacing/header text) so its evidence is never hidden -- it is
simply never summed into `resistance_fdr05_count`,
`resistance_direction_consistency`, or any resistance-specific band.

## Why tolerance-associated KO genes are excluded from the inhibition shortlist

CUX1, SOX2, TFAP2C, DPP9, and LZTR1 (5 of the 7 MULTIMODAL_STRONG genes)
all have a *positive* CRISPR effect size: their knockout was relatively
*favored*, not depleted, under 4-OHT (`tolerance_associated_KO`). This is
the opposite of the functional-fitness direction a knockout-inhibition
strategy aiming to increase tamoxifen sensitivity would need. Their
strong RNA-association evidence (CUX1 and SOX2 in particular have the
best resistance-RNA support of all seven MULTIMODAL_STRONG genes) is
real and valuable -- it is preserved and reported in class B (resistance
biomarker/pathway leaders) -- but expression association is not
functional-perturbation evidence, and RNA direction is never translated
into an inhibition recommendation without CRISPR direction support.

## Near-miss exclusions

The first genes excluded from the freeze, and why:
- **SPRED2** -- `tolerance_associated_KO` (wrong direction entirely).
- **TIPARP** -- `tolerance_associated_KO` (wrong direction).
- **MBTPS2** -- `tolerance_associated_KO` (wrong direction).
- **KDM1A** (and TLK2, TADA2B, USP17L29, HMGB1, EIF4ENIF1, SUPT4H1, ICK)
  -- correct (sensitising_KO) direction and, for KDM1A/TLK2, the
  *strongest* CRISPR evidence of any candidate in the pool (FDR
  0.0004/0.0017), but zero resistance-RNA or human-tumor support
  (0/3 FDR<0.05 everywhere) -- excluded by the multimodal-corroboration
  gate, not by direction. These remain the frozen class-C
  (functional-sensitisation) leaders.

## What was explicitly NOT used

No pathway analysis, no literature/mechanism review, no druggability or
tractability database, no known-inhibitor lookup, no normal-tissue
expression, no essentiality score, and no toxicity/safety information
entered the eligibility gate or the ranking. The freeze rests entirely on
the five already-frozen, already-Codex-reviewed experimental/
computational evidence layers.

## Status

**This shortlist is now frozen.** Any future change to its membership or
order requires an explicitly documented reason (a corrected bug, a newly
discovered value-provenance mismatch, or new experimental evidence) --
not a silent re-run with different parameters.
