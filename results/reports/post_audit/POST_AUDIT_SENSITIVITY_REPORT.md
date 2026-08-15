# Post-Audit Sensitivity Report

**Status:** a separate, new analysis. Does not modify, recompute, or
supersede the frozen `THERAPEUTIC_SHORTLIST_FREEZE.tsv`
(USP34 > VEZF1 > EML5 > CITED2), which remains historically valid as "the
result under the pre-defined multimodal gate." This report tests how that
conclusion changes under alternative, transparent selection rules.

All source code: `src/post_audit_sensitivity_data.py`,
`src/post_audit_sensitivity_visualization.py`. All tables:
`results/tables/post_audit_sensitivity/01`–`09` (plus `06b`, added in the
2026-08-15 structural-tractability re-audit). All figures:
`results/figures/post_audit/A`–`D`. Tests: `tests/test_post_audit_sensitivity.py`.

**2026-08-15 update:** this report's structural-tractability claims were
re-audited and corrected (Section 13A) after a follow-up review found the
first draft's USP34-structure framing was imprecise (it implied USP34
uniquely had a structure; KDM1A and TLK2 also have real experimental
structures). GSE240112, GDSC, and USP34-safety framing were also swept
project-wide for remaining imprecise language; see
`results/reports/post_audit/SCIENCE_FREEZE_REPORT.md` for the
consolidated final state.

---

## 1. What did the external reviewer claim?

That the original candidate-freeze logic required functional CRISPR
sensitisation **plus** transcriptomic corroboration, and that a gene does
not biologically need to become differentially expressed during acquired
resistance for its knockout to sensitise cells to tamoxifen — so the gate
may have excluded stronger functional hits, specifically **KDM1A** and
**TLK2**. Secondary claims: GSE118713 provides only one parental/resistant
model rather than independent resistance backgrounds; GSE111151 (which
does have independent backgrounds) should carry more weight; GSE240112 may
be mischaracterized as "matched" when it is unpaired; the GDSC
pharmacogenomics result ("9 significant rows") may overstate independent
replication; the environment specification has real gaps (h5py, pandas);
and a 2026 human-genetics paper reports USP34 heterozygous
loss-of-function associated with a neurodevelopmental phenotype.

## 2-3. Which claims were verified / which were incorrect or overstated?

Every quantitative claim independently checked was **confirmed correct**
(table `01_external_audit_claim_verification.tsv`, 19 rows): the exact
Hany effect/FDR for all 4 focus genes; the sign convention; that RNA
corroboration is a hard eligibility gate in the implemented code
(`src/evidence_freeze_shortlist_freeze.py`); that KDM1A/TLK2 were excluded
specifically for lacking that corroboration despite `VERY_STRONG` CRISPR
evidence; the GSE118713 single-derivation-event design; GSE111151's
4-independent-background design; the GDSC 9-rows/8-compounds/GDSC1-only/
2-sign-flips pattern; the h5py/pandas/gseapy/lifelines/networkx/requests
environment gaps; and the USP34 neurodevelopmental paper (independently
re-found via web search: "USP34 Haploinsufficiency as a Cause of
Neurodevelopmental Phenotypes," Wigoda et al., *Clinical Genetics* 2026,
DOI 10.1111/cge.70194 — matches the project's own PMID 42315110 citation).

**One claim was found to be a genuine documentation error, not a claim
about the analysis itself**: `docs/DATA_PROVENANCE.md` described GSE240112
as "matched primary/recurrent ER+ tumor pairs." The project's own locked
`GSE240112_PREANALYSIS.md` already correctly documents this dataset as
unpaired (different source biobanks, no pairing statement in GEO
metadata), and the frozen GSE240112 differential-expression analysis was
already run correctly as unpaired — only documentation summary lines were
wrong. **Corrected this phase** (see Section 17 below). A **second**
occurrence of the identical error was later found in `README.md`'s
dataset table during the 2026-08-15 science-freeze repo-wide sweep (the
first pass's fix, and its regression test, checked `docs/DATA_PROVENANCE.md`
specifically and did not catch this second file) — also corrected, and
the regression test broadened to a repo-wide, paragraph-level sweep with
its own test to guard against a third recurrence.

No claim was found to be fabricated or exaggerated.

## 4. What exactly did the original freeze rule do?

`src/evidence_freeze_shortlist_freeze.py::determine_eligibility()`
(quoted): a gene enters the shortlist only if **all** of:
1. `crispr_direction == "sensitising_KO"`
2. CRISPR FDR<0.25 ("real evidence", not just a noisy sign)
3. `resistance_fdr05_count >= 1` **OR** `human_tumor_support == "significant"`
   (`some_multimodal_support`, computed from `RESISTANCE_DATASET_ORDER =
   ["gse118713", "gse240112", "gse111151"]` FDR<0.05, or GSE240112/GSE245601
   FDR<0.05 for the "human_tumor_support" branch)

Criterion 3 is a **hard eligibility gate**, not merely a ranking
criterion, despite the module's own docstring opening with "RNA
significance is explicitly NOT part of the eligibility gate" — that
sentence is inaccurate relative to the actual boolean logic two lines
later in the same docstring and the code itself (a real internal
documentation inconsistency, flagged but not edited, since the file is
part of the frozen original analysis).

Among genes passing eligibility, the ranking hierarchy is: CRISPR-strength
band → resistance-RNA band → direction consistency → human-tumor band →
GSE111151 sample-robustness band → leave-one-out stability → gene symbol
(alphabetical tie-break only). This hierarchy is transparent and
band-based (never a hidden numeric score) — verified against
`rank_eligible()`.

**Classification of the rule's components:**
- **PRE-SPECIFIED RULE** (PREANALYSIS.md, locked before results seen):
  Gate-1 CRISPR FDR<0.1 (§4); the sign convention (§2); the sensitising
  direction requirement follows directly from the project's inhibition
  therapeutic-strategy framing (§1).
- **BIOLOGICAL/DESIGN ASSUMPTION** (introduced later, at the evidence-freeze
  phase, not in the original locked plan): the "some multimodal
  corroboration" hard gate (criterion 3 above). Its own docstring gives an
  explicit, disclosed rationale (distinguishing this list, "List A," from
  the project's separately-computed "List C" functional-only list) — this
  is a reasonable, defensible, but genuinely debatable design choice, not
  a pre-registered rule and not the only reasonable choice.
- **BOOKKEEPING/WORKFLOW RULE**: the CRISPR-FDR<0.25 "real evidence" floor
  (looser than the original Gate-1 FDR<0.1 — meaning EML5 and CITED2, at
  FDR 0.149 and 0.110, enter the freeze process despite not being official
  Gate-1 hits under the pre-registered FDR<0.1 threshold); the band
  boundaries themselves; the tie-break ordering.

A secondary, unresolved documentation tension was found and is reported
(not corrected, since it requires the original author's intent, not a
typo fix): a 2026-08-06 `PREANALYSIS.md` amendment states GSE111151 "will
no longer contribute discovery features" and "does not become a
feature/discovery input," to be judged only by a lineage-consistency rule,
"no ordinary ... FDR ... calculated." A later, separately-dated
`docs/GSE111151_PREANALYSIS.md` (2026-08-12) instead specifies a
conventional blocked edgeR model with real FDR — and it is this real FDR
that the freeze code's `resistance_fdr05_count()` actually consumes as
part of the hard eligibility gate. The two documents were not
cross-referenced to each other; this project's own conclusion should
flag this to the analysis owner rather than have this sensitivity
analysis silently pick a side.

## 5. Why were KDM1A and TLK2 excluded originally?

Verified exactly (`freeze_eligibility_audit.tsv`): both have
`_crispr_band = VERY_STRONG` (the strongest possible CRISPR band) and
`ineligibility_reason = "no resistance-RNA dataset FDR<0.05 and no
significant human-tumor support (functional-only -- see List C
instead)"`. Neither failed on CRISPR strength or direction. Both failed
criterion 3 alone.

This project's **own, separately-computed, pre-existing** List C
("functional sensitisation," RNA not required — `results/tables/
candidate_adjudication/shortlist_C_functional_sensitisation.tsv`) already
ranked KDM1A #1 and TLK2 #2 by CRISPR strength, with USP34 not present in
its top 5 — this was computed before the current external audit, in an
earlier phase of this same project, and independently corroborates the
finding.

## 6. Was that exclusion scientifically necessary, or one reasonable design choice?

**One reasonable, disclosed design choice — not a scientific necessity.**
The rationale (distinguish a "therapeutic, multimodally-corroborated"
list from a "functional-only" list) is legitimate and was never hidden —
List C existed the whole time as the un-gated alternative. But nothing in
the pre-registered PREANALYSIS.md required RNA corroboration as an
eligibility bar; it was introduced later, at the evidence-freeze phase.
Sections 7-19 below quantify exactly what changes when this one design
choice is relaxed.

## 7. What happens when RNA differential expression is removed as an eligibility requirement?

Under Rule 1 (`04_selection_rule_sensitivity.tsv`, `rule=
RULE_1_crispr_only_no_rna_gate`): the eligible universe becomes all 13
significant sensitising Gate-1 hits, ranked by CRISPR-strength band then
FDR. **KDM1A ranks 1st, TLK2 ranks 2nd, VEZF1 ranks 6th, USP34 ranks
10th** (of 13). USP34 falls to the bottom third of its own screen's
significant sensitising hits once RNA is not required for eligibility.

## 8. Functional-first ranking (View A)

Reported by effect and FDR **separately**, since they answer different
questions (Figure A; `02_significant_sensitising_crispr_hits.tsv`):

| Gene | Effect rank (of 13) | FDR rank (of 13) |
|---|---|---|
| KDM1A | 1 | 1 |
| TLK2 | 4 | 2 |
| VEZF1 | 8 | 6 |
| USP34 | 12 | 10 |

KDM1A is the single strongest sensitising hit in the entire screen by
both metrics. TLK2 is 2nd by FDR. USP34 is near the bottom of the
significant-sensitising set by both metrics.

## 9. Multimodal evidence comparison (View B)

Figure B / `03_post_audit_evidence_matrix.tsv`. Of all 13 significant
sensitising hits, **only USP34 and VEZF1 have any non-acute (chronic
-context) RNA dataset significant at FDR<0.05 — and each has exactly
one, of two evidentially different kinds**: USP34 → GSE118713 only
(FDR=0.0073; an established acquired-resistance cell-line model;
GSE240112 FDR=0.228, GSE111151 FDR=0.632, both null). VEZF1 → GSE240112
only (FDR=0.0195; a recurrence-*association* in human tumors, NOT an
experimentally established resistance model; GSE118713 FDR=0.238,
GSE111151 FDR=0.608, both null). Neither is a "multiple independent
models" case — both are "one supporting dataset, of a different
evidential kind from each other" — and Section 18 (leave-one-out) shows
each candidate's entire non-acute RNA support disappears if that one
dataset is removed. KDM1A and TLK2 have **zero** non-acute-dataset
FDR<0.05 hits (all 3 non-acute datasets null for both), and are also null
(not significant) in
GSE245601 — i.e. their transcriptomic signal (if any) is genuinely absent
across every RNA dataset this project has, not merely non-significant in
one.

## 10. GSE111151-specific result (View C)

Figure D / `05_leave_one_dataset_out.tsv` / Rule 3. GSE111151 is this
project's only dataset with multiple independent resistance backgrounds
(4 parental lines, 7 independently-derived TamR sublines, cell-line
-blocked model). **Zero of the 12 GSE111151-testable significant
sensitising Gate-1 hits** (USP17L29 has no GSE111151 count-matrix entry
and is untestable, not tested-and-null -- 12 of the 13, not all 13, are
actually evaluable here; a Codex adversarial review of this analysis
caught this precision gap), including all 4 focus genes, reach GSE111151
FDR<0.05. This is reported as a real, honest null result, not suppressed: the best-designed
resistance dataset in this project does not independently corroborate any
CRISPR sensitising hit at conventional significance, KDM1A/TLK2/USP34/
VEZF1 included.

## 11. Human-evidence comparison (View D)

GSE240112 (real patient tumors, though unpaired/source-confounded) is
significant only for VEZF1 among the 4 focus genes (FDR=0.019); TCGA
paired tumor-vs-normal is significant for none of USP34/VEZF1 at FDR<0.05
in this specific contrast (USP34 FDR=0.21, VEZF1's ER+/ER- contrast is
significant but that is a different comparison than tumor-vs-normal — see
`06_top_candidate_translational_comparison.tsv`). KDM1A/TLK2 have no TCGA
analysis in this project (evidence gap, not a null result) but their only
human-tumor RNA source available (GSE240112) is also null for both. No
candidate has strong, unambiguous human-tumor transcriptomic support;
VEZF1 has the most (one significant human dataset), USP34 has none.

## 12. DepMap comparison

Figure C / `03_post_audit_evidence_matrix.tsv`. **Provenance:** DepMap
Public 26Q1, `CRISPRGeneDependency.csv`. Of 96 breast models in DepMap's
own metadata, 22 are ER+/luminal by this project's `ModelSubtypeFeatures`
rule; of those, 53 breast models (11 ER+/luminal) actually have a
Chronos/dependency-probability value in the release's gene-effect matrix
("dependency-evaluable") — the **n=11 ER+/luminal dependency-evaluable
lines** is the denominator behind every percentage below, not the raw 22.

Baseline ER+/luminal (n=11) dependency, dependency-probability>0.5
fraction: **TLK2 81.8% (9/11)**, HMGB1 81.8% (9/11), SUPT4H1 90.9% (10/11,
the highest across the full 13-gene significant-sensitising universe --
not one of the 4 focus genes and not promoted to a lead here), **VEZF1
27.3% (3/11)**, **KDM1A 0.0% (0/11)**, **USP34 0.0% (0/11)**. TLK2 has,
by a wide margin, the **strongest baseline cancer-cell dependency among
the 4 focus genes specifically** (SUPT4H1 is higher across the broader
13-gene universe, a scope distinction, not a new candidate claim).

**High baseline dependency is not automatically a translational
advantage for a tamoxifen-sensitisation strategy.** It can indicate two
different things, and this project's data do not distinguish between
them: (a) a genuinely strong baseline anticancer vulnerability (a
positive), or (b) a narrower context-specific therapeutic window --
i.e. the gene may already matter broadly for cancer-cell fitness
independent of tamoxifen, making a tamoxifen-*specific* combination
strategy less differentiated and raising normal-tissue-toxicity risk
proportionally. This applies to **both TLK2 (81.8%) and VEZF1 (27.3%)**
and should not be read as a purely favourable signal for either; it is a
real, quantitative "dual-action" hypothesis (functional sensitiser +
baseline vulnerability) that neither USP34 nor KDM1A show, stated as a
hypothesis requiring follow-up (e.g. the normal-tissue comparators in the
proposed experimental design), not a proven advantage.

## 13. Translational/druggability comparison (View E)

`06_top_candidate_translational_comparison.tsv`, `07`/`08` deep dives.
KDM1A: existing clinical-stage selective inhibitor (iadademstat/ORY-1001,
Orphan Drug Designation AML/SCLC, not yet FDA-approved; independent
breast-cancer-stem-cell preclinical literature; independent literature
links LSD1 reduction to preventing resistance in breast cancer models) —
most translationally mature of the 4. TLK2: independent literature
(amplified ~10.5% of ER+ tumors, worse outcomes, mechanistic knockdown
data improving PFS in vivo) but **no potent/selective inhibitor exists**
(TLK1-paralog selectivity is an unsolved medicinal-chemistry problem) —
least translationally mature. USP34: real experimental structures
(7W3R/7W3U), proven covalent reactivity at Cys1903, but no selective
inhibitor and docking explicitly judged not-yet-justified in this
project. VEZF1: a C2H2 zinc-finger transcription factor — a target class
historically much harder to drug directly than any of the 3 enzymes
compared here; no validated selective inhibitor, though a real, weak
(IC50=20 micromolar), structure-unguided preliminary screening hit (T4/
503-1-83) has been published (corrected in Section 13A below — the first
draft of this report understated this to "no inhibitor program
identified" at all).

## 13A. Structural tractability re-audit (2026-08-15 correction)

A follow-up, more precise re-audit was requested after the first draft of
this report stated USP34's advantage as simply "has a structure whereas
KDM1A/TLK2 do not." **That framing was wrong and has been corrected
throughout this report.** `06b_structural_tractability_audit.tsv` records
7 separate facets (A-G) per gene, each independently sourced, rather than
one collapsed boolean:

| | A. exptl structure | B. relevant domain solved | C. ligand/probe-bound structure | D. druggable pocket established | E. validated selective inhibitor | F. clinical-stage pharmacology |
|---|---|---|---|---|---|---|
| KDM1A | YES | YES | YES, several inhibitor co-crystals (e.g. PDB 6NQU) | YES | YES, several (e.g. iadademstat) | YES, 9 compounds |
| TLK2 | **YES (PDB 5O0Y)** | partial (kinase domain only) | partial -- bound to an ATP analog (AGS), NOT an inhibitor | plausible, unconfirmed by an inhibitor co-crystal | NO | NO |
| USP34 | YES (7W3R/7W3U) | partial (catalytic domain only) | YES, a covalent activity-based PROBE (not a drug), proven reactivity at Cys1903 | partial -- large/groove-shaped fpocket hit, not a compact drug pocket | NO | NO |
| VEZF1 | **NO** | NO (homology-modeled only) | NO | NO | NO, but 1 weak (IC50=20uM) published screening hit | NO |

**Corrected conclusion:** KDM1A, TLK2, and USP34 all have real
experimental structures — "has a structure" does not distinguish USP34.
The precise, corrected USP34 distinction is: *"a novel catalytic DUB
target with experimentally resolved covalent chemical addressability at
Cys1903 (proven by a covalent activity-based probe, not an inhibitor) but
no validated selective small-molecule inhibitor identified."* Only VEZF1
genuinely lacks any experimental structure. Sources: RCSB PDB (5O0Y,
6NQU, 2Z5U, web-verified 2026-08-15), PMC6100598 (VEZF1 screening hit),
this project's own frozen `USP34_pocket_analysis.tsv`/
`USP34_docking_decision.tsv`.

## 14. KDM1A deep dive

`07_KDM1A_deep_dive.tsv`. Strongest CRISPR hit in the entire screen
(effect/FDR both rank 1 of 13). Null in every RNA dataset this project
has (GSE118713/GSE111151/GSE240112/GSE245601 all non-significant). No
baseline ER+/luminal DepMap dependency (0%). Was one of the two
pre-registered blind positive controls (never inspected until the
2026-08-10 amendment) — its recovery as the top hit is itself validity
evidence for the screen, exactly as PREANALYSIS.md §10 intended.
Independent literature ties LSD1/KDM1A to ER-complex biology and
endocrine-resistance-prevention in breast cancer models — this is *why*
it was chosen as a blind control, not a new finding. Existing, clinically
advanced, selective inhibitor (iadademstat). Class-level LSD1-inhibitor
toxicity concerns (e.g. hematologic effects) are a known generic class
liability, not independently re-verified against USP34/VEZF1-level
GTEx/HPA rigor in this project. Novelty as a *tamoxifen-response*
modifier is low (its ER-pathway role was already well-established prior
literature); novelty as *this project's specific finding* is that the
unbiased screen recovered it as the single strongest hit.

## 15. TLK2 deep dive

`08_TLK2_deep_dive.tsv`. 2nd-strongest hit by FDR, 4th by effect. Null in
every RNA dataset this project has, same pattern as KDM1A. **Strongest
baseline ER+/luminal DepMap dependency of any of the 4 focus genes
(81.8%)** — a real, distinct signal neither USP34 nor KDM1A shows.
Independent literature: amplified in ~10.5% of ER+ tumors, worse outcomes
regardless of endocrine-therapy status, mechanistic knockdown data
(downregulates ERα/BCL2/SKP2, improves PFS in vivo, blocks metastasis in
one model) — substantive external support, not from this project's
screen. No potent/selective inhibitor exists; TLK1-paralog homology is an
explicit, unsolved selectivity problem in the medicinal-chemistry
literature. Toxicity/normal-tissue liability not independently assessed
in this project. Genuinely higher novelty than KDM1A as a
*tamoxifen-response* finding specifically (TLK2's prior literature is
about amplification/prognosis, not acute drug-context modulation).

## 16. Updated USP34 safety/liability assessment

`09_USP34_updated_liability.tsv`. The Wigoda et al. 2026 neurodevelopmental
paper is real and independently re-verified (Section 2-3 above; already
correctly documented in this project's frozen
`USP34_VEZF1_tissue_liability.tsv`). It reports a **congenital, lifelong,
germline ~50%-dosage** phenotype (6 individuals, 5 confirmed de novo
heterozygous LoF: developmental delay, speech impairment, autism-spectrum
features, craniofacial dysmorphism), driven through USP34's role
stabilizing Axin/canonical Wnt-beta-catenin signalling in development.
**This does not automatically mean adult pharmacological USP34 inhibition
is unsafe** — a time-limited, adult-onset partial inhibition is not the
same exposure as a lifelong germline dosage reduction spanning
neurodevelopmental windows — but it is a genuine, serious, previously
under-weighted translational liability signal, materially stronger than a
normal-tissue-expression finding alone, and should be disclosed
prominently (not just listed) in any USP34 translational-strategy
document going forward.

## 17. Selection-rule sensitivity results

`04_selection_rule_sensitivity.tsv`, Figure D. Summary:

| Rule | Eligible (of 13+original-4 universe) | KDM1A | TLK2 | USP34 | VEZF1 | Lead changes? |
|---|---|---|---|---|---|---|
| 0 original frozen gate | 4 | ineligible | ineligible | **rank 1** | rank 2 | — (baseline) |
| 1 CRISPR only, no RNA gate | 13 | **rank 1** | rank 2 | rank 10 | rank 6 | YES -- KDM1A leads |
| 2 non-acute RNA corroboration | 2 | ineligible | ineligible | rank 2 | **rank 1** | YES -- VEZF1 leads (of the 2 eligible) |
| 3 GSE111151-specific | 0 | ineligible | ineligible | ineligible | ineligible | N/A -- empty set |
| 4 human evidence first | 13 | **rank 1** | rank 2 | rank 4 | rank 3 | YES -- KDM1A/TLK2 lead; VEZF1 outranks USP34 |

The external reviewer's specific claim — that VEZF1 may outrank USP34
once human evidence is ordered before cell-line RNA consistency — is
**confirmed under Rule 4** (`TestRule4HumanEvidenceFirst`, pinned).

**Rule 4 coverage limitation (added 2026-08-16):** Rule 4's "human
evidence" tie-break uses GSE240112 (available genome-wide, all 13 genes)
and TCGA (`tcga_fdr`). TCGA follow-up in this project was originally run
**only on the frozen four-candidate set** (USP34/VEZF1/EML5/CITED2) --
9 of the 13 significant sensitising genes, including KDM1A and TLK2, have
no TCGA column at all (`NaN`, not "not significant"; see
`03_post_audit_evidence_matrix.tsv`). **Rule 4 is therefore a sensitivity
test of one specific, disclosed reordering, not a genome-wide, unbiased
human-evidence ranking** — its result for KDM1A/TLK2 rests on GSE240112
alone (both null there), never on a TCGA comparison that was never run
for them.

**Full transparency on Rule 4's construction** (raised by the Codex
adversarial review, Section "K" below): Rule 4 was deliberately built,
after reading the external audit, specifically to test that reviewer's
named hypothesis — it is a purpose-built sensitivity probe, not a
pre-registered rule discovered independently of the audit. Its "human
evidence first" name is also an approximation: CRISPR-strength band is
still the first sort key (this is disclosed in the function's own
docstring), with human evidence entering only as the tie-break within a
CRISPR band. Within the STRONG band that USP34 and VEZF1 both occupy,
this is exactly the reordering the reviewer proposed and the effect is
real — but the rule was not blind to the question it was built to answer,
and should be read as "does this specific, named reordering produce the
result the reviewer predicted" (yes) rather than as an independently
-discovered, neutral alternative ranking.

## 18. Leave-one-dataset-out results

`05_leave_one_dataset_out.tsv`. USP34's non-acute RNA corroboration
(GSE118713-only, an established resistance model) disappears entirely if
GSE118713 is excluded. VEZF1's non-acute RNA corroboration (GSE240112
-only, a recurrence association, not a resistance model) disappears
entirely if GSE240112 is excluded. Both candidates' RNA support is
therefore **single-dataset-dependent** — neither survives removal of its
one supporting dataset. KDM1A and TLK2 have no non-acute RNA
corroboration to lose under any leave-one-out scenario (they start at
zero).

## 19. Does USP34 remain...?

- **strongest functional sensitiser?** **NO.** 12th of 13 by effect, 10th
  of 13 by FDR (Section 8).
- **strongest multimodal candidate?** **Tied/ambiguous with VEZF1, not a
  clean "yes."** Both have exactly one significant non-acute RNA dataset
  (USP34: GSE118713, a resistance model; VEZF1: GSE240112, a recurrence
  association, not a resistance model -- not directly comparable in
  kind); USP34 additionally has a measured TCGA paired tumor-vs-normal trend
  that does NOT reach FDR<0.05 (FDR=0.21 -- not statistical evidence of a
  real difference, corrected wording after Codex review flagged the
  original draft's "signal" language as disproportionate) and the only
  existing GDSC pharmacogenomic hit; VEZF1 has the only
  significant human-tumor (GSE240112) hit. Neither dominates the other on
  every multimodal axis (Section 9, 11, 20).
- **strongest structurally addressable NOVEL target?** **YES, with a
  corrected, precise basis** — a structural-tractability re-audit (Section
  21A below; `06b_structural_tractability_audit.tsv`) found that KDM1A,
  TLK2, AND USP34 all have real experimental structures — the original
  draft of this report was wrong to imply USP34 was the only one with a
  structure, and that specific sentence has been corrected. The accurate,
  narrower distinction that DOES hold: USP34 is a novel catalytic DUB
  target with experimentally resolved covalent chemical addressability at
  Cys1903 (proven by a covalent activity-based probe, not an inhibitor)
  but no validated selective small-molecule inhibitor identified — KDM1A's
  structures are already extensively inhibitor-bound and clinically
  mature (low remaining novelty), while TLK2's one structure (PDB 5O0Y) is
  bound only to an ATP analog, not an inhibitor, with real open
  medicinal-chemistry opportunity but no chemical starting point as direct
  as USP34's proven covalent reactivity.
- **strongest translational target overall?** **NO** — KDM1A has a more
  translationally mature existing pharmacology (clinical-stage selective
  inhibitor); TLK2 has a stronger baseline dependency signal; USP34's
  translational case rests specifically on structural/covalent
  tractability plus a real (if single-dataset) resistance signal, not on
  functional CRISPR strength or existing pharmacology.

## 20. Does VEZF1 become stronger under any reasonable criterion?

**Yes, under Rule 4 (human-evidence-first)**, VEZF1 outranks USP34 —
directly confirming the external reviewer's specific hypothesis. Under
Rule 2 (non-acute-RNA-only eligibility), VEZF1 also ranks 1st of the 2
eligible genes (ahead of USP34). VEZF1 does not become stronger than
KDM1A or TLK2 under any rule tested (CRISPR band always favors them where
they are eligible).

## 21. Is KDM1A stronger under any reasonable criterion?

**Yes, decisively, under any rule that does not hard-gate on RNA/human
corroboration** (Rules 1 and 4): KDM1A is the single strongest CRISPR hit
in the entire screen and has the most translationally mature existing
pharmacology of the 4 focus genes. It is weaker only under rules that
require RNA/human corroboration (Rules 0, 2, 3), where it is ineligible
because — verified — it has no significant signal in any RNA dataset this
project has tested.

**KDM1A vs. USP34, explicit pairwise interpretation (not a new ranking
rule):** both show 0.0% strong ER+/luminal baseline dependency in the
evaluated DepMap subset (11 dependency-evaluable lines) — indistinguishable
on that specific axis. KDM1A has substantially stronger CRISPR
sensitisation evidence (effect rank 1/13 vs. USP34's 12/13; FDR rank 1/13
vs. 10/13) and existing clinical-stage pharmacology (iadademstat).
USP34's distinguishing advantages are novelty as a target, its
experimentally demonstrated catalytic-cysteine (Cys1903) reactivity to a
covalent activity probe, and its single-model GSE118713 RNA
corroboration (which KDM1A entirely lacks). **This is not a claim that
USP34 is superior overall, nor that KDM1A universally wins** — each
advantage is stated on its own named axis.

## 22. Is TLK2 stronger under any reasonable criterion?

**Yes, on two specific, real axes**: functional CRISPR strength (2nd of
13 by FDR) and baseline ER+/luminal cancer-cell dependency (81.8%,
strongest **among the 4 focus genes** by a wide margin — SUPT4H1 is
higher, 90.9%, across the full 13-gene universe but is not one of the 4
focus genes). This dependency signal is a real "dual-action" hypothesis
(functional sensitiser + baseline vulnerability), **not automatically a
positive** — high baseline dependency can equally indicate a narrower,
less tamoxifen-specific therapeutic window, a caveat that applies
symmetrically to VEZF1's 27.3%. It is weaker on RNA/human corroboration
(identically null to KDM1A) and markedly weaker on existing pharmacology
(no selective
inhibitor exists, an unsolved TLK1-paralog selectivity problem).

## 23. Is there a single scientifically defensible "winner"?

**No.** Section 24 defines the distinct, non-dominated roles the evidence
actually supports.

## 24. Distinct candidate roles supported by the evidence

- **Strongest functional sensitiser:** KDM1A (Section 8, 14).
- **Strongest baseline cancer-cell vulnerability among the four focus
  genes:** TLK2 (Section 12, 15); SUPT4H1 is higher (90.9%) in the full
  13-gene universe but is not one of the four focus genes and is not
  promoted to a lead here. Not automatically an advantage -- see Section
  12/21/22.
- **Strongest single-dataset non-acute RNA corroboration:** USP34
  (GSE118713, an established acquired-resistance cell-line model) and
  VEZF1 (GSE240112, a recurrence *association* in unpaired human tumors —
  a different, weaker evidential kind, not an experimental resistance
  model). These are **not a clean tie** — each is supported by exactly
  one dataset, not multiple, but the two datasets differ substantially in
  design and rigor and should not be read as interchangeable
  (Section 9, 18; Codex final-review finding, corrected here).
- **Strongest human-recurrence-tumor signal:** VEZF1 (Section 11).
- **Most translationally mature existing pharmacology:** KDM1A
  (iadademstat, clinical-stage) (Section 13, 14).
- **Novel target with the most direct experimental evidence of covalent
  chemical reactivity, no existing inhibitor (hedged, not a clean
  "winner"):** USP34 -- precise wording after the Section 13A structural
  re-audit; KDM1A and TLK2 also have real experimental structures, so
  "has a structure" alone does not distinguish USP34. Proven covalent
  reactivity of Cys1903 to a macromolecular activity-based probe is real
  evidence of catalytic-cysteine reactivity, but does **not** by itself
  establish that USP34 is more comparatively small-molecule-tractable
  than TLK2's canonical ATP-competitive kinase pocket -- a Codex
  final-review finding, corrected here (Section 19, 13A).
- **Highest-novelty hypothesis specific to acute tamoxifen-response
  modulation (as opposed to general breast-cancer relevance):** TLK2
  (Section 15).
- **Weakest overall functional evidence among the 4 focus genes:** USP34
  (bottom third of the screen's own significant sensitising hits).

## 25. What exact conclusion should the poster use now?

See Section M of the final response below (Poster Consequence) —
summary: the poster should not claim USP34 is "the strongest CRISPR hit"
or "the best-supported target" in any unqualified sense; it should either
present USP34 specifically as the structurally-tractable lead with an
explicit, disclosed CRISPR-strength caveat, or restructure around multiple
complementary roles (Section 24) rather than a single ranked winner.

---

## Final post-audit conclusion, by explicit category (2026-08-15)

| Category | Answer | Basis |
|---|---|---|
| Strongest functional sensitiser | **KDM1A** | rank 1/13 by both CRISPR effect and FDR (Section 8) |
| Most pharmacologically mature | **KDM1A** | existing clinical-stage selective inhibitor, iadademstat (Section 13, 13A) |
| Strongest baseline ER+/luminal dependency **among the 4 focus genes** | **TLK2** | 81.8% (9/11) DepMap 26Q1 dependency-evaluable lines, dependency-probability>0.5, highest of the 4 focus genes (SUPT4H1 is higher, 90.9%, across the full 13-gene universe). Not automatically an advantage -- may indicate a narrower, less tamoxifen-specific therapeutic window (Section 12) |
| Strongest human recurrence-associated signal | **VEZF1** | only significant GSE240112 hit among the 4 (Section 11) -- an association with recurrence, not proof of causing or reversing tamoxifen resistance |
| Novel target with the most direct evidence of covalent chemical reactivity (hedged -- not a proven comparative-tractability ranking) | **USP34** | corrected, precise basis in Section 13A -- not "the only one with a structure"; proves Cys1903 reactivity to a macromolecular probe, does not by itself prove USP34 is more small-molecule-tractable than TLK2's canonical kinase pocket |
| Best-supported multimodal target | **Not a clean tie -- not clearly decidable** | USP34 (GSE118713, an established resistance model) and VEZF1 (GSE240112, a smaller, unpaired, source-confounded recurrence association) each have exactly one significant non-acute RNA dataset, but the two datasets differ substantially in kind and rigor and should not be read as interchangeable (Section 9, 18) |
| Overall universal winner | **NO** | no gene leads on every axis; see Section 24 for the full non-dominated role breakdown |

This table is a plain restatement of Sections 8-24 above, provided in this
compact form because it was explicitly requested as a standalone
checklist for the science-freeze process.

---

## Reproducibility

**BEFORE** (this session's already-installed but under-documented
environment): full suite 1254 passed, 1 skipped — but `environment.yml`
was missing 5 real, imported, test-covered packages (h5py, gseapy,
lifelines, networkx, requests) and had `pandas` completely unpinned (a
verified, not merely suspected, reproducibility gap — confirmed via a
full AST-level import audit of every `.py` file in `src/`, `tests/`,
`scripts/`, cross-checked against the declared dependency list).
Rebuilding a genuinely fresh conda environment from scratch was judged
impractical within this session's time budget and was not performed;
the static dependency audit is the verification method used instead, and
is exact (not an estimate).

**Fix:** `environment.yml` now declares `gseapy`, `h5py`, `lifelines`,
`networkx`, `requests`, and pins `pandas>=2,<3` with an inline rationale
comment. No test was weakened to make anything pass.

**AFTER:** `tests/test_post_audit_sensitivity.py` — 41 passed (as of the
2026-08-15 science-freeze pass; grew from 39 to 41 during that pass after
two additional regression tests were added for a second GSE240112
documentation occurrence found by the science-freeze Codex review). Full
project suite — see `SCIENCE_FREEZE_REPORT.md` and the final response for
the exact final count.
