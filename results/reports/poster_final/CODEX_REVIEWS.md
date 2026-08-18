---
title: Independent Codex review log — final poster deliverable
status: post_freeze_exploratory
analysis_date: 2026-08-18
branch: poster-final (unmerged)
---

# Codex review log

Three checkpoints, as required by the brief. Each records what was asked,
what came back, and what I actually changed in response. Where I did not
follow the reviewer, the reason is stated.

---

## Checkpoint 1 — before implementation (figure plan + narrative)

**Asked:** whether the 4-beat narrative is supported by the evidence; whether
any planned figure implies more than the data shows; whether the 7-figure set
has a gap; whether presenting four genes drawn from two different selection
rules as one candidate set is honest; and which claim a hostile specialist
attacks first. Codex was given the frozen numbers (screen effect/FDR for all
13 hits, the 2-of-16 corroboration result, the meta-analysis and power
results, the Hallmark enrichment scores, the DepMap counts, the STRING path
structure, and the structural audit) and ran read-only against the repository.

### What came back, and what I did

**1. Beat 4 as briefed is false.** The brief's framing — "only one candidate
has a chemically reachable site" — is contradicted by KDM1A's inhibitor-bound
co-crystals and nine clinical-stage LSD1 inhibitors. Codex further judged
"USP34 leads on novelty of a reachable site" *still* too strong, because
7W3U demonstrates catalytic accessibility to a ubiquitin activity-based probe
in an isolated domain (~12% of the protein), not druggability.

*Action:* beat 4 rewritten as a **maturity spectrum**, not a winner: KDM1A
clinically advanced; TLK2 structurally resolved kinase site, no validated
inhibitor; USP34 early catalytic-site probe evidence, no drug; VEZF1 no
experimental structure. The probe is never called an inhibitor or a drug
anywhere.

*Where I did not follow it:* the brief's Part 1 lists as a **required**
disclosure that "the lead candidate ranks twelfth of thirteen by effect size
in the screen and leads on tractability". I keep that sentence, because the
brief requires it and the author reaffirmed the brief, but I have tightened
"leads on tractability" to name the exact and only sense in which it is true
(a crystallographically demonstrated covalently addressable catalytic
cysteine in a target that is not already mature), and paired it in the same
breath with KDM1A's greater maturity. The tension is logged in
`SCIENTIFIC_REVIEW.md` rather than hidden.

**2. Six specific over-implications**, all accepted:

| Figure | Over-implication | Fix applied |
|---|---|---|
| F1 | "analysis plan predates results" would imply the whole poster was prospectively specified | The diagram now marks *which* parts are pre-registered and *which* are post-freeze (candidate reinterpretation, meta-analysis, power, network) |
| F3 | "replication" overstates comparability of four heterogeneous designs | Word "replicate" removed at gene level; "corroboration" used. Multiplicity family (genome-wide Benjamini–Hochberg within each dataset) now stated on the figure |
| F4 | Hallmark EMT is a transcriptional gene-set score, not a measured "mobility programme"; "dissociation" invites a causal reading | Row label describes the programme, formal gene-set name retained; caption states enrichment only, and that the acute dataset differs in more than duration |
| F5 | A STRING route is not a communication channel; "single DNMT1 bridge" is imprecise | Panel states "four shortest paths, each through DNMT1 and one ubiquitin-encoding node"; no arrows, no mechanism claimed |
| F6 | "independent axes" is not demonstrated by different counts | Title now says the two measurements "need not coincide" |
| F7 | Levels of structural evidence are not interchangeable | Evidence track separates inhibitor-bound precedent / ATP-site availability / activity-based probe / homology-model screening / no experimental structure |

**3. A real gap: the meta-analysis and power results were invisible.**
Without them, beat 2 sounds more decisive than the evidence permits. Codex
also cautioned that power explains only GSE111151's null, not every null.

*Action:* F3 gained a fourth row — a pooled-evidence strip carrying the
random-effects result (no candidate reaches pooled FDR 0.05; smallest 0.100)
and the GSE111151 minimum-detectable-effect result, with the scope limit
stated. Codex's recommendation to show the published USP34 counter-evidence
(PMID 28499884) is taken: it appears on F7 and in the poster text.

**4. Two selection rules must be disclosed, not blended.** Codex was
explicit that the four cannot be presented as one pre-specified panel.

*Action:* both F1 and F2 now mark provenance class per gene — USP34/VEZF1
from the frozen multimodal rule, KDM1A/TLK2 added post-audit as the two
strongest CRISPR hits with zero qualifying RNA corroboration — and record
that EML5/CITED2 were displaced from the original frozen shortlist. The
words "validated", "preregistered panel" and "independently validated" are
never applied to the four.

**5. The claim a hostile reviewer attacks first:** "genes do not replicate
but pathways do", because the datasets are not equivalent replications, EMT
in fact *reverses* in the acute dataset, and broad correlated gene sets can
be significant when candidate effects are unstable.

*Action:* the poster's central claim is now the calibrated form Codex
proposed — sparse candidate-level corroboration with non-significant pooled
evidence, against directionally consistent oestrogen-response and E2F
suppression across four heterogeneous datasets, with EMT differing between
chronic/recurrent and acute settings. "Pathways replicate" as a bare slogan
does not appear.

---

## Checkpoint 2 — after the figures and text existed

**Asked:** to verify every printed number against its source table; to quote
any sentence violating or approaching the prohibited claims; to confirm the
six required disclosures; to judge whether the tractability sentence is now
adequately calibrated given the brief requires it; to name the three sentences
a hostile specialist attacks first; and to say whether a well-informed sceptic
reading only the poster would get an accurate picture.

**A note on how this ran.** The first attempt failed: Codex's read-only
sandbox aborted at startup with `bwrap: pivot_root: Invalid argument`, and
Codex correctly refused to review from prose rather than guess at file
contents — *"No files were accessed or modified, so I cannot responsibly
verify the numbers or claims."* That refusal is the right behaviour and is
recorded here rather than hidden. The checkpoint was re-run with full
filesystem access and completed.

### What came back, and what I did

**1. Numbers: no substantive mismatch.** All fifteen checked values verified,
including 13 of 19,103, USP34 twelfth of thirteen, 2 of 16, the two
significant FDRs, the pooled 0.100, TLK2 9 of 11, the network's 47/147/3 and
four routes, and USP34's ~12% catalytic domain.

*One qualification, accepted:* "could only have detected a 1.5-fold change at
80% power" describes the **median** per-gene threshold, not a universal floor;
the per-gene range is 1.34- to 1.74-fold. *Action:* reworded to "had a median
80%-power detection threshold of about 1.5-fold (per-gene range 1.34- to
1.74-fold)".

**2. Prohibited claims: none found as affirmative claims.** Codex confirmed
the poster explicitly states the negation of most of them. It listed six
near-misses; I acted on five.

| Near-miss | Action |
|---|---|
| "external **validation** datasets" | changed to "corroboration datasets" — these datasets do not validate the genes |
| "works against the **benefit proposed here**" | changed to "the hypothesised sensitisation benefit" |
| reference 5 annotated "counter-evidence to the **mechanism** proposed here" | changed to "**hypothesis**" — the body text says the observation is not a mechanism, so the annotation contradicted it |
| "**chemically reachable** catalytic site" | changed to "**experimentally addressable** catalytic cysteine" |
| "**pre-registered** gate" | changed to "**pre-specified**" everywhere (text, Figure 1, Figure 2), plus a new limitations line stating the plan was dated before results but not lodged in a public registry |
| "KDM1A is the mature target" | **kept** — it is true, refers to pre-existing pharmacology, and the surrounding text makes that clear |

The "pre-registered" correction is the most substantive of these. It is a term
of art, the project does not have a public preregistration, and using it
loosely would have been the first thing a methodologist challenged.

**3. Required disclosures: all six present**, with the power statement the
only one Codex would tighten — done, as above.

**4. Tractability wording.** Codex judged the surrounding wording now
"substantially well calibrated" but flagged "chemically reachable" as still
able to sound like demonstrated small-molecule tractability, and proposed
"experimentally addressable". *Action:* adopted. Codex confirmed the required
final sentence "can stay exactly as written".

**5. The three most attackable sentences** — the categorical "nulls are
uninformative, not negative"; the tractability claim; and "pre-registered".
*Action:* the first was softened to "weak evidence, not strong negatives", the
third was corrected outright, and the second is recorded as a standing
disagreement in `SCIENTIFIC_REVIEW.md` §7 for the supervisor to arbitrate,
since the brief requires it.

**6. Fairness: yes**, with the wording reservations above. Codex's assessment
was that a sceptic "would receive an unusually accurate overall picture", and
that the poster's evidential conclusion — a hypothesis worth testing, not a
supported therapeutic target — is fair.

