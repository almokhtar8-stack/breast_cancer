# Pre-analysis plan: tamoxifen response modifiers in ER+ breast cancer

**Date written:** 2026-08-05
**Status:** locked after first analysis is run. Do not edit the sections below
once analysis begins — append dated amendments only (see bottom).

---

## 1. Question and claim boundary

Which gene knockouts change how ER-positive breast cancer cells respond to
tamoxifen, and which of those candidate targets carry the lowest predicted
risk of contributing to the joint pain, muscle pain and fatigue that drive
patients to discontinue endocrine therapy.

The CRISPR screen underlying the labels (Hany et al. 2023) was run in
**MCF7-V**, a drug-tolerant *parental* clone selected under short-term
4-OHT exposure, not a stably *acquired resistant* derivative (e.g. a TamR
line grown out under months of continuous drug pressure). These are
biologically distinct states — drug tolerance vs. drug resistance — and
results from one do not necessarily transfer to the other.

**Claim boundary:** any gene nominated by this analysis "modulates
tamoxifen response" in the drug-tolerant-persister setting. It must never be
described as "restoring tamoxifen sensitivity" or "reversing resistance" —
that would require validation in an acquired-resistance model, which this
project does not have (GSE118713 and GSE111151 are used only as external
expression features, not as outcome data; see §5).

---

## 2. Label definition

- **Source:** Hany et al. 2023, *Science Advances* 9:eadd3685, Data S1.
  DOI: [10.1126/sciadv.add3685](https://doi.org/10.1126/sciadv.add3685).
- **Contrast:** E2 + 4-OHT versus E2 alone. Oestradiol is present in both
  arms, so the contrast isolates the effect of tamoxifen with hormone
  stimulation held constant.
- **Statistic:** gene-by-treatment interaction term, estimated from raw
  sgRNA counts (not from pre-computed fold-changes or the paper's own
  gene-level summary statistic, so that guide-level filtering in §3 can be
  applied consistently before aggregation).
- **Sign convention:** **negative** = guide/gene more depleted under 4-OHT
  than under E2 alone (candidate knockout sensitizes cells to tamoxifen,
  i.e. cells need that gene to tolerate the drug). Positive = enriched
  under 4-OHT (candidate knockout confers a growth advantage under
  tamoxifen). This sign convention must be carried through every table,
  figure, and model output without renaming or flipping.

---

## 3. Exclusions

Applied before any modelling, on the raw-count matrix:

- **Treatment arms:** db-cAMP and ICI (fulvestrant) arms excluded entirely.
  Only E2 and E2+4-OHT arms are used, per the claim boundary in §1 — this
  is a tamoxifen-response screen, not a general anti-estrogen or cAMP
  screen.
- **Zero-count guides at T0:** any guide with zero raw counts in the T0
  (plasmid/day-0) reference is excluded — a zero baseline makes any
  fold-change or interaction estimate for that guide undefined/unstable.
- **Low-representation genes:** genes represented by fewer than 3
  surviving guides (after the above filters) are excluded from gene-level
  aggregation — insufficient within-gene replication to estimate a
  stable interaction term.

---

## 4. Gate 1 — is there a signal worth modelling?

Applied to the gene-by-treatment interaction results from §2–3, using
FDR < 0.1 as the significance threshold:

| Significant genes (FDR<0.1) | Decision |
|---|---|
| ≥ 30 | Enough for a classifier (discrete hit/non-hit modelling, §7) |
| 10–29 | Continuous modelling only (rank/regression on the interaction statistic; no discrete classifier — too few positives to hold out a meaningful test fold) |
| < 10 | No model. Stop and report a negative/underpowered result. Do not force a model onto a screen that didn't yield a signal. |

This gate is evaluated once, on the full filtered label set, before any
feature-based modelling is attempted. The outcome and the actual gene
count must be recorded in the analysis log at the time this gate is run.

---

## 5. Positive control

**CDK4** was already noticed during exploratory inspection of the Hany
data before this document was finalized. It is therefore **not a blind
control** — its presence in the top hits confirms the pipeline is
producing sane, biologically plausible output, but it cannot be used to
claim the method has predictive power, since it was seen in advance.

**RCOR1** and **KDM1A (LSD1)** are held out as the blind positive
controls: both are LSD1/CoREST-complex members with independent prior
evidence linking them to ER signalling and endocrine therapy response,
and neither was inspected or discussed before this plan was written. Their
behavior in the screen results is not to be checked until the modelling
pipeline (§7) is otherwise finalized, so that seeing them can't retroactively
shape feature or model choices.

> **Dated amendment (2026-08-05):** CDK4 downgraded from blind positive
> control to sanity check for the reason above. RCOR1 and KDM1A substituted
> as the blind controls. This amendment is written into the original
> document rather than appended below because it predates any analysis
> being run — see the note on amendment policy in §12.

---

## 6. Feature sources

| Source | Content | Role |
|---|---|---|
| GSE118713 | MCF7 parental / TAMR / FASR expression, TPM, n=3 per line | Baseline and resistant-state expression features |
| GSE111151 | 4 lines, 7 TamR derivatives, counts, n=1 per line | Additional resistant-derivative expression features |

**Limitations, carried forward into interpretation, not fixed by
preprocessing:**

- GSE118713 is TPM, not raw counts — it is not directly comparable in
  scale or variance structure to GSE111151 (counts) or to the Hany screen
  data. No cross-dataset count-vs-TPM harmonization is attempted; features
  from each source are kept on their native scale and, where combined,
  only via rank/percentile transforms (see §9), never raw values.
- GSE111151 has **no replicates** (n=1 per line) — no within-line variance
  estimate is possible, so any feature derived from it should be treated
  as a point estimate, not a mean, and downstream models should not be
  given false confidence in its precision (e.g. no per-gene error bars
  fabricated from a single sample).

**Leakage rule:** nothing derived from the Hany screen (labels, screen
read counts, screen-derived gene lists, or any statistic computed from
screen data) may enter the feature table. Features are restricted to
GSE118713, GSE111151, and external/public annotation resources (e.g.
pathway membership, used for the split in §7) that were not generated by
this screen. This rule is the single most important safeguard against
inflated performance estimates and must be checked explicitly at feature
table construction time, not just implicitly assumed.

---

## 7. Model

- **Modelling target:** determined by Gate 1 (§4) — classifier if ≥30
  significant genes, continuous target otherwise.
- **Split:** by **pathway**, not by gene. Genes in the same pathway/complex
  are not independently informative (e.g. correlated features, shared
  regulation, shared screen behavior), so a gene-level random split would
  leak pathway-level signal between train and test. Pathway assignment
  uses an external annotation source (not screen-derived, per §6), fixed
  before any train/test split is made.
- **Primary metric:** AUPRC (not AUROC) — appropriate given the expected
  class imbalance between hit and non-hit genes, and more informative
  than AUROC when the positive class is a small minority.
- **Mandatory baselines**, all reported alongside the model on the same
  split and metric:
  1. **Random baseline** — label-shuffled / prevalence-matched random
     classifier, to anchor the floor of the metric.
  2. **Guide-count baseline** — a model using only the number-of-guides-
     per-gene feature (already used as a filter in §3, but its residual
     predictive value on the retained gene set must be checked
     separately) to catch any screen-artifact-driven signal.
  3. **Single-feature baseline** — best individual feature from §6 used
     alone, to establish whether a multi-feature model earns its
     complexity over the strongest univariate predictor.
  4. **Pathway-membership-only baseline** — predicts hit status purely
     from pathway identity (no gene-level expression features), to
     quantify how much of any apparent model performance is actually
     just the pathway-level split in §7 leaking coarse signal.

Model performance is only interpretable relative to these four baselines,
not in absolute terms.

---

## 8. Gate 2

Applied after the model (or continuous analysis) from §7 is run: the
model/analysis must beat the strongest of the four mandatory baselines
(§7) on AUPRC on the held-out pathway split, by a margin large enough to
not plausibly be explained by split variance (exact margin and variance
estimate — e.g. via repeated pathway-grouped resampling — to be fixed and
logged when Gate 2 is actually run, since it depends on the realized
number of pathways and positives from Gate 1).

- **Pass:** proceed to candidate nomination and tolerability scoring (§9).
- **Fail — fallback:** do not force a ranked candidate list out of a model
  that didn't clear its baselines. Fall back to reporting the Gate-1
  significant gene list directly (interaction statistic + FDR, no
  model-based ranking), with the blind controls (§5) called out, and state
  explicitly that no feature-based model outperformed the baselines.

---

## 9. Tolerability layer

Candidate genes that pass Gate 2 (or, under the fallback, the raw
significant-gene list) are additionally scored for plausible contribution
to endocrine-therapy-intolerance symptoms (joint pain, muscle pain,
fatigue), using expression/association evidence external to the screen.

**Scoring rule:** every tolerability-relevant value is converted to a
**percentile rank within its own source dataset** before use or comparison
— never compared as a raw value across datasets or platforms. This is a
direct consequence of the GSE118713/GSE111151 scale mismatch noted in §6
(TPM vs. counts, replicated vs. unreplicated): raw values from different
platforms are not on a common scale and any cross-platform comparison of
raw magnitudes would be spurious. Percentile-within-dataset is the only
comparison this analysis will make or report.

---

## 10. Success criteria

- Gate 1 clears ≥10 significant genes (continuous modelling floor) — a
  hard prerequisite for the rest of the plan to be executable at all.
- The blind positive controls (RCOR1, KDM1A) are recovered among the
  significant genes and, if a model is built, rank favorably relative to
  the random and baseline models in §7.
- If a model is built (Gate 1 ≥30 genes), it clears Gate 2 against all
  four mandatory baselines.
- At least one candidate gene passing Gate 2 (or, under fallback, in the
  Gate-1 significant list) has a tolerability percentile indicating
  comparatively low predicted contribution to joint pain/muscle
  pain/fatigue, giving a concrete, prioritizable candidate for follow-up.
- Every claim in the final report stays within the "modulates tamoxifen
  response" boundary from §1; no claim of resistance reversal is made.

## 11. Failure criteria

- Gate 1 yields <10 significant genes: analysis stops and is reported as
  negative/underpowered, per §4. This is a valid, reportable outcome, not
  a failure of process.
- The blind controls (RCOR1, KDM1A) are **not** recovered as significant
  or do not show the expected direction of effect: this is reported as a
  concern about the label/pipeline's biological validity, and any
  downstream candidate list is presented with that caveat rather than
  suppressed or reinterpreted post hoc to make the controls "fit."
  (CDK4's behavior, per §5, is reported only as a sanity check, never as
  evidence of predictive validity.)
- A model is built but fails Gate 2 against any of the four baselines: no
  model-ranked candidate list is reported; fall back per §8.
- The leakage rule in §6 is found to have been violated at any point after
  results exist: those results are invalidated and rerun, not patched or
  selectively re-reported.
- Any result section is found to rely on cross-dataset comparison of raw
  (non-percentile) values from §6/§9: that section is invalid and must be
  redone using the percentile rule.

---

## 12. Amendment policy

This document is **locked once the first analysis is run**. After that
point, nothing above may be edited, deleted, or silently reworded —
corrections, scope changes, or newly discovered issues are recorded as
**dated amendments appended below this line**, each one explaining what
changed and why. The one exception is the §5 amendment above, which
predates any analysis and was folded into the relevant section for
readability; every amendment from this point forward goes at the bottom.

### Amendments log

**Amendment 2026-08-05:** corrected description of the E2 arm; it is a
hormone-treated comparator, not a vehicle control. The contrast itself is
unchanged.
