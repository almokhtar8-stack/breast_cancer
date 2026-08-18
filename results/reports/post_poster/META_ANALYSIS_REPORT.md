---
title: Random-effects meta-analysis of the resistance effect
status: post_freeze_exploratory
analysis_date: 2026-08-18
branch: post-poster-strengthening
scope: GSE118713, GSE111151, GSE240112 (GSE245601 excluded by design)
---

# Replacing vote counting with a pooled effect

**post_freeze_exploratory.** This analysis describes the existing evidence
differently. It produces **no new shortlist, no new gene order and no new
lead**. Where it disagrees with the frozen ranking, the disagreement is reported
below as a discrepancy for a human to decide on; nothing has been acted on.

## The problem being addressed

The frozen evidence table ranks partly on `resistance_fdr05_count` (column 36 of
`results/tables/evidence_freeze/final_candidate_evidence.tsv`) — the number of
datasets in which a gene cleared FDR 0.05. Counting significant datasets throws
away effect size and direction, and treats FDR 0.049 and 0.051 as categorically
different. A gene significant once and near-significant twice scores below a
gene significant once and null twice.

## What was done

A random-effects meta-analysis of the per-gene resistance effect, at two scopes
(the 13 sensitising candidates, and genome-wide), reported for three arms and
two standard-error definitions.

`GSE245601 is excluded.` It measures acute 12 h response, and
`docs/THERAPEUTIC_SHORTLIST_FREEZE.md` gives the reasoning for why it is never
summed into resistance evidence. That decision is not revisited here.

### Provenance and reproduction

```
Rscript src/post_poster_de_refit.R          # needs the sc245601 R env (see below)
python -m src.post_poster_meta_analysis
pytest tests/test_post_poster_meta_analysis.py
```

`src/post_poster_de_refit.R` refits each frozen model with identical parameters
and **verifies gene-by-gene that it reproduces the committed frozen table**,
aborting if it does not:

| Dataset | Engine | genes in → after filter | max abs Δlog2FC | max abs Δp |
|---|---|---|---|---|
| GSE111151 | edgeR glmQLF | 60,619 → 27,418 | 5.3e-15 | 6.7e-16 |
| GSE240112 | edgeR glmQLF | 27,161 → 18,428 | 5.0e-14 | 6.7e-16 |
| GSE118713 | limma eBayes | 60,155 → 14,838 | 5.1e-15 | 5.6e-16 (Δse 4.9e-15) |

Those row counts are persisted in `refit_validation.tsv`, not merely logged, and
every one of the retained genes joins the frozen table one-to-one.

Environment note: the frozen R analyses run in the `sc245601` micromamba
environment, not the `bc` environment that `environment.yml` describes — `bc`
lacks `r-yaml` and `r-data.table`, so neither the existing frozen R scripts nor
this refit can run there. That is a **pre-existing** packaging gap, recorded
here rather than fixed, because changing `environment.yml` is outside this
branch's scope.

## Three decisions that had to be made honestly

### 1. The three datasets do not estimate the same quantity

GSE118713 is limma on `log2(TPM + 1)`; GSE111151 and GSE240112 are edgeR
`glmQLFit`/`glmQLFTest` on counts. Their log2 fold changes are not automatically
comparable: the `+1` offset attenuates fold changes in an expression-level
dependent way, and TPM normalises for transcript length and composition
differently from TMM.

**Judgement: pooling across that boundary is not defensible as the primary
numeric result.** The primary arm is therefore `edger_only` — same engine, same
negative-binomial log2 fold-change estimand. Even that is only *more* comparable,
not homogeneous: GSE111151 estimates a resistance-status coefficient adjusted for
cell line, GSE240112 an unpaired RT−PT contrast. The `all3` arm is reported as a
sensitivity with the scale caveat attached, not suppressed.

A signed weighted Stouffer combination of the frozen p-values is also carried in
every output (`stouffer_z_equal`, `stouffer_p_equal`, `stouffer_z_units`,
`stouffer_p_units`) as the assumption-light alternative. It is reported as a
descriptive sensitivity, not promoted to primary: it answers a different
question (directional evidence, not a common effect size), and combining
p-values does not repair a p-value from a pseudoreplicated or fully confounded
design.

### 2. edgeR's quasi-F is not a Wald test, so there is no exact SE

`glmQLFTest` computes `F = LR / s2.post` (verified: maximum absolute difference
**0** across all genes), a deviance-based statistic. No standard error inverts
it exactly. Rather than fabricate one, **two exactly-defined quantities are
carried through every table**:

- **`wald`** — the model-based standard error from the fit's own information
  matrix, `Var(β) = s2.post · (X'WX)⁻¹` with `W = μ/(1+φμ)`, converted to log2.
- **`f_calibrated`** — `|log2FC| / √F`, which reproduces the QL test statistic
  by construction.

Because `s2.post` cancels from the ratio `t_wald²/F`, that ratio is a pure
Wald-versus-likelihood-ratio comparison, and it is measured and reported in
`se_variant_diagnostic.tsv`:

| Dataset | median `t_wald²/F` | IQR |
|---|---|---|
| GSE118713 (limma) | 1.000 | identical by construction |
| GSE111151 | 0.724 | 0.681 – 0.792 |
| GSE240112 | 0.896 | 0.885 – 0.918 |

The reconstruction is *correct*: on simulated well-behaved counts it returns
0.9993, dropping the `s2.post` scaling degrades it to 0.965, and among the
candidate dispersions (`fit$dispersion`, trended, tagwise, common) the one the
fit actually used gives the ratio closest to 1. The residual gap is the familiar
Wald/LR divergence for log-link count models at small n.

**`wald` is the primary SE.** That divergence is a reason to *report*
`f_calibrated`, not to pool on it. `|log2FC|/√F` is computed from the very
coefficient and statistic being summarised, so using it as an inverse-variance
weight would make the weight outcome-dependent: it is an algebraic rescaling
that reproduces `F`, not a curvature estimate of `Var(β̂)`, and it is undefined
at exactly zero effect — which would preferentially drop the most null-looking
rows. Inverse-variance pooling requires an information-based SE.
`f_calibrated` is therefore retained as a **labelled sensitivity** showing how
much the Wald/LR gap could matter, and results that move between the two are
flagged below rather than chosen between.

282 genes in GSE240112 and 1 in GSE111151 have a singular `X'WX` and therefore no
`wald` SE; they are counted in `se_variant_diagnostic.tsv` and excluded with a
recorded reason, never dropped silently. Duplicate gene symbols are collapsed by
a deterministic rule (highest `avg_log_cpm`, ties on gene id), and the
rows-in/rows-lost/rows-out counts are persisted in
`gene_symbol_collapse_audit.tsv` — 81, 310 and 0 rows respectively.

### 2b. The pooled p-values are descriptive, not calibrated

Pooling uses a plug-in random-effects estimator: τ² is estimated, then the
pooled effect is referred to a normal critical value as though τ² were known. At
k = 2–3 studies this ignores the large uncertainty in τ² and is
**anti-conservative**. Every pooled p-value and interval below should be read as
indicative. Hartung–Knapp is not applied because it behaves pathologically at
k = 2. In particular, a pooled p near 0.05 is not a finding.

### 3. Pseudoreplication and confounding are carried as data, not prose

- **GSE118713 is pseudoreplicated**: 3 replicates per arm of a single
  MCF7/TAMR lineage. Its declared `n_biological_units` is **1**, not 3. The
  required "with and without" comparison is the `all3` vs `edger_only`
  contrast. A third arm, `all3_de3`, multiplies its SE by √3 — the cluster
  design effect at ICC = 1 — and is labelled in `arm_definitions.tsv` as a
  **heuristic stress test, not a corrected analysis**: variance inflation
  cannot create experimental units that were never run.
- **GSE240112 is confounded**: group is perfectly confounded with biobank at
  n = 3 vs 3. Every output row carries
  `group_confounded_with_biobank` / `any_group_confounded_with_biobank`. No
  weighting scheme separates recurrence from biobank, so its contribution is
  not attributable to recurrence alone.
- **Heterogeneity statistics at k = 2–3 are not identifiable.** τ² (DerSimonian–
  Laird and Paule–Mandel), Q and I² are reported because they were requested,
  and must be read as descriptive. Hartung–Knapp is not applied: it behaves
  pathologically at k = 2.

## Results — the 13 candidates

Primary arm `edger_only`, primary SE `wald`. Positive pooled effect = higher
expression in the resistant/recurrent state. `pooled_z_rank` ranks on signed
evidence strength (the pooled z), which is what the vote count was standing in
for — it is not a ranking on effect magnitude.

| Gene | pooled log2FC | 95% CI | p | BH FDR | I² | vote count | vote rank | pooled z-rank | shift |
|---|---|---|---|---|---|---|---|---|---|
| USP34 | +0.227 | −0.025 – 0.479 | 0.078 | 0.782 | 0 | 1 | 1 | 1 | 0 |
| SUPT4H1 | +0.451 | −0.507 – 1.410 | 0.356 | 0.782 | 79.2 | 0 | 3 | 2 | +1 |
| TADA2B | +0.136 | −0.196 – 0.468 | 0.423 | 0.782 | 0 | 0 | 3 | 3 | 0 |
| KDM1A | +0.118 | −0.198 – 0.434 | 0.463 | 0.782 | 0 | 0 | 3 | 4 | −1 |
| VEZF1 | +0.435 | −0.925 – 1.794 | 0.531 | 0.782 | 92.5 | 1 | 1 | 5 | **−4** |
| TLK2 | +0.069 | −0.169 – 0.308 | 0.570 | 0.782 | 0 | 0 | 3 | 6 | −3 |
| ICK | +0.013 | −0.477 – 0.503 | 0.958 | 0.958 | 44.4 | 0 | 3 | 7 | −4 |
| TSR3 | −0.090 | −0.579 – 0.398 | 0.717 | 0.782 | 32.6 | *none* | – | 8 | – |
| CTDNEP1 | −0.067 | −0.409 – 0.276 | 0.703 | 0.782 | 0 | *none* | – | 9 | – |
| PET117 | −0.088 | −0.484 – 0.307 | 0.662 | 0.782 | 0 | *none* | – | 10 | – |
| EIF4ENIF1 | −0.075 | −0.331 – 0.180 | 0.563 | 0.782 | 0 | 0 | 3 | 11 | **−8** |
| HMGB1 | −0.134 | −0.499 – 0.231 | 0.471 | 0.782 | 0 | 0 | 3 | 12 | **−9** |
| USP17L29 | — | — | — | — | — | 0 | 3 | — | not tested anywhere |

### The headline

**No candidate reaches BH FDR 0.05 for the pooled resistance effect in any arm,
under either SE definition.** The smallest candidate-family FDR anywhere is
0.100 (USP34, `all3`, `wald`) — and given §2b, even that is a descriptive
plug-in value. Pooling does not rescue a resistance association for any of these
genes; it makes explicit that the resistance evidence for all 13 is weak.

### Named discrepancies with the vote count

1. **VEZF1 is the clearest vote-counting failure.** It shares the top vote-count
   rank (count = 1) but falls four places to pooled z-rank 5. Its two
   contributing datasets point in *opposite directions*, with I² = 92.5 and
   Cochran Q p = 0.00088 in the `all3` arm. Its single "vote" is one dataset's
   significance in the face of active disagreement — exactly what a count
   cannot see.
2. **SUPT4H1 rises** from vote-count rank 3 (count = 0, i.e. never significant
   anywhere) to pooled z-rank 2, on the largest pooled effect of any candidate
   (+0.451) with consistent direction. Its vote-count rank of 3 is a nine-way
   tie, so "the count ranks it worst" only means the count cannot separate it
   from eight other genes — which is the point. Its I² is 79.2, so this is not a
   clean win either: it is a gene worth looking at, not a gene to promote.
3. **HMGB1 (−9) and EIF4ENIF1 (−8) fall furthest**, because their pooled effects
   point the *wrong way* for the resistance hypothesis. A count of zero
   significant datasets is direction-blind and hid that.
4. **USP34 is unmoved at rank 1** under both schemes — the one place where vote
   counting and pooling agree, and the candidate with the strongest signed
   pooled evidence.
5. **USP34 is not significant, and is not robustly near-significant either.**
   In the primary arm its pooled p is 0.078 (`wald`); under the `f_calibrated`
   sensitivity it is 0.049. It straddles 0.05, and §2b means neither value is a
   calibrated p at k = 2. USP34 must not be described as "significant in the
   pooled analysis".
6. **Three candidates have no vote count at all.** TSR3, CTDNEP1 and PET117 have
   no row in `final_candidate_evidence.tsv`, so `resistance_fdr05_count` does not
   exist for them. They are pooled here for the first time. All three sit in the
   bottom half.
7. **USP17L29 is absent from all three fitted tables, for two different
   reasons.** It was present before filtering and removed by the expression
   filter in GSE111151 and GSE118713, and is not in GSE240112's gene annotation
   at all (see `not_tested_reason` in the power tables). It has a vote count of
   0, which reads as "tested and null" but means "filtered out for low
   expression" in two datasets and "never measured" in the third. This is a
   reporting defect in the count itself.

### Arm sensitivity

Primary SE (`wald`) throughout, so the three columns are comparable:

| Gene | `edger_only` p | `all3` p | `all3_de3` p |
|---|---|---|---|
| USP34 | 0.078 | 0.0084 (FDR 0.100) | 0.012 |
| SUPT4H1 | 0.356 | 0.318 | 0.315 |
| VEZF1 | 0.531 | 0.287 | 0.366 |

Adding the pseudoreplicated GSE118713 strengthens USP34 by roughly an order of
magnitude in p, and inflating its SE by √3 gives most of that back. USP34's apparent
gain in `all3` is therefore substantially a function of counting three
non-independent replicates as three units.

## Results — genome-wide

Genes measured in every dataset of the arm are pooled; the rest are returned with
`pooled_effect` NaN and an `exclusion_reason`, so nothing disappears.

| Arm | symbols considered | pooled | BH FDR < 0.05 | τ²=0 | Q p<0.05 |
|---|---|---|---|---|---|
| `edger_only` | 31,688 | 13,846 | 452 | 36.3% | 37.3% |
| `all3` | 32,178 | 11,743 | 1,141 | 25.1% | 45.8% |

Roughly 40% of genes show significant between-dataset heterogeneity, which is
the genome-wide version of the VEZF1 finding: these datasets frequently
disagree, and any summary that ignores that will overstate consistency. The
strongest pooled `all3` hit is GREB1 (−7.26), a canonical estrogen-regulated
gene — a sanity check that the pooling is oriented correctly, not a new finding.

## What this does and does not change

**Changes how the evidence should be described.** "Cleared FDR 0.05 in one
dataset" should be replaced by the pooled estimate and its interval, and by the
direction-consistency and heterogeneity statistics. On that basis the resistance
evidence for all 13 candidates is weak, and for VEZF1 it is not merely weak but
internally contradictory.

**Does not contradict the frozen conclusions.** The frozen shortlist was never
built on resistance transcriptomics alone — CRISPR effect is the primary axis —
and this analysis does not touch CRISPR. USP34, the frozen rank-1 gene, remains
rank 1 here. The candidate ordering below rank 1 does change, and that change is
reported for a human to weigh, not applied.
