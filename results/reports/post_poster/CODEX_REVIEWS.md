---
title: Independent review log (Codex)
status: post_freeze_exploratory
branch: post-poster-strengthening
---

# Codex review log

Three checkpoints were planned. Checkpoints 1 and 2 ran with Codex; checkpoint 3
was a structural verification carried out directly. What follows records each
review and exactly what changed in response — including two occasions where
Codex was right and I was wrong.

---

## Checkpoint 1 — plan review, before implementation

**Ran.** Thread `01a0120f-a2e6-72b0-a292-0f0c898cbf40`, read-only sandbox. Codex
was given the full plan for all three tasks and asked seven specific questions.

### Q1 — is the edgeR Wald standard error defensible?

> "Defensible as a secondary, explicitly model-based approximation, but 'exact'
> is too strong. […] Do not derive `df.total` generically as
> `df.residual + df.prior`: edgeR caps total df and robust QL fitting can make
> `df.prior` gene-specific. Use the exact df employed for each gene/test.
> Verify numerically that `t_W² ≈ F`."

**Adopted in full.** `df.total` is now taken straight from the `glmQLFTest`
object rather than re-derived. Reading `glmQLFTest`'s source also revealed that
edgeR 4.4 renamed `var.post` → `s2.post` and `df.residual.zeros` →
`df.residual.adj`; the code accepts either rather than silently producing an
unscaled SE. The `t_W²/F` check is implemented, exported as
`se_variant_diagnostic.tsv`, and test-enforced. The word "exact" was removed
from the module docstrings.

### Q2 — is Stouffer the right primary?

> "Making effect-size pooling secondary is justified. Weighted Stouffer is not
> automatically a better primary analysis. […] Neither √df nor √n is principled
> here. […] Prefer the edgeR-only effect-size synthesis as the least
> heterogeneous quantitative analysis. […] They share an analysis engine but not
> necessarily an estimand — call them 'more comparable', not homogeneous."

**Adopted.** The primary result is the `edger_only` effect-size synthesis, not
Stouffer. Stouffer is carried in every table as a descriptive sensitivity under
two weightings, both stated as unprincipled. The report says "more comparable",
not "homogeneous", and names the estimand difference between the two edgeR
designs.

### Q3 — is the √3 variance inflation acceptable?

> "Report A1 versus A2. A3 may be shown only as a heuristic stress test, not as
> a valid corrected analysis. […] Variance inflation cannot create the missing
> experimental units."

**Adopted.** `all3_de3` is labelled a heuristic stress test in
`arm_definitions.tsv`, in the module, and in the report, with Codex's reasoning
reproduced.

### Q4 — heterogeneity statistics at k = 2–3

> "Report DL, Q, I² and τ² because requested, but label all heterogeneity
> estimates essentially non-identifiable at k = 2–3. Add REML or Paule–Mandel.
> Hartung–Knapp does not rescue k = 2 and can behave pathologically."

**Adopted.** Paule–Mandel added alongside DerSimonian–Laird; Hartung–Knapp not
applied, with the reason stated; the non-identifiability caveat is in the report.

### Q5 — is the BH-implied alpha a legitimate power alpha?

> "Not a legitimate prospective '80% power at FDR 0.05'. It is outcome-dependent
> and circular. […] You may report the realised BH cutoff descriptively, but call
> it 'MDE at the observed rejection threshold'."

**Adopted.** The headline alpha is now prespecified `nominal_0.05`, with
Bonferroni 0.05/13 for the candidate family. The BH threshold is retained but
explicitly marked descriptive-only, and a test asserts the phrase "not a
prospective design alpha" survives in the alpha definitions. Codex's further
point that workflow-level power needs simulation is stated as a limitation.

### Q6 — does MDE escape the observed-power fallacy?

> "'Informative negative' is too strong. […] Prefer
> `observed_abs_log2fc_below_mde80` and
> `mde80_exceeds_external_reference_effect`. […] For evidence of negligible
> effect, prespecify a smallest effect size of interest and use equivalence
> testing."

**Adopted with one deliberate deviation.** Both suggested boolean columns exist.
The brief also required a single interpretation column, so one is provided — but
worded to Codex's constraint: no value asserts a genuine negative, and a test
enumerates every possible label and fails if any contains `genuine_negative`,
`informative_negative`, `no_effect`, `absence` or `proven`. The
prespecified-SESOI / equivalence-testing point is recorded as the analysis that
would be needed and is not available.

### Q7 — does the Task 3 feasibility gate fail?

> "The feasibility gate fails on the evidence presented, and stopping is correct.
> […] Document exact access dates, URLs/accessions, commit hash and archive
> contents; barcode-normalization rules searched; whether notebook-embedded
> objects, serialized R objects, Git LFS, releases and repository history were
> checked."

**Adopted.** The probe was extended to cover releases, tags, `.gitattributes`
(LFS), and the full commit history's file-extension set. All are recorded in
`feasibility_probe.json` and test-asserted.

### Additional traps Codex raised

| Trap | Response |
|---|---|
| Confounding must stay prominent, not just be a joined column | `group_confounded_with_biobank` is on every row **and** the consequence is stated in both reports |
| Do not assign an MDE or null verdict to filtered/missing genes | `not_tested` rows carry no MDE and a reason; test-enforced |
| Confirm candidate counts and name every missing symbol | Done — and this uncovered a correction (below) |
| Harmonise coefficient orientation before signed tests | All three contrasts are resistant/recurrent vs sensitive/primary; stated in the report |
| Gene-symbol duplication before intersections | Deterministic collapse rule (highest `avg_log_cpm`, tie-broken on gene id), with rows-in/rows-out logged: 81, 310 and 0 rows collapsed |
| Separate candidate-family from transcriptome-wide multiplicity | Separate scopes with separate BH corrections |

### A correction I made to my own checkpoint-1 briefing

I told Codex that KDM1A was missing from GSE118713 because of the expression
filter. **That was wrong.** KDM1A (with RCOR1) was withheld from the published
table by the preregistered blinding, retired 2026-08-10 per `CLAUDE.md` and the
`PREANALYSIS.md` amendments log. Both genes were always fitted. The analysis now
reads the unredacted frozen table (14,838 genes), so 12 of 13 candidates are
present in all three datasets and only USP17L29 is absent everywhere. This was
caught by following Codex's "name every dataset-specific missing symbol"
instruction.

---

## Checkpoint 2 — implementation and output review

**Ran**, after a false start. Codex's first attempt failed on every command with
`bwrap: pivot_root: Invalid argument` and it correctly refused to review from
prose alone:

> "Since you explicitly requested an evidence-based review, I won't pretend to
> validate the implementation or headline numbers from the text alone."

Rerun in a working sandbox (thread `01a01228-d337-7d31-8f52-df139c12ecc1`) with
read access to every module, table and report. It confirmed all 68 tests passing
at the time and then rejected my central inferential choice.

### The finding that changed the analysis: `f_calibrated` must not be the primary SE

> "`f_calibrated` is not a defensible inverse-variance SE. […] It is calculated
> from the same observed coefficient and test statistic being summarized. Its
> inverse-variance weight is consequently outcome-dependent. It is not a
> curvature/information estimate of `Var(beta_hat)`. It becomes undefined at
> exactly zero effect/F, potentially excluding precisely null-looking rows. […]
> The cancellation in `t_wald²/F` is sound and supports the conclusion that the
> Wald reconstruction is not simply missing the `s2.post` factor. It does not
> establish that `|beta|/sqrt(F)` is a sampling SE."

**Accepted; this was my error.** I had over-corrected: having shown the Wald
statistic runs below the quasi-F, I promoted the statistic-derived quantity to
primary, trading a calibration mismatch for a circularity. `wald` is now the
primary SE for all pooling, `f_calibrated` is labelled a diagnostic on the
statistic-equivalent scale, and a test pins both the choice and the docstring
reasoning.

Three consequences followed:

1. **A committed inconsistency Codex caught directly.** The power report quoted
   `f_calibrated` medians (0.52 / 0.51 / 0.95) while
   `dataset_sensitivity_summary.tsv` was built from `wald` (0.627 / 0.514 /
   1.007). Reverting the primary fixed the mismatch at the root; the report now
   carries the `wald` numbers, which are the ones in the table.
2. **The MDE was partly outcome-derived under `f_calibrated`**, contradicting
   this module's own stated design-property claim — and the test advertised as
   preventing that only checked the `wald` rows. `compute_power_table()` now
   *raises* on any non-information-based SE, and a test asserts the refusal.
3. **`k = 2` inference is anti-conservative.** Codex: "The code estimates DL τ²
   and then uses a normal critical value and plug-in SE. This ignores
   uncertainty in τ² […] USP34's p=0.0487 should be described as a descriptive
   plug-in result, not merely as an SE-variant flip." Both reports now say so
   explicitly, and a test asserts the module docstring keeps the caveat.

### Hardcoded values (Q2)

| Raised | Response |
|---|---|
| `N_CANDIDATES = 13` duplicated in the power module | now `len(thirteen_candidates())`, derived from config |
| Power module hardcoded its own SE variants and `wald` for the summary | now imports `PRIMARY_SE_VARIANT` from the meta module |
| `n_samples=9` for GSE118713 conflates the three-group fit with the 3-vs-3 contrast | split into `n_samples_in_model` / `n_samples_in_contrast`, test-asserted |
| Probe duplicates the pinned commit literal | left as a pinned constant; the test reads the commit from `docs/gse245601_PREANALYSIS.md` and asserts the probe matches it, so drift fails loudly |

### Drop accounting (Q3)

> "The claim 'every drop is logged' is too strong. […] no candidate silently
> disappears, but not every upstream row loss has a persistent recorded count."

**Correct.** Fixed by persisting what had only been log lines:

- `gene_symbol_collapse_audit.tsv` records rows-in, rows-lost-to-missing-symbol,
  rows-lost-to-duplicate-collapse and rows-out per dataset, with the collapse
  rule stated; a test re-derives every figure from the source table.
- `refit_validation.tsv` now carries `n_genes_before_filter`,
  `n_genes_removed_by_filter`, `n_genes_after_filter` and
  `n_genes_joined_to_frozen` (60,619→27,418; 27,161→18,428; 60,155→14,838).
- The `not_tested_reason` column no longer conflates "expression filter" with
  "annotation". **Acting on this found a factual error in my own report**:
  USP17L29 is present *before* filtering in GSE111151 and GSE118713 — it was
  expression-filtered, not unannotated — and is genuinely absent only from
  GSE240112's annotation. Both reports were corrected.

### Headline claims (Q4)

Codex verified claims 1–6 numerically and flagged three overstatements, all
corrected:

- "best FDR anywhere is 0.10" — the value is `0.100257`, so the report now says
  "smallest candidate-family FDR anywhere is 0.100".
- "the count treats SUPT4H1 worst" — its vote rank of 3 is a nine-way tie, so
  the report now says the count cannot separate it from eight other genes.
- "no data file has ever existed" in the authors' repository — the probe
  classified by extension and globbed only root-level notebooks. The probe now
  uses `rglob` and additionally records every path ever committed whose
  extension is not code-or-image; the report's wording was softened accordingly.

Claim 7 was "partially unsupported as committed-output reporting" for exactly
the report/table mismatch above, now resolved.

### Test criticism (Q5)

> "The MDE inverse test […] is self-consistency, not an independent standard."

**Accepted.** Added `test_power_matches_an_independent_monte_carlo_estimate`,
which simulates `(Z + sqrt(ncp))² / (chi2_df/df)` at finite df and compares the
empirical rejection rate to `power_at` — an independent standard rather than a
second call to the same distribution. Codex's point that the annotation tests
pin the committed JSON rather than reprobe the sources is correct and is
accepted as-is: they validate transformation and integrity, not continuing
external truth, and the probe script exists to be rerun.

### Codex's own conclusion

> "Bottom line: refit fidelity and most table arithmetic are strong. The central
> inferential choice is not: Wald should remain primary for effect-size
> meta-analysis, and the power report/summary must be reconciled before
> release."

Both were done. 77 tests pass after the changes.

---

## Checkpoint 3 — structural constraints, before the final commit

**Ran** (same thread as checkpoint 2). Codex was asked to verify seven
structural constraints directly and report each PASS/FAIL with the command
output it relied on.

| # | Constraint | Verdict |
|---|---|---|
| 1 | No read-only path changed on this branch (`poster/`, root `README.md`, `PREANALYSIS.md`, `docs/*PREANALYSIS*.md`, `results/tables/evidence_freeze/`) | **PASS** — `git diff --quiet 477c992 -- [protected]` exit 0 |
| 2 | Freeze anchors unchanged | **PASS** — commit `9a1b7777…`, shortlist SHA-256 `b6990d7e…` |
| 3 | `git diff science-freeze-2026-08-15 -- results/tables/evidence_freeze/` empty | **PASS** — exit 0 |
| 4 | Every new output labelled post-freeze exploratory | **PASS** — 13 TSV, 9 TSV.GZ, the probe JSON key, and all five report frontmatters |
| 5 | No new result referenced from any poster file or frozen doc | **PASS** — no matches for `post_poster\|meta_analysis\|MDE\|minimum detectable\|f_calibrated\|pooled_effect` |
| 6 | All new work confined to permitted paths | **FAIL** — see below |
| 7 | Nothing committed, no history rewritten, HEAD still `477c992` | **PASS** |

On constraint 1, Codex correctly separated the branch-local diff from the
pre-existing one: `git diff` against the *tag* shows `README.md` and `poster/`
differing, but that difference exists at the branch base `477c992` and was
created by the seven post-freeze audit commits already on `main`. Against
`477c992` there is no branch-local difference at all.

### Constraint 6 — the one failure, and its resolution

Codex found eight modified tracked files outside the permitted paths:

```
results/figures/poster_crispr_discovery_v1/CRISPR_discovery_main.{pdf,svg}
results/figures/poster_depmap_v1/DEPMAP_v1.{pdf,svg}
results/figures/poster_depmap_v2/DEPMAP_v2.{pdf,svg}
results/figures/poster_druggability_v1/DRUGGABILITY_v1.{pdf,svg}
```

These are **not edits**. They are re-renders produced by running the existing
poster figure tests: PDF and SVG output is not byte-reproducible (matplotlib
embeds per-run element ids, and PDFs a creation date), a fact already documented
by commit `477c992` "Document that figure tests dirty non-reproducible
artifacts". They were reverted with `git checkout --` before committing, and no
figure content changed. Codex confirmed that all the checkpoint work itself is
confined to the permitted paths, and that the probe's placement under `scripts/`
rather than `src/` "is appropriate because it performs network calls".

### Two report errors Codex caught, both fixed

1. **An unlabelled arm-sensitivity table mixed SE variants.** For USP34 it
   showed `0.049 / 0.0084 / 0.012` — the first from the *diagnostic* variant and
   the other two from `wald`, a leftover from the pre-checkpoint-2 primary. The
   table now declares `wald` throughout and reads `0.078 / 0.0084 / 0.012`.
2. **The report still said USP17L29 was "never expression-filtered out".** The
   power table shows it was filtered from GSE111151 and GSE118713 and never
   annotated in GSE240112. Corrected in both reports.

Codex also confirmed each checkpoint-2 fix had landed: `PRIMARY_SE_VARIANT =
"wald"`, the power module's refusal of any other SE, the reconciled medians
(0.626918 / 0.513566 / 1.006868), the derived `N_CANDIDATES`, the split sample
counts, both accounting tables, the distinguished missing-gene reasons, the
`rglob` probe, `pooled_z_rank`, and the Monte-Carlo power test — with
`77 passed`.
