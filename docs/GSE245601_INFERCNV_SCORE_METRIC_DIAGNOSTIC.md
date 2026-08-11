# InferCNV whole-genome CNV-score metric diagnostics

**Status: audit-only, diagnostic-only. The frozen 0.01 cutoff is unchanged.
No malignancy label is changed. InferCNV/CopyKAT are not rerun (only
`run.final.infercnv_obj` is reloaded read-only, exactly as in the prior
diagnostic phase). No candidate gene is read. CopyKAT is used only as a
plotting annotation, never to tune or validate a metric.**

Scope: same 12 samples and same three tumor groups as the prior
threshold-diagnostics phase -- disagreement (Tumor_01, Tumor_04, Tumor_08),
good agreement controls (Tumor_02, Tumor_03), and Tumor_10 tracked
separately. Question: is InferCNV's whole-genome mean-squared CNV score
(`cnv_score = mean((expr-1)^2)` across every gene) an appropriate
malignancy metric when abnormalities affect only part of the genome?

All numbers are reproducible from `results/tables/gse245601_infercnv_score_metric_diagnostics.tsv`,
written by `scripts/analysis/gse245601_11_extract_cnv_score_metric_diagnostics.R`
(a read-only reload of `run.final.infercnv_obj`, cross-checked to
reproduce the frozen `cnv_score` for every one of 11,211 cells to within
1e-9 before writing anything -- confirmed, 0 mismatches), and rendered by
`src/gse245601_infercnv_score_metric_diagnostics.py`.

---

## 1. CNV signal-extent metrics (Point 1)

Per epithelial cell: the frozen `cnv_score`; fraction of genes with
`|expr-1|` exceeding three pre-declared, round deviation bands (0.05,
0.10, 0.15 -- chosen to span the dynamic range InferCNV's own output
typically shows in this project, roughly 0.8-1.2 on its heatmap color
scale; not tuned to any outcome); max/p95/p99 deviation across genes; and,
per chromosome, that chromosome's own mean `|expr-1|` for that cell, from
which the count of chromosomes exceeding each deviation band is derived.
These are diagnostic measurements only -- none of them feeds back into any
classifier.

## 2. Disagreement tumors vs. good controls

(Each cell shows the range across the tumor's two arms, Control and
Tamoxifen; medians are per-arm medians across that arm's cells.)

| Tumor | median CNV score | median max `\|dev\|` | median p99 `\|dev\|` | median extent (frac. genes >0.10) | median chromosomes >0.05 | median chromosomes >0.10 |
|---|---:|---:|---:|---:|---:|---:|
| Tumor_01 | 0.0027-0.0033 | 0.278-0.305 | 0.196-0.209 | 0.12-0.15 | 2 | 0 |
| Tumor_04 | 0.0029-0.0031 | 0.246-0.264 | 0.186-0.197 | 0.135-0.136 | 2-3 | 0 |
| Tumor_08 | 0.0060-0.0063 | 0.302-0.303 | 0.228-0.230 | 0.228-0.234 | 7-8 | 1-2 |
| Tumor_02 | 0.0047-0.0126 | 0.400-0.657 | 0.243-0.385 | 0.198-0.295 | 4-9 | 1-3 |
| Tumor_03 | 0.0105-0.0106 | 0.440-0.456 | 0.330-0.336 | 0.291-0.300 | 8 | 4 |
| Tumor_10 | 0.0019-0.0034 | 0.237-0.288 | 0.160-0.195 | 0.067-0.153 | 1-3 | 0 |

**Peak individual-gene deviation overlaps substantially but is not
identical** -- median max deviation ranges 0.25-0.30 for Tumor_01/04/08/10
and 0.40-0.66 for Tumor_02/03; the full per-cell distributions do overlap
(see `score_vs_upper_percentile_grid.png`), but p99 deviation is
systematically somewhat lower in Tumor_01/04 than in Tumor_02/03,
especially Tumor_03. Individual genes reach broadly comparable extremes,
but max/p99 are noisy, single-gene-sensitive statistics and should not be
read as "peak signal is simply equal" -- only as "not obviously smaller by
an order of magnitude."

**Chromosome burden is lower in Tumor_01/04, and Tumor_08 is clearly
intermediate, not a third homogeneous member of one group.** At the
permissive 0.05 band, Tumor_01/04 show fewer chromosomes with elevated
mean deviation (2-3) than Tumor_02/03 (4-9), though the ranges overlap
(Tumor_02_Control alone is as low as 4). At the 0.10 band, **Tumor_01 and
Tumor_04's median cell has zero chromosomes** exceeding it, versus 1-4 for
Tumor_02/03 -- but Tumor_08 sits at 1-2, and at the 0.05 band Tumor_08
(7-8) is comparable to or higher than Tumor_02_Control. The chromosome-
burden heatmaps (`chromosome_burden_heatmap_grid.png`) show Tumor_02/03
with several consistently high-mean-deviation chromosome columns across
most cells, and Tumor_01/04 with only faint, inconsistent ones; Tumor_08
is visually in between. **Important caveat on what "chromosome burden"
here can and cannot show**: it is a per-chromosome mean of *absolute*
deviation, so it does not distinguish a real, contiguous, chromosome-wide
copy-number change from mixed gains and losses averaging up, from a
strong focal sub-region elevating the whole-chromosome mean, or from
noise in a chromosome with few retained genes; InferCNV's own smoothing
also makes neighboring genes non-independent, which can itself produce
sustained-looking structure. This metric shows *higher chromosome-level
mean absolute burden* in Tumor_02/03, not proven chromosome-wide
coherence -- establishing true coherence would need signed deviations,
contiguous-run lengths, or segment-level summaries, none of which were
computed here.

**A key refinement to the prior phase's working hypothesis -- itself
appropriately hedged, not a clean resolution.** The prior document framed
this as "strong localized signal diluted by the whole-genome mean." The
score-vs-extent scatter (`score_vs_extent_grid.png`) shows that within
every single sample, `cnv_score` tracks genomic extent tightly and
similarly across all six tumors. This is expected mathematically --
`cnv_score` is approximately `extent x (typical deviation)^2` by
construction, so a tight score-extent relationship is not independent
evidence that whole-genome averaging is diagnostically *appropriate*; it
only shows that Tumor_01/04/08 are not behaving anomalously relative to
how the same formula treats Tumor_02/03. That means the data support
"these tumors are not a special case of the metric misbehaving" more
strongly than they support "the metric is fine for detecting localized
disease-relevant CNVs in general" -- the latter, broader question about
whole-genome averaging as a metric class is not resolved by this analysis
either way. What differs between groups is that these tumors' cells reach
a **lower measured genomic extent** in the first place (fewer genes,
fewer chromosomes with elevated mean deviation), and because `cnv_score`
is approximately extent-weighted by construction, a lower extent
mechanically produces a lower score, via the same relationship that
governs Tumor_02/03 too. The fixed 0.01 floor does not adapt to that
difference in extent -- but this analysis cannot say whether that lower
measured extent reflects genuinely weaker CNV burden, or a whole-genome
metric failing to register a diagnostically meaningful but focal event;
see the caveat above and Point 3 below.

## 3. InferCNV denoising/centering (Point 3)

**Not available from saved files -- confirmed by directly listing each
sample's InferCNV output directory.** Every one of the 12 samples'
`data/processed/gse245601/infercnv/<sample>/` directories contains exactly
six files (`annotations.tsv`, `gene_order_used.tsv`,
`infercnv.heatmap_thresholds.txt`, `infercnv.observation_groupings.txt`,
`infercnv.png`, `run.final.infercnv_obj`) -- no intermediate per-step
object. This is expected: `scripts/analysis/gse245601_05_infercnv_malignant.R`
explicitly calls `infercnv::run(..., resume_mode = FALSE, save_rds = FALSE)`,
so InferCNV never wrote its internal step-by-step checkpoints to disk. No
pre-denoise expression matrix exists anywhere in this project's frozen
outputs.

**Per the task's explicit instruction, no rerun was performed.** What a
rerun would require, exactly, if this comparison is wanted later: rerun
`infercnv::run()` on the same frozen inputs (same raw counts, same
reference cells, same gene order) for a small number of representative
samples (e.g. one from each of the three groups) with `denoise = FALSE`
(and/or `save_rds = TRUE` to capture intermediate steps), then compare
per-cell score/variance between that run and the frozen `denoise = TRUE`
result. This would be a new, additional InferCNV run -- not a
modification of the frozen one -- and was not undertaken here.

**Boundary on what the frozen object alone can ever establish**, even
with such a rerun: it can characterize what the classifier *received*
(the post-denoise, post-centering signal analyzed throughout this
document). It cannot, on its own, distinguish a tumor with genuinely weak
CNV burden from one whose real signal was attenuated or distorted by
InferCNV's own upstream processing -- that distinction would require the
pre-denoise comparison described above, or an orthogonal ground truth.

## 4. Tumor_10, investigated separately (Point 4)

Tumor_10 does **not** fit the Tumor_01/04/08 pattern. Reusing the
already-verified recomputation from the prior phase
(`src.gse245601_infercnv_threshold_diagnostics`), per-cell failure
categories for Tumor_10:

| Sample | fails CNV-score only | fails correlation only | fails both | passes both |
|---|---:|---:|---:|---:|
| Tumor_10_Control | 175 | 3 | 177 | 61 |
| Tumor_10_Tamoxifen | 59 | 0 | 217 | 2 |

In Tumor_01/04/08, "fails both" was a small minority and "fails
correlation only" was essentially zero -- the story there was almost
purely about CNV score. In Tumor_10, **"fails both" is comparable to or
larger than "fails CNV-score only"** in both arms, and dominates in
Tamoxifen. `tumor10_diagnostic.png` shows why: Tumor_10's cells form two
visually distinct clusters on the score/correlation plane -- a
low-correlation cluster (Kendall correlation ~0-0.2, essentially all
failing) clearly separated from a higher-correlation cluster (~0.4-0.85)
that mostly just sits below the 0.01 score line. Mean correlation is also
markedly lower in the Tamoxifen arm (0.176) than Control (0.307),
consistent with the shift toward "fails both" in that arm -- this is a
descriptive difference between the two arms, not evidence of a Tamoxifen
treatment effect (sample composition or technical differences are equally
plausible and untested here). Genomic extent (median 0.07-0.15, the
lowest of all six tumors) and chromosome burden (median 0-3) are also
low, comparable to or below Tumor_01/04. What can be said precisely: the
frozen classifier's own correlation-to-seed output shows a distinct
low-correlation component in Tumor_10 that Tumor_01/04/08 do not show
(their correlation values are almost never the blocking factor); *why*
that component exists -- biological heterogeneity, an unrepresentative
seed, technical quality, or something else -- is not established here.
Tumor_10 combines a low-extent CNV-score problem **with** this separate
correlation issue -- two compounding issues, not a single explanation --
consistent with the user's instruction not to force it into the
Tumor_01/04/08 story.

## 5. Diagnostic plots

Under `results/figures/gse245601_infercnv_score_metric_diagnostics/`:

- `score_vs_extent_grid.png` -- CNV score vs. genomic extent, one panel
  per (tumor, condition); InferCNV label = color, CopyKAT call = marker
  shape (annotation only); horizontal line at the frozen 0.01 floor.
- `score_vs_upper_percentile_grid.png` -- same layout, x = p99 deviation
  amplitude.
- `per_sample_score_and_extent_distributions.png` -- box plots of CNV
  score and genomic extent for all 12 samples, colored by tumor group.
- `chromosome_burden_heatmap_grid.png` -- per (tumor, condition): cells
  (rows, sorted by CNV score) x chromosomes 1-22 (columns), colored by
  that chromosome's own mean deviation for that cell.
- `tumor10_diagnostic.png` -- Tumor_10 only, colored by recomputed failure
  category: score vs. correlation, and extent vs. correlation.

---

## 6. Answers

**Is the whole-genome mean-squared CNV score systematically penalizing
cells with localized/partial-genome CNVs?** This analysis shows the score
tracks genomic extent tightly and *consistently* across all six tumors --
there is no evidence Tumor_01/04/08 sit in an unusual score-extent regime
relative to Tumor_02/03, i.e. no sign of a tumor-specific distortion. But
that consistency is expected mathematically (`cnv_score` is approximately
extent x deviation², so a tight score-extent relationship is close to
definitional, not an independent validation), so it cannot by itself
answer the broader question of whether whole-genome averaging is an
appropriate way to detect diagnostically meaningful *localized* CNVs in
general -- that broader question remains open. What this analysis does
establish is narrower but still practically important: **the score is
extent-weighted by construction, a fixed absolute threshold (0.01) does
not adapt to a tumor's own extent distribution, and Tumor_01/04/08 have
lower measured extent than Tumor_02/03** (see below) -- whatever the
underlying cause of that lower extent.

**Why do Tumor_01/04/08 have lower scores than Tumor_02/03 despite visible
CNV signal?** Individual-gene peak deviation (max, p99) overlaps
substantially across all six tumors but is somewhat and systematically
lower in Tumor_01/04 than Tumor_02/03 -- not simply "equal." The larger,
more consistent difference is in chromosome-level burden: far fewer
chromosomes show elevated mean absolute deviation in Tumor_01/04 (median 0
at the 0.10 band) than in Tumor_02/03 (1-4). Tumor_08 is genuinely
intermediate, not a third member of one homogeneous "disagreement" group
-- broad at a low deviation bar (comparable to Tumor_02_Control) but
narrow at a higher one, consistent with its previously-found threshold
sensitivity. One important caveat: "chromosome-level mean absolute
deviation" cannot, by itself, distinguish a real contiguous chromosome-wide
copy-number change from mixed gains/losses, a strong focal sub-region, or
noise -- it shows higher measured burden in Tumor_02/03, not proven
chromosome-wide coherence.

**Is denoising/centering contributing materially?** Cannot be determined
from saved files -- no pre-denoise intermediate exists anywhere in this
project's frozen InferCNV outputs, and none was generated for this
diagnostic per the "no rerun" instruction. This remains a genuinely open,
untested possibility, not ruled in or out -- and more fundamentally, the
frozen final object can only characterize what the classifier *received*,
not whether upstream processing attenuated a real signal before this
analysis ever saw it.

**Is Tumor_10 a separate failure mode?** Yes. Unlike Tumor_01/04/08 (where
correlation to seed is essentially never the blocking factor), Tumor_10's
frozen classifier output shows a distinct low-correlation subpopulation
and a markedly lower mean correlation in its Tamoxifen arm specifically.
Its low genomic extent compounds with this separate correlation issue.
Why the low-correlation component exists -- biology, seed
representativeness, technical quality, or something else -- is not
established here, and the Control-vs-Tamoxifen difference is descriptive,
not evidence of a treatment effect. It should not be explained by the
same single mechanism as the other three tumors.

**Keep, modify, or abandon the current metric?** The evidence does not
support silently keeping the current fixed-threshold metric unexamined --
it demonstrably assigns low scores to tumors with lower measured extent
regardless of peak amplitude, and a fixed 0.01 floor cannot adapt to that.
It also does not support concluding InferCNV's underlying signal is
unreliable or should be abandoned -- the signal itself (peak deviations,
chromosome-level differences) looks structured and directionally
consistent with biology, just not proven to be diagnostically appropriate
at the whole-genome-average level. The evidence supports **treating "is a
whole-genome average the right summary for extent-variable CNV signal"
as an open methodological question worth evaluating** -- candidate
directions could include a genome-wide burden statistic paired with a
separate focal/segmental burden measure, chromosome-arm-level summaries,
or an extent-aware rather than fixed threshold -- but no specific
replacement is selected, implemented, or recommended here, and any
candidate would need validation against an independent target, not
selection based on which one recovers the most cells in Tumor_01/04/08.

---

## Files

- `results/tables/gse245601_infercnv_score_metric_diagnostics.tsv` --
  11,211 rows (one per cell), 42 columns: identifiers, frozen labels
  (annotation only), recomputed/cross-checked `cnv_score`, deviation-band
  fractions, max/p95/p99 deviation, chromosome-count-affected columns, and
  22 per-chromosome mean-deviation columns.
- `results/figures/gse245601_infercnv_score_metric_diagnostics/` -- 5
  figures (described above).
- `scripts/analysis/gse245601_11_extract_cnv_score_metric_diagnostics.R`
  (read-only extraction, cross-checks its own output against the frozen
  `cnv_score` before writing) + `src/gse245601_infercnv_score_metric_diagnostics.py`
  + `tests/test_gse245601_infercnv_score_metric_diagnostics.py` (8 tests).

## Codex independent review

Reviewed inline (sandbox filesystem access unavailable in this
environment, as in prior audit phases): metric definitions, calculations,
and every interpretive claim in Section 6. Verdict: the quantitative
results and the R/Python cross-checks (exact `cnv_score` reproduction,
failure-category counts matching the independently-written Point-1 table)
were not challenged. Several interpretive claims were found to overreach
and were revised (reflected in the text above, not just here):

1. The original framing ("the score is NOT penalizing localized CNVs")
   overstated what a tight, consistent score-vs-extent relationship can
   show -- that relationship is close to mathematically definitional
   (`cnv_score ≈ extent x deviation²`), so consistency across tumors rules
   out *tumor-specific distortion* but does not validate whole-genome
   averaging as an appropriate metric for localized signal in general.
   Revised to state the narrower, supported claim and leave the broader
   question explicitly open.
2. "Peak amplitude is comparable" and "coherent chromosome-wide
   elevation" overclaimed precision the underlying statistics (max, p99,
   and an unsigned per-chromosome mean) cannot provide -- p99 is
   systematically somewhat lower in Tumor_01/04, and "chromosome mean
   absolute deviation" cannot distinguish true chromosome-wide coherence
   from mixed-sign changes, a strong focal sub-region, or noise. Revised
   to "higher measured burden," with the caveat stated explicitly, and
   Tumor_08 no longer described as part of one homogeneous group with
   Tumor_01/04.
3. The Tumor_10 "correlation-coherence problem" was reworded to "a
   distinct low seed-correlation component in the frozen classifier
   output" -- avoiding an unsupported claim about biological cause -- and
   the Control-vs-Tamoxifen difference was explicitly marked as
   descriptive, not evidence of a treatment effect.
4. The Point 3 boundary was sharpened: even a future pre/post-denoise
   comparison could only describe what the classifier received, not
   definitively separate weak biological signal from upstream
   attenuation.
5. The Point 6 recommendation was broadened from a specific proposed fix
   (extent-normalized score / extent-adaptive threshold) to a neutral list
   of candidate directions, since a specific proposal risked being
   selected because it would "rescue" the disagreement tumors rather than
   validated independently -- exactly the kind of overreach the task asked
   to avoid.

All five corrections are incorporated into the sections above.
