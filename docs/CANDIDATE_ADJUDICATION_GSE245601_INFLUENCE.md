# Candidate adjudication: why does GSE245601 have the largest leave-one-out impact?

**Date written:** 2026-08-12. This is a sensitivity audit
(`src/candidate_adjudication_gse245601_influence.py`) of an already
-reported observation (`docs/CROSS_DATASET_GENOMEWIDE_ANALYSIS_REPORT.md`:
"removing GSE245601 caused the largest disruption in the global
ranking"). It re-verifies that claim quantitatively and investigates the
mechanism. **No ranking rule is changed as a result of this audit.**

## Re-verification

Top-20 overlap between the main global ranking and each leave-one-out
variant (source: the already-frozen `ranking_stability.tsv`):

| dataset removed | Top-20 overlap with main |
|---|---|
| GSE245601 | **7/20** |
| GSE240112 | 12/20 |
| GSE118713 | 14/20 |
| GSE111151 | 16/20 |
| CRISPR | 17/20 |

Confirmed: GSE245601 removal disrupts the Top 20 more than removing any
other single dataset, including CRISPR.

## Question 1-2: why?

**Not because of coverage.** GSE245601 is testable for 17,987 genes --
neither the largest (GSE111151: 26,982) nor the smallest (GSE118713:
14,701) testable set, and the actual Top-20 genes are all Tier A (5/5
datasets testable), so removing any single dataset can only demote them
to Tier B, never below the >=3-testable eligibility floor. Coverage-tier
dropout is therefore not the mechanism for genes already in the Top 20
(`gse245601_influence_coverage_risk.tsv` shows GSE111151, not GSE245601,
has the most genes at coverage-dropout risk overall -- 2,041 vs. 1,341 --
but that is irrelevant to the already-Tier-A Top 20).

**Not because GSE245601 contributes an unusually large number of FDR<0.05
votes.** Quite the opposite: only 101/17,987 Track A genes (0.56%) and 0/17,987
Track B genes reach FDR<0.05, versus 30.3% for GSE118713 and 28.9% for
GSE240112 (tumor-cell). GSE245601 is one of the *lowest*-power datasets by
this measure, comparable to CRISPR (0.12%) and GSE111151 (0.16%). It is
not inflating the `n_datasets_fdr05` sort key disproportionately.

**The real mechanism is the median/mean-percentile tie-break, at a rank
density where small percentile shifts reorder many genes.** Among the
Top 20 genes themselves, GSE245601's own within-dataset percentile is
often very high (mean 0.81 across the Top 20, several genes >0.99: HRK,
MYBL1, GREB1, LARP6, MSMO1, FGD3, SYTL5, PTGER4) even though almost none
of them individually cross GSE245601's own FDR<0.05 bar -- percentile is
a purely ordinal, rank-based statistic, so a gene can sit in GSE245601's
own top 1% without ever reaching FDR<0.05 in a dataset where FDR<0.05 is
rare overall. Removing GSE245601 from a gene's 5-value percentile set
shifts its median/mean percentile by as much as 0.10-0.22 for several
Top-20 genes (`gse245601_influence_median_shift.tsv`), and at the
extreme, densely-packed top of a 15,255-gene ranking, a percentile shift
of that size is enough to move a gene's rank by hundreds of positions.
This is a real, data-driven pattern (many Top-20 genes happen to rank
well within GSE245601 too), not an artifact of the ranking code.

## Question 3: does GSE245601 receive hidden extra weight?

**No.** It contributes exactly one percentile value and one FDR-count
vote to the hierarchy, identically to every other dataset (verified
during the Phase 27 Codex review of the cross-dataset genome-wide
integration, and re-confirmed here: GSE245601's `n_datasets_fdr05`
contribution is the *smallest* of the five, not the largest). Its
outsized *leave-one-out* impact is a consequence of which specific genes
happen to rank highly in it, not of any structural over-weighting.

## Question 4: would the Track A/B collapsing scheme matter?

Rebuilding the full ranking hierarchy (same frozen sort, `assign_coverage_tier`
+ `build_global_ranking`, unmodified) three times, substituting only the
GSE245601 percentile/FDR columns (`gse245601_influence_track_scheme_comparison.tsv`):

| scheme | Top-20 overlap with the frozen mean-of-A/B scheme |
|---|---|
| frozen mean of Track A + Track B | 20/20 (baseline) |
| **Track A only** | **20/20 -- identical Top 20** |
| **Track B only** | **9/20 -- materially different Top 20** |

**Track A alone reproduces the frozen Top 20 exactly.** This means the
current mean-of-two-tracks scheme is, in practice, almost entirely
carried by Track A (all epithelial cells, the larger and less noisy
population) for the specific question of who is in the Top 20; Track B
(strict malignant cells, much smaller per-patient cell counts) would, if
used alone, produce a substantially different list. This is consistent
with Track B's near-total lack of FDR<0.05 hits (0/17,987) reported
above -- Track B is simply a noisier signal at the sample sizes
available in GSE245601. This asymmetry was already flagged as a
documented (not hidden) limitation in the Phase 27 Codex review
("GSE245601's Track A/Track B mean-of-two-percentiles vs. single-track
genes creates a real, documented missing-data asymmetry"); this audit
adds the specific, quantified consequence -- Track B alone would reorder
more than half the Top 20 -- without changing the frozen ranking, per
Phase 18's explicit instruction to treat this as a sensitivity check only.

## Conclusion

GSE245601's outsized leave-one-out effect on the global Top 20 is real,
mechanistically explained (a percentile-tie-break effect at a densely
-competitive rank density, not a coverage or FDR-count weighting
artifact), and does not indicate hidden extra weight for GSE245601 in the
frozen ranking rule. It does, however, mean the Top 20 is more sensitive
to GSE245601's Track A pseudobulk than any other single evidence source
in the collection -- a fact worth carrying into Phase 20's decision table
as a caveat on Top-20 genes whose strongest support comes primarily from
GSE245601, and worth acknowledging explicitly rather than treating the
global Top 20 as equally robust to every dataset choice.
