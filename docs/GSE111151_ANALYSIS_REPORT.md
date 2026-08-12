# GSE111151 analysis report: independent tamoxifen-resistance validation

**Date written:** 2026-08-12. All numbers below are read directly from
the committed output tables listed in each section, not hand-typed.

## 1. Dataset design

4 isogenic breast cancer cell lines (MCF-7, T-47D, ZR-75-1, BT-474), each
with 1 parental sample and 1-2 independently-derived tamoxifen-resistant
sublines: **4 parental, 7 resistant, 11 samples total**. Design is
**cell-line-blocked** (`~ cell_line + resistance_status`), not a naive
unpaired two-group test -- QC (section 2) shows samples cluster almost
entirely by cell-line identity (PC1+PC2 = 67.8% of variance), so ignoring
cell-line would confound resistance status with baseline cell-line
differences. Within each cell-line block the pairing is genuine (the
resistant sublines are directly derived from that block's own parental
culture), unlike GSE240112's PT/RT comparison. This is **genuine
tamoxifen resistance** (confirmed from the paper title/abstract, PMID
30143015), not a different endocrine-resistance state. Molecular-subtype
caveat: BT-474 is HER2-amplified, unlike the other 3 lines.

## 2. QC

**PASS.** Library sizes 30.9M-75.0M raw counts, no catastrophic outlier.
PCA and the sample-correlation dendrogram both show samples clustering
tightly into 4 groups matching cell-line identity exactly (`results/figures/gse111151/final_review/pca.png`,
`correlation_heatmap.png`), confirming the blocked design is both
necessary and working as intended. **Important correction applied during
this run:** an initial sample-level display used naive (raw-library-size)
CPM, which produced a spurious apparent direction for several candidates
because per-sample RNA composition differs enough that edgeR's TMM norm
factors range 0.82-1.30 across these 11 samples (e.g. MCF-7 parental
TMM factor 1.30 vs. MCF-7_Tam1 0.97). All sample-level values and figures
in this report use TMM-adjusted log2(CPM+1) (`results/tables/gse111151/tmm_norm_factors.tsv`),
consistent with what the statistical model actually estimates -- see
`docs/GSE111151_CODEX_REVIEW.md` for detail on how this was caught.

## 3. USP34

1. **Expressed/testable?** Yes -- passes `filterByExpr` in the
   cell-line-blocked model (`results/tables/gse111151/candidate_table.tsv`).
2. **Effect direction:** up in tamoxifen-resistant sublines.
3. **log2FC:** +0.162.
4. **p-value:** 0.213 (nominal).
5. **Candidate-set BH FDR:** 0.936 (not significant -- and no candidate
   in this panel reaches FDR<0.05; see section 6).
6. **Consistent across replicates?** 3 of 4 cell-line blocks (MCF-7,
   T-47D, BT-474) show the resistant subline(s) higher than that line's
   own parental; only ZR-75-1 is essentially flat/slightly lower
   (`results/figures/gse111151/final_review/usp34_sample_level.png`).
   Not a single-cell-line artifact.
7. **Does GSE111151 independently support the GSE118713 bulk result?**
   Directionally yes (both up), but GSE111151 does not itself reach
   statistical significance at any level (candidate FDR=0.936, nominal
   p=0.213) -- this is directional reinforcement, not independent
   statistical confirmation.
8. **Does it change confidence in USP34?** Modestly, positively: this is
   now the *third* independent human/cell-line RNA dataset (GSE118713,
   GSE240112 tumor-cell, GSE111151) showing the same "up in
   resistance-associated context" direction, with reasonably consistent
   support across 3 of 4 cell-line backgrounds here. It remains
   unconfirmed by any single dataset's own statistics beyond GSE118713.
   See section 9 for the classification.

## 4. VEZF1

1. **Does GSE111151 reproduce the GSE240112 signal?** No -- GSE240112
   found VEZF1 up in recurrent tumors (log2FC=+1.15, candidate
   FDR=0.048, its only candidate-set-significant hit); GSE111151 finds
   VEZF1 **down** in resistant sublines (log2FC=-0.238).
2. **Is the effect direction consistent?** No, it is opposite between
   the two datasets.
3. **Is it statistically supported (in GSE111151)?** No -- nominal
   p=0.188, candidate FDR=0.936.
4. **Is the result driven by one sample?** Partially: MCF-7 (-0.60) and
   T-47D (-0.70) drive the negative direction; ZR-75-1 (+0.07) and
   BT-474 (+0.26) are flat-to-slightly-positive
   (`results/figures/gse111151/final_review/vezf1_sample_level.png`).
   Only 2 of 4 cell-line blocks agree with the overall (negative)
   coefficient -- this is a genuinely mixed, not uniformly consistent,
   picture, and per the fixed classification rule (section 9) does not
   qualify as "discordant" (that label is reserved for a nominally
   significant GSE111151 result opposing a *significant* GSE118713
   result, and GSE118713 itself was not significant for VEZF1 either:
   log2FC=+0.427, FDR=0.238) -- it is reported as **neutral**, with the
   opposite-direction discrepancy from GSE240112 stated plainly rather
   than smoothed over.

## 5. SUPT4H1

1. **Does GSE111151 reproduce the GSE240112 signal?** No meaningful
   signal either way -- GSE240112 showed a large, nominally-interesting
   effect (log2FC=+1.00, nominal p=0.026, candidate FDR=0.157);
   GSE111151 shows essentially nothing (log2FC=+0.018, nominal p=0.916).
2. **Is the effect direction consistent?** Nominally the same sign
   (both positive) but GSE111151's effect is indistinguishable from
   zero.
3. **Is it statistically supported?** No (candidate FDR=0.936).
4. **Is the result driven by one sample?** The only cell line showing a
   sizeable per-sample deviation is MCF-7 (-0.66, actually the opposite
   direction from the tiny overall positive coefficient); T-47D, ZR-75-1,
   and BT-474 are all small and mixed in sign
   (`results/figures/gse111151/final_review/supt4h1_sample_level.png`).
   No consistent signal in either direction.

**Any other candidate that unexpectedly becomes stronger?** ICK: nominal
p=0.288 (still not <0.05) with 3/4 cell lines consistent, classified
`directionally_supportive_but_weak` alongside USP34 -- but its direction
in GSE111151 (log2FC=-0.168, down in resistant) does not match its
GSE240112 direction (log2FC=+0.358, up, ns) or its GSE118713 direction
(log2FC=+0.107, up, ns), so this is not a coherent cross-dataset signal,
just the second gene (after USP34) to clear GSE111151's own weak-support
bar in isolation.

## 6. All 13 candidates

(Cell-line-blocked resistant-vs-parental; sorted by candidate-set BH FDR.
Interpretation = `results/tables/gse111151/candidate_classification.tsv`,
a fixed, documented rule -- see `src/gse111151_evidence_classification.py`
docstring -- not a per-gene judgment call. **No candidate reaches
candidate-set FDR<0.05 in this panel.**)

| Gene | log2FC | p (nominal) | Candidate FDR | Direction | Cell lines consistent | Interpretation |
|---|---|---|---|---|---|---|
| VEZF1 | -0.238 | 0.188 | 0.936 | down | 2/4 | neutral -- no additional support |
| USP34 | +0.162 | 0.213 | 0.936 | up | 3/4 | directionally supportive but weak |
| ICK | -0.168 | 0.288 | 0.936 | down | 3/4 | directionally supportive but weak |
| EIF4ENIF1 | -0.129 | 0.322 | 0.936 | down | 4/4 | neutral -- no additional support |
| TLK2 | +0.073 | 0.530 | 0.936 | up | 4/4 | neutral -- no additional support |
| HMGB1 | -0.131 | 0.454 | 0.936 | down | 2/4 | neutral -- no additional support |
| TSR3 | +0.094 | 0.632 | 0.936 | up | 2/4 | neutral -- no additional support |
| KDM1A | +0.071 | 0.685 | 0.936 | up | 3/4 | neutral -- no additional support |
| TADA2B | +0.056 | 0.729 | 0.936 | up | 2/4 | neutral -- no additional support |
| CTDNEP1 | -0.044 | 0.802 | 0.936 | down | 2/4 | neutral -- no additional support |
| SUPT4H1 | +0.018 | 0.916 | 0.936 | up | 3/4 | neutral -- no additional support |
| PET117 | +0.018 | 0.936 | 0.936 | up | 3/4 | neutral -- no additional support |
| USP17L29 | NA | NA | NA | not tested | NA | untestable |

Every tested candidate's candidate-set BH FDR collapses to the same
value (0.936): with no p-value close to the significance threshold, this
is the mathematically correct BH behavior for an essentially null
12-gene family, not a computation error (independently reproduced during
Codex review; see `docs/GSE111151_CODEX_REVIEW.md`).

PAICS (benchmark, not in the 13-gene family): log2FC=-0.221, p=0.199,
genome-wide FDR=0.618 -- not significant, down in resistant.

## 7. Independently supported candidates

**None** reach candidate-set FDR<0.05 in GSE111151. Two candidates
(USP34, ICK) are classified `directionally_supportive_but_weak` under
the fixed rule (nominal p<0.3 AND >=3/4 cell lines consistent) --
neither is statistically confirmed on its own, and only USP34's
direction coheres with its own prior GSE118713/GSE240112 support (ICK's
does not, section 3/5).

## 8. Cross-dataset shortlist

Full table: `results/tables/gse111151/integrated_evidence_5layer.tsv`
(CRISPR, GSE118713, GSE245601, GSE240112, GSE111151 -- no composite
score computed).

- **Strongest cross-dataset candidate (most datasets pointing the same
  direction, even if not all significant): USP34.** Up-in-resistance
  direction in GSE118713 (significant, FDR=0.0073), GSE240112 tumor-cell
  (ns), and now GSE111151 (ns, but 3/4 cell lines consistent). Only
  GSE245601 (acute 12h) shows a small opposite-sign, itself
  nonsignificant, effect.
- **Strongest functional candidate:** by CRISPR effect size alone, KDM1A
  (-2.167, FDR=0.0004) and TADA2B (-2.064, FDR=0.015) have the largest
  magnitude screen hits; USP34's screen effect (-1.391, FDR=0.042) is
  more modest but is the one with the most consistent downstream RNA
  support.
- **Strongest human-resistance-expression candidate:** two different
  candidates hold the FDR<0.05 record in two different RNA layers, and
  neither should be read as "the" strongest without qualifying which
  layer: USP34 has a candidate-set-significant GSE118713 bulk result
  (FDR=0.0073); VEZF1 has the only candidate-set-significant *single-cell*
  result (GSE240112 tumor-cell, FDR=0.048), which did **not** reproduce
  in GSE111151 (opposite sign, ns) and should not be described as
  validated. USP34's GSE118713 significance, by contrast, is reinforced
  (not contradicted) by directionally-consistent-but-nonsignificant
  results in both GSE240112 and GSE111151.
- **Candidates supported by >=2 genuinely independent datasets (same
  direction, allowing non-significant support):** USP34 (GSE118713 sig +
  GSE240112 ns + GSE111151 ns, all same sign).
- **Contradictory evidence:** VEZF1 (GSE240112 significant up vs.
  GSE111151 ns down); ICK (GSE111151's own weak support does not match
  its GSE240112 or GSE118713 direction).
- **Candidates not prioritized by current RNA evidence** (neutral or no
  support in every RNA layer checked so far, despite real CRISPR screen
  hits -- this reflects the current evidence, not a decision to exclude
  them permanently; per Absolute Rule 8 carried through this project,
  nonsignificance does not refute the CRISPR screen result itself):
  CTDNEP1, HMGB1, TADA2B, TSR3, KDM1A, PET117, TLK2, EIF4ENIF1.
  USP17L29 remains untestable in every
  RNA dataset in this project (GSE118713, GSE240112, GSE111151 all show
  zero/unavailable counts) and cannot be evaluated by transcriptomics at
  all.

## 9. What changed

- **USP34: modestly strengthened.** A third RNA dataset shows the same
  direction as GSE118713, with reasonable (3/4) cell-line consistency,
  though still not itself statistically significant. This does not
  upgrade USP34 to "validated," but it is the most consistent
  multi-dataset story of any candidate in this project to date.
- **VEZF1: did not reproduce.** Its one significant RNA hit (GSE240112)
  goes the opposite direction in GSE111151 (nonsignificantly). This
  weakens confidence that VEZF1's GSE240112 signal reflects a general
  tamoxifen-resistance-associated expression change, though it does not
  disprove it (GSE240112 is a different biological context -- human
  primary-vs-recurrent tumors -- from GSE111151's cell-line acquired
  resistance).
- **No new candidate emerged as strong.** ICK crosses GSE111151's own
  weak-support bar but has no coherent cross-dataset story.

## 10. Figures

`results/figures/gse111151/final_review/`: `pca.png`,
`correlation_heatmap.png`, `candidate_effect_size.png`,
`candidate_heatmap.png`, `usp34_sample_level.png`,
`vezf1_sample_level.png`, `supt4h1_sample_level.png`, `volcano.png`.

## 11. Limitations

- Only 4 cell-line blocks (biological replicates in the blocking sense);
  statistical power is limited, and a nonsignificant result here does
  not refute prior CRISPR/GSE118713/GSE245601/GSE240112 evidence.
- No independent biological replicate of any parental state (1 parental
  sample per cell line); Tam1/Tam2 (where both exist) are independently
  -derived resistant sublines from the same parental culture, not
  technical replicates.
- Long-term (8-12 month) tamoxifen-selected cell-line derivatives may
  acquire cross-resistance or culture-adaptation changes unrelated to
  tamoxifen specifically; this cannot be distinguished from
  tamoxifen-specific resistance biology with this design alone.
- Molecular-subtype heterogeneity (BT-474 is HER2-amplified) is
  addressed by cell-line blocking, not eliminated.
- USP17L29 could not be evaluated (filtered out by `filterByExpr`; zero
  detection is consistent with its behavior in every other RNA dataset
  in this project).
- A real normalization bug (naive vs. TMM-adjusted per-sample CPM) was
  caught and corrected during this run before any figure or
  interpretation was finalized -- see `docs/GSE111151_CODEX_REVIEW.md`.

## 12. Codex verdict

**PASS WITH NOTES.** 3 issues found, all in report wording (an
overclaimed "only candidate with FDR<0.05" statement, an overly strong
"should be dropped" recommendation, and forward references to this
review file before it existed) -- the underlying statistical
implementation (design matrix, TMM propagation, candidate BH, sample
-level consistency numbers, classification rule) was independently
re-verified and matched to numerical precision. Full findings:
`docs/GSE111151_CODEX_REVIEW.md`.

## 13. Tests

47 GSE111151-specific tests passed
(`python3 -m pytest tests/ -k gse111151 -q`); 572 passed for the full
project test suite (`python3 -m pytest tests/ -q`), confirming no
regression to any existing frozen module.

## 14. Git

See final commit hash reported at the end of this run.
