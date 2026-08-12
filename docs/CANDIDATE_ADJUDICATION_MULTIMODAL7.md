# Candidate adjudication: seven-gene head-to-head

**Date written:** 2026-08-12. Every number below is read directly from
`results/tables/candidate_adjudication/multimodal7_exact_evidence.tsv`
and `results/tables/cross_dataset_genomewide/ranking_stability.tsv`
(both re-verified against the original per-dataset files in Phase 3,
0 mismatches out of 98 checked values). No mechanistic claim is made from
expression data alone (Phase 32 wording rules).

## Cross-cutting finding before the per-gene profiles

**CRISPR direction splits the seven genes 2-5, and this split matters more
than the shared MULTIMODAL_STRONG label suggests.** USP34 and VEZF1 have
negative CRISPR effect sizes (`sensitising_KO`: knockout was relatively
depleted under 4-OHT). CUX1, DPP9, LZTR1, SOX2, and TFAP2C have *positive*
effect sizes (`tolerance_associated_KO`: knockout was relatively favored
under 4-OHT) -- the functional-fitness direction a knockout-based
sensitisation strategy would need to *avoid*, not pursue, if the goal is
a gene whose loss re-sensitises resistant cells. This is not evident from
the MULTIMODAL_STRONG label alone and is the single most important fact
for Phase 21-22's therapeutic shortlist.

**The category also rests on two distinct evidentiary paths** (see
`docs/CANDIDATE_ADJUDICATION_CATEGORY_AUDIT.md`): CUX1 and SOX2 have two
independent resistance datasets at FDR<0.05; DPP9, LZTR1, TFAP2C, USP34,
and VEZF1 have exactly one.

## Gene-by-gene profiles

### CUX1

- **CRISPR:** effect +1.27, FDR 0.0450, rank 21/19,103, `tolerance_associated_KO`.
- **Resistance RNA:** GSE118713 log2FC -0.90, FDR 0.00079 (rank 1,366); GSE240112 (tumor) log2FC -0.98, FDR 0.0326 (rank 4,565); GSE111151 log2FC -0.47, FDR 0.353 (not significant, but direction consistent in 3/4 cell lines). Pattern `DOWN_DOWN_DOWN`, consensus `all_down`, 2/3 datasets FDR<0.05.
- **Acute human response (GSE245601):** epithelial FDR 0.80, malignant FDR 0.60 -- no acute signal.
- **Human recurrence (GSE240112):** significant and direction-consistent with the epithelial-sensitivity track (log2FC -0.84).
- **Cross-dataset:** 3/5 datasets FDR<0.05, 3/5 in their own top-10%, global rank **10/15,255**, `MODERATELY_STABLE` (best rank 3, worst rank 428 -- the worst case is specifically leave-CRISPR-out, since CUX1's global-rank advantage partly depends on CRISPR remaining significant; removing GSE245601 instead *improves* its rank to 3, since GSE245601 was a drag, not a support, for this gene).
- **STRENGTHS:** the strongest resistance-RNA support of the seven (2 independent datasets, both FDR<0.05, consistent direction); highest global rank (10th of 15,255).
- **LIMITATIONS:** CRISPR direction is `tolerance_associated_KO` -- knockout was *favored*, not depleted, under 4-OHT; not a sensitising hit by the functional-fitness definition. GSE111151 (the second resistance cell-line panel) does not reach significance.
- **MOST DEFENSIBLE ROLE:** resistance biomarker / pathway candidate. Not a knockout-sensitisation candidate on the CRISPR evidence alone.

### SOX2

- **CRISPR:** effect +1.51, FDR 0.0122, rank 7/19,103, `tolerance_associated_KO`.
- **Resistance RNA:** GSE118713 log2FC +0.69, FDR 0.0093 (rank 2,704); GSE240112 (tumor) log2FC **+7.44**, FDR 5.6e-05 (rank 328); GSE111151 log2FC +0.89, FDR 0.632 (not significant, 3/4 cell lines directionally consistent). Pattern `UP_UP_UP`, consensus `all_up`, 2/3 FDR<0.05.
- **Acute human response:** epithelial FDR 0.69, malignant FDR 0.69 -- no acute signal.
- **Human recurrence:** the largest, most significant effect of any of the seven genes in any dataset. Sample-level inspection (Phase 10) confirms this is broad across all three recurrent-tumor pseudobulk samples, not driven by one sample -- normalized pseudobulk log2(CPM+1) is 0.000/0.000/0.118 in PT1/PT2/PT3 versus 2.500/2.627/2.428 in RT1/RT2/RT3. RT3's raw count (5) is small in absolute terms (reflecting its lower cell count) and is not, by itself, distinguishable from PT3's raw count (also 5) -- but RT3's *normalized* value is comparable to RT1/RT2, confirming the effect is broad across all three recurrent samples and not an RT3-driven artifact.
- **Cross-dataset:** 3/5 FDR<0.05, 2/5 top-10%, global rank **16/15,255**, `MODERATELY_STABLE` (best rank 7, worst rank 1,013 on leave-CRISPR-out).
- **STRENGTHS:** by far the largest and cleanest resistance-recurrence effect among the seven; two independent resistance datasets FDR<0.05; 2nd-highest global rank.
- **LIMITATIONS:** same `tolerance_associated_KO` CRISPR direction issue as CUX1 -- functionally, loss of SOX2 was *not* a sensitising event in the screen. GSE111151 not significant.
- **MOST DEFENSIBLE ROLE:** resistance/recurrence biomarker candidate, not a functional sensitisation candidate.

### TFAP2C

- **CRISPR:** effect +1.54, FDR 0.0482, rank 22/19,103, `tolerance_associated_KO`.
- **Resistance RNA:** GSE118713 log2FC -1.09, FDR 0.00096 (rank 1,430); GSE240112 (tumor) FDR 0.461 (not significant); GSE111151 log2FC -0.54, FDR 0.376 (not significant, only 2/4 cell lines directionally consistent). Pattern `DOWN_DOWN_DOWN`, consensus `all_down`, but only 1/3 FDR<0.05.
- **Acute human response:** epithelial FDR 0.162 -- the closest of the seven to acute significance, still not significant.
- **Cross-dataset:** 2/5 FDR<0.05, global rank **272/15,255**, `DATASET_DEPENDENT` (best rank 13 on leave-GSE111151-out, worst rank 1,803 on leave-CRISPR-out -- its global-rank standing depends heavily on CRISPR remaining in the comparison).
- **STRENGTHS:** one clearly significant resistance dataset (GSE118713); consistent direction sign across all three resistance datasets even where not individually significant.
- **LIMITATIONS:** only 1/3 resistance datasets reaches FDR<0.05; GSE111151 cell-line consistency is the weakest of the seven (2/4); `tolerance_associated_KO` direction.
- **MOST DEFENSIBLE ROLE:** weaker resistance biomarker candidate; not a sensitisation candidate.

### DPP9

- **CRISPR:** effect +1.38, FDR 0.0417, rank 17/19,103, `tolerance_associated_KO`.
- **Resistance RNA:** GSE118713 log2FC -1.04, FDR 0.0021 (rank 1,768) -- the *only* significant resistance dataset. GSE240112 (tumor) FDR 0.802, GSE111151 FDR 0.954 (both essentially null). Pattern `DOWN_UP_UP` (GSE118713 down, the other two nominally up but not significant), consensus `majority_up` -- **note the direction pattern is not fully concordant**: GSE118713 (the one significant dataset) points down, while the two non-significant datasets point up.
- **Acute human response:** no signal (FDR 0.97, 0.76).
- **Cross-dataset:** 2/5 FDR<0.05 (CRISPR + GSE118713 only), global rank **1,090/15,255**, `DATASET_DEPENDENT`, **0/5 leave-one-out Top-20 appearances** -- drops out of the Top 20 under every single-dataset removal.
- **STRENGTHS:** CRISPR and GSE118713 both individually significant.
- **LIMITATIONS:** the weakest evidence base of the seven outside CRISPR+GSE118713; direction pattern is internally inconsistent across resistance datasets; zero robustness to any single dataset's removal.
- **MOST DEFENSIBLE ROLE:** functional-with-single-RNA-support candidate; weakest multimodal case of the seven.

### LZTR1

- **CRISPR:** effect +1.46, FDR 0.0124, rank 8/19,103, `tolerance_associated_KO`.
- **Resistance RNA:** GSE118713 log2FC +1.11, FDR 0.0151 (rank 3,080) -- only significant resistance dataset. GSE240112 FDR 0.608, GSE111151 FDR 0.717 (both null). Pattern `UP_UP_UP`, consensus `all_up`, direction at least fully concordant even where not significant.
- **Acute human response:** no signal.
- **Cross-dataset:** 2/5 FDR<0.05, global rank **1,173/15,255**, `DATASET_DEPENDENT`, 0/5 leave-one-out Top-20 appearances.
- **STRENGTHS:** fully concordant direction across all three resistance datasets; strong CRISPR (rank 8 of the full screen).
- **LIMITATIONS:** only one resistance dataset reaches significance; no human-tumor signal at all.
- **MOST DEFENSIBLE ROLE:** functional candidate with modest, direction-consistent (but not independently significant) resistance-RNA support.

### USP34

- **CRISPR:** effect -1.39, FDR 0.0417, rank 18/19,103, `sensitising_KO` -- **one of only two of the seven with the functionally desirable direction.**
- **Resistance RNA:** GSE118713 log2FC +0.59, FDR 0.0073 (rank 2,530) -- only significant resistance dataset. GSE240112 FDR 0.228, GSE111151 FDR 0.632 (both null). Pattern `UP_UP_UP`, consensus `all_up`.
- **Acute human response:** no signal (FDR 0.90).
- **Cross-dataset:** 2/5 FDR<0.05, global rank **994/15,255**, `DATASET_DEPENDENT`, 0/5 leave-one-out Top-20 appearances.
- **STRENGTHS:** correct (sensitising) CRISPR direction; fully concordant resistance direction across all three datasets.
- **LIMITATIONS:** only GSE118713 reaches significance among the three resistance datasets; no human-tumor support; the weakest global rank of the two sensitising-direction genes.
- **MOST DEFENSIBLE ROLE:** the more defensible of a knockout-sensitisation therapeutic candidate among the seven, on direction grounds, though resistance-RNA support is limited to one dataset.

### VEZF1

- **CRISPR:** effect -1.60, FDR 0.0373, rank 11/19,103, `sensitising_KO` -- the other of the two sensitising-direction genes.
- **Resistance RNA:** GSE240112 (tumor) log2FC +1.15, FDR 0.0195 (rank 3,792) -- only significant resistance dataset (note: for VEZF1 this is GSE240112, not GSE118713 as for the other five). GSE118713 FDR 0.2375 (not significant -- the *secondary* TAMR_vs_FASR contrast is significant, FDR 0.00023, but that contrast is not this dataset's independent vote per the frozen methodology). GSE111151 log2FC -0.24, FDR 0.608, and its direction *disagrees* (down, vs. up in the other two). Pattern `UP_UP_DOWN`, consensus `majority_up` (not full concordance).
- **Acute human response:** no signal.
- **Cross-dataset:** 2/5 FDR<0.05, global rank **1,138/15,255**, `DATASET_DEPENDENT`, 0/5 leave-one-out Top-20 appearances.
- **STRENGTHS:** correct (sensitising) CRISPR direction; the resistance-significant dataset is GSE240112 (human primary-vs-recurrent tumor context) rather than a cell-line panel, adding a different evidence type than USP34's.
- **LIMITATIONS:** GSE111151 direction actively disagrees with the other two resistance datasets; GSE111151 cell-line consistency is the weakest tie of the seven (2/4); GSE118713's own primary contrast is not significant.
- **MOST DEFENSIBLE ROLE:** knockout-sensitisation candidate (correct direction) with a genuine directional inconsistency across resistance datasets that should not be glossed over.

## Head-to-head, ranked ONLY for multimodal therapeutic follow-up

This ranking uses the transparent hierarchy from Phase 22 (CRISPR
sensitising direction first, then CRISPR strength, then resistance-dataset
count/consistency, then human evidence, then stability) and applies **only**
to the question of knockout-sensitisation follow-up -- it is not a
restatement of global rank or of general "importance."

1. **VEZF1** -- sensitising direction, CRISPR FDR 0.037 (band: FDR<0.05, tied with USP34 on criterion 2), one significant resistance dataset (band: tied with USP34 on criterion 3), and -- unlike USP34 -- significant human-tumor support (GSE240112 FDR 0.019, criterion 4). This is what breaks the tie in VEZF1's favor, ahead of criterion 8.
2. **USP34** -- sensitising direction, CRISPR FDR 0.042 (same band as VEZF1), one significant resistance dataset (same band), but no human-tumor significance. Its resistance direction is fully concordant across all three datasets with no disagreement, unlike VEZF1 -- a genuine advantage on criterion 8 ("absence of major contradiction"), but criterion 8 is only reached after criteria 2-7 are exhausted, and criterion 4 already separates the two genes.
3. CUX1, SOX2, TFAP2C, DPP9, LZTR1 -- all `tolerance_associated_KO`: knockout was functionally *favored*, not depleted, under tamoxifen. These do not satisfy the first, non-negotiable criterion (sensitising direction) for this specific question, regardless of how strong their resistance-RNA evidence is (CUX1 and SOX2 in particular have the best resistance-RNA evidence of all seven). They remain excellent resistance-biomarker candidates (see Phase 21 List B).

**Caveat carried forward, not used to reverse the ranking:** VEZF1's
GSE111151 direction (down) disagrees with its other two resistance
datasets (up, up) -- a real inconsistency that does not disqualify it
under this hierarchy (criterion 4 resolves the ranking before criterion 8
is reached) but should be disclosed alongside VEZF1 wherever it is
recommended.
