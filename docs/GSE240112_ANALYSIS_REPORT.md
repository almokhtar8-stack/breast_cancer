# GSE240112 analysis report: primary vs. tamoxifen-treated recurrent ER+ breast tumor scRNA-seq

**Date written:** 2026-08-12. All numbers below are read directly from
the committed output tables listed in each section, not hand-typed.

## 1. Study design

3 primary (PT) and 3 tamoxifen-treated recurrent (RT) ER+ breast tumor
scRNA-seq samples, plus 2 normal breast tissue (NT) samples used only for
optional secondary context. **Unpaired**: PT tissue came from OriGene
Technologies Inc., RT tissue from the Ontario Tumor Bank -- different
source institutions, no pairing statement in GEO metadata or the paper
Methods (`docs/GSE240112_DATA_AUDIT.md` section 2). Important
treatment-history caveat: RT is tamoxifen-treated recurrent disease by
group definition, but other systemic therapy history (chemotherapy,
other endocrine agents) is not available in the public metadata -- a real,
unresolved confounder.

## 2. Data

Public files used: author-processed h5Seurat objects (`NTs.h5seurat`,
`TTs_cancer_060223.h5seurat`, Dropbox-hosted, linked from the author
GitHub repo) plus raw GEO Cell Ranger matrices for the 6 PT/RT samples
(secondary analysis only). Raw UMI counts (not just processed/normalized
values) were available and used throughout. Final cells: 9,942 tumor
cells (`cell.annot == "Breast cancer cells"` for all of them) across the
6 PT/RT samples (PT1=1,029, PT2=1,975, PT3=1,442, RT1=2,721, RT2=2,597,
RT3=178); 0 cells lost in extraction (`results/tables/gse240112/extraction_log.tsv`).

## 3. Malignant/tumor-cell definition

**Author labels (Phase 7 Case A).** The `TTs_cancer_060223.h5seurat`
object is, by construction, the authors' own cancer-cell subset
(`cell.annot` reads `"Breast cancer cells"` for all 9,942 cells),
generated per the paper Methods via `Seurat::FindAllMarkers` against
curated marker panels (PanglaoDB, CellMarker, ProteinAtlas, literature) --
not a CNV-based method. Used directly, no independent re-classification.
**Confidence/limitation:** RT3 contributes only 178 cells (vs. 1,029-2,721
for the other 5 samples); pseudobulk QC shows it clusters with RT1/RT2 on
PC1 (85.7% of variance) and correlates more strongly with them (Spearman
0.92-0.93) than with any PT sample (0.87-0.88), so it is retained as a
smaller, noisier, genuine replicate, not excluded.

## 4. Pseudobulk QC

**PASS.** PC1 (85.7% of variance) cleanly separates all 3 PT samples
(coordinates ~-83 to -84) from all 3 RT samples (~+79 to +92); RT3
groups with RT1/RT2 despite its small size. No sample shows a
catastrophic technical signature; none excluded.
(`results/tables/gse240112_pseudobulk/tumor_cell_pca_coordinates.tsv`,
`tumor_cell_sample_correlation.tsv`, `tumor_cell_qc_summary.tsv`.)

## 5. USP34

- **Detected:** yes. 88.8% of tumor cells (8,826/9,942), in all 6
  samples (`results/tables/gse240112/candidate_detection_audit.tsv`).
- **PT expression** (per-sample log2(CPM+1), tumor-cell pseudobulk):
  PT1=6.87, PT2=6.76, PT3=7.09 (mean 6.90).
- **RT expression:** RT1=7.06, RT2=6.86, RT3=7.66 (mean 7.19).
- **RT-vs-PT log2FC:** +0.400 (up in recurrent tumors).
- **p-value:** 0.118 (nominal, not itself <0.05).
- **Candidate-set BH FDR:** 0.473 (not statistically significant at
  FDR<0.05).
- **Individual-sample consistency (corrected after Codex review -- see
  `docs/GSE240112_CODEX_REVIEW.md`):** PT range is 6.760-7.088. Only
  **RT3** (7.656) clearly exceeds the full PT range; RT1 (7.060) and RT2
  (6.858) both fall inside it. RT3 is also the pseudobulk sample with the
  fewest contributing cells (178, vs. 1,029-2,721 for the other five) and
  the smallest library size (1.14M vs. 43-65M UMIs). The positive
  RT-vs-PT log2FC is therefore **substantially, though not entirely,
  driven by the smallest and noisiest sample** -- excluding RT3
  mentally, RT1 and RT2 both sit within the PT range and would not on
  their own suggest an RT-vs-PT difference. This is a genuine
  fragility in the result, not a detail to gloss over.
- **All-epithelial sensitivity:** log2FC=+0.457 (same direction),
  genome-wide FDR=0.196 (not significant). This is **not an independent
  replication** -- it is the same six sequencing libraries reprocessed
  with a different (broader) cell-selection pipeline, substantially
  overlapping the tumor-cell track's cells, so it is a sensitivity check
  on the cell-population definition, not a second independent dataset.
- **Final interpretation:** directionally consistent with the GSE118713
  bulk finding (up in the resistance-associated context), but the
  GSE240112 tumor-cell result itself is not statistically significant, is
  largely carried by one small, noisy sample, and its all-epithelial
  "sensitivity" companion is not an independent confirmation. See section
  9 for the corrected final status.

## 6. 13 candidates table

(RT-vs-PT tumor-cell pseudobulk; sorted by candidate-set BH FDR. "Epi
sensitivity" = all-epithelial track log2FC, from the same six libraries
reprocessed with a broader cell selection, not an independent dataset --
see section 11; "dir. agree" = same sign as tumor-cell track.)

Caveat on the "Interpretation" column: this is a descriptive, post-hoc
classification (not a preregistered PREANALYSIS.md rule -- only the
FDR<0.05 significance threshold itself was frozen in advance), analogous
to this project's existing GSE245601 `evidence_class` scheme. For genes
where the prior GSE118713 bulk result was not itself significant (which
is most of them, including VEZF1, the one "strengthened" gene here), the
"directionally supportive" tier has no pre-specified expected direction
to be judged against -- it reflects only within-GSE240112 evidence
strength, not cross-dataset concordance (cross-dataset concordance is
handled separately in sections 8-9).

Interpretation column = `results/tables/gse240112/candidate_classification.tsv`
(fixed, documented rule -- see `src/gse240112_evidence_classification.py`
docstring -- not a per-gene judgment call).

| Gene | log2FC | p (nominal) | Candidate FDR | Direction | Epi sensitivity log2FC | Dir. agree | Interpretation |
|---|---|---|---|---|---|---|---|
| VEZF1 | +1.149 | 0.0040 | **0.048** | up in RT | +1.172 | yes | strengthened |
| SUPT4H1 | +1.003 | 0.0261 | 0.157 | up in RT | +0.925 | yes | directionally supportive but weak |
| USP34 | +0.400 | 0.1183 | 0.473 | up in RT | +0.457 | yes | neutral -- no additional support |
| TADA2B | +0.428 | 0.2417 | 0.580 | up in RT | +0.547 | yes | neutral -- no additional support |
| TSR3 | -0.427 | 0.2408 | 0.580 | down in RT | -0.450 | yes | neutral -- no additional support |
| ICK | +0.358 | 0.2978 | 0.596 | up in RT | +0.515 | yes | neutral -- no additional support |
| KDM1A | +0.192 | 0.4514 | 0.677 | up in RT | +0.353 | yes | neutral -- no additional support |
| PET117 | -0.246 | 0.4288 | 0.677 | down in RT | -0.110 | yes | neutral -- no additional support |
| EIF4ENIF1 | +0.063 | 0.7922 | 0.819 | up in RT | +0.144 | yes | neutral -- no additional support |
| CTDNEP1 | -0.121 | 0.6977 | 0.819 | down in RT | -0.050 | yes | neutral -- no additional support |
| HMGB1 | -0.155 | 0.7316 | 0.819 | down in RT | +0.022 | **no** | neutral -- no additional support (both near zero) |
| TLK2 | +0.057 | 0.8191 | 0.819 | up in RT | +0.124 | yes | neutral -- no additional support |
| USP17L29 | NA | NA | NA | not tested | NA | NA | untestable (see section 8 note) |

Note: USP34's log2FC (0.400) falls just short of the classifier's
"sizeable effect" threshold (|log2FC|>=0.5), so within GSE240112 alone it
is classified "neutral" by this fixed rule, distinct from its holistic,
cross-dataset final status in section 9 below (which is explicitly a
different, separately-specified question).

PAICS (benchmark, not in the 13-gene family): log2FC=+0.127, p=0.695,
genome-wide FDR=0.794 -- not significant, up in RT.

## 7. Which candidates strengthened

Classification is per the fixed rule in
`src/gse240112_evidence_classification.py` (candidate-set FDR<0.05 =
strengthened; not FDR-significant but nominal p<0.05, or |log2FC|>=0.5
with all-epithelial-track direction agreement = directionally supportive
but weak; a nominally significant result opposing a significant
GSE118713 bulk direction = discordant; everything else tested = neutral):

- **Statistically supported (candidate-set FDR<0.05):** VEZF1 only.
- **Directionally supportive but statistically weak:** SUPT4H1 (nominal
  p=0.026, |log2FC|=1.00, all-epithelial track agrees).
- **Discordant:** none.
- **Neutral -- no additional support:** USP34, CTDNEP1, EIF4ENIF1, HMGB1,
  KDM1A, PET117, TADA2B, ICK, TLK2, TSR3 (10 genes) -- this includes USP34,
  whose GSE240112 effect (log2FC=+0.40, nominal p=0.118) falls short of
  the fixed rule's thresholds on its own, even though it is directionally
  consistent with GSE118713 (see section 9 for the separate, holistic
  cross-dataset USP34 call).
- **Untestable:** USP17L29 (see section 8).

## 8. Cross-dataset comparison

Full table: `results/tables/gse240112/integrated_evidence_4layer.tsv`
(one row per frozen candidate, four independently-measured layers, no
composite score computed -- CRISPR = functional perturbation under 4-OHT;
GSE118713 = acquired-resistance cell-line state; GSE245601 = acute 12h ex
vivo response; GSE240112 = human primary-vs-recurrent context). USP34
across all four: CRISPR effect_size=-1.391 (FDR=0.042, sensitising);
GSE118713 TAMR-vs-MCF7 log2FC=+0.590 (FDR=0.0073, significant, up);
GSE245601 Track A (epithelial) log2FC=-0.033 (ns), Track B (malignant)
log2FC=-0.181 (ns); GSE240112 tumor-cell log2FC=+0.400 (candidate
FDR=0.473, ns), all-epithelial log2FC=+0.457 (ns).

USP17L29 is untestable in **every** RNA layer examined so far: GSE118713
bulk (`RNA_UNAVAILABLE`, filtered out), GSE245601 (not separately
audited here but consistent with the same paralog-family issue), and now
GSE240112 -- absent from the author-processed feature space (27,161
genes) despite being present in the raw upstream CellRanger GRCh38
reference by exact Ensembl ID match (`ENSG00000231637`), and confirmed
via direct inspection of all 6 raw PT/RT Cell Ranger matrices to have
**zero** counts in every cell in every sample. This is a genuine,
repeated, cross-dataset detectability limitation of this specific
paralog-family gene (chr8p21.3 tandem-repeat region), not a
symbol-mapping error specific to this analysis.

## 9. USP34 final status: **UNCHANGED-NEUTRAL**

**(Revised after Codex review -- the original draft of this section called
"MODESTLY STRENGTHENS" based on two claims that turned out to be
inaccurate; see `docs/GSE240112_CODEX_REVIEW.md` for the full finding and
correction.)**

GSE118713 bulk shows a significant, up-in-resistance USP34 signal
(FDR=0.0073). GSE245601 shows the opposite sign in both tracks (Track A
log2FC=-0.033, Track B log2FC=-0.181), though neither is itself
statistically significant. GSE240112 tumor-cell pseudobulk shows a
positive log2FC (+0.400) in the same direction as GSE118713, but this is
not itself significant (candidate FDR=0.473) and is substantially driven
by RT3, the smallest and noisiest of the six pseudobulk samples (178
cells; see section 5) -- RT1 and RT2 alone fall within the PT range. The
all-epithelial "sensitivity" result is not independent confirmation: it
is the same six libraries reprocessed with a broader, overlapping cell
selection, not a second dataset.

Taken together, GSE240112 does not provide a robust, sample-consistent,
or independent positive signal for USP34 beyond what GSE118713 already
showed, and it does not contradict GSE118713 either (the direction is
still nominally positive, just fragile). It also does not resolve the
GSE245601 discrepancy (which remains negative-but-nonsignificant). The
honest, conservative read is that GSE240112 leaves the overall USP34
picture **unchanged**: still anchored by the significant GSE118713 bulk
result and the CRISPR screen, neither meaningfully strengthened nor
weakened by this dataset. The conceptual hypothesis that USP34 may track
an established/chronic resistant state without acute (12h) induction
remains a plausible reading of the GSE245601-vs-GSE118713 contrast, but
GSE240112 -- given the RT3-driven fragility documented above -- should not
be cited as confirming it.

## 10. Figures

`results/figures/gse240112/final_review/` (see Phase 21 manifest) contains:
tumor-cell PCA/correlation/library QC, broad-compartment UMAP, USP34
tumor-cell UMAP + PT/RT-split UMAP + per-sample distribution + per-sample
detection, and the four 13-candidate overview figures (dot plot,
pseudobulk heatmap, effect-size plot, per-sample heatmap).

## 11. Limitations

- Unpaired, small (n=3 vs 3) design; no evidence of PT/RT patient pairing
  was found (different tissue-source institutions make a shared patient
  unlikely but not provably absent -- Table S1 was not retrieved), so the
  contrast was analyzed unpaired; no causal tamoxifen claim is supported
  by this cross-sectional comparison regardless.
- RT3 contributes only 178 tumor cells; retained per the frozen
  pre-analysis rule (it behaves as a genuine, if noisy, RT replicate in
  pseudobulk QC), but it substantially drives the USP34 tumor-cell
  effect specifically -- see section 5/9 -- and its per-sample values
  carry more noise than the other 5 samples generally.
- The all-epithelial "sensitivity" track (Phase 14) is not an
  independent dataset -- it reprocesses the same six sequencing libraries
  with a broader, overlapping cell selection, so direction agreement
  between it and the primary tumor-cell track should be read as
  robustness to the cell-population definition, not as independent
  replication.
- No PR/HER2 status, other systemic-treatment history, or explicit
  patient-pairing statement was available in the extracted metadata.
- No author-provided all-cell-type PT/RT object exists; the Phase 14
  "all-epithelial" sensitivity population, reconstructed from raw Cell
  Ranger matrices, turned out to be compositionally close to "all
  QC-passed cells" because the raw PT/RT samples themselves show
  minimal immune/stromal/endothelial content (0.08% PTPRC+) -- documented
  in `docs/GSE240112_CELLTYPE_DEFINITIONS.md`, not a code defect.
- USP17L29 remains untestable across every RNA layer in this project to
  date (CRISPR screen supplies its label; no RNA layer can evaluate it).
- scATAC deferred (`docs/GSE240112_SCATAC_DECISION.md`).

## 12. Codex verdict

Initial verdict: **FAIL** (6 issues found, 4 substantive -- inaccurate
USP34 sample-consistency claim, an overstated "independent populations"
claim, one factually false statement, and an overstated patient-pairing
claim; 2 stylistic/clarity). Full findings and corrections applied:
`docs/GSE240112_CODEX_REVIEW.md`. The underlying pseudobulk/edgeR/BH
implementation, sample mapping, HDF5 extraction, and USP17L29
detectability claim were all independently re-verified and found
correct -- every issue was in report interpretation/wording, not
computation. Corrections were applied directly to this report (sections
5, 6, 9, 11) and to `docs/GSE240112_DATA_AUDIT.md`; no source data,
extraction, or statistical computation was rerun.

## 13. Tests

48 GSE240112-specific tests passed
(`python3 -m pytest tests/ -k gse240112 -q`); 536 passed for the full
project test suite (`python3 -m pytest tests/ -q`), confirming no
regression to any existing frozen module.

## 14. Git

See final commit hash reported at the end of this run.

## 15. Additional secondary context (not part of the primary 14-section answer)

**Pathway context (Phase 17):** the authors' own precomputed per-cell
gene-module scores (already present in `TTs_cancer_060223.h5seurat`
metadata, not newly computed here) show a large PT-vs-RT separation in
the `TR` module (almost certainly "tamoxifen resistance"; PT mean=-1.49,
RT mean=+1.03 -- `results/tables/gse240112/pathway_module_scores_by_group.tsv`),
confirming the PT/RT group labels behave as biologically expected before
any candidate gene is considered. No new pathway-enrichment analysis was
run (out of scope for this run).

**Normal-tissue context (Phase 15, descriptive only, secondary):** for
the 12 testable candidates, mean log-normalized expression in NT
epithelial cells (bounded marker rule: >=2 of EPCAM/KRT8/KRT18/KRT19
positive and PTPRC-negative; 7,122/7,529 NT cells qualify) is reported in
`results/tables/gse240112/nt_epithelial_candidate_means.tsv`. No
statistical test is run (n=2 NT samples) and no tumor-specificity/safety
claim is made from it.
