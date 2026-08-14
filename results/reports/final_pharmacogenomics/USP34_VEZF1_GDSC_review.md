# USP34 / VEZF1 GDSC pharmacogenomic review

**Question.** Does baseline USP34 or VEZF1 expression associate with response to
any existing anticancer drug in breast-cancer cell lines?

**Scope.** This is a targeted lookup for the two frozen translational
candidates (USP34 = lead, VEZF1 = second/backup). It is **not** new candidate
discovery, does not touch the frozen four-gene shortlist, and does not rerun
any CRISPR/TCGA/DepMap analysis. All findings below are correlational.

---

## 1. Data provenance

| Field | Value |
|---|---|
| Release | GDSC 8.5 |
| Release date | 2023-10-30 |
| Source | `https://cog.sanger.ac.uk/cancerrxgene/GDSC_release8.5/` (the CancerRxGene bulk-download page returned HTTP 410; this direct Sanger COG storage URL was verified live) |
| Files | `GDSC1_fitted_dose_response_27Oct23.xlsx`, `GDSC2_fitted_dose_response_27Oct23.xlsx`, `screened_compounds_rel_8.5.csv`, `Cell_Lines_Details.xlsx` |
| Raw data location | OUTSIDE git (`/ibex/scratch/aljaroaa/tamoxifen-data/gdsc/`), SHA256-recorded in `PROVENANCE.txt` alongside the raw files, per this project's external-data convention |
| Expression source | DepMap Public 26Q1 -- **already downloaded and verified earlier in this project**, reused (not a new RNA-expression dataset) and joined to GDSC via the exact `SangerModelID`/`COSMIC_ID` identifiers present natively in both GDSC and DepMap's `Model.csv`. No fuzzy cell-line name matching was used. |

Full detail: `results/tables/final_pharmacogenomics/GDSC_data_provenance.tsv`.

### Response-metric direction (verified against the primary GDSC source, not assumed)

Verified directly from `GDSC_Fitted_Data_Description.pdf` (Sanger, v1.0.0, 21 Sep 2017):

- **LN_IC50** -- natural log of the fitted IC50 (uM). **Lower = more sensitive.**
- **AUC** -- fraction of the area under the fitted dose-response curve.
  **Lower = more sensitive** (the curve drops further across the screening range).

Both metrics point the same direction (lower = more drug-sensitive), so a
negative Spearman rho between gene expression and either metric means
"higher expression -> more sensitive."

### Cell-line join

`build_breast_expression_joined()`: GDSC1 = 333,161 rows (970 lines, 402
drugs); GDSC2 = 242,036 rows (969 lines, 295 drugs). Filtering to
`TCGA_DESC == "BRCA"` and joining to DepMap 26Q1 via `SangerModelID` gives
**51 GDSC breast lines -> 46 matched to DepMap 26Q1 expression** (5 lost to
no SangerModelID match, logged, not silently dropped). Of the 46 matched
lines, **14 are ER+/luminal** by this project's own established DepMap
classification rule (`independent_validation_depmap_data.load_model()`,
reused unchanged).

### Endocrine-compound availability (searched explicitly, both by name and by target)

| Compound | Present in GDSC 8.5? | Detail |
|---|---|---|
| Tamoxifen | **Yes** | DRUG_ID 1199, target ESR1, screened in both GDSC1 and GDSC2, n=51 breast lines total |
| Fulvestrant | **Yes** | 3 separate DRUG_IDs (1200 GDSC2, 1414 GDSC1, 1816 GDSC2), target ESR |
| 4-hydroxytamoxifen (4-OHT) | **No** | not found under any name/target search |
| Endoxifen | **No** | not found under any name/target search |

Full detail: `results/tables/final_pharmacogenomics/GDSC_compound_availability.tsv`.

---

## 2. Statistical method

Spearman correlation between DepMap 26Q1 log2(TPM+1) expression and each GDSC
response metric, computed **separately per (drug_id, dataset, metric)** cell
(i.e. GDSC1 and GDSC2 treated as materially different screening campaigns,
never pooled, and grouped by the exact `DRUG_ID` rather than `DRUG_NAME` --
several compounds, e.g. AZD7762, were re-screened under more than one
`DRUG_ID` even within a single GDSC release, and grouping by name alone would
pseudoreplicate). Benjamini-Hochberg FDR correction applied within each
(dataset, metric) family. Minimum N = 15 lines for the main analysis; the
ER+/luminal subset (N=14) is reported separately and explicitly flagged
exploratory (never treated as independent validation).

A direct duplicate check on (SangerModelID, DRUG_ID, DATASET) in the joined
breast dataset found zero duplicate groups (30,394/30,394 groups with exactly
one row), confirming the main analysis has no pseudoreplication.

Full per-drug results: `results/tables/final_pharmacogenomics/USP34_GDSC_drug_associations.tsv`,
`results/tables/final_pharmacogenomics/VEZF1_GDSC_drug_associations.tsv`
(1,278 rows each: every screened drug x dataset x metric combination with
N>=15).

---

## 3. USP34 results

**9 of 1,278 drug-tests reach FDR<0.05**, all in GDSC1, all the same
direction (higher USP34 expression -> more sensitive, i.e. lower
LN_IC50/AUC):

| Drug | Target | Pathway | Metric(s) | FDR |
|---|---|---|---|---|
| AZD7762 | CHEK1, CHEK2 | Cell cycle | LN_IC50 + AUC | 0.0078-0.0084 |
| AZD1332 | NTRK1-3 | RTK signaling | LN_IC50 | 0.0248 |
| AZD6738 | ATR | Genome integrity | LN_IC50 | 0.0248 |
| JAK1_3715 | JAK1 | Other kinases | LN_IC50 | 0.0248 |
| Sphingosine Kinase 1 Inhibitor II | SPHK1 | Other kinases | LN_IC50 | 0.0248 |
| FGFR_3831 | FGFR1-4 | RTK signaling | LN_IC50 | 0.0248 |
| LDN-193189 | BMP | Other | LN_IC50 | 0.0439 |
| Ara-G | Anti-metabolite | Other | LN_IC50 | 0.0458 |

Full table: `results/tables/final_pharmacogenomics/GDSC_top_associations.tsv`;
figure: `results/figures/final_pharmacogenomics/01_USP34_GDSC_drug_response.png`.

**Strongest hit -- AZD7762 (CHEK1/CHEK2).** rho=-0.59 (AUC, FDR=0.0078),
rho=-0.59 (LN_IC50, FDR=0.0084), n=44, GDSC1, DRUG_ID 1402. AZD7762 was
re-screened under a second, distinct DRUG_ID within GDSC1 (1022); that batch
is directionally consistent but **not** FDR-significant (FDR=0.174-0.350) --
only DRUG_ID 1402 reaches FDR<0.05, so this is not "two independently
significant re-screens." GDSC2 (a different DRUG_ID again) shows a nominally
significant, same-direction correlation (rho=-0.37, p=0.011) that does not
survive its own within-screen FDR correction -- partial, not independent,
replication. In the small (N=14) ER+/luminal subset it is markedly stronger
(rho=-0.81, p=0.0008, n=13, using the exact DRUG_ID 1402), but this subset
result is exploratory only, not a validation.

CHEK1 itself is directionally consistent (sensitising_KO) but **not**
FDR-significant in this project's own frozen Hany CRISPR screen (FDR=0.812).
It **is** nominally FDR-significant in this project's frozen GSE118713
transcriptomic tamoxifen-resistance dataset (FDR=0.0497, down-in-TAMR
direction) -- a genuine cross-dataset echo, generated independently of this
GDSC lookup. See Section 7 below for the full crosscheck, which also finds
several OTHER drug-target genes from this section with their own
FDR-significant signals in GSE118713/GSE240112 -- CHEK1 is not unique in
that respect, and Section 7 should be read alongside this section, not in
isolation.

**LDN-193189 (BMP).** FDR=0.0439, coherent with USP34's own independently
established BMP2/Smad1/RUNX2 osteogenesis mechanism (Guo et al. 2018, PMID
30181118, carried forward from the druggability_safety phase) -- but that
mechanism was demonstrated in bone/mesenchymal biology, not breast cancer,
and this correlation alone does not show USP34 regulates BMP signaling in
breast cancer cells. It does **not** replicate in the ER+/luminal subset
(p=0.37) and should not be over-weighted.

**The remaining 6 hits** (AZD1332/NTRK1-3, AZD6738/ATR, JAK1_3715/JAK1,
Sphingosine Kinase 1 Inhibitor II/SPHK1, FGFR_3831/FGFR1-4, Ara-G) have no
known direct, upstream-regulatory, or pathway-level mechanistic connection to
USP34 biology identified anywhere in this project's prior work. They are
reported as real, FDR-significant correlations with no current mechanistic
explanation. As with CHEK1, several of these targets' own genes independently
show FDR-significant signal in the project's frozen resistance datasets
(FGFR2, FGFR4, JAK1 in GSE118713; FGFR2, FGFR3, FGFR4 in GSE240112) -- see
Section 7.

**Estrogen signaling / tamoxifen / fulvestrant.** No FDR-significant
association for USP34 with tamoxifen: GDSC1 AUC is nominally significant on
its own (rho=+0.35, p=0.023) but does not survive FDR correction (FDR=0.318,
n=43); GDSC1 LN_IC50 (p=0.55), GDSC2 LN_IC50 (p=0.52), and GDSC2 AUC (p=0.53)
show no signal at all. No fulvestrant DRUG_ID reaches significance either.
This is consistent with, not contradictory to, the Hany finding: GDSC tests
baseline expression-vs-response correlation across a cross-sectional
cell-line panel, a different question from Hany's CRISPR-perturbation-plus-
tamoxifen design.

**WNT/beta-catenin.** No FDR-significant association for USP34 with any
WNT-pathway compound (tankyrase inhibitors XAV939/MN-64/WIKI4/AZ6102,
mechanistically adjacent to USP34's own AXIN1 axis); best nominal p=0.0076
(TWS119/GSK3), FDR=0.138.

**Hippo/TEAD.** GDSC Release 8.5 contains **zero** TEAD/YAP/Hippo-pathway
compounds under any name or target search -- this hypothesis cannot be tested
pharmacogenomically in GDSC at all (an absence of chemical matter, not a
negative finding).

---

## 4. VEZF1 results

**Zero of 1,278 drug-tests reach FDR<0.05.** Best nominal signal: Paclitaxel
(microtubule stabiliser, mitosis pathway), p=0.0005, FDR=0.132, GDSC2,
n=45 -- does not survive correction.

Full table: `results/tables/final_pharmacogenomics/GDSC_top_associations.tsv`;
figure: `results/figures/final_pharmacogenomics/02_VEZF1_GDSC_drug_response.png`.

No cross-screen consistency check is meaningful with zero FDR-significant
hits. No association with tamoxifen (p=0.29-0.88 across both screens/metrics)
or fulvestrant. As with USP34, **zero TEAD/Hippo compounds exist in GDSC**,
so the VEZF1-TEAD1 hypothesis cannot be tested here either -- this remains an
absence of chemical matter, not evidence against the hypothesis.

This is a genuine negative result, reported as such.

---

## 5. ER+/luminal subset (N=14, exploratory only)

Re-tested using the exact `DRUG_ID` (never `DRUG_NAME` alone, to avoid the
pseudoreplication trap described in Section 2) for every hit carried over
from the full-breast-line analysis. All results here are explicitly tagged
`exploratory_only=True` and are **not independently FDR-corrected within
this subset** -- they are not treated as validation of the full-set findings.

- AZD7762 (CHEK1/CHEK2): strengthens (rho=-0.81, p=0.0008, LN_IC50, n=13).
- LDN-193189 (BMP): does not replicate (rho=-0.26, p=0.37).
- All other USP34 hits and the top VEZF1 nominal hits: weaker or
  non-significant at this N.

Full table: `results/tables/final_pharmacogenomics/GDSC_ER_luminal_subset.tsv`.

---

## 6. Indirect-targeting classification

Every hit above was classified into one of six standard categories
(DIRECT_TARGET / KNOWN_UPSTREAM_REGULATOR / KNOWN_REQUIRED_PARTNER /
PATHWAY_LEVEL_CONNECTION / PHARMACOGENOMIC_ASSOCIATION_ONLY /
NO_KNOWN_CONNECTION):

- AZD7762 (CHEK1/CHEK2) -> **PATHWAY_LEVEL_CONNECTION**
- LDN-193189 (BMP) -> **PATHWAY_LEVEL_CONNECTION**
- The remaining 6 USP34 hits (FGFR, NTRK, ATR, JAK1, SPHK1, Ara-G) -> **PHARMACOGENOMIC_ASSOCIATION_ONLY**
- VEZF1 (no FDR-significant hits) -> **NO_KNOWN_CONNECTION**
- VEZF1/TEAD/Hippo -> **PHARMACOGENOMIC_ASSOCIATION_ONLY** (no compound exists
  to test it; TEAD1 remains **not a validated VEZF1 indirect target**, unchanged
  from the lead_target_deep_dive phase)

No category above is DIRECT_TARGET, KNOWN_UPSTREAM_REGULATOR, or
KNOWN_REQUIRED_PARTNER for either gene. Full table:
`results/tables/final_pharmacogenomics/GDSC_indirect_targeting_classification.tsv`.

---

## 7. Cross-check against frozen project evidence

Every drug-target gene from the hits above (CHEK1, CHEK2, ATR, FGFR1-4,
NTRK1-3, JAK1, SPHK1, NAMPT, ACACA, HDAC1/2) was looked up against this
project's existing frozen cross-dataset table
(`all_genes_cross_dataset_evidence_with_ranking.tsv`), without altering it.
This lookup checks whether the drug-TARGET gene itself (not USP34) shows any
independent signal in the project's own frozen resistance datasets -- it is
a check for coincidental convergence across independent datasets, not a
mechanistic link back to USP34.

None of these target genes reaches FDR<0.05 in the project's own frozen Hany
CRISPR screen. Several DO reach FDR<0.05 in the frozen GSE118713 and/or
GSE240112 transcriptomic resistance datasets:

| Target gene | GSE118713 FDR | GSE118713 direction | GSE240112 FDR | GSE240112 direction |
|---|---|---|---|---|
| CHEK2 | 8.0e-06 | down_in_TAMR | 0.353 (n.s.) | down |
| FGFR4 | 0.0001 | down_in_TAMR | 0.0003 | down |
| FGFR2 | 0.0145 | down_in_TAMR | 0.0114 | down |
| JAK1 | 0.0067 | down_in_TAMR | 0.407 (n.s.) | up |
| NAMPT | 0.0195 | up_in_TAMR | 0.210 (n.s.) | down |
| CHEK1 | 0.0497 | down_in_TAMR | 0.806 (n.s.) | down |
| FGFR3 | 0.289 (n.s.) | up_in_TAMR | 0.0040 | up |

FGFR2 and FGFR4 are the most internally consistent: both reach FDR<0.05 in
BOTH GSE118713 and GSE240112, with the SAME direction (down in resistant
tissue) in both. CHEK2, JAK1, and NAMPT reach FDR<0.05 in GSE118713 only,
with GSE240112 not significant (and, for JAK1 and NAMPT, in the opposite
direction). FGFR3 reaches FDR<0.05 in GSE240112 only.

**These are direction/significance patterns for the target genes' OWN
expression in independent resistance datasets -- they do not, by themselves,
say anything about whether USP34 regulates these genes or these pathways.**
They are reported here as an honest inventory of coincidental cross-dataset
signal, not as validation of the GDSC associations, and not as evidence for
any USP34 mechanism. Full table with all target genes and both directions:
`results/tables/final_pharmacogenomics/GDSC_project_crosscheck.tsv`.

**This lookup does not alter the frozen four-gene shortlist or the frozen
Hany/TCGA/DepMap values in any way.**

---

## 8. Final interpretation

Answered explicitly and separately for each gene (full table:
`results/tables/final_pharmacogenomics/GDSC_final_interpretation.tsv`):

**USP34**
- Does expression associate with any existing drug response? **Yes** -- 9 of
  1,278 drug-tests, FDR<0.05, all GDSC1, all same direction.
- Most biologically coherent pathway? CHEK1/CHEK2 (AZD7762), with an
  independent nominal GSE118713 echo (FDR=0.0497, down-in-TAMR; CHEK1's own
  Hany CRISPR FDR is 0.812, not significant) -- though several other hit
  targets (FGFR2/FGFR4 most consistently) also show their own FDR-significant
  signal in GSE118713 and/or GSE240112 (Section 7), so this is not a uniquely
  strong echo; BMP (LDN-193189), coherent with prior osteogenesis mechanism
  but non-replicating in the ER+/luminal subset.
- Any realistic hypothesis worth future testing? Yes, as a **new, separate**
  hypothesis: whether USP34 perturbation sensitises resistant ER+ cells to
  CHK1/CHK2 inhibition. This is not incorporated into EXP-1's design in this
  phase and is not evidence about the tamoxifen-sensitisation hypothesis
  itself.
- **Does GDSC alter USP34 = LEAD? No.** Tamoxifen/fulvestrant themselves show
  no significant USP34 association in GDSC; nothing here strengthens or
  weakens the frozen lead conclusion.

**VEZF1**
- Does expression associate with any existing drug response? **No**
  FDR-significant hits out of 1,278 drug-tests.
- Any Hippo/TEAD relevance? Cannot be tested -- GDSC contains zero
  TEAD/YAP/Hippo compounds. TEAD1 remains unvalidated as a VEZF1 indirect
  target.
- **Does GDSC alter VEZF1 = BACKUP? No.** This adds a further negative data
  point consistent with VEZF1's already-documented poor druggability from the
  lead_target_deep_dive phase; the frozen SECOND/BACKUP conclusion is
  unchanged.

---

## 9. Causality discipline

All findings above are reported strictly as "X expression is associated with
drug Y response in GDSC breast cancer lines." **No causal claim, no direct
mechanistic claim beyond documented target annotation, and no treatment or
combination recommendation is made anywhere in this report.** The CHK1/CHK2
observation is noted as a hypothesis worth *future* investigation, explicitly
not as a validated finding or recommendation.

---

## 10. Figures

- `01_USP34_GDSC_drug_response.png` -- the 9 FDR-significant USP34 hits plus
  context, Spearman rho with FDR annotation.
- `02_VEZF1_GDSC_drug_response.png` -- top 9 VEZF1 associations by nominal
  p-value, all non-significant (honest negative result).
- `03_USP34_VEZF1_pharmacogenomic_summary.png` -- FDR-significant hit counts,
  endocrine-compound availability, and endocrine/Hippo negative findings,
  side by side for both genes.
