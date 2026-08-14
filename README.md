# Breast Cancer Tamoxifen Resistance Project

Computational reanalysis identifying, validating, and deeply characterizing
genes whose loss sensitizes ER-positive breast cancer to tamoxifen, by
integrating a public genome-wide CRISPR screen with public bulk and
single-cell transcriptomics, TCGA clinical data, DepMap dependency data,
network/pathway analysis, structural biology, and pharmacogenomic
(GDSC) drug-response data.

**This is a computational reanalysis of public data. No wet-lab work was
performed by this project, and no result here should be read as a
validated therapeutic target, a proven mechanism, or a safe/efficacious
drug candidate.**

**Quick links:** [Project workflow](docs/PROJECT_WORKFLOW.md) ·
[Data provenance](docs/DATA_PROVENANCE.md) ·
[Results guide](docs/RESULTS_GUIDE.md) ·
[Reproducibility](docs/REPRODUCIBILITY.md) ·
[Canonical results index](results/final/README.md) ·
[Root pre-analysis plan](PREANALYSIS.md) ·
[Code map](docs/CODE_MAP.md) ·
[Tests](tests/)

---

## A. Overview

This project starts from a public genome-wide CRISPR knockout screen in
tamoxifen-treated ER+ breast cancer cells and, over sixteen analysis phases,
narrows a genome-wide gene list down to a small, evidence-ranked shortlist
of tamoxifen-resistance candidate genes, checks that shortlist against
independent public cohorts (TCGA, DepMap, additional GEO datasets), maps the
mechanistic/pathway context of the top candidates, assesses their structural
druggability and safety liabilities, and designs (but does not run) a
concrete wet-lab validation plan for the two leading candidates. All data
are public; all thresholds were declared in advance in dated pre-analysis
plans; nothing derived from the CRISPR screen is used as a model feature
(the screen supplies labels only).

## B. Main research question

**Which genes, when lost, sensitize tamoxifen-resistant or tamoxifen-treated
ER-positive breast cancer cells to tamoxifen — and, among the strongest
candidates, which is the most promising target for a concrete follow-up
combination-therapy experiment?**

## C. High-level final findings (USP34 and VEZF1)

The project's frozen four-gene therapeutic shortlist, in order, is
**USP34 > VEZF1 > EML5 > CITED2** (see
[`docs/THERAPEUTIC_SHORTLIST_FREEZE.md`](docs/THERAPEUTIC_SHORTLIST_FREEZE.md)).
USP34 and VEZF1 were carried forward into deep mechanistic, structural, and
translational-design work; EML5 and CITED2 were not.

**USP34 — current lead translational target.** Functional CRISPR
sensitisation (Hany screen, FDR=0.042) + low baseline dependency in
ER+/luminal cancer cell lines (DepMap 26Q1: 0.0% strongly dependent) + real,
crystallographically-confirmed catalytic targetability (reactive Cys1903,
PDB 7W3R/7W3U) together make USP34 the lead combination-target *hypothesis*
— not a validated target. No USP34-selective inhibitor currently exists.
Structural pocket analysis (fpocket 4.2.3) found a druggable catalytic
pocket but concluded blind/arbitrary small-molecule docking is **not yet
justified** (`DOCKING_NOT_YET_JUSTIFIED`; see
[Section M](#m-proposed-next-experiment) and
[`docs/RESULTS_GUIDE.md`](docs/RESULTS_GUIDE.md)). A targeted GDSC
Release 8.5 pharmacogenomic lookup found 9 FDR-significant drug-response
correlations for USP34 expression in breast cancer cell lines (strongest:
AZD7762/CHK1-CHK2 inhibitor, FDR=0.008) — all correlational, none causal,
and none involving tamoxifen or fulvestrant directly. A prior mammary
epithelial study (PMID 28499884) found USP34 loss can promote EMT/stem-like
features in some contexts; **this does not invalidate the USP34 tamoxifen
hypothesis** — it is counter-evidence motivating explicit EMT/stemness
monitoring in any USP34-perturbation experiment, not a reason to drop USP34.

**VEZF1 — second / backup translational target.** Strong functional CRISPR
sensitisation (Hany screen, FDR=0.037, nominally the stronger of the two) +
real baseline dependency in ER+/luminal cancer cell lines (DepMap 26Q1:
27.3% strongly dependent) together suggest a dual-action biological
hypothesis, but VEZF1 has poor direct druggability (a zinc-finger
transcription factor). A candidate indirect strategy — inhibiting VEZF1's
proposed downstream partner TEAD1 — remains **explicitly unvalidated**
(`PHARMACOGENOMIC_ASSOCIATION_ONLY` at best; GDSC contains zero TEAD/Hippo
compounds, so this cannot even be tested pharmacogenomically). The GDSC
lookup found **zero** FDR-significant drug-response correlations for VEZF1
in breast cancer cell lines — a genuine negative result. VEZF1 also has a
more directly causal cardiovascular/developmental liability signal (a
postnatal zebrafish cardiac-contractility finding, PMID 31911272) than
USP34's bone-related liability signal.

Full detail: [`results/final/README.md`](results/final/README.md).

## D. Complete analysis workflow (16 phases)

```
 1. CRISPR screen reanalysis (Hany et al. 2023)              -- functional labels, Gate-1 hit calling
 2. GSE118713 bulk RNA resistance analysis                   -- TAMR/FASR vs MCF7 expression
 3. CRISPR x bulk-RNA integration                            -- 13 sensitising candidates prioritized
 4. NEBULA poster figures                                    -- early-phase summary figures
 5. GSE245601 single-cell preprocessing + malignant-cell ID  -- InferCNV (primary) + CopyKAT (sensitivity check)
 6. GSE245601 candidate-level expression / pseudobulk        -- per-candidate acute-treatment signal
 7. GSE240112 primary-vs-recurrent scRNA-seq                 -- 4th independent evidence layer
 8. GSE111151 independent resistance-model validation        -- 5th independent evidence layer
 9. Unbiased genome-wide cross-dataset integration            -- all five evidence layers, every gene
10. Candidate adjudication                                    -- 7 MULTIMODAL_STRONG genes narrowed
11. Evidence freeze                                            -- frozen 4-gene shortlist: USP34>VEZF1>EML5>CITED2
12. Systems / pathway network mapping                          -- shortest paths to resistance nodes
13. Literature mechanism review                                 -- published mechanism check per candidate
14. Independent validation (TCGA-BRCA + DepMap 26Q1)            -- external cohorts, not used for discovery
15. Lead-target deep dive + druggability/safety                 -- structure, pocket, safety liabilities
16. Final translational design + GDSC pharmacogenomics          -- EXP-1..5 wet-lab plan, drug-response lookup
```

Steps 1-4 used only the CRISPR screen and GSE118713; steps 5-9 added
independent transcriptomic evidence layers genome-wide; steps 10-11 froze
the shortlist; steps 12-16 characterize the frozen shortlist's top two
candidates without ever reopening candidate discovery or altering the
frozen ranking. See [`docs/PROJECT_WORKFLOW.md`](docs/PROJECT_WORKFLOW.md)
for the full narrative version of this workflow, with code/output pointers
for every phase.

## E. Datasets

| Dataset | Role | Biological system | Caveats |
|---|---|---|---|
| Hany et al. 2023 CRISPR screen (*Sci Adv* 9:eadd3685) | Functional labels only (never used as a feature) | MCF7-V drug-tolerant parental clone, E2 vs E2+4-OHT | Genome-wide but single cell line/screen; Gate-1 FDR<0.1 |
| GSE118713 | Resistance-associated bulk expression | MCF7 parental / TAMR / FASR cell lines | Cell-line model, not patient tissue |
| GSE245601 | Human ex vivo tamoxifen-response single-cell context | 10 paired primary ER+/HER2- tumors, 10 µM tamoxifen vs control, 12h ex vivo | **Only 3/10 tumor pairs (Tumor_02/03/07) meet the frozen >=50-malignant-cell pseudobulk-eligibility rule**; InferCNV/CopyKAT malignant-cell concordance is highly variable across samples (~56% average agreement) |
| GSE240112 | Primary-vs-recurrent scRNA-seq, 4th evidence layer | Matched primary/recurrent ER+ tumor pairs | Recurrence != tamoxifen resistance specifically; small patient N |
| GSE111151 | Independent post-hoc resistance-model confirmation | 4 parental + 7 TamR derivative cell lines | Used for confirmation only, never for candidate discovery |
| TCGA-BRCA | Independent large-cohort expression/clinical validation | 1,095 primary tumors (deduplicated), bulk RNA-seq + clinical | **Not a tamoxifen-resistance cohort** — a general BRCA cohort; ER+/ER- and treatment status come from clinical annotation, not a resistance phenotype |
| DepMap Public (24Q4 archived; 26Q1 active) | Cancer cell-line dependency (CRISPR) | Pan-cancer + breast + ER+/luminal cell-line panels | Two releases used; 26Q1 is now the active/reported release, 24Q4 is archived for traceability — see [`docs/DATA_PROVENANCE.md`](docs/DATA_PROVENANCE.md) |
| GDSC Release 8.5 (30 Oct 2023) | Pharmacogenomic drug-response correlation, USP34/VEZF1 only | Breast cancer cell line panel, GDSC1 + GDSC2 screens | **Correlational only, never causal**; GDSC1/GDSC2 are separate screening campaigns, never pooled for FDR; some compounds carry multiple DRUG_IDs within one release (pseudoreplication trap, handled explicitly) |

## F. The frozen candidates — three distinct rankings, not one

This project has **three separate rankings**, each answering a different
question, and none of them overrides another:

1. **Frozen therapeutic shortlist ranking** (evidence_freeze phase,
   [`docs/THERAPEUTIC_SHORTLIST_FREEZE.md`](docs/THERAPEUTIC_SHORTLIST_FREEZE.md)):
   USP34 > VEZF1 > EML5 > CITED2, based on combined CRISPR + bulk-RNA +
   single-cell + independent-validation evidence strength. **This ranking is
   frozen and is never altered by any later phase.**
2. **Mechanistic/pathway follow-up order** (systems_network +
   literature_mechanism + independent_validation phases,
   [`four_candidate_followup_rankings.tsv`](results/tables/independent_validation/four_candidate_followup_rankings.tsv)): a separate, explicitly-labeled
   "follow-up order only" signal from TCGA pathway/clinical cross-checks —
   it does **not** alter the frozen therapeutic ranking above, by design.
3. **Final translational lead selection** (final_translational phase,
   [`results/final/README.md`](results/final/README.md)): among the top two
   frozen candidates only (USP34, VEZF1 — EML5 and CITED2 were not carried
   this far), USP34 is designated the current **lead** and VEZF1 the
   **second/backup** target specifically for wet-lab combination-therapy
   follow-up, based on structural/druggability/safety criteria layered on
   top of (not replacing) the frozen ranking. This happens to agree with
   the frozen order here, but is a genuinely separate question (translational
   feasibility, not overall evidence strength).

## G. Repository structure

```
config/            YAML configuration -- every file path used by analysis code lives here
data/
  raw/             External data (symlinked to Ibex scratch storage) -- NOT in git
  processed/       A handful of frozen, force-added intermediate matrices -- in git
  reference/       Gene sets, STRING interactions, TCGA Ensembl ID mapping -- in git
  checksums/       SHA256 manifests for raw downloads
docs/               Narrative documentation: workflow, provenance, results guide, per-phase audits
external_refs/      Cloned reference tool source (InferCNV, CopyKAT) for method audit -- NOT in git
results/
  tables/           One subdirectory per analysis phase, TSV outputs
  figures/          One subdirectory per analysis phase, PNG/PDF outputs
  reports/          One subdirectory per analysis phase, narrative .md reports
  networks/         Cytoscape-importable network files
  final/            Canonical results index -- START HERE for final findings
scripts/
  download/         One-time, deterministic network-download scripts (never run at analysis time)
  analysis/         R scripts (limma, edgeR, Seurat, InferCNV, CopyKAT) invoked from src/ via subprocess
src/                Python analysis modules -- one (or a small family) per phase, deterministic, no network calls
tests/              One pytest module per src/ module, exercising real logic
PREANALYSIS.md       Root pre-analysis plan (CRISPR/bulk phase), with a dated, append-only amendments log
CLAUDE.md            Project-specific engineering rules (data hygiene, hard rules, commit conventions)
```

## H. Key results — start here

- **Final findings for USP34/VEZF1:** [`results/final/README.md`](results/final/README.md)
- **Frozen 4-gene shortlist:** [`docs/THERAPEUTIC_SHORTLIST_FREEZE.md`](docs/THERAPEUTIC_SHORTLIST_FREEZE.md)
- **GDSC drug-response review:** [`results/reports/final_pharmacogenomics/USP34_VEZF1_GDSC_review.md`](results/reports/final_pharmacogenomics/USP34_VEZF1_GDSC_review.md)
- **Translational experimental plan (EXP-1..5):** [`results/reports/final_translational/final_USP34_VEZF1_translational_plan.md`](results/reports/final_translational/final_USP34_VEZF1_translational_plan.md)
- **Independent TCGA/DepMap validation:** [`results/reports/independent_validation/four_candidate_TCGA_DepMap_review.md`](results/reports/independent_validation/four_candidate_TCGA_DepMap_review.md)
- **Full results guide (what every table/figure means and when it's current):** [`docs/RESULTS_GUIDE.md`](docs/RESULTS_GUIDE.md)

## I. Reproducibility

- Analyses ran on KAUST Ibex (SLURM); raw sequencing/screen/TCGA/DepMap/GDSC
  data are **not** committed to this repository (`data/raw` is a symlink to
  external scratch storage, gitignored). Public data are downloaded from
  source by deterministic scripts under [`scripts/download/`](scripts/download/).
- All file paths come from [`config/config.yaml`](config/config.yaml); no
  path is hardcoded in analysis code.
- Every `src/` module has a corresponding pytest module under
  [`tests/`](tests/) exercising its actual logic (recomputing statistics,
  checking real numeric values), not merely that it runs without error.
  **Current test suite: 1,150 passed, 1 skipped** (run `pytest -q`).
- Thresholds are declared in advance in [`PREANALYSIS.md`](PREANALYSIS.md)
  and phase-specific pre-analysis plans under `docs/`, with dated,
  append-only amendment logs — never edited in place after analysis begins.
- Nothing derived from the CRISPR screen enters the feature table; the
  screen supplies labels only (see [`CLAUDE.md`](CLAUDE.md)).
- See [`docs/REPRODUCIBILITY.md`](docs/REPRODUCIBILITY.md) for the full,
  step-by-step reproduction guide, including known environment gaps.

## J. Installation / quick start

```bash
# Create the pinned conda/micromamba environment
micromamba env create -f environment.yml -n bc

# Run the full test suite
micromamba run -n bc python3 -m pytest -q

# Run any single analysis module (paths always come from config/config.yaml)
micromamba run -n bc python3 -m src.final_pharmacogenomics_build_tables
micromamba run -n bc python3 -m src.final_pharmacogenomics_visualization
```

`environment.yml` pins the Python/pandas/scipy/statsmodels stack and the R
packages (edgeR, limma) used directly from Python via subprocess. The
single-cell R pipeline (Seurat, InferCNV, CopyKAT) was run in a separate R
environment on Ibex that is **not** captured in `environment.yml` — see the
known-gaps note in [`docs/REPRODUCIBILITY.md`](docs/REPRODUCIBILITY.md).
Raw data downloads require network access and are not needed to run the
test suite, which operates on already-downloaded/committed tables.

## K. Data availability

All source data are public. Raw files themselves are not committed to this
repository (see [`docs/DATA_PROVENANCE.md`](docs/DATA_PROVENANCE.md) for
exact accessions, release dates, and, where computable, SHA256 checksums):

- Hany et al. 2023 CRISPR screen: *Sci Adv* 9:eadd3685, Data S1
- GSE118713, GSE245601, GSE240112, GSE111151: NCBI GEO
- TCGA-BRCA: GDC (Genomic Data Commons)
- DepMap Public 24Q4 / 26Q1: depmap.org
- GDSC Release 8.5: `cog.sanger.ac.uk/cancerrxgene/GDSC_release8.5/`

Derived tables/figures/reports produced by this project's own code are
committed to this repository under `results/`.

## L. Limitations

- This is a **computational reanalysis of public data only** — no wet-lab
  validation of any finding in this repository has been performed.
- GSE245601 malignant-cell yield varies strongly by tumor; only 3 of 10
  paired tumors meet the frozen pseudobulk-eligibility threshold, and
  InferCNV/CopyKAT malignant-cell-classification concordance is highly
  variable between samples (from near-0% to near-100%, ~56% on average).
- TCGA-BRCA is a general breast-cancer cohort, **not** a tamoxifen-resistance
  cohort; its ER+/ER- and treatment annotations come from clinical records,
  not a resistance phenotype.
- GSE111151 and GSE240112 provide independent but imperfect proxies for
  tamoxifen resistance (established resistant cell-line derivatives; primary
  vs. recurrent tumor pairs, respectively) — neither is a direct replication
  of the Hany screen's specific perturbation.
- DepMap dependency data span two releases (24Q4 archived, 26Q1 active);
  results are cross-checked for reproduction across releases, but the
  underlying cell-line panel composition can shift between releases.
- All GDSC pharmacogenomic findings are **correlational only** — baseline
  expression correlated with drug response across a cell-line panel is not
  evidence of a causal or mechanistic relationship, and none of the GDSC
  findings involve tamoxifen or fulvestrant directly for either USP34 or
  VEZF1.
- No validated USP34-selective small-molecule inhibitor currently exists;
  structural analysis concluded docking is not yet justified
  (`DOCKING_NOT_YET_JUSTIFIED`) pending a real fragment-screening campaign.
- USP34 is the most human-genetically-constrained candidate examined in this
  project (LOEUF=0.152), a caution about tolerability, not a toxicity
  prediction.
- Prior mammary-epithelial literature (PMID 28499884) shows USP34 loss can
  promote EMT/stem-like features in some contexts — counter-evidence
  motivating explicit EMT/stemness monitoring in any USP34 experiment, not
  a reason to exclude USP34.
- VEZF1 has poor direct druggability (zinc-finger transcription factor); its
  proposed indirect target TEAD1 remains **unvalidated**, and GDSC contains
  zero TEAD/Hippo-pathway compounds, so this hypothesis cannot currently be
  tested pharmacogenomically at all.
- VEZF1 has a directly causal cardiovascular/developmental liability signal
  (postnatal zebrafish cardiac-contractility finding, PMID 31911272).
- `environment.yml` does not capture the R single-cell pipeline dependencies
  (Seurat, InferCNV, CopyKAT), which were run in a separate Ibex R
  environment — see [`docs/REPRODUCIBILITY.md`](docs/REPRODUCIBILITY.md).
- The genome-wide cross-dataset integration table anonymizes gene identity
  for part of its ranking-stability analysis; see
  `results/tables/cross_dataset_genomewide/anonymized_gene_mapping.tsv` for
  the mapping and `docs/` for the rationale.
- Sample sizes are small at several points in the pipeline (e.g. 14 ER+/
  luminal GDSC breast lines, 3/10 GSE245601 tumor pairs); results from these
  small-N subsets are always explicitly flagged exploratory, never treated
  as validation.
- No result in this repository should be read as "validated therapeutic
  target," "restores tamoxifen sensitivity," or "safe target."

## M. Proposed next experiment

The final translational phase designed (but did not run) five concrete
wet-lab experiments for the two lead candidates, in an acquired
tamoxifen-resistant ER+ breast cancer cell line:

- **EXP-1 (USP34, primary):** genetic (not pharmacological, since no
  selective inhibitor exists) USP34 perturbation + tamoxifen re-challenge,
  with a full dose-response/Bliss-independence interaction framework (not a
  single-dose comparison) and mandatory EMT/stemness monitoring (CDH1,
  CDH2, SNAI1, AXIN1, active beta-catenin) given the PMID 28499884
  counter-evidence.
- **EXP-2A (USP34 comparator A — lineage/cancer-selectivity):** the same
  perturbation in normal human mammary epithelial cells (e.g. MCF10A), to
  test cancer-selectivity of any sensitising effect.
- **EXP-2B (USP34 comparator B — known liability):** the same perturbation
  in primary human MSCs undergoing osteogenic differentiation, given
  USP34's established BMP2/Smad1/RUNX2 bone-biology role.
- **EXP-3 (VEZF1, secondary):** analogous genetic perturbation + tamoxifen
  re-challenge for VEZF1 in the same resistant line.
- **EXP-4 (VEZF1 comparator — cardiovascular liability):** the same
  perturbation in primary human vascular endothelial cells.
- **EXP-5 (VEZF1-TEAD1 hypothesis test, not a therapeutic experiment):** a
  targeted test of whether VEZF1 acts through TEAD1, designed so that it
  could reject (not just confirm) that indirect-targeting hypothesis before
  it is used for anything therapeutic.

Full design, readouts, outcome categories, and success/failure criteria:
[`results/reports/final_translational/final_USP34_VEZF1_translational_plan.md`](results/reports/final_translational/final_USP34_VEZF1_translational_plan.md).

## N. Citation / acknowledgement

This repository is an independent computational reanalysis project by
Almokhtar Aljarodi (KAUST). It reuses public data from Hany et al. 2023
(*Sci Adv*), NCBI GEO (GSE118713, GSE245601, GSE240112, GSE111151), TCGA
(GDC), DepMap, and GDSC/CancerRxGene — see
[`docs/DATA_PROVENANCE.md`](docs/DATA_PROVENANCE.md) for exact citations of
each source dataset. No formal citation for this repository itself has been
established; if referencing this work, please link to this repository
directly. Licensed under the MIT License (see [`LICENSE`](LICENSE)).
