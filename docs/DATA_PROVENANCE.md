# Data Provenance

Exact source, accession, release version, and storage location for every
external dataset used in this project. Raw data are never committed to
git; only derived tables/figures/code are. See
[`../.gitignore`](../.gitignore) and [`../config/config.yaml`](../config/config.yaml)
for the concrete enforcement of this rule.

All raw data live under `data/raw/`, which is a **symlink** to
`/ibex/scratch/aljaroaa/tamoxifen-data/` (KAUST Ibex scratch storage) and is
gitignored in full.

---

## Hany et al. 2023 CRISPR screen

- **Source:** Hany, D. et al. 2023. *Science Advances* 9, eadd3685.
- **File used:** Data S1 (guide-level counts).
- **Config key:** `data.raw.hany_data_s1`.
- **Role:** functional perturbation labels only (never a model feature).

## GSE118713 (bulk RNA-seq)

- **Source:** NCBI GEO accession GSE118713.
- **System:** MCF7 parental / TAMR / FASR cell lines.
- **Config key:** `data.raw.gse118713_dir`.

## GSE245601 (single-cell RNA-seq)

- **Source:** NCBI GEO accession GSE245601.
- **System:** 10 paired primary ER+/HER2- tumors, ex vivo 10 µM tamoxifen
  vs. control media, 12h.
- **Download:** [`scripts/download/download_gse245601.sh`](../scripts/download/download_gse245601.sh),
  checksum-verified 26-sample manifest at
  [`results/tables/gse245601_sample_manifest.tsv`](../results/tables/gse245601_sample_manifest.tsv).
- **Pre-analysis plan:** [`gse245601_PREANALYSIS.md`](gse245601_PREANALYSIS.md).

## GSE240112 (single-cell RNA-seq)

- **Source:** NCBI GEO accession GSE240112.
- **System:** unpaired primary and recurrent ER+ tumors from different
  patients (3 primary, 3 recurrent) -- **not matched or paired**. Primary
  tissue is sourced from OriGene Technologies Inc. and recurrent tissue
  from the Ontario Tumor Bank, two different source institutions; no
  pairing statement exists anywhere in the GEO metadata or the source
  paper's Methods. Disease state is therefore confounded with
  biobank/source institution and this must be disclosed alongside any
  result from this dataset. See `GSE240112_PREANALYSIS.md` Section C and
  `GSE240112_DATA_AUDIT.md` Section 2 for the full audit this correction
  is based on. (Corrected 2026-08-15, post-audit sensitivity analysis --
  the frozen GSE240112 analysis itself was already run correctly as
  unpaired; only this summary line was wrong.)
- **Pre-analysis plan:** [`GSE240112_PREANALYSIS.md`](GSE240112_PREANALYSIS.md).

## GSE111151 (bulk RNA-seq)

- **Source:** NCBI GEO accession GSE111151.
- **System:** 4 parental + 7 tamoxifen-resistant (TamR) derivative cell lines.
- **Pre-analysis plan:** [`GSE111151_PREANALYSIS.md`](GSE111151_PREANALYSIS.md).

## TCGA-BRCA

- **Source:** NCI Genomic Data Commons (GDC), TCGA-BRCA project.
- **Download:** [`scripts/download/download_tcga_brca.py`](../scripts/download/download_tcga_brca.py).
- **Cohort:** 1,095 primary tumors after deduplicating patients with >1
  RNA-seq aliquot (1,106 raw samples -> 1,095 one-aliquot-per-patient).
- **Gene ID mapping:** [`data/reference/tcga_candidate_ensembl_ids.tsv`](../data/reference/tcga_candidate_ensembl_ids.tsv)
  (independently verified Ensembl IDs for USP34/VEZF1/EML5/CITED2).
- **Caveat:** a general BRCA cohort, **not** a tamoxifen-resistance cohort.

## DepMap Public

- **Source:** depmap.org.
- **Releases used:** **24Q4** (archived, `results/tables/*/archive_24Q4/`)
  and **26Q1** (active/reported release).
- **Download:** [`scripts/download/download_depmap.py`](../scripts/download/download_depmap.py),
  [`download_depmap_26q1.py`](../scripts/download/download_depmap_26q1.py),
  [`verify_depmap_26q1_manual.py`](../scripts/download/verify_depmap_26q1_manual.py).
- **26Q1 access note:** `CRISPRGeneDependency.csv` for 26Q1 required manual
  retrieval; full access-channel documentation (Cloudflare, Figshare, AnVIL
  attempts) is in
  [`results/reports/independent_validation/DEPMAP_26Q1_ACCESS_STATUS.md`](../results/reports/independent_validation/DEPMAP_26Q1_ACCESS_STATUS.md).
- **Excluded file:** an earlier, unconfirmed-provenance `gene_effect_chronos_params.csv`
  from Figshare is explicitly excluded from the active pipeline (confirmed
  genuinely different from the trusted `CRISPRGeneEffect.csv` at the raw-data
  provenance level).
- **Cell-line/subtype classification:** ER+/luminal status is derived
  directly from DepMap's own `Model.csv`/`ModelSubtypeFeatures` fields
  (`src/independent_validation_depmap_data.py::load_model()`), never
  hardcoded or fuzzy-matched.

## GDSC / CancerRxGene Release 8.5

- **Source:** `https://cog.sanger.ac.uk/cancerrxgene/GDSC_release8.5/`
  (direct Sanger COG storage URL; the CancerRxGene bulk-download page
  returned HTTP 410 at the time of download and was not used).
- **Release date:** 2023-10-30.
- **Files:** `GDSC1_fitted_dose_response_27Oct23.xlsx`,
  `GDSC2_fitted_dose_response_27Oct23.xlsx`,
  `screened_compounds_rel_8.5.csv`, `Cell_Lines_Details.xlsx`.
- **Download:** [`scripts/download/download_gdsc.py`](../scripts/download/download_gdsc.py)
  (documents exact URLs and expected SHA256 per file).
- **Raw provenance record:** `PROVENANCE.txt` alongside the raw files
  (SHA256, row/line/drug counts, verified response-metric definitions,
  endocrine-compound search results).
- **Expression source for the GDSC correlation:** DepMap Public 26Q1
  (reused, not a new RNA-expression dataset), joined via the exact
  `SangerModelID`/`COSMIC_ID` identifiers present natively in both files —
  no fuzzy cell-line-name matching.
- **Response-metric direction:** independently verified against
  `GDSC_Fitted_Data_Description.pdf` (Sanger, v1.0.0, 21 Sep 2017) — LOWER
  LN_IC50 and LOWER AUC both mean MORE SENSITIVE.
- **Scope:** used for USP34 and VEZF1 only, per the frozen lead/backup
  translational conclusion — not new candidate discovery.

## Reference gene sets and interaction networks

- **Gene sets:** Hallmark, Reactome, GO Biological Process — downloaded via
  [`scripts/download_genesets.sh`](../scripts/download_genesets.sh),
  timestamped in [`data/reference/genesets/download_timestamp_utc.txt`](../data/reference/genesets/download_timestamp_utc.txt).
- **Protein interactions:** STRING — downloaded via
  [`scripts/download_string_interactions.py`](../scripts/download_string_interactions.py).
- **Transcription-factor targets:** TRRUST — `data/reference/interactions/trrust_human.tsv`.
- **Essential-gene reference:** Hart et al. 2017 CEG2 —
  [`scripts/download/download_ceg2.py`](../scripts/download/download_ceg2.py).

---

## General provenance rules enforced throughout

- Every download happens in a one-time `scripts/download/` script, never at
  analysis runtime (`src/` modules make no network calls).
- Every raw-data directory carries its own `PROVENANCE.txt` with source
  URL, retrieval date, and checksums where computable.
- Config (`config/config.yaml`) is the single source of truth for every
  file path; no path is hardcoded in analysis code.
