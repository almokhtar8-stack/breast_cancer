# Data

All data used by this project are **public**. Nothing here contains patient
identifiers or protected health information — only public GEO/DepMap/PDB sample
and model accessions.

```
data/
  raw/         External downloads. NOT in git (gitignored) -- symlinked/pointed to
               machine-local storage via config/config.yaml
  processed/   A small number of force-added frozen intermediate matrices -- in git
  reference/   Gene sets, STRING interaction tables, ID mappings -- in git
  checksums/   SHA-256 manifests for raw downloads
```

**Where paths come from.** Every path used by analysis code is declared in
[`config/config.yaml`](../config/config.yaml) under `data.raw` / `data.reference`.
The `data.raw` entries point at a **machine-local external data root** (on the
original machine, a cluster scratch directory). A stranger reproducing this work
must repoint those entries at their own download location; no analysis module
hardcodes a raw path.

---

## Datasets

| Dataset | Accession / source | Role in the project | Raw in git? | Downloadable? | Processed output in git? |
|---|---|---|---|---|---|
| **Hany CRISPR screen** | Hany et al., supplementary Data S1 | The primary discovery layer: genome-scale knockout screen, E2+4-OHT vs E2, 19,103 fitted genes → 13 significant sensitising hits | No | Yes, from the publication's supplementary material | Yes — `data/processed/labels.parquet` (force-added) |
| **GSE118713** | GEO | Cell-line tamoxifen-resistance model (MCF7 vs TAMR; FASR present in the dataset but excluded from the poster heatmap) | No | Yes (GEO) | Frozen DE tables in `results/tables/` |
| **GSE111151** | GEO | Independent tamoxifen-resistant sublines across 4 parental backgrounds (4 parental + 7 resistant). Frozen candidate DE was largely null | No | Yes (GEO) | Frozen DE tables in `results/tables/` |
| **GSE240112** | GEO | Primary vs recurrent human tumours — **unpaired**, different patients/biobanks. Recurrence-**associated** context, not a controlled tamoxifen-resistance experiment | No | Yes (GEO) | Frozen DE tables in `results/tables/` |
| **GSE245601** | GEO | Single-cell RNA-seq of acute **12 h, 10 µM ex vivo** tamoxifen exposure. Poster figures use **per-tumour pseudobulk** (3 eligible patients × Control/Tamoxifen), not per-cell rows. Malignant/non-malignant assignment is inferential (inferCNV/copyKAT) | No | Yes (GEO) | Frozen pseudobulk + DE tables in `results/tables/` |
| **DepMap 26Q1** | depmap.org portal | Baseline cancer-cell dependency (`CRISPRGeneEffect.csv`, `CRISPRGeneDependency.csv`, `Model.csv`). 11 ER+/luminal dependency-evaluable lines | No (~440 MB per matrix) | Yes, free registration at depmap.org | Yes — cached 11-line extract in `results/tables/poster_depmap_v1/` |
| **DepMap 24Q4** | depmap.org portal | Retained only for the release-comparison table | No | Yes | Comparison table in `results/tables/independent_validation/` |
| **STRING** | string-db.org REST API | Post-freeze exploratory network layer (species 9606, `required_score=700`) | **Yes** — `data/reference/interactions/*.tsv` (small) | Yes | Graph/Cytoscape exports in `results/tables/` |
| **MSigDB gene sets** | MSigDB (Hallmark, Reactome, GO:BP) | Pathway analysis and program annotation | **Yes** — `data/reference/genesets/*.gmt` | Yes | Frozen GSEA tables in `results/tables/systems_network/` |
| **PDB structures** | RCSB PDB — 6NQU, 5O0Y, 7W3R, 7W3U | Structural tractability figure | No (PDB files live in external storage) | Yes, free from rcsb.org | Yes — PyMOL renders in `results/figures/poster_druggability_v1/renders/` |
| **GDSC** | Sanger GDSC | *Exploratory/supporting only.* USP34 associations were GDSC1-only with **no GDSC2 replication**; not used in the poster | No | Yes | Tables in `results/tables/final_pharmacogenomics/` |
| **TCGA-BRCA** | GDC | *Supporting/orthogonal only*, weak and incomplete for the current focus genes; not used in the poster | No | Yes (GDC) | Tables in `results/tables/independent_validation/` |
| **Hart 2017 CEG2** | Publication supplementary | Core-essential-gene reference list for essentiality context | **Yes** — `data/reference/hart2017_ceg2_684.tsv` | Yes | — |

## Restrictions

- No dataset used here requires controlled access or a data-use agreement.
- Raw matrices are excluded from git purely for **size**, not for privacy.
- Sample identifiers that do appear (GEO GSM/patient codes such as `T02`,
  DepMap `ACH-` model IDs, cell-line names) are **public dataset accessions**,
  not personal identifiers.

## Download helpers

One-time, deterministic download scripts live in
[`scripts/download/`](../scripts/download/) and
[`scripts/download_string_network_v2_four_focus.py`](../scripts/download_string_network_v2_four_focus.py).
No analysis module performs network access at runtime — a hard project rule
(see [`CLAUDE.md`](../CLAUDE.md)). Checksums for downloads are recorded in
`data/checksums/` and in each external directory's `PROVENANCE.txt`.
