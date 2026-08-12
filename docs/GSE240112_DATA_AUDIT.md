# GSE240112 data audit

**Date written:** 2026-08-12
**Status:** written from primary sources (GEO series/GSM records, PubMed/PMC
full text, author GitHub repository, author-hosted processed objects) before
any candidate-gene value was inspected. This document answers the Phase 1
data-audit questions only; it does not authorize any statistical threshold
(see `docs/GSE240112_PREANALYSIS.md` for that).

## 1. Pinned primary sources

| Item | Value |
|---|---|
| GEO series | GSE240112 |
| Publication | Fang K, Ohihoin AG, Liu T, Choppavarapu L, Nosirov B, Wang Y, Yu M, Kamaraju S, Leone G, Jin VX. "Integrated single-cell analysis reveals distinct epigenetic-regulated cancer cell states and a heterogeneity-guided core signature in tamoxifen-resistant breast cancer." *Genome Medicine* 2024;16:139 |
| DOI | 10.1186/s13073-024-01407-3 |
| PMID / PMCID | 39558215 / PMC11572372 |
| Author code repository | https://github.com/KunFang93/BRCA_TR_scRNAscATAC (scripts organized Fig1-Fig7; README explicitly warns the code uses Seurat v4.3.0 and "some might not be compatible with v5") |
| Zenodo archive | 10.5281/zenodo.8247774 -- source-code zip only (52.7 KB), no data |
| Processed objects | Dropbox-hosted (linked from GitHub README), external to GEO and Zenodo -- see section 4 |

## 2. Study design (from GEO series/GSM records + paper Methods)

- Tissue groups: 2 normal breast tissue (NT), 3 primary ER+ breast tumor
  (PT), 3 tamoxifen-treated recurrent tumor (RT). Matched scRNA-seq +
  scATAC-seq per sample (16 libraries total: 8 scRNA, 8 scATAC).
- scRNA platform: GPL18573, Illumina NextSeq 500, 10x Chromium 3' v3,
  Cell Ranger v7.0.1, reference GRCh38.
- Clinical characteristics (paper Methods): all tumors ER+; patients aged
  50-60; tumor grade G1-G3. PR/HER2 status and an explicit
  patient-pairing statement were not found in the extracted Methods text
  (Additional file 1: Table S1 is referenced by the paper but was not
  directly retrieved in this audit).
- **Tissue sourcing (relevant to pairing determination):** NT tissue was
  sourced from NDRI (Philadelphia); PT tissue from OriGene Technologies
  Inc.; RT tissue from the Ontario Tumor Bank. Three different source
  institutions for the three groups. **Corrected framing (after Codex
  review -- see `docs/GSE240112_CODEX_REVIEW.md`):** this makes a shared
  patient between any PT and RT sample unlikely, but does not logically
  prove it never occurred -- Additional file 1: Table S1, which might
  contain an explicit patient-pairing statement, was referenced by the
  paper but not directly retrieved in this audit. No statement anywhere
  in the GEO metadata or the paper Methods text that *was* retrieved
  claims within-patient primary/recurrent pairing. The correct basis for
  the unpaired design decision is the absence of any pairing evidence,
  not a proven absence of a shared patient: **no evidence of pairing was
  found, so the contrast is analyzed as unpaired** (per
  `docs/GSE240112_PREANALYSIS.md` section C).
- Treatment history / confounding: RT tissue is explicitly described as
  tamoxifen-treated recurrent disease; PT tissue is treatment-naive at the
  time of collection (recurrence context). No information was found on
  other systemic therapies (chemotherapy, other endocrine agents) the RT
  patients may have received, which is a real, unresolved confounder that
  cannot be corrected for with the available metadata -- documented as a
  limitation, not something to be estimated or assumed away.

## 3. Sample identity / GSM mapping

Three different naming conventions exist for the same six PT/RT
sequencing libraries, and they do **not** align by simple ordinal
position. Each is separately verified below; do not assume replicate
number = sample number across conventions.

| GSM | GEO `Sample_title` (declared replicate order) | Raw supplementary file prefix (internal Cell Ranger filename) | Processed h5Seurat `orig.ident` (assumed by declared-replicate-order alignment; see caveat) |
|---|---|---|---|
| GSM7681685 | NT replicate1 | NT1 | NT1 |
| GSM7681686 | NT replicate2 | NT2 | NT2 |
| GSM7681687 | PT replicate1 | PT1 | PT1 |
| GSM7681688 | PT replicate2 | PT2 | PT2 |
| GSM7681689 | PT replicate3 | **PT5** | PT3 |
| GSM7681690 | RT replicate1 | **RT3** | RT1 |
| GSM7681691 | RT replicate2 | **RT4** | RT2 |
| GSM7681692 | RT replicate3 | **RT6** | RT3 |

Caveat: the processed h5Seurat objects' `orig.ident`/cell-barcode-prefix
labels (`PT1/PT2/PT3`, `RT1/RT2/RT3`) are the authors' own clean relabeling
and do not literally match the raw GEO supplementary-file internal
prefixes (`PT1/PT2/PT5`, `RT3/RT4/RT6`). The alignment used throughout
this analysis is GEO's own declared `replicate1/2/3` order per group,
which is the only unambiguous, source-documented ordering available. This
mapping is used to interpret the raw Cell Ranger matrices (downloaded
separately for the Phase 6/14 secondary all-epithelial analysis, see
`data/raw/gse240112/cellranger/`); it does not affect the primary
tumor-cell pseudobulk analysis, which uses the processed h5Seurat object's
own `orig.ident` labels directly and never needs to touch the raw
GSM-prefix convention.

scATAC GSMs (GSM7681693-700) are not used in this run; see Phase 16
decision in `docs/GSE240112_ANALYSIS_REPORT.md`.

## 4. Available data representations, in Phase-4 preference order

1. **Author-provided cell-level processed objects with cell-level
   annotation (used as primary; Phase 4 tier 1).** Two h5Seurat objects
   are hosted on Dropbox and linked from the GitHub README (not in GEO or
   Zenodo):
   - `NTs.h5seurat` (671,765,179 bytes) -- all normal-tissue cells,
     7,529 cells x 25,543 genes, raw counts + SCT-normalized data, PCA/
     tSNE/UMAP, standard Seurat QC metadata (`nCount_RNA`, `nFeature_RNA`,
     `percent.mt`, `scDblFinder.class`/score fields), **no cell-type or
     malignancy annotation column**.
   - `TTs_cancer_060223.h5seurat` (4,482,023,978 bytes) -- **cancer cells
     only**, pooled across all 3 PT + 3 RT samples: 9,942 cells x 27,161
     genes, raw counts + SCT-normalized data, PCA/tSNE/UMAP,
     `integrated_snn_res.{0.5,0.7,0.9}` clustering, gene-module scores
     (AURKA/CASP3/ERBB2/ESR1/PLAU/STAT1/TR/VEGF), and a `cell.annot`
     metadata column that reads `"Breast cancer cells"` for all 9,942
     cells (i.e. this object *is* the authors' cancer-cell subset -- it
     was constructed by the authors' own marker-based tumor-cell calling,
     described in the paper Methods as `FindAllMarkers` against curated
     reference marker panels (PanglaoDB, CellMarker, ProteinAtlas,
     literature), explicitly **not** a CNV-based method). Per-sample
     `orig.ident` cell counts: PT1=1,029, PT2=1,975, PT3=1,442, RT1=2,721,
     RT2=2,597, RT3=178. `scDblFinder.class` is `singlet` for all 9,942
     cells (doublets already removed).
   - Both objects were confirmed to load raw counts (CSC sparse
     `data`/`indices`/`indptr` triplet), full metadata, and UMAP
     coordinates via direct low-level HDF5 access (`hdf5r`), working
     around a real, confirmed Seurat v4-vs-v5 API incompatibility in
     `SeuratDisk::LoadH5Seurat()` (fails with
     `GetAssayData(slot=...)` deprecation error on Seurat 5.5.1, exactly
     as the GitHub README's compatibility warning predicts). No object
     data was lost or reconstructed approximately by this workaround --
     every field read matches the file's own declared dimensions.
   - Gene symbols: 12 of the 13 frozen candidates match exactly once in
     the 27,161-gene feature space of the TT object (no duplicates
     anywhere in the feature list); `USP17L29` is absent by exact-symbol
     match (0 occurrences) -- flagged for the Phase 8 alias/mapping audit,
     not assumed untestable without investigation. `PAICS` is present
     exactly once.
2. **Raw Cell Ranger filtered matrices per GSM (Phase 4 tier 3; used only
   for the secondary Phase 6/14 all-epithelial analysis, since no
   author-provided all-cell-type PT/RT object exists).** GEO supplementary
   files for GSM7681687-92 (PT1/PT2/PT5, RT3/RT4/RT6 per the raw internal
   prefix) are standard Cell Ranger triplets
   (`barcodes.tsv.gz`/`features.tsv.gz`/`matrix.mtx.gz`), already
   cell-filtered (barcode counts in the thousands, not millions --
   confirmed by file size), total ~500 MB for all 6 PT/RT samples.
   Downloaded to `data/raw/gse240112/cellranger/`. These are used only to
   reconstruct a broad epithelial/immune/stromal/endothelial compartment
   label (standard canonical markers) for the Phase 6 broad-cell-type
   context and the Phase 14 all-epithelial sensitivity pseudobulk -- not
   for malignancy calling, which remains governed entirely by the
   author-provided `TTs_cancer_060223.h5seurat` object.
3. **No usable all-cell-type PT/RT processed object exists in the
   author's public release.** The only PT/RT-covering processed object
   (`TTs_cancer_060223.h5seurat`) is restricted to cancer cells only by
   construction; there is no `PTs.h5seurat`/`RTs.h5seurat` equivalent to
   the `NTs.h5seurat` normal-tissue object. This is a genuine data-gap,
   not an oversight in this audit -- documented so Phase 6/14 scope
   decisions are traceable to it.
4. Raw FASTQ was not needed and was not downloaded (tier 4, not reached).

## 5. Answers to the Phase 1 checklist

| # | Question | Answer |
|---|---|---|
| 1 | # scRNA samples | 8 (2 NT, 3 PT, 3 RT) |
| 2 | patient/sample identities | NT1-2, PT1-3, RT1-3 (GEO GSM7681685-92); see section 3 for the three-convention mapping |
| 3 | PT samples | 3 |
| 4 | RT samples | 3 |
| 5 | normal samples | 2 |
| 6 | treatment history | RT = tamoxifen-treated recurrent disease; PT = primary, presumed treatment-naive at collection; other systemic therapy history not available |
| 7 | ER/HER2 status | All tumors ER+ (paper Methods); HER2 status not found in extracted text |
| 8 | PT/RT paired or unpaired | Unpaired -- different source institutions per group (OriGene for PT, Ontario Tumor Bank for RT), no pairing statement in metadata |
| 9 | all RT exposed to tamoxifen | Yes, per group definition ("tamoxifen-treated recurrent") |
| 10 | other confounding systemic treatments | Unknown / not available -- documented limitation |
| 11 | scRNA technology | 10x Chromium 3' v3, Cell Ranger v7.0.1, GRCh38 |
| 12 | raw vs processed data availability | Raw (GEO Cell Ranger triplets) and author-processed (Dropbox h5Seurat) both available |
| 13 | author cell-level annotations available | Yes for cancer cells (`cell.annot` in the TT object); no broad cell-type annotation column in either object |
| 14 | malignant/tumor-cell labels available | Yes -- the TT object is itself the author's cancer-cell subset (Phase 7 Case A) |
| 15 | author UMAP/embedding available | Yes, in both h5Seurat objects (`reductions/umap/cell.embeddings`) |
| 16 | raw UMI counts suitable for pseudobulk | Yes, `assays/RNA/counts` in both h5Seurat objects, standard CSC sparse format |

## 6. Table deliverable

See `results/tables/gse240112_data_audit.tsv` (one row per sample).
