---
title: GSE245601 annotation concordance — feasibility report
status: post_freeze_exploratory
analysis_date: 2026-08-18
branch: post-poster-strengthening
verdict: FEASIBILITY GATE FAILED — analysis not performed
---

# GSE245601 annotation concordance: feasibility report

**post_freeze_exploratory.** Nothing here changes, or is derived from, any
frozen result. No frozen label was inspected for the purpose of altering it.

## Verdict

**The feasibility gate failed. No concordance matrix, no Cohen's kappa and no
per-tumour concordance breakdown were produced, and no proxy comparison was
attempted.**

The requested analysis needs the authors' own per-cell `Epi. Tumor` /
`Epi. Nontumor` labels, keyed by cell barcode and sample, so they can be matched
against our inferCNV calls. Those labels are not publicly available. Building
them ourselves in order to compare them against our own reconstruction would
measure agreement between two runs of our pipeline, not agreement with the
paper — which is the circularity the brief specifically ruled out.

Machine-readable outputs:

| File | Contents |
|---|---|
| `results/post_poster/annotation_concordance/feasibility_probe.json` | raw record of every source queried |
| `results/post_poster/annotation_concordance/feasibility_sources_checked.tsv` | one row per source, with the evidence |
| `results/post_poster/annotation_concordance/feasibility_verdict.tsv` | the gate result |
| `results/post_poster/annotation_concordance/our_frozen_label_inventory.tsv` | our side of the matrix that could not be built |

Reproduce with:

```
python scripts/post_poster_probe_gse245601_annotations.py   # network
python -m src.post_poster_annotation_feasibility            # deterministic
```

## What was checked, and what exists

Probe date **2026-08-18**.

### 1. GEO series GSE245601 and all 26 samples

- Series supplementary files: **1** — `GSE245601_RAW.tar`.
- Sample supplementary files: **26**, every one of them a `.h5`.

Every public GEO file is a Cell Ranger filtered feature-barcode matrix. A
per-cell label table would appear as a separate `.csv`/`.tsv`/`.rds`
supplementary file. None exists. This matches what
`docs/gse245601_PREANALYSIS.md` §1 already recorded: *"No cell-level metadata,
cell-type labels, or malignancy calls are included in the public GEO archive."*
That statement is confirmed, not merely repeated.

### 2. The paper's supplementary material (PMC10690085)

Retrieved through the Europe PMC open REST API (`/supplementaryFiles`) — 23
files: 6 supplementary figure PDFs, 7 supplementary tables, and the main figure
images.

**Every spreadsheet was opened and every cell scanned.** 76 worksheets across
the six `.xlsx` files; largest sheet 910 rows; **0 cells matching a 10x cell
barcode**. All sheets are gene-level: differential-expression tables
(`gene`, `p_val`, `avg_log2FC`, `pct.1`, `pct.2`, `p_val_adj`), Hallmark/UNC
enrichment results, and 262–301-gene resistance signatures. Supplementary Table
S1 is a `.docx` describing participant representativeness — 1,107 characters,
no barcodes.

A 40,000-cell annotation cannot hide in a 910-row sheet.

### 3. The authors' code repository at the pinned commit

`hyunsoo77/BC_tamoxifen_response` at
`ceabf3f331c88f464e6a57b0ad1f9c500bedde85` — the commit pinned in
`docs/gse245601_PREANALYSIS.md`, verified still to be the repository's current
`main` HEAD.

- 12 commits over all refs; **0 releases, 0 tags**; no `.gitattributes`, so no
  Git LFS.
- File extensions ever committed across the entire history:
  `.ipynb`, `.md`, `.png`, `.r`. **No data file has ever existed in this
  repository**, at any commit.
- Barcode scan of every notebook: the largest match count in any notebook is
  **48**, in `figure2_02_dge.ipynb`. Inspected directly, these are truncated
  dataframe previews (`<li>'T47D_Control#AAACCCACAGTCAGTT-1'</li>`) of **T47D
  cell-line** barcodes carrying PAM50 subtype calls — a sample that is not even
  in our primary tumour cohort, and not a malignancy label.

The repository publishes the method, not the annotations.

### 4. dbGaP (recorded, not queried)

The paper's Data Availability statement routes processed scRNA-seq data to
**phs003186.v1.p1**, controlled access. Not queried. Controlled access is a
possible future source; it is not a publicly available comparator, so it cannot
close this gate.

## What our side of the matrix would have been

`results/tables/gse245601_malignant_cell_labels.tsv` holds **29,175 epithelial
cells** across 10 tumours × 2 conditions (the frozen 44,140-cell metadata covers
all lineages; only epithelial cells receive a malignancy call). Per-tumour
counts are in `our_frozen_label_inventory.tsv` and were read without
modification.

The inventory does show, without any comparison being made, why the three Track
B eligible tumours are the ones they are: Tumor_02 (134/355 and 183/325
malignant), Tumor_03 (276/513 and 582/1105) and Tumor_07 (80/360 and 146/446)
carry substantial malignant fractions, whereas Tumor_01 calls 3 of 6,108 cells
malignant and Tumor_04 calls 1 of 221. That spread is the sample-dependence
`docs/CNV_METHOD_AUDIT.md` already established. **How it relates to the
authors' calls remains unmeasured**, and this report does not let it be guessed
at.

## Limitation of this report

The 44,140-cell figure is verifiable in this repository
(`docs/PROJECT_STATUS.md`, `docs/CODE_MAP.md`,
`docs/GSE245601_CANDIDATE_DEEPDIVE.md`). The published 40,428 comparison figure
could **not** be traced to any artefact in this repository; it appears in the
task brief but has no repository provenance. Since the gap between the two
counts cannot be attributed to a specific QC stage without the authors'
per-cell data, that decomposition is also not attempted here.

## Why our malignant calls differ from the paper's — the two-sentence answer

We do not know, and with the data the authors have released publicly we cannot
find out: they published Cell Ranger count matrices and their analysis code, but
never the per-cell malignancy labels themselves, which sit behind controlled
access at dbGaP. What we can say is that the difference is not uniform — our
inferCNV calls range from 3 malignant cells in 6,108 for Tumor_01 to 858 in
1,618 for Tumor_03 — so any single number for "how far apart we are" would be
misleading even if we could compute one.
