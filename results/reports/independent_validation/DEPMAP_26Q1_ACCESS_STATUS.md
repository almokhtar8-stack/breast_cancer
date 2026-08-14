# DepMap Public 26Q1 access status -- BLOCKED, awaiting manual download

**Status as of 2026-08-14: the 24Q4 -> 26Q1 DepMap update is STOPPED
pending manual download.** No 26Q1 analysis has been run; no existing
24Q4 output has been changed. This document is the full access-attempt
log and the exact manual-download instructions requested for this
situation.

---

## Why this update is needed (confirmed)

DepMap's announced 2026 release schedule is **26Q1 and 26Q3** (no 26Q2),
per `https://forum.depmap.org/t/depmap-quarterly-release-notes/3560`
(checked 2026-08-14). 26Q1 was announced 2026-04-01 at
`https://forum.depmap.org/t/announcing-the-26q1-release/4606`. This
confirms 26Q1 is the correct target release; the previous phase's use of
24Q4 was correctly disclosed at the time as "the most recent release
actually accessible during that run," not an error.

## What was attempted (official, non-interactive channels)

Not claimed to be an exhaustive enumeration of every conceivable channel
-- these are the specific ones tried, each with its exact result. Row 2's
"Figshare/Figshare+ search API" query phrasings tried, verbatim: "DepMap
Public 26Q1", "DepMap 26Q1 Public", "Model (Public 26Q1)", "CRISPR (Public
26Q1)", "Omics (Public 26Q1)", "26Q1 Public CRISPRGeneEffect", "Models
(Public 26Q1)", "Model.csv 26Q1", "OmicsExpressionProteinCodingGenesTPMLogp1
26Q1", "26Q1 Model", "PublicCellLineMetadata 26Q1", "26Q1 Model.csv",
"DepMap Model 26Q1", "Public 26Q1 Model", "CRISPRGeneEffect Public 26Q1",
"26Q1 metadata cell line" -- all via `POST api.figshare.com/v2/articles/search`.

| # | Channel | Result |
|---|---|---|
| 1 | `depmap.org/portal/api/download/all`, `depmap.org/portal/data_page/?tab=allData` (curl AND the WebFetch tool) | **Cloudflare Turnstile "verifying you're a person" challenge** returned to every non-interactive client. Not bypassed. |
| 2 | Figshare/Figshare+ search API (`api.figshare.com/v2/articles/search`), >15 query phrasings incl. "DepMap Public 26Q1", "DepMap 26Q1 Public", "Model 26Q1", "CRISPRGeneEffect 26Q1", "OmicsExpressionProteinCodingGenesTPMLogp1 26Q1" | Only ONE 26Q1 item found: Figshare+ article **31660582**, *"Chronos parameters (Public 26Q1)"*. **CORRECTION (caught by independent Codex review, re-verified directly against the Figshare API on 2026-08-14): this item's author is "Yejie Yun" (Figshare author id 21559277), NOT the "Broad DepMap" institutional account (id 17476659) that published the verified 24Q4/24Q2/23Q4 bundles.** An earlier version of this document incorrectly claimed the same-account match; that claim was wrong and is retracted. No independent confirmation that this individual is affiliated with the DepMap team could be found. The item's title/description ("DepMap 26Q1 Public Chronos parameters") is the only basis for associating it with DepMap at all -- this is materially weaker provenance than initially stated, and the file downloaded from it (see below) should be treated as **unverified / DepMap-labeled only, not confirmed-official**, on top of the pre-existing corrected/uncorrected ambiguity. Contains `gene_effect.csv`, `t0_offset.csv`, `guide_efficacy.csv`, `replicate_efficacy.csv`, `library_effect.csv` -- **no Model.csv, no CRISPRGeneDependency.csv, no expression file.** |
| 3 | Figshare author-listing page (`figshare.com/authors/Broad_DepMap/...`, both `figshare.com` and `plus.figshare.com`) | Client-rendered SPA shell (HTTP 202/403, empty body) -- cannot be scraped without executing JavaScript, which this project does not do. |
| 4 | Google Cloud Storage bucket referenced in the 26Q1 release notes (`storage.googleapis.com/shared-portal-files/...`) | Individually-known object paths ARE publicly readable (verified: the mutation-pipeline PDF cited in the release notes returns HTTP 200), but **anonymous bucket listing is disabled**, so a data-file object key cannot be discovered without already knowing it. No data-file path is published anywhere found. |
| 5 | AnVIL/Terra (`anvilproject.org/news/2026/03/03/depmap-data-release`) | Requires NIH Researcher Authentication Service (RAS) login via Login.gov/ID.me plus a dbGaP accession (`phs003444.v3.p1`) for **controlled-access** data -- not usable non-interactively, and covers raw sequencing data rather than the processed public files this project needs. |

**What WAS downloaded, with download-integrity (not provenance) verified**
(channel #2 above):
- `gene_effect.csv` from Figshare+ article 31660582
- Size: 413,544,937 bytes; MD5 `c89e3ec7e2c3682e5c3535172177a1ee` -- **matches Figshare's reported MD5 exactly** (locally recomputed and confirmed -- this confirms the download is byte-identical to what Figshare serves, NOT that the underlying data is confirmed-official DepMap data)
- Shape: 1,208 models x 18,531 genes, `SYMBOL (Entrez)` column format -- same schema as the standard `CRISPRGeneEffect.csv`
- All four candidate genes present with Entrez IDs matching the previously-verified mapping: USP34 (9736), VEZF1 (7716), EML5 (161436), CITED2 (10370)
- Saved at `/ibex/scratch/aljaroaa/tamoxifen-data/depmap/26Q1/gene_effect_chronos_params.csv`, with full provenance in that directory's `PROVENANCE.txt`

**Two unresolved caveats on that file, both disclosed, neither glossed
over:**
1. **Uploader identity (downgraded after independent review -- see the
   correction in the table above):** the Figshare item's author is an
   individual account ("Yejie Yun"), not the institutional "Broad DepMap"
   account used for the previously-verified 24Q4/24Q2/23Q4 bundles. The
   item's title/description is DepMap-branded, but this is materially
   weaker provenance than a same-account match would be.
2. It comes from a Figshare item titled *"Chronos parameters"* -- the
   Chronos model's own fitted-parameter bundle. The standard DepMap
   release separately publishes a batch-corrected `CRISPRGeneEffect.csv`
   and an `CRISPRGeneEffectUncorrected.csv`; this item's one-line
   description does not state which of the two `gene_effect.csv` here
   corresponds to.

**This file is NOT used in any analysis for either reason above** --
resolution requires either obtaining the canonical, fully-documented
release bundle from the official portal (which also supplies the still-
missing Model.csv and settles both caveats at once), or independent
confirmation of this specific file's origin and processing status. The
manual-download path below is recommended over trying to further validate
this file.

**Not obtained via any official non-interactive channel:**
- `Model.csv` (cell-line metadata: lineage, `OncotreeSubtype`,
  `ModelSubtypeFeatures` -- **required** to define breast/ER+/luminal
  cohorts; this alone blocks Part 3 of the update regardless of the
  gene-effect file's status)
- `CRISPRGeneDependency.csv` (per-line dependency probability)
- `OmicsExpressionProteinCodingGenesTPMLogp1.csv` (matched CCLE expression)
- The canonical, unambiguous `CRISPRGeneEffect.csv`

## Decision

Per your explicit instruction: **do not scrape, do not bypass the
Cloudflare challenge, do not use an unofficial mirror, do not fall back to
24Q4 for the primary analysis, do not fabricate 26Q1 results.** The
DepMap-26Q1 calculation is stopped here. `config/config.yaml`'s
`independent_validation.depmap.active_release` remains `"24Q4"`; every
existing `results/tables/independent_validation/DepMap_*` /
`four_candidate_independent_validation.tsv` file and figures 3-4 are
byte-identical to before this session (verified by regression diff, see
the main report's provenance note and the git-status/test report below).

---

## Manual download instructions

1. Open **`https://depmap.org/portal/data_page/?tab=allData`** in a
   regular browser (this passes the Cloudflare human-verification
   automatically). If it does not load the downloads table, click through
   the "I'm not a robot" / Turnstile prompt first.
2. In the release selector, choose **`DepMap Public 26Q1`** explicitly
   (do not use "Latest" if the two differ by the time you check -- record
   whichever release string the page shows).
3. Download exactly these four files (search/filter by name on that
   page):
   - `Model.csv`
   - `CRISPRGeneEffect.csv` (**not** `CRISPRGeneEffectUncorrected.csv` --
     use the standard, batch-corrected file, matching what was used for
     24Q4)
   - `CRISPRGeneDependency.csv`
   - `OmicsExpressionProteinCodingGenesTPMLogp1.csv`
4. Note the exact release string, download date, and (if the download
   page shows one) any checksum/version identifier displayed next to each
   file -- paste these into a text file alongside the downloads if
   possible, for the record.
5. Place all four files, **using exactly these filenames**, in:
   ```
   /ibex/scratch/aljaroaa/tamoxifen-data/depmap/26Q1/
   ```
   (this directory already exists and currently contains only
   `gene_effect_chronos_params.csv` and `PROVENANCE.txt` from the partial
   automated attempt above -- leave those two files in place, they will
   not interfere).
6. Run the verification script from the project root:
   ```
   micromamba run -n bc python3 scripts/download/verify_depmap_26q1_manual.py
   ```
   This checks each file's schema, confirms all four candidate genes are
   present, computes a SHA256 for the record, and writes an updated
   `PROVENANCE.txt`. It will fail loudly and tell you exactly what's wrong
   if a file is missing, misnamed, or doesn't look like the expected
   schema.
7. Tell me it's done (or just say "continue" / re-invoke this task) --
   once verification passes, I will set
   `independent_validation.depmap.active_release: "26Q1"` in
   `config/config.yaml` and run the fully-prepared analysis
   (`src/independent_validation_depmap.py`,
   `src/independent_validation_depmap_codependency.py`,
   `src/independent_validation_depmap_comparison.py`, the DepMap-only
   figures, and the report update) immediately -- no further code changes
   should be needed.

## Code already prepared and waiting

- `src/independent_validation_depmap_data.py` -- refactored to take an
  explicit `release` argument on every loader (`load_model`,
  `load_gene_effect`, `load_gene_dependency`, `load_expression`); reads
  `config/config.yaml`'s `independent_validation.depmap.releases.<release>`
  sub-config. Regression-tested: rerunning with `release="24Q4"` reproduces
  the existing tables/figures byte-for-byte identically (aside from a new
  `depmap_release` provenance column).
- `src/independent_validation_depmap.py`,
  `src/independent_validation_depmap_codependency.py`,
  `src/independent_validation_visualization.py`'s `build_figure_03` --
  all take an explicit/optional `release` parameter, defaulting to
  `active_release`.
- `src/independent_validation_depmap_comparison.py` -- builds
  `DepMap_24Q4_vs_26Q1_comparison.tsv` given both releases' data; ready to
  run the instant 26Q1 files are verified.
- `config/config.yaml`'s `independent_validation.depmap` section already
  has a `releases.26Q1` sub-config with the standard 4 filenames pinned
  (matching step 3 above exactly) and `data.raw.depmap_26q1_dir` pointing
  at the directory in step 5.
