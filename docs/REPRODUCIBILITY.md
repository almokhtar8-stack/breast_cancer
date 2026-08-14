# Reproducibility

Step-by-step guide to reproducing this project's analyses, plus a candid
list of the gaps in that reproducibility.

## Environment

```bash
micromamba env create -f environment.yml -n bc
micromamba activate bc
```

`environment.yml` pins: Python 3.11, pandas, numpy, scipy, scikit-learn,
statsmodels, pyarrow, openpyxl, pytest, R 4.5.3 with edgeR 4.8.2, limma
3.66.0, and statmod — the packages invoked directly from Python via
`subprocess` (`scripts/analysis/*_limma.R`, `*_edger.R`).

**Known gap:** the single-cell R pipeline
(`scripts/analysis/gse245601_*.R`, `gse240112_*.R`) uses Seurat, InferCNV,
and CopyKAT, which were installed in a **separate R environment on KAUST
Ibex** and are **not** captured in `environment.yml`. Reproducing the
single-cell preprocessing/malignant-cell-calling steps requires installing
these packages separately (see
[`external_refs/inferCNV`](../external_refs/inferCNV) and
[`external_refs/copykat`](../external_refs/copykat) for the exact reference
tool versions audited in [`CNV_METHOD_AUDIT.md`](CNV_METHOD_AUDIT.md)).
Reproducing the downstream Python analyses (from Phase 6 onward) does
**not** require this — it reads the already-frozen single-cell output
tables.

## Running the test suite

```bash
micromamba run -n bc python3 -m pytest -q
```

Current state: **1,150 passed, 1 skipped.** Every `src/` module has a
corresponding test module that exercises real logic (recomputes statistics
independently, checks specific numeric values, checks pseudoreplication/
duplication invariants) rather than only checking that code runs without
error.

## Running an individual analysis module

Every phase follows the same pattern: a `_data.py` or `_gdsc_data.py`
module (loaders/computation), a `_build_tables.py` module (writes TSVs to
`results/tables/<phase>/`), and a `_visualization.py` module (writes PNGs
to `results/figures/<phase>/`). Example (the most recent phase):

```bash
micromamba run -n bc python3 -m src.final_pharmacogenomics_build_tables
micromamba run -n bc python3 -m src.final_pharmacogenomics_visualization
```

No module makes a network call at runtime — all data are already
downloaded (or, for `results/tables/`, already committed). Downloads are
one-time scripts under `scripts/download/`, run manually and separately.

## Downloading raw data (only needed to fully re-derive from scratch)

```bash
bash scripts/download/download_gse245601.sh
micromamba run -n bc python3 scripts/download/download_tcga_brca.py
micromamba run -n bc python3 scripts/download/download_depmap_26q1.py
micromamba run -n bc python3 scripts/download/download_gdsc.py
# etc. -- see scripts/download/ for the full list; .py scripts run under
# the bc environment, the one .sh script runs directly under bash
```

Each download script documents its exact source URL(s) and expected
checksums; each raw-data directory gets its own `PROVENANCE.txt` recording
what was actually retrieved. Raw data are written under `data/raw/`, a
symlink to external KAUST Ibex scratch storage (`/ibex/scratch/...`) — this
symlink target will not resolve outside that specific Ibex allocation, and
raw data are never committed to git. Downloading raw data from scratch
therefore requires either Ibex access or independently pointing
`config/config.yaml`'s `data.raw` paths at your own copy of each dataset.

## Configuration

All file paths come from [`config/config.yaml`](../config/config.yaml) —
no path is hardcoded anywhere in `src/`. If you relocate the raw-data
directory, only `config/config.yaml` needs to change.

## Determinism

- Every `src/` module is deterministic given its inputs (no random seeds
  left unset where randomness is used, e.g. any resampling/permutation
  step declares its seed).
- Every module logs rows in / rows out / rows lost at each filter or join
  step — nothing is silently dropped.
- Nothing derived from the CRISPR screen enters the feature table; the
  screen supplies labels only (see [`../CLAUDE.md`](../CLAUDE.md)).

## What reproducing this project does NOT require

- No GPU.
- No paid data access — every dataset used is public.
- No wet-lab equipment — this is a fully computational reanalysis.

## What it does require, beyond `environment.yml`

- Access to KAUST Ibex (or equivalent SLURM cluster) is convenient but not
  required — the pipeline is not Ibex-specific beyond path configuration
  and the separate R single-cell environment noted above.
- For the single-cell phases specifically: Seurat, InferCNV, and CopyKAT
  installed separately (see the known gap above).
