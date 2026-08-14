"""Shared DepMap data loaders for the independent-validation phase.

Supports multiple DepMap releases side by side (config/config.yaml's
independent_validation.depmap.releases), since the project keeps both the
completed, archived 24Q4 analysis and the now-active 26Q1 analysis
directly comparable without ever silently mixing files from the two
releases. Every loader function REQUIRES an explicit `release` argument
(no default) so a caller can never accidentally read the wrong release's
files. Exact release/URLs/access rationale are documented in
scripts/download/download_depmap.py (24Q4, programmatic) and
scripts/download/download_depmap_26q1.py + verify_depmap_26q1_manual.py
(26Q1, manually downloaded from the official portal and verified -- see
results/reports/independent_validation/DEPMAP_26Q1_ACCESS_STATUS.md).
26Q1's CRISPRGeneDependency.csv was added in a follow-up manual download
and is now verified and in use (has_dependency_probability() returns True
for both current releases); the flag exists so a future release missing
that file degrades gracefully (E_INSUFFICIENT_DATA) instead of crashing
or being silently skipped.

No network calls happen here.
"""

from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd
import yaml

logger = logging.getLogger(__name__)

# Entrez gene IDs, verified independently via mygene.info and
# rest.ensembl.org/xrefs/symbol on 2026-08-14 (same verification pass as
# data/reference/tcga_candidate_ensembl_ids.tsv). Release-independent.
CANDIDATE_ENTREZ = {"USP34": 9736, "VEZF1": 7716, "EML5": 161436, "CITED2": 10370}


def load_config(config_path: str | Path = "config/config.yaml") -> dict:
    with open(config_path) as f:
        return yaml.safe_load(f)


def _release_cfg(cfg: dict, release: str) -> dict:
    releases = cfg["independent_validation"]["depmap"]["releases"]
    if release not in releases:
        raise KeyError(f"unknown DepMap release {release!r}; known releases: {sorted(releases)}")
    return releases[release]


def raw_dir(cfg: dict, release: str) -> Path:
    rel_cfg = _release_cfg(cfg, release)
    return Path(cfg["data"]["raw"][rel_cfg["raw_dir_key"]])


def _raw_path(cfg: dict, release: str, key: str) -> Path:
    rel_cfg = _release_cfg(cfg, release)
    return raw_dir(cfg, release) / rel_cfg["raw"][key]


def _gene_col(matrix_columns: pd.Index, symbol: str) -> str:
    target = f"{symbol} ({CANDIDATE_ENTREZ[symbol]})"
    if target not in matrix_columns:
        raise KeyError(f"{target} not found in DepMap matrix columns")
    return target


def load_model(cfg: dict, release: str) -> pd.DataFrame:
    """Cell-line metadata with is_breast / is_er_luminal columns.

    ER+/luminal status uses DepMap's own curated ModelSubtypeFeatures free
    -text field (never a manually recalled cell-line list): a breast line
    counts as ER+/luminal if that field contains "ER+" or the "ER," shorthand
    DepMap uses for "ER positive, other marker follows" (e.g. "luminal ER,
    PR+"); "luminal" alone is NOT sufficient (some "luminal" calls are
    explicitly TNBC or HER2-only in this field and must not be pulled in).
    Applied identically regardless of release; if a future release changes
    this field's format, `load_model` will simply select 0 ER+/luminal
    lines rather than silently misclassifying -- callers should sanity
    -check n_er_luminal > 0.
    """
    path = _raw_path(cfg, release, "model_csv")
    m = pd.read_csv(path, index_col="ModelID")
    n_total = len(m)
    if m.index.duplicated().any():
        raise ValueError(f"load_model ({release}): duplicate ModelID rows in {path}")
    m["is_breast"] = m["OncotreeLineage"] == "Breast"
    msf = m["ModelSubtypeFeatures"].fillna("")
    m["is_er_luminal"] = m["is_breast"] & msf.str.contains(r"ER\+|ER,", regex=True)
    logger.info(
        "load_model (%s): %d cell lines total; %d breast; %d ER+/luminal breast (ModelSubtypeFeatures rule)",
        release, n_total, int(m["is_breast"].sum()), int(m["is_er_luminal"].sum()),
    )
    return m


def load_gene_effect(cfg: dict, release: str, genes: list[str]) -> pd.DataFrame:
    path = _raw_path(cfg, release, "crispr_gene_effect_csv")
    df = pd.read_csv(path, index_col=0)
    cols = {_gene_col(df.columns, g): g for g in genes}
    out = df[list(cols)].rename(columns=cols)
    logger.info("load_gene_effect (%s): %d cell lines x %d genes", release, *out.shape)
    return out


def load_gene_dependency(cfg: dict, release: str, genes: list[str]) -> pd.DataFrame:
    path = _raw_path(cfg, release, "crispr_gene_dependency_csv")
    df = pd.read_csv(path, index_col=0)
    cols = {_gene_col(df.columns, g): g for g in genes}
    out = df[list(cols)].rename(columns=cols)
    logger.info("load_gene_dependency (%s): %d cell lines x %d genes", release, *out.shape)
    return out


def has_dependency_probability(cfg: dict, release: str) -> bool:
    return bool(_release_cfg(cfg, release).get("has_dependency_probability", True))


def load_expression(cfg: dict, release: str, genes: list[str]) -> pd.DataFrame:
    """Samples x genes expression matrix, indexed by ModelID.

    Two schemas exist across releases (config's expression_schema field):
    - "simple_modelid_index" (24Q4): one row per model, first column IS
      the ModelID.
    - "profile_level_default_flag" (26Q1): one row per sequencing
      profile, with a `ModelID` column and an `IsDefaultEntryForModel`
      flag DepMap itself sets to mark the canonical profile per model
      (some models have >1 profile). Filtered to that flag, never a
      guessed row-selection rule.
    """
    path = _raw_path(cfg, release, "expression_csv")
    schema = _release_cfg(cfg, release).get("expression_schema", "simple_modelid_index")
    if schema == "profile_level_default_flag":
        df = pd.read_csv(path, index_col=0)
        n_profiles = len(df)
        df = df.loc[df["IsDefaultEntryForModel"] == "Yes"].set_index("ModelID")
        if df.index.duplicated().any():
            raise ValueError(f"load_expression ({release}): duplicate ModelID among IsDefaultEntryForModel=='Yes' rows in {path}")
        logger.info("load_expression (%s): %d profile rows -> %d default-per-model rows", release, n_profiles, len(df))
    elif schema == "simple_modelid_index":
        df = pd.read_csv(path, index_col=0)
    else:
        raise ValueError(f"load_expression ({release}): unknown expression_schema {schema!r}")

    cols = {_gene_col(df.columns, g): g for g in genes}
    out = df[list(cols)].rename(columns=cols)
    logger.info("load_expression (%s): %d cell lines x %d genes", release, *out.shape)
    return out
