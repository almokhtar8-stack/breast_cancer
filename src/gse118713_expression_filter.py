"""Frozen-matrix validation and TPM expression filtering for GSE118713 Phase 2B.

Source: the frozen gene-level TPM matrix written by
``src.gse118713_prep.prepare_gse118713`` (config
``gse118713.output.gene_tpm_parquet``), checksum-pinned in config
(``gse118713.frozen_gene_tpm_sha256``). Sample/group metadata is read from
``gse118713.output.sample_metadata_tsv``, the single source of the
sample_id/group mapping (see PREANALYSIS.md's Phase 2 data-audit
amendment).

Implements PREANALYSIS.md's 2026-08-05 Phase 2B statistical-plan amendment:
genes are retained if TPM >= 1 in at least 3 of the 9 samples; all other
genes are dropped before any modelling. No differential expression,
pathway analysis, modelling, feature construction, or candidate ranking is
performed here -- this module only validates the frozen input and applies
the preregistered filter.
"""

from __future__ import annotations

import hashlib
import logging
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

logger = logging.getLogger(__name__)

ID_COLUMNS: tuple[str, ...] = ("gene_id", "gene_symbol", "symbol_mapping_status")


def _load_config(config_path: str | Path) -> dict:
    with open(config_path) as f:
        return yaml.safe_load(f)


def _sha256_file(path: str | Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


@dataclass(frozen=True)
class FilterConfig:
    """Resolved, config-driven paths and preregistered filtering rule."""

    gene_tpm_parquet: Path
    sample_metadata_tsv: Path
    expected_sha256: str
    expected_n_genes: int
    expected_n_samples: int
    expected_groups: tuple[str, ...]
    expected_replicates_per_group: int
    min_tpm: float
    min_samples: int
    filtered_gene_tpm_tsv: Path
    filtering_summary_tsv: Path

    @classmethod
    def from_config(cls, config: dict) -> "FilterConfig":
        gse = config["gse118713"]
        phase2b = config["gse118713_phase2b"]
        filtering = phase2b["filtering"]
        return cls(
            gene_tpm_parquet=Path(gse["output"]["gene_tpm_parquet"]),
            sample_metadata_tsv=Path(gse["output"]["sample_metadata_tsv"]),
            expected_sha256=gse["frozen_gene_tpm_sha256"],
            expected_n_genes=int(phase2b["expected_n_genes"]),
            expected_n_samples=int(phase2b["expected_n_samples"]),
            expected_groups=tuple(phase2b["expected_groups"]),
            expected_replicates_per_group=int(phase2b["expected_replicates_per_group"]),
            min_tpm=float(filtering["min_tpm"]),
            min_samples=int(filtering["min_samples"]),
            filtered_gene_tpm_tsv=Path(filtering["filtered_gene_tpm_tsv"]),
            filtering_summary_tsv=Path(filtering["filtering_summary_tsv"]),
        )


def load_sample_metadata(cfg: FilterConfig) -> pd.DataFrame:
    """Load and validate the frozen nine-sample metadata table.

    Rejects a metadata table with anything other than exactly the
    expected number of samples, duplicated sample IDs, or a group
    composition other than exactly ``expected_replicates_per_group``
    replicates per expected group.
    """
    meta = pd.read_csv(cfg.sample_metadata_tsv, sep="\t")
    if len(meta) != cfg.expected_n_samples:
        raise ValueError(f"expected exactly {cfg.expected_n_samples} samples in metadata, got {len(meta)}")
    if meta["sample_id"].duplicated().any():
        raise ValueError("duplicated sample_id values in GSE118713 sample metadata")

    group_counts = meta["group"].value_counts()
    missing_groups = sorted(set(cfg.expected_groups) - set(group_counts.index))
    if missing_groups:
        raise ValueError(f"metadata missing expected groups: {missing_groups}")
    extra_groups = sorted(set(group_counts.index) - set(cfg.expected_groups))
    if extra_groups:
        raise ValueError(f"metadata contains unexpected groups: {extra_groups}")
    wrong_counts = {
        g: int(group_counts[g])
        for g in cfg.expected_groups
        if group_counts[g] != cfg.expected_replicates_per_group
    }
    if wrong_counts:
        raise ValueError(
            f"expected exactly {cfg.expected_replicates_per_group} replicates per group, got: {wrong_counts}"
        )
    logger.info(
        "load_sample_metadata: validated %d samples across groups %s", len(meta), sorted(group_counts.index)
    )
    return meta


def validate_frozen_matrix(df: pd.DataFrame, cfg: FilterConfig, sample_ids: list[str]) -> None:
    """Validate the frozen gene-TPM matrix's structural invariants.

    Checks: exactly ``expected_n_genes`` unique nonblank gene IDs; exactly
    the expected sample columns present and no others beyond the ID
    columns; no supplied group-mean columns leaked through; and every TPM
    value is finite and non-negative.
    """
    if df["gene_id"].isna().any() or (df["gene_id"].astype(str).str.strip() == "").any():
        raise ValueError("null or blank gene_id values in frozen GSE118713 matrix")
    if df["gene_id"].duplicated().any():
        raise ValueError(f"{int(df['gene_id'].duplicated().sum())} duplicate gene_id values in frozen matrix")
    if len(df) != cfg.expected_n_genes:
        raise ValueError(f"expected exactly {cfg.expected_n_genes} gene rows, got {len(df)}")

    missing_samples = [c for c in sample_ids if c not in df.columns]
    if missing_samples:
        raise ValueError(f"frozen matrix missing expected sample columns: {missing_samples}")
    if len(sample_ids) != cfg.expected_n_samples:
        raise ValueError(f"expected exactly {cfg.expected_n_samples} sample columns, got {len(sample_ids)}")

    mean_like = [c for c in df.columns if c not in ID_COLUMNS and c not in sample_ids]
    if mean_like:
        raise ValueError(f"frozen matrix contains unexpected non-sample columns (possible mean leakage): {mean_like}")

    values = df[sample_ids].to_numpy(dtype=float)
    if not np.isfinite(values).all():
        raise ValueError("non-finite TPM value in frozen GSE118713 matrix")
    if (values < 0).any():
        raise ValueError("negative TPM value in frozen GSE118713 matrix")

    logger.info(
        "validate_frozen_matrix: %d genes x %d samples, all finite and non-negative", len(df), len(sample_ids)
    )


def verify_checksum(cfg: FilterConfig) -> str:
    """Verify the frozen matrix on disk matches the pinned Phase 2A checksum."""
    actual = _sha256_file(cfg.gene_tpm_parquet)
    if actual != cfg.expected_sha256:
        raise ValueError(
            f"frozen GSE118713 matrix checksum mismatch: expected {cfg.expected_sha256}, got {actual}"
        )
    logger.info("verify_checksum: frozen matrix checksum matches (%s)", actual)
    return actual


def load_frozen_matrix(cfg: FilterConfig, meta: pd.DataFrame) -> pd.DataFrame:
    """Verify checksum, load, and structurally validate the frozen matrix."""
    verify_checksum(cfg)
    df = pd.read_parquet(cfg.gene_tpm_parquet)
    validate_frozen_matrix(df, cfg, sample_ids=list(meta["sample_id"]))
    return df


def detect_gene(df: pd.DataFrame, sample_ids: list[str], min_tpm: float) -> pd.DataFrame:
    """Per-gene, per-sample detection matrix (``TPM >= min_tpm``)."""
    return df[sample_ids] >= min_tpm


def filter_expression(
    df: pd.DataFrame, sample_ids: list[str], min_tpm: float, min_samples: int
) -> tuple[pd.DataFrame, dict[str, object]]:
    """Apply exactly ``TPM >= min_tpm`` in at least ``min_samples`` samples.

    Returns the filtered (row-subset) matrix and a filtering record
    reporting genes in, genes retained, genes removed, the threshold
    applied, and per-sample detection counts -- no silent row loss.
    """
    detection = detect_gene(df, sample_ids, min_tpm)
    n_detected = detection.sum(axis=1)
    keep = n_detected >= min_samples

    genes_in = len(df)
    filtered = df.loc[keep].reset_index(drop=True)
    genes_retained = len(filtered)
    genes_removed = genes_in - genes_retained

    per_sample_detection_counts = {sample_id: int(detection[sample_id].sum()) for sample_id in sample_ids}

    record: dict[str, object] = {
        "genes_in": genes_in,
        "genes_retained": genes_retained,
        "genes_removed": genes_removed,
        "min_tpm_threshold": min_tpm,
        "min_samples_threshold": min_samples,
    }
    for sample_id in sample_ids:
        record[f"{sample_id}_genes_detected"] = per_sample_detection_counts[sample_id]

    logger.info(
        "filter_expression: genes_in=%d genes_retained=%d genes_removed=%d (TPM>=%s in >=%d/%d samples)",
        genes_in,
        genes_retained,
        genes_removed,
        min_tpm,
        min_samples,
        len(sample_ids),
    )
    return filtered, record


def write_filtered_matrix(filtered: pd.DataFrame, cfg: FilterConfig) -> None:
    cfg.filtered_gene_tpm_tsv.parent.mkdir(parents=True, exist_ok=True)
    # mtime=0 keeps the gzip header byte-identical across reruns of identical data.
    filtered.to_csv(
        cfg.filtered_gene_tpm_tsv,
        sep="\t",
        index=False,
        compression={"method": "gzip", "mtime": 0},
    )
    logger.info(
        "write_filtered_matrix: wrote %d filtered genes to %s", len(filtered), cfg.filtered_gene_tpm_tsv
    )


def write_filtering_summary(record: dict[str, object], cfg: FilterConfig) -> pd.DataFrame:
    out = pd.DataFrame([record])
    cfg.filtering_summary_tsv.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(cfg.filtering_summary_tsv, sep="\t", index=False)
    logger.info("write_filtering_summary: wrote filtering record to %s", cfg.filtering_summary_tsv)
    return out


def run_expression_filtering(config_path: str | Path = "config/config.yaml") -> tuple[pd.DataFrame, pd.DataFrame]:
    """Run the full Phase 2B expression-filtering step and write outputs."""
    config = _load_config(config_path)
    cfg = FilterConfig.from_config(config)

    meta = load_sample_metadata(cfg)
    df = load_frozen_matrix(cfg, meta)
    filtered, record = filter_expression(
        df, sample_ids=list(meta["sample_id"]), min_tpm=cfg.min_tpm, min_samples=cfg.min_samples
    )
    write_filtered_matrix(filtered, cfg)
    summary_df = write_filtering_summary(record, cfg)
    return filtered, summary_df


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run_expression_filtering()
