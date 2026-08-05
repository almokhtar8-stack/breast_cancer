"""Gene-level TPM preparation for GSE118713 (MCF7 parental / TAMR / FASR).

Source: GSE118713, "TRANSCRIPT_TAMR_FASR_MCF7_expression_data.tsv.gz",
public 2019-12-11 / last updated 2020-01-19 (config
``data.raw.gse118713_transcript_tpm``; registry entry in
``config/resources.yaml``).

The source file is transcript-level TPM (one row per Ensembl transcript),
not gene-level expression, despite the informal "MCF7 ... TPM" description
in README.md/PREANALYSIS.md §6 -- see PREANALYSIS.md's 2026-08-05 Phase 2
data-audit amendment. This module sums transcript TPM to gene-level TPM,
separately for each of the nine replicate columns (three replicates each
for MCF7, TAMR, FASR), grouped by unversioned Ensembl gene ID. The file's
own group-mean TPM columns (``MCF7.mean.TPM``, ``TAMR.mean.TPM``,
``FASR.mean.TPM``) are never read as sample data: they average over
replicates of one transcript, not over the transcripts of one gene, so
they cannot substitute for the sums computed here.

Unversioned Ensembl gene ID is the primary identifier. Gene symbols are
carried as annotation only and are never used to group or collapse rows,
because one symbol can be attached to more than one Ensembl gene ID.

No differential expression, pathway analysis, modelling, feature
construction, or biological interpretation is performed here -- this
module only loads, aggregates to gene level, and quality-checks.
"""

from __future__ import annotations

import hashlib
import logging
import re
from dataclasses import dataclass, replace
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

logger = logging.getLogger(__name__)

ID_COLUMNS: tuple[str, ...] = ("transcript.id", "gene.id", "gene.symbol")
GENE_ID_COL = "gene.id"
GENE_SYMBOL_COL = "gene.symbol"

SYMBOL_RESOLVED = "resolved"
SYMBOL_MISSING = "missing"
SYMBOL_AMBIGUOUS = "ambiguous"

_VERSION_SUFFIX = re.compile(r"\.\d+$")


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
class Gse118713Config:
    """Resolved, config-driven paths and sample layout. No hardcoded paths."""

    source_path: Path
    gene_tpm_parquet: Path
    sample_metadata_tsv: Path
    qc_tsv: Path
    sample_columns: tuple[str, ...]
    sample_ids: tuple[str, ...]
    groups: tuple[str, ...]
    replicates: tuple[int, ...]
    mean_columns: tuple[str, ...]

    @classmethod
    def from_config(cls, config: dict) -> "Gse118713Config":
        section = config["gse118713"]
        samples = section["samples"]
        return cls(
            source_path=Path(config["data"]["raw"]["gse118713_transcript_tpm"]),
            gene_tpm_parquet=Path(section["output"]["gene_tpm_parquet"]),
            sample_metadata_tsv=Path(section["output"]["sample_metadata_tsv"]),
            qc_tsv=Path(section["output"]["qc_tsv"]),
            sample_columns=tuple(s["column"] for s in samples),
            sample_ids=tuple(s["sample_id"] for s in samples),
            groups=tuple(s["group"] for s in samples),
            replicates=tuple(int(s["replicate"]) for s in samples),
            mean_columns=tuple(section["excluded_mean_columns"]),
        )

    @property
    def column_to_sample_id(self) -> dict[str, str]:
        return dict(zip(self.sample_columns, self.sample_ids))


def write_sample_metadata(cfg: Gse118713Config, config: dict) -> pd.DataFrame:
    """Write the nine-sample metadata table straight from configuration.

    No metadata here is inferred from the source file's column names --
    ``config["gse118713"]["samples"]`` is the single source of the
    sample_id/group/replicate mapping.
    """
    samples = config["gse118713"]["samples"]
    metadata = pd.DataFrame(samples)[["sample_id", "column", "group", "replicate"]]
    if len(metadata) != 9:
        raise ValueError(f"expected exactly 9 configured GSE118713 samples, got {len(metadata)}")

    cfg.sample_metadata_tsv.parent.mkdir(parents=True, exist_ok=True)
    metadata.to_csv(cfg.sample_metadata_tsv, sep="\t", index=False)
    logger.info("write_sample_metadata: wrote %d samples to %s", len(metadata), cfg.sample_metadata_tsv)
    return metadata


def load_transcript_tpm(cfg: Gse118713Config) -> pd.DataFrame:
    """Read the gzip transcript-TPM TSV and validate it before any arithmetic.

    Validates: all expected id columns and sample columns are present; the
    mean columns are present in the file but excluded from
    ``cfg.sample_columns`` (so they cannot accidentally be treated as
    samples downstream); and every replicate TPM value is numeric, finite,
    and non-negative.
    """
    df = pd.read_csv(cfg.source_path, sep="\t")
    logger.info("load_transcript_tpm: read %d transcript rows from %s", len(df), cfg.source_path)

    expected = list(ID_COLUMNS) + list(cfg.sample_columns) + list(cfg.mean_columns)
    missing = [c for c in expected if c not in df.columns]
    if missing:
        raise ValueError(f"GSE118713 source file missing expected columns: {missing}")

    overlap = set(cfg.mean_columns) & set(cfg.sample_columns)
    if overlap:
        raise ValueError(f"configured mean columns overlap configured sample columns: {sorted(overlap)}")
    logger.info(
        "load_transcript_tpm: confirmed %d group-mean columns present in file but excluded from samples: %s",
        len(cfg.mean_columns),
        cfg.mean_columns,
    )

    replicate_values = df[list(cfg.sample_columns)].apply(pd.to_numeric, errors="coerce")
    if replicate_values.isna().to_numpy().any():
        raise ValueError("non-numeric value in a GSE118713 replicate TPM column")
    arr = replicate_values.to_numpy(dtype=float)
    if not np.isfinite(arr).all():
        raise ValueError("non-finite value in a GSE118713 replicate TPM column")
    if (arr < 0).any():
        raise ValueError("negative TPM value in a GSE118713 replicate TPM column")
    df[list(cfg.sample_columns)] = replicate_values

    return df


def strip_ensembl_version(gene_ids: pd.Series) -> pd.Series:
    """Strip a trailing ``.<version>`` suffix from Ensembl gene IDs.

    E.g. ``ENSG00000004059.8`` -> ``ENSG00000004059``. Values with no
    version suffix are returned unchanged.
    """
    return gene_ids.astype(str).str.replace(_VERSION_SUFFIX, "", regex=True)


def validate_gene_ids(gene_ids: pd.Series) -> None:
    """Reject null or blank gene IDs rather than silently dropping them."""
    if gene_ids.isna().any():
        raise ValueError(f"{int(gene_ids.isna().sum())} null gene.id values in GSE118713 source file")
    blank = gene_ids.astype(str).str.strip() == ""
    if blank.any():
        raise ValueError(f"{int(blank.sum())} blank gene.id values in GSE118713 source file")


def _valid_symbol_mask(symbols: pd.Series) -> pd.Series:
    return symbols.notna() & (symbols.astype(str).str.strip() != "")


def reconcile_symbols(df: pd.DataFrame, gene_id_col: str = "gene_id") -> pd.DataFrame:
    """Map each unversioned gene ID to zero, one, or many distinct symbols.

    Per PREANALYSIS.md's Phase 2 data-audit amendment: a gene ID keeps its
    symbol only if exactly one unique nonblank symbol is observed across
    its transcript rows. Zero valid symbols -> missing (flagged, symbol
    left null). More than one distinct nonblank symbol -> ambiguous
    (flagged, symbol left unresolved/null). Symbols are never used to
    group or collapse rows -- this function only annotates the gene IDs
    already produced by grouping on ``gene_id_col``.
    """
    valid_mask = _valid_symbol_mask(df[GENE_SYMBOL_COL])
    valid_symbols = df[GENE_SYMBOL_COL].where(valid_mask)

    def _resolve(group: pd.Series) -> tuple[object, str]:
        unique_symbols = group.dropna().unique()
        if len(unique_symbols) == 0:
            return None, SYMBOL_MISSING
        if len(unique_symbols) == 1:
            return unique_symbols[0], SYMBOL_RESOLVED
        return None, SYMBOL_AMBIGUOUS

    resolved = (
        pd.DataFrame({gene_id_col: df[gene_id_col], "_valid_symbol": valid_symbols})
        .groupby(gene_id_col)["_valid_symbol"]
        .apply(_resolve)
    )
    out = pd.DataFrame(resolved.tolist(), index=resolved.index, columns=[GENE_SYMBOL_COL, "symbol_mapping_status"])
    out.index.name = gene_id_col
    return out.reset_index()


@dataclass(frozen=True)
class PreparationQC:
    """Row/gene bookkeeping and reproducibility metadata for the audit log."""

    transcripts_in: int
    gene_ids_out: int
    transcript_to_gene_collapse_count: int
    missing_symbol_count: int
    ambiguous_symbol_count: int
    resolved_symbol_count: int
    tpm_conserved_all_samples: bool
    per_sample_tpm_totals: dict[str, tuple[float, float]]
    source_path: Path | None = None
    source_sha256: str | None = None
    output_path: Path | None = None
    output_sha256: str | None = None

    def log(self) -> None:
        logger.info(
            "gse118713 QC: transcripts_in=%d gene_ids_out=%d collapse_count=%d "
            "missing_symbols=%d ambiguous_symbols=%d resolved_symbols=%d tpm_conserved=%s",
            self.transcripts_in,
            self.gene_ids_out,
            self.transcript_to_gene_collapse_count,
            self.missing_symbol_count,
            self.ambiguous_symbol_count,
            self.resolved_symbol_count,
            self.tpm_conserved_all_samples,
        )

    def as_record(self, sample_ids: tuple[str, ...]) -> dict[str, object]:
        record: dict[str, object] = {
            "transcripts_in": self.transcripts_in,
            "gene_ids_out": self.gene_ids_out,
            "rows_in": self.transcripts_in,
            "rows_out": self.gene_ids_out,
            "rows_lost": 0,
            "transcript_to_gene_collapse_count": self.transcript_to_gene_collapse_count,
            "missing_symbol_count": self.missing_symbol_count,
            "ambiguous_symbol_count": self.ambiguous_symbol_count,
            "resolved_symbol_count": self.resolved_symbol_count,
            "tpm_conserved_all_samples": self.tpm_conserved_all_samples,
            "transformation": (
                "sum transcript TPM by unversioned Ensembl gene ID, per replicate column; "
                "group-mean columns excluded"
            ),
            "source_path": str(self.source_path) if self.source_path else "",
            "source_sha256": self.source_sha256 or "",
            "output_path": str(self.output_path) if self.output_path else "",
            "output_sha256": self.output_sha256 or "",
        }
        for sample_id in sample_ids:
            before, after = self.per_sample_tpm_totals[sample_id]
            record[f"{sample_id}_tpm_total_before"] = before
            record[f"{sample_id}_tpm_total_after"] = after
        return record


def aggregate_to_gene_level(
    df: pd.DataFrame, cfg: Gse118713Config, tolerance: float = 1e-6
) -> tuple[pd.DataFrame, PreparationQC]:
    """Sum transcript TPM to unversioned-Ensembl-gene-ID level, per sample.

    Verifies total TPM is conserved (rows_in sum == rows_out sum, within
    ``tolerance``) for every one of the nine replicate columns before
    returning -- an unnoticed row loss or an accidental mean-instead-of-sum
    aggregation would break this invariant.
    """
    validate_gene_ids(df[GENE_ID_COL])
    working = df.copy()
    working["gene_id"] = strip_ensembl_version(working[GENE_ID_COL])

    sample_cols = list(cfg.sample_columns)
    gene_tpm = working.groupby("gene_id", sort=True)[sample_cols].sum().reset_index()

    symbol_map = reconcile_symbols(working, gene_id_col="gene_id")
    gene_df = gene_tpm.merge(symbol_map, on="gene_id", how="left", validate="one_to_one")
    gene_df = gene_df.rename(columns=cfg.column_to_sample_id)
    gene_df = gene_df[["gene_id", GENE_SYMBOL_COL, "symbol_mapping_status"] + list(cfg.sample_ids)]
    gene_df = gene_df.rename(columns={GENE_SYMBOL_COL: "gene_symbol"})

    per_sample_totals: dict[str, tuple[float, float]] = {}
    not_conserved = []
    for source_col, sample_id in zip(cfg.sample_columns, cfg.sample_ids):
        total_before = float(working[source_col].sum())
        total_after = float(gene_df[sample_id].sum())
        per_sample_totals[sample_id] = (total_before, total_after)
        if not np.isclose(total_before, total_after, rtol=0, atol=tolerance):
            not_conserved.append(sample_id)
    if not_conserved:
        raise ValueError(f"total TPM not conserved after gene-level aggregation for samples: {not_conserved}")

    status_counts = gene_df["symbol_mapping_status"].value_counts()
    qc = PreparationQC(
        transcripts_in=len(df),
        gene_ids_out=len(gene_df),
        transcript_to_gene_collapse_count=len(df) - len(gene_df),
        missing_symbol_count=int(status_counts.get(SYMBOL_MISSING, 0)),
        ambiguous_symbol_count=int(status_counts.get(SYMBOL_AMBIGUOUS, 0)),
        resolved_symbol_count=int(status_counts.get(SYMBOL_RESOLVED, 0)),
        tpm_conserved_all_samples=True,
        per_sample_tpm_totals=per_sample_totals,
    )
    qc.log()
    return gene_df, qc


def write_outputs(gene_df: pd.DataFrame, qc: PreparationQC, cfg: Gse118713Config) -> PreparationQC:
    cfg.gene_tpm_parquet.parent.mkdir(parents=True, exist_ok=True)
    gene_df.to_parquet(cfg.gene_tpm_parquet, index=False)
    output_sha256 = _sha256_file(cfg.gene_tpm_parquet)
    logger.info(
        "write_outputs: wrote %d gene rows to %s (sha256=%s)",
        len(gene_df),
        cfg.gene_tpm_parquet,
        output_sha256,
    )

    qc_with_output = replace(
        qc,
        source_path=cfg.source_path,
        source_sha256=_sha256_file(cfg.source_path),
        output_path=cfg.gene_tpm_parquet,
        output_sha256=output_sha256,
    )
    cfg.qc_tsv.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame([qc_with_output.as_record(cfg.sample_ids)]).to_csv(cfg.qc_tsv, sep="\t", index=False)
    logger.info("write_outputs: wrote QC summary to %s", cfg.qc_tsv)
    return qc_with_output


def prepare_gse118713(config_path: str | Path = "config/config.yaml") -> pd.DataFrame:
    """Run the full GSE118713 gene-level TPM preparation and write outputs."""
    config = _load_config(config_path)
    cfg = Gse118713Config.from_config(config)

    write_sample_metadata(cfg, config)
    df = load_transcript_tpm(cfg)
    gene_df, qc = aggregate_to_gene_level(df, cfg)
    write_outputs(gene_df, qc, cfg)
    return gene_df


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    prepare_gse118713()
