"""Prepared log2 CPM expression matrix for GSE111151 (4 lines, 7 TamR derivatives).

Source: GSE111151_RAW.tar, public 2018-08-30 / last updated 2019-03-27
(config ``data.raw.gse111151_sample_dir``; registry entry in
``config/resources.yaml``). Each of the 11 per-sample files carries three
value columns: raw integer read counts, published (uncorrected) log2 CPM,
and a separately batch-corrected log2 CPM -- see PREANALYSIS.md's
2026-08-05 Phase 2 data-audit amendment. This module uses the published
uncorrected ``CPM (log2)`` column as the prepared expression matrix. Raw
counts are read and validated for QC only, never written to the prepared
output. ``CPM_batch (log2)`` is read only to confirm it is present and
finite, then discarded -- its batch correction was fit for integrating
this cell-line dataset with a separate patient dataset, not for this
cell-line-only comparison, and it never enters the prepared output.

GSE111151 has no biological replicates (n=1 per sample). Tam1/Tam2, where
both exist for a line, are distinct resistant derivative clones of the
same parental line, not replicates of one derivative. The fixed
parental/resistant pairing is read from ``config["gse111151"]["samples"]``,
never re-derived from filenames.

No differential expression, pathway analysis, modelling, feature
construction, or biological interpretation is performed here.
"""

from __future__ import annotations

import hashlib
import logging
from dataclasses import dataclass, replace
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

logger = logging.getLogger(__name__)

GENE_ID_COL = "EnsEMBL_GenID"
GENE_NAME_COL = "gene_name"
COUNTS_COL = "counts"
LOG2CPM_COL = "CPM (log2)"
LOG2CPM_BATCH_COL = "CPM_batch (log2)"
EXPECTED_COLUMNS: tuple[str, ...] = (GENE_ID_COL, GENE_NAME_COL, COUNTS_COL, LOG2CPM_COL, LOG2CPM_BATCH_COL)


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
class SampleSpec:
    gsm: str
    filename: str
    sample_id: str
    parental_line: str
    status: str
    derivative_id: str | None
    paired_parental_sample_id: str | None


@dataclass(frozen=True)
class Gse111151Config:
    """Resolved, config-driven paths and sample layout. No hardcoded paths."""

    sample_dir: Path
    log2cpm_parquet: Path
    sample_metadata_tsv: Path
    qc_tsv: Path
    samples: tuple[SampleSpec, ...]

    @classmethod
    def from_config(cls, config: dict) -> "Gse111151Config":
        section = config["gse111151"]
        samples = tuple(
            SampleSpec(
                gsm=s["gsm"],
                filename=s["filename"],
                sample_id=s["sample_id"],
                parental_line=s["parental_line"],
                status=s["status"],
                derivative_id=s.get("derivative_id"),
                paired_parental_sample_id=s.get("paired_parental_sample_id"),
            )
            for s in section["samples"]
        )
        return cls(
            sample_dir=Path(config["data"]["raw"]["gse111151_sample_dir"]),
            log2cpm_parquet=Path(section["output"]["log2cpm_parquet"]),
            sample_metadata_tsv=Path(section["output"]["sample_metadata_tsv"]),
            qc_tsv=Path(section["output"]["qc_tsv"]),
            samples=samples,
        )

    @property
    def sample_ids(self) -> tuple[str, ...]:
        return tuple(s.sample_id for s in self.samples)


def write_sample_metadata(cfg: Gse111151Config, config: dict) -> pd.DataFrame:
    """Write the eleven-sample metadata table straight from configuration.

    Validates the configured metadata itself is internally consistent
    (exactly 11 samples, exactly 4 parental, exactly 7 resistant, and every
    resistant sample's paired parental actually exists with matching
    parental_line) before writing, so a typo in config.yaml fails loudly
    rather than producing a silently wrong metadata table.
    """
    samples = config["gse111151"]["samples"]
    metadata = pd.DataFrame(samples)[
        ["gsm", "sample_id", "parental_line", "status", "derivative_id", "paired_parental_sample_id"]
    ]
    if len(metadata) != 11:
        raise ValueError(f"expected exactly 11 configured GSE111151 samples, got {len(metadata)}")

    parental = metadata[metadata["status"] == "parental"]
    resistant = metadata[metadata["status"] == "resistant"]
    if len(parental) != 4:
        raise ValueError(f"expected exactly 4 parental GSE111151 samples, got {len(parental)}")
    if len(resistant) != 7:
        raise ValueError(f"expected exactly 7 resistant GSE111151 samples, got {len(resistant)}")

    parental_by_line = dict(zip(parental["parental_line"], parental["sample_id"]))
    bad_pairs = [
        (row.sample_id, row.paired_parental_sample_id)
        for row in resistant.itertuples()
        if parental_by_line.get(row.parental_line) != row.paired_parental_sample_id
    ]
    if bad_pairs:
        raise ValueError(f"resistant/parental pairing inconsistent in config: {bad_pairs}")

    cfg.sample_metadata_tsv.parent.mkdir(parents=True, exist_ok=True)
    metadata.to_csv(cfg.sample_metadata_tsv, sep="\t", index=False)
    logger.info(
        "write_sample_metadata: wrote %d samples (%d parental, %d resistant, %d pairs) to %s",
        len(metadata),
        len(parental),
        len(resistant),
        len(resistant),
        cfg.sample_metadata_tsv,
    )
    return metadata


def read_one_sample(path: str | Path, sample_id: str) -> pd.DataFrame:
    """Read and validate a single GSE111151 sample file.

    Validates the five expected columns are present; ``EnsEMBL_GenID`` is
    nonblank and unique within the file; ``counts`` are finite, non-negative
    integers; and both ``CPM (log2)`` and ``CPM_batch (log2)`` are finite.
    Returns a frame indexed by ``EnsEMBL_GenID`` with all three value
    columns retained (the caller decides which to keep downstream).
    """
    path = Path(path)
    df = pd.read_csv(path, sep="\t")
    missing = [c for c in EXPECTED_COLUMNS if c not in df.columns]
    if missing:
        raise ValueError(f"{path} missing expected columns: {missing}")

    gene_id = df[GENE_ID_COL]
    if gene_id.isna().any() or (gene_id.astype(str).str.strip() == "").any():
        raise ValueError(f"blank or null {GENE_ID_COL} in {path}")
    if gene_id.duplicated().any():
        n_dup = int(gene_id.duplicated().sum())
        raise ValueError(f"{n_dup} duplicate {GENE_ID_COL} values in {path}")

    counts = pd.to_numeric(df[COUNTS_COL], errors="coerce")
    if counts.isna().any() or not np.isfinite(counts.to_numpy(dtype=float)).all():
        raise ValueError(f"non-finite {COUNTS_COL} in {path}")
    if (counts.to_numpy(dtype=float) < 0).any():
        raise ValueError(f"negative {COUNTS_COL} in {path}")
    if not np.array_equal(counts.to_numpy(dtype=float), np.round(counts.to_numpy(dtype=float))):
        raise ValueError(f"non-integer {COUNTS_COL} in {path}")

    log2cpm = pd.to_numeric(df[LOG2CPM_COL], errors="coerce")
    if log2cpm.isna().any() or not np.isfinite(log2cpm.to_numpy(dtype=float)).all():
        raise ValueError(f"non-finite {LOG2CPM_COL} in {path}")

    log2cpm_batch = pd.to_numeric(df[LOG2CPM_BATCH_COL], errors="coerce")
    if log2cpm_batch.isna().any() or not np.isfinite(log2cpm_batch.to_numpy(dtype=float)).all():
        raise ValueError(f"non-finite {LOG2CPM_BATCH_COL} in {path}")

    out = pd.DataFrame(
        {
            GENE_NAME_COL: df[GENE_NAME_COL].to_numpy(),
            COUNTS_COL: counts.astype("int64").to_numpy(),
            LOG2CPM_COL: log2cpm.to_numpy(),
            LOG2CPM_BATCH_COL: log2cpm_batch.to_numpy(),
        },
        index=pd.Index(gene_id.to_numpy(), name=GENE_ID_COL),
    )
    logger.info("read_one_sample: read %d validated rows from %s (sample_id=%s)", len(out), path, sample_id)
    return out


def load_all_samples(cfg: Gse111151Config) -> dict[str, pd.DataFrame]:
    tables: dict[str, pd.DataFrame] = {}
    for spec in cfg.samples:
        path = cfg.sample_dir / spec.filename
        tables[spec.sample_id] = read_one_sample(path, spec.sample_id)
    return tables


def validate_gene_id_consistency(
    tables: dict[str, pd.DataFrame], sample_order: list[str]
) -> tuple[list[str], bool]:
    """Confirm every file shares exactly the same Ensembl gene ID set.

    Row order equality is reported as a QC fact only -- it is checked here
    for information, but ``build_log2cpm_matrix`` always aligns columns by
    identifier (index reindex), never by row position, so a mismatch in
    order alone would not silently misalign the output even though it is
    still surfaced in the QC log.
    """
    if not sample_order:
        raise ValueError("no samples to validate")
    reference_id = sample_order[0]
    reference_index = tables[reference_id].index
    reference_set = set(reference_index)

    mismatched = []
    order_identical = []
    for sample_id in sample_order:
        idx = tables[sample_id].index
        if set(idx) != reference_set:
            mismatched.append(sample_id)
            order_identical.append(False)
        else:
            order_identical.append(bool(np.array_equal(idx.to_numpy(), reference_index.to_numpy())))

    if mismatched:
        raise ValueError(
            f"Ensembl gene ID set mismatch relative to {reference_id} in samples: {mismatched}"
        )

    all_orders_identical = all(order_identical)
    logger.info(
        "validate_gene_id_consistency: %d shared Ensembl gene IDs across %d samples; "
        "identical row order across all files=%s",
        len(reference_set),
        len(sample_order),
        all_orders_identical,
    )
    return sorted(reference_set), all_orders_identical


def validate_symbol_consistency(tables: dict[str, pd.DataFrame], sample_order: list[str]) -> pd.Series:
    """Confirm the Ensembl-ID-to-symbol annotation agrees across every file."""
    reference_id = sample_order[0]
    reference = tables[reference_id][GENE_NAME_COL]
    for sample_id in sample_order[1:]:
        other = tables[sample_id][GENE_NAME_COL].reindex(reference.index)
        if not reference.equals(other):
            raise ValueError(f"gene_name annotation mismatch between {reference_id} and {sample_id}")
    return reference


@dataclass(frozen=True)
class PreparationQC:
    """Row/gene bookkeeping and reproducibility metadata for the audit log."""

    files_read: int
    rows_in: int
    rows_out: int
    shared_gene_id_count: int
    row_order_identical_across_files: bool
    duplicate_symbol_count: int
    counts_column_used: bool
    log2cpm_column_used: bool
    log2cpm_batch_column_used: bool
    source_dir: Path
    per_file_source_sha256: dict[str, str]
    combined_source_sha256: str
    output_path: Path | None = None
    output_sha256: str | None = None

    def log(self) -> None:
        logger.info(
            "gse111151 QC: files_read=%d rows_in=%d rows_out=%d shared_gene_ids=%d "
            "row_order_identical=%s duplicate_symbols=%d counts_used=%s log2cpm_used=%s "
            "log2cpm_batch_used=%s",
            self.files_read,
            self.rows_in,
            self.rows_out,
            self.shared_gene_id_count,
            self.row_order_identical_across_files,
            self.duplicate_symbol_count,
            self.counts_column_used,
            self.log2cpm_column_used,
            self.log2cpm_batch_column_used,
        )

    def as_record(self) -> dict[str, object]:
        record: dict[str, object] = {
            "files_read": self.files_read,
            "rows_in": self.rows_in,
            "rows_out": self.rows_out,
            "rows_lost": self.rows_in - self.rows_out,
            "shared_gene_id_count": self.shared_gene_id_count,
            "row_order_identical_across_files": self.row_order_identical_across_files,
            "duplicate_symbol_count": self.duplicate_symbol_count,
            "counts_column_used_in_output": self.counts_column_used,
            "log2cpm_column_used_in_output": self.log2cpm_column_used,
            "log2cpm_batch_column_used_in_output": self.log2cpm_batch_column_used,
            "transformation": (
                "published uncorrected CPM (log2) per sample, aligned by Ensembl gene ID; "
                "counts read and validated for QC only; CPM_batch (log2) read for "
                "validation only and excluded from output"
            ),
            "source_dir": str(self.source_dir),
            "combined_source_sha256": self.combined_source_sha256,
            "output_path": str(self.output_path) if self.output_path else "",
            "output_sha256": self.output_sha256 or "",
        }
        for sample_id, sha in self.per_file_source_sha256.items():
            record[f"{sample_id}_source_sha256"] = sha
        return record


def build_log2cpm_matrix(
    tables: dict[str, pd.DataFrame], cfg: Gse111151Config
) -> tuple[pd.DataFrame, PreparationQC]:
    """Build the gene_id x sample published-log2-CPM matrix.

    Alignment is always by Ensembl gene ID (index reindex), not by row
    position. Raises if reindexing introduces any missing value, which
    would mean a gene ID present in the reference set is absent from a
    given sample's table -- this should be unreachable once
    ``validate_gene_id_consistency`` has passed, and is checked again here
    as a second, independent guard against silent row loss.
    """
    sample_order = list(cfg.sample_ids)
    reference_gene_ids, order_identical = validate_gene_id_consistency(tables, sample_order)
    gene_names = validate_symbol_consistency(tables, sample_order)

    gene_index = pd.Index(reference_gene_ids, name="gene_id")
    matrix = pd.DataFrame(index=gene_index)
    for sample_id in sample_order:
        col = tables[sample_id][LOG2CPM_COL].reindex(gene_index)
        if col.isna().any():
            raise ValueError(f"missing gene IDs after identifier-based alignment for sample {sample_id}")
        matrix[sample_id] = col

    aligned_symbols = gene_names.reindex(gene_index)
    duplicate_symbol_count = int(aligned_symbols.duplicated(keep=False).sum())

    out = matrix.reset_index()
    out.insert(1, "gene_symbol", aligned_symbols.to_numpy())

    rows_in = int(np.mean([len(t) for t in tables.values()]))
    if any(len(t) != rows_in for t in tables.values()):
        raise ValueError("input files have differing row counts; cannot summarise a single rows_in")

    per_file_sha256 = {
        spec.sample_id: _sha256_file(cfg.sample_dir / spec.filename) for spec in cfg.samples
    }
    combined_source_sha256 = hashlib.sha256(
        "\n".join(f"{sid}:{sha}" for sid, sha in sorted(per_file_sha256.items())).encode("utf-8")
    ).hexdigest()

    qc = PreparationQC(
        files_read=len(tables),
        rows_in=rows_in,
        rows_out=len(out),
        shared_gene_id_count=len(reference_gene_ids),
        row_order_identical_across_files=order_identical,
        duplicate_symbol_count=duplicate_symbol_count,
        counts_column_used=False,
        log2cpm_column_used=True,
        log2cpm_batch_column_used=False,
        source_dir=cfg.sample_dir,
        per_file_source_sha256=per_file_sha256,
        combined_source_sha256=combined_source_sha256,
    )
    qc.log()
    return out, qc


def write_outputs(log2cpm_df: pd.DataFrame, qc: PreparationQC, cfg: Gse111151Config) -> PreparationQC:
    cfg.log2cpm_parquet.parent.mkdir(parents=True, exist_ok=True)
    log2cpm_df.to_parquet(cfg.log2cpm_parquet, index=False)
    output_sha256 = _sha256_file(cfg.log2cpm_parquet)
    logger.info(
        "write_outputs: wrote %d gene rows to %s (sha256=%s)",
        len(log2cpm_df),
        cfg.log2cpm_parquet,
        output_sha256,
    )

    qc_with_output = replace(qc, output_path=cfg.log2cpm_parquet, output_sha256=output_sha256)
    cfg.qc_tsv.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame([qc_with_output.as_record()]).to_csv(cfg.qc_tsv, sep="\t", index=False)
    logger.info("write_outputs: wrote QC summary to %s", cfg.qc_tsv)
    return qc_with_output


def prepare_gse111151(config_path: str | Path = "config/config.yaml") -> pd.DataFrame:
    """Run the full GSE111151 log2 CPM preparation and write outputs."""
    config = _load_config(config_path)
    cfg = Gse111151Config.from_config(config)

    write_sample_metadata(cfg, config)
    tables = load_all_samples(cfg)
    log2cpm_df, qc = build_log2cpm_matrix(tables, cfg)
    write_outputs(log2cpm_df, qc, cfg)
    return log2cpm_df


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    prepare_gse111151()
