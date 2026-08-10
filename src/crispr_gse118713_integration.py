"""Master integration table: Gate-1 CRISPR hits joined to GSE118713 expression.

Source A: the frozen Hany et al. 2023 gene-by-treatment interaction table
written by ``src.labels.build_labels`` (config ``gate1.labels_path``,
checksum-pinned by ``gate1.frozen_labels_sha256``), restricted here to the
Gate-1 FDR<``gate1.fdr_threshold`` hit set. The full Gate-1 decision
(``gate1.decision_output``) is reproduced from the frozen labels via
``src.gate1_checks.decide_gate1`` and cross-checked field by field, not
just on ``n_passing``.

Source B: the frozen GSE118713 gene-level TPM matrix written by
``src.gse118713_prep.prepare_gse118713`` (config
``gse118713.output.gene_tpm_parquet``, checksum-pinned by
``gse118713.frozen_gene_tpm_sha256``), the Phase 2B expression-filtered
subset written by ``src.gse118713_expression_filter`` (config
``gse118713_phase2b.filtering.filtered_gene_tpm_tsv``, checksum-pinned by
``gse118713_phase2b.filtering.frozen_filtered_gene_tpm_sha256``), and the
post-unblinding TAMR-specificity table written by
``src.gse118713_unredact.build_unredacted_specificity`` (config
``gse118713_phase2b.unredaction.specificity_unredacted_tsv_gz``,
checksum-pinned by
``gse118713_phase2b.unredaction.frozen_specificity_unredacted_sha256``).

This module performs no differential expression, pathway analysis, or
candidate scoring -- it only joins already-computed, already-frozen
results. Per PREANALYSIS.md's 2026-08-10 amendments (blind retirement and
its provenance correction), RCOR1/KDM1A are treated like any other Gate-1
hit here, and their GSE118713 values come from the post-unblinding
specificity table (``src.gse118713_unredact``), which recovers them from
the exact frozen limma fit with only export-stage redaction disabled --
see that module's docstring for the controlled-unredaction procedure and
its row-for-row provenance comparison against the original redacted
output.

Historical-blind compatibility: if this module is ever pointed at an
older, still-redacted specificity table (config
``gse118713_phase2b.specificity.output_tsv_gz``) instead of the current
post-unblinding one, a uniquely-mapped, filter-passing gene absent from
that table is accepted as an expected gap ONLY when its gene ID is one of
the configured ``gse118713_phase2b.blinding.blinded_gene_ids`` -- any
other unexplained absence raises immediately rather than being assumed to
be a blind. This module never silently infers blinding from a merely
missing row.

Mapping rule (no silent alias substitution, per PREANALYSIS.md's Phase 2
data-audit amendment): a CRISPR gene symbol is joined to a GSE118713
Ensembl gene ID only when exactly one DISTINCT gene ID in the full
(pre-filter) GSE118713 annotation has ``symbol_mapping_status ==
"resolved"`` and that resolved symbol equals the CRISPR gene symbol.
Duplicate annotation rows carrying the identical (gene_id, symbol) pair
collapse to one candidate and do not manufacture false ambiguity; a
symbol resolved by more than one DISTINCT gene ID is reported as
``ambiguous`` and excluded from the RNA fields (not silently resolved by
picking one). A symbol resolved by zero gene IDs is reported as
``unmatched``. Every one of the 28 CRISPR hits is kept as one
master-table row regardless of mapping outcome -- missing RNA fields are
``NA`` with an explicit reason column, never a dropped row.
"""

from __future__ import annotations

import hashlib
import logging
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from src.gate1_checks import decide_gate1, load_and_validate_labels

logger = logging.getLogger(__name__)

MAPPING_UNIQUE_FILTERED = "unique_filtered"
MAPPING_UNIQUE_FILTERED_OUT = "unique_filtered_out"
MAPPING_AMBIGUOUS = "ambiguous"
MAPPING_UNMATCHED = "unmatched"

DE_NA_NOT_MAPPED = "not_uniquely_mapped"
DE_NA_FILTERED_OUT = "filtered_out_no_de_fit"
# Only ever assigned when the missing gene ID is one of the configured
# gse118713_phase2b.blinding.blinded_gene_ids -- see the module docstring's
# "Historical-blind compatibility" note. Not expected to be assigned at all
# when specificity_tsv_gz points at the current post-unblinding table,
# since that table has no blinding-related gaps.
DE_NA_HISTORICALLY_BLINDED = "withheld_from_export_before_blind_retirement"

MCF7_REPLICATE_COLUMNS: tuple[str, ...] = ("MCF7_Rep1", "MCF7_Rep2", "MCF7_Rep3")
ANNOTATION_ID_COLUMNS: tuple[str, ...] = ("gene_id", "gene_symbol", "symbol_mapping_status")

MASTER_TABLE_COLUMNS: tuple[str, ...] = (
    "gene_symbol",
    "crispr_effect_size",
    "crispr_se",
    "crispr_p_value",
    "crispr_fdr",
    "crispr_n_guides",
    "gse118713_gene_id",
    "mapping_status",
    "passed_gse118713_expression_filter",
    "mcf7_baseline_log2_tpm_plus1",
    "tamr_vs_mcf7_log2fc",
    "tamr_vs_mcf7_fdr",
    "tamr_vs_fasr_log2fc",
    "tamr_vs_fasr_fdr",
    "gse118713_de_na_reason",
)


def _load_config(config_path: str | Path) -> dict:
    with open(config_path) as f:
        return yaml.safe_load(f)


def _sha256_file(path: str | Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _verify_checksum(path: str | Path, expected_sha256: str, label: str) -> None:
    actual = _sha256_file(path)
    if actual != expected_sha256:
        raise ValueError(f"{label} checksum mismatch at {path}: expected {expected_sha256}, got {actual}")


@dataclass(frozen=True)
class IntegrationConfig:
    """Resolved, config-driven paths. No hardcoded paths."""

    labels_path: Path
    frozen_labels_sha256: str
    fdr_threshold: float
    n_fitted_genes_expected: int
    classifier_min: int
    continuous_min: int
    gate1_decision_tsv: Path
    expected_n_hits: int
    gene_tpm_parquet: Path
    frozen_gene_tpm_sha256: str
    filtered_gene_tpm_tsv: Path
    frozen_filtered_gene_tpm_sha256: str
    specificity_tsv_gz: Path
    frozen_specificity_sha256: str
    blinded_gene_ids: tuple[str, ...]
    master_table_tsv: Path
    master_table_csv: Path
    master_table_parquet: Path
    qc_summary_tsv: Path
    qc_summary_md: Path

    @classmethod
    def from_config(cls, config: dict) -> "IntegrationConfig":
        gate1 = config["gate1"]
        gse = config["gse118713"]
        phase2b = config["gse118713_phase2b"]
        unred = phase2b["unredaction"]
        integ = config["crispr_gse118713_integration"]
        return cls(
            labels_path=Path(gate1["labels_path"]),
            frozen_labels_sha256=gate1["frozen_labels_sha256"],
            fdr_threshold=float(gate1["fdr_threshold"]),
            n_fitted_genes_expected=int(gate1["n_fitted_genes_expected"]),
            classifier_min=int(gate1["branch_thresholds"]["classifier_min"]),
            continuous_min=int(gate1["branch_thresholds"]["continuous_min"]),
            gate1_decision_tsv=Path(gate1["decision_output"]),
            expected_n_hits=int(integ["expected_n_hits"]),
            gene_tpm_parquet=Path(gse["output"]["gene_tpm_parquet"]),
            frozen_gene_tpm_sha256=gse["frozen_gene_tpm_sha256"],
            filtered_gene_tpm_tsv=Path(phase2b["filtering"]["filtered_gene_tpm_tsv"]),
            frozen_filtered_gene_tpm_sha256=phase2b["filtering"]["frozen_filtered_gene_tpm_sha256"],
            specificity_tsv_gz=Path(unred["specificity_unredacted_tsv_gz"]),
            frozen_specificity_sha256=unred["frozen_specificity_unredacted_sha256"],
            blinded_gene_ids=tuple(phase2b["blinding"]["blinded_gene_ids"]),
            master_table_tsv=Path(integ["output"]["master_table_tsv"]),
            master_table_csv=Path(integ["output"]["master_table_csv"]),
            master_table_parquet=Path(integ["output"]["master_table_parquet"]),
            qc_summary_tsv=Path(integ["output"]["qc_summary_tsv"]),
            qc_summary_md=Path(integ["output"]["qc_summary_md"]),
        )


def load_crispr_hits(cfg: IntegrationConfig) -> pd.DataFrame:
    """Load the frozen labels table and extract the Gate-1 FDR<threshold hits.

    Verifies the labels file's checksum, reuses ``load_and_validate_labels``
    for structural validation, and reproduces the FULL Gate-1 decision via
    ``decide_gate1`` -- not just ``n_passing`` -- cross-checking every
    field recorded in ``gate1_decision.tsv`` (total fitted genes, FDR
    threshold, threshold band, branch decision, labels input path, labels
    checksum) against what is actually recomputed from the frozen labels
    file. Any mismatch means the frozen labels table, the gate, or the
    decision record have drifted apart since this integration was
    preregistered, which must not pass silently.
    """
    _verify_checksum(cfg.labels_path, cfg.frozen_labels_sha256, "labels input")

    labels_df = load_and_validate_labels(cfg.labels_path, cfg.n_fitted_genes_expected)
    decision = decide_gate1(
        labels_df,
        fdr_threshold=cfg.fdr_threshold,
        classifier_min=cfg.classifier_min,
        continuous_min=cfg.continuous_min,
    )

    recorded = pd.read_csv(cfg.gate1_decision_tsv, sep="\t")
    if len(recorded) != 1:
        raise ValueError(f"expected exactly one row in {cfg.gate1_decision_tsv}, got {len(recorded)}")
    recorded_row = recorded.iloc[0]

    mismatches: dict[str, tuple[object, object]] = {}
    for field in ("total_fitted_genes", "fdr_threshold", "n_passing", "threshold_band", "branch_decision"):
        recomputed_value = decision[field]
        recorded_value = recorded_row[field]
        if isinstance(recomputed_value, float):
            equal = np.isclose(float(recomputed_value), float(recorded_value), rtol=0, atol=1e-12)
        else:
            equal = recomputed_value == recorded_value
        if not equal:
            mismatches[field] = (recomputed_value, recorded_value)

    recorded_labels_path = str(recorded_row["labels_input_path"])
    if recorded_labels_path != str(cfg.labels_path):
        mismatches["labels_input_path"] = (str(cfg.labels_path), recorded_labels_path)

    recorded_labels_sha256 = str(recorded_row["labels_file_sha256"])
    if recorded_labels_sha256 != cfg.frozen_labels_sha256:
        mismatches["labels_file_sha256"] = (cfg.frozen_labels_sha256, recorded_labels_sha256)

    if mismatches:
        raise ValueError(
            f"recomputed Gate-1 decision does not match the committed decision record "
            f"at {cfg.gate1_decision_tsv}: {mismatches}"
        )
    if decision["n_passing"] != cfg.expected_n_hits:
        raise ValueError(
            f"Gate-1 hit count ({decision['n_passing']}) does not match config "
            f"expected_n_hits ({cfg.expected_n_hits})"
        )

    hits = labels_df.loc[labels_df["fdr"] < cfg.fdr_threshold].copy()
    hits = hits.sort_values(["fdr", "gene"], ascending=[True, True]).reset_index(drop=True)
    logger.info(
        "load_crispr_hits: %d of %d fitted genes pass fdr < %s (full Gate-1 decision record verified)",
        len(hits),
        len(labels_df),
        cfg.fdr_threshold,
    )
    return hits


def load_gse118713_annotation(cfg: IntegrationConfig) -> pd.DataFrame:
    """Load and checksum-verify the full (pre-filter) GSE118713 gene table.

    Validates: required columns present; every TPM sample column finite
    and non-negative (required before ``log2(TPM + 1)`` is ever taken).

    Duplicate row handling: a ``gene_id`` duplicated across rows that are
    otherwise byte-identical is a harmless duplicate (e.g. a repeated
    annotation row) and is silently collapsed to one row. A ``gene_id``
    duplicated across rows that DISAGREE on any other column is a
    genuine data conflict and raises -- it is never guessed at or
    silently resolved by picking one row.
    """
    _verify_checksum(cfg.gene_tpm_parquet, cfg.frozen_gene_tpm_sha256, "GSE118713 full annotation")
    df = pd.read_parquet(cfg.gene_tpm_parquet)

    missing = [c for c in ANNOTATION_ID_COLUMNS if c not in df.columns]
    if missing:
        raise ValueError(f"GSE118713 annotation missing expected columns: {missing}")

    rows_before = len(df)
    df = df.drop_duplicates().reset_index(drop=True)
    if len(df) != rows_before:
        logger.info(
            "load_gse118713_annotation: collapsed %d exact-duplicate row(s)", rows_before - len(df)
        )
    if df["gene_id"].duplicated().any():
        conflicting = sorted(df.loc[df["gene_id"].duplicated(keep=False), "gene_id"].unique())
        raise ValueError(
            f"{len(conflicting)} gene_id(s) have conflicting (non-identical) duplicate rows in "
            f"GSE118713 annotation: {conflicting[:10]}{'...' if len(conflicting) > 10 else ''}"
        )

    missing_mcf7 = [c for c in MCF7_REPLICATE_COLUMNS if c not in df.columns]
    if missing_mcf7:
        raise ValueError(f"GSE118713 annotation missing MCF7 replicate columns: {missing_mcf7}")

    sample_cols = [c for c in df.columns if c not in ANNOTATION_ID_COLUMNS]
    tpm = df[sample_cols].to_numpy(dtype=float)
    if not np.isfinite(tpm).all():
        raise ValueError("non-finite TPM value in GSE118713 annotation")
    if (tpm < 0).any():
        raise ValueError("negative TPM value in GSE118713 annotation")

    logger.info("load_gse118713_annotation: read %d gene rows (checksum + uniqueness + finiteness verified)", len(df))
    return df


def load_gse118713_filtered_ids(cfg: IntegrationConfig) -> set[str]:
    _verify_checksum(cfg.filtered_gene_tpm_tsv, cfg.frozen_filtered_gene_tpm_sha256, "GSE118713 filtered matrix")
    df = pd.read_csv(cfg.filtered_gene_tpm_tsv, sep="\t")
    if "gene_id" not in df.columns:
        raise ValueError(f"GSE118713 filtered matrix missing 'gene_id' column: {cfg.filtered_gene_tpm_tsv}")
    if df["gene_id"].duplicated().any():
        raise ValueError(f"{int(df['gene_id'].duplicated().sum())} duplicate gene_id values in GSE118713 filtered matrix")
    ids = set(df["gene_id"])
    logger.info("load_gse118713_filtered_ids: %d genes passed the Phase 2B expression filter (checksum verified)", len(ids))
    return ids


def load_gse118713_specificity(cfg: IntegrationConfig) -> pd.DataFrame:
    """Load and validate the TAMR-specificity table used as the DE source.

    Validates: required columns; unique ``gene_id``; every present logFC/
    FDR value finite; every present FDR in ``[0, 1]``.
    """
    _verify_checksum(cfg.specificity_tsv_gz, cfg.frozen_specificity_sha256, "GSE118713 specificity table")
    df = pd.read_csv(cfg.specificity_tsv_gz, sep="\t")
    required = {"gene_id", "tamr_vs_mcf7_log2fc", "tamr_vs_mcf7_fdr", "tamr_vs_fasr_log2fc", "tamr_vs_fasr_fdr"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"GSE118713 specificity table missing required columns: {sorted(missing)}")
    if df["gene_id"].duplicated().any():
        raise ValueError(f"{int(df['gene_id'].duplicated().sum())} duplicate gene_id values in GSE118713 specificity table")

    numeric_cols = ["tamr_vs_mcf7_log2fc", "tamr_vs_mcf7_fdr", "tamr_vs_fasr_log2fc", "tamr_vs_fasr_fdr"]
    values = df[numeric_cols].to_numpy(dtype=float)
    if not np.isfinite(values).all():
        raise ValueError("non-finite logFC/FDR value in GSE118713 specificity table")

    fdr_cols = df[["tamr_vs_mcf7_fdr", "tamr_vs_fasr_fdr"]].to_numpy(dtype=float)
    if ((fdr_cols < 0) | (fdr_cols > 1)).any():
        raise ValueError("FDR value outside [0, 1] in GSE118713 specificity table")

    logger.info(
        "load_gse118713_specificity: read %d genes from %s (checksum + uniqueness + range verified)",
        len(df),
        cfg.specificity_tsv_gz,
    )
    return df


def validate_filtered_ids_in_annotation(filtered_ids: set[str], annotation_df: pd.DataFrame) -> None:
    """Every expression-filtered gene ID must exist in the full annotation.

    Checked globally over the whole filtered universe, not only the 28
    CRISPR hits actually used downstream -- a gap here means the two
    frozen GSE118713 tables have drifted apart from each other.
    """
    annotation_ids = set(annotation_df["gene_id"])
    missing = filtered_ids - annotation_ids
    if missing:
        raise ValueError(
            f"{len(missing)} filtered GSE118713 gene(s) absent from the full annotation: "
            f"{sorted(missing)[:10]}{'...' if len(missing) > 10 else ''}"
        )
    logger.info("validate_filtered_ids_in_annotation: all %d filtered gene IDs present in annotation", len(filtered_ids))


def validate_specificity_covers_filtered_universe(specificity_df: pd.DataFrame, filtered_ids: set[str]) -> None:
    """The specificity table's gene set must exactly equal the filtered universe.

    Checked globally, not only the 28 CRISPR hits: a missing gene means a
    filtered gene has no DE result for a reason other than a configured
    historical blind (which ``build_master_table`` also refuses to
    silently assume); an unexpected extra gene means the specificity
    table includes something outside the frozen filtered matrix used to
    fit it.
    """
    specificity_ids = set(specificity_df["gene_id"])
    missing = filtered_ids - specificity_ids
    extra = specificity_ids - filtered_ids
    if missing:
        raise ValueError(
            f"{len(missing)} filtered GSE118713 gene(s) missing from the specificity table: "
            f"{sorted(missing)[:10]}{'...' if len(missing) > 10 else ''}"
        )
    if extra:
        raise ValueError(
            f"{len(extra)} gene(s) in the specificity table are not in the filtered universe: "
            f"{sorted(extra)[:10]}{'...' if len(extra) > 10 else ''}"
        )
    logger.info(
        "validate_specificity_covers_filtered_universe: %d genes exactly match the filtered universe",
        len(specificity_ids),
    )


def build_symbol_index(annotation_df: pd.DataFrame) -> dict[str, list[str]]:
    """Map each resolved gene symbol to every DISTINCT gene ID that resolves to it.

    Only ``symbol_mapping_status == "resolved"`` rows are indexed -- an
    ambiguous or missing gene ID never contributes a candidate symbol, per
    PREANALYSIS.md's Phase 2 data-audit amendment. Candidates are built
    from distinct (gene_id, symbol) pairs, so a duplicated annotation row
    carrying the identical pair twice does not manufacture false
    ambiguity; only genuinely distinct gene IDs sharing a symbol count as
    more than one candidate.
    """
    resolved = annotation_df.loc[annotation_df["symbol_mapping_status"] == "resolved", ["gene_id", "gene_symbol"]]
    distinct_pairs = resolved.drop_duplicates()

    index: dict[str, list[str]] = {}
    for gene_id, symbol in zip(distinct_pairs["gene_id"], distinct_pairs["gene_symbol"]):
        index.setdefault(symbol, []).append(gene_id)
    logger.info(
        "build_symbol_index: %d distinct (gene_id, symbol) pairs indexed under %d distinct symbols "
        "(%d raw resolved rows before dedup)",
        len(distinct_pairs),
        len(index),
        len(resolved),
    )
    return index


def _mcf7_baseline(annotation_df: pd.DataFrame, gene_id: str) -> float:
    row = annotation_df.loc[annotation_df["gene_id"] == gene_id, list(MCF7_REPLICATE_COLUMNS)]
    if len(row) != 1:
        raise ValueError(f"expected exactly one annotation row for gene_id={gene_id!r}, found {len(row)}")
    tpm = row.iloc[0].to_numpy(dtype=float)
    return float(np.mean(np.log2(tpm + 1.0)))


def build_master_table(
    hits_df: pd.DataFrame,
    annotation_df: pd.DataFrame,
    filtered_ids: set[str],
    specificity_df: pd.DataFrame,
    blinded_gene_ids: tuple[str, ...] = (),
) -> pd.DataFrame:
    """Join every Gate-1 hit to GSE118713, keeping exactly one row per hit.

    No CRISPR hit is ever dropped: a mapping failure or filter failure
    becomes an ``NA`` RNA value plus an explicit status/reason column,
    never a removed row. A uniquely mapped, filter-passing gene absent
    from ``specificity_df`` is NOT assumed to be blinded: it raises
    immediately unless its gene ID is in ``blinded_gene_ids`` (see the
    module docstring's "Historical-blind compatibility" note).
    """
    symbol_index = build_symbol_index(annotation_df)
    specificity_by_id = specificity_df.set_index("gene_id")
    blinded_gene_id_set = set(blinded_gene_ids)

    rows: list[dict[str, object]] = []
    for record in hits_df.itertuples(index=False):
        symbol = record.gene
        candidate_ids = symbol_index.get(symbol, [])

        gene_id: str | None
        mapping_status: str
        passed_filter: object = pd.NA
        mcf7_baseline: object = pd.NA

        if len(candidate_ids) == 0:
            gene_id = None
            mapping_status = MAPPING_UNMATCHED
        elif len(candidate_ids) > 1:
            gene_id = None
            mapping_status = MAPPING_AMBIGUOUS
        else:
            gene_id = candidate_ids[0]
            passed_filter = gene_id in filtered_ids
            mapping_status = MAPPING_UNIQUE_FILTERED if passed_filter else MAPPING_UNIQUE_FILTERED_OUT
            mcf7_baseline = _mcf7_baseline(annotation_df, gene_id)

        tamr_vs_mcf7_log2fc: object = pd.NA
        tamr_vs_mcf7_fdr: object = pd.NA
        tamr_vs_fasr_log2fc: object = pd.NA
        tamr_vs_fasr_fdr: object = pd.NA
        de_na_reason: object = pd.NA

        if mapping_status in (MAPPING_UNMATCHED, MAPPING_AMBIGUOUS):
            de_na_reason = DE_NA_NOT_MAPPED
        elif mapping_status == MAPPING_UNIQUE_FILTERED_OUT:
            de_na_reason = DE_NA_FILTERED_OUT
        elif gene_id in specificity_by_id.index:
            de_row = specificity_by_id.loc[gene_id]
            tamr_vs_mcf7_log2fc = float(de_row["tamr_vs_mcf7_log2fc"])
            tamr_vs_mcf7_fdr = float(de_row["tamr_vs_mcf7_fdr"])
            tamr_vs_fasr_log2fc = float(de_row["tamr_vs_fasr_log2fc"])
            tamr_vs_fasr_fdr = float(de_row["tamr_vs_fasr_fdr"])
        elif gene_id in blinded_gene_id_set:
            de_na_reason = DE_NA_HISTORICALLY_BLINDED
        else:
            raise ValueError(
                f"gene_id={gene_id!r} (symbol={symbol!r}) is uniquely mapped and passed the "
                f"expression filter, but is absent from the specificity table and is not a "
                f"configured historical blind ID -- this is an unexplained gap and must not be "
                f"silently assumed to be blinding"
            )

        rows.append(
            {
                "gene_symbol": symbol,
                "crispr_effect_size": float(record.effect_size),
                "crispr_se": float(record.se),
                "crispr_p_value": float(record.p_value),
                "crispr_fdr": float(record.fdr),
                "crispr_n_guides": int(record.n_guides),
                "gse118713_gene_id": gene_id if gene_id is not None else pd.NA,
                "mapping_status": mapping_status,
                "passed_gse118713_expression_filter": passed_filter,
                "mcf7_baseline_log2_tpm_plus1": mcf7_baseline,
                "tamr_vs_mcf7_log2fc": tamr_vs_mcf7_log2fc,
                "tamr_vs_mcf7_fdr": tamr_vs_mcf7_fdr,
                "tamr_vs_fasr_log2fc": tamr_vs_fasr_log2fc,
                "tamr_vs_fasr_fdr": tamr_vs_fasr_fdr,
                "gse118713_de_na_reason": de_na_reason,
            }
        )

    out = pd.DataFrame(rows, columns=list(MASTER_TABLE_COLUMNS))

    mapping_counts = out["mapping_status"].value_counts().to_dict()
    logger.info(
        "build_master_table: %d CRISPR hits in, %d master rows out (0 dropped); mapping_status=%s",
        len(hits_df),
        len(out),
        mapping_counts,
    )
    return out


def write_master_table(df: pd.DataFrame, cfg: IntegrationConfig) -> None:
    for path in (cfg.master_table_tsv, cfg.master_table_csv, cfg.master_table_parquet):
        path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(cfg.master_table_tsv, sep="\t", index=False)
    df.to_csv(cfg.master_table_csv, index=False)
    df.to_parquet(cfg.master_table_parquet, index=False)
    logger.info(
        "write_master_table: wrote %d rows to %s, %s, %s",
        len(df),
        cfg.master_table_tsv,
        cfg.master_table_csv,
        cfg.master_table_parquet,
    )


def compute_qc_summary(df: pd.DataFrame) -> dict[str, object]:
    """Compute the aggregate QC counts required for the integration report.

    Every count here is a summary statistic over the 28-row master table --
    no individual gene identity is required to interpret this summary.
    """
    n_total = len(df)
    n_unique_filtered = int((df["mapping_status"] == MAPPING_UNIQUE_FILTERED).sum())
    n_unique_filtered_out = int((df["mapping_status"] == MAPPING_UNIQUE_FILTERED_OUT).sum())
    n_ambiguous = int((df["mapping_status"] == MAPPING_AMBIGUOUS).sum())
    n_unmatched = int((df["mapping_status"] == MAPPING_UNMATCHED).sum())
    n_historically_blinded_gap = int((df["gse118713_de_na_reason"] == DE_NA_HISTORICALLY_BLINDED).sum())

    mcf7_fdr = pd.to_numeric(df["tamr_vs_mcf7_fdr"], errors="coerce")
    fasr_fdr = pd.to_numeric(df["tamr_vs_fasr_fdr"], errors="coerce")
    has_mcf7_sig = mcf7_fdr < 0.05
    has_fasr_sig = fasr_fdr < 0.05
    n_mcf7_sig = int(has_mcf7_sig.fillna(False).sum())
    n_fasr_sig = int(has_fasr_sig.fillna(False).sum())
    n_both_sig = int((has_mcf7_sig.fillna(False) & has_fasr_sig.fillna(False)).sum())

    de_available = df["gse118713_de_na_reason"].isna()
    n_de_available = int(de_available.sum())
    log2fc_mcf7 = pd.to_numeric(df.loc[de_available, "tamr_vs_mcf7_log2fc"], errors="coerce")
    log2fc_fasr = pd.to_numeric(df.loc[de_available, "tamr_vs_fasr_log2fc"], errors="coerce")
    n_up_mcf7 = int((log2fc_mcf7 > 0).sum())
    n_down_mcf7 = int((log2fc_mcf7 < 0).sum())
    n_up_fasr = int((log2fc_fasr > 0).sum())
    n_down_fasr = int((log2fc_fasr < 0).sum())

    summary = {
        "n_crispr_hits_total": n_total,
        "n_mapped_unique_and_filtered": n_unique_filtered,
        "n_mapped_unique_but_filtered_out": n_unique_filtered_out,
        "n_mapping_ambiguous": n_ambiguous,
        "n_mapping_unmatched": n_unmatched,
        "n_de_gap_historically_blinded": n_historically_blinded_gap,
        "n_de_available": n_de_available,
        "n_tamr_vs_mcf7_fdr_lt_0_05": n_mcf7_sig,
        "n_tamr_vs_fasr_fdr_lt_0_05": n_fasr_sig,
        "n_significant_in_both_comparisons": n_both_sig,
        "n_de_available_up_in_tamr_vs_mcf7": n_up_mcf7,
        "n_de_available_down_in_tamr_vs_mcf7": n_down_mcf7,
        "n_de_available_up_in_tamr_vs_fasr": n_up_fasr,
        "n_de_available_down_in_tamr_vs_fasr": n_down_fasr,
    }
    logger.info("compute_qc_summary: %s", summary)
    return summary


def write_qc_summary(summary: dict[str, object], cfg: IntegrationConfig) -> None:
    cfg.qc_summary_tsv.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame([summary]).to_csv(cfg.qc_summary_tsv, sep="\t", index=False)

    cfg.qc_summary_md.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "# CRISPR-to-GSE118713 master table QC summary",
        "",
        "Descriptive, aggregate counts only. FDR values here are evidence for",
        "candidate interpretation, not a predictive feature or a ranking rule.",
        "",
        f"- CRISPR hits (Gate-1 FDR<0.1): {summary['n_crispr_hits_total']}",
        f"- Uniquely mapped and passed the GSE118713 expression filter: "
        f"{summary['n_mapped_unique_and_filtered']}",
        f"- Uniquely mapped but removed by the GSE118713 expression filter: "
        f"{summary['n_mapped_unique_but_filtered_out']}",
        f"- Ambiguous mapping (excluded from RNA fields): {summary['n_mapping_ambiguous']}",
        f"- Unmatched (no GSE118713 gene ID resolves to this symbol): "
        f"{summary['n_mapping_unmatched']}",
        f"- Mapped and filter-passing but DE unavailable (historical blind gap, "
        f"expected to be zero when reading the post-unblinding specificity table): "
        f"{summary['n_de_gap_historically_blinded']}",
        f"- DE values available: {summary['n_de_available']}",
        f"- TAMR_vs_MCF7 FDR<0.05: {summary['n_tamr_vs_mcf7_fdr_lt_0_05']}",
        f"- TAMR_vs_FASR FDR<0.05: {summary['n_tamr_vs_fasr_fdr_lt_0_05']}",
        f"- Significant in both comparisons: {summary['n_significant_in_both_comparisons']}",
        f"- Directionality where DE is available: "
        f"{summary['n_de_available_up_in_tamr_vs_mcf7']} up / "
        f"{summary['n_de_available_down_in_tamr_vs_mcf7']} down in TAMR_vs_MCF7; "
        f"{summary['n_de_available_up_in_tamr_vs_fasr']} up / "
        f"{summary['n_de_available_down_in_tamr_vs_fasr']} down in TAMR_vs_FASR. "
        "Direction is descriptive only and not a causal claim.",
    ]
    cfg.qc_summary_md.write_text("\n".join(lines) + "\n")
    logger.info("write_qc_summary: wrote %s and %s", cfg.qc_summary_tsv, cfg.qc_summary_md)


def run_integration(config_path: str | Path = "config/config.yaml") -> pd.DataFrame:
    """Run the full CRISPR-to-GSE118713 integration and write all outputs."""
    config = _load_config(config_path)
    cfg = IntegrationConfig.from_config(config)

    hits_df = load_crispr_hits(cfg)
    annotation_df = load_gse118713_annotation(cfg)
    filtered_ids = load_gse118713_filtered_ids(cfg)
    specificity_df = load_gse118713_specificity(cfg)

    validate_filtered_ids_in_annotation(filtered_ids, annotation_df)
    validate_specificity_covers_filtered_universe(specificity_df, filtered_ids)

    master_df = build_master_table(hits_df, annotation_df, filtered_ids, specificity_df, cfg.blinded_gene_ids)
    if len(master_df) != cfg.expected_n_hits:
        raise ValueError(f"expected exactly {cfg.expected_n_hits} master rows, got {len(master_df)}")
    write_master_table(master_df, cfg)

    summary = compute_qc_summary(master_df)
    write_qc_summary(summary, cfg)

    return master_df


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run_integration()
