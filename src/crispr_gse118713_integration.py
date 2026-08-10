"""Master integration table: Gate-1 CRISPR hits joined to GSE118713 expression.

Source A: the frozen Hany et al. 2023 gene-by-treatment interaction table
written by ``src.labels.build_labels`` (config ``gate1.labels_path``),
restricted here to the Gate-1 FDR<``gate1.fdr_threshold`` hit set already
decided by ``src.gate1_checks.decide_gate1`` (config
``gate1.decision_output``).

Source B: the frozen GSE118713 gene-level TPM matrix written by
``src.gse118713_prep.prepare_gse118713`` (config
``gse118713.output.gene_tpm_parquet``, checksum-pinned by
``gse118713.frozen_gene_tpm_sha256``), the Phase 2B expression-filtered
subset written by ``src.gse118713_expression_filter`` (config
``gse118713_phase2b.filtering.filtered_gene_tpm_tsv``), and the
TAMR-specificity table written by ``src.gse118713_tamr_specificity``
(config ``gse118713_phase2b.specificity.output_tsv_gz``).

This module performs no differential expression, pathway analysis, or
candidate scoring -- it only joins already-computed, already-frozen
results. Per PREANALYSIS.md's 2026-08-10 amendment, the RCOR1/KDM1A blind
is retired for this analysis, so both genes are treated like any other
Gate-1 hit here; however, their GSE118713 TAMR_vs_MCF7/TAMR_vs_FASR values
were never computed (the limma pipeline redacted them before writing any
output -- see ``scripts/analysis/gse118713_limma_lib.R``'s
``redact_blinded_genes``) and are reported here as
``gse118713_de_na_reason == "blinded_at_source_not_recomputed"`` rather
than silently rerunning limma.

Mapping rule (no silent alias substitution, per PREANALYSIS.md's Phase 2
data-audit amendment): a CRISPR gene symbol is joined to a GSE118713
Ensembl gene ID only when exactly one gene ID in the full (pre-filter)
GSE118713 annotation has ``symbol_mapping_status == "resolved"`` and that
resolved symbol equals the CRISPR gene symbol. A symbol resolved by more
than one distinct gene ID is reported as ``ambiguous`` and excluded from
the RNA fields (not silently resolved by picking one). A symbol resolved
by zero gene IDs is reported as ``unmatched``. Every one of the 28 CRISPR
hits is kept as one master-table row regardless of mapping outcome --
missing RNA fields are ``NA`` with an explicit reason column, never a
dropped row.
"""

from __future__ import annotations

import hashlib
import logging
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from src.gate1_checks import load_and_validate_labels

logger = logging.getLogger(__name__)

MAPPING_UNIQUE_FILTERED = "unique_filtered"
MAPPING_UNIQUE_FILTERED_OUT = "unique_filtered_out"
MAPPING_AMBIGUOUS = "ambiguous"
MAPPING_UNMATCHED = "unmatched"

DE_NA_NOT_MAPPED = "not_uniquely_mapped"
DE_NA_FILTERED_OUT = "filtered_out_no_de_fit"
DE_NA_BLINDED = "blinded_at_source_not_recomputed"

MCF7_REPLICATE_COLUMNS: tuple[str, ...] = ("MCF7_Rep1", "MCF7_Rep2", "MCF7_Rep3")

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


@dataclass(frozen=True)
class IntegrationConfig:
    """Resolved, config-driven paths. No hardcoded paths."""

    labels_path: Path
    fdr_threshold: float
    n_fitted_genes_expected: int
    gate1_decision_tsv: Path
    expected_n_hits: int
    gene_tpm_parquet: Path
    frozen_gene_tpm_sha256: str
    filtered_gene_tpm_tsv: Path
    specificity_tsv_gz: Path
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
        integ = config["crispr_gse118713_integration"]
        return cls(
            labels_path=Path(gate1["labels_path"]),
            fdr_threshold=float(gate1["fdr_threshold"]),
            n_fitted_genes_expected=int(gate1["n_fitted_genes_expected"]),
            gate1_decision_tsv=Path(gate1["decision_output"]),
            expected_n_hits=int(integ["expected_n_hits"]),
            gene_tpm_parquet=Path(gse["output"]["gene_tpm_parquet"]),
            frozen_gene_tpm_sha256=gse["frozen_gene_tpm_sha256"],
            filtered_gene_tpm_tsv=Path(phase2b["filtering"]["filtered_gene_tpm_tsv"]),
            specificity_tsv_gz=Path(phase2b["specificity"]["output_tsv_gz"]),
            master_table_tsv=Path(integ["output"]["master_table_tsv"]),
            master_table_csv=Path(integ["output"]["master_table_csv"]),
            master_table_parquet=Path(integ["output"]["master_table_parquet"]),
            qc_summary_tsv=Path(integ["output"]["qc_summary_tsv"]),
            qc_summary_md=Path(integ["output"]["qc_summary_md"]),
        )


def load_crispr_hits(cfg: IntegrationConfig) -> pd.DataFrame:
    """Load the frozen labels table and extract the Gate-1 FDR<threshold hits.

    Reuses ``load_and_validate_labels`` rather than reimplementing the
    label-table validation, so the hit count here cannot silently drift
    from a malformed or unexpected labels table. Raises if the recomputed
    FDR<threshold count does not match both ``gate1_decision_tsv`` (on
    disk) and ``expected_n_hits`` (config) -- a mismatch means the frozen
    labels table or the gate has changed since this table was
    preregistered as Phase-3 input, which must not pass silently.
    """
    labels_df = load_and_validate_labels(cfg.labels_path, cfg.n_fitted_genes_expected)
    n_passing = int((labels_df["fdr"] < cfg.fdr_threshold).sum())

    recorded = pd.read_csv(cfg.gate1_decision_tsv, sep="\t")
    recorded_n_passing = int(recorded.loc[0, "n_passing"])
    if n_passing != recorded_n_passing:
        raise ValueError(
            f"recomputed Gate-1 hit count ({n_passing}) does not match the committed "
            f"decision record ({recorded_n_passing}) at {cfg.gate1_decision_tsv}"
        )
    if n_passing != cfg.expected_n_hits:
        raise ValueError(
            f"Gate-1 hit count ({n_passing}) does not match config expected_n_hits "
            f"({cfg.expected_n_hits})"
        )

    hits = labels_df.loc[labels_df["fdr"] < cfg.fdr_threshold].copy()
    hits = hits.sort_values(["fdr", "gene"], ascending=[True, True]).reset_index(drop=True)
    logger.info(
        "load_crispr_hits: %d of %d fitted genes pass fdr < %s (matches gate1_decision.tsv and config)",
        len(hits),
        len(labels_df),
        cfg.fdr_threshold,
    )
    return hits


def load_gse118713_annotation(cfg: IntegrationConfig) -> pd.DataFrame:
    """Load and checksum-verify the full (pre-filter) GSE118713 gene table."""
    actual_sha256 = _sha256_file(cfg.gene_tpm_parquet)
    if actual_sha256 != cfg.frozen_gene_tpm_sha256:
        raise ValueError(
            f"GSE118713 gene TPM checksum mismatch: expected {cfg.frozen_gene_tpm_sha256}, "
            f"got {actual_sha256}"
        )
    df = pd.read_parquet(cfg.gene_tpm_parquet)
    missing = [c for c in ("gene_id", "gene_symbol", "symbol_mapping_status") if c not in df.columns]
    if missing:
        raise ValueError(f"GSE118713 annotation missing expected columns: {missing}")
    missing_mcf7 = [c for c in MCF7_REPLICATE_COLUMNS if c not in df.columns]
    if missing_mcf7:
        raise ValueError(f"GSE118713 annotation missing MCF7 replicate columns: {missing_mcf7}")
    logger.info("load_gse118713_annotation: read %d gene rows (checksum verified)", len(df))
    return df


def load_gse118713_filtered_ids(cfg: IntegrationConfig) -> set[str]:
    df = pd.read_csv(cfg.filtered_gene_tpm_tsv, sep="\t")
    if "gene_id" not in df.columns:
        raise ValueError(f"GSE118713 filtered matrix missing 'gene_id' column: {cfg.filtered_gene_tpm_tsv}")
    ids = set(df["gene_id"])
    logger.info("load_gse118713_filtered_ids: %d genes passed the Phase 2B expression filter", len(ids))
    return ids


def load_gse118713_specificity(cfg: IntegrationConfig) -> pd.DataFrame:
    df = pd.read_csv(cfg.specificity_tsv_gz, sep="\t")
    required = {"gene_id", "tamr_vs_mcf7_log2fc", "tamr_vs_mcf7_fdr", "tamr_vs_fasr_log2fc", "tamr_vs_fasr_fdr"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"GSE118713 specificity table missing required columns: {sorted(missing)}")
    logger.info("load_gse118713_specificity: read %d genes (blinded genes already withheld)", len(df))
    return df


def build_symbol_index(annotation_df: pd.DataFrame) -> dict[str, list[str]]:
    """Map each resolved gene symbol to every distinct gene ID that resolves to it.

    Only ``symbol_mapping_status == "resolved"`` rows are indexed -- an
    ambiguous or missing gene ID never contributes a candidate symbol, per
    PREANALYSIS.md's Phase 2 data-audit amendment.
    """
    resolved = annotation_df.loc[annotation_df["symbol_mapping_status"] == "resolved"]
    index: dict[str, list[str]] = {}
    for gene_id, symbol in zip(resolved["gene_id"], resolved["gene_symbol"]):
        index.setdefault(symbol, []).append(gene_id)
    logger.info(
        "build_symbol_index: %d resolved gene IDs indexed under %d distinct symbols",
        len(resolved),
        len(index),
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
) -> pd.DataFrame:
    """Join every Gate-1 hit to GSE118713, keeping exactly one row per hit.

    No CRISPR hit is ever dropped: a mapping failure, filter failure, or
    blinded-at-source DE gap becomes an ``NA`` RNA value plus an explicit
    status/reason column, never a removed row.
    """
    symbol_index = build_symbol_index(annotation_df)
    specificity_by_id = specificity_df.set_index("gene_id")

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
        else:
            # Uniquely mapped and filter-passing, but absent from the
            # specificity table: only RCOR1/KDM1A can land here, since the
            # limma pipeline redacts blinded genes before writing any DE
            # output (scripts/analysis/gse118713_limma_lib.R).
            de_na_reason = DE_NA_BLINDED

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
    n_blinded_gap = int((df["gse118713_de_na_reason"] == DE_NA_BLINDED).sum())

    mcf7_fdr = pd.to_numeric(df["tamr_vs_mcf7_fdr"], errors="coerce")
    fasr_fdr = pd.to_numeric(df["tamr_vs_fasr_fdr"], errors="coerce")
    has_mcf7_sig = mcf7_fdr < 0.05
    has_fasr_sig = fasr_fdr < 0.05
    n_mcf7_sig = int(has_mcf7_sig.fillna(False).sum())
    n_fasr_sig = int(has_fasr_sig.fillna(False).sum())
    n_both_sig = int((has_mcf7_sig.fillna(False) & has_fasr_sig.fillna(False)).sum())

    de_available = df["gse118713_de_na_reason"].isna()
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
        "n_de_gap_blinded_at_source": n_blinded_gap,
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
        f"- Mapped and filter-passing but DE unavailable (RCOR1/KDM1A blinded "
        f"at source, not recomputed): {summary['n_de_gap_blinded_at_source']}",
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

    master_df = build_master_table(hits_df, annotation_df, filtered_ids, specificity_df)
    if len(master_df) != cfg.expected_n_hits:
        raise ValueError(f"expected exactly {cfg.expected_n_hits} master rows, got {len(master_df)}")
    write_master_table(master_df, cfg)

    summary = compute_qc_summary(master_df)
    write_qc_summary(summary, cfg)

    return master_df


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run_integration()
