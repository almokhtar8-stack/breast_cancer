"""Cross-dataset genome-wide integration, Phases 20-21: a human-readable
directional-pattern string per gene (CRISPR functional direction is never
compared to an RNA up/down sign -- they answer different questions, see
docs/CROSS_DATASET_GENOMEWIDE_DATA_AUDIT.md), and a combined statistical
-support-vs-consistency-support summary (Phase 21 explicitly asks these
be kept as two separate, complementary views, never forced into one
score).

Data source: `results/tables/cross_dataset_genomewide/all_genes_cross_dataset_evidence_with_ranking.tsv`,
`resistance_consensus_all_genes.tsv`, `ranking_stability.tsv`.
"""

from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd
import yaml

logger = logging.getLogger(__name__)

_CRISPR_LABEL = {
    "sensitising_KO": "CRISPR sensitising", "tolerance_associated_KO": "CRISPR tolerance-associated",
    "approximately_neutral": "CRISPR neutral", "not_applicable": "CRISPR untested",
}


def _load_config(config_path: str | Path) -> dict:
    with open(config_path) as f:
        return yaml.safe_load(f)


def _arrow(value: float | None) -> str:
    if pd.isna(value):
        return "·"
    if value > 0:
        return "↑"
    if value < 0:
        return "↓"
    return "="


def build_directional_pattern(row: pd.Series) -> str:
    crispr_part = _CRISPR_LABEL.get(row["crispr_direction"], "CRISPR untested")
    resistance_arrows = " ".join(_arrow(row[c]) for c in ("gse118713_log2fc", "gse240112_tumor_log2fc", "gse111151_log2fc"))
    acute_arrow = _arrow(row["gse245601_epi_log2fc"])
    return f"{crispr_part} | resistance RNA (GSE118713,GSE240112,GSE111151) {resistance_arrows} | acute tamoxifen (GSE245601) {acute_arrow}"


def build_directional_patterns_table(df: pd.DataFrame) -> pd.DataFrame:
    out = df[["gene"]].copy()
    out["pattern"] = df.apply(build_directional_pattern, axis=1)
    out["crispr_direction"] = df["crispr_direction"]
    out["gse118713_arrow"] = df["gse118713_log2fc"].map(_arrow)
    out["gse240112_arrow"] = df["gse240112_tumor_log2fc"].map(_arrow)
    out["gse111151_arrow"] = df["gse111151_log2fc"].map(_arrow)
    out["gse245601_arrow"] = df["gse245601_epi_log2fc"].map(_arrow)
    logger.info("build_directional_patterns_table: %d genes", len(out))
    return out


def build_significance_vs_consistency_summary(df: pd.DataFrame, resistance_consensus: pd.DataFrame, stability: pd.DataFrame | None) -> pd.DataFrame:
    """Two complementary, never-combined views per gene:
    STATISTICAL SUPPORT (n_datasets_fdr05, n_datasets_nominal_p05) and
    CONSISTENCY SUPPORT (resistance direction consensus, within-dataset
    track direction agreement where available, median evidence
    percentile, and leave-one-dataset-out stability label)."""
    out = df[["gene", "n_datasets_fdr05", "n_datasets_nominal_p05", "median_evidence_percentile", "gse245601_track_direction_agreement", "gse240112_track_direction_agreement"]].copy()
    out = out.merge(resistance_consensus[["gene", "resistance_direction_consensus", "resistance_fdr05_count"]], on="gene", how="left")
    if stability is not None:
        out = out.merge(stability[["gene", "stability_label", "n_top20_appearances"]], on="gene", how="left")
    logger.info("build_significance_vs_consistency_summary: %d genes", len(out))
    return out


def run_directional_patterns(config_path: str | Path = "config/config.yaml") -> dict[str, pd.DataFrame]:
    config = _load_config(config_path)
    cfg = config["cross_dataset_genomewide"]
    out = cfg["output"]
    tables_dir = Path(out["wide_matrix_tsv"]).parent

    df = pd.read_csv(tables_dir / "all_genes_cross_dataset_evidence_with_ranking.tsv", sep="\t")
    resistance = pd.read_csv(out["resistance_consensus_tsv"], sep="\t")
    stability_path = Path(out["ranking_stability_tsv"])
    stability = pd.read_csv(stability_path, sep="\t") if stability_path.exists() else None

    patterns = build_directional_patterns_table(df)
    summary = build_significance_vs_consistency_summary(df, resistance, stability)

    patterns.to_csv(out["directional_patterns_tsv"], sep="\t", index=False)
    summary_path = tables_dir / "significance_vs_consistency_summary.tsv"
    summary.to_csv(summary_path, sep="\t", index=False)
    logger.info("wrote %s and %s", out["directional_patterns_tsv"], summary_path)

    return {"patterns": patterns, "summary": summary}


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run_directional_patterns()
