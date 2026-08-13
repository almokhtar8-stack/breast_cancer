"""Systems-network phase 2: signed genome-wide ranked gene vectors for GSEA.

Data sources: results/tables/gse118713_differential_expression_unredacted.tsv.gz
(GEO GSE118713, limma), results/tables/gse240112_pseudobulk/tumor_cell_genomewide_de.tsv.gz
(GEO GSE240112, edgeR), results/tables/gse111151/genomewide_de.tsv.gz (GEO
GSE111151, edgeR), results/tables/gse245601_pseudobulk/track_a_genomewide_de.tsv.gz
and track_b_genomewide_de.tsv.gz (GEO GSE245601, edgeR), data/processed/labels.parquet
(CRISPR screen, Hany-style MCF7-V E2+4-OHT vs E2 model fit). Version: as
frozen by their respective upstream phases (see docs/SYSTEMS_NETWORK_INPUT_AUDIT.md).

Ranking statistic per dataset (documented per docs/SYSTEMS_NETWORK_INPUT_AUDIT.md):
GSE118713 uses the limma moderated t-statistic directly (a true model
statistic). GSE240112/GSE111151/GSE245601 (both tracks) use the transparent
fallback sign(log2fc) * -log10(p_value), because their frozen edgeR exports
do not retain an internal test statistic and rerunning edgeR is not
warranted for this. CRISPR uses the Wald statistic effect_size / se.

Duplicate gene symbols (present in GSE118713 and GSE111151 -- see
docs/SYSTEMS_NETWORK_INPUT_AUDIT.md) are resolved deterministically by
keeping, per symbol, the row with the largest |ranking statistic| (ties
broken by smallest p-value, then first-occurrence row order) -- never
averaged, never arbitrarily dropped.
"""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

logger = logging.getLogger(__name__)

DATASET_LABELS = ["gse118713", "gse240112", "gse111151", "gse245601_track_a", "gse245601_track_b", "crispr"]


def _load_config(config_path: str | Path) -> dict:
    with open(config_path) as f:
        return yaml.safe_load(f)


def _dedup_by_max_abs_stat(df: pd.DataFrame, symbol_col: str, stat_col: str, p_col: str) -> pd.DataFrame:
    n_before = len(df)
    df = df.copy()
    df["_abs_stat"] = df[stat_col].abs()
    df = df.sort_values(["_abs_stat", p_col], ascending=[False, True], kind="mergesort")
    df = df.drop_duplicates(subset=symbol_col, keep="first").drop(columns="_abs_stat")
    n_after = len(df)
    if n_after != n_before:
        logger.info("dedup by max|stat|: %d -> %d rows (%d duplicate-symbol rows resolved)", n_before, n_after, n_before - n_after)
    return df


def build_gse118713_ranking(path: str | Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    df = df.loc[df["contrast"] == "TAMR_vs_MCF7"].copy()
    df = _dedup_by_max_abs_stat(df, "gene_symbol", "moderated_t", "p_value")
    out = df.rename(columns={"gene_symbol": "gene", "moderated_t": "ranking_stat"})[["gene", "ranking_stat", "log2fc", "p_value", "fdr"]]
    out["ranking_stat_method"] = "moderated_t"
    logger.info("build_gse118713_ranking: %d genes ranked", len(out))
    return out


def _fallback_signed_stat(log2fc: pd.Series, p_value: pd.Series) -> pd.Series:
    return np.sign(log2fc) * -np.log10(p_value)


def build_edger_ranking(path: str | Path, gene_col: str, dataset_label: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    df["ranking_stat"] = _fallback_signed_stat(df["log2fc"], df["p_value"])
    df = _dedup_by_max_abs_stat(df, gene_col, "ranking_stat", "p_value")
    out = df.rename(columns={gene_col: "gene"})[["gene", "ranking_stat", "log2fc", "p_value", "fdr"]]
    out["ranking_stat_method"] = "sign(log2fc)*-log10(p_value)"
    logger.info("build_edger_ranking(%s): %d genes ranked", dataset_label, len(out))
    return out


def build_crispr_ranking(path: str | Path) -> pd.DataFrame:
    df = pd.read_parquet(path)
    df = df.copy()
    df["ranking_stat"] = df["effect_size"] / df["se"]
    df = _dedup_by_max_abs_stat(df, "gene", "ranking_stat", "p_value")
    out = df.rename(columns={"effect_size": "log2fc"})[["gene", "ranking_stat", "log2fc", "p_value", "fdr"]]
    out["ranking_stat_method"] = "effect_size/se (Wald z)"
    logger.info("build_crispr_ranking: %d genes ranked", len(out))
    return out


def run_ranking(config_path: str | Path = "config/config.yaml") -> dict[str, pd.DataFrame]:
    config = _load_config(config_path)
    cfg = config["systems_network"]
    out_dir = Path(cfg["output"]["tables_dir"])
    out_dir.mkdir(parents=True, exist_ok=True)
    inputs = cfg["inputs"]

    rankings = {
        "gse118713": build_gse118713_ranking(inputs["gse118713_de_tsv"]),
        "gse240112": build_edger_ranking(inputs["gse240112_tumor_cell_tsv"], "gene", "gse240112"),
        "gse111151": build_edger_ranking(inputs["gse111151_de_tsv"], "gene_name", "gse111151"),
        "gse245601_track_a": build_edger_ranking(inputs["gse245601_track_a_tsv"], "gene", "gse245601_track_a"),
        "gse245601_track_b": build_edger_ranking(inputs["gse245601_track_b_tsv"], "gene", "gse245601_track_b"),
        "crispr": build_crispr_ranking(inputs["crispr_labels_parquet"]),
    }

    # gse245601_track_a is the primary/representative GSE245601 ranking and is
    # written under the plain "gse245601_ranked_genes.tsv" name (matching the
    # other three datasets' naming); track_b is written under its own
    # distinct name as the optional exploratory secondary (Phase 2 spec:
    # "Track A and B [are not] independent datasets").
    output_names = {
        "gse118713": "gse118713_ranked_genes.tsv",
        "gse240112": "gse240112_ranked_genes.tsv",
        "gse111151": "gse111151_ranked_genes.tsv",
        "gse245601_track_a": "gse245601_ranked_genes.tsv",
        "gse245601_track_b": "gse245601_track_b_ranked_genes.tsv",
        "crispr": "crispr_ranked_genes.tsv",
    }
    for label, df in rankings.items():
        df.sort_values("ranking_stat", ascending=False).to_csv(out_dir / output_names[label], sep="\t", index=False)

    return rankings


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run_ranking()
