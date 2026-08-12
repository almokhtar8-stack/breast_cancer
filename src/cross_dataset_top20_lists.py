"""Cross-dataset genome-wide integration, Phases 16-17: the genome-wide, candidate-list-independent
global Top 20 and six additional Top-20 lists answering different
biological questions. None is presented as superior to another. Every
list is built by taking the head of an already-deterministically-sorted
table (Phase 6-15's outputs) -- no list is hand-constructed or
name-filtered.

Data source: the tables written by `src.cross_dataset_ranking` and
`src.cross_dataset_consensus_views`.
"""

from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd
import yaml

logger = logging.getLogger(__name__)


def _load_config(config_path: str | Path) -> dict:
    with open(config_path) as f:
        return yaml.safe_load(f)


def _main_caveat(row: pd.Series) -> str:
    """Programmatic, rule-based caveat string -- never hand-typed per gene."""
    caveats = []
    if row.get("n_datasets_testable", 5) == 3:
        caveats.append("borderline coverage (3/5 datasets)")
    if row.get("gse245601_one_track_only", False):
        caveats.append("GSE245601 evidence from a single track only")
    if row.get("gse240112_outlier_fragility", False):
        caveats.append("GSE240112 signal only in the epithelial sensitivity track, not the primary tumor-cell track")
    if row.get("resistance_direction_consensus") == "mixed":
        caveats.append("conflicting direction across resistance-state RNA datasets")
    if row.get("resistance_direction_consensus") == "insufficient":
        caveats.append("insufficient resistance-dataset coverage for a direction consensus")
    return "; ".join(caveats) if caveats else "none"


def build_top20_global(ranked: pd.DataFrame, resistance_consensus: pd.DataFrame, categories: pd.DataFrame, wide: pd.DataFrame) -> pd.DataFrame:
    merged = ranked.merge(resistance_consensus[["gene", "resistance_direction_consensus"]], on="gene", how="left")
    merged = merged.merge(categories, on="gene", how="left")
    merged = merged.merge(wide[["gene", "gse245601_one_track_only", "gse240112_outlier_fragility"]], on="gene", how="left")
    top = merged.head(20).copy()
    top["main_caveat"] = top.apply(_main_caveat, axis=1)
    out = top[
        [
            "global_rank", "gene", "coverage_tier", "crispr_evidence_percentile", "gse118713_evidence_percentile",
            "gse245601_evidence_percentile", "gse240112_evidence_percentile", "gse111151_evidence_percentile",
            "n_datasets_fdr05", "n_datasets_top10pct", "median_evidence_percentile", "resistance_direction_consensus",
            "evidence_category", "main_caveat",
        ]
    ].rename(
        columns={
            "global_rank": "rank", "coverage_tier": "coverage", "crispr_evidence_percentile": "crispr_percentile",
            "gse118713_evidence_percentile": "gse118713_percentile", "gse245601_evidence_percentile": "gse245601_percentile",
            "gse240112_evidence_percentile": "gse240112_percentile", "gse111151_evidence_percentile": "gse111151_percentile",
            "n_datasets_fdr05": "datasets_fdr05", "n_datasets_top10pct": "datasets_top10pct",
            "median_evidence_percentile": "median_percentile", "resistance_direction_consensus": "resistance_consensus",
        }
    )
    logger.info("build_top20_global: %d rows", len(out))
    return out


def build_top20_multimodal(ranked: pd.DataFrame, categories: pd.DataFrame) -> pd.DataFrame:
    merged = ranked.merge(categories, on="gene", how="left")
    out = merged.loc[merged["evidence_category"] == "MULTIMODAL_STRONG"].head(20).reset_index(drop=True)
    logger.info("build_top20_multimodal: %d rows (category may have fewer than 20 genes)", len(out))
    return out


def build_top20_resistance_consensus(resistance_consensus: pd.DataFrame) -> pd.DataFrame:
    return resistance_consensus.head(20).reset_index(drop=True)


def build_top20_rna_only(rna_only: pd.DataFrame) -> pd.DataFrame:
    return rna_only.head(20).reset_index(drop=True)


def build_top20_crispr_direction(crispr_functional: pd.DataFrame, direction: str) -> pd.DataFrame:
    return crispr_functional.loc[crispr_functional["crispr_direction"] == direction].head(20).reset_index(drop=True)


def build_top20_human_tumor(human_only: pd.DataFrame) -> pd.DataFrame:
    return human_only.head(20).reset_index(drop=True)


def build_top20_crispr_nonsignificant_rna(crispr_ns_rna: pd.DataFrame) -> pd.DataFrame:
    return crispr_ns_rna.head(20).reset_index(drop=True)


def run_top20_lists(config_path: str | Path = "config/config.yaml") -> dict[str, pd.DataFrame]:
    config = _load_config(config_path)
    cfg = config["cross_dataset_genomewide"]
    out = cfg["output"]
    tables_dir = Path(out["wide_matrix_tsv"]).parent

    ranked = pd.read_csv(tables_dir / "global_ranking_eligible.tsv", sep="\t")
    wide = pd.read_csv(tables_dir / "all_genes_cross_dataset_evidence_with_ranking.tsv", sep="\t")
    resistance = pd.read_csv(out["resistance_consensus_tsv"], sep="\t")
    categories = pd.read_csv(tables_dir / "evidence_categories.tsv", sep="\t")
    crispr_functional = pd.read_csv(out["crispr_functional_all_genes_tsv"], sep="\t")
    human_only = pd.read_csv(out["human_only_tsv"], sep="\t")
    rna_only = pd.read_csv(out["rna_only_tsv"], sep="\t")
    crispr_ns_rna = pd.read_csv(out["crispr_nonsignificant_rna_consensus_tsv"], sep="\t")

    results = {
        "top20_global": build_top20_global(ranked, resistance, categories, wide),
        "top20_multimodal": build_top20_multimodal(ranked, categories),
        "top20_resistance_consensus": build_top20_resistance_consensus(resistance),
        "top20_rna_only": build_top20_rna_only(rna_only),
        "top20_crispr_sensitising": build_top20_crispr_direction(crispr_functional, "sensitising_KO"),
        "top20_crispr_tolerance": build_top20_crispr_direction(crispr_functional, "tolerance_associated_KO"),
        "top20_human_tumor": build_top20_human_tumor(human_only),
        "top20_crispr_nonsignificant_rna": build_top20_crispr_nonsignificant_rna(crispr_ns_rna),
    }

    path_keys = {
        "top20_global": "top20_global_tsv", "top20_multimodal": "top20_multimodal_tsv",
        "top20_resistance_consensus": "top20_resistance_consensus_tsv", "top20_rna_only": "top20_rna_only_tsv",
        "top20_crispr_sensitising": "top20_crispr_sensitising_tsv", "top20_crispr_tolerance": "top20_crispr_tolerance_tsv",
        "top20_human_tumor": "top20_human_tumor_tsv", "top20_crispr_nonsignificant_rna": "top20_crispr_nonsignificant_rna_tsv",
    }
    for key, df in results.items():
        df.to_csv(out[path_keys[key]], sep="\t", index=False)
    logger.info("run_top20_lists: wrote %d top-20 tables", len(results))
    return results


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run_top20_lists()
