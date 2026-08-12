"""Candidate adjudication Phases 5-6: the strongest resistance-state RNA
leaders (independent of CRISPR) and the strongest sensitising-CRISPR
leaders (independent of RNA), each enriched with the other axis's
evidence so the two lists can be compared side by side without merging
them into one score.

Data source: `results/tables/cross_dataset_genomewide/resistance_consensus_all_genes.tsv`,
`crispr_functional_all_genes.tsv`, `all_genes_cross_dataset_evidence_with_ranking.tsv`,
`global_ranking_eligible.tsv`, `evidence_categories.tsv` (all frozen, read-only).
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


def _human_support_label(row: pd.Series) -> str:
    hits = []
    if pd.notna(row.get("gse245601_epi_fdr")) and row["gse245601_epi_fdr"] < 0.05:
        hits.append("GSE245601_epi")
    if pd.notna(row.get("gse240112_tumor_fdr")) and row["gse240112_tumor_fdr"] < 0.05:
        hits.append("GSE240112_tumor")
    return "+".join(hits) if hits else "none_FDR<0.05"


def build_top_resistance_genes(resistance: pd.DataFrame, wide: pd.DataFrame, ranked: pd.DataFrame, categories: pd.DataFrame, top_n: int = 20) -> pd.DataFrame:
    """`resistance_consensus_all_genes.tsv` is already sorted by the
    frozen `build_resistance_consensus` hierarchy -- taking its head is
    the resistance-RNA-only leaderboard by construction, not a re-rank."""
    top = resistance.head(top_n).copy()
    top = top.merge(
        wide[["gene", "gse118713_log2fc", "gse118713_fdr", "gse240112_tumor_log2fc", "gse240112_tumor_fdr", "gse111151_log2fc", "gse111151_fdr", "crispr_effect", "crispr_fdr", "crispr_direction", "gse245601_epi_fdr", "gse240112_epi_log2fc"]],
        on="gene", how="left",
    )
    top = top.merge(ranked[["gene", "global_rank"]], on="gene", how="left")
    top = top.merge(categories, on="gene", how="left")
    top["human_support"] = top.apply(_human_support_label, axis=1)
    out = top[
        [
            "gene", "gse118713_log2fc", "gse118713_fdr", "gse240112_tumor_log2fc", "gse240112_tumor_fdr", "gse111151_log2fc", "gse111151_fdr",
            "resistance_direction_consensus", "resistance_fdr05_count", "resistance_nominal_p05_count", "resistance_median_percentile",
            "crispr_effect", "crispr_fdr", "crispr_direction", "human_support", "global_rank", "evidence_category",
        ]
    ]
    logger.info("build_top_resistance_genes: top %d resistance-RNA leaders", len(out))
    return out


def build_top_crispr_sensitising(crispr_functional: pd.DataFrame, wide: pd.DataFrame, resistance: pd.DataFrame, ranked: pd.DataFrame, categories: pd.DataFrame, top_n: int = 20) -> pd.DataFrame:
    """`crispr_functional_all_genes.tsv` is already sorted fdr/p/gene
    within the full 19,103-gene screen; filtering to `sensitising_KO` and
    taking the head is the CRISPR-only leaderboard, not a re-rank."""
    sensitising = crispr_functional.loc[crispr_functional["crispr_direction"] == "sensitising_KO"].head(top_n).copy()
    sensitising = sensitising.merge(
        wide[["gene", "gse118713_log2fc", "gse118713_fdr", "gse245601_epi_log2fc", "gse245601_epi_fdr", "gse245601_malignant_log2fc", "gse245601_malignant_fdr", "gse240112_tumor_log2fc", "gse240112_tumor_fdr", "gse111151_log2fc", "gse111151_fdr"]],
        on="gene", how="left",
    )
    sensitising = sensitising.merge(resistance[["gene", "resistance_direction_consensus", "resistance_fdr05_count"]], on="gene", how="left")
    sensitising = sensitising.merge(ranked[["gene", "global_rank"]], on="gene", how="left")
    sensitising = sensitising.merge(categories, on="gene", how="left")

    def support_tag(row: pd.Series) -> str:
        has_resistance = pd.notna(row["resistance_fdr05_count"]) and row["resistance_fdr05_count"] >= 1
        has_rna_any = any(pd.notna(row[c]) and row[c] < 0.05 for c in ("gse118713_fdr", "gse245601_epi_fdr", "gse245601_malignant_fdr", "gse240112_tumor_fdr", "gse111151_fdr"))
        if has_resistance:
            return "functional_with_RNA_support"
        if has_rna_any:
            return "functional_with_discordant_or_acute_only_RNA"
        return "functional_only_no_RNA_support"

    sensitising["support_class"] = sensitising.apply(support_tag, axis=1)
    out = sensitising[
        [
            "gene", "crispr_effect", "crispr_p", "crispr_fdr", "crispr_evidence_percentile",
            "gse118713_log2fc", "gse118713_fdr", "gse245601_epi_log2fc", "gse245601_epi_fdr", "gse245601_malignant_log2fc", "gse245601_malignant_fdr",
            "gse240112_tumor_log2fc", "gse240112_tumor_fdr", "gse111151_log2fc", "gse111151_fdr",
            "resistance_direction_consensus", "global_rank", "evidence_category", "support_class",
        ]
    ]
    logger.info("build_top_crispr_sensitising: top %d sensitising-CRISPR leaders", len(out))
    return out


def run_leaders(config_path: str | Path = "config/config.yaml") -> dict[str, pd.DataFrame]:
    config = _load_config(config_path)
    cdx_out = config["cross_dataset_genomewide"]["output"]
    tables_dir = Path(cdx_out["wide_matrix_tsv"]).parent

    resistance = pd.read_csv(cdx_out["resistance_consensus_tsv"], sep="\t")
    wide = pd.read_csv(tables_dir / "all_genes_cross_dataset_evidence_with_ranking.tsv", sep="\t")
    ranked = pd.read_csv(tables_dir / "global_ranking_eligible.tsv", sep="\t")
    categories = pd.read_csv(tables_dir / "evidence_categories.tsv", sep="\t")
    crispr_functional = pd.read_csv(cdx_out["crispr_functional_all_genes_tsv"], sep="\t")

    top_resistance = build_top_resistance_genes(resistance, wide, ranked, categories)
    top_crispr = build_top_crispr_sensitising(crispr_functional, wide, resistance, ranked, categories)

    out_dir = Path(config["candidate_adjudication"]["output"]["tables_dir"])
    out_dir.mkdir(parents=True, exist_ok=True)
    top_resistance.to_csv(out_dir / "top_resistance_genes_exact.tsv", sep="\t", index=False)
    top_crispr.to_csv(out_dir / "top_crispr_sensitising_exact.tsv", sep="\t", index=False)
    return {"top_resistance": top_resistance, "top_crispr": top_crispr}


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run_leaders()
