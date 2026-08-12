"""Candidate adjudication Phase 19: quantifies how much the global
ranking is driven by RNA evidence vs. CRISPR functional evidence vs.
resistance-state RNA specifically, via rank correlation and Top-20/Top-100
overlap -- so the global ranking is never mistaken for a functional
-therapeutic ranking. Purely descriptive; does not change any ranking.

Data source: `results/tables/cross_dataset_genomewide/global_ranking_eligible.tsv`,
`rna_only_ranking.tsv`, `crispr_functional_all_genes.tsv`,
`resistance_consensus_all_genes.tsv` (all frozen, read-only).
"""

from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd
import yaml
from scipy.stats import spearmanr

logger = logging.getLogger(__name__)


def _load_config(config_path: str | Path) -> dict:
    with open(config_path) as f:
        return yaml.safe_load(f)


def build_modality_dependency_summary(global_ranked: pd.DataFrame, rna_only: pd.DataFrame, crispr_functional: pd.DataFrame, resistance: pd.DataFrame) -> pd.DataFrame:
    crispr_ranked = crispr_functional.reset_index(drop=True).copy()
    crispr_ranked["crispr_rank"] = crispr_ranked.index + 1
    resistance_ranked = resistance.reset_index(drop=True).copy()
    resistance_ranked["resistance_rank"] = resistance_ranked.index + 1

    merged = global_ranked.merge(rna_only[["gene", "rank"]].rename(columns={"rank": "rna_only_rank"}), on="gene", how="left")
    merged = merged.merge(crispr_ranked[["gene", "crispr_rank"]], on="gene", how="left")
    merged = merged.merge(resistance_ranked[["gene", "resistance_rank"]], on="gene", how="left")

    rows = []
    for other_col, label in [("rna_only_rank", "RNA-only"), ("crispr_rank", "CRISPR-only"), ("resistance_rank", "resistance-consensus-only")]:
        paired = merged.dropna(subset=["global_rank", other_col])
        rho, p = spearmanr(paired["global_rank"], paired[other_col])
        rows.append({"comparison": f"global_rank_vs_{label}", "n_genes_compared": len(paired), "spearman_rho": rho, "p_value": p})
    corr_table = pd.DataFrame(rows)
    logger.info("build_modality_dependency_summary: %s", {r["comparison"]: round(r["spearman_rho"], 3) for _, r in corr_table.iterrows()})
    return corr_table, merged


def build_topn_overlap_table(merged: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for n in (20, 100):
        global_top = set(merged.loc[merged["global_rank"] <= n, "gene"])
        for other_col, label in [("rna_only_rank", "RNA-only"), ("crispr_rank", "CRISPR-sensitising-and-tolerance"), ("resistance_rank", "resistance-consensus")]:
            other_top = set(merged.loc[merged[other_col] <= n, "gene"])
            rows.append({"top_n": n, "comparison_ranking": label, "overlap_with_global": len(global_top & other_top), "overlap_fraction_of_global": len(global_top & other_top) / n})
    out = pd.DataFrame(rows)
    logger.info("build_topn_overlap_table: %s", out.to_dict("records"))
    return out


def run_modality_dependency(config_path: str | Path = "config/config.yaml") -> dict[str, pd.DataFrame]:
    config = _load_config(config_path)
    cdx_out = config["cross_dataset_genomewide"]["output"]
    tables_dir = Path(cdx_out["wide_matrix_tsv"]).parent

    global_ranked = pd.read_csv(tables_dir / "global_ranking_eligible.tsv", sep="\t")[["gene", "global_rank"]]
    rna_only = pd.read_csv(cdx_out["rna_only_tsv"], sep="\t")
    crispr_functional = pd.read_csv(cdx_out["crispr_functional_all_genes_tsv"], sep="\t")
    resistance = pd.read_csv(cdx_out["resistance_consensus_tsv"], sep="\t")

    corr_table, merged = build_modality_dependency_summary(global_ranked, rna_only, crispr_functional, resistance)
    overlap_table = build_topn_overlap_table(merged)

    out_dir = Path(config["candidate_adjudication"]["output"]["tables_dir"])
    out_dir.mkdir(parents=True, exist_ok=True)
    corr_table.to_csv(out_dir / "modality_dependency_correlations.tsv", sep="\t", index=False)
    overlap_table.to_csv(out_dir / "modality_dependency_topn_overlap.tsv", sep="\t", index=False)
    return {"correlations": corr_table, "overlap": overlap_table, "merged": merged}


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run_modality_dependency()
