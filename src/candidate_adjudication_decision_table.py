"""Candidate adjudication Phase 20: the final candidate decision table --
~20-30 genes (all seven MULTIMODAL_STRONG, the top-10 resistance-RNA
leaders, the top-10 CRISPR-sensitising leaders, and the near-miss
multimodal genes, deduplicated), built entirely from already-computed
adjudication tables. No external biology, no druggability, no numeric
composite score.

Data source: `results/tables/candidate_adjudication/*.tsv` (all built by
earlier phases of this adjudication) and the frozen cross-dataset
genome-wide tables.
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


def _main_strength_weakness(row: pd.Series) -> tuple[str, str]:
    strengths, weaknesses = [], []
    if pd.notna(row["crispr_fdr"]) and row["crispr_fdr"] < 0.05:
        strengths.append("CRISPR FDR<0.05")
    elif pd.notna(row["crispr_fdr"]) and row["crispr_fdr"] < 0.10:
        strengths.append("CRISPR FDR<0.10")
    else:
        weaknesses.append("no CRISPR significance")
    if row.get("crispr_direction") == "tolerance_associated_KO":
        weaknesses.append("CRISPR direction is tolerance-associated, not sensitising")
    res_n = row.get("resistance_fdr05_count", 0)
    res_n = 0 if pd.isna(res_n) else res_n
    if res_n >= 2:
        strengths.append(f"{int(res_n)}/3 resistance datasets FDR<0.05")
    elif res_n >= 1:
        strengths.append("1/3 resistance datasets FDR<0.05")
    else:
        weaknesses.append("no resistance-RNA dataset FDR<0.05")
    if row.get("resistance_direction_consensus") == "mixed":
        weaknesses.append("mixed resistance direction")
    human_hit = (pd.notna(row.get("gse245601_epi_fdr")) and row["gse245601_epi_fdr"] < 0.05) or (pd.notna(row.get("gse240112_tumor_fdr")) and row["gse240112_tumor_fdr"] < 0.05)
    if human_hit:
        strengths.append("human-tumor support")
    return "; ".join(strengths) if strengths else "none", "; ".join(weaknesses) if weaknesses else "none"


def build_decision_table(genes: list[str], wide: pd.DataFrame, resistance: pd.DataFrame, categories: pd.DataFrame, ranked: pd.DataFrame, archetypes: pd.DataFrame, stability: pd.DataFrame, crispr_functional: pd.DataFrame) -> pd.DataFrame:
    resistance = resistance.reset_index(drop=True).copy()
    resistance["resistance_rank"] = resistance.index + 1
    crispr_functional = crispr_functional.reset_index(drop=True).copy()
    crispr_functional["crispr_screen_rank"] = crispr_functional.index + 1

    merged = wide.loc[wide["gene"].isin(genes)].merge(resistance, on="gene", how="left").merge(categories, on="gene", how="left").merge(ranked[["gene", "global_rank"]], on="gene", how="left")
    merged = merged.merge(archetypes[["gene", "archetype"]], on="gene", how="left")
    merged = merged.merge(stability[["gene", "stability_label"]], on="gene", how="left")
    merged = merged.merge(crispr_functional[["gene", "crispr_screen_rank"]], on="gene", how="left")

    rows = []
    for _, row in merged.iterrows():
        strength, weakness = _main_strength_weakness(row)
        rows.append(
            {
                "gene": row["gene"], "candidate_archetype": row.get("archetype", "NA"),
                "crispr_effect": row["crispr_effect"], "crispr_fdr": row["crispr_fdr"], "crispr_sensitising": row["crispr_direction"] == "sensitising_KO",
                "gse118713_log2fc": row["gse118713_log2fc"], "gse118713_fdr": row["gse118713_fdr"],
                "gse240112_log2fc": row["gse240112_tumor_log2fc"], "gse240112_fdr": row["gse240112_tumor_fdr"],
                "gse111151_log2fc": row["gse111151_log2fc"], "gse111151_fdr": row["gse111151_fdr"],
                "resistance_consensus": row.get("resistance_direction_consensus", "NA"),
                "gse245601_acute_summary": f"epi_fdr={row['gse245601_epi_fdr']:.3g}, malig_fdr={row['gse245601_malignant_fdr']:.3g}" if pd.notna(row["gse245601_epi_fdr"]) else "not_testable",
                "human_tumor_summary": "significant" if (pd.notna(row.get("gse245601_epi_fdr")) and row["gse245601_epi_fdr"] < 0.05) or (pd.notna(row.get("gse240112_tumor_fdr")) and row["gse240112_tumor_fdr"] < 0.05) else "not_significant",
                "global_rank": row.get("global_rank"), "resistance_rank": row.get("resistance_rank"), "crispr_screen_rank": row.get("crispr_screen_rank"),
                "leave_one_out_stability": row.get("stability_label", "NA"),
                "coverage_tier": row["coverage_tier"], "evidence_category": row["evidence_category"],
                "main_strength": strength, "main_weakness": weakness,
            }
        )
    out = pd.DataFrame(rows).sort_values(by="global_rank", na_position="last").reset_index(drop=True)
    logger.info("build_decision_table: %d genes", len(out))
    return out


def run_decision_table(config_path: str | Path = "config/config.yaml") -> pd.DataFrame:
    config = _load_config(config_path)
    cdx_out = config["cross_dataset_genomewide"]["output"]
    tables_dir = Path(cdx_out["wide_matrix_tsv"]).parent
    adj_tables = Path(config["candidate_adjudication"]["output"]["tables_dir"])

    wide = pd.read_csv(tables_dir / "all_genes_cross_dataset_evidence_with_ranking.tsv", sep="\t")
    resistance = pd.read_csv(cdx_out["resistance_consensus_tsv"], sep="\t")
    categories = pd.read_csv(tables_dir / "evidence_categories.tsv", sep="\t")
    ranked = pd.read_csv(tables_dir / "global_ranking_eligible.tsv", sep="\t")
    archetypes = pd.read_csv(adj_tables / "candidate_archetypes.tsv", sep="\t")
    stability = pd.read_csv(cdx_out["ranking_stability_tsv"], sep="\t")
    crispr_functional = pd.read_csv(cdx_out["crispr_functional_all_genes_tsv"], sep="\t")

    multimodal7 = config["candidate_adjudication"]["multimodal7"]["genes"]
    top_resistance = pd.read_csv(adj_tables / "top_resistance_genes_exact.tsv", sep="\t")["gene"].head(10).tolist()
    top_crispr = pd.read_csv(adj_tables / "top_crispr_sensitising_exact.tsv", sep="\t")["gene"].head(10).tolist()
    near_misses = pd.read_csv(adj_tables / "multimodal_near_misses.tsv", sep="\t")["gene"].tolist()

    genes: list[str] = []
    for group in (multimodal7, top_resistance, top_crispr, near_misses):
        for g in group:
            if g not in genes:
                genes.append(g)

    table = build_decision_table(genes, wide, resistance, categories, ranked, archetypes, stability, crispr_functional)
    table.to_csv(adj_tables / "final_candidate_decision_table.tsv", sep="\t", index=False)
    return table


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run_decision_table()
