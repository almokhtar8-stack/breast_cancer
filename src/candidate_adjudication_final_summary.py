"""Candidate adjudication Phase 28: a compact, manually-scannable summary
table -- one row per adjudicated candidate, ranked within the Phase 21
List A (multimodal therapeutic) first, then by global rank for everything
else. This is the table most likely to be reused directly in discussion,
so every column is a plain, already-computed value (no derived score).

Data source: `results/tables/candidate_adjudication/final_candidate_decision_table.tsv`,
`shortlist_A_multimodal_therapeutic.tsv`, `all_genes_cross_dataset_evidence_with_ranking.tsv`.
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


def build_final_summary(decision_table: pd.DataFrame, wide: pd.DataFrame, list_a_genes: list[str]) -> pd.DataFrame:
    merged = decision_table.merge(wide[["gene", "gse245601_epi_log2fc", "gse245601_epi_fdr", "gse245601_malignant_log2fc", "gse245601_malignant_fdr"]], on="gene", how="left")

    merged["RankWithinTherapeuticShortlist"] = merged["gene"].apply(lambda g: list_a_genes.index(g) + 1 if g in list_a_genes else None)
    out = pd.DataFrame(
        {
            "RankWithinTherapeuticShortlist": merged["RankWithinTherapeuticShortlist"],
            "Gene": merged["gene"],
            "CRISPR_Effect": merged["crispr_effect"].round(3),
            "CRISPR_FDR": merged["crispr_fdr"].round(4),
            "CRISPR_Direction": merged["crispr_sensitising"].map({True: "sensitising_KO", False: "tolerance_associated_KO"}),
            "GSE118713_log2FC": merged["gse118713_log2fc"].round(3),
            "GSE118713_FDR": merged["gse118713_fdr"].round(4),
            "GSE240112_log2FC": merged["gse240112_log2fc"].round(3),
            "GSE240112_FDR": merged["gse240112_fdr"].round(4),
            "GSE111151_log2FC": merged["gse111151_log2fc"].round(3),
            "GSE111151_FDR": merged["gse111151_fdr"].round(4),
            "GSE245601_epi_log2FC": merged["gse245601_epi_log2fc"].round(3),
            "GSE245601_epi_FDR": merged["gse245601_epi_fdr"].round(4),
            "GSE245601_malignant_log2FC": merged["gse245601_malignant_log2fc"].round(3),
            "GSE245601_malignant_FDR": merged["gse245601_malignant_fdr"].round(4),
            "ResistanceDirection": merged["resistance_consensus"],
            "ResistanceSignificantDatasets": merged["gene"].map(lambda g: int((decision_table.loc[decision_table["gene"] == g, ["gse118713_fdr", "gse240112_fdr", "gse111151_fdr"]] < 0.05).sum(axis=1).iloc[0])),
            "HumanSupport": merged["human_tumor_summary"],
            "Stability": merged["leave_one_out_stability"],
            "CandidateArchetype": merged["candidate_archetype"],
            "MainStrength": merged["main_strength"],
            "MainLimitation": merged["main_weakness"],
        }
    )
    out = out.sort_values(by=["RankWithinTherapeuticShortlist", "Gene"], na_position="last").reset_index(drop=True)
    logger.info("build_final_summary: %d rows", len(out))
    return out


def run_final_summary(config_path: str | Path = "config/config.yaml") -> pd.DataFrame:
    config = _load_config(config_path)
    cdx_out = config["cross_dataset_genomewide"]["output"]
    tables_dir = Path(cdx_out["wide_matrix_tsv"]).parent
    adj_tables = Path(config["candidate_adjudication"]["output"]["tables_dir"])

    decision_table = pd.read_csv(adj_tables / "final_candidate_decision_table.tsv", sep="\t")
    wide = pd.read_csv(tables_dir / "all_genes_cross_dataset_evidence_with_ranking.tsv", sep="\t")
    list_a_genes = pd.read_csv(adj_tables / "shortlist_A_multimodal_therapeutic.tsv", sep="\t")["gene"].tolist()

    summary = build_final_summary(decision_table, wide, list_a_genes)
    summary.to_csv(adj_tables / "final_candidate_summary.tsv", sep="\t", index=False)
    return summary


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run_final_summary()
