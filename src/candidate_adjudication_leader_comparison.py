"""Candidate adjudication Phases 15-16: side-by-side comparison of the
seven MULTIMODAL_STRONG genes against (a) the top-10 resistance-RNA
leaders and (b) the top-10 CRISPR-sensitising leaders, using only the
already-computed archetype/axis/evidence tables -- no new numeric score.

Data source: `results/tables/candidate_adjudication/candidate_archetypes.tsv`,
`three_axis_candidate_matrix.tsv`, `adjudication_candidate_pool.tsv`,
`top_resistance_genes_exact.tsv`, `top_crispr_sensitising_exact.tsv`.
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


def _evidence_summary(pool_row: pd.Series) -> tuple[str, str, str]:
    crispr = f"fdr={pool_row['crispr_fdr']:.3g}, {pool_row['crispr_direction']}" if pd.notna(pool_row["crispr_fdr"]) else "not_testable"
    resistance = f"{int(pool_row['resistance_fdr05_count']) if pd.notna(pool_row['resistance_fdr05_count']) else 0}/3 FDR<0.05, {pool_row.get('resistance_direction_consensus', 'NA')}"
    human = "GSE245601_epi_sig" if pd.notna(pool_row.get("gse245601_epi_fdr")) and pool_row["gse245601_epi_fdr"] < 0.05 else ("GSE240112_sig" if pd.notna(pool_row.get("gse240112_tumor_fdr")) and pool_row["gse240112_tumor_fdr"] < 0.05 else "none_FDR<0.05")
    return crispr, resistance, human


def build_comparison_table(genes: list[str], pool: pd.DataFrame, archetypes: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for gene in genes:
        prow = pool.loc[pool["gene"] == gene]
        arow = archetypes.loc[archetypes["gene"] == gene]
        if len(prow) == 0:
            continue
        prow = prow.iloc[0]
        archetype = arow.iloc[0]["archetype"] if len(arow) else "NA"
        crispr, resistance, human = _evidence_summary(prow)
        rows.append(
            {
                "gene": gene, "candidate_archetype": archetype, "crispr_evidence": crispr, "resistance_rna_evidence": resistance,
                "human_evidence": human, "global_evidence": f"rank {int(prow['global_rank']) if pd.notna(prow['global_rank']) else 'NA'}/15,255, {prow.get('coverage_tier', 'NA')}",
                "dataset_stability": arow.iloc[0].get("axis_a_functional", "NA") + "/" + arow.iloc[0].get("axis_b_resistance", "NA") + "/" + arow.iloc[0].get("axis_c_human", "NA") if len(arow) else "NA",
            }
        )
    return pd.DataFrame(rows)


def run_leader_comparison(config_path: str | Path = "config/config.yaml") -> dict[str, pd.DataFrame]:
    config = _load_config(config_path)
    tables_dir = Path(config["candidate_adjudication"]["output"]["tables_dir"])
    multimodal7 = config["candidate_adjudication"]["multimodal7"]["genes"]

    pool = pd.read_csv(tables_dir / "adjudication_candidate_pool.tsv", sep="\t")
    archetypes = pd.read_csv(tables_dir / "candidate_archetypes.tsv", sep="\t")
    top_resistance = pd.read_csv(tables_dir / "top_resistance_genes_exact.tsv", sep="\t")["gene"].head(10).tolist()
    top_crispr = pd.read_csv(tables_dir / "top_crispr_sensitising_exact.tsv", sep="\t")["gene"].head(10).tolist()

    vs_rna = build_comparison_table(multimodal7 + [g for g in top_resistance if g not in multimodal7], pool, archetypes)
    vs_functional = build_comparison_table(multimodal7 + [g for g in top_crispr if g not in multimodal7], pool, archetypes)

    vs_rna.to_csv(tables_dir / "multimodal_vs_rna_leaders.tsv", sep="\t", index=False)
    vs_functional.to_csv(tables_dir / "multimodal_vs_functional_leaders.tsv", sep="\t", index=False)
    logger.info("run_leader_comparison: vs_rna=%d rows, vs_functional=%d rows", len(vs_rna), len(vs_functional))
    return {"vs_rna": vs_rna, "vs_functional": vs_functional}


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run_leader_comparison()
