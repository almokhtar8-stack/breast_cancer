"""Candidate adjudication Phase 17: leave-one-dataset-out rank stability
for the adjudication candidate sets (the seven MULTIMODAL_STRONG genes,
the top-10 resistance-RNA leaders, the top-10 CRISPR-sensitising
leaders), read directly from the already-frozen
`ranking_stability.tsv` -- no new weighting scheme, no re-derivation of
the leave-one-out hierarchy (that was fixed and re-verified during the
cross-dataset genome-wide integration's own Codex review).

Data source: `results/tables/cross_dataset_genomewide/ranking_stability.tsv`
(frozen, read-only).
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


def build_candidate_rank_stability(stability: pd.DataFrame, gene_groups: dict[str, list[str]]) -> pd.DataFrame:
    rows = []
    seen: set[str] = set()
    for group_name, genes in gene_groups.items():
        for gene in genes:
            if gene in seen:
                continue
            seen.add(gene)
            row = stability.loc[stability["gene"] == gene]
            if len(row) == 0:
                rows.append({"gene": gene, "candidate_group": group_name, "global_rank": None, "note": "not in primary ranking (LOW_COVERAGE)"})
                continue
            r = row.iloc[0]
            rows.append(
                {
                    "gene": gene, "candidate_group": group_name, "global_rank": r["rank_main"],
                    "rank_without_crispr": r["rank_without_crispr"], "rank_without_gse118713": r["rank_without_gse118713"],
                    "rank_without_gse245601": r["rank_without_gse245601"], "rank_without_gse240112": r["rank_without_gse240112"],
                    "rank_without_gse111151": r["rank_without_gse111151"], "median_alternate_rank": r["median_rank"],
                    "best_rank": r["best_rank"], "worst_rank": r["worst_rank"], "rank_range": r["worst_rank"] - r["best_rank"],
                    "n_top20_appearances_of_9_schemes": r["n_top20_appearances"], "stability_label": r["stability_label"],
                }
            )
    out = pd.DataFrame(rows).sort_values(by=["global_rank"], na_position="last").reset_index(drop=True)
    logger.info("build_candidate_rank_stability: %d genes across %d candidate groups", len(out), len(gene_groups))
    return out


def run_stability(config_path: str | Path = "config/config.yaml") -> pd.DataFrame:
    config = _load_config(config_path)
    cdx_out = config["cross_dataset_genomewide"]["output"]
    tables_dir = Path(cdx_out["wide_matrix_tsv"]).parent
    adj_tables = Path(config["candidate_adjudication"]["output"]["tables_dir"])

    stability = pd.read_csv(cdx_out["ranking_stability_tsv"], sep="\t")
    top_resistance = pd.read_csv(adj_tables / "top_resistance_genes_exact.tsv", sep="\t")["gene"].head(10).tolist()
    top_crispr = pd.read_csv(adj_tables / "top_crispr_sensitising_exact.tsv", sep="\t")["gene"].head(10).tolist()
    multimodal7 = config["candidate_adjudication"]["multimodal7"]["genes"]

    gene_groups = {"MULTIMODAL_STRONG": multimodal7, "top10_resistance_rna": top_resistance, "top10_crispr_sensitising": top_crispr}
    table = build_candidate_rank_stability(stability, gene_groups)
    table.to_csv(adj_tables / "candidate_rank_stability.tsv", sep="\t", index=False)
    return table


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run_stability()
