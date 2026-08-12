"""Cross-dataset genome-wide integration, Phase 24: surprise discovery.
Identifies, from the already-finalized genome-wide, candidate-list-independent ranking, which genes are
new relative to the project's prior CRISPR-first strategy -- computed
only AFTER the ranking itself was finalized (this module never feeds
back into any ranking decision).

"Previously known" genes = the 28 genes in the already-frozen
`results/tables/crispr_gse118713_master_table.tsv` (Gate-1 FDR<0.1 CRISPR
hits, the project's prior CRISPR-first candidate set). Genes outside that
set that rank strongly in the genome-wide global ranking are the
`is_surprise` genes below.

Data source: `results/tables/cross_dataset_genomewide/global_ranking_eligible.tsv`,
`results/tables/crispr_gse118713_master_table.tsv` (already frozen, read-only).
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


def build_surprise_discovery(ranked: pd.DataFrame, old28_genes: set[str], top_n: int = 20) -> pd.DataFrame:
    """One row per gene in the top ``top_n`` of the global ranking, with
    an explicit ``is_surprise`` flag (True = not among the project's
    prior 28 CRISPR-significant hits)."""
    top = ranked.head(top_n).copy()
    top["is_surprise"] = ~top["gene"].isin(old28_genes)
    top["previously_in_28_crispr_hits"] = top["gene"].isin(old28_genes)
    logger.info("build_surprise_discovery: %d/%d top-%d genes are surprises (outside the prior 28 CRISPR hits)", top["is_surprise"].sum(), len(top), top_n)
    return top


def build_fallen_genes(ranked: pd.DataFrame, old28_genes: set[str], full: pd.DataFrame, rank_threshold: int = 1000) -> pd.DataFrame:
    """Genes that WERE in the prior 28 CRISPR hits but rank far down (or
    are absent from) the genome-wide global ranking -- the reverse question
    to build_surprise_discovery."""
    rank_by_gene = ranked.set_index("gene")["global_rank"]
    rows = []
    for gene in sorted(old28_genes):
        if gene in rank_by_gene.index:
            rank = int(rank_by_gene.loc[gene])
            if rank > rank_threshold:
                rows.append({"gene": gene, "global_rank": rank, "status": "ranked_but_far_down"})
        else:
            in_full = gene in set(full["gene"])
            n_testable = int(full.loc[full["gene"] == gene, "n_datasets_testable"].iloc[0]) if in_full else 0
            rows.append({"gene": gene, "global_rank": pd.NA, "status": f"below_coverage_threshold_{n_testable}_of_5_testable" if in_full else "not_in_universe"})
    out = pd.DataFrame(rows)
    logger.info("build_fallen_genes: %d of the prior 28 CRISPR hits rank below top-%d or lack coverage", len(out), rank_threshold)
    return out


def run_surprise_discovery(config_path: str | Path = "config/config.yaml") -> dict[str, pd.DataFrame]:
    config = _load_config(config_path)
    cfg = config["cross_dataset_genomewide"]
    out = cfg["output"]
    tables_dir = Path(out["wide_matrix_tsv"]).parent

    ranked = pd.read_csv(tables_dir / "global_ranking_eligible.tsv", sep="\t")
    full = pd.read_csv(tables_dir / "all_genes_cross_dataset_evidence_with_ranking.tsv", sep="\t")
    old28 = pd.read_csv("results/tables/crispr_gse118713_master_table.tsv", sep="\t")
    old28_genes = set(old28["gene_symbol"])

    surprises = build_surprise_discovery(ranked, old28_genes)
    fallen = build_fallen_genes(ranked, old28_genes, full)

    surprises_path = tables_dir / "surprise_discovery_top20.tsv"
    fallen_path = tables_dir / "fallen_genes_from_prior_28.tsv"
    surprises.to_csv(surprises_path, sep="\t", index=False)
    fallen.to_csv(fallen_path, sep="\t", index=False)
    logger.info("wrote %s and %s", surprises_path, fallen_path)

    return {"surprises": surprises, "fallen": fallen}


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run_surprise_discovery()
