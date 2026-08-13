"""Systems-network phase 22: USP34 dedicated systems-biology view.

USP34 is ranked #1 in the frozen therapeutic shortlist (docs/
THERAPEUTIC_SHORTLIST_FREEZE.md) -- this phase gives it a dedicated
network table/figure using exactly the same ranking rules and thresholds
already applied to every other gene in this phase (no bespoke relaxed
criteria for USP34).
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


def build_usp34_network_evidence(nodes: pd.DataFrame, edges: pd.DataFrame, candidate_membership: pd.DataFrame) -> pd.DataFrame:
    direct = edges.loc[(edges["source_gene"] == "USP34") | (edges["target_gene"] == "USP34")].copy()
    direct["partner"] = direct.apply(lambda r: r["target_gene"] if r["source_gene"] == "USP34" else r["source_gene"], axis=1)

    usp34_pathways = candidate_membership.loc[
        (candidate_membership["candidate"] == "USP34") & (candidate_membership["candidate_is_member"] | (candidate_membership["candidate_is_leading_edge_datasets"].fillna("") != ""))
    ]

    nodes_idx = nodes.set_index("gene")
    rows = []
    for _, row in direct.iterrows():
        partner = row["partner"]
        n = nodes_idx.loc[partner] if partner in nodes_idx.index else None
        rows.append(
            {
                "partner_gene": partner,
                "hop_distance": 1,
                "interaction_type": row["interaction_type"],
                "database_source": row["database_source"],
                "resistance_pattern": n["resistance_pattern"] if n is not None else "not_in_node_table",
                "crispr_direction": n["crispr_direction"] if n is not None else "not_in_node_table",
                "number_consensus_pathways": int(n["number_consensus_pathways"]) if n is not None else 0,
                "human_tumor_support": bool(n["human_tumor_support"]) if n is not None else False,
                "connects_usp34_to_pathway": ";".join(usp34_pathways["pathway"].tolist()) if partner in ("UBC", "UBA52", "UBB", "RPS27A") else "",
            }
        )
    out = pd.DataFrame(rows).sort_values(["hop_distance", "number_consensus_pathways"], ascending=[True, False])
    logger.info("build_usp34_network_evidence: %d direct partner rows", len(out))
    return out


def run_usp34_deep(config_path: str | Path = "config/config.yaml") -> pd.DataFrame:
    config = _load_config(config_path)
    cfg = config["systems_network"]
    tables_dir = Path(cfg["output"]["tables_dir"])
    networks_dir = Path(cfg["output"]["networks_dir"])

    nodes = pd.read_csv(networks_dir / "usp34_nodes.tsv", sep="\t")
    edges = pd.read_csv(networks_dir / "usp34_edges.tsv", sep="\t")
    candidate_membership = pd.read_csv(tables_dir / "candidate_pathway_membership.tsv", sep="\t")

    out = build_usp34_network_evidence(nodes, edges, candidate_membership)
    out.to_csv(tables_dir / "USP34_network_evidence.tsv", sep="\t", index=False)
    return out


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run_usp34_deep()
