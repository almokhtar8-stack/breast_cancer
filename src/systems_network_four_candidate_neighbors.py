"""Comparative systems-network audit, Part 2: direct (1-hop) network
neighborhood of each of the four frozen candidates (USP34, VEZF1, EML5,
CITED2).

Reads only the already-frozen Cytoscape export
(results/networks/systems_network/cytoscape/network_{nodes,edges}.tsv) --
does not rerun any upstream ranking/GSEA/network-build phase and does not
modify any frozen file. A "direct neighbor" here means literally present as
one row in network_edges.tsv with the candidate as source_gene or
target_gene -- the same STRING required_score>=0.7 / TRRUST / pathway
co-membership thresholds already used to build that file (no threshold is
lowered here to manufacture a neighborhood; see
docs/SYSTEMS_NETWORK_NODE_RULE.md). If a candidate has zero such rows, that
is reported as-is ("NO RESOLVED NETWORK NEIGHBOURHOOD IN CURRENT
ANALYSIS"), never worked around.
"""

from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd

logger = logging.getLogger(__name__)

CANDIDATES = ["USP34", "VEZF1", "EML5", "CITED2"]

NODES_PATH = Path("results/networks/systems_network/cytoscape/network_nodes.tsv")
EDGES_PATH = Path("results/networks/systems_network/cytoscape/network_edges.tsv")
OUT_TABLE = Path("results/tables/systems_network/four_candidate_direct_neighbors.tsv")


def build_direct_neighbors(nodes: pd.DataFrame, edges: pd.DataFrame) -> pd.DataFrame:
    nodes_idx = nodes.set_index("gene")
    rows = []
    for candidate in CANDIDATES:
        direct = edges.loc[(edges["source_gene"] == candidate) | (edges["target_gene"] == candidate)].copy()
        if direct.empty:
            rows.append(
                {
                    "candidate": candidate,
                    "neighbor_gene": "NO_RESOLVED_NETWORK_NEIGHBOURHOOD",
                    "interaction_type": "",
                    "database_source": "",
                    "confidence": pd.NA,
                    "neighbor_node_class": "",
                    "neighbor_crispr_direction": "",
                    "neighbor_crispr_fdr": pd.NA,
                    "neighbor_resistance_pattern": "",
                    "neighbor_number_consensus_pathways": pd.NA,
                    "neighbor_degree_in_network": pd.NA,
                }
            )
            continue
        direct["neighbor"] = direct.apply(lambda r: r["target_gene"] if r["source_gene"] == candidate else r["source_gene"], axis=1)
        for _, r in direct.iterrows():
            neighbor = r["neighbor"]
            n = nodes_idx.loc[neighbor] if neighbor in nodes_idx.index else None
            rows.append(
                {
                    "candidate": candidate,
                    "neighbor_gene": neighbor,
                    "interaction_type": r["interaction_type"],
                    "database_source": r["database_source"],
                    "confidence": r["confidence"],
                    "neighbor_node_class": n["node_class"] if n is not None else "not_in_node_table",
                    "neighbor_crispr_direction": n["crispr_direction"] if n is not None else "",
                    "neighbor_crispr_fdr": n["crispr_fdr"] if n is not None else pd.NA,
                    "neighbor_resistance_pattern": n["resistance_pattern"] if n is not None else "",
                    "neighbor_number_consensus_pathways": n["number_consensus_pathways"] if n is not None else pd.NA,
                    "neighbor_degree_in_network": pd.NA,
                }
            )
    out = pd.DataFrame(rows)

    # attach degree-in-network from network_node_metrics.tsv if available, so
    # a reader can see at a glance whether a neighbor is itself a high-degree
    # generic hub (Part 6) -- purely informational, computed nowhere else here
    metrics_path = Path("results/networks/systems_network/network_node_metrics.tsv")
    if metrics_path.exists():
        metrics = pd.read_csv(metrics_path, sep="\t").set_index("gene")["degree"]
        out["neighbor_degree_in_network"] = out["neighbor_gene"].map(metrics)

    logger.info(
        "build_direct_neighbors: %s",
        {c: int((out["candidate"] == c).sum() - (out.loc[out["candidate"] == c, "neighbor_gene"] == "NO_RESOLVED_NETWORK_NEIGHBOURHOOD").sum()) for c in CANDIDATES},
    )
    return out


def run(nodes_path: Path = NODES_PATH, edges_path: Path = EDGES_PATH, out_table: Path = OUT_TABLE) -> pd.DataFrame:
    nodes = pd.read_csv(nodes_path, sep="\t")
    edges = pd.read_csv(edges_path, sep="\t")
    out = build_direct_neighbors(nodes, edges)
    out_table.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(out_table, sep="\t", index=False)
    logger.info("wrote %s (%d rows)", out_table, len(out))
    return out


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run()
