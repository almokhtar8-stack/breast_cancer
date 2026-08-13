"""Systems-network phase 17: candidate-centered subnetworks.

For each frozen candidate: the candidate itself, its direct (1-hop) curated
neighbors in the focused network (edges.tsv), and second-degree connectors
only when justified -- a 2-hop node is included only if it connects to >=2
of the candidate's distinct 1-hop neighbors (a real bridging pattern, not
indiscriminate radius expansion). This threshold is fixed and applied
identically to all four candidates, including EML5 (which, per Phase 14/15,
has zero curated neighbors at all -- its subnetwork is correctly just the
single isolated node, not an error).
"""

from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd
import yaml

logger = logging.getLogger(__name__)

CANDIDATES = ["USP34", "VEZF1", "EML5", "CITED2"]
MIN_SHARED_NEIGHBORS_FOR_2HOP = 2


def _load_config(config_path: str | Path) -> dict:
    with open(config_path) as f:
        return yaml.safe_load(f)


def _neighbors(edges: pd.DataFrame, gene: str) -> set[str]:
    a = set(edges.loc[edges["source_gene"] == gene, "target_gene"])
    b = set(edges.loc[edges["target_gene"] == gene, "source_gene"])
    return a | b


def build_candidate_subnetwork(candidate: str, edges: pd.DataFrame, nodes: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    one_hop = _neighbors(edges, candidate)

    two_hop_counts: dict[str, int] = {}
    for n in one_hop:
        for m in _neighbors(edges, n):
            if m == candidate or m in one_hop:
                continue
            two_hop_counts[m] = two_hop_counts.get(m, 0) + 1
    two_hop = {g for g, c in two_hop_counts.items() if c >= MIN_SHARED_NEIGHBORS_FOR_2HOP}

    subnetwork_genes = {candidate} | one_hop | two_hop
    sub_nodes = nodes.loc[nodes["gene"].isin(subnetwork_genes)].copy()
    sub_nodes["subnetwork_role"] = sub_nodes["gene"].apply(
        lambda g: "candidate" if g == candidate else ("direct_partner" if g in one_hop else "second_degree_connector")
    )
    sub_edges = edges.loc[edges["source_gene"].isin(subnetwork_genes) & edges["target_gene"].isin(subnetwork_genes)].copy()

    logger.info("build_candidate_subnetwork(%s): %d nodes (%d direct, %d second-degree), %d edges", candidate, len(sub_nodes), len(one_hop), len(two_hop), len(sub_edges))
    return sub_nodes, sub_edges


def run_subnetworks(config_path: str | Path = "config/config.yaml") -> None:
    config = _load_config(config_path)
    cfg = config["systems_network"]
    networks_dir = Path(cfg["output"]["networks_dir"])

    nodes = pd.read_csv(networks_dir / "nodes.tsv", sep="\t")
    edges = pd.read_csv(networks_dir / "edges.tsv", sep="\t")

    for candidate in CANDIDATES:
        sub_nodes, sub_edges = build_candidate_subnetwork(candidate, edges, nodes)
        prefix = candidate.lower()
        sub_nodes.to_csv(networks_dir / f"{prefix}_nodes.tsv", sep="\t", index=False)
        sub_edges.to_csv(networks_dir / f"{prefix}_edges.tsv", sep="\t", index=False)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run_subnetworks()
