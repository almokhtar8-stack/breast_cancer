"""Systems-network phase 20: hub / connector analysis.

Degree, betweenness, and closeness centrality on the focused network.
Centrality is descriptive of network position only -- a high-degree node is
often simply a well-studied, promiscuously-annotated protein (e.g. ubiquitin
core machinery), not evidence of therapeutic importance by itself (task
spec, explicit caution). candidate_connector / bridge_gene / multi_pathway_node
flags are separate, more specific structural roles computed alongside raw
centrality, not derived from centrality rank alone.
"""

from __future__ import annotations

import logging
from pathlib import Path

import networkx as nx
import pandas as pd
import yaml

logger = logging.getLogger(__name__)

CANDIDATES = ["USP34", "VEZF1", "EML5", "CITED2"]


def _load_config(config_path: str | Path) -> dict:
    with open(config_path) as f:
        return yaml.safe_load(f)


def build_graph(nodes: pd.DataFrame, edges: pd.DataFrame) -> nx.Graph:
    g = nx.Graph()
    g.add_nodes_from(nodes["gene"])
    for _, row in edges.iterrows():
        g.add_edge(row["source_gene"], row["target_gene"])
    return g


def build_node_metrics(g: nx.Graph, nodes: pd.DataFrame) -> pd.DataFrame:
    degree = dict(g.degree())
    betweenness = nx.betweenness_centrality(g)
    largest_cc = max(nx.connected_components(g), key=len)
    closeness = nx.closeness_centrality(g.subgraph(largest_cc))

    def candidate_neighbor_count(gene: str) -> int:
        neigh = set(g.neighbors(gene)) if gene in g else set()
        return len(neigh & set(CANDIDATES))

    rows = []
    for gene in nodes["gene"]:
        deg = degree.get(gene, 0)
        btw = betweenness.get(gene, 0.0)
        clo = closeness.get(gene, float("nan"))
        n_cand_neigh = candidate_neighbor_count(gene)
        rows.append(
            {
                "gene": gene,
                "degree": deg,
                "betweenness": btw,
                "closeness": clo,
                "is_candidate_connector": n_cand_neigh > 0 and gene not in CANDIDATES,
                "n_candidates_directly_connected": n_cand_neigh,
                "is_bridge_gene_high_betweenness": btw > 0 and btw >= pd.Series(list(betweenness.values())).quantile(0.9),
            }
        )
    out = pd.DataFrame(rows).sort_values("betweenness", ascending=False)
    logger.info("build_node_metrics: %d nodes; top-degree=%s", len(out), out.sort_values("degree", ascending=False)["gene"].head(5).tolist())
    return out


def run_hubs(config_path: str | Path = "config/config.yaml") -> pd.DataFrame:
    config = _load_config(config_path)
    cfg = config["systems_network"]
    networks_dir = Path(cfg["output"]["networks_dir"])

    nodes = pd.read_csv(networks_dir / "nodes.tsv", sep="\t")
    edges = pd.read_csv(networks_dir / "edges.tsv", sep="\t")
    g = build_graph(nodes, edges)

    out = build_node_metrics(g, nodes)
    out.to_csv(networks_dir / "network_node_metrics.tsv", sep="\t", index=False)
    return out


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run_hubs()
