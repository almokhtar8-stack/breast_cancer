"""Systems-network phase 19: network community/module detection.

Louvain community detection (networkx.community.louvain_communities) on the
focused network (undirected, edge-type-agnostic -- physical_PPI, regulatory,
and pathway_co_membership edges are all treated as connectivity for this
structural step only). An algorithmic community is NOT interpreted as a
biological pathway by itself -- each community is annotated with which
candidate genes, CRISPR-sensitising genes, and resistance-RNA genes (Phase
15 node attributes) it contains and left at that; no new "this community IS
pathway X" claim is manufactured.
"""

from __future__ import annotations

import logging
from pathlib import Path

import networkx as nx
import pandas as pd
import yaml

logger = logging.getLogger(__name__)

CANDIDATES = ["USP34", "VEZF1", "EML5", "CITED2"]
LOUVAIN_SEED = 20260813


def _load_config(config_path: str | Path) -> dict:
    with open(config_path) as f:
        return yaml.safe_load(f)


def build_graph(nodes: pd.DataFrame, edges: pd.DataFrame) -> nx.Graph:
    g = nx.Graph()
    g.add_nodes_from(nodes["gene"])
    for _, row in edges.iterrows():
        g.add_edge(row["source_gene"], row["target_gene"])
    return g


def build_network_modules(g: nx.Graph, nodes: pd.DataFrame) -> pd.DataFrame:
    node_class = nodes.set_index("gene")["node_class"]
    crispr_dir = nodes.set_index("gene")["crispr_direction"]

    non_isolated = [n for n in g.nodes if g.degree(n) > 0]
    subgraph = g.subgraph(non_isolated)
    communities = nx.community.louvain_communities(subgraph, seed=LOUVAIN_SEED, weight=None)

    rows = []
    for i, comm in enumerate(sorted(communities, key=len, reverse=True)):
        genes = sorted(comm)
        candidates_present = [g_ for g_ in genes if g_ in CANDIDATES]
        crispr_sensitising = [g_ for g_ in genes if crispr_dir.get(g_) == "sensitising_KO"]
        resistance_genes = [g_ for g_ in genes if node_class.get(g_) in ("resistance_gene", "leading_edge")]
        rows.append(
            {
                "module_id": i,
                "n_genes": len(genes),
                "genes": ";".join(genes),
                "candidate_genes_present": ",".join(candidates_present),
                "crispr_sensitising_genes": ",".join(crispr_sensitising),
                "resistance_rna_genes_count": len(resistance_genes),
                "note": "algorithmic community (Louvain, structural only) -- not itself a claimed biological pathway",
            }
        )
    isolated = [n for n in g.nodes if g.degree(n) == 0]
    if isolated:
        rows.append(
            {
                "module_id": -1,
                "n_genes": len(isolated),
                "genes": ";".join(sorted(isolated)),
                "candidate_genes_present": ",".join(sorted(set(isolated) & set(CANDIDATES))),
                "crispr_sensitising_genes": "",
                "resistance_rna_genes_count": 0,
                "note": "isolated nodes (zero edges in the focused network), not part of any community",
            }
        )
    out = pd.DataFrame(rows)
    logger.info("build_network_modules: %d communities (+isolated group), sizes=%s", len(communities), [len(c) for c in communities])
    return out


def run_modules(config_path: str | Path = "config/config.yaml") -> pd.DataFrame:
    config = _load_config(config_path)
    cfg = config["systems_network"]
    networks_dir = Path(cfg["output"]["networks_dir"])

    nodes = pd.read_csv(networks_dir / "nodes.tsv", sep="\t")
    edges = pd.read_csv(networks_dir / "edges.tsv", sep="\t")
    g = build_graph(nodes, edges)

    out = build_network_modules(g, nodes)
    out.to_csv(networks_dir / "network_modules.tsv", sep="\t", index=False)
    return out


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run_modules()
