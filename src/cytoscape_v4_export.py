"""Exports Cytoscape-ready node/edge/shortest-path tables from the EXACT
post-freeze exploratory STRING graph used by the v2/v3/v4 network figures.

Data source: `poster_network_mechanism_v2.build_network()`, imported and
called unmodified (species 9606, STRING interaction_partners,
required_score >= 0.7, local frozen TSVs -- no STRING requery, no network
call, no graph modification). This module only reshapes that graph into
three TSVs for Cytoscape import; it computes nothing new beyond
deterministic component ids and the already-defined shortest-path
enumeration. Every verification below is asserted at export time.

Outputs:
  results/tables/cytoscape_v4_network_edges.tsv
  results/tables/cytoscape_v4_network_nodes.tsv
  results/tables/cytoscape_v4_candidate_shortest_paths.tsv
"""

from __future__ import annotations

import logging
from pathlib import Path

import networkx as nx
import pandas as pd

from src.poster_network_mechanism_v2 import CANDIDATES, build_network

logger = logging.getLogger(__name__)

OUT_DIR = Path("results/tables")
EDGES_TSV = OUT_DIR / "cytoscape_v4_network_edges.tsv"
NODES_TSV = OUT_DIR / "cytoscape_v4_network_nodes.tsv"
PATHS_TSV = OUT_DIR / "cytoscape_v4_candidate_shortest_paths.tsv"

KIND_TO_LEVEL = {"candidate": 0, "level1_partner": 1, "level2_bridge": 2}


def _component_ids(graph: nx.Graph) -> dict[str, int]:
    """Deterministic component id: components ordered by (size desc,
    alphabetically-first member) and numbered 0, 1, 2, ..."""
    components = sorted(nx.connected_components(graph), key=lambda c: (-len(c), min(c)))
    return {node: idx for idx, comp in enumerate(components) for node in comp}


def build_edge_table(graph: nx.Graph) -> pd.DataFrame:
    """One row per unique undirected edge (endpoints alphabetized so A-B
    and B-A can never both appear), sorted deterministically."""
    rows = []
    for u, v, d in graph.edges(data=True):
        a, b = sorted((u, v))
        rows.append({
            "source": a,
            "target": b,
            "string_score": d["score"],
            "interaction_type": d["interaction_type"],
            "source_candidate": a in CANDIDATES,
            "target_candidate": b in CANDIDATES,
        })
    df = pd.DataFrame(rows).sort_values(["source", "target"]).reset_index(drop=True)
    assert len(df) == graph.number_of_edges()
    assert not df.duplicated(subset=["source", "target"]).any()
    logger.info("edge table: %d edges in graph -> %d rows out (0 lost)", graph.number_of_edges(), len(df))
    return df


def build_node_table(graph: nx.Graph) -> pd.DataFrame:
    component_of = _component_ids(graph)
    rows = []
    for node, data in graph.nodes(data=True):
        rows.append({
            "gene": node,
            "node_kind": data["kind"],
            "candidate": node in CANDIDATES,
            "level": KIND_TO_LEVEL[data["kind"]],
            "component_id": component_of[node],
            "degree": data["degree"],
            "betweenness": data["betweenness"],
        })
    df = pd.DataFrame(rows).sort_values("gene").reset_index(drop=True)
    assert len(df) == graph.number_of_nodes()
    logger.info("node table: %d nodes in graph -> %d rows out (0 lost)", graph.number_of_nodes(), len(df))
    return df


def build_shortest_path_table(graph: nx.Graph) -> pd.DataFrame:
    """ALL equally short paths for every candidate pair, computed from the
    exact graph. Disconnected pairs are recorded explicitly as NO_PATH
    rows rather than silently dropped."""
    rows = []
    for i, a in enumerate(CANDIDATES):
        for b in CANDIDATES[i + 1:]:
            try:
                all_paths = sorted(nx.all_shortest_paths(graph, a, b))
                for number, path in enumerate(all_paths, start=1):
                    rows.append({
                        "candidate_1": a,
                        "candidate_2": b,
                        "path_length": len(path) - 1,
                        "path_number": number,
                        "path": " -- ".join(path),
                    })
            except nx.NetworkXNoPath:
                rows.append({
                    "candidate_1": a,
                    "candidate_2": b,
                    "path_length": pd.NA,
                    "path_number": pd.NA,
                    "path": "NO_PATH",
                })
    df = pd.DataFrame(rows)
    n_pairs = len(CANDIDATES) * (len(CANDIDATES) - 1) // 2
    assert df.groupby(["candidate_1", "candidate_2"]).ngroups == n_pairs
    logger.info("shortest-path table: %d candidate pairs in -> %d rows out (disconnected pairs kept as NO_PATH)",
                n_pairs, len(df))
    return df


def verify_graph_invariants(graph: nx.Graph) -> None:
    assert graph.number_of_nodes() == 47
    assert graph.number_of_edges() == 147
    assert nx.number_connected_components(graph) == 3
    assert nx.has_path(graph, "KDM1A", "USP34")
    assert not nx.has_path(graph, "TLK2", "KDM1A")
    assert not nx.has_path(graph, "TLK2", "USP34")
    assert graph.degree("VEZF1") == 0
    kdm1a_usp34 = list(nx.all_shortest_paths(graph, "KDM1A", "USP34"))
    assert all(len(p) - 1 == 3 for p in kdm1a_usp34)
    assert len(kdm1a_usp34) == 4


def export_all() -> None:
    graph = build_network()
    verify_graph_invariants(graph)

    edges = build_edge_table(graph)
    nodes = build_node_table(graph)
    paths = build_shortest_path_table(graph)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    edges.to_csv(EDGES_TSV, sep="\t", index=False)
    nodes.to_csv(NODES_TSV, sep="\t", index=False)
    paths.to_csv(PATHS_TSV, sep="\t", index=False)
    logger.info("wrote %s (%d rows), %s (%d rows), %s (%d rows)",
                EDGES_TSV, len(edges), NODES_TSV, len(nodes), PATHS_TSV, len(paths))


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    export_all()
