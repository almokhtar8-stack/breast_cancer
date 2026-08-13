"""Ad-hoc query on the already-built systems-network graph: shortest undirected
paths from USP34 to CTNNB1, PTEN, EP300, SOX2.

Reads ONLY the existing Cytoscape export tables -- does not rerun any upstream
systems-network phase and does not alter results/networks/systems_network/.
The graph is built exactly as it was exported: every edge keeps its own
interaction_type/database_source/confidence, and path length is unweighted
(edge count), consistent with how the rest of the systems-network phase treats
this graph (no confidence-weighted shortest-path claim is made anywhere else
in that phase either).

A 2+ edge path is never described as a "direct interaction" -- only an edge
that appears literally in network_edges.tsv between two genes is direct.
"""

from __future__ import annotations

import logging
from pathlib import Path

import networkx as nx
import pandas as pd

logger = logging.getLogger(__name__)

TARGET = "USP34"
DESTINATIONS = ["CTNNB1", "PTEN", "EP300", "SOX2"]
MAX_PATHS_REPORTED = 10

EDGES_PATH = Path("results/networks/systems_network/cytoscape/network_edges.tsv")
OUT_TABLE = Path("results/tables/systems_network/USP34_shortest_paths.tsv")


def load_graph(edges_path: Path = EDGES_PATH) -> tuple[nx.Graph, pd.DataFrame]:
    edges = pd.read_csv(edges_path, sep="\t")
    g = nx.Graph()
    for _, row in edges.iterrows():
        # Keep all parallel evidence for a pair -- a gene pair can have more
        # than one edge row (e.g. STRING + pathway_co_membership). Store them
        # as a list on the Graph edge so no evidence is silently dropped when
        # networkx collapses parallel edges into one undirected edge.
        u, v = row["source_gene"], row["target_gene"]
        rec = {
            "interaction_type": row["interaction_type"],
            "database_source": row["database_source"],
            "confidence": row["confidence"],
            "pathway": row["pathway"],
            "evidence_notes": row["evidence_notes"],
        }
        if g.has_edge(u, v):
            g[u][v]["records"].append(rec)
        else:
            g.add_edge(u, v, records=[rec])
    return g, edges


def _edge_records(g: nx.Graph, u: str, v: str) -> list[dict]:
    return g[u][v]["records"]


def shortest_paths_report(g: nx.Graph) -> pd.DataFrame:
    rows = []
    for dest in DESTINATIONS:
        if dest not in g or TARGET not in g:
            logger.warning("gene missing from graph: %s", dest)
            continue
        if not nx.has_path(g, TARGET, dest):
            rows.append(
                {
                    "target": dest,
                    "path_length_edges": None,
                    "path_rank": None,
                    "n_equally_short_paths": 0,
                    "path": "NO_PATH_IN_NETWORK",
                    "edge_index_in_path": None,
                    "source_gene": None,
                    "target_gene": None,
                    "interaction_type": None,
                    "database_source": None,
                    "confidence": None,
                    "pathway": None,
                    "evidence_notes": None,
                }
            )
            continue

        all_paths = list(nx.all_shortest_paths(g, TARGET, dest))
        path_len = len(all_paths[0]) - 1
        n_paths = len(all_paths)
        reported = all_paths[:MAX_PATHS_REPORTED]

        for rank, path in enumerate(reported, start=1):
            for i in range(len(path) - 1):
                u, v = path[i], path[i + 1]
                for rec in _edge_records(g, u, v):
                    rows.append(
                        {
                            "target": dest,
                            "path_length_edges": path_len,
                            "path_rank": rank,
                            "n_equally_short_paths": n_paths,
                            "path": " -> ".join(path),
                            "edge_index_in_path": i + 1,
                            "source_gene": u,
                            "target_gene": v,
                            "interaction_type": rec["interaction_type"],
                            "database_source": rec["database_source"],
                            "confidence": rec["confidence"],
                            "pathway": rec["pathway"],
                            "evidence_notes": rec["evidence_notes"],
                        }
                    )
    return pd.DataFrame(rows)


def shared_intermediate_nodes(g: nx.Graph) -> dict[str, list[str]]:
    """For each destination, the set of intermediate genes (excluding USP34
    and the destination itself) appearing in ANY of its shortest paths; then
    report which intermediates are shared across >=2 destinations' path sets."""
    per_dest: dict[str, set[str]] = {}
    for dest in DESTINATIONS:
        if dest not in g or TARGET not in g or not nx.has_path(g, TARGET, dest):
            per_dest[dest] = set()
            continue
        all_paths = list(nx.all_shortest_paths(g, TARGET, dest))
        intermediates: set[str] = set()
        for path in all_paths:
            intermediates.update(path[1:-1])
        per_dest[dest] = intermediates

    gene_to_dests: dict[str, list[str]] = {}
    for dest, genes in per_dest.items():
        for gene in genes:
            gene_to_dests.setdefault(gene, []).append(dest)

    shared = {gene: dests for gene, dests in gene_to_dests.items() if len(dests) >= 2}
    logger.info("shared intermediate nodes across >=2 USP34 target paths: %s", shared)
    return shared


def run(edges_path: Path = EDGES_PATH, out_table: Path = OUT_TABLE) -> pd.DataFrame:
    g, _ = load_graph(edges_path)
    report = shortest_paths_report(g)
    out_table.parent.mkdir(parents=True, exist_ok=True)
    report.to_csv(out_table, sep="\t", index=False)
    logger.info("wrote %s (%d rows)", out_table, len(report))
    shared_intermediate_nodes(g)
    return report


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run()
