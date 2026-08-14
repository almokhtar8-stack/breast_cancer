"""Comparative systems-network audit, Part 4: shortest undirected paths from
each of the four frozen candidates to a small, candidate-specific set of
biologically relevant resistance-system nodes.

Targets are NOT the same across candidates (the task explicitly forbids
forcing USP34's WNT-effector target list onto the others). Each list is
derived from that candidate's own existing systems-network evidence, read
directly from already-frozen tables -- never chosen post hoc to make a path
"work":

  - USP34: reuses the existing, already-computed
    results/tables/systems_network/USP34_shortest_paths.tsv (CTNNB1, PTEN,
    EP300, SOX2) verbatim -- not recomputed.
  - CITED2: its five direct STRING interaction partners at
    required_score>=0.7 that are also independently biologically named in
    the task (EP300, CREBBP, TFAP2C, HIF1A, TP53) -- all five are already
    literally 1-hop neighbors in network_edges.tsv (Part 2 output), so this
    section reports and confirms that directness rather than discovering a
    longer path.
  - VEZF1: VEZF1's connected component in the frozen network is exactly
    {VEZF1, DMTN} (network_edges.tsv has exactly one row touching VEZF1).
    DMTN is therefore the only node any path can reach; no other vascular/
    developmental resistance gene is wired to VEZF1 in this frozen network,
    even though several are VEZF1 pathway co-members (Part 3) -- pathway
    co-membership only becomes a pairwise edge when the pathway's
    node-universe leading edge is <=10 genes (MAX_PATHWAY_GENES_FOR_EDGES,
    src/systems_network_build.py), and VEZF1's larger pathways
    (e.g. GOBP_BLOOD_VESSEL_MORPHOGENESIS) exceed that cap, so they do not
    produce edges here. That cap is an existing, documented project design
    decision and is not altered by this audit.
  - EML5: absent from the network entirely (zero edges) -- no path of any
    length exists to any node, and none is manufactured.

No path here exceeds 2 hops.
"""

from __future__ import annotations

import logging
from pathlib import Path

import networkx as nx
import pandas as pd

logger = logging.getLogger(__name__)

EDGES_PATH = Path("results/networks/systems_network/cytoscape/network_edges.tsv")
USP34_PATHS_PATH = Path("results/tables/systems_network/USP34_shortest_paths.tsv")
OUT_TABLE = Path("results/tables/systems_network/four_candidate_shortest_paths.tsv")

MAX_PATHS_REPORTED = 10

CITED2_TARGETS = ["EP300", "CREBBP", "TFAP2C", "HIF1A", "TP53"]
VEZF1_TARGETS = ["DMTN"]


def load_graph(edges_path: Path = EDGES_PATH) -> nx.Graph:
    edges = pd.read_csv(edges_path, sep="\t")
    g = nx.Graph()
    for _, row in edges.iterrows():
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
    return g


def _paths_for(g: nx.Graph, source: str, dest: str) -> list[list[str]]:
    if source not in g or dest not in g or not nx.has_path(g, source, dest):
        return []
    return list(nx.all_shortest_paths(g, source, dest))


def build_candidate_target_rows(g: nx.Graph, candidate: str, dest: str) -> list[dict]:
    rows = []
    all_paths = _paths_for(g, candidate, dest)
    if not all_paths:
        rows.append(
            {
                "candidate": candidate,
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
                "evidence_notes": None,
                "flag_weak_remote": False,
            }
        )
        return rows

    path_len = len(all_paths[0]) - 1
    n_paths = len(all_paths)
    weak_remote = path_len > 2
    for rank, path in enumerate(all_paths[:MAX_PATHS_REPORTED], start=1):
        for i in range(len(path) - 1):
            u, v = path[i], path[i + 1]
            for rec in g[u][v]["records"]:
                rows.append(
                    {
                        "candidate": candidate,
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
                        "evidence_notes": rec["evidence_notes"],
                        "flag_weak_remote": weak_remote,
                    }
                )
    return rows


def build_four_candidate_shortest_paths(g: nx.Graph, usp34_paths: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict] = []

    # USP34: reuse verbatim, tag with candidate column
    usp = usp34_paths.copy()
    usp.insert(0, "candidate", "USP34")
    usp["flag_weak_remote"] = usp["path_length_edges"].apply(lambda x: bool(x) and x > 2)
    rows.extend(usp.to_dict("records"))

    for target in CITED2_TARGETS:
        rows.extend(build_candidate_target_rows(g, "CITED2", target))

    for target in VEZF1_TARGETS:
        rows.extend(build_candidate_target_rows(g, "VEZF1", target))

    # EML5: not in the graph at all -- report explicitly, no path attempted
    rows.append(
        {
            "candidate": "EML5",
            "target": "NONE",
            "path_length_edges": None,
            "path_rank": None,
            "n_equally_short_paths": 0,
            "path": "NO_RESOLVED_NETWORK_NEIGHBOURHOOD_IN_CURRENT_ANALYSIS",
            "edge_index_in_path": None,
            "source_gene": None,
            "target_gene": None,
            "interaction_type": None,
            "database_source": None,
            "confidence": None,
            "evidence_notes": "EML5 has zero edges in the frozen network; no shortest-path analysis is possible.",
            "flag_weak_remote": False,
        }
    )

    out = pd.DataFrame(rows)
    cols = [
        "candidate",
        "target",
        "path_length_edges",
        "path_rank",
        "n_equally_short_paths",
        "path",
        "edge_index_in_path",
        "source_gene",
        "target_gene",
        "interaction_type",
        "database_source",
        "confidence",
        "evidence_notes",
        "flag_weak_remote",
    ]
    return out[cols]


def run(edges_path: Path = EDGES_PATH, usp34_paths_path: Path = USP34_PATHS_PATH, out_table: Path = OUT_TABLE) -> pd.DataFrame:
    g = load_graph(edges_path)
    usp34_paths = pd.read_csv(usp34_paths_path, sep="\t")
    out = build_four_candidate_shortest_paths(g, usp34_paths)
    out_table.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(out_table, sep="\t", index=False)
    logger.info("wrote %s (%d rows)", out_table, len(out))
    return out


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run()
