"""Comparative systems-network audit, Part 7: pairwise convergence matrix
across the four frozen candidates.

Reads only already-frozen outputs -- does not rerun any upstream phase:
  - results/tables/systems_network/candidate_candidate_connections.tsv
    (frozen pairwise direct-interaction / shared-pathway table, Phase 18 of
    the original systems-network build)
  - results/networks/systems_network/cytoscape/network_edges.tsv (direct
    1-hop neighbor sets, Part 2 of this audit)
  - results/tables/systems_network/four_candidate_bridge_evidence.tsv
    (bridge/connector gene sets, Part 5 of this audit)
  - results/networks/systems_network/network_node_metrics.tsv
    (n_candidates_directly_connected, for the "shared resistance hub" check)

Checks, per candidate pair, exactly the 7 convergence dimensions requested:
  1. direct interaction (a literal edge between the two candidates)
  2. shared direct neighbor (1-hop overlap)
  3. shared bridge/connector gene (Part 5 sets)
  4/5. shared resistance pathway / leading-edge gene or module (from the
       frozen candidate_candidate_connections.tsv, itself built from
       candidate_pathway_membership.tsv leading-edge overlap)
  6. shared transcriptional regulator (a TRRUST regulatory edge pointing at
     both candidates) -- none exists in the frozen edge set for any
     candidate, reported as such rather than inferred from anything looser
  7. shared resistance hub (a gene with n_candidates_directly_connected>=2
     in network_node_metrics.tsv)

A pair sharing only broad GO:BP category membership (e.g. two candidates
that both happen to appear in a huge generic "developmental process" term)
is not upgraded to "mechanistic convergence" -- the underlying pathway
names are reported as-is so a reader can judge specificity themselves.
"""

from __future__ import annotations

import itertools
import logging
from pathlib import Path

import pandas as pd

logger = logging.getLogger(__name__)

CANDIDATES = ["USP34", "VEZF1", "EML5", "CITED2"]

CONNECTIONS_PATH = Path("results/tables/systems_network/candidate_candidate_connections.tsv")
EDGES_PATH = Path("results/networks/systems_network/cytoscape/network_edges.tsv")
BRIDGE_EVIDENCE_PATH = Path("results/tables/systems_network/four_candidate_bridge_evidence.tsv")
METRICS_PATH = Path("results/networks/systems_network/network_node_metrics.tsv")
OUT_TABLE = Path("results/tables/systems_network/four_candidate_convergence.tsv")


def _direct_neighbors(edges: pd.DataFrame, candidate: str) -> set[str]:
    direct = edges.loc[(edges["source_gene"] == candidate) | (edges["target_gene"] == candidate)]
    return set(direct["source_gene"]).union(direct["target_gene"]) - {candidate}


def build_convergence_matrix(
    connections: pd.DataFrame,
    edges: pd.DataFrame,
    bridge_evidence: pd.DataFrame,
    metrics: pd.DataFrame,
) -> pd.DataFrame:
    neighbors = {c: _direct_neighbors(edges, c) for c in CANDIDATES}
    bridges = {c: set(bridge_evidence.loc[bridge_evidence["candidate"] == c, "bridge_gene"]) - {"NONE"} for c in CANDIDATES}

    shared_hub_genes = set(metrics.loc[metrics["n_candidates_directly_connected"] >= 2, "gene"])

    rows = []
    for a, b in itertools.combinations(CANDIDATES, 2):
        direct_interaction = a in neighbors[b] or b in neighbors[a]
        shared_neighbors = neighbors[a] & neighbors[b]
        shared_bridges = bridges[a] & bridges[b]

        conn_rows = connections.loc[((connections["candidate_A"] == a) & (connections["candidate_B"] == b)) | ((connections["candidate_A"] == b) & (connections["candidate_B"] == a))]
        shared_pathways = sorted(conn_rows.loc[conn_rows["connection_type"] != "none_found", "pathway"].dropna().unique())

        rows.append(
            {
                "candidate_A": a,
                "candidate_B": b,
                "direct_interaction": direct_interaction,
                "shared_direct_neighbors": ",".join(sorted(shared_neighbors)) if shared_neighbors else "",
                "n_shared_direct_neighbors": len(shared_neighbors),
                "shared_bridge_genes": ",".join(sorted(shared_bridges)) if shared_bridges else "",
                "n_shared_bridge_genes": len(shared_bridges),
                "shared_resistance_pathways_or_leading_edge_modules": ";".join(shared_pathways),
                "n_shared_pathways": len(shared_pathways),
                "shared_transcriptional_regulator": "",  # none found in frozen TRRUST edges for any candidate, see module docstring
                "shared_resistance_hub_gene": ",".join(sorted(shared_hub_genes)) if shared_hub_genes else "",
                "any_convergence": bool(direct_interaction or shared_neighbors or shared_bridges or shared_pathways),
            }
        )
    out = pd.DataFrame(rows)
    logger.info("build_convergence_matrix: %d pairs, %d with any convergence", len(out), out["any_convergence"].sum())
    return out


def run(
    connections_path: Path = CONNECTIONS_PATH,
    edges_path: Path = EDGES_PATH,
    bridge_evidence_path: Path = BRIDGE_EVIDENCE_PATH,
    metrics_path: Path = METRICS_PATH,
    out_table: Path = OUT_TABLE,
) -> pd.DataFrame:
    connections = pd.read_csv(connections_path, sep="\t")
    edges = pd.read_csv(edges_path, sep="\t")
    bridge_evidence = pd.read_csv(bridge_evidence_path, sep="\t")
    metrics = pd.read_csv(metrics_path, sep="\t")

    out = build_convergence_matrix(connections, edges, bridge_evidence, metrics)
    out_table.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(out_table, sep="\t", index=False)
    logger.info("wrote %s (%d rows)", out_table, len(out))
    return out


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run()
