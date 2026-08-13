"""Systems-network phase 18: candidate-candidate connections.

Checks, for every pair of the four frozen candidates, whether they are
connected by: a direct edge in the focused curated network (any type --
physical PPI, TF_target, or pathway co-membership), a shared 1-hop
neighbor in that network, or shared direct/leading-edge membership in the
same pathway per candidate_pathway_membership.tsv (Phase 7) -- the last
check is NOT restricted to the small, size-capped network edge set, so it
can surface pathway-level convergence the network's edge-count caps would
otherwise hide (see docs/SYSTEMS_NETWORK_NODE_RULE.md).
"""

from __future__ import annotations

import itertools
import logging
from pathlib import Path

import pandas as pd
import yaml

logger = logging.getLogger(__name__)

CANDIDATES = ["USP34", "VEZF1", "EML5", "CITED2"]


def _load_config(config_path: str | Path) -> dict:
    with open(config_path) as f:
        return yaml.safe_load(f)


def _neighbors(edges: pd.DataFrame, gene: str) -> set[str]:
    a = set(edges.loc[edges["source_gene"] == gene, "target_gene"])
    b = set(edges.loc[edges["target_gene"] == gene, "source_gene"])
    return a | b


def build_candidate_candidate_connections(edges: pd.DataFrame, candidate_membership: pd.DataFrame) -> pd.DataFrame:
    direct_evidence = candidate_membership.loc[
        candidate_membership["candidate_is_member"] | (candidate_membership["candidate_is_leading_edge_datasets"].fillna("") != "")
    ]
    pathway_to_candidates: dict[tuple[str, str], set[str]] = {}
    for _, row in direct_evidence.iterrows():
        pathway_to_candidates.setdefault((row["collection"], row["pathway"]), set()).add(row["candidate"])

    rows = []
    for a, b in itertools.combinations(CANDIDATES, 2):
        direct_edge = edges.loc[
            ((edges["source_gene"] == a) & (edges["target_gene"] == b)) | ((edges["source_gene"] == b) & (edges["target_gene"] == a))
        ]
        neigh_a, neigh_b = _neighbors(edges, a), _neighbors(edges, b)
        shared_neighbors = neigh_a & neigh_b

        shared_pathways = [(coll, pw, coll) for (coll, pw), cands in pathway_to_candidates.items() if {a, b} <= cands]

        if len(direct_edge):
            for _, e in direct_edge.iterrows():
                rows.append(
                    {
                        "candidate_A": a,
                        "candidate_B": b,
                        "connection_type": e["interaction_type"],
                        "intermediate_gene_if_any": "",
                        "pathway": e["pathway"],
                        "evidence_source": e["database_source"],
                        "direct_or_indirect": "direct",
                    }
                )
        if shared_neighbors:
            rows.append(
                {
                    "candidate_A": a,
                    "candidate_B": b,
                    "connection_type": "shared_1hop_network_neighbor",
                    "intermediate_gene_if_any": ";".join(sorted(shared_neighbors)),
                    "pathway": "",
                    "evidence_source": "focused_network (STRING/TRRUST/pathway_co_membership)",
                    "direct_or_indirect": "indirect",
                }
            )
        for coll, pw, _ in shared_pathways:
            rows.append(
                {
                    "candidate_A": a,
                    "candidate_B": b,
                    "connection_type": "shared_pathway_membership_or_leading_edge",
                    "intermediate_gene_if_any": "",
                    "pathway": f"{coll}:{pw}",
                    "evidence_source": "MSigDB / GSEA leading edge",
                    "direct_or_indirect": "indirect",
                }
            )
        if not len(direct_edge) and not shared_neighbors and not shared_pathways:
            rows.append(
                {
                    "candidate_A": a,
                    "candidate_B": b,
                    "connection_type": "none_found",
                    "intermediate_gene_if_any": "",
                    "pathway": "",
                    "evidence_source": "",
                    "direct_or_indirect": "none",
                }
            )

    out = pd.DataFrame(rows)
    logger.info("build_candidate_candidate_connections: %d rows across %d candidate pairs", len(out), len(list(itertools.combinations(CANDIDATES, 2))))
    return out


def run_candidate_connections(config_path: str | Path = "config/config.yaml") -> pd.DataFrame:
    config = _load_config(config_path)
    cfg = config["systems_network"]
    tables_dir = Path(cfg["output"]["tables_dir"])
    networks_dir = Path(cfg["output"]["networks_dir"])

    edges = pd.read_csv(networks_dir / "edges.tsv", sep="\t")
    candidate_membership = pd.read_csv(tables_dir / "candidate_pathway_membership.tsv", sep="\t")

    out = build_candidate_candidate_connections(edges, candidate_membership)
    out.to_csv(tables_dir / "candidate_candidate_connections.tsv", sep="\t", index=False)
    return out


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run_candidate_connections()
