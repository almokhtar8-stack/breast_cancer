"""Systems-network phase 21: biological intervention-point classes.

Biology only -- no drug names, no structures, no toxicity/essentiality
normal-tissue data (that is explicitly the next phase, not this one).
Candidates for this table: the four frozen candidates themselves, plus
every node directly connected to >=1 candidate in the focused network
(node_metrics.is_candidate_connector), plus the top resistance_gene /
multimodal-pathway-driver nodes by number_consensus_pathways (Phase 15).
"""

from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd
import yaml

logger = logging.getLogger(__name__)

CANDIDATES = ["USP34", "VEZF1", "EML5", "CITED2"]
TOP_N_RESISTANCE_DRIVERS = 15


def _load_config(config_path: str | Path) -> dict:
    with open(config_path) as f:
        return yaml.safe_load(f)


def _network_role(gene: str, metrics: pd.Series) -> str:
    if gene in CANDIDATES:
        return "candidate"
    if metrics["is_bridge_gene_high_betweenness"]:
        return "bridge_gene_high_betweenness"
    if metrics["degree"] >= 20:
        return "high_degree_hub"
    if metrics["degree"] > 0:
        return "peripheral_connected_node"
    return "isolated"


def build_intervention_points(nodes: pd.DataFrame, metrics: pd.DataFrame, edges: pd.DataFrame) -> pd.DataFrame:
    metrics_idx = metrics.set_index("gene")
    nodes_idx = nodes.set_index("gene")

    connectors = set(metrics.loc[metrics["is_candidate_connector"], "gene"])
    top_drivers = set(
        nodes.loc[nodes["node_class"].isin(["resistance_gene", "leading_edge"])]
        .sort_values("number_consensus_pathways", ascending=False)
        .head(TOP_N_RESISTANCE_DRIVERS)["gene"]
    )
    genes = set(CANDIDATES) | connectors | top_drivers

    def candidates_connected(gene: str) -> list[str]:
        a = set(edges.loc[edges["source_gene"] == gene, "target_gene"])
        b = set(edges.loc[edges["target_gene"] == gene, "source_gene"])
        return sorted((a | b) & set(CANDIDATES))

    rows = []
    for gene in sorted(genes):
        m = metrics_idx.loc[gene]
        n = nodes_idx.loc[gene]
        connected_cands = candidates_connected(gene)

        if gene in CANDIDATES:
            relationship = "candidate_itself"
        elif len(connected_cands) >= 2:
            relationship = "shared_connector_between_candidates"
        elif connected_cands:
            relationship = f"direct_network_neighbor_of_{connected_cands[0]}"
        else:
            relationship = "resistance_pathway_driver_not_candidate_connected"

        reasons = []
        if gene in CANDIDATES:
            reasons.append("frozen therapeutic shortlist candidate")
        if connected_cands:
            reasons.append(f"direct curated network connection to {','.join(connected_cands)}")
        if n["node_class"] in ("crispr_sensitiser",):
            reasons.append("genome-wide CRISPR sensitising hit (FDR<0.05)")
        if n["number_consensus_pathways"] >= 5:
            reasons.append(f"leading-edge in {n['number_consensus_pathways']} STRONG_CONSENSUS resistance pathways")
        if m["is_bridge_gene_high_betweenness"]:
            reasons.append("high network betweenness (top decile) -- potential bridge between modules")

        rows.append(
            {
                "gene": gene,
                "relationship_to_candidate": relationship,
                "network_role": _network_role(gene, m),
                "CRISPR_support": n["crispr_direction"],
                "RNA_resistance_support": n["resistance_pattern"],
                "human_support": bool(n["human_tumor_support"]),
                "pathway": n["pathways_supported"] if pd.notna(n["pathways_supported"]) else "",
                "reason_for_interest": "; ".join(reasons) if reasons else "included as a candidate's direct network neighbor",
            }
        )
    out = pd.DataFrame(rows)
    logger.info("build_intervention_points: %d genes classified", len(out))
    return out


def run_intervention_points(config_path: str | Path = "config/config.yaml") -> pd.DataFrame:
    config = _load_config(config_path)
    cfg = config["systems_network"]
    tables_dir = Path(cfg["output"]["tables_dir"])
    networks_dir = Path(cfg["output"]["networks_dir"])

    nodes = pd.read_csv(networks_dir / "nodes.tsv", sep="\t")
    metrics = pd.read_csv(networks_dir / "network_node_metrics.tsv", sep="\t")
    edges = pd.read_csv(networks_dir / "edges.tsv", sep="\t")

    out = build_intervention_points(nodes, metrics, edges)
    out.to_csv(tables_dir / "potential_intervention_nodes_biology_only.tsv", sep="\t", index=False)
    return out


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run_intervention_points()
