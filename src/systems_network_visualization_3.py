"""Systems-network phases 22-23: USP34 dedicated network figure and the
four-candidate network figure. Both restricted to 1-hop direct partners
only, for readability (the full 2-hop data remains in
usp34_nodes.tsv/edges.tsv etc, Phase 17) -- avoids the dense-hairball
anti-pattern explicitly called out in the task spec.
"""

from __future__ import annotations

import logging
from pathlib import Path

import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd
import yaml

logger = logging.getLogger(__name__)

CANDIDATES = ["USP34", "VEZF1", "EML5", "CITED2"]
CANDIDATE_COLOR = "#8b1a1a"
SENSITISER_COLOR = "#c0392b"
TOLERANCE_COLOR = "#2874a6"
NEUTRAL_COLOR = "#7f8c8d"
EDGE_STYLE = {"physical_PPI": "solid", "regulatory": "dashed", "pathway_co_membership": "dotted", "functional_association": "solid"}
EDGE_WIDTH = {"physical_PPI": 1.6, "regulatory": 1.6, "pathway_co_membership": 1.0, "functional_association": 0.8}


def _load_config(config_path: str | Path) -> dict:
    with open(config_path) as f:
        return yaml.safe_load(f)


def _node_color(gene: str, nodes_idx: pd.DataFrame) -> str:
    if gene in CANDIDATES:
        return CANDIDATE_COLOR
    if gene not in nodes_idx.index:
        return NEUTRAL_COLOR
    direction = nodes_idx.loc[gene, "crispr_direction"]
    if direction == "sensitising_KO":
        return SENSITISER_COLOR
    if direction == "tolerance_associated_KO":
        return TOLERANCE_COLOR
    return NEUTRAL_COLOR


def plot_usp34_network(nodes: pd.DataFrame, edges: pd.DataFrame, out_path: Path) -> None:
    one_hop_edges = edges.loc[(edges["source_gene"] == "USP34") | (edges["target_gene"] == "USP34")]
    g = nx.Graph()
    g.add_node("USP34")
    for _, row in one_hop_edges.iterrows():
        g.add_edge(row["source_gene"], row["target_gene"], interaction_type=row["interaction_type"])

    nodes_idx = nodes.set_index("gene")
    pos = nx.spring_layout(g, seed=20260813, k=0.9)

    fig, ax = plt.subplots(figsize=(8, 7))
    for etype, style in EDGE_STYLE.items():
        elist = [(u, v) for u, v, d in g.edges(data=True) if d.get("interaction_type") == etype]
        nx.draw_networkx_edges(g, pos, edgelist=elist, style=style, width=EDGE_WIDTH[etype], edge_color="#999999", ax=ax)

    colors = [_node_color(n, nodes_idx) for n in g.nodes]
    sizes = [1400 if n == "USP34" else 700 for n in g.nodes]
    nx.draw_networkx_nodes(g, pos, node_color=colors, node_size=sizes, edgecolors="black", linewidths=0.8, ax=ax)
    nx.draw_networkx_labels(g, pos, font_size=9, font_weight="bold", ax=ax)

    legend_elems = [
        plt.Line2D([0], [0], marker="o", color="w", markerfacecolor=CANDIDATE_COLOR, markersize=12, label="USP34 (candidate)"),
        plt.Line2D([0], [0], marker="o", color="w", markerfacecolor=SENSITISER_COLOR, markersize=10, label="CRISPR sensitising_KO (FDR<0.05)"),
        plt.Line2D([0], [0], marker="o", color="w", markerfacecolor=TOLERANCE_COLOR, markersize=10, label="CRISPR tolerance_associated_KO (FDR<0.05)"),
        plt.Line2D([0], [0], marker="o", color="w", markerfacecolor=NEUTRAL_COLOR, markersize=10, label="no significant CRISPR effect"),
        plt.Line2D([0], [0], color="#999999", linestyle="solid", label="physical_PPI (STRING, escore/dscore>0)"),
        plt.Line2D([0], [0], color="#999999", linestyle="solid", linewidth=0.8, label="functional_association (STRING, other channels)"),
    ]
    ax.legend(handles=legend_elems, loc="upper center", bbox_to_anchor=(0.5, -0.02), fontsize=8, ncol=1, frameon=False)
    ax.set_title("USP34: direct (1-hop) curated network partners\n(all high-confidence, score>=0.7 STRING edges; ubiquitin/proteasome-adjacent, no direct WNT-effector partner at this threshold)", fontsize=10)
    ax.axis("off")
    fig.subplots_adjust(bottom=0.32)
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    logger.info("plot_usp34_network: saved %s (%d nodes)", out_path, g.number_of_nodes())


def plot_four_candidate_network(nodes: pd.DataFrame, edges: pd.DataFrame, out_path: Path, top_n_per_candidate: int = 6) -> None:
    nodes_idx = nodes.set_index("gene")
    g = nx.Graph()
    g.add_nodes_from(CANDIDATES)

    for candidate in CANDIDATES:
        one_hop = edges.loc[(edges["source_gene"] == candidate) | (edges["target_gene"] == candidate)].copy()
        if not len(one_hop):
            continue
        one_hop["partner"] = one_hop.apply(lambda r: r["target_gene"] if r["source_gene"] == candidate else r["source_gene"], axis=1)
        one_hop["n_pathways"] = one_hop["partner"].map(lambda g_: nodes_idx.loc[g_, "number_consensus_pathways"] if g_ in nodes_idx.index else 0)
        top = one_hop.sort_values("n_pathways", ascending=False).head(top_n_per_candidate)
        for _, row in top.iterrows():
            g.add_edge(candidate, row["partner"], interaction_type=row["interaction_type"])

    pos = nx.spring_layout(g, seed=20260813, k=1.1)
    fig, ax = plt.subplots(figsize=(10, 9))
    for etype, style in EDGE_STYLE.items():
        elist = [(u, v) for u, v, d in g.edges(data=True) if d.get("interaction_type") == etype]
        nx.draw_networkx_edges(g, pos, edgelist=elist, style=style, width=EDGE_WIDTH[etype], edge_color="#999999", ax=ax)

    colors = [_node_color(n, nodes_idx) for n in g.nodes]
    sizes = [1600 if n in CANDIDATES else 600 for n in g.nodes]
    nx.draw_networkx_nodes(g, pos, node_color=colors, node_size=sizes, edgecolors="black", linewidths=0.8, ax=ax)
    nx.draw_networkx_labels(g, pos, font_size=8, font_weight="bold", ax=ax)

    legend_elems = [
        plt.Line2D([0], [0], marker="o", color="w", markerfacecolor=CANDIDATE_COLOR, markersize=12, label="frozen candidate"),
        plt.Line2D([0], [0], marker="o", color="w", markerfacecolor=SENSITISER_COLOR, markersize=10, label="CRISPR sensitising_KO"),
        plt.Line2D([0], [0], marker="o", color="w", markerfacecolor=TOLERANCE_COLOR, markersize=10, label="CRISPR tolerance_associated_KO"),
        plt.Line2D([0], [0], marker="o", color="w", markerfacecolor=NEUTRAL_COLOR, markersize=10, label="no significant CRISPR effect"),
        plt.Line2D([0], [0], color="#999999", linestyle="solid", label="physical_PPI"),
        plt.Line2D([0], [0], color="#999999", linestyle="dashed", label="regulatory (TF_target, TRRUST)"),
        plt.Line2D([0], [0], color="#999999", linestyle="dotted", label="pathway_co_membership"),
    ]
    ax.legend(handles=legend_elems, loc="upper center", bbox_to_anchor=(0.5, -0.02), fontsize=8, ncol=2, frameon=False)
    ax.set_title(f"Four frozen candidates + top {top_n_per_candidate} strongest-supported direct network neighbors each\n(no direct or shared 1-hop edge exists between any two candidates -- see candidate_candidate_connections.tsv)", fontsize=10)
    ax.axis("off")
    fig.subplots_adjust(bottom=0.22)
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    logger.info("plot_four_candidate_network: saved %s (%d nodes, %d edges)", out_path, g.number_of_nodes(), g.number_of_edges())


def run_visualization_3(config_path: str | Path = "config/config.yaml") -> None:
    config = _load_config(config_path)
    cfg = config["systems_network"]
    networks_dir = Path(cfg["output"]["networks_dir"])
    final_review_dir = Path(cfg["output"]["final_review_dir"])

    nodes = pd.read_csv(networks_dir / "nodes.tsv", sep="\t")
    edges = pd.read_csv(networks_dir / "edges.tsv", sep="\t")

    plot_usp34_network(nodes, edges, final_review_dir / "03_USP34_network.png")
    plot_four_candidate_network(nodes, edges, final_review_dir / "04_four_candidate_network.png")


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run_visualization_3()
