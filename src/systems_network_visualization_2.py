"""Systems-network phase 8: candidate x pathway matrix figure.

Restricted to STRONG_CONSENSUS resistance pathways from Hallmark + Reactome
only (GO:BP excluded from this figure for readability -- GO:BP contributes
161 of 161 unique STRONG_CONSENSUS candidate-connected pathways almost
entirely redundant near-synonym terms; the full GO:BP evidence remains in
results/tables/systems_network/candidate_pathway_membership.tsv, Phase 7).
"""

from __future__ import annotations

import logging
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import yaml

logger = logging.getLogger(__name__)

CANDIDATES = ["USP34", "VEZF1", "EML5", "CITED2"]


def _load_config(config_path: str | Path) -> dict:
    with open(config_path) as f:
        return yaml.safe_load(f)


def select_matrix_pathways(membership: pd.DataFrame) -> pd.DataFrame:
    sub = membership.loc[
        (membership["resistance_consensus_class"] == "STRONG_CONSENSUS") & (membership["collection"].isin(["hallmark", "reactome"]))
    ].copy()
    return sub


def plot_candidate_pathway_matrix(sub: pd.DataFrame, consensus: pd.DataFrame, out_path: Path) -> None:
    pathway_keys = sorted(sub[["collection", "pathway"]].drop_duplicates().itertuples(index=False, name=None))
    consensus_idx = consensus.set_index(["collection", "pathway"])["median_NES"]

    labels = [f"{c.upper()}: {p.replace('HALLMARK_', '').replace('REACTOME_', '').replace('_', ' ').title()}" for c, p in pathway_keys]
    n_rows = len(pathway_keys)
    n_cols = len(CANDIDATES)

    fig, ax = plt.subplots(figsize=(7.5, max(5, n_rows * 0.35)))

    # background shading by resistance direction
    for i, key in enumerate(pathway_keys):
        nes = consensus_idx.get(key, 0)
        color = "#fbe3e0" if nes > 0 else "#dde8f5"
        ax.axhspan(i - 0.5, i + 0.5, color=color, zorder=0)

    lookup = sub.set_index(["candidate", "collection", "pathway"])
    for j, candidate in enumerate(CANDIDATES):
        for i, (collection, pathway) in enumerate(pathway_keys):
            key = (candidate, collection, pathway)
            if key not in lookup.index:
                continue
            row = lookup.loc[key]
            if isinstance(row, pd.DataFrame):
                row = row.iloc[0]
            is_member = bool(row["candidate_is_member"])
            is_le = isinstance(row["candidate_is_leading_edge_datasets"], str) and len(row["candidate_is_leading_edge_datasets"]) > 0
            has_interactor = isinstance(row["interactor_genes_in_pathway"], str) and len(row["interactor_genes_in_pathway"]) > 0

            if is_member:
                ax.scatter(j, i, marker="s", s=140, color="#1a1a1a", zorder=3)
            if is_le:
                ax.scatter(j, i, marker="^", s=100, facecolors="none", edgecolors="#1a1a1a", linewidths=1.6, zorder=4)
            if has_interactor and not is_member:
                ax.scatter(j, i, marker="o", s=60, facecolors="none", edgecolors="#555555", linewidths=1.3, zorder=2)

    ax.set_xlim(-0.5, n_cols - 0.5)
    ax.set_ylim(n_rows - 0.5, -0.5)
    ax.set_xticks(range(n_cols))
    ax.set_xticklabels(CANDIDATES, fontsize=10, fontweight="bold")
    ax.set_yticks(range(n_rows))
    ax.set_yticklabels(labels, fontsize=7.5)
    ax.set_title("Candidate x resistance-consensus pathway matrix\n(Hallmark + Reactome, STRONG_CONSENSUS pathways only)", fontsize=10)

    legend_elems = [
        plt.Line2D([0], [0], marker="s", color="w", markerfacecolor="#1a1a1a", markersize=10, label="direct curated pathway membership"),
        plt.Line2D([0], [0], marker="^", color="w", markerfacecolor="none", markeredgecolor="#1a1a1a", markersize=10, label="in GSEA leading edge (>=1 dataset)"),
        plt.Line2D([0], [0], marker="o", color="w", markerfacecolor="none", markeredgecolor="#555555", markersize=9, label="pathway contains a direct STRING interactor (score>=0.7)"),
        plt.Rectangle((0, 0), 1, 1, color="#fbe3e0", label="pathway median NES > 0 (up in resistance)"),
        plt.Rectangle((0, 0), 1, 1, color="#dde8f5", label="pathway median NES < 0 (down in resistance)"),
    ]
    ax.legend(handles=legend_elems, loc="upper center", bbox_to_anchor=(0.5, -0.10), fontsize=7.5, ncol=1, frameon=False)
    fig.subplots_adjust(bottom=0.28)
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    logger.info("plot_candidate_pathway_matrix: saved %s (%d pathways)", out_path, n_rows)


def run_visualization_2(config_path: str | Path = "config/config.yaml") -> None:
    config = _load_config(config_path)
    cfg = config["systems_network"]
    tables_dir = Path(cfg["output"]["tables_dir"])
    final_review_dir = Path(cfg["output"]["final_review_dir"])

    membership = pd.read_csv(tables_dir / "candidate_pathway_membership.tsv", sep="\t")
    consensus = pd.read_csv(tables_dir / "resistance_pathway_consensus.tsv", sep="\t")

    sub = select_matrix_pathways(membership)
    plot_candidate_pathway_matrix(sub, consensus, final_review_dir / "02_candidate_pathway_matrix.png")


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run_visualization_2()
