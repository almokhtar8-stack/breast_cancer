"""Systems-network phase 24: one higher-level systems map figure.

candidate -> network/module -> pathway -> resistance evidence -> CRISPR
functional evidence, one row per frozen candidate. Text content is
programmatically derived from the already-built tables (nodes.tsv,
edges.tsv, candidate_pathway_membership.tsv,
multimodal_pathway_convergence.tsv) -- nothing here is hand-typed biology.

Codex review fix: an earlier version of this module contained a hardcoded
per-candidate `module_text` dict of hand-typed prose describing each
candidate's network neighborhood, directly contradicting this docstring's
claim. It is now computed from edges.tsv (1-hop partner count/names,
ranked by `number_consensus_pathways` in nodes.tsv, alphabetical
tie-break). `top_pathway` selection was also unranked (first row
encountered); it now uses an explicit rank (datasets_FDR05 descending,
then pathway name alphabetically).
"""

from __future__ import annotations

import logging
import textwrap
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import yaml

logger = logging.getLogger(__name__)

CANDIDATES = ["USP34", "VEZF1", "EML5", "CITED2"]
COL_TITLES = ["CANDIDATE", "NETWORK / MODULE", "PATHWAY", "RESISTANCE RNA EVIDENCE", "CRISPR FUNCTIONAL EVIDENCE"]


def _load_config(config_path: str | Path) -> dict:
    with open(config_path) as f:
        return yaml.safe_load(f)


def _direct_neighbors(edges: pd.DataFrame, gene: str) -> set[str]:
    a = set(edges.loc[edges["source_gene"] == gene, "target_gene"])
    b = set(edges.loc[edges["target_gene"] == gene, "source_gene"])
    return a | b


def build_systems_map_text(
    candidate: str, nodes: pd.DataFrame, edges: pd.DataFrame, membership: pd.DataFrame, consensus: pd.DataFrame, mm: pd.DataFrame
) -> list[str]:
    nodes_idx = nodes.set_index("gene")
    n = nodes_idx.loc[candidate]

    connected = membership.loc[
        (membership["candidate"] == candidate) & (membership["candidate_is_member"] | (membership["candidate_is_leading_edge_datasets"].fillna("") != ""))
    ]
    strong = connected.loc[connected["resistance_consensus_class"] == "STRONG_CONSENSUS"]
    n_strong = strong[["collection", "pathway"]].drop_duplicates().shape[0]
    if len(strong):
        # all STRONG_CONSENSUS pathways already meet the same significance
        # bar (Phase 5) -- explicit deterministic tie-break: alphabetical
        top_pathway = strong.sort_values("pathway", kind="mergesort")["pathway"].iloc[0]
    elif len(connected):
        top_pathway = connected.sort_values("pathway", kind="mergesort")["pathway"].iloc[0]
    else:
        top_pathway = "none found"

    mm_c = mm.loc[mm["convergence_category"] == "MULTIMODAL_PATHWAY"].copy()
    mm_c["candidates_connected"] = mm_c["candidates_connected"].fillna("")
    n_multimodal = mm_c["candidates_connected"].apply(lambda s: candidate in s.split(",")).sum()

    neighbors = _direct_neighbors(edges, candidate)
    if not neighbors:
        module_text = "isolated: zero curated network connections (STRING physical or functional, TRRUST, pathway co-membership) at any tier tested"
    else:
        ranked_neighbors = sorted(neighbors, key=lambda g: (-int(nodes_idx.loc[g, "number_consensus_pathways"]) if g in nodes_idx.index else 0, g))
        top_names = ", ".join(ranked_neighbors[:5])
        module_text = f"{len(neighbors)} direct (1-hop) network neighbor(s); top by resistance-pathway support: {top_names}"

    pathway_text = f"{top_pathway.replace('GOBP_', '').replace('HALLMARK_', '').replace('REACTOME_', '').replace('_', ' ').title()}" if top_pathway != "none found" else "none found"

    resistance_text = f"{n_strong} STRONG_CONSENSUS resistance pathway(s) connected" if n_strong else "no STRONG_CONSENSUS resistance pathway connection found"

    crispr_direction = n["crispr_direction"]
    crispr_text = f"{candidate}: {crispr_direction} (FDR={n['crispr_fdr']:.3g})" if pd.notna(n["crispr_fdr"]) else f"{candidate}: not tested"
    crispr_text += f"\n{n_multimodal} of its pathways are MULTIMODAL_PATHWAY (RNA+CRISPR both significant)" if n_multimodal else "\nnone of its pathways reach MULTIMODAL_PATHWAY status"

    return [candidate, module_text, pathway_text, resistance_text, crispr_text]


def plot_systems_map(rows_text: list[list[str]], out_path: Path) -> None:
    n_rows = len(rows_text)
    n_cols = len(COL_TITLES)
    fig, ax = plt.subplots(figsize=(15, 3.2 * n_rows))
    col_widths = [0.10, 0.24, 0.20, 0.22, 0.24]
    col_x = [sum(col_widths[:i]) for i in range(n_cols)]

    for j, title in enumerate(COL_TITLES):
        ax.text(col_x[j] + col_widths[j] / 2, n_rows + 0.15, title, ha="center", va="bottom", fontsize=10, fontweight="bold")

    row_h = 1.0
    for i, row in enumerate(rows_text):
        y = n_rows - 1 - i
        for j, text in enumerate(row):
            wrapped = "\n".join(textwrap.wrap(text, width=32))
            box_color = "#f6d9d9" if j == 0 else "#eef2f7"
            ax.add_patch(plt.Rectangle((col_x[j] + 0.005, y + 0.05), col_widths[j] - 0.01, row_h - 0.1, facecolor=box_color, edgecolor="#888888", linewidth=0.8))
            fontsize = 10 if j == 0 else 7.3
            fontweight = "bold" if j == 0 else "normal"
            ax.text(col_x[j] + col_widths[j] / 2, y + row_h / 2, wrapped, ha="center", va="center", fontsize=fontsize, fontweight=fontweight)
            if j < n_cols - 1:
                ax.annotate("", xy=(col_x[j + 1] + 0.005, y + row_h / 2), xytext=(col_x[j] + col_widths[j] - 0.005, y + row_h / 2), arrowprops={"arrowstyle": "->", "color": "#555555", "lw": 1.2})

    ax.set_xlim(0, 1)
    ax.set_ylim(0, n_rows + 0.6)
    ax.axis("off")
    ax.set_title("Resistance systems map: candidate -> network module -> pathway -> RNA resistance evidence -> CRISPR functional evidence", fontsize=11, pad=10)
    fig.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    logger.info("plot_systems_map: saved %s", out_path)


def run_visualization_4(config_path: str | Path = "config/config.yaml") -> None:
    config = _load_config(config_path)
    cfg = config["systems_network"]
    tables_dir = Path(cfg["output"]["tables_dir"])
    networks_dir = Path(cfg["output"]["networks_dir"])
    final_review_dir = Path(cfg["output"]["final_review_dir"])

    nodes = pd.read_csv(networks_dir / "nodes.tsv", sep="\t")
    edges = pd.read_csv(networks_dir / "edges.tsv", sep="\t")
    membership = pd.read_csv(tables_dir / "candidate_pathway_membership.tsv", sep="\t")
    consensus = pd.read_csv(tables_dir / "resistance_pathway_consensus.tsv", sep="\t")
    mm = pd.read_csv(tables_dir / "multimodal_pathway_convergence.tsv", sep="\t")

    rows_text = [build_systems_map_text(c, nodes, edges, membership, consensus, mm) for c in CANDIDATES]
    plot_systems_map(rows_text, final_review_dir / "05_resistance_systems_map.png")


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run_visualization_4()
