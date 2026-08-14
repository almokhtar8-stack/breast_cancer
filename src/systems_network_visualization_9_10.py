"""Comparative systems-network audit figures 09 and 10.

Figure 09: a clean head-to-head scorecard comparing the four frozen
candidates (own five-layer evidence heatmap + direct-neighbor counts).
Figure 10: a deliberately small, hand-laid-out mechanism map --
candidate -> supported neighbor/bridge -> resistance pathway/module --
using only connections already established as real in this audit (Parts
2-5). Reads only already-written tables from this audit and the frozen
systems-network node table; does not touch any other figure.
"""

from __future__ import annotations

import logging
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

logger = logging.getLogger(__name__)

CANDIDATES = ["USP34", "VEZF1", "EML5", "CITED2"]

NODES_PATH = Path("results/networks/systems_network/cytoscape/network_nodes.tsv")
NEIGHBORS_PATH = Path("results/tables/systems_network/four_candidate_direct_neighbors.tsv")
AUDIT_PATH = Path("results/tables/systems_network/four_candidate_network_audit.tsv")

OUT_FIG_09 = Path("results/figures/systems_network/final_review/09_four_candidate_network_comparison.png")
OUT_FIG_10 = Path("results/figures/systems_network/final_review/10_four_candidate_mechanism_map.png")

DIVERGING_CMAP = LinearSegmentedColormap.from_list("audit_diverging", ["#1c5cab", "#f0efec", "#e34948"])
VABS = 1.75

EVIDENCE_COLUMNS = [
    ("crispr_effect", "crispr_fdr", "CRISPR\neffect_size"),
    ("gse118713_log2fc", "gse118713_fdr", "GSE118713\n(TAMR vs MCF7)"),
    ("gse240112_log2fc", "gse240112_fdr", "GSE240112\n(recurrent vs\nprimary)"),
    ("gse111151_log2fc", "gse111151_fdr", "GSE111151\n(resistant vs\nparental)"),
    ("gse245601_acute_log2fc", "gse245601_acute_fdr", "GSE245601\nACUTE 12h\n(not resistance)"),
]


def build_figure_09(nodes: pd.DataFrame, neighbors: pd.DataFrame, audit: pd.DataFrame, out_fig: Path = OUT_FIG_09) -> None:
    nodes_idx = nodes.set_index("gene")
    n_neighbors = {c: int(audit.set_index("candidate").loc[c, "n_direct_neighbors"]) for c in CANDIDATES}
    classification = {c: audit.set_index("candidate").loc[c, "systems_mechanism_classification"] for c in CANDIDATES}

    fig = plt.figure(figsize=(12, 8.5), dpi=200)
    fig.patch.set_facecolor("white")
    gs = fig.add_gridspec(2, 2, height_ratios=[3, 1.3], width_ratios=[4, 1.7], wspace=0.08, hspace=0.55)
    ax_heat = fig.add_subplot(gs[0, 0])
    ax_text = fig.add_subplot(gs[0, 1], sharey=ax_heat)
    ax_bar = fig.add_subplot(gs[1, :])

    values = np.zeros((len(CANDIDATES), len(EVIDENCE_COLUMNS)))
    for j, (val_col, _, _) in enumerate(EVIDENCE_COLUMNS):
        for i, c in enumerate(CANDIDATES):
            values[i, j] = nodes_idx.loc[c, val_col]

    im = ax_heat.imshow(values, cmap=DIVERGING_CMAP, vmin=-VABS, vmax=VABS, aspect="auto")
    ax_heat.axvline(x=0.5, color="white", linewidth=3)
    ax_heat.axvline(x=3.5, color="#111111", linewidth=2.5)

    for i, c in enumerate(CANDIDATES):
        for j, (val_col, fdr_col, _) in enumerate(EVIDENCE_COLUMNS):
            val = nodes_idx.loc[c, val_col]
            fdr = nodes_idx.loc[c, fdr_col]
            marker = "**" if pd.notna(fdr) and fdr < 0.05 else ""
            text_color = "white" if abs(val) > VABS * 0.4 else "#111111"
            ax_heat.text(j, i, f"{val:+.2f}{marker}", ha="center", va="center", fontsize=11, color=text_color, fontweight="bold" if marker else "normal")

    ax_heat.set_xticks(range(len(EVIDENCE_COLUMNS)))
    ax_heat.set_xticklabels([c[2] for c in EVIDENCE_COLUMNS], fontsize=9)
    ax_heat.set_yticks(range(len(CANDIDATES)))
    ax_heat.set_yticklabels(CANDIDATES, fontsize=13, fontweight="bold")
    ax_heat.tick_params(length=0)
    for spine in ax_heat.spines.values():
        spine.set_visible(False)
    ax_heat.set_title("Own-candidate evidence across five layers (** FDR<0.05)", fontsize=11, loc="left", fontweight="bold")

    cbar = fig.colorbar(im, ax=ax_heat, fraction=0.046, pad=0.03)
    cbar.ax.set_title("effect /\nlog2FC", fontsize=7.5, pad=6)
    cbar.ax.tick_params(labelsize=8)

    ax_text.set_xlim(0, 1)
    ax_text.axis("off")
    for i, c in enumerate(CANDIDATES):
        label = classification[c].split("_", 1)[1].replace("_", " ")
        wrapped = "\n".join([label[:len(label) // 2].rsplit(" ", 1)[0], label[len(label[: len(label) // 2].rsplit(" ", 1)[0]) + 1 :]]) if len(label) > 16 else label
        ax_text.text(0.28, i, wrapped, ha="left", va="center", fontsize=9, color="#111111", fontweight="bold")

    bar_colors = ["#256abf", "#6da7ec", "#cde2fb", "#1c5cab"]
    bars = ax_bar.bar(CANDIDATES, [n_neighbors[c] for c in CANDIDATES], color=bar_colors, width=0.55)
    for rect, c in zip(bars, CANDIDATES):
        ax_bar.text(rect.get_x() + rect.get_width() / 2, rect.get_height() + 0.4, str(n_neighbors[c]), ha="center", va="bottom", fontsize=12, fontweight="bold")
    ax_bar.set_ylabel("direct (1-hop)\nnetwork neighbors", fontsize=9)
    ax_bar.set_ylim(0, max(n_neighbors.values()) + 3)
    for spine in ["top", "right"]:
        ax_bar.spines[spine].set_visible(False)
    ax_bar.set_title("Direct network neighborhood size (frozen network, required_score>=0.7 / TRRUST / pathway co-membership)", fontsize=10, loc="left")

    fig.suptitle("Four-candidate systems-network head-to-head", fontsize=14, fontweight="bold", y=0.99)
    out_fig.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_fig, facecolor="white")
    plt.close(fig)
    logger.info("wrote %s", out_fig)


# ---------------------------------------------------------------------------
# Figure 10: small mechanism map
# ---------------------------------------------------------------------------

NODE_POS = {
    "USP34": (-8.0, 0.0, "candidate"),
    "CTNNB1": (-5.3, -3.3, "gene"),
    "GOBP_CANONICAL_WNT_\nSIGNALING_PATHWAY": (-8.6, 3.6, "pathway"),
    "HALLMARK_WNT_BETA_\nCATENIN_SIGNALING": (-2.6, -5.8, "pathway"),
    "VEZF1": (-1.2, 5.6, "candidate"),
    "DMTN": (-6.0, 7.6, "gene"),
    "HALLMARK_HEME_\nMETABOLISM": (-6.3, 10.4, "pathway"),
    "GOBP_BLOOD_VESSEL_\nMORPHOGENESIS": (2.8, 8.4, "pathway"),
    "CITED2": (6.8, 0.0, "candidate"),
    "HALLMARK_UV_\nRESPONSE_DN": (9.5, -3.6, "pathway"),
    "TFAP2C": (9.8, 2.6, "gene"),
    "HALLMARK_ESTROGEN_\nRESPONSE_EARLY": (14.6, 3.6, "pathway"),
    "EP300": (6.2, 4.6, "gene_hub_caution"),
    "EML5": (-8.0, -7.0, "candidate_isolated"),
}

# (source, target, edge_type, label, label_t, label_offset) -- label_t in
# [0,1] controls where along the edge the annotation sits (0=at source, 1=at
# target, 0.5=midpoint); label_offset is an additional (dx, dy) nudge so a
# label can be steered clear of a node's own name-label bounding box.
EDGES = [
    ("USP34", "GOBP_CANONICAL_WNT_\nSIGNALING_PATHWAY", "pathway_module", "direct member +\nleading-edge", 0.5, (0, 0)),
    ("USP34", "CTNNB1", "indirect_2hop", "via ubiquitin bridge\n(RPS27A/UBB/UBC -- generic,\nno independent support)", 0.5, (0, 0)),
    ("CTNNB1", "HALLMARK_WNT_BETA_\nCATENIN_SIGNALING", "pathway_module", "leading-edge", 0.5, (0, 0)),
    ("VEZF1", "GOBP_BLOOD_VESSEL_\nMORPHOGENESIS", "pathway_module", "direct member +\nleading-edge", 0.55, (0, 0)),
    ("VEZF1", "HALLMARK_HEME_\nMETABOLISM", "pathway_module", "direct member +\nleading-edge", 0.5, (0, 0)),
    ("VEZF1", "DMTN", "pathway_co_membership", "A-tier: DMTN independently\nFDR<0.001 (2 datasets)", 0.5, (-2.0, -0.6)),
    ("DMTN", "HALLMARK_HEME_\nMETABOLISM", "pathway_module", "leading-edge", 0.5, (0, 0)),
    ("CITED2", "GOBP_BLOOD_VESSEL_\nMORPHOGENESIS", "pathway_module", "direct member +\nleading-edge (CONVERGENCE)", 0.62, (0, 0)),
    ("CITED2", "HALLMARK_UV_\nRESPONSE_DN", "pathway_module", "direct + leading-edge +\nMULTIMODAL (RNA+CRISPR)", 0.5, (0, 0)),
    ("CITED2", "TFAP2C", "physical_PPI", "A-tier: CRISPR FDR=0.048\n+ RNA FDR=0.001", 0.42, (0, 0)),
    ("TFAP2C", "HALLMARK_ESTROGEN_\nRESPONSE_EARLY", "pathway_module", "leading-edge", 0.5, (0.1, 0.9)),
    ("CITED2", "EP300", "physical_PPI", "A-tier RNA support, but\n#3 highest-degree hub\nin whole network", 0.62, (0, 0)),
]

EDGE_STYLE = {
    "physical_PPI": dict(color="#111111", linestyle="-", linewidth=2.2),
    "functional_association": dict(color="#555555", linestyle="--", linewidth=1.8),
    "pathway_co_membership": dict(color="#7a4fa3", linestyle="-.", linewidth=1.8),
    "indirect_2hop": dict(color="#b0392f", linestyle=":", linewidth=2.2),
    "pathway_module": dict(color="#1c5cab", linestyle="-", linewidth=1.4),
}

NODE_STYLE = {
    "candidate": dict(color="#0072B2", size=2600, shape="o"),
    "candidate_isolated": dict(color="#9a9a9a", size=1800, shape="o"),
    "gene": dict(color="#E69F00", size=1500, shape="o"),
    "gene_hub_caution": dict(color="#E69F00", size=1500, shape="o"),
    "pathway": dict(color="#009E73", size=0, shape="box"),
}


def _draw_pathway_box(ax, x, y, label):
    ax.add_patch(FancyBboxPatch((x - 1.35, y - 0.42), 2.7, 0.84, boxstyle="round,pad=0.08", linewidth=1.6, edgecolor="#009E73", facecolor="#e7f7f2", zorder=3))
    ax.text(x, y, label, fontsize=7.6, ha="center", va="center", color="#0b5f45", fontweight="bold", zorder=4)


def build_figure_10(out_fig: Path = OUT_FIG_10) -> None:
    fig, ax = plt.subplots(figsize=(15, 13), dpi=200)
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")

    for u, v, etype, label, label_t, (odx, ody) in EDGES:
        x1, y1, _ = NODE_POS[u]
        x2, y2, _ = NODE_POS[v]
        style = EDGE_STYLE[etype]
        arrowstyle = "-" if etype == "pathway_module" else "-|>"
        ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle=arrowstyle, mutation_scale=16, color=style["color"], linestyle=style["linestyle"], linewidth=style["linewidth"], zorder=1, shrinkA=28, shrinkB=28, capstyle="round"))
        xm, ym = x1 + (x2 - x1) * label_t + odx, y1 + (y2 - y1) * label_t + ody
        ax.text(xm, ym, label, fontsize=7.2, color=style["color"], ha="center", va="center", backgroundcolor="white", zorder=2, style="italic")

    for name, (x, y, kind) in NODE_POS.items():
        if kind == "pathway":
            _draw_pathway_box(ax, x, y, name)
            continue
        style = NODE_STYLE[kind]
        edgecolor = "#b0392f" if kind == "gene_hub_caution" else "white"
        lw = 3 if kind == "gene_hub_caution" else 1.5
        linestyle = "--" if kind == "candidate_isolated" else "-"
        ax.scatter([x], [y], s=style["size"], c=style["color"], edgecolors=edgecolor, linewidths=lw, linestyle=linestyle, zorder=3)
        label = name.replace("\n", " ")
        fontsize = 12.5 if kind.startswith("candidate") else 10
        r = math.hypot(x, y) or 1
        lx, ly = x + (x / r) * 0.95, y + (y / r) * 0.95
        ax.text(lx, ly, label, fontsize=fontsize, fontweight="bold", ha="center", va="center", color="#111111", zorder=5, bbox=dict(boxstyle="round,pad=0.12", facecolor="white", edgecolor="none", alpha=0.9))

    ax.annotate("NO RESOLVED NETWORK\nNEIGHBOURHOOD IN\nCURRENT ANALYSIS", xy=NODE_POS["EML5"][:2], xytext=(NODE_POS["EML5"][0], NODE_POS["EML5"][1] - 1.7), fontsize=8.5, ha="center", color="#555555", style="italic")

    legend_elements = [
        Line2D([0], [0], color=EDGE_STYLE["physical_PPI"]["color"], lw=2.2, linestyle="-", label="direct physical_PPI (STRING)"),
        Line2D([0], [0], color=EDGE_STYLE["pathway_co_membership"]["color"], lw=1.8, linestyle="-.", label="direct pathway_co_membership"),
        Line2D([0], [0], color=EDGE_STYLE["indirect_2hop"]["color"], lw=2.2, linestyle=":", label="indirect 2-hop (via generic bridge)"),
        Line2D([0], [0], color=EDGE_STYLE["pathway_module"]["color"], lw=1.4, linestyle="-", label="curated pathway / module membership"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor="#0072B2", markersize=15, label="candidate"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor="#E69F00", markersize=12, label="connector gene (A-tier data-supported)"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor="#E69F00", markeredgecolor="#b0392f", markeredgewidth=3, markersize=12, label="connector gene -- generic-hub caution"),
        Line2D([0], [0], marker="s", color="w", markerfacecolor="#e7f7f2", markeredgecolor="#009E73", markersize=13, label="pathway / module (STRONG_CONSENSUS)"),
    ]
    fig.legend(handles=legend_elements, loc="lower center", bbox_to_anchor=(0.5, 0.0), fontsize=9, frameon=False, ncol=3)

    fig.suptitle("Four-candidate mechanism map -- only connections established in this audit (EML5 has none; not forced onto the map)", fontsize=13, fontweight="bold", x=0.02, ha="left", y=0.995)
    ax.set_xlim(-11, 15.5)
    ax.set_ylim(-9.5, 11.5)
    ax.set_aspect("equal")
    ax.axis("off")

    fig.subplots_adjust(left=0.02, right=0.98, top=0.94, bottom=0.14)
    out_fig.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_fig, facecolor="white")
    plt.close(fig)
    logger.info("wrote %s", out_fig)


def run() -> None:
    nodes = pd.read_csv(NODES_PATH, sep="\t")
    neighbors = pd.read_csv(NEIGHBORS_PATH, sep="\t")
    audit = pd.read_csv(AUDIT_PATH, sep="\t")
    build_figure_09(nodes, neighbors, audit, OUT_FIG_09)
    build_figure_10(OUT_FIG_10)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run()
