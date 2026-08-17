"""POST-FREEZE EXPLORATORY network/mechanism figure (v3) -- VISUAL
REFINEMENT ONLY of v2. This module imports `build_network()` and
`network_stats()` directly from `poster_network_mechanism_v2` and calls
them unmodified: identical nodes, identical edges, identical pathway
annotations, identical STRING source. Nothing here recomputes, filters,
or reinterprets any part of the network. Only the rendering changes:
component-aware compact layout, larger candidate/hub emphasis, reduced
label set, quieter edges, one consolidated legend, and a shorter title.

Reused unmodified from v2: `CANDIDATES`, `FOCUS_COLORS`, `PATHWAY_COLORS`,
`CANONICAL_BIOLOGY_GENES`, `PARTNER_COLOR`, `BRIDGE_COLOR`, `DGRAY`,
`MGRAY`.
"""

from __future__ import annotations

import logging
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
from adjustText import adjust_text
from matplotlib.lines import Line2D
from matplotlib.patches import Wedge

from src.poster_network_mechanism_v2 import (
    CANDIDATES,
    CANONICAL_BIOLOGY_GENES,
    DGRAY,
    FOCUS_COLORS,
    MGRAY,
    PARTNER_COLOR,
    BRIDGE_COLOR,
    PATHWAY_COLORS,
    build_network,
    network_stats,
)

logger = logging.getLogger(__name__)

OUT_DIR = Path("results/figures/poster_network_mechanism_v3")


# ---------------------------------------------------------------------------
# Compact, component-aware layout -- the graph topology within each
# connected component is a genuine networkx Kamada-Kawai layout; only the
# TRANSLATION and SCALE of each component as a whole is chosen manually,
# to pack the real three-component structure onto the canvas instead of
# letting Kamada-Kawai (which has no basis for placing disconnected
# components relative to each other) scatter them arbitrarily. No edge is
# added, removed, or rerouted by this step.
# ---------------------------------------------------------------------------

def _normalize(pos: dict) -> dict:
    arr = np.array(list(pos.values()))
    center = arr.mean(axis=0)
    scale = np.abs(arr - center).max() or 1.0
    return {k: (v - center) / scale for k, v in pos.items()}


def compact_layout(graph: nx.Graph) -> dict:
    components = sorted(nx.connected_components(graph), key=len, reverse=True)
    main_nodes = components[0]
    other_components = components[1:]
    # Identify the TLK2-containing component and the true singleton(s) by
    # membership, not by an assumed order -- robust to any future change
    # in component sizes.
    tlk2_nodes = next((c for c in other_components if "TLK2" in c), other_components[0])
    singleton_components = [c for c in other_components if c is not tlk2_nodes]

    sub_main = graph.subgraph(main_nodes)
    pos_main = _normalize(nx.kamada_kawai_layout(sub_main, weight="layout_weight"))

    sub_tlk2 = graph.subgraph(tlk2_nodes)
    pos_tlk2 = _normalize(nx.kamada_kawai_layout(sub_tlk2, weight="layout_weight"))

    pos: dict = {}
    for node, (x, y) in pos_main.items():
        pos[node] = np.array([x * 1.0, y * 0.78])

    tlk2_scale, tlk2_center = 0.32, np.array([-0.32, 0.82])
    for node, (x, y) in pos_tlk2.items():
        pos[node] = tlk2_center + np.array([x, y]) * tlk2_scale

    singleton_anchors = [np.array([0.98, -0.10]), np.array([1.28, -0.62])]
    for comp, anchor in zip(singleton_components, singleton_anchors):
        for node in comp:
            pos[node] = anchor

    return pos


# ---------------------------------------------------------------------------
# Node sizing / color -- same underlying degree/betweenness metrics v2
# already computed and stored on the graph; only the visual scale
# constants change, tuned for stronger candidate/hub dominance at poster
# distance.
# ---------------------------------------------------------------------------

def _node_size(kind: str, degree: int, betweenness: float) -> float:
    base = {"candidate": 6600, "level1_partner": 480, "level2_bridge": 190}[kind]
    hub_boost = 1.0 + min(degree / 12.0, 1.6) + min(betweenness * 7.0, 1.6)
    return base * (hub_boost if kind != "candidate" else 1.0 + min(betweenness * 2.5, 0.45))


def _node_color(node: str, kind: str) -> str:
    if kind == "candidate":
        return FOCUS_COLORS[node]
    if kind == "level1_partner":
        return PARTNER_COLOR
    return BRIDGE_COLOR


def _should_label(node: str, data: dict) -> bool:
    """Identical rule to v2: candidates, every Level-1 partner, and
    Level-2 nodes that are either a canonical biology marker gene or a
    real hub (degree >= 6) -- unchanged from v2, reused here to keep the
    same set of genes considered "informative" across both figures."""
    if data["kind"] in ("candidate", "level1_partner"):
        return True
    return node in CANONICAL_BIOLOGY_GENES or data["degree"] >= 6


def build_network_mechanism_main(stub: Path) -> nx.Graph:
    graph = build_network()
    weight = {(u, v): max(0.15, 1.0 - d["score"]) for u, v, d in graph.edges(data=True)}
    nx.set_edge_attributes(graph, weight, "layout_weight")
    pos = compact_layout(graph)

    fig, ax = plt.subplots(figsize=(17.5, 13.0), dpi=300)
    ax.set_facecolor("white")
    fig.patch.set_facecolor("white")
    fig.subplots_adjust(left=0.005, right=0.995, top=0.885, bottom=0.01)

    # Edges deliberately quiet: thin, light, low-alpha functional
    # associations; only physical_PPI edges (real STRING evidence, not a
    # hand-picked "mechanistic path") are drawn a touch stronger.
    for u, v, d in graph.edges(data=True):
        lw = (1.3 if d["interaction_type"] == "physical_PPI" else 0.55) * (0.7 + 0.5 * d["score"])
        alpha = 0.38 if d["hop"] == 1 else 0.22
        ax.plot([pos[u][0], pos[v][0]], [pos[u][1], pos[v][1]], color="#C7CFD8",
                 linewidth=lw, alpha=alpha, zorder=1, solid_capstyle="round")

    xs = [pos[n][0] for n in graph.nodes]
    ys = [pos[n][1] for n in graph.nodes]
    ax.set_xlim(min(xs) - 0.16, max(xs) + 0.16)
    ax.set_ylim(min(ys) - 0.14, max(ys) + 0.16)
    ax.set_aspect("equal")
    ax.axis("off")

    fig.canvas.draw()
    bbox_px = ax.get_window_extent()
    x0, x1 = ax.get_xlim()
    data_per_point = (x1 - x0) / bbox_px.width * (fig.dpi / 72.0)

    def node_radius_data(kind: str, degree: int, betweenness: float) -> float:
        return np.sqrt(_node_size(kind, degree, betweenness) / np.pi) * data_per_point

    # Pathway-membership halo: kept subtle (thin ring, muted via alpha) so
    # candidate identity colors and hub sizing stay dominant.
    for node, data in graph.nodes(data=True):
        pathways = data["pathways"]
        if not pathways:
            continue
        r = node_radius_data(data["kind"], data["degree"], data["betweenness"]) * 1.28
        n = len(pathways)
        for i, pw in enumerate(pathways):
            theta1, theta2 = 360 * i / n, 360 * (i + 1) / n
            wedge = Wedge(pos[node], r, theta1, theta2, width=r * 0.22,
                           facecolor=PATHWAY_COLORS[pw], edgecolor="none", alpha=0.85, zorder=2)
            ax.add_patch(wedge)

    sizes = [_node_size(d["kind"], d["degree"], d["betweenness"]) for _, d in graph.nodes(data=True)]
    colors = [_node_color(n, d["kind"]) for n, d in graph.nodes(data=True)]
    linewidths = [2.6 if d["kind"] == "candidate" else 0.8 for _, d in graph.nodes(data=True)]
    ax.scatter(xs, ys, s=sizes, c=colors, edgecolors="white", linewidths=linewidths, zorder=3)

    texts = []
    for node, data in graph.nodes(data=True):
        if not _should_label(node, data):
            continue
        x, y = pos[node]
        if data["kind"] == "candidate":
            ax.text(x, y, node, fontsize=18, fontweight="bold", color="white", ha="center", va="center", zorder=5)
        else:
            is_hub = data["degree"] >= 10 or node in CANONICAL_BIOLOGY_GENES
            fontsize = 12.5 if is_hub else 9.0
            fontweight = "bold" if is_hub else "normal"
            r = node_radius_data(data["kind"], data["degree"], data["betweenness"])
            texts.append(ax.text(x, y - r * 1.25 - 0.008, node, fontsize=fontsize, fontweight=fontweight,
                                  color=DGRAY, ha="center", va="top", zorder=5))

    adjust_text(texts, x=xs, y=ys, ax=ax, expand=(1.25, 1.5), force_text=(0.3, 0.55), force_static=(0.25, 0.45),
                arrowprops=dict(arrowstyle="-", color=MGRAY, lw=0.6, alpha=0.8, shrinkA=0, shrinkB=1))

    # Short micro-annotations: one line per component, placed in the
    # nearby whitespace, not paragraphs.
    ax.text(-0.05, -0.98, "KDM1A–USP34: connected component",
             fontsize=10.5, color=MGRAY, ha="center", va="top", style="italic")
    tlk2_center = np.array([-0.32, 0.82])
    ax.text(tlk2_center[0], tlk2_center[1] - 0.50, "TLK2: separate chromatin neighborhood",
             fontsize=10.5, color=MGRAY, ha="center", va="top", style="italic")
    singleton_positions = [pos[n] for n in graph.nodes if graph.degree(n) == 0]
    if singleton_positions:
        vx, vy = singleton_positions[0]
        ax.text(vx, vy - 0.20, "VEZF1: no high-confidence\nSTRING partners",
                 fontsize=10.5, color=MGRAY, ha="center", va="top", style="italic")

    fig.text(0.015, 0.975, "Molecular networks reveal distinct candidate neighborhoods",
              fontsize=24, fontweight="bold", color=DGRAY, ha="left", va="top")
    fig.text(0.015, 0.947, "High-confidence STRING associations (score ≥ 0.7); post-freeze exploratory analysis.",
              fontsize=12, color="#555555", ha="left", va="top")

    # One consolidated legend block, tucked into the (real, unused)
    # whitespace to the right of the singleton anchors.
    legend_handles = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor=FOCUS_COLORS[c], markeredgecolor="white",
               markersize=15, label=c) for c in CANDIDATES
    ]
    legend_handles += [
        Line2D([0], [0], marker="o", color="none", markerfacecolor=PARTNER_COLOR, markeredgecolor="white",
               markersize=9, label="Level-1 partner"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor=BRIDGE_COLOR, markeredgecolor="white",
               markersize=6, label="Level-2 bridge"),
        Line2D([0], [0], color="#C7CFD8", linewidth=1.3, label="physical PPI"),
        Line2D([0], [0], color="#C7CFD8", linewidth=0.55, alpha=0.6, label="functional association"),
    ]
    legend_handles += [
        Line2D([0], [0], marker="o", color="none", markerfacecolor=color, markeredgecolor="none",
               markersize=9, label=label) for label, color in PATHWAY_COLORS.items()
    ]
    ax.legend(handles=legend_handles, loc="upper right", frameon=False, fontsize=9.6, labelcolor=DGRAY,
              ncol=1, bbox_to_anchor=(1.005, 1.0), title="network  /  program membership", title_fontsize=9.6,
              handletextpad=0.6, labelspacing=0.55)

    stub.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(stub.with_suffix(".png"), facecolor="white", bbox_inches="tight", dpi=300)
    fig.savefig(stub.with_suffix(".pdf"), facecolor="white", bbox_inches="tight")
    fig.savefig(stub.with_suffix(".svg"), facecolor="white", bbox_inches="tight")
    plt.close(fig)
    logger.info("wrote %s.png/.pdf/.svg (%d nodes, %d edges)", stub, graph.number_of_nodes(), graph.number_of_edges())
    return graph


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    build_network_mechanism_main(OUT_DIR / "NETWORK_mechanism_v3")
