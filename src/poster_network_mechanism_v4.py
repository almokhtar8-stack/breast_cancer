"""POST-FREEZE EXPLORATORY network/mechanism figure (v4) -- a VISUAL /
PRESENTATION REBUILD of the v2 network as a five-panel publication-style
figure: one local neighborhood panel per candidate (A KDM1A, B TLK2,
C USP34, D VEZF1) plus one candidate-connectivity summary panel (E).

No STRING requery, no new interaction data, no change to any network-
generation rule. This module imports `build_network()` from
`poster_network_mechanism_v2` and calls it unmodified -- the same 47-node
/ 147-edge / 3-component graph. Each panel's node set is derived
PROGRAMMATICALLY from that graph:

  local subgraph of candidate X =
      {X}
    | {direct neighbors of X in the v2 graph}                    (Level 1)
    | {v2 nodes of kind "level2_bridge" adjacent to >=1 of
       those direct neighbors}                                   (Level 2)
  with ALL v2 edges induced on that node set.

No gene list is hand-typed; the graph decides membership. VEZF1 has zero
neighbors at the frozen score>=0.7 threshold, so its panel honestly shows
a single node (the old VEZF1--DMTN pathway-co-membership relation is NOT
a STRING interaction and is not drawn). Panel E's shortest path is
computed at render time with networkx -- all STRING edges are UNDIRECTED
functional associations, so no arrows are drawn anywhere and no
activation/inhibition is implied.
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

from src.poster_network_mechanism_v2 import (
    BRIDGE_COLOR,
    CANDIDATES,
    DGRAY,
    FOCUS_COLORS,
    MGRAY,
    PARTNER_COLOR,
    build_network,
)

logger = logging.getLogger(__name__)

OUT_DIR = Path("results/figures/poster_network_mechanism_v4")

EDGE_COLOR = "#C4CCD5"

# KNOWN, ACCEPTED NONDETERMINISM: adjustText's collision solver does not
# produce identical label pixel positions across runs (verified by rendering
# repeatedly in one process with a fixed numpy seed and a fixed PYTHONHASHSEED
# -- the PNG still differs; with adjust_text disabled it is bit-identical every
# time). Removing the solver WAS tested and rejected: it reintroduces real label
# collisions (HDAC1 over the KDM1A node, HMG20B/CDH1, CTBP1/SNAI1), so
# readability is preferred over byte-reproducibility here.
#
# This affects LABEL POSITIONS ONLY. The graph itself -- nodes, edges, scores,
# components, shortest paths -- is fully deterministic and independently
# asserted in tests/test_poster_network_mechanism_v4.py. See
# tests/test_poster_release_integrity.py and poster/README.md.
PANEL_BORDER = "#DDDDDD"


def local_subgraph(graph: nx.Graph, candidate: str) -> nx.Graph:
    """The candidate's local neighborhood, derived purely from the v2
    graph: the candidate, its direct (Level-1) neighbors, and any v2
    Level-2 bridge node adjacent to at least one of those neighbors, with
    all v2 edges induced on that node set. Deterministic; nothing is
    hand-listed."""
    partners = set(graph.neighbors(candidate))
    level2 = {
        node for node in graph.nodes
        if graph.nodes[node]["kind"] == "level2_bridge"
        and any(graph.has_edge(node, partner) for partner in partners)
    }
    nodes = {candidate} | partners | level2
    return graph.subgraph(nodes).copy()


def candidate_shortest_paths(graph: nx.Graph) -> dict:
    """Candidate-pair connectivity computed directly from the v2 graph at
    render time: the networkx shortest path where one exists (plus the
    count of equally short alternatives, reported rather than hidden),
    None where the pair is disconnected."""
    out: dict = {}
    for i, a in enumerate(CANDIDATES):
        for b in CANDIDATES[i + 1:]:
            try:
                path = nx.shortest_path(graph, a, b)
                n_alternatives = len(list(nx.all_shortest_paths(graph, a, b)))
                out[(a, b)] = {"path": path, "n_edges": len(path) - 1, "n_equally_short": n_alternatives}
            except nx.NetworkXNoPath:
                out[(a, b)] = None
    return out


# ---------------------------------------------------------------------------
# Rendering
# ---------------------------------------------------------------------------

def _node_size(kind: str, degree_local: int, density_factor: float = 1.0) -> float:
    """Marker area; non-candidate nodes shrink with `density_factor` in
    crowded panels so a dense neighborhood (KDM1A) stays readable without
    dropping any node."""
    base = {"candidate": 1900, "level1_partner": 330, "level2_bridge": 150}[kind]
    if kind == "candidate":
        return base
    return base * density_factor * (1.0 + min(degree_local / 18.0, 0.6))


def _panel_frame(ax, letter: str, title: str) -> None:
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_color(PANEL_BORDER)
        spine.set_linewidth(0.8)
    ax.set_title(f"{letter}   {title}", loc="left", fontsize=13, fontweight="bold", color=DGRAY, pad=8)


def _draw_local_network(ax, sub: nx.Graph, candidate: str, label_all: bool) -> None:
    if sub.number_of_edges() > 0:
        weight = {(u, v): max(0.15, 1.0 - d["score"]) for u, v, d in sub.edges(data=True)}
        nx.set_edge_attributes(sub, weight, "layout_weight")
        pos = nx.kamada_kawai_layout(sub, weight="layout_weight")
        # Recenter so the candidate sits at the panel middle.
        shift = np.array(pos[candidate])
        pos = {n: np.array(p) - shift for n, p in pos.items()}
    else:
        pos = {candidate: np.array([0.0, 0.0])}

    xs = np.array([pos[n][0] for n in sub.nodes])
    ys = np.array([pos[n][1] for n in sub.nodes])
    span = max(xs.max() - xs.min(), ys.max() - ys.min(), 0.5)
    pad = span * 0.16
    ax.set_xlim(xs.min() - pad, xs.max() + pad)
    ax.set_ylim(ys.min() - pad, ys.max() + pad)
    ax.set_aspect("equal")

    # Dense neighborhoods get smaller non-candidate markers so no node or
    # edge is hidden -- panel size differences stay honest.
    density_factor = min(1.0, np.sqrt(14.0 / sub.number_of_nodes()))

    for u, v, d in sub.edges(data=True):
        lw = 1.1 if d["interaction_type"] == "physical_PPI" else 0.55
        ax.plot([pos[u][0], pos[v][0]], [pos[u][1], pos[v][1]], color=EDGE_COLOR,
                 linewidth=lw, alpha=0.65, zorder=1, solid_capstyle="round")

    local_degree = dict(sub.degree())
    others = [n for n in sub.nodes if sub.nodes[n]["kind"] != "candidate"]
    sizes = [_node_size(sub.nodes[n]["kind"], local_degree[n], density_factor) for n in others]
    colors = [PARTNER_COLOR if sub.nodes[n]["kind"] == "level1_partner" else BRIDGE_COLOR for n in others]
    ax.scatter([pos[n][0] for n in others], [pos[n][1] for n in others], s=sizes, c=colors,
                edgecolors="white", linewidths=0.7, zorder=3)
    # Candidate drawn last so it always sits on top of its dense neighborhood.
    ax.scatter([pos[candidate][0]], [pos[candidate][1]], s=_node_size("candidate", 0),
                c=FOCUS_COLORS[candidate], edgecolors="white", linewidths=1.8, zorder=4)
    ax.text(pos[candidate][0], pos[candidate][1], candidate, fontsize=11.5, fontweight="bold",
             color="white", ha="center", va="center", zorder=6)

    # Seed each label on the far side of its node from the candidate, so
    # labels of nodes hugging the candidate never start underneath it;
    # adjust_text then resolves any remaining collisions.
    cand_xy = np.array(pos[candidate])
    texts = []
    for node in others:
        kind = sub.nodes[node]["kind"]
        x, y = pos[node]
        if not label_all and kind == "level2_bridge" and local_degree[node] <= 1:
            continue
        fontsize = 8.2 if kind == "level1_partner" else 7.2
        r = np.sqrt(_node_size(kind, local_degree[node], density_factor) / np.pi) * span / 300
        away = np.array([x, y]) - cand_xy
        norm = np.linalg.norm(away)
        direction = away / norm if norm > 1e-9 else np.array([0.0, -1.0])
        seed = np.array([x, y]) + direction * (r + span * 0.03)
        va = "bottom" if direction[1] > 0.3 else ("top" if direction[1] < -0.3 else "center")
        texts.append(ax.text(seed[0], seed[1], node, fontsize=fontsize, color=DGRAY,
                              ha="center", va=va, zorder=5))
    if texts:
        adjust_text(texts, x=list(xs), y=list(ys), ax=ax, expand=(1.2, 1.45),
                    force_text=(0.3, 0.5), force_static=(0.25, 0.4),
                    arrowprops=dict(arrowstyle="-", color=MGRAY, lw=0.5, alpha=0.7, shrinkA=0, shrinkB=1))


def _draw_vezf1_panel(ax, sub: nx.Graph) -> None:
    assert sub.number_of_nodes() == 1 and sub.number_of_edges() == 0
    ax.set_xlim(-1, 1)
    ax.set_ylim(-1, 1)
    ax.set_aspect("equal")
    ax.scatter([0], [0.18], s=1900, c=FOCUS_COLORS["VEZF1"], edgecolors="white", linewidths=1.8, zorder=3)
    ax.text(0, 0.18, "VEZF1", fontsize=11.5, fontweight="bold", color="white", ha="center", va="center", zorder=5)
    ax.text(0, -0.42, "No high-confidence STRING partners\ndetected (score ≥ 0.7)",
             fontsize=9.5, color=MGRAY, ha="center", va="center", style="italic")


def _draw_connectivity_panel(ax, graph: nx.Graph, paths: dict) -> None:
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)

    result = paths[("KDM1A", "USP34")]
    path = result["path"]
    n_nodes = len(path)
    x_positions = np.linspace(0.9, 6.1, n_nodes)
    y_path = 6.4

    for x0, x1 in zip(x_positions[:-1], x_positions[1:]):
        ax.plot([x0, x1], [y_path, y_path], color="#9AA6B2", linewidth=2.2, zorder=1, solid_capstyle="round")

    for x, node in zip(x_positions, path):
        is_candidate = graph.nodes[node]["kind"] == "candidate"
        color = FOCUS_COLORS[node] if is_candidate else PARTNER_COLOR
        size = 2000 if is_candidate else 900
        ax.scatter([x], [y_path], s=size, c=color, edgecolors="white",
                    linewidths=1.6 if is_candidate else 0.9, zorder=3)
        if is_candidate:
            ax.text(x, y_path, node, fontsize=10.5, fontweight="bold", color="white",
                     ha="center", va="center", zorder=5)
        else:
            ax.text(x, y_path - 0.95, node, fontsize=9.5, color=DGRAY, ha="center", va="top", zorder=5)

    ax.text((x_positions[0] + x_positions[-1]) / 2, y_path - 2.5,
             f"Shortest path: {result['n_edges']} edges\n"
             f"({result['n_equally_short']} equally short routes, all via {path[1]})",
             fontsize=9.5, color=MGRAY, ha="center", va="top", style="italic")

    ax.scatter([7.9], [7.1], s=1500, c=FOCUS_COLORS["TLK2"], edgecolors="white", linewidths=1.6, zorder=3)
    ax.text(7.9, 7.1, "TLK2", fontsize=10, fontweight="bold", color="white", ha="center", va="center", zorder=5)
    ax.text(7.9, 5.9, "Separate connected\ncomponent", fontsize=9.5, color=MGRAY,
             ha="center", va="top", style="italic")

    ax.scatter([7.9], [3.0], s=1500, c=FOCUS_COLORS["VEZF1"], edgecolors="white", linewidths=1.6, zorder=3)
    ax.text(7.9, 3.0, "VEZF1", fontsize=10, fontweight="bold", color="white", ha="center", va="center", zorder=5)
    ax.text(7.9, 1.8, "Isolated at this\nthreshold", fontsize=9.5, color=MGRAY,
             ha="center", va="top", style="italic")


def build_network_mechanism_main(stub: Path) -> nx.Graph:
    graph = build_network()
    paths = candidate_shortest_paths(graph)

    fig = plt.figure(figsize=(16.5, 10.5), dpi=300)
    gs = fig.add_gridspec(2, 6, hspace=0.28, wspace=0.35,
                          left=0.035, right=0.975, top=0.865, bottom=0.05)
    ax_a = fig.add_subplot(gs[0, 0:2])
    ax_b = fig.add_subplot(gs[0, 2:4])
    ax_c = fig.add_subplot(gs[0, 4:6])
    ax_d = fig.add_subplot(gs[1, 0:2])
    ax_e = fig.add_subplot(gs[1, 2:6])

    sub_kdm1a = local_subgraph(graph, "KDM1A")
    sub_tlk2 = local_subgraph(graph, "TLK2")
    sub_usp34 = local_subgraph(graph, "USP34")
    sub_vezf1 = local_subgraph(graph, "VEZF1")

    _panel_frame(ax_a, "A", "KDM1A local neighborhood")
    _draw_local_network(ax_a, sub_kdm1a, "KDM1A", label_all=False)

    _panel_frame(ax_b, "B", "TLK2 local neighborhood")
    _draw_local_network(ax_b, sub_tlk2, "TLK2", label_all=True)

    _panel_frame(ax_c, "C", "USP34 local neighborhood")
    _draw_local_network(ax_c, sub_usp34, "USP34", label_all=True)

    _panel_frame(ax_d, "D", "VEZF1 local neighborhood")
    _draw_vezf1_panel(ax_d, sub_vezf1)

    _panel_frame(ax_e, "E", "Candidate connectivity")
    _draw_connectivity_panel(ax_e, graph, paths)

    fig.text(0.035, 0.975, "Candidate-specific molecular neighborhoods and connectivity",
              fontsize=21, fontweight="bold", color=DGRAY, ha="left", va="top")
    fig.text(0.035, 0.938, "High-confidence STRING associations (score ≥ 0.7); post-freeze exploratory analysis.",
              fontsize=11, color="#555555", ha="left", va="top")

    legend_handles = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor=DGRAY, markeredgecolor="white",
               markersize=12, label="Candidate"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor=PARTNER_COLOR, markeredgecolor="white",
               markersize=8, label="Direct partner"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor=BRIDGE_COLOR, markeredgecolor="white",
               markersize=6, label="Level-2 / bridge"),
    ]
    fig.legend(handles=legend_handles, loc="upper right", frameon=False, fontsize=9.5,
               labelcolor=DGRAY, ncol=3, bbox_to_anchor=(0.975, 0.985), handletextpad=0.5, columnspacing=1.2)

    stub.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(stub.with_suffix(".png"), facecolor="white", bbox_inches="tight", dpi=300)
    fig.savefig(stub.with_suffix(".pdf"), facecolor="white", bbox_inches="tight")
    fig.savefig(stub.with_suffix(".svg"), facecolor="white", bbox_inches="tight")
    plt.close(fig)
    logger.info("wrote %s.png/.pdf/.svg", stub)
    return graph


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    build_network_mechanism_main(OUT_DIR / "NETWORK_mechanism_v4")
