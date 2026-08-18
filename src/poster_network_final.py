"""Poster figure 5 (final): what connects the candidates, and how weakly.

post_freeze_exploratory. No new STRING query and no new interaction data: the
graph is built by `poster_network_mechanism_v2.build_network()`, called
unmodified -- the same pinned query (species 9606, combined score >= 0.70),
the same 47 nodes, 147 edges and 3 components, from the same committed
`data/reference/interactions/string_v2_*` files. Values are verified before
anything is drawn.

WHAT THIS FIXES. The frozen figure 03 is crowded and its label placement uses
adjustText, a collision solver whose output moves between runs -- which is why
that figure's PNG is documented as not byte-reproducible. Here the layout is
seeded, every label sits at a fixed offset, no solver is used, and the PNG is
byte-identical across renders (test-enforced). Label density is cut by naming
only the candidates, the bridge, and the ubiquitin nodes.

WHAT PANEL B IS FOR. It is the honest reading of the connection. KDM1A and
USP34 are three edges apart, and there are four equally short routes -- but
every one of them passes through DNMT1 and then through a ubiquitin-encoding
gene (UBB, UBC, UBA52, RPS27A). The four routes are therefore not four
independent lines of support; they are one route whose middle position is
filled by four members of the same small gene family. Drawing that is what
stops "connected in a network" from sounding like evidence of a mechanism.

WHAT THIS FIGURE DOES NOT CLAIM. STRING edges are undirected functional
associations -- not activation, not inhibition, not necessarily physical
binding, and not a communication channel. No arrows are drawn anywhere. A
shortest path is not a mechanism. TLK2 sits in a separate component and VEZF1
has no partner at all at this threshold; both facts are shown, not hidden.
"""

from __future__ import annotations

import logging
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import networkx as nx

from src.poster_final_common import (
    FONT,
    OUT_DIR,
    figure_footer,
    headline,
    pin_reproducibility,
    save,
    verify,
)
from src.poster_network_mechanism_v2 import build_network
from src.poster_palette import GENE_COLOURS, NEUTRAL, WHITE

logger = logging.getLogger(__name__)

FIGURE = "F5_network_connectivity"
CANDIDATES = ("KDM1A", "TLK2", "USP34", "VEZF1")
LAYOUT_SEED = 20260818          # fixed: the layout must not move between runs
BRIDGE = "DNMT1"
UBIQUITIN = ("UBB", "UBC", "UBA52", "RPS27A")

REFERENCE = {"n_nodes": 47, "n_edges": 147, "n_components": 3,
             "path_length": 3, "n_shortest_paths": 4, "vezf1_degree": 0}


def component_layout(G):
    """Deterministic layout: each connected component is spring-laid with a
    fixed seed, rescaled into its own slot, then the slots are packed. A plain
    whole-graph spring layout pushes the isolated VEZF1 node so far from the
    main component that most of the panel is empty; this keeps every component
    readable at the same node size without any randomness between runs."""
    comps = sorted(nx.connected_components(G), key=len, reverse=True)
    # slots: (centre_x, centre_y, half_width, half_height), largest first
    slots = [(0.36, 0.50, 0.36, 0.46), (0.86, 0.74, 0.11, 0.11), (0.86, 0.30, 0.05, 0.05)]
    pos: dict = {}
    for comp, (cx, cy, hw, hh) in zip(comps, slots):
        sub = G.subgraph(comp)
        if len(comp) == 1:
            pos[next(iter(comp))] = (cx, cy)
            continue
        raw = nx.spring_layout(sub, seed=LAYOUT_SEED, k=0.85, iterations=600)
        xs = [p[0] for p in raw.values()]
        ys = [p[1] for p in raw.values()]
        sx = (max(xs) - min(xs)) or 1.0
        sy = (max(ys) - min(ys)) or 1.0
        for n, (x, y) in raw.items():
            pos[n] = (cx + (2 * (x - min(xs)) / sx - 1) * hw,
                      cy + (2 * (y - min(ys)) / sy - 1) * hh)
    return pos


def analyse(G):
    paths = list(nx.all_shortest_paths(G, "KDM1A", "USP34"))
    return {"n_nodes": G.number_of_nodes(), "n_edges": G.number_of_edges(),
            "n_components": nx.number_connected_components(G),
            "path_length": len(paths[0]) - 1, "n_shortest_paths": len(paths),
            "vezf1_degree": G.degree("VEZF1"), "paths": paths}


def gate(stats):
    checks = [(k, stats[k], v, 0) for k, v in REFERENCE.items()]
    # every shortest route must pass through the same bridge, and its middle
    # node must be a ubiquitin gene -- the claim panel B makes
    all_via_bridge = all(p[1] == BRIDGE for p in stats["paths"])
    all_mid_ubiquitin = all(p[2] in UBIQUITIN for p in stats["paths"])
    checks.append(("every_route_via_single_bridge", 1.0 if all_via_bridge else 0.0, 1.0, 0))
    checks.append(("every_route_middle_is_ubiquitin", 1.0 if all_mid_ubiquitin else 0.0, 1.0, 0))
    checks.append(("n_distinct_middle_nodes", len({p[2] for p in stats["paths"]}), 4, 0))
    return verify(FIGURE, checks)


def build(stub: Path):
    pin_reproducibility(FIGURE)
    G = build_network()
    stats = analyse(G)
    verification = gate(stats)

    fig = plt.figure(figsize=(15.5, 8.6))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.05, 1.0], left=0.02, right=0.985,
                          top=0.885, bottom=0.070, wspace=0.06)
    axl, axr = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])

    # ---- panel A: the whole standardised neighbourhood -----------------------
    pos = component_layout(G)
    named = set(CANDIDATES) | {BRIDGE} | set(UBIQUITIN)
    nx.draw_networkx_edges(G, pos, ax=axl, edge_color=NEUTRAL["grid"], width=1.1)
    other = [n for n in G if n not in named]
    nx.draw_networkx_nodes(G, pos, nodelist=other, ax=axl, node_size=130,
                           node_color=NEUTRAL["backdrop"], edgecolors=WHITE, linewidths=1.0)
    nx.draw_networkx_nodes(G, pos, nodelist=[BRIDGE], ax=axl, node_size=430,
                           node_color=NEUTRAL["ink_2"], edgecolors=WHITE, linewidths=1.6)
    nx.draw_networkx_nodes(G, pos, nodelist=list(UBIQUITIN), ax=axl, node_size=250,
                           node_color=NEUTRAL["ink_muted"], edgecolors=WHITE, linewidths=1.4)
    for gene in CANDIDATES:
        nx.draw_networkx_nodes(G, pos, nodelist=[gene], ax=axl, node_size=760,
                               node_color=GENE_COLOURS[gene], edgecolors=WHITE, linewidths=2.2)
    # fixed label offsets -- no collision solver anywhere
    OFFSET = {"KDM1A": (0.0, 0.085), "TLK2": (0.0, 0.115), "USP34": (0.0, 0.085),
              "VEZF1": (0.0, 0.085), BRIDGE: (0.0, -0.080)}
    for node, (dx, dy) in OFFSET.items():
        x, y = pos[node]
        axl.text(x + dx, y + dy, node, ha="center",
                 va="bottom" if dy > 0 else "top", fontsize=FONT["annot"] + 1, fontweight="bold",
                 color=GENE_COLOURS.get(node, NEUTRAL["ink"]), zorder=20)
    # slot annotations in DATA coordinates, matching the layout slots above
    axl.text(0.86, 0.905, "separate\ncomponent", ha="center", va="bottom", fontsize=FONT["note"] - 1,
             color=NEUTRAL["ink_muted"])
    axl.text(0.86, 0.215, "no partners", ha="center", va="top", fontsize=FONT["note"] - 1,
             color=NEUTRAL["ink_muted"])
    axl.set_title(f"A · One standardised query, applied to all four candidates\n"
                  f"{stats['n_nodes']} proteins, {stats['n_edges']} associations, {stats['n_components']} components",
                  fontsize=FONT["panel"] - 3, color=NEUTRAL["ink"], loc="left", pad=14)
    axl.axis("off")

    # ---- panel B: the four routes, drawn as what they are --------------------
    axr.set_xlim(0, 10)
    axr.set_ylim(0, 10)
    axr.axis("off")
    col_x = [1.1, 3.7, 6.2, 8.6]
    mid_y = [7.3, 6.2, 5.1, 4.0]
    axr.scatter([col_x[0]], [5.8], s=1500, c=GENE_COLOURS["KDM1A"], edgecolors=WHITE,
                linewidths=2.4, zorder=5)
    axr.text(col_x[0], 5.8 - 0.95, "KDM1A", ha="center", va="top", fontsize=FONT["annot"] + 1,
             fontweight="bold", color=GENE_COLOURS["KDM1A"])
    axr.scatter([col_x[3]], [5.8], s=1500, c=GENE_COLOURS["USP34"], edgecolors=WHITE,
                linewidths=2.4, zorder=5)
    axr.text(col_x[3], 5.8 - 0.95, "USP34", ha="center", va="top", fontsize=FONT["annot"] + 1,
             fontweight="bold", color=GENE_COLOURS["USP34"])
    axr.scatter([col_x[1]], [5.8], s=1100, c=NEUTRAL["ink_2"], edgecolors=WHITE,
                linewidths=2.0, zorder=5)
    axr.text(col_x[1], 5.8 - 0.90, BRIDGE, ha="center", va="top", fontsize=FONT["annot"],
             fontweight="bold", color=NEUTRAL["ink"])
    for y, mid in zip(mid_y, [p[2] for p in stats["paths"]]):
        axr.plot([col_x[1], col_x[2]], [5.8, y], color=NEUTRAL["rule"], lw=1.4, zorder=1)
        axr.plot([col_x[2], col_x[3]], [y, 5.8], color=NEUTRAL["rule"], lw=1.4, zorder=1)
        axr.scatter([col_x[2]], [y], s=520, c=NEUTRAL["ink_muted"], edgecolors=WHITE,
                    linewidths=1.6, zorder=5)
        axr.text(col_x[2] + 0.45, y, mid, ha="left", va="center", fontsize=FONT["note"] - 1,
                 color=NEUTRAL["ink_2"])
    axr.plot([col_x[0], col_x[1]], [5.8, 5.8], color=NEUTRAL["rule"], lw=1.4, zorder=1)
    axr.text(5.0, 9.4,
             f"B · All {stats['n_shortest_paths']} routes are the same route",
             ha="center", va="top", fontsize=FONT["panel"] - 3, color=NEUTRAL["ink"], fontweight="bold")
    axr.text(5.0, 8.75,
             f"{stats['path_length']} associations apart. Every route\npasses through {BRIDGE}, then a ubiquitin gene.",
             ha="center", va="top", fontsize=FONT["note"], color=NEUTRAL["ink_2"], linespacing=1.4)
    axr.text(5.0, 2.15,
             "One connection, filled four ways —\nnot four independent lines of support.",
             ha="center", va="top", fontsize=FONT["note"], color=NEUTRAL["ink"], linespacing=1.4)
    axr.text(5.0, 0.85,
             "TLK2 is a separate component;\nVEZF1 has no partner at all.",
             ha="center", va="top", fontsize=FONT["note"] - 1, color=NEUTRAL["ink_muted"], linespacing=1.4)

    headline(
        fig,
        "The two connected candidates are joined by one route, not four",
        "Associations are undirected: they are not activation, not inhibition, not necessarily physical binding, and no arrows are drawn. A short\n"
        "path between two proteins is not a mechanism. Panel B is the reason this connection is weak evidence rather than strong: the four\n"
        "equally short routes differ only in which member of one small gene family occupies the middle position.",
        key=FIGURE)
    figure_footer(fig, "Pinned STRING query (score ≥ 0.70). Layout seeded; no collision solver.")
    return save(fig, stub), verification


def main(out_dir: Path = OUT_DIR):
    return build(out_dir / FIGURE)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    main()
