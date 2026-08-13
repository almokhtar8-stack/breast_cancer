"""Ad-hoc figure: USP34 shortest undirected paths to CTNNB1, PTEN, EP300, SOX2,
using ONLY the edges that are actually part of one of those shortest paths
(results/tables/systems_network/USP34_shortest_paths.tsv). Other edges exist
between these same genes in the full network (e.g. EP300-CTNNB1 is itself a
direct physical_PPI edge) but are out of scope for this figure and are
deliberately omitted so the picture answers exactly the question asked --
"how does USP34 reach each target" -- without implying those incidental edges
are part of a USP34 shortest path.

Reads only existing systems-network output tables; does not alter or rerun
any upstream phase.
"""

from __future__ import annotations

import logging
import math
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch

logger = logging.getLogger(__name__)

PATHS_TABLE = Path("results/tables/systems_network/USP34_shortest_paths.tsv")
OUT_FIG = Path("results/figures/systems_network/final_review/07_USP34_shortest_paths.png")

# Okabe-Ito colorblind-safe categorical set
COLOR_SOURCE = "#0072B2"
COLOR_INTERMEDIATE = "#E69F00"
COLOR_TARGET = "#009E73"
COLOR_SHARED_INTERMEDIATE_RING = "#D55E00"

EDGE_STYLE = {
    "physical_PPI": {"color": "#222222", "linestyle": "-", "linewidth": 2.4},
    "functional_association": {"color": "#888888", "linestyle": "--", "linewidth": 1.6},
}

# Manual layout: USP34 at the center, intermediates on an inner ring,
# targets on an outer ring, angularly placed near the intermediate(s) they
# connect to -- deterministic and chosen only for label/edge readability, not
# derived from any centrality or force-directed algorithm.
INTERMEDIATE_ANGLES_DEG = {
    "RPS27A": 150,
    "USP9X": 45,
    "UBB": 210,
    "UBC": 315,
}
TARGET_ANGLES_DEG = {
    "CTNNB1": 90,
    "SOX2": 15,
    "EP300": 340,
    "PTEN": 255,
}
R_INTERMEDIATE = 1.6
R_TARGET = 3.15


def _xy(angle_deg: float, r: float) -> tuple[float, float]:
    a = math.radians(angle_deg)
    return r * math.cos(a), r * math.sin(a)


def build_figure(paths_df: pd.DataFrame, out_fig: Path = OUT_FIG) -> None:
    edge_cols = ["source_gene", "target_gene", "interaction_type", "database_source", "confidence"]
    edges = paths_df.loc[paths_df["path"] != "NO_PATH_IN_NETWORK", edge_cols].drop_duplicates()

    intermediates = sorted(set(INTERMEDIATE_ANGLES_DEG) & (set(edges["source_gene"]) | set(edges["target_gene"])))
    targets = sorted(set(TARGET_ANGLES_DEG) & (set(edges["source_gene"]) | set(edges["target_gene"])))

    pos = {"USP34": (0.0, 0.0)}
    for g in intermediates:
        pos[g] = _xy(INTERMEDIATE_ANGLES_DEG[g], R_INTERMEDIATE)
    for g in targets:
        pos[g] = _xy(TARGET_ANGLES_DEG[g], R_TARGET)

    # shared intermediates: connect to >=2 of the 4 targets
    dest_per_intermediate: dict[str, set[str]] = {g: set() for g in intermediates}
    for _, row in edges.iterrows():
        u, v = row["source_gene"], row["target_gene"]
        for a, b in [(u, v), (v, u)]:
            if a in dest_per_intermediate and b in TARGET_ANGLES_DEG:
                dest_per_intermediate[a].add(b)
    shared_intermediates = {g for g, dests in dest_per_intermediate.items() if len(dests) >= 2}

    fig, ax = plt.subplots(figsize=(9, 9), dpi=200)
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")

    for _, row in edges.iterrows():
        u, v = row["source_gene"], row["target_gene"]
        style = EDGE_STYLE[row["interaction_type"]]
        x1, y1 = pos[u]
        x2, y2 = pos[v]
        ax.add_patch(
            FancyArrowPatch(
                (x1, y1),
                (x2, y2),
                arrowstyle="-",
                color=style["color"],
                linestyle=style["linestyle"],
                linewidth=style["linewidth"],
                zorder=1,
                shrinkA=18,
                shrinkB=18,
            )
        )
        xm, ym = (x1 + x2) / 2, (y1 + y2) / 2
        ax.text(xm, ym, f"{row['confidence']:.2f}", fontsize=7, color=style["color"], ha="center", va="center", backgroundcolor="white", zorder=2)

    # Labels are placed OUTSIDE each node marker (not centered inside it):
    # gene symbols vary in length (USP34 vs CTNNB1 vs UBB) and a fixed-size
    # circle cannot fit all of them without either shrinking font to
    # illegibility or letting the outer characters spill onto the white
    # background, where white-on-white text silently disappears.
    node_specs = [("USP34", COLOR_SOURCE, 950, 11, (0, 0))]
    node_specs += [(g, COLOR_INTERMEDIATE, 650, 10, None) for g in intermediates]
    node_specs += [(g, COLOR_TARGET, 750, 10.5, None) for g in targets]

    for gene, color, size, fontsize, label_offset in node_specs:
        x, y = pos[gene]
        edgecolor = COLOR_SHARED_INTERMEDIATE_RING if gene in shared_intermediates else "white"
        lw = 3.5 if gene in shared_intermediates else 1.5
        ax.scatter([x], [y], s=size, c=color, edgecolors=edgecolor, linewidths=lw, zorder=3)
        if label_offset is None:
            # push the label radially outward from the origin so it clears
            # the marker regardless of which ring the node sits on
            r = math.hypot(x, y)
            dx, dy = (x / r, y / r) if r > 1e-6 else (0, 1)
            lx, ly = x + dx * 0.42, y + dy * 0.42
            ha = "left" if dx > 0.15 else ("right" if dx < -0.15 else "center")
        else:
            lx, ly = x, y - 0.42
            ha = "center"
        ax.text(
            lx,
            ly,
            gene,
            fontsize=fontsize,
            fontweight="bold",
            ha=ha,
            va="center",
            color="#111111",
            zorder=5,
            bbox=dict(boxstyle="round,pad=0.15", facecolor="white", edgecolor="none", alpha=0.85),
        )

    legend_elements = [
        Line2D([0], [0], color=EDGE_STYLE["physical_PPI"]["color"], lw=2.4, linestyle="-", label="physical_PPI edge (STRING)"),
        Line2D([0], [0], color=EDGE_STYLE["functional_association"]["color"], lw=1.6, linestyle="--", label="functional_association edge (STRING)"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor=COLOR_SOURCE, markersize=14, label="USP34 (source)"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor=COLOR_INTERMEDIATE, markersize=12, label="intermediate gene"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor=COLOR_TARGET, markersize=13, label="target gene"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor=COLOR_INTERMEDIATE, markeredgecolor=COLOR_SHARED_INTERMEDIATE_RING, markeredgewidth=3, markersize=12, label="intermediate shared by >=2 targets"),
    ]
    ax.legend(handles=legend_elements, loc="upper left", bbox_to_anchor=(-0.15, 1.05), fontsize=9, frameon=False)

    ax.set_title("USP34 shortest undirected paths (all length = 2 edges)\nlabels on edges = STRING combined_score", fontsize=12, fontweight="bold")
    ax.set_xlim(-4.2, 4.2)
    ax.set_ylim(-4.2, 4.2)
    ax.set_aspect("equal")
    ax.axis("off")

    fig.tight_layout()
    out_fig.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_fig, facecolor="white")
    plt.close(fig)
    logger.info("wrote %s", out_fig)


def run(paths_table: Path = PATHS_TABLE, out_fig: Path = OUT_FIG) -> None:
    df = pd.read_csv(paths_table, sep="\t")
    build_figure(df, out_fig)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run()
