"""Literature/mechanism review figure: candidate -> experimentally
supported regulator/pathway -> proposed resistance relevance.

Deliberately small (11 nodes, 3 tiers) -- reads only the already-built
claim-evidence table for the exact PMIDs shown in edge labels; no new
literature search or project-data computation happens here.
"""

from __future__ import annotations

import logging
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

logger = logging.getLogger(__name__)

OUT_FIG = Path("results/figures/literature_mechanism/01_candidate_mechanism_evidence_map.png")

EDGE_STYLE = {
    "OUR_DATA": dict(color="#0072B2", linestyle="-", linewidth=2.6),
    "LITERATURE_DIRECT": dict(color="#111111", linestyle="-", linewidth=2.0),
    "LITERATURE_INDIRECT": dict(color="#555555", linestyle="--", linewidth=1.8),
    "HYPOTHESIS_INFERENCE": dict(color="#b0392f", linestyle=":", linewidth=2.0),
    "LITERATURE_CONTRADICTORY": dict(color="#b0392f", linestyle="-", linewidth=2.0),
}

NODE_KIND_COLOR = {
    "candidate": "#0072B2",
    "candidate_isolated": "#9a9a9a",
    "regulator": "#E69F00",
    "outcome": "#009E73",
    "outcome_contested": "#b0392f",
    "outcome_unresolved": "#9a9a9a",
}

# x,y positions: col1 (candidates) x~-8, col2 (regulator/pathway) x~0, col3 (outcome) x~8
NODES = {
    "USP34": (-8, 6, "candidate"),
    "AXIN1 / Wnt-beta-catenin\ntranscriptional output": (0, 6, "regulator"),
    "Possible endocrine-resistance\nphenotype (untested in breast)": (8, 6.6, "outcome"),

    "VEZF1": (-8, 1.5, "candidate"),
    "Angiogenic gene program\n(VEGFR2, endothelin-1,\nrepression of CITED2)": (0, 1.5, "regulator"),
    "Possible vascular/tumor-niche\nadaptation (unestablished\nin drug resistance)": (8, 1.9, "outcome"),

    "CITED2": (-8, -3, "candidate"),
    "ER transcriptional activity /\np53 suppression": (0, -3, "regulator"),
    "Endocrine resistance\n(direction CONTESTED --\nclinical data opposite)": (8, -2.6, "outcome_contested"),

    "EML5": (-8, -7, "candidate_isolated"),
    "LITERATURE-MECHANISM\nUNRESOLVED": (0, -7, "outcome_unresolved"),
}

EDGES = [
    ("USP34", "AXIN1 / Wnt-beta-catenin\ntranscriptional output", "LITERATURE_DIRECT", "PMID 21383061\n(HEK293T/colorectal,\nnot breast)", 0.5, (0, 0.55)),
    ("AXIN1 / Wnt-beta-catenin\ntranscriptional output", "Possible endocrine-resistance\nphenotype (untested in breast)", "HYPOTHESIS_INFERENCE", "inference only --\nnever tested in breast", 0.5, (0, 0.55)),
    ("USP34", "Possible endocrine-resistance\nphenotype (untested in breast)", "OUR_DATA", "own CRISPR FDR=0.042\n+ WNT leading-edge", 0.5, (0, 2.8)),

    ("VEZF1", "Angiogenic gene program\n(VEGFR2, endothelin-1,\nrepression of CITED2)", "LITERATURE_DIRECT", "PMID 15031128,\n11504723, 29794136\n(non-cancer systems)", 0.5, (0, 0.6)),
    ("Angiogenic gene program\n(VEGFR2, endothelin-1,\nrepression of CITED2)", "Possible vascular/tumor-niche\nadaptation (unestablished\nin drug resistance)", "HYPOTHESIS_INFERENCE", "no paper links this\nto tumor resistance", 0.5, (0, 0.6)),
    ("VEZF1", "Possible vascular/tumor-niche\nadaptation (unestablished\nin drug resistance)", "OUR_DATA", "own CRISPR+RNA +\n2 STRONG_CONSENSUS pathways", 0.5, (0, 2.8)),

    ("CITED2", "ER transcriptional activity /\np53 suppression", "LITERATURE_DIRECT", "PMID 23811274,\n27627783 (breast, ER+)", 0.5, (0, 0.6)),
    ("ER transcriptional activity /\np53 suppression", "Endocrine resistance\n(direction CONTESTED --\nclinical data opposite)", "LITERATURE_CONTRADICTORY", "PMID 19904269:\nresistance-selection vs\nclinical outcome disagree", 0.5, (0, 0.6)),
    ("CITED2", "Endocrine resistance\n(direction CONTESTED --\nclinical data opposite)", "OUR_DATA", "own RNA resistance\nsupport (1/3 datasets)", 0.5, (0, 2.8)),

    ("VEZF1", "CITED2", "LITERATURE_DIRECT", "represses\n(PMID 29794136)", 0.5, (-2.6, 0)),
]

# OUR_DATA edges bow well below the row so they never cross the col2 box
# sitting at the same y as the candidate node.
ARC_RAD = {
    ("USP34", "Possible endocrine-resistance\nphenotype (untested in breast)"): -0.35,
    ("VEZF1", "Possible vascular/tumor-niche\nadaptation (unestablished\nin drug resistance)"): -0.35,
    ("CITED2", "Endocrine resistance\n(direction CONTESTED --\nclinical data opposite)"): -0.35,
    ("VEZF1", "CITED2"): 0.15,
}


def build_figure(out_fig: Path = OUT_FIG) -> None:
    fig, ax = plt.subplots(figsize=(14, 11), dpi=200)
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")

    for u, v, etype, label, t, (odx, ody) in EDGES:
        x1, y1, _ = NODES[u]
        x2, y2, _ = NODES[v]
        style = EDGE_STYLE[etype]
        rad = ARC_RAD.get((u, v), 0.0)
        ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>", mutation_scale=16, color=style["color"], linestyle=style["linestyle"], linewidth=style["linewidth"], zorder=1, shrinkA=45, shrinkB=45, connectionstyle=f"arc3,rad={rad}"))
        xm, ym = x1 + (x2 - x1) * t + odx, y1 + (y2 - y1) * t + ody
        ax.text(xm, ym, label, fontsize=7.3, color=style["color"], ha="center", va="center", backgroundcolor="white", zorder=2, style="italic")

    for name, (x, y, kind) in NODES.items():
        color = NODE_KIND_COLOR[kind]
        is_box = kind in ("regulator", "outcome", "outcome_contested", "outcome_unresolved")
        if is_box:
            w, h = 3.6, 1.5
            edgecolor = "#b0392f" if kind == "outcome_contested" else ("#666666" if kind == "outcome_unresolved" else color)
            facecolor = "#fbeceb" if kind == "outcome_contested" else ("#f0f0ef" if kind == "outcome_unresolved" else "#fff6e8" if kind == "regulator" else "#e7f7f2")
            ax.add_patch(FancyBboxPatch((x - w / 2, y - h / 2), w, h, boxstyle="round,pad=0.08", linewidth=1.8, edgecolor=edgecolor, facecolor=facecolor, zorder=3))
            ax.text(x, y, name, fontsize=8.3, ha="center", va="center", color="#111111", fontweight="bold", zorder=4)
        else:
            linestyle = "--" if kind == "candidate_isolated" else "-"
            ax.scatter([x], [y], s=2400, c=color, edgecolors="white", linewidths=2, linestyle=linestyle, zorder=3)
            ax.text(x, y, name, fontsize=13, fontweight="bold", ha="center", va="center", color="white", zorder=4)

    legend_elements = [
        Line2D([0], [0], color=EDGE_STYLE["OUR_DATA"]["color"], lw=2.6, linestyle="-", label="our project evidence (CRISPR/RNA)"),
        Line2D([0], [0], color=EDGE_STYLE["LITERATURE_DIRECT"]["color"], lw=2.0, linestyle="-", label="direct published experimental evidence"),
        Line2D([0], [0], color=EDGE_STYLE["LITERATURE_INDIRECT"]["color"], lw=1.8, linestyle="--", label="indirect literature-supported mechanism"),
        Line2D([0], [0], color=EDGE_STYLE["HYPOTHESIS_INFERENCE"]["color"], lw=2.0, linestyle=":", label="hypothesis / inference (untested link)"),
        Line2D([0], [0], color=EDGE_STYLE["LITERATURE_CONTRADICTORY"]["color"], lw=2.0, linestyle="-", label="literature evidence, but internally contradictory"),
    ]
    fig.legend(handles=legend_elements, loc="lower center", bbox_to_anchor=(0.5, 0.0), fontsize=9, frameon=False, ncol=2)

    fig.suptitle("Candidate -> experimentally supported regulator/pathway -> proposed resistance relevance\n(EML5: no literature or network mechanism found; not forced)", fontsize=12.5, fontweight="bold", x=0.02, ha="left", y=0.995)

    ax.set_xlim(-11, 11.5)
    ax.set_ylim(-9.5, 8.5)
    ax.set_aspect("equal")
    ax.axis("off")

    fig.subplots_adjust(left=0.02, right=0.98, top=0.92, bottom=0.16)
    out_fig.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_fig, facecolor="white")
    plt.close(fig)
    logger.info("wrote %s", out_fig)


def run(out_fig: Path = OUT_FIG) -> None:
    build_figure(out_fig)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run()
