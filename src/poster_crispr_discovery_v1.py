"""ONE poster-grade CRISPR discovery figure -- the entry figure for the
poster story, answering: which genes sensitize ER+ breast-cancer cells to
tamoxifen, and where do the four focus candidates sit within that
discovery result?

Data source: `src.post_audit_sensitivity_data.load_significant_sensitising_hits()`
and `.load_genomewide_crispr()`, called unmodified. Both are already-frozen
loaders used elsewhere in this project (the post-audit sensitivity
analysis) -- this module performs NO new discovery, NO re-ranking, and NO
recomputation of the CRISPR effect size, FDR, or the pre-specified
significance gate (FDR < 0.10, PREANALYSIS.md Section 4). Every number
that appears on the figure (hit count, genome-wide gene count, effect
sizes, ranks) is read from these frozen tables at render time -- nothing
is hand-typed.
"""

from __future__ import annotations

import logging
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from src import post_audit_sensitivity_data as pad
from src import poster_story_v1_data as sv1

logger = logging.getLogger(__name__)

OUT_DIR = Path("results/figures/poster_crispr_discovery_v1")

FOCUS_FOUR = sv1.FOCUS_FOUR  # ["KDM1A", "TLK2", "USP34", "VEZF1"]
FOCUS_COLORS = sv1.FOCUS_COLORS

DGRAY = "#262626"
MGRAY = "#8c8c8c"
MUTED = "#b3b3b3"
MUTED_LABEL = "#6e6e6e"


def build_crispr_discovery_main(stub: Path) -> None:
    genomewide = pad.load_genomewide_crispr()
    hits = pad.load_significant_sensitising_hits()  # already sorted by rank_by_effect ascending
    n_genomewide = len(genomewide)
    n_hits = len(hits)

    genes = hits["gene"].tolist()
    effects = hits["effect_size"].to_numpy()
    y_pos = np.arange(n_hits)  # 0 = strongest (rank 1), inverted below so it plots at the top
    is_focus = [g in FOCUS_FOUR for g in genes]

    fig, ax = plt.subplots(figsize=(12.5, 7.6), dpi=300)

    for i, gene in enumerate(genes):
        focus = is_focus[i]
        color = FOCUS_COLORS[gene] if focus else MUTED
        lw = 3.2 if focus else 1.6
        ms = 190 if focus else 90
        zorder = 5 if focus else 3
        ax.hlines(y=y_pos[i], xmin=0, xmax=effects[i], color=color, linewidth=lw, zorder=zorder)
        ax.scatter([effects[i]], [y_pos[i]], s=ms, color=color, zorder=zorder + 1,
                   edgecolor="white", linewidth=1.0)

    ax.axvline(0, color="#d0d0d0", linewidth=1.1, zorder=1)

    ax.set_yticks(y_pos)
    ax.set_yticklabels(genes, fontsize=12.5)
    for i, tick_label in enumerate(ax.get_yticklabels()):
        if is_focus[i]:
            tick_label.set_color(FOCUS_COLORS[genes[i]])
            tick_label.set_fontweight("bold")
            tick_label.set_fontsize(15)
        else:
            tick_label.set_color(MUTED_LABEL)
            tick_label.set_fontweight("normal")
    ax.invert_yaxis()

    ax.set_xlabel("CRISPR effect size  (more negative = stronger sensitising knockout)",
                  fontsize=11.5, color=DGRAY, labelpad=10)
    ax.set_ylabel("")
    xmin = effects.min() * 1.12
    ax.set_xlim(xmin, 0.32)
    ax.set_ylim(n_hits - 0.4, -0.9)

    for spine in ("top", "right", "left"):
        ax.spines[spine].set_visible(False)
    ax.spines["bottom"].set_color("#c9c9c9")
    ax.tick_params(axis="y", length=0)
    ax.tick_params(axis="x", labelsize=10, colors=DGRAY)
    ax.grid(axis="x", color="#eeeeee", linewidth=0.8, zorder=0)
    ax.set_axisbelow(True)

    fig.text(0.06, 0.975, "Genome-scale CRISPR screen identifies the tamoxifen-sensitising hits",
              fontsize=21, fontweight="bold", color=DGRAY, ha="left", va="top")
    fig.text(0.06, 0.925,
              f"{n_hits} genes met the pre-specified significance threshold (FDR < 0.10) among "
              f"{n_genomewide:,} genome-wide fitted genes.",
              fontsize=12, color="#555555", ha="left", va="top")
    fig.text(0.06, 0.895,
              "Bold, colored labels: the four candidates carried forward into subsequent transcriptomic analysis. "
              "Gray: the other significant sensitising hits.",
              fontsize=10, color=MGRAY, ha="left", va="top", style="italic")

    fig.subplots_adjust(left=0.14, right=0.97, top=0.85, bottom=0.11)

    stub.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(stub.with_suffix(".png"), facecolor="white", bbox_inches="tight", dpi=300)
    fig.savefig(stub.with_suffix(".pdf"), facecolor="white", bbox_inches="tight")
    fig.savefig(stub.with_suffix(".svg"), facecolor="white", bbox_inches="tight")
    plt.close(fig)
    logger.info("wrote %s.png/.pdf/.svg (%d significant sensitising hits of %d fitted genes)",
                stub, n_hits, n_genomewide)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    build_crispr_discovery_main(OUT_DIR / "CRISPR_discovery_main")
