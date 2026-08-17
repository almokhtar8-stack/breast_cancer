"""ONE poster-grade hero heatmap: cross-context transcriptomic behavior of
the 4 focus genes (KDM1A, TLK2, USP34, VEZF1).

No new analysis, no recomputation of any frozen result. Reuses
`poster_story_v1_data.build_hero_heatmap_pairs()` unchanged -- the same
already-verified log2 fold-change values (see
results/reports/poster_hero_heatmap/NOTE.md for the short data note).
"""

from __future__ import annotations

import logging
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.cm import ScalarMappable
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.patches import FancyBboxPatch, Rectangle

from src import poster_story_v1_data as sv1

logger = logging.getLogger(__name__)

OUT_DIR = Path("results/figures/poster_hero_heatmap")

FOCUS_FOUR = sv1.FOCUS_FOUR
FOCUS_COLORS = sv1.FOCUS_COLORS
DGRAY = "#2a2a2a"
MGRAY = "#8c8c8c"
REF_GRAY = "#e2e0dc"

plt.rcParams["font.family"] = "DejaVu Sans"

# A refined, poster-grade diverging palette -- cool slate blue (decrease)
# through a warm ivory midpoint to a warm terracotta (increase). Deliberately
# less saturated / less "textbook" than default RdBu_r.
DIVERGING_CMAP = LinearSegmentedColormap.from_list(
    "poster_diverging",
    ["#2E6C8E", "#7FAFC4", "#F4EFE6", "#E3A180", "#C1543A"],
    N=256,
)


def build_hero_heatmap(stub: Path) -> None:
    pairs = sv1.build_hero_heatmap_pairs()
    genes = FOCUS_FOUR
    blocks = []
    for dataset in sv1.DATASET_ORDER:
        sub = pairs[pairs["dataset"] == dataset]
        blocks.append(dict(
            dataset=dataset,
            ref_label=sub["ref_label"].iloc[0],
            cmp_label=sub["cmp_label"].iloc[0],
            log2fc=np.array([sub[sub["gene"] == g]["log2fc"].iloc[0] for g in genes], dtype=float),
        ))

    vmax = np.nanmax(np.abs(np.concatenate([b["log2fc"] for b in blocks])))
    norm = TwoSlopeNorm(vcenter=0, vmin=-vmax, vmax=vmax)

    # ---- geometry ----
    row_h = 1.0
    pair_gap = 0.06          # gap between reference and comparison row within a pair
    block_gap = 0.85         # whitespace between dataset blocks
    col_w = 1.0
    n_genes = len(genes)

    fig, ax = plt.subplots(figsize=(9.5, 11.2), dpi=300)

    y_cursor = 0.0
    block_positions = []
    for b in blocks:
        y_ref = y_cursor
        y_cmp = y_cursor - (row_h + pair_gap)
        block_positions.append((b, y_ref, y_cmp))
        y_cursor = y_cmp - (row_h + block_gap)

    top_y = row_h * 0.5 + 0.35
    bottom_y = y_cursor + (row_h + block_gap) - row_h * 0.5 - 0.1

    for b, y_ref, y_cmp in block_positions:
        # reference row -- flat neutral gray, deliberately off the color scale
        for j in range(n_genes):
            ax.add_patch(Rectangle((j * col_w - 0.44, y_ref - row_h / 2), 0.88, row_h,
                                    facecolor=REF_GRAY, edgecolor="white", linewidth=2.2))
        # comparison row -- the real signal
        for j in range(n_genes):
            fc = b["log2fc"][j]
            color = DIVERGING_CMAP(norm(fc)) if np.isfinite(fc) else "#f2f2f2"
            ax.add_patch(Rectangle((j * col_w - 0.44, y_cmp - row_h / 2), 0.88, row_h,
                                    facecolor=color, edgecolor="white", linewidth=2.2))

        # pairing bracket, left of the two rows
        bx0, bx1 = -0.72, -0.58
        ax.plot([bx0, bx0], [y_cmp, y_ref], color=MGRAY, linewidth=1.4, solid_capstyle="round", zorder=1)
        ax.plot([bx0, bx1], [y_ref, y_ref], color=MGRAY, linewidth=1.4, solid_capstyle="round", zorder=1)
        ax.plot([bx0, bx1], [y_cmp, y_cmp], color=MGRAY, linewidth=1.4, solid_capstyle="round", zorder=1)

        # condition labels
        ax.text(-0.82, y_ref, b["ref_label"], ha="right", va="center", fontsize=11, color=MGRAY)
        ax.text(-0.82, y_cmp, b["cmp_label"], ha="right", va="center", fontsize=11.5, color=DGRAY, fontweight="bold")

        # dataset name, left margin
        y_mid = (y_ref + y_cmp) / 2
        ax.text(-2.55, y_mid, b["dataset"], ha="left", va="center", fontsize=15, fontweight="bold", color=DGRAY)

    # gene headers
    header_y = top_y - 0.05
    for j, gene in enumerate(genes):
        ax.text(j * col_w, header_y, gene, ha="center", va="bottom", fontsize=19, fontweight="bold",
                 color=FOCUS_COLORS.get(gene, DGRAY))

    # reference-row swatch, drawn in dedicated whitespace clearly BELOW the
    # last block's cells (not just below the tight bottom_y bound)
    legend_y = bottom_y - 0.4
    ax.add_patch(Rectangle((0.15, legend_y - 0.11), 0.22, 0.22, facecolor=REF_GRAY, edgecolor="none", clip_on=False))
    ax.text(0.48, legend_y, "reference condition", fontsize=10.5, color=MGRAY, va="center", ha="left")

    ax.set_xlim(-2.75, (n_genes - 1) * col_w + 0.6)
    ax.set_ylim(legend_y - 0.35, top_y + 0.55)
    ax.axis("off")

    # ---- colorbar (the only "legend", replaces per-cell numbers) --
    # placed relative to `ax` itself so matplotlib shrinks the main axes to
    # make room automatically -- never overlaps the heatmap cells.
    sm = ScalarMappable(norm=norm, cmap=DIVERGING_CMAP)
    sm.set_array([])
    cb = fig.colorbar(sm, ax=ax, fraction=0.045, pad=0.1, shrink=0.32, anchor=(0.0, 0.02))
    cb.set_label("log$_2$ fold-change\nvs. reference", fontsize=10.5, color=DGRAY)
    cb.ax.tick_params(labelsize=9, color=DGRAY, labelcolor=DGRAY)
    cb.outline.set_visible(False)

    fig.text(0.01, 1.045, "Cross-context transcriptomic behavior of the four focus candidates",
              fontsize=23, fontweight="bold", color=DGRAY, ha="left", va="bottom")
    fig.text(0.01, -0.03,
              "Gray = reference condition. Color = log$_2$ fold-change vs. that reference -- resistance-model, "
              "human-recurrence, and acute ex vivo contexts.",
              fontsize=11.5, color="#555555", ha="left", va="top")

    stub.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(stub.with_suffix(".png"), facecolor="white", bbox_inches="tight", dpi=300)
    fig.savefig(stub.with_suffix(".pdf"), facecolor="white", bbox_inches="tight")
    fig.savefig(stub.with_suffix(".svg"), facecolor="white", bbox_inches="tight")
    plt.close(fig)
    logger.info("wrote %s.png/.pdf/.svg", stub)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    build_hero_heatmap(OUT_DIR / "HERO_main_heatmap")
