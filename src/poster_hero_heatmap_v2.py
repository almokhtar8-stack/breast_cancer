"""ONE poster-grade hero heatmap, v2: fixes the v1 blank-looking reference
rows by giving them a real (within-dataset-only) sequential encoding, and
tightens the layout so the heatmap dominates the canvas.

No new analysis, no recomputation of any frozen result. Reuses
`poster_story_v1_data.build_hero_heatmap_pairs()` unchanged -- the exact
same already-verified log2 fold-change AND reference-condition mean values
already used in v1 (see results/reports/poster_hero_heatmap_v2/NOTE.md).

Two-layer honest encoding per dataset block:
  - reference row: SEQUENTIAL palette, min-max scaled WITHIN that dataset
    block's 4 genes only (never claimed comparable across blocks)
  - comparison row: DIVERGING palette, log2 fold-change vs. that same
    dataset's reference (comparable across all 4 blocks -- log2FC is a
    ratio, not a raw scale)
"""

from __future__ import annotations

import logging
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.cm import ScalarMappable
from matplotlib.colors import LinearSegmentedColormap, Normalize, TwoSlopeNorm
from matplotlib.patches import Rectangle

from src import poster_story_v1_data as sv1

logger = logging.getLogger(__name__)

OUT_DIR = Path("results/figures/poster_hero_heatmap_v2")

FOCUS_FOUR = sv1.FOCUS_FOUR
FOCUS_COLORS = sv1.FOCUS_COLORS
DGRAY = "#262626"
MGRAY = "#8c8c8c"

# Sequential palette for the reference row -- warm neutral (mist -> deep
# taupe), deliberately a different hue family from the diverging palette
# below so a viewer never confuses "high baseline" with "positive fold
# change".
SEQUENTIAL_CMAP = LinearSegmentedColormap.from_list(
    "poster_sequential", ["#FAF7F1", "#E7DECD", "#C7B79C", "#A08F6C", "#6E5F44"], N=256,
)

# Diverging palette for the comparison row -- slate blue (decrease) through
# warm stone (near zero) to terracotta (increase).
DIVERGING_CMAP = LinearSegmentedColormap.from_list(
    "poster_diverging", ["#2E6C8E", "#7FAFC4", "#F3EEE4", "#E3A180", "#C1543A"], N=256,
)


def build_hero_heatmap_v2(stub: Path) -> None:
    pairs = sv1.build_hero_heatmap_pairs()
    genes = FOCUS_FOUR
    blocks = []
    for dataset in sv1.DATASET_ORDER:
        sub = pairs[pairs["dataset"] == dataset]
        ref_vals = np.array([sub[sub["gene"] == g]["ref_value"].iloc[0] for g in genes], dtype=float)
        blocks.append(dict(
            dataset=dataset,
            ref_label=sub["ref_label"].iloc[0],
            cmp_label=sub["cmp_label"].iloc[0],
            ref_vals=ref_vals,
            log2fc=np.array([sub[sub["gene"] == g]["log2fc"].iloc[0] for g in genes], dtype=float),
        ))

    vmax = np.nanmax(np.abs(np.concatenate([b["log2fc"] for b in blocks])))
    div_norm = TwoSlopeNorm(vcenter=0, vmin=-vmax, vmax=vmax)

    # ---- tight geometry: the heatmap body should dominate the canvas ----
    row_h = 1.0
    pair_gap = 0.03
    block_gap = 0.42
    col_w = 1.0
    n_genes = len(genes)

    fig, ax = plt.subplots(figsize=(8.6, 8.2), dpi=300)

    y_cursor = 0.0
    block_positions = []
    for b in blocks:
        y_ref = y_cursor
        y_cmp = y_cursor - (row_h + pair_gap)
        block_positions.append((b, y_ref, y_cmp))
        y_cursor = y_cmp - (row_h + block_gap)

    top_y = row_h * 0.5 + 0.42
    bottom_y = y_cursor + (row_h + block_gap) - row_h * 0.5

    for b, y_ref, y_cmp in block_positions:
        # within-block sequential scaling for the reference row -- real
        # data-bearing color, explicitly NOT cross-dataset comparable
        rv = b["ref_vals"]
        rmin, rmax = np.nanmin(rv), np.nanmax(rv)
        seq_norm = Normalize(vmin=rmin, vmax=rmax) if rmax > rmin else Normalize(vmin=rmin - 1, vmax=rmax + 1)

        for j in range(n_genes):
            ref_color = SEQUENTIAL_CMAP(seq_norm(rv[j]))
            ax.add_patch(Rectangle((j * col_w - 0.46, y_ref - row_h / 2), 0.92, row_h,
                                    facecolor=ref_color, edgecolor="white", linewidth=2.4))
        for j in range(n_genes):
            fc = b["log2fc"][j]
            cmp_color = DIVERGING_CMAP(div_norm(fc)) if np.isfinite(fc) else "#f0f0f0"
            ax.add_patch(Rectangle((j * col_w - 0.46, y_cmp - row_h / 2), 0.92, row_h,
                                    facecolor=cmp_color, edgecolor="white", linewidth=2.4))

        # pairing bracket, left of the two rows
        bx0, bx1 = -0.70, -0.58
        ax.plot([bx0, bx0], [y_cmp, y_ref], color=MGRAY, linewidth=1.5, solid_capstyle="round", zorder=3)
        ax.plot([bx0, bx1], [y_ref, y_ref], color=MGRAY, linewidth=1.5, solid_capstyle="round", zorder=3)
        ax.plot([bx0, bx1], [y_cmp, y_cmp], color=MGRAY, linewidth=1.5, solid_capstyle="round", zorder=3)

        ax.text(-0.80, y_ref, b["ref_label"], ha="right", va="center", fontsize=10.5, color=MGRAY)
        ax.text(-0.80, y_cmp, b["cmp_label"], ha="right", va="center", fontsize=11.5, color=DGRAY, fontweight="bold")

        y_mid = (y_ref + y_cmp) / 2
        ax.text(-2.55, y_mid, b["dataset"], ha="left", va="center", fontsize=16, fontweight="bold", color=DGRAY)

    # gene headers
    header_y = top_y - 0.08
    for j, gene in enumerate(genes):
        ax.text(j * col_w, header_y, gene, ha="center", va="bottom", fontsize=20, fontweight="bold",
                 color=FOCUS_COLORS.get(gene, DGRAY))

    ax.set_xlim(-2.75, (n_genes - 1) * col_w + 0.62)
    ax.set_ylim(bottom_y - 0.15, top_y + 0.5)
    ax.axis("off")

    # ---- two compact legends, side by side, below the heatmap body ----
    seq_sm = ScalarMappable(norm=Normalize(vmin=0, vmax=1), cmap=SEQUENTIAL_CMAP)
    seq_sm.set_array([])
    div_sm = ScalarMappable(norm=div_norm, cmap=DIVERGING_CMAP)
    div_sm.set_array([])

    cax1 = fig.add_axes([0.17, -0.045, 0.28, 0.028])
    cb1 = fig.colorbar(seq_sm, cax=cax1, orientation="horizontal")
    cb1.set_ticks([0, 1])
    cb1.set_ticklabels(["low", "high"])
    cb1.ax.tick_params(labelsize=8.5, length=0, color=MGRAY, labelcolor=MGRAY)
    cb1.outline.set_visible(False)
    fig.text(0.17, -0.065, "reference row: relative baseline (within dataset only)", fontsize=8.7, color=MGRAY, va="top")

    cax2 = fig.add_axes([0.56, -0.045, 0.28, 0.028])
    cb2 = fig.colorbar(div_sm, cax=cax2, orientation="horizontal")
    cb2.set_ticks([-vmax, 0, vmax])
    cb2.set_ticklabels([f"-{vmax:.1f}", "0", f"+{vmax:.1f}"])
    cb2.ax.tick_params(labelsize=8.5, length=0, color=DGRAY, labelcolor=DGRAY)
    cb2.outline.set_visible(False)
    fig.text(0.56, -0.065, "comparison row: log$_2$ fold-change vs. reference (cross-context comparable)",
              fontsize=8.7, color=DGRAY, va="top")

    fig.text(0.0, 1.05, "Cross-context transcriptomic behavior of the four focus candidates",
              fontsize=23, fontweight="bold", color=DGRAY, ha="left", va="bottom")
    fig.text(0.0, 1.005,
              "Reference rows show baseline expression within each dataset; comparison rows show log$_2$ "
              "fold-change vs. that reference.",
              fontsize=11.5, color="#555555", ha="left", va="bottom")

    stub.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(stub.with_suffix(".png"), facecolor="white", bbox_inches="tight", dpi=300)
    fig.savefig(stub.with_suffix(".pdf"), facecolor="white", bbox_inches="tight")
    fig.savefig(stub.with_suffix(".svg"), facecolor="white", bbox_inches="tight")
    plt.close(fig)
    logger.info("wrote %s.png/.pdf/.svg", stub)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    build_hero_heatmap_v2(OUT_DIR / "HERO_main_heatmap_v2")
