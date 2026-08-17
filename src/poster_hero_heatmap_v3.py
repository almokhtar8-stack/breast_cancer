"""ONE poster-grade hero heatmap, v3: a TRUE paired-condition heatmap.

Fixes v2's mixed-semantics problem (reference row and comparison row used
DIFFERENT color meanings, even though they looked like one heatmap). v3
gives each dataset block three rows:
  1. reference condition   \  SAME sequential expression scale,
  2. comparison condition  /  normalized together, within this block only
  3. Delta (log2FC) strip  -- ONE shared diverging scale across all blocks

No new analysis, no recomputation of any frozen result. Reuses
`poster_story_v1_data.build_hero_heatmap_pairs()` unchanged -- the exact
same already-verified reference/comparison means and log2 fold-change
values already used in v1 and v2 (see
results/reports/poster_hero_heatmap_v3/NOTE.md).
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

OUT_DIR = Path("results/figures/poster_hero_heatmap_v3")

FOCUS_FOUR = sv1.FOCUS_FOUR
FOCUS_COLORS = sv1.FOCUS_COLORS
DGRAY = "#262626"
MGRAY = "#8c8c8c"

# Sequential palette for expression rows (reference + comparison, SAME
# scale within a dataset block) -- warm neutral mist -> deep taupe.
SEQUENTIAL_CMAP = LinearSegmentedColormap.from_list(
    "poster_sequential", ["#FAF7F1", "#E7DECD", "#C7B79C", "#A08F6C", "#6E5F44"], N=256,
)

# Diverging palette for the Delta strip -- ONE shared scale across all 4
# dataset blocks. Slate blue (decrease) -> warm stone (~0) -> terracotta
# (increase).
DIVERGING_CMAP = LinearSegmentedColormap.from_list(
    "poster_diverging", ["#2E6C8E", "#7FAFC4", "#F3EEE4", "#E3A180", "#C1543A"], N=256,
)


def build_hero_heatmap_v3(stub: Path) -> None:
    pairs = sv1.build_hero_heatmap_pairs()
    genes = FOCUS_FOUR
    blocks = []
    for dataset in sv1.DATASET_ORDER:
        sub = pairs[pairs["dataset"] == dataset]
        ref_vals = np.array([sub[sub["gene"] == g]["ref_value"].iloc[0] for g in genes], dtype=float)
        cmp_vals = np.array([sub[sub["gene"] == g]["cmp_value"].iloc[0] for g in genes], dtype=float)
        blocks.append(dict(
            dataset=dataset,
            ref_label=sub["ref_label"].iloc[0],
            cmp_label=sub["cmp_label"].iloc[0],
            ref_vals=ref_vals,
            cmp_vals=cmp_vals,
            log2fc=np.array([sub[sub["gene"] == g]["log2fc"].iloc[0] for g in genes], dtype=float),
        ))

    vmax = np.nanmax(np.abs(np.concatenate([b["log2fc"] for b in blocks])))
    div_norm = TwoSlopeNorm(vcenter=0, vmin=-vmax, vmax=vmax)

    # ---- geometry ----
    row_h = 1.0
    pair_gap = 0.03
    delta_h = 0.42
    delta_top_gap = 0.22
    block_gap = 0.55
    col_w = 1.0
    n_genes = len(genes)

    fig, ax = plt.subplots(figsize=(8.8, 9.6), dpi=300)

    y_cursor = 0.0
    block_positions = []
    for b in blocks:
        y_ref = y_cursor
        y_cmp = y_ref - (row_h + pair_gap)
        y_delta = y_cmp - row_h / 2 - delta_top_gap - delta_h / 2
        block_positions.append((b, y_ref, y_cmp, y_delta))
        y_cursor = y_delta - delta_h / 2 - block_gap

    top_y = row_h * 0.5 + 0.42
    bottom_y = y_cursor + (delta_h / 2 + block_gap)

    for b, y_ref, y_cmp, y_delta in block_positions:
        # SAME sequential scale for reference + comparison rows, computed
        # together across both rows' 8 values within this dataset block
        # only (never claimed comparable to any other block).
        block_vals = np.concatenate([b["ref_vals"], b["cmp_vals"]])
        vmin_b, vmax_b = np.nanmin(block_vals), np.nanmax(block_vals)
        seq_norm = Normalize(vmin=vmin_b, vmax=vmax_b) if vmax_b > vmin_b else Normalize(vmin=vmin_b - 1, vmax=vmax_b + 1)

        for j in range(n_genes):
            ref_color = SEQUENTIAL_CMAP(seq_norm(b["ref_vals"][j]))
            ax.add_patch(Rectangle((j * col_w - 0.46, y_ref - row_h / 2), 0.92, row_h,
                                    facecolor=ref_color, edgecolor="white", linewidth=2.2, zorder=2))
        for j in range(n_genes):
            cmp_color = SEQUENTIAL_CMAP(seq_norm(b["cmp_vals"][j]))
            ax.add_patch(Rectangle((j * col_w - 0.46, y_cmp - row_h / 2), 0.92, row_h,
                                    facecolor=cmp_color, edgecolor="white", linewidth=2.2, zorder=2))
        for j in range(n_genes):
            fc = b["log2fc"][j]
            delta_color = DIVERGING_CMAP(div_norm(fc)) if np.isfinite(fc) else "#f0f0f0"
            ax.add_patch(Rectangle((j * col_w - 0.46, y_delta - delta_h / 2), 0.92, delta_h,
                                    facecolor=delta_color, edgecolor="white", linewidth=1.6, zorder=2))

        # pairing bracket -- links ONLY the reference+comparison pair, the
        # Delta strip sits visibly separated below it (its own meaning)
        bx0, bx1 = -0.70, -0.58
        ax.plot([bx0, bx0], [y_cmp, y_ref], color=MGRAY, linewidth=1.5, solid_capstyle="round", zorder=3)
        ax.plot([bx0, bx1], [y_ref, y_ref], color=MGRAY, linewidth=1.5, solid_capstyle="round", zorder=3)
        ax.plot([bx0, bx1], [y_cmp, y_cmp], color=MGRAY, linewidth=1.5, solid_capstyle="round", zorder=3)

        ax.text(-0.80, y_ref, b["ref_label"], ha="right", va="center", fontsize=10.5, color=MGRAY)
        ax.text(-0.80, y_cmp, b["cmp_label"], ha="right", va="center", fontsize=11, color=DGRAY, fontweight="bold")
        ax.text(-0.80, y_delta, "Δ", ha="right", va="center", fontsize=12, color=DGRAY, fontweight="bold",
                 style="italic")

        y_mid = (y_ref + y_cmp) / 2
        ax.text(-2.55, y_mid, b["dataset"], ha="left", va="center", fontsize=16, fontweight="bold", color=DGRAY)

    # gene headers
    header_y = top_y - 0.08
    for j, gene in enumerate(genes):
        ax.text(j * col_w, header_y, gene, ha="center", va="bottom", fontsize=20, fontweight="bold",
                 color=FOCUS_COLORS.get(gene, DGRAY))

    ax.set_xlim(-2.75, (n_genes - 1) * col_w + 0.62)
    ax.set_ylim(bottom_y - 0.05, top_y + 0.5)
    ax.axis("off")

    # ---- two compact, unambiguous legends ----
    seq_sm = ScalarMappable(norm=Normalize(vmin=0, vmax=1), cmap=SEQUENTIAL_CMAP)
    seq_sm.set_array([])
    div_sm = ScalarMappable(norm=div_norm, cmap=DIVERGING_CMAP)
    div_sm.set_array([])

    cax1 = fig.add_axes([0.15, -0.05, 0.30, 0.026])
    cb1 = fig.colorbar(seq_sm, cax=cax1, orientation="horizontal")
    cb1.set_ticks([0, 1])
    cb1.set_ticklabels(["low", "high"])
    cb1.ax.tick_params(labelsize=8.6, length=0, color=MGRAY, labelcolor=MGRAY)
    cb1.outline.set_visible(False)
    fig.text(0.15, -0.075, "Expression -- low → high, scaled WITHIN each dataset block",
              fontsize=8.8, color=MGRAY, va="top")

    cax2 = fig.add_axes([0.57, -0.05, 0.30, 0.026])
    cb2 = fig.colorbar(div_sm, cax=cax2, orientation="horizontal")
    cb2.set_ticks([-vmax, 0, vmax])
    cb2.set_ticklabels(["down", "0", "up"])
    cb2.ax.tick_params(labelsize=8.6, length=0, color=DGRAY, labelcolor=DGRAY)
    cb2.outline.set_visible(False)
    fig.text(0.57, -0.075, "Δ log$_2$FC -- down → up, ONE shared scale across all datasets",
              fontsize=8.8, color=DGRAY, va="top")

    fig.text(0.0, 1.06, "Cross-context expression dynamics of four tamoxifen-sensitisation candidates",
              fontsize=21.5, fontweight="bold", color=DGRAY, ha="left", va="bottom")
    fig.text(0.0, 1.015,
              "Expression states are scaled within dataset; Δ rows show cross-context-comparable log$_2$ fold-change.",
              fontsize=11.2, color="#555555", ha="left", va="bottom")

    stub.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(stub.with_suffix(".png"), facecolor="white", bbox_inches="tight", dpi=300)
    fig.savefig(stub.with_suffix(".pdf"), facecolor="white", bbox_inches="tight")
    fig.savefig(stub.with_suffix(".svg"), facecolor="white", bbox_inches="tight")
    plt.close(fig)
    logger.info("wrote %s.png/.pdf/.svg", stub)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    build_hero_heatmap_v3(OUT_DIR / "HERO_main_heatmap_v3")
