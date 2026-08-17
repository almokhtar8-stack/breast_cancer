"""ONE poster-grade hero heatmap, v6: v5 with the dead white space
squeezed out so the colored heatmap body dominates the canvas.

v6 changes LAYOUT/GEOMETRY ONLY. It imports the exact same row-building
functions, gene-wise within-dataset z-score values, cmap, and state
colors from `poster_hero_heatmap_v4`/`poster_hero_heatmap_v5` -- nothing
about the data, the transform, row order, or gene order is touched. What
changes relative to v5:

  - wider gene columns and taller rows (bigger, more substantial cells)
  - a tighter left-side label/bracket/annotation column
  - shrunk vertical gaps between rows, sub-groups, and dataset blocks
  - a smaller top title band and a single-line, compact bottom legend

See results/reports/poster_hero_heatmap_v6/NOTE.md for the full
rationale and the exact set of geometry constants changed from v5.
"""

from __future__ import annotations

import logging
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import TwoSlopeNorm
from matplotlib.patches import Rectangle

from src import poster_hero_heatmap_v4 as v4
from src import poster_hero_heatmap_v5 as v5

logger = logging.getLogger(__name__)

OUT_DIR = Path("results/figures/poster_hero_heatmap_v6")

FOCUS_FOUR = v4.FOCUS_FOUR
FOCUS_COLORS = v4.FOCUS_COLORS
DIVERGING_CMAP = v4.DIVERGING_CMAP
STATE_COLORS = v4.STATE_COLORS
DATASET_BUILDERS = v4.DATASET_BUILDERS
# Public-audit wording correction (2026-08-17): v5's GSE245601 headline read
# "Single-cell: before vs 12 h after tamoxifen". The dataset IS single-cell
# RNA-seq, but the six rows plotted here are PER-TUMOUR PSEUDOBULK (3
# patients x Control/Tamoxifen, see v4._build_gse245601_rows), so the old
# headline risked implying cell-level resolution / pseudo-replication. Only
# the displayed label changes; no plotted value, row, or z-score is affected.
CONTEXT_TITLE = dict(v5.CONTEXT_TITLE)
CONTEXT_TITLE["GSE245601"] = "Acute 12 h tamoxifen — per-tumour pseudobulk"
STATE_LABELS = v5.STATE_LABELS
DIVIDER_DATASETS = v5.DIVIDER_DATASETS
BRACKET_DATASETS = v5.BRACKET_DATASETS

DGRAY = "#262626"
MGRAY = "#8c8c8c"


def build_hero_heatmap_v6(stub: Path) -> None:
    blocks = [(name, builder()) for name, builder in DATASET_BUILDERS]
    all_rows = [r for _, rows in blocks for r in rows]

    # Same display limits as v5: no clipping (max |z| = 1.90 < 2), see
    # results/reports/poster_hero_heatmap_v5/NOTE.md.
    zmax = max(abs(r.values[g]) for r in all_rows for g in FOCUS_FOUR)
    div_norm = TwoSlopeNorm(vcenter=0, vmin=-zmax, vmax=zmax)

    # ---- bigger cells, tighter gaps everywhere -- the only geometry
    # difference from v5; the row/gene order and every colored value are
    # untouched. ----
    row_h = 0.40
    row_gap = 0.018
    subgroup_gap = 0.050
    block_gap = 0.22
    header_top_pad = 0.05
    # Header band heights are sized from the actual rendered text extent
    # at this figure's size/fontsize (measured directly, not guessed) so
    # the accession line never collides with the first data row below it.
    title_band = 0.58
    accession_band = 0.38
    prerow_gap = 0.15
    col_w = 1.28
    n_genes = len(FOCUS_FOUR)
    n_rows = len(all_rows)

    cell_left = -0.46 * col_w
    cell_right_edge = (n_genes - 1) * col_w + 0.46 * col_w

    ann_x0, ann_x1 = cell_left - 0.16, cell_left - 0.03
    label_x = ann_x0 - 0.05
    # gap sized to fit the longest row label ("ZR-75-1_Tam1") plus the
    # bracket itself without ever crossing the label text
    bracket_x1 = label_x - 0.95
    bracket_x0 = bracket_x1 - 0.12
    header_x = bracket_x0 - 0.35
    grouplabel_x = header_x  # reuse the header column -- ample clearance from row labels

    y = 0.0
    row_y: dict[tuple[str, str], float] = {}
    header_positions: dict[str, tuple[float, float, float]] = {}
    divider_ys: dict[str, list[float]] = {}
    group_labels: list[tuple[float, str]] = []

    for dataset, rows in blocks:
        y -= header_top_pad
        title_y = y
        y -= title_band
        accession_y = y
        y -= accession_band
        rule_y = y + prerow_gap * 0.45
        y -= prerow_gap
        prev_subgroup = None
        first_of_subgroup_y: dict[int, float] = {}
        for r in rows:
            if prev_subgroup is not None and r.subgroup != prev_subgroup:
                y -= subgroup_gap
                if dataset in DIVIDER_DATASETS:
                    divider_ys.setdefault(dataset, []).append(y + subgroup_gap / 2)
            row_y[(dataset, r.label)] = y
            if r.subgroup not in first_of_subgroup_y:
                first_of_subgroup_y[r.subgroup] = y
            prev_subgroup = r.subgroup
            y -= (row_h + row_gap)
        header_positions[dataset] = (title_y, accession_y, rule_y)
        y -= block_gap

        if dataset == "GSE240112":
            group_labels.append((first_of_subgroup_y[0], "PRIMARY"))
            group_labels.append((first_of_subgroup_y[1], "RECURRENT"))

    top_y = -header_top_pad + title_band + 0.06
    bottom_y = y + block_gap - row_h * 0.5

    fig_w = 13.2
    fig_h = 8.0
    fig, ax = plt.subplots(figsize=(fig_w, fig_h), dpi=300)
    # Reclaim matplotlib's default axes margins so the heatmap -- not the
    # default ~12% empty border -- defines the usable canvas width.
    fig.subplots_adjust(left=0.005, right=0.995, top=0.995, bottom=0.005)

    for dataset, rows in blocks:
        for r in rows:
            yy = row_y[(dataset, r.label)]
            for j, gene in enumerate(FOCUS_FOUR):
                color = DIVERGING_CMAP(div_norm(r.values[gene]))
                ax.add_patch(Rectangle((j * col_w - 0.46 * col_w, yy - row_h / 2), 0.92 * col_w, row_h,
                                        facecolor=color, edgecolor="white", linewidth=1.0, zorder=2))
            ax.add_patch(Rectangle((ann_x0, yy - row_h / 2), ann_x1 - ann_x0, row_h,
                                    facecolor=STATE_COLORS[r.state], edgecolor="white", linewidth=0.6, zorder=2))
            ax.text(label_x, yy, r.label, ha="right", va="center", fontsize=8.6, color=DGRAY)

        if dataset in BRACKET_DATASETS:
            anchors: dict[str, list[str]] = {}
            for r in rows:
                if r.pair_anchor is not None:
                    anchors.setdefault(r.pair_anchor, []).append(r.label)
            for anchor_label, child_labels in anchors.items():
                y_anchor = row_y[(dataset, anchor_label)]
                child_ys = [row_y[(dataset, c)] for c in child_labels]
                ax.plot([bracket_x0, bracket_x0], [min(child_ys), y_anchor], color=MGRAY, linewidth=1.0,
                        solid_capstyle="round", zorder=3)
                ax.plot([bracket_x0, bracket_x1], [y_anchor, y_anchor], color=MGRAY, linewidth=1.0,
                        solid_capstyle="round", zorder=3)
                for cy in child_ys:
                    ax.plot([bracket_x0, bracket_x1], [cy, cy], color=MGRAY, linewidth=1.0,
                            solid_capstyle="round", zorder=3)

        if dataset in divider_ys:
            for yy in divider_ys[dataset]:
                ax.plot([label_x - 0.02, cell_right_edge], [yy, yy], color="#e6e6e6",
                        linewidth=0.9, linestyle=(0, (2, 2)), zorder=1)

    for yy, text in group_labels:
        ax.text(grouplabel_x, yy, text, ha="left", va="center", fontsize=8.0, fontweight="bold", color=MGRAY)

    for dataset, (title_y, accession_y, rule_y) in header_positions.items():
        ax.text(header_x, title_y, CONTEXT_TITLE[dataset], ha="left", va="top", fontsize=14.5,
                 fontweight="bold", color=DGRAY)
        ax.text(header_x, accession_y, dataset, ha="left", va="top", fontsize=9.0, color=MGRAY,
                 style="italic")
        ax.plot([header_x - 0.04, cell_right_edge], [rule_y] * 2, color="#e6e6e6", linewidth=1.0, zorder=1)

    header_y = top_y - 0.06
    for j, gene in enumerate(FOCUS_FOUR):
        ax.text(j * col_w, header_y, gene, ha="center", va="bottom", fontsize=23, fontweight="bold",
                 color=FOCUS_COLORS.get(gene, DGRAY))

    ax.set_xlim(header_x - 0.10, cell_right_edge + 0.08)
    ax.set_ylim(bottom_y - 0.01, top_y + 0.68)
    ax.axis("off")

    div_sm = ScalarMappable(norm=div_norm, cmap=DIVERGING_CMAP)
    div_sm.set_array([])
    cax = fig.add_axes([0.085, -0.020, 0.19, 0.020])
    cb = fig.colorbar(div_sm, cax=cax, orientation="horizontal")
    cb.set_ticks([-zmax, 0, zmax])
    cb.set_ticklabels(["Low", "0", "High"])
    cb.ax.tick_params(labelsize=8.4, length=0, color=DGRAY, labelcolor=DGRAY)
    cb.outline.set_visible(False)
    fig.text(0.085, 0.010, "Relative expression within dataset", fontsize=9.0, fontweight="bold",
              color=DGRAY, va="bottom")
    fig.text(0.285, -0.025, "Gene-wise z-score", fontsize=7.6, color=MGRAY, va="center", style="italic")

    legend_x = 0.46
    for i, state in enumerate(["baseline", "resistant", "recurrent", "acute_tam"]):
        sx = legend_x + i * 0.135
        sy = -0.020
        fig.add_artist(plt.Rectangle((sx, sy), 0.013, 0.014, transform=fig.transFigure,
                                      facecolor=STATE_COLORS[state], edgecolor="none"))
        fig.text(sx + 0.020, sy + 0.007, STATE_LABELS[state], fontsize=7.5, color=DGRAY, va="center")

    # Title fontsize is capped so its rendered width fits within the
    # reclaimed axes width (measured directly, not guessed) -- otherwise
    # bbox_inches='tight' expands the canvas to fit the title and leaves
    # a large dead zone to the right of the actual heatmap.
    fig.text(0.0, 1.028, "Candidate expression across resistance, recurrence and acute tamoxifen response",
              fontsize=18.5, fontweight="bold", color=DGRAY, ha="left", va="bottom")
    fig.text(0.0, 1.006,
              "Rows are real biological observations; color shows gene-wise standardized expression within each dataset.",
              fontsize=10.2, color="#555555", ha="left", va="bottom")

    stub.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(stub.with_suffix(".png"), facecolor="white", bbox_inches="tight", pad_inches=0.08, dpi=300)
    fig.savefig(stub.with_suffix(".pdf"), facecolor="white", bbox_inches="tight", pad_inches=0.08)
    fig.savefig(stub.with_suffix(".svg"), facecolor="white", bbox_inches="tight", pad_inches=0.08)
    plt.close(fig)
    logger.info("wrote %s.png/.pdf/.svg (%d real biological rows)", stub, n_rows)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    build_hero_heatmap_v6(OUT_DIR / "HERO_sample_level_heatmap_v6")
