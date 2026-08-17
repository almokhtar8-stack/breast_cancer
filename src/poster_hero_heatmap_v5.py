"""ONE poster-grade hero heatmap, v5: polishes the v4 sample-level
heatmap into a landscape, poster-grade central figure.

v5 changes PRESENTATION ONLY. It imports the exact same row-building
functions (`DATASET_BUILDERS`) and the exact same gene-wise,
within-dataset z-score transform from `poster_hero_heatmap_v4` --
nothing about the data or the visualization-only transform is
recomputed, re-derived, or altered. What changes:

  - biological context names ("Cell-line resistance model", ...) replace
    accession codes as the primary block headline; the accession is kept
    as small, muted secondary text directly beneath it
  - block headers are attached directly above their rows (not a
    vertically-centered column off to the side)
  - a wide, landscape canvas so the heatmap body dominates the figure
  - a single compact legend instead of two separate colorbars
  - simplified, non-crossing bracket/divider annotations per block

No CRISPR, pathway, DepMap, structural, or clinical content is included
here -- this figure answers exactly one question: how do the four focus
candidates behave across independent real transcriptomic contexts. See
results/reports/poster_hero_heatmap_v5/NOTE.md for the full rationale.
"""

from __future__ import annotations

import logging
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.cm import ScalarMappable
from matplotlib.colors import TwoSlopeNorm
from matplotlib.patches import Rectangle

from src import poster_hero_heatmap_v4 as v4

logger = logging.getLogger(__name__)

OUT_DIR = Path("results/figures/poster_hero_heatmap_v5")

FOCUS_FOUR = v4.FOCUS_FOUR
FOCUS_COLORS = v4.FOCUS_COLORS
DIVERGING_CMAP = v4.DIVERGING_CMAP
STATE_COLORS = v4.STATE_COLORS
DATASET_BUILDERS = v4.DATASET_BUILDERS

DGRAY = "#262626"
MGRAY = "#8c8c8c"

# Biological context reads first; the accession is secondary/technical.
CONTEXT_TITLE = {
    "GSE118713": "Cell-line resistance model",
    "GSE111151": "Independent tamoxifen-resistant sublines",
    "GSE240112": "Primary vs recurrent tumours",
    "GSE245601": "Single-cell: before vs 12 h after tamoxifen",
}

STATE_LABELS = {
    "baseline": "Baseline / untreated",
    "resistant": "Tamoxifen-resistant",
    "recurrent": "Recurrent tumour",
    "acute_tam": "Acute tamoxifen (12 h)",
}

# datasets whose two conditions are grouped by a plain divider (no real
# sample-to-sample pairing exists to bracket)
DIVIDER_DATASETS = {"GSE118713", "GSE240112"}
# datasets with a genuine parent/patient -> derivative pairing to bracket
BRACKET_DATASETS = {"GSE111151", "GSE245601"}


def build_hero_heatmap_v5(stub: Path) -> None:
    blocks = [(name, builder()) for name, builder in DATASET_BUILDERS]
    all_rows = [r for _, rows in blocks for r in rows]

    # Display limits: no clipping applied -- max |z| across all 116
    # gene x sample observations is ~1.90, already inside +/-2, so a
    # +/-2 clip would not move a single color and is omitted (documented
    # in NOTE.md rather than silently applied).
    zmax = max(abs(r.values[g]) for r in all_rows for g in FOCUS_FOUR)
    div_norm = TwoSlopeNorm(vcenter=0, vmin=-zmax, vmax=zmax)

    row_h = 0.30
    row_gap = 0.024
    subgroup_gap = 0.075
    block_gap = 0.34
    # header region above each block's rows, laid out top-to-bottom as
    # explicit, non-overlapping bands (title text, then accession text,
    # then a fixed clearance gap before the first data row)
    header_top_pad = 0.06
    title_band = 0.46
    accession_band = 0.32
    prerow_gap = 0.14
    header_pad = header_top_pad + title_band + accession_band + prerow_gap
    col_w = 1.0
    n_genes = len(FOCUS_FOUR)
    n_rows = len(all_rows)

    ann_x0, ann_x1 = -0.60, -0.49
    label_x = -0.70
    bracket_x0, bracket_x1 = -1.52, -1.39
    header_x = -1.90
    grouplabel_x = -1.50

    y = 0.0
    row_y: dict[tuple[str, str], float] = {}
    header_positions: dict[str, tuple[float, float, float]] = {}  # dataset -> (title_y, accession_y, rule_y)
    divider_ys: dict[str, list[float]] = {}
    group_labels: list[tuple[float, str, str]] = []

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
        block_bottom = y + row_gap
        header_positions[dataset] = (title_y, accession_y, rule_y)
        y -= block_gap

        if dataset == "GSE240112":
            group_labels.append((first_of_subgroup_y[0], "PRIMARY", "baseline"))
            group_labels.append((first_of_subgroup_y[1], "RECURRENT", "recurrent"))

    top_y = -header_top_pad + title_band + 0.10
    bottom_y = y + block_gap - row_h * 0.5

    fig_h = 8.6
    fig_w = 13.6
    fig, ax = plt.subplots(figsize=(fig_w, fig_h), dpi=300)

    for dataset, rows in blocks:
        for r in rows:
            yy = row_y[(dataset, r.label)]
            for j, gene in enumerate(FOCUS_FOUR):
                color = DIVERGING_CMAP(div_norm(r.values[gene]))
                ax.add_patch(Rectangle((j * col_w - 0.46, yy - row_h / 2), 0.92, row_h,
                                        facecolor=color, edgecolor="white", linewidth=0.9, zorder=2))
            ax.add_patch(Rectangle((ann_x0, yy - row_h / 2), ann_x1 - ann_x0, row_h,
                                    facecolor=STATE_COLORS[r.state], edgecolor="white", linewidth=0.6, zorder=2))
            ax.text(label_x, yy, r.label, ha="right", va="center", fontsize=7.6, color=DGRAY)

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
                ax.plot([label_x - 0.02, (n_genes - 1) * col_w + 0.46], [yy, yy], color="#e6e6e6",
                        linewidth=0.9, linestyle=(0, (2, 2)), zorder=1)

    for yy, text, state in group_labels:
        ax.text(grouplabel_x, yy, text, ha="left", va="center", fontsize=7.6, fontweight="bold",
                 color=MGRAY)

    for dataset, (title_y, accession_y, rule_y) in header_positions.items():
        ax.text(header_x, title_y, CONTEXT_TITLE[dataset], ha="left", va="top", fontsize=13.5,
                 fontweight="bold", color=DGRAY)
        ax.text(header_x, accession_y, dataset, ha="left", va="top", fontsize=8.6, color=MGRAY,
                 style="italic")
        ax.plot([header_x - 0.05, (n_genes - 1) * col_w + 0.46], [rule_y] * 2,
                color="#e6e6e6", linewidth=1.0, zorder=1)

    header_y = top_y - 0.10
    for j, gene in enumerate(FOCUS_FOUR):
        ax.text(j * col_w, header_y, gene, ha="center", va="bottom", fontsize=21, fontweight="bold",
                 color=FOCUS_COLORS.get(gene, DGRAY))

    ax.set_xlim(-2.35, (n_genes - 1) * col_w + 0.62)
    ax.set_ylim(bottom_y - 0.04, top_y + 0.36)
    ax.axis("off")

    # ---- one compact legend ----
    div_sm = ScalarMappable(norm=div_norm, cmap=DIVERGING_CMAP)
    div_sm.set_array([])
    cax = fig.add_axes([0.10, -0.045, 0.20, 0.020])
    cb = fig.colorbar(div_sm, cax=cax, orientation="horizontal")
    cb.set_ticks([-zmax, 0, zmax])
    cb.set_ticklabels(["Low", "0", "High"])
    cb.ax.tick_params(labelsize=8.6, length=0, color=DGRAY, labelcolor=DGRAY)
    cb.outline.set_visible(False)
    fig.text(0.10, -0.010, "Relative expression within dataset", fontsize=9.2, fontweight="bold",
              color=DGRAY, va="bottom")
    fig.text(0.10, -0.062, "Gene-wise z-score", fontsize=7.8, color=MGRAY, va="top", style="italic")

    legend_x = 0.40
    for i, state in enumerate(["baseline", "resistant", "recurrent", "acute_tam"]):
        sy = -0.005 - i * 0.016
        fig.add_artist(plt.Rectangle((legend_x, sy), 0.014, 0.012, transform=fig.transFigure,
                                      facecolor=STATE_COLORS[state], edgecolor="none"))
        fig.text(legend_x + 0.022, sy + 0.006, STATE_LABELS[state], fontsize=8.0, color=DGRAY, va="center")

    fig.text(0.0, 1.035, "Candidate expression across resistance, recurrence and acute tamoxifen response",
              fontsize=21, fontweight="bold", color=DGRAY, ha="left", va="bottom")
    fig.text(0.0, 1.010,
              "Rows are real biological observations; color shows gene-wise standardized expression within each dataset.",
              fontsize=11.3, color="#555555", ha="left", va="bottom")

    stub.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(stub.with_suffix(".png"), facecolor="white", bbox_inches="tight", dpi=300)
    fig.savefig(stub.with_suffix(".pdf"), facecolor="white", bbox_inches="tight")
    fig.savefig(stub.with_suffix(".svg"), facecolor="white", bbox_inches="tight")
    plt.close(fig)
    logger.info("wrote %s.png/.pdf/.svg (%d real biological rows)", stub, n_rows)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    build_hero_heatmap_v5(OUT_DIR / "HERO_sample_level_heatmap_v5")
