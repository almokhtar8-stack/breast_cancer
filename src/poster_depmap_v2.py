"""ONE poster-grade DepMap figure (v2) -- a communication rebuild of v1.

v1 showed the full 11 x 4 Chronos heatmap, which is scientifically
complete but forces the viewer through Chronos, colour intensity,
dependency probability, ring markers, 11 cell-line rows and a second bar
plot before reaching the conclusion. v2 shows the conclusion directly:
one dot per candidate in a 2D comparison of

    x = tamoxifen sensitisation strength  (frozen CRISPR)
    y = % of ER+/luminal lines with strong baseline dependency  (DepMap)

answering "does tamoxifen sensitisation occur in genes that cancer cells
already depend on at baseline?"

Both axes come from the exact frozen sources already validated in v1 --
no DepMap re-run, no new release, no new threshold, no changed ranking:
  - x: `post_audit_sensitivity_data.load_significant_sensitising_hits()`
    `effect_size`, displayed as `-1 * effect_size`. This sign flip is a
    DISPLAY transform only (a more negative CRISPR effect means a stronger
    sensitising knockout, so negating it makes "further right = stronger"
    intuitive); the frozen values themselves are never modified.
  - y: `poster_depmap_v1.load_cellline_table()` /
    `dependency_summary()`, i.e. the frozen ER+/luminal 11-line subset and
    the frozen strong-dependency criterion (DepMap 26Q1 dependency
    probability > 0.5 from config).

Interpretation guard rails: nothing here says safe/superior/bad/toxic, and
baseline dependency in cancer cell lines is NOT normal-tissue safety and
NOT proof of tamoxifen specificity. The figure only separates
treatment-associated sensitisation from general cancer-cell dependency.
"""

from __future__ import annotations

import logging
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

from src import post_audit_sensitivity_data as pad
from src import poster_depmap_v1 as dm1

logger = logging.getLogger(__name__)

OUT_DIR = Path("results/figures/poster_depmap_v2")

FOCUS_FOUR = dm1.FOCUS_FOUR            # ["KDM1A", "TLK2", "USP34", "VEZF1"]
FOCUS_COLORS = dm1.FOCUS_COLORS        # frozen Okabe-Ito poster identity colors
RELEASE = dm1.RELEASE

DGRAY = "#262626"
MGRAY = "#8c8c8c"
LOW_BAND = "#F5F7F4"
HIGH_BAND = "#F3EFF6"

# Label offsets (points) per candidate -- pure text placement so the four
# labels never overlap each other or the axes; the DOT positions are always
# the unmodified data coordinates.
LABEL_OFFSETS = {
    "KDM1A": (0, 34),
    "TLK2": (0, -46),
    "USP34": (0, 34),
    "VEZF1": (0, 34),
}


def build_plot_table() -> pd.DataFrame:
    """One row per candidate with the exact plotted coordinates, both
    derived from frozen sources."""
    hits = pad.load_significant_sensitising_hits().set_index("gene")
    summary = dm1.dependency_summary(dm1.load_cellline_table()).set_index("gene")

    rows = []
    for gene in FOCUS_FOUR:
        effect = float(hits.loc[gene, "effect_size"])
        rows.append({
            "gene": gene,
            "crispr_effect": effect,
            "sensitisation_strength": -1.0 * effect,   # display transform only
            "n_lines": int(summary.loc[gene, "n_lines"]),
            "n_strongly_dependent": int(summary.loc[gene, "n_strongly_dependent"]),
            "pct_strongly_dependent": float(summary.loc[gene, "pct_strongly_dependent"]),
            "median_chronos": float(summary.loc[gene, "median_chronos"]),
        })
    table = pd.DataFrame(rows)
    logger.info("plot table: %d candidates; x from frozen CRISPR effect, y from frozen DepMap %s",
                len(table), RELEASE)
    return table


def build_depmap_v2(stub: Path) -> None:
    table = build_plot_table()
    n_lines = int(table["n_lines"].iloc[0])
    assert (table["n_lines"] == n_lines).all(), "all candidates must share the same evaluable subset"

    fig, ax = plt.subplots(figsize=(11.2, 8.4), dpi=300)
    fig.subplots_adjust(left=0.115, right=0.965, top=0.815, bottom=0.195)

    x_min = table["sensitisation_strength"].min() - 0.22
    x_max = table["sensitisation_strength"].max() + 0.22
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(-9, 100)

    # Subtle interpretation bands -- no invented cutoff on either axis, just
    # a visual reminder of what "up" means.
    ax.axhspan(-9, 33, facecolor=LOW_BAND, edgecolor="none", zorder=0)
    ax.axhspan(66, 100, facecolor=HIGH_BAND, edgecolor="none", zorder=0)
    ax.text(x_min + 0.04, 30.5, "Low baseline dependency", fontsize=11, color="#7f8c7f",
             ha="left", va="top", style="italic")
    ax.text(x_min + 0.04, 96.5, "High baseline dependency", fontsize=11, color="#8b7f96",
             ha="left", va="top", style="italic")

    for y in (25, 50, 75, 100):
        ax.axhline(y, color="#e3e3e3", linewidth=0.8, zorder=1)

    for row in table.itertuples():
        color = FOCUS_COLORS[row.gene]
        ax.scatter([row.sensitisation_strength], [row.pct_strongly_dependent], s=1500,
                    color=color, edgecolor="white", linewidth=2.6, zorder=4)
        dx, dy = LABEL_OFFSETS[row.gene]
        ax.annotate(
            f"{row.gene}\n{row.n_strongly_dependent}/{row.n_lines} dependent",
            xy=(row.sensitisation_strength, row.pct_strongly_dependent),
            xytext=(dx, dy), textcoords="offset points",
            fontsize=15, fontweight="bold", color=color, ha="center",
            va="bottom" if dy > 0 else "top", zorder=5, linespacing=1.35,
        )

    ax.set_xlabel("Tamoxifen sensitisation strength", fontsize=14, color=DGRAY, labelpad=36)
    ax.set_ylabel(f"ER+/luminal cell lines with strong\nbaseline dependency (%)",
                   fontsize=13.5, color=DGRAY, labelpad=12)
    ax.set_yticks([0, 25, 50, 75, 100])
    ax.set_yticklabels(["0%", "25%", "50%", "75%", "100%"], fontsize=11.5, color=DGRAY)
    ax.set_xticks([])
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    ax.spines["left"].set_color("#c9c9c9")
    ax.spines["bottom"].set_color("#c9c9c9")
    ax.tick_params(axis="y", length=0)

    # Directional cue instead of raw numbers on x -- the viewer never has to
    # reason about negative CRISPR effect sizes.
    ax.annotate("", xy=(x_max - 0.06, -16.5), xytext=(x_min + 0.06, -16.5),
                 xycoords=("data", "data"), textcoords=("data", "data"),
                 annotation_clip=False,
                 arrowprops=dict(arrowstyle="-|>", color="#b6b6b6", linewidth=1.4, shrinkA=0, shrinkB=0))
    ax.text(x_min + 0.06, -14.5, "Weaker", fontsize=11, color=MGRAY, ha="left", va="bottom",
             clip_on=False)
    ax.text(x_max - 0.06, -14.5, "Stronger", fontsize=11, color=MGRAY, ha="right", va="bottom",
             clip_on=False)

    fig.text(0.035, 0.975, "Tamoxifen sensitisation is distinct from baseline cancer dependency",
              fontsize=20, fontweight="bold", color=DGRAY, ha="left", va="top")
    fig.text(0.035, 0.928,
              "Functional CRISPR sensitisation versus DepMap dependency across ER+/luminal "
              "breast-cancer models.",
              fontsize=11.5, color="#555555", ha="left", va="top")

    fig.text(0.035, 0.030,
              f"Baseline dependency = fraction of {n_lines} ER+/luminal cell lines with DepMap "
              f"dependency probability > 0.5.",
              fontsize=9.6, color=MGRAY, ha="left", va="center")

    stub.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(stub.with_suffix(".png"), facecolor="white", bbox_inches="tight", dpi=300)
    fig.savefig(stub.with_suffix(".pdf"), facecolor="white", bbox_inches="tight")
    fig.savefig(stub.with_suffix(".svg"), facecolor="white", bbox_inches="tight")
    plt.close(fig)
    logger.info("wrote %s.png/.pdf/.svg", stub)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    build_depmap_v2(OUT_DIR / "DEPMAP_v2")
