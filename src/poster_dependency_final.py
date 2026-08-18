"""Poster figure 6 (final): baseline dependency, as counts.

post_freeze_exploratory. No new dependency analysis: the per-cell-line table
and the strong-dependency criterion (DepMap Public 26Q1 dependency
probability > 0.5, from config) are read unchanged through
`poster_depmap_v1.load_cellline_table()` / `.dependency_summary()`, the same
frozen loaders the frozen poster figure 05 uses, and are verified before
anything is drawn.

WHAT THIS FIXES. The frozen figure 05 is a scatter whose x-axis is the CRISPR
effect size -- restating figure 1 -- against a percentage, carrying one number
per gene on two axes with no numeric x scale a reader can use. Counts out of a
stated denominator are the honest form for 11 cell lines: a percentage of 11
implies a precision the sample size does not support, and 81.8% is simply 9 of
11. Bars are sorted so the outlier is immediate.

WHAT THE SCATTER CARRIED THAT BARS CANNOT. Its geometry asserted that
sensitisation and baseline requirement are separate measurements that need not
coincide. That claim is recovered in the panel title and in the annotation on
KDM1A, which is the strongest sensitiser in the screen and is required by none
of the 11 lines. Note the wording: these are two distinct measurements that
need not coincide. This figure does not establish statistical independence,
and does not claim it.

WHAT HIGH DEPENDENCY MEANS HERE. It is not an advantage. A gene most cell
lines already need at baseline offers a narrower, less tamoxifen-specific
therapeutic window, and that is why TLK2's 9 of 11 is presented as a
limitation rather than a strength. Baseline dependency in cancer cell lines is
also not normal-tissue safety and says nothing about toxicity.
"""

from __future__ import annotations

import logging
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from src import poster_depmap_v1 as d1
from src.poster_final_common import (
    FONT,
    OUT_DIR,
    figure_footer,
    headline,
    pin_reproducibility,
    save,
    style_axes,
    verify,
)
from src.poster_palette import GENE_COLOURS, NEUTRAL, WHITE

logger = logging.getLogger(__name__)

FIGURE = "F6_baseline_dependency"
N_LINES = 11
REFERENCE_COUNTS = {"TLK2": 9, "VEZF1": 3, "KDM1A": 0, "USP34": 0}


def load_dependency():
    table = d1.load_cellline_table()
    summary = d1.dependency_summary(table)
    logger.info("dependency: %d rows over %d genes", len(table), summary["gene"].nunique())
    return summary


def gate(summary):
    idx = summary.set_index("gene")
    checks = []
    for gene, n in REFERENCE_COUNTS.items():
        checks.append((f"{gene}_n_strongly_dependent", int(idx.loc[gene, "n_strongly_dependent"]), n, 0))
        checks.append((f"{gene}_n_lines", int(idx.loc[gene, "n_lines"]), N_LINES, 0))
    checks.append(("all_four_genes_present", len(summary), 4, 0))
    return verify(FIGURE, checks)


def build(stub: Path):
    pin_reproducibility(FIGURE)
    summary = load_dependency()
    verification = gate(summary)

    ordered = summary.sort_values("n_strongly_dependent", ascending=True)
    genes = ordered["gene"].tolist()
    counts = ordered["n_strongly_dependent"].astype(int).tolist()

    fig, ax = plt.subplots(figsize=(15.5, 8.0))
    fig.subplots_adjust(left=0.125, right=0.975, top=0.965, bottom=0.275)
    style_axes(ax, grid_axis="x")

    y = range(len(genes))
    # A faint full-width bar shows the denominator behind every gene, so a
    # count of 0 is visibly "0 out of 11" rather than an absent bar.
    ax.barh(list(y), [N_LINES] * len(genes), height=0.62, color=NEUTRAL["grid"], zorder=1)
    ax.barh(list(y), counts, height=0.62,
            color=[GENE_COLOURS[g] for g in genes], zorder=2)

    for i, (gene, n) in enumerate(zip(genes, counts)):
        ax.text(n + 0.22, i, f"{n} of {N_LINES}", va="center", ha="left",
                fontsize=FONT["panel"], fontweight="bold", color=NEUTRAL["ink"], zorder=4)
    ax.set_yticks(list(y))
    ax.set_yticklabels(genes, fontsize=FONT["panel"] + 2, fontweight="bold")
    for tick, gene in zip(ax.get_yticklabels(), genes):
        tick.set_color(GENE_COLOURS[gene])
    ax.set_xlim(0, N_LINES + 1.6)
    ax.set_xticks(range(0, N_LINES + 1, 2))
    ax.tick_params(axis="x", labelsize=FONT["tick"])
    ax.set_xlabel(f"Cell lines strongly dependent at baseline, of {N_LINES} luminal lines",
                  fontsize=FONT["axis"], color=NEUTRAL["ink_2"], labelpad=12)

    # Annotations are attached by GENE, not by a hardcoded row index, so they
    # cannot drift onto the wrong bar if the sort order ever changes.
    row = {g: i for i, g in enumerate(genes)}
    ax.annotate("strongest sensitiser in the screen,\nyet required by none of the 11 lines",
                xy=(2.10, row["KDM1A"] + 0.26), xytext=(3.2, row["KDM1A"] + 0.60),
                fontsize=FONT["annot"], color=NEUTRAL["ink_2"], va="center", ha="left",
                arrowprops=dict(arrowstyle="-", color=NEUTRAL["rule"], lw=1.0,
                                shrinkA=2, shrinkB=6))
    ax.annotate("needed by most lines with or without the drug:\na narrower, less tamoxifen-specific window",
                xy=(8.4, row["TLK2"]), xytext=(3.2, row["TLK2"] - 0.62),
                fontsize=FONT["annot"], color=NEUTRAL["ink_2"], va="center", ha="left",
                arrowprops=dict(arrowstyle="-", color=NEUTRAL["rule"], lw=1.0,
                                shrinkA=2, shrinkB=6))

    headline(
        fig,
        "Needing a gene at baseline and being sensitised by losing it are two\ndifferent measurements, and they need not coincide",
        "TLK2 is required for survival by 9 of 11 oestrogen-receptor-positive cell lines whether or not tamoxifen is present. That is a limitation,\n"
        "not a strength: a gene most cells already need offers a narrower, less tamoxifen-specific therapeutic window. KDM1A is the opposite —\n"
        "the strongest sensitiser in the screen, required by none of the 11. Counts, not percentages, out of the 11 lines with a dependency value.\n"
        "Baseline dependency in cancer cell lines is not normal-tissue safety and says nothing about toxicity.",
        key=FIGURE)
    figure_footer(fig, "DepMap 26Q1, dependency probability > 0.5; 11 evaluable lines.")
    return save(fig, stub), verification


def main(out_dir: Path = OUT_DIR):
    return build(out_dir / FIGURE)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    main()
