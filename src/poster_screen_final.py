"""Poster figure 2 (final): the genome-wide screen, with CERTAINTY visible.

post_freeze_exploratory. No new statistics: every value is read at render
time from the frozen CRISPR loaders
(`post_audit_sensitivity_data.load_significant_sensitising_hits()` and
`.load_genomewide_crispr()`), the same loaders the frozen poster figure 01
uses, and is verified against them before anything is drawn.

WHAT THIS FIXES. The frozen figure 01 sorts the 13 hits by effect size
alone. A reader cannot see why USP34 (rank 12 of 13 by effect) and VEZF1
(rank 8) were carried forward while USP17L29 (rank 2) and TADA2B (rank 3)
were not -- and a careful reader was in fact confused by exactly this.
Two things were invisible and are now drawn:

  1. CERTAINTY. Effect size and false discovery rate are plotted as two
     axes, so a large effect with weak certainty is visibly distinct from a
     small effect with strong certainty. KDM1A and TLK2 sit far to the
     certain side; USP17L29 and TADA2B have large effects but sit with the
     rest of the pack on certainty.
  2. PROVENANCE. The four poster genes did NOT come from one selection
     rule, and this figure does not imply they did. USP34 and VEZF1 entered
     through the frozen multimodal rule (screen hit + at least one
     corroborating RNA/human dataset at FDR<0.05); KDM1A and TLK2 were added
     after an external audit challenged that rule, as the two strongest
     screen hits, and have zero qualifying corroboration. Marker shape
     carries that distinction; EML5 and CITED2, displaced from the original
     frozen shortlist, are named in the note.

Nothing here re-ranks, re-tests, or changes the pre-specified gate
(FDR < 0.10 with a sensitising direction, PREANALYSIS.md Section 4 --
written and dated before results existed, but not a public preregistration).
"""

from __future__ import annotations

import logging
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from src import post_audit_sensitivity_data as pad
from src.poster_final_common import (
    FONT,
    OUT_DIR,
    PROVENANCE,
    PROVENANCE_MARKER,
    SCREEN_FDR,
    figure_footer,
    headline,
    pin_reproducibility,
    save,
    style_axes,
    verify,
)
from src.poster_palette import GENE_COLOURS, NEUTRAL, WHITE

logger = logging.getLogger(__name__)

FIGURE = "F2_screen_certainty"
FOCUS = ("KDM1A", "TLK2", "USP34", "VEZF1")

# Frozen reference values (docs/THERAPEUTIC_SHORTLIST_FREEZE.md and the
# post-audit report). The gate below checks the loaded table against these
# before the figure is drawn.
REFERENCE = {
    "KDM1A": (-2.167336, 0.000385),
    "TLK2": (-1.848198, 0.001676),
    "VEZF1": (-1.602445, 0.037258),
    "USP34": (-1.391298, 0.041685),
}
REFERENCE_N_HITS = 13
REFERENCE_N_GENOMEWIDE = 19103

# Fixed label anchors in DATA coordinates, hand-placed once against the
# frozen values. No collision solver is used anywhere, so label positions are
# byte-identical on every render. Five hits share FDR 0.037258 and therefore
# sit on one horizontal line; their labels alternate above and below it.
LABEL_POS: dict[str, tuple[float, float, str]] = {
    "KDM1A": (2.10, 3.62, "right"), "TLK2": (1.79, 2.98, "right"),
    "TADA2B": (2.05, 2.13, "center"), "USP17L29": (2.15, 1.24, "center"),
    "PET117": (1.83, 1.06, "center"), "SUPT4H1": (1.71, 1.72, "left"),
    "TSR3": (1.62, 0.92, "center"), "VEZF1": (1.63, 1.80, "left"),
    "EIF4ENIF1": (1.52, 1.63, "right"), "HMGB1": (1.44, 1.79, "right"),
    "CTDNEP1": (1.47, 0.95, "left"), "USP34": (1.30, 1.66, "right"),
    "ICK": (1.225, 1.22, "center"),
}


def load_screen():
    """Frozen screen hits + genome-wide gene count, unmodified."""
    hits = pad.load_significant_sensitising_hits().copy()
    genomewide = pad.load_genomewide_crispr()
    logger.info("screen: %d genome-wide fitted genes, %d hits at FDR<%.2f",
                len(genomewide), len(hits), SCREEN_FDR)
    return hits, len(genomewide)


def gate(hits, n_genomewide):
    checks = [("n_significant_hits", len(hits), REFERENCE_N_HITS, 0),
              ("n_genomewide_fitted", n_genomewide, REFERENCE_N_GENOMEWIDE, 0)]
    idx = hits.set_index("gene")
    for gene, (eff, fdr) in REFERENCE.items():
        checks.append((f"{gene}_effect_size", float(idx.loc[gene, "effect_size"]), eff, 1e-6))
        checks.append((f"{gene}_fdr", float(idx.loc[gene, "fdr"]), fdr, 1e-6))
    # the pre-specified gate itself must hold for every plotted hit
    checks.append(("max_fdr_among_hits", float(hits["fdr"].max()), 0.0, SCREEN_FDR))
    checks.append(("max_effect_among_hits", float(hits["effect_size"].max()), -1.0, 1.0))
    # USP34's rank by effect is a claim made on the poster
    ranks = hits.sort_values("effect_size").reset_index(drop=True)
    usp34_rank = int(ranks.index[ranks["gene"] == "USP34"][0]) + 1
    checks.append(("usp34_rank_by_effect", usp34_rank, 12, 0))
    return verify(FIGURE, checks)


def build(stub: Path):
    pin_reproducibility(FIGURE)
    hits, n_genomewide = load_screen()
    verification = gate(hits, n_genomewide)

    fig, ax = plt.subplots(figsize=(15.5, 9.4))
    fig.subplots_adjust(left=0.095, right=0.985, top=0.955, bottom=0.245)

    x = -hits["effect_size"].to_numpy()          # display flip: right = stronger
    y = -np.log10(hits["fdr"].to_numpy())
    genes = hits["gene"].tolist()

    style_axes(ax, grid_axis="both")
    gate_y = -np.log10(SCREEN_FDR)
    ax.axhline(gate_y, color=NEUTRAL["rule"], linestyle="--", linewidth=1.4, zorder=2)
    ax.text(0.985, gate_y + 0.05, f"pre-specified gate: FDR < {SCREEN_FDR:.2f}",
            fontsize=FONT["annot"], color=NEUTRAL["ink_2"], va="bottom", ha="right",
            transform=ax.get_yaxis_transform())

    for gene, gx, gy in zip(genes, x, y):
        focus = gene in FOCUS
        if focus:
            marker = PROVENANCE_MARKER[PROVENANCE[gene]]
            ax.scatter([gx], [gy], s=430, c=GENE_COLOURS[gene], marker=marker,
                       edgecolors=WHITE, linewidths=2.0, zorder=12)
        else:
            ax.scatter([gx], [gy], s=190, c=NEUTRAL["backdrop"], marker="o",
                       edgecolors=WHITE, linewidths=1.2, zorder=6)
        lx, ly, ha = LABEL_POS[gene]
        ax.annotate(gene, (gx, gy), xytext=(lx, ly),
                    ha=ha, va="center", fontsize=FONT["annot"] + (3 if focus else 0),
                    fontweight="bold" if focus else "normal",
                    color=GENE_COLOURS[gene] if focus else NEUTRAL["ink_2"], zorder=15,
                    arrowprops=dict(arrowstyle="-", color=NEUTRAL["rule"], lw=0.9,
                                    shrinkA=2, shrinkB=7))

    ax.set_xlabel("Sensitisation strength   (screen effect size, sign-flipped: further right = stronger)",
                  fontsize=FONT["axis"], color=NEUTRAL["ink_2"], labelpad=10)
    ax.set_ylabel("Certainty   $-\\log_{10}$(false discovery rate)", fontsize=20,
                  color=NEUTRAL["ink_2"], labelpad=10)
    ax.tick_params(labelsize=FONT["tick"])
    ax.set_xlim(1.13, 2.28)
    ax.set_ylim(0.72, 3.95)

    # provenance legend, inside the plot area
    from matplotlib.lines import Line2D
    handles = [
        Line2D([], [], marker="o", linestyle="none", markersize=15, markerfacecolor=NEUTRAL["ink"],
               markeredgecolor=WHITE, label="frozen rule"),
        Line2D([], [], marker="s", linestyle="none", markersize=14, markerfacecolor=NEUTRAL["ink"],
               markeredgecolor=WHITE, label="added post-audit"),
        Line2D([], [], marker="o", linestyle="none", markersize=11, markerfacecolor=NEUTRAL["backdrop"],
               markeredgecolor=WHITE, label="other screen hit"),
    ]
    leg = fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.55, 0.070),
                     ncol=3, frameon=False, fontsize=FONT["legend"], handletextpad=0.8,
                     columnspacing=4.0, borderpad=0.3)
    for t in leg.get_texts():
        t.set_color(NEUTRAL["ink_2"])

    headline(
        fig,
        "Effect size and certainty are different things, and the four\ncandidates were not chosen on effect size alone",
        f"{len(hits)} of {n_genomewide:,} genes met the pre-specified gate. KDM1A and TLK2 are the most certain hits; USP34 sits 12th of 13 by effect size.\n"
        "Marker shape shows how each candidate entered the set — the four did not come from one selection rule. EML5 and CITED2 were\n"
        "displaced from the original frozen shortlist by the same reinterpretation. Colour identifies the gene, never its rank.",
        key=FIGURE)
    figure_footer(fig, "Frozen CRISPR screen (Hany et al. 2023); FDR across 19,103 genes.")
    return save(fig, stub), verification


def main(out_dir: Path = OUT_DIR):
    written, verification = build(out_dir / FIGURE)
    return written, verification


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    main()
