"""ONE poster-grade CRISPR discovery figure, v2 -- COLOUR CHANGE ONLY.

post_freeze_exploratory. This is `poster_crispr_discovery_v1` with exactly one
difference: the four candidate colours are imported from `src.poster_palette`
(`GENE_COLOURS`) instead of from `post_audit_sensitivity_visualization`
(`FOCUS_COLORS`). Nothing else changes -- same sort order, same bar and marker
geometry, same title and subtitle strings, same typography, same grey values
for the nine non-candidate genes, same axis treatment, same figure size, same
margins, same save calls.

v1 is NOT edited and NOT deleted. It and the shared `FOCUS_COLORS` dict are
still used by other committed figures, so changing them in place would have
silently recoloured those figures; this module exists to avoid exactly that.

The three changed values (KDM1A alone is unchanged):

    TLK2   #CC79A7 -> #009E73
    USP34  #0072B2 -> #6A3D9A
    VEZF1  #E69F00 -> #56B4E9

Rationale is recorded in `src/poster_palette.py` and in
`results/reports/figure1_palette/NOTE.md`.

Data source: `src.post_audit_sensitivity_data.load_significant_sensitising_hits()`
and `.load_genomewide_crispr()`, called unmodified. Both are already-frozen
loaders used elsewhere in this project (the post-audit sensitivity
analysis) -- this module performs NO new discovery, NO re-ranking, and NO
recomputation of the CRISPR effect size, FDR, or the pre-specified
significance gate (FDR < 0.10, PREANALYSIS.md Section 4). Every number
that appears on the figure (hit count, genome-wide gene count, effect
sizes, ranks) is read from these frozen tables at render time -- nothing
is hand-typed. Before anything is drawn, those values are asserted against
the frozen source; a mismatch raises rather than substituting.
"""

from __future__ import annotations

import logging
import os
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from src import post_audit_sensitivity_data as pad
from src import poster_story_v1_data as sv1
from src.poster_palette import GENE_COLOURS

logger = logging.getLogger(__name__)

OUT_DIR = Path("results/figures/poster_crispr_v2")
STATUS_LABEL = "post_freeze_exploratory"

FOCUS_FOUR = sv1.FOCUS_FOUR  # ["KDM1A", "TLK2", "USP34", "VEZF1"]
# THE ONE DIFFERENCE FROM v1: colours come from the palette module, not from
# the shared post_audit_sensitivity_visualization dict.
FOCUS_COLORS = GENE_COLOURS

DGRAY = "#262626"
MGRAY = "#8c8c8c"
MUTED = "#b3b3b3"
MUTED_LABEL = "#6e6e6e"


class VerificationError(AssertionError):
    """Plotted values disagree with the frozen source table."""


# The thirteen genes, in the order the figure plots them (ascending effect
# size, i.e. strongest sensitiser first), with their frozen effect sizes.
# Checked before rendering; a mismatch raises rather than substituting.
FROZEN_ORDER: tuple[tuple[str, float], ...] = (
    ("KDM1A", -2.167336), ("USP17L29", -2.144824), ("TADA2B", -2.064365),
    ("TLK2", -1.848198), ("PET117", -1.827951), ("SUPT4H1", -1.683788),
    ("TSR3", -1.618753), ("VEZF1", -1.602445), ("EIF4ENIF1", -1.550654),
    ("HMGB1", -1.474065), ("CTDNEP1", -1.392970), ("USP34", -1.391298),
    ("ICK", -1.224554),
)
FROZEN_N_GENOMEWIDE = 19103


def verify_against_frozen(hits, n_genomewide: int) -> None:
    """Assert the thirteen genes, their effect sizes and their order match the
    frozen source exactly. Raises listing every mismatch."""
    problems: list[str] = []
    genes = hits["gene"].tolist()
    effects = hits["effect_size"].to_numpy()

    if len(genes) != len(FROZEN_ORDER):
        problems.append(f"hit count: got {len(genes)}, expected {len(FROZEN_ORDER)}")
    if n_genomewide != FROZEN_N_GENOMEWIDE:
        problems.append(f"genome-wide gene count: got {n_genomewide}, expected {FROZEN_N_GENOMEWIDE}")
    for i, (gene, effect) in enumerate(FROZEN_ORDER):
        if i >= len(genes):
            problems.append(f"position {i + 1}: missing, expected {gene}")
            continue
        if genes[i] != gene:
            problems.append(f"position {i + 1}: got {genes[i]!r}, expected {gene!r}")
        elif abs(float(effects[i]) - effect) > 1e-6:
            problems.append(f"{gene} effect size: got {float(effects[i])!r}, expected {effect}")
    if problems:
        raise VerificationError(
            "CRISPR discovery v2: plotted values disagree with the frozen source; "
            "refusing to plot:\n  " + "\n  ".join(problems))
    logger.info("verification passed: %d genes, effect sizes and order match the frozen source",
                len(genes))


def build_crispr_discovery_main(stub: Path) -> None:
    genomewide = pad.load_genomewide_crispr()
    hits = pad.load_significant_sensitising_hits()  # already sorted by rank_by_effect ascending
    n_genomewide = len(genomewide)
    n_hits = len(hits)
    verify_against_frozen(hits, n_genomewide)

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
    # pinned so the PDF is byte-reproducible (matplotlib otherwise stamps a
    # wall-clock /CreationDate)
    os.environ.setdefault("SOURCE_DATE_EPOCH", "1600000000")
    meta = {"Title": "CRISPR discovery v2", "Description": STATUS_LABEL}
    fig.savefig(stub.with_suffix(".png"), facecolor="white", bbox_inches="tight", dpi=300,
                metadata=meta)
    fig.savefig(stub.with_suffix(".pdf"), facecolor="white", bbox_inches="tight",
                metadata={"Title": "CRISPR discovery v2", "Subject": STATUS_LABEL})
    fig.savefig(stub.with_suffix(".svg"), facecolor="white", bbox_inches="tight",
                metadata={"Title": "CRISPR discovery v2", "Description": STATUS_LABEL})
    plt.close(fig)
    logger.info("wrote %s.png/.pdf/.svg (%d significant sensitising hits of %d fitted genes)",
                stub, n_hits, n_genomewide)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    build_crispr_discovery_main(OUT_DIR / "CRISPR_discovery_v2")
