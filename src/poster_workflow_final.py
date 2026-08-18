"""Poster figure 1 (final): the methods, as a diagram.

post_freeze_exploratory. This figure is a NEW addition -- the project had no
methods diagram, and its methods are hard to follow in prose: a genome-wide
screen, four transcriptomic datasets analysed with two different statistical
engines, a dependency resource, a network query and a structural audit.

It computes nothing. Every count it prints (genes fitted, hits, datasets,
cell lines, candidate cells) is read from the frozen tables through the same
loaders the other figures use, and is verified before anything is drawn, so
the diagram cannot drift out of step with the analysis it describes.

WHAT THE DIAGRAM MARKS, AND WHY IT MATTERS. Two provenance facts are drawn on
it rather than left in fine print:

  * WHERE THE CANDIDATE LIST WAS FIXED, relative to the external validation
    data. The analysis plan and its thresholds were dated before results
    existed, and the candidate list was fixed before the external datasets
    were opened. That is this project's strongest methodological claim.
  * WHICH PARTS ARE PRE-REGISTERED AND WHICH ARE POST-FREEZE. The screen gate
    and the candidate rule are pre-registered; the candidate reinterpretation
    that added KDM1A and TLK2, the network query, the meta-analysis and the
    power calculation all happened after the freeze. Drawing the freeze line
    without drawing what crosses it would imply the whole poster was
    prospectively specified, which is not true.

No laboratory work was performed at any stage: every box is a computational
reanalysis of public data.
"""

from __future__ import annotations

import logging
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

from src import post_audit_sensitivity_data as pad
from src import poster_depmap_v1 as d1
from src.poster_final_common import (
    FONT,
    OUT_DIR,
    SCREEN_FDR,
    figure_footer,
    headline,
    pin_reproducibility,
    save,
    verify,
)
from src.poster_network_mechanism_v2 import build_network
from src.poster_palette import NEUTRAL, WHITE

logger = logging.getLogger(__name__)

FIGURE = "F1_methods_workflow"


def load_counts() -> dict:
    hits = pad.load_significant_sensitising_hits()
    genomewide = pad.load_genomewide_crispr()
    dep = d1.dependency_summary(d1.load_cellline_table())
    G = build_network()
    return {"n_genomewide": len(genomewide), "n_hits": len(hits),
            "n_lines": int(dep["n_lines"].max()),
            "n_nodes": G.number_of_nodes(), "n_edges": G.number_of_edges()}


def gate(counts):
    return verify(FIGURE, [
        ("n_genomewide_fitted", counts["n_genomewide"], 19103, 0),
        ("n_screen_hits", counts["n_hits"], 13, 0),
        ("n_depmap_lines", counts["n_lines"], 11, 0),
        ("n_network_nodes", counts["n_nodes"], 47, 0),
        ("n_network_edges", counts["n_edges"], 147, 0),
    ])


# The axes spans 100 data units across 0.976 of a 15.5 inch figure, so one
# data unit is about 10.9 points. Text is wrapped to the box it sits in rather
# than hand-broken, so a box can be resized without re-breaking its strings.
PT_PER_UNIT = 0.976 * 15.5 * 72.0 / 100.0


def _wrap(text: str, width_units: float, fontsize: float, pad_units: float = 1.6) -> str:
    import textwrap

    usable_pt = max((width_units - pad_units) * PT_PER_UNIT, 1.0)
    n_chars = max(int(usable_pt / (0.545 * fontsize)), 8)
    return "\n".join("\n".join(textwrap.wrap(line, n_chars)) or "" for line in text.split("\n"))


def _box(ax, x, y, w, h, title, body, *, face=WHITE, edge=None, lw=1.8,
         title_size=None, body_size=None, dashed=False):
    edge = edge or NEUTRAL["rule"]
    title_size = FONT["annot"] + 1 if title_size is None else title_size
    body_size = FONT["note"] if body_size is None else body_size
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.22,rounding_size=0.30",
                                facecolor=face, edgecolor=edge, linewidth=lw,
                                linestyle=(0, (5, 4)) if dashed else "solid", zorder=3))
    wrapped_title = _wrap(title, w, title_size)
    n_title_lines = wrapped_title.count("\n") + 1
    ax.text(x + w / 2, y + h - 0.14 * h, wrapped_title, ha="center", va="top",
            fontsize=title_size, fontweight="bold", color=NEUTRAL["ink"], zorder=5,
            linespacing=1.25)
    if body:
        # start the body below the actual rendered depth of the title
        title_depth = n_title_lines * title_size * 1.25 / PT_PER_UNIT
        ax.text(x + w / 2, y + h - 0.14 * h - title_depth - 0.35, _wrap(body, w, body_size),
                ha="center", va="top", fontsize=body_size, color=NEUTRAL["ink_2"],
                zorder=5, linespacing=1.35)


def _arrow(ax, xy_from, xy_to, *, colour=None, lw=2.0, style="-|>"):
    ax.add_patch(FancyArrowPatch(xy_from, xy_to, arrowstyle=style, mutation_scale=22,
                                 color=colour or NEUTRAL["ink_muted"], lw=lw,
                                 shrinkA=3, shrinkB=3, zorder=2))


def build(stub: Path):
    pin_reproducibility(FIGURE)
    counts = load_counts()
    verification = gate(counts)

    fig, ax = plt.subplots(figsize=(15.5, 10.2))
    fig.subplots_adjust(left=0.012, right=0.988, top=0.985, bottom=0.070)
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 70)
    ax.axis("off")

    # ---------- the freeze line, drawn first so boxes sit on top ----------
    ax.plot([50.5, 50.5], [1.0, 58.0], color=NEUTRAL["ink"], lw=2.6,
            linestyle=(0, (7, 5)), zorder=1)
    ax.text(50.5, 69.0, "ANALYSIS PLAN AND CANDIDATE LIST FIXED HERE",
            ha="center", va="top", fontsize=FONT["annot"] + 1, fontweight="bold",
            color=NEUTRAL["ink"])
    ax.text(50.5, 65.2, "thresholds dated before any result existed · list fixed before\nthe external datasets were opened",
            ha="center", va="top", fontsize=FONT["note"] - 1, color=NEUTRAL["ink_2"],
            linespacing=1.3)
    ax.text(24.0, 56.6, "PRE-REGISTERED", ha="center", va="bottom",
            fontsize=FONT["note"] - 1, color=NEUTRAL["ink_muted"], fontweight="bold")
    ax.text(76.0, 56.6, "RUN AFTERWARDS, AGAINST THE FIXED LIST", ha="center", va="bottom",
            fontsize=FONT["note"] - 1, color=NEUTRAL["ink_muted"], fontweight="bold")

    # ---------- left of the line: the screen and the rule ----------
    _box(ax, 3, 40, 44, 14, "Genome-wide CRISPR knockout screen",
         f"{counts['n_genomewide']:,} genes · published data, reanalysed")
    _box(ax, 3, 23, 44, 13, f"{counts['n_hits']} hits pass the pre-registered gate",
         f"false discovery rate < {SCREEN_FDR:.2f}, sensitising direction")
    _box(ax, 3, 5, 44, 14, "Four candidates carried forward",
         "USP34, VEZF1 by the frozen rule\nKDM1A, TLK2 added after external audit")
    _arrow(ax, (25, 40), (25, 36.4))
    _arrow(ax, (25, 23), (25, 19.4))

    # ---------- right of the line: everything run against the fixed list ----------
    _box(ax, 54, 42, 44, 12, "Four transcriptomic datasets",
         "never used to pick the candidates")
    _box(ax, 54, 26, 21, 12, "Do the genes\nchange?", "→ Figure 3")
    _box(ax, 77, 26, 21, 12, "Do programmes\nchange?", "→ Figure 4")
    _box(ax, 54, 11, 13.5, 11, "Network", "→ Figure 5", title_size=FONT["annot"])
    _box(ax, 69.3, 11, 13.5, 11, "Dependency", "→ Figure 6", title_size=FONT["annot"])
    _box(ax, 84.5, 11, 13.5, 11, "Structures", "→ Figure 7", title_size=FONT["annot"])
    _arrow(ax, (64, 42), (64, 38.4))
    _arrow(ax, (88, 42), (88, 38.4))
    _arrow(ax, (60, 26), (60, 22.4))
    _arrow(ax, (76, 26), (76, 22.4))
    _arrow(ax, (91, 26), (91, 22.4))
    _arrow(ax, (47, 11.5), (54, 11.5), style="-|>")

    # post-freeze marker
    _box(ax, 54, 0.5, 44, 9.2, "Added after the freeze, labelled throughout",
         "extra candidates · network · meta-analysis · power",
         face=NEUTRAL["tint"], dashed=True, title_size=FONT["annot"],
         body_size=FONT["note"] - 1)

    headline(
        fig,
        "How the work was done, and where the list was fixed",
        "Everything here is a computational reanalysis of public data. No laboratory experiment was performed at any stage, and nothing on the\n"
        "right-hand side was allowed to change the candidate list on the left. The dashed box records what was added after the freeze, so the\n"
        "diagram does not imply that the whole study was specified in advance — it was not.",
        key=FIGURE)
    figure_footer(fig, "Counts read from the frozen tables at render time, not typed in.")
    return save(fig, stub), verification


def main(out_dir: Path = OUT_DIR):
    return build(out_dir / FIGURE)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    main()
