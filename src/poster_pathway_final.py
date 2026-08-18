"""Poster figure 4 (final): programme-level signal across the four datasets.

post_freeze_exploratory. No new enrichment: every normalised enrichment score
and false discovery rate is read at render time from the frozen Hallmark
gene-set tables via `poster_exploration_v2_data.load_pathway_trajectories()`,
the same loader the frozen poster figure 04 uses, and verified against the
frozen values before anything is drawn.

WHAT THIS FIXES. The frozen figure 04 labels its rows with raw gene-set
identifiers (HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION) that a non-specialist
cannot parse, and it gives the project's one genuinely novel observation no
more prominence than the rest. Here each row is labelled by the job the
programme does, with the formal gene-set name kept small underneath for
traceability, and the adhesion/motility row is the visual centre.

WHAT THE ADHESION/MOTILITY ROW SAYS, AND WHAT IT DOES NOT. Hallmark
"epithelial-mesenchymal transition" enrichment is positive in the three
long-term settings (chronic resistance models and recurrent tumours) and
negative after 12 hours of tamoxifen. That is a difference in a
transcriptional gene-set score between settings. It is NOT a measured
migration phenotype, and it is not evidence that acute response converts into
long-term biology: the acute dataset differs from the other three in far more
than duration (different tissue, different design, different statistical
unit). The figure states the observation and stops there.
"""

from __future__ import annotations

import logging
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm

from src import poster_exploration_v2_data as ed
from src.poster_final_common import (
    FONT,
    OUT_DIR,
    figure_footer,
    headline,
    pin_reproducibility,
    save,
    verify,
)
from src.poster_palette import DIVERGING, NEUTRAL, WHITE

logger = logging.getLogger(__name__)

FIGURE = "F4_programme_signal"
FDR_SIG = 0.05

# (formal gene-set name, plain-language job, is this the highlighted row)
ROWS: list[tuple[str, str, bool]] = [
    ("HALLMARK_ESTROGEN_RESPONSE_EARLY", "Oestrogen response (early)", False),
    ("HALLMARK_ESTROGEN_RESPONSE_LATE", "Oestrogen response (late)", False),
    ("HALLMARK_E2F_TARGETS", "Cell-cycle entry", False),
    ("HALLMARK_WNT_BETA_CATENIN_SIGNALING", "WNT / β-catenin signalling", False),
    ("HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "Cell adhesion and motility", True),
]
DATASETS = ["gse118713", "gse111151", "gse240112", "gse245601"]
DATASET_LABEL = {
    "gse118713": "Cell-line\nresistance",
    "gse111151": "Independent\nsublines",
    "gse240112": "Recurrent\ntumours",
    "gse245601": "Acute 12 h\ntamoxifen",
}
LONG_TERM = ("gse118713", "gse111151", "gse240112")

# Frozen reference values checked before plotting (label, dataset, NES).
REFERENCE = {
    ("Cell adhesion and motility", "gse118713"): 1.977988,
    ("Cell adhesion and motility", "gse111151"): 2.166115,
    ("Cell adhesion and motility", "gse240112"): 1.506975,
    ("Cell adhesion and motility", "gse245601"): -1.559312,
    ("Oestrogen response (early)", "gse118713"): -2.881807,
    ("Oestrogen response (early)", "gse245601"): -2.070304,
}


def load_rows():
    """Frozen Hallmark trajectories, relabelled but never re-valued."""
    spec = [("hallmark", formal, plain) for formal, plain, _ in ROWS]
    long = ed.load_pathway_trajectories(spec)
    logger.info("pathway: %d (programme, dataset) values loaded", len(long))
    return long


def gate(long):
    idx = long.set_index(["pathway_label", "dataset"])
    checks = [("n_values", len(long), len(ROWS) * len(DATASETS), 0)]
    for (label, ds), nes in REFERENCE.items():
        checks.append((f"NES[{label}|{ds}]", float(idx.loc[(label, ds), "NES"]), nes, 1e-6))
    # the highlighted observation, asserted as a directional fact
    emt = idx.xs("Cell adhesion and motility", level="pathway_label")
    for ds in LONG_TERM:
        checks.append((f"EMT_positive_in_{ds}", 1.0 if float(emt.loc[ds, "NES"]) > 0 else 0.0, 1.0, 0))
        checks.append((f"EMT_significant_in_{ds}", 1.0 if float(emt.loc[ds, "fdr"]) < FDR_SIG else 0.0, 1.0, 0))
    checks.append(("EMT_negative_in_acute", 1.0 if float(emt.loc["gse245601", "NES"]) < 0 else 0.0, 1.0, 0))
    checks.append(("EMT_significant_in_acute", 1.0 if float(emt.loc["gse245601", "fdr"]) < FDR_SIG else 0.0, 1.0, 0))
    # oestrogen response suppressed everywhere -- the internal positive control
    er = idx.xs("Oestrogen response (early)", level="pathway_label")
    checks.append(("ER_early_negative_in_all_four", float((er["NES"] < 0).sum()), 4.0, 0))
    return verify(FIGURE, checks)


def build(stub: Path):
    pin_reproducibility(FIGURE)
    long = load_rows()
    verification = gate(long)
    idx = long.set_index(["pathway_label", "dataset"])

    cmap = LinearSegmentedColormap.from_list("poster_div", DIVERGING)
    vmax = float(np.abs(long["NES"]).max())
    norm = TwoSlopeNorm(vcenter=0.0, vmin=-vmax, vmax=vmax)

    fig = plt.figure(figsize=(15.5, 9.4))
    ax = fig.add_axes([0.310, 0.235, 0.560, 0.665])
    n_rows, n_cols = len(ROWS), len(DATASETS)

    for r, (formal, plain, highlight) in enumerate(ROWS):
        y = n_rows - 1 - r
        for c, ds in enumerate(DATASETS):
            nes = float(idx.loc[(plain, ds), "NES"])
            fdr = float(idx.loc[(plain, ds), "fdr"])
            ax.add_patch(plt.Rectangle((c - 0.46, y - 0.42), 0.92, 0.84,
                                       facecolor=cmap(norm(nes)), edgecolor=WHITE,
                                       linewidth=2.0, zorder=2))
            # significance by fill, exactly as everywhere else: a filled cell
            # is significant, a hollow cell is not
            if fdr >= FDR_SIG:
                ax.add_patch(plt.Rectangle((c - 0.46, y - 0.42), 0.92, 0.84,
                                           facecolor=WHITE, edgecolor=NEUTRAL["rule"],
                                           linewidth=2.0, zorder=3))
                ax.text(c, y, "n.s.", ha="center", va="center", fontsize=FONT["annot"],
                        color=NEUTRAL["ink_muted"], zorder=4)
            else:
                light = abs(nes) > 0.55 * vmax
                ax.text(c, y, f"{nes:+.2f}", ha="center", va="center", fontsize=FONT["annot"] + 2,
                        fontweight="bold" if highlight else "normal",
                        color=WHITE if light else NEUTRAL["ink"], zorder=4)
        # row labels: the job first, the formal name small beneath
        ax.text(-0.62, y + 0.10, plain, ha="right", va="center", fontsize=FONT["panel"] - 3,
                fontweight="bold" if highlight else "normal", color=NEUTRAL["ink"])
        ax.text(-0.62, y - 0.20, formal.replace("HALLMARK_", "").replace("_", " ").lower(),
                ha="right", va="center", fontsize=FONT["note"] - 2, color=NEUTRAL["ink_muted"], style="italic")
        if highlight:
            ax.add_patch(plt.Rectangle((-0.52, y - 0.48), n_cols - 0.04 + 0.06, 0.96,
                                       facecolor="none", edgecolor=NEUTRAL["ink"],
                                       linewidth=2.4, zorder=6))

    ax.set_xlim(-0.55, n_cols - 0.45)
    ax.set_ylim(-0.6, n_rows - 0.10)
    ax.set_xticks(range(n_cols))
    ax.set_xticklabels([DATASET_LABEL[d] for d in DATASETS], fontsize=FONT["tick"], color=NEUTRAL["ink"])
    ax.tick_params(axis="x", length=0, pad=10)
    ax.set_yticks([])
    for side in ax.spines:
        ax.spines[side].set_visible(False)
    ax.xaxis.set_ticks_position("top")
    ax.xaxis.set_label_position("top")

    # a rule separating the three long-term settings from the acute one
    ax.vlines(2.5, -0.52, n_rows - 0.42, color=NEUTRAL["ink"], linewidth=1.6, zorder=7)
    ax.text(1.0, n_rows - 0.30, "three long-term settings", ha="center", va="bottom",
            fontsize=FONT["note"], color=NEUTRAL["ink_2"])
    ax.text(3.0, n_rows - 0.30, "acute", ha="center", va="bottom",
            fontsize=FONT["note"], color=NEUTRAL["ink_2"])

    # colour bar
    cax = fig.add_axes([0.310, 0.118, 0.215, 0.020])
    sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    cb = fig.colorbar(sm, cax=cax, orientation="horizontal")
    cb.set_label("normalised enrichment score  (− suppressed, + enriched)",
                 fontsize=FONT["note"], color=NEUTRAL["ink_2"], labelpad=8)
    cb.ax.tick_params(labelsize=FONT["note"] - 2, colors=NEUTRAL["ink_2"], length=0)
    cb.outline.set_visible(False)
    fig.text(0.580, 0.128, "“n.s.” = not significant at FDR 0.05",
             fontsize=FONT["note"], color=NEUTRAL["ink_2"], va="center")

    headline(
        fig,
        "Programme-level signal is present in all four datasets, and the\nadhesion/motility programme points the other way after 12 hours",
        "Oestrogen response is suppressed in all four datasets and cell-cycle entry in all four — the expected pharmacology, and an internal\n"
        "check that the datasets behave sensibly. The boxed row is the observation this project adds: Hallmark epithelial-to-mesenchymal\n"
        "transition enrichment is positive in the three long-term settings and negative after 12 hours of tamoxifen. This is a gene-set score,\n"
        "not a measured migration phenotype, and the acute dataset differs from the other three in more than duration.",
        key=FIGURE)
    figure_footer(fig, "Frozen Hallmark gene-set enrichment tables; no new enrichment run.")
    return save(fig, stub), verification


def main(out_dir: Path = OUT_DIR):
    return build(out_dir / FIGURE)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    main()
