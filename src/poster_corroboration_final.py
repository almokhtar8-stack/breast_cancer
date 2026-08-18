"""Poster figure 3 (final): candidate-level corroboration across four datasets.

post_freeze_exploratory. Replaces the frozen figure 02 within-dataset z-score
heatmap, whose colour standardisation stretched noise across the full colour
range -- a reader misread genes at false discovery rate 0.49 as strongly
changed. Significance is now on an axis, where it cannot be misread.

The volcano panels reuse, unchanged, the loaders and the verification gate of
`src.poster_candidate_volcano_v2` (which in turn reuses v1's): the same four
frozen genome-wide differential-expression tables, the same 16 reference
values checked before anything is drawn, and the same zero-displacement rule
-- coincident candidates are drawn as concentric rings at their true shared
position, never displaced. Only presentation changes here: the project
palette, poster-scale type, and a fourth row.

WHY THE FOURTH ROW EXISTS. The three volcano rows alone would make "the
candidates do not corroborate" sound more decisive than the evidence permits.
The pooled-evidence strip carries the two post-freeze analyses that qualify
that statement: a random-effects meta-analysis (no candidate reaches pooled
FDR 0.05 under the primary estimator; smallest 0.100) and a minimum-detectable
-effect calculation showing that GSE111151's nulls are uninformative rather
than negative. Both are read from their committed tables, not retyped.

Wording is deliberate. These four datasets are not four equivalent
replications -- they are a chronic resistance model, an independent resistant
subline panel, an unpaired recurrent-versus-primary tissue comparison and an
acute 12-hour exposure. The figure says "corroboration", never "replication",
and never calls the acute dataset a resistance measurement.
"""

from __future__ import annotations

import logging
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from src.poster_candidate_volcano_v1 import CANDIDATES, load_panels, verify_against_frozen
from src.poster_candidate_volcano_v2 import (
    CONCENTRIC_SIZES,
    ZOOM_XLIM,
    ZOOM_YLIM,
    candidate_points,
    coincident_clusters,
)
from src.poster_final_common import (
    FONT,
    OUT_DIR,
    SIG_FDR,
    figure_footer,
    headline,
    pin_reproducibility,
    save,
    significance_marker,
    style_axes,
    verify,
)
from src.poster_palette import GENE_COLOURS, NEUTRAL, WHITE

logger = logging.getLogger(__name__)

FIGURE = "F3_candidate_corroboration"

PANEL_TITLE = {
    "GSE118713": "Cell-line\nresistance model",
    "GSE111151": "Independent\nresistant sublines",
    "GSE240112": "Recurrent vs\nprimary tumours",
    "GSE245601": "Acute 12 h\ntamoxifen",
}
PANEL_SUB = {
    "GSE118713": "GSE118713 · bulk",
    "GSE111151": "GSE111151 · bulk",
    "GSE240112": "GSE240112 · unpaired",
    "GSE245601": "GSE245601 · not resistance",
}

META_TSV = Path("results/post_poster/meta_analysis/candidates_meta_analysis.tsv")
POWER_TSV = Path("results/post_poster/power/dataset_sensitivity_summary.tsv")

# Fixed label anchors in data coordinates; no collision solver anywhere.
LABEL_POS: dict[tuple[str, str], tuple[float, float, str]] = {
    ("GSE118713", "USP34"): (0.30, 2.13, "right"), ("GSE118713", "VEZF1"): (1.05, 0.85, "left"),
    ("GSE118713", "KDM1A"): (0.90, 0.40, "left"), ("GSE118713", "TLK2"): (-0.70, 0.52, "right"),
    ("GSE111151", "USP34"): (0.95, 0.72, "left"), ("GSE111151", "VEZF1"): (-0.85, 0.72, "right"),
    ("GSE111151", "KDM1A"): (-0.80, 0.30, "right"), ("GSE111151", "TLK2"): (0.95, 0.30, "left"),
    ("GSE240112", "USP34"): (0.95, 0.98, "left"), ("GSE240112", "VEZF1"): (0.85, 2.18, "right"),
    ("GSE240112", "KDM1A"): (0.90, 0.40, "left"), ("GSE240112", "TLK2"): (-0.70, 0.45, "right"),
    ("GSE245601", "USP34"): (0.80, 0.62, "left"), ("GSE245601", "VEZF1"): (0.80, 0.24, "left"),
    ("GSE245601", "KDM1A"): (-0.85, 0.88, "right"), ("GSE245601", "TLK2"): (-0.85, 0.36, "right"),
}
DRAW_ORDER = ("KDM1A", "TLK2", "USP34", "VEZF1")


def load_pooled_evidence():
    """The two post-freeze qualifiers, read from their committed tables."""
    meta = pd.read_csv(META_TSV, sep="\t")
    primary = meta[(meta["arm"] == "all3") & (meta["se_variant"] == "wald")]
    best = primary.nsmallest(1, "fdr").iloc[0]
    power = pd.read_csv(POWER_TSV, sep="\t")
    p111 = power[(power["dataset"] == "GSE111151") & (power["alpha_name"] == "nominal_0.05")].iloc[0]
    return {"meta_min_fdr": float(best["fdr"]), "meta_min_gene": str(best["gene_symbol"]),
            "meta_effect": float(best["pooled_effect"]),
            "meta_n_below_005": int((primary["fdr"] < SIG_FDR).sum()),
            "mde_median": float(p111["median_mde80"]),
            "mde_n_tested": int(p111["n_candidates_tested"])}


def gate(panels, pooled):
    """Volcano gate (16 frozen values, unchanged) plus the pooled-evidence
    values this figure prints."""
    verification = verify_against_frozen(panels)      # raises on any mismatch
    n_sig = sum(int(float(p.candidate_rows.loc[g, "fdr"]) < SIG_FDR)
                for p in panels for g in CANDIDATES)
    extra = verify(FIGURE, [
        ("n_significant_of_16", n_sig, 2, 0),
        ("n_candidate_dataset_cells", len(panels) * len(CANDIDATES), 16, 0),
        ("meta_primary_min_fdr", pooled["meta_min_fdr"], 0.100257, 1e-5),
        ("meta_n_candidates_below_fdr_005", pooled["meta_n_below_005"], 0, 0),
        ("gse111151_median_mde80", pooled["mde_median"], 0.626918, 1e-5),
        ("gse111151_n_candidates_tested", pooled["mde_n_tested"], 12, 0),
    ])
    return pd.concat([verification.assign(figure=FIGURE), extra], ignore_index=True)


def build(stub: Path):
    pin_reproducibility(FIGURE)
    panels = load_panels()
    pooled = load_pooled_evidence()
    verification = gate(panels, pooled)

    fig = plt.figure(figsize=(15.5, 11.0))
    gs = fig.add_gridspec(2, 4, height_ratios=[1.0, 0.46], left=0.062, right=0.986,
                          top=0.935, bottom=0.075, wspace=0.09, hspace=0.30)
    axes = [fig.add_subplot(gs[0, i]) for i in range(4)]
    strip = fig.add_subplot(gs[1, :])
    sig_y = -np.log10(SIG_FDR)
    n_filled = 0

    for i, (ax, panel) in enumerate(zip(axes, panels)):
        lfc = panel.df["log2fc"].to_numpy()
        y = -np.log10(panel.df["fdr"].to_numpy())
        cand = panel.df["gene"].isin(CANDIDATES).to_numpy()
        style_axes(ax)
        ax.set_xlim(*ZOOM_XLIM)
        ax.set_ylim(*ZOOM_YLIM)
        ax.set_xticks((-1.0, -0.5, 0.0, 0.5, 1.0))
        inz = (lfc >= ZOOM_XLIM[0]) & (lfc <= ZOOM_XLIM[1]) & (y <= ZOOM_YLIM[1]) & ~cand
        ax.scatter(lfc[inz], y[inz], s=5, c=NEUTRAL["backdrop"], alpha=0.35,
                   linewidths=0, zorder=1)
        ax.axhline(sig_y, color=NEUTRAL["rule"], linestyle="--", linewidth=1.3, zorder=3)
        ax.axvline(0.0, color=NEUTRAL["rule"], linewidth=0.8, zorder=3)

        measured = candidate_points(panel)
        clusters = coincident_clusters(measured)
        for z, gene in enumerate(DRAW_ORDER):
            gx, gy = measured[gene]                     # TRUE position; zero displacement
            fdr = float(panel.candidate_rows.loc[gene, "fdr"])
            passes = fdr < SIG_FDR
            members = clusters[gene]
            rank = members.index(gene)
            size = 240.0 if len(members) == 1 else CONCENTRIC_SIZES[rank] * 1.25
            zo = 10 + z if len(members) == 1 else 10 + (len(members) - rank)
            significance_marker(ax, gx, gy, gene, passes, size=size, lw=3.0,
                                zorder=zo, clip_on=False)
            if passes:
                n_filled += 1
            lx, ly, ha = LABEL_POS[(panel.accession, gene)]
            txt = f"{gene}\nFDR {fdr:.4f}" if passes else gene
            ax.annotate(txt, (gx, gy), xytext=(lx, ly), ha=ha, va="center",
                        fontsize=FONT["annot"], color=GENE_COLOURS[gene] if passes else NEUTRAL["ink"],
                        fontweight="bold" if passes else "normal", zorder=20,
                        arrowprops=dict(arrowstyle="-", color=NEUTRAL["rule"], lw=0.9,
                                        shrinkA=2, shrinkB=7))
        ax.text(0.5, 1.105, PANEL_TITLE[panel.accession], transform=ax.transAxes,
                ha="center", va="bottom", fontsize=FONT["panel"] - 3, fontweight="bold",
                color=NEUTRAL["ink"], linespacing=1.2)
        ax.text(0.5, 1.022, PANEL_SUB[panel.accession], transform=ax.transAxes,
                ha="center", va="bottom", fontsize=FONT["note"], color=NEUTRAL["ink_muted"])
        ax.set_xlabel("log$_2$ fold change", fontsize=FONT["axis"] - 2, color=NEUTRAL["ink_2"])
        ax.tick_params(labelsize=FONT["tick"])
        if i > 0:
            ax.tick_params(labelleft=False)
    if n_filled != 2:
        raise ValueError(f"drew {n_filled} filled candidate points, expected 2")

    axes[0].set_ylabel("$-\\log_{10}$(false discovery rate)", fontsize=FONT["axis"] - 2, color=NEUTRAL["ink_2"])
    axes[0].text(ZOOM_XLIM[0] + 0.04, sig_y + 0.05, "FDR 0.05", fontsize=FONT["note"],
                 color=NEUTRAL["ink_muted"], va="bottom")

    # --- pooled-evidence strip -------------------------------------------------
    strip.axis("off")
    strip.add_patch(plt.Rectangle((0, 0), 1, 1, transform=strip.transAxes,
                                  facecolor=NEUTRAL["tint"], edgecolor="none", zorder=0))
    strip.text(0.014, 0.90, "Why these nulls are weak evidence, not strong negatives",
               fontsize=FONT["panel"], fontweight="bold", color=NEUTRAL["ink"], va="top",
               transform=strip.transAxes)
    strip.text(0.014, 0.63,
               f"Pooling the three resistance-context datasets (random-effects meta-analysis), {pooled['meta_n_below_005']} of 13 screen hits reach pooled FDR 0.05;\n"
               f"the smallest is {pooled['meta_min_fdr']:.3f} ({pooled['meta_min_gene']}). GSE111151 — the only dataset with several independent resistance backgrounds —\n"
               f"could detect a {2 ** pooled['mde_median']:.1f}-fold change at 80% power, and all {pooled['mde_n_tested']} of its testable candidate effects were smaller than that.",
               fontsize=FONT["note"], color=NEUTRAL["ink_2"], va="top", transform=strip.transAxes,
               linespacing=1.55)
    strip.text(0.014, 0.07,
               "Its nulls are uninformative, not negative. The power argument applies to GSE111151; the other nulls are weak for other reasons.",
               fontsize=FONT["note"], color=NEUTRAL["ink_muted"], va="top", transform=strip.transAxes, style="italic")

    headline(
        fig,
        "Two of sixteen candidate gene-and-dataset combinations reach\nfalse discovery rate 0.05",
        "Filled = reaches FDR 0.05; hollow ring = does not. Each candidate is corroborated in at most one dataset. Grey = all\n"
        "other genes in that dataset (Benjamini–Hochberg across the transcriptome, within each dataset). These are not four\n"
        "equivalent replications: panel 4 is an acute 12-hour response, not resistance; panel 3 compares different, unpaired\n"
        "patients from two tissue banks. Coincident candidates are concentric rings at their true position — nothing is displaced.",
        key=FIGURE)
    figure_footer(fig)
    return save(fig, stub), verification


def main(out_dir: Path = OUT_DIR):
    return build(out_dir / FIGURE)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    main()
