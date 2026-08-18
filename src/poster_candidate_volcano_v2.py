"""Candidate volcano figure, v2 -- proposed Figure 2 replacement, two variants.

post_freeze_exploratory. v2 changes LAYOUT and COLOUR only. It reuses v1's data
loaders and verification gate unchanged (`src.poster_candidate_volcano_v1`), so
every plotted value is still read from a frozen genome-wide DE table, is still
verified against the frozen candidate reference numbers before anything is
plotted, and the significant count is still asserted to be exactly 2. v1 is
retained alongside for provenance and is not imported for colours.

What changed from v1 (see results/reports/post_poster/FIGURE2_VOLCANO_V2_NOTE.md):

  * Gene colours are the project's data-tier palette, defined ONCE here in
    ``GENE_COLOURS`` so later passes over other figures can import the same
    source of truth. Everything else uses the repository's neutral ladder.
  * Shared axes: one common x range and one common y range across all four
    panels of a variant, and identical panel geometry, so the same fold
    change sits in the same place in every panel.
  * A shared zoom row beneath the four main panels (variant A) covers the
    region where the candidates actually sit; a rectangle on each main panel
    marks it.
  * Coincident candidates: in variant A's zoom row they are spread
    horizontally by a fixed, documented amount (``COINCIDENCE_TOL``,
    ``COINCIDENCE_STEP``, disclosed on the figure); y is never moved. In
    variant B NO point is displaced: coincident candidates are drawn as
    concentric rings at their true shared position, so every plotted x
    equals the source log2FC exactly (test-enforced).
  * Labels sit at fixed, hand-placed positions with leader lines -- no
    collision solver -- so output is reproducible.
  * Two variants from the same verified data in one run: A with the
    genome-wide backdrop, B candidates only.
  * A deuteranopia / protanopia simulation of the palette and of the
    rendered PNGs (Machado, Oliveira & Fernandes 2009 matrices, computed
    directly -- no CVD library is installed and none is added).

Data sources: identical to v1 (config/config.yaml
``cross_dataset_genomewide.inputs``: unredacted GSE118713 TAMR_vs_MCF7 limma
table; GSE111151 edgeR; GSE240112 tumour-cell-track edgeR pseudobulk;
GSE245601 Track A edgeR pseudobulk). No statistic is recomputed.

Determinism: no network, no randomness. PDF bytes reproducible via
``SOURCE_DATE_EPOCH``; SVG ids pinned via ``svg.hashsalt``.
"""

from __future__ import annotations

import hashlib
import itertools
import logging
import os
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle

from src.poster_candidate_volcano_v1 import (  # loaders + gate, unchanged
    CANDIDATES,
    EXPECTED_SIGNIFICANT,
    PANEL_ORDER,
    SIG_FDR,
    STATUS_LABEL,
    Panel,
    load_panels,
    verify_against_frozen,
)

logger = logging.getLogger(__name__)

OUT_DIR = Path("results/figures/poster_candidate_volcano_v2")

# ---------------------------------------------------------------------------
# Colour: the project's data-tier gene colours (single source of truth) and
# the repository's neutral ladder. Nothing else is used inside a plot.
# ---------------------------------------------------------------------------
GENE_COLOURS: dict[str, str] = {
    "USP34": "#6A3D9A",
    "KDM1A": "#D55E00",
    "TLK2": "#009E73",
    "VEZF1": "#56B4E9",
}
NEUTRAL: dict[str, str] = {
    "text": "#262626",        # figure text
    "subtitle": "#555555",    # subtitles and axis names
    "tick": "#8c8c8c",        # tick labels
    "line": "#b0b0b0",        # threshold and reference lines
    "background": "#c9c9c9",  # genome-wide background points
    "grid": "#e6e6e6",        # gridlines
}

# Two-line panel titles (so all four panels have identical header height),
# biological description first; the accession is drawn separately beneath.
PANEL_TITLES: dict[str, str] = {
    "GSE118713": "Cell-line resistance\nmodel",
    "GSE111151": "Independent tamoxifen-\nresistant sublines",
    "GSE240112": "Primary vs recurrent\ntumours",
    "GSE245601": "Acute 12 h tamoxifen\nper-tumour pseudobulk",
}

CAPTION = ("2 of 16 candidate gene-dataset combinations reach FDR 0.05. "
           "Each candidate is supported by at most one dataset.")
TITLE = "Candidate significance across the four transcriptomic contexts"

# ---------------------------------------------------------------------------
# Shared axis ranges. Variant A must contain every plotted gene in every
# panel; the zoom region (and variant B) must contain every candidate.
# Values are derived from the data at run time and asserted against these
# declared bounds so nothing is ever clipped silently.
# ---------------------------------------------------------------------------
ZOOM_XLIM: tuple[float, float] = (-1.5, 1.5)   # log2FC
ZOOM_YLIM: tuple[float, float] = (0.0, 2.5)    # -log10(FDR)
VARIANT_B_XLIM = ZOOM_XLIM
VARIANT_B_YLIM = ZOOM_YLIM
RANGE_PAD = 1.06  # variant A limits = data extreme x pad

# Coincidence rule: candidates whose pairwise |dx| and |dy| are both below
# COINCIDENCE_TOL (in the zoom / variant-B data units) are spread
# horizontally, in fixed alphabetical order, at COINCIDENCE_STEP intervals
# centred on the cluster's mean x. y is never changed.
COINCIDENCE_TOL: float = 0.06
COINCIDENCE_STEP: float = 0.14

# Fixed label anchors in zoom / variant-B data coordinates: (x, y, ha).
# Hand-placed once against the frozen values; a leader line joins each label
# to its point. No collision solver is used anywhere.
LABEL_POS: dict[tuple[str, str], tuple[float, float, str]] = {
    ("GSE118713", "USP34"): (0.28, 2.14, "right"), ("GSE118713", "VEZF1"): (0.95, 0.90, "left"),
    ("GSE118713", "KDM1A"): (0.80, 0.42, "left"), ("GSE118713", "TLK2"): (-0.60, 0.55, "right"),
    ("GSE111151", "USP34"): (0.80, 0.75, "left"), ("GSE111151", "VEZF1"): (-0.75, 0.75, "right"),
    ("GSE111151", "KDM1A"): (-0.65, 0.28, "right"), ("GSE111151", "TLK2"): (0.80, 0.32, "left"),
    ("GSE240112", "USP34"): (0.90, 0.95, "left"), ("GSE240112", "VEZF1"): (0.90, 2.20, "right"),
    ("GSE240112", "KDM1A"): (0.80, 0.42, "left"), ("GSE240112", "TLK2"): (-0.60, 0.45, "right"),
    ("GSE245601", "USP34"): (0.65, 0.60, "left"), ("GSE245601", "VEZF1"): (0.65, 0.22, "left"),
    ("GSE245601", "KDM1A"): (-0.75, 0.85, "right"), ("GSE245601", "TLK2"): (-0.75, 0.32, "right"),
}

# Explicit draw order (later = on top) so no ring is ever fully hidden.
DRAW_ORDER: tuple[str, ...] = ("KDM1A", "TLK2", "USP34", "VEZF1")

# Variant B: coincident candidates share their TRUE centre and differ in
# radius. Marker areas (points^2) for the 1st, 2nd, 3rd member of a cluster
# in alphabetical order; the largest is drawn first so none is hidden.
CONCENTRIC_SIZES: tuple[float, ...] = (190.0, 460.0, 820.0)

OFFSET_DISCLOSURE = ("Zoom row: candidates within 0.06 log2FC of each other are drawn 0.14 apart "
                     "horizontally so each can be counted (5 points; y unchanged; measured values in the manifest).")


# ---------------------------------------------------------------------------
# Geometry helpers
# ---------------------------------------------------------------------------
def candidate_points(panel: Panel) -> dict[str, tuple[float, float]]:
    """(x, y) = (log2FC, -log10 FDR) per candidate, as measured."""
    return {g: (float(panel.candidate_rows.loc[g, "log2fc"]),
                float(-np.log10(panel.candidate_rows.loc[g, "fdr"]))) for g in CANDIDATES}


def spread_coincident(points: dict[str, tuple[float, float]]) -> tuple[dict[str, tuple[float, float]], list[dict]]:
    """Apply the coincidence rule. Returns (display points, offset records).
    Only x moves, by a fixed step in fixed alphabetical order; y is untouched
    so significance side is preserved by construction."""
    genes = sorted(points)
    # union-find clusters on the pairwise tolerance
    parent = {g: g for g in genes}

    def find(g):
        while parent[g] != g:
            g = parent[g]
        return g

    for a, b in itertools.combinations(genes, 2):
        if abs(points[a][0] - points[b][0]) < COINCIDENCE_TOL and abs(points[a][1] - points[b][1]) < COINCIDENCE_TOL:
            parent[find(a)] = find(b)
    clusters: dict[str, list[str]] = {}
    for g in genes:
        clusters.setdefault(find(g), []).append(g)

    display = dict(points)
    records: list[dict] = []
    for members in clusters.values():
        if len(members) < 2:
            continue
        members = sorted(members)  # fixed alphabetical order
        mean_x = float(np.mean([points[g][0] for g in members]))
        for i, g in enumerate(members):
            new_x = mean_x + (i - (len(members) - 1) / 2.0) * COINCIDENCE_STEP
            display[g] = (new_x, points[g][1])
            records.append({"gene": g, "measured_x": points[g][0], "display_x": new_x,
                            "x_offset": new_x - points[g][0], "y_unchanged": points[g][1],
                            "cluster": "+".join(members)})
    return display, records


def variant_a_limits(panels: list[Panel]) -> tuple[tuple[float, float], tuple[float, float]]:
    xmax = max(float(np.abs(p.df["log2fc"]).max()) for p in panels) * RANGE_PAD
    ymax = max(float((-np.log10(p.df["fdr"])).max()) for p in panels) * RANGE_PAD
    return (-xmax, xmax), (0.0, ymax)


def assert_candidates_inside(panels: list[Panel], xlim, ylim, where: str) -> None:
    for p in panels:
        for g, (x, y) in candidate_points(p).items():
            if not (xlim[0] <= x <= xlim[1] and ylim[0] <= y <= ylim[1]):
                raise ValueError(f"{where}: candidate {g}/{p.accession} at ({x:.3f},{y:.3f}) "
                                 f"is outside x{xlim} y{ylim}; refusing to clip a candidate")


# ---------------------------------------------------------------------------
# Drawing primitives
# ---------------------------------------------------------------------------
ZOOM_XTICKS = (-1.0, -0.5, 0.0, 0.5, 1.0)


def _style_axes(ax, xlim, ylim, sig_y):
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    if tuple(xlim) == ZOOM_XLIM:
        # edge labels of neighbouring panels would otherwise touch ("1.5-1.5")
        ax.set_xticks(ZOOM_XTICKS)
    ax.axhline(sig_y, color=NEUTRAL["line"], linestyle="--", linewidth=1.0, zorder=3)
    ax.axvline(0.0, color=NEUTRAL["line"], linestyle="-", linewidth=0.7, zorder=3)
    ax.tick_params(labelsize=8.5, colors=NEUTRAL["tick"], length=3)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    for spine in ("left", "bottom"):
        ax.spines[spine].set_color(NEUTRAL["tick"])


def coincident_clusters(points: dict[str, tuple[float, float]]) -> dict[str, list[str]]:
    """gene -> alphabetically sorted members of its coincidence cluster
    (singletons included), using the same tolerance as spread_coincident."""
    _, records = spread_coincident(points)
    members: dict[str, list[str]] = {g: [g] for g in points}
    for r in records:
        members[r["gene"]] = r["cluster"].split("+")
    return members


def _draw_candidates(ax, panel: Panel, ring_size: float, fill_size: float, lw: float,
                     labels: bool, coincidence: str, label_fontsize: float) -> tuple[int, list[dict]]:
    """Draw the four candidates on ``ax`` in fixed order.

    ``coincidence`` is one of
      "measured"   -- every point at its measured position, plain rings;
      "spread"     -- coincident points displaced horizontally (variant A zoom row);
      "concentric" -- coincident points at their TRUE position as concentric
                      rings of increasing radius (variant B); zero displacement.
    Returns (n_filled, offset records) -- records are empty unless "spread"."""
    if coincidence not in ("measured", "spread", "concentric"):
        raise ValueError(coincidence)
    measured = candidate_points(panel)
    display, records = spread_coincident(measured) if coincidence == "spread" else (dict(measured), [])
    clusters = coincident_clusters(measured) if coincidence == "concentric" else {g: [g] for g in measured}
    n_filled = 0
    for z, gene in enumerate(DRAW_ORDER):
        x, y = display[gene]
        fdr = float(panel.candidate_rows.loc[gene, "fdr"])
        significant = fdr < SIG_FDR
        colour = GENE_COLOURS[gene]
        members = clusters[gene]
        rank = members.index(gene)              # 0 for singletons
        size = ring_size if len(members) == 1 else CONCENTRIC_SIZES[rank]
        # concentric: larger rings drawn beneath smaller ones so none is hidden
        zorder = 10 + z if len(members) == 1 else 10 + (len(members) - rank)
        if significant:
            ax.scatter([x], [y], s=fill_size, c=colour, edgecolors="white",
                       linewidths=1.2, zorder=10 + z)
            n_filled += 1
        else:
            # clip_on=False: a ring centred just above y=0 must not lose its
            # lower arc to the axis edge, or it stops being countable
            ax.scatter([x], [y], s=size, facecolors="none", edgecolors=colour,
                       linewidths=lw, zorder=zorder, clip_on=False)
        if labels:
            lx, ly, ha = LABEL_POS[(panel.accession, gene)]
            text = f"{gene}\nFDR {fdr:.4f}" if significant else gene
            ax.annotate(text, (x, y), xytext=(lx, ly), ha=ha, va="center",
                        fontsize=label_fontsize, color=colour if significant else NEUTRAL["text"],
                        fontweight="bold" if significant else "normal", zorder=20,
                        arrowprops=dict(arrowstyle="-", color=NEUTRAL["line"], lw=0.7,
                                        shrinkA=2, shrinkB=5))
    return n_filled, records


def _panel_header(ax, panel: Panel):
    ax.text(0.5, 1.115, PANEL_TITLES[panel.accession], transform=ax.transAxes, ha="center",
            va="bottom", fontsize=10.8, color=NEUTRAL["text"], fontweight="bold", linespacing=1.15)
    ax.text(0.5, 1.035, panel.accession, transform=ax.transAxes, ha="center", va="bottom",
            fontsize=9.0, color=NEUTRAL["tick"])


def _legend_strip(fig, y: float, with_background: bool, disclosure: str | None = None):
    handles = [
        Line2D([], [], marker="o", linestyle="none", markersize=9, markerfacecolor="none",
               markeredgecolor=NEUTRAL["text"], markeredgewidth=1.8, label="candidate, FDR ≥ 0.05"),
        Line2D([], [], marker="o", linestyle="none", markersize=9.5, markerfacecolor=NEUTRAL["text"],
               markeredgecolor="white", label="candidate, FDR < 0.05"),
        Line2D([], [], linestyle="--", color=NEUTRAL["line"], linewidth=1.2, label="FDR 0.05"),
    ]
    if with_background:
        handles.append(Line2D([], [], marker="o", linestyle="none", markersize=4.5,
                              markerfacecolor=NEUTRAL["background"], markeredgecolor="none",
                              label="all genes in the dataset"))
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, y), ncol=len(handles),
               frameon=False, fontsize=9.0, handletextpad=0.5, columnspacing=1.8,
               labelcolor=NEUTRAL["subtitle"])
    if disclosure:
        # an alteration of what is drawn belongs ON the figure, not only in the note
        fig.text(0.5, y - 0.006, disclosure, ha="center", va="top", fontsize=8.2,
                 color=NEUTRAL["subtitle"])


def _figure_text(fig, top_y: float, caption_gap: float = 0.045):
    fig.suptitle(TITLE, x=0.045, y=top_y, ha="left", fontsize=17, fontweight="bold",
                 color=NEUTRAL["text"])
    fig.text(0.045, top_y - caption_gap, CAPTION, fontsize=10.2, color=NEUTRAL["subtitle"], va="top")
    fig.text(0.985, 0.008, STATUS_LABEL + " — candidate figure, not the frozen poster Figure 2",
             ha="right", va="bottom", fontsize=7.0, color=NEUTRAL["tick"], style="italic")


def _pin_reproducibility():
    os.environ.setdefault("SOURCE_DATE_EPOCH", "1600000000")
    plt.rcParams["svg.hashsalt"] = "poster_candidate_volcano_v2"


def _save(fig, stub: Path) -> dict[str, Path]:
    stub.parent.mkdir(parents=True, exist_ok=True)
    meta_common = {"Title": "Candidate volcano v2 (proposed Figure 2 replacement)"}
    written: dict[str, Path] = {}
    for ext, meta in (("png", {**meta_common, "Description": STATUS_LABEL}),
                      ("pdf", {**meta_common, "Subject": STATUS_LABEL, "Keywords": STATUS_LABEL}),
                      ("svg", {**meta_common, "Description": STATUS_LABEL})):
        path = stub.with_suffix(f".{ext}")
        fig.savefig(path, facecolor="white", dpi=300, metadata=meta)
        written[ext] = path
    plt.close(fig)
    return written


# ---------------------------------------------------------------------------
# Variant A: genome-wide backdrop + shared zoom row
# ---------------------------------------------------------------------------
def build_variant_a(panels: list[Panel], stub: Path) -> tuple[dict[str, Path], dict]:
    _pin_reproducibility()
    xlim, ylim = variant_a_limits(panels)
    assert_candidates_inside(panels, xlim, ylim, "variant A main")
    assert_candidates_inside(panels, ZOOM_XLIM, ZOOM_YLIM, "variant A zoom")
    sig_y = -np.log10(SIG_FDR)

    fig = plt.figure(figsize=(13.3, 10.2))
    gs = fig.add_gridspec(2, 4, height_ratios=[2.0, 1.15], left=0.055, right=0.985,
                          top=0.80, bottom=0.135, wspace=0.10, hspace=0.38)
    main_axes = [fig.add_subplot(gs[0, i]) for i in range(4)]
    zoom_axes = [fig.add_subplot(gs[1, i]) for i in range(4)]

    n_filled_main = 0
    n_filled_zoom = 0
    offset_records: list[dict] = []
    for i, (ax, zax, panel) in enumerate(zip(main_axes, zoom_axes, panels)):
        lfc = panel.df["log2fc"].to_numpy()
        y = -np.log10(panel.df["fdr"].to_numpy())
        cand_mask = panel.df["gene"].isin(CANDIDATES).to_numpy()
        # main panel: genome-wide backdrop, then candidates (no labels here --
        # at this scale labels would sit on the clouds; the zoom row carries them)
        _style_axes(ax, xlim, ylim, sig_y)
        ax.scatter(lfc[~cand_mask], y[~cand_mask], s=4, c=NEUTRAL["background"], alpha=0.45,
                   linewidths=0, zorder=1)
        nf, _ = _draw_candidates(ax, panel, ring_size=70, fill_size=85, lw=1.6,
                                 labels=False, coincidence="measured", label_fontsize=9)
        n_filled_main += nf
        # zoom rectangle on the main panel
        ax.add_patch(Rectangle((ZOOM_XLIM[0], ZOOM_YLIM[0]), ZOOM_XLIM[1] - ZOOM_XLIM[0],
                               ZOOM_YLIM[1] - ZOOM_YLIM[0], fill=False,
                               edgecolor=NEUTRAL["subtitle"], linewidth=0.9, zorder=8))
        _panel_header(ax, panel)
        ax.set_xlabel("log2 fold change", fontsize=9.5, color=NEUTRAL["subtitle"])
        if i > 0:
            ax.tick_params(labelleft=False)

        # zoom panel: local backdrop + candidates with labels and offsets
        _style_axes(zax, ZOOM_XLIM, ZOOM_YLIM, sig_y)
        inz = (lfc >= ZOOM_XLIM[0]) & (lfc <= ZOOM_XLIM[1]) & (y <= ZOOM_YLIM[1]) & ~cand_mask
        zax.scatter(lfc[inz], y[inz], s=5, c=NEUTRAL["background"], alpha=0.35, linewidths=0, zorder=1)
        nz, recs = _draw_candidates(zax, panel, ring_size=110, fill_size=130, lw=2.2,
                                    labels=True, coincidence="spread", label_fontsize=9.2)
        n_filled_zoom += nz
        for r in recs:
            r["dataset"] = panel.accession
        offset_records += recs
        zax.set_xlabel("log2 fold change", fontsize=9.5, color=NEUTRAL["subtitle"])
        zax.set_title("zoom: boxed region above", fontsize=8.6, color=NEUTRAL["tick"], pad=5)
        if i > 0:
            zax.tick_params(labelleft=False)

    if n_filled_main != EXPECTED_SIGNIFICANT or n_filled_zoom != EXPECTED_SIGNIFICANT:
        raise ValueError(f"variant A drew {n_filled_main}/{n_filled_zoom} filled points, "
                         f"expected {EXPECTED_SIGNIFICANT}")
    main_axes[0].set_ylabel(r"$-\log_{10}$(FDR)", fontsize=10, color=NEUTRAL["subtitle"])
    zoom_axes[0].set_ylabel(r"$-\log_{10}$(FDR)", fontsize=10, color=NEUTRAL["subtitle"])
    main_axes[0].text(xlim[0] + 0.02 * (xlim[1] - xlim[0]), sig_y + 0.10, "FDR 0.05",
                      fontsize=8, color=NEUTRAL["tick"], va="bottom")
    _figure_text(fig, 0.975)
    _legend_strip(fig, 0.040, with_background=True, disclosure=OFFSET_DISCLOSURE)
    written = _save(fig, stub)
    info = {"xlim": xlim, "ylim": ylim, "zoom_xlim": ZOOM_XLIM, "zoom_ylim": ZOOM_YLIM,
            "offsets": offset_records}
    return written, info


# ---------------------------------------------------------------------------
# Variant B: candidates only, tighter shared range, no inset, ZERO displacement
# (coincident candidates are concentric rings at their true position)
# ---------------------------------------------------------------------------
def build_variant_b(panels: list[Panel], stub: Path) -> tuple[dict[str, Path], dict]:
    _pin_reproducibility()
    xlim, ylim = VARIANT_B_XLIM, VARIANT_B_YLIM
    assert_candidates_inside(panels, xlim, ylim, "variant B")
    sig_y = -np.log10(SIG_FDR)

    fig = plt.figure(figsize=(13.3, 6.6))
    gs = fig.add_gridspec(1, 4, left=0.055, right=0.985, top=0.735, bottom=0.165, wspace=0.10)
    axes = [fig.add_subplot(gs[0, i]) for i in range(4)]

    n_filled = 0
    offset_records: list[dict] = []
    for i, (ax, panel) in enumerate(zip(axes, panels)):
        _style_axes(ax, xlim, ylim, sig_y)
        nf, recs = _draw_candidates(ax, panel, ring_size=CONCENTRIC_SIZES[0], fill_size=220, lw=2.6,
                                    labels=True, coincidence="concentric", label_fontsize=10.2)
        n_filled += nf
        for r in recs:
            r["dataset"] = panel.accession
        offset_records += recs
        _panel_header(ax, panel)
        ax.set_xlabel("log2 fold change", fontsize=9.5, color=NEUTRAL["subtitle"])
        if i > 0:
            ax.tick_params(labelleft=False)
    if n_filled != EXPECTED_SIGNIFICANT:
        raise ValueError(f"variant B drew {n_filled} filled points, expected {EXPECTED_SIGNIFICANT}")
    axes[0].set_ylabel(r"$-\log_{10}$(FDR)", fontsize=10, color=NEUTRAL["subtitle"])
    axes[0].text(xlim[0] + 0.02 * (xlim[1] - xlim[0]), sig_y + 0.04, "FDR 0.05",
                 fontsize=8, color=NEUTRAL["tick"], va="bottom")
    _figure_text(fig, 0.965, caption_gap=0.062)
    _legend_strip(fig, 0.02, with_background=False)
    written = _save(fig, stub)
    return written, {"xlim": xlim, "ylim": ylim, "offsets": offset_records}


# ---------------------------------------------------------------------------
# Colour-vision-deficiency simulation (computed directly; no dependency)
# Machado, Oliveira & Fernandes (2009), IEEE TVCG 15(6): severity-1.0
# matrices, applied in LINEAR sRGB.
# ---------------------------------------------------------------------------
CVD_MATRICES: dict[str, np.ndarray] = {
    "protanopia": np.array([[0.152286, 1.052583, -0.204868],
                            [0.114503, 0.786281, 0.099216],
                            [-0.003882, -0.048116, 1.051998]]),
    "deuteranopia": np.array([[0.367322, 0.860646, -0.227968],
                              [0.280085, 0.672501, 0.047413],
                              [-0.011820, 0.042940, 0.968881]]),
}


def _hex_to_rgb01(h: str) -> np.ndarray:
    h = h.lstrip("#")
    return np.array([int(h[i:i + 2], 16) for i in (0, 2, 4)], dtype=float) / 255.0


def _srgb_to_linear(c: np.ndarray) -> np.ndarray:
    return np.where(c <= 0.04045, c / 12.92, ((c + 0.055) / 1.055) ** 2.4)


def _linear_to_srgb(c: np.ndarray) -> np.ndarray:
    c = np.clip(c, 0.0, 1.0)
    return np.where(c <= 0.0031308, 12.92 * c, 1.055 * np.power(c, 1 / 2.4) - 0.055)


def simulate_cvd(rgb01: np.ndarray, kind: str) -> np.ndarray:
    """rgb01: (..., 3) sRGB in [0,1] -> simulated sRGB in [0,1]."""
    lin = _srgb_to_linear(rgb01)
    sim = lin @ CVD_MATRICES[kind].T
    return _linear_to_srgb(sim)


def _rgb01_to_lab(rgb01: np.ndarray) -> np.ndarray:
    lin = _srgb_to_linear(rgb01)
    m = np.array([[0.4124564, 0.3575761, 0.1804375],
                  [0.2126729, 0.7151522, 0.0721750],
                  [0.0193339, 0.1191920, 0.9503041]])
    xyz = lin @ m.T
    xyz = xyz / np.array([0.95047, 1.0, 1.08883])  # D65 white
    eps, kappa = 216 / 24389, 24389 / 27
    f = np.where(xyz > eps, np.cbrt(xyz), (kappa * xyz + 16) / 116)
    L = 116 * f[..., 1] - 16
    a = 500 * (f[..., 0] - f[..., 1])
    b = 200 * (f[..., 1] - f[..., 2])
    return np.stack([L, a, b], axis=-1)


def palette_cvd_report(colours: dict[str, str] = GENE_COLOURS) -> pd.DataFrame:
    """Pairwise CIE76 ΔE and ΔL* between the gene colours under normal vision,
    protanopia and deuteranopia."""
    rows = []
    for kind in ("normal", "protanopia", "deuteranopia"):
        rgb = {g: _hex_to_rgb01(h) for g, h in colours.items()}
        if kind != "normal":
            rgb = {g: simulate_cvd(v, kind) for g, v in rgb.items()}
        lab = {g: _rgb01_to_lab(v) for g, v in rgb.items()}
        for a, b in itertools.combinations(sorted(colours), 2):
            d = lab[a] - lab[b]
            rows.append({"vision": kind, "gene_a": a, "gene_b": b,
                         "delta_e76": float(np.sqrt((d ** 2).sum())),
                         "delta_L": float(abs(d[0]))})
    return pd.DataFrame(rows)


def simulate_png(src: Path, dst_stub: Path) -> dict[str, Path]:
    """Write deuteranopia and protanopia simulations of a rendered PNG."""
    from PIL import Image

    im = np.asarray(Image.open(src).convert("RGB"), dtype=float) / 255.0
    out: dict[str, Path] = {}
    for kind in ("deuteranopia", "protanopia"):
        sim = (simulate_cvd(im, kind) * 255.0).round().astype(np.uint8)
        path = dst_stub.with_name(f"{dst_stub.name}_{kind}.png")
        Image.fromarray(sim).save(path)
        out[kind] = path
    return out


# ---------------------------------------------------------------------------
# Manifest + entry point
# ---------------------------------------------------------------------------
def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_outputs(panels: list[Panel], verification: pd.DataFrame, out_dir: Path,
                  written: dict[str, dict[str, Path]], info: dict[str, dict],
                  cvd_pngs: dict[str, dict[str, Path]]) -> None:
    rows = []
    for variant, files in written.items():
        inf = info[variant]
        row = {
            "figure_name": f"poster_candidate_volcano_v2_{variant}",
            "variant": variant,
            "scientific_question": "How many of the 16 candidate gene-dataset combinations are actually significant?",
            "analysis_status": STATUS_LABEL,
            "post_freeze": "yes",
            "render_module": "src/poster_candidate_volcano_v2.py",
            "wrapper_script": "scripts/poster/02b_candidate_volcano_v2.py",
            "gse118713_contrast": next(p.contrast for p in panels if p.accession == "GSE118713"),
            "significance_threshold_fdr": SIG_FDR,
            "n_significant_candidate_points": EXPECTED_SIGNIFICANT,
            "shared_xlim": f"{inf['xlim'][0]:.3f},{inf['xlim'][1]:.3f}",
            "shared_ylim": f"{inf['ylim'][0]:.3f},{inf['ylim'][1]:.3f}",
            "zoom_xlim": f"{ZOOM_XLIM[0]},{ZOOM_XLIM[1]}" if variant == "variant_a_genomewide" else "",
            "zoom_ylim": f"{ZOOM_YLIM[0]},{ZOOM_YLIM[1]}" if variant == "variant_a_genomewide" else "",
            "coincidence_tol": COINCIDENCE_TOL, "coincidence_step": COINCIDENCE_STEP,
            "n_points_offset": len(inf["offsets"]),
            **{f"colour_{g}": c for g, c in GENE_COLOURS.items()},
            **{f"source_{p.accession.lower()}": p.source_path for p in panels},
            **{f"n_genes_{p.accession.lower()}": p.n_genes for p in panels},
            **{f"sha256_{ext}": _sha256(path) for ext, path in files.items()},
            **{f"sha256_png_{kind}_simulation": _sha256(path) for kind, path in cvd_pngs[variant].items()},
        }
        rows.append(row)
    pd.DataFrame(rows).to_csv(out_dir / "volcano_v2_manifest.tsv", sep="\t", index=False)

    values = pd.concat(
        [p.candidate_rows.reset_index().assign(dataset=p.accession, significant=lambda d: d["fdr"] < SIG_FDR)
         for p in panels], ignore_index=True)[["dataset", "gene", "log2fc", "fdr", "significant"]]
    values.insert(0, "analysis_status", STATUS_LABEL)
    values.to_csv(out_dir / "candidate_values_plotted.tsv", sep="\t", index=False)

    verification = verification.copy()
    verification.insert(0, "analysis_status", STATUS_LABEL)
    verification.to_csv(out_dir / "verification_against_frozen.tsv", sep="\t", index=False)

    offsets = pd.DataFrame(info["variant_a_genomewide"]["offsets"])
    offsets["applies_to"] = "variant_a_genomewide zoom row only; variant B draws concentric rings at measured x"
    offsets.insert(0, "analysis_status", STATUS_LABEL)
    offsets.to_csv(out_dir / "cosmetic_offsets.tsv", sep="\t", index=False)

    cvd = palette_cvd_report()
    cvd.insert(0, "analysis_status", STATUS_LABEL)
    cvd.to_csv(out_dir / "cvd_palette_simulation.tsv", sep="\t", index=False)
    logger.info("wrote manifest, values, verification, offsets and CVD tables to %s", out_dir)


def main(out_dir: Path = OUT_DIR) -> dict[str, dict[str, Path]]:
    panels = load_panels()
    verification = verify_against_frozen(panels)
    out_dir.mkdir(parents=True, exist_ok=True)
    written: dict[str, dict[str, Path]] = {}
    info: dict[str, dict] = {}
    cvd_pngs: dict[str, dict[str, Path]] = {}
    for variant, builder in (("variant_a_genomewide", build_variant_a),
                             ("variant_b_candidates_only", build_variant_b)):
        stub = out_dir / f"candidate_volcano_v2_{variant}"
        written[variant], info[variant] = builder(panels, stub)
        cvd_pngs[variant] = simulate_png(written[variant]["png"], out_dir / "cvd_simulation" / stub.name) \
            if (out_dir / "cvd_simulation").mkdir(parents=True, exist_ok=True) is None else {}
    write_outputs(panels, verification, out_dir, written, info, cvd_pngs)
    return written


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    main()
