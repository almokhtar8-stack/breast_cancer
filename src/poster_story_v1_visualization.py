"""POSTER-STORY-V1 figure bank -- visualization only.

Complete visual-story reset around ONE central heatmap (see
results/reports/poster_story_v1/STORY_PLAN.md). No new analysis, no
recomputation of any frozen p-value/FDR/effect-size, no re-ranking
anywhere in this module. See src/poster_story_v1_data.py for the exact
source of every plotted number.

Style system (identical philosophy to poster_exploration_v3, reused
deliberately -- it already met the poster-grade bar):
  - white background always; large type; one short takeaway line per
    figure; few panels; minimal gridlines; direct labeling
  - candidate identity colors (FOCUS_COLORS): KDM1A=orange, TLK2=light
    blue, USP34=strong blue, VEZF1=gold
"""

from __future__ import annotations

import logging
import shutil
import subprocess
import tempfile
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import TwoSlopeNorm
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch, Rectangle

from src import poster_story_v1_data as sv1
from src.poster_exploration_v2_visualization import _render_kdm1a_tlk2_panels
from src.poster_figures_bank_visualization import _render_structure_bank_panels
from src.poster_figures_visualization import _autocrop_white, rng

logger = logging.getLogger(__name__)

FIGURES = Path("results/figures/poster_story_v1")

FOCUS_FOUR = sv1.FOCUS_FOUR
FOCUS_COLORS = sv1.FOCUS_COLORS
GRAY = sv1.GRAY
LGRAY = sv1.LGRAY
DGRAY = sv1.DGRAY

plt.rcParams["font.family"] = "DejaVu Sans"
plt.rcParams["text.color"] = DGRAY
plt.rcParams["axes.edgecolor"] = DGRAY
plt.rcParams["axes.labelcolor"] = DGRAY
plt.rcParams["xtick.color"] = DGRAY
plt.rcParams["ytick.color"] = DGRAY


def _save(fig, stub: Path, vector: bool = True) -> None:
    stub.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(stub.with_suffix(".png"), facecolor="white", bbox_inches="tight", dpi=300)
    if vector:
        fig.savefig(stub.with_suffix(".pdf"), facecolor="white", bbox_inches="tight")
        fig.savefig(stub.with_suffix(".svg"), facecolor="white", bbox_inches="tight")
    plt.close(fig)
    logger.info("wrote %s.png%s", stub, " + .pdf/.svg" if vector else "")


def _clean_axes(ax, left: bool = True) -> None:
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    if not left:
        ax.spines["left"].set_visible(False)


def _hero_title(fig, text: str, y: float = 1.02, fontsize: float = 22) -> None:
    fig.text(0.01, y, text, fontsize=fontsize, fontweight="bold", color=DGRAY, ha="left", va="bottom")


def _takeaway(fig, text: str, y: float = -0.02, fontsize: float = 12.5) -> None:
    fig.text(0.01, y, text, fontsize=fontsize, color="#4a4a4a", ha="left", va="top")


DATASET_LABELS_SHORT = {
    "GSE118713": "GSE118713",
    "GSE111151": "GSE111151",
    "GSE240112": "GSE240112",
    "GSE245601": "GSE245601",
}


# ===========================================================================
# HERO HEATMAP CANDIDATES (Phase 2 comparison)
# ===========================================================================

def _dataset_blocks(pairs: pd.DataFrame, genes: list[str]) -> list[dict]:
    """One block per dataset, in a fixed row order, each carrying its own
    ref/cmp labels and a genes x 2 value grid."""
    blocks = []
    for dataset in sv1.DATASET_ORDER:
        sub = pairs[pairs["dataset"] == dataset]
        if len(sub) == 0:
            continue
        ref_label = sub["ref_label"].iloc[0]
        cmp_label = sub["cmp_label"].iloc[0]
        ref_vals = [sub[sub["gene"] == g]["ref_value"].iloc[0] if (sub["gene"] == g).any() else np.nan for g in genes]
        cmp_vals = [sub[sub["gene"] == g]["cmp_value"].iloc[0] if (sub["gene"] == g).any() else np.nan for g in genes]
        log2fc = [sub[sub["gene"] == g]["log2fc"].iloc[0] if (sub["gene"] == g).any() else np.nan for g in genes]
        blocks.append(dict(dataset=dataset, ref_label=ref_label, cmp_label=cmp_label,
                            ref_vals=np.array(ref_vals, dtype=float), cmp_vals=np.array(cmp_vals, dtype=float),
                            log2fc=np.array(log2fc, dtype=float)))
    return blocks


def _draw_gene_header(ax, genes: list[str], y: float = 1.0) -> None:
    n = len(genes)
    for j, gene in enumerate(genes):
        color = FOCUS_COLORS.get(gene, GRAY)
        weight = "bold" if gene in FOCUS_FOUR else "normal"
        fontsize = 15 if gene in FOCUS_FOUR else 9
        ax.text(j, y, gene, ha="center", va="bottom", fontsize=fontsize, fontweight=weight, color=color,
                 transform=ax.get_xaxis_transform())


def build_candidate1_raw_paired(stub: Path) -> None:
    """Candidate 1: raw-ish paired rows, colored by WITHIN-DATASET-BLOCK
    min-max scaling (sequential colormap) -- explicitly NOT claiming
    cross-dataset color comparability (TPM vs log2CPM), disclosed in the
    caption."""
    genes = FOCUS_FOUR
    pairs = sv1.build_hero_heatmap_pairs()
    blocks = _dataset_blocks(pairs, genes)

    n_rows = sum(2 for _ in blocks) + (len(blocks) - 1)  # +1 blank spacer row between blocks
    fig, ax = plt.subplots(figsize=(9, 10), dpi=300)
    cmap = plt.get_cmap("Blues")

    y = n_rows - 1
    row_labels = []
    for bi, b in enumerate(blocks):
        grid = np.vstack([b["ref_vals"], b["cmp_vals"]])
        vmin, vmax = np.nanmin(grid), np.nanmax(grid)
        norm_grid = (grid - vmin) / (vmax - vmin + 1e-9)
        for r, label in enumerate([b["ref_label"], b["cmp_label"]]):
            for j in range(len(genes)):
                v = norm_grid[r, j]
                ax.add_patch(Rectangle((j - 0.45, y - 0.4), 0.9, 0.8, facecolor=cmap(v), edgecolor="white", linewidth=1.5))
                ax.text(j, y, f"{grid[r, j]:.2f}", ha="center", va="center", fontsize=8.5,
                         color="white" if v > 0.6 else DGRAY)
            row_labels.append((y, f"{b['dataset']}\n{label}"))
            y -= 1
        if bi < len(blocks) - 1:
            y -= 1  # spacer

    for yy, lbl in row_labels:
        ax.text(-0.7, yy, lbl, ha="right", va="center", fontsize=9.5, color=DGRAY)
    _draw_gene_header(ax, genes)
    ax.set_xlim(-2.6, len(genes) - 0.3)
    ax.set_ylim(-1, n_rows + 0.6)
    ax.axis("off")
    _hero_title(fig, "Candidate 1: raw paired values (within-dataset color scale)", fontsize=15, y=1.03)
    _takeaway(fig, "Color is scaled independently within each dataset block -- NOT comparable across datasets "
                    "(TPM vs. log2CPM). Honest but requires this caveat.", y=-0.03, fontsize=10.5)
    _save(fig, stub)


def build_candidate2_delta_only(stub: Path) -> None:
    """Candidate 2: pure delta-from-control heatmap, ONE row per dataset
    (no paired rows shown) -- the cleanest scientifically, but drops the
    explicit visual pairing the brief asked for."""
    genes = FOCUS_FOUR
    pairs = sv1.build_hero_heatmap_pairs()
    blocks = _dataset_blocks(pairs, genes)

    fig, ax = plt.subplots(figsize=(8, 6), dpi=300)
    grid = np.vstack([b["log2fc"] for b in blocks])
    vmax = np.nanmax(np.abs(grid))
    norm = TwoSlopeNorm(vcenter=0, vmin=-vmax, vmax=vmax)
    im = ax.imshow(grid, cmap="RdBu_r", norm=norm, aspect="auto")
    for i, b in enumerate(blocks):
        for j in range(len(genes)):
            v = grid[i, j]
            ax.text(j, i, f"{v:+.2f}", ha="center", va="center", fontsize=11,
                     color="white" if abs(v) > vmax * 0.55 else DGRAY)
    ax.set_yticks(range(len(blocks)))
    ax.set_yticklabels([f"{b['dataset']}\n({b['cmp_label']} vs. {b['ref_label']})" for b in blocks], fontsize=9.5)
    ax.set_xticks(range(len(genes)))
    ax.set_xticklabels(genes, fontsize=13, fontweight="bold")
    for tick, gene in zip(ax.get_xticklabels(), genes):
        tick.set_color(FOCUS_COLORS[gene])
    ax.tick_params(length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)
    cb = fig.colorbar(im, ax=ax, fraction=0.05, pad=0.03)
    cb.set_label("log2 fold-change", fontsize=10)
    _hero_title(fig, "Candidate 2: delta-only (single row per dataset)", fontsize=15, y=1.03)
    _takeaway(fig, "Cleanest scientifically (log2FC is comparable across studies) but does not show the paired "
                    "conditions themselves.", y=-0.04, fontsize=10.5)
    _save(fig, stub)


def build_candidate3_paired_plus_delta_annotation(stub: Path) -> None:
    """Candidate 3: paired raw-ish rows (as candidate 1) PLUS a small
    delta arrow/value annotation on the right margin -- both the raw
    pattern and the honest quantitative delta are visible."""
    genes = FOCUS_FOUR
    pairs = sv1.build_hero_heatmap_pairs()
    blocks = _dataset_blocks(pairs, genes)

    n_rows = sum(2 for _ in blocks) + (len(blocks) - 1)
    fig, ax = plt.subplots(figsize=(11, 10), dpi=300)
    cmap = plt.get_cmap("Purples")

    y = n_rows - 1
    row_labels = []
    for bi, b in enumerate(blocks):
        grid = np.vstack([b["ref_vals"], b["cmp_vals"]])
        vmin, vmax = np.nanmin(grid), np.nanmax(grid)
        norm_grid = (grid - vmin) / (vmax - vmin + 1e-9)
        y_ref, y_cmp = y, y - 1
        for r, (label, yy) in enumerate([(b["ref_label"], y_ref), (b["cmp_label"], y_cmp)]):
            for j in range(len(genes)):
                v = norm_grid[r, j]
                ax.add_patch(Rectangle((j - 0.45, yy - 0.4), 0.9, 0.8, facecolor=cmap(v), edgecolor="white", linewidth=1.5))
                ax.text(j, yy, f"{grid[r, j]:.2f}", ha="center", va="center", fontsize=8.5,
                         color="white" if v > 0.6 else DGRAY)
            row_labels.append((yy, f"{b['dataset']}\n{label}"))
        ax.plot([len(genes) - 0.3, len(genes) - 0.1], [y_ref, y_ref], color=LGRAY, linewidth=1)
        ax.plot([len(genes) - 0.3, len(genes) - 0.1], [y_cmp, y_cmp], color=LGRAY, linewidth=1)
        for j in range(len(genes)):
            fc = b["log2fc"][j]
            color = "#b0392f" if fc > 0 else "#2b6cb0"
            ax.annotate("", xy=(j, y_cmp - 0.55), xytext=(j, y_ref + 0.55),
                        arrowprops=dict(arrowstyle="-|>", color=color, lw=1.6, alpha=0.85))
        y -= 3
        if bi < len(blocks) - 1:
            y -= 0  # already spaced by the 3-step block advance

    for yy, lbl in row_labels:
        ax.text(-0.7, yy, lbl, ha="right", va="center", fontsize=9.5, color=DGRAY)
    _draw_gene_header(ax, genes)
    ax.set_xlim(-2.9, len(genes) + 0.3)
    ax.set_ylim(min(y for y, _ in row_labels) - 1.2, n_rows + 0.6)
    ax.axis("off")
    _hero_title(fig, "Candidate 3: paired rows + delta-direction annotation", fontsize=15, y=1.03)
    _takeaway(fig, "Arrows show direction of change; color remains within-block only (same caveat as Candidate 1).",
              y=-0.03, fontsize=10.5)
    _save(fig, stub)


def _build_hybrid_heatmap(pairs: pd.DataFrame, genes: list[str], figsize: tuple, title: str,
                           gene_fontsize_focus: float = 15, gene_fontsize_other: float = 9) -> plt.Figure:
    """The Candidate 4/5 hybrid design: reference row always neutral gray
    (it IS the reference point, not zero signal); comparison row colored
    by log2 fold-change vs. that same-dataset reference -- honest (log2FC
    is a comparable unit across studies), visually paired (bracket +
    shared block), and never mixes TPM/log2CPM colors on one scale."""
    blocks = _dataset_blocks(pairs, genes)
    vmax = np.nanmax(np.abs(np.vstack([b["log2fc"] for b in blocks])))
    norm = TwoSlopeNorm(vcenter=0, vmin=-vmax, vmax=vmax)
    cmap = plt.get_cmap("RdBu_r")

    n_blocks = len(blocks)
    fig = plt.figure(figsize=figsize, dpi=300)
    gs = fig.add_gridspec(n_blocks, 1, hspace=0.55)

    for bi, b in enumerate(blocks):
        ax = fig.add_subplot(gs[bi])
        for j in range(len(genes)):
            fc = b["log2fc"][j]
            not_testable = not np.isfinite(fc)
            ax.add_patch(Rectangle((j - 0.47, 0.55), 0.94, 0.9, facecolor="#e8e8e8", edgecolor="white", linewidth=1.2, zorder=2))
            cmp_color = "#f5f5f5" if not_testable else cmap(norm(fc))
            ax.add_patch(Rectangle((j - 0.47, -0.45), 0.94, 0.9, facecolor=cmp_color, edgecolor="white", linewidth=1.2,
                                    hatch="////" if not_testable else None, zorder=2))
        ax.plot([-0.75, -0.75], [0.15, 0.95], color=DGRAY, linewidth=1.6, zorder=1)
        ax.plot([-0.75, -0.6], [0.55, 0.55], color=DGRAY, linewidth=1.6, zorder=1)
        ax.plot([-0.75, -0.75], [-0.9, -0.1], color=DGRAY, linewidth=1.6, zorder=1)
        ax.plot([-0.75, -0.6], [-0.5, -0.5], color=DGRAY, linewidth=1.6, zorder=1)
        ax.text(-0.85, 1.0, b["ref_label"], ha="right", va="center", fontsize=10.5, color=DGRAY)
        ax.text(-0.85, 0.1, b["cmp_label"], ha="right", va="center", fontsize=10.5, color=DGRAY, fontweight="bold")
        ax.text(-2.7, 1.75, b["dataset"], ha="left", va="bottom", fontsize=13, fontweight="bold", color=DGRAY)
        if bi == 0:
            for j, gene in enumerate(genes):
                color = FOCUS_COLORS.get(gene, GRAY)
                weight = "bold" if gene in FOCUS_FOUR else "normal"
                fs = gene_fontsize_focus if gene in FOCUS_FOUR else gene_fontsize_other
                ax.text(j, 1.75, gene, ha="center", va="bottom", fontsize=fs, fontweight=weight, color=color)
        for j in range(len(genes)):
            fc = b["log2fc"][j]
            if not np.isfinite(fc):
                ax.text(j, 0.0, "n/t", ha="center", va="center", fontsize=8.2, color=GRAY, style="italic", zorder=3)
            else:
                ax.text(j, 0.0, f"{fc:+.2f}", ha="center", va="center", fontsize=8.6,
                         color="white" if abs(fc) > vmax * 0.6 else DGRAY, zorder=3)
        ax.set_xlim(-2.9, len(genes) - 0.2)
        ax.set_ylim(-1.05, 2.05)
        ax.axis("off")

    fig.suptitle(title, fontsize=15.5, fontweight="bold", color=DGRAY, x=0.02, ha="left", y=1.01)
    return fig


def build_candidate4_hybrid_hero(stub: Path) -> None:
    genes = FOCUS_FOUR
    pairs = sv1.build_hero_heatmap_pairs()
    fig = _build_hybrid_heatmap(pairs, genes, (9.5, 10.5), "Candidate 4: paired-row hybrid (4 focus genes) -- HERO")
    _takeaway(fig, "Gray = reference condition (baseline). Colored = log2 fold-change vs. that same-dataset "
                    "reference. Comparable across datasets because log2FC is a ratio, not a raw scale.", y=-0.02)
    _save(fig, stub)


def build_hero_heatmap(stub: Path) -> None:
    """The selected hero figure for the main figure bank -- identical
    design/data to Candidate 4 (selected after comparing all 5
    candidates, see FINAL_FIGURE_RECOMMENDATION.md), with poster-final
    framing text instead of "Candidate 4" comparison-sheet labeling."""
    genes = FOCUS_FOUR
    pairs = sv1.build_hero_heatmap_pairs()
    fig = _build_hybrid_heatmap(pairs, genes, (9.5, 10.5),
                                 "How do the 4 focus genes behave across every real transcriptomic context?")
    _takeaway(fig, "Gray = reference condition. Colored = log2 fold-change vs. that same-dataset reference "
                    "(comparable across studies). Two resistance-model cell-line panels, one human recurrence "
                    "panel, one acute human ex-vivo panel.", y=-0.02)
    _save(fig, stub)


def build_candidate5_hybrid_13hit(stub: Path) -> None:
    sens = sv1.load_significant_sensitising_hits().sort_values("rank_by_effect")
    genes = sens["gene"].tolist()
    m = sv1.load_13hit_log2fc_matrix().set_index("gene").loc[genes]

    pairs_rows = []
    label_map = {
        "GSE118713": ("MCF7 (baseline)", "TAMR (resistant)"),
        "GSE111151": ("Parental", "Resistant"),
        "GSE240112": ("Primary", "Recurrent"),
        "GSE245601": ("Control", "Tamoxifen (12h)"),
    }
    for dataset in sv1.DATASET_ORDER:
        ref_label, cmp_label = label_map[dataset]
        for gene in genes:
            pairs_rows.append(dict(dataset=dataset, ref_label=ref_label, cmp_label=cmp_label, gene=gene,
                                    ref_value=np.nan, cmp_value=np.nan, log2fc=m.loc[gene, dataset]))
    pairs = pd.DataFrame(pairs_rows)

    fig = _build_hybrid_heatmap(pairs, genes, (15, 10.5),
                                 "Candidate 5: paired-row hybrid, all 13 significant sensitising hits",
                                 gene_fontsize_focus=13, gene_fontsize_other=8.5)
    _takeaway(fig, "Same design and log2FC encoding as Candidate 4, extended to all 13 significant sensitising "
                    "hits -- shows the 4 focus genes (bold, colored) are not cherry-picked in isolation.", y=-0.02)
    _save(fig, stub)


HEATMAP_CANDIDATE_ITEMS = [
    ("CANDIDATE1_raw_paired", "Candidate 1: raw paired (within-block scale)"),
    ("CANDIDATE2_delta_only", "Candidate 2: delta-only (unpaired rows)"),
    ("CANDIDATE3_paired_plus_delta_annotation", "Candidate 3: paired + delta arrows"),
    ("CANDIDATE4_hybrid_hero", "Candidate 4: paired-row hybrid (HERO, 4 genes)"),
    ("CANDIDATE5_hybrid_13hit", "Candidate 5: paired-row hybrid (13 hits)"),
]


def build_heatmap_contact_sheet(out_dir: Path = FIGURES) -> None:
    n = len(HEATMAP_CANDIDATE_ITEMS)
    fig, axes = plt.subplots(3, 2, figsize=(16, 22), dpi=150)
    axes = axes.flatten()
    for ax, (stem, label) in zip(axes, HEATMAP_CANDIDATE_ITEMS):
        img = plt.imread(out_dir / f"{stem}.png")
        ax.imshow(img)
        ax.axis("off")
        ax.set_title(label, fontsize=13, fontweight="bold", color=DGRAY, pad=8)
    for ax in axes[n:]:
        ax.axis("off")
    fig.suptitle("Hero heatmap candidates -- side-by-side comparison", fontsize=19, fontweight="bold", y=1.0)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(out_dir / "CONTACT_HEATMAP_CANDIDATES.png", facecolor="white", bbox_inches="tight", dpi=150)
    plt.close(fig)
    logger.info("wrote %s", out_dir / "CONTACT_HEATMAP_CANDIDATES.png")


# ===========================================================================
# Remaining figure bank (Phase 3)
# ===========================================================================

from src import poster_exploration_v3_visualization as pv3v  # noqa: E402


def build_discovery(stub: Path) -> None:
    """CRISPR discovery -- reused design from poster_exploration_v3 (F1),
    already independently judged poster-grade; not redesigned here."""
    pv3v.build_fig1_crispr_discovery(stub)


def build_pathway_convergence(stub: Path) -> None:
    """Pathway convergence -- reused design from poster_exploration_v3
    (F3), already independently judged the strongest pathway-biology
    visual hero in that phase; not redesigned here."""
    pv3v.build_fig3_pathway_convergence(stub)


def build_structural_comparison(stub: Path) -> None:
    """Structural / pharmacological comparison -- reused design from
    poster_exploration_v3 (F6), already independently judged the
    strongest visual hero in that phase; not redesigned here."""
    pv3v.build_fig6_structural_comparison(stub)


def build_depmap_context(stub: Path) -> None:
    """Baseline DepMap dependency -- reused design from
    poster_exploration_v3 (F5); not redesigned here."""
    pv3v.build_fig5_depmap_context(stub)


def build_final_synthesis(stub: Path) -> None:
    """Final synthesis / why-these-genes -- reused design from
    poster_exploration_v3 (F4), already independently judged elegant and
    non-dashboard; not redesigned here."""
    pv3v.build_fig4_postaudit_interpretation(stub)


def build_disease_clinical_context(stub: Path) -> None:
    """NEW this phase: real human recurrence signal (GSE240112) combined
    with real malignant-vs-non-malignant tumour-cell specificity
    (GSE245601, copyKAT-derived malignancy calls) -- a real data layer
    never visualized in any prior poster phase (see DATA_AUDIT.md)."""
    d_recur = sv1.build_gse240112_recurrence_delta()
    d_malig = sv1.build_malignant_vs_nonmalignant_paired_delta()

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(16, 8), dpi=300, gridspec_kw=dict(wspace=0.35))
    pv3v._combined_delta_panel(axA, d_recur, "gene", "delta_log2cpm",
                                "GSE240112 -- recurrent vs. primary (unpaired)")
    axA.set_ylabel(r"$\Delta$log$_2$(CPM) vs. primary mean", fontsize=14)

    pv3v._combined_delta_panel(axB, d_malig, "gene", "delta_log2",
                                "GSE245601 -- malignant vs. non-malignant (paired)")
    axB.set_ylabel(r"$\Delta$log$_2$(CPM+1), malignant - non-malignant", fontsize=13.5)

    handles_a = [Line2D([0], [0], marker="o", color="none", markerfacecolor=GRAY, markeredgecolor="white",
                         markersize=10, label="primary (n=3) / recurrent (n=3) tumour")]
    axA.legend(handles=handles_a, loc="upper left", fontsize=10, frameon=False)
    handles_b = [Line2D([0], [0], marker="o", color="none", markerfacecolor=GRAY, markeredgecolor="white",
                         markersize=10, label="patient (n=5), malignant vs. own non-malignant cells")]
    axB.legend(handles=handles_b, loc="upper left", fontsize=9.5, frameon=False)

    _hero_title(fig, "Is the signal specific to recurrence, and to the tumour cells themselves?", y=1.06)
    _takeaway(fig, "Right panel: malignancy calls are copyKAT-derived (real, already-frozen); a positive delta "
                    "means higher expression in malignant cells from the same patient's own tissue.", y=-0.05)
    _save(fig, stub)


def build_network_backup(stub: Path) -> None:
    """Backup/appendix only -- explicitly NOT part of the main sequence.
    USP34's real 1-hop STRING neighbors (the richest coverage among the
    original 4 candidates); KDM1A and TLK2 have zero frozen network rows
    at all (see DATA_AUDIT.md) so this can never carry a 4-focus-gene
    story -- shown honestly as a single-gene backup panel only."""
    nb = sv1.load_direct_neighbors("USP34").sort_values("confidence", ascending=False).reset_index(drop=True)
    n = len(nb)
    angles = np.linspace(0, 2 * np.pi, n, endpoint=False)
    radius = 1.0

    fig, ax = plt.subplots(figsize=(8, 8), dpi=300)
    ax.scatter([0], [0], s=380, color=FOCUS_COLORS["USP34"], edgecolor="white", linewidth=1.6, zorder=5)
    ax.text(0, -0.16, "USP34", ha="center", va="top", fontsize=12, fontweight="bold", color=FOCUS_COLORS["USP34"], zorder=6)
    for angle, (_, r) in zip(angles, nb.iterrows()):
        px, py = radius * np.cos(angle), radius * np.sin(angle)
        conf = float(r["confidence"]) if pd.notna(r["confidence"]) else 0.5
        ax.plot([0, px], [0, py], color=GRAY, linewidth=0.5 + 2.2 * conf, alpha=0.35 + 0.5 * conf, zorder=1)
        ax.scatter([px], [py], s=50, color=GRAY, edgecolor="white", linewidth=0.6, zorder=3)
        ha = "left" if px >= 0 else "right"
        ax.text(px * 1.13, py * 1.13, r["neighbor_gene"], fontsize=7.2, color=DGRAY, ha=ha, va="center")
    ax.set_xlim(-1.6, 1.6)
    ax.set_ylim(-1.6, 1.6)
    ax.set_aspect("equal")
    ax.axis("off")
    _hero_title(fig, "Backup: USP34's real STRING network (KDM1A/TLK2 have none)", y=1.02, fontsize=15)
    _takeaway(fig, "Explicitly NOT a main-sequence figure -- KDM1A and TLK2 have zero frozen network rows in this "
                    "project. Kept here only as an honest appendix item.", y=-0.02)
    _save(fig, stub, vector=False)


# ===========================================================================
# Main-set contact sheet + runner
# ===========================================================================

MAIN_STORY_ITEMS = [
    ("FIG1_crispr_discovery", "1. CRISPR discovery"),
    ("HERO_cross_context_heatmap", "2. HERO: cross-context heatmap"),
    ("FIG2_structural_comparison", "3. Structural comparison"),
    ("FIG3_pathway_convergence", "4. Pathway convergence"),
    ("FIG4_disease_clinical_context", "5. Disease / clinical context"),
    ("FIG5_depmap_context", "6. Baseline DepMap dependency"),
    ("FIG6_final_synthesis", "7. Final synthesis"),
]


def build_main_story_contact_sheet(out_dir: Path = FIGURES) -> None:
    n = len(MAIN_STORY_ITEMS)
    ncols = 2
    nrows = -(-n // ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(9 * ncols, 6.4 * nrows), dpi=160)
    axes = axes.flatten()
    for ax, (stem, label) in zip(axes, MAIN_STORY_ITEMS):
        img = plt.imread(out_dir / f"{stem}.png")
        ax.imshow(img)
        ax.axis("off")
        weight = "bold"
        fontsize = 15 if stem.startswith("HERO") else 13
        ax.set_title(label, fontsize=fontsize, fontweight=weight, color=DGRAY, pad=8)
    for ax in axes[n:]:
        ax.axis("off")
    fig.suptitle("Poster story v1 -- main figure sequence (7 figures, heatmap as hero)", fontsize=19,
                 fontweight="bold", y=1.0)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(out_dir / "CONTACT_MAIN_STORY.png", facecolor="white", bbox_inches="tight", dpi=160)
    plt.close(fig)
    logger.info("wrote %s", out_dir / "CONTACT_MAIN_STORY.png")


def run(out_dir: str | Path = FIGURES) -> None:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    build_candidate1_raw_paired(out_dir / "CANDIDATE1_raw_paired")
    build_candidate2_delta_only(out_dir / "CANDIDATE2_delta_only")
    build_candidate3_paired_plus_delta_annotation(out_dir / "CANDIDATE3_paired_plus_delta_annotation")
    build_candidate4_hybrid_hero(out_dir / "CANDIDATE4_hybrid_hero")
    build_candidate5_hybrid_13hit(out_dir / "CANDIDATE5_hybrid_13hit")
    build_heatmap_contact_sheet(out_dir)

    build_hero_heatmap(out_dir / "HERO_cross_context_heatmap")
    build_discovery(out_dir / "FIG1_crispr_discovery")
    build_structural_comparison(out_dir / "FIG2_structural_comparison")
    build_pathway_convergence(out_dir / "FIG3_pathway_convergence")
    build_disease_clinical_context(out_dir / "FIG4_disease_clinical_context")
    build_depmap_context(out_dir / "FIG5_depmap_context")
    build_final_synthesis(out_dir / "FIG6_final_synthesis")
    build_network_backup(out_dir / "BACKUP_network_usp34")
    build_main_story_contact_sheet(out_dir)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run()
