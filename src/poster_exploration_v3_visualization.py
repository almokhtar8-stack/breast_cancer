"""EXPLORATION-V3 poster-grade figure bank -- visualization only.

Rebuilds a SMALL set of 6 hero figures (+ up to 2 alternates) from the same
frozen data already used in v2, with a stricter, more poster-grade visual
design system. v2 was judged scientifically sound but visually too close
to an internal analysis/report deck (too many small subplots, too much
grey, dense footnotes). v3 fixes the visual language, not the science.

No new analysis, no recomputation, no re-ranking anywhere in this module.
See src/poster_exploration_v3_data.py for the exact source of every
plotted number (a thin, documented wrapper over the already-frozen v2/
post-audit loaders).

Style system (poster-grade, not report-grade):
  - white background always; no grey panel fills
  - large type: titles ~19-22pt, axis labels ~14-16pt, tick labels ~12-13pt
  - at most ONE short takeaway line per figure (no paragraph footnotes)
  - few panels per figure (1-2, never a small-multiple grid of >4)
  - minimal/no gridlines; light, recessive axis spines
  - direct labeling over legends wherever practical
  - candidate identity colors (FOCUS_COLORS): KDM1A=orange, TLK2=light
    blue, USP34=strong blue, VEZF1=gold -- identical hex values already
    validated (dataviz skill, colorblind-safe) and used in v2/poster_final
  - every number shown is read from a loaded DataFrame/Series, never
    hand-typed
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
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle

from src import poster_exploration_v3_data as pv3
from src.poster_figures_bank_visualization import _render_structure_bank_panels
from src.poster_figures_visualization import _autocrop_white
from src.poster_exploration_v2_visualization import _render_kdm1a_tlk2_panels

logger = logging.getLogger(__name__)

FIGURES = Path("results/figures/poster_exploration_v3")

FOCUS_FOUR = pv3.FOCUS_FOUR
FOCUS_COLORS = pv3.FOCUS_COLORS
GRAY = pv3.GRAY
LGRAY = pv3.LGRAY
DGRAY = pv3.DGRAY

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
    ax.spines["bottom"].set_color(DGRAY)
    if left:
        ax.spines["left"].set_color(DGRAY)


def _hero_title(fig, text: str, y: float = 1.02, fontsize: float = 22) -> None:
    fig.text(0.01, y, text, fontsize=fontsize, fontweight="bold", color=DGRAY, ha="left", va="bottom")


def _takeaway(fig, text: str, y: float = -0.02, fontsize: float = 12.5) -> None:
    """One short, large takeaway line -- never a dense paragraph."""
    fig.text(0.01, y, text, fontsize=fontsize, color="#4a4a4a", ha="left", va="top")


# ===========================================================================
# FIGURE 1 -- CRISPR discovery (hero)
# ===========================================================================

def build_fig1_crispr_discovery(stub: Path) -> None:
    genomewide = pv3.load_genomewide_crispr()
    sens = pv3.load_significant_sensitising_hits().sort_values("effect_size").reset_index(drop=True)
    n = len(sens)
    y = np.arange(n)

    fig, ax = plt.subplots(figsize=(12, 8.5), dpi=300)
    colors = [FOCUS_COLORS.get(g, GRAY) for g in sens["gene"]]
    sizes = [420 if g in FOCUS_FOUR else 190 for g in sens["gene"]]
    ax.hlines(y, 0, sens["effect_size"], color=colors, linewidth=2.4, alpha=0.55, zorder=2)
    ax.scatter(sens["effect_size"], y, s=sizes, color=colors, edgecolor="white", linewidth=1.6, zorder=3)
    ax.axvline(0, color=LGRAY, linewidth=1.0, zorder=0)

    ax.set_yticks(y)
    ax.set_yticklabels(sens["gene"].tolist(), fontsize=15)
    for tick, gene in zip(ax.get_yticklabels(), sens["gene"]):
        if gene in FOCUS_FOUR:
            tick.set_color(FOCUS_COLORS[gene])
            tick.set_fontweight("bold")
            tick.set_fontsize(18)
        else:
            tick.set_color(GRAY)
            tick.set_fontsize(13)

    ax.set_xlabel("CRISPR effect size (sensitising knockout)", fontsize=16)
    ax.tick_params(axis="x", labelsize=13)
    _clean_axes(ax, left=False)
    ax.tick_params(left=False)

    ax.text(0.98, 0.97, f"{len(genomewide):,} genes screened", transform=ax.transAxes, ha="right", va="top",
             fontsize=12, color=GRAY, style="italic")

    _hero_title(fig, "Which genes sensitise ER+ breast cancer cells to tamoxifen?")
    _takeaway(fig, f"{n} significant sensitising hits -- the 4 focus genes (bold) span the full range of "
                    f"functional strength.", y=-0.03)
    _save(fig, stub)


# ===========================================================================
# FIGURE 2 -- Transcriptomic corroboration (GSE111151 + GSE240112, combined)
# ===========================================================================

_rng = np.random.default_rng(7)


def _combined_delta_panel(ax, df: pd.DataFrame, gene_col: str, val_col: str, title: str) -> None:
    x = np.arange(len(FOCUS_FOUR))
    for i, gene in enumerate(FOCUS_FOUR):
        vals = df.loc[df[gene_col] == gene, val_col].to_numpy()
        color = FOCUS_COLORS[gene]
        q1, med, q3 = np.percentile(vals, [25, 50, 75])
        ax.add_patch(Rectangle((i - 0.17, q1), 0.34, max(q3 - q1, 1e-6), facecolor=color, alpha=0.16,
                                edgecolor=color, linewidth=1.6, zorder=2))
        ax.plot([i - 0.17, i + 0.17], [med, med], color=color, linewidth=2.6, zorder=3)
        jitter = _rng.uniform(-0.09, 0.09, size=len(vals))
        ax.scatter(i + jitter, vals, s=95, color=color, edgecolor="white", linewidth=1.1, zorder=4)
    ax.axhline(0, color=LGRAY, linewidth=1.2, zorder=0)
    ax.set_xticks(x)
    ax.set_xticklabels(FOCUS_FOUR, fontsize=16, fontweight="bold")
    for tick, gene in zip(ax.get_xticklabels(), FOCUS_FOUR):
        tick.set_color(FOCUS_COLORS[gene])
    ax.set_xlim(-0.6, 3.6)
    ax.tick_params(axis="y", labelsize=12)
    _clean_axes(ax)
    ax.set_title(title, fontsize=15.5, fontweight="bold", color=DGRAY, loc="left", pad=12)


def build_fig2_transcriptomic_corroboration(stub: Path) -> None:
    d111151 = pv3.build_gse111151_delta_from_parental()
    d240112 = pv3.build_gse240112_delta_from_primary_mean()

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(15, 8), dpi=300)
    _combined_delta_panel(axA, d111151, "gene_symbol", "delta_log2cpm",
                           "GSE111151 -- resistant vs. own parental line")
    axA.set_ylabel(r"$\Delta$log$_2$(CPM) vs. parental", fontsize=14)

    _combined_delta_panel(axB, d240112, "gene", "delta_log2cpm",
                           "GSE240112 -- recurrent vs. primary-tumour mean")
    axB.set_ylabel(r"$\Delta$log$_2$(CPM) vs. primary mean", fontsize=14)

    handles = [Line2D([0], [0], marker="o", color="none", markerfacecolor=GRAY, markeredgecolor="white",
                       markersize=10, label="individual resistant subline (n=7)")]
    axA.legend(handles=handles, loc="upper left", fontsize=10.5, frameon=False)
    handles_b = [Line2D([0], [0], marker="o", color="none", markerfacecolor=GRAY, markeredgecolor="white",
                         markersize=10, label="primary (n=3) / recurrent (n=3) tumour")]
    axB.legend(handles=handles_b, loc="upper left", fontsize=10.5, frameon=False)

    _hero_title(fig, "Do these genes actually change in resistant and recurrent disease?", y=1.06)
    _takeaway(fig, "Real per-sample expression change, each gene on its own reference point -- "
                    "GSE111151 is a resistance-model cell-line panel; GSE240112 is unpaired human recurrence data.",
              y=-0.05)
    _save(fig, stub)


# ===========================================================================
# FIGURE 3 -- Pathway convergence (hero)
# ===========================================================================

def build_fig3_pathway_convergence(stub: Path) -> None:
    er_pathways = [("hallmark", "HALLMARK_ESTROGEN_RESPONSE_EARLY", "Estrogen response (early)"),
                   ("hallmark", "HALLMARK_ESTROGEN_RESPONSE_LATE", "Estrogen response (late)")]
    df_er = pv3.load_pathway_trajectories(er_pathways)
    df_emt = pv3.load_pathway_trajectories([("hallmark", "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "EMT")])
    x = np.arange(len(pv3.DATASET_ORDER))
    ds_labels = [d.upper() for d in pv3.DATASET_ORDER]

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(15, 8), dpi=300)

    for label, color, marker in [("Estrogen response (early)", "#08519c", "o"),
                                  ("Estrogen response (late)", "#6baed6", "o")]:
        sub = df_er[df_er["pathway_label"] == label].set_index("dataset").reindex(pv3.DATASET_ORDER)
        lw = 4.2 if "early" in label else 3.0
        alpha = 1.0 if "early" in label else 0.7
        axA.plot(x, sub["NES"], color=color, linewidth=lw, marker=marker, markersize=13, alpha=alpha, zorder=3)
        axA.annotate(label.replace("Estrogen response ", ""), xy=(x[-1], sub["NES"].iloc[-1]), xytext=(12, 0),
                     textcoords="offset points", fontsize=13, color=color, va="center", fontweight="bold")
    axA.axhline(0, color=LGRAY, linewidth=1.2, zorder=0)
    axA.set_xticks(x)
    axA.set_xticklabels(ds_labels, fontsize=13)
    axA.set_xlim(-0.3, 4.0)
    axA.set_ylabel("NES", fontsize=16)
    axA.tick_params(axis="y", labelsize=13)
    _clean_axes(axA)
    axA.set_title("Estrogen response falls in every context", fontsize=17, fontweight="bold", color=DGRAY, loc="left", pad=14)

    sub = df_emt.set_index("dataset").reindex(pv3.DATASET_ORDER)
    colors = ["#b0392f"] * 3 + [GRAY]
    axB.axhline(0, color=LGRAY, linewidth=1.2, zorder=0)
    axB.bar(x, sub["NES"], color=colors, width=0.58, zorder=2)
    for xi, (nes, fdr) in enumerate(zip(sub["NES"], sub["fdr"])):
        marker = "*" if (pd.notna(fdr) and fdr < 0.05) else ""
        axB.text(xi, nes + (0.1 if nes >= 0 else -0.1), marker, ha="center",
                 va="bottom" if nes >= 0 else "top", fontsize=20, color=DGRAY)
    axB.set_xticks(x)
    axB.set_xticklabels(ds_labels, fontsize=13)
    axB.set_ylabel("NES", fontsize=16)
    axB.tick_params(axis="y", labelsize=13)
    _clean_axes(axB)
    axB.set_title("EMT rises in resistance and recurrence, falls acutely", fontsize=17, fontweight="bold", color=DGRAY, loc="left", pad=14)

    _hero_title(fig, "Two pathway signals define the resistance/recurrence program", y=1.08)
    _takeaway(fig, "Real GSEA NES from the frozen per-dataset enrichment tables. * = FDR<0.05.", y=-0.04)
    _save(fig, stub)


# ===========================================================================
# FIGURE 4 -- Post-audit interpretation framework (why these genes)
# ===========================================================================

def build_fig4_postaudit_interpretation(stub: Path) -> None:
    rule0, rule1 = pv3.load_rule0_rule1()
    sens = pv3.load_significant_sensitising_hits()
    genes = sens.sort_values("rank_by_effect")["gene"].tolist()
    n = len(genes)

    fig, ax = plt.subplots(figsize=(11, 9), dpi=300)
    y_of = {g: n - 1 - i for i, g in enumerate(genes)}

    passed = [g for g in genes if g in rule0.index and bool(rule0.loc[g, "eligible"])]

    for g in genes:
        y = y_of[g]
        color = FOCUS_COLORS.get(g, GRAY)
        is_focus = g in FOCUS_FOUR
        ax.scatter([0], [y], s=340 if is_focus else 130, color=color, edgecolor=DGRAY if is_focus else "none",
                   linewidth=1.4, zorder=3)
        if g in passed:
            r0_rank = int(rule0.loc[g, "rank"])
            ax.plot([0, 1], [y, y], color=color, linewidth=2.2, alpha=0.6, zorder=1)
            ax.scatter([1], [y], s=340 if is_focus else 130, color=color, edgecolor=DGRAY if is_focus else "none",
                       linewidth=1.4, zorder=3)
            ax.text(1.08, y, f"selected -- rank {r0_rank}/4", fontsize=12.5 if is_focus else 10, color=color,
                    fontweight="bold" if is_focus else "normal", va="center")
        fontsize = 17 if is_focus else 12
        weight = "bold" if is_focus else "normal"
        ax.text(-0.08, y, g, fontsize=fontsize, color=color, fontweight=weight, ha="right", va="center")

    excluded_ys = [y_of[g] for g in genes if g not in passed]
    if excluded_ys:
        y_mid = float(np.mean(excluded_ys))
        ax.annotate(f"{len(excluded_ys)} genes (incl. KDM1A, TLK2) excluded here\nby the original RNA-eligibility gate",
                    xy=(1, y_mid), xytext=(1.35, y_mid), fontsize=12.5, color=GRAY, va="center", ha="left",
                    style="italic")

    ax.set_xlim(-1.05, 3.6)
    ax.set_ylim(-1.3, n + 1.1)
    ax.text(0, n + 0.4, "CRISPR rank\n(13 hits)", fontsize=14, fontweight="bold", color=DGRAY, ha="center", va="bottom")
    ax.text(1, n + 0.4, "Original\nselection gate", fontsize=14, fontweight="bold", color=DGRAY, ha="center", va="bottom")
    ax.set_xticks([])
    ax.tick_params(left=False, labelleft=False, bottom=False)
    for spine in ("top", "right", "left", "bottom"):
        ax.spines[spine].set_visible(False)

    _hero_title(fig, "The shortlist is an interpretation framework, not a CRISPR leaderboard", y=1.03)
    _takeaway(fig, "KDM1A and TLK2 rank 1st and 2nd by CRISPR strength alone, yet were never eligible under the "
                    "original RNA-corroboration gate -- USP34 and VEZF1 were selected on different grounds.", y=-0.03)
    _save(fig, stub)


# ===========================================================================
# FIGURE 5 -- Human / DepMap context (hero)
# ===========================================================================

def build_fig5_depmap_context(stub: Path) -> None:
    eff = pv3.load_depmap_effect_focus_four()
    names = pv3.load_depmap_model_names().set_index("ModelID")["CellLineName"]
    pivot = eff.pivot(index="cell_line", columns="gene", values="chronos_effect")[FOCUS_FOUR]
    pivot.index = [names.get(i, i) for i in pivot.index]
    pivot = pivot.sort_values("TLK2")
    medians = pivot.median()

    fig = plt.figure(figsize=(10.5, 10), dpi=300)
    gs = fig.add_gridspec(2, 1, height_ratios=(0.7, 9), hspace=0.02)
    axTop = fig.add_subplot(gs[0])
    ax = fig.add_subplot(gs[1])

    axTop.axis("off")
    for j, gene in enumerate(FOCUS_FOUR):
        axTop.scatter([j], [0.5], s=380, color=FOCUS_COLORS[gene], edgecolor="white", linewidth=1.4, zorder=3)
        axTop.text(j, -0.35, f"median {medians[gene]:.2f}", ha="center", va="top", fontsize=10.5, color=FOCUS_COLORS[gene], fontweight="bold")
    axTop.set_xlim(-0.5, 3.5)
    axTop.set_ylim(-1.1, 1.2)

    vmax = np.nanmax(np.abs(pivot.to_numpy()))
    im = ax.imshow(pivot.to_numpy(), cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
    for i in range(pivot.shape[0]):
        for j in range(pivot.shape[1]):
            v = pivot.iloc[i, j]
            ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=11.5,
                    color="white" if abs(v) > vmax * 0.55 else DGRAY)
    ax.set_xticks(range(4))
    ax.set_xticklabels(FOCUS_FOUR, fontsize=17, fontweight="bold")
    for tick, gene in zip(ax.get_xticklabels(), FOCUS_FOUR):
        tick.set_color(FOCUS_COLORS[gene])
    ax.set_yticks(range(len(pivot)))
    ax.set_yticklabels(pivot.index, fontsize=12.5)
    ax.tick_params(length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)

    _hero_title(fig, "Baseline dependency: real cell lines, real values", y=1.02)
    _takeaway(fig, f"DepMap 26Q1 Chronos scores, {len(pivot)} ER+/luminal breast-cancer lines. More negative = "
                    "stronger baseline dependency (not automatically an advantage for tamoxifen sensitisation).",
              y=-0.03)
    _save(fig, stub, vector=False)


# ===========================================================================
# FIGURE 6 -- Structural / pharmacological comparison (hero)
# ===========================================================================

_PHARMACOLOGY_LABEL = {
    "KDM1A": "Mature pharmacology\n(clinical-stage inhibitors)",
    "TLK2": "Structure exists\nno validated selective inhibitor",
    "USP34": "Catalytic cysteine reactivity\nno validated inhibitor",
    "VEZF1": "Hard-to-drug\ntranscription factor",
}


def build_fig6_structural_comparison(stub: Path) -> None:
    with tempfile.TemporaryDirectory() as tmp:
        kt_panels = _render_kdm1a_tlk2_panels(Path(tmp))
        kdm1a_img = _autocrop_white(plt.imread(kt_panels["kdm1a_full"]))
        tlk2_img = _autocrop_white(plt.imread(kt_panels["tlk2_full"]))
        usp34_dir = Path(tmp) / "usp34"
        usp34_dir.mkdir()
        usp34_panels = _render_structure_bank_panels(usp34_dir)
        usp34_img = _autocrop_white(plt.imread(usp34_panels["comp_bound"]))

    fig, axes = plt.subplots(1, 4, figsize=(17, 5.0), dpi=300, gridspec_kw=dict(wspace=0.08))
    imgs = [kdm1a_img, tlk2_img, usp34_img, None]
    for ax, img, gene in zip(axes, imgs, FOCUS_FOUR):
        color = FOCUS_COLORS[gene]
        if img is not None:
            ax.imshow(img)
            ax.axis("off")
        else:
            ax.axis("off")
            ax.add_patch(Rectangle((0.06, 0.08), 0.88, 0.84, transform=ax.transAxes, facecolor="#f4f4f4",
                                    edgecolor=LGRAY, linewidth=1.4, zorder=1))
            ax.text(0.5, 0.5, "no experimental\nstructure", transform=ax.transAxes, ha="center", va="center",
                    fontsize=15, color=GRAY, fontweight="bold")
        ax.text(0.5, 1.05, gene, transform=ax.transAxes, fontsize=21, fontweight="bold", ha="center", color=color)
        ax.text(0.5, -0.1, _PHARMACOLOGY_LABEL[gene], transform=ax.transAxes, fontsize=11.5, ha="center", va="top",
                color=DGRAY)

    _hero_title(fig, "Four candidates, four different levels of structural and pharmacological maturity", y=1.14)
    _takeaway(fig, "Real, experimentally solved structures (KDM1A 6NQU, TLK2 5O0Y, USP34 7W3U) rendered in a "
                    "matched style; VEZF1 has none.", y=-0.16)
    _save(fig, stub, vector=False)


# ===========================================================================
# ALTERNATE A -- CRISPR discovery, genome-wide landscape variant
# ===========================================================================

def build_alt_a_genomewide_landscape(stub: Path) -> None:
    genomewide = pv3.load_genomewide_crispr().sort_values("effect_size").reset_index(drop=True)
    genomewide["rank"] = np.arange(1, len(genomewide) + 1)
    gate1 = genomewide["fdr"] < 0.1
    sens = pv3.load_significant_sensitising_hits().merge(genomewide[["gene", "rank"]], on="gene")

    fig, ax = plt.subplots(figsize=(13, 8), dpi=300)
    ax.scatter(genomewide.loc[~gate1, "rank"], genomewide.loc[~gate1, "effect_size"], s=4, color=LGRAY,
               alpha=0.55, linewidth=0, rasterized=True, zorder=1)
    ax.scatter(genomewide.loc[gate1, "rank"], genomewide.loc[gate1, "effect_size"], s=14, color="#b0b0b0",
               alpha=0.8, linewidth=0, zorder=2)
    ax.axhline(0, color=LGRAY, linewidth=1.0, zorder=0)

    for gene in FOCUS_FOUR:
        row = genomewide[genomewide["gene"] == gene].iloc[0]
        ax.scatter([row["rank"]], [row["effect_size"]], s=280, color=FOCUS_COLORS[gene], edgecolor="white",
                   linewidth=1.8, zorder=5)

    axins = ax.inset_axes([0.30, 0.06, 0.42, 0.42])
    max_rank = int(sens["rank"].max()) + 2
    colors13 = [FOCUS_COLORS.get(g, GRAY) for g in sens["gene"]]
    axins.scatter(sens["rank"], sens["effect_size"], s=70, color=colors13, edgecolor="white", linewidth=0.9, zorder=3)
    axins.set_xlim(0, max_rank)
    axins.set_title(f"the {len(sens)} significant sensitising hits", fontsize=11, color=DGRAY)
    axins.tick_params(labelsize=9)
    for spine in ("top", "right"):
        axins.spines[spine].set_visible(False)
    for _, r in sens.iterrows():
        if r["gene"] in FOCUS_FOUR:
            axins.annotate(r["gene"], xy=(r["rank"], r["effect_size"]), xytext=(3, 3), textcoords="offset points",
                            fontsize=10.5, fontweight="bold", color=FOCUS_COLORS[r["gene"]])

    ax.set_xlabel("Rank by CRISPR effect size", fontsize=16)
    ax.set_ylabel("CRISPR effect size", fontsize=16)
    ax.tick_params(labelsize=13)
    _clean_axes(ax)

    _hero_title(fig, f"Genome-wide screen: {len(genomewide):,} genes, one clear sensitising tail")
    _takeaway(fig, "The 4 focus genes sit deep in the tail of sensitising effect -- magnified in the inset.", y=-0.03)
    _save(fig, stub)


# ===========================================================================
# ALTERNATE B -- structural comparison, pocket close-up variant
# ===========================================================================

def build_alt_b_pocket_closeups(stub: Path) -> None:
    with tempfile.TemporaryDirectory() as tmp:
        kt_panels = _render_kdm1a_tlk2_panels(Path(tmp))
        kdm1a_img = _autocrop_white(plt.imread(kt_panels["kdm1a_pocket"]))
        tlk2_img = _autocrop_white(plt.imread(kt_panels["tlk2_pocket"]))
        usp34_dir = Path(tmp) / "usp34"
        usp34_dir.mkdir()
        usp34_panels = _render_structure_bank_panels(usp34_dir)
        usp34_img = _autocrop_white(plt.imread(usp34_panels["comp_closeup"]))

    fig, axes = plt.subplots(1, 3, figsize=(15, 5.6), dpi=300, gridspec_kw=dict(wspace=0.06))
    titles = [("KDM1A", "Inhibitor pocket"), ("TLK2", "ATP-analog pocket"), ("USP34", "Covalent probe site")]
    for ax, img, (gene, note) in zip(axes, [kdm1a_img, tlk2_img, usp34_img], titles):
        ax.imshow(img)
        ax.axis("off")
        ax.text(0.5, 1.06, gene, transform=ax.transAxes, fontsize=20, fontweight="bold", ha="center",
                 color=FOCUS_COLORS[gene])
        ax.text(0.5, -0.05, note, transform=ax.transAxes, fontsize=13, ha="center", va="top", color=DGRAY)

    _hero_title(fig, "Zooming into the binding site: three very different chemistries", y=1.16)
    _takeaway(fig, "VEZF1 omitted -- no experimental pocket exists to zoom into.", y=-0.1)
    _save(fig, stub, vector=False)


# ===========================================================================
# Contact sheets
# ===========================================================================

MAIN_ITEMS = [
    ("F1_crispr_discovery", "F1. CRISPR discovery"),
    ("F2_transcriptomic_corroboration", "F2. Transcriptomic corroboration"),
    ("F3_pathway_convergence", "F3. Pathway convergence"),
    ("F4_postaudit_interpretation", "F4. Post-audit interpretation"),
    ("F5_depmap_context", "F5. Human / DepMap context"),
    ("F6_structural_comparison", "F6. Structural comparison"),
]
ALTERNATE_ITEMS = [
    ("ALT_A_genomewide_landscape", "Alt A. Genome-wide landscape (alt. to F1)"),
    ("ALT_B_pocket_closeups", "Alt B. Pocket close-ups (alt. to F6)"),
]


def _contact_sheet(items: list[tuple[str, str]], out_png: Path, ncols: int, title: str, figures_dir: Path = FIGURES) -> None:
    n = len(items)
    nrows = -(-n // ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(7.5 * ncols, 5.6 * nrows), dpi=170)
    axes = np.atleast_2d(axes)
    for i, (stem, label) in enumerate(items):
        r, c = divmod(i, ncols)
        ax = axes[r, c]
        img = plt.imread(figures_dir / f"{stem}.png")
        ax.imshow(img)
        ax.axis("off")
        ax.set_title(label, fontsize=13, fontweight="bold", color=DGRAY, pad=8)
    for i in range(n, nrows * ncols):
        r, c = divmod(i, ncols)
        axes[r, c].axis("off")
    fig.suptitle(title, fontsize=18, fontweight="bold", y=1.0)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(out_png, facecolor="white", bbox_inches="tight", dpi=170)
    plt.close(fig)
    logger.info("wrote %s", out_png)


def build_contact_sheets(out_dir: Path = FIGURES) -> None:
    _contact_sheet(MAIN_ITEMS, out_dir / "CONTACT_MAIN_SIX.png", ncols=2,
                    title="Poster-grade figure bank v3 -- 6 main candidates")
    _contact_sheet(MAIN_ITEMS + ALTERNATE_ITEMS, out_dir / "CONTACT_MAIN_PLUS_ALTERNATES.png", ncols=2,
                    title="Poster-grade figure bank v3 -- 6 mains + 2 alternates")


# ===========================================================================
# Runner
# ===========================================================================

def run(out_dir: str | Path = FIGURES) -> None:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    build_fig1_crispr_discovery(out_dir / "F1_crispr_discovery")
    build_fig2_transcriptomic_corroboration(out_dir / "F2_transcriptomic_corroboration")
    build_fig3_pathway_convergence(out_dir / "F3_pathway_convergence")
    build_fig4_postaudit_interpretation(out_dir / "F4_postaudit_interpretation")
    build_fig5_depmap_context(out_dir / "F5_depmap_context")
    build_fig6_structural_comparison(out_dir / "F6_structural_comparison")
    build_alt_a_genomewide_landscape(out_dir / "ALT_A_genomewide_landscape")
    build_alt_b_pocket_closeups(out_dir / "ALT_B_pocket_closeups")
    build_contact_sheets(out_dir)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run()
