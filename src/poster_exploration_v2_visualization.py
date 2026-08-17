"""EXPLORATION-V2 poster figure bank -- visualization only.

Complete reset of the visualization strategy (2026-08-16). Builds ~28
genuinely different exploratory figures across 8 sections (A-H) into
results/figures/poster_exploration_v2/, for later visual selection -- NOT
a final 6-figure poster. See results/reports/poster_exploration_v2/
FIGURE_CANDIDATE_GUIDE.md for the per-figure rationale and
POSTER_LAYOUT... (not built this phase -- no poster is assembled here).

Style system (restrained, "show the data"):
  - white background; near-black/dark-gray text and axes
  - non-focal data (background genes, other cell lines, other candidates)
    is always light neutral gray -- never colored
  - candidate accent color (FOCUS_COLORS) used ONLY where gene identity is
    the point of a specific panel
  - no colored card/box backgrounds, no scorecards, no dashboards
  - bold panel letters (A, B, C) for multi-panel figures; no "F1."/"A1."
    poster-section numbering baked into any image (that lives in the
    filename and the contact sheets only)
  - every figure exports .png (300dpi) + .pdf; .svg where content is
    vector-friendly (not a rasterized PyMOL render)

No new analysis, no recomputation, no re-ranking anywhere in this module.
See src/poster_exploration_v2_data.py for the exact source of every
plotted number.
"""

from __future__ import annotations

import logging
import shutil
import subprocess
import tempfile
import textwrap
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch, Rectangle

from src import poster_exploration_v2_data as pv2
from src.poster_figures_bank_visualization import _render_structure_bank_panels
from src.poster_figures_visualization import _autocrop_white, rng

logger = logging.getLogger(__name__)

FIGURES = Path("results/figures/poster_exploration_v2")

FOCUS_FOUR = pv2.FOCUS_FOUR
FOCUS_COLORS = pv2.FOCUS_COLORS
GRAY = pv2.GRAY
LGRAY = pv2.LGRAY
VLGRAY = pv2.VLGRAY
DGRAY = pv2.DGRAY

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


def _clean_axes(ax) -> None:
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)


def _panel_letter(ax, letter: str, fontsize: float = 15) -> None:
    ax.text(-0.02, 1.06, letter, transform=ax.transAxes, fontsize=fontsize, fontweight="bold", color=DGRAY, va="bottom")


def _title(ax, text: str, fontsize: float = 12) -> None:
    ax.set_title(text, fontsize=fontsize, color=DGRAY, loc="left", pad=8)


# VEZF1 and USP34 sit almost on top of each other on the volcano/hexbin
# axes (effect -1.60 vs -1.39, FDR-log10 ~1.43 vs ~1.38) -- explicit
# per-gene label offsets avoid the text overlap a uniform offset causes.
_VOLCANO_LABEL_OFFSETS = {
    "KDM1A": (6, 5), "TLK2": (6, 5), "VEZF1": (-8, 12), "USP34": (8, -14),
}


def _footnote(fig, text: str, y: float = -0.02, fontsize: float = 7.8, width: int = 130) -> None:
    """Wraps at a fixed character width so long footnotes cannot force
    bbox_inches='tight' to stretch the whole canvas width to fit one
    unwrapped line (observed as a real bug on first render of Figure B4)."""
    wrapped = textwrap.fill(text, width=width)
    fig.text(0.01, y, wrapped, fontsize=fontsize, color=GRAY, style="italic")


# ===========================================================================
# SECTION A -- CRISPR discovery
# ===========================================================================

def build_a1_genomewide_ranked_landscape(stub: Path) -> None:
    df = pv2.load_genomewide_crispr().sort_values("effect_size").reset_index(drop=True)
    df["rank"] = np.arange(1, len(df) + 1)
    gate1 = df["fdr"] < 0.1
    sens13 = pv2.load_significant_sensitising_hits().merge(df[["gene", "rank"]], on="gene")

    fig, ax = plt.subplots(figsize=(11, 7.5), dpi=300)
    ax.scatter(df.loc[~gate1, "rank"], df.loc[~gate1, "effect_size"], s=3, color=LGRAY, alpha=0.5,
               linewidth=0, rasterized=True, zorder=1)
    ax.scatter(df.loc[gate1, "rank"], df.loc[gate1, "effect_size"], s=10, color=GRAY, alpha=0.7,
               linewidth=0, zorder=2)
    ax.axhline(0, color=LGRAY, linewidth=0.8, zorder=0)

    for gene in FOCUS_FOUR:
        row = df[df["gene"] == gene].iloc[0]
        ax.scatter([row["rank"]], [row["effect_size"]], s=55, color=FOCUS_COLORS[gene],
                   edgecolor="white", linewidth=0.8, zorder=5)

    ax.set_xlabel("Rank by CRISPR effect size (most sensitising -> most tolerance-associated)", fontsize=11)
    ax.set_ylabel("CRISPR effect size", fontsize=11)
    ax.set_xlim(-300, len(df) + 300)
    _clean_axes(ax)
    _title(ax, f"Genome-wide CRISPR screen, ranked by effect (n={len(df):,} genes)")

    # zoomed inset -- the significant sensitising tail (rank 1-13), where
    # the 4 focus genes actually sit, magnified since they are otherwise
    # visually indistinguishable at genome-wide scale
    axins = ax.inset_axes([0.32, 0.08, 0.36, 0.38])
    max_rank = int(sens13["rank"].max()) + 2
    colors13 = [FOCUS_COLORS.get(g, GRAY) for g in sens13["gene"]]
    axins.scatter(sens13["rank"], sens13["effect_size"], s=45, color=colors13, edgecolor="white", linewidth=0.6, zorder=3)
    axins.set_xlim(0, max_rank)
    axins.set_title(f"zoom: rank 1-{max_rank}\n(13 sensitising hits)", fontsize=7.6, color=DGRAY)
    axins.tick_params(labelsize=6.5)
    for spine in ("top", "right"):
        axins.spines[spine].set_visible(False)
    for _, r in sens13.iterrows():
        if r["gene"] in FOCUS_FOUR:
            axins.annotate(r["gene"], xy=(r["rank"], r["effect_size"]), xytext=(2, 2), textcoords="offset points",
                            fontsize=6.3, fontweight="bold", color=FOCUS_COLORS[r["gene"]])

    _footnote(fig, f"{len(df):,} fitted genes -> {int(gate1.sum())} Gate-1 significant hits (FDR<0.1) -> "
                    f"{len(sens13)} sensitising hits -> 4 focus genes highlighted (gray = other 9).", y=-0.04)
    _save(fig, stub)


def build_a2_lollipop_13_hits(stub: Path) -> None:
    sens = pv2.load_significant_sensitising_hits().sort_values("effect_size").reset_index(drop=True)
    y = np.arange(len(sens))
    colors = [FOCUS_COLORS.get(g, GRAY) for g in sens["gene"]]
    sizes = 40 + 260 * (1 - sens["fdr"].clip(upper=0.1) / 0.1)

    fig, ax = plt.subplots(figsize=(8, 7.5), dpi=300)
    ax.hlines(y, 0, sens["effect_size"], color=colors, linewidth=1.6, alpha=0.55, zorder=2)
    ax.scatter(sens["effect_size"], y, s=sizes, color=colors, edgecolor="white", linewidth=0.8, zorder=3)
    ax.axvline(0, color=LGRAY, linewidth=0.8, zorder=0)
    ax.set_yticks(y)
    labels = [f"$\\bf{{{g}}}$" if g in FOCUS_FOUR else g for g in sens["gene"]]
    ax.set_yticklabels(labels, fontsize=10.5)
    for tick, gene in zip(ax.get_yticklabels(), sens["gene"]):
        tick.set_color(FOCUS_COLORS.get(gene, DGRAY))
    ax.set_xlabel("CRISPR effect size (all 13 significant, FDR<0.1, sensitising)", fontsize=10.5)
    _clean_axes(ax)
    _title(ax, "13 significant sensitising CRISPR hits")

    handles = [Line2D([0], [0], marker="o", color="none", markerfacecolor=GRAY, markersize=6, label="FDR near 0.1 (weaker)"),
               Line2D([0], [0], marker="o", color="none", markerfacecolor=GRAY, markersize=11, label="FDR near 0 (stronger)")]
    ax.legend(handles=handles, loc="lower right", fontsize=8, frameon=False, title="point size = FDR", title_fontsize=8)
    _footnote(fig, "KDM1A and TLK2 are functionally stronger than USP34 and VEZF1 by this measure alone; "
                    "USP34/VEZF1 were retained for other evidence (see evidence-intersection figures).", y=-0.04)
    _save(fig, stub)


def build_a3_redesigned_volcano(stub: Path) -> None:
    df = pv2.load_genomewide_crispr()
    x = df["effect_size"].to_numpy()
    y = -np.log10(df["fdr"].to_numpy())
    gate1 = df["fdr"] < 0.1

    fig, ax = plt.subplots(figsize=(8, 7), dpi=300)
    ax.scatter(x[~gate1], y[~gate1], s=3, color=LGRAY, alpha=0.45, linewidth=0, rasterized=True, zorder=1)
    ax.scatter(x[gate1], y[gate1], s=9, color=GRAY, alpha=0.7, linewidth=0, zorder=2)
    ax.axhline(-np.log10(0.1), color=LGRAY, linewidth=0.7, linestyle="--", zorder=0)
    for gene in FOCUS_FOUR:
        row = df[df["gene"] == gene].iloc[0]
        gx, gy = float(row["effect_size"]), -np.log10(float(row["fdr"]))
        ax.scatter([gx], [gy], s=95, color=FOCUS_COLORS[gene], edgecolor="white", linewidth=0.9, zorder=5)
        ax.annotate(gene, xy=(gx, gy), xytext=_VOLCANO_LABEL_OFFSETS[gene], textcoords="offset points", fontsize=9.5,
                    fontweight="bold", color=FOCUS_COLORS[gene],
                    ha="right" if _VOLCANO_LABEL_OFFSETS[gene][0] < 0 else "left")
    ax.set_xlabel("CRISPR effect size", fontsize=10.5)
    ax.set_ylabel(r"$-\log_{10}$(FDR)", fontsize=10.5)
    _clean_axes(ax)
    _title(ax, "Genome-wide volcano")
    _save(fig, stub)


def build_a4_hexbin_density(stub: Path) -> None:
    """Density view of the genome-wide screen -- a hexbin of the full
    19,103-gene cloud with the 4 focus genes and RCOR1 overlaid, showing
    how far into the tail the focus genes sit relative to the bulk
    distribution (a different visual than the scatter-based A1/A3)."""
    df = pv2.load_genomewide_crispr()
    x = df["effect_size"].to_numpy()
    y = -np.log10(df["fdr"].to_numpy())

    fig, ax = plt.subplots(figsize=(8, 7), dpi=300)
    hb = ax.hexbin(x, y, gridsize=55, cmap="Greys", mincnt=1, linewidths=0.15, edgecolors=VLGRAY)
    cb = fig.colorbar(hb, ax=ax, fraction=0.045, pad=0.02)
    cb.set_label("gene density", fontsize=9)
    cb.ax.tick_params(labelsize=8)

    for gene in FOCUS_FOUR:
        row = df[df["gene"] == gene].iloc[0]
        gx, gy = float(row["effect_size"]), -np.log10(float(row["fdr"]))
        ax.scatter([gx], [gy], s=110, color=FOCUS_COLORS[gene], edgecolor="white", linewidth=1.1, zorder=5)
        ax.annotate(gene, xy=(gx, gy), xytext=_VOLCANO_LABEL_OFFSETS[gene], textcoords="offset points", fontsize=9.5,
                    fontweight="bold", color=FOCUS_COLORS[gene],
                    ha="right" if _VOLCANO_LABEL_OFFSETS[gene][0] < 0 else "left")
    rcor1 = pv2.load_blind_control_row()
    if rcor1 is not None:
        rx, ry = float(rcor1["effect_size"]), -np.log10(float(rcor1["fdr"]))
        ax.scatter([rx], [ry], s=80, facecolor="none", edgecolor=DGRAY, linewidth=1.3, zorder=5)
        ax.annotate("RCOR1", xy=(rx, ry), xytext=(6, -12), textcoords="offset points", fontsize=8.5, color=DGRAY)

    ax.set_xlabel("CRISPR effect size", fontsize=10.5)
    ax.set_ylabel(r"$-\log_{10}$(FDR)", fontsize=10.5)
    _clean_axes(ax)
    _title(ax, "Screen density with focus genes in the tail")
    _save(fig, stub, vector=False)


def run_section_a(out_dir: Path = FIGURES) -> None:
    build_a1_genomewide_ranked_landscape(out_dir / "A1_genomewide_ranked_landscape")
    build_a2_lollipop_13_hits(out_dir / "A2_lollipop_13_hits")
    build_a3_redesigned_volcano(out_dir / "A3_redesigned_volcano")
    build_a4_hexbin_density(out_dir / "A4_hexbin_density")


# ===========================================================================
# SECTION B -- bulk / single-cell transcriptomics (real sample-level data)
# ===========================================================================

def _strip_panel(ax, groups: list[str], values_by_group: dict[str, np.ndarray], color: str, width: float = 0.32) -> None:
    for i, g in enumerate(groups):
        vals = values_by_group[g]
        if len(vals) == 0:
            continue
        mean = float(np.mean(vals))
        ax.plot([i - width / 2, i + width / 2], [mean, mean], color=color, linewidth=2.0, zorder=3)
        jitter = rng.uniform(-width * 0.35, width * 0.35, size=len(vals))
        ax.scatter(np.full(len(vals), i) + jitter, vals, s=42, color=color, edgecolor="white",
                   linewidth=0.5, zorder=4, alpha=0.9)
    ax.set_xticks(range(len(groups)))
    ax.set_xticklabels(groups, fontsize=8.7)


def build_b1_gse118713_sample_panel(stub: Path) -> None:
    df = pv2.load_gse118713_focus_gene_samples()
    order = ["MCF7", "TAMR", "FASR"]
    fig, axes = plt.subplots(2, 2, figsize=(9, 7.5), dpi=300, sharex=True)
    for ax, gene in zip(axes.flat, FOCUS_FOUR):
        sub = df[df["gene_symbol"] == gene]
        groups = {c: sub.loc[sub["condition"] == c, "tpm"].to_numpy() for c in order}
        _strip_panel(ax, order, groups, FOCUS_COLORS[gene])
        ax.set_ylabel("TPM", fontsize=9)
        _clean_axes(ax)
        _title(ax, gene, fontsize=11)
    fig.suptitle("GSE118713 -- real per-sample expression (3 replicates per condition)", fontsize=12.5, x=0.02, ha="left", y=1.01)
    _footnote(fig, "One parental line -> one TAMR-derivation + one FASR-derivation event; 3 replicate aliquots "
                    "each, not 3 independent derivations.", y=-0.03)
    _save(fig, stub)


def build_b2_gse111151_trajectories(stub: Path) -> None:
    df = pv2.load_gse111151_focus_gene_samples()
    parental_order = ["MCF-7", "T-47D", "ZR-75-1", "BT-474"]
    fig, axes = plt.subplots(2, 2, figsize=(9.5, 8), dpi=300, sharex=True)
    for ax, gene in zip(axes.flat, FOCUS_FOUR):
        sub = df[df["gene_symbol"] == gene]
        for i, pline in enumerate(parental_order):
            prow = sub[(sub["parental_line"] == pline) & (sub["status"] == "parental")]
            if len(prow) == 0:
                continue
            pval = float(prow["log2cpm"].iloc[0])
            ax.scatter([i], [pval], s=60, facecolor="white", edgecolor=DGRAY, linewidth=1.4, zorder=4)
            derivs = sub[(sub["parental_line"] == pline) & (sub["status"] == "resistant")]
            for _, drow in derivs.iterrows():
                dval = float(drow["log2cpm"])
                ax.plot([i, i], [pval, dval], color=FOCUS_COLORS[gene], linewidth=1.3, alpha=0.75, zorder=2)
                ax.scatter([i], [dval], s=55, color=FOCUS_COLORS[gene], edgecolor="white", linewidth=0.6, zorder=3)
        ax.set_xticks(range(len(parental_order)))
        ax.set_xticklabels(parental_order, fontsize=8.3, rotation=20, ha="right")
        ax.set_ylabel("log2(CPM)", fontsize=9)
        _clean_axes(ax)
        _title(ax, gene, fontsize=11)
    handles = [Line2D([0], [0], marker="o", color="none", markerfacecolor="white", markeredgecolor=DGRAY, markersize=8, label="parental"),
               Line2D([0], [0], marker="o", color="none", markerfacecolor=GRAY, markersize=8, label="TamR derivative")]
    fig.legend(handles=handles, loc="lower center", ncol=2, fontsize=8.5, frameon=False, bbox_to_anchor=(0.5, -0.02))
    fig.suptitle("GSE111151 -- real parental -> TamR-derivative trajectories (4 backgrounds, 7 sublines)", fontsize=12.5, x=0.02, ha="left", y=1.02)
    _footnote(fig, "BT-474 is HER2-amplified, a different molecular subtype from the other three parental lines. "
                    "Vertical connectors reflect the real recorded parental->derivative pairing, not an invented one.", y=-0.06)
    _save(fig, stub)


def build_b3_gse240112_tumours(stub: Path) -> None:
    df = pv2.load_gse240112_focus_gene_tumours()
    order = ["PT", "RT"]
    order_labels = {"PT": "primary (n=3)", "RT": "recurrent (n=3)"}
    fig, axes = plt.subplots(2, 2, figsize=(8.5, 7.5), dpi=300, sharex=True)
    for ax, gene in zip(axes.flat, FOCUS_FOUR):
        sub = df[df["gene"] == gene]
        for i, grp in enumerate(order):
            vals = sub.loc[sub["group"] == grp, "log2cpm"].to_numpy()
            mean = float(np.mean(vals))
            ax.plot([i - 0.18, i + 0.18], [mean, mean], color=FOCUS_COLORS[gene], linewidth=2.0, zorder=3)
            jitter = rng.uniform(-0.1, 0.1, size=len(vals))
            ax.scatter(np.full(len(vals), i) + jitter, vals, s=70, color=FOCUS_COLORS[gene],
                       edgecolor="white", linewidth=0.7, zorder=4)
        ax.set_xticks(range(len(order)))
        ax.set_xticklabels([order_labels[o] for o in order], fontsize=8.7)
        ax.set_xlim(-0.5, 1.5)
        ax.set_ylabel("log2(CPM)", fontsize=9)
        _clean_axes(ax)
        weight = "bold" if gene == "VEZF1" else "normal"
        ax.set_title(gene, fontsize=11.5 if gene == "VEZF1" else 11, color=DGRAY, loc="left", pad=8, fontweight=weight)
    fig.suptitle("GSE240112 -- real per-tumour values, 3 primary vs 3 recurrent (UNPAIRED)", fontsize=12.5, x=0.02, ha="left", y=1.01)
    _footnote(fig, "Different patients, different biobanks -- no primary observation corresponds to any recurrent "
                    "observation. No connecting line is drawn. VEZF1 (bold) has this project's strongest recurrence-"
                    "associated RNA signal among the 4 focus genes.", y=-0.04)
    _save(fig, stub)


def build_b4_gse245601_acute_paired(stub: Path) -> None:
    df = pv2.load_gse245601_paired_focus_genes()
    patients = sorted(df["patient"].unique())
    conditions = ["Control", "Tamoxifen"]
    fig, axes = plt.subplots(2, 2, figsize=(8.5, 7.5), dpi=300, sharex=True)
    for ax, gene in zip(axes.flat, FOCUS_FOUR):
        sub = df[df["gene"] == gene]
        for patient in patients:
            prow = sub[sub["patient"] == patient].set_index("condition")
            vals = [float(prow.loc[c, "log2_expr"]) for c in conditions]
            ax.plot([0, 1], vals, color=FOCUS_COLORS[gene], linewidth=1.5, alpha=0.8, marker="o",
                    markersize=6.5, markeredgecolor="white", markeredgewidth=0.6, zorder=3)
        ax.set_xticks([0, 1])
        ax.set_xticklabels(conditions, fontsize=9)
        ax.set_xlim(-0.35, 1.35)
        ax.set_ylabel("log2(CPM+1)", fontsize=9)
        _clean_axes(ax)
        _title(ax, gene, fontsize=11)
    fig.suptitle(f"GSE245601 -- real patient-matched acute 12h response ({len(patients)} eligible patients)", fontsize=12.5, x=0.02, ha="left", y=1.01)
    _footnote(fig, "ACUTE 12h ex vivo tamoxifen exposure, NOT chronic resistance. Each line = one patient's own "
                    "tumour, Control -> Tamoxifen (real pairing, connecting lines appropriate here unlike GSE240112). "
                    "KDM1A/TLK2 values computed from frozen raw counts using the same formula as the frozen USP34/"
                    "VEZF1 values (see DATA_FOR_VISUALIZATION_AUDIT.md).", y=-0.05)
    _save(fig, stub)


def build_b5a_integrated_dataset_centric(stub: Path) -> None:
    """Dataset-centric integrated overview: one row per dataset, one small
    panel per focus gene within that row, all real values (condensed
    mean+spread strips, not binary yes/no marks)."""
    g118713 = pv2.load_gse118713_focus_gene_samples()
    g111151 = pv2.load_gse111151_focus_gene_samples()
    g240112 = pv2.load_gse240112_focus_gene_tumours()
    g245601 = pv2.load_gse245601_paired_focus_genes()

    fig, axes = plt.subplots(4, 4, figsize=(13, 10), dpi=300)
    row_specs = [
        ("GSE118713\n(resistance)", g118713, "condition", "tpm", ["MCF7", "TAMR", "FASR"], "gene_symbol"),
        ("GSE111151\n(resistance)", g111151, "status", "log2cpm", ["parental", "resistant"], "gene_symbol"),
        ("GSE240112\n(recurrence)", g240112, "group", "log2cpm", ["PT", "RT"], "gene"),
        ("GSE245601\n(acute)", g245601, "condition", "log2_expr", ["Control", "Tamoxifen"], "gene"),
    ]
    for r, (row_label, data, group_col, val_col, groups, gene_col) in enumerate(row_specs):
        for c, gene in enumerate(FOCUS_FOUR):
            ax = axes[r, c]
            sub = data[data[gene_col] == gene]
            vals_by_group = {g: sub.loc[sub[group_col] == g, val_col].to_numpy() for g in groups}
            _strip_panel(ax, groups, vals_by_group, FOCUS_COLORS[gene], width=0.4)
            ax.tick_params(labelsize=7)
            for spine in ("top", "right"):
                ax.spines[spine].set_visible(False)
            if c == 0:
                ax.set_ylabel(row_label, fontsize=9, color=DGRAY)
            if r == 0:
                ax.set_title(gene, fontsize=11, color=FOCUS_COLORS[gene], fontweight="bold")
    fig.suptitle("Dataset-centric: real values, 4 datasets (rows) x 4 focus genes (columns)", fontsize=13, x=0.02, ha="left", y=1.0)
    _footnote(fig, "Same underlying observations as Figures B1-B4, reorganized to compare across datasets within "
                    "one row at a time. GSE118713/GSE111151 = resistance models; GSE240112 = recurrence-associated "
                    "unpaired human tumours; GSE245601 = acute 12h, patient-paired.", y=-0.02)
    _save(fig, stub)


def build_b5b_integrated_gene_centric(stub: Path) -> None:
    """Gene-centric integrated overview: one row per focus gene, one panel
    per dataset -- the same data as B5a, transposed, to make each gene's
    story legible in a single row."""
    g118713 = pv2.load_gse118713_focus_gene_samples()
    g111151 = pv2.load_gse111151_focus_gene_samples()
    g240112 = pv2.load_gse240112_focus_gene_tumours()
    g245601 = pv2.load_gse245601_paired_focus_genes()

    col_specs = [
        ("GSE118713", g118713, "condition", "tpm", ["MCF7", "TAMR", "FASR"], "gene_symbol"),
        ("GSE111151", g111151, "status", "log2cpm", ["parental", "resistant"], "gene_symbol"),
        ("GSE240112", g240112, "group", "log2cpm", ["PT", "RT"], "gene"),
        ("GSE245601", g245601, "condition", "log2_expr", ["Control", "Tamoxifen"], "gene"),
    ]
    fig, axes = plt.subplots(4, 4, figsize=(13, 10), dpi=300)
    for r, gene in enumerate(FOCUS_FOUR):
        for c, (col_label, data, group_col, val_col, groups, gene_col) in enumerate(col_specs):
            ax = axes[r, c]
            sub = data[data[gene_col] == gene]
            vals_by_group = {g: sub.loc[sub[group_col] == g, val_col].to_numpy() for g in groups}
            _strip_panel(ax, groups, vals_by_group, FOCUS_COLORS[gene], width=0.4)
            ax.tick_params(labelsize=7)
            for spine in ("top", "right"):
                ax.spines[spine].set_visible(False)
            if c == 0:
                ax.set_ylabel(gene, fontsize=11, color=FOCUS_COLORS[gene], fontweight="bold")
            if r == 0:
                ax.set_title(col_label, fontsize=9.5, color=DGRAY)
    fig.suptitle("Gene-centric: same real values, reorganized one row per focus gene", fontsize=13, x=0.02, ha="left", y=1.0)
    _footnote(fig, "USP34's only significant RNA corroboration is GSE118713; VEZF1's is GSE240112; KDM1A and TLK2 "
                    "have no significant RNA corroboration in any of these 4 datasets despite the strongest CRISPR "
                    "signal (Section A) -- an honest, visible null, not hidden.", y=-0.02)
    _save(fig, stub)


def build_c1a_pathway_trajectories_lines(stub: Path) -> None:
    pathways = pv2.HERO_PATHWAYS + pv2.EXTRA_PATHWAYS
    df = pv2.load_pathway_trajectories(pathways)
    labels = [p[2] for p in pathways]
    x = np.arange(len(pv2.DATASET_ORDER))
    cmap = plt.get_cmap("tab10")

    fig, ax = plt.subplots(figsize=(9.5, 7), dpi=300)
    ax.axvspan(2.5, 3.5, color=VLGRAY, zorder=0)
    ax.axhline(0, color=LGRAY, linewidth=0.8, zorder=0)
    endpoints = []
    for i, label in enumerate(labels):
        sub = df[df["pathway_label"] == label].set_index("dataset").reindex(pv2.DATASET_ORDER)
        color = cmap(i % 10)
        ax.plot(x, sub["NES"], color=color, linewidth=1.8, marker="o", markersize=5.5, alpha=0.9, zorder=3)
        last_valid = sub["NES"].last_valid_index()
        if last_valid is not None:
            xi = pv2.DATASET_ORDER.index(last_valid)
            endpoints.append([label, color, xi, float(sub.loc[last_valid, "NES"])])

    # greedy vertical de-collision: enforce a minimum label gap in data
    # units among endpoints that land close together (observed collision
    # among EMT / E2F targets / Estrogen response (late), all ~-1.5 to -1.6).
    # The REAL endpoint (real_y) is kept separate from the dodged label
    # position (label_y) and a thin leader line is drawn between them
    # whenever they diverge, so a label can never appear to belong to the
    # wrong trajectory (a real bug in an earlier version of this figure,
    # caught by adversarial review 2026-08-16).
    for e in endpoints:
        e.append(e[3])  # label_y, initialized to real_y
    endpoints.sort(key=lambda e: e[3])
    min_gap = 0.22
    for j in range(1, len(endpoints)):
        if endpoints[j][4] - endpoints[j - 1][4] < min_gap:
            endpoints[j][4] = endpoints[j - 1][4] + min_gap
    for label, color, xi, real_y, label_y in endpoints:
        if abs(label_y - real_y) > 1e-9:
            ax.plot([xi, xi + 0.18], [real_y, label_y], color=color, linewidth=0.7, alpha=0.6, zorder=2)
        ax.annotate(label, xy=(xi + 0.18, label_y), xytext=(8, 0), textcoords="offset points",
                    fontsize=8.3, color=color, va="center")
    ax.set_xticks(x)
    ax.set_xticklabels([pv2.DATASET_LABELS[d] for d in pv2.DATASET_ORDER], fontsize=9)
    ax.set_ylabel("NES", fontsize=10.5)
    ax.set_xlim(-0.3, 4.3)
    _clean_axes(ax)
    _title(ax, "Pathway NES trajectories across transcriptomic contexts")
    _footnote(fig, "Shaded column = acute 12h context (GSE245601), kept visually separate from the two resistance-"
                    "model and one recurrence-associated contexts. NES magnitude is not directly comparable across "
                    "pathways of different gene-set size -- compare each pathway's own trend across datasets, not "
                    "NES values between different pathways.", y=-0.05)
    _save(fig, stub)


def build_c1b_pathway_small_multiples(stub: Path) -> None:
    pathways = pv2.HERO_PATHWAYS + pv2.EXTRA_PATHWAYS
    df = pv2.load_pathway_trajectories(pathways)
    labels = [p[2] for p in pathways]
    x = np.arange(len(pv2.DATASET_ORDER))

    fig, axes = plt.subplots(len(labels), 1, figsize=(6.5, 11), dpi=300, sharex=True)
    for ax, label in zip(axes, labels):
        sub = df[df["pathway_label"] == label].set_index("dataset").reindex(pv2.DATASET_ORDER)
        colors = ["#4C4C9D"] * 2 + ["#B0392F", GRAY]
        ax.axhline(0, color=LGRAY, linewidth=0.7, zorder=0)
        ax.bar(x, sub["NES"], color=colors, width=0.6, zorder=2)
        for xi, (nes, fdr) in enumerate(zip(sub["NES"], sub["fdr"])):
            if pd.notna(fdr) and fdr < 0.05:
                ax.text(xi, nes + (0.15 if nes >= 0 else -0.15), "*", ha="center",
                        va="bottom" if nes >= 0 else "top", fontsize=11, color=DGRAY)
        ax.set_ylabel(label, fontsize=8.7, rotation=0, ha="right", va="center")
        for spine in ("top", "right", "left"):
            ax.spines[spine].set_visible(False)
        ax.tick_params(left=False, labelleft=False)
    axes[-1].set_xticks(x)
    axes[-1].set_xticklabels([pv2.DATASET_LABELS[d] for d in pv2.DATASET_ORDER], fontsize=8.7)
    fig.suptitle("Pathway NES, one row per pathway (blue=resistance, red=recurrence, gray=acute)", fontsize=12, x=0.02, ha="left", y=0.995)
    _footnote(fig, "Each row's y-axis autoscales independently -- compare a pathway's own trend across the 4 "
                    "columns, never magnitudes between different rows (NES is not comparable across differently "
                    "sized gene sets).", y=-0.01)
    _save(fig, stub)


def build_c2_estrogen_emt_hero(stub: Path) -> None:
    pathways = [("hallmark", "HALLMARK_ESTROGEN_RESPONSE_EARLY", "Estrogen response (early)"),
                ("hallmark", "HALLMARK_ESTROGEN_RESPONSE_LATE", "Estrogen response (late)")]
    df_er = pv2.load_pathway_trajectories(pathways)
    df_emt = pv2.load_pathway_trajectories([("hallmark", "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "EMT")])
    x = np.arange(len(pv2.DATASET_ORDER))
    ds_labels = [d.upper() for d in pv2.DATASET_ORDER]

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(12, 6), dpi=300)
    for label, color, marker in [("Estrogen response (early)", "#08519c", "o"), ("Estrogen response (late)", "#6baed6", "s")]:
        sub = df_er[df_er["pathway_label"] == label].set_index("dataset").reindex(pv2.DATASET_ORDER)
        axA.plot(x, sub["NES"], color=color, linewidth=2.4, marker=marker, markersize=9, label=label, zorder=3)
    axA.axhline(0, color=LGRAY, linewidth=0.8, zorder=0)
    axA.set_xticks(x)
    axA.set_xticklabels(ds_labels, fontsize=9)
    axA.set_ylabel("NES", fontsize=11)
    axA.legend(loc="lower right", fontsize=9, frameon=False)
    _clean_axes(axA)
    _title(axA, "Estrogen response: DOWN in every context", fontsize=13)

    sub = df_emt.set_index("dataset").reindex(pv2.DATASET_ORDER)
    colors = ["#b0392f"] * 3 + [GRAY]
    axB.axhline(0, color=LGRAY, linewidth=0.8, zorder=0)
    axB.bar(x, sub["NES"], color=colors, width=0.55, zorder=2)
    for xi, (nes, fdr) in enumerate(zip(sub["NES"], sub["fdr"])):
        marker = "*" if (pd.notna(fdr) and fdr < 0.05) else ""
        axB.text(xi, nes + (0.08 if nes >= 0 else -0.08), marker, ha="center",
                 va="bottom" if nes >= 0 else "top", fontsize=13, color=DGRAY)
    axB.set_xticks(x)
    axB.set_xticklabels(ds_labels, fontsize=9)
    axB.set_ylabel("NES", fontsize=11)
    _clean_axes(axB)
    _title(axB, "EMT: UP in resistance/recurrence, DOWN in acute 12h", fontsize=13)

    _footnote(fig, "Both panels: real NES from the already-frozen per-dataset GSEA tables. * = FDR<0.05. "
                    "GSE118713/GSE111151 = resistance models; GSE240112 = recurrence-associated (human); "
                    "GSE245601 = acute 12h (gray). NES is not directly comparable across gene sets of different "
                    "size (e.g. early vs. late estrogen response) -- compare each pathway's own trend across "
                    "datasets.", y=-0.04)
    _save(fig, stub)


def build_c3_enrichment_curves(stub: Path) -> None:
    combos = [
        ("gse118713", "HALLMARK_ESTROGEN_RESPONSE_EARLY", "GSE118713: Estrogen response (early)"),
        ("gse245601", "HALLMARK_ESTROGEN_RESPONSE_EARLY", "GSE245601 (acute): Estrogen response (early)"),
        ("gse111151", "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "GSE111151: EMT"),
        ("gse245601", "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "GSE245601 (acute): EMT"),
    ]
    fig, axes = plt.subplots(2, 2, figsize=(11, 9), dpi=300, gridspec_kw=dict(hspace=0.5, wspace=0.3))
    for ax, (dataset, pathway, title) in zip(axes.flat, combos):
        curve = pv2.build_enrichment_curve(dataset, pathway)
        ax.plot(curve["rank"], curve["running_es"], color="#08519c", linewidth=1.8, zorder=3)
        ax.axhline(0, color=LGRAY, linewidth=0.7, zorder=0)
        hit_ranks = curve.loc[curve["in_gene_set"], "rank"].to_numpy()
        ax.vlines(hit_ranks, ymin=-0.02, ymax=0.02, color=DGRAY, linewidth=0.35, alpha=0.5, zorder=2)
        peak_idx = curve["running_es"].abs().idxmax()
        ax.scatter([curve.loc[peak_idx, "rank"]], [curve.loc[peak_idx, "running_es"]], s=35, color="#b0392f", zorder=4)
        ax.set_xlabel("gene rank", fontsize=9)
        ax.set_ylabel("running ES", fontsize=9)
        _clean_axes(ax)
        _title(ax, title, fontsize=10.5)
    _footnote(fig, "Real running-sum GSEA enrichment-score curves, reconstructed from the already-frozen ranked gene "
                    "list + Hallmark gene-set membership (same statistic the frozen NES/FDR already summarize -- see "
                    "DATA_FOR_VISUALIZATION_AUDIT.md). Tick marks = gene-set members' rank positions. Red dot = peak "
                    "deviation (the enrichment score).", y=-0.03)
    _save(fig, stub)


# ===========================================================================
# SECTION D -- network / systems biology (original-4-candidate coverage only)
# ===========================================================================

_ORIGINAL_FOUR_COLORS = {
    "USP34": FOCUS_COLORS["USP34"], "VEZF1": FOCUS_COLORS["VEZF1"],
    "EML5": "#6b6b6b", "CITED2": "#a35b9e",
}


def build_d1_candidate_program_network(stub: Path) -> None:
    """Bipartite candidate -> STRONG_CONSENSUS-pathway network. Only the
    ORIGINAL 4 candidates have frozen network data (KDM1A/TLK2 were never
    analyzed in the systems-network phase -- see
    DATA_FOR_VISUALIZATION_AUDIT.md). CITED2's 49 real STRONG_CONSENSUS
    memberships are aggregated into one labeled summary node (its true
    count is stated, not hidden) to keep the figure legible; USP34/VEZF1's
    small real counts (4 and 2) are shown as individual pathway nodes;
    EML5's real zero-membership is shown explicitly, not omitted."""
    membership = pv2.load_candidate_pathway_membership()
    candidates = ["USP34", "VEZF1", "EML5", "CITED2"]

    def short_label(p: str) -> str:
        return p.split("_", 1)[1].replace("_", " ").title() if "_" in p else p

    fig, ax = plt.subplots(figsize=(11, 8), dpi=300)
    cand_y = {c: i for i, c in enumerate(candidates)}
    cand_x = 0.0
    pathway_x = 1.0

    pathway_nodes = []  # (label, y, candidate)
    y_cursor = 0.0
    row_h = 0.62
    for c in candidates:
        sub = membership[membership["candidate"] == c]
        if c == "CITED2":
            label = f"{len(sub)} STRONG_CONSENSUS programs\n(embryonic / vascular / skeletal\ndevelopment themes)"
            pathway_nodes.append((label, y_cursor, c))
            y_cursor -= row_h * 2.2
        elif len(sub) == 0:
            pathway_nodes.append(("(no STRONG_CONSENSUS\npathway membership)", y_cursor, c))
            y_cursor -= row_h * 1.3
        else:
            for _, r in sub.iterrows():
                pathway_nodes.append((short_label(r["pathway"]), y_cursor, c))
                y_cursor -= row_h

    total_h = -y_cursor
    cand_ys = {c: total_h - (i + 0.5) * (total_h / len(candidates)) for i, c in enumerate(candidates)}

    for label, y, c in pathway_nodes:
        color = _ORIGINAL_FOUR_COLORS[c]
        is_summary = "STRONG_CONSENSUS" in label or "no STRONG" in label
        ax.plot([cand_x, pathway_x], [cand_ys[c], y], color=color, linewidth=1.1, alpha=0.55, zorder=1)
        node_color = "white" if "no STRONG" in label else color
        edge_color = GRAY if "no STRONG" in label else color
        ax.scatter([pathway_x], [y], s=170 if is_summary else 90, color=node_color, edgecolor=edge_color,
                   linewidth=1.3, zorder=3)
        ax.text(pathway_x + 0.04, y, label, fontsize=7.6 if not is_summary else 8.3, color=DGRAY,
                va="center", ha="left")

    for c, y in cand_ys.items():
        ax.scatter([cand_x], [y], s=280, color=_ORIGINAL_FOUR_COLORS[c], edgecolor="white", linewidth=1.4, zorder=4)
        ax.text(cand_x - 0.05, y, c, fontsize=12, fontweight="bold", color=_ORIGINAL_FOUR_COLORS[c],
                va="center", ha="right")

    ax.set_xlim(-0.55, 2.1)
    ax.set_ylim(y_cursor - row_h, total_h + row_h)
    ax.axis("off")
    _title(ax, "Candidate -> pathway-program network (original 4 candidates only)")
    _footnote(fig, "KDM1A and TLK2 have no frozen network analysis in this project (the systems-network phase ran "
                    "only on the original 4-candidate set) and cannot be shown here -- a real coverage gap, not a "
                    "redesigned network. USP34-VEZF1 share no direct interaction and no STRONG_CONSENSUS pathway "
                    "(any_convergence=False in the frozen table).", y=-0.02)
    _save(fig, stub)


def build_d2_usp34_local_neighborhood(stub: Path) -> None:
    """USP34's real 1-hop STRING direct-neighbor table (10 rows) -- USP34
    and CITED2 have the richest neighbor coverage of the original 4
    candidates (VEZF1/EML5 have only 1 each); shown in full, not a hand-
    picked top-N, to avoid silently discarding real evidence;
    low-confidence edges are drawn lighter rather than omitted."""
    nb = pv2.load_direct_neighbors("USP34").sort_values("confidence", ascending=False).reset_index(drop=True)
    n = len(nb)
    angles = np.linspace(0, 2 * np.pi, n, endpoint=False)
    radius = 1.0

    fig, ax = plt.subplots(figsize=(9, 9), dpi=300)
    ax.scatter([0], [0], s=420, color=FOCUS_COLORS["USP34"], edgecolor="white", linewidth=1.6, zorder=5)
    ax.text(0, -0.16, "USP34", ha="center", va="top", fontsize=12, fontweight="bold",
            color=FOCUS_COLORS["USP34"], zorder=6)

    for angle, (_, r) in zip(angles, nb.iterrows()):
        px, py = radius * np.cos(angle), radius * np.sin(angle)
        conf = float(r["confidence"]) if pd.notna(r["confidence"]) else 0.5
        ax.plot([0, px], [0, py], color=GRAY, linewidth=0.5 + 2.2 * conf, alpha=0.35 + 0.5 * conf, zorder=1)
        is_sig = pd.notna(r["neighbor_crispr_fdr"]) and float(r["neighbor_crispr_fdr"]) < 0.1
        color = "#b0392f" if is_sig else GRAY
        ax.scatter([px], [py], s=55, color=color, edgecolor="white", linewidth=0.6, zorder=3)
        ha = "left" if px >= 0 else "right"
        ax.text(px * 1.13, py * 1.13, r["neighbor_gene"], fontsize=7.3, color=DGRAY, ha=ha, va="center")

    ax.set_xlim(-1.6, 1.6)
    ax.set_ylim(-1.6, 1.6)
    ax.set_aspect("equal")
    ax.axis("off")
    n_sig = int((pd.to_numeric(nb["neighbor_crispr_fdr"], errors="coerce") < 0.1).sum())
    _title(ax, f"USP34's real 1-hop STRING neighbors (n={n}, line width/opacity = confidence)")
    _footnote(fig, f"Red nodes = neighbor gene itself reaches Gate-1 CRISPR significance (FDR<0.1) in this screen "
                    f"({n_sig} of {n} do, by this frozen table). KDM1A/TLK2 have no frozen direct-neighbor table at "
                    f"all; VEZF1 and EML5 each have only 1 frozen neighbor row (too sparse for their own network "
                    f"figure).", y=-0.02)
    _save(fig, stub, vector=False)


def build_e1_upset_evidence_intersection(stub: Path) -> None:
    """Manual UpSet-style plot (matplotlib, no upsetplot dependency) for
    the 13 significant sensitising hits x 6 evidence-set booleans.
    'high_baseline_dependency' and 'low_baseline_dependency' are kept as
    two SEPARATE, non-overlapping sets deliberately -- collapsing them
    into one 'dependency' axis would treat an advantage and a liability
    as the same evidence type."""
    sets_df = pv2.build_evidence_sets_13()
    cols = ["resistance_model_rna_support", "recurrence_associated_rna_support", "high_baseline_dependency",
            "low_baseline_dependency", "experimental_structure", "validated_inhibitor"]
    col_labels = ["Resistance-model\nRNA support", "Recurrence-associated\nRNA support", "High baseline\ndependency",
                  "Low baseline\ndependency", "Experimental\nstructure", "Validated\ninhibitor"]

    combos = sets_df.groupby(cols)["gene"].apply(list).reset_index()
    combos["n_sets"] = combos[cols].sum(axis=1)
    combos["size"] = combos["gene"].apply(len)
    combos = combos[combos["n_sets"] > 0].sort_values(["n_sets", "size"], ascending=[False, False]).reset_index(drop=True)

    n_combos = len(combos)
    fig = plt.figure(figsize=(13, 8), dpi=300)
    gs = fig.add_gridspec(2, 1, height_ratios=(2.4, 1.6), hspace=0.05)
    axBar = fig.add_subplot(gs[0])
    axMat = fig.add_subplot(gs[1], sharex=axBar)

    for i, r in combos.iterrows():
        colors = [FOCUS_COLORS.get(g, GRAY) for g in r["gene"]]
        bottom = 0
        for g, c in zip(r["gene"], colors):
            axBar.bar([i], [1], bottom=bottom, color=c, width=0.65, edgecolor="white", linewidth=0.4)
            bottom += 1
    axBar.set_ylabel("genes in intersection", fontsize=9.5)
    axBar.set_xticks([])
    _clean_axes(axBar)
    _title(axBar, "Evidence-set intersections across the 13 significant sensitising hits (colored segments = genes; 4 focus genes in accent color)", fontsize=11.5)

    for i, r in combos.iterrows():
        active = [j for j, c in enumerate(cols) if r[c]]
        axMat.scatter([i] * len(cols), range(len(cols)), s=90, color=VLGRAY, zorder=2)
        axMat.scatter([i] * len(active), active, s=90, color=DGRAY, zorder=3)
        if len(active) > 1:
            axMat.plot([i, i], [min(active), max(active)], color=DGRAY, linewidth=1.6, zorder=1)
    axMat.set_yticks(range(len(cols)))
    axMat.set_yticklabels(col_labels, fontsize=8.7)
    axMat.set_xlim(-0.7, n_combos - 0.3)
    axMat.set_xlabel(f"{n_combos} distinct evidence combinations observed among the 13 hits", fontsize=9.5)
    axMat.set_xticks([])
    for spine in ("top", "right", "bottom"):
        axMat.spines[spine].set_visible(False)

    focus_gene_x = {}
    for i, r in combos.iterrows():
        for g in r["gene"]:
            if g in FOCUS_FOUR:
                focus_gene_x.setdefault(g, i)
    for g, i in focus_gene_x.items():
        axBar.annotate(g, xy=(i, 0.5), xytext=(0, 12), textcoords="offset points", fontsize=7.6,
                        fontweight="bold", color=FOCUS_COLORS[g], ha="center")

    _footnote(fig, "'High baseline dependency' (>=50% of 11 ER+/luminal lines strongly dependent) and 'low baseline "
                    "dependency' (<10%) are kept as separate, non-overlapping sets -- never merged into one "
                    "'dependency' axis, since one is context and the other may be a liability.", y=-0.02)
    _save(fig, stub)


def build_e2_quantitative_evidence_map(stub: Path) -> None:
    em = pv2.load_evidence_matrix_13()
    struct = pv2.load_structural_tractability_audit().set_index("gene")
    x = -np.log10(em["crispr_fdr"])
    y = em["frac_strongly_dependent_er_luminal"] * 100

    fig, ax = plt.subplots(figsize=(10, 7.5), dpi=300)
    for _, r in em.iterrows():
        gene = r["gene"]
        gx, gy = -np.log10(r["crispr_fdr"]), r["frac_strongly_dependent_er_luminal"] * 100
        has_rna = (pd.notna(r["gse118713_fdr"]) and r["gse118713_fdr"] < 0.05) or \
                  (pd.notna(r["gse111151_fdr"]) and r["gse111151_fdr"] < 0.05) or \
                  (pd.notna(r["gse240112_fdr"]) and r["gse240112_fdr"] < 0.05)
        marker = "^" if has_rna else "o"
        has_structure = gene in struct.index and bool(struct.loc[gene, "A_experimental_human_structure_exists"])
        edge_color = "#b0392f" if has_structure else "white"
        edge_width = 2.2 if has_structure else 0.8
        color = FOCUS_COLORS.get(gene, GRAY)
        size = 230 if gene in FOCUS_FOUR else 110
        ax.scatter([gx], [gy], s=size, marker=marker, color=color, edgecolor=edge_color, linewidth=edge_width, zorder=3)
        has_inhibitor = gene in struct.index and str(struct.loc[gene, "E_validated_selective_small_molecule_inhibitor"]).upper().startswith("YES")
        if has_inhibitor:
            ax.scatter([gx], [gy], s=size * 2.6, marker="o", facecolor="none", edgecolor="#2a7a2a", linewidth=1.3, zorder=2)
        label = gene if gene in FOCUS_FOUR else None
        if label:
            ax.annotate(label, xy=(gx, gy), xytext=(8, 6), textcoords="offset points", fontsize=10,
                        fontweight="bold", color=color)
    ax.set_xlabel(r"CRISPR sensitisation strength ($-\log_{10}$ FDR)", fontsize=11)
    ax.set_ylabel("Baseline ER+/luminal DepMap dependency\n(% of 11 dependency-evaluable lines)", fontsize=10.5)
    _clean_axes(ax)
    _title(ax, "Quantitative evidence map: 13 significant sensitising hits")

    handles = [
        Line2D([0], [0], marker="^", color="none", markerfacecolor=GRAY, markersize=9, label="has non-acute RNA support (FDR<0.05)"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor=GRAY, markersize=9, label="no non-acute RNA support"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor="none", markeredgecolor="#b0392f", markeredgewidth=2, markersize=9, label="experimental structure exists"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor="none", markeredgecolor="#2a7a2a", markersize=13, label="validated selective inhibitor"),
    ]
    ax.legend(handles=handles, loc="upper left", fontsize=8, frameon=False)
    _footnote(fig, "All 13 significant sensitising hits; 4 focus genes labeled and enlarged. No composite score -- "
                    "each visual channel (position, marker shape, ring color, halo) encodes one distinct, real "
                    "evidence dimension.", y=-0.04)
    _save(fig, stub)


def build_e3_historical_vs_postaudit(stub: Path) -> None:
    rule_sens = pv2.load_selection_rule_sensitivity()
    rule0 = rule_sens[rule_sens["rule"] == "RULE_0_original_frozen_gate"].set_index("gene")
    rule1 = rule_sens[rule_sens["rule"] == "RULE_1_crispr_only_no_rna_gate"].set_index("gene")
    genes = sorted(set(rule1.index) & set(pv2.load_significant_sensitising_hits()["gene"]))
    genes = sorted(genes, key=lambda g: rule1.loc[g, "rank"])

    fig, ax = plt.subplots(figsize=(9, 8.5), dpi=300)
    y = np.arange(len(genes))
    for i, gene in enumerate(genes):
        r1_rank = rule1.loc[gene, "rank"]
        r0_eligible = bool(rule0.loc[gene, "eligible"]) if gene in rule0.index else False
        color = FOCUS_COLORS.get(gene, GRAY)
        if r0_eligible:
            r0_rank = rule0.loc[gene, "rank"]
            ax.plot([0, 1], [i, i], color=LGRAY, linewidth=1, zorder=1)
            ax.scatter([1], [i], s=90, color=color, edgecolor=DGRAY, linewidth=1, zorder=3)
            ax.text(1.08, i, f"eligible, rank {int(r0_rank)}/4", fontsize=8, color=DGRAY, va="center")
        else:
            ax.text(1.08, i, "excluded by RNA-eligibility gate", fontsize=8, color=GRAY, va="center", style="italic")
        ax.scatter([0], [i], s=90, color=color, edgecolor=DGRAY, linewidth=1, zorder=3)
        weight = "bold" if gene in FOCUS_FOUR else "normal"
        ax.text(-0.08, i, f"{gene}  (CRISPR-only rank {int(r1_rank)}/13)", fontsize=9, color=color,
                ha="right", va="center", fontweight=weight)
    ax.set_xlim(-1.7, 2.6)
    ax.set_ylim(-1, len(genes))
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["Rule 1: CRISPR-only\n(no RNA gate)", "Rule 0: original frozen gate\n(RNA eligibility required)"], fontsize=9.5)
    ax.invert_yaxis()
    ax.tick_params(left=False, labelleft=False)
    for spine in ("top", "right", "left"):
        ax.spines[spine].set_visible(False)
    _title(ax, "Historical selection gate vs. CRISPR-only reordering")
    _footnote(fig, "Post-audit interpretation no longer treats absence of RNA change as disqualifying -- Rule 1 "
                    "shows the same 13 genes' CRISPR-only rank; Rule 0 shows which of those genes the ORIGINAL "
                    "frozen gate accepted (only the 4 original candidates could ever pass it by construction).", y=-0.03)
    _save(fig, stub)


def build_f1_depmap_heatmap(stub: Path) -> None:
    eff = pv2.load_depmap_effect_focus_four()
    names = pv2.load_depmap_model_names().set_index("ModelID")["CellLineName"]
    pivot = eff.pivot(index="cell_line", columns="gene", values="chronos_effect")[FOCUS_FOUR]
    pivot.index = [names.get(i, i) for i in pivot.index]
    pivot = pivot.sort_values("TLK2")

    fig, ax = plt.subplots(figsize=(6.5, 8), dpi=300)
    vmax = np.nanmax(np.abs(pivot.to_numpy()))
    im = ax.imshow(pivot.to_numpy(), cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
    for i in range(pivot.shape[0]):
        for j in range(pivot.shape[1]):
            v = pivot.iloc[i, j]
            ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=8,
                    color="white" if abs(v) > vmax * 0.55 else DGRAY)
    ax.set_xticks(range(4))
    ax.set_xticklabels(FOCUS_FOUR, fontsize=10.5)
    for tick, gene in zip(ax.get_xticklabels(), FOCUS_FOUR):
        tick.set_color(FOCUS_COLORS[gene])
        tick.set_fontweight("bold")
    ax.set_yticks(range(len(pivot)))
    ax.set_yticklabels(pivot.index, fontsize=9)
    cb = fig.colorbar(im, ax=ax, fraction=0.05, pad=0.03)
    cb.set_label("Chronos gene effect", fontsize=9)
    _title(ax, f"DepMap 26Q1: real per-line values, {len(pivot)} ER+/luminal models")
    _footnote(fig, "Rows sorted by TLK2 effect. More negative = stronger dependency. Real cell-line names, not "
                    "anonymized IDs.", y=-0.04)
    _save(fig, stub)


def build_f2_cellline_fingerprint(stub: Path) -> None:
    eff = pv2.load_depmap_effect_focus_four()
    names = pv2.load_depmap_model_names().set_index("ModelID")["CellLineName"]
    pivot = eff.pivot(index="cell_line", columns="gene", values="chronos_effect")[FOCUS_FOUR]

    fig, ax = plt.subplots(figsize=(9.5, 7.5), dpi=300)
    x = np.arange(len(FOCUS_FOUR))
    cmap = plt.get_cmap("tab20")
    end_labels = []
    for i, (cell_line, row) in enumerate(pivot.iterrows()):
        name = names.get(cell_line, cell_line)
        color = cmap(i % 20)
        ax.plot(x, row.values, color=color, linewidth=1.4, marker="o", markersize=5, alpha=0.85, zorder=2)
        real_y = float(row.values[-1])
        end_labels.append([name, color, real_y, real_y])  # [name, color, real_y, label_y]

    # greedy vertical de-collision for the 11 end-of-line labels (several
    # VEZF1 endpoints cluster within ~0.05 of each other). The REAL
    # endpoint (real_y) is kept separate from the dodged label position
    # (label_y), with a thin leader line drawn between them whenever they
    # diverge, so a label can never appear to belong to the wrong cell
    # line (a real bug in an earlier version of this figure, caught by
    # adversarial review 2026-08-16).
    end_labels.sort(key=lambda e: e[3])
    min_gap = 0.055
    for j in range(1, len(end_labels)):
        if end_labels[j][3] - end_labels[j - 1][3] < min_gap:
            end_labels[j][3] = end_labels[j - 1][3] + min_gap
    for name, color, real_y, label_y in end_labels:
        if abs(label_y - real_y) > 1e-9:
            ax.plot([3, 3.15], [real_y, label_y], color=color, linewidth=0.6, alpha=0.6, zorder=2)
        ax.annotate(name, xy=(3.15, label_y), xytext=(6, 0), textcoords="offset points", fontsize=7, color=color, va="center")
    ax.axhline(-1, color=GRAY, linewidth=0.8, linestyle="--", zorder=0)
    ax.text(0.02, -1, "strong-dependency reference", fontsize=7.5, color=GRAY, va="bottom", transform=ax.get_yaxis_transform())
    ax.set_xticks(x)
    ax.set_xticklabels(FOCUS_FOUR, fontsize=11)
    for tick, gene in zip(ax.get_xticklabels(), FOCUS_FOUR):
        tick.set_color(FOCUS_COLORS[gene])
        tick.set_fontweight("bold")
    ax.set_ylabel("Chronos gene effect", fontsize=10.5)
    ax.set_xlim(-0.3, 4.3)
    _clean_axes(ax)
    _title(ax, "Per-cell-line fingerprint across the 4 focus genes (11 real ER+/luminal lines)")
    _footnote(fig, "Each line is one real DepMap cell-line model, connected across its own 4 gene-effect values -- "
                    "shows heterogeneity within each gene's distribution, not just an aggregate percentage.", y=-0.03)
    _save(fig, stub)


def build_f3_human_depmap_combo(stub: Path) -> None:
    g240112 = pv2.load_gse240112_focus_gene_tumours()
    g118713 = pv2.load_gse118713_focus_gene_samples()
    eff = pv2.load_depmap_effect_focus_four()

    fig = plt.figure(figsize=(13, 6), dpi=300)
    gs = fig.add_gridspec(1, 3, width_ratios=(1, 1, 1.6), wspace=0.35)

    axA = fig.add_subplot(gs[0])
    sub = g240112[g240112["gene"] == "VEZF1"]
    for i, grp in enumerate(["PT", "RT"]):
        vals = sub.loc[sub["group"] == grp, "log2cpm"].to_numpy()
        mean = float(np.mean(vals))
        axA.plot([i - 0.18, i + 0.18], [mean, mean], color=FOCUS_COLORS["VEZF1"], linewidth=2.2, zorder=3)
        jitter = rng.uniform(-0.1, 0.1, size=len(vals))
        axA.scatter(np.full(len(vals), i) + jitter, vals, s=75, color=FOCUS_COLORS["VEZF1"], edgecolor="white", linewidth=0.7, zorder=4)
    axA.set_xticks([0, 1])
    axA.set_xticklabels(["primary", "recurrent"], fontsize=9.5)
    axA.set_ylabel("log2(CPM)", fontsize=10)
    _clean_axes(axA)
    _title(axA, "VEZF1: GSE240112\n(recurrence, unpaired)", fontsize=10.5)

    axB = fig.add_subplot(gs[1])
    sub = g118713[g118713["gene_symbol"] == "USP34"]
    for i, cond in enumerate(["MCF7", "TAMR", "FASR"]):
        vals = sub.loc[sub["condition"] == cond, "tpm"].to_numpy()
        mean = float(np.mean(vals))
        axB.plot([i - 0.18, i + 0.18], [mean, mean], color=FOCUS_COLORS["USP34"], linewidth=2.2, zorder=3)
        jitter = rng.uniform(-0.1, 0.1, size=len(vals))
        axB.scatter(np.full(len(vals), i) + jitter, vals, s=60, color=FOCUS_COLORS["USP34"], edgecolor="white", linewidth=0.6, zorder=4)
    axB.set_xticks(range(3))
    axB.set_xticklabels(["MCF7", "TAMR", "FASR"], fontsize=9)
    axB.set_ylabel("TPM", fontsize=10)
    _clean_axes(axB)
    _title(axB, "USP34: GSE118713\n(resistance model)", fontsize=10.5)

    axC = fig.add_subplot(gs[2])
    names = pv2.load_depmap_model_names().set_index("ModelID")["CellLineName"]
    pivot = eff.pivot(index="cell_line", columns="gene", values="chronos_effect")[FOCUS_FOUR]
    pivot.index = [names.get(i, i) for i in pivot.index]
    pivot = pivot.sort_values("TLK2")
    vmax = np.nanmax(np.abs(pivot.to_numpy()))
    im = axC.imshow(pivot.to_numpy(), cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
    axC.set_xticks(range(4))
    axC.set_xticklabels(FOCUS_FOUR, fontsize=9)
    for tick, gene in zip(axC.get_xticklabels(), FOCUS_FOUR):
        tick.set_color(FOCUS_COLORS[gene])
        tick.set_fontweight("bold")
    axC.set_yticks(range(len(pivot)))
    axC.set_yticklabels(pivot.index, fontsize=7.3)
    cb = fig.colorbar(im, ax=axC, fraction=0.05, pad=0.03)
    cb.set_label("Chronos", fontsize=8)
    _title(axC, "Baseline DepMap dependency\n(11 ER+/luminal lines)", fontsize=10.5)

    fig.suptitle("Human tumour / resistance-model context vs. baseline cancer-cell dependency", fontsize=13, x=0.02, ha="left", y=1.03)
    _footnote(fig, "Panels A-B: real human/cell-line evidence supporting VEZF1 and USP34 respectively. Panel C: "
                    "baseline dependency for all 4 focus genes side by side -- a different, orthogonal quantity, "
                    "not more of the same evidence.", y=-0.04)
    _save(fig, stub)


def build_f4_tcga_secondary(stub: Path) -> None:
    """TCGA kept deliberately small/secondary -- both USP34 and VEZF1
    associations are weak/non-significant (FDR 0.21 and 0.90), so this is
    NOT built as a hero figure, per the task's explicit instruction."""
    forest = pv2.load_tcga_forest_original_four()
    sub = forest[(forest["candidate"].isin(["USP34", "VEZF1"])) & (forest["comparison"] == "tumor_vs_normal_PAIRED")]

    fig, ax = plt.subplots(figsize=(6, 3), dpi=300)
    y = np.arange(len(sub))
    for i, (_, r) in enumerate(sub.iterrows()):
        color = FOCUS_COLORS[r["candidate"]]
        ax.plot([r["ci_low"], r["ci_high"]], [i, i], color=color, linewidth=2, zorder=2)
        ax.scatter([r["mean_diff"]], [i], s=70, color=color, edgecolor="white", linewidth=0.8, zorder=3)
        ax.text(r["ci_high"] + 0.02, i, f"FDR={r['fdr']:.2f}", fontsize=8, color=DGRAY, va="center")
    ax.axvline(0, color=LGRAY, linewidth=0.8, zorder=0)
    ax.set_yticks(y)
    ax.set_yticklabels(sub["candidate"], fontsize=10)
    for tick, gene in zip(ax.get_yticklabels(), sub["candidate"]):
        tick.set_color(FOCUS_COLORS[gene])
    ax.set_xlabel("TCGA-BRCA paired tumour vs normal, log2FC (95% CI)", fontsize=8.7)
    _clean_axes(ax)
    _title(ax, "TCGA (secondary -- both weak/non-significant)", fontsize=10.5)
    _footnote(fig, "Deliberately small: TCGA was only assessed for the original 4 candidates, and neither USP34 "
                    "nor VEZF1 reaches FDR<0.05 here -- not a hero-figure result.", y=-0.14)
    _save(fig, stub)


# ===========================================================================
# SECTION G -- structural / pharmacological comparison
# ===========================================================================

def _render_kdm1a_tlk2_panels(out_dir: Path) -> dict[str, Path]:
    paths = pv2.kdm1a_tlk2_structure_paths()
    script = Path(__file__).resolve().parent / "poster_exploration_v2_structure.py"
    pymol_bin = shutil.which("pymol")
    if pymol_bin is None:
        raise RuntimeError("pymol not found on PATH -- install pymol-open-source (see environment.yml)")
    cmd = [pymol_bin, "-cq", str(script), "--", str(paths["6NQU"]), str(paths["5O0Y"]), str(out_dir)]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"pymol render failed:\n{result.stdout}\n{result.stderr}")
    names = ["kdm1a_full", "kdm1a_pocket", "tlk2_full", "tlk2_pocket"]
    panels = {n: out_dir / f"{n}.png" for n in names}
    for p in panels.values():
        assert p.exists(), f"expected pymol output missing: {p}"
    return panels


def build_g1_three_structures_plus_vezf1(stub: Path) -> None:
    with tempfile.TemporaryDirectory() as tmp:
        kt_panels = _render_kdm1a_tlk2_panels(Path(tmp))
        kdm1a_img = _autocrop_white(plt.imread(kt_panels["kdm1a_full"]))
        tlk2_img = _autocrop_white(plt.imread(kt_panels["tlk2_full"]))
        usp34_dir = Path(tmp) / "usp34"
        usp34_dir.mkdir()
        usp34_panels = _render_structure_bank_panels(usp34_dir)
        usp34_img = _autocrop_white(plt.imread(usp34_panels["comp_bound"]))

    fig, axes = plt.subplots(1, 4, figsize=(15, 4.6), dpi=300)
    for ax, img, gene, note in zip(axes[:3], [kdm1a_img, tlk2_img, usp34_img], FOCUS_FOUR[:3],
                                    ["6NQU, inhibitor-bound\n(GSK2879552)", "5O0Y, ATP-analog-bound\n(AGS, not an inhibitor)",
                                     "7W3U, covalent probe-bound\n(ubiquitin-propargylamide)"]):
        ax.imshow(img)
        ax.axis("off")
        ax.text(0.5, 1.03, gene, transform=ax.transAxes, fontsize=13, fontweight="bold", ha="center",
                color=FOCUS_COLORS[gene])
        ax.text(0.5, -0.04, note, transform=ax.transAxes, fontsize=7.8, ha="center", va="top", color=DGRAY)

    axV = axes[3]
    axV.axis("off")
    axV.add_patch(Rectangle((0.05, 0.05), 0.9, 0.9, transform=axV.transAxes, facecolor=VLGRAY,
                             edgecolor=GRAY, linewidth=1, zorder=1))
    axV.text(0.5, 1.03, "VEZF1", transform=axV.transAxes, fontsize=13, fontweight="bold", ha="center",
             color=FOCUS_COLORS["VEZF1"])
    axV.text(0.5, 0.5, "No experimental\nstructure exists.\n\nOnly a homology model\n(built on the unrelated\nZif268 zinc finger,\nPDB 1AAY) has been\npublished.", transform=axV.transAxes,
             fontsize=9, ha="center", va="center", color=DGRAY)

    fig.suptitle("Structural comparison, matched rendering style", fontsize=14, x=0.02, ha="left", y=1.08)
    _footnote(fig, "KDM1A and TLK2 structures fetched this phase from RCSB (real, published, experimentally solved); "
                    "USP34 structures already frozen from the final_translational phase. Same background, cartoon "
                    "style, and ray-tracing settings across all three -- VEZF1 is shown as a genuine absence, not a "
                    "fabricated structure.", y=-0.05)
    _save(fig, stub, vector=False)


def build_g2_structure_pharmacology_maturity(stub: Path) -> None:
    struct = pv2.load_structural_tractability_audit().set_index("gene")
    with tempfile.TemporaryDirectory() as tmp:
        kt_panels = _render_kdm1a_tlk2_panels(Path(tmp))
        kdm1a_img = _autocrop_white(plt.imread(kt_panels["kdm1a_full"]))
        tlk2_img = _autocrop_white(plt.imread(kt_panels["tlk2_full"]))
        usp34_dir = Path(tmp) / "usp34"
        usp34_dir.mkdir()
        usp34_panels = _render_structure_bank_panels(usp34_dir)
        usp34_img = _autocrop_white(plt.imread(usp34_panels["comp_bound"]))

    labels = {
        "KDM1A": "Clinical-stage selective pharmacology\n(e.g. iadademstat)",
        "TLK2": "Structure exists; no validated\nselective clinical inhibitor",
        "USP34": "Activity-probe reactivity;\nno validated selective inhibitor",
        "VEZF1": "No experimental structure;\ndifficult TF pharmacology",
    }
    fig, axes = plt.subplots(1, 4, figsize=(15, 5), dpi=300, gridspec_kw=dict(hspace=0.05))
    for ax, img, gene in zip(axes[:3], [kdm1a_img, tlk2_img, usp34_img], FOCUS_FOUR[:3]):
        ax.imshow(img)
        ax.axis("off")
    axes[3].axis("off")
    axes[3].add_patch(Rectangle((0.05, 0.05), 0.9, 0.9, transform=axes[3].transAxes, facecolor=VLGRAY,
                                 edgecolor=GRAY, linewidth=1))
    for ax, gene in zip(axes, FOCUS_FOUR):
        ax.text(0.5, 1.05, gene, transform=ax.transAxes, fontsize=12.5, fontweight="bold", ha="center", color=FOCUS_COLORS[gene])
        ax.text(0.5, -0.08, labels[gene], transform=ax.transAxes, fontsize=8.3, ha="center", va="top", color=DGRAY)
    fig.suptitle("Structure and pharmacological maturity, side by side", fontsize=14, x=0.02, ha="left", y=1.1)
    _save(fig, stub, vector=False)


def build_g3_pocket_closeups(stub: Path) -> None:
    with tempfile.TemporaryDirectory() as tmp:
        kt_panels = _render_kdm1a_tlk2_panels(Path(tmp))
        kdm1a_img = _autocrop_white(plt.imread(kt_panels["kdm1a_pocket"]))
        tlk2_img = _autocrop_white(plt.imread(kt_panels["tlk2_pocket"]))
        usp34_dir = Path(tmp) / "usp34"
        usp34_dir.mkdir()
        usp34_panels = _render_structure_bank_panels(usp34_dir)
        usp34_img = _autocrop_white(plt.imread(usp34_panels["comp_closeup"]))

    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.8), dpi=300)
    titles = [("KDM1A", "Inhibitor pocket (KWM)"), ("TLK2", "ATP pocket (AGS)"),
              ("USP34", "Catalytic cleft, covalent probe bound\n(Cys1903, ubiquitin-propargylamide)")]
    for ax, img, (gene, note) in zip(axes, [kdm1a_img, tlk2_img, usp34_img], titles):
        ax.imshow(img)
        ax.axis("off")
        ax.text(0.5, 1.05, gene, transform=ax.transAxes, fontsize=13, fontweight="bold", ha="center", color=FOCUS_COLORS[gene])
        ax.text(0.5, -0.04, note, transform=ax.transAxes, fontsize=9, ha="center", va="top", color=DGRAY)
    fig.suptitle("Pocket close-ups (VEZF1 omitted -- no experimental pocket exists)", fontsize=13, x=0.02, ha="left", y=1.08)
    _save(fig, stub, vector=False)


# ===========================================================================
# SECTION H -- translational / follow-up (data-grounded, no role cards)
# ===========================================================================

_H2_ARMS = {
    "KDM1A": ["control", "4-OHT", "+ iadademstat", "4-OHT + iadademstat"],
    "TLK2": ["control", "4-OHT", "TLK2 knockout", "4-OHT + TLK2 KO"],
    "USP34": ["control", "4-OHT", "USP34 knockout", "4-OHT + USP34 KO"],
    "VEZF1": ["control", "4-OHT", "VEZF1 knockdown", "4-OHT + VEZF1 KD"],
}


def _draw_cell(ax, x, y, r=0.16, face="white", edge=DGRAY, nucleus=None):
    ax.add_patch(plt.Circle((x, y), r, facecolor=face, edgecolor=edge, linewidth=1.1, zorder=3))
    if nucleus:
        ax.add_patch(plt.Circle((x, y), r * 0.45, facecolor=nucleus, edgecolor="none", zorder=4))


def build_h2_experimental_schematic(stub: Path) -> None:
    fig, axes = plt.subplots(4, 1, figsize=(9, 9), dpi=300)
    for ax, gene in zip(axes, FOCUS_FOUR):
        color = FOCUS_COLORS[gene]
        arms = _H2_ARMS[gene]
        for i, arm in enumerate(arms):
            cx = 0.6 + i * 1.5
            perturbed = "KO" in arm or "KD" in arm or "iadademstat" in arm
            treated = "4-OHT" in arm
            face = color if perturbed else "white"
            _draw_cell(ax, cx, 0.5, face=face, edge=DGRAY, nucleus="white" if perturbed else None)
            if treated:
                ax.add_patch(plt.Circle((cx, 0.5), 0.24, facecolor="none", edgecolor="#b0392f", linewidth=1.6,
                                         linestyle="--", zorder=2))
            ax.text(cx, 0.08, arm, ha="center", va="top", fontsize=7.6, color=DGRAY)
        ax.set_xlim(0, 0.6 + len(arms) * 1.5)
        ax.set_ylim(-0.15, 1.0)
        ax.axis("off")
        ax.text(-0.02, 0.5, gene, transform=ax.transAxes, fontsize=12, fontweight="bold", color=color, ha="right", va="center")

    handles = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor="white", markeredgecolor=DGRAY, markersize=9, label="unperturbed cell"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor=GRAY, markeredgecolor=DGRAY, markersize=9, label="genetic/pharmacological perturbation"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor="none", markeredgecolor="#b0392f", markersize=11, label="4-OHT / tamoxifen exposure"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=3, fontsize=8.5, frameon=False, bbox_to_anchor=(0.5, -0.02))
    fig.suptitle("Proposed follow-up experimental arms, per candidate", fontsize=13, x=0.02, ha="left", y=1.0)
    _footnote(fig, "KDM1A: existing clinical-stage inhibitor enables a pharmacological arm. TLK2/USP34/VEZF1: no "
                    "validated inhibitor exists, so genetic perturbation (knockout/knockdown) is the only current "
                    "option -- a real, disclosed asymmetry, not an oversight.", y=-0.1)
    _save(fig, stub)


def build_h3_translational_maturity_landscape(stub: Path) -> None:
    """A different axis pairing than E2/F-series: CRISPR sensitisation
    strength vs. an ORDINAL pharmacological-maturity ladder (not a
    composite score -- a single, named, ordered categorical axis read
    directly from Table 06b's own facets)."""
    em = pv2.load_evidence_matrix_13()
    struct = pv2.load_structural_tractability_audit().set_index("gene")

    def maturity(gene: str) -> int:
        if gene not in struct.index:
            return 0
        s = struct.loc[gene]
        if str(s["F_clinical_stage_pharmacology"]).upper().startswith("YES"):
            return 4
        if str(s["E_validated_selective_small_molecule_inhibitor"]).upper().startswith("YES"):
            return 3
        if str(s["C_ligand_or_probe_bound_structure"]).upper().startswith(("YES", "PARTIAL")):
            return 2
        if bool(s["A_experimental_human_structure_exists"]):
            return 1
        return 0

    fig, ax = plt.subplots(figsize=(9.5, 7), dpi=300)
    for _, r in em.iterrows():
        gene = r["gene"]
        gx = -np.log10(r["crispr_fdr"])
        gy = maturity(gene)
        jitter = rng.uniform(-0.08, 0.08)
        color = FOCUS_COLORS.get(gene, GRAY)
        size = 220 if gene in FOCUS_FOUR else 90
        ax.scatter([gx], [gy + jitter], s=size, color=color, edgecolor="white", linewidth=0.9, zorder=3)
        if gene in FOCUS_FOUR:
            ax.annotate(gene, xy=(gx, gy + jitter), xytext=(8, 0), textcoords="offset points", fontsize=10,
                        fontweight="bold", color=color, va="center")
    ax.set_yticks([0, 1, 2, 3, 4])
    ax.set_yticklabels(["no structure", "structure only", "probe/analog-\nbound structure",
                         "validated selective\ninhibitor", "clinical-stage\npharmacology"], fontsize=8.7)
    ax.set_xlabel(r"CRISPR sensitisation strength ($-\log_{10}$ FDR)", fontsize=10.5)
    ax.set_ylim(-0.5, 4.5)
    _clean_axes(ax)
    _title(ax, "Translational maturity ladder vs. CRISPR sensitisation strength")
    _footnote(fig, "Y-axis is a single ordered categorical read directly from Table 06b's own facets -- not a "
                    "composite or weighted score. All 13 significant sensitising hits shown; 4 focus genes labeled.", y=-0.03)
    _save(fig, stub)


def run_section_h(out_dir: Path = FIGURES) -> None:
    build_h2_experimental_schematic(out_dir / "H2_experimental_schematic")
    build_h3_translational_maturity_landscape(out_dir / "H3_translational_maturity_landscape")


def run_section_g(out_dir: Path = FIGURES) -> None:
    build_g1_three_structures_plus_vezf1(out_dir / "G1_three_structures_plus_vezf1")
    build_g2_structure_pharmacology_maturity(out_dir / "G2_structure_pharmacology_maturity")
    build_g3_pocket_closeups(out_dir / "G3_pocket_closeups")


def run_section_f(out_dir: Path = FIGURES) -> None:
    build_f1_depmap_heatmap(out_dir / "F1_depmap_heatmap")
    build_f2_cellline_fingerprint(out_dir / "F2_cellline_fingerprint")
    build_f3_human_depmap_combo(out_dir / "F3_human_depmap_combo")
    build_f4_tcga_secondary(out_dir / "F4_tcga_secondary")


def run_section_e(out_dir: Path = FIGURES) -> None:
    build_e1_upset_evidence_intersection(out_dir / "E1_upset_evidence_intersection")
    build_e2_quantitative_evidence_map(out_dir / "E2_quantitative_evidence_map")
    build_e3_historical_vs_postaudit(out_dir / "E3_historical_vs_postaudit")


def run_section_d(out_dir: Path = FIGURES) -> None:
    build_d1_candidate_program_network(out_dir / "D1_candidate_program_network")
    build_d2_usp34_local_neighborhood(out_dir / "D2_usp34_local_neighborhood")


def run_section_c(out_dir: Path = FIGURES) -> None:
    build_c1a_pathway_trajectories_lines(out_dir / "C1a_pathway_trajectories_lines")
    build_c1b_pathway_small_multiples(out_dir / "C1b_pathway_small_multiples")
    build_c2_estrogen_emt_hero(out_dir / "C2_estrogen_emt_hero")
    build_c3_enrichment_curves(out_dir / "C3_enrichment_curves")


def run_section_b(out_dir: Path = FIGURES) -> None:
    build_b1_gse118713_sample_panel(out_dir / "B1_gse118713_sample_panel")
    build_b2_gse111151_trajectories(out_dir / "B2_gse111151_trajectories")
    build_b3_gse240112_tumours(out_dir / "B3_gse240112_tumours")
    build_b4_gse245601_acute_paired(out_dir / "B4_gse245601_acute_paired")
    build_b5a_integrated_dataset_centric(out_dir / "B5a_integrated_dataset_centric")
    build_b5b_integrated_gene_centric(out_dir / "B5b_integrated_gene_centric")


# ===========================================================================
# Contact sheets
# ===========================================================================

def _contact_sheet(items: list[tuple[str, str]], out_png: Path, ncols: int, figures_dir: Path = FIGURES,
                    title: str = "") -> None:
    n = len(items)
    nrows = -(-n // ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(6.2 * ncols, 5.2 * nrows), dpi=160)
    axes = np.atleast_2d(axes)
    for i, (stem, label) in enumerate(items):
        r, c = divmod(i, ncols)
        ax = axes[r, c]
        img = plt.imread(figures_dir / f"{stem}.png")
        ax.imshow(img)
        ax.axis("off")
        ax.set_title(label, fontsize=11, fontweight="bold", color=DGRAY, pad=6)
    for i in range(n, nrows * ncols):
        r, c = divmod(i, ncols)
        axes[r, c].axis("off")
    if title:
        fig.suptitle(title, fontsize=15, fontweight="bold", y=1.0)
    fig.tight_layout(rect=[0, 0, 1, 0.97 if title else 1])
    fig.savefig(out_png, facecolor="white", bbox_inches="tight", dpi=160)
    plt.close(fig)
    logger.info("wrote %s", out_png)


CONTACT_A_ITEMS = [
    ("A1_genomewide_ranked_landscape", "A1. Genome-wide ranked landscape"),
    ("A2_lollipop_13_hits", "A2. 13-hit lollipop"),
    ("A3_redesigned_volcano", "A3. Redesigned volcano"),
    ("A4_hexbin_density", "A4. Hexbin density"),
]
CONTACT_B_ITEMS = [
    ("B1_gse118713_sample_panel", "B1. GSE118713 samples"),
    ("B2_gse111151_trajectories", "B2. GSE111151 trajectories"),
    ("B3_gse240112_tumours", "B3. GSE240112 tumours"),
    ("B4_gse245601_acute_paired", "B4. GSE245601 acute paired"),
    ("B5a_integrated_dataset_centric", "B5a. Dataset-centric"),
    ("B5b_integrated_gene_centric", "B5b. Gene-centric"),
]
CONTACT_C_ITEMS = [
    ("C1a_pathway_trajectories_lines", "C1a. Pathway trajectories (lines)"),
    ("C1b_pathway_small_multiples", "C1b. Pathway small multiples"),
    ("C2_estrogen_emt_hero", "C2. Estrogen/EMT hero"),
    ("C3_enrichment_curves", "C3. Enrichment curves"),
    ("D1_candidate_program_network", "D1. Candidate-program network"),
    ("D2_usp34_local_neighborhood", "D2. USP34 local neighborhood"),
]
CONTACT_D_ITEMS = [
    ("E1_upset_evidence_intersection", "E1. UpSet evidence intersection"),
    ("E2_quantitative_evidence_map", "E2. Quantitative evidence map"),
    ("E3_historical_vs_postaudit", "E3. Historical vs. post-audit"),
    ("F1_depmap_heatmap", "F1. DepMap heatmap"),
    ("F2_cellline_fingerprint", "F2. Cell-line fingerprint"),
    ("F3_human_depmap_combo", "F3. Human + DepMap combo"),
    ("F4_tcga_secondary", "F4. TCGA (secondary)"),
]
CONTACT_E_ITEMS = [
    ("G1_three_structures_plus_vezf1", "G1. 3 structures + VEZF1"),
    ("G2_structure_pharmacology_maturity", "G2. Structure + pharmacology"),
    ("G3_pocket_closeups", "G3. Pocket close-ups"),
    ("H2_experimental_schematic", "H2. Experimental schematic"),
    ("H3_translational_maturity_landscape", "H3. Translational maturity"),
]
SHORTLIST_ITEMS = [
    ("A1_genomewide_ranked_landscape", "A1. Genome-wide landscape"),
    ("A2_lollipop_13_hits", "A2. 13-hit lollipop"),
    ("B1_gse118713_sample_panel", "B1. GSE118713 samples"),
    ("B2_gse111151_trajectories", "B2. GSE111151 trajectories"),
    ("B3_gse240112_tumours", "B3. GSE240112 tumours"),
    ("B5a_integrated_dataset_centric", "B5a. Dataset-centric overview"),
    ("C2_estrogen_emt_hero", "C2. Estrogen/EMT hero"),
    ("E1_upset_evidence_intersection", "E1. UpSet evidence intersection"),
    ("E2_quantitative_evidence_map", "E2. Quantitative evidence map"),
    ("E3_historical_vs_postaudit", "E3. Historical vs. post-audit"),
    ("F1_depmap_heatmap", "F1. DepMap heatmap"),
    ("G1_three_structures_plus_vezf1", "G1. 3 structures + VEZF1"),
]


def build_contact_sheets(out_dir: Path = FIGURES) -> None:
    _contact_sheet(CONTACT_A_ITEMS, out_dir / "CONTACT_A_CRISPR.png", ncols=2, title="Section A -- CRISPR discovery")
    _contact_sheet(CONTACT_B_ITEMS, out_dir / "CONTACT_B_TRANSCRIPTOMICS.png", ncols=3, title="Section B -- Transcriptomics")
    _contact_sheet(CONTACT_C_ITEMS, out_dir / "CONTACT_C_PATHWAYS_NETWORK.png", ncols=3, title="Sections C+D -- Pathways & network")
    _contact_sheet(CONTACT_D_ITEMS, out_dir / "CONTACT_D_EVIDENCE_HUMAN_DEPMAP.png", ncols=3, title="Sections E+F -- Evidence intersection & human/DepMap")
    _contact_sheet(CONTACT_E_ITEMS, out_dir / "CONTACT_E_STRUCTURE_TRANSLATIONAL.png", ncols=3, title="Sections G+H -- Structure & translational")
    _contact_sheet(SHORTLIST_ITEMS, out_dir / "CONTACT_ALL_SHORTLIST.png", ncols=3,
                   title="Shortlist -- 12 strongest exploratory candidates")


def run(out_dir: str | Path = FIGURES) -> None:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    run_section_a(out_dir)
    run_section_b(out_dir)
    run_section_c(out_dir)
    run_section_d(out_dir)
    run_section_e(out_dir)
    run_section_f(out_dir)
    run_section_g(out_dir)
    run_section_h(out_dir)
    build_contact_sheets(out_dir)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run()
