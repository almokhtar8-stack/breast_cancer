"""Poster FIGURE BANK -- candidate scientific figures (visualization only).

Builds ~12 publication-style candidate figures into
results/figures/poster_candidates/, one per genuinely distinct scientific
question this repository has already answered. No poster layout is decided
here; see results/reports/poster/FIGURE_BANK_REVIEW.md for the CORE /
SUPPORTING / SUPPLEMENTARY judgment and results/figures/poster_candidates/
POSTER_FIGURE_CONTACT_SHEET.png for a visual index.

Design system (Step 15/17 of the request this module implements):
  - no giant titles inside any panel -- small panel letters + a short
    descriptive label only; a figure's headline claim is a poster section
    heading, not text baked into the PNG
  - candidate identity colors are fixed and consistent across every figure
    in this bank (GENE_COLORS, imported from poster_figures_visualization,
    the exact palette already used in the frozen poster/ set)
  - background/non-candidate data stays neutral grey; candidate color is
    never used for anything that isn't candidate identity
  - every figure exports .png (300dpi), .svg, and .pdf

No new analysis, no recomputation, no re-ranking. See
src/poster_figures_bank_data.py for the exact source table for every
plotted number.
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
from matplotlib.colors import TwoSlopeNorm
from matplotlib.lines import Line2D
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch

from src import poster_figures_bank_data as pbd
from src import poster_figures_data as pfd
from src.poster_figures_visualization import (
    BAD,
    DGRAY,
    GENE_COLORS,
    GENE_ORDER,
    GOOD,
    GRAY,
    LGRAY,
    _autocrop_white,
    _strip_box,
)

logger = logging.getLogger(__name__)

FIGURES = Path("results/figures/poster_candidates")


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


def _panel_letter(ax, letter: str, label: str, fontsize=12) -> None:
    ax.set_title(f"{letter}. {label}", fontsize=fontsize, fontweight="bold", loc="left", pad=8)


# ---------------------------------------------------------------------------
# 01a -- genome-wide CRISPR ranked-effect plot
# ---------------------------------------------------------------------------

def build_01a_crispr_ranked_effect(stub: Path) -> None:
    df = pfd.load_genomewide_crispr().sort_values("effect_size").reset_index(drop=True)
    df["rank"] = np.arange(1, len(df) + 1)
    gate1 = df["fdr"] < 0.1

    fig, ax = plt.subplots(figsize=(11, 7), dpi=300)
    ax.scatter(df.loc[~gate1, "rank"], df.loc[~gate1, "effect_size"], s=4, color=GRAY, alpha=0.22,
               linewidth=0, rasterized=True, zorder=1)
    ax.scatter(df.loc[gate1, "rank"], df.loc[gate1, "effect_size"], s=18, color=DGRAY, alpha=0.6,
               linewidth=0, zorder=2, label=f"Gate-1 hits (FDR<0.1, n={int(gate1.sum())})")
    ax.axhline(0, color=GRAY, linewidth=0.8, zorder=0)

    ordered = sorted(GENE_ORDER, key=lambda g: df.loc[df["gene"] == g, "rank"].iloc[0])
    label_x0 = df["rank"].min() + 60
    label_ys = np.linspace(1.35, -0.35, len(ordered))
    for gene, ly in zip(ordered, label_ys):
        row = df[df["gene"] == gene].iloc[0]
        ax.scatter([row["rank"]], [row["effect_size"]], s=200, color=GENE_COLORS[gene],
                   edgecolor=DGRAY, linewidth=1.5, zorder=5)
        ax.plot([label_x0, row["rank"]], [ly, row["effect_size"]], color=GENE_COLORS[gene], linewidth=1.0, zorder=4)
        ax.text(label_x0 + 30, ly, f"{gene}  (rank {int(row['rank'])}/{len(df)}, effect {row['effect_size']:.2f})",
                ha="left", va="center", fontsize=10, fontweight="bold", color=GENE_COLORS[gene])

    ax.set_xlabel(f"Gene rank, sorted by CRISPR effect size (n={len(df):,} fitted genes)\n" +
                  r"$\bf{negative}$ = sensitising knockout under 4-OHT   |   positive = tolerance-associated knockout", fontsize=11)
    ax.set_ylabel("CRISPR effect size (Hany et al. screen)", fontsize=12)
    _clean_axes(ax)
    ax.legend(loc="lower right", fontsize=10, frameon=False)
    _panel_letter(ax, "1a", "Genome-wide CRISPR screen, ranked by effect", fontsize=13)
    _save(fig, stub)


# ---------------------------------------------------------------------------
# 01c -- refined volcano (stripped down, USP34/VEZF1 dominant)
# ---------------------------------------------------------------------------

def build_01c_crispr_volcano_refined(stub: Path) -> None:
    df = pfd.load_genomewide_crispr()
    neglog10fdr = -np.log10(df["fdr"].clip(lower=1e-300))
    gate1 = df["fdr"] < 0.1

    fig, ax = plt.subplots(figsize=(9.5, 7.5), dpi=300)
    ax.scatter(df.loc[~gate1, "effect_size"], neglog10fdr[~gate1], s=3.5, color=GRAY, alpha=0.18,
               linewidth=0, rasterized=True, zorder=1)
    ax.scatter(df.loc[gate1, "effect_size"], neglog10fdr[gate1], s=10, color=DGRAY, alpha=0.5,
               linewidth=0, zorder=2)
    ax.axhline(-np.log10(0.1), color=GRAY, linewidth=0.9, linestyle="--", zorder=0)
    ax.axvline(0, color=GRAY, linewidth=0.7, zorder=0)

    primary = ["USP34", "VEZF1"]
    secondary_offsets = {"EML5": (10, -6), "CITED2": (-14, 14)}
    for gene, (dx, dy) in secondary_offsets.items():
        row = df[df["gene"] == gene].iloc[0]
        y = -np.log10(row["fdr"])
        ax.scatter([row["effect_size"]], [y], s=110, color=GENE_COLORS[gene], edgecolor=DGRAY,
                   linewidth=1.0, alpha=0.85, zorder=4)
        ha = "left" if dx > 0 else "right"
        ax.annotate(gene, (row["effect_size"], y), xytext=(dx, dy), textcoords="offset points",
                    fontsize=9.5, color=GENE_COLORS[gene], fontweight="bold", ha=ha)
    for gene in primary:
        row = df[df["gene"] == gene].iloc[0]
        y = -np.log10(row["fdr"])
        ax.scatter([row["effect_size"]], [y], s=260, color=GENE_COLORS[gene], edgecolor=DGRAY,
                   linewidth=1.8, zorder=5)
        dx = 8 if gene == "USP34" else -8
        ha = "left" if gene == "USP34" else "right"
        ax.annotate(gene, (row["effect_size"], y), xytext=(dx, 10), textcoords="offset points",
                    fontsize=13, color=GENE_COLORS[gene], fontweight="bold", ha=ha)

    ax.set_xlabel("CRISPR effect size  (negative = sensitising KO under 4-OHT)", fontsize=11.5)
    ax.set_ylabel(r"$-\log_{10}$(FDR)", fontsize=12.5)
    _clean_axes(ax)
    _panel_letter(ax, "1c", "CRISPR screen -- refined volcano", fontsize=13)
    _save(fig, stub)


# ---------------------------------------------------------------------------
# 03 -- GSE118713 resistance landscape (PCA + volcano + candidate inset)
# ---------------------------------------------------------------------------

GROUP_COLORS_118713 = {"MCF7": GRAY, "TAMR": "#0072B2", "FASR": "#CC79A7"}


def build_03_gse118713_landscape(stub: Path) -> None:
    pca = pbd.load_gse118713_pca()
    volc = pbd.load_gse118713_volcano()

    fig = plt.figure(figsize=(14, 6.6), dpi=300)
    gs = fig.add_gridspec(2, 3, width_ratios=(1.0, 1.35, 0.75), height_ratios=(1, 1), wspace=0.42, hspace=0.15)
    axA = fig.add_subplot(gs[:, 0])
    axB = fig.add_subplot(gs[:, 1])
    axC = fig.add_subplot(gs[:, 2])

    for cond, color in GROUP_COLORS_118713.items():
        sub = pca[pca["group"] == cond]
        axA.scatter(sub["PC1"], sub["PC2"], s=140, color=color, edgecolor=DGRAY, linewidth=1.1, zorder=3, label=cond)
    v1 = pca["pc1_variance_explained_pct"].iloc[0]
    v2 = pca["pc2_variance_explained_pct"].iloc[0]
    axA.set_xlabel(f"PC1 ({v1:.1f}% variance)", fontsize=11)
    axA.set_ylabel(f"PC2 ({v2:.1f}% variance)", fontsize=11)
    axA.axhline(0, color=LGRAY, linewidth=0.8, zorder=0)
    axA.axvline(0, color=LGRAY, linewidth=0.8, zorder=0)
    _clean_axes(axA)
    axA.legend(loc="best", fontsize=9.5, frameon=False)
    _panel_letter(axA, "A", "GSE118713 sample PCA (n=9)")

    gate1 = volc["is_gate1_hit"].fillna(False).astype(bool)
    neglog10fdr = -np.log10(volc["fdr"].clip(lower=1e-300))
    axB.scatter(volc.loc[~gate1, "log2fc"], neglog10fdr[~gate1], s=4, color=GRAY, alpha=0.2, linewidth=0,
                rasterized=True, zorder=1)
    axB.scatter(volc.loc[gate1, "log2fc"], neglog10fdr[gate1], s=14, color=DGRAY, alpha=0.55, linewidth=0, zorder=2)
    for gene in GENE_ORDER:
        hit = volc[volc["gene_symbol"] == gene]
        if len(hit) == 0:
            continue
        row = hit.iloc[0]
        y = -np.log10(max(row["fdr"], 1e-300))
        axB.scatter([row["log2fc"]], [y], s=170, color=GENE_COLORS[gene], edgecolor=DGRAY, linewidth=1.4, zorder=5)
        axB.annotate(gene, (row["log2fc"], y), xytext=(6, 6), textcoords="offset points",
                     fontsize=10, color=GENE_COLORS[gene], fontweight="bold")
    axB.axvline(0, color=GRAY, linewidth=0.7, zorder=0)
    axB.set_xlabel("log$_2$FC (TAMR vs MCF7)", fontsize=11)
    axB.set_ylabel(r"$-\log_{10}$(FDR)", fontsize=11)
    _clean_axes(axB)
    _panel_letter(axB, "B", f"TAMR vs MCF7 -- genome-wide (n={len(volc):,})")

    samples = pfd.load_gse118713_usp34_samples()
    order = ["MCF7", "TAMR", "FASR"]
    for i, cond in enumerate(order):
        vals = samples.loc[samples["condition"] == cond, "tpm"].to_numpy()
        _strip_box(axC, i, vals, GROUP_COLORS_118713[cond])
    axC.set_xticks(range(len(order)))
    axC.set_xticklabels(order, fontsize=10, fontweight="bold")
    axC.set_ylabel("USP34 TPM", fontsize=10.5)
    _clean_axes(axC)
    _panel_letter(axC, "C", "USP34 (inset)", fontsize=11)

    fig.text(0.01, -0.02, "GSE118713: 3 MCF7 (parental) + 3 TAMR + 3 FASR (acquired-resistance derivatives) bulk RNA-seq replicates.",
              fontsize=8.8, color=GRAY, style="italic")
    _save(fig, stub)


# ---------------------------------------------------------------------------
# 04 -- GSE240112 recurrence (genome-wide volcano + VEZF1 per-tumor)
# ---------------------------------------------------------------------------

CONDITION_COLORS_240112 = {"PT": GRAY, "RT": "#E69F00"}


def build_04_gse240112_recurrence(stub: Path) -> None:
    volc = pbd.load_gse240112_volcano()
    samples = pfd.load_gse240112_vezf1_samples()
    stats = pfd.load_gse240112_vezf1_stats()

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(12.5, 6), dpi=300, gridspec_kw={"width_ratios": (1.5, 1)})

    valid = volc["fdr"].notna()
    neglog10fdr = -np.log10(volc.loc[valid, "fdr"].clip(lower=1e-300))
    other = valid & ~volc["is_candidate"]
    axA.scatter(volc.loc[other, "log2fc"], -np.log10(volc.loc[other, "fdr"].clip(lower=1e-300)), s=4,
                color=GRAY, alpha=0.2, linewidth=0, rasterized=True, zorder=1)
    label_offsets = {"USP34": (6, 6), "VEZF1": (10, -14), "EML5": (10, -4), "CITED2": (10, 12)}
    for gene in GENE_ORDER:
        hit = volc[volc["gene"] == gene]
        if len(hit) == 0 or hit.iloc[0][["log2fc", "fdr"]].isna().any():
            continue
        row = hit.iloc[0]
        y = -np.log10(max(row["fdr"], 1e-300))
        axA.scatter([row["log2fc"]], [y], s=190, color=GENE_COLORS[gene], edgecolor=DGRAY, linewidth=1.5, zorder=5)
        dx, dy = label_offsets[gene]
        axA.annotate(gene, (row["log2fc"], y), xytext=(dx, dy), textcoords="offset points",
                     fontsize=11, color=GENE_COLORS[gene], fontweight="bold")
    axA.axvline(0, color=GRAY, linewidth=0.7, zorder=0)
    axA.set_xlabel("log$_2$FC (recurrent vs primary tumor)", fontsize=11.5)
    axA.set_ylabel(r"$-\log_{10}$(FDR)", fontsize=12)
    _clean_axes(axA)
    _panel_letter(axA, "A", f"GSE240112 -- genome-wide (n={valid.sum():,})")

    order = ["PT", "RT"]
    labels = ["Primary\n(n=3)", "Recurrent\n(n=3)"]
    for i, cond in enumerate(order):
        vals = samples.loc[samples["group"] == cond, "log2cpm"].to_numpy()
        _strip_box(axB, i, vals, CONDITION_COLORS_240112[cond], width=0.4)
    axB.set_xticks(range(len(order)))
    axB.set_xticklabels(labels, fontsize=11, fontweight="bold")
    axB.set_ylabel("VEZF1 log2(CPM)", fontsize=11.5)
    _clean_axes(axB)
    _panel_letter(axB, "B", f"VEZF1 per tumor  (log2FC {stats['log2fc']:.2f}, FDR {stats['genomewide_fdr']:.3f})", fontsize=11)

    fig.text(0.01, -0.03, "3 primary + 3 recurrent tumors, unpaired (different patients). Association with recurrence, not a causal resistance model.",
              fontsize=8.8, color=GRAY, style="italic")
    _save(fig, stub)


# ---------------------------------------------------------------------------
# 05 -- cross-dataset candidate effect map (improved forest)
# ---------------------------------------------------------------------------

def build_05_cross_dataset_effects(stub: Path) -> None:
    forest = pfd.build_forest_table()
    dataset_order = ["gse118713", "gse240112_tumor", "gse111151", "gse245601_epi"]
    n_ds = len(dataset_order)

    fig, ax = plt.subplots(figsize=(9.5, 8), dpi=300)
    row = 0
    yticks, yticklabels = [], []
    for gene in GENE_ORDER:
        for dkey in dataset_order:
            r = forest[(forest["gene"] == gene) & (forest["dataset_key"] == dkey)].iloc[0]
            sig = r["fdr"] < 0.05
            acute = dkey == "gse245601_epi"
            facecolor = GENE_COLORS[gene] if sig else "white"
            ax.scatter([r["log2fc"]], [row], s=140, facecolor=facecolor, edgecolor=GENE_COLORS[gene],
                       linewidth=2.0, zorder=3, marker="s" if acute else "o")
            yticks.append(row)
            label = pfd.DATASET_LABELS[dkey].split(" (")[0]
            yticklabels.append(label + (" [ACUTE]" if acute else ""))
            row -= 1
        row -= 0.7

    ax.axvline(0, color=GRAY, linewidth=1.0, zorder=0)
    ax.set_yticks(yticks)
    ax.set_yticklabels(yticklabels, fontsize=9)
    for lbl in ax.get_yticklabels():
        if "[ACUTE]" in lbl.get_text():
            lbl.set_color(BAD)
            lbl.set_fontstyle("italic")
    ax.set_xlabel("log$_2$ fold-change  (resistant/recurrent/tumor vs comparator; for the [ACUTE] row, 12h 4-OHT vs vehicle -- a different comparator)", fontsize=10.5)
    _clean_axes(ax)

    gene_row_centers, r = [], 0
    for gene in GENE_ORDER:
        gene_row_centers.append(r - (n_ds - 1) / 2)
        r -= n_ds + 0.7
    for gene, yc in zip(GENE_ORDER, gene_row_centers):
        ax.text(1.03, yc, gene, transform=ax.get_yaxis_transform(), fontsize=12, fontweight="bold",
                color=GENE_COLORS[gene], va="center", ha="left")

    handles = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor=DGRAY, markeredgecolor=DGRAY, markersize=9, label="FDR<0.05 (filled)"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor="white", markeredgecolor=DGRAY, markersize=9, label="not significant (open)"),
        Line2D([0], [0], marker="s", color="none", markerfacecolor=DGRAY, markeredgecolor=DGRAY, markersize=8, label="acute 12h context (square)"),
    ]
    ax.legend(handles=handles, loc="upper right", fontsize=8.8, frameon=False)
    _panel_letter(ax, "05", "Cross-dataset effect consistency, 4 candidates x 4 datasets", fontsize=13)
    _save(fig, stub)


# ---------------------------------------------------------------------------
# 06 -- resistance pathway landscape (pathway x dataset NES heatmap)
# ---------------------------------------------------------------------------

def build_06_pathway_landscape(stub: Path) -> None:
    pl = pbd.load_pathway_landscape()
    order = [lbl for _, _, lbl in pbd.PATHWAY_LANDSCAPE]
    ds_order = list(pbd.PATHWAY_RNA_DATASETS.keys())
    ds_labels = list(pbd.PATHWAY_RNA_DATASETS.values())

    mat = pl.pivot(index="pathway_label", columns="dataset_key", values="NES").loc[order, ds_order]
    fdr_mat = pl.pivot(index="pathway_label", columns="dataset_key", values="fdr").loc[order, ds_order]

    fig = plt.figure(figsize=(9.5, 8.6), dpi=300)
    gs = fig.add_gridspec(2, 1, height_ratios=(len(order), 1.15), hspace=0.55)
    ax = fig.add_subplot(gs[0])

    vmax = np.nanmax(np.abs(mat.to_numpy()))
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)
    im = ax.imshow(mat.to_numpy(), cmap="RdBu_r", norm=norm, aspect="auto")

    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            v = mat.iloc[i, j]
            if np.isnan(v):
                ax.text(j, i, "n/a", ha="center", va="center", fontsize=8, color=GRAY)
                continue
            sig = fdr_mat.iloc[i, j] < 0.05
            txt_color = "white" if abs(v) > vmax * 0.55 else DGRAY
            ax.text(j, i, f"{v:.1f}{'*' if sig else ''}", ha="center", va="center", fontsize=9,
                    color=txt_color, fontweight="bold" if sig else "normal")

    ax.set_xticks(range(len(ds_order)))
    ax.set_xticklabels(ds_labels, fontsize=10, rotation=20, ha="right")
    # the acute-12h column uses a different biological context than the 3
    # resistance datasets -- flagged in red/italic, consistent with figure 05
    acute_idx = ds_order.index("gse245601")
    ax.get_xticklabels()[acute_idx].set_color(BAD)
    ax.get_xticklabels()[acute_idx].set_fontstyle("italic")
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels(order, fontsize=10)
    ax.set_xticks(np.arange(-0.5, len(ds_order), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(order), 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=2)
    ax.tick_params(which="minor", length=0)
    # thick divider separating the 3 resistance datasets from the acute
    # 12h column -- same NES color scale is used for both (a genuine GSEA
    # NES is a genuine GSEA NES regardless of dataset), but the biological
    # context must not be visually conflated, hence the divider + red label
    ax.axvline(acute_idx - 0.5, color=DGRAY, linewidth=2.5, zorder=5)
    _panel_letter(ax, "06", "Resistance pathway landscape (real GSEA NES, RNA datasets)", fontsize=13)

    cbar = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.02)
    cbar.set_label("NES  (blue = down, red = up, vs each dataset's own comparator)", fontsize=9)
    cbar.ax.tick_params(labelsize=8)
    ax.text(0, -1.15, "* FDR<0.05   |   red column label = acute 12h context, not resistance", fontsize=8.5, color=GRAY, transform=ax.transData)

    axc = fig.add_subplot(gs[1])
    crispr = pbd.load_pathway_landscape_crispr().set_index("pathway_label").loc[order]
    cvals = crispr["NES"].to_numpy().reshape(1, -1)
    cvmax = np.nanmax(np.abs(cvals))
    cnorm = TwoSlopeNorm(vmin=-cvmax, vcenter=0, vmax=cvmax)
    axc.imshow(cvals, cmap="PuOr_r", norm=cnorm, aspect="auto")
    for j, v in enumerate(cvals[0]):
        if np.isnan(v):
            axc.text(j, 0, "n/a", ha="center", va="center", fontsize=8, color=GRAY)
        else:
            axc.text(j, 0, f"{v:.1f}", ha="center", va="center", fontsize=9, color=DGRAY)
    axc.set_xticks(range(len(order)))
    axc.set_xticklabels(order, fontsize=8.5, rotation=30, ha="right")
    axc.set_yticks([0])
    axc.set_yticklabels(["CRISPR\ncontext"], fontsize=9)
    axc.tick_params(length=0)
    for spine in axc.spines.values():
        spine.set_visible(False)
    axc.text(0.0, 1.45, "CRISPR pathway context uses a different metric (gene-dependency ranking, not expression fold-change) -- separate color scale, not directly comparable to the NES panel above.",
             transform=axc.transAxes, fontsize=8.5, color=GRAY, style="italic")

    _save(fig, stub)


# ---------------------------------------------------------------------------
# 07 -- candidate mechanism map (curated, small, real connectivity only)
# ---------------------------------------------------------------------------

def build_07_candidate_mechanism_map(stub: Path) -> None:
    # Every count/label in this figure is computed from the real
    # STRONG_CONSENSUS pathway-membership and pairwise-convergence tables --
    # nothing is hand-typed (a hardcoded "~40" was caught and fixed here
    # during review; the real CITED2 count is 49).
    memberships = pbd.load_candidate_strong_consensus_pathways()
    counts = memberships.groupby("candidate").size()
    convergence = pbd.load_candidate_convergence().set_index(["candidate_A", "candidate_B"])

    def conv(a, b):
        key = (a, b) if (a, b) in convergence.index else (b, a)
        return int(convergence.loc[key, "n_shared_pathways"])

    fig, ax = plt.subplots(figsize=(12, 7.5), dpi=300)
    ax.set_xlim(0, 12)
    ax.set_ylim(0, 8)
    ax.axis("off")

    candidate_pos = {
        "USP34": (1.6, 6.4),
        "VEZF1": (1.6, 1.6),
        "CITED2": (5.2, 4.0),
        "EML5": (10.4, 4.0),
    }
    label_wnt = f"Wnt / beta-catenin signalling\n({counts.get('USP34', 0)} STRONG_CONSENSUS GO terms)"
    label_heme = f"Heme metabolism +\nblood vessel morphogenesis\n({counts.get('VEZF1', 0)} STRONG_CONSENSUS terms)"
    label_cardio = f"Cardiovascular / embryonic\nmorphogenesis programs\n({counts.get('CITED2', 0)} STRONG_CONSENSUS terms)"
    label_uv = "UV_RESPONSE_DN /\nstress response"
    program_pos = {
        label_wnt: (5.2, 7.05),
        label_heme: (5.2, 0.95),
        label_cardio: (8.2, 6.2),
        label_uv: (8.2, 1.8),
    }
    program_of = {"USP34": label_wnt, "VEZF1": label_heme}

    for label, (x, y) in program_pos.items():
        ax.add_patch(FancyBboxPatch((x - 1.35, y - 0.55), 2.7, 1.1, boxstyle="round,pad=0.05,rounding_size=0.12",
                                     facecolor=LGRAY, edgecolor=GRAY, linewidth=1.0, zorder=2))
        ax.text(x, y, label, ha="center", va="center", fontsize=8.5, color=DGRAY, zorder=3, linespacing=1.25)

    for gene, (x, y) in candidate_pos.items():
        ax.add_patch(Circle((x, y), 0.62, facecolor=GENE_COLORS[gene], edgecolor=DGRAY, linewidth=1.6, zorder=4))
        ax.text(x, y, gene, ha="center", va="center", fontsize=11.5, fontweight="bold", color="white", zorder=5)

    for gene, prog_label in program_of.items():
        gx, gy = candidate_pos[gene]
        px, py = program_pos[prog_label]
        ax.add_patch(FancyArrowPatch((gx, gy), (px, py), arrowstyle="-", color=GENE_COLORS[gene], linewidth=2.0, zorder=1))

    # CITED2 -- broad, dual footprint (two programs); the umbrella program
    # names are interpretive groupings of the real GO-term list, labeled as
    # such rather than presented as a single GSEA entity.
    for prog_label in [label_cardio, label_uv]:
        px, py = program_pos[prog_label]
        cx, cy = candidate_pos["CITED2"]
        ax.add_patch(FancyArrowPatch((cx, cy), (px, py), arrowstyle="-", color=GENE_COLORS["CITED2"], linewidth=2.0, zorder=1))

    # real candidate-candidate convergence (verified counts, not inferred)
    n_usp34_cited2 = conv("USP34", "CITED2")
    n_vezf1_cited2 = conv("VEZF1", "CITED2")

    x1, y1 = candidate_pos["USP34"]
    x2, y2 = candidate_pos["CITED2"]
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-", color=GRAY, linewidth=1.0,
                                  linestyle="--", zorder=1, connectionstyle="arc3,rad=0.15"))
    ax.text((x1 + x2) / 2 - 0.3, (y1 + y2) / 2 + 0.55, f"{n_usp34_cited2} shared GO term(s)", fontsize=7.8, color=GRAY, style="italic")

    x1, y1 = candidate_pos["VEZF1"]
    x2, y2 = candidate_pos["CITED2"]
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-", color=GRAY, linewidth=1.6,
                                  linestyle="--", zorder=1, connectionstyle="arc3,rad=-0.15"))
    ax.text((x1 + x2) / 2 - 0.2, (y1 + y2) / 2 - 0.55, f"{n_vezf1_cited2} shared GO terms\n(vascular/morphogenesis)", fontsize=7.8, color=GRAY, style="italic", ha="center")

    ex, ey = candidate_pos["EML5"]
    n_eml5 = counts.get("EML5", 0)
    ax.text(ex, ey - 1.05, f"{n_eml5} STRONG_CONSENSUS\npathway memberships\n(0 of 5,847 pathways scored)", ha="center", va="top",
            fontsize=8.5, color=GENE_COLORS["EML5"], style="italic", linespacing=1.3)

    ax.text(0.0, 0.99, "07. Candidate -> pathway-program map", transform=ax.transAxes, fontsize=13,
            fontweight="bold", color=DGRAY, ha="left", va="top")
    ax.text(0.0, 0.02, "Edges = real STRONG_CONSENSUS GSEA membership / verified shared pathway terms only -- no generic network hairball; USP34-VEZF1, USP34-EML5, VEZF1-EML5, EML5-CITED2 have zero verified connectivity and are correctly shown as unconnected. Grey box labels are interpretive groupings of the underlying GO-term lists, not single GSEA entities.",
            transform=ax.transAxes, fontsize=8.3, color=GRAY, style="italic", wrap=True)

    _save(fig, stub)


# ---------------------------------------------------------------------------
# 08 -- TCGA human tumor validation
# ---------------------------------------------------------------------------

def build_08_tcga_validation(stub: Path) -> None:
    forest = pbd.load_tcga_expression_forest()
    clinical = pbd.load_tcga_clinical_er_adjusted().set_index("candidate")

    contrasts = [("tumor_vs_normal_PAIRED", "Tumor vs matched normal\n(paired, n=113)"),
                 ("ER+ vs ER- (clinical IHC)", "ER+ vs ER-\n(n=600 vs 178)")]

    genes_top_to_bottom = list(GENE_ORDER)
    genes_ascending = genes_top_to_bottom[::-1]  # tick position 0 = bottom

    fig, axes = plt.subplots(1, 2, figsize=(11, 6), dpi=300)
    for panel_letter, ax, (ckey, clabel) in zip(("A", "B"), axes, contrasts):
        sub = forest[forest["comparison"] == ckey].set_index("candidate")
        ys = np.arange(len(genes_ascending))
        for y, gene in zip(ys, genes_ascending):
            r = sub.loc[gene]
            sig = r["fdr"] < 0.05
            ax.plot([r["ci_low"], r["ci_high"]], [y, y], color=GENE_COLORS[gene], linewidth=2.2, zorder=2)
            ax.scatter([r["mean_diff"]], [y], s=130 if sig else 90, color=GENE_COLORS[gene] if sig else "white",
                       edgecolor=GENE_COLORS[gene], linewidth=2.0, zorder=3)
        ax.axvline(0, color=GRAY, linewidth=0.9, zorder=0)
        ax.set_yticks(ys)
        ax.set_yticklabels(genes_ascending, fontsize=11, fontweight="bold")
        for label, gene in zip(ax.get_yticklabels(), genes_ascending):
            label.set_color(GENE_COLORS[gene])
        ax.set_xlabel("mean difference (log2 TPM), 95% CI", fontsize=10)
        _clean_axes(ax)
        _panel_letter(ax, panel_letter, clabel, fontsize=10.5)

    ph_flags = (clinical["ph_assumption_p"] <= 0.05).sum()
    fig.suptitle("08. TCGA-BRCA human tumor expression validation", fontsize=13, fontweight="bold", x=0.02, ha="left", y=1.03)
    fig.text(0.02, -0.02, f"Cox survival HRs (ER+, age/stage-adjusted) exist for all 4 candidates but are not shown here -- "
                          f"{ph_flags}/4 violate the proportional-hazards assumption and are not presented as a clean clinical finding.",
             fontsize=8.5, color=GRAY, style="italic")
    _save(fig, stub)


# ---------------------------------------------------------------------------
# 09 -- Hany CRISPR effect vs DepMap baseline dependency (28 Gate-1 genes)
# ---------------------------------------------------------------------------

def build_09_hany_vs_depmap(stub: Path) -> None:
    df = pbd.load_depmap_gate1_dependency_summary()
    df = df[df["in_depmap"]].copy()

    fig, ax = plt.subplots(figsize=(10, 8), dpi=300)
    others = ~df["gene"].isin(GENE_ORDER)
    sizes_other = 30 + 260 * df.loc[others, "frac_strongly_dependent_er_luminal"].fillna(0)
    ax.scatter(df.loc[others, "crispr_effect_size"], df.loc[others, "median_chronos_er_luminal"],
               s=sizes_other, color=GRAY, alpha=0.45, edgecolor="white", linewidth=0.5, zorder=2)

    for gene in GENE_ORDER:
        hit = df[df["gene"] == gene]
        if len(hit) == 0:
            continue
        r = hit.iloc[0]
        size = 60 + 420 * (r["frac_strongly_dependent_er_luminal"] or 0)
        ax.scatter([r["crispr_effect_size"]], [r["median_chronos_er_luminal"]], s=size, color=GENE_COLORS[gene],
                   edgecolor=DGRAY, linewidth=1.8, zorder=5)
        ax.annotate(f"{gene}\n({r['frac_strongly_dependent_er_luminal']*100:.0f}% dependency-probability>0.5)",
                    (r["crispr_effect_size"], r["median_chronos_er_luminal"]), xytext=(10, 8),
                    textcoords="offset points", fontsize=9.5, color=GENE_COLORS[gene], fontweight="bold", linespacing=1.2)

    ax.axvline(0, color=GRAY, linewidth=0.8, zorder=0)
    ax.axhline(0, color=GRAY, linewidth=0.8, zorder=0)
    ax.set_xlabel("Hany CRISPR effect size (negative = sensitising KO under 4-OHT)", fontsize=11.5)
    ax.set_ylabel("Median DepMap 26Q1 Chronos gene effect\n(ER+/luminal breast lines, n=11) -- more negative = more dependent", fontsize=11)
    _clean_axes(ax)

    # y=low (more negative Chronos) = MORE dependent; y=high (near/above 0) =
    # LESS dependent. Labels placed accordingly (corrected after review --
    # an earlier draft had these two swapped).
    ax.text(0.02, 0.03, "sensitising + dependent\n(dual-action hypothesis)", transform=ax.transAxes, fontsize=9,
            color=GRAY, style="italic")
    ax.text(0.02, 0.97, "context-specific sensitiser\n(low baseline dependency)", transform=ax.transAxes, fontsize=9,
            color=GRAY, style="italic", va="top")

    n_other = int(others.sum())
    handles = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor=GRAY, alpha=0.5, markersize=6, label=f"other Gate-1 genes (n={n_other})"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor=GRAY, markeredgecolor=DGRAY, markersize=12, label="marker size = % ER+/luminal lines with dependency probability>0.5"),
    ]
    ax.legend(handles=handles, loc="lower right", fontsize=8.5, frameon=False)
    _panel_letter(ax, "09", "Drug-context sensitisation vs baseline cancer dependency (28 Gate-1 genes)", fontsize=12.5)
    fig.text(0.01, -0.03, "Quadrant labels are interpretive hypotheses, not proven biological classes.", fontsize=8.3, color=GRAY, style="italic")
    _save(fig, stub)


def build_09b_usp34_vezf1_line_dependencies(stub: Path) -> None:
    long = pfd.load_depmap_gene_effect_er_luminal()
    wide = long.pivot(index="cell_line", columns="gene", values="chronos_effect")[["USP34", "VEZF1"]]

    fig, ax = plt.subplots(figsize=(6.5, 7.5), dpi=300)
    for line, row in wide.iterrows():
        color = BAD if row["VEZF1"] < -0.5 else GRAY
        ax.plot([0, 1], [row["USP34"], row["VEZF1"]], color=color, linewidth=1.3, alpha=0.75, zorder=2)
    ax.scatter(np.zeros(len(wide)), wide["USP34"], s=100, color=GENE_COLORS["USP34"], edgecolor=DGRAY, linewidth=1.2, zorder=3)
    ax.scatter(np.ones(len(wide)), wide["VEZF1"], s=100, color=GENE_COLORS["VEZF1"], edgecolor=DGRAY, linewidth=1.2, zorder=3)

    for line, row in wide.iterrows():
        if row["VEZF1"] < -0.5:
            ax.annotate(line, (1, row["VEZF1"]), xytext=(8, 0), textcoords="offset points", fontsize=8, color=BAD, va="center")

    ax.axhline(-0.5, color=GRAY, linewidth=1.0, linestyle="--", zorder=0)
    ax.text(1.05, -0.5, "informal Chronos<-0.5\nreference line", fontsize=8.5, color=GRAY, va="center", transform=ax.get_yaxis_transform())
    fig.text(0.02, -0.03, "This -0.5 Chronos line is an informal visual reference, distinct from the dependency-probability>0.5 metric used for the % figures in the companion Hany-vs-DepMap map.",
             fontsize=8.3, color=GRAY, style="italic")
    ax.set_xlim(-0.35, 1.65)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["USP34", "VEZF1"], fontsize=13, fontweight="bold")
    ax.get_xticklabels()[0].set_color(GENE_COLORS["USP34"])
    ax.get_xticklabels()[1].set_color(GENE_COLORS["VEZF1"])
    ax.set_ylabel("Chronos gene effect (DepMap 26Q1)", fontsize=11.5)
    _clean_axes(ax)
    _panel_letter(ax, "09b", f"Matched ER+/luminal lines (n={len(wide)})", fontsize=12)
    _save(fig, stub)


# ---------------------------------------------------------------------------
# 10 -- tissue liability context (expression magnitude vs functional liability)
# ---------------------------------------------------------------------------

def build_10_tissue_liability(stub: Path) -> None:
    ctx = pbd.load_normal_tissue_context().set_index("candidate")
    liab = pbd.load_tissue_liability()

    tissues = [("gtex_breast_tpm", "Breast"), ("gtex_blood_tpm", "Whole blood"),
               ("gtex_liver_tpm", "Liver"), ("gtex_heart_tpm", "Heart"), ("gtex_kidney_tpm", "Kidney")]

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(14.5, 6.2), dpi=300, gridspec_kw={"width_ratios": (1, 1.5), "wspace": 1.5})

    y = np.arange(len(tissues))
    for gi, gene in enumerate(["USP34", "VEZF1"]):
        vals = [ctx.loc[gene, col] for col, _ in tissues]
        offset = -0.15 if gene == "USP34" else 0.15
        axL.hlines(y + offset, 0, vals, color=GENE_COLORS[gene], linewidth=2, alpha=0.6, zorder=2)
        axL.scatter(vals, y + offset, s=90, color=GENE_COLORS[gene], edgecolor=DGRAY, linewidth=1.0, zorder=3, label=gene)
    axL.set_yticks(y)
    axL.set_yticklabels([t[1] for t in tissues], fontsize=11)
    axL.set_xlabel("GTEx v8 expression (TPM)", fontsize=10.5)
    _clean_axes(axL)
    axL.legend(loc="lower right", fontsize=10, frameon=False)
    _panel_letter(axL, "A", "Normal-tissue expression (not a liability score)", fontsize=11)

    tier_rank = {"DOCUMENTED_DEVELOPMENTAL": 0, "DOCUMENTED_POSTNATAL_CAUSAL": 0, "DOCUMENTED_HUMAN": 0,
                 "EXPRESSION_ONLY": 1, "INFERRED_ONLY": 1, "INSUFFICIENT_DATA": 2, "NONE_IDENTIFIED": 2}
    tier_color = {0: BAD, 1: "#E69F00", 2: LGRAY}
    systems = sorted(liab["organ_system"].unique())
    yr = np.arange(len(systems))
    for gi, gene in enumerate(["USP34", "VEZF1"]):
        offset = -0.17 if gene == "USP34" else 0.17
        for si, sys_ in enumerate(systems):
            row = liab[(liab["candidate"] == gene) & (liab["organ_system"] == sys_)]
            if len(row) == 0:
                continue
            cls = row.iloc[0]["classification"]
            rank = tier_rank.get(cls, 2)
            axR.scatter([gi], [si + offset], s=170, color=tier_color[rank], edgecolor=GENE_COLORS[gene], linewidth=2.2, zorder=3)
    axR.set_xlim(-0.6, 1.6)
    axR.set_xticks([0, 1])
    axR.set_xticklabels(["USP34", "VEZF1"], fontsize=12, fontweight="bold")
    axR.set_yticks(range(len(systems)))
    axR.set_yticklabels(systems, fontsize=9)
    _clean_axes(axR)
    handles = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor=BAD, markersize=10, label="documented liability (developmental and/or postnatal -- see source table for which)"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor="#E69F00", markersize=10, label="expression-only / inferred"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor=LGRAY, markeredgecolor=GRAY, markersize=10, label="insufficient data / none identified"),
    ]
    axR.legend(handles=handles, loc="upper left", bbox_to_anchor=(1.0, 1.0), fontsize=8.3, frameon=False)
    _panel_letter(axR, "B", "Documented functional-liability evidence (literature)", fontsize=11)

    fig.text(0.5, -0.03, "A developmental/germline phenotype is not automatically a liability of partial adult pharmacological inhibition -- see results/tables/lead_target_deep_dive/USP34_VEZF1_tissue_liability.tsv for the reversibility/zygosity detail behind each red marker.",
             fontsize=8.3, color=GRAY, style="italic", ha="center")
    fig.suptitle("10. Tissue expression vs functional liability -- kept as two separate axes", fontsize=13, fontweight="bold", x=0.02, ha="left", y=1.02)
    _save(fig, stub)


# ---------------------------------------------------------------------------
# 11 -- GDSC pharmacogenomics (top associations + AZD7762 scatter)
# ---------------------------------------------------------------------------

def build_11_gdsc_pharmacogenomics(stub: Path) -> None:
    import scipy.stats as ss

    top = pbd.load_gdsc_top_associations()
    # "top associations" also includes a TOP_EFFECT_SIZE_NOT_NECESSARILY_SIGNIFICANT
    # tier (largest rho regardless of FDR) -- restrict to the genuinely
    # FDR-significant tier so the panel matches its own title (this is where
    # every VEZF1 GDSC row was found to be non-significant during review).
    top = top[(top["metric"] == "LN_IC50") & (top["tier"] == "FDR_SIGNIFICANT")].sort_values("spearman_rho")

    fig = plt.figure(figsize=(13, 6.4), dpi=300)
    gs = fig.add_gridspec(1, 2, width_ratios=(1.3, 1), wspace=0.4)
    axA = fig.add_subplot(gs[0])
    axB = fig.add_subplot(gs[1])

    y = np.arange(len(top))
    colors = [GENE_COLORS[g] for g in top["gene"]]
    axA.hlines(y, 0, top["spearman_rho"], color=colors, linewidth=2.2, alpha=0.6, zorder=2)
    axA.scatter(top["spearman_rho"], y, s=90, color=colors, edgecolor=DGRAY, linewidth=1.0, zorder=3)
    axA.axvline(0, color=GRAY, linewidth=0.8, zorder=0)
    axA.set_yticks(y)
    labels = [f"{r.drug_name} ({r.gene})" for r in top.itertuples()]
    axA.set_yticklabels(labels, fontsize=9)
    axA.set_xlabel("Spearman rho (expression vs LN_IC50, FDR<0.05)", fontsize=10.5)
    _clean_axes(axA)
    _panel_letter(axA, "A", "Top FDR-significant GDSC drug associations", fontsize=11.5)

    df = pfd.load_gdsc_usp34_azd7762()
    stats_row = pfd.load_gdsc_usp34_azd7762_stats()
    x = df["USP34"].to_numpy()
    yv = df["LN_IC50"].to_numpy()
    slope, intercept, r, _, se = ss.linregress(x, yv)
    xs = np.linspace(x.min(), x.max(), 100)
    yhat = slope * xs + intercept
    n = len(x)
    resid = yv - (slope * x + intercept)
    dof = n - 2
    s_err = np.sqrt(np.sum(resid ** 2) / dof)
    conf = ss.t.ppf(0.975, dof) * s_err * np.sqrt(1 / n + (xs - x.mean()) ** 2 / np.sum((x - x.mean()) ** 2))
    axB.fill_between(xs, yhat - conf, yhat + conf, color=GENE_COLORS["USP34"], alpha=0.15, zorder=1)
    axB.plot(xs, yhat, color=GENE_COLORS["USP34"], linewidth=2, zorder=2)
    er = df["is_er_luminal"]
    axB.scatter(x[~er], yv[~er], s=55, color=GRAY, edgecolor=DGRAY, linewidth=0.6, zorder=3, label="other breast lines")
    axB.scatter(x[er], yv[er], s=75, color=GENE_COLORS["USP34"], edgecolor=DGRAY, linewidth=0.8, zorder=4, label="ER+/luminal")
    axB.set_xlabel("USP34 expression (log2 TPM)", fontsize=10.5)
    axB.set_ylabel("AZD7762 LN_IC50", fontsize=10.5)
    _clean_axes(axB)
    axB.legend(loc="upper right", fontsize=8.5, frameon=False)
    _panel_letter(axB, "B", f"USP34 x AZD7762 (rho={stats_row['spearman_rho']:.2f}, FDR={stats_row['fdr']:.3f}, n={n})", fontsize=10.5)

    fig.text(0.01, -0.03, "Pharmacogenomic association only -- not evidence of direct drug-target inhibition. Reported statistic is Spearman rank correlation; the OLS line/band in panel B is a visual trend guide, not the fitted test.",
             fontsize=8.5, color=GRAY, style="italic")
    _save(fig, stub)


# ---------------------------------------------------------------------------
# 12 / 12b -- USP34 structural biology (real PyMOL renders, two variants)
# ---------------------------------------------------------------------------

def _render_structure_bank_panels(out_dir: Path) -> dict[str, Path]:
    paths = pfd.usp34_structure_paths()
    script = Path(__file__).resolve().parent / "poster_structure_render_bank.py"
    pymol_bin = shutil.which("pymol")
    if pymol_bin is None:
        raise RuntimeError("pymol not found on PATH -- install pymol-open-source (see environment.yml)")
    cmd = [pymol_bin, "-cq", str(script), "--", str(paths["7W3R"]), str(paths["7W3U"]), str(out_dir)]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"pymol render failed:\n{result.stdout}\n{result.stderr}")
    logger.info("pymol structure-bank render stdout: %s", result.stdout.strip()[-500:])
    names = ["hero_main", "hero_inset", "comp_apo", "comp_bound", "comp_closeup"]
    panels = {n: out_dir / f"{n}.png" for n in names}
    for p in panels.values():
        assert p.exists(), f"expected pymol output missing: {p}"
    return panels


def build_12_structure_surface(stub: Path) -> None:
    with tempfile.TemporaryDirectory() as tmp:
        panels = _render_structure_bank_panels(Path(tmp))
        main_img = _autocrop_white(plt.imread(panels["hero_main"]))
        inset_img = _autocrop_white(plt.imread(panels["hero_inset"]))

    fig = plt.figure(figsize=(11, 12), dpi=300)
    gs = fig.add_gridspec(2, 1, height_ratios=(5.2, 1.5), hspace=0.05)
    ax = fig.add_subplot(gs[0])
    ax.imshow(main_img)
    ax.axis("off")
    _panel_letter(ax, "12", "USP34 catalytic domain (7W3R) -- surface, region around the catalytic dyad highlighted", fontsize=13)

    # inset placed inside the main axes, lower-left corner (empty white
    # space in the surface render -- does not overlap the molecule or the
    # highlighted cleft patch)
    axIns = ax.inset_axes([0.03, 0.03, 0.36, 0.38])
    axIns.imshow(inset_img)
    axIns.axis("off")
    for spine_pos in ("top", "bottom", "left", "right"):
        axIns.spines[spine_pos].set_visible(True)
        axIns.spines[spine_pos].set_color(DGRAY)
        axIns.spines[spine_pos].set_linewidth(1.2)
    axIns.set_title("covalent linker close-up (7W3U)", fontsize=8.5, color=DGRAY, pad=3)

    axKey = fig.add_subplot(gs[1])
    axKey.axis("off")
    key_items = [("Cys1903 (catalytic)", BAD), ("His2164 (catalytic)", GENE_COLORS["VEZF1"]),
                 ("AYE (covalent probe linker, 7W3U)", GOOD), ("surface within 9A of the dyad (illustrative, not a computed pocket)", "#AFC6E0")]
    # 2x2 grid (not a single row) -- the longest label no longer collides
    # with its neighbor regardless of text length
    col_x = [0.02, 0.52]
    row_y = [0.75, 0.35]
    for i, (label, color) in enumerate(key_items):
        kx, ky = col_x[i % 2], row_y[i // 2]
        axKey.scatter([kx], [ky], s=110, color=color, edgecolor=DGRAY, linewidth=0.8, transform=axKey.transAxes, clip_on=False)
        axKey.text(kx + 0.025, ky, label, fontsize=9, color=DGRAY, va="center", ha="left", transform=axKey.transAxes)
    axKey.text(0.5, -0.1, "Selective USP34 inhibitor: not yet identified. Docking not pursued (no validated ligand set for calibration). The real fpocket-scored pocket (score 0.845) is documented in USP34_pocket_analysis.tsv, not re-derived here.",
               fontsize=8.3, color=GRAY, style="italic", ha="center", transform=axKey.transAxes)

    _save(fig, stub, vector=False)


def build_12b_structure_comparison(stub: Path) -> None:
    with tempfile.TemporaryDirectory() as tmp:
        panels = _render_structure_bank_panels(Path(tmp))
        apo_img = _autocrop_white(plt.imread(panels["comp_apo"]))
        bound_img = _autocrop_white(plt.imread(panels["comp_bound"]))
        closeup_img = _autocrop_white(plt.imread(panels["comp_closeup"]))

    fig = plt.figure(figsize=(12.5, 9.5), dpi=300)
    gs = fig.add_gridspec(2, 2, height_ratios=(1, 1), width_ratios=(1.4, 1), hspace=0.12, wspace=0.05)

    axA = fig.add_subplot(gs[0, 0])
    axA.imshow(apo_img)
    axA.axis("off")
    _panel_letter(axA, "A", "7W3R -- apo (1.92 Å)", fontsize=12)

    axB = fig.add_subplot(gs[1, 0])
    axB.imshow(bound_img)
    axB.axis("off")
    _panel_letter(axB, "B", "7W3U -- covalent ubiquitin-probe complex (3.13 Å)", fontsize=12)

    axC = fig.add_subplot(gs[:, 1])
    axC.imshow(closeup_img)
    axC.axis("off")
    _panel_letter(axC, "C", "Catalytic-cleft close-up (bound)", fontsize=12)

    fig.text(0.02, -0.01, "Matched camera orientation (aligned on the USP34 catalytic domain, chain A). Both structures already frozen/verified in the final_translational phase; re-rendered here for visualization only.",
             fontsize=8.3, color=GRAY, style="italic")
    _save(fig, stub, vector=False)


# ---------------------------------------------------------------------------
# 13 -- validation experiment, biological-glyph schematic (compact)
# ---------------------------------------------------------------------------

def _cell_glyph(ax, cx, cy, r, membrane_color, ko=False, treated=False):
    ax.add_patch(Circle((cx, cy), r, facecolor="white", edgecolor=membrane_color, linewidth=2.6, zorder=2))
    nucleus_color = GRAY if ko else DGRAY
    ax.add_patch(Circle((cx, cy), r * 0.42, facecolor=nucleus_color, alpha=0.25 if ko else 0.4, edgecolor=nucleus_color, linewidth=1.2, zorder=3))
    if ko:
        d = r * 0.3
        ax.plot([cx - d, cx + d], [cy - d, cy + d], color=BAD, linewidth=2.4, zorder=4)
        ax.plot([cx - d, cx + d], [cy + d, cy - d], color=BAD, linewidth=2.4, zorder=4)
    if treated:
        for ang in np.linspace(0, 2 * np.pi, 8, endpoint=False):
            ax.plot([cx + r * np.cos(ang), cx + (r + 0.13) * np.cos(ang)],
                    [cy + r * np.sin(ang), cy + (r + 0.13) * np.sin(ang)], color=GENE_COLORS["VEZF1"], linewidth=1.6, zorder=1)


def build_13_validation_experiment(stub: Path) -> None:
    fig, ax = plt.subplots(figsize=(13.5, 7.3), dpi=300)
    ax.set_xlim(0, 13.5)
    ax.set_ylim(0, 7.5)
    ax.axis("off")

    ax.text(0.0, 1.0, "13. Lead validation experiment -- USP34 (proposed)", transform=ax.transAxes,
            fontsize=13.5, fontweight="bold", color=DGRAY, ha="left", va="top")

    arms = [
        ("Control", False, False), ("4-OHT", False, True),
        ("USP34 KO", True, False), ("USP34 KO\n+ 4-OHT", True, True),
    ]
    axs = [1.6, 4.3, 7.0, 9.7]
    for x, (label, ko, treated) in zip(axs, arms):
        _cell_glyph(ax, x, 6.1, 0.62, GENE_COLORS["USP34"], ko=ko, treated=treated)
        ax.text(x, 5.3, label, ha="center", va="top", fontsize=10.5, fontweight="bold", color=DGRAY, linespacing=1.3)
        ax.add_patch(FancyArrowPatch((x, 4.85), (x, 4.25), arrowstyle="-|>", mutation_scale=14, color=DGRAY, linewidth=1.4))

    readouts = [
        (1.6, "EFFICACY", "viability * clonogenic\nsurvival * dose-response", GOOD),
        (5.35, "MECHANISM", "USP34 loss * AXIN1/\nactive beta-catenin * ESR1 targets", GENE_COLORS["USP34"]),
        (9.7, "COUNTER-\nLIABILITY", "E/N-cadherin * SNAIL *\nEMT / stemness markers", BAD),
    ]
    for x, title, body, color in readouts:
        ax.add_patch(FancyBboxPatch((x - 1.55, 2.55), 3.1, 1.55, boxstyle="round,pad=0.06,rounding_size=0.14",
                                     facecolor=color, alpha=0.14, edgecolor=color, linewidth=1.6, zorder=2))
        ax.text(x, 3.75, title, ha="center", fontsize=9.3, fontweight="bold", color=color, zorder=3)
        ax.text(x, 3.15, body, ha="center", fontsize=7.8, color=DGRAY, zorder=3, linespacing=1.3)

    ax.text(6.65, 1.9, "GOAL: separate direct dependency  vs.  tamoxifen sensitisation  vs.  dual action", ha="center",
            fontsize=10.5, fontweight="bold", color=DGRAY)

    # normal-cell comparators, smaller, off to the side
    _cell_glyph(ax, 12.4, 6.1, 0.42, GRAY, ko=True, treated=False)
    ax.text(12.4, 5.35, "Normal mammary\nepithelial", ha="center", fontsize=8, color=DGRAY, linespacing=1.2)
    _cell_glyph(ax, 12.4, 4.2, 0.42, GRAY, ko=True, treated=False)
    ax.text(12.4, 3.45, "MSC\nosteogenic", ha="center", fontsize=8, color=DGRAY, linespacing=1.2)
    ax.text(12.4, 2.7, "-> selectivity /\nbone-liability\ncomparators", ha="center", fontsize=7.3, color=GRAY, style="italic", linespacing=1.2)

    ax.add_patch(FancyBboxPatch((0.6, 0.35), 8.6, 0.85, boxstyle="round,pad=0.06,rounding_size=0.16",
                                 facecolor="white", edgecolor=GENE_COLORS["VEZF1"], linewidth=1.4, zorder=2))
    ax.text(4.9, 0.775, "secondary follow-up: VEZF1 CRISPRi +/- 4-OHT (smaller, later priority)", ha="center", va="center",
            fontsize=8.5, fontweight="bold", color=GENE_COLORS["VEZF1"], zorder=3)

    _save(fig, stub)


CONTACT_SHEET_ITEMS = [
    ("01a_crispr_ranked_effect", "1a. CRISPR ranked effect"),
    ("01c_crispr_volcano_refined", "1c. CRISPR volcano (refined)"),
    ("03_GSE118713_resistance_landscape", "03. GSE118713 landscape"),
    ("04_GSE240112_recurrence", "04. GSE240112 recurrence"),
    ("05_cross_dataset_candidate_effects", "05. Cross-dataset effects"),
    ("06_resistance_pathway_landscape", "06. Pathway landscape"),
    ("07_candidate_mechanism_map", "07. Mechanism map"),
    ("08_TCGA_human_validation", "08. TCGA validation"),
    ("09_Hany_vs_DepMap_context_map", "09. Hany vs DepMap map"),
    ("09b_USP34_VEZF1_line_dependencies", "09b. USP34/VEZF1 lines"),
    ("10_tissue_liability_context", "10. Tissue liability"),
    ("11_GDSC_USP34_pharmacogenomics", "11. GDSC pharmacogenomics"),
    ("12_USP34_structure_surface", "12. Structure (surface)"),
    ("12b_USP34_structure_comparison", "12b. Structure (comparison)"),
    ("13_validation_experiment", "13. Validation experiment"),
]


def build_contact_sheet(figures_dir: Path = FIGURES, out_png: Path | None = None, ncols: int = 4) -> None:
    """Read-only visual index of every candidate figure -- for selection
    only, not a poster panel itself."""
    out_png = out_png or figures_dir / "POSTER_FIGURE_CONTACT_SHEET.png"
    n = len(CONTACT_SHEET_ITEMS)
    nrows = -(-n // ncols)

    fig, axes = plt.subplots(nrows, ncols, figsize=(5.2 * ncols, 4.6 * nrows), dpi=170)
    axes = np.atleast_2d(axes)
    for i, (stem, label) in enumerate(CONTACT_SHEET_ITEMS):
        r, c = divmod(i, ncols)
        ax = axes[r, c]
        img_path = figures_dir / f"{stem}.png"
        img = plt.imread(img_path)
        ax.imshow(img)
        ax.axis("off")
        ax.set_title(label, fontsize=12, fontweight="bold", color=DGRAY, pad=6)
    for i in range(n, nrows * ncols):
        r, c = divmod(i, ncols)
        axes[r, c].axis("off")

    fig.suptitle("Poster figure bank -- candidate contact sheet (for visual selection only)",
                 fontsize=15, fontweight="bold", y=1.0)
    fig.tight_layout(rect=[0, 0, 1, 0.98])
    fig.savefig(out_png, facecolor="white", bbox_inches="tight", dpi=170)
    plt.close(fig)
    logger.info("wrote %s", out_png)


def run(figures_dir: Path = FIGURES) -> None:
    figures_dir.mkdir(parents=True, exist_ok=True)
    build_01a_crispr_ranked_effect(figures_dir / "01a_crispr_ranked_effect")
    build_01c_crispr_volcano_refined(figures_dir / "01c_crispr_volcano_refined")
    build_03_gse118713_landscape(figures_dir / "03_GSE118713_resistance_landscape")
    build_04_gse240112_recurrence(figures_dir / "04_GSE240112_recurrence")
    build_05_cross_dataset_effects(figures_dir / "05_cross_dataset_candidate_effects")
    build_06_pathway_landscape(figures_dir / "06_resistance_pathway_landscape")
    build_07_candidate_mechanism_map(figures_dir / "07_candidate_mechanism_map")
    build_08_tcga_validation(figures_dir / "08_TCGA_human_validation")
    build_09_hany_vs_depmap(figures_dir / "09_Hany_vs_DepMap_context_map")
    build_09b_usp34_vezf1_line_dependencies(figures_dir / "09b_USP34_VEZF1_line_dependencies")
    build_10_tissue_liability(figures_dir / "10_tissue_liability_context")
    build_11_gdsc_pharmacogenomics(figures_dir / "11_GDSC_USP34_pharmacogenomics")
    build_12_structure_surface(figures_dir / "12_USP34_structure_surface")
    build_12b_structure_comparison(figures_dir / "12b_USP34_structure_comparison")
    build_13_validation_experiment(figures_dir / "13_validation_experiment")
    build_contact_sheet(figures_dir)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run()
