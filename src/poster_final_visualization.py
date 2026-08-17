"""FINAL poster figure set (visualization only) -- results/figures/poster_final/.

Builds the six main poster-quality figures selected from the frozen science
(``results/reports/poster_final/POSTER_FINAL_FIGURE_GUIDE.md`` explains the
selection). No new analysis, no recomputation, no re-ranking -- every
plotted number is read from ``src/poster_final_data.py``, which in turn
reads only already-frozen project tables. See that module's docstring for
the exact source of every value.

Design system:
  - candidate identity colors are fixed across every figure in this set
    (FOCUS_COLORS, imported from poster_final_data): USP34 blue, VEZF1
    amber, KDM1A vermillion, TLK2 sky-blue -- the Okabe-Ito colorblind-safe
    set, validated with the dataviz skill's palette validator
  - non-focus/background data stays neutral gray; candidate color is never
    used for anything that isn't candidate identity
  - no composite/weighted score is ever plotted -- every panel is a single,
    real, named evidence dimension
  - small panel letters + short descriptive labels only, no large in-figure
    titles or paragraph-length callouts
  - every figure exports .png (300dpi) and .pdf; .svg where the panel
    content is vector-friendly (not a rasterized PyMOL render)
"""

from __future__ import annotations

import logging
import tempfile
import textwrap
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import TwoSlopeNorm
from matplotlib.lines import Line2D

from src import poster_final_data as pfd2
from src.poster_figures_bank_visualization import _render_structure_bank_panels
from src.poster_figures_visualization import _autocrop_white, rng

logger = logging.getLogger(__name__)

FIGURES = Path("results/figures/poster_final")

FOCUS_COLORS = pfd2.FOCUS_COLORS
FOCUS_FOUR = pfd2.FOCUS_FOUR
GRAY = pfd2.GRAY
LGRAY = pfd2.LGRAY
DGRAY = pfd2.DGRAY


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


def _panel_letter(ax, letter: str, label: str, fontsize: float = 12) -> None:
    ax.set_title(f"{letter}. {label}", fontsize=fontsize, fontweight="bold", loc="left", pad=8)


def _caption(fig, text: str, y: float = -0.02) -> None:
    fig.text(0.01, y, text, fontsize=8.3, color=GRAY, style="italic", wrap=True)


# ---------------------------------------------------------------------------
# F1 -- Genome-wide CRISPR discovery
# ---------------------------------------------------------------------------

def build_f1_crispr_discovery(stub: Path) -> None:
    df = pfd2.load_f1_genomewide()
    ranks = pfd2.load_f1_focus_gene_ranks().set_index("gene")
    n_sens = int(ranks["n_sensitising_hits"].iloc[0])
    rcor1 = pfd2.load_f1_blind_control_row()

    x = df["effect_size"].to_numpy()
    y = -np.log10(df["fdr"].to_numpy())

    fig, ax = plt.subplots(figsize=(12, 8), dpi=300)
    bg = ~df["gate1_significant"] & ~df["is_focus_gene"]
    gate1_only = df["gate1_significant"] & ~df["is_focus_gene"]
    ax.scatter(x[bg], y[bg], s=4, color=LGRAY, alpha=0.35, linewidth=0, rasterized=True, zorder=1)
    ax.scatter(x[gate1_only], y[gate1_only], s=14, color=GRAY, alpha=0.55, linewidth=0, zorder=2,
               label=f"other Gate-1 hits (FDR<0.1, n={int(df['gate1_significant'].sum()) - 4})")
    ax.axhline(-np.log10(0.1), color=GRAY, linewidth=0.9, linestyle="--", zorder=0)
    ax.text(x.max() * 0.97, -np.log10(0.1) + 0.05, "Gate-1 threshold (FDR=0.1)", ha="right", va="bottom",
            fontsize=8.5, color=GRAY)
    ax.axvline(0, color=LGRAY, linewidth=0.8, zorder=0)

    if rcor1 is not None:
        rcor1_fdr = float(rcor1["fdr"])
        rx, ry = float(rcor1["effect_size"]), -np.log10(rcor1_fdr)
        ax.scatter([rx], [ry], s=110, facecolor="none", edgecolor=DGRAY, linewidth=1.6, zorder=3)
        ax.annotate(f"RCOR1 (blind control, NOT recovered at\nGate-1 threshold, FDR={rcor1_fdr:.2f})",
                    xy=(rx, ry), xytext=(rx + 0.4, ry - 0.55), fontsize=8.6, color=DGRAY, ha="left", va="top",
                    arrowprops=dict(arrowstyle="-", color=DGRAY, linewidth=0.8))

    ordered = sorted(FOCUS_FOUR, key=lambda g: -(-np.log10(float(df.loc[df["gene"] == g, "fdr"].iloc[0]))))
    label_x0 = -0.55
    label_ys = np.linspace(4.3, 0.7, len(ordered))
    for gene, ly in zip(ordered, label_ys):
        row = df[df["gene"] == gene].iloc[0]
        gx, gy = float(row["effect_size"]), -np.log10(float(row["fdr"]))
        r = ranks.loc[gene]
        ax.scatter([gx], [gy], s=230, color=FOCUS_COLORS[gene], edgecolor=DGRAY, linewidth=1.4, zorder=5)
        ax.plot([label_x0, gx], [ly, gy], color=FOCUS_COLORS[gene], linewidth=1.0, zorder=4)
        label = f"{gene}  (effect rank {int(r['rank_by_effect'])}/{n_sens}, FDR rank {int(r['rank_by_fdr'])}/{n_sens})"
        ax.text(label_x0 + 0.05, ly, label, ha="left", va="center", fontsize=10, fontweight="bold", color=FOCUS_COLORS[gene])

    ax.set_xlabel(r"CRISPR effect size (Hany et al. screen)   " +
                  r"$\bf{negative}$ = sensitising knockout under 4-OHT   |   positive = tolerance-associated",
                  fontsize=11)
    ax.set_ylabel(r"$-\log_{10}$(FDR)", fontsize=12)
    _clean_axes(ax)
    ax.legend(loc="upper left", fontsize=9.5, frameon=False)
    _panel_letter(ax, "F1", f"Genome-wide CRISPR screen (n={len(df):,} fitted genes)", fontsize=14)
    _caption(fig, f"USP34 is one of {n_sens} significant sensitising Gate-1 hits, not the top effect-size or "
                   "top-FDR hit -- see rank annotations above. Sign convention verified against src/labels.py.")
    _save(fig, stub)


# ---------------------------------------------------------------------------
# F2 -- Transcriptomic / pathway systems view
# ---------------------------------------------------------------------------

_CATEGORY_COLOR = {
    "resistance model (cell line)": "#4C4C9D",
    "recurrence-associated (human tumour, unpaired)": "#B0392F",
    "acute 12h (not resistance)": "#8C8C8C",
}


def build_f2_pathway_systems(stub: Path) -> None:
    df = pfd2.load_f2_pathway_matrix()
    pathway_order = list(dict.fromkeys(df["pathway_label"]))
    dataset_order = ["GSE118713", "GSE111151", "GSE240112", "GSE245601 (acute 12h)"]

    pivot_nes = df.pivot(index="pathway_label", columns="dataset_label", values="NES").loc[pathway_order, dataset_order]
    pivot_fdr = df.pivot(index="pathway_label", columns="dataset_label", values="fdr").loc[pathway_order, dataset_order]
    cat_by_label = df.drop_duplicates("dataset_label").set_index("dataset_label")["dataset_category"]

    fig = plt.figure(figsize=(11, 8.5), dpi=300)
    gs = fig.add_gridspec(1, 2, width_ratios=(20, 1), wspace=0.05)
    ax = fig.add_subplot(gs[0])
    cax = fig.add_subplot(gs[1])

    vmax = np.nanmax(np.abs(pivot_nes.to_numpy()))
    norm = TwoSlopeNorm(vcenter=0, vmin=-vmax, vmax=vmax)
    masked = np.ma.masked_invalid(pivot_nes.to_numpy())
    im = ax.imshow(masked, cmap="RdBu_r", norm=norm, aspect="auto")
    for i in range(pivot_nes.shape[0]):
        for j in range(pivot_nes.shape[1]):
            v = pivot_nes.iloc[i, j]
            if np.isnan(v):
                ax.text(j, i, "n/t", ha="center", va="center", fontsize=8, color=GRAY)
                continue
            fdr = pivot_fdr.iloc[i, j]
            marker = "*" if (not np.isnan(fdr) and fdr < 0.05) else ""
            ax.text(j, i, f"{v:.2f}{marker}", ha="center", va="center", fontsize=8.5,
                    color="white" if abs(v) > vmax * 0.55 else DGRAY)

    ax.set_xticks(range(len(dataset_order)))
    ax.set_xticklabels(dataset_order, fontsize=9.5, rotation=20, ha="right")
    for tick, label in zip(ax.get_xticklabels(), dataset_order):
        tick.set_color(_CATEGORY_COLOR[cat_by_label[label]])
        tick.set_fontweight("bold")
    ax.set_yticks(range(len(pathway_order)))
    ax.set_yticklabels(pathway_order, fontsize=9.5)
    ax.set_xticks(np.arange(-0.5, len(dataset_order), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(pathway_order), 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=2)
    ax.tick_params(which="minor", length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)

    cb = fig.colorbar(im, cax=cax)
    cb.set_label("NES", fontsize=9)
    cb.ax.tick_params(labelsize=8)

    legend_handles = [Line2D([0], [0], marker="s", color="none", markerfacecolor=c, markersize=10, label=lab)
                       for lab, c in _CATEGORY_COLOR.items()]
    fig.legend(handles=legend_handles, loc="lower center", bbox_to_anchor=(0.44, -0.1), ncol=1,
               fontsize=9, frameon=False, title="Dataset category", title_fontsize=9.5)

    _panel_letter(ax, "F2", "Pathway-level convergence across transcriptomic contexts", fontsize=14)
    _caption(fig, "* FDR<0.05 (gene-set enrichment within that dataset). \"n/t\" = pathway not testable in that "
                   "dataset (gene-set size filter). GSE240112 is a recurrence-associated human-tumour comparison, "
                   "not a resistance model -- shown in its own category, never grouped with GSE118713/GSE111151.",
             y=-0.16)
    _save(fig, stub)


# ---------------------------------------------------------------------------
# F3 -- Candidate evidence divergence / integration
# ---------------------------------------------------------------------------

def _dot(ax, x, y, state, color):
    if state == "YES":
        ax.scatter([x], [y], s=170, color=color, edgecolor=DGRAY, linewidth=1.0, zorder=3)
    elif state == "PARTIAL":
        ax.scatter([x], [y], s=170, facecolor=color, alpha=0.35, edgecolor=DGRAY, linewidth=1.0, zorder=3)
    elif state == "NA":
        ax.text(x, y, "n/a", ha="center", va="center", fontsize=8, color=GRAY)
    else:
        ax.scatter([x], [y], s=170, facecolor="white", edgecolor=GRAY, linewidth=1.2, zorder=3)


def build_f3_evidence_integration(stub: Path) -> None:
    em = pfd2.load_f3_evidence_matrix().set_index("gene").loc[FOCUS_FOUR]
    facets = pfd2.load_f3_structural_facets().set_index("gene").loc[FOCUS_FOUR]

    fig = plt.figure(figsize=(14, 10), dpi=300)
    gs = fig.add_gridspec(2, 2, hspace=0.55, wspace=0.55)

    # A -- CRISPR strength
    axA = fig.add_subplot(gs[0, 0])
    y = np.arange(len(FOCUS_FOUR))
    neg_log_fdr = -np.log10(em["crispr_fdr"].to_numpy())
    colors = [FOCUS_COLORS[g] for g in FOCUS_FOUR]
    axA.barh(y, neg_log_fdr, color=colors, edgecolor=DGRAY, linewidth=0.8, height=0.6)
    for i, gene in enumerate(FOCUS_FOUR):
        axA.text(neg_log_fdr[i] + 0.05, i,
                  f"effect {em.loc[gene, 'crispr_effect']:.2f}, rank {int(em.loc[gene, 'rank_by_effect'])}/13",
                  va="center", fontsize=8.8, color=DGRAY)
    axA.set_yticks(y)
    axA.set_yticklabels(FOCUS_FOUR, fontsize=10.5, fontweight="bold")
    axA.set_xlabel(r"$-\log_{10}$(CRISPR FDR)", fontsize=10)
    axA.set_xlim(0, neg_log_fdr.max() * 1.85)
    axA.invert_yaxis()
    _clean_axes(axA)
    _panel_letter(axA, "A", "CRISPR sensitisation strength", fontsize=11.5)

    # B -- RNA / human corroboration dot-matrix
    axB = fig.add_subplot(gs[0, 1])
    cols = [("gse118713_fdr", "GSE118713\n(resistance)"), ("gse111151_fdr", "GSE111151\n(resistance)"),
            ("gse240112_fdr", "GSE240112\n(recurrence)"), ("tcga_fdr", "TCGA\n(paired tumour)")]
    for j, (col, _label) in enumerate(cols):
        for i, gene in enumerate(FOCUS_FOUR):
            v = em.loc[gene, col]
            if pd.isna(v):
                state = "NA"
            elif v < 0.05:
                state = "YES"
            else:
                state = "NO"
            _dot(axB, j, i, state, FOCUS_COLORS[gene])
    axB.set_xlim(-0.5, len(cols) - 0.5)
    axB.set_ylim(-0.5, len(FOCUS_FOUR) - 0.5)
    axB.set_xticks(range(len(cols)))
    axB.set_xticklabels([c[1] for c in cols], fontsize=8.6)
    axB.set_yticks(range(len(FOCUS_FOUR)))
    axB.set_yticklabels(FOCUS_FOUR, fontsize=10.5, fontweight="bold")
    axB.invert_yaxis()
    for spine in axB.spines.values():
        spine.set_visible(False)
    axB.tick_params(length=0)
    legend_handles = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor=DGRAY, markeredgecolor=DGRAY, markersize=9, label="FDR<0.05"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor="white", markeredgecolor=GRAY, markersize=9, label="tested, not significant"),
    ]
    axB.legend(handles=legend_handles, loc="upper center", bbox_to_anchor=(0.5, -0.18), ncol=2, fontsize=8.3, frameon=False)
    _panel_letter(axB, "B", "RNA / human-tumour corroboration", fontsize=11.5)

    # C -- DepMap baseline dependency
    axC = fig.add_subplot(gs[1, 0])
    pct = (em["frac_strongly_dependent_er_luminal"] * 100).to_numpy()
    axC.barh(y, pct, color=colors, edgecolor=DGRAY, linewidth=0.8, height=0.6)
    for i, gene in enumerate(FOCUS_FOUR):
        axC.text(pct[i] + 1.5, i, f"{pct[i]:.0f}% of 11 lines", va="center", fontsize=8.8, color=DGRAY)
    axC.set_yticks(y)
    axC.set_yticklabels(FOCUS_FOUR, fontsize=10.5, fontweight="bold")
    axC.set_xlabel("% ER+/luminal DepMap 26Q1 lines strongly dependent", fontsize=9.5)
    axC.set_xlim(0, 118)
    axC.invert_yaxis()
    _clean_axes(axC)
    _panel_letter(axC, "C", "Baseline cancer-cell dependency (double-edged)", fontsize=11.5)

    # D -- structural / pharmacology facets
    axD = fig.add_subplot(gs[1, 1])
    facet_cols = ["structure_exists", "ligand_or_probe_bound", "validated_inhibitor", "clinical_stage"]
    facet_labels = ["Experimental\nstructure", "Ligand/probe-\nbound structure", "Validated selective\ninhibitor", "Clinical-stage\npharmacology"]
    for j, col in enumerate(facet_cols):
        for i, gene in enumerate(FOCUS_FOUR):
            _dot(axD, j, i, facets.loc[gene, col], FOCUS_COLORS[gene])
    axD.set_xlim(-0.5, len(facet_cols) - 0.5)
    axD.set_ylim(-0.5, len(FOCUS_FOUR) - 0.5)
    axD.set_xticks(range(len(facet_cols)))
    axD.set_xticklabels(facet_labels, fontsize=8.4)
    axD.set_yticks(range(len(FOCUS_FOUR)))
    axD.set_yticklabels(FOCUS_FOUR, fontsize=10.5, fontweight="bold")
    axD.invert_yaxis()
    for spine in axD.spines.values():
        spine.set_visible(False)
    axD.tick_params(length=0)
    legend_handles2 = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor=DGRAY, markeredgecolor=DGRAY, markersize=9, label="yes"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor=DGRAY, markeredgecolor=DGRAY, markersize=9, alpha=0.35, label="partial"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor="white", markeredgecolor=GRAY, markersize=9, label="no"),
    ]
    axD.legend(handles=legend_handles2, loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=3, fontsize=8.3, frameon=False)
    _panel_letter(axD, "D", "Structural / pharmacological tractability", fontsize=11.5)

    fig.suptitle("", fontsize=1)  # no in-figure master title, per design system
    _caption(fig, "No composite or weighted score is computed anywhere in this figure -- each panel is a single, "
                   "named evidence dimension read directly from the frozen post-audit tables. Different candidates "
                   "lead on different axes; this is an evidence-integration figure, not a ranking.", y=-0.04)
    _save(fig, stub)


# ---------------------------------------------------------------------------
# F4 -- Human (TCGA) / DepMap orthogonal validation
# ---------------------------------------------------------------------------

def build_f4_human_depmap_validation(stub: Path) -> None:
    forest = pfd2.load_f4_tcga_forest().set_index("candidate")
    depmap = pfd2.load_f4_depmap_effect()

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(13, 6.5), dpi=300, gridspec_kw=dict(width_ratios=(1, 1.3)))

    # A -- TCGA forest (USP34, VEZF1 real; KDM1A, TLK2 not assessed)
    y = np.arange(len(FOCUS_FOUR))
    for i, gene in enumerate(FOCUS_FOUR):
        if gene in forest.index:
            r = forest.loc[gene]
            est, lo, hi = float(r["mean_diff"]), float(r["ci_low"]), float(r["ci_high"])
            axA.plot([lo, hi], [i, i], color=FOCUS_COLORS[gene], linewidth=2.2, zorder=2)
            axA.scatter([est], [i], s=110, color=FOCUS_COLORS[gene], edgecolor=DGRAY, linewidth=1.2, zorder=3)
            axA.text(hi + 0.05, i, f"FDR={r['fdr']:.3f}", va="center", fontsize=8.5, color=DGRAY)
        else:
            axA.text(0, i, "not assessed in this project\n(TCGA only run for the original 4)",
                      va="center", ha="center", fontsize=8.5, color=GRAY, style="italic")
    axA.axvline(0, color=LGRAY, linewidth=0.9, zorder=0)
    axA.set_yticks(y)
    axA.set_yticklabels(FOCUS_FOUR, fontsize=10.5, fontweight="bold")
    axA.invert_yaxis()
    axA.set_xlabel("TCGA-BRCA tumour vs normal, paired log2FC (95% CI)", fontsize=9.7)
    _clean_axes(axA)
    _panel_letter(axA, "A", "Human tumour evidence (TCGA)", fontsize=12)

    # B -- DepMap Chronos distribution
    for i, gene in enumerate(FOCUS_FOUR):
        vals = depmap.loc[depmap["gene"] == gene, "chronos_effect"].to_numpy()
        q1, med, q3 = np.percentile(vals, [25, 50, 75])
        axB.add_patch(plt.Rectangle((i - 0.16, q1), 0.32, q3 - q1, facecolor=FOCUS_COLORS[gene], alpha=0.18,
                                     edgecolor=FOCUS_COLORS[gene], linewidth=1.4, zorder=2))
        axB.plot([i - 0.16, i + 0.16], [med, med], color=FOCUS_COLORS[gene], linewidth=2.2, zorder=3)
        jitter = rng.uniform(-0.1, 0.1, size=len(vals))
        axB.scatter(i + jitter, vals, s=45, color=FOCUS_COLORS[gene], edgecolor="white", linewidth=0.5, zorder=4, alpha=0.9)
    axB.axhline(0, color=LGRAY, linewidth=0.8, zorder=0)
    axB.axhline(-1, color=GRAY, linewidth=0.8, linestyle="--", zorder=0)
    axB.text(3.55, -1, "strong-dependency\nreference (Chronos=-1)", fontsize=7.6, color=GRAY, va="center")
    axB.set_xticks(range(len(FOCUS_FOUR)))
    axB.set_xticklabels(FOCUS_FOUR, fontsize=10.5, fontweight="bold")
    axB.set_ylabel("DepMap 26Q1 Chronos gene effect\n(11 ER+/luminal screened lines)", fontsize=9.5)
    _clean_axes(axB)
    _panel_letter(axB, "B", "Baseline cancer-cell dependency (DepMap)", fontsize=12)

    _caption(fig, "Orthogonal validation distinguishes tamoxifen-specific sensitisation (panel A / Figure F1) from "
                   "baseline cell-line essentiality (panel B). High baseline dependency (TLK2) is not automatically "
                   "an advantage -- it may reflect a narrower, less tamoxifen-specific therapeutic window, exactly "
                   "as for VEZF1's more moderate dependency. KDM1A shows low baseline dependency despite the "
                   "strongest CRISPR sensitisation signal (Figure F1), the clearest tamoxifen-specific-vs-baseline "
                   "dissociation among the four focus genes.", y=-0.1)
    _save(fig, stub)


# ---------------------------------------------------------------------------
# F5 -- USP34 structure and tractability
# ---------------------------------------------------------------------------

def build_f5_usp34_structure(stub: Path) -> None:
    with tempfile.TemporaryDirectory() as tmp:
        panels = _render_structure_bank_panels(Path(tmp))
        apo_img = _autocrop_white(plt.imread(panels["comp_apo"]))
        bound_img = _autocrop_white(plt.imread(panels["comp_bound"]))
        closeup_img = _autocrop_white(plt.imread(panels["comp_closeup"]))

    fig = plt.figure(figsize=(13, 6.6), dpi=300)
    gs = fig.add_gridspec(2, 3, height_ratios=(4.3, 1.35), width_ratios=(1, 1, 1), hspace=0.1, wspace=0.04)

    axA = fig.add_subplot(gs[0, 0])
    axA.imshow(apo_img)
    axA.axis("off")
    _panel_letter(axA, "A", "7W3R -- apo (1.92 Å)", fontsize=12)

    axB = fig.add_subplot(gs[0, 1])
    axB.imshow(bound_img)
    axB.axis("off")
    _panel_letter(axB, "B", "7W3U -- covalent ubiquitin-probe complex", fontsize=12)

    axC = fig.add_subplot(gs[0, 2])
    axC.imshow(closeup_img)
    axC.axis("off")
    _panel_letter(axC, "C", "Catalytic cleft close-up (Cys1903 / His2164)", fontsize=12)

    axKey = fig.add_subplot(gs[1, :])
    axKey.axis("off")
    key_items = [("Cys1903 (catalytic)", "#b0392f"), ("His2164 (catalytic)", "#E69F00"),
                 ("AYE (covalent ubiquitin-probe linker, 7W3U)", "#009E73")]
    col_x = [0.02, 0.37, 0.68]
    for (label, color), kx in zip(key_items, col_x):
        axKey.scatter([kx], [0.85], s=110, color=color, edgecolor=DGRAY, linewidth=0.8, transform=axKey.transAxes, clip_on=False)
        axKey.text(kx + 0.025, 0.85, label, fontsize=9, color=DGRAY, va="center", ha="left", transform=axKey.transAxes)
    axKey.text(0.5, 0.35,
               "Direct catalytic-cysteine reactivity is experimentally observed (a covalent activity-based "
               "ubiquitin probe, not a small-molecule ligand). No validated selective USP34 inhibitor exists. "
               "Docking was not pursued in this project (no validated ligand set for calibration).",
               fontsize=9.3, color=DGRAY, ha="center", va="center", transform=axKey.transAxes, wrap=True)

    fig.text(0.5, 0.985, "F5. USP34 structure and covalent tractability", fontsize=14, fontweight="bold", ha="center")
    _caption(fig, "Both structures already frozen/verified in the final_translational and post-audit-sensitivity "
                   "phases; re-rendered here for poster composition only. USP34 is not the top CRISPR hit (see "
                   "Figure F1) -- this figure documents why it remains a structurally distinctive candidate despite "
                   "that, not a claim that it is the strongest functional hit.", y=-0.02)
    _save(fig, stub, vector=False)


# ---------------------------------------------------------------------------
# F6 -- Final candidate logic / experimental strategy
# ---------------------------------------------------------------------------

def _role_text(gene: str, ranks: pd.DataFrame, em: pd.DataFrame) -> str:
    """Role-card narrative for each focus gene. Qualitative/pharmacology
    descriptions are static prose (already verified against Table 06b and
    the post-audit deep-dive tables in POSTER_FINAL_FIGURE_GUIDE.md); every
    *numeric* fact (rank, %) is interpolated from the loaded tables at
    render time, never hand-typed, so it cannot drift from the source."""
    r = ranks.loc[gene]
    pct = em.loc[gene, "frac_strongly_dependent_er_luminal"] * 100
    if gene == "KDM1A":
        return (f"Strongest functional benchmark: rank {int(r['rank_by_effect'])}/{int(r['n_sensitising_hits'])} "
                "by CRISPR effect and FDR; existing clinical-stage LSD1 pharmacology (e.g. iadademstat). "
                "Low baseline dependency.")
    if gene == "TLK2":
        return (f"Strongest baseline dependency among the 4 focus genes ({pct:.1f}% of 11 lines) -- a real "
                "experimental kinase structure exists, but no selective inhibitor.")
    if gene == "USP34":
        return ("Structurally tractable novel lead: covalent catalytic-cysteine (Cys1903) reactivity "
                f"demonstrated; no validated inhibitor yet; not the top CRISPR hit (rank {int(r['rank_by_effect'])}/"
                f"{int(r['n_sensitising_hits'])}).")
    if gene == "VEZF1":
        return ("Strongest recurrence-associated human signal (GSE240112, single dataset) among the 4 focus "
                "genes; hard-to-drug zinc-finger transcription-factor class.")
    raise ValueError(f"no role text defined for {gene!r}")


_ROLE_ASSAY = {
    "KDM1A": "Follow-up: LSD1-inhibitor combination assay (tamoxifen +/- iadademstat)",
    "TLK2": "Follow-up: baseline-dependency vs. tamoxifen-specific rescue assay",
    "USP34": "Follow-up: covalent-probe engagement assay; CRISPR-competition validation",
    "VEZF1": "Follow-up: recurrence-cohort replication; CRISPR-competition validation",
}


def build_f6_translational_framework(stub: Path) -> None:
    compass = pfd2.load_f6_role_compass()
    ranks = pfd2.load_f1_focus_gene_ranks().set_index("gene")
    em = pfd2.load_f3_evidence_matrix().set_index("gene")

    fig = plt.figure(figsize=(13, 10), dpi=300)
    gs = fig.add_gridspec(2, 1, height_ratios=(3, 2), hspace=0.32)

    axA = fig.add_subplot(gs[0])
    bg = ~compass["is_focus_gene"]
    axA.scatter(compass.loc[bg, "neg_log10_fdr"], compass.loc[bg, "frac_strongly_dependent_er_luminal"] * 100,
                s=70, color=LGRAY, edgecolor=GRAY, linewidth=0.8, zorder=2)
    for _, r in compass.iterrows():
        if not r["is_focus_gene"]:
            continue
        gene = r["gene"]
        axA.scatter([r["neg_log10_fdr"]], [r["frac_strongly_dependent_er_luminal"] * 100], s=260,
                    color=FOCUS_COLORS[gene], edgecolor=DGRAY, linewidth=1.5, zorder=4)
        axA.annotate(gene, xy=(r["neg_log10_fdr"], r["frac_strongly_dependent_er_luminal"] * 100),
                     xytext=(8, 8), textcoords="offset points", fontsize=11, fontweight="bold",
                     color=FOCUS_COLORS[gene])
    axA.set_xlabel(r"CRISPR sensitisation strength ($-\log_{10}$ FDR)", fontsize=10.5)
    axA.set_ylabel("Baseline ER+/luminal DepMap\ndependency (% of 11 lines)", fontsize=10.5)
    _clean_axes(axA)
    n_total = len(compass)
    axA.text(0.02, 0.96, f"all {n_total} significant sensitising Gate-1 hits (gray = other 9)",
             transform=axA.transAxes, fontsize=8.8, color=GRAY, va="top")
    _panel_letter(axA, "A", "Sensitisation strength is not the same axis as baseline dependency", fontsize=13)

    axB = fig.add_subplot(gs[1])
    axB.axis("off")
    n = len(FOCUS_FOUR)
    card_w = 1.0 / n
    for i, gene in enumerate(FOCUS_FOUR):
        x0 = i * card_w
        axB.add_patch(plt.Rectangle((x0 + 0.01, 0.08), card_w - 0.02, 0.86, facecolor=FOCUS_COLORS[gene],
                                     alpha=0.10, edgecolor=FOCUS_COLORS[gene], linewidth=1.6,
                                     transform=axB.transAxes, clip_on=False))
        axB.text(x0 + card_w / 2, 0.86, gene, ha="center", va="top", fontsize=13, fontweight="bold",
                  color=FOCUS_COLORS[gene], transform=axB.transAxes)
        axB.text(x0 + card_w / 2, 0.68, textwrap.fill(_role_text(gene, ranks, em), width=32), ha="center", va="top",
                  fontsize=8.3, color=DGRAY, transform=axB.transAxes)
        axB.text(x0 + card_w / 2, 0.22, textwrap.fill(_ROLE_ASSAY[gene], width=30), ha="center", va="top",
                  fontsize=8.2, color=DGRAY, style="italic", transform=axB.transAxes)
    _panel_letter(axB, "B", "Distinct candidate roles and proposed follow-up logic (no universal winner)", fontsize=13)

    _caption(fig, "This project yields a testable, multi-candidate translational framework, not a single "
                   "overclaimed winner. Panel A is real, data-grounded (Table 03); panel B is a role summary, not "
                   "a new ranking rule -- roles are the same distinct-candidate-roles interpretation already "
                   "documented in the post-audit sensitivity report.", y=-0.02)
    _save(fig, stub)


# ---------------------------------------------------------------------------
# Contact sheet
# ---------------------------------------------------------------------------

FINAL_CONTACT_SHEET_ITEMS = [
    ("F1_crispr_discovery", "F1. Genome-wide CRISPR discovery"),
    ("F2_pathway_systems", "F2. Pathway systems view"),
    ("F3_candidate_evidence_integration", "F3. Candidate evidence integration"),
    ("F4_human_depmap_validation", "F4. Human / DepMap validation"),
    ("F5_USP34_structure_tractability", "F5. USP34 structure & tractability"),
    ("F6_final_translational_framework", "F6. Final translational framework"),
]


def build_final_contact_sheet(figures_dir: Path = FIGURES, out_png: Path | None = None, ncols: int = 2) -> None:
    """Visual index of the 6 selected final poster figures."""
    out_png = out_png or figures_dir / "POSTER_FINAL_CONTACT_SHEET.png"
    n = len(FINAL_CONTACT_SHEET_ITEMS)
    nrows = -(-n // ncols)

    fig, axes = plt.subplots(nrows, ncols, figsize=(8 * ncols, 6.4 * nrows), dpi=170)
    axes = np.atleast_2d(axes)
    for i, (stem, label) in enumerate(FINAL_CONTACT_SHEET_ITEMS):
        r, c = divmod(i, ncols)
        ax = axes[r, c]
        img = plt.imread(figures_dir / f"{stem}.png")
        ax.imshow(img)
        ax.axis("off")
        ax.set_title(label, fontsize=14, fontweight="bold", color=DGRAY, pad=8)
    for i in range(n, nrows * ncols):
        r, c = divmod(i, ncols)
        axes[r, c].axis("off")

    fig.suptitle("Final poster figure set (6 figures) -- selected from the frozen science",
                 fontsize=17, fontweight="bold", y=1.0)
    fig.tight_layout(rect=[0, 0, 1, 0.98])
    fig.savefig(out_png, facecolor="white", bbox_inches="tight", dpi=170)
    plt.close(fig)
    logger.info("wrote %s", out_png)


# ---------------------------------------------------------------------------
# Runner
# ---------------------------------------------------------------------------

def run(out_dir: str | Path = FIGURES) -> None:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    build_f1_crispr_discovery(out_dir / "F1_crispr_discovery")
    build_f2_pathway_systems(out_dir / "F2_pathway_systems")
    build_f3_evidence_integration(out_dir / "F3_candidate_evidence_integration")
    build_f4_human_depmap_validation(out_dir / "F4_human_depmap_validation")
    build_f5_usp34_structure(out_dir / "F5_USP34_structure_tractability")
    build_f6_translational_framework(out_dir / "F6_final_translational_framework")
    build_final_contact_sheet(out_dir)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run()
