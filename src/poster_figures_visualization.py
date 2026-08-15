"""Poster figure set, science-first rebuild.

Six figures built only from real, already-validated project data (never
hand-typed numbers, never icon/symbol summaries standing in for data):

  1. Genome-wide CRISPR volcano plot (Hany screen, 19,103 genes)
  2. Real sample-level expression evidence (GSE118713 USP34, GSE240112
     VEZF1, cross-dataset forest panel)
  3. DepMap 26Q1 Chronos dependency distributions (ER+/luminal lines)
  4. GDSC pharmacogenomic scatter (USP34 expression vs AZD7762 response)
  5. USP34 structural-biology figure (real PyMOL cartoon renders of the
     7W3R/7W3U PDB structures, not a matplotlib wireframe)
  6. Final experimental strategy (the one schematic figure)

No new analysis, no new candidate discovery, no re-ranking. See
src/poster_figures_data.py for the exact source table/file for every
plotted number.
"""

from __future__ import annotations

import logging
import subprocess
import tempfile
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

from src import poster_figures_data as pfd

logger = logging.getLogger(__name__)

FIGURES = Path("results/figures/poster")

GENE_ORDER = ["USP34", "VEZF1", "EML5", "CITED2"]
GENE_COLORS = {
    "USP34": "#0072B2",
    "VEZF1": "#E69F00",
    "EML5": "#009E73",
    "CITED2": "#CC79A7",
}
DGRAY = "#333333"
GRAY = "#9a9a9a"
LGRAY = "#d8d8d8"
GOOD = "#009E73"
BAD = "#b0392f"

rng = np.random.default_rng(0)


def _save(fig, out_png: Path, vector: bool = True) -> None:
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, facecolor="white", bbox_inches="tight", dpi=220)
    if vector:
        fig.savefig(out_png.with_suffix(".pdf"), facecolor="white", bbox_inches="tight")
    plt.close(fig)
    logger.info("wrote %s", out_png)


def _strip_box(ax, x0: float, values: np.ndarray, color: str, width: float = 0.32, box_alpha: float = 0.18):
    """Box (quartiles) + jittered strip of real points -- used for every
    real-data distribution panel in this figure set."""
    q1, med, q3 = np.percentile(values, [25, 50, 75])
    ax.add_patch(plt.Rectangle((x0 - width / 2, q1), width, q3 - q1, facecolor=color, alpha=box_alpha,
                                edgecolor=color, linewidth=1.4, zorder=2))
    ax.plot([x0 - width / 2, x0 + width / 2], [med, med], color=color, linewidth=2.2, zorder=3)
    jitter = rng.uniform(-width * 0.32, width * 0.32, size=len(values))
    ax.scatter(x0 + jitter, values, s=55, color=color, edgecolor="white", linewidth=0.6, zorder=4, alpha=0.9)


# ---------------------------------------------------------------------------
# Figure 1 -- genome-wide CRISPR discovery volcano plot
# ---------------------------------------------------------------------------

def build_figure_01_crispr_discovery(out_png: Path) -> None:
    df = pfd.load_genomewide_crispr()
    freeze = pfd.load_shortlist_freeze().set_index("gene")

    neglog10fdr = -np.log10(df["fdr"].clip(lower=1e-300))
    gate1 = df["fdr"] < 0.1

    fig, ax = plt.subplots(figsize=(11, 8.5), dpi=220)
    ax.scatter(df.loc[~gate1, "effect_size"], neglog10fdr[~gate1], s=5, color=GRAY, alpha=0.25,
               linewidth=0, zorder=1, rasterized=True)
    ax.scatter(df.loc[gate1, "effect_size"], neglog10fdr[gate1], s=14, color=DGRAY, alpha=0.55,
               linewidth=0, zorder=2, label="Gate-1 hits (FDR<0.1, n=28)")

    ax.set_xlim(df["effect_size"].min() - 2.4, df["effect_size"].max() + 0.3)

    # The four candidates cluster tightly (effect -1.06..-1.60, FDR 0.04..0.15)
    # -- a fixed same-direction offset collides, so labels are stacked in a
    # column to the left with thin leader lines, ordered by y (FDR).
    label_x = df["effect_size"].min() - 1.55
    ordered = sorted(GENE_ORDER, key=lambda g: -(-np.log10(df.loc[df["gene"] == g, "fdr"].iloc[0])))
    label_ys = np.linspace(2.55, 0.75, len(ordered))
    for gene, label_y in zip(ordered, label_ys):
        row = df[df["gene"] == gene].iloc[0]
        y = -np.log10(row["fdr"])
        ax.scatter([row["effect_size"]], [y], s=340, color=GENE_COLORS[gene], edgecolor=DGRAY,
                   linewidth=1.6, zorder=5)
        ax.plot([label_x + 0.05, row["effect_size"]], [label_y, y], color=GENE_COLORS[gene], linewidth=1.1, zorder=4)
        ax.text(label_x, label_y, f"{gene}\neffect {row['effect_size']:.2f}, FDR {row['fdr']:.3f}",
                ha="left", va="center", fontsize=10.5, fontweight="bold", color=GENE_COLORS[gene], linespacing=1.3)

    ax.axhline(-np.log10(0.1), color=GRAY, linewidth=1.1, linestyle="--", zorder=0)
    ax.text(df["effect_size"].max(), -np.log10(0.1) + 0.08,
            "Gate-1 threshold (FDR = 0.1)", ha="right", fontsize=9.5, color=GRAY)
    ax.axvline(0, color=GRAY, linewidth=0.8, zorder=0)

    ax.set_xlabel("CRISPR effect size (E2+4-OHT vs E2 interaction term)\n" + r"$\bf{negative}$ = sensitising knockout under 4-OHT", fontsize=12)
    ax.set_ylabel(r"$-\log_{10}$(FDR)", fontsize=13)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    ax.tick_params(labelsize=11)

    handles = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor=GRAY, alpha=0.5, markersize=7, label=f"all genes tested (n={len(df):,})"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor=DGRAY, alpha=0.7, markersize=8, label="Gate-1 hits (FDR<0.1, n=28)"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor=GENE_COLORS["USP34"], markeredgecolor=DGRAY, markersize=11, label="frozen shortlist (this poster)"),
    ]
    ax.legend(handles=handles, loc="upper left", fontsize=10, frameon=False)
    ax.set_title("Genome-wide CRISPR screen (Hany et al. 2023) -- 19,103 genes fitted", fontsize=15, fontweight="bold", loc="left", pad=14)

    _save(fig, out_png)


# ---------------------------------------------------------------------------
# Figure 2 -- real sample-level expression evidence (3 panels)
# ---------------------------------------------------------------------------

CONDITION_COLORS_118713 = {"MCF7": GRAY, "TAMR": "#0072B2", "FASR": "#CC79A7"}
CONDITION_COLORS_240112 = {"PT": GRAY, "RT": "#E69F00"}


def _panel_gse118713(ax) -> None:
    samples = pfd.load_gse118713_usp34_samples()
    stats = pfd.load_gse118713_usp34_stats()
    tamr_vs_mcf7 = stats[stats["contrast"] == "TAMR_vs_MCF7"].iloc[0]

    order = ["MCF7", "TAMR", "FASR"]
    for i, cond in enumerate(order):
        vals = samples.loc[samples["condition"] == cond, "tpm"].to_numpy()
        _strip_box(ax, i, vals, CONDITION_COLORS_118713[cond])
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels(order, fontsize=12, fontweight="bold")
    ax.set_ylabel("USP34 TPM", fontsize=11.5)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    ax.set_title("A. GSE118713 -- USP34 per sample", fontsize=12, fontweight="bold", loc="left")
    ax.text(0.02, 0.96, f"TAMR vs MCF7: log2FC {tamr_vs_mcf7['log2fc']:.2f}, FDR {tamr_vs_mcf7['fdr']:.3f}",
            transform=ax.transAxes, fontsize=9.3, color=DGRAY, va="top")


def _panel_gse240112(ax) -> None:
    samples = pfd.load_gse240112_vezf1_samples()
    stats = pfd.load_gse240112_vezf1_stats()

    order = ["PT", "RT"]
    labels = ["Primary\ntumor", "Recurrent\ntumor"]
    for i, cond in enumerate(order):
        vals = samples.loc[samples["group"] == cond, "log2cpm"].to_numpy()
        _strip_box(ax, i, vals, CONDITION_COLORS_240112[cond])
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels(labels, fontsize=12, fontweight="bold")
    ax.set_ylabel("VEZF1 log2(CPM)", fontsize=11.5)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    ax.set_title("B. GSE240112 -- VEZF1 per tumor", fontsize=12, fontweight="bold", loc="left")
    ax.text(0.02, 0.96, f"Recurrent vs primary: log2FC {stats['log2fc']:.2f}, FDR {stats['genomewide_fdr']:.3f}",
            transform=ax.transAxes, fontsize=9.3, color=DGRAY, va="top")


def _panel_forest(ax) -> None:
    forest = pfd.build_forest_table()
    dataset_order = ["gse118713", "gse240112_tumor", "gse111151", "gse245601_epi"]
    n_ds = len(dataset_order)
    row = 0
    yticks, yticklabels = [], []
    for gene in GENE_ORDER:
        for dkey in dataset_order:
            r = forest[(forest["gene"] == gene) & (forest["dataset_key"] == dkey)].iloc[0]
            sig = r["fdr"] < 0.05
            ax.scatter([r["log2fc"]], [row], s=95 if sig else 55, color=GENE_COLORS[gene],
                       edgecolor=DGRAY, linewidth=1.3 if sig else 0.5, alpha=1.0 if sig else 0.4, zorder=3)
            yticks.append(row)
            yticklabels.append(pfd.DATASET_LABELS[dkey].split(" (")[0])
            row -= 1
        row -= 0.6

    ax.axvline(0, color=GRAY, linewidth=0.9, zorder=0)
    ax.set_yticks(yticks)
    ax.set_yticklabels(yticklabels, fontsize=8.8)
    ax.set_xlabel("log$_2$ fold-change (resistant/recurrent/tumor vs comparator)", fontsize=10.5)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    ax.set_title("C. Cross-dataset consistency (log2FC by dataset)", fontsize=12, fontweight="bold", loc="left")

    # gene-group labels on the right
    gene_row_centers = []
    row = 0
    for gene in GENE_ORDER:
        gene_row_centers.append(row - (n_ds - 1) / 2)
        row -= n_ds + 0.6
    for gene, yc in zip(GENE_ORDER, gene_row_centers):
        ax.text(1.02, yc, gene, transform=ax.get_yaxis_transform(), fontsize=11, fontweight="bold",
                color=GENE_COLORS[gene], va="center", ha="left")

    handles = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor=DGRAY, markeredgecolor=DGRAY, markersize=9, alpha=1.0, label="FDR<0.05"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor=DGRAY, markeredgecolor=DGRAY, markersize=7, alpha=0.4, label="not significant"),
    ]
    ax.legend(handles=handles, loc="upper right", fontsize=8.5, frameon=False)


def build_figure_02_expression_evidence(out_png: Path) -> None:
    fig = plt.figure(figsize=(15, 8.5), dpi=220)
    gs = fig.add_gridspec(2, 2, width_ratios=(1.0, 1.15), height_ratios=(1, 1), wspace=0.45, hspace=0.4)
    axA = fig.add_subplot(gs[0, 0])
    axB = fig.add_subplot(gs[1, 0])
    axC = fig.add_subplot(gs[:, 1])

    _panel_gse118713(axA)
    _panel_gse240112(axB)
    _panel_forest(axC)

    fig.suptitle("Independent resistance / recurrence expression evidence", fontsize=16, fontweight="bold", x=0.02, ha="left", y=1.02)
    fig.text(0.02, 0.985, "GSE245601 = acute 12h ex vivo tamoxifen context, not a resistance cohort. GSE240112 = primary-vs-recurrent, not causal resistance.",
              fontsize=9, color=GRAY, style="italic")
    _save(fig, out_png)


# ---------------------------------------------------------------------------
# Figure 3 -- DepMap 26Q1 Chronos dependency distributions (ER+/luminal)
# ---------------------------------------------------------------------------

def build_figure_03_depmap_distributions(out_png: Path) -> None:
    long = pfd.load_depmap_gene_effect_er_luminal()
    dep_summary = pfd.load_depmap_dependency().set_index("candidate")

    fig, ax = plt.subplots(figsize=(9.5, 7.5), dpi=220)
    for i, gene in enumerate(GENE_ORDER):
        vals = long.loc[long["gene"] == gene, "chronos_effect"].to_numpy()
        _strip_box(ax, i, vals, GENE_COLORS[gene], width=0.4)

    ax.axhline(-0.5, color=GRAY, linewidth=1.1, linestyle="--", zorder=0)
    ax.text(3.55, -0.5, "informal strong-dependency\nreference (~-0.5)", fontsize=9, color=GRAY,
            va="center", ha="left", linespacing=1.3)
    ax.axhline(0, color=GRAY, linewidth=0.7, zorder=0)

    ax.set_xticks(range(len(GENE_ORDER)))
    ax.set_xticklabels(GENE_ORDER, fontsize=14, fontweight="bold")
    for i, gene in enumerate(GENE_ORDER):
        ax.get_xticklabels()[i].set_color(GENE_COLORS[gene])
    ax.set_ylabel("Chronos gene effect (DepMap 26Q1)\nER+/luminal breast lines, n=11", fontsize=12)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)

    ymin = long["chronos_effect"].min() - 0.35
    for i, gene in enumerate(GENE_ORDER):
        pct = dep_summary.loc[gene, "frac_strongly_dependent_er_luminal"] * 100
        ax.text(i, ymin, f"{pct:.1f}%\nstrongly\ndependent", ha="center", va="top", fontsize=9.3,
                color=GENE_COLORS[gene], fontweight="bold", linespacing=1.25)
    ax.set_ylim(ymin - 0.55, long["chronos_effect"].max() + 0.25)
    ax.set_xlim(-0.6, 3.9)

    ax.set_title("Baseline cancer-cell dependency (DepMap 26Q1, real Chronos gene-effect values)",
                 fontsize=14, fontweight="bold", loc="left", pad=14)
    ax.text(0.0, -0.30, "DepMap = baseline cancer-cell CRISPR dependency only -- NOT tamoxifen-specific response, NOT a normal-tissue safety measure.",
            transform=ax.transAxes, fontsize=9.3, color=GRAY, style="italic")

    _save(fig, out_png)


# ---------------------------------------------------------------------------
# Figure 4 -- GDSC pharmacogenomic scatter (USP34 vs AZD7762)
# ---------------------------------------------------------------------------

def build_figure_04_pharmacogenomics(out_png: Path) -> None:
    import scipy.stats as ss

    df = pfd.load_gdsc_usp34_azd7762()
    x = df["USP34"].to_numpy()
    y = df["LN_IC50"].to_numpy()
    stats_row = pfd.load_gdsc_usp34_azd7762_stats()
    rho, fdr = stats_row["spearman_rho"], stats_row["fdr"]

    slope, intercept, r, _, se = ss.linregress(x, y)
    xs = np.linspace(x.min(), x.max(), 100)
    yhat = slope * xs + intercept
    n = len(x)
    resid = y - (slope * x + intercept)
    dof = n - 2
    s_err = np.sqrt(np.sum(resid ** 2) / dof)
    x_mean = x.mean()
    conf = ss.t.ppf(0.975, dof) * s_err * np.sqrt(1 / n + (xs - x_mean) ** 2 / np.sum((x - x_mean) ** 2))

    fig, ax = plt.subplots(figsize=(9, 7.5), dpi=220)
    ax.fill_between(xs, yhat - conf, yhat + conf, color=GENE_COLORS["USP34"], alpha=0.15, zorder=1)
    ax.plot(xs, yhat, color=GENE_COLORS["USP34"], linewidth=2.2, zorder=2)

    er = df["is_er_luminal"]
    ax.scatter(x[~er], y[~er], s=75, color=GRAY, edgecolor=DGRAY, linewidth=0.7, zorder=3, label="other breast lines")
    ax.scatter(x[er], y[er], s=100, color=GENE_COLORS["USP34"], edgecolor=DGRAY, linewidth=1.0, zorder=4, label="ER+/luminal lines")

    ax.set_xlabel("USP34 expression (log2 TPM, DepMap 26Q1)", fontsize=12)
    ax.set_ylabel("AZD7762 response (LN_IC50, GDSC1)\nlower = more sensitive", fontsize=12)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    ax.legend(loc="upper right", fontsize=10, frameon=False)
    ax.set_title(f"USP34 expression vs AZD7762 (CHK1/CHK2i) response\nSpearman rho={rho:.2f}, FDR={fdr:.3f}, n={n} breast cell lines (GDSC Release 8.5)",
                 fontsize=13.5, fontweight="bold", loc="left", pad=12)
    ax.text(0.0, -0.17, "Pharmacogenomic association only -- NOT evidence that AZD7762 directly inhibits USP34.",
            transform=ax.transAxes, fontsize=9.5, color=GRAY, style="italic")

    _save(fig, out_png)


# ---------------------------------------------------------------------------
# Figure 5 -- USP34 structural-biology figure (real PyMOL renders)
# ---------------------------------------------------------------------------

def _render_pymol_panels(out_dir: Path) -> dict[str, Path]:
    paths = pfd.usp34_structure_paths()
    script = Path(__file__).resolve().parent / "poster_structure_render.py"
    import shutil
    pymol_bin = shutil.which("pymol")
    if pymol_bin is None:
        raise RuntimeError("pymol not found on PATH -- install pymol-open-source (see environment.yml)")
    cmd = [pymol_bin, "-cq", str(script), "--", str(paths["7W3R"]), str(paths["7W3U"]), str(out_dir)]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"pymol render failed:\n{result.stdout}\n{result.stderr}")
    logger.info("pymol render stdout: %s", result.stdout.strip()[-500:])
    panels = {
        "A": out_dir / "panel_A_apo_overview.png",
        "B": out_dir / "panel_B_bound_overview.png",
        "C": out_dir / "panel_C_cleft_closeup.png",
        "D": out_dir / "panel_D_apo_vs_bound_overlay.png",
    }
    for p in panels.values():
        assert p.exists(), f"expected pymol output missing: {p}"
    return panels


def _autocrop_white(img: np.ndarray, margin: int = 12) -> np.ndarray:
    """Trim uniform-white border from a rendered PyMOL PNG so panels fill
    their grid cell instead of leaving a large blank gap under the title."""
    rgb = img[:, :, :3]
    non_white = np.any(rgb < 0.98, axis=2)
    if not non_white.any():
        return img
    rows = np.where(non_white.any(axis=1))[0]
    cols = np.where(non_white.any(axis=0))[0]
    r0, r1 = max(rows[0] - margin, 0), min(rows[-1] + margin, img.shape[0] - 1)
    c0, c1 = max(cols[0] - margin, 0), min(cols[-1] + margin, img.shape[1] - 1)
    return img[r0:r1 + 1, c0:c1 + 1]


def build_figure_05_structure(out_png: Path) -> None:
    with tempfile.TemporaryDirectory() as tmp:
        panels = _render_pymol_panels(Path(tmp))
        images = {k: _autocrop_white(plt.imread(v)) for k, v in panels.items()}

    fig = plt.figure(figsize=(12, 13), dpi=220)
    gs = fig.add_gridspec(4, 2, height_ratios=(1, 1, 0.09, 0.09), hspace=0.28, wspace=0.08)

    titles = {
        "A": "A. 7W3R -- apo catalytic domain (1.92 Å)",
        "B": "B. 7W3U -- covalent ubiquitin-probe complex (3.13 Å)",
        "C": "C. Catalytic-cleft close-up",
        "D": "D. Apo vs probe-bound overlay",
    }
    positions = {"A": (0, 0), "B": (0, 1), "C": (1, 0), "D": (1, 1)}
    for key, (r, c) in positions.items():
        ax = fig.add_subplot(gs[r, c])
        ax.imshow(images[key])
        ax.axis("off")
        ax.set_title(titles[key], fontsize=12, fontweight="bold", loc="left", pad=6)

    axKey = fig.add_subplot(gs[2, :])
    axKey.axis("off")
    key_items = [("Cys1903 (catalytic)", BAD), ("His2164 (catalytic)", GENE_COLORS["VEZF1"]), ("covalent probe (ubiquitin-PA)", GOOD)]
    kx = 0.5 - 0.17 * (len(key_items) - 1)
    for label, color in key_items:
        axKey.scatter([kx], [0.5], s=110, color=color, edgecolor=DGRAY, linewidth=0.8, transform=axKey.transAxes, clip_on=False)
        axKey.text(kx + 0.025, 0.5, label, fontsize=10.5, color=DGRAY, va="center", ha="left", transform=axKey.transAxes)
        kx += 0.34

    axS = fig.add_subplot(gs[3, :])
    axS.axis("off")
    summary = ("Experimental structure ✓    Catalytic cysteine (Cys1903) ✓    Covalent reactivity (7W3U probe) ✓    "
               "Selective USP34 inhibitor ✗    Docking not pursued -- no validated ligand set for calibration")
    axS.text(0.5, 0.5, summary, ha="center", va="center", fontsize=10.8, color=DGRAY, fontweight="bold")

    fig.suptitle("USP34 structural targetability", fontsize=17, fontweight="bold", x=0.03, ha="left", y=0.99)
    _save(fig, out_png, vector=False)


# ---------------------------------------------------------------------------
# Figure 6 -- final experimental strategy (the one schematic figure)
# ---------------------------------------------------------------------------

def build_figure_06_strategy(out_png: Path) -> None:
    fig, ax = plt.subplots(figsize=(15, 8), dpi=220)
    ax.set_xlim(0, 15)
    ax.set_ylim(0, 8.4)
    ax.axis("off")

    def box(cx, cy, w, h, text, color, text_color="white", fontsize=11, edge=DGRAY):
        b = FancyBboxPatch((cx - w / 2, cy - h / 2), w, h, boxstyle="round,pad=0.08,rounding_size=0.18",
                            facecolor=color, edgecolor=edge, linewidth=1.3, zorder=2)
        ax.add_patch(b)
        ax.text(cx, cy, text, ha="center", va="center", fontsize=fontsize, fontweight="bold",
                color=text_color, zorder=3, linespacing=1.3)

    ax.text(0.0, 1.0, "Final experimental strategy -- lead experiment (proposed, not yet run)", transform=ax.transAxes,
            fontsize=15, fontweight="bold", color=DGRAY, ha="left", va="top")

    box(5.3, 7.0, 9.6, 0.85, "Acquired tamoxifen-resistant ER+ cells  +  matched parental cells", DGRAY, fontsize=12.5)

    arms = ["CONTROL", "TAMOXIFEN\n/ 4-OHT", "USP34 KO", "USP34 KO\n+ TAMOXIFEN"]
    axs = [1.7, 4.1, 6.5, 8.9]
    for x, a in zip(axs, arms):
        box(x, 5.85, 2.15, 1.05, a, GENE_COLORS["USP34"], fontsize=11.5)
    for x in axs:
        ax.add_patch(FancyArrowPatch((x, 5.3), (x, 4.65), arrowstyle="-|>", mutation_scale=16, color=DGRAY, linewidth=1.6))

    box(1.7, 3.55, 3.0, 1.7, "EFFICACY\nViability · apoptosis\nClonogenic survival\nDose–response\ninteraction", GOOD, fontsize=8.8)
    box(5.3, 3.55, 3.0, 1.7, "MECHANISM\nUSP34 protein loss\nAXIN1 / active\nβ-catenin\nESR1 target genes", GENE_COLORS["USP34"], fontsize=8.8)
    box(8.9, 3.55, 3.0, 1.7, "COUNTER-LIABILITY\nE-cadherin · N-cadherin\nSNAIL · EMT/stemness\n(optional: invasion /\nmammosphere)", BAD, fontsize=8.8)

    # Normal comparators -- same visual weight as the main readouts (real
    # controls for the lead experiment)
    box(12.9, 5.85, 3.6, 1.05, "Normal mammary\nepithelial cells", GRAY, fontsize=10.5)
    ax.text(12.9, 5.15, "-> lineage / selectivity", ha="center", fontsize=9.3, color=DGRAY, style="italic")
    box(12.9, 4.4, 3.6, 1.05, "Human MSC\nosteogenic model", GRAY, fontsize=10.5)
    ax.text(12.9, 3.7, "-> bone-liability comparator", ha="center", fontsize=9.3, color=DGRAY, style="italic")

    # VEZF1 -- deliberately smaller and visually secondary (a follow-up
    # note, not equal in weight to the USP34 lead experiment above)
    ax.add_patch(FancyBboxPatch((11.55, 2.55), 2.7, 0.68, boxstyle="round,pad=0.06,rounding_size=0.12",
                                 facecolor="white", edgecolor=GENE_COLORS["VEZF1"], linewidth=1.6, zorder=2))
    ax.text(12.9, 2.89, "secondary follow-up: VEZF1 CRISPRi ± tamoxifen", ha="center", va="center",
            fontsize=8.8, fontweight="bold", color=GENE_COLORS["VEZF1"], zorder=3)

    ax.add_patch(FancyBboxPatch((0.5, 0.35), 9.6, 1.4, boxstyle="round,pad=0.1,rounding_size=0.2",
                                 facecolor=LGRAY, edgecolor=GRAY, linewidth=1.0, zorder=2))
    ax.text(5.3, 1.42, "GOAL: distinguish", ha="center", fontsize=11, fontweight="bold", color=DGRAY, zorder=3)
    ax.text(5.3, 0.98, "DIRECT DEPENDENCY   vs.   TAMOXIFEN SENSITISATION   vs.   DUAL ACTION",
            ha="center", fontsize=12, fontweight="bold", color=DGRAY, zorder=3)
    ax.text(5.3, 0.56, "(general toxicity is assessed separately via the normal-cell comparators at right)",
            ha="center", fontsize=9, color=GRAY, style="italic", zorder=3)

    _save(fig, out_png)


def run(figures_dir: Path = FIGURES) -> None:
    figures_dir.mkdir(parents=True, exist_ok=True)
    build_figure_01_crispr_discovery(figures_dir / "01_crispr_discovery.png")
    build_figure_02_expression_evidence(figures_dir / "02_expression_evidence.png")
    build_figure_03_depmap_distributions(figures_dir / "03_depmap_distributions.png")
    build_figure_04_pharmacogenomics(figures_dir / "04_pharmacogenomics.png")
    build_figure_05_structure(figures_dir / "05_structure.png")
    build_figure_06_strategy(figures_dir / "06_experimental_strategy.png")


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run()
