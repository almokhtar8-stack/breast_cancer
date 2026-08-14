"""Independent-validation Part 13 figures. Reads only the already-built
tables under results/tables/independent_validation/ (and, for figure 1,
the underlying per-sample expression it needs that isn't pre-aggregated
in any table) -- no new statistics are computed here beyond what the
analysis modules already produced. Palette follows the CVD-safe
categorical set already used by this project's other figures
(systems_network_visualization_9_10.py): blue #0072B2, orange #E69F00,
green #009E73, red #b0392f, neutral grays.
"""

from __future__ import annotations

import logging
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D

from src.independent_validation_depmap_data import load_config as load_depmap_config
from src.independent_validation_depmap_data import load_gene_effect, load_model
from src.independent_validation_tcga_data import build_cohort_table
from src.independent_validation_tcga_data import load_config as load_tcga_config
from src.independent_validation_tcga_data import load_expression

logger = logging.getLogger(__name__)

CANDIDATES = ["USP34", "VEZF1", "EML5", "CITED2"]
TABLES = Path("results/tables/independent_validation")
FIGURES = Path("results/figures/independent_validation")

BLUE, ORANGE, GREEN, RED, GRAY, DGRAY = "#0072B2", "#E69F00", "#009E73", "#b0392f", "#9a9a9a", "#555555"
DIVERGING = LinearSegmentedColormap.from_list("bwr_diverging", ["#0072B2", "#f2f2f2", "#b0392f"])


def build_figure_01(out_fig: Path) -> None:
    cfg = load_tcga_config()
    cohort = build_cohort_table(cfg)
    expr = load_expression(cfg, genes=CANDIDATES)
    df = cohort.join(expr)
    primary = df.loc[df["is_primary_tumor"]]

    fig, axes = plt.subplots(2, 2, figsize=(11, 9), dpi=200)
    fig.patch.set_facecolor("white")
    group_defs = [
        ("ER-", primary.loc[primary["ER_STATUS"] == "Negative"], GRAY),
        ("ER+", primary.loc[primary["ER_STATUS"] == "Positive"], BLUE),
        ("Luminal A", primary.loc[primary["PAM50_SUBTYPE"] == "Luminal A"], GREEN),
        ("Luminal B", primary.loc[primary["PAM50_SUBTYPE"] == "Luminal B"], ORANGE),
        ("Normal", df.loc[df["is_normal"]], DGRAY),
        ("Tumor", primary, RED),
    ]
    for ax, candidate in zip(axes.flat, CANDIDATES):
        data = [g[candidate].dropna().values for _, g, _ in group_defs]
        labels = [f"{name}\n(n={len(v)})" for (name, _, _), v in zip(group_defs, data)]
        colors = [c for _, _, c in group_defs]
        bp = ax.boxplot(data, patch_artist=True, showfliers=False, widths=0.6)
        for patch, color in zip(bp["boxes"], colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.55)
            patch.set_edgecolor(DGRAY)
        for median in bp["medians"]:
            median.set_color("black")
        ax.set_xticks(range(1, len(labels) + 1))
        ax.set_xticklabels(labels, fontsize=8)
        ax.set_ylabel("log2(TPM+1)")
        ax.set_title(candidate, fontsize=12, fontweight="bold")
        ax.spines[["top", "right"]].set_visible(False)
    fig.suptitle("TCGA-BRCA candidate expression across clinical/subtype groups", fontsize=13, fontweight="bold", y=0.995)
    fig.text(0.5, 0.955, "ER status = clinical IHC (primary); PAM50 = molecular-subtype proxy (secondary); Normal/Tumor = unpaired, see TCGA_candidate_expression.tsv for the paired test", ha="center", fontsize=8, style="italic", color=DGRAY)
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    out_fig.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_fig, facecolor="white")
    plt.close(fig)
    logger.info("wrote %s", out_fig)


def build_figure_02(out_fig: Path) -> None:
    df = pd.read_csv(TABLES / "TCGA_candidate_pathway_associations.tsv", sep="\t")
    df = df.loc[df["pathway"] != "NONE"]
    mat = df.pivot_table(index="candidate", columns="pathway", values="spearman_rho")
    fdr = df.pivot_table(index="candidate", columns="pathway", values="fdr")
    mat = mat.reindex(CANDIDATES)
    fdr = fdr.reindex(CANDIDATES)

    fig, ax = plt.subplots(figsize=(max(8, 0.9 * mat.shape[1] + 3), 4.2), dpi=200)
    fig.patch.set_facecolor("white")
    im = ax.imshow(mat.values, cmap=DIVERGING, vmin=-0.3, vmax=0.3, aspect="auto")
    ax.set_xticks(range(mat.shape[1]))
    ax.set_xticklabels([c.replace("HALLMARK_", "").replace("GOBP_", "").replace("_", " ") for c in mat.columns], rotation=40, ha="right", fontsize=8)
    ax.set_yticks(range(mat.shape[0]))
    ax.set_yticklabels(mat.index, fontsize=10, fontweight="bold")
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            v = mat.values[i, j]
            f = fdr.values[i, j]
            if np.isnan(v):
                ax.text(j, i, "n/a", ha="center", va="center", fontsize=7, color=DGRAY)
                continue
            star = "*" if f < 0.05 else ""
            ax.text(j, i, f"{v:.2f}{star}", ha="center", va="center", fontsize=8, color="black")
    cbar = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02)
    cbar.ax.set_title("Spearman\nrho", fontsize=8)
    ax.set_title("TCGA-BRCA (ER+ tumors): candidate expression vs pathway ssGSEA score\n(* = FDR<0.05; EML5 has no candidate-specific pathway declared)", fontsize=10.5, fontweight="bold")
    fig.tight_layout()
    out_fig.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_fig, facecolor="white")
    plt.close(fig)
    logger.info("wrote %s", out_fig)


def build_figure_03(out_fig: Path, release: str | None = None) -> None:
    cfg = load_depmap_config()
    release = release or cfg["independent_validation"]["depmap"]["active_release"]
    model = load_model(cfg, release)
    effect = load_gene_effect(cfg, release, CANDIDATES)
    breast_ids = model.index[model["is_breast"]]
    luminal_ids = model.index[model["is_er_luminal"]]

    fig, axes = plt.subplots(1, 4, figsize=(14, 4.2), dpi=200, sharey=True)
    fig.patch.set_facecolor("white")
    group_defs = [("All cancer\nlines", effect.index, GRAY), ("Breast\nlines", breast_ids, ORANGE), ("ER+/luminal\nbreast lines", luminal_ids, BLUE)]
    for ax, candidate in zip(axes, CANDIDATES):
        data = [effect.loc[effect.index.isin(ids), candidate].dropna().values for _, ids, _ in group_defs]
        colors = [c for _, _, c in group_defs]
        bp = ax.boxplot(data, patch_artist=True, showfliers=True, widths=0.6, flierprops=dict(marker="o", markersize=2, alpha=0.4))
        for patch, color in zip(bp["boxes"], colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.55)
            patch.set_edgecolor(DGRAY)
        ax.axhline(0, color="black", lw=0.8, ls=":")
        ax.axhline(-0.5, color=RED, lw=0.8, ls="--", alpha=0.6)
        ax.set_xticks([1, 2, 3])
        ax.set_xticklabels([f"{n}\n(n={len(v)})" for (n, _, _), v in zip(group_defs, data)], fontsize=7.5)
        ax.set_title(candidate, fontsize=12, fontweight="bold")
        ax.spines[["top", "right"]].set_visible(False)
    axes[0].set_ylabel(f"DepMap {release} Chronos gene effect\n(more negative = greater dependency)")
    fig.legend(handles=[Line2D([0], [0], color=RED, lw=1.2, ls="--", label="common informal 'dependent' reference line (-0.5)")], loc="lower center", bbox_to_anchor=(0.5, -0.02), fontsize=8, frameon=False)
    fig.suptitle(f"DepMap Public {release} baseline CRISPR dependency (NOT a tamoxifen-context screen)", fontsize=12.5, fontweight="bold")
    fig.tight_layout(rect=(0, 0.04, 1, 0.93))
    out_fig.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_fig, facecolor="white")
    plt.close(fig)
    logger.info("wrote %s", out_fig)


def build_figure_04(out_fig: Path) -> None:
    integ = pd.read_csv(TABLES / "four_candidate_independent_validation.tsv", sep="\t").set_index("candidate").loc[CANDIDATES]
    depmap_release = integ["depmap_release"].iloc[0]

    cols = [
        ("Hany CRISPR\n(drug-context)", lambda r: f"sensitising_KO\n(FDR={r.frozen_crispr_fdr:.3f})", lambda r: GREEN if r.frozen_crispr_fdr < 0.05 else ORANGE),
        ("TCGA ER+/ER-\nexpression", lambda r: "significant" if "significant)" in r.tcga_er_pos_vs_neg_result and "not significant" not in r.tcga_er_pos_vs_neg_result else "n.s.", lambda r: GREEN if "significant)" in r.tcga_er_pos_vs_neg_result and "not significant" not in r.tcga_er_pos_vs_neg_result else GRAY),
        ("TCGA pathway\nassociation", lambda r: "significant" if "(significant)" in r.tcga_strongest_pathway_association else "n.s./none", lambda r: GREEN if "(significant)" in r.tcga_strongest_pathway_association else GRAY),
        ("TCGA clinical\n(OS, ER+ adj.)", lambda r: "significant" if r.integration_clinically_relevant else "n.s.", lambda r: GREEN if r.integration_clinically_relevant else GRAY),
        (f"DepMap {depmap_release}\nbaseline dependency", lambda r: r.depmap_essentiality_concern.split("_", 1)[1].replace("_", " ").title(), lambda r: {"A": RED, "B": ORANGE, "C": ORANGE, "D": GREEN, "E": GRAY}[r.depmap_essentiality_concern[0]]),
        ("Independent-validation\nstrength", lambda r: r.integration_validation_strength.split("_", 1)[1].replace("_", " ").title(), lambda r: {"1": GREEN, "2": BLUE, "3": ORANGE, "4": GRAY, "5": RED}[r.integration_validation_strength[0]]),
    ]

    col_widths = [1.3, 1.3, 1.3, 1.3, 1.6, 1.8]
    col_edges = [0.0]
    for w in col_widths:
        col_edges.append(col_edges[-1] + w)
    n_cols_total = col_edges[-1]

    fig, ax = plt.subplots(figsize=(14.5, 1.4 + 1.15 * len(CANDIDATES)), dpi=200)
    fig.patch.set_facecolor("white")
    n_rows = len(CANDIDATES)
    for i, candidate in enumerate(CANDIDATES):
        r = integ.loc[candidate]
        for j, (_, valuefn, colorfn) in enumerate(cols):
            val = valuefn(r)
            color = colorfn(r)
            x0, w = col_edges[j], col_widths[j]
            ax.add_patch(plt.Rectangle((x0, n_rows - 1 - i), w, 1, facecolor=color, alpha=0.35, edgecolor="white"))
            ax.text(x0 + w / 2, n_rows - 1 - i + 0.5, val, ha="center", va="center", fontsize=8.3, wrap=True)
        ax.text(-0.15, n_rows - 1 - i + 0.5, candidate, ha="right", va="center", fontsize=11, fontweight="bold")
    for j, (label, _, _) in enumerate(cols):
        x0, w = col_edges[j], col_widths[j]
        ax.text(x0 + w / 2, n_rows + 0.15, label, ha="center", va="bottom", fontsize=8.3, fontweight="bold")
    ax.set_xlim(-2.6, n_cols_total)
    ax.set_ylim(0, n_rows + 1.0)
    ax.axis("off")
    fig.suptitle(f"Integrated candidate validation: project function + TCGA human relevance + DepMap Public {depmap_release} baseline dependency", fontsize=12, fontweight="bold", x=0.5, y=0.995)
    fig.text(0.5, 0.02, "Frozen therapeutic ranking (USP34 > VEZF1 > EML5 > CITED2) is unchanged by this figure -- see report for the separate follow-up rankings.", ha="center", fontsize=8, style="italic", color=DGRAY)
    fig.tight_layout(rect=(0.02, 0.04, 1, 0.90))
    out_fig.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_fig, facecolor="white")
    plt.close(fig)
    logger.info("wrote %s", out_fig)


def build_figure_05(out_fig: Path) -> None:
    df = pd.read_csv(TABLES / "TCGA_candidate_clinical.tsv", sep="\t")
    df = df.loc[df["model"] == "adjusted_age_stage"]
    df = df.set_index(["candidate", "cohort"])

    fig, ax = plt.subplots(figsize=(8, 5), dpi=200)
    fig.patch.set_facecolor("white")
    y_labels, y_pos = [], []
    y = 0
    for candidate in CANDIDATES:
        for cohort, color in [("all_primary_tumors", GRAY), ("ER_positive", BLUE)]:
            row = df.loc[(candidate, cohort)]
            if pd.isna(row.hr_per_sd):
                y += 1
                continue
            sig = row.fdr < 0.05
            ax.errorbar(row.hr_per_sd, y, xerr=[[row.hr_per_sd - row.ci_low], [row.ci_high - row.hr_per_sd]], fmt="o", color=color if sig else GRAY, markersize=7 if sig else 5, capsize=3, elinewidth=1.5 if sig else 1)
            y_labels.append(f"{candidate} ({'all' if cohort=='all_primary_tumors' else 'ER+'}){'  *FDR<0.05' if sig else ''}")
            y_pos.append(y)
            y += 1
    ax.axvline(1.0, color="black", lw=0.8, ls=":")
    ax.set_yticks(y_pos)
    ax.set_yticklabels(y_labels, fontsize=8.5)
    ax.set_xlabel("HR per SD of expression (age+stage-adjusted Cox, overall survival)")
    ax.invert_yaxis()
    ax.spines[["top", "right"]].set_visible(False)
    fig.suptitle("TCGA-BRCA overall-survival association (NOT a tamoxifen-response result)", fontsize=11.5, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    out_fig.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_fig, facecolor="white")
    plt.close(fig)
    logger.info("wrote %s", out_fig)


def build_figure_06(out_fig: Path) -> None:
    df = pd.read_csv(TABLES / "DepMap_24Q4_vs_26Q1_comparison.tsv", sep="\t").set_index("candidate").loc[CANDIDATES]

    panels = [
        ("All cancer lines\n(median gene effect)", "24Q4_all_median", "26Q1_all_median", False),
        ("Breast lines\n(median gene effect)", "24Q4_breast_median", "26Q1_breast_median", False),
        ("ER+/luminal breast lines\n(median gene effect)", "24Q4_ERluminal_median", "26Q1_ERluminal_median", False),
        ("ER+/luminal breast lines\n(% strongly dependent, p>0.5)", "24Q4_strong_dependency_fraction_ERluminal", "26Q1_strong_dependency_fraction_ERluminal", True),
    ]
    fig, axes = plt.subplots(1, 4, figsize=(16.5, 4.4), dpi=200, sharey=True)
    fig.patch.set_facecolor("white")
    y_pos = list(range(len(CANDIDATES)))
    for ax, (title, old_col, new_col, is_pct) in zip(axes, panels):
        for y, candidate in zip(y_pos, CANDIDATES):
            o, n = df.loc[candidate, old_col], df.loc[candidate, new_col]
            if is_pct:
                o, n = o * 100, n * 100
            ax.plot([o, n], [y, y], color=DGRAY, lw=1.3, zorder=1)
            ax.scatter([o], [y], color=GRAY, s=70, zorder=2, label="24Q4" if y == 0 else None)
            ax.scatter([n], [y], color=BLUE, s=70, zorder=2, label="26Q1" if y == 0 else None)
        ax.axvline(0, color="black", lw=0.7, ls=":")
        ax.set_title(title, fontsize=9.5, fontweight="bold")
        ax.spines[["top", "right"]].set_visible(False)
    axes[0].set_yticks(y_pos)
    axes[0].set_yticklabels(CANDIDATES, fontsize=10.5, fontweight="bold")
    axes[0].invert_yaxis()
    axes[1].set_xlabel("Median Chronos gene effect (more negative = greater dependency)")
    axes[3].set_xlabel("% strongly dependent")
    fig.legend(handles=[Line2D([0], [0], marker="o", color="w", markerfacecolor=GRAY, markersize=8, label="24Q4"), Line2D([0], [0], marker="o", color="w", markerfacecolor=BLUE, markersize=8, label="26Q1")], loc="lower center", bbox_to_anchor=(0.5, -0.02), ncol=2, fontsize=9, frameon=False)
    fig.suptitle("DepMap 24Q4 -> 26Q1: baseline dependency comparison", fontsize=13, fontweight="bold", x=0.5, y=0.99)
    fig.text(0.5, 0.925, "Both continuous gene-effect and probability-based strong-dependency % are now available for 26Q1 (CRISPRGeneDependency.csv obtained)", ha="center", fontsize=8.5, style="italic", color=DGRAY)
    fig.tight_layout(rect=(0, 0.06, 1, 0.88))
    out_fig.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_fig, facecolor="white")
    plt.close(fig)
    logger.info("wrote %s", out_fig)


def run() -> None:
    build_figure_01(FIGURES / "01_TCGA_four_candidate_expression.png")
    build_figure_02(FIGURES / "02_TCGA_candidate_pathway_associations.png")
    build_figure_03(FIGURES / "03_DepMap_four_candidate_dependency.png")
    build_figure_04(FIGURES / "04_integrated_candidate_validation.png")
    build_figure_05(FIGURES / "05_TCGA_candidate_survival.png")
    build_figure_06(FIGURES / "06_DepMap_24Q4_vs_26Q1_comparison.png")


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run()
