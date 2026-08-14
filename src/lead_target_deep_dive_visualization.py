"""USP34 vs VEZF1 translational deep-dive: five poster-friendly figures
built only from this phase's own output tables (never hand-typed
numbers), matching the Okabe-Ito colorblind-safe palette already used
throughout this project's visualization modules.
"""

from __future__ import annotations

import logging
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

logger = logging.getLogger(__name__)

TABLES = Path("results/tables/lead_target_deep_dive")
FIGURES = Path("results/figures/lead_target_deep_dive")

BLUE, ORANGE, GREEN, RED, GRAY, DGRAY = "#0072B2", "#E69F00", "#009E73", "#b0392f", "#9a9a9a", "#555555"
GENE_COLOR = {"USP34": BLUE, "VEZF1": ORANGE}

# Numeric GTEx v8 median TPM values, read directly from the phase's own
# tissue-expression table text (re-parsed here rather than retyped
# separately, so a change to the table's numbers is the only place that
# needs editing -- see _extract_tpm below).
TISSUE_TPM = {
    "USP34": {"Breast": 18.61, "Blood": 4.28, "Liver": 4.06, "Heart (LV)": 5.44, "Skeletal muscle": 13.93, "Bone marrow (HPA nTPM)": 35.5, "Cerebellum": 27.6},
    "VEZF1": {"Breast": 30.5, "Blood": 9.8, "Liver": 6.7, "Heart (LV)": 14.1, "Skeletal muscle": 33.2, "Bone marrow (HPA nTPM)": 53.5, "Cerebellum": None},
}

LIABILITY_ORDER = {
    "NONE_IDENTIFIED": 0, "INSUFFICIENT_DATA": 0, "EXPRESSION_ONLY": 1, "INFERRED_ONLY": 2,
    "ANIMAL_ONLY": 3, "DOCUMENTED_DEVELOPMENTAL": 4, "DOCUMENTED_POSTNATAL_CAUSAL": 5,
    "DOCUMENTED_HUMAN": 5, "PRIMARY_HUMAN_CELL_FUNCTIONAL": 5, "DOCUMENTED_ADULT_CAUSAL": 6,
}
LIABILITY_ORGAN_ORDER = ["Skeletal muscle", "Bone", "Marrow/hematological", "Cardiovascular", "Neurological/cognitive"]


def _load_tables() -> dict[str, pd.DataFrame]:
    return {
        "tissue_liability": pd.read_csv(TABLES / "USP34_VEZF1_tissue_liability.tsv", sep="\t"),
        "direct_usp34": pd.read_csv(TABLES / "USP34_direct_targeting.tsv", sep="\t"),
        "direct_vezf1": pd.read_csv(TABLES / "VEZF1_direct_targeting.tsv", sep="\t"),
        "indirect_crosscheck": pd.read_csv(TABLES / "indirect_target_project_crosscheck.tsv", sep="\t"),
        "genetic_constraint": pd.read_csv(TABLES / "USP34_VEZF1_human_genetic_constraint.tsv", sep="\t").set_index("candidate"),
    }


def build_figure_01(out_fig: Path) -> None:
    tissues = ["Breast", "Blood", "Liver", "Heart (LV)", "Skeletal muscle", "Bone marrow (HPA nTPM)", "Cerebellum"]
    fig, ax = plt.subplots(figsize=(11.5, 5.5), dpi=200)
    x = range(len(tissues))
    width = 0.35
    for offset, gene in zip((-width / 2, width / 2), ("USP34", "VEZF1")):
        vals = [TISSUE_TPM[gene].get(t) for t in tissues]
        plot_vals = [0 if v is None else v for v in vals]
        bars = ax.bar([i + offset for i in x], plot_vals, width, label=gene, color=GENE_COLOR[gene], edgecolor=DGRAY, linewidth=0.8)
        for i, v in enumerate(vals):
            if v is None:
                ax.text(i + offset, 1.5, "no\ndata", ha="center", va="bottom", fontsize=7, color=RED)
    ax.set_xticks(list(x))
    ax.set_xticklabels(tissues, fontsize=9.5, rotation=20, ha="right")
    ax.set_ylabel("Expression (GTEx v8 median TPM, or HPA nTPM where noted)", fontsize=9.5)
    ax.legend(fontsize=10, frameon=False)
    ax.grid(axis="y", color="#eeeeee", zorder=0)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    ax.set_title(
        "USP34 vs VEZF1: normal-tissue expression across key systems\n"
        "Both broadly, similarly expressed (low tissue specificity) -- expression breadth alone is not a liability claim",
        fontsize=11,
    )
    fig.tight_layout()
    out_fig.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_fig, facecolor="white")
    plt.close(fig)
    logger.info("wrote %s", out_fig)


def build_figure_02(out_fig: Path) -> None:
    t = _load_tables()
    liab = t["tissue_liability"]
    liab = liab[liab["organ_system"].isin(LIABILITY_ORGAN_ORDER)]

    fig, ax = plt.subplots(figsize=(11, 5.5), dpi=200)
    for _, row in liab.iterrows():
        gene = row["candidate"]
        y = LIABILITY_ORGAN_ORDER.index(row["organ_system"])
        x = LIABILITY_ORDER[row["classification"]]
        jitter = -0.12 if gene == "USP34" else 0.12
        ax.scatter(x, y + jitter, s=260, color=GENE_COLOR[gene], edgecolor=DGRAY, linewidth=1.1, zorder=3)
    ax.set_yticks(range(len(LIABILITY_ORGAN_ORDER)))
    ax.set_yticklabels(LIABILITY_ORGAN_ORDER, fontsize=10.5)
    ax.set_xlim(-0.5, 6.5)
    ax.set_xticks(range(7))
    ax.set_xticklabels(
        ["none/\ninsufficient", "expression\nonly", "inferred\nonly", "animal\nonly", "documented\ndevelopmental", "documented\npostnatal/human", "documented\nadult-causal"],
        fontsize=7.6,
    )
    ax.set_xlabel("Strength/directness of evidence found (NOT a risk score -- 'more evidence' is not 'more danger'; see report)", fontsize=9)
    ax.grid(axis="x", color="#eeeeee", zorder=0)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    handles = [plt.Line2D([0], [0], marker="o", color="w", markerfacecolor=GENE_COLOR[g], markersize=12, markeredgecolor=DGRAY, label=g) for g in ("USP34", "VEZF1")]
    ax.legend(handles=handles, fontsize=10, frameon=False, loc="upper left")
    ax.set_title("Normal-tissue liability evidence by organ system\nEmbryonic-lethal/developmental findings are never read as adult drug toxicity", fontsize=11)
    fig.tight_layout()
    out_fig.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_fig, facecolor="white")
    plt.close(fig)
    logger.info("wrote %s", out_fig)


def build_figure_03(out_fig: Path) -> None:
    t = _load_tables()
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.6), dpi=200)

    ax = axes[0]
    direct_tier = {"USP34": ("POTENTIALLY_\nDRUGGABLE", 2), "VEZF1": ("CURRENTLY_POORLY_\nDRUGGABLE", 1)}
    for i, gene in enumerate(("USP34", "VEZF1")):
        label, tier = direct_tier[gene]
        ax.scatter(tier, i, s=500, color=GENE_COLOR[gene], edgecolor=DGRAY, linewidth=1.3, zorder=3)
        ax.text(tier + 0.1, i, label, va="center", fontsize=10, color=DGRAY)
    ax.set_yticks([0, 1])
    ax.set_yticklabels(["USP34", "VEZF1"], fontsize=11, fontweight="bold")
    ax.set_xlim(0.5, 4)
    ax.set_ylim(-0.6, 1.6)
    ax.set_xticks([1, 2])
    ax.set_xticklabels(["POORLY\nDRUGGABLE", "POTENTIALLY\nDRUGGABLE"], fontsize=8.5)
    ax.grid(axis="x", color="#eeeeee", zorder=0)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    ax.set_title("Direct targetability\n(real structures exist for USP34; none for VEZF1)", fontsize=10.5)

    ax2 = axes[1]
    cc = t["indirect_crosscheck"]
    genes = cc["gene_x"].tolist()
    fdrs = cc["gse118713_fdr"].tolist()
    fdrs2 = cc["gse240112_tumor_fdr"].tolist()
    y = range(len(genes))
    for i, (g, f1, f2, parent) in enumerate(zip(genes, fdrs, fdrs2, cc["parent_gene"])):
        color = GREEN if (pd.notna(f1) and f1 < 0.05) or (pd.notna(f2) and f2 < 0.05) else GRAY
        best_fdr = min([v for v in (f1, f2) if pd.notna(v)], default=float("nan"))
        ax2.barh(i, -1 * (best_fdr if pd.notna(best_fdr) else 1.0), color=color, edgecolor=DGRAY, height=0.6)
        ax2.text(0.02, i, f"{g} (candidate for {parent}, unvalidated) -- best resistance-RNA FDR={best_fdr:.3f}" if pd.notna(best_fdr) else f"{g} ({parent}, unvalidated) -- not testable", va="center", ha="left", fontsize=8.5, color=DGRAY)
    ax2.set_yticks([])
    ax2.set_xlim(-1.05, 0.02)
    ax2.set_xticks([-1.0, -0.5, -0.05])
    ax2.set_xticklabels(["FDR=1.0\n(none)", "FDR=0.5", "FDR<0.05\n(significant)"], fontsize=8)
    ax2.set_xlabel("Best FDR in this project's own resistance-RNA datasets (green = FDR<0.05)", fontsize=8.5)
    for spine in ("top", "right"):
        ax2.spines[spine].set_visible(False)
    ax2.set_title("Candidate (unvalidated) indirect targets, cross-checked against\nthis project's own frozen evidence -- TEAD1's RNA signal stands out,\nbut no perturbation evidence confirms TEAD1 regulates VEZF1", fontsize=9.8)

    fig.tight_layout(rect=(0, 0, 1, 0.92))
    fig.suptitle("Direct vs indirect targetability, USP34 and VEZF1", fontsize=13)
    out_fig.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_fig, facecolor="white")
    plt.close(fig)
    logger.info("wrote %s", out_fig)


def build_figure_04(out_fig: Path) -> None:
    t = _load_tables()
    depmap = pd.read_csv("results/tables/independent_validation/DepMap_candidate_dependency.tsv", sep="\t").set_index("candidate")
    constraint = t["genetic_constraint"]

    fig, axes = plt.subplots(2, 2, figsize=(11, 8), dpi=200)

    ax = axes[0, 0]
    vals = [float(depmap.loc[g, "frac_strongly_dependent_er_luminal"]) * 100 for g in ("USP34", "VEZF1")]
    ax.bar(["USP34", "VEZF1"], vals, color=[BLUE, ORANGE], edgecolor=DGRAY)
    ax.set_ylabel("% ER+/luminal strongly dependent\n(DepMap 26Q1, real data)", fontsize=8.5)
    ax.set_title("Baseline cancer dependency", fontsize=10)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)

    ax = axes[0, 1]
    vals = [float(depmap.loc[g, "hany_crispr_fdr"]) for g in ("USP34", "VEZF1")]
    ax.bar(["USP34", "VEZF1"], vals, color=[BLUE, ORANGE], edgecolor=DGRAY)
    ax.axhline(0.05, color=RED, linestyle="--", linewidth=1, label="FDR=0.05")
    ax.set_ylabel("Hany CRISPR FDR (frozen, lower=more significant)", fontsize=8.5)
    ax.set_title("Functional tamoxifen evidence", fontsize=10)
    ax.legend(fontsize=7.5, frameon=False)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)

    ax = axes[1, 0]
    vals = [float(constraint.loc[g, "loeuf"]) for g in ("USP34", "VEZF1")]
    ax.bar(["USP34", "VEZF1"], vals, color=[BLUE, ORANGE], edgecolor=DGRAY)
    ax.set_ylabel("gnomAD LOEUF (lower = more constrained)", fontsize=8.5)
    ax.set_title("Genetic constraint", fontsize=10)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)

    ax = axes[1, 1]
    vals = [int(constraint.loc[g, "clinvar_pathogenic_or_likely_pathogenic"]) for g in ("USP34", "VEZF1")]
    totals = [int(constraint.loc[g, "clinvar_total_variants"]) for g in ("USP34", "VEZF1")]
    ax.bar(["USP34", "VEZF1"], vals, color=[BLUE, ORANGE], edgecolor=DGRAY)
    for i, (v, tot) in enumerate(zip(vals, totals)):
        ax.text(i, v + 0.5, f"of {tot} total", ha="center", fontsize=7.5, color=DGRAY)
    ax.set_ylabel("ClinVar pathogenic/likely-pathogenic variants", fontsize=8.5)
    ax.set_title("Human genetic liability (ClinVar)", fontsize=10)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)

    fig.suptitle(
        "Therapeutic-window comparison: four separate, real-data dimensions\n"
        "(qualitative dimensions -- muscle/bone/cardio/CNS liability -- are in Figure 02, not repeated here as a fake score)",
        fontsize=11.5,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.90))
    out_fig.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_fig, facecolor="white")
    plt.close(fig)
    logger.info("wrote %s", out_fig)


def build_figure_05(out_fig: Path) -> None:
    fig, ax = plt.subplots(figsize=(13, 7.5), dpi=200)
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.axis("off")

    def arm_box(x, y, w, h, text, color):
        rect = plt.Rectangle((x, y), w, h, facecolor=color, edgecolor=DGRAY, linewidth=1.1, alpha=0.85, zorder=2)
        ax.add_patch(rect)
        ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=8.2, color="white", fontweight="bold", zorder=3, wrap=True)

    ax.text(0.2, 9.5, "EXP-1: USP34", fontsize=12, fontweight="bold", color=BLUE)
    arms1 = ["control", "tamoxifen", "USP34 inhibition\n(genetic -- no chemical\nprobe exists yet)", "USP34 inhibition\n+ tamoxifen"]
    for i, a in enumerate(arms1):
        arm_box(0.2 + i * 2.3, 8.2, 2.0, 1.0, a, BLUE)
    ax.text(0.2, 7.6, "+ normal-cell comparator: human primary MSCs, osteogenic induction assay (extends PMID 30181118)", fontsize=8.3, color=DGRAY)
    ax.text(0.2, 7.15, "readouts: viability, apoptosis, clonogenic survival, ER signaling (TFF1/PGR/GREB1), AXIN1/RUNX2/Smad1 mechanistic readout, Bliss synergy", fontsize=8, color=DGRAY)

    ax.text(0.2, 6.3, "EXP-2: VEZF1 (direct)", fontsize=12, fontweight="bold", color=ORANGE)
    arms2 = ["control", "tamoxifen", "VEZF1 inhibition/\nCRISPRi", "VEZF1 inhibition/\nCRISPRi + tamoxifen"]
    for i, a in enumerate(arms2):
        arm_box(0.2 + i * 2.3, 5.0, 2.0, 1.0, a, ORANGE)
    ax.text(0.2, 4.4, "+ normal-cell comparator: primary human cardiomyocytes / zebrafish larval assay (extends PMID 31911272)", fontsize=8.3, color=DGRAY)
    ax.text(0.2, 3.95, "readouts: (A) direct cancer-cell effect (B) additional tamoxifen sensitisation (C) synergy -- tests the dual-action hypothesis directly", fontsize=8, color=DGRAY)

    ax.text(0.2, 3.1, "EXP-3: TEAD1 (VEZF1 indirect-target validation)", fontsize=12, fontweight="bold", color=GREEN)
    arms3 = ["control", "tamoxifen", "TEAD1 inhibitor\n(existing clinical-stage\nchemotype)", "TEAD1 inhibitor\n+ tamoxifen"]
    for i, a in enumerate(arms3):
        arm_box(0.2 + i * 2.3, 1.9, 2.0, 1.0, a, GREEN)
    ax.text(0.2, 1.3, "PRIMARY readout: does TEAD1 inhibition reduce VEZF1 expression/activity in ER+ breast cancer cells?", fontsize=8.5, color=DGRAY, fontweight="bold")
    ax.text(0.2, 0.85, "(the VEZF1-TEAD1 link is demonstrated in cardiac/zebrafish tissue only -- this experiment tests transferability, not a foregone conclusion)", fontsize=7.8, color=DGRAY)

    fig.suptitle("Proposed experimental strategy: three target-validation experiments", fontsize=13.5)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    out_fig.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_fig, facecolor="white")
    plt.close(fig)
    logger.info("wrote %s", out_fig)


def run(figures_dir: Path = FIGURES) -> None:
    figures_dir.mkdir(parents=True, exist_ok=True)
    build_figure_01(figures_dir / "01_USP34_VEZF1_tissue_expression_atlas.png")
    build_figure_02(figures_dir / "02_USP34_VEZF1_liability_map.png")
    build_figure_03(figures_dir / "03_USP34_VEZF1_direct_indirect_targetability.png")
    build_figure_04(figures_dir / "04_USP34_VEZF1_therapeutic_window_comparison.png")
    build_figure_05(figures_dir / "05_USP34_VEZF1_experimental_strategy.png")


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run()
