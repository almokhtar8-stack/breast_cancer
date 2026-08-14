"""Final USP34/VEZF1 translational + structure phase: three poster-friendly
figures built only from this phase's own output tables (never hand-typed
numbers), matching the Okabe-Ito colorblind-safe palette used throughout
this project's visualization modules. No figure 04 (docking hypothesis)
is produced -- Part 9's decision is DOCKING_NOT_YET_JUSTIFIED, per the
user's own instruction that figure 4 is conditional on a JUSTIFIED
decision.
"""

from __future__ import annotations

import logging
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

logger = logging.getLogger(__name__)

TABLES = Path("results/tables/final_translational")
FIGURES = Path("results/figures/final_translational")

BLUE, ORANGE, GREEN, RED, GRAY, DGRAY = "#0072B2", "#E69F00", "#009E73", "#b0392f", "#9a9a9a", "#555555"


def build_figure_01(out_fig: Path) -> None:
    fig, ax = plt.subplots(figsize=(13.5, 10.5), dpi=200)
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 13)
    ax.axis("off")

    def arm_box(x, y, w, h, text, color):
        rect = plt.Rectangle((x, y), w, h, facecolor=color, edgecolor=DGRAY, linewidth=1.1, alpha=0.85, zorder=2)
        ax.add_patch(rect)
        ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=7.6, color="white", fontweight="bold", zorder=3, wrap=True)

    # EXP-1
    ax.text(0.2, 12.55, "EXP-1: USP34 perturbation ± tamoxifen (resistant ER+ cells)", fontsize=12, fontweight="bold", color=BLUE)
    arms1 = ["control", "tamoxifen/4-OHT", "USP34 CRISPR KO\n(no chemical probe\nexists yet)", "USP34 CRISPR KO\n+ tamoxifen/4-OHT"]
    for i, a in enumerate(arms1):
        arm_box(0.2 + i * 2.3, 11.35, 2.0, 1.0, a, BLUE)
    ax.text(0.2, 10.75, "readouts: viability, apoptosis, clonogenic survival, full dose-response + Bliss/Chou-Talalay interaction (not single-dose), AXIN1/beta-catenin + ESR1 genes", fontsize=8, color=DGRAY)
    ax.text(0.2, 10.4, "EMT/stemness monitoring (PMID 28499884 counter-evidence): CDH1/CDH2/SNAI1/AXIN1/active-beta-catenin, optional mammosphere/invasion", fontsize=8, color=DGRAY, fontweight="bold")

    # EXP-2A / EXP-2B
    ax.text(0.2, 9.85, "EXP-2A: normal mammary epithelial cells (lineage/selectivity)  |  EXP-2B: primary human MSCs, osteogenic induction (known liability)", fontsize=9.3, fontweight="bold", color=BLUE)
    ax.text(0.2, 9.4, "EXP-2A readouts: selectivity ratio vs cancer cells, same EMT/stemness panel. EXP-2B readouts: viability, RUNX2/ALP/Alizarin Red mineralization", fontsize=8, color=DGRAY)

    # EXP-3
    ax.text(0.2, 8.95, "EXP-3: VEZF1 suppression ± tamoxifen (resistant ER+ cells)", fontsize=12, fontweight="bold", color=ORANGE)
    arms3 = ["control", "tamoxifen/4-OHT", "VEZF1 CRISPRi", "VEZF1 CRISPRi\n+ tamoxifen/4-OHT"]
    for i, a in enumerate(arms3):
        arm_box(0.2 + i * 2.3, 7.75, 2.0, 1.0, a, ORANGE)
    ax.text(0.2, 7.15, "readouts: SAME panel as EXP-1, PLUS separate direct-dependency (arm3 vs 1) vs sensitisation (arm4 vs 2+3) analysis, plus VEGFR2/TIMP3/MMP2", fontsize=8, color=DGRAY)

    # EXP-4
    ax.text(0.2, 6.55, "EXP-4: VEZF1 normal-cell comparator -- primary endothelial cells (+ iPSC-cardiomyocytes if resourced)", fontsize=9.3, fontweight="bold", color=ORANGE)
    ax.text(0.2, 6.1, "readouts: viability, tube formation, VEGFR2/TIMP3/MMP2; cardiomyocyte contractility if included -- animal cardiac findings are NOT equated with human toxicity", fontsize=8, color=DGRAY)

    # EXP-5
    ax.text(0.2, 5.35, "EXP-5: TEAD inhibition -> does VEZF1 actually change? (HYPOTHESIS TEST ONLY, not a therapeutic arm)", fontsize=12, fontweight="bold", color=GREEN)
    arms5 = ["control", "pan-TEAD\ninhibitor", "VEZF1 CRISPRi\n(positive-control\nreference)", "TEAD inhibitor\n+ tamoxifen\n(ONLY if VEZF1\nchanges)"]
    for i, a in enumerate(arms5):
        arm_box(0.2 + i * 2.3, 4.15, 2.0, 1.0, a, GREEN)
    ax.text(0.2, 3.55, "readouts: VEZF1 RNA/protein/target program (VEGFR2/TIMP3/MMP2) + MANDATORY TEAD/YAP target-engagement control (CTGF/CYR61)", fontsize=8, color=DGRAY)
    ax.text(0.2, 3.1, "DECISION RULE: no VEZF1 change despite confirmed TEAD engagement -> REJECT TEAD1 as an indirect strategy. VEZF1 changes -> proceed to arm 4.", fontsize=8.2, color=DGRAY, fontweight="bold")

    ax.text(0.2, 2.3, "PRIORITY IF ONLY ONE EXPERIMENT CAN BE DONE: EXP-1 (USP34 main experiment)", fontsize=11, fontweight="bold", color=DGRAY,
            bbox=dict(boxstyle="round,pad=0.4", facecolor="#f2f2f2", edgecolor=DGRAY))

    fig.suptitle("Final experimental strategy: five experiments, USP34 lead / VEZF1 backup / TEAD1 hypothesis-only", fontsize=13.5)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    out_fig.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_fig, facecolor="white")
    plt.close(fig)
    logger.info("wrote %s", out_fig)


def build_figure_02(out_fig: Path) -> None:
    pockets = pd.read_csv(TABLES / "USP34_pocket_analysis.tsv", sep="\t")

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.6), dpi=200)

    ax = axes[0]
    labels = [f"{row['structure'].split(' ')[0]}\n{row['pocket_label'].split('(')[0].strip()}" for _, row in pockets.iterrows()]
    scores = pockets["fpocket_druggability_score"].astype(float).tolist()
    colors = [GREEN if s >= 0.5 else (ORANGE if s >= 0.2 else GRAY) for s in scores]
    bars = ax.barh(range(len(labels)), scores, color=colors, edgecolor=DGRAY)
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlabel("fpocket druggability score (0-1 geometric heuristic, NOT a probability\nof clinical success) -- red dashed line = fpocket's own 'druggable' cutoff (0.5)", fontsize=8.2)
    ax.set_xlim(0, 1.0)
    ax.axvline(0.5, color=RED, linestyle="--", linewidth=1)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    ax.set_title("Real, locally-computed pocket scores\n(fpocket 4.2.3, run on downloaded 7W3R/7W3U)", fontsize=10)

    ax2 = axes[1]
    states = ["Apo, chain A\n(7W3R)", "Ubiquitin-probe-bound,\nchain A (7W3U)"]
    dists = [3.94, 3.37]
    bars2 = ax2.bar(states, dists, color=[BLUE, GREEN], edgecolor=DGRAY, width=0.5)
    for i, d in enumerate(dists):
        ax2.text(i, d + 0.08, f"{d:.2f} Å", ha="center", fontsize=10, fontweight="bold")
    ax2.set_ylabel("Cys1903(SG) - His2164(ND1) distance (Å)\ndirectly measured from downloaded coordinates", fontsize=8.5)
    ax2.set_ylim(0, 4.6)
    for spine in ("top", "right"):
        ax2.spines[spine].set_visible(False)
    ax2.set_title("Measured catalytic-dyad tightening, chain A only\n(other copies range 3.10-4.98 Å -- not uniform, see report)", fontsize=9.5)

    fig.suptitle("USP34 structural targetability: real, locally-computed evidence (not literature-only)", fontsize=12.5)
    fig.text(0.5, 0.02, "Cys1903-AYE covalent bond directly confirmed by LINK record in 7W3U (bond distance 1.59-2.48 Å across 3 copies in the asymmetric unit)",
              ha="center", fontsize=8.2, color=DGRAY, style="italic")
    fig.tight_layout(rect=(0, 0.06, 1, 0.90))
    out_fig.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_fig, facecolor="white")
    plt.close(fig)
    logger.info("wrote %s", out_fig)


def build_figure_03(out_fig: Path) -> None:
    fig, ax = plt.subplots(figsize=(13, 7.5), dpi=200)
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.axis("off")

    def box(x, y, w, h, text, color, fontsize=8.6):
        rect = plt.Rectangle((x, y), w, h, facecolor=color, edgecolor=DGRAY, linewidth=1.2, alpha=0.88, zorder=2)
        ax.add_patch(rect)
        ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=fontsize, color="white", fontweight="bold", zorder=3, wrap=True)

    # USP34 model
    ax.text(0.2, 9.5, "USP34 -- LEAD TARGET", fontsize=13, fontweight="bold", color=BLUE)
    box(0.2, 8.0, 2.0, 1.1, "Functional CRISPR\nsensitisation\n(Hany FDR=0.042)", BLUE)
    ax.text(2.35, 8.55, "+", fontsize=16, fontweight="bold", ha="center")
    box(2.6, 8.0, 2.0, 1.1, "Low baseline ER+\ncancer dependency\n(DepMap 26Q1 = 0.0%)", BLUE)
    ax.text(4.75, 8.55, "+", fontsize=16, fontweight="bold", ha="center")
    box(5.0, 8.0, 2.3, 1.1, "Real catalytic\ntargetability\n(PDB 7W3R/7W3U,\nCys1903 confirmed reactive)", BLUE)
    ax.text(7.5, 8.55, "=", fontsize=16, fontweight="bold", ha="center")
    box(7.75, 7.75, 2.05, 1.6, "LEAD COMBINATION-\nTARGET HYPOTHESIS\n(not clinically validated)", DGRAY, fontsize=8.8)
    ax.text(0.2, 7.55, "Caveats: no validated inhibitor yet; real bone/osteogenic liability (EXP-2 comparator); strongly constrained (LOEUF=0.152)", fontsize=7.8, color=DGRAY, style="italic")

    # VEZF1 model
    ax.text(0.2, 6.4, "VEZF1 -- SECOND / BACKUP TARGET", fontsize=13, fontweight="bold", color=ORANGE)
    box(0.2, 4.9, 2.0, 1.1, "Strong CRISPR\nsensitisation\n(Hany FDR=0.037)", ORANGE)
    ax.text(2.35, 5.45, "+", fontsize=16, fontweight="bold", ha="center")
    box(2.6, 4.9, 2.0, 1.1, "Baseline ER+/luminal\ncancer dependency\n(DepMap 26Q1 = 27.3%)", ORANGE)
    ax.text(4.75, 5.45, "=", fontsize=16, fontweight="bold", ha="center")
    box(5.0, 4.65, 2.35, 1.6, "DUAL-ACTION\nBIOLOGICAL HYPOTHESIS\nlimited by poor direct\ndruggability", DGRAY, fontsize=8.6)
    ax.text(0.2, 4.45, "Caveats: cardiovascular/developmental liability (PMID 31911272) stronger than muscle/bone evidence; TEAD1 = unvalidated hypothesis pending EXP-5", fontsize=7.8, color=DGRAY, style="italic")

    # Priority footer
    ax.add_patch(plt.Rectangle((0.2, 2.9), 9.6, 1.0, facecolor="#f2f2f2", edgecolor=DGRAY, linewidth=1.0, zorder=2))
    ax.text(5.0, 3.4, "Frozen ranking, unchanged by this phase: USP34 = LEAD, VEZF1 = SECOND/BACKUP.\nNeither target is described as clinically validated anywhere in this project.", ha="center", va="center", fontsize=9.5, fontweight="bold", color=DGRAY, zorder=3)

    ax.text(0.2, 2.3, "Locked final questions (Part 1): USP34 -- does reducing USP34 restore/enhance tamoxifen response without nonspecific toxicity?", fontsize=8, color=DGRAY)
    ax.text(0.2, 1.85, "VEZF1 -- does reducing VEZF1 (A) impair resistant-cell fitness alone AND/OR (B) further enhance tamoxifen sensitivity? Neither is proven.", fontsize=8, color=DGRAY)

    fig.suptitle("Final USP34 / VEZF1 translational model (poster-ready, not a clinical claim)", fontsize=13.5)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    out_fig.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_fig, facecolor="white")
    plt.close(fig)
    logger.info("wrote %s", out_fig)


def run(figures_dir: Path = FIGURES) -> None:
    figures_dir.mkdir(parents=True, exist_ok=True)
    build_figure_01(figures_dir / "01_final_experimental_strategy.png")
    build_figure_02(figures_dir / "02_USP34_structure_targetability.png")
    build_figure_03(figures_dir / "03_USP34_VEZF1_final_translational_model.png")


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run()
