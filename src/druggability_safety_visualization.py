"""Druggability + normal-tissue/selectivity review: four poster-friendly
figures built only from this phase's own output tables (never hand-typed
numbers), matching the Okabe-Ito colorblind-safe palette already used by
src/independent_validation_visualization.py.

USP34, VEZF1, EML5, CITED2 -- frozen candidates and frozen ranking,
unchanged; nothing here reorders or rescoring the therapeutic shortlist.
"""

from __future__ import annotations

import logging
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

logger = logging.getLogger(__name__)

CANDIDATES = ["USP34", "VEZF1", "EML5", "CITED2"]
TABLES = Path("results/tables/druggability_safety")
FIGURES = Path("results/figures/druggability_safety")

BLUE, ORANGE, GREEN, RED, GRAY, DGRAY = "#0072B2", "#E69F00", "#009E73", "#b0392f", "#9a9a9a", "#555555"
CANDIDATE_COLOR = {"USP34": BLUE, "VEZF1": ORANGE, "EML5": GREEN, "CITED2": RED}

DRUGGABILITY_TIER = {
    "CURRENTLY_POORLY_DRUGGABLE": 1,
    "INDIRECT_OR_MODALITY_DEPENDENT": 2,
    "POTENTIALLY_DRUGGABLE": 3,
    "DIRECTLY_DRUGGABLE": 4,
}
DRUGGABILITY_TIER_COLOR = {1: GRAY, 2: ORANGE, 3: BLUE, 4: GREEN}

NORMAL_TISSUE_TIER = {"LOWER": 1, "MODERATE": 2, "MODERATE-HIGH": 3}
# DOCUMENTED_CAUSAL_POSTNATAL = mouse conditional-KO active during
# skeletal development/differentiation (e.g. USP34's MSC/pre-osteoblast
# -Cre); DOCUMENTED_ADULT_CAUSAL = conditional deletion genuinely induced
# in already-mature adult tissue (e.g. CITED2's Mx1-Cre HSC studies) --
# kept as distinct tiers per the Codex phase review, which correctly
# flagged that conflating the two overstates USP34's evidence.
BONE_TIER = {"NONE_IDENTIFIED": 0, "INFERRED_ONLY": 1, "DOCUMENTED_CAUSAL_POSTNATAL": 2, "DOCUMENTED_ADULT_CAUSAL": 3}


def _load_tables() -> dict[str, pd.DataFrame]:
    return {
        "druggability": pd.read_csv(TABLES / "candidate_druggability.tsv", sep="\t").set_index("candidate"),
        "normal_tissue": pd.read_csv(TABLES / "candidate_normal_tissue_context.tsv", sep="\t").set_index("candidate"),
        "genetic_constraint": pd.read_csv(TABLES / "candidate_genetic_constraint.tsv", sep="\t").set_index("candidate"),
        "bone": pd.read_csv(TABLES / "candidate_bone_musculoskeletal_context.tsv", sep="\t").set_index("candidate"),
        "window": pd.read_csv(TABLES / "candidate_therapeutic_window_summary.tsv", sep="\t").set_index("candidate"),
    }


def build_figure_01(out_fig: Path) -> None:
    """Direct-druggability classification per candidate: an ordinal dot
    plot (never a fabricated numeric druggability score), annotated with
    the concrete structural/tool-compound evidence each tier rests on."""
    t = _load_tables()
    drug = t["druggability"]

    fig, ax = plt.subplots(figsize=(11.5, 5.2), dpi=200)
    for i, c in enumerate(CANDIDATES):
        row = drug.loc[c]
        tier = DRUGGABILITY_TIER[row["druggability_classification"]]
        ax.scatter(tier, i, s=420, color=CANDIDATE_COLOR[c], edgecolor=DGRAY, linewidth=1.2, zorder=3)
        has_pdb = "PDB" in row["structural_data"] and "No experimental PDB" not in row["structural_data"] and "no experimental pdb" not in row["structural_data"].lower()
        has_tool = not (row["known_ligands_tools_probes_degraders"].strip().upper().startswith("NOT FOUND"))
        annot = f"PDB structure: {'yes' if has_pdb else 'no'}  |  tool compound/probe: {'yes (unoptimized)' if has_tool else 'none verified'}"
        ax.text(tier + 0.12, i, annot, va="center", ha="left", fontsize=8.8, color=DGRAY)
    ax.set_yticks(range(len(CANDIDATES)))
    ax.set_yticklabels(CANDIDATES, fontsize=11, fontweight="bold")
    ax.set_ylim(-0.7, len(CANDIDATES) - 0.3)
    ax.set_xlim(0.5, 4.6)
    ax.set_xticks([1, 2, 3, 4])
    ax.set_xticklabels(
        ["CURRENTLY_POORLY\nDRUGGABLE", "INDIRECT_OR_\nMODALITY_DEPENDENT", "POTENTIALLY_\nDRUGGABLE", "DIRECTLY_\nDRUGGABLE"],
        fontsize=8.3,
    )
    ax.grid(axis="x", color="#e5e5e5", zorder=0)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    ax.set_title(
        "Conservative druggability classification (no candidate reaches DIRECTLY_DRUGGABLE)\n"
        "Enzyme/scaffold vs. zinc-finger TF vs. intrinsically disordered co-regulator require different modalities -- see candidate_druggability.tsv",
        fontsize=10.5,
    )
    fig.tight_layout()
    out_fig.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_fig, facecolor="white")
    plt.close(fig)
    logger.info("wrote %s", out_fig)


def build_figure_02(out_fig: Path) -> None:
    """Normal-tissue expression (GTEx v8 median TPM) across five key organ
    systems, small multiples per candidate, log-scaled (values span four
    orders of magnitude across candidates). Expression breadth is not
    essentiality and is not toxicity -- see report caveats."""
    t = _load_tables()
    tissue = t["normal_tissue"]
    cols = ["gtex_breast_tpm", "gtex_blood_tpm", "gtex_liver_tpm", "gtex_heart_tpm", "gtex_kidney_tpm"]
    labels = ["Breast", "Blood", "Liver", "Heart", "Kidney"]

    fig, axes = plt.subplots(1, 4, figsize=(15, 4.4), dpi=200, sharey=True)
    for ax, c in zip(axes, CANDIDATES):
        vals = tissue.loc[c, cols].astype(float).values
        ax.bar(labels, vals, color=CANDIDATE_COLOR[c], edgecolor=DGRAY, linewidth=0.8)
        ax.set_yscale("log")
        ax.set_title(c, fontsize=12, fontweight="bold", color=CANDIDATE_COLOR[c])
        ax.tick_params(axis="x", labelrotation=35, labelsize=9)
        ax.grid(axis="y", which="major", color="#e5e5e5", zorder=0)
        for spine in ("top", "right"):
            ax.spines[spine].set_visible(False)
        spec = tissue.loc[c, "hpa_tissue_specificity_category"]
        ax.set_xlabel(f"HPA: {spec[:32]}{'...' if len(spec) > 32 else ''}", fontsize=7.6, color=DGRAY)
    axes[0].set_ylabel("GTEx v8 median TPM (log scale)", fontsize=10)
    fig.suptitle(
        "Normal-tissue expression breadth (GTEx v8 + Human Protein Atlas) -- NOT a measure of essentiality or toxicity",
        fontsize=12.5,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.90))
    out_fig.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_fig, facecolor="white")
    plt.close(fig)
    logger.info("wrote %s", out_fig)


def build_figure_03(out_fig: Path) -> None:
    """Therapeutic-window map: cancer-vulnerability axis (real DepMap
    ER+/luminal dependency %) vs. an EXPLORATORY, clearly-labeled
    normal-tissue+bone concern index (sum of two qualitative tiers from
    this phase's own tables -- not a hidden master score, and never
    DepMap-derived). Point size = druggability tier; a star marks
    Hany-CRISPR FDR<0.05."""
    t = _load_tables()
    depmap = pd.read_csv("results/tables/independent_validation/DepMap_candidate_dependency.tsv", sep="\t").set_index("candidate")
    constraint = t["genetic_constraint"]
    window = t["window"]

    def normal_tier(c: str) -> int:
        txt = window.loc[c, "normal_tissue_concern"]
        for label, val in sorted(NORMAL_TISSUE_TIER.items(), key=lambda kv: -len(kv[0])):
            if txt.startswith(label):
                return val
        raise ValueError(f"unrecognized normal-tissue concern label for {c}: {txt[:30]}")

    def bone_tier(c: str) -> int:
        cat = t["bone"].loc[c, "bone_concern_category"]
        return BONE_TIER[cat]

    fig, ax = plt.subplots(figsize=(10.5, 7.2), dpi=200)
    jitter = {"USP34": -1.3, "VEZF1": 0.0, "EML5": 0.0, "CITED2": 1.3}
    label_offset = {"USP34": (18, -58), "VEZF1": (0, 32), "EML5": (0, 32), "CITED2": (80, 32)}
    label_ha = {"USP34": "left", "VEZF1": "center", "EML5": "center", "CITED2": "left"}
    for c in CANDIDATES:
        x = float(depmap.loc[c, "frac_strongly_dependent_er_luminal"]) * 100 + jitter[c]
        y = normal_tier(c) + bone_tier(c)
        tier = DRUGGABILITY_TIER[t["druggability"].loc[c, "druggability_classification"]]
        size = 260 + 220 * tier
        fdr = float(depmap.loc[c, "hany_crispr_fdr"])
        ax.scatter(x, y, s=size, color=CANDIDATE_COLOR[c], edgecolor=DGRAY, linewidth=1.3, zorder=3, alpha=0.88)
        marker = " *" if fdr < 0.05 else ""
        ax.annotate(
            f"{c}{marker}\nDepMap ER+/lum={x - jitter[c]:.1f}%\nHany FDR={fdr:.3f}",
            (x, y), textcoords="offset points", xytext=label_offset[c], ha=label_ha[c], va="center",
            fontsize=8.6, color=DGRAY,
            bbox=dict(boxstyle="round,pad=0.25", facecolor="white", edgecolor="none", alpha=0.75),
            arrowprops=dict(arrowstyle="-", color="#bbbbbb", lw=0.8),
        )
    ax.set_xlim(-6, 36)
    ax.set_ylim(-0.5, 7.0)
    ax.set_xlabel("Cancer-cell vulnerability: DepMap ER+/luminal strongly-dependent fraction (%, real data)", fontsize=10)
    ax.set_ylabel("Normal-tissue + bone/musculoskeletal concern index\n(exploratory sum of two curated qualitative tiers, 0-6; see report)", fontsize=9.5)
    ax.grid(color="#eeeeee", zorder=0)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    ax.set_title(
        "Therapeutic-window map: real cancer-dependency data (x) vs. exploratory normal-tissue/bone concern (y)\n"
        "* = Hany CRISPR FDR<0.05 (frozen). Point size = druggability tier. An explicit, labeled composite -- not a hidden score.",
        fontsize=10.3,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.90))
    out_fig.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_fig, facecolor="white")
    plt.close(fig)
    logger.info("wrote %s", out_fig)


def build_figure_04(out_fig: Path) -> None:
    """Bone/musculoskeletal context: expression (left, HPA bulk-tissue
    nTPM, missing data shown as an explicit gap) kept strictly separate
    from functional evidence (right, categorical, species-labeled --
    never implying low expression means safety, never mixing species
    without labeling)."""
    t = _load_tables()
    bone = t["bone"]

    fig, axes = plt.subplots(1, 2, figsize=(14.5, 5.6), dpi=200)

    ax = axes[0]
    x = range(len(CANDIDATES))
    width = 0.35
    bm = [bone.loc[c, "bone_marrow_expression_hpa_ntpm"] for c in CANDIDATES]
    sm = [bone.loc[c, "skeletal_muscle_expression_hpa_ntpm"] for c in CANDIDATES]
    bm_missing = [pd.isna(v) for v in bm]
    sm_missing = [pd.isna(v) for v in sm]
    bm_plot = [0 if pd.isna(v) else v for v in bm]
    sm_plot = [0 if pd.isna(v) else v for v in sm]
    bars1 = ax.bar([i - width / 2 for i in x], bm_plot, width, label="Bone marrow (HPA nTPM)", color=BLUE, edgecolor=DGRAY)
    bars2 = ax.bar([i + width / 2 for i in x], sm_plot, width, label="Skeletal muscle (HPA nTPM)", color=ORANGE, edgecolor=DGRAY)
    for i, (miss, bar) in enumerate(zip(bm_missing, bars1)):
        if miss:
            ax.text(bar.get_x() + bar.get_width() / 2, 1.5, "no data", ha="center", va="bottom", fontsize=8, color=RED, rotation=90)
    for i, (miss, bar) in enumerate(zip(sm_missing, bars2)):
        if miss:
            ax.text(bar.get_x() + bar.get_width() / 2, 1.5, "no data", ha="center", va="bottom", fontsize=8, color=RED, rotation=90)
    ax.set_xticks(list(x))
    ax.set_xticklabels(CANDIDATES, fontsize=10.5, fontweight="bold")
    ax.set_ylabel("HPA consensus nTPM (bulk tissue)", fontsize=10)
    ax.legend(fontsize=8.5, frameon=False)
    ax.grid(axis="y", color="#eeeeee", zorder=0)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    ax.set_title("EXPRESSION only (bulk tissue; no osteoblast/osteoclast/chondrocyte\ncell-type data exists in HPA for any candidate)", fontsize=9.8)

    ax2 = axes[1]
    order = BONE_TIER
    species_label = {
        "USP34": "human MSC + mouse (causal, postnatal dev.)",
        "VEZF1": "mouse only (embryonic lineage, inferred)",
        "EML5": "no evidence in any species",
        "CITED2": "human + mouse + rat (causal, adult, 4 axes)",
    }
    for i, c in enumerate(CANDIDATES):
        cat = bone.loc[c, "bone_concern_category"]
        y = order[cat]
        ax2.scatter(y, i, s=420, color=CANDIDATE_COLOR[c], edgecolor=DGRAY, linewidth=1.2, zorder=3)
        ax2.text(y + 0.08, i, species_label[c], va="center", fontsize=8.8, color=DGRAY)
    ax2.set_yticks(range(len(CANDIDATES)))
    ax2.set_yticklabels(CANDIDATES, fontsize=10.5, fontweight="bold")
    ax2.set_xlim(-0.4, 4.4)
    ax2.set_ylim(-0.7, len(CANDIDATES) - 0.3)
    ax2.set_xticks([0, 1, 2, 3])
    ax2.set_xticklabels(["NONE\nIDENTIFIED", "INFERRED\nONLY", "DOCUMENTED\n(postnatal dev.)", "DOCUMENTED\n(adult-onset)"], fontsize=8.3)
    ax2.grid(axis="x", color="#eeeeee", zorder=0)
    for spine in ("top", "right"):
        ax2.spines[spine].set_visible(False)
    ax2.set_title("FUNCTIONAL bone/skeletal role (published, species labeled)\nNever extrapolated across species; embryonic-lethal != adult toxicity", fontsize=9.8)

    fig.suptitle("Bone / musculoskeletal context: expression and function are shown separately, on purpose", fontsize=12.5)
    fig.tight_layout(rect=(0, 0, 1, 0.90))
    out_fig.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_fig, facecolor="white")
    plt.close(fig)
    logger.info("wrote %s", out_fig)


def run(figures_dir: Path = FIGURES) -> None:
    figures_dir.mkdir(parents=True, exist_ok=True)
    build_figure_01(figures_dir / "01_four_candidate_druggability_summary.png")
    build_figure_02(figures_dir / "02_four_candidate_normal_tissue_context.png")
    build_figure_03(figures_dir / "03_four_candidate_therapeutic_window_map.png")
    build_figure_04(figures_dir / "04_four_candidate_bone_musculoskeletal_context.png")


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run()
