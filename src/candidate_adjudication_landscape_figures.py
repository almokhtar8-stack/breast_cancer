"""Candidate adjudication Phases 24-26: the three summary figures for the
full ~20-30 gene adjudicated candidate pool -- an evidence landscape grid,
a function-vs-resistance scatter map, and a seven-gene head-to-head panel.
No overall-score bar anywhere (Phase 26's explicit instruction); no
favorite genes hand-highlighted outside what the adjudicated tables
already single out.

Data source: `results/tables/candidate_adjudication/final_candidate_decision_table.tsv`,
`three_axis_candidate_matrix.tsv`, `multimodal7_exact_evidence.tsv`.
"""

from __future__ import annotations

import logging
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import yaml  # noqa: E402

logger = logging.getLogger(__name__)

ARCHETYPE_COLORS = {
    "A_FUNCTIONAL_RESISTANCE_CONVERGENCE": "#C44E52", "B_FUNCTIONAL_HUMAN_CONTEXT": "#DD8452",
    "C_FUNCTIONAL_ONLY": "#8172B2", "D_RESISTANCE_BIOMARKER_PATHWAY": "#4C72B0",
    "E_HUMAN_RECURRENCE_DOMINANT": "#55A868", "F_ACUTE_RESPONSE_DOMINANT": "#64B5CD",
    "G_CONTEXT_DEPENDENT": "#CCB974", "H_LOW_INSUFFICIENT_EVIDENCE": "#8C8C8C",
}


def _load_config(config_path: str | Path) -> dict:
    with open(config_path) as f:
        return yaml.safe_load(f)


def plot_evidence_landscape(decision_table: pd.DataFrame, out_path: Path) -> None:
    df = decision_table.sort_values("global_rank", na_position="last").reset_index(drop=True)
    genes = df["gene"].tolist()
    cols = [
        ("crispr_fdr", "CRISPR", lambda v: v < 0.10),
        ("gse118713_fdr", "GSE118713", lambda v: v < 0.05),
        ("gse240112_fdr", "GSE240112", lambda v: v < 0.05),
        ("gse111151_fdr", "GSE111151", lambda v: v < 0.05),
    ]
    grid = np.zeros((len(genes), len(cols)))
    for i, (col, _, sig_fn) in enumerate(cols):
        for j, gene in enumerate(genes):
            v = df.loc[j, col]
            grid[j, i] = 1 if pd.notna(v) and sig_fn(v) else 0

    fig, (ax_grid, ax_arch) = plt.subplots(1, 2, figsize=(9, 0.4 * len(genes) + 2), gridspec_kw={"width_ratios": [3, 1]})
    ax_grid.imshow(grid, cmap="Greens", vmin=0, vmax=1.3, aspect="auto")
    ax_grid.set_xticks(range(len(cols)))
    ax_grid.set_xticklabels([c[1] for c in cols], fontsize=8)
    ax_grid.set_yticks(range(len(genes)))
    ax_grid.set_yticklabels(genes, fontsize=7.5)
    for j in range(len(genes)):
        for i in range(len(cols)):
            ax_grid.text(i, j, "Y" if grid[j, i] else "-", ha="center", va="center", fontsize=7)
    ax_grid.set_title("FDR<0.05 (FDR<0.10 for CRISPR)", fontsize=9)

    ax_arch.axis("off")
    for j, gene in enumerate(genes):
        arch = df.loc[j, "candidate_archetype"]
        color = ARCHETYPE_COLORS.get(arch, "#8C8C8C")
        ax_arch.barh(j, 1, color=color, height=0.8)
        ax_arch.text(1.05, j, arch.split("_", 1)[0], fontsize=7, va="center")
    ax_arch.set_ylim(-0.5, len(genes) - 0.5)
    ax_arch.invert_yaxis()
    ax_arch.set_title("archetype", fontsize=9)

    fig.suptitle("Candidate evidence landscape (sorted by global rank; no overall score)", fontsize=11, y=1.01)
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info("wrote %s", out_path)


def plot_function_vs_resistance_map(axes_matrix: pd.DataFrame, wide: pd.DataFrame, multimodal7: list[str], out_path: Path) -> None:
    band_to_pct = {"VERY_STRONG": 1.0, "STRONG": 0.75, "MODERATE": 0.5, "WEAK": 0.25, "DISCORDANT": 0.25, "NO_EVIDENCE": 0.0}
    df = axes_matrix.merge(wide[["gene", "gse245601_evidence_percentile", "gse240112_evidence_percentile"]], on="gene", how="left")
    x = df["axis_a_functional"].map(band_to_pct) + np.random.default_rng(20260812).uniform(-0.05, 0.05, len(df))
    y = df["axis_b_resistance"].map(band_to_pct) + np.random.default_rng(20260813).uniform(-0.05, 0.05, len(df))
    size = 30 + 200 * df[["gse245601_evidence_percentile", "gse240112_evidence_percentile"]].max(axis=1).fillna(0)
    is_multimodal7 = df["gene"].isin(multimodal7)

    fig, ax = plt.subplots(figsize=(7.5, 6.5))
    ax.scatter(x[~is_multimodal7], y[~is_multimodal7], s=size[~is_multimodal7], color="#8172B2", alpha=0.6, edgecolor="white", linewidth=0.5, label="other adjudicated candidates")
    ax.scatter(x[is_multimodal7], y[is_multimodal7], s=size[is_multimodal7], color="#C44E52", alpha=0.85, edgecolor="black", linewidth=0.8, label="MULTIMODAL_STRONG (7)")
    for xi, yi, gene in zip(x, y, df["gene"]):
        ax.annotate(gene, (xi, yi), fontsize=6.5, xytext=(3, 3), textcoords="offset points")
    ax.set_xlabel("functional (CRISPR) evidence band ->")
    ax.set_ylabel("resistance-state RNA evidence band ->")
    ax.set_xlim(-0.15, 1.15)
    ax.set_ylim(-0.15, 1.15)
    ax.axvline(0.5, color="gray", linewidth=0.5, linestyle="--")
    ax.axhline(0.5, color="gray", linewidth=0.5, linestyle="--")
    ax.text(0.98, 0.98, "functional +\nresistance\nconvergence", ha="right", va="top", fontsize=8, color="gray", transform=ax.transAxes)
    ax.text(0.98, 0.02, "functional only", ha="right", va="bottom", fontsize=8, color="gray", transform=ax.transAxes)
    ax.text(0.02, 0.98, "resistance-RNA\nonly", ha="left", va="top", fontsize=8, color="gray", transform=ax.transAxes)
    ax.legend(loc="lower left", fontsize=8, frameon=False)
    ax.set_title("Marker size = strongest human-tumor percentile (GSE245601/GSE240112)", fontsize=9)
    fig.suptitle("Function vs. resistance-state RNA map", fontsize=11)
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info("wrote %s", out_path)


def plot_multimodal7_head_to_head(master: pd.DataFrame, out_path: Path) -> None:
    genes = master["gene"].tolist()
    fig, axes = plt.subplots(1, 4, figsize=(13, 0.5 * len(genes) + 1.5))

    ax = axes[0]
    colors = ["#C44E52" if d == "sensitising_KO" else "#4C72B0" for d in master["crispr_direction_class"]]
    ax.barh(genes, -np.log10(master["crispr_fdr"]), color=colors)
    ax.axvline(-np.log10(0.10), color="black", linestyle="--", linewidth=0.7)
    ax.set_title("CRISPR\n-log10(FDR)\nred=sensitising", fontsize=8.5)
    ax.invert_yaxis()

    ax = axes[1]
    ax.barh(genes, master["resistance_fdr05_count"], color="#55A868")
    ax.set_xlim(0, 3)
    ax.set_title("resistance datasets\nFDR<0.05 (of 3)", fontsize=8.5)
    ax.invert_yaxis()
    ax.set_yticklabels([])

    ax = axes[2]
    human_sig = ((master["gse245601_epi_fdr"] < 0.05) | (master["gse240112_tumor_fdr"] < 0.05)).astype(int)
    ax.barh(genes, human_sig, color="#DD8452")
    ax.set_xlim(0, 1.2)
    ax.set_xticks([0, 1])
    ax.set_title("any human-tumor\ndataset FDR<0.05", fontsize=8.5)
    ax.invert_yaxis()
    ax.set_yticklabels([])

    ax = axes[3]
    stability_colors = {"ROBUST": "#2E7D32", "MODERATELY_STABLE": "#F9A825", "DATASET_DEPENDENT": "#C62828"}
    colors = [stability_colors.get(s, "gray") for s in master["global_stability_class"]]
    ax.barh(genes, [1] * len(genes), color=colors)
    ax.set_xticks([])
    ax.set_title("leave-one-out\nstability", fontsize=8.5)
    ax.invert_yaxis()
    ax.set_yticklabels([])
    handles = [plt.Rectangle((0, 0), 1, 1, color=c) for c in stability_colors.values()]
    ax.legend(handles, stability_colors.keys(), loc="upper center", bbox_to_anchor=(0.5, -0.08), fontsize=6.5, ncol=1, frameon=False)

    fig.suptitle("Seven MULTIMODAL_STRONG genes head-to-head (no overall-score bar)", fontsize=11, y=1.02)
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info("wrote %s", out_path)


def run_landscape_figures(config_path: str | Path = "config/config.yaml") -> None:
    config = _load_config(config_path)
    adj = config["candidate_adjudication"]["output"]
    tables_dir = Path(adj["tables_dir"])
    final_review = Path(adj["final_review_dir"])
    final_review.mkdir(parents=True, exist_ok=True)

    cdx_out = config["cross_dataset_genomewide"]["output"]
    wide = pd.read_csv(Path(cdx_out["wide_matrix_tsv"]).parent / "all_genes_cross_dataset_evidence_with_ranking.tsv", sep="\t")

    decision_table = pd.read_csv(tables_dir / "final_candidate_decision_table.tsv", sep="\t")
    axes_matrix = pd.read_csv(tables_dir / "three_axis_candidate_matrix.tsv", sep="\t")
    master = pd.read_csv(tables_dir / "multimodal7_exact_evidence.tsv", sep="\t")
    multimodal7 = config["candidate_adjudication"]["multimodal7"]["genes"]

    plot_evidence_landscape(decision_table, final_review / "05_candidate_evidence_landscape.png")
    plot_function_vs_resistance_map(axes_matrix, wide, multimodal7, final_review / "06_function_vs_resistance_map.png")
    plot_multimodal7_head_to_head(master, final_review / "07_multimodal7_head_to_head.png")


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run_landscape_figures()
