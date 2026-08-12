"""Candidate adjudication Phases 12-14, 19, 24-26: the final_review figure
set. CRISPR and RNA effect sizes are never placed on one shared
continuous axis (different units/scales) -- separate panels or
standardized within-dataset percentiles are used instead, per the task's
explicit instruction.

Data source: `results/tables/candidate_adjudication/*.tsv` (built by the
other `src/candidate_adjudication_*.py` modules) and the frozen
cross-dataset genome-wide tables they were built from.
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

DATASET_PCT_COLS = ["crispr_evidence_percentile", "gse118713_evidence_percentile", "gse240112_evidence_percentile", "gse111151_evidence_percentile", "gse245601_epi_percentile", "gse245601_malignant_percentile"]
DATASET_PCT_LABELS = ["CRISPR\n(functional)", "GSE118713\n(resistance)", "GSE240112\n(recurrence)", "GSE111151\n(resistance)", "GSE245601 epi\n(acute)", "GSE245601 malig.\n(acute)"]
DATASET_FDR_COLS = ["crispr_fdr", "gse118713_fdr", "gse240112_tumor_fdr", "gse111151_fdr", "gse245601_epi_fdr", "gse245601_malignant_fdr"]


def _load_config(config_path: str | Path) -> dict:
    with open(config_path) as f:
        return yaml.safe_load(f)


def _prep_multimodal7_percentiles(master: pd.DataFrame, wide: pd.DataFrame) -> pd.DataFrame:
    """Track-A/B percentiles aren't in the master table (only the
    dataset-level mean is) -- pulled fresh from the frozen wide-with
    -ranking table for the heatmap/effect figures, which want each track
    shown separately."""
    extra = wide.loc[wide["gene"].isin(master["gene"]), ["gene", "gse245601_track_a_percentile", "gse245601_track_b_percentile"]]
    extra = extra.rename(columns={"gse245601_track_a_percentile": "gse245601_epi_percentile", "gse245601_track_b_percentile": "gse245601_malignant_percentile"})
    return master.merge(extra, on="gene", how="left")


def plot_multimodal7_cross_dataset_effects(master: pd.DataFrame, out_path: Path) -> None:
    """Figure 01: two side-by-side panels (CRISPR effect size; RNA log2FC
    across the four RNA datasets), never on one shared axis."""
    genes = master["gene"].tolist()
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 0.5 * len(genes) + 2), gridspec_kw={"width_ratios": [1, 2.2]})

    y = np.arange(len(genes))
    crispr_colors = ["#C44E52" if f < 0.05 else "#4C72B0" for f in master["crispr_fdr"]]
    ax1.barh(y, master["crispr_effect"], color=crispr_colors)
    ax1.set_yticks(y)
    ax1.set_yticklabels(genes)
    ax1.axvline(0, color="black", linewidth=0.8)
    ax1.set_xlabel("CRISPR effect size")
    ax1.set_title("Functional (CRISPR)\nred = FDR<0.05", fontsize=9)
    ax1.invert_yaxis()

    rna_datasets = [("gse118713_log2fc", "gse118713_fdr", "GSE118713"), ("gse240112_tumor_log2fc", "gse240112_tumor_fdr", "GSE240112"), ("gse111151_log2fc", "gse111151_fdr", "GSE111151")]
    n_d = len(rna_datasets)
    bar_h = 0.8 / n_d
    for i, (eff_col, fdr_col, label) in enumerate(rna_datasets):
        offsets = y - 0.4 + bar_h * (i + 0.5)
        colors = ["#C44E52" if f is not None and pd.notna(f) and f < 0.05 else "#8CA0C6" for f in master[fdr_col]]
        ax2.barh(offsets, master[eff_col], height=bar_h, color=colors, label=label)
    ax2.set_yticks(y)
    ax2.set_yticklabels(genes)
    ax2.axvline(0, color="black", linewidth=0.8)
    ax2.set_xlabel("log2FC (resistant/recurrent vs. parental/primary)")
    ax2.set_title("Resistance-state RNA (GSE118713 / GSE240112 / GSE111151)\nred = FDR<0.05, blue = not significant", fontsize=9)
    ax2.invert_yaxis()
    dataset_labels_text = "  |  ".join(f"row {i+1}={l}" for i, (_, _, l) in enumerate(rna_datasets))
    ax2.text(0.5, -0.12, f"within each gene's 3-bar group, top-to-bottom: {dataset_labels_text}", transform=ax2.transAxes, fontsize=7, ha="center")

    fig.suptitle("Seven MULTIMODAL_STRONG genes: CRISPR and RNA effects (separate panels -- not a shared scale)", fontsize=11, y=1.02)
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info("wrote %s", out_path)


def plot_multimodal7_evidence_heatmap(master_with_tracks: pd.DataFrame, out_path: Path) -> None:
    """Figure 02: percentile heatmap, FDR<0.05 marked with an asterisk;
    GSE245601's two tracks shown as separate columns (never double-counted
    in any summary statistic elsewhere -- this is a display of the two
    tracks' own percentiles, not a re-vote)."""
    mat = master_with_tracks.set_index("gene")[DATASET_PCT_COLS]
    fdr = master_with_tracks.set_index("gene")[DATASET_FDR_COLS]
    fig, ax = plt.subplots(figsize=(8, 0.5 * len(mat) + 1.5))
    im = ax.imshow(mat.to_numpy(dtype=float), cmap="YlOrRd", vmin=0, vmax=1, aspect="auto")
    ax.set_xticks(range(len(DATASET_PCT_LABELS)))
    ax.set_xticklabels(DATASET_PCT_LABELS, fontsize=8)
    ax.set_yticks(range(len(mat)))
    ax.set_yticklabels(mat.index)
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            v = mat.iloc[i, j]
            f = fdr.iloc[i, j]
            if pd.isna(v):
                ax.text(j, i, "NA", ha="center", va="center", fontsize=7, color="gray")
            else:
                star = "*" if pd.notna(f) and f < 0.05 else ""
                ax.text(j, i, f"{v:.3f}{star}", ha="center", va="center", fontsize=7, color="black" if v < 0.7 else "white")
    fig.colorbar(im, ax=ax, label="within-dataset evidence percentile")
    ax.set_title("* = FDR<0.05 in that dataset/track", fontsize=9)
    fig.suptitle("Seven MULTIMODAL_STRONG genes: evidence percentile by dataset/track", fontsize=11, y=1.02)
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info("wrote %s", out_path)


def plot_multimodal7_support_matrix(master: pd.DataFrame, out_path: Path) -> None:
    """Figure 03: a checkmark/cross significance-and-consistency matrix."""
    genes = master["gene"].tolist()
    columns = [
        ("crispr_fdr_lt_010", "CRISPR\nFDR<0.10"),
        ("gse118713_col", "GSE118713\nFDR<0.05"),
        ("gse240112_col", "GSE240112\nFDR<0.05"),
        ("gse111151_col", "GSE111151\nFDR<0.05"),
        ("gse245601_col", "GSE245601\nFDR<0.05 (either track)"),
        ("resist3of3", "resistance\ndirection 3/3"),
        ("resist2of3", "resistance\ndirection >=2/3"),
        ("human_col", "human tumor\nsupport"),
        ("robust_col", "leave-one-out\nMODERATELY_STABLE+"),
    ]
    m = master.copy()
    m["gse118713_col"] = m["gse118713_fdr"] < 0.05
    m["gse240112_col"] = m["gse240112_tumor_fdr"] < 0.05
    m["gse111151_col"] = m["gse111151_fdr"] < 0.05
    m["gse245601_col"] = (m["gse245601_epi_fdr"] < 0.05) | (m["gse245601_malignant_fdr"] < 0.05)
    m["resist3of3"] = m["resistance_fdr05_count"] >= 3
    m["resist2of3"] = m["resistance_fdr05_count"] >= 2
    m["human_col"] = m["human_tumor_datasets_testable"] > 0
    m["robust_col"] = m["global_stability_class"].isin(["ROBUST", "MODERATELY_STABLE"])

    grid = np.array([[bool(m.loc[m["gene"] == g, col].iloc[0]) for col, _ in columns] for g in genes])
    fig, ax = plt.subplots(figsize=(9, 0.55 * len(genes) + 1.5))
    ax.imshow(grid.astype(int), cmap="Greens", vmin=0, vmax=1.4, aspect="auto")
    ax.set_xticks(range(len(columns)))
    ax.set_xticklabels([lbl for _, lbl in columns], fontsize=7.5)
    ax.set_yticks(range(len(genes)))
    ax.set_yticklabels(genes)
    for i in range(grid.shape[0]):
        for j in range(grid.shape[1]):
            ax.text(j, i, "Y" if grid[i, j] else "-", ha="center", va="center", fontsize=9, color="black")
    fig.suptitle("Seven MULTIMODAL_STRONG genes: significance and consistency matrix", fontsize=11, y=1.02)
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info("wrote %s", out_path)


def plot_modality_dependency(merged: pd.DataFrame, out_path: Path) -> None:
    """Figure 04: global rank vs. RNA-only rank and vs. CRISPR rank, log
    -log scatter (both axes are ranks, a comparable scale by construction,
    unlike raw effect sizes)."""
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.2))
    for ax, col, label, color in [(axes[0], "rna_only_rank", "RNA-only rank", "#4C72B0"), (axes[1], "crispr_rank", "CRISPR rank", "#C44E52")]:
        sub = merged.dropna(subset=["global_rank", col])
        ax.scatter(sub["global_rank"], sub[col], s=4, alpha=0.15, color=color)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("global rank (log scale)")
        ax.set_ylabel(f"{label} (log scale)")
        ax.set_title(label, fontsize=10)
    fig.suptitle("Global ranking is far more correlated with RNA-only rank than with CRISPR rank", fontsize=11, y=1.03)
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info("wrote %s", out_path)


def run_visualization(config_path: str | Path = "config/config.yaml") -> None:
    config = _load_config(config_path)
    adj = config["candidate_adjudication"]["output"]
    tables_dir = Path(adj["tables_dir"])
    final_review = Path(adj["final_review_dir"])
    final_review.mkdir(parents=True, exist_ok=True)

    cdx_out = config["cross_dataset_genomewide"]["output"]
    wide = pd.read_csv(Path(cdx_out["wide_matrix_tsv"]).parent / "all_genes_cross_dataset_evidence_with_ranking.tsv", sep="\t")

    master = pd.read_csv(tables_dir / "multimodal7_exact_evidence.tsv", sep="\t")
    master_with_tracks = _prep_multimodal7_percentiles(master, wide)

    plot_multimodal7_cross_dataset_effects(master, final_review / "01_multimodal7_cross_dataset_effects.png")
    plot_multimodal7_evidence_heatmap(master_with_tracks, final_review / "02_multimodal7_evidence_heatmap.png")
    plot_multimodal7_support_matrix(master, final_review / "03_multimodal7_support_matrix.png")

    from src.candidate_adjudication_modality_dependency import run_modality_dependency

    modality = run_modality_dependency(config_path)
    plot_modality_dependency(modality["merged"], final_review / "04_global_ranking_modality_dependency.png")


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run_visualization()
