"""GSE245601 candidate deep-dive Phases 4-9 figures: patient pseudobulk
(all-epithelial and malignant), patient direction matrix, malignancy x
treatment response, and prevalence-vs-intensity decomposition.

Every "line" in every patient-level plot connects that ONE tumor's own
Control pseudobulk point to that SAME tumor's own Tamoxifen pseudobulk
point -- never individual cells. This is stated in every figure's title/
caption, never omitted.
"""

from __future__ import annotations

import logging
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.patches as mpatches  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import yaml  # noqa: E402

from src.gse245601_candidate_deepdive_data import GENES

logger = logging.getLogger(__name__)

PATIENT_ORDER = [f"Tumor_{i:02d}" for i in range(1, 11)]
PATIENT_COLORS = dict(zip(PATIENT_ORDER, plt.cm.tab10(np.linspace(0, 1, 10))))


def _load_config(config_path: str | Path) -> dict:
    with open(config_path) as f:
        return yaml.safe_load(f)


def _frozen_annotation(gene: str, frozen: pd.DataFrame) -> str:
    row = frozen.set_index("gene").loc[gene]
    return f"frozen Track A edgeR: log2FC={row['gse245601_epi_log2fc']:.3f}, FDR={row['gse245601_epi_fdr']:.3g}"


def plot_patient_pseudobulk(pb: pd.DataFrame, gene: str, frozen: pd.DataFrame, direction_summary: pd.DataFrame, out_path: Path, low_cell_flag: bool = False) -> None:
    sub = pb.loc[pb["gene"] == gene].set_index("patient").reindex([p for p in PATIENT_ORDER if p in pb["patient"].values])
    fig, ax = plt.subplots(figsize=(6, 4.5))
    for patient, row in sub.iterrows():
        if pd.isna(row.get("normalized_expression_Control")) or pd.isna(row.get("normalized_expression_Tamoxifen")):
            continue
        color = PATIENT_COLORS[patient]
        low = low_cell_flag and (row.get("low_cell_count_warning_Control", False) or row.get("low_cell_count_warning_Tamoxifen", False))
        ls = "--" if low else "-"
        alpha = 0.5 if low else 0.9
        ax.plot([0, 1], [row["normalized_expression_Control"], row["normalized_expression_Tamoxifen"]], color=color, marker="o", markersize=6, linewidth=1.5, alpha=alpha, linestyle=ls, label=patient)

    ax.set_xticks([0, 1])
    ax.set_xticklabels(["Control", "Tamoxifen\n(12h)"])
    ax.set_xlim(-0.3, 1.3)
    ax.set_ylabel("pseudobulk log2(CPM+1)")

    d = direction_summary.set_index("gene").loc[gene]
    n_tot = int(d["n_patients_total"])
    title = f"{gene}: patient-level pseudobulk, Control→Tamoxifen\n(each line = ONE tumor's matched pseudobulk pair, NOT individual cells)\n{int(d['n_patients_increase'])}/{n_tot} tumors increase, {int(d['n_patients_decrease'])}/{n_tot} decrease\n{_frozen_annotation(gene, frozen)}"
    ax.set_title(title, fontsize=9)
    ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), fontsize=7, frameon=False, title="tumor" + (" (dashed = <50 malignant cells in that arm)" if low_cell_flag else ""), title_fontsize=6.5)
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info("wrote %s", out_path)


def plot_patient_direction_matrix(all_epi_pb: pd.DataFrame, out_path: Path) -> None:
    genes = GENES
    fig, ax = plt.subplots(figsize=(8, 0.6 * len(genes) + 1.5))
    grid = np.zeros((len(genes), len(PATIENT_ORDER)))
    for gi, gene in enumerate(genes):
        sub = all_epi_pb.loc[all_epi_pb["gene"] == gene].set_index("patient")
        for pi, patient in enumerate(PATIENT_ORDER):
            if patient not in sub.index or pd.isna(sub.loc[patient, "log_fold_change_descriptive"]):
                grid[gi, pi] = np.nan
                continue
            grid[gi, pi] = sub.loc[patient, "log_fold_change_descriptive"]

    vmax = np.nanmax(np.abs(grid))
    im = ax.imshow(grid, cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
    for gi in range(len(genes)):
        for pi in range(len(PATIENT_ORDER)):
            v = grid[gi, pi]
            if np.isnan(v):
                ax.text(pi, gi, "NA", ha="center", va="center", fontsize=7, color="gray")
            else:
                arrow = "↑" if v > 0 else ("↓" if v < 0 else "→")
                ax.text(pi, gi, arrow, ha="center", va="center", fontsize=11, color="black")
    ax.set_xticks(range(len(PATIENT_ORDER)))
    ax.set_xticklabels(PATIENT_ORDER, rotation=45, ha="right", fontsize=8)
    ax.set_yticks(range(len(genes)))
    ax.set_yticklabels(genes, fontsize=9)
    fig.colorbar(im, ax=ax, label="descriptive log2FC (Tam − Control pseudobulk)", shrink=0.7)
    ax.set_title("Patient direction matrix (all-epithelial pseudobulk)\narrow = direction, color = magnitude; NEVER a per-cell comparison", fontsize=9.5)
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info("wrote %s", out_path)


def plot_malignancy_vs_nonmalignancy_response(malig_cond: pd.DataFrame, gene: str, out_path: Path) -> None:
    sub = malig_cond.loc[(malig_cond["gene"] == gene) & (malig_cond["malignancy_status"] != "all_epithelial")]
    fig, axes = plt.subplots(1, 2, figsize=(9, 4.5), sharey=True)
    for ax, malignancy, title in [(axes[0], "non-malignant epithelial", "Non-malignant epithelial"), (axes[1], "malignant", "Malignant")]:
        m = sub.loc[sub["malignancy_status"] == malignancy]
        wide = m.pivot_table(index="patient", columns="condition", values=["pseudobulk_normalized_expression", "n_cells"])
        for patient in PATIENT_ORDER:
            if patient not in wide.index:
                continue
            row = wide.loc[patient]
            c_val, t_val = row.get(("pseudobulk_normalized_expression", "Control")), row.get(("pseudobulk_normalized_expression", "Tamoxifen"))
            c_n, t_n = row.get(("n_cells", "Control"), 0), row.get(("n_cells", "Tamoxifen"), 0)
            if pd.isna(c_val) or pd.isna(t_val):
                continue
            low = (malignancy == "malignant") and (c_n < 50 or t_n < 50)
            ax.plot([0, 1], [c_val, t_val], color=PATIENT_COLORS[patient], marker="o", markersize=5, linewidth=1.3, alpha=0.5 if low else 0.9, linestyle="--" if low else "-", label=patient)
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["Control", "Tam (12h)"])
        ax.set_xlim(-0.3, 1.3)
        ax.set_title(title, fontsize=9.5)
        if ax is axes[0]:
            ax.set_ylabel("pseudobulk log2(CPM+1)")
    axes[1].legend(loc="center left", bbox_to_anchor=(1.02, 0.5), fontsize=6.5, frameon=False, title="tumor\n(dashed=<50 malignant\ncells in that arm)", title_fontsize=6)
    fig.suptitle(f"{gene}: malignant vs. non-malignant epithelial response (each line = one tumor's pseudobulk pair)", fontsize=10, y=1.03)
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info("wrote %s", out_path)


def plot_prevalence_vs_intensity(prevalence_intensity: pd.DataFrame, gene: str, out_path: Path) -> None:
    sub = prevalence_intensity.loc[(prevalence_intensity["gene"] == gene) & (prevalence_intensity["malignancy_status"] == "all_epithelial")]
    fig, axes = plt.subplots(1, 2, figsize=(9, 4.2))

    ax = axes[0]
    wide = sub.pivot_table(index="patient", columns="condition", values="fraction_expressing")
    for patient in PATIENT_ORDER:
        if patient not in wide.index:
            continue
        row = wide.loc[patient]
        if pd.isna(row.get("Control")) or pd.isna(row.get("Tamoxifen")):
            continue
        ax.plot([0, 1], [100 * row["Control"], 100 * row["Tamoxifen"]], color=PATIENT_COLORS[patient], marker="o", markersize=5, linewidth=1.3, alpha=0.85)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["Control", "Tam"])
    ax.set_ylabel("% cells expressing")
    ax.set_title("Prevalence\n(% cells with raw count>0)", fontsize=9)
    ax.set_xlim(-0.3, 1.3)

    ax = axes[1]
    wide2 = sub.pivot_table(index="patient", columns="condition", values="mean_normalized_positive_cells_only")
    for patient in PATIENT_ORDER:
        if patient not in wide2.index:
            continue
        row = wide2.loc[patient]
        if pd.isna(row.get("Control")) or pd.isna(row.get("Tamoxifen")):
            continue
        ax.plot([0, 1], [row["Control"], row["Tamoxifen"]], color=PATIENT_COLORS[patient], marker="o", markersize=5, linewidth=1.3, alpha=0.85, label=patient)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["Control", "Tam"])
    ax.set_ylabel("mean log-norm expression\n(positive cells only)")
    ax.set_title("Intensity among expressing cells", fontsize=9)
    ax.set_xlim(-0.3, 1.3)
    ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), fontsize=6.5, frameon=False, title="tumor", title_fontsize=6.5)

    fig.suptitle(f"{gene}: prevalence vs. intensity (all epithelial, patient-level; each line = one tumor)", fontsize=10.5, y=1.03)
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info("wrote %s", out_path)


def run_visualization_1(config_path: str | Path = "config/config.yaml") -> None:
    config = _load_config(config_path)
    tables_dir = Path(config["gse245601_candidate_deepdive"]["output"]["tables_dir"])
    figures_dir = Path(config["gse245601_candidate_deepdive"]["output"]["figures_dir"])
    figures_dir.mkdir(parents=True, exist_ok=True)

    all_epi_pb = pd.read_csv(tables_dir / "patient_all_epithelial_pseudobulk.tsv", sep="\t")
    malignant_pb = pd.read_csv(tables_dir / "patient_malignant_pseudobulk.tsv", sep="\t")
    direction_summary = pd.read_csv(tables_dir / "patient_direction_summary.tsv", sep="\t")
    malig_cond = pd.read_csv(tables_dir / "malignancy_condition_patient_summary.tsv", sep="\t")
    prevalence_intensity = pd.read_csv(tables_dir / "expression_prevalence_intensity.tsv", sep="\t")
    frozen = pd.read_csv("results/tables/evidence_freeze/final_candidate_evidence.tsv", sep="\t")

    for i, gene in enumerate(GENES, start=1):
        plot_patient_pseudobulk(all_epi_pb, gene, frozen, direction_summary, figures_dir / f"{i:02d}_{gene}_patient_pseudobulk_all_epithelial.png")

    plot_patient_direction_matrix(all_epi_pb, figures_dir / "05_patient_direction_matrix.png")

    for i, gene in enumerate(GENES, start=6):
        plot_patient_pseudobulk(malignant_pb, gene, frozen, direction_summary, figures_dir / f"{i:02d}_{gene}_patient_pseudobulk_malignant.png", low_cell_flag=True)

    for i, gene in enumerate(GENES, start=10):
        plot_malignancy_vs_nonmalignancy_response(malig_cond, gene, figures_dir / f"{i:02d}_{gene}_malignant_vs_nonmalignant_response.png")

    for i, gene in enumerate(GENES, start=14):
        plot_prevalence_vs_intensity(prevalence_intensity, gene, figures_dir / f"{i:02d}_{gene}_prevalence_vs_intensity.png")


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run_visualization_1()
