"""GSE245601 candidate deep-dive Phases 10-12: descriptive cell-level
distribution plots, the USP34 special decomposition figure, and UMAP
expression maps. All descriptive -- no per-cell p-value anywhere. UMAP
coordinates are the frozen ones from the pseudobulk cell-level summary
table; never recomputed.
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

from src.gse245601_candidate_deepdive_data import GENES, load_per_cell_table

logger = logging.getLogger(__name__)

CONDITION_COLORS = {"Control": "#4C72B0", "Tamoxifen": "#C44E52"}


def _load_config(config_path: str | Path) -> dict:
    with open(config_path) as f:
        return yaml.safe_load(f)


def plot_cell_distribution(per_cell: pd.DataFrame, gene: str, out_path: Path) -> None:
    col = f"{gene}_log_norm"
    groups = [
        ("all epithelial", per_cell),
        ("malignant", per_cell.loc[per_cell["malignancy_status"] == "malignant"]),
        ("non-malignant epithelial", per_cell.loc[per_cell["malignancy_status"] == "non-malignant epithelial"]),
    ]
    fig, axes = plt.subplots(1, 3, figsize=(11, 4.2), sharey=True)
    for ax, (label, grp) in zip(axes, groups):
        data = [grp.loc[grp["condition"] == c, col].to_numpy() for c in ("Control", "Tamoxifen")]
        parts = ax.violinplot(data, showmedians=True, widths=0.8)
        for pc, c in zip(parts["bodies"], ("Control", "Tamoxifen")):
            pc.set_facecolor(CONDITION_COLORS[c])
            pc.set_alpha(0.6)
        ax.set_xticks([1, 2])
        ax.set_xticklabels(["Control", "Tam"])
        n_c, n_t = len(data[0]), len(data[1])
        ax.set_title(f"{label}\n(n={n_c} Control, {n_t} Tam cells)", fontsize=8.5)
    axes[0].set_ylabel(f"{gene} log-normalized expression")
    fig.suptitle(f"{gene}: cell-level expression distributions (DESCRIPTIVE ONLY -- distributions use ALL cells, not a subsample; no per-cell p-value)", fontsize=9.5, y=1.03)
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info("wrote %s", out_path)


def plot_umap_expression(per_cell: pd.DataFrame, gene: str, out_path: Path) -> None:
    col = f"{gene}_log_norm"
    vmax = per_cell[col].quantile(0.99)
    fig, axes = plt.subplots(2, 2, figsize=(9, 8.4))
    panels = [
        ("Control, all epithelial", per_cell.loc[per_cell["condition"] == "Control"]),
        ("Tamoxifen (12h), all epithelial", per_cell.loc[per_cell["condition"] == "Tamoxifen"]),
        ("Control, malignant only", per_cell.loc[(per_cell["condition"] == "Control") & (per_cell["malignancy_status"] == "malignant")]),
        ("Tamoxifen (12h), malignant only", per_cell.loc[(per_cell["condition"] == "Tamoxifen") & (per_cell["malignancy_status"] == "malignant")]),
    ]
    for ax, (title, grp) in zip(axes.flat, panels):
        ax.scatter(per_cell["umap_1"], per_cell["umap_2"], s=2, c="#E8E8E8", linewidths=0)
        order = grp[col].sort_values().index
        sc = ax.scatter(grp.loc[order, "umap_1"], grp.loc[order, "umap_2"], s=3, c=grp.loc[order, col], cmap="viridis", vmin=0, vmax=vmax, linewidths=0)
        ax.set_title(f"{title} (n={len(grp)})", fontsize=8.5)
        ax.set_xticks([])
        ax.set_yticks([])
    fig.colorbar(sc, ax=axes, label=f"{gene} log-norm expression", shrink=0.6)
    fig.suptitle(f"{gene}: frozen UMAP coordinates, gene expression overlay (embedding not recomputed)", fontsize=10.5)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info("wrote %s", out_path)


def plot_usp34_decomposition(per_cell: pd.DataFrame, all_epi_pb: pd.DataFrame, malignant_pb: pd.DataFrame, malig_cond: pd.DataFrame, frozen: pd.DataFrame, out_path: Path) -> None:
    from src.gse245601_candidate_deepdive_visualization_1 import PATIENT_COLORS, PATIENT_ORDER

    gene = "USP34"
    fig = plt.figure(figsize=(16, 11))
    gs = fig.add_gridspec(3, 4, hspace=0.55, wspace=0.4)

    def line_panel(ax, pb, title, low_flag=False):
        sub = pb.loc[pb["gene"] == gene].set_index("patient")
        for patient in PATIENT_ORDER:
            if patient not in sub.index:
                continue
            row = sub.loc[patient]
            if pd.isna(row.get("normalized_expression_Control")) or pd.isna(row.get("normalized_expression_Tamoxifen")):
                continue
            low = low_flag and (row.get("low_cell_count_warning_Control", False) or row.get("low_cell_count_warning_Tamoxifen", False))
            ax.plot([0, 1], [row["normalized_expression_Control"], row["normalized_expression_Tamoxifen"]], color=PATIENT_COLORS[patient], marker="o", markersize=4, linewidth=1.2, alpha=0.4 if low else 0.9, linestyle="--" if low else "-")
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["Ctrl", "Tam"], fontsize=7)
        ax.set_title(title, fontsize=8.5)

    ax = fig.add_subplot(gs[0, 0])
    line_panel(ax, all_epi_pb, "A. All-epithelial\npatient pseudobulk")
    ax.set_ylabel("log2(CPM+1)", fontsize=7)

    ax = fig.add_subplot(gs[0, 1])
    line_panel(ax, malignant_pb, "B. Malignant\npatient pseudobulk", low_flag=True)

    nonmal_pb = malig_cond.loc[(malig_cond["gene"] == gene) & (malig_cond["malignancy_status"] == "non-malignant epithelial")]
    ax = fig.add_subplot(gs[0, 2])
    wide = nonmal_pb.pivot_table(index="patient", columns="condition", values="pseudobulk_normalized_expression")
    for patient in PATIENT_ORDER:
        if patient not in wide.index:
            continue
        row = wide.loc[patient]
        if pd.isna(row.get("Control")) or pd.isna(row.get("Tamoxifen")):
            continue
        ax.plot([0, 1], [row["Control"], row["Tamoxifen"]], color=PATIENT_COLORS[patient], marker="o", markersize=4, linewidth=1.2, alpha=0.9)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["Ctrl", "Tam"], fontsize=7)
    ax.set_title("C. Non-malignant\npatient pseudobulk", fontsize=8.5)

    prevalence = malig_cond.loc[(malig_cond["gene"] == gene) & (malig_cond["malignancy_status"] == "all_epithelial")]
    ax = fig.add_subplot(gs[0, 3])
    wide = prevalence.pivot_table(index="patient", columns="condition", values="fraction_expressing")
    for patient in PATIENT_ORDER:
        if patient not in wide.index:
            continue
        row = wide.loc[patient]
        if pd.isna(row.get("Control")) or pd.isna(row.get("Tamoxifen")):
            continue
        ax.plot([0, 1], [100 * row["Control"], 100 * row["Tamoxifen"]], color=PATIENT_COLORS[patient], marker="o", markersize=4, linewidth=1.2, alpha=0.9)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["Ctrl", "Tam"], fontsize=7)
    ax.set_title("D. %% USP34+ cells\n(all epithelial)", fontsize=8.5)
    ax.set_ylabel("%% expressing", fontsize=7)

    ax = fig.add_subplot(gs[1, 0])
    wide = prevalence.pivot_table(index="patient", columns="condition", values="mean_normalized_positive_cells_only")
    for patient in PATIENT_ORDER:
        if patient not in wide.index:
            continue
        row = wide.loc[patient]
        if pd.isna(row.get("Control")) or pd.isna(row.get("Tamoxifen")):
            continue
        ax.plot([0, 1], [row["Control"], row["Tamoxifen"]], color=PATIENT_COLORS[patient], marker="o", markersize=4, linewidth=1.2, alpha=0.9)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["Ctrl", "Tam"], fontsize=7)
    ax.set_title("E. Expression among\nUSP34+ cells only", fontsize=8.5)

    ax = fig.add_subplot(gs[1, 1])
    mal = malig_cond.loc[(malig_cond["gene"] == gene) & (malig_cond["malignancy_status"].isin(["malignant", "non-malignant epithelial"]))]
    for cond, marker in [("Control", "o"), ("Tamoxifen", "s")]:
        m = mal.loc[mal["condition"] == cond].pivot_table(index="patient", columns="malignancy_status", values="pseudobulk_normalized_expression")
        if "malignant" in m.columns and "non-malignant epithelial" in m.columns:
            ax.scatter(m["non-malignant epithelial"], m["malignant"], c=CONDITION_COLORS[cond], marker=marker, s=25, label=cond, alpha=0.8)
    lims = ax.get_xlim()
    ax.plot(lims, lims, color="gray", linewidth=0.7, linestyle=":")
    ax.set_xlabel("non-malignant pseudobulk", fontsize=6.5)
    ax.set_ylabel("malignant pseudobulk", fontsize=6.5)
    ax.set_title("F. Malignant vs.\nnon-malignant expression", fontsize=8.5)
    ax.legend(fontsize=6, frameon=False)

    per_cell_gene = per_cell[["umap_1", "umap_2", "condition", "malignancy_status", f"{gene}_log_norm"]]
    vmax = per_cell_gene[f"{gene}_log_norm"].quantile(0.99)
    # G/H: UMAP panels
    ax_g = fig.add_subplot(gs[2, 0])
    ax_h = fig.add_subplot(gs[2, 1])
    for ax, cond, label in [(ax_g, "Control", "G. UMAP: Control"), (ax_h, "Tamoxifen", "H. UMAP: Tamoxifen (12h)")]:
        ax.scatter(per_cell_gene["umap_1"], per_cell_gene["umap_2"], s=1.5, c="#E8E8E8", linewidths=0)
        grp = per_cell_gene.loc[per_cell_gene["condition"] == cond]
        order = grp[f"{gene}_log_norm"].sort_values().index
        sc = ax.scatter(grp.loc[order, "umap_1"], grp.loc[order, "umap_2"], s=2.5, c=grp.loc[order, f"{gene}_log_norm"], cmap="viridis", vmin=0, vmax=vmax, linewidths=0)
        ax.set_title(label, fontsize=8.5)
        ax.set_xticks([])
        ax.set_yticks([])
    fig.colorbar(sc, ax=[ax_g, ax_h], label=f"{gene} log-norm", shrink=0.7)

    ax_text = fig.add_subplot(gs[2, 2:])
    ax_text.axis("off")
    f = frozen.set_index("gene").loc[gene]
    summary_text = (
        f"Frozen Track A (all epithelial): log2FC={f['gse245601_epi_log2fc']:.3f}, FDR={f['gse245601_epi_fdr']:.3g}\n"
        f"Frozen Track B (strict malignant): log2FC={f['gse245601_malignant_log2fc']:.3f}, FDR={f['gse245601_malignant_fdr']:.3g}\n\n"
        "Patient pseudobulk direction (all epithelial): see panel A / figure 05.\n"
        "Reliable malignant samples (>=50 cells both arms) = Tumor_02/03/07 only\n"
        "(dashed lines in panel B = <50 malignant cells in that arm -- descriptive,\n"
        "not part of the frozen inferential comparison).\n\n"
        "No individual cell is ever compared before/after treatment. Each line\n"
        "connects one tumor's own Control pseudobulk point to its own Tamoxifen\n"
        "pseudobulk point."
    )
    ax_text.text(0.0, 1.0, summary_text, va="top", ha="left", fontsize=8.5, family="monospace", transform=ax_text.transAxes)

    fig.suptitle("USP34 acute 12h decomposition: what underlies the apparent pseudobulk change?", fontsize=13, y=1.0)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info("wrote %s", out_path)


def run_visualization_2(config_path: str | Path = "config/config.yaml") -> None:
    config = _load_config(config_path)
    tables_dir = Path(config["gse245601_candidate_deepdive"]["output"]["tables_dir"])
    figures_dir = Path(config["gse245601_candidate_deepdive"]["output"]["figures_dir"])
    figures_dir.mkdir(parents=True, exist_ok=True)

    per_cell = load_per_cell_table(config)

    for i, gene in enumerate(GENES, start=18):
        plot_cell_distribution(per_cell, gene, figures_dir / f"{i:02d}_{gene}_cell_distribution.png")

    all_epi_pb = pd.read_csv(tables_dir / "patient_all_epithelial_pseudobulk.tsv", sep="\t")
    malignant_pb = pd.read_csv(tables_dir / "patient_malignant_pseudobulk.tsv", sep="\t")
    malig_cond = pd.read_csv(tables_dir / "malignancy_condition_patient_summary.tsv", sep="\t")
    frozen = pd.read_csv("results/tables/evidence_freeze/final_candidate_evidence.tsv", sep="\t")
    plot_usp34_decomposition(per_cell, all_epi_pb, malignant_pb, malig_cond, frozen, figures_dir / "22_USP34_decomposition.png")

    for i, gene in enumerate(GENES, start=23):
        plot_umap_expression(per_cell, gene, figures_dir / f"{i:02d}_{gene}_umap_control_tam.png")


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run_visualization_2()
