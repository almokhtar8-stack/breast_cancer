"""GSE245601 candidate deep-dive Phases 13-14: epithelial cluster
-specific acute response, using the project's only existing clustering
(`seurat_clusters`, computed once on the full multi-lineage object;
"epithelial clusters" here means the subset of cluster IDs populated by
epithelial cells -- there is no separate epithelial-only subclustering
and none is computed here).

Minimum-representation rule (declared BEFORE inspecting any gene
-specific cluster result, based only on cell-count structure): a cluster
is "sufficiently represented" if at least 3 tumors have >=10 cells in
BOTH Control and Tamoxifen within that cluster. This is a conservative,
descriptive-support threshold, not a statistical power calculation --
clusters failing it are still reported (for transparency) but flagged
LOW_SUPPORT and excluded from any interpretive claim.
"""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from src.gse245601_candidate_deepdive_data import GENES, load_per_cell_table

logger = logging.getLogger(__name__)

MIN_CELLS_PER_PATIENT_CONDITION = 10
MIN_TUMORS_SUPPORTED = 3


def _load_config(config_path: str | Path) -> dict:
    with open(config_path) as f:
        return yaml.safe_load(f)


def determine_cluster_support(per_cell: pd.DataFrame) -> pd.DataFrame:
    tab = per_cell.groupby(["seurat_clusters", "patient", "condition"], observed=True).size().reset_index(name="n_cells")
    piv = tab.pivot_table(index=["seurat_clusters", "patient"], columns="condition", values="n_cells", fill_value=0)
    for c in ("Control", "Tamoxifen"):
        if c not in piv.columns:
            piv[c] = 0
    piv["min_cells_both_arms"] = piv[["Control", "Tamoxifen"]].min(axis=1)
    piv["patient_supported"] = piv["min_cells_both_arms"] >= MIN_CELLS_PER_PATIENT_CONDITION
    support = piv.reset_index().groupby("seurat_clusters")["patient_supported"].sum().rename("n_tumors_supported").reset_index()
    support["sufficiently_represented"] = support["n_tumors_supported"] >= MIN_TUMORS_SUPPORTED
    total_cells = per_cell.groupby("seurat_clusters", observed=True).size().rename("total_cells_in_cluster")
    support = support.merge(total_cells, on="seurat_clusters")
    logger.info("determine_cluster_support: %d/%d epithelial clusters sufficiently represented (rule: >=%d tumors with >=%d cells/arm)", int(support["sufficiently_represented"].sum()), len(support), MIN_TUMORS_SUPPORTED, MIN_CELLS_PER_PATIENT_CONDITION)
    return support


def build_cluster_candidate_response(per_cell: pd.DataFrame, genes: list[str], support: pd.DataFrame) -> pd.DataFrame:
    supported_clusters = set(support.loc[support["sufficiently_represented"], "seurat_clusters"])
    rows = []
    for cluster, cgrp in per_cell.groupby("seurat_clusters", observed=True):
        is_supported = cluster in supported_clusters
        for gene in genes:
            pb = cgrp.groupby(["patient", "condition"], observed=True).agg(raw=(f"{gene}_raw_count", "sum"), lib=("nCount_RNA", "sum")).reset_index()
            pb["norm"] = np.log2(pb["raw"] / pb["lib"].replace(0, np.nan) * 1e6 + 1)
            wide = pb.pivot_table(index="patient", columns="condition", values="norm")
            if "Control" not in wide.columns or "Tamoxifen" not in wide.columns:
                continue
            valid = wide.dropna(subset=["Control", "Tamoxifen"])
            delta = valid["Tamoxifen"] - valid["Control"]
            n_up = int((delta > 0).sum())
            n_down = int((delta < 0).sum())
            rows.append(
                {
                    "seurat_clusters": cluster, "gene": gene, "sufficiently_represented": is_supported,
                    "n_tumors_with_matched_pair": len(valid), "n_tumors_increase": n_up, "n_tumors_decrease": n_down,
                    "median_descriptive_log2fc": float(delta.median()) if len(delta) else float("nan"),
                    "mean_descriptive_log2fc": float(delta.mean()) if len(delta) else float("nan"),
                }
            )
    out = pd.DataFrame(rows)
    logger.info("build_cluster_candidate_response: %d cluster x gene rows", len(out))
    return out


def plot_cluster_heatmap(cluster_response: pd.DataFrame, out_path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    clusters = sorted(cluster_response["seurat_clusters"].unique(), key=int)
    genes = GENES
    grid = np.full((len(clusters), len(genes)), np.nan)
    supported = np.zeros((len(clusters), len(genes)), dtype=bool)
    for ci, cluster in enumerate(clusters):
        for gi, gene in enumerate(genes):
            row = cluster_response.loc[(cluster_response["seurat_clusters"] == cluster) & (cluster_response["gene"] == gene)]
            if len(row) == 0:
                continue
            grid[ci, gi] = row["median_descriptive_log2fc"].iloc[0]
            supported[ci, gi] = bool(row["sufficiently_represented"].iloc[0])

    fig, ax = plt.subplots(figsize=(6, 0.4 * len(clusters) + 2))
    vmax = np.nanmax(np.abs(grid))
    im = ax.imshow(grid, cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
    for ci in range(len(clusters)):
        for gi in range(len(genes)):
            v = grid[ci, gi]
            if np.isnan(v):
                continue
            marker = "" if supported[ci, gi] else " (low support)"
            ax.text(gi, ci, f"{v:.2f}{marker}", ha="center", va="center", fontsize=6 if supported[ci, gi] else 5, color="black" if supported[ci, gi] else "gray")
    ax.set_xticks(range(len(genes)))
    ax.set_xticklabels(genes)
    ax.set_yticks(range(len(clusters)))
    ax.set_yticklabels([f"cluster {c}" for c in clusters], fontsize=8)
    fig.colorbar(im, ax=ax, label="median descriptive log2FC across matched tumors")
    ax.set_title("Cluster-specific candidate response\n(gray/small text = below the >=3-tumor support threshold)", fontsize=9)
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info("wrote %s", out_path)


def plot_usp34_cluster_response(per_cell: pd.DataFrame, cluster_response: pd.DataFrame, out_path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from src.gse245601_candidate_deepdive_visualization_1 import PATIENT_COLORS, PATIENT_ORDER

    gene = "USP34"
    sub = cluster_response.loc[(cluster_response["gene"] == gene) & (cluster_response["sufficiently_represented"])].sort_values("median_descriptive_log2fc")
    clusters = sub["seurat_clusters"].tolist()
    if len(clusters) == 0:
        logger.warning("plot_usp34_cluster_response: no sufficiently-represented clusters for %s", gene)
        return

    fig, axes = plt.subplots(1, len(clusters), figsize=(3.2 * len(clusters), 4), sharey=True)
    if len(clusters) == 1:
        axes = [axes]
    for ax, cluster in zip(axes, clusters):
        cgrp = per_cell.loc[per_cell["seurat_clusters"] == cluster]
        pb = cgrp.groupby(["patient", "condition"], observed=True).agg(raw=(f"{gene}_raw_count", "sum"), lib=("nCount_RNA", "sum")).reset_index()
        pb["norm"] = np.log2(pb["raw"] / pb["lib"] * 1e6 + 1)
        wide = pb.pivot_table(index="patient", columns="condition", values="norm")
        for patient in PATIENT_ORDER:
            if patient not in wide.index or wide.loc[patient].isna().any():
                continue
            row = wide.loc[patient]
            ax.plot([0, 1], [row["Control"], row["Tamoxifen"]], color=PATIENT_COLORS[patient], marker="o", markersize=5, linewidth=1.3, alpha=0.85, label=patient)
        row_summary = sub.loc[sub["seurat_clusters"] == cluster].iloc[0]
        ax.set_title(f"cluster {cluster}\n{int(row_summary['n_tumors_increase'])}up/{int(row_summary['n_tumors_decrease'])}down of {int(row_summary['n_tumors_with_matched_pair'])}", fontsize=8.5)
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["Ctrl", "Tam"])
    axes[0].set_ylabel("log2(CPM+1)")
    axes[-1].legend(loc="center left", bbox_to_anchor=(1.02, 0.5), fontsize=6.5, frameon=False, title="tumor", title_fontsize=6.5)
    fig.suptitle(f"{gene}: patient-level response within sufficiently-represented epithelial clusters", fontsize=10.5, y=1.03)
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info("wrote %s", out_path)


def run_clusters(config_path: str | Path = "config/config.yaml") -> dict[str, pd.DataFrame]:
    config = _load_config(config_path)
    out_dir = Path(config["gse245601_candidate_deepdive"]["output"]["tables_dir"])
    figures_dir = Path(config["gse245601_candidate_deepdive"]["output"]["figures_dir"])
    out_dir.mkdir(parents=True, exist_ok=True)

    per_cell = load_per_cell_table(config)
    support = determine_cluster_support(per_cell)
    support.to_csv(out_dir / "cluster_support.tsv", sep="\t", index=False)

    cluster_response = build_cluster_candidate_response(per_cell, GENES, support)
    cluster_response.to_csv(out_dir / "cluster_candidate_response.tsv", sep="\t", index=False)

    plot_cluster_heatmap(cluster_response, figures_dir / "27_cluster_specific_candidate_response.png")
    plot_usp34_cluster_response(per_cell, cluster_response, figures_dir / "28_USP34_cluster_response.png")

    return {"support": support, "cluster_response": cluster_response}


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run_clusters()
