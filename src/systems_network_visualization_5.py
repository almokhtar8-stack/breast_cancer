"""Systems-network phase 25: candidate x pathway x CRISPR matrix figure.

Same pathway selection as Phase 8 (STRONG_CONSENSUS, Hallmark+Reactome only,
for readability) with two additional column groups: RNA NES per dataset
(three resistance + acute GSE245601) and CRISPR pathway-level GSEA NES/FDR.
"""

from __future__ import annotations

import logging
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml
from matplotlib.colors import TwoSlopeNorm

logger = logging.getLogger(__name__)

CANDIDATES = ["USP34", "VEZF1", "EML5", "CITED2"]
RNA_DATASETS = ["gse118713", "gse240112", "gse111151", "gse245601"]
RNA_LABELS = ["GSE118713", "GSE240112", "GSE111151", "GSE245601\n(acute)"]


def _load_config(config_path: str | Path) -> dict:
    with open(config_path) as f:
        return yaml.safe_load(f)


def select_pathways(membership: pd.DataFrame) -> list[tuple[str, str]]:
    sub = membership.loc[(membership["resistance_consensus_class"] == "STRONG_CONSENSUS") & (membership["collection"].isin(["hallmark", "reactome"]))]
    return sorted(sub[["collection", "pathway"]].drop_duplicates().itertuples(index=False, name=None))


def plot_matrix(pathway_keys: list[tuple[str, str]], membership: pd.DataFrame, gsea_tables: dict[str, pd.DataFrame], gsea_crispr: pd.DataFrame, out_path: Path) -> None:
    n_rows = len(pathway_keys)
    labels = [f"{c.upper()}: {p.replace('HALLMARK_', '').replace('REACTOME_', '').replace('_', ' ').title()}" for c, p in pathway_keys]

    membership_lookup = membership.set_index(["collection", "pathway", "candidate"])
    cand_matrix = np.zeros((n_rows, len(CANDIDATES)))
    for i, key in enumerate(pathway_keys):
        for j, cand in enumerate(CANDIDATES):
            idx = key + (cand,)
            if idx in membership_lookup.index:
                row = membership_lookup.loc[idx]
                if isinstance(row, pd.DataFrame):
                    row = row.iloc[0]
                if row["candidate_is_member"]:
                    cand_matrix[i, j] = 2
                elif isinstance(row["candidate_is_leading_edge_datasets"], str) and row["candidate_is_leading_edge_datasets"]:
                    cand_matrix[i, j] = 1

    rna_matrix = np.full((n_rows, len(RNA_DATASETS)), np.nan)
    for i, key in enumerate(pathway_keys):
        for j, d in enumerate(RNA_DATASETS):
            t = gsea_tables[d].set_index(["collection", "pathway"])
            if key in t.index:
                rna_matrix[i, j] = t.loc[key, "NES"]

    crispr_matrix = np.full((n_rows, 1), np.nan)
    crispr_fdr = np.full((n_rows, 1), np.nan)
    ct = gsea_crispr.set_index(["collection", "pathway"])
    for i, key in enumerate(pathway_keys):
        if key in ct.index:
            crispr_matrix[i, 0] = ct.loc[key, "NES"]
            crispr_fdr[i, 0] = ct.loc[key, "fdr"]

    fig, axes = plt.subplots(1, 3, figsize=(13, max(6, n_rows * 0.32)), gridspec_kw={"width_ratios": [len(CANDIDATES), len(RNA_DATASETS), 1.3], "wspace": 0.15})

    ax0 = axes[0]
    cmap0 = plt.matplotlib.colors.ListedColormap(["#f5f5f5", "#a6bddb", "#08519c"])
    ax0.imshow(cand_matrix, cmap=cmap0, vmin=0, vmax=2, aspect="auto")
    ax0.set_xticks(range(len(CANDIDATES)))
    ax0.set_xticklabels(CANDIDATES, fontsize=9, fontweight="bold")
    ax0.set_yticks(range(n_rows))
    ax0.set_yticklabels(labels, fontsize=7.5)
    ax0.set_title("Candidate\nmembership", fontsize=9)

    ax1 = axes[1]
    vmax = np.nanmax(np.abs(rna_matrix)) if np.any(~np.isnan(rna_matrix)) else 1
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)
    cmap1 = plt.cm.RdBu_r.copy()
    cmap1.set_bad("#e5e5e5")
    im1 = ax1.imshow(np.ma.masked_invalid(rna_matrix), cmap=cmap1, norm=norm, aspect="auto")
    ax1.set_xticks(range(len(RNA_DATASETS)))
    ax1.set_xticklabels(RNA_LABELS, fontsize=7.5)
    ax1.set_yticks([])
    ax1.set_title("RNA (NES)", fontsize=9)
    ax1.axvline(2.5, color="black", linewidth=0.8, linestyle="--")

    ax2 = axes[2]
    im2 = ax2.imshow(np.ma.masked_invalid(crispr_matrix), cmap=cmap1, norm=norm, aspect="auto")
    for i in range(n_rows):
        if crispr_fdr[i, 0] < 0.05:
            ax2.text(0, i, "*", ha="center", va="center", fontsize=11, fontweight="bold")
    ax2.set_xticks([0])
    ax2.set_xticklabels(["CRISPR"], fontsize=8)
    ax2.set_yticks([])
    ax2.set_title("Function\n(NES)", fontsize=9)

    legend_elems = [
        plt.Rectangle((0, 0), 1, 1, facecolor="#f5f5f5", edgecolor="black", label="no evidence"),
        plt.Rectangle((0, 0), 1, 1, facecolor="#a6bddb", edgecolor="black", label="leading-edge (>=1 dataset)"),
        plt.Rectangle((0, 0), 1, 1, facecolor="#08519c", edgecolor="black", label="direct curated member"),
    ]
    fig.legend(handles=legend_elems, loc="lower center", bbox_to_anchor=(0.22, -0.02), fontsize=7.5, ncol=1, frameon=False)
    cbar = fig.colorbar(im1, ax=[ax1, ax2], shrink=0.5, pad=0.02, location="right")
    cbar.set_label("NES", fontsize=8)
    fig.suptitle("Candidate x resistance-consensus pathway x CRISPR functional matrix\n(* = CRISPR pathway FDR<0.05; STRONG_CONSENSUS Hallmark+Reactome pathways only)", fontsize=10)
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    logger.info("plot_matrix: saved %s (%d pathways)", out_path, n_rows)


def run_visualization_5(config_path: str | Path = "config/config.yaml") -> None:
    config = _load_config(config_path)
    cfg = config["systems_network"]
    tables_dir = Path(cfg["output"]["tables_dir"])
    final_review_dir = Path(cfg["output"]["final_review_dir"])

    membership = pd.read_csv(tables_dir / "candidate_pathway_membership.tsv", sep="\t")
    gsea_tables = {d: pd.read_csv(tables_dir / f"gsea_{d}.tsv", sep="\t") for d in RNA_DATASETS}
    gsea_crispr = pd.read_csv(tables_dir / "gsea_crispr.tsv", sep="\t")

    pathway_keys = select_pathways(membership)
    plot_matrix(pathway_keys, membership, gsea_tables, gsea_crispr, final_review_dir / "06_candidate_pathway_crispr_matrix.png")


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run_visualization_5()
