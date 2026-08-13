"""Systems-network phase 6: resistance pathway consensus heatmap."""

from __future__ import annotations

import logging
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml
from matplotlib.colors import TwoSlopeNorm

logger = logging.getLogger(__name__)

RESISTANCE_COLS = ["gse118713", "gse240112", "gse111151"]
COL_LABELS = ["GSE118713", "GSE240112", "GSE111151", "GSE245601\n(acute 12h)"]


def _load_config(config_path: str | Path) -> dict:
    with open(config_path) as f:
        return yaml.safe_load(f)


def select_top_consensus_pathways(consensus: pd.DataFrame, n_each_direction: int = 15) -> pd.DataFrame:
    strong = consensus.loc[consensus["consensus_category"] == "STRONG_CONSENSUS"]
    top_pos = strong.loc[strong["median_NES"] > 0].sort_values("median_NES", ascending=False).head(n_each_direction)
    top_neg = strong.loc[strong["median_NES"] < 0].sort_values("median_NES", ascending=True).head(n_each_direction)
    return pd.concat([top_pos, top_neg], ignore_index=True)


def build_consensus_matrix(top_pathways: pd.DataFrame, gsea_tables: dict[str, pd.DataFrame]) -> tuple[pd.DataFrame, pd.DataFrame, list[str]]:
    keys = list(zip(top_pathways["collection"], top_pathways["pathway"]))
    labels = [f"{c.upper()}: {p.replace('HALLMARK_', '').replace('GOBP_', '').replace('REACTOME_', '').replace('_', ' ').title()[:55]}" for c, p in keys]

    nes = pd.DataFrame(np.nan, index=range(len(keys)), columns=RESISTANCE_COLS + ["gse245601"], dtype=float)
    fdr_sig = pd.DataFrame(False, index=range(len(keys)), columns=RESISTANCE_COLS + ["gse245601"], dtype=bool)
    for label in RESISTANCE_COLS + ["gse245601"]:
        table = gsea_tables[label].set_index(["collection", "pathway"])
        for i, key in enumerate(keys):
            if key in table.index:
                row = table.loc[key]
                nes.loc[i, label] = row["NES"]
                fdr_sig.loc[i, label] = bool(row["fdr"] < 0.05)
    return nes, fdr_sig, labels


def plot_resistance_consensus_heatmap(nes: pd.DataFrame, fdr_sig: pd.DataFrame, labels: list[str], out_path: Path) -> None:
    n_rows = len(labels)
    fig, ax = plt.subplots(figsize=(9, max(6, n_rows * 0.28)))

    # insert a visual gap column between resistance (3 cols) and acute (1 col)
    plot_matrix = np.full((n_rows, 5), np.nan)
    plot_matrix[:, 0:3] = nes[RESISTANCE_COLS].values
    plot_matrix[:, 4] = nes["gse245601"].values

    vmax = np.nanmax(np.abs(plot_matrix))
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)
    masked = np.ma.masked_invalid(plot_matrix)
    cmap = plt.cm.RdBu_r.copy()
    cmap.set_bad(color="#e5e5e5")

    im = ax.imshow(masked, cmap=cmap, norm=norm, aspect="auto")

    for i in range(n_rows):
        for j, label in enumerate(RESISTANCE_COLS):
            if fdr_sig.iloc[i][label]:
                ax.text(j, i, "*", ha="center", va="center", color="black", fontsize=11, fontweight="bold")
        if fdr_sig.iloc[i]["gse245601"]:
            ax.text(4, i, "*", ha="center", va="center", color="black", fontsize=11, fontweight="bold")

    ax.set_xticks([0, 1, 2, 4])
    ax.set_xticklabels(COL_LABELS, fontsize=9)
    ax.axvline(3, color="white", linewidth=6)
    ax.axvline(3, color="black", linewidth=0.8, linestyle="--")
    ax.set_yticks(range(n_rows))
    ax.set_yticklabels(labels, fontsize=7)
    ax.set_title("Resistance pathway consensus (GSE118713/GSE240112/GSE111151)\nvs. GSE245601 acute 12h tamoxifen response (shown separately, not part of consensus)", fontsize=10)
    cbar = fig.colorbar(im, ax=ax, shrink=0.6, pad=0.02)
    cbar.set_label("NES (preranked GSEA)", fontsize=9)
    ax.text(0.01, -0.06, "* FDR < 0.05    gray = pathway not tested in that dataset (below min gene-set size after overlap)", transform=ax.transAxes, fontsize=7.5)
    fig.tight_layout()
    fig.savefig(out_path, dpi=200)
    plt.close(fig)
    logger.info("plot_resistance_consensus_heatmap: saved %s (%d rows)", out_path, n_rows)


def run_visualization_1(config_path: str | Path = "config/config.yaml") -> None:
    config = _load_config(config_path)
    cfg = config["systems_network"]
    tables_dir = Path(cfg["output"]["tables_dir"])
    final_review_dir = Path(cfg["output"]["final_review_dir"])
    final_review_dir.mkdir(parents=True, exist_ok=True)

    consensus = pd.read_csv(tables_dir / "resistance_pathway_consensus.tsv", sep="\t")
    gsea_tables = {
        "gse118713": pd.read_csv(tables_dir / "gsea_gse118713.tsv", sep="\t"),
        "gse240112": pd.read_csv(tables_dir / "gsea_gse240112.tsv", sep="\t"),
        "gse111151": pd.read_csv(tables_dir / "gsea_gse111151.tsv", sep="\t"),
        "gse245601": pd.read_csv(tables_dir / "gsea_gse245601.tsv", sep="\t"),
    }

    top_pathways = select_top_consensus_pathways(consensus)
    nes, fdr_sig, labels = build_consensus_matrix(top_pathways, gsea_tables)
    plot_resistance_consensus_heatmap(nes, fdr_sig, labels, final_review_dir / "01_resistance_pathway_consensus.png")


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run_visualization_1()
