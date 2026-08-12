"""GSE111151 Phase 8 compact candidate figure set: candidate effect-size
plot (all 13, frozen-list order, never sorted by p-value), candidate
pseudobulk-style expression heatmap across all 11 samples, per-gene
sample-level plots for USP34/VEZF1/SUPT4H1 (and any other candidate that
stands out), and an optional volcano plot with the 13 candidates
highlighted. PCA/correlation figures are produced by
``src/gse111151_qc.py``, not duplicated here.

Sample-level plots connect each cell line's parental sample to its own
resistant subline(s) with a line -- unlike the GSE240112 primary/
recurrent comparison, this pairing is genuine (isogenic derivation from
the same parental culture, not independent patients), so a connecting
line is appropriate here (docs/GSE111151_PREANALYSIS.md section C).

Data source: GSE111151 (Hultsch et al., BMC Cancer 2018, PMID 30143015),
version as downloaded 2026-08-12.
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

from src.gse111151_qc import load_tmm_norm_factors  # noqa: E402

logger = logging.getLogger(__name__)

CELL_LINE_COLORS = {"MCF-7": "#4C72B0", "T-47D": "#DD8452", "ZR-75-1": "#55A868", "BT-474": "#C44E52"}
CELL_LINE_ORDER = ["MCF-7", "T-47D", "ZR-75-1", "BT-474"]


def _load_config(config_path: str | Path) -> dict:
    with open(config_path) as f:
        return yaml.safe_load(f)


def build_candidate_log2cpm_matrix(
    counts_path: str | Path, candidate_genes: list[str], candidate_ensembl_ids: dict[str, str], effective_lib_sizes: pd.Series | None = None
) -> pd.DataFrame:
    """genes x samples TMM-adjusted log2(CPM+1) matrix for the candidates
    present in the count matrix (by Ensembl ID). See
    ``src.gse111151_qc.load_tmm_norm_factors`` for why TMM adjustment
    matters here (naive CPM was found to be materially misleading for
    several samples)."""
    counts = pd.read_csv(counts_path, sep="\t").set_index("gene_id")
    counts = counts.drop(columns="gene_name")
    lib_sizes = effective_lib_sizes if effective_lib_sizes is not None else counts.sum(axis=0)
    present_ids = {g: candidate_ensembl_ids[g] for g in candidate_genes if candidate_ensembl_ids[g] in counts.index}
    sub = counts.loc[list(present_ids.values())]
    cpm = sub.div(lib_sizes, axis=1) * 1e6
    log2cpm = np.log2(cpm + 1)
    log2cpm.index = list(present_ids.keys())
    return log2cpm


def plot_candidate_heatmap(log2cpm: pd.DataFrame, sample_order: list[str], out_path: str | Path) -> None:
    mat = log2cpm[sample_order]
    zscored = mat.sub(mat.mean(axis=1), axis=0).div(mat.std(axis=1, ddof=1).replace(0, np.nan), axis=0)
    fig, ax = plt.subplots(figsize=(7, 0.4 * len(zscored) + 1.5))
    im = ax.imshow(zscored.to_numpy(), cmap="RdBu_r", vmin=-2, vmax=2, aspect="auto")
    ax.set_xticks(range(len(zscored.columns)))
    ax.set_xticklabels(zscored.columns, rotation=90, fontsize=8)
    ax.set_yticks(range(len(zscored.index)))
    ax.set_yticklabels(zscored.index)
    fig.colorbar(im, ax=ax, label="row z-score (log2 CPM+1)")
    ax.set_title("GSE111151 candidate expression across all 11 samples", fontsize=10)
    fig.tight_layout()
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info("wrote %s", out_path)


def plot_effect_size(candidate_table: pd.DataFrame, candidates: list[str], out_path: str | Path) -> None:
    fig, ax = plt.subplots(figsize=(8, 0.45 * len(candidates) + 1.5))
    y_positions = range(len(candidates))
    max_abs = 0.1
    for gene in candidates:
        row = candidate_table.loc[candidate_table["gene"] == gene]
        if len(row) > 0 and bool(row["tested"].iloc[0]):
            max_abs = max(max_abs, abs(float(row["log2fc"].iloc[0])))
    for y, gene in zip(y_positions, candidates):
        row = candidate_table.loc[candidate_table["gene"] == gene]
        if len(row) == 0 or not bool(row["tested"].iloc[0]):
            ax.text(0, y, "untested", va="center", ha="center", fontsize=8, color="gray")
            continue
        log2fc = float(row["log2fc"].iloc[0])
        fdr = float(row["candidate_set_bh_fdr"].iloc[0])
        color = "#C44E52" if fdr < 0.05 else "#4C72B0"
        ax.barh(y, log2fc, color=color)
        marker = "*" if fdr < 0.05 else ""
        offset = 0.06 * max_abs
        ax.text(log2fc + (offset if log2fc >= 0 else -offset), y, f"FDR={fdr:.3f}{marker}", va="center", ha="left" if log2fc >= 0 else "right", fontsize=7)
    ax.axvline(0, color="black", linewidth=0.8)
    ax.set_xlim(-max_abs * 1.7, max_abs * 1.7)
    ax.set_yticks(list(y_positions))
    ax.set_yticklabels(candidates)
    ax.invert_yaxis()
    ax.set_xlabel("log2FC (resistant vs parental, cell-line-blocked)")
    ax.set_title("Frozen 13 candidates: GSE111151 effect size (red = candidate-set BH FDR<0.05)", fontsize=9)
    fig.tight_layout()
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info("wrote %s", out_path)


def plot_sample_level(sample_level: pd.DataFrame, gene: str, out_path: str | Path) -> None:
    """Per cell line: parental point connected by a line to each of its
    own resistant subline(s) -- a genuine isogenic pairing, not an
    invented one."""
    sub = sample_level.loc[sample_level["gene"] == gene]
    fig, ax = plt.subplots(figsize=(6.5, 5))
    for i, cell_line in enumerate(CELL_LINE_ORDER):
        cl_sub = sub.loc[sub["cell_line"] == cell_line]
        parental = cl_sub.loc[cl_sub["resistance_status"] == "parental"]
        resistant = cl_sub.loc[cl_sub["resistance_status"] == "resistant"]
        if len(parental) == 0:
            continue
        p_val = parental["log2cpm"].iloc[0]
        color = CELL_LINE_COLORS[cell_line]
        ax.scatter([i], [p_val], color=color, marker="o", s=100, zorder=3, label=cell_line if i == 0 else None)
        for _, r_row in resistant.iterrows():
            ax.plot([i, i], [p_val, r_row["log2cpm"]], color=color, linestyle="-", linewidth=1.2, zorder=2)
            ax.scatter([i], [r_row["log2cpm"]], color=color, marker="^", s=100, zorder=3)
    ax.set_xticks(range(len(CELL_LINE_ORDER)))
    ax.set_xticklabels(CELL_LINE_ORDER)
    ax.set_ylabel(f"{gene} log2(CPM+1)")
    ax.set_title(f"{gene}: parental (circle) vs resistant subline(s) (triangle), by cell line", fontsize=10)
    fig.tight_layout()
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info("wrote %s", out_path)


def plot_volcano(de_df: pd.DataFrame, candidate_table: pd.DataFrame, out_path: str | Path) -> None:
    fig, ax = plt.subplots(figsize=(7, 6))
    ax.scatter(de_df["log2fc"], -np.log10(de_df["p_value"]), s=3, color="lightgray", linewidths=0)
    tested = candidate_table.loc[candidate_table["tested"]]
    ax.scatter(tested["log2fc"], -np.log10(tested["p_value"]), s=40, color="#C44E52", zorder=3)
    for _, row in tested.iterrows():
        ax.annotate(row["gene"], (row["log2fc"], -np.log10(row["p_value"])), fontsize=8, xytext=(4, 4), textcoords="offset points")
    ax.axhline(-np.log10(0.05), color="gray", linestyle="--", linewidth=0.8)
    ax.set_xlabel("log2FC (resistant vs parental)")
    ax.set_ylabel("-log10(nominal p-value)")
    ax.set_title("GSE111151 genome-wide volcano, frozen 13 candidates highlighted", fontsize=10)
    fig.tight_layout()
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info("wrote %s", out_path)


def run_candidate_visualization(config_path: str | Path = "config/config.yaml") -> None:
    config = _load_config(config_path)
    cfg = config["gse111151"]
    candidates = cfg["candidates"]["thirteen"]
    ensembl_ids = cfg["candidate_ensembl_ids"]
    figures_dir = Path(cfg["output"]["candidate_viz"]["figures_dir"])

    metadata = pd.read_csv(cfg["output"]["metadata_tsv"], sep="\t")
    sample_order = metadata.sort_values(["cell_line", "resistance_status"])["sample_id"].tolist()
    effective_lib_sizes = load_tmm_norm_factors(cfg["output"]["de"]["tmm_norm_factors_tsv"])

    log2cpm = build_candidate_log2cpm_matrix(cfg["output"]["counts_tsv"], candidates, ensembl_ids, effective_lib_sizes)
    plot_candidate_heatmap(log2cpm, sample_order, figures_dir / "candidate_heatmap.png")

    candidate_table = pd.read_csv(cfg["output"]["candidate_table_tsv"], sep="\t")
    plot_effect_size(candidate_table, candidates, figures_dir / "candidate_effect_size.png")

    sample_level = pd.read_csv(cfg["output"]["sample_level_tsv"], sep="\t")
    for gene in ["USP34", "VEZF1", "SUPT4H1"]:
        plot_sample_level(sample_level, gene, figures_dir / f"{gene.lower()}_sample_level.png")

    de_df = pd.read_csv(cfg["output"]["de"]["genomewide_tsv"], sep="\t")
    plot_volcano(de_df, candidate_table, figures_dir / "volcano.png")

    logger.info("run_candidate_visualization: wrote figures to %s", figures_dir)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run_candidate_visualization()
