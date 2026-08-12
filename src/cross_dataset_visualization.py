"""Cross-dataset genome-wide integration, Phase 22: the 12-figure final
review set. Every figure is built from the already-committed tables
produced by the earlier phases of this module family -- no hand-typed
numbers, no figure built before its underlying table.

Data source: `results/tables/cross_dataset_genomewide/*.tsv`.
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

PCT_COLS = ["crispr_evidence_percentile", "gse118713_evidence_percentile", "gse245601_evidence_percentile", "gse240112_evidence_percentile", "gse111151_evidence_percentile"]
PCT_LABELS = ["CRISPR", "GSE118713", "GSE245601", "GSE240112", "GSE111151"]


def _load_config(config_path: str | Path) -> dict:
    with open(config_path) as f:
        return yaml.safe_load(f)


def _heatmap(matrix: pd.DataFrame, out_path: Path, title: str, cmap: str = "viridis", vmin: float = 0, vmax: float = 1, cbar_label: str = "evidence percentile") -> None:
    fig, ax = plt.subplots(figsize=(7, 0.35 * len(matrix) + 1.5))
    im = ax.imshow(matrix.to_numpy(dtype=float), cmap=cmap, vmin=vmin, vmax=vmax, aspect="auto")
    ax.set_xticks(range(len(matrix.columns)))
    ax.set_xticklabels(matrix.columns, rotation=45, ha="right", fontsize=8)
    ax.set_yticks(range(len(matrix.index)))
    ax.set_yticklabels(matrix.index, fontsize=8)
    fig.colorbar(im, ax=ax, label=cbar_label)
    ax.set_title(title, fontsize=10)
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info("wrote %s", out_path)


def plot_all_gene_evidence_heatmap(ranked: pd.DataFrame, out_path: Path, n_genes: int = 100) -> None:
    top = ranked.head(n_genes).set_index("gene")
    matrix = top[PCT_COLS]
    matrix.columns = PCT_LABELS
    _heatmap(matrix, out_path, f"Cross-dataset evidence percentile, top {n_genes} genes by global rank")


def plot_top20_heatmap(top20: pd.DataFrame, out_path: Path) -> None:
    m = top20.set_index("gene")[["crispr_percentile", "gse118713_percentile", "gse245601_percentile", "gse240112_percentile", "gse111151_percentile"]]
    m.columns = PCT_LABELS
    _heatmap(m, out_path, "Global Top 20: evidence percentile heatmap")


def plot_top20_evidence_matrix(top20: pd.DataFrame, out_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(9, 0.4 * len(top20) + 1.5))
    y = np.arange(len(top20))
    ax.barh(y, top20["datasets_fdr05"], color="#4C72B0", label="datasets FDR<0.05")
    ax.barh(y, top20["datasets_top10pct"], color="#C44E52", alpha=0.5, label="datasets top-10%", left=0)
    ax.set_yticks(y)
    ax.set_yticklabels(top20["gene"])
    ax.invert_yaxis()
    ax.set_xlabel("count (of 5 independent datasets)")
    ax.legend(fontsize=8)
    ax.set_title("Global Top 20: significance and top-10% dataset counts", fontsize=10)
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info("wrote %s", out_path)


def plot_resistance_consensus_heatmap(top20_resistance: pd.DataFrame, out_path: Path) -> None:
    m = top20_resistance.set_index("gene")[["resistance_up_count", "resistance_down_count", "resistance_fdr05_count"]]
    _heatmap(m, out_path, "Resistance-state RNA consensus Top 20 (GSE118713+GSE240112+GSE111151)", cmap="RdBu_r", vmin=0, vmax=3, cbar_label="count (of 3 datasets)")


def plot_rna_only_top20(top20_rna: pd.DataFrame, out_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(7, 0.4 * len(top20_rna) + 1.5))
    order = top20_rna.sort_values("rank")
    ax.barh(np.arange(len(order)), order["median_evidence_percentile"], color="#55A868")
    ax.set_yticks(np.arange(len(order)))
    ax.set_yticklabels(order["gene"])
    ax.invert_yaxis()
    ax.set_xlabel("median evidence percentile (RNA datasets only, CRISPR excluded)")
    ax.set_title("RNA-only Top 20 (GSE118713, GSE245601, GSE240112, GSE111151)", fontsize=10)
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info("wrote %s", out_path)


def plot_crispr_vs_resistance_scatter(ranked: pd.DataFrame, resistance_consensus: pd.DataFrame, out_path: Path) -> None:
    merged = ranked.merge(resistance_consensus[["gene", "resistance_median_percentile"]], on="gene", how="left")
    fig, ax = plt.subplots(figsize=(6.5, 6))
    ax.scatter(merged["crispr_evidence_percentile"], merged["resistance_median_percentile"], s=6, alpha=0.25, color="gray", linewidths=0)
    top = merged.head(20)
    ax.scatter(top["crispr_evidence_percentile"], top["resistance_median_percentile"], s=40, color="#C44E52", zorder=3)
    for _, row in top.iterrows():
        ax.annotate(row["gene"], (row["crispr_evidence_percentile"], row["resistance_median_percentile"]), fontsize=7, xytext=(3, 3), textcoords="offset points")
    ax.set_xlabel("CRISPR evidence percentile")
    ax.set_ylabel("resistance-RNA median evidence percentile (GSE118713+GSE240112+GSE111151)")
    ax.set_title("CRISPR functional evidence vs. resistance-RNA evidence (global Top 20 highlighted)", fontsize=10)
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info("wrote %s", out_path)


def plot_ranking_stability(stability: pd.DataFrame, top20_genes: list[str], out_path: Path) -> None:
    sub = stability.loc[stability["gene"].isin(top20_genes)].set_index("gene").loc[top20_genes]
    fig, ax = plt.subplots(figsize=(7, 0.4 * len(sub) + 1.5))
    y = np.arange(len(sub))
    ax.hlines(y, sub["best_rank"], sub["worst_rank"], color="lightgray", linewidth=4)
    ax.scatter(sub["rank_main"], y, color="#4C72B0", s=50, zorder=3, label="main global rank")
    ax.scatter(sub["median_rank"], y, color="#C44E52", marker="x", s=40, zorder=3, label="median alternate rank")
    ax.set_yticks(y)
    ax.set_yticklabels(sub.index)
    ax.invert_yaxis()
    ax.set_xscale("log")
    ax.set_xlabel("rank across all ranking-scheme variants (log scale; best=1)")
    ax.legend(fontsize=8)
    ax.set_title("Ranking stability: best-worst rank range, global Top 20", fontsize=10)
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info("wrote %s", out_path)


def plot_leave_one_out(stability: pd.DataFrame, top20_genes: list[str], out_path: Path) -> None:
    cols = ["rank_without_crispr", "rank_without_gse118713", "rank_without_gse245601", "rank_without_gse240112", "rank_without_gse111151"]
    labels = ["w/o CRISPR", "w/o GSE118713", "w/o GSE245601", "w/o GSE240112", "w/o GSE111151"]
    sub = stability.loc[stability["gene"].isin(top20_genes)].set_index("gene").loc[top20_genes]
    fig, ax = plt.subplots(figsize=(8, 5.5))
    for gene in sub.index:
        ax.plot(range(len(cols)), sub.loc[gene, cols], marker="o", markersize=4, linewidth=1, label=gene, alpha=0.8)
    ax.axhline(20, color="black", linestyle="--", linewidth=0.8, label="Top-20 cutoff")
    ax.set_xticks(range(len(cols)))
    ax.set_xticklabels(labels, rotation=30, ha="right")
    ax.set_yscale("log")
    ax.set_ylabel("global rank when the named dataset is removed (log scale)")
    ax.set_title("Leave-one-dataset-out ranking, global Top 20", fontsize=10)
    ax.legend(fontsize=6, ncol=2, loc="upper left", bbox_to_anchor=(1.0, 1.0))
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info("wrote %s", out_path)


def plot_evidence_category_counts(categories: pd.DataFrame, out_path: Path) -> None:
    counts = categories["evidence_category"].value_counts()
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.bar(counts.index, counts.values, color="#4C72B0")
    ax.set_yscale("log")
    ax.set_ylabel("number of genes (log scale)")
    ax.tick_params(axis="x", rotation=30)
    ax.set_title("Multimodal evidence-category counts, full gene universe", fontsize=10)
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info("wrote %s", out_path)


def plot_directional_pattern_matrix(patterns: pd.DataFrame, top20_genes: list[str], out_path: Path) -> None:
    arrow_map = {"↑": 1, "↓": -1, "=": 0, "·": np.nan}
    sub = patterns.loc[patterns["gene"].isin(top20_genes)].set_index("gene").loc[top20_genes]
    m = sub[["gse118713_arrow", "gse240112_arrow", "gse111151_arrow", "gse245601_arrow"]].apply(lambda col: col.map(arrow_map))
    m.columns = ["GSE118713", "GSE240112", "GSE111151", "GSE245601 (acute)"]
    fig, ax = plt.subplots(figsize=(6, 0.4 * len(m) + 1.5))
    im = ax.imshow(m.to_numpy(dtype=float), cmap="RdBu_r", vmin=-1, vmax=1, aspect="auto")
    ax.set_xticks(range(len(m.columns)))
    ax.set_xticklabels(m.columns, rotation=30, ha="right", fontsize=8)
    ax.set_yticks(range(len(m.index)))
    ax.set_yticklabels(m.index, fontsize=8)
    fig.colorbar(im, ax=ax, label="direction (blue=down, red=up)")
    ax.set_title("Directional pattern, global Top 20 (CRISPR direction shown separately in the table)", fontsize=9)
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info("wrote %s", out_path)


def plot_coverage_vs_evidence(full: pd.DataFrame, out_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(6.5, 5.5))
    tier_order = ["E", "D", "C", "B", "A"]
    positions = {t: i for i, t in enumerate(tier_order)}
    x = full["coverage_tier"].map(positions) + np.random.default_rng(0).uniform(-0.3, 0.3, len(full))
    ax.scatter(x, full["equal_dataset_mean_percentile"], s=3, alpha=0.15, color="gray", linewidths=0)
    ax.set_xticks(range(len(tier_order)))
    ax.set_xticklabels([f"Tier {t}\n({i + 1}/5)" for i, t in enumerate(tier_order)])
    ax.set_xlabel("coverage tier (E=1/5 datasets testable ... A=5/5)")
    ax.set_ylabel("equal_dataset_mean_percentile")
    ax.set_title("Coverage vs. evidence strength, full gene universe", fontsize=10)
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info("wrote %s", out_path)


def plot_surprise_candidates(ranked: pd.DataFrame, old_28_genes: set[str], out_path: Path, n: int = 20) -> None:
    surprises = ranked.loc[~ranked["gene"].isin(old_28_genes)].head(n)
    fig, ax = plt.subplots(figsize=(7, 0.4 * len(surprises) + 1.5))
    y = np.arange(len(surprises))
    ax.barh(y, surprises["median_evidence_percentile"], color="#DD8452")
    ax.set_yticks(y)
    ax.set_yticklabels(surprises["gene"])
    ax.invert_yaxis()
    ax.set_xlabel("median evidence percentile")
    ax.set_title(f"Surprise candidates: top {n} global-ranked genes outside the original 28 CRISPR-significant hits", fontsize=9)
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info("wrote %s", out_path)


def run_visualization(config_path: str | Path = "config/config.yaml") -> None:
    config = _load_config(config_path)
    cfg = config["cross_dataset_genomewide"]
    out = cfg["output"]
    tables_dir = Path(out["wide_matrix_tsv"]).parent
    final_review = Path(out["final_review_dir"])

    ranked = pd.read_csv(tables_dir / "global_ranking_eligible.tsv", sep="\t")
    full = pd.read_csv(tables_dir / "all_genes_cross_dataset_evidence_with_ranking.tsv", sep="\t")
    top20 = pd.read_csv(out["top20_global_tsv"], sep="\t")
    top20_resistance = pd.read_csv(out["top20_resistance_consensus_tsv"], sep="\t")
    top20_rna = pd.read_csv(out["top20_rna_only_tsv"], sep="\t")
    resistance_consensus = pd.read_csv(out["resistance_consensus_tsv"], sep="\t")
    stability = pd.read_csv(out["ranking_stability_tsv"], sep="\t")
    categories = pd.read_csv(tables_dir / "evidence_categories.tsv", sep="\t")
    patterns = pd.read_csv(out["directional_patterns_tsv"], sep="\t")

    top20_genes = top20["gene"].tolist()

    plot_all_gene_evidence_heatmap(ranked, final_review / "01_all_gene_evidence_heatmap_top100.png")
    plot_top20_heatmap(top20, final_review / "02_top20_global_heatmap.png")
    plot_top20_evidence_matrix(top20, final_review / "03_top20_global_evidence_matrix.png")
    plot_resistance_consensus_heatmap(top20_resistance, final_review / "04_top20_resistance_consensus_heatmap.png")
    plot_rna_only_top20(top20_rna, final_review / "05_top20_rna_only.png")
    plot_crispr_vs_resistance_scatter(ranked, resistance_consensus, final_review / "06_crispr_vs_resistance_scatter.png")
    plot_ranking_stability(stability, top20_genes, final_review / "07_ranking_stability.png")
    plot_leave_one_out(stability, top20_genes, final_review / "08_leave_one_dataset_out.png")
    plot_evidence_category_counts(categories, final_review / "09_evidence_category_counts.png")
    plot_directional_pattern_matrix(patterns, top20_genes, final_review / "10_directional_pattern_matrix.png")
    plot_coverage_vs_evidence(full, final_review / "11_coverage_vs_evidence.png")

    # old 28 CRISPR-significant hits (Gate-1 FDR<0.1) -- read from the already-frozen master table, not redefined here
    old28 = pd.read_csv("results/tables/crispr_gse118713_master_table.tsv", sep="\t")
    old28_genes = set(old28["gene_symbol"]) if "gene_symbol" in old28.columns else set(old28.iloc[:, 0])
    plot_surprise_candidates(ranked, old28_genes, final_review / "12_surprise_candidates.png")

    logger.info("run_visualization: wrote 12 figures to %s", final_review)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run_visualization()
