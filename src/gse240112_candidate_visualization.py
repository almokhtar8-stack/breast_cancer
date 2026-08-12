"""GSE240112 Phase 13 compact 13-candidate overview: dot plot (PT/RT
average expression + % detected), pseudobulk heatmap (genes x biological
samples), effect-size plot (log2FC with candidate-set BH status), and an
individual-sample heatmap to check whether an effect is consistent across
samples or driven by one. Genes are never ranked purely by p-value here --
the effect-size plot orders by the frozen candidate list, not by
significance.

Data source: GSE240112 (Fang et al., Genome Medicine 2024, PMID 39558215),
tumor-cell pseudobulk RT vs PT and per-cell candidate expression, version
as downloaded 2026-08-12.
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

SAMPLE_ORDER = ["PT1", "PT2", "PT3", "RT1", "RT2", "RT3"]


def _load_config(config_path: str | Path) -> dict:
    with open(config_path) as f:
        return yaml.safe_load(f)


def build_dotplot_data(raw_counts_tsv: str | Path, metadata_tsv: str | Path, candidates: list[str]) -> pd.DataFrame:
    """One row per (gene, group): mean raw count and % of tumor cells
    detecting the gene. Genes absent from the raw-counts table (untested
    candidates) are skipped -- the caller reports their untested status
    separately (candidate_table.tsv), this function never invents a value."""
    raw = pd.read_csv(raw_counts_tsv, sep="\t")
    meta = pd.read_csv(metadata_tsv, sep="\t")[["cell_id", "group"]]
    merged = raw.merge(meta, on="cell_id", how="inner", validate="one_to_one")
    if len(merged) != len(raw):
        raise ValueError("join between candidate raw-counts table and metadata lost or duplicated rows")

    rows = []
    for gene in candidates:
        if gene not in raw.columns:
            continue
        for group in ["PT", "RT"]:
            sub = merged.loc[merged["group"] == group, gene]
            rows.append(
                {
                    "gene": gene,
                    "group": group,
                    "mean_raw_count": float(sub.mean()),
                    "pct_detected": 100.0 * (sub > 0).mean(),
                }
            )
    out = pd.DataFrame(rows)
    logger.info("build_dotplot_data: %d gene x group rows", len(out))
    return out


def plot_dotplot(dotplot_df: pd.DataFrame, candidates: list[str], out_path: str | Path) -> None:
    present_genes = [g for g in candidates if g in set(dotplot_df["gene"])]
    fig, ax = plt.subplots(figsize=(6, 0.4 * len(present_genes) + 1.5))
    groups = ["PT", "RT"]
    for gi, gene in enumerate(present_genes):
        for xi, group in enumerate(groups):
            row = dotplot_df.loc[(dotplot_df["gene"] == gene) & (dotplot_df["group"] == group)]
            if len(row) == 0:
                continue
            size = max(row["pct_detected"].iloc[0] * 4, 5)
            color = np.log1p(row["mean_raw_count"].iloc[0])
            ax.scatter(xi, gi, s=size, c=[color], cmap="viridis", vmin=0, vmax=np.log1p(dotplot_df["mean_raw_count"].max()))
    ax.set_xticks(range(len(groups)))
    ax.set_xticklabels(groups)
    ax.set_yticks(range(len(present_genes)))
    ax.set_yticklabels(present_genes)
    ax.set_title("Candidate detection: dot size = % tumor cells detected,\ncolor = log1p(mean raw count)", fontsize=9)
    ax.set_xlim(-0.5, 1.5)
    fig.tight_layout()
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info("wrote %s", out_path)


def build_pseudobulk_log2cpm(counts_tsv: str | Path, candidates: list[str]) -> pd.DataFrame:
    """genes x samples log2(CPM+1) matrix for the candidates present in
    the pseudobulk count table."""
    counts = pd.read_csv(counts_tsv, sep="\t").set_index("gene")
    lib_sizes = counts.sum(axis=0)
    present = [g for g in candidates if g in counts.index]
    cpm = counts.loc[present].div(lib_sizes, axis=1) * 1e6
    log2cpm = np.log2(cpm + 1)
    return log2cpm[SAMPLE_ORDER]


def plot_pseudobulk_heatmap(log2cpm: pd.DataFrame, out_path: str | Path) -> None:
    zscored = log2cpm.sub(log2cpm.mean(axis=1), axis=0).div(log2cpm.std(axis=1, ddof=1).replace(0, np.nan), axis=0)
    fig, ax = plt.subplots(figsize=(6, 0.4 * len(zscored) + 1.5))
    im = ax.imshow(zscored.to_numpy(), cmap="RdBu_r", vmin=-2, vmax=2, aspect="auto")
    ax.set_xticks(range(len(zscored.columns)))
    ax.set_xticklabels(zscored.columns, rotation=45)
    ax.set_yticks(range(len(zscored.index)))
    ax.set_yticklabels(zscored.index)
    fig.colorbar(im, ax=ax, label="row z-score (log2 CPM+1)")
    ax.set_title("Candidate pseudobulk expression (biological samples)", fontsize=10)
    fig.tight_layout()
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info("wrote %s", out_path)


def plot_effect_size(candidate_table: pd.DataFrame, candidates: list[str], out_path: str | Path) -> None:
    """log2FC (RT vs PT) for all 13 candidates, in frozen list order (not
    sorted by significance). Bar color/marker indicates candidate-set BH
    status; untested genes shown as an explicit gap with a text label."""
    fig, ax = plt.subplots(figsize=(8.5, 0.45 * len(candidates) + 1.5))
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
    ax.set_xlim(-max_abs * 1.6, max_abs * 1.6)
    ax.set_yticks(list(y_positions))
    ax.set_yticklabels(candidates)
    ax.invert_yaxis()
    ax.set_xlabel("log2FC (RT vs PT, tumor-cell pseudobulk)")
    ax.set_title("Frozen 13 candidates: effect size (red = candidate-set BH FDR<0.05)", fontsize=9)
    fig.tight_layout()
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info("wrote %s", out_path)


def plot_sample_level_heatmap(sample_level: pd.DataFrame, candidates: list[str], out_path: str | Path) -> None:
    """genes x individual samples, z-scored per gene, to check whether an
    effect is consistent across all samples in a group or driven by one."""
    present = [g for g in candidates if g in set(sample_level["gene"])]
    pivot = sample_level.loc[sample_level["gene"].isin(present)].pivot(index="gene", columns="sample_id", values="log2cpm")
    pivot = pivot.loc[present, SAMPLE_ORDER]
    zscored = pivot.sub(pivot.mean(axis=1), axis=0).div(pivot.std(axis=1, ddof=1).replace(0, np.nan), axis=0)

    fig, ax = plt.subplots(figsize=(6, 0.4 * len(zscored) + 1.5))
    im = ax.imshow(zscored.to_numpy(), cmap="RdBu_r", vmin=-2, vmax=2, aspect="auto")
    ax.set_xticks(range(len(zscored.columns)))
    ax.set_xticklabels(zscored.columns, rotation=45)
    ax.set_yticks(range(len(zscored.index)))
    ax.set_yticklabels(zscored.index)
    fig.colorbar(im, ax=ax, label="row z-score (log2 CPM+1)")
    ax.set_title("Candidate expression by individual sample", fontsize=10)
    fig.tight_layout()
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info("wrote %s", out_path)


def run_candidate_visualization(config_path: str | Path = "config/config.yaml") -> None:
    config = _load_config(config_path)
    cfg = config["gse240112"]
    candidates = cfg["candidates"]["thirteen"]
    figures_dir = Path(cfg["output"]["candidate_viz"]["figures_dir"])

    dotplot_df = build_dotplot_data(cfg["output"]["tt_cancer_candidate_raw_tsv"], cfg["output"]["tt_cancer_metadata_tsv"], candidates)
    plot_dotplot(dotplot_df, candidates, figures_dir / "candidate_dotplot.png")

    log2cpm = build_pseudobulk_log2cpm(cfg["output"]["tumor_cell"]["counts_tsv"], candidates)
    plot_pseudobulk_heatmap(log2cpm, figures_dir / "candidate_pseudobulk_heatmap.png")

    candidate_table = pd.read_csv(cfg["output"]["candidate_table_tsv"], sep="\t")
    plot_effect_size(candidate_table, candidates, figures_dir / "candidate_effect_size.png")

    sample_level = pd.read_csv(Path(cfg["output"]["candidate_table_tsv"]).parent / "candidate_sample_level_log2cpm.tsv", sep="\t")
    plot_sample_level_heatmap(sample_level, candidates, figures_dir / "candidate_sample_level_heatmap.png")

    logger.info("run_candidate_visualization: wrote 4 figures to %s", figures_dir)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run_candidate_visualization()
