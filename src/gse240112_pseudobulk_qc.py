"""Pseudobulk QC for the GSE240112 tumor-cell primary-vs-recurrent
comparison -- library size, detected genes, sample-sample correlation,
and PCA on the pseudobulk count matrix written by
``scripts/analysis/gse240112_02_build_pseudobulk.R``.

Unlike the project's GSE245601 pseudobulk QC module, PT and RT samples
here are unpaired (different, unrelated patients -- see
docs/GSE240112_DATA_AUDIT.md section 2), so no patient-blocking or
paired-line visualization is used; samples are grouped and colored by
PT/RT only. No sample is excluded here regardless of outlier appearance
(docs/GSE240112_PREANALYSIS.md section G on RT3's low cell count) --
a catastrophic technical problem would be reported, not silently dropped.

Data source: GSE240112 (Fang et al., Genome Medicine 2024, PMID 39558215),
tumor-cell pseudobulk, version as downloaded 2026-08-12.
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
from scipy.cluster.hierarchy import dendrogram, linkage  # noqa: E402

logger = logging.getLogger(__name__)

GROUP_COLORS = {"PT": "#4C72B0", "RT": "#C44E52"}


def _load_config(config_path: str | Path) -> dict:
    with open(config_path) as f:
        return yaml.safe_load(f)


def load_pseudobulk(counts_path: str | Path, metadata_path: str | Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Load the tumor-cell pseudobulk count matrix (genes x samples,
    first column 'gene') and its metadata table. Read-only."""
    counts = pd.read_csv(counts_path, sep="\t").set_index("gene")
    metadata = pd.read_csv(metadata_path, sep="\t")
    missing = set(metadata["sample_id"]) - set(counts.columns)
    if missing:
        raise ValueError(f"metadata sample_id(s) not present as count matrix columns: {missing}")
    logger.info("load_pseudobulk: %d genes x %d samples", counts.shape[0], len(metadata))
    return counts, metadata


def compute_log2cpm(counts: pd.DataFrame, sample_ids: list[str]) -> pd.DataFrame:
    """log2(CPM+1), library size from each sample's own column sum."""
    lib_sizes = counts[sample_ids].sum(axis=0)
    cpm = counts[sample_ids].div(lib_sizes, axis=1) * 1e6
    return np.log2(cpm + 1)


def compute_sample_correlations(counts: pd.DataFrame, sample_ids: list[str]) -> pd.DataFrame:
    """Pairwise Spearman correlation of log2(CPM+1), long format."""
    log2cpm = compute_log2cpm(counts, sample_ids)
    corr = log2cpm.corr(method="spearman")
    corr.index.name = "sample_id_1"
    long = corr.reset_index().melt(id_vars="sample_id_1", var_name="sample_id_2", value_name="correlation")
    logger.info("compute_sample_correlations: %dx%d Spearman correlation matrix", len(sample_ids), len(sample_ids))
    return long


def select_top_variable_genes(counts: pd.DataFrame, sample_ids: list[str], n_genes: int) -> pd.DataFrame:
    """Top-N most variable genes by variance of log2(CPM+1)."""
    log2cpm = compute_log2cpm(counts, sample_ids)
    variance = log2cpm.var(axis=1, ddof=1)
    ordered = variance.sort_values(ascending=False, kind="stable")
    top_genes = ordered.index[:n_genes]
    logger.info("select_top_variable_genes: selected %d of %d genes", len(top_genes), len(counts))
    return counts.loc[top_genes]


def compute_pca(top_variable_counts: pd.DataFrame, metadata: pd.DataFrame, sample_ids: list[str]) -> pd.DataFrame:
    """PCA of pseudobulk samples. Genes centered before SVD; PC sign fixed
    deterministically. Long format: one row per (sample_id, pc)."""
    log2cpm = compute_log2cpm(top_variable_counts, sample_ids)
    matrix = log2cpm[sample_ids].to_numpy(dtype=float).T
    centered = matrix - matrix.mean(axis=0, keepdims=True)

    u, s, _vt = np.linalg.svd(centered, full_matrices=False)
    scores = u * s
    variance_explained = (s**2) / np.sum(s**2)

    for k in range(scores.shape[1]):
        col = scores[:, k]
        sign = 1.0 if col[np.argmax(np.abs(col))] >= 0 else -1.0
        scores[:, k] = col * sign

    meta_by_sample = metadata.set_index("sample_id")
    records = []
    for i, sample_id in enumerate(sample_ids):
        for k in range(scores.shape[1]):
            records.append(
                {
                    "sample_id": sample_id,
                    "group": meta_by_sample.loc[sample_id, "group"],
                    "pc": f"PC{k + 1}",
                    "coordinate": float(scores[i, k]),
                    "variance_explained_fraction": float(variance_explained[k]),
                }
            )
    logger.info("compute_pca: %d PCs for %d samples", scores.shape[1], len(sample_ids))
    return pd.DataFrame(records)


def build_qc_summary(metadata: pd.DataFrame) -> pd.DataFrame:
    """Per-sample QC summary -- passthrough of the pseudobulk construction
    metadata. No sample is flagged for exclusion here."""
    out = metadata.sort_values(["group", "sample_id"]).reset_index(drop=True)
    logger.info(
        "build_qc_summary: library size range %.2fM-%.2fM, detected genes range %d-%d",
        out["total_library_size"].min() / 1e6,
        out["total_library_size"].max() / 1e6,
        out["n_detected_genes"].min(),
        out["n_detected_genes"].max(),
    )
    return out


def plot_library_and_detected_genes(metadata: pd.DataFrame, out_path: str | Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(9, 4.5))
    order = metadata.sort_values(["group", "sample_id"])
    for ax, col, title in zip(axes, ("total_library_size", "n_detected_genes"), ("total pseudobulk library size", "detected genes")):
        colors = [GROUP_COLORS[g] for g in order["group"]]
        ax.bar(order["sample_id"], order[col], color=colors)
        ax.set_title(title, fontsize=10)
        ax.set_ylabel(col)
        ax.tick_params(axis="x", rotation=45)
    from matplotlib.patches import Patch

    handles = [Patch(color=c, label=g) for g, c in GROUP_COLORS.items()]
    axes[1].legend(handles=handles, fontsize=8)
    fig.suptitle("GSE240112 tumor-cell pseudobulk: library size / detected genes", fontsize=11)
    fig.tight_layout()
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info("wrote %s", out_path)


def plot_sample_correlation_heatmap(corr_long: pd.DataFrame, out_path: str | Path) -> None:
    corr = corr_long.pivot(index="sample_id_1", columns="sample_id_2", values="correlation")
    dist = 1 - corr
    link = linkage(dist.to_numpy()[np.triu_indices(len(dist), k=1)], method="average")

    fig, (ax_dendro, ax_heat) = plt.subplots(1, 2, figsize=(9, 4.5), gridspec_kw={"width_ratios": [1, 2]})
    dend = dendrogram(link, labels=corr.index.tolist(), orientation="left", ax=ax_dendro)
    ax_dendro.set_xticks([])

    order = dend["ivl"]
    ordered_corr = corr.loc[order, order]
    im = ax_heat.imshow(ordered_corr.to_numpy(), cmap="viridis", vmin=ordered_corr.to_numpy().min(), vmax=1.0)
    ax_heat.set_xticks(range(len(order)))
    ax_heat.set_xticklabels(order, rotation=90, fontsize=8)
    ax_heat.set_yticks(range(len(order)))
    ax_heat.set_yticklabels(order, fontsize=8)
    fig.colorbar(im, ax=ax_heat, label="Spearman correlation (log2 CPM+1)")
    fig.suptitle("GSE240112 tumor-cell pseudobulk: sample-sample correlation", fontsize=11)
    fig.tight_layout()
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info("wrote %s", out_path)


def plot_pca(pca_df: pd.DataFrame, out_path: str | Path) -> None:
    pc1 = pca_df.loc[pca_df["pc"] == "PC1"].set_index("sample_id")
    pc2 = pca_df.loc[pca_df["pc"] == "PC2"].set_index("sample_id")
    var1 = pc1["variance_explained_fraction"].iloc[0] * 100
    var2 = pc2["variance_explained_fraction"].iloc[0] * 100

    fig, ax = plt.subplots(figsize=(6, 5.5))
    for sample_id in pc1.index:
        group = pc1.loc[sample_id, "group"]
        ax.scatter(pc1.loc[sample_id, "coordinate"], pc2.loc[sample_id, "coordinate"], color=GROUP_COLORS[group], s=80, label=group)
        ax.annotate(sample_id, (pc1.loc[sample_id, "coordinate"], pc2.loc[sample_id, "coordinate"]), fontsize=8, xytext=(4, 4), textcoords="offset points")

    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    ax.legend(by_label.values(), by_label.keys(), loc="best", fontsize=9)
    ax.set_xlabel(f"PC1 ({var1:.1f}% variance)")
    ax.set_ylabel(f"PC2 ({var2:.1f}% variance)")
    ax.set_title("GSE240112 tumor-cell pseudobulk PCA (PT vs RT)", fontsize=10)
    fig.tight_layout()
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info("wrote %s", out_path)


def run_pseudobulk_qc(config_path: str | Path = "config/config.yaml", n_pca_genes: int = 2000) -> dict[str, pd.DataFrame]:
    config = _load_config(config_path)
    cfg = config["gse240112"]
    out = cfg["output"]

    counts, metadata = load_pseudobulk(out["tumor_cell"]["counts_tsv"], out["tumor_cell"]["metadata_tsv"])
    sample_ids = metadata["sample_id"].tolist()

    qc_summary = build_qc_summary(metadata)
    corr_long = compute_sample_correlations(counts, sample_ids)
    n_pca = min(n_pca_genes, len(counts))
    top_var = select_top_variable_genes(counts, sample_ids, n_pca)
    pca_df = compute_pca(top_var, metadata, sample_ids)

    Path(out["qc"]["qc_tsv"]).parent.mkdir(parents=True, exist_ok=True)
    qc_summary.to_csv(out["qc"]["qc_tsv"], sep="\t", index=False)
    corr_long.to_csv(out["qc"]["correlation_tsv"], sep="\t", index=False)
    pca_df.to_csv(out["qc"]["pca_tsv"], sep="\t", index=False)

    figures_dir = Path(out["qc"]["figures_dir"])
    plot_library_and_detected_genes(metadata, figures_dir / "tumor_cell_library_and_genes.png")
    plot_sample_correlation_heatmap(corr_long, figures_dir / "tumor_cell_correlation_heatmap.png")
    plot_pca(pca_df, figures_dir / "tumor_cell_pca.png")

    return {"qc_summary": qc_summary, "correlation": corr_long, "pca": pca_df}


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run_pseudobulk_qc()
