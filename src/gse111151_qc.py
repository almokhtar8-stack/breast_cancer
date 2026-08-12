"""QC for the GSE111151 independent tamoxifen-resistance validation panel
-- library size, detected genes, sample-sample correlation, and PCA on
the 11-sample raw count matrix written by
``scripts/analysis/gse111151_01_build_count_matrix.R``.

Samples are labeled by cell_line (4 levels) and resistance_status
(parental/resistant); strong clustering by cell line is expected (see
docs/GSE111151_PREANALYSIS.md section J) and is not itself a QC failure
-- it is exactly why the edgeR model blocks on cell_line rather than
ignoring it. No sample is excluded here regardless of how it looks; a
catastrophic technical problem would be reported, not silently dropped.

Data source: GSE111151 (Hultsch et al., BMC Cancer 2018, PMID 30143015),
version as downloaded 2026-08-12 (data/raw/gse111151/MANIFEST.tsv).
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

CELL_LINE_COLORS = {"MCF-7": "#4C72B0", "T-47D": "#DD8452", "ZR-75-1": "#55A868", "BT-474": "#C44E52"}
STATUS_MARKERS = {"parental": "o", "resistant": "^"}


def _load_config(config_path: str | Path) -> dict:
    with open(config_path) as f:
        return yaml.safe_load(f)


def load_counts(counts_path: str | Path, metadata_path: str | Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Load the genes x samples raw count matrix (first two columns
    'gene_id', 'gene_name') and its sample metadata. Read-only."""
    counts = pd.read_csv(counts_path, sep="\t").set_index("gene_id")
    gene_names = counts.pop("gene_name")
    metadata = pd.read_csv(metadata_path, sep="\t")
    missing = set(metadata["sample_id"]) - set(counts.columns)
    if missing:
        raise ValueError(f"metadata sample_id(s) not present as count matrix columns: {missing}")
    logger.info("load_counts: %d genes x %d samples", counts.shape[0], len(metadata))
    return counts, metadata, gene_names


def load_tmm_norm_factors(path: str | Path) -> pd.Series:
    """Per-sample effective library size (raw library size x edgeR TMM
    norm.factor), written by ``scripts/analysis/gse111151_02_edger.R``.
    Using these (rather than a raw-library-size-only CPM) matters here:
    TMM factors range 0.82-1.30 across these 11 samples (e.g. MCF-7
    parental=1.30 vs MCF-7_Tam1=0.97), so naive CPM materially
    misrepresents several samples' relative expression -- confirmed by a
    spurious apparent per-cell-line direction for several candidates
    before this correction (docs/GSE111151_ANALYSIS_REPORT.md limitations)."""
    df = pd.read_csv(path, sep="\t")
    return df.set_index("sample_id")["effective_library_size"]


def compute_log2cpm(counts: pd.DataFrame, sample_ids: list[str], effective_lib_sizes: pd.Series | None = None) -> pd.DataFrame:
    """log2(CPM+1). If ``effective_lib_sizes`` (TMM-adjusted, from
    ``load_tmm_norm_factors``) is provided, it is used instead of each
    sample's raw column sum -- see ``load_tmm_norm_factors`` docstring."""
    if effective_lib_sizes is not None:
        lib_sizes = effective_lib_sizes.loc[sample_ids]
    else:
        lib_sizes = counts[sample_ids].sum(axis=0)
    cpm = counts[sample_ids].div(lib_sizes, axis=1) * 1e6
    return np.log2(cpm + 1)


def compute_sample_correlations(counts: pd.DataFrame, sample_ids: list[str], effective_lib_sizes: pd.Series | None = None) -> pd.DataFrame:
    log2cpm = compute_log2cpm(counts, sample_ids, effective_lib_sizes)
    corr = log2cpm.corr(method="spearman")
    corr.index.name = "sample_id_1"
    long = corr.reset_index().melt(id_vars="sample_id_1", var_name="sample_id_2", value_name="correlation")
    logger.info("compute_sample_correlations: %dx%d Spearman correlation matrix", len(sample_ids), len(sample_ids))
    return long


def select_top_variable_genes(counts: pd.DataFrame, sample_ids: list[str], n_genes: int, effective_lib_sizes: pd.Series | None = None) -> pd.DataFrame:
    log2cpm = compute_log2cpm(counts, sample_ids, effective_lib_sizes)
    variance = log2cpm.var(axis=1, ddof=1)
    ordered = variance.sort_values(ascending=False, kind="stable")
    top_genes = ordered.index[:n_genes]
    logger.info("select_top_variable_genes: selected %d of %d genes", len(top_genes), len(counts))
    return counts.loc[top_genes]


def compute_pca(top_variable_counts: pd.DataFrame, metadata: pd.DataFrame, sample_ids: list[str], effective_lib_sizes: pd.Series | None = None) -> pd.DataFrame:
    log2cpm = compute_log2cpm(top_variable_counts, sample_ids, effective_lib_sizes)
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
                    "cell_line": meta_by_sample.loc[sample_id, "cell_line"],
                    "resistance_status": meta_by_sample.loc[sample_id, "resistance_status"],
                    "pc": f"PC{k + 1}",
                    "coordinate": float(scores[i, k]),
                    "variance_explained_fraction": float(variance_explained[k]),
                }
            )
    logger.info("compute_pca: %d PCs for %d samples", scores.shape[1], len(sample_ids))
    return pd.DataFrame(records)


def build_qc_summary(metadata: pd.DataFrame) -> pd.DataFrame:
    out = metadata.sort_values(["cell_line", "resistance_status"]).reset_index(drop=True)
    logger.info(
        "build_qc_summary: library size range %.1fM-%.1fM, detected genes range %d-%d",
        out["library_size"].min() / 1e6,
        out["library_size"].max() / 1e6,
        out["n_detected_genes"].min(),
        out["n_detected_genes"].max(),
    )
    return out


def plot_library_and_detected_genes(metadata: pd.DataFrame, out_path: str | Path) -> None:
    order = metadata.sort_values(["cell_line", "resistance_status"])
    fig, axes = plt.subplots(1, 2, figsize=(9, 4.5))
    for ax, col, title in zip(axes, ("library_size", "n_detected_genes"), ("library size (raw counts)", "detected genes")):
        colors = [CELL_LINE_COLORS[c] for c in order["cell_line"]]
        ax.bar(order["sample_id"], order[col], color=colors)
        ax.set_title(title, fontsize=10)
        ax.tick_params(axis="x", rotation=90)
    from matplotlib.patches import Patch

    handles = [Patch(color=c, label=g) for g, c in CELL_LINE_COLORS.items()]
    axes[1].legend(handles=handles, fontsize=7)
    fig.suptitle("GSE111151: library size / detected genes by sample", fontsize=11)
    fig.tight_layout()
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info("wrote %s", out_path)


def plot_sample_correlation_heatmap(corr_long: pd.DataFrame, out_path: str | Path) -> None:
    corr = corr_long.pivot(index="sample_id_1", columns="sample_id_2", values="correlation")
    dist = 1 - corr
    link = linkage(dist.to_numpy()[np.triu_indices(len(dist), k=1)], method="average")

    fig, (ax_dendro, ax_heat) = plt.subplots(1, 2, figsize=(10, 5), gridspec_kw={"width_ratios": [1, 2]})
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
    fig.suptitle("GSE111151: sample-sample correlation", fontsize=11)
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

    fig, ax = plt.subplots(figsize=(6.5, 5.5))
    for sample_id in pc1.index:
        cell_line = pc1.loc[sample_id, "cell_line"]
        status = pc1.loc[sample_id, "resistance_status"]
        ax.scatter(pc1.loc[sample_id, "coordinate"], pc2.loc[sample_id, "coordinate"], color=CELL_LINE_COLORS[cell_line], marker=STATUS_MARKERS[status], s=90)
        ax.annotate(sample_id, (pc1.loc[sample_id, "coordinate"], pc2.loc[sample_id, "coordinate"]), fontsize=7, xytext=(4, 4), textcoords="offset points")

    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch

    color_handles = [Patch(color=c, label=g) for g, c in CELL_LINE_COLORS.items()]
    marker_handles = [Line2D([0], [0], marker=m, color="gray", linestyle="", label=s, markersize=8) for s, m in STATUS_MARKERS.items()]
    ax.legend(handles=color_handles + marker_handles, loc="best", fontsize=7)
    ax.set_xlabel(f"PC1 ({var1:.1f}% variance)")
    ax.set_ylabel(f"PC2 ({var2:.1f}% variance)")
    ax.set_title("GSE111151 PCA (color=cell line, shape=resistance status)", fontsize=10)
    fig.tight_layout()
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info("wrote %s", out_path)


def run_qc(config_path: str | Path = "config/config.yaml", n_pca_genes: int = 2000) -> dict[str, pd.DataFrame]:
    config = _load_config(config_path)
    cfg = config["gse111151"]
    out = cfg["output"]

    counts, metadata, _gene_names = load_counts(out["counts_tsv"], out["metadata_tsv"])
    sample_ids = metadata["sample_id"].tolist()
    effective_lib_sizes = load_tmm_norm_factors(out["de"]["tmm_norm_factors_tsv"])

    qc_summary = build_qc_summary(metadata)
    corr_long = compute_sample_correlations(counts, sample_ids, effective_lib_sizes)
    n_pca = min(n_pca_genes, len(counts))
    top_var = select_top_variable_genes(counts, sample_ids, n_pca, effective_lib_sizes)
    pca_df = compute_pca(top_var, metadata, sample_ids, effective_lib_sizes)

    Path(out["qc"]["qc_tsv"]).parent.mkdir(parents=True, exist_ok=True)
    qc_summary.to_csv(out["qc"]["qc_tsv"], sep="\t", index=False)
    corr_long.to_csv(out["qc"]["correlation_tsv"], sep="\t", index=False)
    pca_df.to_csv(out["qc"]["pca_tsv"], sep="\t", index=False)

    figures_dir = Path(out["qc"]["figures_dir"])
    plot_library_and_detected_genes(metadata, figures_dir / "library_and_genes.png")
    plot_sample_correlation_heatmap(corr_long, figures_dir / "correlation_heatmap.png")
    plot_pca(pca_df, figures_dir / "pca.png")

    return {"qc_summary": qc_summary, "correlation": corr_long, "pca": pca_df}


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run_qc()
