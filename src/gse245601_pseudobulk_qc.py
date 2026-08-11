"""Pseudobulk QC for GSE245601 Track A (epithelial) and Track B (strict
malignant) -- library size, detected genes, sample-sample correlation, and
PCA, all computed on the pseudobulk count matrices written by
``scripts/analysis/gse245601_13_build_pseudobulk.R``.

No batch-correction method (Harmony, Seurat integration, ComBat, or
otherwise) is applied anywhere in this module -- patient pairing is
handled statistically in the edgeR design, not here. This module answers
one question only: do pseudobulk samples cluster predominantly by
patient, by treatment, or by an obvious technical feature? It does not
exclude any sample; if a catastrophic quality problem is found, it is
reported, not silently dropped (gse245601_PREANALYSIS.md section 13 /
Phase 4 instruction).

Mirrors the established QC pattern in ``src/gse118713_qc.py`` (log2(x+1)
transform, Spearman sample-sample correlation, PCA via numpy SVD on the
most variable genes, deterministic PC sign-fixing) adapted for CPM-based
pseudobulk counts and a (patient, condition) design instead of a single
group label.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import yaml  # noqa: E402
from scipy.cluster.hierarchy import dendrogram, linkage  # noqa: E402

logger = logging.getLogger(__name__)

PATIENT_COLORS = plt.get_cmap("tab10").colors
CONDITION_MARKERS = {"Control": "o", "Tamoxifen": "^"}


def _load_config(config_path: str | Path) -> dict:
    with open(config_path) as f:
        return yaml.safe_load(f)


def load_pseudobulk(counts_path: str | Path, metadata_path: str | Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Load a pseudobulk count matrix (genes x samples, first column
    'gene') and its metadata table. Read-only."""
    counts = pd.read_csv(counts_path, sep="\t").set_index("gene")
    metadata = pd.read_csv(metadata_path, sep="\t")
    missing = set(metadata["sample_id"]) - set(counts.columns)
    if missing:
        raise ValueError(f"metadata sample_id(s) not present as count matrix columns: {missing}")
    logger.info("load_pseudobulk: %d genes x %d samples", counts.shape[0], len(metadata))
    return counts, metadata


def compute_log2cpm(counts: pd.DataFrame, sample_ids: list[str]) -> pd.DataFrame:
    """log2(CPM+1), library size from each sample's own column sum (not
    the frozen metadata total_library_size column, though the two must
    agree -- checked by the caller/tests)."""
    lib_sizes = counts[sample_ids].sum(axis=0)
    cpm = counts[sample_ids].div(lib_sizes, axis=1) * 1e6
    return np.log2(cpm + 1)


def compute_sample_correlations(counts: pd.DataFrame, sample_ids: list[str]) -> pd.DataFrame:
    """Pairwise Spearman correlation of log2(CPM+1) across every gene,
    long format (sample_id_1, sample_id_2, correlation)."""
    log2cpm = compute_log2cpm(counts, sample_ids)
    corr = log2cpm.corr(method="spearman")
    corr.index.name = "sample_id_1"
    long = corr.reset_index().melt(id_vars="sample_id_1", var_name="sample_id_2", value_name="correlation")
    logger.info("compute_sample_correlations: %dx%d Spearman correlation matrix", len(sample_ids), len(sample_ids))
    return long


def select_top_variable_genes(counts: pd.DataFrame, sample_ids: list[str], n_genes: int) -> pd.DataFrame:
    """Top-N most variable genes by variance of log2(CPM+1) across
    samples (ties broken by gene name for determinism)."""
    log2cpm = compute_log2cpm(counts, sample_ids)
    variance = log2cpm.var(axis=1, ddof=1)
    ordered = variance.sort_values(ascending=False, kind="stable")
    top_genes = ordered.index[:n_genes]
    logger.info("select_top_variable_genes: selected %d of %d genes", len(top_genes), len(counts))
    return counts.loc[top_genes]


def compute_pca(top_variable_counts: pd.DataFrame, metadata: pd.DataFrame, sample_ids: list[str]) -> pd.DataFrame:
    """PCA of pseudobulk samples over the supplied (already gene-selected)
    matrix. Genes centered (mean-subtracted) before SVD; PC sign fixed
    deterministically (largest-magnitude score positive). Long format: one
    row per (sample_id, pc)."""
    log2cpm = compute_log2cpm(top_variable_counts, sample_ids)
    matrix = log2cpm[sample_ids].to_numpy(dtype=float).T  # samples x genes
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
                    "patient": meta_by_sample.loc[sample_id, "patient"],
                    "condition": meta_by_sample.loc[sample_id, "condition"],
                    "pc": f"PC{k + 1}",
                    "coordinate": float(scores[i, k]),
                    "variance_explained_fraction": float(variance_explained[k]),
                }
            )
    logger.info("compute_pca: %d PCs for %d samples", scores.shape[1], len(sample_ids))
    return pd.DataFrame(records)


def build_qc_summary(metadata: pd.DataFrame) -> pd.DataFrame:
    """Per-sample QC summary -- passthrough of the frozen pseudobulk
    construction metadata (n_contributing_cells, total_library_size,
    n_detected_genes), sorted for readability. No sample is flagged for
    exclusion here; catastrophic problems are reported in the pipeline
    log/report text, not auto-excluded."""
    out = metadata.sort_values(["patient", "condition"]).reset_index(drop=True)
    logger.info(
        "build_qc_summary: library size range %.2fM-%.2fM, detected genes range %d-%d",
        out["total_library_size"].min() / 1e6,
        out["total_library_size"].max() / 1e6,
        out["n_detected_genes"].min(),
        out["n_detected_genes"].max(),
    )
    return out


def _patient_color(patient: str, patients_sorted: list[str]) -> str:
    return PATIENT_COLORS[patients_sorted.index(patient) % len(PATIENT_COLORS)]


def plot_library_and_detected_genes(metadata: pd.DataFrame, track_label: str, out_path: str | Path) -> None:
    """Two panels: library size and detected-gene count, one paired
    Control->Tamoxifen line per patient."""
    patients_sorted = sorted(metadata["patient"].unique())
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    for ax, col, title in zip(axes, ("total_library_size", "n_detected_genes"), ("total pseudobulk library size", "detected genes")):
        for patient in patients_sorted:
            sub = metadata.loc[metadata["patient"] == patient].set_index("condition")
            if not {"Control", "Tamoxifen"}.issubset(sub.index):
                continue
            color = _patient_color(patient, patients_sorted)
            ax.plot([0, 1], [sub.loc["Control", col], sub.loc["Tamoxifen", col]], color=color, marker="o", label=patient)
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["Control", "Tamoxifen"])
        ax.set_title(title, fontsize=10)
        ax.set_ylabel(col)
    axes[1].legend(fontsize=6, loc="center left", bbox_to_anchor=(1.0, 0.5))
    fig.suptitle(f"{track_label}: pseudobulk library size / detected genes by patient", fontsize=11)
    fig.tight_layout()
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info("wrote %s", out_path)


def plot_sample_correlation_heatmap(corr_long: pd.DataFrame, track_label: str, out_path: str | Path) -> None:
    """Correlation heatmap with a dendrogram, samples labeled by
    patient_condition."""
    corr = corr_long.pivot(index="sample_id_1", columns="sample_id_2", values="correlation")
    dist = 1 - corr
    link = linkage(dist.to_numpy()[np.triu_indices(len(dist), k=1)], method="average")

    fig, (ax_dendro, ax_heat) = plt.subplots(1, 2, figsize=(12, 6), gridspec_kw={"width_ratios": [1, 3]})
    dend = dendrogram(link, labels=corr.index.tolist(), orientation="left", ax=ax_dendro)
    ax_dendro.set_xticks([])

    order = dend["ivl"]
    ordered_corr = corr.loc[order, order]
    im = ax_heat.imshow(ordered_corr.to_numpy(), cmap="viridis", vmin=ordered_corr.to_numpy().min(), vmax=1.0)
    ax_heat.set_xticks(range(len(order)))
    ax_heat.set_xticklabels(order, rotation=90, fontsize=7)
    ax_heat.set_yticks(range(len(order)))
    ax_heat.set_yticklabels(order, fontsize=7)
    fig.colorbar(im, ax=ax_heat, label="Spearman correlation (log2 CPM+1)")
    fig.suptitle(f"{track_label}: pseudobulk sample-sample correlation", fontsize=11)
    fig.tight_layout()
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info("wrote %s", out_path)


def plot_pca(pca_df: pd.DataFrame, track_label: str, out_path: str | Path) -> None:
    """PC1 vs PC2, color = patient, marker shape = condition."""
    patients_sorted = sorted(pca_df["patient"].unique())
    pc1 = pca_df.loc[pca_df["pc"] == "PC1"].set_index("sample_id")
    pc2 = pca_df.loc[pca_df["pc"] == "PC2"].set_index("sample_id")
    var1 = pc1["variance_explained_fraction"].iloc[0] * 100
    var2 = pc2["variance_explained_fraction"].iloc[0] * 100

    fig, ax = plt.subplots(figsize=(6.5, 5.5))
    for sample_id in pc1.index:
        patient = pc1.loc[sample_id, "patient"]
        condition = pc1.loc[sample_id, "condition"]
        color = _patient_color(patient, patients_sorted)
        marker = CONDITION_MARKERS[condition]
        ax.scatter(pc1.loc[sample_id, "coordinate"], pc2.loc[sample_id, "coordinate"], color=color, marker=marker, s=60)
        ax.annotate(patient, (pc1.loc[sample_id, "coordinate"], pc2.loc[sample_id, "coordinate"]), fontsize=6, xytext=(3, 3), textcoords="offset points")

    from matplotlib.lines import Line2D

    marker_handles = [Line2D([0], [0], marker=m, color="gray", linestyle="", label=c, markersize=8) for c, m in CONDITION_MARKERS.items()]
    ax.legend(handles=marker_handles, loc="best", fontsize=8)
    ax.set_xlabel(f"PC1 ({var1:.1f}% variance)")
    ax.set_ylabel(f"PC2 ({var2:.1f}% variance)")
    ax.set_title(f"{track_label}: pseudobulk PCA (color=patient, shape=condition)", fontsize=10)
    fig.tight_layout()
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info("wrote %s", out_path)


@dataclass(frozen=True)
class PseudobulkQcPaths:
    counts_tsv: Path
    metadata_tsv: Path
    qc_tsv: Path
    correlation_tsv: Path
    pca_tsv: Path
    figures_dir: Path
    track_label: str


def run_track_qc(paths: PseudobulkQcPaths, n_pca_genes: int = 2000) -> dict[str, pd.DataFrame]:
    counts, metadata = load_pseudobulk(paths.counts_tsv, paths.metadata_tsv)
    sample_ids = metadata["sample_id"].tolist()

    qc_summary = build_qc_summary(metadata)
    corr_long = compute_sample_correlations(counts, sample_ids)
    n_pca = min(n_pca_genes, len(counts))
    top_var = select_top_variable_genes(counts, sample_ids, n_pca)
    pca_df = compute_pca(top_var, metadata, sample_ids)

    paths.qc_tsv.parent.mkdir(parents=True, exist_ok=True)
    qc_summary.to_csv(paths.qc_tsv, sep="\t", index=False)
    corr_long.to_csv(paths.correlation_tsv, sep="\t", index=False)
    pca_df.to_csv(paths.pca_tsv, sep="\t", index=False)

    plot_library_and_detected_genes(metadata, paths.track_label, paths.figures_dir / f"{paths.track_label}_library_and_genes.png")
    plot_sample_correlation_heatmap(corr_long, paths.track_label, paths.figures_dir / f"{paths.track_label}_correlation_heatmap.png")
    plot_pca(pca_df, paths.track_label, paths.figures_dir / f"{paths.track_label}_pca.png")

    return {"qc_summary": qc_summary, "correlation": corr_long, "pca": pca_df}


def run_pseudobulk_qc(config_path: str | Path = "config/config.yaml") -> dict[str, dict[str, pd.DataFrame]]:
    config = _load_config(config_path)
    cfg = config["gse245601_pseudobulk"]
    out = cfg["output"]

    track_a_paths = PseudobulkQcPaths(
        counts_tsv=Path(out["track_a"]["counts_tsv"]),
        metadata_tsv=Path(out["track_a"]["metadata_tsv"]),
        qc_tsv=Path(out["qc"]["track_a_qc_tsv"]),
        correlation_tsv=Path(out["qc"]["track_a_correlation_tsv"]),
        pca_tsv=Path(out["qc"]["track_a_pca_tsv"]),
        figures_dir=Path(out["qc"]["figures_dir"]),
        track_label="track_a",
    )
    track_b_paths = PseudobulkQcPaths(
        counts_tsv=Path(out["track_b"]["counts_tsv"]),
        metadata_tsv=Path(out["track_b"]["metadata_tsv"]),
        qc_tsv=Path(out["qc"]["track_b_qc_tsv"]),
        correlation_tsv=Path(out["qc"]["track_b_correlation_tsv"]),
        pca_tsv=Path(out["qc"]["track_b_pca_tsv"]),
        figures_dir=Path(out["qc"]["figures_dir"]),
        track_label="track_b",
    )

    return {"track_a": run_track_qc(track_a_paths), "track_b": run_track_qc(track_b_paths, n_pca_genes=2000)}


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run_pseudobulk_qc()
