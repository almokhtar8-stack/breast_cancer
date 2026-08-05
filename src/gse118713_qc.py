"""Sample QC (totals, correlations, PCA, clustering) for GSE118713 Phase 2B.

Source: the frozen gene-level TPM matrix validated and loaded by
``src.gse118713_expression_filter`` (config
``gse118713.output.gene_tpm_parquet``, checksum-pinned). Sample/group
metadata comes from ``gse118713.output.sample_metadata_tsv``.

Per-sample descriptive statistics (total TPM, genes detected, median TPM,
IQR) and sample-sample Spearman correlations are computed on every gene in
the frozen matrix. PCA is restricted to the 2,000 most variable genes
*after* the PREANALYSIS.md-preregistered TPM filter (see
``src.gse118713_expression_filter.filter_expression``), per the Phase 2B
statistical-plan amendment. No sample is ever removed automatically here --
this module only describes the data; removal decisions are out of scope.

Figures are generated only from the TSV tables this module writes, never
from in-memory objects, so a figure can always be regenerated from the
saved outputs alone.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import squareform

logger = logging.getLogger(__name__)


def _load_config(config_path: str | Path) -> dict:
    with open(config_path) as f:
        return yaml.safe_load(f)


@dataclass(frozen=True)
class QcConfig:
    n_pca_genes: int
    sample_qc_tsv: Path
    sample_correlations_tsv: Path
    pca_coordinates_tsv: Path
    pca_figure_pdf: Path
    correlation_figure_pdf: Path

    @classmethod
    def from_config(cls, config: dict) -> "QcConfig":
        qc = config["gse118713_phase2b"]["qc"]
        return cls(
            n_pca_genes=int(qc["n_pca_genes"]),
            sample_qc_tsv=Path(qc["sample_qc_tsv"]),
            sample_correlations_tsv=Path(qc["sample_correlations_tsv"]),
            pca_coordinates_tsv=Path(qc["pca_coordinates_tsv"]),
            pca_figure_pdf=Path(qc["pca_figure_pdf"]),
            correlation_figure_pdf=Path(qc["correlation_figure_pdf"]),
        )


def compute_sample_summary(df: pd.DataFrame, meta: pd.DataFrame, sample_ids: list[str]) -> pd.DataFrame:
    """Per-sample total TPM, genes>0, genes>=1, median TPM, and IQR.

    Computed over every gene in the (unfiltered) frozen matrix -- these are
    global descriptive statistics, not restricted to the TPM filter used
    for PCA gene selection.
    """
    rows = []
    group_by_sample = dict(zip(meta["sample_id"], meta["group"]))
    for sample_id in sample_ids:
        values = df[sample_id].to_numpy(dtype=float)
        q1, q3 = np.percentile(values, [25, 75])
        rows.append(
            {
                "sample_id": sample_id,
                "group": group_by_sample[sample_id],
                "total_tpm": float(values.sum()),
                "n_genes_tpm_gt_0": int((values > 0).sum()),
                "n_genes_tpm_ge_1": int((values >= 1).sum()),
                "median_tpm": float(np.median(values)),
                "iqr_tpm": float(q3 - q1),
            }
        )
    logger.info("compute_sample_summary: computed QC summary for %d samples", len(rows))
    return pd.DataFrame(rows)


def compute_sample_correlations(df: pd.DataFrame, sample_ids: list[str]) -> pd.DataFrame:
    """Pairwise Spearman correlation of log2(TPM+1) across every gene, long format."""
    log2_df = np.log2(df[sample_ids].to_numpy(dtype=float) + 1)
    corr = pd.DataFrame(log2_df, columns=sample_ids).corr(method="spearman")

    records = []
    for s1 in sample_ids:
        for s2 in sample_ids:
            records.append({"sample_id_1": s1, "sample_id_2": s2, "spearman_r": float(corr.loc[s1, s2])})
    logger.info("compute_sample_correlations: computed %dx%d Spearman correlation matrix", len(sample_ids), len(sample_ids))
    return pd.DataFrame(records)


def select_top_variable_genes(filtered_df: pd.DataFrame, sample_ids: list[str], n_genes: int) -> pd.DataFrame:
    """Deterministically select the ``n_genes`` most variable filtered genes.

    Variance is computed on log2(TPM+1) across samples. Ties are broken by
    ascending ``gene_id`` so the selection does not depend on row order or
    the sort stability of the underlying implementation.
    """
    log2_vals = np.log2(filtered_df[sample_ids].to_numpy(dtype=float) + 1)
    variance = log2_vals.var(axis=1, ddof=1)
    ranked = filtered_df.assign(_variance=variance).sort_values(
        ["_variance", "gene_id"], ascending=[False, True]
    )
    n_take = min(n_genes, len(ranked))
    top = ranked.head(n_take).drop(columns="_variance").reset_index(drop=True)
    logger.info("select_top_variable_genes: selected %d of %d filtered genes", len(top), len(filtered_df))
    return top


def compute_pca(top_variable_df: pd.DataFrame, meta: pd.DataFrame, sample_ids: list[str]) -> pd.DataFrame:
    """PCA of samples over the supplied (already gene-selected) matrix.

    Genes are centered (mean-subtracted across samples) before SVD.
    Component sign is fixed deterministically: for each PC, the sign is
    chosen so the sample with the largest-magnitude score is positive.
    Returned in long format: one row per (sample_id, pc).
    """
    log2_vals = np.log2(top_variable_df[sample_ids].to_numpy(dtype=float) + 1)
    matrix = log2_vals.T  # samples x genes
    centered = matrix - matrix.mean(axis=0, keepdims=True)

    u, s, _vt = np.linalg.svd(centered, full_matrices=False)
    scores = u * s  # samples x k
    variance_explained = (s**2) / np.sum(s**2)

    n_components = scores.shape[1]
    for k in range(n_components):
        col = scores[:, k]
        sign = 1.0 if col[np.argmax(np.abs(col))] >= 0 else -1.0
        scores[:, k] = col * sign

    group_by_sample = dict(zip(meta["sample_id"], meta["group"]))
    records = []
    for i, sample_id in enumerate(sample_ids):
        for k in range(n_components):
            records.append(
                {
                    "sample_id": sample_id,
                    "group": group_by_sample[sample_id],
                    "pc": f"PC{k + 1}",
                    "coordinate": float(scores[i, k]),
                    "variance_explained_fraction": float(variance_explained[k]),
                }
            )
    logger.info("compute_pca: computed %d principal components for %d samples", n_components, len(sample_ids))
    return pd.DataFrame(records)


def write_sample_summary(summary_df: pd.DataFrame, cfg: QcConfig) -> None:
    cfg.sample_qc_tsv.parent.mkdir(parents=True, exist_ok=True)
    summary_df.to_csv(cfg.sample_qc_tsv, sep="\t", index=False)
    logger.info("write_sample_summary: wrote %s", cfg.sample_qc_tsv)


def write_sample_correlations(corr_df: pd.DataFrame, cfg: QcConfig) -> None:
    cfg.sample_correlations_tsv.parent.mkdir(parents=True, exist_ok=True)
    corr_df.to_csv(cfg.sample_correlations_tsv, sep="\t", index=False)
    logger.info("write_sample_correlations: wrote %s", cfg.sample_correlations_tsv)


def write_pca_coordinates(pca_df: pd.DataFrame, cfg: QcConfig) -> None:
    cfg.pca_coordinates_tsv.parent.mkdir(parents=True, exist_ok=True)
    pca_df.to_csv(cfg.pca_coordinates_tsv, sep="\t", index=False)
    logger.info("write_pca_coordinates: wrote %s", cfg.pca_coordinates_tsv)


_GROUP_COLORS = {"MCF7": "#4C72B0", "TAMR": "#DD8452", "FASR": "#55A868"}


def plot_pca(pca_coordinates_tsv: str | Path, output_pdf: str | Path) -> None:
    """Regenerate the PCA scatter (PC1 vs PC2) purely from the saved PCA table."""
    pca_df = pd.read_csv(pca_coordinates_tsv, sep="\t")
    wide = pca_df.pivot(index=["sample_id", "group"], columns="pc", values="coordinate").reset_index()
    variance = pca_df.drop_duplicates("pc").set_index("pc")["variance_explained_fraction"]

    fig, ax = plt.subplots(figsize=(6, 5))
    for group, sub in wide.groupby("group"):
        ax.scatter(sub["PC1"], sub["PC2"], label=group, color=_GROUP_COLORS.get(group), s=60)
    for _, row in wide.iterrows():
        ax.annotate(row["sample_id"], (row["PC1"], row["PC2"]), fontsize=6, xytext=(3, 3), textcoords="offset points")
    ax.set_xlabel(f"PC1 ({variance['PC1'] * 100:.1f}%)")
    ax.set_ylabel(f"PC2 ({variance['PC2'] * 100:.1f}%)")
    ax.set_title("GSE118713: PCA of top variable genes")
    ax.legend()
    fig.tight_layout()

    output_pdf = Path(output_pdf)
    output_pdf.parent.mkdir(parents=True, exist_ok=True)
    # CreationDate omitted so identical inputs always produce a byte-identical PDF.
    fig.savefig(output_pdf, metadata={"CreationDate": None})
    plt.close(fig)
    logger.info("plot_pca: wrote %s from %s", output_pdf, pca_coordinates_tsv)


def plot_sample_correlation(sample_correlations_tsv: str | Path, output_pdf: str | Path) -> None:
    """Regenerate the correlation heatmap + dendrogram purely from the saved correlation table."""
    corr_long = pd.read_csv(sample_correlations_tsv, sep="\t")
    corr_wide = corr_long.pivot(index="sample_id_1", columns="sample_id_2", values="spearman_r")
    sample_ids = list(corr_wide.index)

    distance = 1.0 - corr_wide.to_numpy()
    np.fill_diagonal(distance, 0.0)
    condensed = squareform(distance, checks=False)
    link = linkage(condensed, method="average")

    fig, (ax_dendro, ax_heat) = plt.subplots(
        2, 1, figsize=(7, 8), gridspec_kw={"height_ratios": [1, 3]}, sharex=False
    )
    dendro = dendrogram(link, labels=sample_ids, ax=ax_dendro)
    ax_dendro.set_title("GSE118713: hierarchical clustering (1 - Spearman r, average linkage)")

    ordered = dendro["ivl"]
    ordered_corr = corr_wide.loc[ordered, ordered]
    im = ax_heat.imshow(ordered_corr.to_numpy(), vmin=ordered_corr.to_numpy().min(), vmax=1.0, cmap="viridis")
    ax_heat.set_xticks(range(len(ordered)))
    ax_heat.set_xticklabels(ordered, rotation=90, fontsize=7)
    ax_heat.set_yticks(range(len(ordered)))
    ax_heat.set_yticklabels(ordered, fontsize=7)
    fig.colorbar(im, ax=ax_heat, label="Spearman r")
    fig.tight_layout()

    output_pdf = Path(output_pdf)
    output_pdf.parent.mkdir(parents=True, exist_ok=True)
    # CreationDate omitted so identical inputs always produce a byte-identical PDF.
    fig.savefig(output_pdf, metadata={"CreationDate": None})
    plt.close(fig)
    logger.info("plot_sample_correlation: wrote %s from %s", output_pdf, sample_correlations_tsv)


def run_sample_qc(
    df: pd.DataFrame,
    filtered_df: pd.DataFrame,
    meta: pd.DataFrame,
    sample_ids: list[str],
    config_path: str | Path = "config/config.yaml",
) -> dict[str, pd.DataFrame]:
    """Run the full Phase 2B sample-QC step and write tables + figures."""
    config = _load_config(config_path)
    cfg = QcConfig.from_config(config)

    summary_df = compute_sample_summary(df, meta, sample_ids)
    write_sample_summary(summary_df, cfg)

    corr_df = compute_sample_correlations(df, sample_ids)
    write_sample_correlations(corr_df, cfg)

    top_variable = select_top_variable_genes(filtered_df, sample_ids, cfg.n_pca_genes)
    pca_df = compute_pca(top_variable, meta, sample_ids)
    write_pca_coordinates(pca_df, cfg)

    plot_pca(cfg.pca_coordinates_tsv, cfg.pca_figure_pdf)
    plot_sample_correlation(cfg.sample_correlations_tsv, cfg.correlation_figure_pdf)

    return {"sample_qc": summary_df, "sample_correlations": corr_df, "pca_coordinates": pca_df}
