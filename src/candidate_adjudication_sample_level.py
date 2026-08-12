"""Candidate adjudication Phase 10: sample-level inspection of the seven
MULTIMODAL_STRONG genes in each of the four RNA datasets, using each
dataset's own already-validated normalized representation (never a new
statistic): GSE118713 log2(TPM+1) per replicate; GSE111151 TMM-adjusted
log2(CPM+1) per sample, parental->resistant connected within each
isogenic cell-line block; GSE240112 pseudobulk log2(CPM+1) per
PT1-3/RT1-3 sample, RT3 always visibly labeled (never removed); GSE245601
patient-level pseudobulk log2(CPM+1) deltas, Track A/B kept separate,
never treated as more than one observation per patient.

Reuses the already-frozen, generic (gene-list-parametrized) plotting
primitives from `src.gse111151_candidate_visualization`,
`src.gse240112_candidate_visualization`, and
`src.gse245601_candidate_visualization` wherever they already accept an
arbitrary gene list -- no new statistic invented for genes outside the
original 13-candidate set.
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

from src.gse111151_qc import compute_log2cpm as g111151_log2cpm
from src.gse111151_qc import load_counts as g111151_load_counts
from src.gse111151_qc import load_tmm_norm_factors
from src.gse240112_candidate_visualization import SAMPLE_ORDER as GSE240112_SAMPLE_ORDER
from src.gse240112_candidate_visualization import build_pseudobulk_log2cpm as g240112_log2cpm
from src.gse245601_candidate_visualization import plot_paired_candidates_grid
from src.gse245601_pseudobulk_qc import load_pseudobulk as g245601_load_pseudobulk

logger = logging.getLogger(__name__)

CELL_LINE_COLORS = {"MCF-7": "#4C72B0", "T-47D": "#DD8452", "ZR-75-1": "#55A868", "BT-474": "#C44E52"}
CELL_LINE_ORDER = ["MCF-7", "T-47D", "ZR-75-1", "BT-474"]


def _load_config(config_path: str | Path) -> dict:
    with open(config_path) as f:
        return yaml.safe_load(f)


def plot_gse118713_sample_level(genes: list[str], config: dict, out_dir: Path) -> None:
    g = config["gse118713"]
    tpm = pd.read_parquet(g["output"]["gene_tpm_parquet"])
    meta = pd.read_csv(g["output"]["sample_metadata_tsv"], sep="\t")
    present = [gene for gene in genes if gene in set(tpm["gene_symbol"])]
    n = len(genes)
    fig, axes = plt.subplots(1, n, figsize=(2.6 * n, 3.2), squeeze=False)
    group_order = ["MCF7", "TAMR", "FASR"]
    group_colors = {"MCF7": "#4C72B0", "TAMR": "#C44E52", "FASR": "#55A868"}
    for i, gene in enumerate(genes):
        ax = axes[0, i]
        if gene not in present:
            ax.text(0.5, 0.5, "not tested", ha="center", va="center", fontsize=8, transform=ax.transAxes)
            ax.set_title(gene, fontsize=9)
            continue
        row = tpm.loc[tpm["gene_symbol"] == gene].iloc[0]
        for xi, group in enumerate(group_order):
            sample_ids = meta.loc[meta["group"] == group, "sample_id"].tolist()
            vals = np.log2(row[sample_ids].astype(float) + 1)
            ax.scatter([xi] * len(vals), vals, color=group_colors[group], s=30, zorder=3)
            ax.scatter([xi], [vals.mean()], color="black", marker="_", s=200, zorder=4)
        ax.set_xticks(range(len(group_order)))
        ax.set_xticklabels(group_order, fontsize=8)
        ax.set_title(gene, fontsize=9)
        if i == 0:
            ax.set_ylabel("log2(TPM+1)", fontsize=8)
    fig.suptitle("GSE118713: individual replicates (n=3/group), black bar = mean", fontsize=10, y=1.05)
    fig.tight_layout()
    out_path = out_dir / "gse118713_multimodal7_sample_level.png"
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info("wrote %s", out_path)


def plot_gse111151_sample_level(genes: list[str], config: dict, out_dir: Path) -> None:
    g = config["gse111151"]
    counts, metadata, gene_names = g111151_load_counts(g["output"]["counts_tsv"], g["output"]["metadata_tsv"])
    effective_lib = load_tmm_norm_factors(g["output"]["de"]["tmm_norm_factors_tsv"])
    log2cpm = g111151_log2cpm(counts, metadata["sample_id"].tolist(), effective_lib)
    de = pd.read_csv(g["output"]["de"]["genomewide_tsv"], sep="\t")
    symbol_to_id = de.set_index("gene_name")["gene_id"].to_dict()

    n = len(genes)
    fig, axes = plt.subplots(1, n, figsize=(2.6 * n, 3.4), squeeze=False)
    meta_idx = metadata.set_index("sample_id")
    for i, gene in enumerate(genes):
        ax = axes[0, i]
        gene_id = symbol_to_id.get(gene)
        if gene_id is None or gene_id not in log2cpm.index:
            ax.text(0.5, 0.5, "not tested", ha="center", va="center", fontsize=8, transform=ax.transAxes)
            ax.set_title(gene, fontsize=9)
            continue
        vals = log2cpm.loc[gene_id]
        for cl_i, cl in enumerate(CELL_LINE_ORDER):
            cl_samples = meta_idx.loc[meta_idx["cell_line"] == cl]
            parental = cl_samples.loc[cl_samples["resistance_status"] == "parental"].index
            resistant = cl_samples.loc[cl_samples["resistance_status"] == "resistant"].index
            color = CELL_LINE_COLORS[cl]
            if len(parental) > 0:
                ax.scatter([0], [vals.loc[parental[0]]], color=color, s=35, zorder=3)
            for r in resistant:
                ax.scatter([1], [vals.loc[r]], color=color, s=35, zorder=3)
                if len(parental) > 0:
                    ax.plot([0, 1], [vals.loc[parental[0]], vals.loc[r]], color=color, alpha=0.6, linewidth=1, zorder=2)
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["parental", "resistant"], fontsize=8)
        ax.set_title(gene, fontsize=9)
        ax.set_xlim(-0.3, 1.3)
        if i == 0:
            ax.set_ylabel("TMM log2(CPM+1)", fontsize=8)
    handles = [plt.Line2D([0], [0], color=CELL_LINE_COLORS[cl], marker="o", linestyle="-", label=cl) for cl in CELL_LINE_ORDER]
    fig.legend(handles=handles, loc="upper center", ncol=4, fontsize=8, frameon=False, bbox_to_anchor=(0.5, 1.1))
    fig.suptitle("GSE111151: parental->resistant within each isogenic cell line (genuine pairing)", fontsize=10, y=1.18)
    fig.tight_layout()
    out_path = out_dir / "gse111151_multimodal7_sample_level.png"
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info("wrote %s", out_path)


def plot_gse240112_sample_level(genes: list[str], config: dict, out_dir: Path) -> None:
    g = config["gse240112"]["output"]["tumor_cell"]
    log2cpm = g240112_log2cpm(g["counts_tsv"], genes)
    n = len(genes)
    fig, axes = plt.subplots(1, n, figsize=(2.6 * n, 3.2), squeeze=False)
    colors = {"PT1": "#4C72B0", "PT2": "#4C72B0", "PT3": "#4C72B0", "RT1": "#C44E52", "RT2": "#C44E52", "RT3": "#8B0000"}
    for i, gene in enumerate(genes):
        ax = axes[0, i]
        if gene not in log2cpm.index:
            ax.text(0.5, 0.5, "not tested", ha="center", va="center", fontsize=8, transform=ax.transAxes)
            ax.set_title(gene, fontsize=9)
            continue
        vals = log2cpm.loc[gene]
        for xi, s in enumerate(GSE240112_SAMPLE_ORDER):
            marker = "^" if s == "RT3" else "o"
            ax.scatter([xi], [vals[s]], color=colors[s], s=60 if s == "RT3" else 40, marker=marker, zorder=3)
        ax.set_xticks(range(len(GSE240112_SAMPLE_ORDER)))
        ax.set_xticklabels(GSE240112_SAMPLE_ORDER, fontsize=7, rotation=45)
        ax.set_title(gene, fontsize=9)
        if i == 0:
            ax.set_ylabel("pseudobulk log2(CPM+1)", fontsize=8)
    fig.suptitle("GSE240112 tumor-cell pseudobulk: PT1-3 (blue) vs RT1-3 (red); RT3 shown as a triangle (low cell count)", fontsize=9.5, y=1.05)
    fig.tight_layout()
    out_path = out_dir / "gse240112_multimodal7_sample_level.png"
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info("wrote %s", out_path)


def plot_gse245601_sample_level(genes: list[str], config: dict, out_dir: Path) -> None:
    pb = config["gse245601_pseudobulk"]
    track_a_counts, track_a_meta = g245601_load_pseudobulk(pb["output"]["track_a"]["counts_tsv"], pb["output"]["track_a"]["metadata_tsv"])
    track_b_counts, track_b_meta = g245601_load_pseudobulk(pb["output"]["track_b"]["counts_tsv"], pb["output"]["track_b"]["metadata_tsv"])
    track_a_de = pd.read_csv(pb["output"]["de"]["track_a_genomewide_tsv"], sep="\t")
    track_b_de = pd.read_csv(pb["output"]["de"]["track_b_genomewide_tsv"], sep="\t")

    plot_paired_candidates_grid(track_a_counts, track_a_meta, genes, "GSE245601 Track A (all epithelial, 10 patients) -- acute 12h response", out_dir / "gse245601_track_a_multimodal7_sample_level.png", tested_genes=set(track_a_de["gene"]))
    n_track_b_patients = track_b_meta["patient"].nunique()
    plot_paired_candidates_grid(track_b_counts, track_b_meta, genes, f"GSE245601 Track B (strict malignant, {n_track_b_patients} patients -- small n) -- acute 12h response", out_dir / "gse245601_track_b_multimodal7_sample_level.png", tested_genes=set(track_b_de["gene"]))
    logger.info("wrote GSE245601 track A/B sample-level grids (%d Track B patients)", n_track_b_patients)


def run_sample_level(config_path: str | Path = "config/config.yaml") -> None:
    config = _load_config(config_path)
    genes = config["candidate_adjudication"]["multimodal7"]["genes"]
    out_dir = Path(config["candidate_adjudication"]["output"]["sample_level_dir"])
    out_dir.mkdir(parents=True, exist_ok=True)

    plot_gse118713_sample_level(genes, config, out_dir)
    plot_gse111151_sample_level(genes, config, out_dir)
    plot_gse240112_sample_level(genes, config, out_dir)
    plot_gse245601_sample_level(genes, config, out_dir)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run_sample_level()
