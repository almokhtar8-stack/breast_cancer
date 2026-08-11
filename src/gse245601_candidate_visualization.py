"""Patient-level candidate visualizations for GSE245601 pseudobulk Track A
/ Track B (gse245601_PREANALYSIS.md section 13). All figures are built
directly from the already-written pseudobulk counts and candidate tables
-- no new statistic is computed here beyond simple per-patient log2(CPM+1)
deltas already used for the patient-direction-consistency columns in
``src.gse245601_candidate_extraction``. PAICS is always drawn as a
visually separate, clearly labeled benchmark row -- never merged into the
13-candidate ranking or ordering.
"""

from __future__ import annotations

import logging
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from src.gse245601_candidate_extraction import compute_patient_direction  # noqa: E402
from src.gse245601_pseudobulk_qc import PATIENT_COLORS  # noqa: E402

logger = logging.getLogger(__name__)


def build_patient_delta_matrix(
    counts: pd.DataFrame, metadata: pd.DataFrame, genes: list[str], tested_genes: set[str] | None = None
) -> pd.DataFrame:
    """genes x patients matrix of per-patient (Tamoxifen - Control)
    log2(CPM+1) delta. Genes not present in the count matrix, OR not in
    ``tested_genes`` (e.g. filtered out by edgeR's filterByExpr for being
    too lowly expressed to trust), get an all-NaN row -- a near-zero delta
    for an untested, unreliable gene must never render identically to a
    genuinely-measured near-zero delta."""
    all_patients = sorted(metadata["patient"].unique())
    rows = {}
    for gene in genes:
        if tested_genes is not None and gene not in tested_genes:
            rows[gene] = [np.nan] * len(all_patients)
            continue
        patient_dir = compute_patient_direction(counts, metadata, gene)
        delta_by_patient = dict(zip(patient_dir["patient"], patient_dir["delta"]))
        rows[gene] = [delta_by_patient.get(p, np.nan) for p in all_patients]
    return pd.DataFrame(rows, index=all_patients).T


def plot_paired_candidates_grid(
    counts: pd.DataFrame,
    metadata: pd.DataFrame,
    candidate_genes: list[str],
    track_label: str,
    out_path: str | Path,
    tested_genes: set[str] | None = None,
) -> None:
    """Compact small-multiples grid: one panel per candidate gene,
    Control->Tamoxifen paired line per patient (color = patient). Genes
    not in ``tested_genes`` (filtered by edgeR as too lowly expressed) are
    drawn as "not tested" rather than plotting an unreliable raw-count
    trend line."""
    from src.gse245601_pseudobulk_qc import compute_log2cpm

    sample_ids = metadata["sample_id"].tolist()
    log2cpm = compute_log2cpm(counts, sample_ids)
    meta_indexed = metadata.set_index("sample_id")
    patients_sorted = sorted(metadata["patient"].unique())

    n = len(candidate_genes)
    n_cols = 4
    n_rows = int(np.ceil(n / n_cols))
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(3.0 * n_cols, 2.6 * n_rows), squeeze=False)

    for i, gene in enumerate(candidate_genes):
        ax = axes[i // n_cols, i % n_cols]
        gene_untested = gene not in log2cpm.index or (tested_genes is not None and gene not in tested_genes)
        if gene_untested:
            ax.text(0.5, 0.5, "not tested\n(filtered)", ha="center", va="center", fontsize=8, transform=ax.transAxes)
            ax.set_title(gene, fontsize=9)
            ax.set_xticks([])
            ax.set_yticks([])
            continue
        gene_vals = log2cpm.loc[gene]
        for patient in patients_sorted:
            patient_samples = meta_indexed.loc[meta_indexed["patient"] == patient]
            if not {"Control", "Tamoxifen"}.issubset(set(patient_samples["condition"])):
                continue
            control_id = patient_samples.loc[patient_samples["condition"] == "Control"].index[0]
            tam_id = patient_samples.loc[patient_samples["condition"] == "Tamoxifen"].index[0]
            color = PATIENT_COLORS[patients_sorted.index(patient) % len(PATIENT_COLORS)]
            ax.plot([0, 1], [gene_vals[control_id], gene_vals[tam_id]], color=color, marker="o", markersize=3, linewidth=1, alpha=0.8)
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["Ctrl", "Tam"], fontsize=7)
        ax.set_title(gene, fontsize=9)
        ax.tick_params(axis="y", labelsize=6)

    for j in range(n, n_rows * n_cols):
        axes[j // n_cols, j % n_cols].axis("off")

    handles = [plt.Line2D([0], [0], color=PATIENT_COLORS[i % len(PATIENT_COLORS)], label=p) for i, p in enumerate(patients_sorted)]
    fig.legend(handles=handles, loc="upper center", ncol=len(patients_sorted), fontsize=7, frameon=False, bbox_to_anchor=(0.5, 1.03))
    fig.suptitle(f"{track_label}: paired candidate log2(CPM+1), Control vs Tamoxifen", fontsize=11, y=1.06)
    fig.tight_layout()
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info("wrote %s", out_path)


def plot_candidate_effect_heatmap(delta_matrix: pd.DataFrame, track_label: str, out_path: str | Path) -> None:
    """Rows = candidates, columns = patients, value = paired
    Tamoxifen-Control log2(CPM+1) delta. Diverging colormap centered at 0."""
    arr = delta_matrix.to_numpy(dtype=float)
    vmax = np.nanmax(np.abs(arr)) if np.isfinite(arr).any() else 1.0
    fig, ax = plt.subplots(figsize=(0.9 * delta_matrix.shape[1] + 2, 0.35 * delta_matrix.shape[0] + 1.5))
    im = ax.imshow(arr, cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
    ax.set_xticks(range(delta_matrix.shape[1]))
    ax.set_xticklabels(delta_matrix.columns, rotation=90, fontsize=8)
    ax.set_yticks(range(delta_matrix.shape[0]))
    ax.set_yticklabels(delta_matrix.index, fontsize=8)
    for i in range(arr.shape[0]):
        for j in range(arr.shape[1]):
            if not np.isnan(arr[i, j]):
                ax.text(j, i, f"{arr[i, j]:.2f}", ha="center", va="center", fontsize=6)
            else:
                ax.text(j, i, "N/A", ha="center", va="center", fontsize=6, color="gray")
    fig.colorbar(im, ax=ax, label="paired delta, log2(CPM+1)\n(Tamoxifen - Control)")
    ax.set_title(f"{track_label}: candidate response heatmap", fontsize=10)
    fig.tight_layout()
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info("wrote %s", out_path)


def plot_track_comparison_summary(
    track_a_candidates: pd.DataFrame,
    track_b_candidates: pd.DataFrame,
    track_a_delta: pd.DataFrame,
    track_b_delta: pd.DataFrame,
    paics_row_a: pd.DataFrame,
    paics_row_b: pd.DataFrame,
    candidate_genes: list[str],
    out_path: str | Path,
) -> None:
    """Point-range ("forest-style") plot: one row per gene (13 candidates,
    ordered by Track A p-value, plus PAICS drawn last and visually
    separated by a divider), Track A and Track B log2FC shown as two
    offset points with error bars = SEM of the per-patient deltas that
    feed each track's direction-consistency check (an empirical,
    transparent uncertainty proxy -- not edgeR's internal GLM SE)."""
    order = track_a_candidates.sort_values("p_value")["gene"].tolist()
    all_rows = order + ["PAICS"]

    def _stats(cand_df, delta_df, gene):
        row = cand_df.loc[cand_df["gene"] == gene]
        if len(row) == 0 or not bool(row.iloc[0]["tested"]):
            return np.nan, np.nan
        log2fc = float(row.iloc[0]["log2fc"])
        if gene in delta_df.index:
            deltas = delta_df.loc[gene].dropna().to_numpy()
            sem = deltas.std(ddof=1) / np.sqrt(len(deltas)) if len(deltas) > 1 else np.nan
        else:
            sem = np.nan
        return log2fc, sem

    fig, ax = plt.subplots(figsize=(6.5, 0.4 * len(all_rows) + 1.5))
    y_positions = np.arange(len(all_rows))[::-1]
    for y, gene in zip(y_positions, all_rows):
        if gene == "PAICS":
            a_fc = float(paics_row_a["log2fc"].iloc[0]) if paics_row_a["tested"].iloc[0] else np.nan
            b_fc = float(paics_row_b["log2fc"].iloc[0]) if paics_row_b["tested"].iloc[0] else np.nan
            a_sem = b_sem = np.nan
        else:
            a_fc, a_sem = _stats(track_a_candidates, track_a_delta, gene)
            b_fc, b_sem = _stats(track_b_candidates, track_b_delta, gene)

        if not np.isnan(a_fc):
            ax.errorbar(a_fc, y + 0.15, xerr=a_sem, fmt="o", color="#2166AC", markersize=5, capsize=3, label="Track A" if y == y_positions[0] else None)
        if not np.isnan(b_fc):
            ax.errorbar(b_fc, y - 0.15, xerr=b_sem, fmt="s", color="#B2182B", markersize=5, capsize=3, label="Track B" if y == y_positions[0] else None)

    ax.axvline(0, color="black", linewidth=0.6)
    ax.axhline(y_positions[-1] + 0.5, color="gray", linewidth=0.8, linestyle="--")  # separates PAICS
    ax.set_yticks(y_positions)
    ax.set_yticklabels(all_rows)
    ax.set_xlabel("log2FC (Tamoxifen vs Control)")
    ax.set_title("Track A vs Track B candidate effects (PAICS = benchmark, below divider)", fontsize=10)
    ax.legend(loc="best", fontsize=8)
    fig.tight_layout()
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info("wrote %s", out_path)


def run_candidate_visualization(config_path: str | Path = "config/config.yaml") -> None:
    import yaml

    from src.gse245601_pseudobulk_qc import load_pseudobulk

    with open(config_path) as f:
        config = yaml.safe_load(f)
    pb_cfg = config["gse245601_pseudobulk"]
    ci_cfg = config["gse245601_candidate_integration"]
    candidates_13 = pb_cfg["candidates"]["thirteen"]
    paics_gene = pb_cfg["candidates"]["paics"]
    figures_dir = Path(ci_cfg["output"]["figures_dir"])

    track_a_counts, track_a_meta = load_pseudobulk(pb_cfg["output"]["track_a"]["counts_tsv"], pb_cfg["output"]["track_a"]["metadata_tsv"])
    track_b_counts, track_b_meta = load_pseudobulk(pb_cfg["output"]["track_b"]["counts_tsv"], pb_cfg["output"]["track_b"]["metadata_tsv"])
    track_a_candidates = pd.read_csv(ci_cfg["output"]["track_a_candidates_tsv"], sep="\t")
    track_b_candidates = pd.read_csv(ci_cfg["output"]["track_b_candidates_tsv"], sep="\t")
    paics = pd.read_csv(ci_cfg["output"]["paics_singlecell_tsv"], sep="\t")
    paics_row_a = paics.loc[paics["track"] == "track_a_epithelial"]
    paics_row_b = paics.loc[paics["track"] == "track_b_malignant"]

    track_a_tested = set(track_a_candidates.loc[track_a_candidates["tested"], "gene"])
    track_b_tested = set(track_b_candidates.loc[track_b_candidates["tested"], "gene"])

    plot_paired_candidates_grid(
        track_a_counts, track_a_meta, candidates_13, "track_a_epithelial", figures_dir / "track_a_paired_candidates_grid.png", tested_genes=track_a_tested
    )
    plot_paired_candidates_grid(
        track_b_counts, track_b_meta, candidates_13, "track_b_malignant", figures_dir / "track_b_paired_candidates_grid.png", tested_genes=track_b_tested
    )

    track_a_delta = build_patient_delta_matrix(track_a_counts, track_a_meta, candidates_13, tested_genes=track_a_tested)
    track_b_delta = build_patient_delta_matrix(track_b_counts, track_b_meta, candidates_13, tested_genes=track_b_tested)
    plot_candidate_effect_heatmap(track_a_delta, "track_a_epithelial", figures_dir / "track_a_candidate_heatmap.png")
    plot_candidate_effect_heatmap(track_b_delta, "track_b_malignant", figures_dir / "track_b_candidate_heatmap.png")

    plot_track_comparison_summary(
        track_a_candidates, track_b_candidates, track_a_delta, track_b_delta, paics_row_a, paics_row_b, candidates_13,
        figures_dir / "track_a_vs_track_b_summary.png"
    )


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run_candidate_visualization()
