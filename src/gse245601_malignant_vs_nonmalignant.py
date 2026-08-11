"""Malignant vs non-malignant epithelial context check for GSE245601
(gse245601_PREANALYSIS.md section 13, Phase 8). Validation/context only --
never the Control-vs-Tamoxifen treatment question, and candidate-gene
results are never used to modify the frozen InferCNV labels anywhere in
this module. Purpose: do the frozen malignant and non-malignant epithelial
groups show biologically plausible differences? This module does not
claim "messier expression = cancer" -- dispersion/heterogeneity metrics
are reported descriptively, not interpreted causally.

Reads the per-cell descriptive table and the malignant-vs-nonmalignant
pseudobulk written by ``scripts/analysis/gse245601_16_malignant_vs_nonmalignant.R``.
The secondary candidate table (13 genes + PAICS, malignant vs
non-malignant) uses patient as the biological unit (a paired Wilcoxon
signed-rank test across patients, never per-cell) -- with only 5 eligible
patients, most candidates will not reach significance; that is reported
honestly, not worked around.
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
from scipy.stats import wilcoxon  # noqa: E402

from src.gse245601_candidate_extraction import benjamini_hochberg  # noqa: E402
from src.gse245601_pseudobulk_qc import compute_log2cpm  # noqa: E402

logger = logging.getLogger(__name__)

MALIGNANT = "malignant"
NONMALIGNANT = "non-malignant epithelial"
STATUS_COLORS = {MALIGNANT: "#B2182B", NONMALIGNANT: "#2166AC"}


def _load_config(config_path: str | Path) -> dict:
    with open(config_path) as f:
        return yaml.safe_load(f)


def load_cell_level_summary(path: str | Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    required = ("cell_id", "patient", "condition", "malignancy_status", "cnv_score")
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"cell-level summary missing required columns: {missing}")
    bad = set(df["malignancy_status"].unique()) - {MALIGNANT, NONMALIGNANT}
    if bad:
        raise ValueError(f"unexpected malignancy_status values: {bad}")
    logger.info("load_cell_level_summary: %d cells (%d malignant, %d non-malignant)", len(df), (df["malignancy_status"] == MALIGNANT).sum(), (df["malignancy_status"] == NONMALIGNANT).sum())
    return df


def plot_cnv_score_distribution(df: pd.DataFrame, out_path: str | Path) -> None:
    """CNV score distribution by malignancy status -- an internal-
    consistency check (malignant cells were classified USING this score,
    so separation here is expected by construction, not independent
    evidence) rather than a validation of the classifier."""
    fig, ax = plt.subplots(figsize=(6, 4.5))
    for status, color in STATUS_COLORS.items():
        vals = df.loc[df["malignancy_status"] == status, "cnv_score"]
        ax.hist(vals, bins=60, alpha=0.6, color=color, label=f"{status} (n={len(vals)})", density=True)
    ax.set_xlabel("CNV score")
    ax.set_ylabel("density")
    ax.set_title("CNV score distribution by frozen malignancy status\n(expected by construction -- not independent validation)", fontsize=9)
    ax.legend(fontsize=8)
    fig.tight_layout()
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info("wrote %s", out_path)


def plot_umap_separation(df: pd.DataFrame, out_path: str | Path) -> None:
    """UMAP (already computed during the original candidate-blind
    clustering, not recomputed here), epithelial cells only, colored by
    malignancy status."""
    fig, ax = plt.subplots(figsize=(6.5, 5.5))
    for status, color in STATUS_COLORS.items():
        sub = df.loc[df["malignancy_status"] == status]
        ax.scatter(sub["umap_1"], sub["umap_2"], s=3, alpha=0.3, color=color, label=f"{status} (n={len(sub)})", linewidths=0)
    ax.set_xlabel("UMAP_1")
    ax.set_ylabel("UMAP_2")
    ax.set_title("Epithelial cells in existing UMAP space, by frozen malignancy status", fontsize=10)
    legend = ax.legend(fontsize=8, markerscale=3)
    for lh in legend.legend_handles:
        lh.set_alpha(1)
    fig.tight_layout()
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info("wrote %s", out_path)


def plot_broad_programs_and_dispersion(df: pd.DataFrame, out_path: str | Path) -> None:
    """Cell-cycle scores (Seurat's standard CellCycleScoring, established
    field-standard gene sets, not a candidate gene), the existing
    epithelial-program module score, and basic QC-derived dispersion
    metrics (nFeature_RNA, percent_mt), by malignancy status. Described
    descriptively only -- higher/lower dispersion is not interpreted as
    "more cancer-like"."""
    panels = [
        ("S_Score", "S-phase score"),
        ("G2M_Score", "G2M-phase score"),
        ("score_epithelial_program", "epithelial program score"),
        ("nFeature_RNA", "detected genes per cell"),
        ("percent_mt", "% mitochondrial reads"),
    ]
    fig, axes = plt.subplots(1, len(panels), figsize=(3.0 * len(panels), 4.2))
    for ax, (col, title) in zip(axes, panels):
        data = [df.loc[df["malignancy_status"] == s, col].dropna().to_numpy() for s in (NONMALIGNANT, MALIGNANT)]
        bp = ax.boxplot(data, tick_labels=["non-mal.", "malignant"], showfliers=False, patch_artist=True)
        for patch, status in zip(bp["boxes"], (NONMALIGNANT, MALIGNANT)):
            patch.set_facecolor(STATUS_COLORS[status])
        ax.set_title(title, fontsize=8)
        ax.tick_params(axis="x", labelsize=7)
    fig.suptitle("Broad epithelial programs and dispersion metrics by frozen malignancy status", fontsize=10)
    fig.tight_layout()
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info("wrote %s", out_path)


def build_malignant_vs_nonmalignant_candidate_table(
    counts: pd.DataFrame, metadata: pd.DataFrame, candidate_genes: list[str], paics_gene: str
) -> pd.DataFrame:
    """Secondary, exploratory table: paired (by patient) malignant vs
    non-malignant epithelial log2(CPM+1) delta for the 13 candidates +
    PAICS, using the patients present in the malignant-vs-nonmalignant
    pseudobulk (patient is the unit; a paired Wilcoxon signed-rank test is
    used, not a per-cell test). PAICS is flagged as a benchmark and
    excluded from the candidate-set BH family, same convention as the
    Track A/B tables."""
    sample_ids = metadata["sample_id"].tolist()
    log2cpm = compute_log2cpm(counts, sample_ids)
    meta_indexed = metadata.set_index("sample_id")
    patients = sorted(metadata["patient"].unique())

    rows = []
    for gene in candidate_genes + [paics_gene]:
        if gene not in log2cpm.index:
            rows.append(
                {
                    "gene": gene,
                    "is_paics_benchmark": gene == paics_gene,
                    "tested": False,
                    "mean_delta_malignant_minus_nonmalignant": np.nan,
                    "p_value": np.nan,
                    "n_patients": 0,
                    "direction": "not_tested",
                }
            )
            continue
        gene_vals = log2cpm.loc[gene]
        deltas = []
        for patient in patients:
            psamples = meta_indexed.loc[meta_indexed["patient"] == patient]
            if not {"malignant", "nonmalignant"}.issubset(set(psamples["malignancy_status"])):
                continue
            mal_id = psamples.loc[psamples["malignancy_status"] == "malignant"].index[0]
            nonmal_id = psamples.loc[psamples["malignancy_status"] == "nonmalignant"].index[0]
            deltas.append(gene_vals[mal_id] - gene_vals[nonmal_id])
        deltas = np.array(deltas)
        if len(deltas) >= 2 and np.any(deltas != 0):
            try:
                # Parameters pinned explicitly (matching scipy's current
                # defaults) so behavior does not silently drift with a
                # future scipy version: two-sided exact/auto signed-rank
                # test, zeros dropped before ranking (Wilcoxon's original
                # convention).
                _, p_value = wilcoxon(deltas, zero_method="wilcox", alternative="two-sided", method="auto")
            except ValueError:
                p_value = np.nan
        else:
            p_value = np.nan
        mean_delta = deltas.mean() if len(deltas) > 0 else np.nan
        rows.append(
            {
                "gene": gene,
                "is_paics_benchmark": gene == paics_gene,
                "tested": True,
                "mean_delta_malignant_minus_nonmalignant": mean_delta,
                "p_value": p_value,
                "n_patients": len(deltas),
                "direction": "up_in_malignant" if mean_delta > 0 else ("down_in_malignant" if mean_delta < 0 else "unchanged"),
            }
        )

    out = pd.DataFrame(rows)
    candidate_mask = ~out["is_paics_benchmark"] & out["tested"] & out["p_value"].notna()
    out["candidate_set_bh_fdr"] = np.nan
    if candidate_mask.any():
        out.loc[candidate_mask, "candidate_set_bh_fdr"] = benjamini_hochberg(out.loc[candidate_mask, "p_value"])
    logger.info(
        "build_malignant_vs_nonmalignant_candidate_table: %d/%d candidates tested (n_patients=%s)",
        int((out["tested"] & ~out["is_paics_benchmark"]).sum()),
        len(candidate_genes),
        sorted(out.loc[~out["is_paics_benchmark"], "n_patients"].unique().tolist()),
    )
    return out


def run_malignant_vs_nonmalignant(config_path: str | Path = "config/config.yaml") -> dict[str, pd.DataFrame]:
    config = _load_config(config_path)
    pb_cfg = config["gse245601_pseudobulk"]
    ci_cfg = config["gse245601_candidate_integration"]
    candidates_13 = pb_cfg["candidates"]["thirteen"]
    paics_gene = pb_cfg["candidates"]["paics"]

    cell_level = load_cell_level_summary(pb_cfg["output"]["malignant_vs_nonmalignant"]["cell_level_summary_tsv"])
    figures_dir = Path(pb_cfg["output"]["malignant_vs_nonmalignant"]["figures_dir"])

    plot_cnv_score_distribution(cell_level, figures_dir / "malignant_vs_nonmalignant_cnv_score.png")
    plot_umap_separation(cell_level, figures_dir / "malignant_vs_nonmalignant_umap.png")
    plot_broad_programs_and_dispersion(cell_level, figures_dir / "malignant_vs_nonmalignant_programs_dispersion.png")

    counts = pd.read_csv(pb_cfg["output"]["malignant_vs_nonmalignant"]["counts_tsv"], sep="\t").set_index("gene")
    metadata = pd.read_csv(pb_cfg["output"]["malignant_vs_nonmalignant"]["metadata_tsv"], sep="\t")
    candidate_table = build_malignant_vs_nonmalignant_candidate_table(counts, metadata, candidates_13, paics_gene)

    out_path = Path(ci_cfg["output"]["malignant_vs_nonmalignant_candidates_tsv"])
    out_path.parent.mkdir(parents=True, exist_ok=True)
    candidate_table.to_csv(out_path, sep="\t", index=False)
    logger.info("run_malignant_vs_nonmalignant: wrote %s", out_path)

    return {"cell_level": cell_level, "candidate_table": candidate_table}


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run_malignant_vs_nonmalignant()
