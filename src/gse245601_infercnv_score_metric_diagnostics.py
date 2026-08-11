"""Diagnostics for whether InferCNV's whole-genome mean-squared CNV score
(``cnv_score = mean((expr-1)^2)`` across every gene InferCNV output) is an
appropriate malignancy metric when abnormalities affect only part of the
genome (docs/CNV_METHOD_AUDIT.md Point 1 follow-up).

This module reruns nothing and changes no label, threshold, or classifier.
It reads exactly ``results/tables/gse245601_infercnv_score_metric_diagnostics.tsv``,
written by ``scripts/analysis/gse245601_11_extract_cnv_score_metric_diagnostics.R``
(itself a read-only reload of the already-frozen ``run.final.infercnv_obj``,
not a rerun of InferCNV -- that script cross-checks its recomputed
``cnv_score`` against the frozen one for every cell before writing
anything, so this module can trust the two agree).

Every metric here is a diagnostic MEASUREMENT, never a proposed
replacement classifier: fraction of genes exceeding a deviation level,
upper-percentile deviation amplitude, and chromosome-level burden. CopyKAT
calls are used only as a plotting annotation (marker shape) -- never to
choose, tune, or validate any of these metrics or any threshold.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.patches as mpatches  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import yaml  # noqa: E402

from src.gse245601_infercnv_threshold_diagnostics import (  # noqa: E402
    COLOR_MALIGNANT,
    COLOR_NONMALIGNANT,
    COPYKAT_MARKERS,
    FAIL_BOTH,
    FAIL_CNV_ONLY,
    FAIL_CORR_ONLY,
    PASS_BOTH,
    ThresholdDiagnosticsConfig,
    classify_failure_category,
    load_cnv_score_table,
    recompute_group_thresholds,
)
from src.gse245601_cnv_method_comparison import LABEL_MALIGNANT, LABEL_NONMALIGNANT  # noqa: E402

logger = logging.getLogger(__name__)

GROUP_DISAGREEMENT = "disagreement"
GROUP_GOOD_CONTROL = "good_control"
GROUP_SIGNAL_QUALITY = "signal_quality"

REQUIRED_COLUMNS = (
    "cell_id",
    "sample_id",
    "patient",
    "condition",
    "primary_malignancy_label",
    "cnv_correlation_to_seed",
    "cnv_score",
    "cnv_score_recomputed",
    "max_abs_deviation",
    "p95_abs_deviation",
    "p99_abs_deviation",
)


def _load_config(config_path: str | Path) -> dict:
    with open(config_path) as f:
        return yaml.safe_load(f)


@dataclass(frozen=True)
class ScoreMetricConfig:
    """Resolved, config-driven paths and constants. No hardcoded values."""

    selected_samples: tuple[str, ...]
    disagreement_tumors: tuple[str, ...]
    good_control_tumors: tuple[str, ...]
    signal_quality_tumors: tuple[str, ...]
    infercnv_labels_tsv: Path
    copykat_labels_tsv: Path
    deviation_levels: tuple[float, ...]
    diagnostics_tsv: Path
    figures_dir: Path

    @classmethod
    def from_config(cls, config: dict) -> "ScoreMetricConfig":
        cfg = config["gse245601_infercnv_score_metric_diagnostics"]
        inputs = cfg["inputs"]
        out = cfg["output"]
        return cls(
            selected_samples=tuple(cfg["selected_samples"]),
            disagreement_tumors=tuple(cfg["disagreement_tumors"]),
            good_control_tumors=tuple(cfg["good_control_tumors"]),
            signal_quality_tumors=tuple(cfg["signal_quality_tumors"]),
            infercnv_labels_tsv=Path(inputs["infercnv_labels_tsv"]),
            copykat_labels_tsv=Path(inputs["copykat_labels_tsv"]),
            deviation_levels=tuple(float(x) for x in cfg["deviation_levels"]),
            diagnostics_tsv=Path(out["diagnostics_tsv"]),
            figures_dir=Path(out["figures_dir"]),
        )

    def threshold_cfg(self) -> ThresholdDiagnosticsConfig:
        """Reconstructs the ThresholdDiagnosticsConfig this module needs
        (only the pieces used by ``recompute_group_thresholds``) so
        failure categories can be recomputed here without duplicating the
        rule constants a second time."""
        return ThresholdDiagnosticsConfig(
            selected_samples=self.selected_samples,
            disagreement_tumors=self.disagreement_tumors,
            good_control_tumors=self.good_control_tumors,
            signal_quality_tumors=self.signal_quality_tumors,
            infercnv_labels_tsv=self.infercnv_labels_tsv,
            copykat_labels_tsv=self.copykat_labels_tsv,
            seed_top_fraction=0.05,
            min_seed_cells=2,
            cnv_score_sd_multiplier=2.0,
            cnv_score_clamp=(0.01, 0.05),
            correlation_sd_multiplier=1.5,
            correlation_clamp=(0.2, 0.4),
            cnv_score_sweep_range=(0.0, 0.05),
            cnv_score_grid_points=1,
            correlation_sweep_range=(0.2, 0.4),
            correlation_grid_points=1,
            diagnostics_tsv=Path("unused"),
            sensitivity_grid_tsv=Path("unused"),
            local_score_sensitivity_tsv=Path("unused"),
            figures_dir=Path("unused"),
        )


def tumor_group_for(patient: str, cfg: ScoreMetricConfig) -> str:
    """Classifies a tumor into exactly one of the three pre-declared
    groups from the prior threshold-diagnostics phase. Raises on an
    unrecognized patient rather than silently leaving it ungrouped."""
    if patient in cfg.disagreement_tumors:
        return GROUP_DISAGREEMENT
    if patient in cfg.good_control_tumors:
        return GROUP_GOOD_CONTROL
    if patient in cfg.signal_quality_tumors:
        return GROUP_SIGNAL_QUALITY
    raise ValueError(f"tumor_group_for: {patient!r} is not in any configured tumor group")


def load_score_metric_table(path: str | Path) -> pd.DataFrame:
    """Load the frozen per-cell CNV signal-extent table written by
    gse245601_11_extract_cnv_score_metric_diagnostics.R. Read-only -- that
    script already cross-checked cnv_score_recomputed against the frozen
    cnv_score for every cell; this loader re-checks that invariant rather
    than assuming the file on disk still satisfies it."""
    df = pd.read_csv(path, sep="\t")
    missing = [c for c in REQUIRED_COLUMNS if c not in df.columns]
    if missing:
        raise ValueError(f"score metric table missing required columns: {missing}")
    diff = (df["cnv_score"] - df["cnv_score_recomputed"]).abs()
    if (diff > 1e-6).any():
        raise ValueError(
            f"cnv_score and cnv_score_recomputed disagree for {(diff > 1e-6).sum()} cell(s) "
            f"(max diff={diff.max():.3g}) -- this table should not be trusted"
        )
    logger.info("load_score_metric_table: read %d cells across %d samples", len(df), df["sample_id"].nunique())
    return df


def build_group_comparison_summary(metric_df: pd.DataFrame, cfg: ScoreMetricConfig) -> pd.DataFrame:
    """One row per (patient, condition): median/range of cnv_score, the
    genomic-extent fraction at each deviation level, upper-percentile
    amplitude, and chromosome burden -- the numbers behind Point 2's
    magnitude/extent/chromosome-count comparison. Adds the tumor's
    pre-declared group (disagreement/good_control/signal_quality) as a
    plain label, not a derived classification."""
    frac_cols = [c for c in metric_df.columns if c.startswith("fraction_genes_dev_gt_")]
    n_chr_cols = [c for c in metric_df.columns if c.startswith("n_chromosomes_affected_")]

    rows = []
    for (patient, condition), sub in metric_df.groupby(["patient", "condition"]):
        row = {
            "patient": patient,
            "condition": condition,
            "tumor_group": tumor_group_for(patient, cfg),
            "n_cells": len(sub),
            "median_cnv_score": sub["cnv_score"].median(),
            "median_max_abs_deviation": sub["max_abs_deviation"].median(),
            "median_p95_abs_deviation": sub["p95_abs_deviation"].median(),
            "median_p99_abs_deviation": sub["p99_abs_deviation"].median(),
        }
        for c in frac_cols:
            row[f"median_{c}"] = sub[c].median()
        for c in n_chr_cols:
            row[f"median_{c}"] = sub[c].median()
        rows.append(row)

    out = pd.DataFrame(rows)
    group_order = {GROUP_DISAGREEMENT: 0, GROUP_GOOD_CONTROL: 1, GROUP_SIGNAL_QUALITY: 2}
    out["_group_order"] = out["tumor_group"].map(group_order)
    out = out.sort_values(["_group_order", "patient", "condition"]).drop(columns="_group_order").reset_index(drop=True)
    logger.info("build_group_comparison_summary: %d (patient, condition) rows", len(out))
    return out


def compute_per_cell_failure_category(cfg: ScoreMetricConfig) -> pd.DataFrame:
    """Recomputes the frozen classifier's per-cell failure category
    (reusing the exact, already-verified logic from
    ``src.gse245601_infercnv_threshold_diagnostics`` -- not reimplemented
    here) for every selected cell. Used only to annotate WHY a cell was
    rejected when investigating Tumor_10 (Point 4); never used to
    redefine or re-derive the frozen labels."""
    tcfg = cfg.threshold_cfg()
    cell_df = load_cnv_score_table(tcfg.infercnv_labels_tsv)
    sub = cell_df.loc[cell_df["sample_id"].isin(cfg.selected_samples)].copy()

    categories = np.empty(len(sub), dtype=object)
    for (sample_id, group), grp_df in sub.groupby(["sample_id", "threshold_group"]):
        th_value, th_corr = recompute_group_thresholds(grp_df["cnv_score"], grp_df["cnv_correlation_to_seed"], tcfg)
        idx = grp_df.index
        cats = [
            classify_failure_category(s, c, th_value, th_corr)
            for s, c in zip(grp_df["cnv_score"], grp_df["cnv_correlation_to_seed"])
        ]
        categories[sub.index.get_indexer(idx)] = cats

    sub["failure_category"] = categories
    logger.info("compute_per_cell_failure_category: %d cells categorized", len(sub))
    return sub[["cell_id", "sample_id", "failure_category"]]


def build_tumor10_diagnostic_table(metric_df: pd.DataFrame, cfg: ScoreMetricConfig) -> pd.DataFrame:
    """Merges CNV-extent metrics with the recomputed failure category for
    Tumor_10 cells only -- to determine whether Tumor_10's low agreement
    is driven by weak CNV score, poor correlation to seed, or a broader
    (genuinely noisier) CNV structure, without assuming it matches the
    Tumor_01/04/08 pattern."""
    failure_df = compute_per_cell_failure_category(cfg)
    sub = metric_df.loc[metric_df["patient"].isin(cfg.signal_quality_tumors)].copy()
    out = sub.merge(failure_df[["cell_id", "failure_category"]], on="cell_id", how="left", validate="one_to_one")
    if out["failure_category"].isna().any():
        raise ValueError("build_tumor10_diagnostic_table: failure_category missing for some Tumor_10 cells")
    logger.info("build_tumor10_diagnostic_table: %d Tumor_10 cells", len(out))
    return out


def _annotated_copykat(metric_df: pd.DataFrame, copykat_df: pd.DataFrame | None) -> pd.DataFrame:
    if copykat_df is None:
        out = metric_df.copy()
        out["sensitivity_malignancy_label"] = "not_defined"
        return out
    ck = copykat_df[["cell_id", "sensitivity_malignancy_label"]]
    return metric_df.merge(ck, on="cell_id", how="left")


def _scatter_grid(
    metric_df: pd.DataFrame, x_col: str, x_label: str, cfg: ScoreMetricConfig, out_path: str | Path, hline: float | None = 0.01
) -> None:
    fig, axes = plt.subplots(6, 2, figsize=(9, 20), squeeze=False)
    for i, patient in enumerate(("Tumor_01", "Tumor_04", "Tumor_08", "Tumor_02", "Tumor_03", "Tumor_10")):
        for j, condition in enumerate(("Control", "Tamoxifen")):
            ax = axes[i, j]
            sub = metric_df.loc[(metric_df["patient"] == patient) & (metric_df["condition"] == condition)]
            for label, color in ((LABEL_MALIGNANT, COLOR_MALIGNANT), (LABEL_NONMALIGNANT, COLOR_NONMALIGNANT)):
                for ck_label, marker in COPYKAT_MARKERS.items():
                    pts = sub.loc[(sub["primary_malignancy_label"] == label) & (sub["sensitivity_malignancy_label"] == ck_label)]
                    if len(pts) == 0:
                        continue
                    ax.scatter(pts[x_col], pts["cnv_score"], c=color, marker=marker, s=10, alpha=0.5, linewidths=0)
            if hline is not None:
                ax.axhline(hline, color="black", linewidth=0.7, linestyle="--")
            ax.set_title(f"{patient}_{condition} (n={len(sub)})", fontsize=8)
            ax.set_xlabel(x_label, fontsize=7)
            ax.set_ylabel("CNV score", fontsize=7)
            ax.tick_params(labelsize=6)
    fig.tight_layout()
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info("wrote %s", out_path)


def plot_score_vs_extent_grid(metric_df: pd.DataFrame, copykat_df: pd.DataFrame | None, cfg: ScoreMetricConfig, out_path: str | Path) -> None:
    """CNV score vs. genomic extent (fraction of genes with |dev|>0.10),
    one panel per (tumor, condition). If cells with strong but localized
    CNVs are being diluted, they should show up as points with LOW extent
    but comparatively HIGH max/p99 deviation (see the companion plot) --
    this panel shows the score/extent relationship itself."""
    annotated = _annotated_copykat(metric_df, copykat_df)
    _scatter_grid(annotated, "fraction_genes_dev_gt_0.1", "fraction of genes with |dev|>0.10 (genomic extent)", cfg, out_path)


def plot_score_vs_upper_percentile_grid(metric_df: pd.DataFrame, copykat_df: pd.DataFrame | None, cfg: ScoreMetricConfig, out_path: str | Path) -> None:
    """CNV score vs. p99 deviation amplitude, one panel per (tumor,
    condition). A cell with a high p99 (its most-affected genes deviate
    strongly) but a low whole-genome score is a direct, per-cell instance
    of dilution: strong localized signal, low whole-genome average."""
    annotated = _annotated_copykat(metric_df, copykat_df)
    _scatter_grid(annotated, "p99_abs_deviation", "p99 |deviation| across genes (upper-percentile amplitude)", cfg, out_path)


def plot_per_sample_distributions(metric_df: pd.DataFrame, cfg: ScoreMetricConfig, out_path: str | Path) -> None:
    """Two rows of box plots (CNV score; genomic extent), one box per
    sample, ordered and colored by tumor group -- a compact view of
    whether disagreement tumors differ from good controls in overall
    distribution, not just in a handful of example cells."""
    group_colors = {GROUP_DISAGREEMENT: "#B2182B", GROUP_GOOD_CONTROL: "#2166AC", GROUP_SIGNAL_QUALITY: "#969696"}
    order = []
    for tumor_list in (cfg.disagreement_tumors, cfg.good_control_tumors, cfg.signal_quality_tumors):
        for t in tumor_list:
            for cond in ("Control", "Tamoxifen"):
                order.append(f"{t}_{cond}")

    fig, axes = plt.subplots(2, 1, figsize=(11, 7))
    for ax, col, title in zip(axes, ("cnv_score", "fraction_genes_dev_gt_0.1"), ("CNV score", "genomic extent (frac. genes |dev|>0.10)")):
        data = [metric_df.loc[metric_df["sample_id"] == s, col].to_numpy() for s in order]
        bp = ax.boxplot(data, tick_labels=order, showfliers=False, patch_artist=True)
        for patch, sample_id in zip(bp["boxes"], order):
            patient = sample_id.rsplit("_", 1)[0]
            patch.set_facecolor(group_colors[tumor_group_for(patient, cfg)])
        ax.set_title(title, fontsize=10)
        ax.tick_params(axis="x", labelrotation=90, labelsize=7)
        if col == "cnv_score":
            ax.axhline(0.01, color="black", linewidth=0.7, linestyle="--")

    legend_handles = [mpatches.Patch(color=c, label=g) for g, c in group_colors.items()]
    fig.legend(handles=legend_handles, loc="upper center", ncol=3, fontsize=8, frameon=False, bbox_to_anchor=(0.5, 1.02))
    fig.tight_layout()
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info("wrote %s", out_path)


def plot_chromosome_burden_heatmap_grid(metric_df: pd.DataFrame, cfg: ScoreMetricConfig, out_path: str | Path) -> None:
    """One small heatmap per (tumor, condition): rows = cells (sorted by
    CNV score descending), columns = chr1..chr22, color = that
    chromosome's own mean |deviation| for that cell. Shows how many
    chromosomes carry elevated signal per cell -- the direct visual for
    "genomic extent" and "number of affected chromosomes"."""
    chr_cols = [f"chr{k}_mean_abs_dev" for k in range(1, 23)]
    fig, axes = plt.subplots(6, 2, figsize=(9, 24), squeeze=False, gridspec_kw={"hspace": 0.6})
    vmax = float(np.nanpercentile(metric_df[chr_cols].to_numpy(), 99))
    for i, patient in enumerate(("Tumor_01", "Tumor_04", "Tumor_08", "Tumor_02", "Tumor_03", "Tumor_10")):
        for j, condition in enumerate(("Control", "Tamoxifen")):
            ax = axes[i, j]
            sub = metric_df.loc[(metric_df["patient"] == patient) & (metric_df["condition"] == condition)].sort_values(
                "cnv_score", ascending=False
            )
            arr = sub[chr_cols].to_numpy()
            im = ax.imshow(arr, aspect="auto", cmap="viridis", vmin=0, vmax=vmax, interpolation="none")
            ax.set_yticks([])
            ax.set_xticks(range(0, 22, 3))
            ax.set_xticklabels(range(1, 23, 3), fontsize=6)
            ax.set_title(f"{patient}_{condition} (n={len(sub)})", fontsize=8)
            ax.set_xlabel("chromosome", fontsize=7)
    fig.colorbar(im, ax=axes, fraction=0.02, pad=0.02, label="mean |deviation| per chromosome")
    fig.subplots_adjust(hspace=0.6)
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info("wrote %s", out_path)


def plot_tumor10_diagnostic(tumor10_df: pd.DataFrame, copykat_df: pd.DataFrame | None, out_path: str | Path) -> None:
    """Two panels for Tumor_10 cells only, colored by the recomputed
    failure category (not tumor group): CNV score vs. correlation to seed
    (the two axes the classifier actually uses), and genomic extent vs.
    correlation to seed (to check whether Tumor_10's problem is more about
    a noisier/less-coherent CNV profile than about score magnitude)."""
    annotated = _annotated_copykat(tumor10_df, copykat_df)
    cat_colors = {PASS_BOTH: COLOR_MALIGNANT, FAIL_CNV_ONLY: "#F1A340", FAIL_CORR_ONLY: "#66C2A5", FAIL_BOTH: COLOR_NONMALIGNANT}

    fig, axes = plt.subplots(1, 2, figsize=(11, 5))
    for cat, color in cat_colors.items():
        pts = annotated.loc[annotated["failure_category"] == cat]
        axes[0].scatter(pts["cnv_score"], pts["cnv_correlation_to_seed"], c=color, s=10, alpha=0.6, label=cat, linewidths=0)
        axes[1].scatter(pts["fraction_genes_dev_gt_0.1"], pts["cnv_correlation_to_seed"], c=color, s=10, alpha=0.6, label=cat, linewidths=0)
    axes[0].axvline(0.01, color="black", linewidth=0.7, linestyle="--")
    axes[0].set_xlabel("CNV score")
    axes[0].set_ylabel("Kendall correlation to seed")
    axes[0].set_title("Tumor_10: score vs. correlation")
    axes[1].set_xlabel("genomic extent (frac. genes |dev|>0.10)")
    axes[1].set_ylabel("Kendall correlation to seed")
    axes[1].set_title("Tumor_10: extent vs. correlation")
    axes[1].legend(fontsize=7, loc="best")
    fig.tight_layout()
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info("wrote %s", out_path)


def run_score_metric_diagnostics(config_path: str | Path = "config/config.yaml") -> dict[str, object]:
    """Full pipeline: load the frozen metric table (built by the R
    extraction script -- not run automatically here, since it requires the
    sc245601 env and is comparatively slow; run it manually first if the
    table is stale), build the group comparison summary and the Tumor_10
    diagnostic table, and render every figure."""
    config = _load_config(config_path)
    cfg = ScoreMetricConfig.from_config(config)

    metric_df = load_score_metric_table(cfg.diagnostics_tsv)
    copykat_df = pd.read_csv(cfg.copykat_labels_tsv, sep="\t") if cfg.copykat_labels_tsv.exists() else None

    group_summary = build_group_comparison_summary(metric_df, cfg)
    tumor10_df = build_tumor10_diagnostic_table(metric_df, cfg)

    plot_score_vs_extent_grid(metric_df, copykat_df, cfg, cfg.figures_dir / "score_vs_extent_grid.png")
    plot_score_vs_upper_percentile_grid(metric_df, copykat_df, cfg, cfg.figures_dir / "score_vs_upper_percentile_grid.png")
    plot_per_sample_distributions(metric_df, cfg, cfg.figures_dir / "per_sample_score_and_extent_distributions.png")
    plot_chromosome_burden_heatmap_grid(metric_df, cfg, cfg.figures_dir / "chromosome_burden_heatmap_grid.png")
    plot_tumor10_diagnostic(tumor10_df, copykat_df, cfg.figures_dir / "tumor10_diagnostic.png")

    return {"metric_df": metric_df, "group_summary": group_summary, "tumor10_df": tumor10_df}


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run_score_metric_diagnostics()
