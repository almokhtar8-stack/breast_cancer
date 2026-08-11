"""Diagnostics for the InferCNV malignant-cell classification rule
(PREANALYSIS.md section 9 / amendment 3), for the 12 selected GSE245601
samples (docs/CNV_METHOD_AUDIT.md Point 1 follow-up).

This module reruns nothing and changes no label. It reads exactly
``results/tables/gse245601_malignant_cell_labels.tsv`` -- already frozen,
already containing each epithelial cell's ``cnv_score``,
``cnv_correlation_to_seed``, ``threshold_group`` (the (sample, Seurat
cluster) or ``whole_sample_pooled`` group actually used for thresholding),
and ``primary_malignancy_label`` -- and recomputes the adaptive per-group
thresholds from those already-frozen per-cell values using the exact
formula in ``scripts/analysis/gse245601_05_infercnv_malignant.R``:

    th_value = clamp(mean(cnv_score) - 2 * sd(cnv_score), 0.01, 0.05)
    th_corr  = clamp(mean(correlation) - 1.5 * sd(correlation), 0.2, 0.4)
    malignant iff cnv_score > th_value AND correlation > th_corr

Recomputing this from already-frozen per-cell values is pure arithmetic
audit, not a rerun: ``verify_reconstruction`` checks that reapplying this
formula to the frozen data reproduces the frozen ``primary_malignancy_label``
for every single cell, which is the central correctness gate of this
module -- if it does not reproduce exactly, something about this
reconstruction (not the frozen label) is wrong, and the module raises
rather than reporting a diagnosis built on a mismatched premise.

Provenance ([A]/[B]/[C] tags per PREANALYSIS.md section 3):
- CNV score, top-5%-per-group seed selection, Kendall-tau-to-seed
  correlation, and the mean-k*SD threshold form are [A+B] (paper
  description + author code, ``ng_2021_and_thresholding`` variant).
- The exact clamp bounds ([0.01, 0.05], [0.2, 0.4]) are part of that same
  [A+B] code variant.
- The choice of the clamped ``ng_2021_and_thresholding`` variant itself
  (over the simpler unclamped ``ng_2021`` variant) is [C] -- documented in
  PREANALYSIS.md section 9 as a frozen decision, not a literal recoverable
  default.
- The small-cluster-pooling rule (<10 cells -> ``whole_sample_pooled``) is
  [C] (PREANALYSIS.md 2026-08-11 amendment 3) -- already baked into the
  frozen ``threshold_group`` column consumed here, not re-derived.
- ``MIN_REFERENCE_CELLS=20`` (sample-level exclusion) is [C] but irrelevant
  here: all 12 selected samples have status "ok".

Point 5's sensitivity sweep is diagnostic only: it evaluates malignant
yield on a fixed grid spanning exactly the pre-declared clamp bounds above
(never a search for a better cutoff, never informed by CopyKAT agreement).
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
from matplotlib.lines import Line2D  # noqa: E402

from src.gse245601_cnv_method_comparison import LABEL_MALIGNANT, LABEL_NONMALIGNANT  # noqa: E402

logger = logging.getLogger(__name__)

REQUIRED_COLUMNS = (
    "cell_id",
    "sample_id",
    "patient",
    "condition",
    "threshold_group",
    "cnv_score",
    "cnv_correlation_to_seed",
    "primary_malignancy_label",
)

PASS_BOTH = "passes_both_malignant"
FAIL_CNV_ONLY = "fails_cnv_score_only"
FAIL_CORR_ONLY = "fails_correlation_only"
FAIL_BOTH = "fails_both"
FAILURE_CATEGORIES = (PASS_BOTH, FAIL_CNV_ONLY, FAIL_CORR_ONLY, FAIL_BOTH)

COLOR_MALIGNANT = "#B2182B"
COLOR_NONMALIGNANT = "#2166AC"
COPYKAT_MARKERS = {"malignant": "o", "non-malignant epithelial": "^", "not_defined": "x"}


def _load_config(config_path: str | Path) -> dict:
    with open(config_path) as f:
        return yaml.safe_load(f)


@dataclass(frozen=True)
class ThresholdDiagnosticsConfig:
    """Resolved, config-driven paths and rule constants. No hardcoded
    values -- every constant here is copied from PREANALYSIS.md section 9
    / amendment 3 via config, not reinvented."""

    selected_samples: tuple[str, ...]
    disagreement_tumors: tuple[str, ...]
    good_control_tumors: tuple[str, ...]
    signal_quality_tumors: tuple[str, ...]
    infercnv_labels_tsv: Path
    copykat_labels_tsv: Path
    seed_top_fraction: float
    min_seed_cells: int
    cnv_score_sd_multiplier: float
    cnv_score_clamp: tuple[float, float]
    correlation_sd_multiplier: float
    correlation_clamp: tuple[float, float]
    cnv_score_sweep_range: tuple[float, float]
    cnv_score_grid_points: int
    correlation_sweep_range: tuple[float, float]
    correlation_grid_points: int
    diagnostics_tsv: Path
    sensitivity_grid_tsv: Path
    local_score_sensitivity_tsv: Path
    figures_dir: Path

    @classmethod
    def from_config(cls, config: dict) -> "ThresholdDiagnosticsConfig":
        cfg = config["gse245601_infercnv_threshold_diagnostics"]
        inputs = cfg["inputs"]
        rule = cfg["rule"]
        grid = cfg["sensitivity_grid"]
        out = cfg["output"]
        return cls(
            selected_samples=tuple(cfg["selected_samples"]),
            disagreement_tumors=tuple(cfg["disagreement_tumors"]),
            good_control_tumors=tuple(cfg["good_control_tumors"]),
            signal_quality_tumors=tuple(cfg["signal_quality_tumors"]),
            infercnv_labels_tsv=Path(inputs["infercnv_labels_tsv"]),
            copykat_labels_tsv=Path(inputs["copykat_labels_tsv"]),
            seed_top_fraction=float(rule["seed_top_fraction"]),
            min_seed_cells=int(rule["min_seed_cells"]),
            cnv_score_sd_multiplier=float(rule["cnv_score_sd_multiplier"]),
            cnv_score_clamp=tuple(rule["cnv_score_clamp"]),
            correlation_sd_multiplier=float(rule["correlation_sd_multiplier"]),
            correlation_clamp=tuple(rule["correlation_clamp"]),
            cnv_score_sweep_range=tuple(grid["cnv_score_sweep_range"]),
            cnv_score_grid_points=int(grid["cnv_score_grid_points"]),
            correlation_sweep_range=tuple(grid["correlation_sweep_range"]),
            correlation_grid_points=int(grid["correlation_grid_points"]),
            diagnostics_tsv=Path(out["diagnostics_tsv"]),
            sensitivity_grid_tsv=Path(out["sensitivity_grid_tsv"]),
            local_score_sensitivity_tsv=Path(out["local_score_sensitivity_tsv"]),
            figures_dir=Path(out["figures_dir"]),
        )


def load_cnv_score_table(path: str | Path) -> pd.DataFrame:
    """Load the frozen per-cell InferCNV malignancy table. Read-only --
    this is the same file gse245601_05_infercnv_malignant.R wrote; nothing
    is recomputed by rerunning InferCNV here."""
    df = pd.read_csv(path, sep="\t")
    missing = [c for c in REQUIRED_COLUMNS if c not in df.columns]
    if missing:
        raise ValueError(f"cnv score table missing required columns: {missing}")
    bad_labels = set(df["primary_malignancy_label"].unique()) - {LABEL_MALIGNANT, LABEL_NONMALIGNANT}
    if bad_labels:
        raise ValueError(f"unexpected primary_malignancy_label values: {bad_labels}")
    if df["cnv_score"].isna().any() or df["cnv_correlation_to_seed"].isna().any():
        # R's threshold formula does not pass na.rm=TRUE for cnv_score (only for correlation);
        # this module's pandas-based recomputation skips NaN for both, which only matches R's
        # behavior when neither column actually contains a NaN -- verified here, not assumed.
        raise ValueError(
            "cnv score table contains NaN cnv_score or cnv_correlation_to_seed values -- this module's "
            "threshold recomputation is only verified equivalent to the frozen R formula when neither "
            "column has missing values"
        )
    logger.info("load_cnv_score_table: read %d epithelial cells across %d samples", len(df), df["sample_id"].nunique())
    return df


def _clamp(x: float, lo: float, hi: float) -> float:
    return min(max(x, lo), hi)


def recompute_group_thresholds(
    scores: pd.Series, correlations: pd.Series, cfg: ThresholdDiagnosticsConfig
) -> tuple[float, float]:
    """Exact reimplementation of the frozen formula:
    th_value = clamp(mean(score) - 2*SD(score), 0.01, 0.05);
    th_corr = clamp(mean(corr, na.rm=TRUE) - 1.5*SD(corr, na.rm=TRUE), 0.2,
    0.4). ``pandas`` uses the same ddof=1 (sample) standard deviation as
    R's ``sd()``. Note an asymmetry versus the R source: R's correlation
    threshold explicitly passes ``na.rm=TRUE``, but its CNV-score
    threshold does not (plain ``mean(grp_scores)``/``sd(grp_scores)``) --
    pandas ``.mean()``/``.std()`` skip NaN by default for *both* here,
    which would only diverge from R if ``cnv_score`` itself contained a
    NaN. ``load_cnv_score_table`` asserts neither ``cnv_score`` nor
    ``cnv_correlation_to_seed`` contains any NaN in the loaded table, so
    this asymmetry is verified inert rather than merely assumed inert."""
    th_value_raw = scores.mean() - cfg.cnv_score_sd_multiplier * scores.std()
    th_value = _clamp(th_value_raw, *cfg.cnv_score_clamp)
    th_corr_raw = correlations.mean() - cfg.correlation_sd_multiplier * correlations.std()
    th_corr = _clamp(th_corr_raw, *cfg.correlation_clamp)
    return th_value, th_corr


def identify_seed_cells(group_df: pd.DataFrame, cfg: ThresholdDiagnosticsConfig) -> list[str]:
    """Exact reimplementation of seed selection: top
    ``max(round(n*seed_top_fraction), min_seed_cells)`` cells by cnv_score
    within the group (matching R's ``max(round(n*0.05), 2)``; Python's and
    R's ``round()`` both round half-to-even, and both round a cell COUNT
    here, so the two never disagree on the seed size). This assumes no
    exact tie in cnv_score at the n_seed/n_seed+1 selection boundary --
    reasonable for continuous floating-point scores and not otherwise
    checked here, since this module's threshold values come directly from
    the frozen ``cnv_score``/``cnv_correlation_to_seed`` columns, not from
    reconstructed seed identities; ``identify_seed_cells`` is used only
    for the seed-margin diagnostic (how far the seed itself sits above the
    threshold), not for recomputing thresholds or labels."""
    n = len(group_df)
    n_seed = max(round(n * cfg.seed_top_fraction), cfg.min_seed_cells)
    n_seed = min(n_seed, n)
    top = group_df.sort_values("cnv_score", ascending=False).head(n_seed)
    return top["cell_id"].tolist()


def classify_failure_category(score: float, corr: float, th_value: float, th_corr: float) -> str:
    """Four mutually exclusive categories matching the frozen rule's
    strict ``>`` comparisons exactly."""
    if pd.isna(score) or pd.isna(corr):
        raise ValueError(f"classify_failure_category: NaN input (score={score!r}, corr={corr!r}) not expected in this dataset")
    passes_cnv = score > th_value
    passes_corr = corr > th_corr
    if passes_cnv and passes_corr:
        return PASS_BOTH
    if passes_cnv and not passes_corr:
        return FAIL_CORR_ONLY
    if not passes_cnv and passes_corr:
        return FAIL_CNV_ONLY
    return FAIL_BOTH


def verify_reconstruction(cell_df: pd.DataFrame, cfg: ThresholdDiagnosticsConfig) -> None:
    """Recomputes th_value/th_corr per (sample_id, threshold_group) from
    the frozen per-cell values and checks that the recomputed malignant
    call reproduces ``primary_malignancy_label`` for EVERY cell in the 12
    selected samples, exactly. Raises on any mismatch -- this is the gate
    that proves the reconstruction is faithful before anything downstream
    is computed from it."""
    sub = cell_df.loc[cell_df["sample_id"].isin(cfg.selected_samples)]
    mismatches = []
    for (sample_id, group), grp_df in sub.groupby(["sample_id", "threshold_group"]):
        th_value, th_corr = recompute_group_thresholds(grp_df["cnv_score"], grp_df["cnv_correlation_to_seed"], cfg)
        predicted = (grp_df["cnv_score"] > th_value) & (grp_df["cnv_correlation_to_seed"] > th_corr)
        predicted_label = np.where(predicted, LABEL_MALIGNANT, LABEL_NONMALIGNANT)
        wrong = grp_df.loc[predicted_label != grp_df["primary_malignancy_label"].to_numpy()]
        if len(wrong) > 0:
            mismatches.append((sample_id, group, len(wrong)))
    if mismatches:
        raise ValueError(f"verify_reconstruction: recomputed thresholds do not reproduce frozen labels exactly: {mismatches}")
    logger.info("verify_reconstruction: recomputed thresholds reproduce frozen labels exactly for all %d selected cells", len(sub))


def build_group_diagnostics_table(cell_df: pd.DataFrame, cfg: ThresholdDiagnosticsConfig) -> pd.DataFrame:
    """One row per (sample_id, threshold_group): recomputed thresholds,
    score/correlation distribution, the four-way failure-category
    breakdown, seed-cell diagnostics, and gap-to-cutoff statistics for
    failing cells (Point 4's "how far from the threshold" question)."""
    sub = cell_df.loc[cell_df["sample_id"].isin(cfg.selected_samples)]
    rows = []
    for (sample_id, group), grp_df in sub.groupby(["sample_id", "threshold_group"]):
        th_value, th_corr = recompute_group_thresholds(grp_df["cnv_score"], grp_df["cnv_correlation_to_seed"], cfg)
        categories = [
            classify_failure_category(s, c, th_value, th_corr)
            for s, c in zip(grp_df["cnv_score"], grp_df["cnv_correlation_to_seed"])
        ]
        cat_counts = pd.Series(categories).value_counts()

        seed_ids = identify_seed_cells(grp_df, cfg)
        seed_scores = grp_df.loc[grp_df["cell_id"].isin(seed_ids), "cnv_score"]

        fails_cnv = grp_df.loc[[c in (FAIL_CNV_ONLY, FAIL_BOTH) for c in categories], "cnv_score"]
        fails_corr = grp_df.loc[[c in (FAIL_CORR_ONLY, FAIL_BOTH) for c in categories], "cnv_correlation_to_seed"]
        cnv_gap = (th_value - fails_cnv) if len(fails_cnv) > 0 else pd.Series(dtype=float)
        corr_gap = (th_corr - fails_corr) if len(fails_corr) > 0 else pd.Series(dtype=float)

        patient, condition = grp_df["patient"].iloc[0], grp_df["condition"].iloc[0]
        rows.append(
            {
                "sample_id": sample_id,
                "patient": patient,
                "condition": condition,
                "threshold_group": group,
                "n_cells": len(grp_df),
                "cnv_score_threshold": th_value,
                "correlation_threshold": th_corr,
                "median_cnv_score": grp_df["cnv_score"].median(),
                "min_cnv_score": grp_df["cnv_score"].min(),
                "max_cnv_score": grp_df["cnv_score"].max(),
                "median_correlation": grp_df["cnv_correlation_to_seed"].median(),
                "min_correlation": grp_df["cnv_correlation_to_seed"].min(),
                "max_correlation": grp_df["cnv_correlation_to_seed"].max(),
                "malignant_count": int(cat_counts.get(PASS_BOTH, 0)),
                "fails_cnv_score_only_count": int(cat_counts.get(FAIL_CNV_ONLY, 0)),
                "fails_correlation_only_count": int(cat_counts.get(FAIL_CORR_ONLY, 0)),
                "fails_both_count": int(cat_counts.get(FAIL_BOTH, 0)),
                "passes_both_count": int(cat_counts.get(PASS_BOTH, 0)),
                "n_seed_cells": len(seed_ids),
                "seed_mean_cnv_score": seed_scores.mean(),
                "seed_margin_above_cnv_threshold": seed_scores.mean() - th_value,
                "median_cnv_score_gap_below_threshold": cnv_gap.median() if len(cnv_gap) > 0 else float("nan"),
                "median_correlation_gap_below_threshold": corr_gap.median() if len(corr_gap) > 0 else float("nan"),
            }
        )

    out = pd.DataFrame(rows).sort_values(["sample_id", "threshold_group"]).reset_index(drop=True)
    total = out[["malignant_count", "fails_cnv_score_only_count", "fails_correlation_only_count", "fails_both_count"]].sum(axis=1)
    if not (total == out["n_cells"]).all():
        raise ValueError("build_group_diagnostics_table: failure-category counts do not sum to n_cells for every group")
    logger.info("build_group_diagnostics_table: %d (sample, threshold_group) rows built", len(out))
    return out


def build_sensitivity_grid(cell_df: pd.DataFrame, cfg: ThresholdDiagnosticsConfig) -> pd.DataFrame:
    """Diagnostic-only sweep (Point 5): for each selected sample, pools
    ALL of that sample's epithelial cells (across threshold_group) and
    evaluates malignant yield on a fixed grid (correlation: the rule's own
    clamp bounds; CNV score: a slightly wider range down to the metric's
    natural 0.0 floor, see config comment) -- never a search for a better
    cutoff, never informed by CopyKAT. This intentionally ignores the real
    rule's per-group adaptivity (it asks a different, simpler question:
    how sensitive is the overall malignant count to the threshold pair, in
    the abstract) -- the per-group actual thresholds are reported
    separately in the diagnostics table and overlaid on the figures."""
    value_grid = np.linspace(cfg.cnv_score_sweep_range[0], cfg.cnv_score_sweep_range[1], cfg.cnv_score_grid_points)
    corr_grid = np.linspace(cfg.correlation_sweep_range[0], cfg.correlation_sweep_range[1], cfg.correlation_grid_points)

    rows = []
    for sample_id in cfg.selected_samples:
        grp_df = cell_df.loc[cell_df["sample_id"] == sample_id]
        n_cells = len(grp_df)
        scores = grp_df["cnv_score"].to_numpy()
        corrs = grp_df["cnv_correlation_to_seed"].to_numpy()
        for tv in value_grid:
            for tc in corr_grid:
                n_malignant = int(((scores > tv) & (corrs > tc)).sum())
                rows.append(
                    {
                        "sample_id": sample_id,
                        "cnv_score_threshold": tv,
                        "correlation_threshold": tc,
                        "n_cells": n_cells,
                        "malignant_count": n_malignant,
                        "malignant_fraction": n_malignant / n_cells if n_cells > 0 else float("nan"),
                    }
                )
    out = pd.DataFrame(rows)
    logger.info(
        "build_sensitivity_grid: %d samples x %d x %d grid = %d rows",
        len(cfg.selected_samples),
        cfg.cnv_score_grid_points,
        cfg.correlation_grid_points,
        len(out),
    )
    return out


def build_local_score_sensitivity(
    cell_df: pd.DataFrame, cfg: ThresholdDiagnosticsConfig, band_half_width: float = 0.004, n_points: int = 17
) -> pd.DataFrame:
    """Dense, CNV-score-only sweep in a narrow band directly around the
    rule's actual clamp floor (0.01) -- added specifically to test whether
    the low malignant yield in Tumor_01/04/08 reflects a knife-edge right
    at the operating threshold, or a broader shortfall (the coarse
    ``build_sensitivity_grid`` sweep, at 0.005 spacing, cannot resolve
    this). Correlation is deliberately ignored here (fraction is computed
    from ``cnv_score`` alone, pooled across all of a sample's cells) since
    ``build_group_diagnostics_table`` already establishes correlation is
    not the binding constraint for these samples -- isolating the
    CNV-score axis alone gives the cleanest read on its local steepness."""
    grid = np.linspace(0.01 - band_half_width, 0.01 + band_half_width, n_points)
    rows = []
    for sample_id in cfg.selected_samples:
        scores = cell_df.loc[cell_df["sample_id"] == sample_id, "cnv_score"].to_numpy()
        n_cells = len(scores)
        for tv in grid:
            n_pass = int((scores > tv).sum())
            rows.append(
                {
                    "sample_id": sample_id,
                    "cnv_score_threshold": tv,
                    "n_cells": n_cells,
                    "n_pass_score_only": n_pass,
                    "fraction_pass_score_only": n_pass / n_cells if n_cells > 0 else float("nan"),
                }
            )
    out = pd.DataFrame(rows)
    logger.info("build_local_score_sensitivity: %d samples x %d grid points around the 0.01 floor", len(cfg.selected_samples), n_points)
    return out


def build_local_sensitivity_figure(local_df: pd.DataFrame, cfg: ThresholdDiagnosticsConfig, out_path: str | Path) -> None:
    """One line per sample: fraction of that sample's cells passing the
    CNV-score criterion alone, as a dense function of the threshold value
    right around the actual 0.01 floor -- steep near 0.01 means the
    result is sensitive to exactly where the floor is set; flat near 0.01
    means the operating point is not special (the population is already
    far from the boundary in that region)."""
    fig, ax = plt.subplots(figsize=(7, 5))
    colors = plt.get_cmap("tab20").colors
    for i, sample_id in enumerate(cfg.selected_samples):
        sub = local_df.loc[local_df["sample_id"] == sample_id].sort_values("cnv_score_threshold")
        ax.plot(sub["cnv_score_threshold"], sub["fraction_pass_score_only"], color=colors[i % len(colors)], label=sample_id, linewidth=1.5)
    ax.axvline(0.01, color="black", linewidth=1, linestyle="--", label="actual floor (0.01)")
    ax.set_xlabel("CNV score threshold (dense local sweep)")
    ax.set_ylabel("fraction of cells passing CNV-score criterion alone")
    ax.set_title("Local sensitivity of CNV-score criterion around the actual 0.01 floor")
    ax.legend(fontsize=6, ncol=2, loc="upper right")
    fig.tight_layout()
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info("build_local_sensitivity_figure: wrote %s", out_path)


def build_scatter_figure(
    sample_id: str, cell_df: pd.DataFrame, group_table: pd.DataFrame, copykat_df: pd.DataFrame | None, out_path: str | Path
) -> None:
    """One figure per sample: one panel per threshold_group actually
    present in that sample, x=cnv_score, y=correlation, colored by the
    frozen InferCNV call, CopyKAT status shown as marker shape (visual
    overlay only -- never used to draw or move the threshold lines), and
    the two ACTUAL recomputed thresholds for that group drawn as lines."""
    sample_cells = cell_df.loc[cell_df["sample_id"] == sample_id]
    sample_groups = group_table.loc[group_table["sample_id"] == sample_id].sort_values("n_cells", ascending=False)

    if copykat_df is not None:
        ck = copykat_df.loc[copykat_df["sample_id"] == sample_id, ["cell_id", "sensitivity_malignancy_label"]]
        sample_cells = sample_cells.merge(ck, on="cell_id", how="left")
    else:
        sample_cells = sample_cells.copy()
        sample_cells["sensitivity_malignancy_label"] = "not_defined"

    n_groups = len(sample_groups)
    n_cols = min(4, n_groups)
    n_rows = int(np.ceil(n_groups / n_cols))
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(4.2 * n_cols, 3.8 * n_rows), squeeze=False)

    for i, group_row in enumerate(sample_groups.itertuples(index=False)):
        ax = axes[i // n_cols, i % n_cols]
        grp_cells = sample_cells.loc[sample_cells["threshold_group"] == group_row.threshold_group]
        for infercnv_label, color in ((LABEL_MALIGNANT, COLOR_MALIGNANT), (LABEL_NONMALIGNANT, COLOR_NONMALIGNANT)):
            for ck_label, marker in COPYKAT_MARKERS.items():
                pts = grp_cells.loc[
                    (grp_cells["primary_malignancy_label"] == infercnv_label)
                    & (grp_cells["sensitivity_malignancy_label"] == ck_label)
                ]
                if len(pts) == 0:
                    continue
                ax.scatter(
                    pts["cnv_score"], pts["cnv_correlation_to_seed"], c=color, marker=marker, s=14, alpha=0.6, linewidths=0
                )
        ax.axvline(group_row.cnv_score_threshold, color="black", linewidth=0.8, linestyle="--")
        ax.axhline(group_row.correlation_threshold, color="black", linewidth=0.8, linestyle="--")
        ax.set_title(f"group={group_row.threshold_group} (n={group_row.n_cells})", fontsize=8)
        ax.set_xlabel("CNV score", fontsize=7)
        ax.set_ylabel("Kendall correlation to seed", fontsize=7)
        ax.tick_params(labelsize=6)

    for j in range(n_groups, n_rows * n_cols):
        axes[j // n_cols, j % n_cols].axis("off")

    label_handles = [
        mpatches.Patch(color=COLOR_MALIGNANT, label="InferCNV malignant"),
        mpatches.Patch(color=COLOR_NONMALIGNANT, label="InferCNV non-malignant"),
    ]
    marker_handles = [
        Line2D([0], [0], marker=m, color="gray", linestyle="", label=f"CopyKAT {lbl}", markersize=6)
        for lbl, m in COPYKAT_MARKERS.items()
    ]
    fig.legend(handles=label_handles + marker_handles, loc="upper center", ncol=5, fontsize=7, frameon=False, bbox_to_anchor=(0.5, 1.02))
    fig.suptitle(f"{sample_id} -- CNV score vs. correlation-to-seed, by threshold group", fontsize=10, y=1.08)

    fig.tight_layout()
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info("build_scatter_figure: wrote %s", out_path)


def build_sensitivity_figure(sample_id: str, grid_df: pd.DataFrame, group_table: pd.DataFrame, out_path: str | Path) -> None:
    """Heatmap of malignant fraction over the (cnv_score_threshold,
    correlation_threshold) grid for one sample, with the sample's actual
    per-group recomputed thresholds overlaid as markers."""
    sample_grid = grid_df.loc[grid_df["sample_id"] == sample_id]
    value_grid = sorted(sample_grid["cnv_score_threshold"].unique())
    corr_grid = sorted(sample_grid["correlation_threshold"].unique())
    pivot = sample_grid.pivot(index="correlation_threshold", columns="cnv_score_threshold", values="malignant_fraction")
    pivot = pivot.loc[corr_grid, value_grid]

    fig, ax = plt.subplots(figsize=(5.5, 4.5))
    im = ax.imshow(pivot.to_numpy(), aspect="auto", origin="lower", cmap="viridis", vmin=0, vmax=1)
    ax.set_xticks(range(len(value_grid)))
    ax.set_xticklabels([f"{v:.3f}" for v in value_grid], rotation=90, fontsize=7)
    ax.set_yticks(range(len(corr_grid)))
    ax.set_yticklabels([f"{v:.3f}" for v in corr_grid], fontsize=7)
    ax.set_xlabel("CNV score threshold (grid)")
    ax.set_ylabel("Correlation threshold (grid)")
    fig.colorbar(im, ax=ax, label="malignant fraction (pooled, whole sample)")

    sample_groups = group_table.loc[group_table["sample_id"] == sample_id]
    for row in sample_groups.itertuples(index=False):
        x = np.interp(row.cnv_score_threshold, value_grid, range(len(value_grid)))
        y = np.interp(row.correlation_threshold, corr_grid, range(len(corr_grid)))
        ax.scatter([x], [y], color="red", marker="+", s=80, linewidths=1.5)

    ax.set_title(f"{sample_id}\nred + = actual per-group threshold(s) used", fontsize=9)
    fig.tight_layout()
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info("build_sensitivity_figure: wrote %s", out_path)


def run_threshold_diagnostics(
    config_path: str | Path = "config/config.yaml", samples: tuple[str, ...] | None = None
) -> dict[str, object]:
    """Full pipeline: load, verify exact reconstruction, build the group
    diagnostics table, the sensitivity grid, and every figure. ``samples``
    restricts figure generation to a subset (used by tests) -- the
    diagnostics table and grid are still always computed over all of
    ``cfg.selected_samples`` found in ``samples`` (or all of them if
    ``samples`` is None)."""
    config = _load_config(config_path)
    cfg = ThresholdDiagnosticsConfig.from_config(config)

    cell_df = load_cnv_score_table(cfg.infercnv_labels_tsv)
    copykat_df = pd.read_csv(cfg.copykat_labels_tsv, sep="\t") if cfg.copykat_labels_tsv.exists() else None

    verify_reconstruction(cell_df, cfg)

    group_table = build_group_diagnostics_table(cell_df, cfg)
    grid_table = build_sensitivity_grid(cell_df, cfg)
    local_table = build_local_score_sensitivity(cell_df, cfg)

    Path(cfg.diagnostics_tsv).parent.mkdir(parents=True, exist_ok=True)
    group_table.to_csv(cfg.diagnostics_tsv, sep="\t", index=False)
    grid_table.to_csv(cfg.sensitivity_grid_tsv, sep="\t", index=False)
    local_table.to_csv(cfg.local_score_sensitivity_tsv, sep="\t", index=False)
    logger.info(
        "wrote %s (%d rows), %s (%d rows), and %s (%d rows)",
        cfg.diagnostics_tsv,
        len(group_table),
        cfg.sensitivity_grid_tsv,
        len(grid_table),
        cfg.local_score_sensitivity_tsv,
        len(local_table),
    )

    samples_to_plot = samples if samples is not None else cfg.selected_samples
    figure_paths = {}
    for sample_id in samples_to_plot:
        scatter_path = cfg.figures_dir / f"{sample_id}_score_vs_correlation.png"
        build_scatter_figure(sample_id, cell_df, group_table, copykat_df, scatter_path)
        sensitivity_path = cfg.figures_dir / f"{sample_id}_sensitivity_grid.png"
        build_sensitivity_figure(sample_id, grid_table, group_table, sensitivity_path)
        figure_paths[sample_id] = {"scatter": scatter_path, "sensitivity": sensitivity_path}

    local_sensitivity_fig_path = cfg.figures_dir / "local_score_sensitivity_all_samples.png"
    build_local_sensitivity_figure(local_table, cfg, local_sensitivity_fig_path)

    return {
        "group_table": group_table,
        "grid_table": grid_table,
        "local_table": local_table,
        "figure_paths": figure_paths,
        "local_sensitivity_figure": local_sensitivity_fig_path,
    }


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run_threshold_diagnostics()
