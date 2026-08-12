"""USP34 expression visualization on the existing GSE245601 UMAP
(read-only follow-up to docs/gse245601_PREANALYSIS.md section 13).

Descriptive only -- no differential-expression test, no candidate
ranking, no new statistic beyond simple per-group summaries. Reuses the
already-frozen per-cell table from the pseudobulk phase
(``gse245601_pseudobulk_..._cell_level_summary.tsv``: condition,
malignancy_status, and UMAP coordinates for every epithelial cell,
already computed during the original candidate-blind clustering) and
joins in USP34's own log-normalized expression, extracted read-only by
``scripts/analysis/gse245601_17_extract_usp34_expression.R``. Clustering,
annotation, UMAP, malignancy labels, and pseudobulk eligibility rules are
none of them recomputed or modified here.
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

MALIGNANT = "malignant"
NONMALIGNANT = "non-malignant epithelial"
CMAP = "viridis"


def _load_config(config_path: str | Path) -> dict:
    with open(config_path) as f:
        return yaml.safe_load(f)


def load_usp34_data(cell_level_summary_tsv: str | Path, expression_tsv: str | Path) -> pd.DataFrame:
    """Join the frozen per-cell metadata/UMAP table with USP34's
    per-cell expression, by cell_id. Read-only; raises if the two tables
    do not cover exactly the same cell set (a mismatch would mean the
    expression extraction ran against a different cell population)."""
    meta = pd.read_csv(cell_level_summary_tsv, sep="\t")
    expr = pd.read_csv(expression_tsv, sep="\t")
    if set(meta["cell_id"]) != set(expr["cell_id"]):
        raise ValueError("cell_id sets differ between the frozen cell-level summary and the USP34 expression table")
    out = meta.merge(expr, on="cell_id", how="inner", validate="one_to_one")
    if len(out) != len(meta):
        raise ValueError("join between cell-level summary and USP34 expression table lost or duplicated rows")
    bad = set(out["malignancy_status"].unique()) - {MALIGNANT, NONMALIGNANT}
    if bad:
        raise ValueError(f"unexpected malignancy_status values: {bad}")
    logger.info("load_usp34_data: %d epithelial cells (%d expressing USP34 > 0)", len(out), (out["usp34_log_norm_expression"] > 0).sum())
    return out


def build_group_summary(df: pd.DataFrame) -> pd.DataFrame:
    """One row per requested group: n_cells, pct_expressing (>0),
    mean_expression, median_expression -- log-normalized expression
    throughout, not raw counts."""
    groups = [
        ("all_epithelial", "Control", df["condition"] == "Control"),
        ("all_epithelial", "Tamoxifen", df["condition"] == "Tamoxifen"),
        ("malignant", "Control", (df["malignancy_status"] == MALIGNANT) & (df["condition"] == "Control")),
        ("malignant", "Tamoxifen", (df["malignancy_status"] == MALIGNANT) & (df["condition"] == "Tamoxifen")),
        ("nonmalignant", "Control", (df["malignancy_status"] == NONMALIGNANT) & (df["condition"] == "Control")),
        ("nonmalignant", "Tamoxifen", (df["malignancy_status"] == NONMALIGNANT) & (df["condition"] == "Tamoxifen")),
    ]
    rows = []
    for population, condition, mask in groups:
        sub = df.loc[mask, "usp34_log_norm_expression"]
        rows.append(
            {
                "population": population,
                "condition": condition,
                "n_cells": int(mask.sum()),
                "pct_expressing": 100.0 * (sub > 0).mean() if len(sub) > 0 else float("nan"),
                "mean_expression": sub.mean() if len(sub) > 0 else float("nan"),
                "median_expression": sub.median() if len(sub) > 0 else float("nan"),
            }
        )
    out = pd.DataFrame(rows)
    logger.info("build_group_summary: %d groups", len(out))
    return out


def _scatter_panel(ax, df: pd.DataFrame, vmax: float, title: str):
    order = df["usp34_log_norm_expression"].sort_values().index  # draw higher-expressing cells last (on top)
    sc = ax.scatter(
        df.loc[order, "umap_1"], df.loc[order, "umap_2"], c=df.loc[order, "usp34_log_norm_expression"],
        cmap=CMAP, vmin=0, vmax=vmax, s=3, linewidths=0,
    )
    ax.set_xlabel("UMAP_1")
    ax.set_ylabel("UMAP_2")
    ax.set_title(f"{title} (n={len(df)})", fontsize=9)
    return sc


def plot_all_epithelial(df: pd.DataFrame, vmax: float, out_path: str | Path) -> None:
    fig, ax = plt.subplots(figsize=(6, 5.5))
    sc = _scatter_panel(ax, df, vmax, "USP34 -- all epithelial cells")
    fig.colorbar(sc, ax=ax, label="log-normalized expression")
    fig.tight_layout()
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info("wrote %s", out_path)


def plot_split_panels(df: pd.DataFrame, split_col: str, split_values: list[str], vmax: float, suptitle: str, out_path: str | Path) -> None:
    fig, axes = plt.subplots(1, len(split_values), figsize=(5.5 * len(split_values), 5), sharex=True, sharey=True)
    if len(split_values) == 1:
        axes = [axes]
    sc = None
    for ax, value in zip(axes, split_values):
        sub = df.loc[df[split_col] == value]
        sc = _scatter_panel(ax, sub, vmax, value)
    fig.colorbar(sc, ax=axes, label="log-normalized expression", fraction=0.025, pad=0.02)
    fig.suptitle(suptitle, fontsize=11)
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info("wrote %s", out_path)


def run_usp34_expression(config_path: str | Path = "config/config.yaml") -> dict[str, pd.DataFrame]:
    config = _load_config(config_path)
    cfg = config["gse245601_usp34_expression"]
    figures_dir = Path(cfg["output"]["figures_dir"])

    df = load_usp34_data(cfg["inputs"]["cell_level_summary_tsv"], cfg["output"]["expression_tsv"])
    summary = build_group_summary(df)

    summary_path = Path(cfg["output"]["summary_tsv"])
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(summary_path, sep="\t", index=False)
    logger.info("wrote %s", summary_path)

    vmax = float(np.percentile(df["usp34_log_norm_expression"], 99))
    vmax = max(vmax, 1e-6)

    plot_all_epithelial(df, vmax, figures_dir / "usp34_all_epithelial_umap.png")
    plot_split_panels(df, "condition", ["Control", "Tamoxifen"], vmax, "USP34 by treatment (all epithelial cells)", figures_dir / "usp34_by_treatment_umap.png")
    plot_split_panels(
        df, "malignancy_status", [NONMALIGNANT, MALIGNANT], vmax, "USP34 by frozen malignancy label", figures_dir / "usp34_by_malignancy_umap.png"
    )
    malignant_df = df.loc[df["malignancy_status"] == MALIGNANT]
    plot_split_panels(
        malignant_df, "condition", ["Control", "Tamoxifen"], vmax, "USP34 by treatment -- malignant cells only", figures_dir / "usp34_malignant_by_treatment_umap.png"
    )

    return {"cell_level": df, "summary": summary}


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run_usp34_expression()
