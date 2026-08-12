"""GSE240112 Phase 12 USP34-focused visualization: UMAP feature plots,
primary-vs-recurrent split, and per-sample distribution/detection, on the
tumor-cell compartment (the only PT/RT-covering single-cell object the
authors publicly released -- see docs/GSE240112_DATA_AUDIT.md section 4;
there is no separate "all cells" object to plot alongside it, so the
tumor-cell UMAP serves both the Phase 12 "all cells" and "epithelial/
tumor compartment" panels, labeled accordingly).

Descriptive only -- no per-cell significance test is computed or implied
here; the only inferential (p-value/FDR) result for USP34 comes from the
sample-level pseudobulk edgeR test (gse240112_candidate_extraction.py).

Data source: GSE240112 (Fang et al., Genome Medicine 2024, PMID 39558215),
tumor-cell object, version as downloaded 2026-08-12.
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

GROUP_COLORS = {"PT": "#4C72B0", "RT": "#C44E52"}
SAMPLE_ORDER = ["PT1", "PT2", "PT3", "RT1", "RT2", "RT3"]
CMAP = "viridis"


def _load_config(config_path: str | Path) -> dict:
    with open(config_path) as f:
        return yaml.safe_load(f)


def load_usp34_cell_data(metadata_tsv: str | Path, lognorm_tsv: str | Path, gene: str = "USP34") -> pd.DataFrame:
    """Join tumor-cell metadata (sample, group, UMAP coords) with USP34's
    per-cell log-normalized expression. Raises on any cell_id mismatch."""
    meta = pd.read_csv(metadata_tsv, sep="\t")
    expr = pd.read_csv(lognorm_tsv, sep="\t")
    if gene not in expr.columns:
        raise ValueError(f"{gene} not present in the per-cell log-normalized expression table")
    expr = expr[["cell_id", gene]].rename(columns={gene: "log_norm_expression"})
    if set(meta["cell_id"]) != set(expr["cell_id"]):
        raise ValueError("cell_id sets differ between tumor-cell metadata and candidate expression table")
    out = meta.merge(expr, on="cell_id", how="inner", validate="one_to_one")
    if len(out) != len(meta):
        raise ValueError("join lost or duplicated rows")
    logger.info("load_usp34_cell_data: %d tumor cells (%d expressing %s > 0)", len(out), (out["log_norm_expression"] > 0).sum(), gene)
    return out


def build_per_sample_summary(df: pd.DataFrame) -> pd.DataFrame:
    """One row per sample: n_cells, pct_expressing (>0), mean/median
    log-normalized expression."""
    rows = []
    for sample_id in SAMPLE_ORDER:
        sub = df.loc[df["orig.ident"] == sample_id, "log_norm_expression"]
        rows.append(
            {
                "sample_id": sample_id,
                "group": sample_id[:2],
                "n_cells": int(len(sub)),
                "pct_expressing": 100.0 * (sub > 0).mean() if len(sub) > 0 else float("nan"),
                "mean_expression": sub.mean() if len(sub) > 0 else float("nan"),
                "median_expression": sub.median() if len(sub) > 0 else float("nan"),
            }
        )
    out = pd.DataFrame(rows)
    logger.info("build_per_sample_summary: %d samples", len(out))
    return out


def plot_umap_all(df: pd.DataFrame, vmax: float, out_path: str | Path) -> None:
    fig, ax = plt.subplots(figsize=(6, 5.5))
    order = df["log_norm_expression"].sort_values().index
    sc = ax.scatter(df.loc[order, "umap_1"], df.loc[order, "umap_2"], c=df.loc[order, "log_norm_expression"], cmap=CMAP, vmin=0, vmax=vmax, s=4, linewidths=0)
    fig.colorbar(sc, ax=ax, label="log-normalized expression")
    ax.set_xlabel("UMAP_1")
    ax.set_ylabel("UMAP_2")
    ax.set_title(f"USP34 -- tumor-cell compartment (n={len(df)})", fontsize=10)
    fig.tight_layout()
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info("wrote %s", out_path)


def plot_umap_by_group(df: pd.DataFrame, vmax: float, out_path: str | Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(11, 5), sharex=True, sharey=True)
    sc = None
    for ax, group in zip(axes, ["PT", "RT"]):
        sub = df.loc[df["group"] == group]
        order = sub["log_norm_expression"].sort_values().index
        sc = ax.scatter(sub.loc[order, "umap_1"], sub.loc[order, "umap_2"], c=sub.loc[order, "log_norm_expression"], cmap=CMAP, vmin=0, vmax=vmax, s=4, linewidths=0)
        ax.set_title(f"{group} (n={len(sub)})", fontsize=10)
        ax.set_xlabel("UMAP_1")
    axes[0].set_ylabel("UMAP_2")
    fig.colorbar(sc, ax=axes, label="log-normalized expression", fraction=0.025, pad=0.02)
    fig.suptitle("USP34 -- tumor cells, primary (PT) vs recurrent (RT)", fontsize=11)
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info("wrote %s", out_path)


def plot_per_sample_distribution(df: pd.DataFrame, out_path: str | Path) -> None:
    fig, ax = plt.subplots(figsize=(7, 5))
    data = [df.loc[df["orig.ident"] == s, "log_norm_expression"].to_numpy() for s in SAMPLE_ORDER]
    colors = [GROUP_COLORS[s[:2]] for s in SAMPLE_ORDER]
    parts = ax.violinplot(data, showmeans=True, showextrema=False)
    for pc, color in zip(parts["bodies"], colors):
        pc.set_facecolor(color)
        pc.set_alpha(0.7)
    ax.set_xticks(range(1, len(SAMPLE_ORDER) + 1))
    ax.set_xticklabels(SAMPLE_ORDER)
    ax.set_ylabel("USP34 log-normalized expression")
    ax.set_title("USP34 tumor-cell expression distribution by sample", fontsize=10)
    from matplotlib.patches import Patch

    handles = [Patch(color=c, label=g) for g, c in GROUP_COLORS.items()]
    ax.legend(handles=handles, fontsize=8)
    fig.tight_layout()
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info("wrote %s", out_path)


def plot_per_sample_detection(summary: pd.DataFrame, out_path: str | Path) -> None:
    fig, ax = plt.subplots(figsize=(7, 4.5))
    colors = [GROUP_COLORS[g] for g in summary["group"]]
    ax.bar(summary["sample_id"], summary["pct_expressing"], color=colors)
    ax.set_ylabel("% tumor cells expressing USP34")
    ax.set_title("USP34 detection rate by sample", fontsize=10)
    from matplotlib.patches import Patch

    handles = [Patch(color=c, label=g) for g, c in GROUP_COLORS.items()]
    ax.legend(handles=handles, fontsize=8)
    fig.tight_layout()
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info("wrote %s", out_path)


def run_usp34_visualization(config_path: str | Path = "config/config.yaml") -> dict[str, pd.DataFrame]:
    config = _load_config(config_path)
    cfg = config["gse240112"]
    figures_dir = Path(cfg["output"]["usp34"]["figures_dir"])

    df = load_usp34_cell_data(cfg["output"]["tt_cancer_metadata_tsv"], cfg["output"]["tt_cancer_candidate_lognorm_tsv"])
    summary = build_per_sample_summary(df)

    summary_path = Path(cfg["output"]["usp34"]["per_sample_tsv"])
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(summary_path, sep="\t", index=False)
    logger.info("wrote %s", summary_path)

    vmax = max(float(np.percentile(df["log_norm_expression"], 99)), 1e-6)

    plot_umap_all(df, vmax, figures_dir / "usp34_tumor_cell_umap.png")
    plot_umap_by_group(df, vmax, figures_dir / "usp34_umap_by_group.png")
    plot_per_sample_distribution(df, figures_dir / "usp34_per_sample_distribution.png")
    plot_per_sample_detection(summary, figures_dir / "usp34_per_sample_detection.png")

    return {"cell_level": df, "summary": summary}


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run_usp34_visualization()
