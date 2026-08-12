"""GSE245601 candidate deep-dive Phase 18: characterizes USP34-high vs
USP34-low cells WITHIN the 12h-Tamoxifen epithelial cells only, using a
predeclared, data-driven, non-manufactured threshold (top vs bottom
quartile of log-normalized USP34 among USP34-detectable Tamoxifen
epithelial cells). Purely descriptive -- no claim of future clonal
resistance, no temporal claim.
"""

from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd
import yaml

logger = logging.getLogger(__name__)


def _load_config(config_path: str | Path) -> dict:
    with open(config_path) as f:
        return yaml.safe_load(f)


def build_usp34_high_distribution(per_cell: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    tam = per_cell.loc[per_cell["condition"] == "Tamoxifen"].copy()
    detectable = tam.loc[tam["USP34_raw_count"] > 0].copy()
    q75 = detectable["USP34_log_norm"].quantile(0.75)
    q25 = detectable["USP34_log_norm"].quantile(0.25)
    detectable["usp34_group"] = pd.cut(detectable["USP34_log_norm"], bins=[-1, q25, q75, detectable["USP34_log_norm"].max() + 1], labels=["low_quartile", "mid", "high_quartile"], include_lowest=True)

    high = detectable.loc[detectable["usp34_group"] == "high_quartile"]

    by_cluster = high.groupby("seurat_clusters", observed=True).size().rename("n_usp34_high_cells").reset_index()
    by_cluster["pct_of_all_usp34_high"] = 100 * by_cluster["n_usp34_high_cells"] / len(high)
    cluster_totals = tam.groupby("seurat_clusters", observed=True).size().rename("n_total_tam_cells_in_cluster")
    by_cluster = by_cluster.merge(cluster_totals, on="seurat_clusters", how="left")
    by_cluster["pct_of_cluster_that_is_usp34_high"] = 100 * by_cluster["n_usp34_high_cells"] / by_cluster["n_total_tam_cells_in_cluster"]

    by_patient = high.groupby("patient", observed=True).size().rename("n_usp34_high_cells").reset_index()
    by_patient["pct_of_all_usp34_high"] = 100 * by_patient["n_usp34_high_cells"] / len(high)
    patient_totals = tam.groupby("patient", observed=True).size().rename("n_total_tam_cells")
    by_patient = by_patient.merge(patient_totals, on="patient", how="left")
    by_patient["pct_of_patient_tam_cells_that_are_usp34_high"] = 100 * by_patient["n_usp34_high_cells"] / by_patient["n_total_tam_cells"]

    malignancy_breakdown = high["malignancy_status"].value_counts(normalize=True).mul(100).rename("pct").reset_index().rename(columns={"index": "malignancy_status"})
    tam_malignancy_breakdown = tam["malignancy_status"].value_counts(normalize=True).mul(100)
    n_malignant_tam = int((tam["malignancy_status"] == "malignant").sum())
    n_malignant_high = int((high["malignancy_status"] == "malignant").sum())
    pct_malignant_that_are_high = 100 * n_malignant_high / n_malignant_tam if n_malignant_tam else float("nan")

    summary_rows = [
        {"metric": "n_tam_epithelial_cells", "value": len(tam)},
        {"metric": "n_tam_usp34_detectable_cells", "value": len(detectable)},
        {"metric": "usp34_log_norm_q25_threshold", "value": q25},
        {"metric": "usp34_log_norm_q75_threshold", "value": q75},
        {"metric": "n_usp34_high_cells", "value": len(high)},
        {"metric": "pct_high_cells_malignant", "value": float(malignancy_breakdown.set_index("malignancy_status").loc["malignant", "pct"]) if "malignant" in malignancy_breakdown["malignancy_status"].values else 0.0},
        {"metric": "pct_high_cells_nonmalignant", "value": float(malignancy_breakdown.set_index("malignancy_status").loc["non-malignant epithelial", "pct"]) if "non-malignant epithelial" in malignancy_breakdown["malignancy_status"].values else 0.0},
        {"metric": "pct_all_tam_cells_malignant_for_reference", "value": float(tam_malignancy_breakdown.get("malignant", 0.0))},
        {"metric": "pct_malignant_tam_cells_that_are_usp34_high", "value": pct_malignant_that_are_high},
        {"metric": "n_patients_contributing_usp34_high_cells", "value": high["patient"].nunique()},
        {"metric": "max_single_patient_pct_of_usp34_high_cells", "value": float(by_patient["pct_of_all_usp34_high"].max()) if len(by_patient) else float("nan")},
        {"metric": "max_contributing_patient", "value": by_patient.sort_values("pct_of_all_usp34_high", ascending=False)["patient"].iloc[0] if len(by_patient) else "NA"},
    ]

    detail = pd.concat([by_cluster.assign(breakdown_by="cluster").rename(columns={"seurat_clusters": "group"}), by_patient.assign(breakdown_by="patient").rename(columns={"patient": "group"})], ignore_index=True)
    summary = pd.DataFrame(summary_rows).set_index("metric")["value"]
    logger.info(
        "build_usp34_high_distribution: %d high cells, %.1f%% malignant (vs %.1f%% baseline), top patient=%s (%.1f%% of high cells)",
        len(high), summary["pct_high_cells_malignant"], summary["pct_all_tam_cells_malignant_for_reference"], summary["max_contributing_patient"], summary["max_single_patient_pct_of_usp34_high_cells"],
    )
    return detail, summary.reset_index().rename(columns={"index": "metric"})


def plot_usp34_high_cells(per_cell: pd.DataFrame, detail: pd.DataFrame, out_path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 3, figsize=(14, 4.5))

    tam = per_cell.loc[per_cell["condition"] == "Tamoxifen"]
    detectable = tam.loc[tam["USP34_raw_count"] > 0]
    ax = axes[0]
    ax.hist(detectable["USP34_log_norm"], bins=40, color="#4C72B0")
    q25, q75 = detectable["USP34_log_norm"].quantile([0.25, 0.75])
    ax.axvline(q25, color="gray", linestyle="--", linewidth=1)
    ax.axvline(q75, color="#C44E52", linestyle="--", linewidth=1.5, label="top-quartile threshold")
    ax.set_xlabel("USP34 log-norm expression")
    ax.set_ylabel("n cells (Tam, USP34-detectable)")
    ax.set_title("A. Threshold definition\n(top vs bottom quartile)", fontsize=9)
    ax.legend(fontsize=7, frameon=False)

    ax = axes[1]
    by_cluster = detail.loc[detail["breakdown_by"] == "cluster"].sort_values("pct_of_cluster_that_is_usp34_high", ascending=False)
    ax.barh(by_cluster["group"].astype(str), by_cluster["pct_of_cluster_that_is_usp34_high"], color="#55A868")
    ax.set_xlabel("% of cluster's Tam cells that are USP34-high")
    ax.set_ylabel("seurat_clusters")
    ax.set_title("B. Which clusters contain\nUSP34-high Tam cells?", fontsize=9)

    ax = axes[2]
    by_patient = detail.loc[detail["breakdown_by"] == "patient"].sort_values("pct_of_all_usp34_high", ascending=False)
    from src.gse245601_candidate_deepdive_visualization_1 import PATIENT_COLORS

    colors = [PATIENT_COLORS.get(p, "gray") for p in by_patient["group"]]
    ax.bar(range(len(by_patient)), by_patient["pct_of_all_usp34_high"], color=colors)
    ax.set_ylabel("% of all USP34-high Tam cells")
    ax.set_xticks(range(len(by_patient)))
    ax.set_xticklabels(by_patient["group"], rotation=45, ha="right", fontsize=7)
    ax.set_title("C. Which patients contribute\nUSP34-high Tam cells?", fontsize=9)

    fig.suptitle("USP34-high cells within 12h-Tamoxifen epithelial cells (DESCRIPTIVE -- not a claim about future resistant clones)", fontsize=10, y=1.03)
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info("wrote %s", out_path)


def run_usp34_high(config_path: str | Path = "config/config.yaml") -> dict[str, pd.DataFrame]:
    from src.gse245601_candidate_deepdive_data import load_per_cell_table

    config = _load_config(config_path)
    out_dir = Path(config["gse245601_candidate_deepdive"]["output"]["tables_dir"])
    figures_dir = Path(config["gse245601_candidate_deepdive"]["output"]["figures_dir"])

    per_cell = load_per_cell_table(config)
    detail, summary = build_usp34_high_distribution(per_cell)
    detail.to_csv(out_dir / "usp34_high_cell_distribution.tsv", sep="\t", index=False)
    summary.to_csv(out_dir / "usp34_high_cell_summary.tsv", sep="\t", index=False)

    plot_usp34_high_cells(per_cell, detail, figures_dir / "32_USP34_high_treated_cells.png")
    return {"detail": detail, "summary": summary}


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run_usp34_high()
