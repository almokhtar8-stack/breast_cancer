"""GSE245601 candidate deep-dive Phase 20: quantifies per-tumor
compositional shift (malignant fraction, cluster proportions) between
Control and Tamoxifen, then checks whether candidate pseudobulk effects
correlate with compositional change -- descriptive correlation only, no
causal claim.
"""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from src.gse245601_candidate_deepdive_data import GENES

logger = logging.getLogger(__name__)


def _load_config(config_path: str | Path) -> dict:
    with open(config_path) as f:
        return yaml.safe_load(f)


def build_composition_change(per_cell: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for patient, pgrp in per_cell.groupby("patient", observed=True):
        malignant_frac = pgrp.groupby("condition", observed=True)["malignancy_status"].apply(lambda s: (s == "malignant").mean())
        counts = pgrp.groupby(["seurat_clusters", "condition"], observed=True).size().unstack("condition", fill_value=0)
        for c in ("Control", "Tamoxifen"):
            if c not in counts.columns:
                counts[c] = 0
        control_props = counts["Control"] / counts["Control"].sum() if counts["Control"].sum() > 0 else counts["Control"] * 0.0
        tam_props = counts["Tamoxifen"] / counts["Tamoxifen"].sum() if counts["Tamoxifen"].sum() > 0 else counts["Tamoxifen"] * 0.0
        cluster_shift = float((tam_props - control_props).abs().sum()) / 2  # total variation distance
        rows.append(
            {
                "patient": patient,
                "malignant_fraction_control": malignant_frac.get("Control", np.nan),
                "malignant_fraction_tamoxifen": malignant_frac.get("Tamoxifen", np.nan),
                "malignant_fraction_change": malignant_frac.get("Tamoxifen", np.nan) - malignant_frac.get("Control", np.nan),
                "cluster_composition_total_variation_distance": cluster_shift,
            }
        )
    out = pd.DataFrame(rows)
    logger.info("build_composition_change: %d patients, mean cluster TV distance=%.3f", len(out), out["cluster_composition_total_variation_distance"].mean())
    return out


def correlate_composition_with_effect(composition: pd.DataFrame, all_epi_pb: pd.DataFrame, genes: list[str]) -> pd.DataFrame:
    from scipy.stats import spearmanr

    rows = []
    comp = composition.set_index("patient")
    for gene in genes:
        pb = all_epi_pb.loc[all_epi_pb["gene"] == gene].set_index("patient")
        merged = pb.join(comp, how="inner").dropna(subset=["log_fold_change_descriptive", "cluster_composition_total_variation_distance"])
        if len(merged) < 4:
            rows.append({"gene": gene, "n_patients": len(merged), "rho_vs_cluster_shift": np.nan, "rho_vs_malignant_fraction_change": np.nan, "note": "too few patients for a meaningful correlation"})
            continue
        rho_cluster, _ = spearmanr(merged["log_fold_change_descriptive"], merged["cluster_composition_total_variation_distance"])
        merged2 = merged.dropna(subset=["malignant_fraction_change"])
        rho_mal, _ = spearmanr(merged2["log_fold_change_descriptive"], merged2["malignant_fraction_change"]) if len(merged2) >= 4 else (np.nan, np.nan)
        rows.append({"gene": gene, "n_patients": len(merged), "rho_vs_cluster_shift": rho_cluster, "rho_vs_malignant_fraction_change": rho_mal, "note": "descriptive Spearman correlation, not a causal claim"})
    out = pd.DataFrame(rows)
    logger.info("correlate_composition_with_effect: %s", {r["gene"]: round(r["rho_vs_cluster_shift"], 2) if pd.notna(r["rho_vs_cluster_shift"]) else None for _, r in out.iterrows()})
    return out


def plot_composition_vs_effect(composition: pd.DataFrame, all_epi_pb: pd.DataFrame, genes: list[str], out_path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from src.gse245601_candidate_deepdive_visualization_1 import PATIENT_COLORS

    comp = composition.set_index("patient")
    fig, axes = plt.subplots(1, len(genes), figsize=(4 * len(genes), 4), sharey=True)
    for ax, gene in zip(axes, genes):
        pb = all_epi_pb.loc[all_epi_pb["gene"] == gene].set_index("patient")
        merged = pb.join(comp, how="inner").dropna(subset=["log_fold_change_descriptive", "cluster_composition_total_variation_distance"])
        for patient, row in merged.iterrows():
            ax.scatter(row["cluster_composition_total_variation_distance"], row["log_fold_change_descriptive"], color=PATIENT_COLORS.get(patient, "gray"), s=40)
        ax.axhline(0, color="gray", linewidth=0.6)
        ax.set_xlabel("cluster composition shift\n(total variation distance)", fontsize=8)
        ax.set_title(gene, fontsize=9)
    axes[0].set_ylabel("descriptive pseudobulk log2FC (all epithelial)")
    fig.suptitle("Composition shift vs. candidate pseudobulk effect, per tumor (descriptive correlation, not causal)", fontsize=10.5, y=1.03)
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info("wrote %s", out_path)


def run_composition(config_path: str | Path = "config/config.yaml") -> dict[str, pd.DataFrame]:
    from src.gse245601_candidate_deepdive_data import load_per_cell_table

    config = _load_config(config_path)
    out_dir = Path(config["gse245601_candidate_deepdive"]["output"]["tables_dir"])
    figures_dir = Path(config["gse245601_candidate_deepdive"]["output"]["figures_dir"])

    per_cell = load_per_cell_table(config)
    composition = build_composition_change(per_cell)
    composition.to_csv(out_dir / "composition_change.tsv", sep="\t", index=False)

    all_epi_pb = pd.read_csv(out_dir / "patient_all_epithelial_pseudobulk.tsv", sep="\t")
    correlation = correlate_composition_with_effect(composition, all_epi_pb, GENES)
    correlation.to_csv(out_dir / "composition_effect_correlation.tsv", sep="\t", index=False)

    plot_composition_vs_effect(composition, all_epi_pb, GENES, figures_dir / "35_composition_vs_candidate_effect.png")
    return {"composition": composition, "correlation": correlation}


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run_composition()
