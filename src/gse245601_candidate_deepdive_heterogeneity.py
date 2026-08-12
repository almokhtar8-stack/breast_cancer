"""GSE245601 candidate deep-dive Phase 19: patient heterogeneity heatmaps
-- all-epithelial (reuses the Phase 5 direction-matrix figure/table
unchanged, no duplicate computation) and a second, malignant-only
version that additionally flags any patient/gene cell with <50 malignant
cells in either arm (below the frozen Track B eligibility threshold).
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


def plot_heterogeneity_heatmap(pb: pd.DataFrame, title: str, out_path: Path, flag_low_cell: bool = False) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from src.gse245601_candidate_deepdive_visualization_1 import PATIENT_ORDER

    genes = GENES
    grid = np.full((len(genes), len(PATIENT_ORDER)), np.nan)
    low_flag = np.zeros((len(genes), len(PATIENT_ORDER)), dtype=bool)
    for gi, gene in enumerate(genes):
        sub = pb.loc[pb["gene"] == gene].set_index("patient")
        for pi, patient in enumerate(PATIENT_ORDER):
            if patient not in sub.index or pd.isna(sub.loc[patient].get("log_fold_change_descriptive")):
                continue
            grid[gi, pi] = sub.loc[patient, "log_fold_change_descriptive"]
            if flag_low_cell:
                low_flag[gi, pi] = bool(sub.loc[patient].get("low_cell_count_warning_Control", False)) or bool(sub.loc[patient].get("low_cell_count_warning_Tamoxifen", False))

    fig, ax = plt.subplots(figsize=(8, 0.6 * len(genes) + 1.5))
    vmax = np.nanmax(np.abs(grid)) if not np.all(np.isnan(grid)) else 1.0
    im = ax.imshow(grid, cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
    for gi in range(len(genes)):
        for pi in range(len(PATIENT_ORDER)):
            v = grid[gi, pi]
            if np.isnan(v):
                ax.text(pi, gi, "NA", ha="center", va="center", fontsize=6.5, color="gray")
                continue
            arrow = "↑" if v > 0 else ("↓" if v < 0 else "→")
            text = arrow + ("†" if low_flag[gi, pi] else "")
            ax.text(pi, gi, text, ha="center", va="center", fontsize=10 if not low_flag[gi, pi] else 8, color="black")
    ax.set_xticks(range(len(PATIENT_ORDER)))
    ax.set_xticklabels(PATIENT_ORDER, rotation=45, ha="right", fontsize=8)
    ax.set_yticks(range(len(genes)))
    ax.set_yticklabels(genes, fontsize=9)
    fig.colorbar(im, ax=ax, label="descriptive log2FC (Tam − Control pseudobulk)", shrink=0.7)
    subtitle = "† = <50 malignant cells in Control or Tamoxifen arm (descriptive only)" if flag_low_cell else ""
    ax.set_title(f"{title}\n{subtitle}", fontsize=9.5)
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info("wrote %s", out_path)


def run_heterogeneity(config_path: str | Path = "config/config.yaml") -> None:
    config = _load_config(config_path)
    tables_dir = Path(config["gse245601_candidate_deepdive"]["output"]["tables_dir"])
    figures_dir = Path(config["gse245601_candidate_deepdive"]["output"]["figures_dir"])

    all_epi_pb = pd.read_csv(tables_dir / "patient_all_epithelial_pseudobulk.tsv", sep="\t")
    malignant_pb = pd.read_csv(tables_dir / "patient_malignant_pseudobulk.tsv", sep="\t")

    plot_heterogeneity_heatmap(all_epi_pb, "Patient heterogeneity: all-epithelial pseudobulk", figures_dir / "33_candidate_patient_heterogeneity.png")
    plot_heterogeneity_heatmap(malignant_pb, "Patient heterogeneity: strict-malignant pseudobulk", figures_dir / "34_candidate_patient_heterogeneity_malignant.png", flag_low_cell=True)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run_heterogeneity()
