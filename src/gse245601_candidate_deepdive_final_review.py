"""GSE245601 candidate deep-dive Phase 23: curates the final_review
directory with only the most useful figures -- copies the best single
-gene USP34 figures already built, and builds three new compact
all-4-gene summary figures (malignancy response, prevalence/intensity,
acute summary) so a reviewer does not have to open all ~35 exploratory
figures.
"""

from __future__ import annotations

import logging
import shutil
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402
import yaml  # noqa: E402

from src.gse245601_candidate_deepdive_data import GENES

logger = logging.getLogger(__name__)


def _load_config(config_path: str | Path) -> dict:
    with open(config_path) as f:
        return yaml.safe_load(f)


COPY_MAP = {
    "01_USP34_patient_pseudobulk.png": "01_USP34_patient_pseudobulk_all_epithelial.png",
    "02_USP34_malignant_vs_nonmalignant.png": "10_USP34_malignant_vs_nonmalignant_response.png",
    "03_USP34_prevalence_vs_intensity.png": "14_USP34_prevalence_vs_intensity.png",
    "04_USP34_decomposition.png": "22_USP34_decomposition.png",
    "05_USP34_cluster_response.png": "28_USP34_cluster_response.png",
    "06_candidate_patient_heterogeneity.png": "33_candidate_patient_heterogeneity.png",
    "07_candidate_patient_heterogeneity_malignant.png": "34_candidate_patient_heterogeneity_malignant.png",
}


def build_all4_malignancy_response(malig_cond: pd.DataFrame, out_path: Path) -> None:
    from src.gse245601_candidate_deepdive_visualization_1 import PATIENT_COLORS, PATIENT_ORDER

    fig, axes = plt.subplots(2, len(GENES), figsize=(4 * len(GENES), 7), sharey="row")
    for gi, gene in enumerate(GENES):
        for ri, (malignancy, label) in enumerate([("non-malignant epithelial", "Non-malignant"), ("malignant", "Malignant")]):
            ax = axes[ri, gi]
            m = malig_cond.loc[(malig_cond["gene"] == gene) & (malig_cond["malignancy_status"] == malignancy)]
            wide = m.pivot_table(index="patient", columns="condition", values="pseudobulk_normalized_expression")
            n_cells_wide = m.pivot_table(index="patient", columns="condition", values="n_cells")
            for patient in PATIENT_ORDER:
                if patient not in wide.index or wide.loc[patient].isna().any():
                    continue
                # same <50-cell reliability threshold used everywhere else in this deep-dive
                # (config gse245601_candidate_deepdive.malignant_cell_min_count) -- only meaningful
                # for the malignant row; the non-malignant compartment always has ample cells
                low = malignancy == "malignant" and n_cells_wide.loc[patient].min() < 50
                ax.plot(
                    [0, 1], [wide.loc[patient, "Control"], wide.loc[patient, "Tamoxifen"]], color=PATIENT_COLORS[patient],
                    marker="o", markersize=4, linewidth=1.1, alpha=0.4 if low else 0.85, linestyle="--" if low else "-",
                )
            ax.set_xticks([0, 1])
            ax.set_xticklabels(["Ctrl", "Tam"], fontsize=7)
            if ri == 0:
                ax.set_title(gene, fontsize=10)
            if gi == 0:
                label_suffix = "\n(dashed = <50 malignant cells\nin that arm)" if malignancy == "malignant" else ""
                ax.set_ylabel(f"{label}\nlog2(CPM+1){label_suffix}", fontsize=7.5)
    fig.suptitle("All four candidates: malignant vs. non-malignant response (each line = one tumor's pseudobulk pair)", fontsize=11, y=1.02)
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info("wrote %s", out_path)


def build_all4_prevalence_intensity(prevalence: pd.DataFrame, out_path: Path) -> None:
    from src.gse245601_candidate_deepdive_visualization_1 import PATIENT_COLORS, PATIENT_ORDER

    fig, axes = plt.subplots(2, len(GENES), figsize=(4 * len(GENES), 7))
    for gi, gene in enumerate(GENES):
        sub = prevalence.loc[(prevalence["gene"] == gene) & (prevalence["malignancy_status"] == "all_epithelial")]
        ax = axes[0, gi]
        wide = sub.pivot_table(index="patient", columns="condition", values="fraction_expressing")
        for patient in PATIENT_ORDER:
            if patient not in wide.index or wide.loc[patient].isna().any():
                continue
            ax.plot([0, 1], [100 * wide.loc[patient, "Control"], 100 * wide.loc[patient, "Tamoxifen"]], color=PATIENT_COLORS[patient], marker="o", markersize=4, linewidth=1.1, alpha=0.85)
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["Ctrl", "Tam"], fontsize=7)
        ax.set_title(gene, fontsize=10)
        if gi == 0:
            ax.set_ylabel("%% expressing\n(prevalence)", fontsize=8)

        ax = axes[1, gi]
        wide2 = sub.pivot_table(index="patient", columns="condition", values="mean_normalized_positive_cells_only")
        for patient in PATIENT_ORDER:
            if patient not in wide2.index or wide2.loc[patient].isna().any():
                continue
            ax.plot([0, 1], [wide2.loc[patient, "Control"], wide2.loc[patient, "Tamoxifen"]], color=PATIENT_COLORS[patient], marker="o", markersize=4, linewidth=1.1, alpha=0.85)
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["Ctrl", "Tam"], fontsize=7)
        if gi == 0:
            ax.set_ylabel("mean expr., positive\ncells only (intensity)", fontsize=8)
    fig.suptitle("All four candidates: prevalence vs. intensity (all epithelial, patient-level)", fontsize=11, y=1.02)
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info("wrote %s", out_path)


def build_all4_acute_summary(frozen: pd.DataFrame, direction: pd.DataFrame, out_path: Path) -> None:
    import numpy as np

    fig, axes = plt.subplots(1, 3, figsize=(13, 4))
    f = frozen.set_index("gene").loc[GENES]

    ax = axes[0]
    colors = ["#C44E52" if v < 0.05 else "#8CA0C6" for v in f["gse245601_epi_fdr"]]
    ax.barh(GENES, f["gse245601_epi_log2fc"], color=colors)
    ax.axvline(0, color="black", linewidth=0.8)
    for i, gene in enumerate(GENES):
        ax.text(f.loc[gene, "gse245601_epi_log2fc"], i, f"  FDR={f.loc[gene, 'gse245601_epi_fdr']:.2g}", va="center", fontsize=7)
    ax.set_title("Frozen Track A\n(all epithelial) log2FC", fontsize=9)
    ax.set_xlabel("log2FC (Tam vs Control)")

    ax = axes[1]
    up = direction.set_index("gene")["n_patients_increase"]
    down = direction.set_index("gene")["n_patients_decrease"]
    y = np.arange(len(GENES))
    ax.barh(y, [up[g] for g in GENES], color="#B2182B", label="increase")
    ax.barh(y, [-down[g] for g in GENES], color="#2166AC", label="decrease")
    ax.set_yticks(y)
    ax.set_yticklabels(GENES)
    ax.axvline(0, color="black", linewidth=0.8)
    ax.set_xlabel("n tumors (of 10)")
    ax.set_title("Patient-level direction\n(all-epithelial pseudobulk)", fontsize=9)
    ax.legend(fontsize=7, frameon=False)

    ax = axes[2]
    colors_b = ["#C44E52" if pd.notna(v) and v < 0.05 else "#8CA0C6" for v in f["gse245601_malignant_fdr"]]
    ax.barh(GENES, f["gse245601_malignant_log2fc"], color=colors_b)
    ax.axvline(0, color="black", linewidth=0.8)
    for i, gene in enumerate(GENES):
        fdr = f.loc[gene, "gse245601_malignant_fdr"]
        ax.text(f.loc[gene, "gse245601_malignant_log2fc"], i, f"  FDR={fdr:.2g}" if pd.notna(fdr) else "  NA", va="center", fontsize=7)
    ax.set_title("Frozen Track B\n(strict malignant) log2FC", fontsize=9)
    ax.set_xlabel("log2FC (Tam vs Control)")

    fig.suptitle("All four candidates: acute 12h summary (red = FDR<0.05)", fontsize=11, y=1.03)
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info("wrote %s", out_path)


def run_final_review(config_path: str | Path = "config/config.yaml") -> None:
    config = _load_config(config_path)
    tables_dir = Path(config["gse245601_candidate_deepdive"]["output"]["tables_dir"])
    figures_dir = Path(config["gse245601_candidate_deepdive"]["output"]["figures_dir"])
    final_review = Path(config["gse245601_candidate_deepdive"]["output"]["final_review_dir"])
    final_review.mkdir(parents=True, exist_ok=True)

    for dest_name, src_name in COPY_MAP.items():
        src = figures_dir / src_name
        if not src.exists():
            raise FileNotFoundError(f"expected source figure missing: {src}")
        shutil.copyfile(src, final_review / dest_name)
        logger.info("copied %s -> final_review/%s", src_name, dest_name)

    malig_cond = pd.read_csv(tables_dir / "malignancy_condition_patient_summary.tsv", sep="\t")
    prevalence = pd.read_csv(tables_dir / "expression_prevalence_intensity.tsv", sep="\t")
    frozen = pd.read_csv("results/tables/evidence_freeze/final_candidate_evidence.tsv", sep="\t")
    direction = pd.read_csv(tables_dir / "patient_direction_summary.tsv", sep="\t")

    build_all4_malignancy_response(malig_cond, final_review / "08_all4_malignancy_response.png")
    build_all4_prevalence_intensity(prevalence, final_review / "09_all4_prevalence_intensity.png")
    build_all4_acute_summary(frozen, direction, final_review / "10_all4_acute_summary.png")


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run_final_review()
