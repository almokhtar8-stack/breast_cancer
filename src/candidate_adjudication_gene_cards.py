"""Candidate adjudication Phase 11: one self-contained, factual summary
figure per MULTIMODAL_STRONG gene (CRISPR effect/FDR/rank, all four RNA
datasets' effect/FDR, dataset significance status, resistance direction
pattern, stability class, and a short factual interpretation sentence
built from the already-computed evidence -- no pathway/druggability
content, no mechanistic claim beyond what the data show).

Data source: `results/tables/candidate_adjudication/multimodal7_exact_evidence.tsv`
(built by `src.candidate_adjudication_master_table`, itself read from the
frozen cross-dataset genome-wide tables).
"""

from __future__ import annotations

import logging
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402
import yaml  # noqa: E402

logger = logging.getLogger(__name__)

RNA_DATASET_LABELS = [
    ("gse118713_log2fc", "gse118713_fdr", "GSE118713 (resistance)"),
    ("gse240112_tumor_log2fc", "gse240112_tumor_fdr", "GSE240112 (recurrence)"),
    ("gse111151_log2fc", "gse111151_fdr", "GSE111151 (resistance)"),
    ("gse245601_epi_log2fc", "gse245601_epi_fdr", "GSE245601 epi (acute)"),
    ("gse245601_malignant_log2fc", "gse245601_malignant_fdr", "GSE245601 malig. (acute)"),
]


def _load_config(config_path: str | Path) -> dict:
    with open(config_path) as f:
        return yaml.safe_load(f)


def _interpretation_sentence(row: pd.Series) -> str:
    direction_word = "sensitised (depleted)" if row["crispr_direction_class"] == "sensitising_KO" else "made more tolerant" if row["crispr_direction_class"] == "tolerance_associated_KO" else "not significantly affected"
    n_resistance_sig = int(row["resistance_fdr05_count"])
    return (
        f"Knockout {direction_word} cells under 4-OHT (CRISPR FDR={row['crispr_fdr']:.3f}). "
        f"{n_resistance_sig}/3 resistance-state RNA datasets reach FDR<0.05, direction pattern "
        f"{row['resistance_direction_pattern']} (consensus: {row['resistance_direction_consensus']}). "
        f"Global rank {int(row['global_rank'])}/15,255, {row['global_stability_class']}."
    )


def plot_gene_card(row: pd.Series, out_path: Path) -> None:
    """Three genuinely separate axes: CRISPR effect (its own numeric
    scale), RNA log2FC across the four RNA datasets (a different numeric
    scale -- CRISPR effect size and RNA log2FC are never drawn on one
    shared axis, per the task's explicit instruction), and a text summary."""
    fig, (ax_crispr, ax_rna, ax_text) = plt.subplots(1, 3, figsize=(11.5, 3.6), gridspec_kw={"width_ratios": [0.7, 1.3, 1]})

    ax_crispr.barh([0], [row["crispr_effect"]], color="#C44E52" if row["crispr_fdr"] < 0.05 else "#4C72B0")
    ax_crispr.set_yticks([0])
    ax_crispr.set_yticklabels(["CRISPR\n(functional)"], fontsize=8)
    ax_crispr.axvline(0, color="black", linewidth=0.8)
    ax_crispr.set_xlabel("CRISPR effect size", fontsize=7.5)
    ax_crispr.set_title("Functional", fontsize=9.5)
    ax_crispr.text(row["crispr_effect"], 0.35, f"fdr={row['crispr_fdr']:.3g}", ha="center", fontsize=6.5)

    labels = [lbl for _, _, lbl in RNA_DATASET_LABELS]
    effects = [row[e] for e, _, _ in RNA_DATASET_LABELS]
    fdrs = [row[f] for _, f, _ in RNA_DATASET_LABELS]
    colors = ["#4C72B0" if pd.isna(v) or v >= 0.05 else "#C44E52" for v in fdrs]
    y_pos = range(len(labels))
    ax_rna.barh(y_pos, [0 if pd.isna(e) else e for e in effects], color=colors)
    ax_rna.set_yticks(y_pos)
    ax_rna.set_yticklabels(labels, fontsize=8)
    ax_rna.axvline(0, color="black", linewidth=0.8)
    ax_rna.set_xlabel("log2FC", fontsize=7.5)
    ax_rna.set_title("RNA (resistance/recurrence/acute)\nred = FDR<0.05", fontsize=9.5)
    for yi, (e, f) in enumerate(zip(effects, fdrs)):
        if pd.notna(e):
            ax_rna.text(e + (0.05 if e >= 0 else -0.05), yi, f"fdr={f:.3g}" if pd.notna(f) else "NA", va="center", ha="left" if e >= 0 else "right", fontsize=6.5)

    fig.suptitle(f"{row['gene']}", fontsize=12, y=1.03)

    ax = ax_text
    ax.axis("off")
    summary = (
        f"Evidence category: {row['evidence_category']}\n"
        f"CRISPR direction: {row['crispr_direction_class']}\n"
        f"CRISPR rank: {int(row['crispr_rank_full_screen'])}/19,103\n"
        f"Global rank: {int(row['global_rank'])}/15,255 (coverage {row['coverage_tier']})\n"
        f"Datasets FDR<0.05: {int(row['n_datasets_fdr05'])}/5\n"
        f"Datasets top-10%%: {int(row['n_datasets_top10pct'])}/5\n"
        f"Stability: {row['global_stability_class']}\n"
        f"GSE111151 cell-line consistency: {int(row['gse111151_consistent_cell_lines_n'])}/{int(row['gse111151_consistent_cell_lines_total'])}\n"
        f"Resistance consensus: {row['resistance_direction_consensus']}\n"
    )
    ax.text(0.0, 1.0, summary, va="top", ha="left", fontsize=8.5, family="monospace", transform=ax.transAxes)
    ax.text(0.0, 0.02, _interpretation_sentence(row), va="bottom", ha="left", fontsize=7.5, wrap=True, transform=ax.transAxes, style="italic")

    fig.tight_layout()
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info("wrote %s", out_path)


def run_gene_cards(config_path: str | Path = "config/config.yaml") -> None:
    config = _load_config(config_path)
    tables_dir = Path(config["candidate_adjudication"]["output"]["tables_dir"])
    out_dir = Path(config["candidate_adjudication"]["output"]["gene_cards_dir"])
    df = pd.read_csv(tables_dir / "multimodal7_exact_evidence.tsv", sep="\t")
    for _, row in df.iterrows():
        plot_gene_card(row, out_dir / f"{row['gene']}.png")


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run_gene_cards()
