"""Evidence freeze Phases 7-9, 15: the four final_review figures. RNA
columns are always shown in the canonical order GSE118713 | GSE240112 |
GSE111151 || GSE245601 (acute), with a visible divider separating the
three resistance/recurrence-context datasets from the acute 12h
tamoxifen-response dataset -- never implying they measure the same
biological state. CRISPR direction and RNA direction are always shown in
separate columns/panels, never merged into one symbol.

Data source: `results/tables/evidence_freeze/final_candidate_evidence.tsv`
(built by `src.evidence_freeze_tables`, itself read from the frozen
candidate-adjudication and cross-dataset genome-wide tables).
"""

from __future__ import annotations

import logging
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.patches as mpatches  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import yaml  # noqa: E402

logger = logging.getLogger(__name__)

# direction+significance -> (display color, category label)
_UP_SIG, _UP_NS = "#B2182B", "#F4A6A6"
_DOWN_SIG, _DOWN_NS = "#2166AC", "#A6C8E8"
_FLAT, _NA = "#E8E8E8", "#FFFFFF"


def _load_config(config_path: str | Path) -> dict:
    with open(config_path) as f:
        return yaml.safe_load(f)


def _cell_color(log2fc: float, fdr: float) -> str:
    if pd.isna(log2fc):
        return _NA
    if log2fc > 0:
        return _UP_SIG if (pd.notna(fdr) and fdr < 0.05) else _UP_NS
    if log2fc < 0:
        return _DOWN_SIG if (pd.notna(fdr) and fdr < 0.05) else _DOWN_NS
    return _FLAT


def _cell_text(log2fc: float, fdr: float) -> str:
    if pd.isna(log2fc):
        return "NA"
    arrow = "↑" if log2fc > 0 else ("↓" if log2fc < 0 else "→")
    return arrow + ("*" if (pd.notna(fdr) and fdr < 0.05) else "")


def plot_four_rna_direction_matrix(df: pd.DataFrame, out_path: Path) -> None:
    df = df.sort_values(by=["freeze_shortlisted", "global_rank"], ascending=[False, True], na_position="last").reset_index(drop=True)
    genes = df["gene"].tolist()
    cols = [("gse118713_log2fc", "gse118713_fdr"), ("gse240112_log2fc", "gse240112_fdr"), ("gse111151_log2fc", "gse111151_fdr"), ("gse245601_epi_log2fc", "gse245601_epi_fdr")]
    col_labels = ["GSE118713", "GSE240112", "GSE111151", "GSE245601\n(acute 12h)"]

    fig, ax = plt.subplots(figsize=(6.5, 0.32 * len(genes) + 2.2))
    for gi, gene in enumerate(genes):
        row = df.iloc[gi]
        for ci, (lc, fc) in enumerate(cols):
            x = ci if ci < 3 else ci + 0.6  # extra gap before the 4th (acute) column
            color = _cell_color(row[lc], row[fc])
            ax.add_patch(mpatches.Rectangle((x - 0.45, gi - 0.45), 0.9, 0.9, facecolor=color, edgecolor="black", linewidth=0.4))
            ax.text(x, gi, _cell_text(row[lc], row[fc]), ha="center", va="center", fontsize=8)

    divider_x = 2.5 + 0.3
    ax.axvline(divider_x, color="black", linewidth=1.5, linestyle="--")
    ax.text(1.0, len(genes) + 0.3, "RESISTANCE / RECURRENCE CONTEXT", ha="center", fontsize=9, fontweight="bold")
    ax.text(3.6, len(genes) + 0.3, "ACUTE 12h TAMOXIFEN", ha="center", fontsize=9, fontweight="bold", color="#7A3B00")

    xticks = [0, 1, 2, 3.6]
    ax.set_xticks(xticks)
    ax.set_xticklabels(col_labels, fontsize=8.5)
    ax.set_yticks(range(len(genes)))
    ax.set_yticklabels(genes, fontsize=8)
    ax.set_xlim(-0.6, 4.3)
    ax.set_ylim(-0.6, len(genes) + 0.8)
    ax.invert_yaxis()
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.tick_params(length=0)

    legend_handles = [
        mpatches.Patch(color=_UP_SIG, label="up, FDR<0.05"), mpatches.Patch(color=_UP_NS, label="up, not significant"),
        mpatches.Patch(color=_DOWN_SIG, label="down, FDR<0.05"), mpatches.Patch(color=_DOWN_NS, label="down, not significant"),
        mpatches.Patch(color=_NA, label="not testable (NA)", ec="black"),
    ]
    ax.legend(handles=legend_handles, loc="upper center", bbox_to_anchor=(0.5, -0.02 / max(len(genes), 1) - 0.06), ncol=3, fontsize=7.5, frameon=False)
    fig.suptitle("Four RNA datasets: resistance/recurrence context (1-3) vs. acute 12h tamoxifen response (4) -- never the same biological state", fontsize=9.5, y=1.02)
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info("wrote %s", out_path)


def plot_five_layer_support_matrix(df: pd.DataFrame, out_path: Path) -> None:
    df = df.sort_values(by=["freeze_shortlisted", "global_rank"], ascending=[False, True], na_position="last").reset_index(drop=True)
    genes = df["gene"].tolist()

    # CRISPR uses a deliberately different color pair (purple/orange) from
    # the RNA red/blue up-down pair -- CRISPR direction and RNA direction
    # are different quantities (functional-fitness phenotype vs. expression
    # association) and must never look like the same scale.
    _SENS_SIG, _SENS_NS = "#5E3C99", "#C4B3E0"
    _TOL_SIG, _TOL_NS = "#E08214", "#F5CB8C"

    fig, ax = plt.subplots(figsize=(9.2, 0.32 * len(genes) + 3.0))
    col_defs = [
        ("CRISPR", None), ("GSE118713", ("gse118713_log2fc", "gse118713_fdr")),
        ("GSE240112", ("gse240112_log2fc", "gse240112_fdr")), ("GSE111151", ("gse111151_log2fc", "gse111151_fdr")),
        ("GSE245601\n(acute)", ("gse245601_epi_log2fc", "gse245601_epi_fdr")),
    ]
    # positions: CRISPR | (gap) | GSE118713 GSE240112 GSE111151 | (bigger gap) | GSE245601 acute
    x_positions = [0, 1.5, 2.5, 3.5, 5.1]
    for gi, gene in enumerate(genes):
        row = df.iloc[gi]
        for ci, (label, cols) in enumerate(col_defs):
            x = x_positions[ci]
            if label == "CRISPR":
                if pd.isna(row["crispr_fdr"]):
                    color, text = _NA, "NA"
                elif row["crispr_direction"] == "sensitising_KO":
                    color = _SENS_SIG if row["crispr_fdr"] < 0.05 else _SENS_NS
                    text = "sens.KO" + ("*" if row["crispr_fdr"] < 0.05 else "")
                else:
                    color = _TOL_SIG if row["crispr_fdr"] < 0.05 else _TOL_NS
                    text = "tol.KO" + ("*" if row["crispr_fdr"] < 0.05 else "")
            else:
                lc, fc = cols
                color, text = _cell_color(row[lc], row[fc]), _cell_text(row[lc], row[fc])
            ax.add_patch(mpatches.Rectangle((x - 0.45, gi - 0.45), 0.9, 0.9, facecolor=color, edgecolor="black", linewidth=0.4))
            ax.text(x, gi, text, ha="center", va="center", fontsize=7.5)

    ax.axvline(0.75, color="black", linewidth=1.2, linestyle="--")  # CRISPR vs RNA
    divider_x = 4.3
    ax.axvline(divider_x, color="black", linewidth=1.5, linestyle="--")  # resistance vs acute
    ax.text(2.5, -1.3, "RESISTANCE / RECURRENCE CONTEXT", ha="center", fontsize=8.5, fontweight="bold")
    ax.text(5.1, -1.3, "ACUTE 12h\nTAMOXIFEN", ha="center", fontsize=8.5, fontweight="bold", color="#7A3B00")

    ax.set_xticks(x_positions)
    ax.set_xticklabels([c[0] for c in col_defs], fontsize=8.5)
    ax.set_yticks(range(len(genes)))
    ax.set_yticklabels(genes, fontsize=8)
    ax.set_xlim(-0.6, 5.7)
    ax.set_ylim(-1.9, len(genes) + 0.6)
    ax.invert_yaxis()
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.tick_params(length=0)

    legend_handles = [
        mpatches.Patch(color=_SENS_SIG, label="CRISPR sensitising_KO, FDR<0.05"), mpatches.Patch(color=_SENS_NS, label="CRISPR sensitising_KO, not sig."),
        mpatches.Patch(color=_TOL_SIG, label="CRISPR tolerance_KO, FDR<0.05"), mpatches.Patch(color=_TOL_NS, label="CRISPR tolerance_KO, not sig."),
        mpatches.Patch(color=_UP_SIG, label="RNA up, FDR<0.05"), mpatches.Patch(color=_UP_NS, label="RNA up, not sig."),
        mpatches.Patch(color=_DOWN_SIG, label="RNA down, FDR<0.05"), mpatches.Patch(color=_DOWN_NS, label="RNA down, not sig."),
        mpatches.Patch(color=_NA, label="not testable (NA)", ec="black"),
    ]
    ax.legend(handles=legend_handles, loc="upper center", bbox_to_anchor=(0.5, -0.1), ncol=3, fontsize=7, frameon=False)
    ax.set_title("Five-layer support matrix: CRISPR (distinct purple/orange scale, never the same as RNA red/blue) + 4 RNA datasets", fontsize=9.5, pad=14)
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info("wrote %s", out_path)


def plot_therapeutic_shortlist_head_to_head(df: pd.DataFrame, shortlist_genes: list[str], out_path: Path) -> None:
    sub = df.loc[df["gene"].isin(shortlist_genes)].copy()
    sub["_order"] = sub["gene"].map({g: i for i, g in enumerate(shortlist_genes)})
    sub = sub.sort_values("_order")
    if not (sub["crispr_direction"] == "sensitising_KO").all():
        raise ValueError("a non-sensitising gene was passed to the therapeutic shortlist head-to-head figure")

    genes = sub["gene"].tolist()
    fig, axes = plt.subplots(1, 5, figsize=(14, 0.55 * len(genes) + 1.5))

    ax = axes[0]
    ax.barh(genes, -np.log10(sub["crispr_fdr"]), color="#2166AC")
    ax.axvline(-np.log10(0.25), color="black", linestyle="--", linewidth=0.9)
    ax.set_title("CRISPR\n-log10(FDR)\n(all sensitising_KO)\ndashed line = eligibility gate, FDR=0.25", fontsize=7)
    ax.invert_yaxis()

    ax = axes[1]
    ax.barh(genes, sub["resistance_fdr05_count"], color="#55A868")
    ax.set_xlim(0, 3)
    ax.set_title("resistance datasets\nFDR<0.05 (of 3)", fontsize=8)
    ax.invert_yaxis()
    ax.set_yticklabels([])

    ax = axes[2]
    ax.barh(genes, (sub["human_tumor_support"] == "significant").astype(int), color="#DD8452")
    ax.set_xlim(0, 1.2)
    ax.set_xticks([0, 1])
    ax.set_title("human-tumor\nsupport", fontsize=8)
    ax.invert_yaxis()
    ax.set_yticklabels([])

    ax = axes[3]
    consistency_map = {"all_up": 3, "all_down": 3, "majority_up": 2, "majority_down": 2, "mixed": 1, "insufficient": 0}
    ax.barh(genes, sub["resistance_direction_consistency"].map(consistency_map), color="#8172B2")
    ax.set_xlim(0, 3.5)
    ax.set_title("resistance direction\nconsistency", fontsize=8)
    ax.invert_yaxis()
    ax.set_yticklabels([])

    ax = axes[4]
    stability_colors = {"ROBUST": "#2E7D32", "MODERATELY_STABLE": "#F9A825", "DATASET_DEPENDENT": "#C62828"}
    colors = [stability_colors.get(s, "gray") for s in sub["ranking_stability"]]
    ax.barh(genes, [1] * len(genes), color=colors)
    ax.set_xticks([])
    ax.set_title("leave-one-out\nstability", fontsize=8)
    ax.invert_yaxis()
    ax.set_yticklabels([])
    handles = [mpatches.Patch(color=c, label=k) for k, c in stability_colors.items()]
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.08), fontsize=6.5, ncol=1, frameon=False)

    fig.suptitle("Therapeutic (inhibition/sensitisation) shortlist head-to-head -- no overall-score bar", fontsize=11, y=1.02)
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info("wrote %s", out_path)


def plot_frozen_shortlist_summary(df: pd.DataFrame, shortlist_genes: list[str], out_path: Path) -> None:
    sub = df.loc[df["gene"].isin(shortlist_genes)].copy()
    sub["_order"] = sub["gene"].map({g: i for i, g in enumerate(shortlist_genes)})
    sub = sub.sort_values("_order").reset_index(drop=True)

    fig, ax = plt.subplots(figsize=(9, 0.9 * len(sub) + 1.5))
    ax.axis("off")
    y0 = 1.0
    dy = 1.0 / (len(sub) + 0.5)
    for i, row in sub.iterrows():
        y = y0 - (i + 0.5) * dy
        text = (
            f"{i + 1}. {row['gene']}   CRISPR: {row['crispr_direction']} (effect={row['crispr_effect']:.2f}, FDR={row['crispr_fdr']:.3g})\n"
            f"     RNA: {row['full_rna_pattern_4']}   (* = FDR<0.05; first 3 = resistance/recurrence, 4th = acute 12h Tam)\n"
            f"     Resistance sig count: {row['resistance_fdr05_count']}/3   Human support: {row['human_tumor_support']}   Stability: {row['ranking_stability']}"
        )
        ax.text(0.02, y, text, fontsize=9.5, va="center", family="monospace", transform=ax.transAxes)
        if i < len(sub) - 1:
            ax.axhline(y - dy / 2, xmin=0.0, xmax=1.0, color="#DDDDDD", linewidth=0.8)
    fig.suptitle("Frozen therapeutic (inhibition/sensitisation) shortlist -- a nonsignificant arrow is NOT confirmed differential expression", fontsize=10.5, y=1.0)
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info("wrote %s", out_path)


def run_visualization(config_path: str | Path = "config/config.yaml", shortlist_genes: list[str] | None = None) -> None:
    config = _load_config(config_path)
    ef = config["evidence_freeze"]["output"]
    tables_dir = Path(ef["tables_dir"])
    final_review = Path(ef["final_review_dir"])
    final_review.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(tables_dir / "final_candidate_evidence.tsv", sep="\t")
    if shortlist_genes is None:
        # THERAPEUTIC_SHORTLIST_FREEZE.tsv (built by src.evidence_freeze_shortlist_freeze)
        # is the single source of truth for shortlist membership AND order --
        # never re-derived from df's own row order, which is sorted for display,
        # not necessarily by freeze_shortlist_rank
        freeze_manifest = pd.read_csv(tables_dir / "THERAPEUTIC_SHORTLIST_FREEZE.tsv", sep="\t").sort_values("freeze_rank")
        shortlist_genes = freeze_manifest["gene"].tolist()

    plot_four_rna_direction_matrix(df, final_review / "01_four_rna_direction_matrix.png")
    plot_five_layer_support_matrix(df, final_review / "02_five_layer_support_matrix.png")
    plot_therapeutic_shortlist_head_to_head(df, shortlist_genes, final_review / "03_therapeutic_shortlist_head_to_head.png")
    plot_frozen_shortlist_summary(df, shortlist_genes, final_review / "04_frozen_shortlist_summary.png")


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run_visualization()
