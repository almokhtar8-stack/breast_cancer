"""Figure comparing the four USP34 shortest-path bridge genes (USP9X,
RPS27A, UBC, UBB) across the five evidence layers audited in
results/tables/systems_network/USP34_bridge_gene_evidence.tsv.

Reads only that already-written table; does not touch or regenerate any
other systems-network figure.
"""

from __future__ import annotations

import logging
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap

logger = logging.getLogger(__name__)

TABLE_PATH = Path("results/tables/systems_network/USP34_bridge_gene_evidence.tsv")
OUT_FIG = Path("results/figures/systems_network/final_review/08_USP34_bridge_gene_evidence.png")

GENE_ORDER = ["USP9X", "UBB", "RPS27A", "UBC"]  # B (partial support) first, then C (network-only)

COLUMNS = [
    ("crispr_effect", "crispr_p", "crispr_fdr", "CRISPR\neffect_size"),
    ("gse118713_log2fc", "gse118713_p", "gse118713_fdr", "GSE118713\n(TAMR vs MCF7)"),
    ("gse240112_log2fc", "gse240112_p", "gse240112_fdr", "GSE240112\n(recurrent vs\nprimary)"),
    ("gse111151_log2fc", "gse111151_p", "gse111151_fdr", "GSE111151\n(resistant vs\nparental)"),
    ("gse245601_track_a_log2fc", "gse245601_track_a_p", "gse245601_track_a_fdr", "GSE245601\nACUTE 12h\n(not resistance)"),
]

# diverging blue<->red, gray midpoint (dataviz skill reference palette)
DIVERGING_CMAP = LinearSegmentedColormap.from_list(
    "bridge_evidence_diverging",
    ["#1c5cab", "#f0efec", "#e34948"],
)
VABS = 1.25

CLASSIFICATION_LABEL = {
    "A_DATA_SUPPORTED_BRIDGE": "A: DATA-SUPPORTED",
    "B_PARTIAL_SUPPORT": "B: PARTIAL SUPPORT",
    "C_NETWORK_ONLY_GENERIC_BRIDGE": "C: NETWORK-ONLY",
}


def build_figure(df: pd.DataFrame, out_fig: Path = OUT_FIG) -> None:
    df = df.set_index("gene").loc[GENE_ORDER]

    n_rows, n_cols = len(GENE_ORDER), len(COLUMNS)
    values = np.zeros((n_rows, n_cols))
    for j, (val_col, _, _, _) in enumerate(COLUMNS):
        values[:, j] = df[val_col].values

    fig, ax = plt.subplots(figsize=(10.5, 5.2), dpi=200)
    fig.patch.set_facecolor("white")

    im = ax.imshow(values, cmap=DIVERGING_CMAP, vmin=-VABS, vmax=VABS, aspect="auto")

    # visual separator between the 3 resistance-dataset columns and CRISPR /
    # the acute-response column, mirroring the acute/resistance separation
    # used in 01_resistance_pathway_consensus.png
    ax.axvline(x=0.5, color="white", linewidth=3)
    ax.axvline(x=3.5, color="#111111", linewidth=2.5)

    for i in range(n_rows):
        for j, (val_col, p_col, fdr_col, _) in enumerate(COLUMNS):
            val = df[val_col].iloc[i]
            p = df[p_col].iloc[i]
            fdr = df[fdr_col].iloc[i]
            marker = "**" if fdr < 0.05 else ("*" if p < 0.05 else "")
            text_color = "white" if abs(val) > VABS * 0.45 else "#111111"
            ax.text(j, i, f"{val:+.2f}{marker}", ha="center", va="center", fontsize=10.5, color=text_color, fontweight="bold" if marker else "normal")

    ax.set_xticks(range(n_cols))
    ax.set_xticklabels([c[3] for c in COLUMNS], fontsize=9)
    ax.set_yticks(range(n_rows))
    ax.set_yticklabels(GENE_ORDER, fontsize=12, fontweight="bold")

    for i, gene in enumerate(GENE_ORDER):
        label = CLASSIFICATION_LABEL[df.loc[gene, "classification"]]
        ax.text(n_cols - 0.5 + 0.35, i, label, ha="left", va="center", fontsize=9.5, color="#111111", fontweight="bold")

    ax.set_xlim(-0.5, n_cols - 0.5)
    ax.tick_params(length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)

    cbar = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.28)
    cbar.set_label("effect / log2FC (blue = down or sensitising,\nred = up or tolerance-associated)", fontsize=8)
    cbar.ax.tick_params(labelsize=8)

    ax.set_title(
        "USP34 bridge-gene evidence audit -- CRISPR + 3 resistance datasets + GSE245601 acute (Track A)\n"
        "* nominal p<0.05   ** FDR<0.05 (none reach FDR<0.05 here)   |   right column = conservative classification",
        fontsize=10.5,
        loc="left",
    )

    fig.tight_layout()
    out_fig.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_fig, facecolor="white")
    plt.close(fig)
    logger.info("wrote %s", out_fig)


def run(table_path: Path = TABLE_PATH, out_fig: Path = OUT_FIG) -> None:
    df = pd.read_csv(table_path, sep="\t")
    build_figure(df, out_fig)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run()
