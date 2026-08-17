"""ONE poster-grade pathway figure (v2) -- an orthogonal transcriptomic
check of the biological programs implicated by the candidate-network
section: "do network-relevant programs actually change in the
transcriptomic resistance/recurrence/acute models?"

Data source: `poster_exploration_v2_data.load_pathway_trajectories()`,
called unmodified -- it reads the already-frozen
`results/tables/systems_network/gsea_{dataset}.tsv` GSEA tables. This
module performs NO new enrichment, NO re-ranking, and NO recomputation of
any NES or FDR; every plotted value is read from the frozen tables at
render time, never hand-typed.

Pathway selection is theme-first, not aesthetics-first: the displayed
gene sets are the pre-specified Hallmark readouts of the biological
themes raised by the network/mechanism section (estrogen response, EMT,
WNT/beta-catenin, E2F/cell cycle), kept regardless of effect size --
WNT is included even where weak/non-significant because it was a
pre-specified network theme. HALLMARK_G2M_CHECKPOINT (the optional
cell-cycle companion) is omitted as directionally redundant with E2F
targets and non-significant in two contexts; the "chromatin regulation"
network theme has no single comparable Hallmark readout across the four
contexts and none was invented (both documented in NOTE.md).

The figure asserts nothing causal: network + pathway agreement supports
a hypothesis, not candidate -> pathway -> resistance causality.
"""

from __future__ import annotations

import logging
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.lines import Line2D

from src import poster_exploration_v2_data as ed

logger = logging.getLogger(__name__)

OUT_DIR = Path("results/figures/poster_pathway_v2")

DATASET_ORDER = ed.DATASET_ORDER  # ["gse118713", "gse111151", "gse240112", "gse245601"]

CONTEXT_LABEL = {
    "gse118713": "Cell-line\nresistance model",
    "gse111151": "Independent\nresistant sublines",
    "gse240112": "Primary vs\nrecurrent tumours",
    "gse245601": "Acute 12 h\ntamoxifen",
}
CONTEXT_ACCESSION = {
    "gse118713": "GSE118713", "gse111151": "GSE111151",
    "gse240112": "GSE240112", "gse245601": "GSE245601",
}

DGRAY = "#262626"
MGRAY = "#8c8c8c"
BLUE = "#2E6C8E"        # suppressed / negative NES -- heatmap palette endpoint
TERRACOTTA = "#C1543A"  # enriched / positive NES -- heatmap palette endpoint
ACUTE_BAND = "#EAF2F5"

FDR_THRESHOLD = 0.05  # project convention, same as poster_pathway_v1

# Pre-specified network themes -> the exact comparable frozen Hallmark
# readout for each. Fixed from the network-section biology BEFORE
# inspecting plot aesthetics; each row of the figure must map to one of
# these themes (tested).
NETWORK_THEME_PATHWAYS = {
    "Estrogen response": [
        ("hallmark", "HALLMARK_ESTROGEN_RESPONSE_EARLY", "Estrogen response — early"),
        ("hallmark", "HALLMARK_ESTROGEN_RESPONSE_LATE", "Estrogen response — late"),
    ],
    "EMT": [
        ("hallmark", "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "EMT"),
    ],
    "WNT / beta-catenin": [
        ("hallmark", "HALLMARK_WNT_BETA_CATENIN_SIGNALING", "WNT / β-catenin"),
    ],
    "E2F / cell cycle": [
        ("hallmark", "HALLMARK_E2F_TARGETS", "E2F targets"),
    ],
    # "Chromatin regulation": no single comparable Hallmark gene set exists
    # across all four frozen tables -- intentionally no entry (NOTE.md #9).
}
PATHWAYS = [p for theme in NETWORK_THEME_PATHWAYS.values() for p in theme]

NES_CMAP = LinearSegmentedColormap.from_list("nes_div", [BLUE, "#F4F1EC", TERRACOTTA])
NES_NORM = Normalize(vmin=-3.0, vmax=3.0)


def _size(nes: float) -> float:
    """Marker area scales with |NES|."""
    return 130.0 + 240.0 * abs(nes)


def load_matrix() -> dict[str, dict[str, tuple[float, float]]]:
    """pathway_label -> dataset -> (NES, fdr), read from the frozen GSEA
    tables via the already-tested loader."""
    long = ed.load_pathway_trajectories(PATHWAYS)
    out: dict[str, dict[str, tuple[float, float]]] = {}
    for _, row in long.iterrows():
        out.setdefault(row["pathway_label"], {})[row["dataset"]] = (float(row["NES"]), float(row["fdr"]))
    return out


def build_pathway_v2(stub: Path) -> None:
    matrix = load_matrix()
    row_labels = [label for _, _, label in PATHWAYS]

    x = np.arange(len(DATASET_ORDER))
    y = np.arange(len(row_labels))
    acute_idx = DATASET_ORDER.index("gse245601")

    fig, ax = plt.subplots(figsize=(11.8, 7.6), dpi=300)

    ax.axvspan(acute_idx - 0.5, acute_idx + 0.5, facecolor=ACUTE_BAND, edgecolor="none", zorder=0)
    ax.axvline(acute_idx - 0.5, color="#b9c6cd", linewidth=1.2, zorder=1)

    for i, label in enumerate(row_labels):
        ax.axhline(i, color="#eeeeee", linewidth=0.8, zorder=1)
        for j, ds in enumerate(DATASET_ORDER):
            nes, fdr = matrix[label][ds]
            color = NES_CMAP(NES_NORM(nes))
            significant = fdr < FDR_THRESHOLD
            ax.scatter([j], [i], s=_size(nes), zorder=3,
                        facecolor=color if significant else "white",
                        edgecolor=color if significant else NES_CMAP(NES_NORM(nes)),
                        linewidth=1.4 if significant else 2.4)

    ax.set_xlim(-0.55, len(DATASET_ORDER) - 0.45)
    ax.set_ylim(len(row_labels) - 0.4, -0.6)  # first pathway on top
    ax.set_xticks(x)
    ax.set_xticklabels([CONTEXT_LABEL[d] for d in DATASET_ORDER], fontsize=12.5, color=DGRAY)
    ax.set_yticks(y)
    ax.set_yticklabels(row_labels, fontsize=13.5, color=DGRAY)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.tick_params(length=0)

    # Small gray accession IDs BELOW the context tick labels (blended
    # transform: x in data coords, y in axes fraction -- never collides
    # with the tick text).
    import matplotlib.transforms as mtransforms
    blend = mtransforms.blended_transform_factory(ax.transData, ax.transAxes)
    for j, ds in enumerate(DATASET_ORDER):
        ax.text(j, -0.102, CONTEXT_ACCESSION[ds], transform=blend, fontsize=8.2, color=MGRAY,
                 ha="center", va="top", style="italic")

    # "Resistance / recurrence" vs "Acute response" grouping cue above columns.
    ax.text(1.0, -0.52, "Resistance / recurrence-associated", fontsize=10.5, color=MGRAY,
             ha="center", va="bottom")
    ax.text(acute_idx, -0.52, "Acute response", fontsize=10.5, color=MGRAY, ha="center", va="bottom")

    fig.text(0.055, 0.975, "Network-relevant programs show distinct pathway remodeling",
              fontsize=19.5, fontweight="bold", color=DGRAY, ha="left", va="top")
    fig.text(0.055, 0.928, "Frozen GSEA results across resistance, recurrence and acute treatment-response contexts.",
              fontsize=11, color="#555555", ha="left", va="top")

    # --- compact legend row under the matrix -------------------------------
    # Color meaning: Suppressed <- 0 -> Enriched (NES), plus size + fill keys.
    cax = fig.add_axes([0.10, 0.075, 0.20, 0.025])
    gradient = np.linspace(-3, 3, 256).reshape(1, -1)
    cax.imshow(gradient, aspect="auto", cmap=NES_CMAP, norm=NES_NORM, extent=[-3, 3, 0, 1])
    cax.set_yticks([])
    cax.set_xticks([-3, 0, 3])
    cax.set_xticklabels(["-3", "0", "+3"], fontsize=8.5, color=DGRAY)
    for spine in cax.spines.values():
        spine.set_color("#cccccc")
    cax.text(-3.25, 0.5, "Suppressed", fontsize=9.5, color=BLUE, ha="right", va="center", fontweight="bold")
    cax.text(3.25, 0.5, "Enriched", fontsize=9.5, color=TERRACOTTA, ha="left", va="center", fontweight="bold")
    cax.set_title("Normalized Enrichment Score (NES)", fontsize=9, color=MGRAY, pad=4)

    size_handles = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor="#b0b0b0", markeredgecolor="none",
               markersize=np.sqrt(_size(nes)) * 0.9, label=f"|NES| = {nes:g}")
        for nes in (1.0, 2.0, 3.0)
    ]
    fill_handles = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor="#b0b0b0", markeredgecolor="none",
               markersize=11, label=f"FDR < {FDR_THRESHOLD:g}"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor="white", markeredgecolor="#b0b0b0",
               markersize=11, label=f"FDR ≥ {FDR_THRESHOLD:g} (open)"),
    ]
    fig.legend(handles=size_handles + fill_handles, loc="lower center", frameon=False, fontsize=9,
               labelcolor=DGRAY, ncol=5, bbox_to_anchor=(0.68, 0.045), handletextpad=0.4, columnspacing=1.3)

    fig.subplots_adjust(left=0.24, right=0.97, top=0.845, bottom=0.19)

    stub.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(stub.with_suffix(".png"), facecolor="white", bbox_inches="tight", dpi=300)
    fig.savefig(stub.with_suffix(".pdf"), facecolor="white", bbox_inches="tight")
    fig.savefig(stub.with_suffix(".svg"), facecolor="white", bbox_inches="tight")
    plt.close(fig)
    logger.info("wrote %s.png/.pdf/.svg", stub)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    build_pathway_v2(OUT_DIR / "PATHWAY_v2")
