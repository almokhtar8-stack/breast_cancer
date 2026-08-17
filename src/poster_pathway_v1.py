"""ONE poster-grade pathway-biology figure -- answers "what biological
programs distinguish endocrine resistance/recurrence from acute
tamoxifen response?"

Data source: `poster_exploration_v2_data.load_pathway_trajectories()`,
called unmodified. That function reads the already-frozen
`results/tables/systems_network/gsea_{dataset}.tsv` GSEA output tables
(one row per Hallmark/Reactome/GO_BP pathway per dataset, produced by the
project's existing, already-frozen GSEA run) and returns exactly the
already-computed NES/FDR for the requested pathways -- this module
performs NO new enrichment analysis, NO re-ranking, and NO recomputation
of any NES or FDR value. Every number on the figure (NES positions,
significance markers, pathway counts) is read from that table at render
time, never hand-typed.
"""

from __future__ import annotations

import logging
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from src import poster_exploration_v2_data as ed

logger = logging.getLogger(__name__)

OUT_DIR = Path("results/figures/poster_pathway_v1")

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
# Resistance-model contexts get related warm tones, recurrence a distinct
# muted purple, acute a distinct cool teal -- the same hues already used
# for the equivalent biological states in the sample-level hero heatmap
# (poster_hero_heatmap_v4.STATE_COLORS), for visual continuity across the
# poster.
CONTEXT_COLORS = {
    "gse118713": "#C2A24F",
    "gse111151": "#B5793A",
    "gse240112": "#7C6F9E",
    "gse245601": "#4F8FA6",
}

DGRAY = "#262626"
MGRAY = "#8c8c8c"
BLUE_DARK = "#2E6C8E"    # Estrogen Response -- Early; also EMT-negative (acute)
BLUE_LIGHT = "#7FAFC4"   # Estrogen Response -- Late
TERRACOTTA = "#C1543A"   # EMT-positive (resistance / recurrence)

FDR_THRESHOLD = 0.05  # project convention, e.g. src/systems_network_node_attributes.py

PATHWAYS_A = [
    ("hallmark", "HALLMARK_ESTROGEN_RESPONSE_EARLY", "Estrogen Response -- Early"),
    ("hallmark", "HALLMARK_ESTROGEN_RESPONSE_LATE", "Estrogen Response -- Late"),
]
PATHWAY_B = ("hallmark", "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "EMT")
PATHWAYS_C = [
    ("hallmark", "HALLMARK_TGF_BETA_SIGNALING", "TGF-β signaling"),
    ("hallmark", "HALLMARK_APICAL_JUNCTION", "Apical junction"),
    ("hallmark", "HALLMARK_E2F_TARGETS", "E2F targets"),
]


def _wide(pathways: list[tuple[str, str, str]]) -> dict[str, dict[str, tuple[float, float]]]:
    """pathway_label -> dataset -> (NES, fdr), read from the frozen GSEA
    tables via the already-tested `load_pathway_trajectories` loader."""
    long = ed.load_pathway_trajectories(pathways)
    out: dict[str, dict[str, tuple[float, float]]] = {}
    for _, row in long.iterrows():
        out.setdefault(row["pathway_label"], {})[row["dataset"]] = (row["NES"], row["fdr"])
    return out


def build_pathway_main(stub: Path) -> None:
    wide_a = _wide(PATHWAYS_A)
    wide_b = _wide([PATHWAY_B])
    wide_c = _wide(PATHWAYS_C)

    x = np.arange(len(DATASET_ORDER))
    acute_idx = DATASET_ORDER.index("gse245601")

    fig = plt.figure(figsize=(13.5, 12.6), dpi=300)
    gs = fig.add_gridspec(3, 1, height_ratios=[1.05, 1.05, 0.62], hspace=0.55)
    axA = fig.add_subplot(gs[0])
    axB = fig.add_subplot(gs[1])
    axC = fig.add_subplot(gs[2])

    # ---- Panel A: Estrogen Response Early/Late trajectories ----
    for label, color in zip(["Estrogen Response -- Early", "Estrogen Response -- Late"], [BLUE_DARK, BLUE_LIGHT]):
        vals = np.array([wide_a[label][d][0] for d in DATASET_ORDER])
        axA.plot(x, vals, color=color, linewidth=3.2, zorder=3, solid_capstyle="round")
        axA.scatter(x, vals, s=230, color=color, zorder=4, edgecolor="white", linewidth=1.4, label=label)

    axA.axhline(0, color="#c9c9c9", linewidth=1.2, zorder=1)
    axA.set_ylabel("NES", fontsize=12, color=DGRAY)
    axA.set_title("A   Estrogen response is suppressed across every context", loc="left", fontsize=14.5,
                  fontweight="bold", color=DGRAY, pad=14)
    axA.legend(loc="lower left", frameon=False, fontsize=10.5, labelcolor=DGRAY)

    # ---- Panel B: EMT lollipop, colored by sign ----
    b_vals = np.array([wide_b["EMT"][d][0] for d in DATASET_ORDER])
    b_colors = [TERRACOTTA if v > 0 else BLUE_DARK for v in b_vals]
    axB.vlines(x, 0, b_vals, color=b_colors, linewidth=6, zorder=3)
    axB.scatter(x, b_vals, s=340, color=b_colors, zorder=4, edgecolor="white", linewidth=1.6)
    axB.axhline(0, color="#c9c9c9", linewidth=1.2, zorder=1)
    axB.set_ylabel("NES", fontsize=12, color=DGRAY)
    axB.set_title("B   EMT reverses direction between resistance/recurrence and acute response", loc="left",
                  fontsize=14.5, fontweight="bold", color=DGRAY, pad=14)

    # ---- shared "acute is different" background band on A and B --
    # axvspan blends x-data with y-axes-fraction, so it never distorts
    # y autoscale (unlike a Rectangle drawn in raw data coordinates).
    for ax in (axA, axB):
        ax.axvspan(acute_idx - 0.5, acute_idx + 0.5, facecolor="#EAF2F5", edgecolor="none", zorder=0)

    for ax in (axA, axB):
        ax.set_xlim(-0.55, len(DATASET_ORDER) - 0.45)
        ax.set_xticks(x)
        ax.set_xticklabels([CONTEXT_LABEL[d] for d in DATASET_ORDER], fontsize=11.5, color=DGRAY)
        for spine in ("top", "right", "bottom"):
            ax.spines[spine].set_visible(False)
        ax.spines["left"].set_color("#c9c9c9")
        ax.tick_params(axis="x", length=0)
        ax.tick_params(axis="y", labelsize=10, colors=DGRAY)

    axB.set_xlabel(" · ".join(CONTEXT_ACCESSION[d] for d in DATASET_ORDER),
                   fontsize=8.6, color=MGRAY, style="italic", labelpad=14)

    # ---- Panel C: compact supporting-pathway dot strip, context-colored ----
    # Small fixed per-context vertical offset (never applied to NES/x) so
    # near-identical NES values (e.g. TGF-beta: gse111151 vs gse240112,
    # 1.46 vs 1.47) don't fully occlude one another -- a standard
    # beeswarm-style separation, not a data change.
    n_ctx = len(DATASET_ORDER)
    y_jitter = (np.arange(n_ctx) - (n_ctx - 1) / 2) * 0.09
    c_labels = [label for _, _, label in PATHWAYS_C]
    yC = np.arange(len(c_labels))
    for i, label in enumerate(c_labels):
        vals = np.array([wide_c[label][d][0] for d in DATASET_ORDER])
        fdrs = np.array([wide_c[label][d][1] for d in DATASET_ORDER])
        axC.plot(vals, [yC[i]] * len(vals), color="#dddddd", linewidth=1.3, zorder=1)
        for j, d in enumerate(DATASET_ORDER):
            sig = fdrs[j] < FDR_THRESHOLD
            axC.scatter([vals[j]], [yC[i] + y_jitter[j]], s=120, zorder=3,
                        facecolor=CONTEXT_COLORS[d] if sig else "white",
                        edgecolor=CONTEXT_COLORS[d], linewidth=1.8)

    axC.axvline(0, color="#c9c9c9", linewidth=1.2, zorder=0)
    axC.set_yticks(yC)
    axC.set_yticklabels(c_labels, fontsize=11.5, color=DGRAY)
    axC.invert_yaxis()
    axC.set_xlabel("NES", fontsize=11, color=DGRAY)
    axC.set_title("C   Supporting programs reinforce the same pattern", loc="left", fontsize=13,
                  fontweight="bold", color=DGRAY, pad=12)
    for spine in ("top", "right", "left"):
        axC.spines[spine].set_visible(False)
    axC.spines["bottom"].set_color("#c9c9c9")
    axC.tick_params(axis="y", length=0)
    axC.tick_params(axis="x", labelsize=9.5, colors=DGRAY)

    handles = [plt.Line2D([0], [0], marker="o", color="none", markerfacecolor=CONTEXT_COLORS[d],
                           markeredgecolor=CONTEXT_COLORS[d], markersize=9,
                           label=CONTEXT_LABEL[d].replace("\n", " ")) for d in DATASET_ORDER]
    open_handle = plt.Line2D([0], [0], marker="o", color="none", markerfacecolor="white",
                              markeredgecolor=MGRAY, markersize=9, label=f"open = FDR ≥ {FDR_THRESHOLD:g}")
    axC.legend(handles=handles + [open_handle], loc="upper center", bbox_to_anchor=(0.5, -0.38),
               ncol=3, frameon=False, fontsize=9, labelcolor=DGRAY, handletextpad=0.5, columnspacing=1.3)

    fig.text(0.03, 0.985, "Pathway remodeling separates resistance from acute tamoxifen response",
              fontsize=23, fontweight="bold", color=DGRAY, ha="left", va="top")
    fig.text(0.03, 0.958, "Normalized enrichment scores from the frozen pathway analyses.",
              fontsize=12, color="#555555", ha="left", va="top")

    fig.subplots_adjust(left=0.075, right=0.98, top=0.88, bottom=0.06)

    stub.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(stub.with_suffix(".png"), facecolor="white", bbox_inches="tight", dpi=300)
    fig.savefig(stub.with_suffix(".pdf"), facecolor="white", bbox_inches="tight")
    fig.savefig(stub.with_suffix(".svg"), facecolor="white", bbox_inches="tight")
    plt.close(fig)
    logger.info("wrote %s.png/.pdf/.svg", stub)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    build_pathway_main(OUT_DIR / "PATHWAY_main")
