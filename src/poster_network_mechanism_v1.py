"""ONE poster-grade network/mechanism figure -- answers "which local
gene/protein mechanism context supports each of the four focus genes, and
where do plausible indirect targeting routes exist?"

Data sources, all read unmodified (no new network build, no new GSEA,
no recomputation of any edge, score, NES or FDR):
  - `results/tables/systems_network/four_candidate_direct_neighbors.tsv`
    and `four_candidate_shortest_paths.tsv` -- the project's existing,
    already-frozen STRING-derived local network for USP34 (and, for
    VEZF1, its one frozen pathway-co-membership edge). Built for an
    earlier four-candidate shortlist (USP34, VEZF1, EML5, CITED2, see
    `docs/SYSTEMS_NETWORK_NODE_RULE.md` Category A) -- KDM1A and TLK2 were
    never part of that network build (KDM1A was a blind positive control
    at the time; TLK2 was not on the shortlist at all), so NO frozen
    PPI/network table exists for either gene. This module does not
    invent one.
  - `results/tables/systems_network/gsea_crispr.tsv` -- the frozen GSEA
    run on the genome-wide CRISPR ranking. Read only to check leading-edge
    *pathway membership* for KDM1A and TLK2 (real, frozen, already-computed
    values), which is a different and weaker kind of evidence than a
    network edge and is presented as such (rounded-box nodes, never drawn
    as a gene-gene interaction).
  - `src/post_audit_sensitivity_data.load_significant_sensitising_hits()`
    -- the same already-tested CRISPR loader used by the CRISPR discovery
    figure, for each candidate's rank/effect/FDR.

No candidate is given a network it does not have in the frozen tables.
"""

from __future__ import annotations

import logging
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Circle, FancyBboxPatch

from src import poster_exploration_v2_data as ed
from src import post_audit_sensitivity_data as pad

logger = logging.getLogger(__name__)

OUT_DIR = Path("results/figures/poster_network_mechanism_v1")
NET_DIR = Path("results/tables/systems_network")

FOCUS_FOUR = ["USP34", "VEZF1", "KDM1A", "TLK2"]  # panel order = evidence-strength order
FOCUS_COLORS = ed.FOCUS_COLORS  # frozen Okabe-Ito identity colors, reused unmodified

DGRAY = "#262626"
MGRAY = "#8c8c8c"
LGRAY = "#d8d8d8"
PARTNER_COLOR = "#4C6B8A"   # neutral network-partner blue-gray -- never a candidate identity color
PATHWAY_COLOR = "#B5793A"   # neutral warm tone for pathway-membership boxes


def _load_crispr_stats() -> dict[str, dict]:
    """Rank/effect/FDR per focus gene, from the already-tested sensitising-
    hit loader -- identical source used by poster_crispr_discovery_v1."""
    hits = pad.load_significant_sensitising_hits()
    out = {}
    for _, row in hits[hits["gene"].isin(FOCUS_FOUR)].iterrows():
        out[row["gene"]] = {
            "rank": int(row["rank_by_effect"]),
            "effect": float(row["effect_size"]),
            "fdr": float(row["fdr"]),
        }
    return out


def _load_usp34_network() -> tuple[pd.DataFrame, pd.DataFrame]:
    """USP34's frozen direct STRING partners and its frozen 2-hop shortest
    paths, read unmodified from the existing systems-network tables."""
    neighbors = pd.read_csv(NET_DIR / "four_candidate_direct_neighbors.tsv", sep="\t")
    neighbors = neighbors[neighbors["candidate"] == "USP34"].reset_index(drop=True)
    paths = pd.read_csv(NET_DIR / "four_candidate_shortest_paths.tsv", sep="\t")
    paths = paths[(paths["candidate"] == "USP34") & (paths["path_rank"] == 1.0)].reset_index(drop=True)
    return neighbors, paths


def _load_vezf1_network() -> pd.DataFrame:
    neighbors = pd.read_csv(NET_DIR / "four_candidate_direct_neighbors.tsv", sep="\t")
    return neighbors[neighbors["candidate"] == "VEZF1"].reset_index(drop=True)


# KDM1A and TLK2 have no frozen network table, so the only real evidence
# available is GSEA leading-edge pathway membership in the CRISPR-ranked
# enrichment. Rather than a blind lowest-FDR cut (which surfaces generic,
# non-interpretable hits such as REACTOME_SARS_COV_INFECTIONS -- a real
# but mechanistically uninformative leading-edge membership, present only
# because of broadly-shared transcriptional-machinery genes), pathways are
# hand-selected for mechanistic interpretability from the full leading-edge
# survey below, matching the same curation approach already used for
# Panel C of poster_pathway_v1. NES/FDR values themselves are always read
# from the frozen table, never hand-typed.
PATHWAYS_KDM1A = ["REACTOME_CHROMATIN_MODIFYING_ENZYMES", "REACTOME_ESTROGEN_DEPENDENT_GENE_EXPRESSION"]
PATHWAYS_TLK2 = ["GOBP_REGULATION_OF_CHROMATIN_ORGANIZATION"]

PATHWAY_SHORT_LABEL = {
    "REACTOME_CHROMATIN_MODIFYING_ENZYMES": "Chromatin-\nmodifying\nenzymes",
    "REACTOME_ESTROGEN_DEPENDENT_GENE_EXPRESSION": "Estrogen-\ndependent gene\nexpression",
    "GOBP_REGULATION_OF_CHROMATIN_ORGANIZATION": "Regulation of\nchromatin organization",
}


def _full_leading_edge_membership(gene: str) -> pd.DataFrame:
    """Every frozen CRISPR-ranked-GSEA pathway where `gene` is leading-edge,
    sorted by FDR -- used only to document/verify the curated shortlist
    above, never to pick nodes automatically."""
    gsea = pd.read_csv(NET_DIR / "gsea_crispr.tsv", sep="\t")
    pattern = re.compile(rf"(^|;){re.escape(gene)}(;|$)")
    mask = gsea["leading_edge_genes"].fillna("").apply(lambda s: bool(pattern.search(s)))
    return gsea[mask].sort_values("fdr").reset_index(drop=True)


def _load_pathway_membership(gene: str, pathway_names: list[str]) -> pd.DataFrame:
    """The curated, named pathways above, with NES/FDR read unmodified
    from the frozen CRISPR-ranked GSEA table."""
    full = _full_leading_edge_membership(gene)
    hits = full[full["pathway"].isin(pathway_names)].copy()
    hits["_order"] = hits["pathway"].apply(pathway_names.index)
    return hits.sort_values("_order").reset_index(drop=True)


def _fmt_fdr(fdr: float) -> str:
    """Avoids a misleading rounded 'FDR=0.000' for small-but-nonzero FDRs."""
    return "FDR<0.001" if fdr < 0.001 else f"FDR={fdr:.3f}"


# ---------------------------------------------------------------------------
# Generic radial node-graph renderer
# ---------------------------------------------------------------------------

def _draw_node(ax, xy, kind: str, label: str, color: str, radius: float, fontsize: float, textcolor: str = "white"):
    x, y = xy
    if kind == "box":
        n_lines = label.count("\n") + 1
        max_chars = max(len(line) for line in label.split("\n"))
        w, h = max(1.05, 0.088 * max_chars), 0.28 + 0.26 * n_lines
        patch = FancyBboxPatch((x - w / 2, y - h / 2), w, h, boxstyle="round,pad=0.02,rounding_size=0.08",
                                facecolor=color, edgecolor="white", linewidth=1.3, zorder=4)
    else:
        patch = Circle(xy, radius=radius, facecolor=color, edgecolor="white", linewidth=1.3, zorder=4)
    ax.add_patch(patch)
    ax.text(x, y, label, ha="center", va="center", fontsize=fontsize, fontweight="bold",
             color=textcolor, zorder=5)
    return patch


def _draw_edge(ax, xy0, xy1, style: str, linewidth: float, color: str = LGRAY):
    dashes = (0, (1,)) if style == "solid" else (0, (4, 3))
    ax.plot([xy0[0], xy1[0]], [xy0[1], xy1[1]], color=color, linewidth=linewidth,
             linestyle="-" if style == "solid" else "--", zorder=2, solid_capstyle="round")


def _panel_frame(ax, title: str, subtitle: str, title_color: str, ylim_bottom: float = -1.9):
    ax.set_xlim(-1.75, 1.75)
    ax.set_ylim(ylim_bottom, 0.65)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title(title, loc="left", fontsize=15.5, fontweight="bold", color=title_color, pad=6)
    ax.text(0.0, 1.22, subtitle, transform=ax.transAxes, fontsize=10.5, color=MGRAY,
             ha="left", va="top", style="italic")


def _place(angle_deg: float, radius: float) -> tuple[float, float]:
    a = np.deg2rad(angle_deg)
    return radius * np.cos(a), radius * np.sin(a)


def _panel_usp34(ax, crispr_stats: dict) -> None:
    neighbors, paths = _load_usp34_network()

    c_stat = crispr_stats["USP34"]
    subtitle = f"CRISPR sensitiser, rank {c_stat['rank']}/13 · {_fmt_fdr(c_stat['fdr'])}"
    _panel_frame(ax, "USP34", subtitle, FOCUS_COLORS["USP34"])

    _draw_node(ax, (0, 0), "circle", "USP34", FOCUS_COLORS["USP34"], 0.30, 12.5)

    # Direct, project-supported STRING partners (solid edges, real).
    direct = {
        "RPS27A": {"angle": 200, "radius": 0.95},
        "UBC": {"angle": 270, "radius": 0.95},
        "USP9X": {"angle": 340, "radius": 0.95},
    }
    for gene, pos in direct.items():
        row = neighbors[neighbors["neighbor_gene"] == gene].iloc[0]
        xy = _place(pos["angle"], pos["radius"])
        width = 4.5 if row["confidence"] >= 0.9 else 2.6
        _draw_edge(ax, (0, 0), xy, "solid", width, PARTNER_COLOR)
        _draw_node(ax, xy, "circle", gene, PARTNER_COLOR, 0.24, 9.5)

    # Indirect, 2-hop bridge targets (dashed edges) -- the plausible
    # indirect route into WNT/stemness biology, read from the frozen
    # rank-1 shortest paths.
    bridges = {"CTNNB1": ("RPS27A", 220, 1.5), "SOX2": ("USP9X", 320, 1.5)}
    for target, (via, angle, radius) in bridges.items():
        via_xy = _place(direct[via]["angle"], direct[via]["radius"])
        xy = _place(angle, radius)
        row = paths[(paths["source_gene"] == via) & (paths["target_gene"] == target)]
        if row.empty:
            continue
        _draw_edge(ax, via_xy, xy, "dashed", 2.0, MGRAY)
        patch = _draw_node(ax, xy, "circle", target, "white", 0.22, 9.0, textcolor=DGRAY)
        patch.set_edgecolor(PARTNER_COLOR)
        patch.set_linewidth(1.6)

    ax.text(0, -1.62, "solid = direct STRING partner · dashed = 2-hop bridge to a WNT/stemness effector",
             ha="center", va="top", fontsize=8.3, color=MGRAY)


def _panel_vezf1(ax, crispr_stats: dict) -> None:
    neighbors = _load_vezf1_network()
    row = neighbors.iloc[0]  # DMTN, the single frozen edge

    c_stat = crispr_stats["VEZF1"]
    subtitle = f"CRISPR sensitiser, rank {c_stat['rank']}/13 · {_fmt_fdr(c_stat['fdr'])}"
    _panel_frame(ax, "VEZF1", subtitle, FOCUS_COLORS["VEZF1"])

    _draw_node(ax, (0, 0), "circle", "VEZF1", FOCUS_COLORS["VEZF1"], 0.30, 12.5)

    xy = _place(0, 1.0)
    _draw_edge(ax, (0, 0), xy, "solid", 4.5, PARTNER_COLOR)
    _draw_node(ax, xy, "circle", row["neighbor_gene"], PARTNER_COLOR, 0.24, 9.5)

    ax.text(0, -0.62, "shared leading-edge pathway:\nHALLMARK_HEME_METABOLISM /\nBLOOD_VESSEL_MORPHOGENESIS",
             ha="center", va="top", fontsize=8.6, color=DGRAY)
    ax.text(0, -1.16,
            "DMTN is itself human recurrence-associated\n(up in gse118713 & gse240112, both FDR<0.001)",
             ha="center", va="top", fontsize=8.3, color=MGRAY, style="italic")
    ax.text(0, -1.62, "single frozen network edge -- honestly minimal, not padded",
             ha="center", va="top", fontsize=8.3, color=MGRAY)


def _panel_pathway_concept(ax, gene: str, crispr_stats: dict, pathway_names: list[str], shared_note: str) -> None:
    pw = _load_pathway_membership(gene, pathway_names)
    c_stat = crispr_stats[gene]
    subtitle = f"CRISPR sensitiser, rank {c_stat['rank']}/13 · {_fmt_fdr(c_stat['fdr'])}"
    _panel_frame(ax, gene, subtitle, FOCUS_COLORS[gene], ylim_bottom=-2.15)

    _draw_node(ax, (0, 0), "circle", gene, FOCUS_COLORS[gene], 0.30, 12.5)

    n = len(pw)
    angles, radius = ([240, 300], 1.2) if n == 2 else ([270], 1.1)
    for i, (_, prow) in enumerate(pw.iterrows()):
        xy = _place(angles[i], radius)
        _draw_edge(ax, (0, 0), xy, "dashed", 2.0, MGRAY)
        label = PATHWAY_SHORT_LABEL[prow["pathway"]]
        _draw_node(ax, xy, "box", label, PATHWAY_COLOR, 0.34, 7.3)

    ax.text(0, -1.78, "box = frozen GSEA leading-edge pathway membership (CRISPR ranking), not a network edge",
             ha="center", va="top", fontsize=7.8, color=MGRAY)
    ax.text(0, -1.96, shared_note, ha="center", va="top", fontsize=8.3, color=MGRAY)


def build_network_mechanism_main(stub: Path) -> None:
    crispr_stats = _load_crispr_stats()

    fig = plt.figure(figsize=(14.5, 8.5), dpi=300)
    gs = fig.add_gridspec(2, 2, hspace=0.3, wspace=0.1)
    ax_usp34 = fig.add_subplot(gs[0, 0])
    ax_vezf1 = fig.add_subplot(gs[0, 1])
    ax_kdm1a = fig.add_subplot(gs[1, 0])
    ax_tlk2 = fig.add_subplot(gs[1, 1])

    _panel_usp34(ax_usp34, crispr_stats)
    _panel_vezf1(ax_vezf1, crispr_stats)
    _panel_pathway_concept(ax_kdm1a, "KDM1A", crispr_stats, PATHWAYS_KDM1A,
                            shared_note="no PPI/network partners in frozen analysis -- pathway membership only")
    _panel_pathway_concept(ax_tlk2, "TLK2", crispr_stats, PATHWAYS_TLK2,
                            shared_note="no PPI/network partners in frozen analysis -- shares its one pathway hit with KDM1A")

    fig.text(0.03, 0.985, "Local mechanism context is not equal across the four candidates",
              fontsize=23, fontweight="bold", color=DGRAY, ha="left", va="top")
    fig.text(0.03, 0.958,
              "Frozen STRING network evidence exists only for USP34 and VEZF1; KDM1A and TLK2 have real CRISPR-ranked pathway leading-edge membership but no frozen PPI network.",
              fontsize=11.3, color="#555555", ha="left", va="top")

    legend_handles = [
        Line2D([0], [0], color=PARTNER_COLOR, linewidth=4.5, label="direct STRING partner"),
        Line2D([0], [0], color=MGRAY, linewidth=2.0, linestyle="--", label="indirect (2-hop bridge / pathway membership)"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor=PARTNER_COLOR, markeredgecolor="white",
               markersize=11, label="gene node"),
        Line2D([0], [0], marker="s", color="none", markerfacecolor=PATHWAY_COLOR, markeredgecolor="white",
               markersize=11, label="pathway-membership node"),
    ]
    fig.legend(handles=legend_handles, loc="lower center", ncol=4, frameon=False, fontsize=10,
               labelcolor=DGRAY, bbox_to_anchor=(0.5, 0.0))

    fig.subplots_adjust(left=0.04, right=0.98, top=0.81, bottom=0.075)

    stub.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(stub.with_suffix(".png"), facecolor="white", bbox_inches="tight", dpi=300)
    fig.savefig(stub.with_suffix(".pdf"), facecolor="white", bbox_inches="tight")
    fig.savefig(stub.with_suffix(".svg"), facecolor="white", bbox_inches="tight")
    plt.close(fig)
    logger.info("wrote %s.png/.pdf/.svg", stub)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    build_network_mechanism_main(OUT_DIR / "NETWORK_mechanism_main")
