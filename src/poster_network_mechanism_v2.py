"""POST-FREEZE EXPLORATORY network/mechanism analysis and figure (v2) --
answers "what molecular neighborhood surrounds each of the four current
poster candidates (KDM1A, TLK2, USP34, VEZF1), and do their local networks
converge on biological programs relevant to tamoxifen response?"

This is a NEW, post-freeze exploratory analysis, not a revision of any
frozen result. It replaces the asymmetric v1 figure (which reused an
older network built for a different candidate shortlist) with a single
consistently-built network: the SAME STRING interaction source, organism,
confidence threshold, query depth and filtering rule applied identically
to all four current candidates.

Data sources, all read unmodified from local files written by
`scripts/download_string_network_v2_four_focus.py` (STRING string-db.org,
`interaction_partners` endpoint, species=9606, required_score=700 -- the
same "high confidence" band and endpoint already used by the project's
earlier network build in `src/systems_network_build.py`):
  - `data/reference/interactions/string_v2_level1_functional.tsv` /
    `..._physical.tsv` -- Level-1 direct partners of the four candidates.
  - `data/reference/interactions/string_v2_level2_functional.tsv` /
    `..._physical.tsv` -- partners of the (uncapped) Level-1 partner set,
    used only to discover real Level-2 bridges, never to invent edges.

No STRING edge, score, or partner list is fabricated or hand-typed. No
frozen candidate ranking, CRISPR result, or prior figure is modified.
"""

from __future__ import annotations

import logging
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
from adjustText import adjust_text
from matplotlib.lines import Line2D
from matplotlib.patches import Wedge

from src import poster_exploration_v2_data as ed

logger = logging.getLogger(__name__)

OUT_DIR = Path("results/figures/poster_network_mechanism_v2")
STRING_DIR = Path("data/reference/interactions")
GENESET_DIR = Path("data/reference/genesets")

CANDIDATES = ["KDM1A", "TLK2", "USP34", "VEZF1"]
FOCUS_COLORS = ed.FOCUS_COLORS  # frozen Okabe-Ito identity colors, reused unmodified

REQUIRED_SCORE = 0.7   # STRING "high confidence" band, 0-1 scale (required_score=700 in the API call)
LEVEL1_CAP = 12        # objective cap on displayed Level-1 partners per candidate, same K for all four

# Low-information / uncharacterized STRING entries excluded from display
# (not a real, interpretable gene symbol).
EXCLUDED_NODES = {"H7C0V5_HUMAN"}

# Pre-specified (before the Level-2 query was filtered) canonical marker
# genes for the biological programs named in the task: estrogen/endocrine
# signaling, EMT, WNT/beta-catenin, cell-cycle/E2F. A Level-2 node is kept
# only if it lands in this list or already exists in the Level-1 node set
# (a real bridge/convergence hit) -- nothing here was chosen after seeing
# which genes would result.
CANONICAL_BIOLOGY_GENES = {
    "ESR1", "ESR2", "GREB1", "FOXA1", "PGR",
    "CTNNB1", "APC", "AXIN1", "GSK3B", "TCF7L2",
    "CDH1", "VIM", "SNAI1", "SNAI2", "ZEB1", "TWIST1",
    "RB1", "E2F1", "CCND1", "CCNE1", "CDK4", "CDK6",
}

# Pathway/program membership sources -- already-frozen, locally-cached
# MSigDB GMT files used elsewhere in the project (e.g.
# `poster_exploration_v2_data.load_gene_set`), read unmodified.
PATHWAY_SETS = {
    "Estrogen response": [("hallmark.gmt", "HALLMARK_ESTROGEN_RESPONSE_EARLY"), ("hallmark.gmt", "HALLMARK_ESTROGEN_RESPONSE_LATE")],
    "EMT": [("hallmark.gmt", "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION")],
    "WNT / beta-catenin": [("hallmark.gmt", "HALLMARK_WNT_BETA_CATENIN_SIGNALING")],
    "E2F / cell cycle": [("hallmark.gmt", "HALLMARK_E2F_TARGETS")],
    "Chromatin regulation": [("reactome.gmt", "REACTOME_CHROMATIN_MODIFYING_ENZYMES")],
}
PATHWAY_COLORS = {
    "Estrogen response": "#2E6C8E",
    "EMT": "#C1543A",
    "WNT / beta-catenin": "#6A5AA8",
    "E2F / cell cycle": "#3C8C5B",
    "Chromatin regulation": "#B5793A",
}

DGRAY = "#262626"
MGRAY = "#8c8c8c"
PARTNER_COLOR = "#8296AC"
BRIDGE_COLOR = "#B7C4D2"


def _load_gmt_set(fname: str, pathway: str) -> set[str]:
    with open(GENESET_DIR / fname) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if parts[0] == pathway:
                return set(parts[2:])
    raise KeyError(f"{pathway} not found in {fname}")


def _read_string_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t")


def _label_physical(edges: pd.DataFrame, physical_path: Path) -> pd.DataFrame:
    """Labels each functional-network edge as physical_PPI if the exact
    pair also appears in STRING's own physical-network query result for
    the same node set -- identical labeling rule already used in
    `src/systems_network_build.py`."""
    physical = _read_string_tsv(physical_path)
    physical_pairs = set(zip(physical["preferredName_A"], physical["preferredName_B"])) | \
        set(zip(physical["preferredName_B"], physical["preferredName_A"]))
    edges = edges.copy()
    edges["interaction_type"] = [
        "physical_PPI" if (a, b) in physical_pairs else "functional_association"
        for a, b in zip(edges["preferredName_A"], edges["preferredName_B"])
    ]
    return edges


def load_level1_edges() -> pd.DataFrame:
    """Each candidate's direct STRING partners, capped at the top
    LEVEL1_CAP by score (deterministic alphabetical tie-break), same rule
    for all four candidates. Candidates with fewer than LEVEL1_CAP
    partners (TLK2: 5, USP34: 10, VEZF1: 0) show all of them -- the cap
    only ever removes low-ranking partners, never pads a sparse
    candidate."""
    functional = _read_string_tsv(STRING_DIR / "string_v2_level1_functional.tsv")
    functional = functional[~functional["preferredName_B"].isin(EXCLUDED_NODES)]
    capped = []
    for candidate in CANDIDATES:
        sub = functional[functional["preferredName_A"] == candidate]
        sub = sub.sort_values(["score", "preferredName_B"], ascending=[False, True]).head(LEVEL1_CAP)
        capped.append(sub)
    capped = pd.concat(capped, ignore_index=True) if capped else functional.iloc[0:0]
    return _label_physical(capped, STRING_DIR / "string_v2_level1_physical.tsv")


def load_level2_edges(level1_edges: pd.DataFrame) -> pd.DataFrame:
    """Real STRING edges from the capped Level-1 partner set to (a) any
    node already in the graph -- a genuine bridge/convergence edge -- or
    (b) a pre-specified canonical biology gene. Source genes are
    restricted to the capped Level-1 partner set so every displayed
    bridge originates from a node that is itself displayed."""
    existing_nodes = set(CANDIDATES) | set(level1_edges["preferredName_B"])
    level1_partner_set = set(level1_edges["preferredName_B"])

    functional = _read_string_tsv(STRING_DIR / "string_v2_level2_functional.tsv")
    functional = functional[functional["preferredName_A"].isin(level1_partner_set)]
    functional = functional[~functional["preferredName_B"].isin(EXCLUDED_NODES)]
    functional = functional[functional["preferredName_A"] != functional["preferredName_B"]]

    is_convergence = functional["preferredName_B"].isin(existing_nodes)
    is_canonical = functional["preferredName_B"].isin(CANONICAL_BIOLOGY_GENES)
    bridges = functional[is_convergence | is_canonical].copy()
    bridges["bridge_kind"] = np.where(
        bridges["preferredName_B"].isin(existing_nodes), "convergence_within_existing_nodes", "canonical_biology_gene"
    )

    # Deduplicate by unordered pair, keep the highest-scoring direction.
    bridges["_pair"] = bridges.apply(lambda r: tuple(sorted([r["preferredName_A"], r["preferredName_B"]])), axis=1)
    bridges = bridges.sort_values("score", ascending=False).drop_duplicates("_pair").drop(columns="_pair")

    return _label_physical(bridges, STRING_DIR / "string_v2_level2_physical.tsv")


def _pathway_membership(gene: str) -> list[str]:
    hits = []
    for label, sources in PATHWAY_SETS.items():
        for fname, pathway in sources:
            if gene in _load_gmt_set(fname, pathway):
                hits.append(label)
                break
    return hits


def build_network() -> nx.Graph:
    level1 = load_level1_edges()
    level2 = load_level2_edges(level1)

    graph = nx.Graph()
    for candidate in CANDIDATES:
        graph.add_node(candidate, kind="candidate")
    for _, row in level1.iterrows():
        graph.add_node(row["preferredName_B"], kind="level1_partner")
        graph.add_edge(row["preferredName_A"], row["preferredName_B"], score=row["score"],
                        interaction_type=row["interaction_type"], hop=1)
    for _, row in level2.iterrows():
        a, b = row["preferredName_A"], row["preferredName_B"]
        if b not in graph:
            graph.add_node(b, kind="level2_bridge")
        if not graph.has_edge(a, b):
            graph.add_edge(a, b, score=row["score"], interaction_type=row["interaction_type"], hop=2)

    for node in graph.nodes:
        graph.nodes[node]["pathways"] = _pathway_membership(node)

    degree = dict(graph.degree())
    betweenness = nx.betweenness_centrality(graph, weight=None)
    nx.set_node_attributes(graph, degree, "degree")
    nx.set_node_attributes(graph, betweenness, "betweenness")

    return graph


def network_stats(graph: nx.Graph) -> dict:
    components = list(nx.connected_components(graph))
    candidate_pairs = {}
    for i, a in enumerate(CANDIDATES):
        for b in CANDIDATES[i + 1:]:
            if a not in graph or b not in graph:
                candidate_pairs[(a, b)] = None
                continue
            try:
                candidate_pairs[(a, b)] = nx.shortest_path_length(graph, a, b)
            except nx.NetworkXNoPath:
                candidate_pairs[(a, b)] = None

    shared_neighbors = {}
    for i, a in enumerate(CANDIDATES):
        for b in CANDIDATES[i + 1:]:
            if a not in graph or b not in graph:
                continue
            shared = set(graph.neighbors(a)) & set(graph.neighbors(b))
            if shared:
                shared_neighbors[(a, b)] = shared

    return {
        "n_nodes": graph.number_of_nodes(),
        "n_edges": graph.number_of_edges(),
        "n_components": len(components),
        "components": components,
        "candidate_shortest_paths": candidate_pairs,
        "shared_direct_neighbors": shared_neighbors,
        "top_hubs_by_degree": sorted(graph.degree, key=lambda x: -x[1])[:8],
        "top_bridges_by_betweenness": sorted(nx.get_node_attributes(graph, "betweenness").items(), key=lambda x: -x[1])[:8],
    }


# ---------------------------------------------------------------------------
# Rendering -- one genuine graph-layout canvas, not a manually-drawn diagram
# ---------------------------------------------------------------------------

def _node_size(kind: str, degree: int, betweenness: float) -> float:
    base = {"candidate": 2600, "level1_partner": 420, "level2_bridge": 190}[kind]
    hub_boost = 1.0 + min(degree / 12.0, 1.4) + min(betweenness * 6.0, 1.2)
    return base * (hub_boost if kind != "candidate" else 1.0 + min(betweenness * 2.0, 0.3))


def _node_color(node: str, kind: str) -> str:
    if kind == "candidate":
        return FOCUS_COLORS[node]
    if kind == "level1_partner":
        return PARTNER_COLOR
    return BRIDGE_COLOR


def _should_label(node: str, data: dict) -> bool:
    if data["kind"] == "candidate":
        return True
    if data["kind"] == "level1_partner":
        return True
    return node in CANONICAL_BIOLOGY_GENES or data["degree"] >= 6


def build_network_mechanism_main(stub: Path) -> nx.Graph:
    graph = build_network()
    weight = {(u, v): max(0.15, 1.0 - d["score"]) for u, v, d in graph.edges(data=True)}
    nx.set_edge_attributes(graph, weight, "layout_weight")
    pos = nx.kamada_kawai_layout(graph, weight="layout_weight", scale=1.0)

    fig, ax = plt.subplots(figsize=(24, 17.5), dpi=300)
    ax.set_facecolor("white")
    fig.patch.set_facecolor("white")
    fig.subplots_adjust(left=0.01, right=0.99, top=0.94, bottom=0.02)

    # Edges: physical_PPI solid + slightly thicker, functional_association
    # thin solid -- no dashed "fake shortcut" edges; a Level-2 relationship
    # is simply drawn via its real intermediate node.
    for u, v, d in graph.edges(data=True):
        lw = (1.7 if d["interaction_type"] == "physical_PPI" else 0.8) * (1.0 + d["score"])
        alpha = 0.75 if d["hop"] == 1 else 0.45
        ax.plot([pos[u][0], pos[v][0]], [pos[u][1], pos[v][1]], color="#B9C2CC",
                 linewidth=lw, alpha=alpha, zorder=1, solid_capstyle="round")

    xs = [pos[n][0] for n in graph.nodes]
    ys = [pos[n][1] for n in graph.nodes]
    ax.set_xlim(min(xs) - 0.18, max(xs) + 0.18)
    ax.set_ylim(min(ys) - 0.18, max(ys) + 0.14)
    ax.set_aspect("equal")
    ax.axis("off")

    # Points-per-data-unit, computed from the REAL axes transform (not a
    # hardcoded constant) so marker/label geometry stays correct at any
    # figsize -- aspect is equal, so a single scalar applies to x and y.
    fig.canvas.draw()
    bbox_px = ax.get_window_extent()
    x0, x1 = ax.get_xlim()
    data_per_point = (x1 - x0) / bbox_px.width * (fig.dpi / 72.0)

    def node_radius_data(kind: str, degree: int, betweenness: float) -> float:
        size = _node_size(kind, degree, betweenness)
        return np.sqrt(size / np.pi) * data_per_point

    # Pathway-membership halo: a thin colored ring behind the node, one
    # color per program the node belongs to (node annotation, not a
    # separate drawn pathway box/edge).
    halo_scale = 1.35
    for node, data in graph.nodes(data=True):
        pathways = data["pathways"]
        if not pathways:
            continue
        r = node_radius_data(data["kind"], data["degree"], data["betweenness"]) * halo_scale
        n = len(pathways)
        for i, pw in enumerate(pathways):
            theta1, theta2 = 360 * i / n, 360 * (i + 1) / n
            wedge = Wedge(pos[node], r, theta1, theta2, width=r * 0.28,
                           facecolor=PATHWAY_COLORS[pw], edgecolor="none", zorder=2)
            ax.add_patch(wedge)

    sizes = [_node_size(d["kind"], d["degree"], d["betweenness"]) for _, d in graph.nodes(data=True)]
    colors = [_node_color(n, d["kind"]) for n, d in graph.nodes(data=True)]
    edgecolors = ["white" for _ in graph.nodes]
    linewidths = [2.2 if d["kind"] == "candidate" else 0.9 for _, d in graph.nodes(data=True)]
    ax.scatter(xs, ys, s=sizes, c=colors, edgecolors=edgecolors, linewidths=linewidths, zorder=3)

    # Candidate labels sit fixed on their (large) node. All other labels
    # are placed just below their node, then resolved against every other
    # label AND every node position by adjustText's real collision-
    # avoidance solver (used here only to move overlapping TEXT -- node
    # positions and edges are never touched).
    texts = []
    for node, data in graph.nodes(data=True):
        if not _should_label(node, data):
            continue
        x, y = pos[node]
        if data["kind"] == "candidate":
            ax.text(x, y, node, fontsize=12.5, fontweight="bold", color="white", ha="center", va="center", zorder=5)
        else:
            fontsize = 9.5 if data["kind"] == "level1_partner" else 8.3
            fontweight = "bold" if node in CANONICAL_BIOLOGY_GENES else "normal"
            r = node_radius_data(data["kind"], data["degree"], data["betweenness"])
            texts.append(ax.text(x, y - r * 1.3 - 0.01, node, fontsize=fontsize, fontweight=fontweight,
                                  color=DGRAY, ha="center", va="top", zorder=5))

    adjust_text(texts, x=xs, y=ys, ax=ax, expand=(1.3, 1.6), force_text=(0.35, 0.6), force_static=(0.3, 0.5),
                arrowprops=dict(arrowstyle="-", color=MGRAY, lw=0.6, alpha=0.8, shrinkA=0, shrinkB=1))

    fig.text(0.02, 0.975, "Molecular interaction networks reveal candidate-specific mechanistic neighborhoods",
              fontsize=19.5, fontweight="bold", color=DGRAY, ha="left", va="top")
    fig.text(0.02, 0.952,
              "STRING interaction_partners, required score ≥ 0.7, species Homo sapiens — post-freeze exploratory analysis, same rule for all four candidates.",
              fontsize=10.5, color="#555555", ha="left", va="top")

    legend_handles = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor=FOCUS_COLORS[c], markeredgecolor="white",
               markersize=13, label=c) for c in CANDIDATES
    ]
    legend_handles += [
        Line2D([0], [0], marker="o", color="none", markerfacecolor=PARTNER_COLOR, markeredgecolor="white",
               markersize=9, label="Level-1 direct partner"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor=BRIDGE_COLOR, markeredgecolor="white",
               markersize=6, label="Level-2 bridge node"),
        Line2D([0], [0], color="#B9C2CC", linewidth=1.7, label="physical PPI"),
        Line2D([0], [0], color="#B9C2CC", linewidth=0.8, alpha=0.6, label="functional association"),
    ]
    leg1 = ax.legend(handles=legend_handles, loc="upper right", frameon=False, fontsize=9.3,
                      labelcolor=DGRAY, bbox_to_anchor=(0.995, 0.98))
    ax.add_artist(leg1)

    pathway_handles = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor=color, markeredgecolor="none", markersize=9, label=label)
        for label, color in PATHWAY_COLORS.items()
    ]
    ax.legend(handles=pathway_handles, loc="lower right", frameon=False, fontsize=9.3, labelcolor=DGRAY,
              title="node halo = program membership", title_fontsize=9.3, bbox_to_anchor=(0.995, 0.02))

    stub.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(stub.with_suffix(".png"), facecolor="white", bbox_inches="tight", dpi=300)
    fig.savefig(stub.with_suffix(".pdf"), facecolor="white", bbox_inches="tight")
    fig.savefig(stub.with_suffix(".svg"), facecolor="white", bbox_inches="tight")
    plt.close(fig)
    logger.info("wrote %s.png/.pdf/.svg (%d nodes, %d edges)", stub, graph.number_of_nodes(), graph.number_of_edges())
    return graph


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    build_network_mechanism_main(OUT_DIR / "NETWORK_mechanism_v2")
