"""Small dedicated test for the v3 visual-refinement-only network figure --
confirms the graph is IDENTICAL to v2 (same nodes, same edges, same
annotations) and only rendering changed. No new testing framework."""

from __future__ import annotations

import hashlib
from pathlib import Path

import networkx as nx
from PIL import Image

from src import poster_network_mechanism_v2 as nm2
from src import poster_network_mechanism_v3 as nm3

FIGURES = Path("results/figures/poster_network_mechanism_v3")
SHORTLIST = Path("results/tables/evidence_freeze/THERAPEUTIC_SHORTLIST_FREEZE.tsv")
FROZEN_SHA256 = "b6990d7e20f1191df7ebbca18be7542cc312e6e8489f05af6e4b4d31e66fb9dc"


def test_v3_reuses_v2_build_network_function_unmodified():
    """v3 must call v2's own build_network, not a re-implementation."""
    assert nm3.build_network is nm2.build_network
    assert nm3.network_stats is nm2.network_stats


def test_node_count_edge_count_and_components_match_v2():
    g2 = nm2.build_network()
    g3 = nm3.build_network()
    assert g3.number_of_nodes() == g2.number_of_nodes() == 47
    assert g3.number_of_edges() == g2.number_of_edges() == 147
    assert nx.number_connected_components(g3) == nx.number_connected_components(g2) == 3


def test_node_identity_matches_v2_exactly():
    g2, g3 = nm2.build_network(), nm3.build_network()
    assert set(g2.nodes) == set(g3.nodes)


def test_edge_identity_matches_v2_exactly():
    g2, g3 = nm2.build_network(), nm3.build_network()
    edges2 = {frozenset(e) for e in g2.edges()}
    edges3 = {frozenset(e) for e in g3.edges()}
    assert edges2 == edges3


def test_pathway_annotations_match_v2_exactly():
    g2, g3 = nm2.build_network(), nm3.build_network()
    for node in g2.nodes:
        assert g2.nodes[node]["pathways"] == g3.nodes[node]["pathways"]


def test_edge_scores_and_interaction_types_match_v2_exactly():
    g2, g3 = nm2.build_network(), nm3.build_network()
    for u, v, d2 in g2.edges(data=True):
        d3 = g3.get_edge_data(u, v)
        assert d3["score"] == d2["score"]
        assert d3["interaction_type"] == d2["interaction_type"]


def test_all_four_candidates_present():
    graph = nm3.build_network()
    for candidate in nm3.CANDIDATES:
        assert candidate in graph
        assert graph.nodes[candidate]["kind"] == "candidate"


def test_vezf1_degree_remains_zero():
    graph = nm3.build_network()
    assert graph.degree("VEZF1") == 0


def test_tlk2_outside_kdm1a_usp34_component():
    graph = nm3.build_network()
    components = list(nx.connected_components(graph))
    main = next(c for c in components if "KDM1A" in c and "USP34" in c)
    assert "TLK2" not in main
    tlk2_component = next(c for c in components if "TLK2" in c)
    assert "KDM1A" not in tlk2_component and "USP34" not in tlk2_component


def test_kdm1a_and_usp34_share_a_component():
    graph = nm3.build_network()
    assert nx.has_path(graph, "KDM1A", "USP34")


def test_figure_generation_succeeds_and_outputs_exist():
    stub = FIGURES / "NETWORK_mechanism_v3"
    nm3.build_network_mechanism_main(stub)
    for ext in ("png", "pdf", "svg"):
        path = stub.with_suffix(f".{ext}")
        assert path.exists()
        assert path.stat().st_size > 0


def test_png_has_a_sane_minimum_resolution():
    with Image.open(FIGURES / "NETWORK_mechanism_v3.png") as image:
        width, _ = image.size
        assert width >= 3000


def test_v2_outputs_not_overwritten():
    v2 = Path("results/figures/poster_network_mechanism_v2/NETWORK_mechanism_v2.png")
    assert v2.exists()


def test_note_documents_visual_only_status():
    note = Path("results/reports/poster_network_mechanism_v3/NOTE.md")
    assert note.exists()
    text = note.read_text().lower()
    assert "no science changed" in text or "rendering-only" in text or "rendering only" in text
    assert "identical" in text


def test_frozen_shortlist_sha256_unchanged():
    digest = hashlib.sha256(SHORTLIST.read_bytes()).hexdigest()
    assert digest == FROZEN_SHA256
