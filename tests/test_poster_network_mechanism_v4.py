"""Small dedicated test for the v4 five-panel network figure -- confirms
v4 is a pure presentation rebuild of the v2 graph: same source data, no
STRING requery, deterministic panel derivation, honest connectivity."""

from __future__ import annotations

import hashlib
from pathlib import Path

import networkx as nx
from PIL import Image

from src import poster_network_mechanism_v2 as nm2
from src import poster_network_mechanism_v4 as nm4

FIGURES = Path("results/figures/poster_network_mechanism_v4")
SHORTLIST = Path("results/tables/evidence_freeze/THERAPEUTIC_SHORTLIST_FREEZE.tsv")
FROZEN_SHA256 = "b6990d7e20f1191df7ebbca18be7542cc312e6e8489f05af6e4b4d31e66fb9dc"


def test_no_string_download_or_api_call_in_v4():
    src_text = Path("src/poster_network_mechanism_v4.py").read_text()
    assert "requests" not in src_text
    assert "string-db.org" not in src_text
    assert "urllib" not in src_text


def test_v4_reuses_exact_v2_build_network():
    assert nm4.build_network is nm2.build_network


def test_every_local_panel_edge_exists_in_v2_graph():
    graph = nm4.build_network()
    v2_edges = {frozenset(e) for e in graph.edges()}
    for candidate in nm4.CANDIDATES:
        sub = nm4.local_subgraph(graph, candidate)
        for u, v in sub.edges():
            assert frozenset((u, v)) in v2_edges, f"{candidate} panel edge {u}-{v} not in v2 graph"


def test_local_panel_node_sets_are_deterministically_derived():
    graph = nm4.build_network()
    for candidate in nm4.CANDIDATES:
        sub1 = nm4.local_subgraph(graph, candidate)
        sub2 = nm4.local_subgraph(nm4.build_network(), candidate)
        assert set(sub1.nodes) == set(sub2.nodes)
        assert {frozenset(e) for e in sub1.edges()} == {frozenset(e) for e in sub2.edges()}
        partners = set(graph.neighbors(candidate))
        for node in sub1.nodes:
            if node == candidate or node in partners:
                continue
            assert graph.nodes[node]["kind"] == "level2_bridge"
            assert any(graph.has_edge(node, p) for p in partners)


def test_expected_panel_sizes_from_current_graph():
    graph = nm4.build_network()
    sizes = {c: nm4.local_subgraph(graph, c).number_of_nodes() for c in nm4.CANDIDATES}
    assert sizes == {"KDM1A": 29, "TLK2": 6, "USP34": 13, "VEZF1": 1}


def test_no_candidate_appears_in_another_candidates_panel():
    graph = nm4.build_network()
    for candidate in nm4.CANDIDATES:
        sub = nm4.local_subgraph(graph, candidate)
        others = set(nm4.CANDIDATES) - {candidate}
        assert others.isdisjoint(set(sub.nodes))


def test_kdm1a_and_usp34_share_a_global_component():
    graph = nm4.build_network()
    assert nx.has_path(graph, "KDM1A", "USP34")


def test_tlk2_in_separate_global_component():
    graph = nm4.build_network()
    assert not nx.has_path(graph, "TLK2", "KDM1A")
    assert not nx.has_path(graph, "TLK2", "USP34")


def test_vezf1_degree_zero():
    graph = nm4.build_network()
    assert graph.degree("VEZF1") == 0
    sub = nm4.local_subgraph(graph, "VEZF1")
    assert sub.number_of_nodes() == 1 and sub.number_of_edges() == 0


def test_shortest_path_computed_programmatically_matches_expected():
    graph = nm4.build_network()
    paths = nm4.candidate_shortest_paths(graph)
    result = paths[("KDM1A", "USP34")]
    assert result is not None
    assert result["n_edges"] == 3
    assert result["path"] == ["KDM1A", "DNMT1", "UBC", "USP34"]
    all_paths = sorted(nx.all_shortest_paths(graph, "KDM1A", "USP34"))
    assert result["path"] in all_paths
    assert result["n_equally_short"] == len(all_paths) == 4
    assert all(p[1] == "DNMT1" for p in all_paths)
    # no fabricated paths for the disconnected pairs
    assert paths[("KDM1A", "TLK2")] is None
    assert paths[("KDM1A", "VEZF1")] is None
    assert paths[("TLK2", "USP34")] is None
    assert paths[("TLK2", "VEZF1")] is None
    assert paths[("USP34", "VEZF1")] is None


def test_no_arrows_rendered_on_edges():
    """All connectors are plain undirected lines: the module must never
    use matplotlib arrow primitives (adjustText's internal leader-line
    arrowprops uses arrowstyle '-', a plain line, and is allowed)."""
    src_text = Path("src/poster_network_mechanism_v4.py").read_text()
    assert "FancyArrow" not in src_text
    assert "ax.annotate" not in src_text
    assert "ax.arrow" not in src_text
    assert 'arrowstyle="->"' not in src_text and "arrowstyle='->'" not in src_text


def test_figure_generation_succeeds_and_outputs_exist():
    stub = FIGURES / "NETWORK_mechanism_v4"
    nm4.build_network_mechanism_main(stub)
    for ext in ("png", "pdf", "svg"):
        path = stub.with_suffix(f".{ext}")
        assert path.exists()
        assert path.stat().st_size > 0


def test_png_has_a_sane_minimum_resolution():
    with Image.open(FIGURES / "NETWORK_mechanism_v4.png") as image:
        width, height = image.size
        assert width >= 3000
        assert width > height  # landscape


def test_earlier_versions_not_overwritten():
    for v in ("v1/NETWORK_mechanism_main", "v2/NETWORK_mechanism_v2", "v3/NETWORK_mechanism_v3"):
        assert Path(f"results/figures/poster_network_mechanism_{v}.png").exists()


def test_note_documents_key_facts():
    note = Path("results/reports/poster_network_mechanism_v4/NOTE.md")
    assert note.exists()
    text = note.read_text().lower()
    assert "no strin" in text.replace("requery", "requery") and "requery" in text
    assert "post-freeze exploratory" in text
    assert "undirected" in text
    assert "isolated" in text


def test_frozen_shortlist_sha256_unchanged():
    digest = hashlib.sha256(SHORTLIST.read_bytes()).hexdigest()
    assert digest == FROZEN_SHA256
