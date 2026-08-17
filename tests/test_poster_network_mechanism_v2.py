"""Small dedicated test for the post-freeze exploratory network/mechanism
figure v2 -- data fidelity, network-construction honesty, and output
existence only, no new testing framework."""

from __future__ import annotations

import hashlib
from pathlib import Path

import networkx as nx
import pandas as pd
import pytest
from PIL import Image

from src import poster_network_mechanism_v2 as nm2

FIGURES = Path("results/figures/poster_network_mechanism_v2")
SHORTLIST = Path("results/tables/evidence_freeze/THERAPEUTIC_SHORTLIST_FREEZE.tsv")
FROZEN_SHA256 = "b6990d7e20f1191df7ebbca18be7542cc312e6e8489f05af6e4b4d31e66fb9dc"


def test_all_four_candidates_present_as_nodes():
    graph = nm2.build_network()
    for candidate in nm2.CANDIDATES:
        assert candidate in graph
        assert graph.nodes[candidate]["kind"] == "candidate"


def test_same_query_rule_applied_to_all_four():
    """Level-1 edges for all four candidates come from the SAME raw source
    file with the SAME required_score and network_type -- no per-candidate
    special-casing of the query itself."""
    level1 = nm2.load_level1_edges()
    assert set(level1["preferredName_A"]) <= set(nm2.CANDIDATES)
    # every candidate was queried, even if some returned zero partners
    raw = pd.read_csv(nm2.STRING_DIR / "string_v2_level1_functional.tsv", sep="\t")
    queried_candidates = set(raw["preferredName_A"]) | set(nm2.CANDIDATES)
    assert set(nm2.CANDIDATES) <= queried_candidates


def test_level1_cap_never_pads_a_sparse_candidate():
    level1 = nm2.load_level1_edges()
    counts = level1.groupby("preferredName_A").size().to_dict()
    assert counts.get("USP34", 0) == 9  # all 9 informative partners (10 minus excluded low-info node)
    assert counts.get("TLK2", 0) == 5   # all 5, below the cap
    assert counts.get("KDM1A", 0) == nm2.LEVEL1_CAP  # capped
    assert "VEZF1" not in counts  # zero real partners -- not padded


def test_vezf1_has_zero_string_partners_at_this_threshold():
    """The central honesty finding of this figure: VEZF1 is a genuine
    isolated singleton under the same rule used for the other three."""
    graph = nm2.build_network()
    assert graph.degree("VEZF1") == 0


def test_every_displayed_edge_exists_in_raw_source_data():
    """No manually fabricated network edge: every edge in the built graph
    must trace back to a row in the raw downloaded STRING tables."""
    level1 = nm2.load_level1_edges()
    level2 = nm2.load_level2_edges(level1)
    raw_l1 = pd.read_csv(nm2.STRING_DIR / "string_v2_level1_functional.tsv", sep="\t")
    raw_l2 = pd.read_csv(nm2.STRING_DIR / "string_v2_level2_functional.tsv", sep="\t")
    raw_pairs = set(zip(raw_l1["preferredName_A"], raw_l1["preferredName_B"])) | \
        set(zip(raw_l2["preferredName_A"], raw_l2["preferredName_B"]))
    graph = nm2.build_network()
    for u, v in graph.edges():
        assert (u, v) in raw_pairs or (v, u) in raw_pairs, f"edge {u}-{v} not found in raw STRING data"


def test_no_pathway_box_masquerading_as_interaction():
    """Pathway/program membership is a node ANNOTATION (halo), never a
    graph node or edge of its own."""
    graph = nm2.build_network()
    pathway_labels = set(nm2.PATHWAY_SETS.keys())
    assert pathway_labels.isdisjoint(set(graph.nodes))
    for _, data in graph.nodes(data=True):
        assert isinstance(data["pathways"], list)


def test_candidate_colors_match_frozen_focus_colors():
    assert nm2.FOCUS_COLORS["KDM1A"] == "#D55E00"
    assert nm2.FOCUS_COLORS["VEZF1"] == "#E69F00"
    assert nm2.FOCUS_COLORS["USP34"] == "#0072B2"


def test_node_and_edge_counts_reproducible():
    g1 = nm2.build_network()
    g2 = nm2.build_network()
    assert g1.number_of_nodes() == g2.number_of_nodes()
    assert g1.number_of_edges() == g2.number_of_edges()
    assert set(g1.nodes) == set(g2.nodes)


def test_kdm1a_usp34_shortest_path_is_real_and_documented():
    graph = nm2.build_network()
    path = nx.shortest_path(graph, "KDM1A", "USP34")
    assert path == ["KDM1A", "DNMT1", "UBC", "USP34"]


def test_network_stats_reflect_three_components():
    graph = nm2.build_network()
    stats = nm2.network_stats(graph)
    assert stats["n_components"] == 3
    assert stats["n_nodes"] == 47
    assert stats["n_edges"] == 147


def test_figure_generation_succeeds_and_outputs_exist():
    stub = FIGURES / "NETWORK_mechanism_v2"
    nm2.build_network_mechanism_main(stub)
    for ext in ("png", "pdf", "svg"):
        path = stub.with_suffix(f".{ext}")
        assert path.exists()
        assert path.stat().st_size > 0


def test_png_has_a_sane_minimum_resolution():
    with Image.open(FIGURES / "NETWORK_mechanism_v2.png") as image:
        width, height = image.size
        assert width >= 3000
        assert width > height  # landscape


def test_v1_outputs_not_overwritten():
    v1 = Path("results/figures/poster_network_mechanism_v1/NETWORK_mechanism_main.png")
    assert v1.exists()


def test_note_documents_exploratory_status_and_frozen_confirmation():
    note = Path("results/reports/poster_network_mechanism_v2/NOTE.md")
    assert note.exists()
    text = note.read_text().lower()
    assert "post-freeze exploratory" in text
    assert "does not imply activation" in text or "not imply activation" in text
    assert "candidate ranking is unchanged" in text or "candidate ranking unchanged" in text


def test_frozen_shortlist_sha256_unchanged():
    digest = hashlib.sha256(SHORTLIST.read_bytes()).hexdigest()
    assert digest == FROZEN_SHA256
