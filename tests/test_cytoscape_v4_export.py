"""Test for the Cytoscape export of the v2/v4 exploratory STRING graph --
verifies table fidelity against the exact source graph, not merely that
the export runs."""

from __future__ import annotations

import hashlib
from pathlib import Path

import networkx as nx
import pandas as pd

from src import cytoscape_v4_export as cx
from src.poster_network_mechanism_v2 import build_network

SHORTLIST = Path("results/tables/evidence_freeze/THERAPEUTIC_SHORTLIST_FREEZE.tsv")
FROZEN_SHA256 = "b6990d7e20f1191df7ebbca18be7542cc312e6e8489f05af6e4b4d31e66fb9dc"


def test_no_string_requery_in_export_module():
    src_text = Path("src/cytoscape_v4_export.py").read_text()
    assert "requests" not in src_text
    assert "string-db.org" not in src_text


def test_graph_invariants_hold():
    cx.verify_graph_invariants(build_network())


def test_edge_table_matches_graph_exactly_no_duplicates():
    graph = build_network()
    edges = cx.build_edge_table(graph)
    assert len(edges) == 147
    table_pairs = {frozenset((r.source, r.target)) for r in edges.itertuples()}
    graph_pairs = {frozenset(e) for e in graph.edges()}
    assert table_pairs == graph_pairs
    assert (edges["source"] < edges["target"]).all()  # alphabetized => no A-B/B-A duplication
    for row in edges.itertuples():
        d = graph.get_edge_data(row.source, row.target)
        assert d["score"] == row.string_score
        assert d["interaction_type"] == row.interaction_type


def test_node_table_matches_graph_exactly():
    graph = build_network()
    nodes = cx.build_node_table(graph)
    assert len(nodes) == 47
    assert set(nodes["gene"]) == set(graph.nodes)
    assert nodes["component_id"].nunique() == 3
    indexed = nodes.set_index("gene")
    for gene in graph.nodes:
        assert indexed.loc[gene, "degree"] == graph.degree(gene)
        assert indexed.loc[gene, "node_kind"] == graph.nodes[gene]["kind"]
    assert indexed.loc["VEZF1", "degree"] == 0
    # KDM1A/USP34 share a component; TLK2 and VEZF1 each sit elsewhere
    assert indexed.loc["KDM1A", "component_id"] == indexed.loc["USP34", "component_id"]
    assert indexed.loc["TLK2", "component_id"] != indexed.loc["KDM1A", "component_id"]
    assert indexed.loc["VEZF1", "component_id"] != indexed.loc["KDM1A", "component_id"]


def test_shortest_path_table_reports_all_four_kdm1a_usp34_paths():
    graph = build_network()
    paths = cx.build_shortest_path_table(graph)
    ku = paths[(paths["candidate_1"] == "KDM1A") & (paths["candidate_2"] == "USP34")]
    assert len(ku) == 4
    assert (ku["path_length"] == 3).all()
    expected = {
        "KDM1A -- DNMT1 -- RPS27A -- USP34",
        "KDM1A -- DNMT1 -- UBA52 -- USP34",
        "KDM1A -- DNMT1 -- UBB -- USP34",
        "KDM1A -- DNMT1 -- UBC -- USP34",
    }
    assert set(ku["path"]) == expected
    no_path = paths[paths["path"] == "NO_PATH"]
    assert len(no_path) == 5  # every pair involving TLK2 or VEZF1


def test_exported_files_exist_and_match_fresh_rebuild():
    cx.export_all()
    for path in (cx.EDGES_TSV, cx.NODES_TSV, cx.PATHS_TSV):
        assert path.exists()
    on_disk = pd.read_csv(cx.EDGES_TSV, sep="\t")
    fresh = cx.build_edge_table(build_network())
    pd.testing.assert_frame_equal(on_disk, fresh, check_dtype=False)


def test_frozen_shortlist_sha256_unchanged():
    digest = hashlib.sha256(SHORTLIST.read_bytes()).hexdigest()
    assert digest == FROZEN_SHA256
