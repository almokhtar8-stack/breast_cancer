"""Small dedicated test for the network/mechanism figure -- data fidelity,
evidence-honesty guarantees, and output existence only, no new testing
framework."""

from __future__ import annotations

import hashlib
from pathlib import Path

import pytest
from PIL import Image

from src import poster_network_mechanism_v1 as nm

FIGURES = Path("results/figures/poster_network_mechanism_v1")
SHORTLIST = Path("results/tables/evidence_freeze/THERAPEUTIC_SHORTLIST_FREEZE.tsv")
FROZEN_SHA256 = "b6990d7e20f1191df7ebbca18be7542cc312e6e8489f05af6e4b4d31e66fb9dc"


def test_focus_four_and_colors():
    assert set(nm.FOCUS_FOUR) == {"KDM1A", "TLK2", "USP34", "VEZF1"}
    assert set(nm.FOCUS_COLORS) >= set(nm.FOCUS_FOUR)


def test_crispr_stats_match_frozen_table():
    stats = nm._load_crispr_stats()
    assert stats["KDM1A"]["rank"] == 1
    assert stats["TLK2"]["rank"] == 4
    assert stats["VEZF1"]["rank"] == 8
    assert stats["USP34"]["rank"] == 12
    assert stats["KDM1A"]["fdr"] == pytest.approx(0.000385, abs=1e-5)


def test_usp34_network_is_real_frozen_data():
    neighbors, paths = nm._load_usp34_network()
    assert set(neighbors["neighbor_gene"]) >= {"RPS27A", "UBC", "USP9X"}
    assert (neighbors["database_source"] == "STRING").all()
    bridge_targets = set(paths["target_gene"])
    assert {"CTNNB1", "SOX2"} <= bridge_targets


def test_vezf1_has_exactly_one_frozen_edge():
    neighbors = nm._load_vezf1_network()
    assert len(neighbors) == 1
    assert neighbors.iloc[0]["neighbor_gene"] == "DMTN"
    assert neighbors.iloc[0]["database_source"] == "MSigDB_hallmark"


def test_kdm1a_and_tlk2_absent_from_every_frozen_network_table():
    """The central honesty guarantee of this figure: KDM1A/TLK2 must not
    appear in any frozen PPI/network table, confirming no network edge is
    being fabricated for them."""
    direct = __import__("pandas").read_csv(nm.NET_DIR / "four_candidate_direct_neighbors.tsv", sep="\t")
    assert not direct["candidate"].isin(["KDM1A", "TLK2"]).any()
    assert not direct["neighbor_gene"].isin(["KDM1A", "TLK2"]).any()
    membership = __import__("pandas").read_csv(nm.NET_DIR / "candidate_pathway_membership.tsv", sep="\t")
    assert not membership["candidate"].isin(["KDM1A", "TLK2"]).any()


def test_kdm1a_and_tlk2_pathway_membership_is_real_and_leading_edge():
    for gene, pathways in [("KDM1A", nm.PATHWAYS_KDM1A), ("TLK2", nm.PATHWAYS_TLK2)]:
        hits = nm._load_pathway_membership(gene, pathways)
        assert len(hits) == len(pathways)
        assert set(hits["pathway"]) == set(pathways)
        assert hits["NES"].notna().all()
        assert hits["fdr"].notna().all()


def test_kdm1a_and_tlk2_share_the_chromatin_organization_pathway():
    shared = "GOBP_REGULATION_OF_CHROMATIN_ORGANIZATION"
    kdm1a_full = nm._full_leading_edge_membership("KDM1A")
    tlk2_full = nm._full_leading_edge_membership("TLK2")
    assert shared in set(kdm1a_full["pathway"])
    assert shared in set(tlk2_full["pathway"])


def test_no_hand_typed_nes_or_fdr_values_in_rendering_logic():
    src_text = Path("src/poster_network_mechanism_v1.py").read_text()
    # spot-check that specific frozen numeric values never appear as literals
    assert "0.000385" not in src_text
    assert "-2.881807" not in src_text


def test_fmt_fdr_avoids_misleading_zero():
    assert nm._fmt_fdr(0.000385) == "FDR<0.001"
    assert nm._fmt_fdr(0.042) == "FDR=0.042"


def test_figure_generation_succeeds_and_outputs_exist():
    stub = FIGURES / "NETWORK_mechanism_main"
    nm.build_network_mechanism_main(stub)
    for ext in ("png", "pdf", "svg"):
        path = stub.with_suffix(f".{ext}")
        assert path.exists()
        assert path.stat().st_size > 0


def test_png_has_a_sane_minimum_resolution():
    with Image.open(FIGURES / "NETWORK_mechanism_main.png") as image:
        width, height = image.size
        assert width >= 2000
        assert width > height  # landscape


def test_note_exists_and_documents_the_evidence_gap():
    note = Path("results/reports/poster_network_mechanism_v1/NOTE.md")
    assert note.exists()
    text = note.read_text().lower()
    assert "no pp" in text or "no network edge" in text or "no direct" in text
    assert "blind" in text  # documents why KDM1A had no frozen network build


def test_frozen_shortlist_sha256_unchanged():
    digest = hashlib.sha256(SHORTLIST.read_bytes()).hexdigest()
    assert digest == FROZEN_SHA256
