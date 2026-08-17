"""Small dedicated test for the pathway figure v2 -- data fidelity,
theme-rule selection, and output existence only."""

from __future__ import annotations

import hashlib
from pathlib import Path

import pytest
from PIL import Image

from src import poster_exploration_v2_data as ed
from src import poster_pathway_v2 as pw2

FIGURES = Path("results/figures/poster_pathway_v2")
SHORTLIST = Path("results/tables/evidence_freeze/THERAPEUTIC_SHORTLIST_FREEZE.tsv")
FROZEN_SHA256 = "b6990d7e20f1191df7ebbca18be7542cc312e6e8489f05af6e4b4d31e66fb9dc"


def test_exact_four_contexts():
    assert pw2.DATASET_ORDER == ["gse118713", "gse111151", "gse240112", "gse245601"]
    assert set(pw2.CONTEXT_LABEL) == set(pw2.DATASET_ORDER)


def test_all_values_loaded_from_frozen_source_tables():
    matrix = pw2.load_matrix()
    long = ed.load_pathway_trajectories(pw2.PATHWAYS)
    assert len(long) == len(pw2.PATHWAYS) * len(pw2.DATASET_ORDER)
    for _, row in long.iterrows():
        nes, fdr = matrix[row["pathway_label"]][row["dataset"]]
        assert nes == pytest.approx(row["NES"])
        assert fdr == pytest.approx(row["fdr"])
        assert not (pytest.approx(nes) == 0)  # every value real, none defaulted


def test_pathway_selection_follows_predefined_network_theme_rule():
    """Every displayed pathway must map to a pre-specified network theme,
    and every theme entry must be a real Hallmark set present in all four
    frozen tables."""
    themed = [p for theme in pw2.NETWORK_THEME_PATHWAYS.values() for p in theme]
    assert pw2.PATHWAYS == themed
    expected_names = {
        "HALLMARK_ESTROGEN_RESPONSE_EARLY",
        "HALLMARK_ESTROGEN_RESPONSE_LATE",
        "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
        "HALLMARK_WNT_BETA_CATENIN_SIGNALING",
        "HALLMARK_E2F_TARGETS",
    }
    assert {name for _, name, _ in pw2.PATHWAYS} == expected_names
    long = ed.load_pathway_trajectories(pw2.PATHWAYS)
    for _, name, label in pw2.PATHWAYS:
        sub = long[long["pathway_label"] == label]
        assert len(sub) == 4
        assert sub["NES"].notna().all() and sub["fdr"].notna().all()


def test_no_hand_typed_nes_or_fdr_in_rendering_logic():
    src_text = Path("src/poster_pathway_v2.py").read_text()
    for literal in ("-2.88", "1.978", "-1.559", "2.166", "1.767", "6.66e-27", "0.168"):
        assert literal not in src_text


def test_context_labels_never_mislabel_recurrence_or_acute_as_resistance():
    assert "resistan" not in pw2.CONTEXT_LABEL["gse240112"].lower()
    assert "recurrent" in pw2.CONTEXT_LABEL["gse240112"].lower()
    assert "resistan" not in pw2.CONTEXT_LABEL["gse245601"].lower()
    assert "acute" in pw2.CONTEXT_LABEL["gse245601"].lower()


def test_emt_direction_derived_correctly_from_source():
    matrix = pw2.load_matrix()
    for ds in ("gse118713", "gse111151", "gse240112"):
        assert matrix["EMT"][ds][0] > 0
    assert matrix["EMT"]["gse245601"][0] < 0


def test_estrogen_direction_derived_correctly_from_source():
    matrix = pw2.load_matrix()
    for label in ("Estrogen response — early", "Estrogen response — late"):
        for ds in pw2.DATASET_ORDER:
            nes, fdr = matrix[label][ds]
            assert nes < 0
            assert fdr < pw2.FDR_THRESHOLD


def test_figure_generation_succeeds_and_outputs_exist():
    stub = FIGURES / "PATHWAY_v2"
    pw2.build_pathway_v2(stub)
    for ext in ("png", "pdf", "svg"):
        path = stub.with_suffix(f".{ext}")
        assert path.exists()
        assert path.stat().st_size > 0


def test_png_has_a_sane_minimum_resolution():
    with Image.open(FIGURES / "PATHWAY_v2.png") as image:
        width, height = image.size
        assert width >= 2000
        assert width > height


def test_pathway_v1_not_overwritten():
    assert Path("results/figures/poster_pathway_v1/PATHWAY_main.png").exists()


def test_note_documents_caveats_and_theme_gaps():
    note = Path("results/reports/poster_pathway_v2/NOTE.md")
    assert note.exists()
    text = note.read_text().lower()
    assert "recurrence-associated" in text
    assert "acute 12 h" in text
    assert "chromatin regulation" in text  # theme-gap documented
    assert "g2m" in text  # companion-set decision documented


def test_frozen_shortlist_sha256_unchanged():
    digest = hashlib.sha256(SHORTLIST.read_bytes()).hexdigest()
    assert digest == FROZEN_SHA256
