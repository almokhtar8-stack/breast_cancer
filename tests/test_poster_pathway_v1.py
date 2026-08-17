"""Small dedicated test for the pathway biology figure -- data fidelity
and output existence only, no new testing framework."""

from __future__ import annotations

import hashlib
import re
from pathlib import Path

import pytest
from PIL import Image

from src import poster_exploration_v2_data as ed
from src import poster_pathway_v1 as pw

FIGURES = Path("results/figures/poster_pathway_v1")
SHORTLIST = Path("results/tables/evidence_freeze/THERAPEUTIC_SHORTLIST_FREEZE.tsv")
FROZEN_SHA256 = "b6990d7e20f1191df7ebbca18be7542cc312e6e8489f05af6e4b4d31e66fb9dc"


def test_exact_four_contexts():
    assert pw.DATASET_ORDER == ["gse118713", "gse111151", "gse240112", "gse245601"]
    assert set(pw.CONTEXT_LABEL) == set(pw.DATASET_ORDER)
    assert set(pw.CONTEXT_COLORS) == set(pw.DATASET_ORDER)


def test_frozen_nes_values_match_source_tables():
    """Locks in the real frozen NES values the figure depends on --
    fails if the frozen GSEA tables are ever altered."""
    hits = ed.load_pathway_trajectories(pw.PATHWAYS_A + [pw.PATHWAY_B])
    row = hits[(hits["pathway_label"] == "Estrogen Response -- Early") & (hits["dataset"] == "gse118713")].iloc[0]
    assert row["NES"] == pytest.approx(-2.881807, abs=1e-4)
    row = hits[(hits["pathway_label"] == "EMT") & (hits["dataset"] == "gse245601")].iloc[0]
    assert row["NES"] == pytest.approx(-1.559312, abs=1e-4)
    row = hits[(hits["pathway_label"] == "EMT") & (hits["dataset"] == "gse118713")].iloc[0]
    assert row["NES"] == pytest.approx(1.977988, abs=1e-4)


def test_pathway_names_are_real_hallmark_sets_present_in_source():
    all_pathways = pw.PATHWAYS_A + [pw.PATHWAY_B] + pw.PATHWAYS_C
    long = ed.load_pathway_trajectories(all_pathways)
    for _, _, label in all_pathways:
        sub = long[long["pathway_label"] == label]
        assert len(sub) == len(pw.DATASET_ORDER)
        assert sub["NES"].notna().all()


def test_emt_direction_matches_the_resistance_vs_acute_story():
    """The central biological claim this figure makes: EMT is positive in
    every resistance/recurrence context and negative in the acute context."""
    long = ed.load_pathway_trajectories([pw.PATHWAY_B])
    wide = long.set_index("dataset")["NES"]
    for ds in ("gse118713", "gse111151", "gse240112"):
        assert wide[ds] > 0
    assert wide["gse245601"] < 0


def test_acute_context_labeled_correctly_and_gse240112_never_called_resistance():
    assert "acute" in pw.CONTEXT_LABEL["gse245601"].lower()
    assert "resistan" not in pw.CONTEXT_LABEL["gse245601"].lower()
    assert "resistan" not in pw.CONTEXT_LABEL["gse240112"].lower()
    assert "recurrent" in pw.CONTEXT_LABEL["gse240112"].lower()


def test_no_hand_typed_nes_values_in_rendering_logic():
    src_text = Path("src/poster_pathway_v1.py").read_text()
    # spot-check a couple of the real frozen NES values never appear as
    # literals in the render source
    assert "-2.88" not in src_text
    assert "1.977" not in src_text
    assert "-1.559" not in src_text


def test_figure_generation_succeeds_and_outputs_exist():
    stub = FIGURES / "PATHWAY_main"
    pw.build_pathway_main(stub)
    for ext in ("png", "pdf", "svg"):
        path = stub.with_suffix(f".{ext}")
        assert path.exists()
        assert path.stat().st_size > 0


def test_png_has_a_sane_minimum_resolution():
    with Image.open(FIGURES / "PATHWAY_main.png") as image:
        width, height = image.size
        assert width >= 2000


def test_note_exists_and_documents_caveats():
    note = Path("results/reports/poster_pathway_v1/NOTE.md")
    assert note.exists()
    text = note.read_text().lower()
    assert "recurrence-associated" in text
    assert "acute 12h ex vivo" in text or "acute 12 h" in text


def test_frozen_shortlist_sha256_unchanged():
    digest = hashlib.sha256(SHORTLIST.read_bytes()).hexdigest()
    assert digest == FROZEN_SHA256
