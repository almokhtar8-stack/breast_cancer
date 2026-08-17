"""Small dedicated test for the v2 hero heatmap -- data fidelity and
output existence only, no new testing framework."""

from __future__ import annotations

from pathlib import Path

import numpy as np
from PIL import Image

from src import poster_story_v1_data as sv1

FIGURES = Path("results/figures/poster_hero_heatmap_v2")
FOCUS_FOUR = ["KDM1A", "TLK2", "USP34", "VEZF1"]


def test_hero_pairs_unchanged_source_used():
    """v2 reuses the exact same frozen pairs builder as v1 -- no new
    science introduced for this figure."""
    pairs = sv1.build_hero_heatmap_pairs()
    assert len(pairs) == 16
    assert set(pairs["gene"]) == set(FOCUS_FOUR)
    assert set(pairs["dataset"]) == set(sv1.DATASET_ORDER)


def test_tlk2_is_negative_in_all_four_contexts():
    """Locks in the real, already-established finding the design brief
    explicitly required to remain visible (TLK2 consistently down) --
    this test fails if the biology is ever accidentally changed, not
    because a threshold was chosen to make it pass."""
    pairs = sv1.build_hero_heatmap_pairs()
    tlk2 = pairs[pairs["gene"] == "TLK2"]
    assert (tlk2["log2fc"] < 0).all()


def test_acute_context_negative_across_all_four_genes():
    pairs = sv1.build_hero_heatmap_pairs()
    acute = pairs[pairs["dataset"] == "GSE245601"]
    assert (acute["log2fc"] < 0).all()


def test_output_files_exist_and_are_non_degenerate():
    for ext in ("png", "pdf", "svg"):
        path = FIGURES / f"HERO_main_heatmap_v2.{ext}"
        assert path.exists(), f"{path} was not written"
        assert path.stat().st_size > 0


def test_png_has_a_sane_minimum_resolution():
    with Image.open(FIGURES / "HERO_main_heatmap_v2.png") as image:
        width, height = image.size
        assert width >= 1200
        assert height >= 800


def test_note_exists():
    note = Path("results/reports/poster_hero_heatmap_v2/NOTE.md")
    assert note.exists()
    text = note.read_text().lower()
    assert "reference row" in text or "reference/row" in text
    assert "comparison row" in text
