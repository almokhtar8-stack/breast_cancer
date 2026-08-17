"""Small dedicated test for the CRISPR discovery figure -- data fidelity
and output existence only, no new testing framework."""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np
from PIL import Image

from src import post_audit_sensitivity_data as pad
from src import poster_story_v1_data as sv1

FIGURES = Path("results/figures/poster_crispr_discovery_v1")
FOCUS_FOUR = ["KDM1A", "TLK2", "USP34", "VEZF1"]


def test_expected_number_of_significant_hits():
    hits = pad.load_significant_sensitising_hits()
    assert len(hits) == 13


def test_ordering_is_by_effect_size_ascending():
    """Rank 1 must be the most negative (strongest sensitising) effect,
    and the plotted order must follow rank_by_effect exactly."""
    hits = pad.load_significant_sensitising_hits()
    assert list(hits["rank_by_effect"]) == list(range(1, len(hits) + 1))
    assert (np.diff(hits["effect_size"].to_numpy()) >= 0).all()
    assert hits.iloc[0]["gene"] == "KDM1A"  # strongest hit in the frozen table


def test_four_focus_genes_present_among_hits():
    hits = pad.load_significant_sensitising_hits()
    present = set(hits["gene"]) & set(FOCUS_FOUR)
    assert present == set(FOCUS_FOUR)
    assert sv1.FOCUS_FOUR == FOCUS_FOUR


def test_no_hand_typed_scientific_values_in_rendering_logic():
    src_text = Path("src/poster_crispr_discovery_v1.py").read_text()
    assert "19103" not in src_text
    assert "19,103" not in src_text
    # no bare literal "13" used as a hit count (the real count must come
    # from len(hits)/n_hits, not a typed constant)
    assert re.search(r"\bn_hits\s*=\s*13\b", src_text) is None
    assert re.search(r"\brange\(1,\s*14\)", src_text) is None


def test_genomewide_and_hit_counts_are_computed_not_typed():
    genomewide = pad.load_genomewide_crispr()
    hits = pad.load_significant_sensitising_hits()
    assert len(genomewide) == 19103
    assert len(hits) == 13


def test_figure_generation_succeeds_and_outputs_exist():
    from src.poster_crispr_discovery_v1 import build_crispr_discovery_main
    stub = FIGURES / "CRISPR_discovery_main"
    build_crispr_discovery_main(stub)
    for ext in ("png", "pdf", "svg"):
        path = stub.with_suffix(f".{ext}")
        assert path.exists()
        assert path.stat().st_size > 0


def test_png_has_a_sane_minimum_resolution():
    with Image.open(FIGURES / "CRISPR_discovery_main.png") as image:
        width, height = image.size
        assert width >= 1800
        assert width > height  # landscape


def test_note_exists_and_documents_source():
    note = Path("results/reports/poster_crispr_discovery_v1/NOTE.md")
    assert note.exists()
    text = note.read_text().lower()
    assert "load_significant_sensitising_hits" in text
    assert "fdr < 0.10" in text or "fdr<0.10" in text or "0.10" in text
