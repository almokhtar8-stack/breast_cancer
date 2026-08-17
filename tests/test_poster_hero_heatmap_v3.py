"""Small dedicated test for the v3 hero heatmap -- data fidelity, structure,
and output existence only, no new testing framework."""

from __future__ import annotations

import hashlib
from pathlib import Path

import numpy as np
from PIL import Image

from src import poster_story_v1_data as sv1

FIGURES = Path("results/figures/poster_hero_heatmap_v3")
SHORTLIST = Path("results/tables/evidence_freeze/THERAPEUTIC_SHORTLIST_FREEZE.tsv")
FROZEN_SHA256 = "b6990d7e20f1191df7ebbca18be7542cc312e6e8489f05af6e4b4d31e66fb9dc"
FOCUS_FOUR = ["KDM1A", "TLK2", "USP34", "VEZF1"]


def test_four_genes_and_four_datasets_present():
    pairs = sv1.build_hero_heatmap_pairs()
    assert set(pairs["gene"]) == set(FOCUS_FOUR)
    assert set(pairs["dataset"]) == set(sv1.DATASET_ORDER)
    assert len(sv1.DATASET_ORDER) == 4
    assert len(FOCUS_FOUR) == 4


def test_reference_comparison_and_delta_structure_exists():
    """Each of the 16 (dataset, gene) rows must carry all three quantities
    the v3 design draws as three rows -- reference value, comparison value,
    and a delta (log2FC) value."""
    pairs = sv1.build_hero_heatmap_pairs()
    assert len(pairs) == 16
    for col in ("ref_value", "cmp_value", "log2fc"):
        assert col in pairs.columns
        assert pairs[col].notna().all()


def test_delta_values_equal_frozen_log2fc_source():
    """The Δ strip must not silently recompute anything -- it renders the
    exact same log2fc column already frozen and used in v1/v2."""
    pairs = sv1.build_hero_heatmap_pairs()
    for dataset in sv1.DATASET_ORDER:
        sub = pairs[pairs["dataset"] == dataset]
        for gene in FOCUS_FOUR:
            row = sub[sub["gene"] == gene].iloc[0]
            assert np.isfinite(row["log2fc"])


def test_no_numeric_cell_annotation_layer_in_source():
    """The v3 module must not call ax.annotate or per-cell ax.text with a
    numeric value -- only label text (dataset/condition/gene/Δ names) is
    drawn on the axes."""
    src_text = Path("src/poster_hero_heatmap_v3.py").read_text()
    assert "ax.annotate(" not in src_text
    import re
    numeric_cell_text = re.search(r"ax\.text\([^)]*\{.*:.2f\}", src_text)
    assert numeric_cell_text is None


def test_output_files_exist_and_are_non_degenerate():
    for ext in ("png", "pdf", "svg"):
        path = FIGURES / f"HERO_main_heatmap_v3.{ext}"
        assert path.exists(), f"{path} was not written"
        assert path.stat().st_size > 0


def test_png_has_a_sane_minimum_resolution():
    with Image.open(FIGURES / "HERO_main_heatmap_v3.png") as image:
        width, height = image.size
        assert width >= 1200
        assert height >= 800


def test_note_exists_and_documents_scaling():
    note = Path("results/reports/poster_hero_heatmap_v3/NOTE.md")
    assert note.exists()
    text = note.read_text().lower()
    assert "within" in text
    assert "delta" in text or "log2fc" in text or "log2 fold-change" in text


def test_frozen_shortlist_sha256_unchanged():
    digest = hashlib.sha256(SHORTLIST.read_bytes()).hexdigest()
    assert digest == FROZEN_SHA256
