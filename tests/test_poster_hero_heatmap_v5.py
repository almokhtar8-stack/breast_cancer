"""Small dedicated test for the v5 hero heatmap -- confirms v5 is a pure
presentation change over v4 (identical rows/genes/z-scores), plus output
existence. No new testing framework."""

from __future__ import annotations

import hashlib
import re
from pathlib import Path

from PIL import Image

from src import poster_hero_heatmap_v4 as v4
from src import poster_hero_heatmap_v5 as v5

FIGURES = Path("results/figures/poster_hero_heatmap_v5")
SHORTLIST = Path("results/tables/evidence_freeze/THERAPEUTIC_SHORTLIST_FREEZE.tsv")
FROZEN_SHA256 = "b6990d7e20f1191df7ebbca18be7542cc312e6e8489f05af6e4b4d31e66fb9dc"
FOCUS_FOUR = ["KDM1A", "TLK2", "USP34", "VEZF1"]


def test_same_29_biological_rows_as_v4():
    v5_rows = [r for _, rows in [(n, b()) for n, b in v5.DATASET_BUILDERS] for r in rows]
    v4_rows = [r for _, rows in [(n, b()) for n, b in v4.DATASET_BUILDERS] for r in rows]
    assert len(v5_rows) == len(v4_rows) == 29


def test_exact_four_genes():
    assert v5.FOCUS_FOUR == v4.FOCUS_FOUR == FOCUS_FOUR


def test_gse240112_remains_unpaired():
    rows = dict(v5.DATASET_BUILDERS)["GSE240112"]()
    assert all(r.pair_anchor is None for r in rows)


def test_gse245601_has_exactly_three_genuine_matched_pairs():
    rows = dict(v5.DATASET_BUILDERS)["GSE245601"]()
    paired = [r for r in rows if r.pair_anchor is not None]
    assert len(paired) == 3
    for r in paired:
        assert r.pair_anchor.split()[0] == r.label.split()[0]  # same patient


def test_gse111151_group_relationships_preserved():
    rows = dict(v5.DATASET_BUILDERS)["GSE111151"]()
    derivatives = [r for r in rows if r.pair_anchor is not None]
    assert len(derivatives) == 7
    parentals = {r.label for r in rows if r.pair_anchor is None}
    assert parentals == {"MCF-7", "T-47D", "ZR-75-1", "BT-474"}


def test_zscore_matrix_identical_to_v4():
    """v5 must not recompute, re-derive, or perturb a single z-score --
    it reuses v4's DATASET_BUILDERS directly."""
    v5_blocks = dict((n, b()) for n, b in v5.DATASET_BUILDERS)
    v4_blocks = dict((n, b()) for n, b in v4.DATASET_BUILDERS)
    for dataset in v4_blocks:
        v5_by_label = {r.label: r.values for r in v5_blocks[dataset]}
        v4_by_label = {r.label: r.values for r in v4_blocks[dataset]}
        assert v5_by_label.keys() == v4_by_label.keys()
        for label, values in v4_by_label.items():
            assert v5_by_label[label] == values


def test_no_numeric_cell_annotation_layer_in_source():
    src_text = Path("src/poster_hero_heatmap_v5.py").read_text()
    assert "ax.annotate(" not in src_text
    numeric_cell_text = re.search(r"ax\.text\([^)]*\{.*:.2f\}", src_text)
    assert numeric_cell_text is None


def test_output_files_exist_and_are_non_degenerate():
    for ext in ("png", "pdf", "svg"):
        path = FIGURES / f"HERO_sample_level_heatmap_v5.{ext}"
        assert path.exists(), f"{path} was not written"
        assert path.stat().st_size > 0


def test_png_is_landscape_and_high_resolution():
    with Image.open(FIGURES / "HERO_sample_level_heatmap_v5.png") as image:
        width, height = image.size
        assert width > height, "poster hero figure should be landscape"
        assert width >= 2400


def test_note_exists_and_documents_no_clipping_rationale():
    note = Path("results/reports/poster_hero_heatmap_v5/NOTE.md")
    assert note.exists()
    text = note.read_text().lower()
    assert "clip" in text
    assert "identical to v4" in text


def test_frozen_shortlist_sha256_unchanged():
    digest = hashlib.sha256(SHORTLIST.read_bytes()).hexdigest()
    assert digest == FROZEN_SHA256
