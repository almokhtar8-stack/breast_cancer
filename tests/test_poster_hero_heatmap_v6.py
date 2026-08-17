"""Small dedicated test for the v6 hero heatmap -- confirms v6 is a pure
geometry/layout change over v5 (identical rows/genes/order/z-scores),
plus output existence. No new testing framework."""

from __future__ import annotations

import hashlib
import re
from pathlib import Path

from PIL import Image

from src import poster_hero_heatmap_v4 as v4
from src import poster_hero_heatmap_v6 as v6

FIGURES = Path("results/figures/poster_hero_heatmap_v6")
SHORTLIST = Path("results/tables/evidence_freeze/THERAPEUTIC_SHORTLIST_FREEZE.tsv")
FROZEN_SHA256 = "b6990d7e20f1191df7ebbca18be7542cc312e6e8489f05af6e4b4d31e66fb9dc"
FOCUS_FOUR = ["KDM1A", "TLK2", "USP34", "VEZF1"]


def test_29_biological_rows_unchanged():
    rows = [r for _, rows in [(n, b()) for n, b in v6.DATASET_BUILDERS] for r in rows]
    assert len(rows) == 29


def test_four_genes_unchanged():
    assert v6.FOCUS_FOUR == v4.FOCUS_FOUR == FOCUS_FOUR


def test_same_row_order_as_v4():
    v6_labels = [r.label for _, rows in [(n, b()) for n, b in v6.DATASET_BUILDERS] for r in rows]
    v4_labels = [r.label for _, rows in [(n, b()) for n, b in v4.DATASET_BUILDERS] for r in rows]
    assert v6_labels == v4_labels


def test_zscore_matrix_identical_to_v5():
    """v6 must not recompute, re-derive, reorder, or perturb a single
    z-score -- it reuses v4's DATASET_BUILDERS directly, same as v5."""
    v6_blocks = dict((n, b()) for n, b in v6.DATASET_BUILDERS)
    v4_blocks = dict((n, b()) for n, b in v4.DATASET_BUILDERS)
    for dataset in v4_blocks:
        v6_by_label = {r.label: r.values for r in v6_blocks[dataset]}
        v4_by_label = {r.label: r.values for r in v4_blocks[dataset]}
        assert v6_by_label.keys() == v4_by_label.keys()
        for label, values in v4_by_label.items():
            assert v6_by_label[label] == values


def test_gse245601_pairing_preserved():
    rows = dict(v6.DATASET_BUILDERS)["GSE245601"]()
    paired = [r for r in rows if r.pair_anchor is not None]
    assert len(paired) == 3
    for r in paired:
        assert r.pair_anchor.split()[0] == r.label.split()[0]


def test_gse240112_remains_unpaired():
    rows = dict(v6.DATASET_BUILDERS)["GSE240112"]()
    assert all(r.pair_anchor is None for r in rows)


def test_no_numeric_cell_annotation_layer_in_source():
    src_text = Path("src/poster_hero_heatmap_v6.py").read_text()
    assert "ax.annotate(" not in src_text
    numeric_cell_text = re.search(r"ax\.text\([^)]*\{.*:.2f\}", src_text)
    assert numeric_cell_text is None


def test_output_files_exist_and_are_non_degenerate():
    for ext in ("png", "pdf", "svg"):
        path = FIGURES / f"HERO_sample_level_heatmap_v6.{ext}"
        assert path.exists(), f"{path} was not written"
        assert path.stat().st_size > 0


def test_png_is_landscape_dense_and_high_resolution():
    with Image.open(FIGURES / "HERO_sample_level_heatmap_v6.png") as image:
        width, height = image.size
        assert width > height
        ratio = width / height
        assert 1.2 <= ratio <= 1.9  # roughly 4:3 to slightly wider than 3:2
        assert width >= 2400


def test_frozen_shortlist_sha256_unchanged():
    digest = hashlib.sha256(SHORTLIST.read_bytes()).hexdigest()
    assert digest == FROZEN_SHA256
