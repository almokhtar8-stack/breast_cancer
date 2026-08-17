"""Small dedicated test for the v4 sample-level hero heatmap -- data
fidelity, real-sample structure, and output existence only, no new
testing framework."""

from __future__ import annotations

import hashlib
import re
from pathlib import Path

import numpy as np
from PIL import Image

from src import poster_hero_heatmap_v4 as v4

FIGURES = Path("results/figures/poster_hero_heatmap_v4")
SHORTLIST = Path("results/tables/evidence_freeze/THERAPEUTIC_SHORTLIST_FREEZE.tsv")
FROZEN_SHA256 = "b6990d7e20f1191df7ebbca18be7542cc312e6e8489f05af6e4b4d31e66fb9dc"
FOCUS_FOUR = ["KDM1A", "TLK2", "USP34", "VEZF1"]
DATASETS = ["GSE118713", "GSE111151", "GSE240112", "GSE245601"]

EXPECTED_N = {"GSE118713": 6, "GSE111151": 11, "GSE240112": 6, "GSE245601": 6}


def test_four_genes_and_four_datasets_present():
    assert set(v4.FOCUS_FOUR) == set(FOCUS_FOUR)
    assert {name for name, _ in v4.DATASET_BUILDERS} == set(DATASETS)


def test_biological_row_counts_match_real_sample_n():
    """Every row is one real biological observation -- no group means, no
    invented pseudoreplicates."""
    for dataset, builder in v4.DATASET_BUILDERS:
        rows = builder()
        assert len(rows) == EXPECTED_N[dataset], f"{dataset}: expected {EXPECTED_N[dataset]} rows, got {len(rows)}"
        labels = [r.label for r in rows]
        assert len(labels) == len(set(labels)), f"{dataset}: duplicate row labels imply pseudoreplication"


def test_real_sample_identifiers_used():
    gse118713 = dict(v4.DATASET_BUILDERS)["GSE118713"]()
    assert {r.label for r in gse118713} == {"MCF7-1", "MCF7-2", "MCF7-3", "TAMR-1", "TAMR-2", "TAMR-3"}

    gse111151 = dict(v4.DATASET_BUILDERS)["GSE111151"]()
    labels = {r.label for r in gse111151}
    assert "MCF-7" in labels and "MCF-7_Tam1" in labels
    assert "T-47D_Tam2" in labels  # T-47D has two independently-derived sublines

    gse245601 = dict(v4.DATASET_BUILDERS)["GSE245601"]()
    assert {r.label for r in gse245601} == {"T02 Ctrl", "T02 TAM", "T03 Ctrl", "T03 TAM", "T07 Ctrl", "T07 TAM"}


def test_gse240112_is_unpaired():
    rows = dict(v4.DATASET_BUILDERS)["GSE240112"]()
    assert all(r.pair_anchor is None for r in rows), "GSE240112 primary/recurrent tumours must never be bracketed as pairs"
    assert {r.state for r in rows} == {"baseline", "recurrent"}


def test_gse245601_pairing_is_genuine_within_patient():
    rows = dict(v4.DATASET_BUILDERS)["GSE245601"]()
    tam_rows = [r for r in rows if r.pair_anchor is not None]
    assert len(tam_rows) == 3
    for r in tam_rows:
        patient_prefix = r.label.split()[0]
        assert r.pair_anchor == f"{patient_prefix} Ctrl"


def test_gse111151_relationships_match_frozen_metadata():
    rows = dict(v4.DATASET_BUILDERS)["GSE111151"]()
    derivatives = [r for r in rows if r.pair_anchor is not None]
    assert len(derivatives) == 7  # 7 independently-derived TamR sublines
    parentals = {r.label for r in rows if r.pair_anchor is None}
    for r in derivatives:
        assert r.pair_anchor in parentals


def test_zscore_is_visualization_only_transform():
    """Within each dataset block, each gene's z-scored column has mean
    ~0 and, for N>=3 with nonzero variance, sample SD ~1 -- confirming
    this is a plain z-score and nothing else was silently altered."""
    for dataset, builder in v4.DATASET_BUILDERS:
        rows = builder()
        for gene in FOCUS_FOUR:
            z = np.array([r.values[gene] for r in rows])
            assert abs(z.mean()) < 1e-6
            assert abs(z.std(ddof=1) - 1.0) < 1e-6


def test_no_numeric_cell_annotation_layer_in_source():
    src_text = Path("src/poster_hero_heatmap_v4.py").read_text()
    assert "ax.annotate(" not in src_text
    numeric_cell_text = re.search(r"ax\.text\([^)]*\{.*:.2f\}", src_text)
    assert numeric_cell_text is None


def test_output_files_exist_and_are_non_degenerate():
    for ext in ("png", "pdf", "svg"):
        path = FIGURES / f"HERO_sample_level_heatmap_v4.{ext}"
        assert path.exists(), f"{path} was not written"
        assert path.stat().st_size > 0


def test_png_has_a_sane_minimum_resolution():
    with Image.open(FIGURES / "HERO_sample_level_heatmap_v4.png") as image:
        width, height = image.size
        assert width >= 1200
        assert height >= 1600  # tall, many-row figure


def test_note_exists_and_documents_transform():
    note = Path("results/reports/poster_hero_heatmap_v4/NOTE.md")
    assert note.exists()
    text = note.read_text().lower()
    assert "z-score" in text or "z =" in text
    assert "unpaired" in text


def test_frozen_shortlist_sha256_unchanged():
    digest = hashlib.sha256(SHORTLIST.read_bytes()).hexdigest()
    assert digest == FROZEN_SHA256
