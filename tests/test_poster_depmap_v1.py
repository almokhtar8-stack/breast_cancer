"""Small dedicated test for the DepMap figure -- data fidelity against the
frozen DepMap rules/threshold, correct scope wording, and output
existence."""

from __future__ import annotations

import hashlib
from pathlib import Path

import pytest
from PIL import Image

from src import post_audit_sensitivity_data as pad
from src import poster_depmap_v1 as dm

FIGURES = Path("results/figures/poster_depmap_v1")
SHORTLIST = Path("results/tables/evidence_freeze/THERAPEUTIC_SHORTLIST_FREEZE.tsv")
FROZEN_SHA256 = "b6990d7e20f1191df7ebbca18be7542cc312e6e8489f05af6e4b4d31e66fb9dc"


def test_exact_four_focus_genes():
    assert dm.FOCUS_FOUR == ["KDM1A", "TLK2", "USP34", "VEZF1"]
    # the retired shortlist members must not appear in this figure
    assert "EML5" not in dm.FOCUS_FOUR and "CITED2" not in dm.FOCUS_FOUR


def test_release_is_the_frozen_active_release():
    cfg = pad.load_config()
    assert dm.RELEASE == cfg["independent_validation"]["depmap"]["active_release"]


def test_strong_dependency_threshold_recovered_from_config_not_invented():
    cfg = pad.load_config()
    threshold = cfg["independent_validation"]["depmap"]["strong_dependency_probability_threshold"]
    assert threshold == 0.5
    table = dm.load_cellline_table()
    # the boolean column must be exactly "dependency probability > threshold"
    recomputed = table["dependency_probability"] > threshold
    assert (table["strongly_dependent"] == recomputed).all()


def test_exact_frozen_er_luminal_subset_of_eleven_lines():
    table = dm.load_cellline_table()
    lines = sorted(table["cell_line"].unique())
    assert len(lines) == 11
    assert lines == sorted([
        "CAMA1", "EFM19", "HCC1428", "KPL1", "MCF7", "MDAMB361",
        "MDAMB415", "MFM223", "T47D", "UACC3133", "ZR751",
    ])
    assert len(table) == 11 * 4


def test_every_chronos_value_loaded_from_frozen_data_matches_frozen_summary():
    """Per-gene medians/percentages derived here must equal the frozen
    post-audit DepMap summary computed by the project's own loader."""
    summary = dm.dependency_summary(dm.load_cellline_table()).set_index("gene")
    frozen = pad.load_depmap_summary_for_genes(dm.FOCUS_FOUR).set_index("gene")
    for gene in dm.FOCUS_FOUR:
        assert summary.loc[gene, "median_chronos"] == pytest.approx(
            frozen.loc[gene, "median_chronos_er_luminal"], abs=1e-9)
        assert summary.loc[gene, "pct_strongly_dependent"] / 100.0 == pytest.approx(
            frozen.loc[gene, "frac_strongly_dependent_er_luminal"], abs=1e-9)
        assert summary.loc[gene, "n_lines"] == 11


def test_no_hand_typed_chronos_or_percentage_values_in_rendering_logic():
    src_text = Path("src/poster_depmap_v1.py").read_text()
    for literal in ("-0.808", "-0.137", "-0.199", "-0.062", "81.8", "27.3", "0.8182"):
        assert literal not in src_text


def test_percentages_calculated_programmatically():
    summary = dm.dependency_summary(dm.load_cellline_table()).set_index("gene")
    for gene in dm.FOCUS_FOUR:
        n_strong = summary.loc[gene, "n_strongly_dependent"]
        n_lines = summary.loc[gene, "n_lines"]
        assert summary.loc[gene, "pct_strongly_dependent"] == pytest.approx(100.0 * n_strong / n_lines)


def test_expected_dependency_profile_per_candidate():
    summary = dm.dependency_summary(dm.load_cellline_table()).set_index("gene")
    pct = summary["pct_strongly_dependent"]
    assert pct["KDM1A"] == 0.0                  # low baseline dependency
    assert pct["USP34"] == 0.0                  # low baseline dependency
    assert pct["TLK2"] > 75.0                   # high baseline dependency
    assert 0.0 < pct["VEZF1"] < 50.0            # intermediate
    # TLK2 is the strongest of the FOUR focus genes by median Chronos
    assert summary["median_chronos"].idxmin() == "TLK2"


def test_scope_wording_does_not_overclaim_tlk2_beyond_the_four_focus_genes():
    text = (Path("src/poster_depmap_v1.py").read_text()
            + Path("results/reports/poster_depmap_v1/NOTE.md").read_text()).lower()
    for overclaim in ("strongest baseline dependency of all 13",
                      "strongest dependency gene overall",
                      "strongest of the 13"):
        assert overclaim not in text
    # the correct scoping must be stated
    assert "among the four" in text


def test_no_normal_tissue_safety_or_toxicity_claim():
    note = Path("results/reports/poster_depmap_v1/NOTE.md").read_text().lower()
    src_text = Path("src/poster_depmap_v1.py").read_text().lower()
    assert "not evidence of normal-tissue safety" in src_text
    assert "not normal-tissue safety" in note or "not normal tissue safety" in note
    assert "proves normal-tissue safety" not in note
    # baseline dependency must never be framed as toxicity
    assert "dependency = toxicity" not in note and "baseline toxicity" not in note


def test_note_states_depmap_is_not_tamoxifen_sensitisation():
    note = Path("results/reports/poster_depmap_v1/NOTE.md").read_text().lower()
    assert "not tamoxifen-specific sensitisation" in note
    assert "26q1" in note
    assert "0.5" in note


def test_figure_generation_succeeds_and_outputs_exist():
    stub = FIGURES / "DEPMAP_v1"
    dm.build_depmap_v1(stub)
    for ext in ("png", "pdf", "svg"):
        path = stub.with_suffix(f".{ext}")
        assert path.exists()
        assert path.stat().st_size > 0


def test_png_has_a_sane_minimum_resolution():
    with Image.open(FIGURES / "DEPMAP_v1.png") as image:
        width, _ = image.size
        assert width >= 2000


def test_frozen_shortlist_sha256_unchanged():
    digest = hashlib.sha256(SHORTLIST.read_bytes()).hexdigest()
    assert digest == FROZEN_SHA256
