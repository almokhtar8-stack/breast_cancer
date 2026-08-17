"""Small dedicated test for the DepMap v2 scatter -- verifies both plotted
coordinates come from the exact frozen sources reused from v1, and that the
sign-flip is a display transform only."""

from __future__ import annotations

import hashlib
from pathlib import Path

import pytest
from PIL import Image

from src import post_audit_sensitivity_data as pad
from src import poster_depmap_v1 as dm1
from src import poster_depmap_v2 as dm2

FIGURES = Path("results/figures/poster_depmap_v2")
SHORTLIST = Path("results/tables/evidence_freeze/THERAPEUTIC_SHORTLIST_FREEZE.tsv")
FROZEN_SHA256 = "b6990d7e20f1191df7ebbca18be7542cc312e6e8489f05af6e4b4d31e66fb9dc"


def test_exact_four_candidates():
    assert dm2.FOCUS_FOUR == ["KDM1A", "TLK2", "USP34", "VEZF1"]
    assert dm2.FOCUS_FOUR is dm1.FOCUS_FOUR  # reuses v1's definition, not a copy


def test_reuses_v1_frozen_loaders_not_a_reimplementation():
    assert dm2.RELEASE == dm1.RELEASE
    assert dm2.FOCUS_COLORS is dm1.FOCUS_COLORS


def test_crispr_effects_loaded_from_frozen_source():
    table = dm2.build_plot_table().set_index("gene")
    frozen = pad.load_significant_sensitising_hits().set_index("gene")
    for gene in dm2.FOCUS_FOUR:
        assert table.loc[gene, "crispr_effect"] == pytest.approx(frozen.loc[gene, "effect_size"])
        assert table.loc[gene, "crispr_effect"] < 0  # sensitising direction preserved in the data


def test_x_coordinate_is_exactly_minus_one_times_frozen_effect():
    table = dm2.build_plot_table().set_index("gene")
    for gene in dm2.FOCUS_FOUR:
        assert table.loc[gene, "sensitisation_strength"] == pytest.approx(
            -1.0 * table.loc[gene, "crispr_effect"], abs=1e-12)
        assert table.loc[gene, "sensitisation_strength"] > 0


def test_y_percentages_derived_from_frozen_eleven_line_subset():
    table = dm2.build_plot_table().set_index("gene")
    frozen = dm1.dependency_summary(dm1.load_cellline_table()).set_index("gene")
    for gene in dm2.FOCUS_FOUR:
        assert table.loc[gene, "n_lines"] == 11
        assert table.loc[gene, "n_strongly_dependent"] == frozen.loc[gene, "n_strongly_dependent"]
        assert table.loc[gene, "pct_strongly_dependent"] == pytest.approx(
            frozen.loc[gene, "pct_strongly_dependent"])
        # percentage must be the computed ratio, not an independent number
        assert table.loc[gene, "pct_strongly_dependent"] == pytest.approx(
            100.0 * table.loc[gene, "n_strongly_dependent"] / table.loc[gene, "n_lines"])


def test_expected_dependency_counts_per_candidate():
    table = dm2.build_plot_table().set_index("gene")
    assert table.loc["KDM1A", "n_strongly_dependent"] == 0
    assert table.loc["TLK2", "n_strongly_dependent"] == 9
    assert table.loc["USP34", "n_strongly_dependent"] == 0
    assert table.loc["VEZF1", "n_strongly_dependent"] == 3


def test_key_contrast_kdm1a_versus_tlk2_holds_in_the_data():
    table = dm2.build_plot_table().set_index("gene")
    # both are strong sensitisers (top of the frozen effect ordering)
    assert table.loc["KDM1A", "sensitisation_strength"] > table.loc["USP34", "sensitisation_strength"]
    assert table.loc["TLK2", "sensitisation_strength"] > table.loc["USP34", "sensitisation_strength"]
    # but they separate sharply on baseline dependency
    assert table.loc["KDM1A", "pct_strongly_dependent"] == 0.0
    assert table.loc["TLK2", "pct_strongly_dependent"] > 75.0


def test_no_scientific_coordinates_hand_typed():
    src_text = Path("src/poster_depmap_v2.py").read_text()
    for literal in ("2.167", "1.848", "1.602", "1.391", "81.8", "27.3", "-0.808", "-0.137"):
        assert literal not in src_text


def _drawn_figure_text() -> str:
    """Only the string literals the figure actually renders -- checking the
    whole module/NOTE would false-positive on sentences that explicitly
    NEGATE an overclaim (e.g. "no region is labelled a therapeutic window")."""
    import re
    src = Path("src/poster_depmap_v2.py").read_text()
    drawn = re.findall(r'(?:ax\.text|fig\.text|ax\.annotate|ax\.set_xlabel|ax\.set_ylabel)\((.*?)\)',
                        src, flags=re.S)
    return " ".join(drawn).lower()


def test_figure_draws_no_evaluative_or_overclaiming_labels():
    drawn = _drawn_figure_text()
    for word in ("safe", "superior", "toxic", "therapeutic window", "specific", "good", "bad"):
        assert word not in drawn, f"figure draws evaluative wording: {word!r}"


def test_note_states_the_required_scientific_cautions():
    note = Path("results/reports/poster_depmap_v2/NOTE.md").read_text().lower()
    for affirmative_overclaim in ("kdm1a is safe", "clinically superior", "tlk2 is bad", "tlk2 is toxic"):
        assert affirmative_overclaim not in note
    assert "not normal-tissue safety" in note
    assert "does not prove tamoxifen specificity" in note
    assert "among these\nfour focus genes" in note or "among these four focus genes" in note


def test_figure_generation_succeeds_and_outputs_exist():
    stub = FIGURES / "DEPMAP_v2"
    dm2.build_depmap_v2(stub)
    for ext in ("png", "pdf", "svg"):
        path = stub.with_suffix(f".{ext}")
        assert path.exists()
        assert path.stat().st_size > 0


def test_png_has_a_sane_minimum_resolution():
    with Image.open(FIGURES / "DEPMAP_v2.png") as image:
        width, height = image.size
        assert width >= 2000
        assert width > height


def test_v1_outputs_not_overwritten():
    assert Path("results/figures/poster_depmap_v1/DEPMAP_v1.png").exists()
    assert Path("results/reports/poster_depmap_v1/NOTE.md").exists()


def test_frozen_shortlist_sha256_unchanged():
    digest = hashlib.sha256(SHORTLIST.read_bytes()).hexdigest()
    assert digest == FROZEN_SHA256
