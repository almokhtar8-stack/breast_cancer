"""Palette isolation and colour-vision safety for the final poster figure set.

post_freeze_exploratory."""

from __future__ import annotations

import inspect
import re
from pathlib import Path

import numpy as np
import pytest

from src import (
    poster_corroboration_final,
    poster_dependency_final,
    poster_final_common,
    poster_network_final,
    poster_palette,
    poster_pathway_final,
    poster_screen_final,
    poster_structure_final,
    poster_workflow_final,
)
from src.poster_palette import (
    CHROME,
    DIVERGING,
    FIGURE_ALLOWED,
    GENE_COLOURS,
    NEUTRAL,
    palette_cvd_report,
    print_point_size,
    simulate_cvd,
)

RENDERERS = [poster_workflow_final, poster_screen_final, poster_corroboration_final,
             poster_pathway_final, poster_network_final, poster_dependency_final,
             poster_structure_final]
HEX = re.compile(r"#[0-9A-Fa-f]{6}")


def test_palette_matches_the_specification_exactly():
    assert GENE_COLOURS == {"USP34": "#6A3D9A", "KDM1A": "#D55E00",
                            "TLK2": "#009E73", "VEZF1": "#56B4E9"}
    assert NEUTRAL["ink"] == "#262626" and NEUTRAL["backdrop"] == "#c9c9c9"
    assert DIVERGING == ["#2E6C8E", "#7FAFC4", "#F3EEE4", "#E3A180", "#C1543A"]


@pytest.mark.parametrize("module", RENDERERS + [poster_final_common], ids=lambda m: m.__name__)
def test_no_hex_literal_outside_the_palette_module(module):
    """Colour is defined once. A renderer that hard-codes a hex has drifted."""
    found = set(HEX.findall(inspect.getsource(module)))
    assert not found, f"{module.__name__} defines colours directly: {sorted(found)}"


@pytest.mark.parametrize("module", RENDERERS, ids=lambda m: m.__name__)
def test_no_chrome_colour_reaches_a_figure(module):
    src = inspect.getsource(module)
    assert "CHROME" not in src, f"{module.__name__} imports poster chrome into a figure"


@pytest.mark.parametrize("module", RENDERERS, ids=lambda m: m.__name__)
def test_no_colour_imported_from_a_v1_or_v2_renderer(module):
    src = inspect.getsource(module)
    for banned in ("FOCUS_COLORS", "STATE_COLORS", "poster_hero_heatmap"):
        assert banned not in src, f"{module.__name__} still imports {banned}"


def test_figure_allowed_set_covers_exactly_the_data_tier_colours():
    """Every data-tier colour is drawable; every chrome-only colour is not.

    `#6A3D9A` is deliberately in both dictionaries -- it is the poster's violet
    AND USP34's gene colour -- so it is excluded from the chrome-only set
    rather than treated as a leak."""
    assert {c.upper() for c in GENE_COLOURS.values()} <= FIGURE_ALLOWED

    # Exactly two colours appear in both tiers, and both are deliberate in the
    # specification: the violet is the poster's chrome violet AND USP34's gene
    # colour, and the warm white is the chrome warm white AND the neutral
    # midpoint of the diverging ramp. Pinning them here means a THIRD overlap
    # -- an actual leak of poster chrome into a figure -- fails the test.
    overlap = {c.upper() for c in CHROME.values()} & FIGURE_ALLOWED
    assert overlap == {"#6A3D9A", "#F3EEE4"}, sorted(overlap)


def test_chrome_is_never_exported_into_the_common_figure_helpers():
    assert "CHROME" not in inspect.getsource(poster_final_common)


# --- colour-vision -----------------------------------------------------------
def test_all_four_genes_stay_distinguishable_under_both_simulations():
    rep = palette_cvd_report()
    for vision in ("deuteranopia", "protanopia"):
        sub = rep[rep.vision == vision]
        assert (sub["delta_e76"] > 20).all(), sub.to_string()


def test_purple_and_light_blue_separate_on_lightness_not_hue():
    rep = palette_cvd_report()
    pv = rep[(rep.gene_a == "USP34") & (rep.gene_b == "VEZF1")]
    assert (pv[pv.vision != "normal"]["delta_L"] > 25).all()


def test_cvd_simulation_preserves_white_and_black():
    for kind in ("deuteranopia", "protanopia"):
        assert np.allclose(simulate_cvd(np.array([1.0, 1.0, 1.0]), kind), 1.0, atol=0.02)
        assert np.allclose(simulate_cvd(np.array([0.0, 0.0, 0.0]), kind), 0.0, atol=1e-6)


def test_every_rendered_figure_has_both_cvd_simulations():
    sim_dir = Path("results/figures/poster_final/cvd_simulation")
    assert sim_dir.is_dir()
    for key in ("F1_methods_workflow", "F2_screen_certainty", "F3_candidate_corroboration",
                "F4_programme_signal", "F5_network_connectivity", "F6_baseline_dependency",
                "F7_reachability"):
        for kind in ("deuteranopia", "protanopia"):
            assert (sim_dir / f"{key}_{kind}.png").exists(), f"missing {key} {kind}"


# --- print legibility ---------------------------------------------------------
def test_print_point_size_maths():
    # a 15.5 inch figure placed 380 mm wide shrinks slightly
    assert print_point_size(21.0, 15.5, 380) == pytest.approx(21.0 * (380 / 25.4) / 15.5)


def test_every_figure_clears_the_20pt_print_floor():
    import pandas as pd

    m = pd.read_csv("results/figures/poster_final/figure_manifest.tsv", sep="\t")
    failing = m.loc[~m["meets_20pt_floor"], ["figure", "smallest_printed_pt"]]
    assert failing.empty, f"below the 20 pt floor:\n{failing.to_string(index=False)}"
    assert (m["smallest_printed_pt"] >= poster_palette.MIN_PRINT_PT).all()
