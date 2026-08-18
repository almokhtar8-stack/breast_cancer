"""Tests for the v2 candidate volcano figure (two variants). post_freeze_exploratory."""

from __future__ import annotations

import copy
import itertools

import numpy as np
import pandas as pd
import pytest

import src.poster_candidate_volcano_v2 as v2
from src.poster_candidate_volcano_v2 import (
    CANDIDATES,
    COINCIDENCE_STEP,
    COINCIDENCE_TOL,
    EXPECTED_SIGNIFICANT,
    GENE_COLOURS,
    NEUTRAL,
    SIG_FDR,
    STATUS_LABEL,
    ZOOM_XLIM,
    ZOOM_YLIM,
    build_variant_a,
    build_variant_b,
    candidate_points,
    load_panels,
    palette_cvd_report,
    simulate_cvd,
    spread_coincident,
    variant_a_limits,
    verify_against_frozen,
)


@pytest.fixture(scope="module")
def panels():
    return load_panels()


@pytest.fixture(scope="module")
def verification(panels):
    return verify_against_frozen(panels)


# --- gate ------------------------------------------------------------------
def test_gate_still_passes_on_all_16_values(verification):
    assert len(verification) == 16
    assert verification["fdr_match"].all()
    lfc = verification[verification["reference_log2fc"] != ""]
    assert len(lfc) == 2 and lfc["log2fc_match"].all()


def test_gate_still_fails_loudly_on_a_corrupted_value(panels):
    broken = copy.deepcopy(panels)
    broken[2].candidate_rows.loc["VEZF1", "fdr"] = 0.5
    with pytest.raises(ValueError, match="verification gate FAILED"):
        verify_against_frozen(broken)


def test_exactly_two_points_are_significant(panels):
    sig = [(p.accession, g) for p in panels for g in CANDIDATES
           if float(p.candidate_rows.loc[g, "fdr"]) < SIG_FDR]
    assert len(sig) == EXPECTED_SIGNIFICANT == 2
    assert set(sig) == {("GSE118713", "USP34"), ("GSE240112", "VEZF1")}
    assert SIG_FDR == 0.05


# --- colour ----------------------------------------------------------------
def test_gene_colours_match_the_project_palette_exactly():
    assert GENE_COLOURS == {"USP34": "#6A3D9A", "KDM1A": "#D55E00",
                            "TLK2": "#009E73", "VEZF1": "#56B4E9"}
    assert set(GENE_COLOURS) == set(CANDIDATES)


def test_v2_does_not_take_colours_from_v1_or_the_heatmap_renderer():
    import inspect
    source = inspect.getsource(v2)
    assert "FOCUS_COLORS" not in source
    assert "poster_hero_heatmap" not in source


def test_only_palette_and_neutral_ladder_colours_are_named_in_the_module():
    """No magenta, no dark violet, nothing outside the declared colours
    (white is the page and marker-edge ground, not a plot colour)."""
    import inspect, re
    hexes = {h.upper() for h in re.findall(r"#[0-9A-Fa-f]{6}", inspect.getsource(v2))}
    allowed = {c.upper() for c in GENE_COLOURS.values()} | {c.upper() for c in NEUTRAL.values()}
    assert hexes <= allowed, hexes - allowed


def test_cvd_simulation_keeps_all_four_genes_distinguishable():
    rep = palette_cvd_report()
    for vision in ("deuteranopia", "protanopia"):
        sub = rep[rep.vision == vision]
        assert (sub["delta_e76"] > 20).all(), sub
    # the pair the brief flags: purple vs light blue separate on lightness
    pv = rep[(rep.gene_a == "USP34") & (rep.gene_b == "VEZF1")]
    assert (pv["delta_L"] > 25).all()


def test_cvd_simulation_matrices_preserve_white_and_black():
    for kind in ("deuteranopia", "protanopia"):
        assert np.allclose(simulate_cvd(np.array([1.0, 1.0, 1.0]), kind), 1.0, atol=0.02)
        assert np.allclose(simulate_cvd(np.array([0.0, 0.0, 0.0]), kind), 0.0, atol=1e-6)


# --- shared axes / visibility ------------------------------------------------
def test_variant_a_limits_contain_every_gene_in_every_panel(panels):
    xlim, ylim = variant_a_limits(panels)
    for p in panels:
        x = p.df["log2fc"].to_numpy(); y = -np.log10(p.df["fdr"].to_numpy())
        assert x.min() >= xlim[0] and x.max() <= xlim[1]
        assert y.min() >= ylim[0] and y.max() <= ylim[1]


def test_no_candidate_falls_outside_the_visible_range_in_either_variant(panels):
    xlim, ylim = variant_a_limits(panels)
    for p in panels:
        for g, (x, y) in candidate_points(p).items():
            assert xlim[0] <= x <= xlim[1] and ylim[0] <= y <= ylim[1]          # A main
            assert ZOOM_XLIM[0] <= x <= ZOOM_XLIM[1] and ZOOM_YLIM[0] <= y <= ZOOM_YLIM[1]  # A zoom, B


def test_shared_axis_limits_are_identical_across_the_four_panels(panels, tmp_path):
    import matplotlib.pyplot as plt

    for builder in (build_variant_a, build_variant_b):
        plt.close("all")
        # render, then inspect every axes of the (still-open) figure before it is closed
        captured = {}
        orig_save = v2._save

        def spy(fig, stub):
            captured["axes"] = [(tuple(a.get_xlim()), tuple(a.get_ylim())) for a in fig.axes
                                if a.get_xlabel() == "log2 fold change"]
            return orig_save(fig, stub)

        v2._save = spy
        try:
            builder(panels, tmp_path / builder.__name__)
        finally:
            v2._save = orig_save
        lims = captured["axes"]
        assert len(lims) in (4, 8)
        # group by y-range: main row and zoom row are each internally identical
        for ylim in {l[1] for l in lims}:
            group = [l for l in lims if l[1] == ylim]
            assert len(group) == 4 and len(set(group)) == 1


# --- coincidence rule ----------------------------------------------------------
def test_coincident_points_are_spread_horizontally_and_never_vertically():
    pts = {"A": (0.00, 0.10), "B": (0.01, 0.11), "C": (0.02, 0.12), "D": (1.0, 1.0)}
    disp, recs = spread_coincident(pts)
    assert {r["gene"] for r in recs} == {"A", "B", "C"}
    assert disp["D"] == pts["D"]
    xs = [disp[g][0] for g in ("A", "B", "C")]
    assert np.allclose(np.diff(xs), COINCIDENCE_STEP)
    assert np.isclose(np.mean(xs), np.mean([pts[g][0] for g in ("A", "B", "C")]))
    for g in ("A", "B", "C"):
        assert disp[g][1] == pts[g][1]


def test_offset_never_changes_threshold_side_on_the_real_data(panels):
    sig_y = -np.log10(SIG_FDR)
    for p in panels:
        m = candidate_points(p); d, _ = spread_coincident(m)
        for g in CANDIDATES:
            assert (m[g][1] >= sig_y) == (d[g][1] >= sig_y)
            assert d[g][1] == m[g][1]


def test_real_data_offsets_are_the_documented_clusters(panels):
    clusters = {}
    for p in panels:
        _, recs = spread_coincident(candidate_points(p))
        for r in recs:
            clusters.setdefault(p.accession, set()).add(r["cluster"])
    assert clusters == {"GSE111151": {"KDM1A+TLK2"}, "GSE245601": {"TLK2+USP34+VEZF1"}}
    assert COINCIDENCE_TOL == 0.06 and COINCIDENCE_STEP == 0.14


# --- build + outputs ---------------------------------------------------------
def test_both_variants_build_in_all_three_formats(panels, tmp_path):
    for builder in (build_variant_a, build_variant_b):
        written, info = builder(panels, tmp_path / builder.__name__)
        assert set(written) == {"png", "pdf", "svg"}
        for path in written.values():
            assert path.exists() and path.stat().st_size > 0


def test_pdf_and_png_bytes_are_reproducible(panels, tmp_path):
    for builder in (build_variant_a, build_variant_b):
        a, _ = builder(panels, tmp_path / "one")
        b, _ = builder(panels, tmp_path / "two")
        assert a["pdf"].read_bytes() == b["pdf"].read_bytes()
        assert a["png"].read_bytes() == b["png"].read_bytes()


def test_main_writes_manifest_covering_both_variants(tmp_path):
    written = v2.main(out_dir=tmp_path)
    assert set(written) == {"variant_a_genomewide", "variant_b_candidates_only"}
    m = pd.read_csv(tmp_path / "volcano_v2_manifest.tsv", sep="\t")
    assert list(m["variant"]) == ["variant_a_genomewide", "variant_b_candidates_only"]
    assert (m["analysis_status"] == STATUS_LABEL).all()
    assert (m["post_freeze"] == "yes").all()
    assert (m["gse118713_contrast"] == "TAMR_vs_MCF7").all()
    for ext in ("png", "pdf", "svg"):
        assert m[f"sha256_{ext}"].str.len().eq(64).all()
        assert m[f"sha256_{ext}"].nunique() == 2  # the two variants differ
    for g, c in GENE_COLOURS.items():
        assert (m[f"colour_{g}"] == c).all()
    assert m.loc[0, "zoom_xlim"] == "-1.5,1.5" and m.loc[0, "zoom_ylim"] == "0.0,2.5"
    assert pd.isna(m.loc[1, "zoom_xlim"])  # variant B has no zoom row
    assert m.loc[1, "shared_xlim"] == "-1.500,1.500" and m.loc[1, "shared_ylim"] == "0.000,2.500"
    for name in ("candidate_values_plotted.tsv", "verification_against_frozen.tsv",
                 "cosmetic_offsets.tsv", "cvd_palette_simulation.tsv"):
        t = pd.read_csv(tmp_path / name, sep="\t")
        assert (t["analysis_status"] == STATUS_LABEL).all()
    assert (tmp_path / "cvd_simulation").is_dir()
    assert len(list((tmp_path / "cvd_simulation").glob("*.png"))) == 4
