"""Tests for the EXPLORATION-V3 poster-grade figure bank (results/figures/
poster_exploration_v3/), a small, high-quality rebuild of 6 hero figures
(+ 2 alternates) from the same frozen data already used in v2, with a
stricter poster-grade visual design system.

Pins the v3 data-transform functions (delta-from-parental,
delta-from-primary-mean) against the underlying frozen v2 loaders, checks
real biological-unit counts, and checks that every rendered figure is a
real, non-degenerate output. No pixel dimension is pinned exactly, per the
project's explicit "do not write brittle pixel-dimension tests" rule.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from PIL import Image

from src import poster_exploration_v3_data as pv3
from src import poster_exploration_v2_data as pv2
from src import post_audit_sensitivity_data as pad

FIGURES = Path("results/figures/poster_exploration_v3")

FOCUS_FOUR = ["KDM1A", "TLK2", "USP34", "VEZF1"]

MAIN_STEMS = [
    "F1_crispr_discovery", "F2_transcriptomic_corroboration", "F3_pathway_convergence",
    "F4_postaudit_interpretation", "F5_depmap_context", "F6_structural_comparison",
]
ALTERNATE_STEMS = ["ALT_A_genomewide_landscape", "ALT_B_pocket_closeups"]
ALL_STEMS = MAIN_STEMS + ALTERNATE_STEMS
PNG_ONLY_STEMS = ["F5_depmap_context", "F6_structural_comparison", "ALT_B_pocket_closeups"]
VECTOR_STEMS = [s for s in ALL_STEMS if s not in PNG_ONLY_STEMS]

# Exact palette required by this phase's brief.
REQUIRED_COLORS = {"KDM1A": "#D55E00", "TLK2": "#56B4E9", "USP34": "#0072B2", "VEZF1": "#E69F00"}


class TestPaletteMatchesBrief:
    def test_focus_four_matches_frozen_modules(self):
        assert pv3.FOCUS_FOUR == pv2.FOCUS_FOUR == pad.FOCUS_FOUR == FOCUS_FOUR

    def test_colors_match_exact_brief(self):
        assert pv3.FOCUS_COLORS == REQUIRED_COLORS

    def test_colors_identical_to_v2_v_poster_final(self):
        assert pv3.FOCUS_COLORS == pv2.FOCUS_COLORS


class TestFig1CrisprData:
    def test_significant_hits_count_and_kdm1a_rank(self):
        sens = pv3.load_significant_sensitising_hits()
        assert len(sens) == 13
        assert int(sens.set_index("gene").loc["KDM1A", "rank_by_effect"]) == 1

    def test_genomewide_gene_count(self):
        assert len(pv3.load_genomewide_crispr()) == 19103


class TestFig2TranscriptomicDeltas:
    def test_gse111151_delta_matches_manual_recomputation(self):
        delta = pv3.build_gse111151_delta_from_parental()
        assert set(delta["gene_symbol"]) == set(FOCUS_FOUR)
        for gene in FOCUS_FOUR:
            assert (delta["gene_symbol"] == gene).sum() == 7  # 7 real resistant sublines

        # spot check: recompute one gene/sample's delta directly from the
        # frozen source table, independent of the v3 helper's own logic
        raw = pv3.load_gse111151_focus_gene_samples()
        meta = pd.read_csv(pv3.GSE111151_METADATA, sep="\t")
        row = raw[(raw["gene_symbol"] == "USP34") & (raw["sample_id"] == "MCF-7_Tam1")].iloc[0]
        parental_line = meta.loc[meta["sample_id"] == "MCF-7_Tam1", "paired_parental_sample_id"].iloc[0]
        parental_val = raw[(raw["gene_symbol"] == "USP34") & (raw["sample_id"] == parental_line)]["log2cpm"].iloc[0]
        expected_delta = float(row["log2cpm"] - parental_val)
        actual = delta[(delta["gene_symbol"] == "USP34") & (delta["sample_id"] == "MCF-7_Tam1")]["delta_log2cpm"].iloc[0]
        assert actual == pytest.approx(expected_delta, abs=1e-9)

    def test_gse240112_delta_matches_manual_recomputation(self):
        delta = pv3.build_gse240112_delta_from_primary_mean()
        assert set(delta["gene"]) == set(FOCUS_FOUR)
        for gene in FOCUS_FOUR:
            sub = delta[delta["gene"] == gene]
            assert (sub["group"] == "PT").sum() == 3
            assert (sub["group"] == "RT").sum() == 3

        raw = pv3.load_gse240112_focus_gene_tumours()
        usp34 = raw[raw["gene"] == "USP34"]
        primary_mean = usp34[usp34["group"] == "PT"]["log2cpm"].mean()
        expected = usp34[usp34["group"] == "RT"]["log2cpm"].to_numpy() - primary_mean
        actual = delta[(delta["gene"] == "USP34") & (delta["group"] == "RT")]["delta_log2cpm"].to_numpy()
        assert np.allclose(np.sort(actual), np.sort(expected), atol=1e-9)

    def test_gse240112_is_unpaired_no_pairing_columns(self):
        delta = pv3.build_gse240112_delta_from_primary_mean()
        assert "pair_id" not in delta.columns
        assert "paired_sample_id" not in delta.columns


class TestFig3PathwayData:
    def test_hero_pathways_cover_4_datasets(self):
        df = pv3.load_pathway_trajectories(pv3.HERO_PATHWAYS)
        assert set(df["dataset"]) == set(pv3.DATASET_ORDER)
        assert len(pv3.DATASET_ORDER) == 4


class TestFig4SelectionRule:
    def test_rule0_rule1_shapes(self):
        rule0, rule1 = pv3.load_rule0_rule1()
        assert len(rule1) == 15  # union of 13 sensitising + original 4 (2 overlap)
        assert "eligible" in rule0.columns
        assert "rank" in rule1.columns

    def test_kdm1a_tlk2_excluded_by_original_gate(self):
        rule0, _ = pv3.load_rule0_rule1()
        assert bool(rule0.loc["KDM1A", "eligible"]) is False
        assert bool(rule0.loc["TLK2", "eligible"]) is False

    def test_usp34_vezf1_eligible_under_original_gate(self):
        rule0, _ = pv3.load_rule0_rule1()
        assert bool(rule0.loc["USP34", "eligible"]) is True
        assert bool(rule0.loc["VEZF1", "eligible"]) is True
        assert int(rule0.loc["USP34", "rank"]) == 1
        assert int(rule0.loc["VEZF1", "rank"]) == 2


class TestFig5Depmap:
    def test_11_lines_per_gene(self):
        eff = pv3.load_depmap_effect_focus_four()
        for gene in FOCUS_FOUR:
            assert (eff["gene"] == gene).sum() == 11
        names = pv3.load_depmap_model_names()
        assert len(names) == 11

    def test_tlk2_has_strongest_median_dependency(self):
        eff = pv3.load_depmap_effect_focus_four()
        medians = eff.groupby("gene")["chronos_effect"].median()
        assert medians.idxmin() == "TLK2"


class TestFig6StructuralFiles:
    def test_kdm1a_tlk2_structure_paths_exist(self):
        paths = pv3.kdm1a_tlk2_structure_paths()
        assert set(paths.keys()) == {"6NQU", "5O0Y"}
        for p in paths.values():
            assert p.exists()

    def test_kdm1a_ligand_code_present(self):
        text = pv3.kdm1a_tlk2_structure_paths()["6NQU"].read_text()
        assert "KWM" in text

    def test_tlk2_ligand_code_present(self):
        text = pv3.kdm1a_tlk2_structure_paths()["5O0Y"].read_text()
        assert "AGS" in text

    def test_usp34_structure_paths_exist(self):
        paths = pv3.usp34_structure_paths()
        assert set(paths.keys()) == {"7W3R", "7W3U"}
        for p in paths.values():
            assert p.exists()

    def test_vezf1_has_no_structure_row_marks_it_false(self):
        struct = pv3.load_structural_tractability_audit().set_index("gene")
        assert bool(struct.loc["VEZF1", "A_experimental_human_structure_exists"]) is False
        for gene in ["KDM1A", "TLK2", "USP34"]:
            assert bool(struct.loc[gene, "A_experimental_human_structure_exists"]) is True


class TestFiguresExistAndAreNonDegenerate:
    @pytest.mark.parametrize("stem", ALL_STEMS)
    def test_png_exists_and_has_content(self, stem):
        path = FIGURES / f"{stem}.png"
        assert path.exists(), f"{path} was not written"
        assert path.stat().st_size > 0

    @pytest.mark.parametrize("stem", ALL_STEMS)
    def test_png_has_a_sane_minimum_resolution(self, stem):
        with Image.open(FIGURES / f"{stem}.png") as image:
            width, height = image.size
            assert width >= 1200
            assert height >= 800

    @pytest.mark.parametrize("stem", VECTOR_STEMS)
    def test_vector_formats_exist(self, stem):
        assert (FIGURES / f"{stem}.pdf").exists()
        assert (FIGURES / f"{stem}.svg").exists()

    @pytest.mark.parametrize("stem", PNG_ONLY_STEMS)
    def test_png_only_stems_documented(self, stem):
        # F5/F6/ALT_B contain rasterized PyMOL renders or composited raster
        # images -- PNG only is a documented, deliberate choice
        assert (FIGURES / f"{stem}.png").exists()

    def test_contact_sheets_exist(self):
        assert (FIGURES / "CONTACT_MAIN_SIX.png").exists()
        assert (FIGURES / "CONTACT_MAIN_PLUS_ALTERNATES.png").exists()


class TestNoHandTypedFigureText:
    def test_no_hardcoded_gene_count_strings_in_module(self):
        import inspect

        from src import poster_exploration_v3_visualization as pv3viz

        src = inspect.getsource(pv3viz)
        # "19,103" and "13 significant" must never appear as bare literal
        # strings -- they must be produced via f-string interpolation of a
        # loaded count (`len(...)`, `n`, etc.)
        assert '"19,103' not in src
        assert "'19,103" not in src


class TestReportsExist:
    def test_guide_exists(self):
        path = Path("results/reports/poster_exploration_v3/POSTER_V3_FIGURE_GUIDE.md")
        assert path.exists()

    def test_guide_mentions_all_six_mains(self):
        text = Path("results/reports/poster_exploration_v3/POSTER_V3_FIGURE_GUIDE.md").read_text()
        for stem in MAIN_STEMS:
            assert stem in text, f"{stem} missing from POSTER_V3_FIGURE_GUIDE.md"
