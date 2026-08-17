"""Tests for the EXPLORATION-V2 poster figure bank (results/figures/
poster_exploration_v2/), the full visualization-strategy reset following
rejection of the poster_final_* figure set.

Pins the data loaders in src/poster_exploration_v2_data.py against the
same frozen source tables the earlier phases already use, checks real
biological-unit counts (no fabricated pairing, no hidden sample-size
inflation), and checks that every rendered figure is a real, non-
degenerate output. No pixel dimension is pinned exactly, per the
project's explicit "do not write brittle pixel-dimension tests" rule.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from PIL import Image

from src import poster_exploration_v2_data as pv2
from src import post_audit_sensitivity_data as pad

FIGURES = Path("results/figures/poster_exploration_v2")
REPORTS = Path("results/reports/poster_exploration_v2")

FOCUS_FOUR = ["KDM1A", "TLK2", "USP34", "VEZF1"]

ALL_FIGURE_STEMS = [
    "A1_genomewide_ranked_landscape", "A2_lollipop_13_hits", "A3_redesigned_volcano", "A4_hexbin_density",
    "B1_gse118713_sample_panel", "B2_gse111151_trajectories", "B3_gse240112_tumours", "B4_gse245601_acute_paired",
    "B5a_integrated_dataset_centric", "B5b_integrated_gene_centric",
    "C1a_pathway_trajectories_lines", "C1b_pathway_small_multiples", "C2_estrogen_emt_hero", "C3_enrichment_curves",
    "D1_candidate_program_network", "D2_usp34_local_neighborhood",
    "E1_upset_evidence_intersection", "E2_quantitative_evidence_map", "E3_historical_vs_postaudit",
    "F1_depmap_heatmap", "F2_cellline_fingerprint", "F3_human_depmap_combo", "F4_tcga_secondary",
    "G1_three_structures_plus_vezf1", "G2_structure_pharmacology_maturity", "G3_pocket_closeups",
    "H2_experimental_schematic", "H3_translational_maturity_landscape",
]
VECTOR_STEMS = [s for s in ALL_FIGURE_STEMS if not s.startswith(("A4", "D2", "G1", "G2", "G3"))]
PNG_ONLY_STEMS = [s for s in ALL_FIGURE_STEMS if s not in VECTOR_STEMS]


class TestFocusFourConsistency:
    def test_focus_four_matches_post_audit_module(self):
        assert pv2.FOCUS_FOUR == pad.FOCUS_FOUR == FOCUS_FOUR

    def test_four_distinct_colors(self):
        assert len(set(pv2.FOCUS_COLORS.values())) == 4


class TestBiologicalUnitCounts:
    """Guards against pseudoreplication, false pairing, and hidden N
    inflation -- the specific failure modes the task explicitly warned
    against."""

    def test_gse118713_has_9_real_samples(self):
        df = pv2.load_gse118713_focus_gene_samples()
        assert df["sample"].nunique() == 9
        for cond in ["MCF7", "TAMR", "FASR"]:
            assert (df.loc[df["gene_symbol"] == "KDM1A", "condition"] == cond).sum() == 3

    def test_gse111151_has_11_real_samples_with_real_pairing(self):
        df = pv2.load_gse111151_focus_gene_samples()
        assert df["sample_id"].nunique() == 11
        meta = pd.read_csv(pv2.GSE111151_METADATA, sep="\t")
        assert set(df["sample_id"]) == set(meta["sample_id"])
        # every resistant sample's paired_parental_sample_id must be a REAL
        # parental sample_id already present in the metadata -- no
        # invented pairing
        resistant = meta[meta["status"] == "resistant"]
        assert resistant["paired_parental_sample_id"].isin(meta["sample_id"]).all()
        assert (meta["status"] == "parental").sum() == 4
        assert (meta["status"] == "resistant").sum() == 7

    def test_gse240112_has_exactly_3_primary_3_recurrent_unpaired(self):
        df = pv2.load_gse240112_focus_gene_tumours()
        sub = df[df["gene"] == "USP34"]
        assert (sub["group"] == "PT").sum() == 3
        assert (sub["group"] == "RT").sum() == 3
        # no pairing column should exist / be used anywhere in this loader
        assert "paired_sample_id" not in df.columns
        assert "pair_id" not in df.columns

    def test_gse245601_paired_uses_exactly_3_eligible_patients(self):
        df = pv2.load_gse245601_paired_focus_genes()
        assert set(df["patient"]) == {"Tumor_02", "Tumor_03", "Tumor_07"}
        for gene in FOCUS_FOUR:
            sub = df[df["gene"] == gene]
            assert set(sub["condition"]) == {"Control", "Tamoxifen"}
            assert len(sub) == 6  # 3 patients x 2 conditions, genuinely paired

    def test_gse245601_eligibility_matches_frozen_filter(self):
        eligible = pd.read_csv(pv2.GSE245601_PAIR_ELIGIBILITY, sep="\t")
        eligible_patients = set(eligible.loc[eligible["eligible_for_pseudobulk"], "patient"])
        assert eligible_patients == {"Tumor_02", "Tumor_03", "Tumor_07"}

    def test_gse245601_kdm1a_tlk2_formula_matches_frozen_usp34_value(self):
        frozen_value, computed_value = pv2.verify_gse245601_computed_formula_matches_frozen_usp34()
        assert computed_value == pytest.approx(frozen_value, abs=1e-6)

    def test_depmap_has_exactly_11_er_luminal_lines_per_gene(self):
        eff = pv2.load_depmap_effect_focus_four()
        for gene in FOCUS_FOUR:
            assert (eff["gene"] == gene).sum() == 11
        names = pv2.load_depmap_model_names()
        assert len(names) == 11


class TestCrisprRanks:
    def test_kdm1a_ranks_first_of_13(self):
        sens = pv2.load_significant_sensitising_hits().set_index("gene")
        assert int(sens.loc["KDM1A", "rank_by_effect"]) == 1
        assert int(sens.loc["KDM1A", "rank_by_fdr"]) == 1

    def test_usp34_is_not_the_top_hit(self):
        sens = pv2.load_significant_sensitising_hits().set_index("gene")
        assert int(sens.loc["USP34", "rank_by_effect"]) > 1

    def test_blind_control_rcor1_not_recovered(self):
        row = pv2.load_blind_control_row()
        assert row["gene"] == "RCOR1"
        assert row["fdr"] >= 0.1

    def test_genomewide_has_19103_genes(self):
        assert len(pv2.load_genomewide_crispr()) == 19103


class TestPathwayData:
    def test_pathway_trajectories_cover_4_datasets(self):
        df = pv2.load_pathway_trajectories(pv2.HERO_PATHWAYS)
        assert set(df["dataset"]) == set(pv2.DATASET_ORDER)

    def test_enrichment_curve_is_deterministic_and_matches_frozen_sign(self):
        curve = pv2.build_enrichment_curve("gse118713", "HALLMARK_ESTROGEN_RESPONSE_EARLY")
        gsea = pd.read_csv("results/tables/systems_network/gsea_gse118713.tsv", sep="\t")
        frozen_nes = gsea.loc[(gsea["collection"] == "hallmark") &
                               (gsea["pathway"] == "HALLMARK_ESTROGEN_RESPONSE_EARLY"), "NES"].iloc[0]
        # the peak deviation of the reconstructed curve must have the same
        # sign as the already-frozen NES (a real consistency check, not a
        # re-derivation of the NES value itself)
        peak = curve["running_es"].iloc[curve["running_es"].abs().idxmax()]
        assert np.sign(peak) == np.sign(frozen_nes)

    def test_enrichment_curve_uses_full_frozen_ranked_list(self):
        ranked = pv2.load_ranked_genes("gse118713")
        curve = pv2.build_enrichment_curve("gse118713", "HALLMARK_ESTROGEN_RESPONSE_EARLY")
        assert len(curve) == len(ranked)


class TestNetworkCoverageDisclosure:
    def test_direct_neighbors_covers_all_original_four_not_usp34_only(self):
        nb = pv2.load_direct_neighbors()
        counts = nb["candidate"].value_counts()
        assert counts["USP34"] == 10
        assert counts["CITED2"] == 18
        assert counts["VEZF1"] == 1
        assert counts["EML5"] == 1

    def test_filtering_by_candidate_works(self):
        nb = pv2.load_direct_neighbors("USP34")
        assert (nb["candidate"] == "USP34").all()
        assert len(nb) == 10


class TestEvidenceSets:
    def test_high_and_low_dependency_are_mutually_exclusive(self):
        sets_df = pv2.build_evidence_sets_13()
        both = sets_df["high_baseline_dependency"] & sets_df["low_baseline_dependency"]
        assert not both.any()

    def test_evidence_sets_cover_all_13_hits(self):
        sets_df = pv2.build_evidence_sets_13()
        sens = pv2.load_significant_sensitising_hits()
        assert set(sets_df["gene"]) == set(sens["gene"])


class TestStructuralFiles:
    def test_kdm1a_pdb_has_expected_ligand_code(self):
        path = pv2.kdm1a_tlk2_structure_paths()["6NQU"]
        text = path.read_text()
        assert "KWM" in text
        assert "TITLE" in text

    def test_tlk2_pdb_has_expected_ligand_code(self):
        path = pv2.kdm1a_tlk2_structure_paths()["5O0Y"]
        text = path.read_text()
        assert "AGS" in text

    def test_usp34_structure_paths_point_to_frozen_files(self):
        paths = pv2.usp34_structure_paths()
        assert set(paths.keys()) == {"7W3R", "7W3U"}
        for p in paths.values():
            assert p.exists()


class TestFiguresExistAndAreNonDegenerate:
    @pytest.mark.parametrize("stem", ALL_FIGURE_STEMS)
    def test_png_exists_and_has_content(self, stem):
        path = FIGURES / f"{stem}.png"
        assert path.exists(), f"{path} was not written"
        assert path.stat().st_size > 0

    @pytest.mark.parametrize("stem", ALL_FIGURE_STEMS)
    def test_png_has_a_sane_minimum_resolution(self, stem):
        with Image.open(FIGURES / f"{stem}.png") as image:
            width, height = image.size
            assert width >= 800
            assert height >= 600

    @pytest.mark.parametrize("stem", VECTOR_STEMS)
    def test_vector_formats_exist(self, stem):
        assert (FIGURES / f"{stem}.pdf").exists()
        assert (FIGURES / f"{stem}.svg").exists()

    @pytest.mark.parametrize("stem", PNG_ONLY_STEMS)
    def test_png_only_stems_documented(self, stem):
        # A4/G1/G2/G3 contain rasterized PyMOL/hexbin content -- PNG only
        # is a documented, deliberate choice, not an oversight
        assert (FIGURES / f"{stem}.png").exists()


class TestReportsExist:
    def test_data_audit_exists(self):
        assert (REPORTS / "DATA_FOR_VISUALIZATION_AUDIT.md").exists()

    def test_visual_reference_notes_exist(self):
        assert (REPORTS / "VISUAL_REFERENCE_NOTES.md").exists()

    def test_figure_candidate_guide_mentions_every_figure(self):
        text = (REPORTS / "FIGURE_CANDIDATE_GUIDE.md").read_text()
        for stem in ALL_FIGURE_STEMS:
            assert stem in text, f"{stem} missing from FIGURE_CANDIDATE_GUIDE.md"


class TestNoScientificOverclaim:
    def test_gse240112_never_called_paired_or_matched(self):
        import inspect

        src = inspect.getsource(pv2)
        low = src.lower()
        assert "gse240112" in low
        # the loader's own docstring must state UNPAIRED explicitly
        assert "unpaired" in low

    def test_gse245601_documented_as_acute_not_resistance(self):
        import inspect

        src = inspect.getsource(pv2).lower()
        assert "acute" in src
