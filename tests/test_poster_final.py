"""Tests for the FINAL poster figure set (results/figures/poster_final/).

Pins the data loaders in src/poster_final_data.py against the same frozen
source tables the post-audit sensitivity analysis and figure bank already
use, and checks that every final figure is a real, non-degenerate render.
No pixel-dimension is pinned exactly (aspect-ratio / minimum-resolution
checks only), per the 2026-08-16 explicit instruction not to write brittle
pixel-dimension tests.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from PIL import Image

from src import poster_final_data as pfd2
from src import post_audit_sensitivity_data as pad

FIGURES = Path("results/figures/poster_final")
REPORTS = Path("results/reports/poster_final")

FOCUS_FOUR = ["KDM1A", "TLK2", "USP34", "VEZF1"]

MAIN_FIGURE_STEMS = [
    "F1_crispr_discovery",
    "F2_pathway_systems",
    "F3_candidate_evidence_integration",
    "F4_human_depmap_validation",
    "F5_USP34_structure_tractability",
    "F6_final_translational_framework",
]


class TestFocusFourAndColors:
    def test_focus_four_matches_post_audit_module(self):
        assert pfd2.FOCUS_FOUR == pad.FOCUS_FOUR == FOCUS_FOUR

    def test_four_distinct_colors_assigned(self):
        colors = [pfd2.FOCUS_COLORS[g] for g in FOCUS_FOUR]
        assert len(set(colors)) == 4

    def test_colors_are_valid_hex(self):
        for c in pfd2.FOCUS_COLORS.values():
            assert c.startswith("#") and len(c) == 7
            int(c[1:], 16)  # raises if not valid hex


class TestF1GenomewideData:
    def test_genomewide_has_19103_fitted_genes(self):
        df = pfd2.load_f1_genomewide()
        assert len(df) == 19103

    def test_gate1_significant_flag_matches_threshold(self):
        df = pfd2.load_f1_genomewide()
        assert (df.loc[df["gate1_significant"], "fdr"] < 0.1).all()
        assert (df.loc[~df["gate1_significant"], "fdr"] >= 0.1).all()

    def test_focus_gene_flag_covers_exactly_the_four(self):
        df = pfd2.load_f1_genomewide()
        flagged = set(df.loc[df["is_focus_gene"], "gene"])
        assert flagged == set(FOCUS_FOUR)

    def test_focus_gene_ranks_match_significant_sensitising_hits_table(self):
        ranks = pfd2.load_f1_focus_gene_ranks().set_index("gene")
        sens = pad.load_significant_sensitising_hits().set_index("gene")
        for gene in FOCUS_FOUR:
            assert ranks.loc[gene, "rank_by_effect"] == sens.loc[gene, "rank_by_effect"]
            assert ranks.loc[gene, "rank_by_fdr"] == sens.loc[gene, "rank_by_fdr"]
            assert ranks.loc[gene, "n_sensitising_hits"] == len(sens)

    def test_kdm1a_ranks_first_by_effect_and_fdr(self):
        ranks = pfd2.load_f1_focus_gene_ranks().set_index("gene")
        assert int(ranks.loc["KDM1A", "rank_by_effect"]) == 1
        assert int(ranks.loc["KDM1A", "rank_by_fdr"]) == 1

    def test_usp34_is_not_the_top_hit(self):
        ranks = pfd2.load_f1_focus_gene_ranks().set_index("gene")
        assert int(ranks.loc["USP34", "rank_by_effect"]) > 1
        assert int(ranks.loc["USP34", "rank_by_fdr"]) > 1

    def test_blind_control_row_is_rcor1_and_not_gate1_significant(self):
        row = pfd2.load_f1_blind_control_row()
        assert row is not None
        assert row["gene"] == "RCOR1"
        assert row["fdr"] >= 0.1  # real finding: RCOR1 was NOT recovered at Gate-1 threshold


class TestF2PathwayMatrix:
    def test_pathway_matrix_has_10_pathways_x_4_datasets(self):
        df = pfd2.load_f2_pathway_matrix()
        assert df["pathway_label"].nunique() == 10
        assert df["dataset_key"].nunique() == 4

    def test_every_row_has_a_dataset_category(self):
        df = pfd2.load_f2_pathway_matrix()
        assert df["dataset_category"].notna().all()

    def test_gse240112_is_never_in_the_resistance_model_category(self):
        df = pfd2.load_f2_pathway_matrix()
        gse240112_rows = df[df["dataset_key"] == "gse240112"]
        assert (gse240112_rows["dataset_category"] == "recurrence-associated (human tumour, unpaired)").all()

    def test_gse245601_is_never_in_the_resistance_model_category(self):
        df = pfd2.load_f2_pathway_matrix()
        acute_rows = df[df["dataset_key"] == "gse245601"]
        assert (acute_rows["dataset_category"] == "acute 12h (not resistance)").all()


class TestF3EvidenceMatrix:
    def test_evidence_matrix_covers_exactly_the_four_focus_genes(self):
        em = pfd2.load_f3_evidence_matrix()
        assert list(em["gene"]) == FOCUS_FOUR

    def test_evidence_matrix_matches_post_audit_build_evidence_matrix(self):
        em = pfd2.load_f3_evidence_matrix().set_index("gene")
        full = pad.build_evidence_matrix().set_index("gene")
        for gene in FOCUS_FOUR:
            assert em.loc[gene, "crispr_effect"] == pytest.approx(full.loc[gene, "crispr_effect"])
            assert em.loc[gene, "crispr_fdr"] == pytest.approx(full.loc[gene, "crispr_fdr"])

    def test_kdm1a_and_tlk2_have_no_tcga_value(self):
        em = pfd2.load_f3_evidence_matrix().set_index("gene")
        assert pd.isna(em.loc["KDM1A", "tcga_fdr"])
        assert pd.isna(em.loc["TLK2", "tcga_fdr"])

    def test_usp34_and_vezf1_have_a_tcga_value(self):
        em = pfd2.load_f3_evidence_matrix().set_index("gene")
        assert pd.notna(em.loc["USP34", "tcga_fdr"])
        assert pd.notna(em.loc["VEZF1", "tcga_fdr"])

    def test_structural_facets_cover_exactly_the_four_focus_genes(self):
        facets = pfd2.load_f3_structural_facets()
        assert list(facets["gene"]) == FOCUS_FOUR

    def test_structural_facet_states_are_valid_tokens(self):
        facets = pfd2.load_f3_structural_facets()
        valid = {"YES", "NO", "PARTIAL"}
        for col in ["structure_exists", "ligand_or_probe_bound", "validated_inhibitor", "clinical_stage"]:
            assert set(facets[col]).issubset(valid)

    def test_kdm1a_has_validated_inhibitor_and_clinical_stage_pharmacology(self):
        facets = pfd2.load_f3_structural_facets().set_index("gene")
        assert facets.loc["KDM1A", "validated_inhibitor"] == "YES"
        assert facets.loc["KDM1A", "clinical_stage"] == "YES"

    def test_vezf1_has_no_experimental_structure(self):
        facets = pfd2.load_f3_structural_facets().set_index("gene")
        assert facets.loc["VEZF1", "structure_exists"] == "NO"

    def test_no_focus_gene_has_a_validated_inhibitor_except_kdm1a(self):
        facets = pfd2.load_f3_structural_facets().set_index("gene")
        for gene in ["TLK2", "USP34", "VEZF1"]:
            assert facets.loc[gene, "validated_inhibitor"] == "NO"


class TestF4HumanDepmapValidation:
    def test_tcga_forest_only_covers_usp34_and_vezf1(self):
        forest = pfd2.load_f4_tcga_forest()
        assert set(forest["candidate"]) == {"USP34", "VEZF1"}

    def test_depmap_effect_covers_all_four_focus_genes_and_11_lines(self):
        d4 = pfd2.load_f4_depmap_effect()
        assert set(d4["gene"]) == set(FOCUS_FOUR)
        for gene in FOCUS_FOUR:
            assert (d4["gene"] == gene).sum() == 11

    def test_tlk2_has_the_strongest_median_dependency_among_focus_four(self):
        d4 = pfd2.load_f4_depmap_effect()
        medians = d4.groupby("gene")["chronos_effect"].median()
        assert medians.idxmin() == "TLK2"  # most negative Chronos = strongest dependency

    def test_depmap_effect_matches_post_audit_evidence_matrix_direction(self):
        d4 = pfd2.load_f4_depmap_effect()
        em = pad.build_evidence_matrix().set_index("gene")
        medians = d4.groupby("gene")["chronos_effect"].median()
        for gene in FOCUS_FOUR:
            assert medians[gene] == pytest.approx(em.loc[gene, "median_chronos_er_luminal"], abs=1e-6)


class TestF5StructuralRow:
    def test_usp34_structural_row_is_usp34(self):
        row = pfd2.load_f5_usp34_structural_row()
        assert row["gene"] == "USP34"

    def test_usp34_has_no_validated_inhibitor(self):
        row = pfd2.load_f5_usp34_structural_row()
        assert row["E_validated_selective_small_molecule_inhibitor"].startswith("NO")

    def test_structure_paths_point_to_the_two_frozen_pdb_files(self):
        paths = pfd2.usp34_structure_paths()
        assert set(paths.keys()) == {"7W3R", "7W3U"}


class TestF6RoleCompass:
    def test_role_compass_has_13_significant_sensitising_genes(self):
        compass = pfd2.load_f6_role_compass()
        assert len(compass) == 13

    def test_role_compass_focus_flag_matches_the_four(self):
        compass = pfd2.load_f6_role_compass()
        assert set(compass.loc[compass["is_focus_gene"], "gene"]) == set(FOCUS_FOUR)

    def test_neg_log10_fdr_is_positive_and_finite(self):
        compass = pfd2.load_f6_role_compass()
        assert np.isfinite(compass["neg_log10_fdr"]).all()
        assert (compass["neg_log10_fdr"] > 0).all()


class TestFiguresExistAndAreNonDegenerate:
    """Robustness checks only -- existence, non-zero dimensions, a sane
    minimum resolution, and no fully-transparent output. No exact pixel
    dimension is pinned (that pattern was deliberately retired in the
    2026-08-16 cleanup patch for tests/test_nebula_plots_final.py and must
    not be reintroduced here)."""

    @pytest.mark.parametrize("stem", MAIN_FIGURE_STEMS)
    def test_png_exists_and_has_content(self, stem):
        path = FIGURES / f"{stem}.png"
        assert path.exists(), f"{path} was not written"
        assert path.stat().st_size > 0

    @pytest.mark.parametrize("stem", MAIN_FIGURE_STEMS)
    def test_png_has_a_sane_minimum_resolution(self, stem):
        with Image.open(FIGURES / f"{stem}.png") as image:
            width, height = image.size
            assert width >= 1200
            assert height >= 800

    @pytest.mark.parametrize("stem", [s for s in MAIN_FIGURE_STEMS if s != "F5_USP34_structure_tractability"])
    def test_vector_formats_exist_for_non_structural_figures(self, stem):
        assert (FIGURES / f"{stem}.pdf").exists()
        assert (FIGURES / f"{stem}.svg").exists()

    def test_f5_is_png_only_by_documented_design(self):
        # F5's content is a rasterized ray-traced PyMOL render -- documented
        # in POSTER_FINAL_FIGURE_GUIDE.md's F5 "Format note" as PNG-only,
        # matching the same constraint already accepted for the earlier
        # poster/05_structure.png and figure-bank structural figures.
        assert (FIGURES / "F5_USP34_structure_tractability.png").exists()
        assert not (FIGURES / "F5_USP34_structure_tractability.pdf").exists()

    def test_contact_sheets_exist(self):
        assert (FIGURES / "POSTER_FINAL_CONTACT_SHEET.png").exists()
        assert (FIGURES / "RETIRED_CANDIDATES_CONTACT_SHEET.png").exists()


class TestReportsExist:
    def test_figure_guide_exists_and_mentions_all_six_figures(self):
        text = (REPORTS / "POSTER_FINAL_FIGURE_GUIDE.md").read_text()
        for label in ["F1.", "F2.", "F3.", "F4.", "F5.", "F6."]:
            assert label in text

    def test_figure_guide_has_the_external_audit_section(self):
        text = (REPORTS / "POSTER_FINAL_FIGURE_GUIDE.md").read_text()
        assert "How the external audits changed the final figure strategy" in text

    def test_layout_recommendation_exists_and_does_not_claim_assembly(self):
        text = (REPORTS / "POSTER_LAYOUT_RECOMMENDATION.md").read_text()
        assert "no poster has been assembled" in text.lower()


class TestNoScientificOverclaim:
    """Guards against the specific overclaiming failure modes the task
    explicitly warned against."""

    def test_no_composite_score_language_in_guide(self):
        text = " ".join((REPORTS / "POSTER_FINAL_FIGURE_GUIDE.md").read_text().lower().split())
        assert "composite score" not in text or "no composite" in text
        assert "weighted score" not in text or "no composite or weighted score" in text

    def test_gse240112_never_called_chronic_resistance_in_final_data_module(self):
        import inspect
        import re

        src = " ".join(inspect.getsource(pfd2).lower().split())
        for m in re.finditer(r"chronic resistance", src):
            window = src[max(0, m.start() - 60):m.start()]
            assert "never" in window or "not " in window, f"unguarded 'chronic resistance' near: {window!r}"

    def test_guide_states_no_universal_winner(self):
        text = (REPORTS / "POSTER_FINAL_FIGURE_GUIDE.md").read_text().lower()
        assert "no universal winner" in text
