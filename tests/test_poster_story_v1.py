"""Tests for the POSTER-STORY-V1 figure bank (results/figures/poster_story_v1/),
a complete visual-story reset around one central cross-context heatmap.

Pins the v1 data-transform functions (hero heatmap pairs, malignant-vs-
non-malignant paired delta, 13-hit log2FC matrix) against the underlying
frozen tables, checks the honest handling of untestable genes (no
fabricated "+nan" display), and checks that every rendered figure is a
real, non-degenerate output. No pixel dimension is pinned exactly.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from PIL import Image

from src import poster_story_v1_data as sv1
from src import poster_exploration_v3_data as pv3
from src import post_audit_sensitivity_data as pad

FIGURES = Path("results/figures/poster_story_v1")
REPORTS = Path("results/reports/poster_story_v1")

FOCUS_FOUR = ["KDM1A", "TLK2", "USP34", "VEZF1"]

HEATMAP_CANDIDATE_STEMS = [
    "CANDIDATE1_raw_paired", "CANDIDATE2_delta_only", "CANDIDATE3_paired_plus_delta_annotation",
    "CANDIDATE4_hybrid_hero", "CANDIDATE5_hybrid_13hit",
]
MAIN_STEMS = [
    "FIG1_crispr_discovery", "HERO_cross_context_heatmap", "FIG2_structural_comparison",
    "FIG3_pathway_convergence", "FIG4_disease_clinical_context", "FIG5_depmap_context",
    "FIG6_final_synthesis",
]
BACKUP_STEMS = ["BACKUP_network_usp34"]
ALL_FIGURE_STEMS = HEATMAP_CANDIDATE_STEMS + MAIN_STEMS + BACKUP_STEMS
PNG_ONLY_STEMS = ["FIG2_structural_comparison", "FIG5_depmap_context", "BACKUP_network_usp34"]
VECTOR_STEMS = [s for s in ALL_FIGURE_STEMS if s not in PNG_ONLY_STEMS]


class TestPaletteConsistency:
    def test_focus_four_matches_frozen_modules(self):
        assert sv1.FOCUS_FOUR == pv3.FOCUS_FOUR == pad.FOCUS_FOUR == FOCUS_FOUR

    def test_colors_identical_to_prior_phases(self):
        assert sv1.FOCUS_COLORS == pv3.FOCUS_COLORS


class TestHeroHeatmapData:
    def test_hero_pairs_has_16_rows(self):
        pairs = sv1.build_hero_heatmap_pairs()
        assert len(pairs) == 16  # 4 datasets x 4 genes

    def test_hero_pairs_covers_all_4_datasets(self):
        pairs = sv1.build_hero_heatmap_pairs()
        assert set(pairs["dataset"]) == set(sv1.DATASET_ORDER)
        assert len(sv1.DATASET_ORDER) == 4

    def test_gse118713_log2fc_matches_manual_recomputation(self):
        """Independent recomputation from the raw per-sample loader, not
        using the same code path as build_hero_heatmap_pairs."""
        pairs = sv1.build_hero_heatmap_pairs()
        row = pairs[(pairs["dataset"] == "GSE118713") & (pairs["gene"] == "USP34")].iloc[0]

        raw = pv3.load_gse118713_focus_gene_samples if hasattr(pv3, "load_gse118713_focus_gene_samples") else None
        from src import poster_exploration_v2_data as pv2
        raw_df = pv2.load_gse118713_focus_gene_samples()
        sub = raw_df[raw_df["gene_symbol"] == "USP34"]
        mcf7 = sub.loc[sub["condition"] == "MCF7", "tpm"].mean()
        tamr = sub.loc[sub["condition"] == "TAMR", "tpm"].mean()
        expected = float(np.log2(tamr + 1) - np.log2(mcf7 + 1))
        assert row["log2fc"] == pytest.approx(expected, abs=1e-9)

    def test_gse245601_log2fc_matches_frozen_master_table_sign(self):
        """The acute-context delta's sign should agree with the already-
        frozen gse245601_epi_log2fc column in the master cross-dataset
        table (a real, independent cross-check, not a re-derivation)."""
        pairs = sv1.build_hero_heatmap_pairs()
        cross = sv1.load_cross_dataset_raw().set_index("gene")
        for gene in FOCUS_FOUR:
            row = pairs[(pairs["dataset"] == "GSE245601") & (pairs["gene"] == gene)].iloc[0]
            frozen_sign = np.sign(cross.loc[gene, "gse245601_epi_log2fc"])
            assert np.sign(row["log2fc"]) == frozen_sign or abs(row["log2fc"]) < 0.01


class Test13HitMatrix:
    def test_13hit_matrix_has_13_genes_4_datasets(self):
        m = sv1.load_13hit_log2fc_matrix()
        assert len(m) == 13
        assert set(m.columns) == {"gene", "GSE118713", "GSE111151", "GSE240112", "GSE245601"}

    def test_usp17l29_is_nan_not_zero(self):
        """USP17L29 is untestable in GSE111151/GSE240112/GSE245601 (no
        count-matrix entry) -- must be a real NaN, never silently zeroed."""
        m = sv1.load_13hit_log2fc_matrix().set_index("gene")
        assert pd.isna(m.loc["USP17L29", "GSE111151"])
        assert pd.isna(m.loc["USP17L29", "GSE240112"])
        assert pd.isna(m.loc["USP17L29", "GSE245601"])

    def test_kdm1a_and_tlk2_are_testable_in_all_4(self):
        m = sv1.load_13hit_log2fc_matrix().set_index("gene")
        for gene in ["KDM1A", "TLK2"]:
            assert m.loc[gene].notna().all()


class TestMalignantVsNonMalignant:
    def test_paired_delta_matches_frozen_value(self):
        """Re-verifies the same formula-consistency check already
        established in poster_story_v1_data.py against the independent
        frozen `malignant_vs_nonmalignant_candidates.tsv` mean delta."""
        computed = sv1.build_malignant_vs_nonmalignant_paired_delta()
        frozen = sv1.load_malignant_vs_nonmalignant_frozen_delta().set_index("gene")
        for gene in FOCUS_FOUR:
            sub = computed[computed["gene"] == gene]
            assert len(sub) == 5  # 5 real patients
            mean_delta = sub["delta_log2"].mean()
            assert mean_delta == pytest.approx(frozen.loc[gene, "mean_delta_malignant_minus_nonmalignant"], abs=1e-6)

    def test_all_4_focus_genes_present(self):
        computed = sv1.build_malignant_vs_nonmalignant_paired_delta()
        assert set(computed["gene"]) == set(FOCUS_FOUR)

    def test_copykat_labels_file_exists_and_is_real(self):
        path = Path("results/tables/gse245601_copykat_sensitivity_labels.tsv")
        assert path.exists()
        df = pd.read_csv(path, sep="\t", nrows=5)
        assert "sensitivity_malignancy_label" in df.columns


class TestNetworkBackupHonesty:
    def test_kdm1a_tlk2_have_zero_network_rows(self):
        from src import poster_exploration_v2_data as pv2
        nb = pv2.load_direct_neighbors()
        assert "KDM1A" not in set(nb["candidate"])
        assert "TLK2" not in set(nb["candidate"])

    def test_usp34_network_backup_uses_real_neighbor_count(self):
        from src import poster_exploration_v2_data as pv2
        nb = pv2.load_direct_neighbors("USP34")
        assert len(nb) == 10


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
            assert width >= 1200
            assert height >= 800

    @pytest.mark.parametrize("stem", VECTOR_STEMS)
    def test_vector_formats_exist(self, stem):
        assert (FIGURES / f"{stem}.pdf").exists()
        assert (FIGURES / f"{stem}.svg").exists()

    def test_contact_sheets_exist(self):
        assert (FIGURES / "CONTACT_HEATMAP_CANDIDATES.png").exists()
        assert (FIGURES / "CONTACT_MAIN_STORY.png").exists()


class TestNoHandTypedNumbers:
    def test_no_literal_19103_string_in_module(self):
        import inspect

        from src import poster_story_v1_visualization as sv1v

        src = inspect.getsource(sv1v)
        assert '"19,103' not in src
        assert "'19,103" not in src

    def test_heatmap_cell_text_computed_not_literal(self):
        """The 'n/t' handling for untestable cells must be triggered by a
        real np.isfinite check, not a hardcoded gene-name exclusion list."""
        import inspect

        from src import poster_story_v1_visualization as sv1v

        src = inspect.getsource(sv1v.build_candidate5_hybrid_13hit) + inspect.getsource(sv1v._build_hybrid_heatmap)
        assert "USP17L29" not in src
        assert "isfinite" in src


class TestReportsExist:
    def test_story_plan_exists(self):
        assert (REPORTS / "STORY_PLAN.md").exists()

    def test_data_audit_exists(self):
        assert (REPORTS / "DATA_AUDIT.md").exists()

    def test_final_recommendation_exists(self):
        path = REPORTS / "FINAL_FIGURE_RECOMMENDATION.md"
        assert path.exists()
        text = path.read_text()
        for stem in MAIN_STEMS:
            assert stem in text, f"{stem} missing from FINAL_FIGURE_RECOMMENDATION.md"
