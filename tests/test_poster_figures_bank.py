"""Targeted tests for the poster FIGURE BANK phase (candidate figures).

Pins critical values from the new src/poster_figures_bank_data.py loaders
against already-frozen source tables, verifies every loader returns real
per-sample/per-line/per-pathway data (not a summary standing in for it),
and checks that every candidate figure in results/figures/poster_candidates/
is a real, non-trivial render.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from src import poster_figures_bank_data as pbd
from src import poster_figures_data as pfd

FIGURES = Path("results/figures/poster_candidates")
REPORTS = Path("results/reports/poster")
GENE_ORDER = ["USP34", "VEZF1", "EML5", "CITED2"]


class TestCrisprGate1:
    def test_returns_exactly_28_gate1_hits(self):
        df = pbd.load_crispr_gate1()
        assert len(df) == 28

    # Only USP34 and VEZF1 actually clear the Gate-1 FDR<0.1 threshold among
    # the 4 frozen candidates (EML5 FDR=0.149, CITED2 FDR=0.110) -- both are
    # real, already-pinned facts (see EXPECTED_HANY in test_poster_figures.py).
    @pytest.mark.parametrize("gene", ["USP34", "VEZF1"])
    def test_gate1_table_agrees_with_the_frozen_shortlist_for_the_2_gate1_candidates(self, gene):
        gate1 = pbd.load_crispr_gate1().set_index("gene_symbol")
        freeze = pfd.load_shortlist_freeze().set_index("gene")
        assert gate1.loc[gene, "crispr_effect_size"] == pytest.approx(freeze.loc[gene, "crispr_effect"], abs=1e-6)
        assert gate1.loc[gene, "crispr_fdr"] == pytest.approx(freeze.loc[gene, "crispr_fdr"], abs=1e-6)

    def test_eml5_and_cited2_are_not_gate1_hits(self):
        # a real, verified fact -- not a bug: both FDRs are >0.1
        gate1 = pbd.load_crispr_gate1()
        assert "EML5" not in set(gate1["gene_symbol"])
        assert "CITED2" not in set(gate1["gene_symbol"])


class TestGSE118713Bank:
    def test_pca_has_9_real_samples(self):
        pca = pbd.load_gse118713_pca()
        assert len(pca) == 9
        assert set(pca["group"]) == {"MCF7", "TAMR", "FASR"}

    def test_pca_variance_fractions_are_a_real_percentage_not_placeholder(self):
        pca = pbd.load_gse118713_pca()
        v1 = pca["pc1_variance_explained_pct"].iloc[0]
        v2 = pca["pc2_variance_explained_pct"].iloc[0]
        assert 0 < v1 <= 100
        assert 0 < v2 <= 100
        assert v1 + v2 <= 100

    def test_volcano_is_genomewide_not_a_candidate_only_subset(self):
        volc = pbd.load_gse118713_volcano()
        assert len(volc) > 10000

    def test_volcano_usp34_row_agrees_with_the_frozen_de_table(self):
        volc = pbd.load_gse118713_volcano()
        row = volc[volc["gene_symbol"] == "USP34"].iloc[0]
        stats = pfd.load_gse118713_usp34_stats()
        tamr_vs_mcf7 = stats[stats["contrast"] == "TAMR_vs_MCF7"].iloc[0]
        assert row["log2fc"] == pytest.approx(tamr_vs_mcf7["log2fc"], abs=1e-6)
        assert row["fdr"] == pytest.approx(tamr_vs_mcf7["fdr"], abs=1e-6)


class TestGSE240112Bank:
    def test_volcano_is_genomewide(self):
        volc = pbd.load_gse240112_volcano()
        assert len(volc) > 10000
        assert "is_candidate" in volc.columns
        assert volc["is_candidate"].sum() == 4

    def test_vezf1_row_agrees_with_the_frozen_candidate_table(self):
        volc = pbd.load_gse240112_volcano()
        row = volc[volc["gene"] == "VEZF1"].iloc[0]
        stats = pfd.load_gse240112_vezf1_stats()
        assert row["log2fc"] == pytest.approx(stats["log2fc"], abs=1e-6)
        assert row["fdr"] == pytest.approx(stats["genomewide_fdr"], abs=1e-6)


class TestPathwayLandscape:
    def test_returns_10_pathways_x_4_datasets(self):
        pl = pbd.load_pathway_landscape()
        assert pl["pathway_label"].nunique() == 10
        assert pl["dataset_key"].nunique() == 4

    def test_estrogen_response_early_nes_matches_gsea_source_table(self):
        pl = pbd.load_pathway_landscape()
        row = pl[(pl["pathway_label"] == "Estrogen response (early)") & (pl["dataset_key"] == "gse118713")].iloc[0]
        gsea = pd.read_csv("results/tables/systems_network/gsea_gse118713.tsv", sep="\t")
        source = gsea[(gsea["collection"] == "hallmark") & (gsea["pathway"] == "HALLMARK_ESTROGEN_RESPONSE_EARLY")].iloc[0]
        assert row["NES"] == pytest.approx(source["NES"], abs=1e-9)
        assert row["fdr"] == pytest.approx(source["fdr"], abs=1e-9)

    def test_gse245601_acute_context_is_a_distinct_dataset_key(self):
        assert "gse245601" in pbd.PATHWAY_RNA_DATASETS
        assert "acute" in pbd.PATHWAY_RNA_DATASETS["gse245601"].lower()

    def test_crispr_pathway_context_is_a_separate_loader_not_merged_with_rna(self):
        rna = pbd.load_pathway_landscape()
        crispr = pbd.load_pathway_landscape_crispr()
        assert "NES" in crispr.columns
        assert set(crispr.columns) != set(rna.columns) or "dataset_key" not in crispr.columns


class TestCandidateStrongConsensusPathways:
    def test_membership_counts_match_the_audited_values(self):
        df = pbd.load_candidate_strong_consensus_pathways()
        counts = df.groupby("candidate").size()
        assert counts.get("USP34", 0) == 4
        assert counts.get("VEZF1", 0) == 2
        assert counts.get("CITED2", 0) == 49

    def test_eml5_has_zero_strong_consensus_pathways(self):
        df = pbd.load_candidate_strong_consensus_pathways()
        assert "EML5" not in set(df["candidate"])


class TestTcgaBank:
    def test_forest_table_has_exactly_two_contrasts_per_candidate(self):
        forest = pbd.load_tcga_expression_forest()
        assert len(forest) == 2 * len(GENE_ORDER)
        assert set(forest["comparison"]) == set(pbd.TCGA_FOREST_CONTRASTS)

    def test_usp34_paired_tumor_vs_normal_matches_frozen_table(self):
        forest = pbd.load_tcga_expression_forest()
        row = forest[(forest["candidate"] == "USP34") & (forest["comparison"] == "tumor_vs_normal_PAIRED")].iloc[0]
        source = pd.read_csv("results/tables/independent_validation/TCGA_candidate_expression.tsv", sep="\t")
        expected = source[(source["candidate"] == "USP34") & (source["comparison"] == "tumor_vs_normal_PAIRED")].iloc[0]
        assert row["mean_diff"] == pytest.approx(expected["mean_diff"], abs=1e-9)
        assert row["n_a"] == expected["n_a"] == 113

    def test_clinical_er_adjusted_returns_one_row_per_candidate(self):
        clinical = pbd.load_tcga_clinical_er_adjusted()
        assert len(clinical) == len(GENE_ORDER)
        assert (clinical["cohort"] == "ER_positive").all()
        assert (clinical["model"] == "adjusted_age_stage").all()


class TestDepMapGate1Summary:
    def test_returns_28_rows_one_per_gate1_gene(self):
        df = pbd.load_depmap_gate1_dependency_summary()
        assert len(df) == 28

    # EML5 and CITED2 are not Gate-1 hits (see TestCrisprGate1), so they are
    # not rows in this Gate-1-only DepMap summary -- only USP34/VEZF1 overlap
    # with the frozen 4-candidate dependency table.
    @pytest.mark.parametrize("gene", ["USP34", "VEZF1"])
    def test_usp34_and_vezf1_match_the_frozen_4_candidate_dependency_table(self, gene):
        df = pbd.load_depmap_gate1_dependency_summary().set_index("gene")
        frozen = pfd.load_depmap_dependency().set_index("candidate")
        assert df.loc[gene, "frac_strongly_dependent_er_luminal"] == pytest.approx(
            frozen.loc[gene, "frac_strongly_dependent_er_luminal"], abs=1e-6
        )

    def test_genes_absent_from_depmap_are_flagged_not_silently_dropped(self):
        df = pbd.load_depmap_gate1_dependency_summary()
        assert (~df["in_depmap"]).any()
        missing = df[~df["in_depmap"]]
        assert missing["median_chronos_er_luminal"].isna().all()


class TestGdscTopAssociations:
    def test_table_has_both_a_significant_and_a_top_effect_tier(self):
        # the raw table intentionally includes a
        # TOP_EFFECT_SIZE_NOT_NECESSARILY_SIGNIFICANT tier alongside
        # FDR_SIGNIFICANT -- figures must filter to the significant tier
        # explicitly rather than assuming every row is significant (a real
        # bug caught during this phase's review: every VEZF1 row here is in
        # the non-significant tier).
        df = pbd.load_gdsc_top_associations()
        assert set(df["tier"]) == {"FDR_SIGNIFICANT", "TOP_EFFECT_SIZE_NOT_NECESSARILY_SIGNIFICANT"}

    def test_fdr_significant_tier_rows_are_actually_fdr_significant(self):
        df = pbd.load_gdsc_top_associations()
        sig = df[df["tier"] == "FDR_SIGNIFICANT"]
        assert len(sig) > 0
        assert (sig["fdr"] < 0.05).all()

    def test_every_vezf1_row_is_in_the_non_significant_tier(self):
        df = pbd.load_gdsc_top_associations()
        vezf1 = df[df["gene"] == "VEZF1"]
        assert len(vezf1) > 0
        assert (vezf1["tier"] == "TOP_EFFECT_SIZE_NOT_NECESSARILY_SIGNIFICANT").all()

    def test_both_candidates_with_gdsc_hits_are_present(self):
        df = pbd.load_gdsc_top_associations()
        assert set(df["gene"]).issubset({"USP34", "VEZF1"})
        assert "USP34" in set(df["gene"])

    def test_figure_11_panel_a_only_plots_the_significant_tier(self):
        src = Path("src/poster_figures_bank_visualization.py").read_text()
        assert 'tier"] == "FDR_SIGNIFICANT"' in src


class TestTissueLiabilityAndContext:
    def test_liability_table_covers_only_usp34_and_vezf1(self):
        df = pbd.load_tissue_liability()
        assert set(df["candidate"]) == {"USP34", "VEZF1"}

    def test_normal_tissue_context_covers_all_4_candidates(self):
        df = pbd.load_normal_tissue_context()
        assert set(df["candidate"]) == set(GENE_ORDER)


class TestNoTeadOnePromotion:
    def test_tead1_does_not_appear_anywhere_in_bank_data_module(self):
        src = Path("src/poster_figures_bank_data.py").read_text()
        assert "TEAD1" not in src

    def test_tead1_does_not_appear_anywhere_in_bank_visualization_module(self):
        src = Path("src/poster_figures_bank_visualization.py").read_text()
        assert "TEAD1" not in src


class TestFiguresExist:
    EXPECTED_STEMS = [stem for stem, _ in __import__(
        "src.poster_figures_bank_visualization", fromlist=["CONTACT_SHEET_ITEMS"]
    ).CONTACT_SHEET_ITEMS]

    @pytest.mark.parametrize("stem", EXPECTED_STEMS)
    def test_figure_png_exists_and_is_nontrivial(self, stem):
        path = FIGURES / f"{stem}.png"
        assert path.exists(), f"missing {path}"
        assert path.stat().st_size > 20_000, f"{path} looks too small to be a real render"

    def test_contact_sheet_exists(self):
        assert (FIGURES / "POSTER_FIGURE_CONTACT_SHEET.png").exists()

    def test_no_figure_is_absurdly_tall_or_wide(self):
        from PIL import Image
        for path in FIGURES.glob("*.png"):
            w, h = Image.open(path).size
            aspect = max(w, h) / min(w, h)
            assert aspect < 4.0, f"{path.name} has an extreme aspect ratio ({w}x{h}) -- likely a layout bug"

    def test_vector_formats_exist_for_conventional_plots(self):
        # structural-biology figures (12, 12b) are PNG-only by design (they
        # embed large raster PyMOL renders); every other figure should have
        # both a pdf and an svg alongside the png.
        for stem in self.EXPECTED_STEMS:
            if stem.startswith("12"):
                continue
            assert (FIGURES / f"{stem}.pdf").exists(), f"missing {stem}.pdf"
            assert (FIGURES / f"{stem}.svg").exists(), f"missing {stem}.svg"


class TestAuditAndReviewDocs:
    def test_scientific_figure_audit_exists(self):
        assert (REPORTS / "SCIENTIFIC_FIGURE_AUDIT.md").exists()

    def test_figure_bank_review_exists_and_covers_every_candidate(self):
        path = REPORTS / "FIGURE_BANK_REVIEW.md"
        assert path.exists()
        text = path.read_text()
        for stem, _ in __import__(
            "src.poster_figures_bank_visualization", fromlist=["CONTACT_SHEET_ITEMS"]
        ).CONTACT_SHEET_ITEMS:
            assert stem in text, f"{stem} not mentioned in FIGURE_BANK_REVIEW.md"
