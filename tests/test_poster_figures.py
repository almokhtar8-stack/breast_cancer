"""Targeted tests for the science-first poster figure rebuild.

These tests protect the specific numbers the user pinned for the poster
(Hany CRISPR effect/FDR and DepMap ER+/luminal dependency, for all four
frozen candidates) against silent drift, verify every real-data loader
returns genuine per-sample/per-cell-line data (not a symbol/icon summary),
and check that the six poster figures are real, non-trivial renders.
"""

from __future__ import annotations

import shutil
from pathlib import Path

import pandas as pd
import pytest

from src import poster_figures_data as pfd

FIGURES = Path("results/figures/poster")
REPORTS = Path("results/reports/poster")
GENE_ORDER = ["USP34", "VEZF1", "EML5", "CITED2"]

EXPECTED_HANY = {
    "USP34": dict(effect=-1.391298, fdr=0.041685),
    "VEZF1": dict(effect=-1.602445, fdr=0.037258),
    "EML5": dict(effect=-1.058423, fdr=0.148773),
    "CITED2": dict(effect=-1.495356, fdr=0.109955),
}

EXPECTED_DEPMAP_ER_LUMINAL_PCT = {
    "USP34": 0.0,
    "VEZF1": 27.3,
    "EML5": 0.0,
    "CITED2": 0.0,
}


class TestPinnedHanyValues:
    @pytest.mark.parametrize("gene", GENE_ORDER)
    def test_hany_effect_and_fdr_match_pinned_values(self, gene):
        df = pfd.load_shortlist_freeze().set_index("gene")
        assert df.loc[gene, "crispr_effect"] == pytest.approx(EXPECTED_HANY[gene]["effect"], abs=1e-6)
        assert df.loc[gene, "crispr_fdr"] == pytest.approx(EXPECTED_HANY[gene]["fdr"], abs=1e-6)

    @pytest.mark.parametrize("gene", GENE_ORDER)
    def test_genomewide_table_reproduces_the_same_pinned_values(self, gene):
        # the genome-wide discovery plot (Figure 1) reads a DIFFERENT file
        # (data/processed/labels.parquet) from the frozen shortlist table --
        # both must agree exactly, since they describe the same screen
        df = pfd.load_genomewide_crispr()
        row = df[df["gene"] == gene].iloc[0]
        assert row["effect_size"] == pytest.approx(EXPECTED_HANY[gene]["effect"], abs=1e-6)
        assert row["fdr"] == pytest.approx(EXPECTED_HANY[gene]["fdr"], abs=1e-6)

    def test_genomewide_table_is_the_real_full_screen_not_a_subset(self):
        df = pfd.load_genomewide_crispr()
        assert len(df) == 19103

    def test_no_fdr_value_is_shared_between_two_candidates(self):
        fdrs = [v["fdr"] for v in EXPECTED_HANY.values()]
        assert len(set(fdrs)) == len(fdrs)


class TestPinnedDepMapValues:
    @pytest.mark.parametrize("gene", GENE_ORDER)
    def test_er_luminal_strong_dependency_matches_pinned_value(self, gene):
        dep = pfd.load_depmap_dependency().set_index("candidate")
        pct = dep.loc[gene, "frac_strongly_dependent_er_luminal"] * 100
        assert pct == pytest.approx(EXPECTED_DEPMAP_ER_LUMINAL_PCT[gene], abs=0.05)

    def test_denominator_is_eleven_screened_lines_for_every_candidate(self):
        dep = pfd.load_depmap_dependency().set_index("candidate")
        assert (dep["n_er_luminal"] == 11).all()

    def test_depmap_table_is_current_26q1_never_archived_24q4(self):
        dep = pfd.load_depmap_dependency()
        assert (dep["depmap_release"] == "26Q1").all()

    def test_real_per_line_chronos_values_match_the_11_line_denominator(self):
        long = pfd.load_depmap_gene_effect_er_luminal()
        assert long["cell_line"].nunique() == 11
        assert set(long["gene"]) == set(GENE_ORDER)
        assert len(long) == 11 * len(GENE_ORDER)

    def test_real_per_line_chronos_medians_are_consistent_with_frozen_summary(self):
        # Figure 3 plots the REAL per-line points -- their median must
        # reproduce the already-frozen summary table's median exactly
        long = pfd.load_depmap_gene_effect_er_luminal()
        dep = pfd.load_depmap_dependency().set_index("candidate")
        for gene in GENE_ORDER:
            vals = long.loc[long["gene"] == gene, "chronos_effect"]
            assert vals.median() == pytest.approx(dep.loc[gene, "median_gene_effect_er_luminal"], abs=1e-9)


class TestGSE118713RealSampleData:
    def test_nine_real_samples_across_three_conditions(self):
        samples = pfd.load_gse118713_usp34_samples()
        assert len(samples) == 9
        assert set(samples["condition"]) == {"MCF7", "TAMR", "FASR"}
        assert (samples.groupby("condition").size() == 3).all()

    def test_tamr_vs_mcf7_stats_match_frozen_cross_dataset_table(self):
        stats = pfd.load_gse118713_usp34_stats()
        row = stats[stats["contrast"] == "TAMR_vs_MCF7"].iloc[0]
        cross = pfd.load_cross_dataset_evidence().set_index("gene")
        assert row["log2fc"] == pytest.approx(cross.loc["USP34", "gse118713_log2fc"], abs=1e-6)
        assert row["fdr"] == pytest.approx(cross.loc["USP34", "gse118713_fdr"], abs=1e-6)

    def test_tamr_samples_have_higher_mean_tpm_than_mcf7(self):
        # sanity check on the real data direction shown in Figure 2 Panel A
        samples = pfd.load_gse118713_usp34_samples()
        tamr_mean = samples.loc[samples["condition"] == "TAMR", "tpm"].mean()
        mcf7_mean = samples.loc[samples["condition"] == "MCF7", "tpm"].mean()
        assert tamr_mean > mcf7_mean


class TestGSE240112RealSampleData:
    def test_six_real_tumor_samples(self):
        samples = pfd.load_gse240112_vezf1_samples()
        assert len(samples) == 6
        assert set(samples["group"]) == {"PT", "RT"}
        assert (samples.groupby("group").size() == 3).all()

    def test_stats_match_frozen_cross_dataset_table(self):
        stats = pfd.load_gse240112_vezf1_stats()
        cross = pfd.load_cross_dataset_evidence().set_index("gene")
        assert stats["log2fc"] == pytest.approx(cross.loc["VEZF1", "gse240112_tumor_log2fc"], abs=1e-6)
        assert stats["genomewide_fdr"] == pytest.approx(cross.loc["VEZF1", "gse240112_tumor_fdr"], abs=1e-6)

    def test_recurrent_samples_have_higher_mean_log2cpm_than_primary(self):
        samples = pfd.load_gse240112_vezf1_samples()
        rt_mean = samples.loc[samples["group"] == "RT", "log2cpm"].mean()
        pt_mean = samples.loc[samples["group"] == "PT", "log2cpm"].mean()
        assert rt_mean > pt_mean


class TestForestTable:
    def test_covers_all_four_genes_and_four_datasets(self):
        forest = pfd.build_forest_table()
        assert set(forest["gene"]) == set(GENE_ORDER)
        assert set(forest["dataset_key"]) == set(pfd.DATASET_LABELS.keys())
        assert len(forest) == len(GENE_ORDER) * len(pfd.DATASET_LABELS)

    def test_gse245601_column_is_labeled_acute_not_resistance(self):
        assert "acute" in pfd.DATASET_LABELS["gse245601_epi"].lower()
        assert "resistan" not in pfd.DATASET_LABELS["gse245601_epi"].lower()

    def test_gse240112_column_is_not_labeled_causal_resistance(self):
        assert "resistan" not in pfd.DATASET_LABELS["gse240112_tumor"].lower()


class TestGDSCPharmacogenomics:
    def test_real_per_cell_line_data_44_breast_lines(self):
        df = pfd.load_gdsc_usp34_azd7762()
        assert len(df) == 44
        assert df["SANGER_MODEL_ID"].nunique() == 44

    def test_association_is_negative_consistent_with_frozen_finding(self):
        import scipy.stats as ss
        df = pfd.load_gdsc_usp34_azd7762()
        rho, p = ss.spearmanr(df["USP34"], df["LN_IC50"])
        assert rho < 0
        assert p < 0.05

    def test_er_luminal_flag_is_a_real_subset_not_all_or_none(self):
        df = pfd.load_gdsc_usp34_azd7762()
        assert 0 < df["is_er_luminal"].sum() < len(df)

    def test_stats_are_read_from_the_frozen_table_not_hand_typed(self):
        # regression test: an earlier draft hardcoded "FDR=0.008" as a
        # literal string in the figure title instead of loading it -- this
        # loader must return the real, already-computed value
        row = pfd.load_gdsc_usp34_azd7762_stats()
        assert row["drug_id"] == 1402
        assert row["dataset"] == "GDSC1"
        assert row["metric"] == "LN_IC50"
        assert row["fdr"] == pytest.approx(0.0078, abs=2e-3)
        assert row["spearman_rho"] < 0

    def test_visualization_module_does_not_hardcode_the_fdr_string(self):
        text = Path("src/poster_figures_visualization.py").read_text()
        assert "FDR=0.008" not in text


class TestStructureData:
    def test_both_pdb_files_exist_and_are_real_structures(self):
        paths = pfd.usp34_structure_paths()
        for pdb_id, path in paths.items():
            assert path.exists(), f"{pdb_id} PDB file missing -- run scripts/download/download_usp34_structures.py"
            assert path.stat().st_size > 100_000, f"{pdb_id} PDB file suspiciously small"

    def test_catalytic_residues_are_real_cys_and_his_in_both_structures(self):
        from Bio.PDB import PDBParser
        import warnings
        warnings.filterwarnings("ignore")
        parser = PDBParser(QUIET=True)
        paths = pfd.usp34_structure_paths()
        for pdb_id, path in paths.items():
            structure = parser.get_structure(pdb_id, str(path))
            chain_a = structure[0]["A"]
            assert chain_a[(" ", 1903, " ")].resname == "CYS"
            assert chain_a[(" ", 2164, " ")].resname == "HIS"

    def test_covalent_probe_residue_present_only_in_bound_structure(self):
        from Bio.PDB import PDBParser
        import warnings
        warnings.filterwarnings("ignore")
        parser = PDBParser(QUIET=True)
        paths = pfd.usp34_structure_paths()
        apo = parser.get_structure("7W3R", str(paths["7W3R"]))
        bound = parser.get_structure("7W3U", str(paths["7W3U"]))
        apo_resnames = {r.resname for r in apo[0]["A"] if r.id[0] != " "}
        bound_resnames = {r.resname for r in bound[0]["A"] if r.id[0] != " "}
        assert "AYE" not in apo_resnames
        assert "AYE" in bound_resnames

    def test_pymol_is_available_for_structural_rendering(self):
        assert shutil.which("pymol") is not None, "pymol-open-source not on PATH -- required for Figure 5"


class TestNoHardcodedDuplication:
    def test_data_module_never_hardcodes_the_pinned_numbers(self):
        text = Path("src/poster_figures_data.py").read_text()
        for gene, vals in EXPECTED_HANY.items():
            assert str(vals["effect"]) not in text, f"{gene} effect hardcoded in poster_figures_data.py"
            assert str(vals["fdr"]) not in text, f"{gene} fdr hardcoded in poster_figures_data.py"

    def test_visualization_module_never_hardcodes_the_pinned_numbers(self):
        text = Path("src/poster_figures_visualization.py").read_text()
        for gene, vals in EXPECTED_HANY.items():
            assert str(vals["effect"]) not in text, f"{gene} effect hardcoded in poster_figures_visualization.py"
            assert str(vals["fdr"]) not in text, f"{gene} fdr hardcoded in poster_figures_visualization.py"


class TestFiguresExist:
    @pytest.mark.parametrize("name", [
        "01_crispr_discovery.png",
        "02_expression_evidence.png",
        "03_depmap_distributions.png",
        "04_pharmacogenomics.png",
        "05_structure.png",
        "06_experimental_strategy.png",
    ])
    def test_png_exists_and_is_nontrivial(self, name):
        path = FIGURES / name
        assert path.exists()
        assert path.stat().st_size > 20_000, f"{name} suspiciously small -- likely a blank/failed render"

    @pytest.mark.parametrize("name", [
        "01_crispr_discovery.pdf",
        "02_expression_evidence.pdf",
        "03_depmap_distributions.pdf",
        "04_pharmacogenomics.pdf",
        "06_experimental_strategy.pdf",
    ])
    def test_vector_pdf_exists(self, name):
        # Figure 5 is a PyMOL-composed raster montage and is PNG-only by
        # design (a PDF re-export would not add fidelity over the ray-traced
        # raster panels).
        assert (FIGURES / name).exists()

    def test_no_figure_is_absurdly_tall_or_wide(self):
        # regression guard for the tight_layout-with-mixed-content bug
        # caught during this rebuild, which silently produced a broken
        # multi-thousand-pixel image
        from PIL import Image
        for path in FIGURES.glob("*.png"):
            w, h = Image.open(path).size
            aspect = max(w, h) / min(w, h)
            assert aspect < 4.0, f"{path.name} has an extreme aspect ratio ({w}x{h}) -- likely a layout bug"


class TestPosterGuide:
    def test_guide_exists_and_covers_all_six_figures(self):
        path = REPORTS / "poster_figure_guide.md"
        assert path.exists()
        text = path.read_text()
        for name in [
            "01_crispr_discovery.png", "02_expression_evidence.png",
            "03_depmap_distributions.png", "04_pharmacogenomics.png",
            "05_structure.png", "06_experimental_strategy.png",
        ]:
            assert name in text

    def test_guide_documents_what_changed_from_the_prior_version(self):
        text = REPORTS.joinpath("poster_figure_guide.md").read_text()
        assert "changed" in text.lower() or "previous version" in text.lower() or "rebuilt" in text.lower()


if __name__ == "__main__":
    import sys

    sys.exit(pytest.main([__file__, "-v"]))
