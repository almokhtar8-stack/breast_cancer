from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from src.cross_dataset_evidence_tables import build_evidence_long, build_wide_matrix

REPO_ROOT = Path(__file__).parent.parent


def _synthetic_inputs():
    universe = pd.DataFrame({"gene": ["G1", "G2", "G3"]})
    crispr = pd.DataFrame({"symbol": ["G1", "G2"], "effect": [-1.5, 1.2], "p_value": [0.01, 0.02], "fdr": [0.03, 0.04], "n_guides": [4, 4]})
    gse118713 = pd.DataFrame({"symbol": ["G1"], "effect": [0.5], "p_value": [0.01], "fdr": [0.02]})
    gse245601_a = pd.DataFrame({"symbol": ["G1", "G2"], "effect": [0.3, -0.2], "p_value": [0.1, 0.2], "fdr": [0.3, 0.4]})
    gse245601_b = pd.DataFrame({"symbol": ["G1"], "effect": [-0.1], "p_value": [0.5], "fdr": [0.6]})
    gse240112_tumor = pd.DataFrame({"symbol": ["G2"], "effect": [0.8], "p_value": [0.02], "fdr": [0.05]})
    gse240112_epi = pd.DataFrame({"symbol": ["G2", "G3"], "effect": [0.7, 0.1], "p_value": [0.03, 0.5], "fdr": [0.06, 0.6]})
    gse111151 = pd.DataFrame({"symbol": ["G3"], "effect": [-0.4], "p_value": [0.2], "fdr": [0.5]})
    return universe, crispr, gse118713, gse245601_a, gse245601_b, gse240112_tumor, gse240112_epi, gse111151


class TestBuildWideMatrix:
    def test_all_genes_present_with_na_for_missing(self):
        args = _synthetic_inputs()
        out = build_wide_matrix(*args).set_index("gene")
        assert set(out.index) == {"G1", "G2", "G3"}
        assert pd.isna(out.loc["G3", "crispr_effect"])
        assert not out.loc["G3", "crispr_testable"]

    def test_crispr_direction_never_conflated_with_rna_sign(self):
        args = _synthetic_inputs()
        out = build_wide_matrix(*args).set_index("gene")
        assert out.loc["G1", "crispr_direction"] == "sensitising_KO"
        assert out.loc["G2", "crispr_direction"] == "tolerance_associated_KO"
        # RNA direction columns use their own vocabulary, never "sensitising"/"tolerance"
        assert "sensitising" not in str(out.loc["G1", "gse118713_direction"])

    def test_gse240112_tumor_and_epi_both_retained_as_separate_columns(self):
        args = _synthetic_inputs()
        out = build_wide_matrix(*args).set_index("gene")
        assert out.loc["G2", "gse240112_tumor_log2fc"] == pytest.approx(0.8)
        assert out.loc["G2", "gse240112_epi_log2fc"] == pytest.approx(0.7)

    def test_outlier_fragility_flag(self):
        args = _synthetic_inputs()
        out = build_wide_matrix(*args).set_index("gene")
        # G3 tested only in epithelial sensitivity track, not tumor-cell primary track
        assert bool(out.loc["G3", "gse240112_outlier_fragility"])
        assert not bool(out.loc["G2", "gse240112_outlier_fragility"])

    def test_testable_flag_requires_tumor_cell_track_not_epithelial_alone(self):
        # regression test for a real bug caught by Codex review: a gene tested ONLY in the
        # epithelial sensitivity track (no tumor-cell value at all) must NOT be marked
        # gse240112_testable=True, since it receives no gse240112 ranking percentile
        args = _synthetic_inputs()
        out = build_wide_matrix(*args).set_index("gene")
        assert pd.isna(out.loc["G3", "gse240112_tumor_log2fc"])
        assert pd.notna(out.loc["G3", "gse240112_epi_log2fc"])
        assert not bool(out.loc["G3", "gse240112_testable"])
        # G2 IS tumor-cell tested -> must remain testable=True
        assert bool(out.loc["G2", "gse240112_testable"])

    def test_gse245601_track_direction_agreement(self):
        args = _synthetic_inputs()
        out = build_wide_matrix(*args).set_index("gene")
        # G1: epi=+0.3, malignant=-0.1 -> disagree
        assert not bool(out.loc["G1", "gse245601_track_direction_agreement"])


class TestBuildEvidenceLong:
    def test_exactly_five_dataset_rows_per_gene(self):
        args = _synthetic_inputs()
        wide = build_wide_matrix(*args)
        long_df = build_evidence_long(wide)
        counts = long_df.groupby("gene").size()
        assert (counts == 5).all()

    def test_five_datasets_named_no_track_double_counted(self):
        args = _synthetic_inputs()
        wide = build_wide_matrix(*args)
        long_df = build_evidence_long(wide)
        assert set(long_df["dataset"].unique()) == {"crispr", "gse118713", "gse245601", "gse240112", "gse111151"}
        # no "track_a"/"track_b"/"tumor"/"epithelial" appearing as a distinct dataset value
        assert not any("track" in d.lower() for d in long_df["dataset"].unique())

    def test_untested_gene_marked_not_testable_not_dropped(self):
        args = _synthetic_inputs()
        wide = build_wide_matrix(*args)
        long_df = build_evidence_long(wide)
        row = long_df.loc[(long_df["gene"] == "G3") & (long_df["dataset"] == "crispr")].iloc[0]
        assert not row["testable"]
        assert pd.isna(row["effect"])


class TestRealData:
    def test_real_wide_and_long_tables_if_present(self):
        wide_path = REPO_ROOT / "results" / "tables" / "cross_dataset_genomewide" / "all_genes_cross_dataset_evidence.tsv"
        long_path = REPO_ROOT / "results" / "tables" / "cross_dataset_genomewide" / "evidence_long.tsv"
        if not wide_path.exists() or not long_path.exists():
            pytest.skip("cross-dataset evidence tables not present in this checkout")
        wide = pd.read_csv(wide_path, sep="\t")
        long_df = pd.read_csv(long_path, sep="\t")
        assert not wide["gene"].duplicated().any()
        assert len(long_df) == len(wide) * 5
        assert set(long_df["dataset"].unique()) == {"crispr", "gse118713", "gse245601", "gse240112", "gse111151"}
