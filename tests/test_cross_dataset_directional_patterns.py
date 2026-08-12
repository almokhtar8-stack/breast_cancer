from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from src.cross_dataset_directional_patterns import (
    build_directional_pattern,
    build_directional_patterns_table,
    build_significance_vs_consistency_summary,
)

REPO_ROOT = Path(__file__).parent.parent


class TestBuildDirectionalPattern:
    def test_sensitising_ko_labeled_correctly(self):
        row = pd.Series({"crispr_direction": "sensitising_KO", "gse118713_log2fc": 0.5, "gse240112_tumor_log2fc": 0.3, "gse111151_log2fc": np.nan, "gse245601_epi_log2fc": -0.2})
        out = build_directional_pattern(row)
        assert "CRISPR sensitising" in out
        assert "↑" in out
        assert "·" in out  # untestable gse111151
        assert "↓" in out

    def test_never_says_rna_sensitising_or_tolerance(self):
        # RNA direction vocabulary must never borrow CRISPR's sensitising/tolerance words
        row = pd.Series({"crispr_direction": "sensitising_KO", "gse118713_log2fc": 0.5, "gse240112_tumor_log2fc": 0.3, "gse111151_log2fc": 0.1, "gse245601_epi_log2fc": 0.2})
        out = build_directional_pattern(row)
        rna_part = out.split("|", 1)[1]
        assert "sensitising" not in rna_part.lower()
        assert "tolerance" not in rna_part.lower()

    def test_untested_crispr_labeled_explicitly(self):
        row = pd.Series({"crispr_direction": "not_applicable", "gse118713_log2fc": np.nan, "gse240112_tumor_log2fc": np.nan, "gse111151_log2fc": np.nan, "gse245601_epi_log2fc": np.nan})
        out = build_directional_pattern(row)
        assert "CRISPR untested" in out


class TestBuildDirectionalPatternsTable:
    def test_one_row_per_gene(self):
        df = pd.DataFrame(
            {
                "gene": ["A", "B"], "crispr_direction": ["sensitising_KO", "not_applicable"],
                "gse118713_log2fc": [0.5, np.nan], "gse240112_tumor_log2fc": [0.3, np.nan],
                "gse111151_log2fc": [0.1, np.nan], "gse245601_epi_log2fc": [-0.1, np.nan],
            }
        )
        out = build_directional_patterns_table(df)
        assert len(out) == 2
        assert "pattern" in out.columns


class TestBuildSignificanceVsConsistencySummary:
    def test_merges_resistance_consensus(self):
        df = pd.DataFrame(
            {"gene": ["A"], "n_datasets_fdr05": [2], "n_datasets_nominal_p05": [3], "median_evidence_percentile": [0.8],
             "gse245601_track_direction_agreement": [True], "gse240112_track_direction_agreement": [False]}
        )
        resistance = pd.DataFrame({"gene": ["A"], "resistance_direction_consensus": ["all_up"], "resistance_fdr05_count": [2]})
        out = build_significance_vs_consistency_summary(df, resistance, None)
        assert out.loc[0, "resistance_direction_consensus"] == "all_up"

    def test_stability_optional(self):
        df = pd.DataFrame(
            {"gene": ["A"], "n_datasets_fdr05": [1], "n_datasets_nominal_p05": [1], "median_evidence_percentile": [0.5],
             "gse245601_track_direction_agreement": [None], "gse240112_track_direction_agreement": [None]}
        )
        resistance = pd.DataFrame({"gene": ["A"], "resistance_direction_consensus": ["mixed"], "resistance_fdr05_count": [0]})
        out = build_significance_vs_consistency_summary(df, resistance, None)
        assert "stability_label" not in out.columns


class TestRealData:
    def test_real_directional_patterns_if_present(self):
        path = REPO_ROOT / "results" / "tables" / "cross_dataset_genomewide" / "directional_patterns.tsv"
        if not path.exists():
            pytest.skip("directional patterns table not present in this checkout")
        out = pd.read_csv(path, sep="\t")
        assert not out["gene"].duplicated().any()
        assert out["pattern"].str.contains("CRISPR").all()
