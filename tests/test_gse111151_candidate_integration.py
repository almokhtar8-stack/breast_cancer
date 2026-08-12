from pathlib import Path

import pandas as pd
import pytest

from src.gse111151_candidate_integration import build_integrated_table

REPO_ROOT = Path(__file__).parent.parent


def _sources():
    crispr_bulk = pd.DataFrame(
        {
            "gene_symbol": ["USP34", "VEZF1"],
            "crispr_effect_size": [-1.39, -1.60],
            "crispr_fdr": [0.042, 0.037],
            "tamr_vs_mcf7_log2fc": [0.59, 0.43],
            "tamr_vs_mcf7_fdr": [0.0073, 0.24],
        }
    )
    track_a = pd.DataFrame({"gene": ["USP34"], "tested": [True], "log2fc": [-0.03], "candidate_set_bh_fdr": [0.76]})
    track_b = pd.DataFrame({"gene": ["USP34"], "tested": [True], "log2fc": [-0.18], "candidate_set_bh_fdr": [0.41]})
    g240_candidates = pd.DataFrame(
        {"gene": ["USP34", "VEZF1", "USP17L29"], "tested": [True, True, False], "log2fc": [0.40, 1.15, float("nan")], "candidate_set_bh_fdr": [0.473, 0.048, float("nan")]}
    )
    g111_candidates = pd.DataFrame(
        {"gene": ["USP34", "VEZF1", "USP17L29"], "tested": [True, True, False], "log2fc": [0.16, -0.24, float("nan")], "p_value": [0.21, 0.19, float("nan")], "candidate_set_bh_fdr": [0.94, 0.94, float("nan")]}
    )
    g111_classification = pd.DataFrame(
        {
            "gene": ["USP34", "VEZF1", "USP17L29"],
            "classification": ["directionally_supportive_but_weak", "neutral_no_additional_support", "untestable"],
            "n_cell_lines_consistent": [3, 2, float("nan")],
            "n_cell_lines_with_both_arms": [4, 4, float("nan")],
        }
    )
    return crispr_bulk, track_a, track_b, g240_candidates, g111_candidates, g111_classification


class TestBuildIntegratedTable:
    def test_one_row_per_candidate_exactly(self):
        crispr_bulk, track_a, track_b, g240c, g111c, g111cls = _sources()
        out = build_integrated_table(crispr_bulk, track_a, track_b, g240c, g111c, g111cls, ["USP34", "VEZF1", "USP17L29"])
        assert len(out) == 3
        assert list(out["gene"]) == ["USP34", "VEZF1", "USP17L29"]

    def test_gene_absent_from_a_source_reported_as_na_not_dropped(self):
        crispr_bulk, track_a, track_b, g240c, g111c, g111cls = _sources()
        out = build_integrated_table(crispr_bulk, track_a, track_b, g240c, g111c, g111cls, ["USP34", "USP17L29"]).set_index("gene")
        assert pd.isna(out.loc["USP17L29", "crispr_effect_size"])
        assert not out.loc["USP17L29", "gse111151_tested"]
        assert out.loc["USP17L29", "gse111151_interpretation"] == "untestable"

    def test_values_pulled_through_correctly(self):
        crispr_bulk, track_a, track_b, g240c, g111c, g111cls = _sources()
        out = build_integrated_table(crispr_bulk, track_a, track_b, g240c, g111c, g111cls, ["USP34"]).set_index("gene")
        assert out.loc["USP34", "crispr_effect_size"] == pytest.approx(-1.39)
        assert out.loc["USP34", "gse111151_log2fc"] == pytest.approx(0.16)
        assert out.loc["USP34", "gse111151_n_cell_lines_consistent"] == 3
        assert out.loc["USP34", "gse111151_interpretation"] == "directionally_supportive_but_weak"

    def test_no_composite_score_column(self):
        crispr_bulk, track_a, track_b, g240c, g111c, g111cls = _sources()
        out = build_integrated_table(crispr_bulk, track_a, track_b, g240c, g111c, g111cls, ["USP34"])
        forbidden = {"score", "composite_score", "weighted_score", "rank"}
        assert forbidden.isdisjoint(set(out.columns))


class TestRealData:
    def test_real_integration_table_if_present(self):
        path = REPO_ROOT / "results" / "tables" / "gse111151" / "integrated_evidence_5layer.tsv"
        if not path.exists():
            pytest.skip("GSE111151 integration table not present in this checkout")
        out = pd.read_csv(path, sep="\t")
        assert len(out) == 13
        expected = {
            "USP34", "CTDNEP1", "EIF4ENIF1", "HMGB1", "KDM1A", "PET117", "TADA2B", "VEZF1", "ICK", "SUPT4H1", "TLK2", "TSR3", "USP17L29",
        }
        assert set(out["gene"]) == expected
