from pathlib import Path

import pandas as pd
import pytest

from src.gse240112_candidate_integration import build_integrated_table

REPO_ROOT = Path(__file__).parent.parent


def _sources():
    crispr_bulk = pd.DataFrame(
        {
            "gene_symbol": ["USP34", "VEZF1"],
            "crispr_effect_size": [-1.39, -1.60],
            "crispr_fdr": [0.042, 0.037],
            "tamr_vs_mcf7_log2fc": [0.59, 0.43],
            "tamr_vs_mcf7_fdr": [0.0073, 0.24],
            "evidence_class": ["PRIMARY_RESISTANCE_SUPPORT", "SECONDARY_CONTEXT_SUPPORT"],
        }
    )
    track_a = pd.DataFrame({"gene": ["USP34"], "tested": [True], "log2fc": [-0.03], "candidate_set_bh_fdr": [0.76]})
    track_b = pd.DataFrame({"gene": ["USP34"], "tested": [True], "log2fc": [-0.18], "candidate_set_bh_fdr": [0.41]})
    g240_candidates = pd.DataFrame(
        {"gene": ["USP34", "VEZF1", "USP17L29"], "tested": [True, True, False], "log2fc": [0.40, 1.15, float("nan")],
         "p_value": [0.118, 0.004, float("nan")], "candidate_set_bh_fdr": [0.473, 0.048, float("nan")]}
    )
    g240_sensitivity = pd.DataFrame(
        {"gene": ["USP34", "VEZF1"], "all_epithelial_log2fc": [0.1, -0.2], "all_epithelial_genomewide_fdr": [0.9, 0.8], "direction_agreement": [True, False]}
    )
    return crispr_bulk, track_a, track_b, g240_candidates, g240_sensitivity


class TestBuildIntegratedTable:
    def test_one_row_per_candidate_exactly(self):
        crispr_bulk, track_a, track_b, g240c, g240s = _sources()
        out = build_integrated_table(crispr_bulk, track_a, track_b, g240c, g240s, ["USP34", "VEZF1", "USP17L29"])
        assert len(out) == 3
        assert list(out["gene"]) == ["USP34", "VEZF1", "USP17L29"]

    def test_gene_absent_from_a_source_reported_as_na_not_dropped(self):
        crispr_bulk, track_a, track_b, g240c, g240s = _sources()
        out = build_integrated_table(crispr_bulk, track_a, track_b, g240c, g240s, ["USP34", "USP17L29"]).set_index("gene")
        assert pd.isna(out.loc["USP17L29", "crispr_effect_size"])
        assert not out.loc["USP17L29", "gse240112_tumor_cell_tested"]
        assert pd.isna(out.loc["USP17L29", "gse240112_tumor_cell_log2fc"])

    def test_values_pulled_through_correctly(self):
        crispr_bulk, track_a, track_b, g240c, g240s = _sources()
        out = build_integrated_table(crispr_bulk, track_a, track_b, g240c, g240s, ["USP34"]).set_index("gene")
        assert out.loc["USP34", "crispr_effect_size"] == pytest.approx(-1.39)
        assert out.loc["USP34", "gse118713_tamr_vs_mcf7_log2fc"] == pytest.approx(0.59)
        assert out.loc["USP34", "gse245601_track_a_epithelial_log2fc"] == pytest.approx(-0.03)
        assert out.loc["USP34", "gse240112_tumor_cell_log2fc"] == pytest.approx(0.40)
        assert out.loc["USP34", "gse240112_tumor_cell_candidate_bh_fdr"] == pytest.approx(0.473)

    def test_no_composite_score_column(self):
        crispr_bulk, track_a, track_b, g240c, g240s = _sources()
        out = build_integrated_table(crispr_bulk, track_a, track_b, g240c, g240s, ["USP34"])
        forbidden = {"score", "composite_score", "weighted_score", "rank"}
        assert forbidden.isdisjoint(set(out.columns))


class TestRealData:
    def test_real_integration_table_if_present(self):
        path = REPO_ROOT / "results" / "tables" / "gse240112" / "integrated_evidence_4layer.tsv"
        if not path.exists():
            pytest.skip("GSE240112 integration table not present in this checkout")
        out = pd.read_csv(path, sep="\t")
        assert len(out) == 13
        expected = {
            "USP34", "CTDNEP1", "EIF4ENIF1", "HMGB1", "KDM1A", "PET117", "TADA2B", "VEZF1", "ICK", "SUPT4H1", "TLK2", "TSR3", "USP17L29",
        }
        assert set(out["gene"]) == expected
