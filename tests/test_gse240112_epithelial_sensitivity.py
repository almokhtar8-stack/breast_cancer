from pathlib import Path

import pandas as pd
import pytest

from src.gse240112_epithelial_sensitivity import build_sensitivity_comparison

REPO_ROOT = Path(__file__).parent.parent


def _tumor_table():
    return pd.DataFrame(
        {
            "gene": ["USP34", "VEZF1", "USP17L29"],
            "tested": [True, True, False],
            "log2fc": [0.4, 1.15, float("nan")],
            "candidate_set_bh_fdr": [0.47, 0.05, float("nan")],
        }
    )


def _epi_de():
    return pd.DataFrame(
        {
            "gene": ["USP34", "VEZF1"],
            "log2fc": [0.2, -0.3],
            "fdr": [0.5, 0.6],
        }
    )


class TestBuildSensitivityComparison:
    def test_direction_agreement_computed_when_both_tested(self):
        out = build_sensitivity_comparison(_tumor_table(), _epi_de(), ["USP34", "VEZF1"]).set_index("gene")
        assert out.loc["USP34", "direction_agreement"] == True  # noqa: E712 -- both positive
        assert out.loc["VEZF1", "direction_agreement"] == False  # noqa: E712 -- tumor positive, epi negative

    def test_untested_gene_reported_as_na_not_dropped(self):
        out = build_sensitivity_comparison(_tumor_table(), _epi_de(), ["USP34", "USP17L29"]).set_index("gene")
        assert not out.loc["USP17L29", "tumor_cell_tested"]
        assert pd.isna(out.loc["USP17L29", "tumor_cell_log2fc"])
        assert out.loc["USP17L29", "direction_agreement"] is None

    def test_gene_missing_from_epithelial_track_reported_as_untested_there(self):
        tumor = _tumor_table()
        epi = pd.DataFrame({"gene": ["USP34"], "log2fc": [0.2], "fdr": [0.5]})
        out = build_sensitivity_comparison(tumor, epi, ["USP34", "VEZF1"]).set_index("gene")
        assert not out.loc["VEZF1", "all_epithelial_tested"]
        assert pd.isna(out.loc["VEZF1", "all_epithelial_log2fc"])


class TestRealData:
    def test_real_sensitivity_table_if_present(self):
        path = REPO_ROOT / "results" / "tables" / "gse240112" / "tumor_vs_epithelial_sensitivity.tsv"
        if not path.exists():
            pytest.skip("GSE240112 epithelial sensitivity table not present in this checkout")
        out = pd.read_csv(path, sep="\t")
        assert len(out) == 13
        assert "USP34" in set(out["gene"])
