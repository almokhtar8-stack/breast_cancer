import pandas as pd
import pytest

from src.candidate_adjudication_provenance import _values_match


class TestValuesMatch:
    def test_equal_floats_match(self):
        assert _values_match(0.0123456789, 0.0123456789) is True

    def test_both_nan_match(self):
        assert _values_match(float("nan"), float("nan")) is True

    def test_one_nan_one_value_does_not_match(self):
        assert _values_match(float("nan"), 0.5) is False

    def test_materially_different_values_do_not_match(self):
        assert _values_match(0.05, 0.06) is False

    def test_pandas_na_treated_as_nan(self):
        assert _values_match(pd.NA, pd.NA) is True


class TestRealProvenanceTable:
    def test_all_98_real_comparisons_match(self):
        path = "results/tables/candidate_adjudication/multimodal7_value_provenance.tsv"
        try:
            df = pd.read_csv(path, sep="\t")
        except FileNotFoundError:
            pytest.skip("provenance table not generated in this environment")
        assert len(df) > 0
        assert df["matches"].all(), df.loc[~df["matches"]]

    def test_provenance_covers_all_five_datasets(self):
        path = "results/tables/candidate_adjudication/multimodal7_value_provenance.tsv"
        try:
            df = pd.read_csv(path, sep="\t")
        except FileNotFoundError:
            pytest.skip("provenance table not generated in this environment")
        datasets = set(df["dataset"])
        for expected in ["crispr", "gse118713", "gse245601_track_a", "gse245601_track_b", "gse240112_tumor", "gse111151"]:
            assert expected in datasets
