from pathlib import Path

import pandas as pd
import pytest

from src.evidence_freeze_source_verification import _values_match

REPO_ROOT = Path(__file__).parent.parent


class TestValuesMatch:
    def test_equal_floats_match(self):
        assert _values_match(-1.3912980676838431, -1.3912980676838431) is True

    def test_both_nan_match(self):
        assert _values_match(float("nan"), float("nan")) is True

    def test_one_nan_does_not_match(self):
        assert _values_match(float("nan"), 0.5) is False

    def test_different_values_do_not_match(self):
        assert _values_match(0.04, 0.05) is False


@pytest.fixture(scope="module")
def verification():
    path = REPO_ROOT / "results" / "tables" / "evidence_freeze" / "source_value_verification.tsv"
    if not path.exists():
        pytest.skip("source-value verification table not generated in this environment")
    return pd.read_csv(path, sep="\t")


class TestRealVerificationTable:
    def test_all_comparisons_match(self, verification):
        assert len(verification) > 0
        assert verification["match"].all(), verification.loc[~verification["match"]]

    def test_covers_all_four_frozen_genes(self, verification):
        assert set(verification["gene"]) == {"USP34", "VEZF1", "EML5", "CITED2"}

    def test_covers_all_five_datasets(self, verification):
        datasets = set(verification["dataset"])
        for expected in ["crispr", "gse118713", "gse240112_tumor", "gse111151", "gse245601_track_a", "gse245601_track_b"]:
            assert expected in datasets

    def test_covers_effect_p_and_fdr_for_every_dataset(self, verification):
        metrics_per_dataset = verification.groupby("dataset")["metric"].apply(set)
        for dataset, metrics in metrics_per_dataset.items():
            assert {"log2fc", "p_value", "fdr"}.issubset(metrics) or {"effect", "p_value", "fdr"}.issubset(metrics), f"{dataset} missing a metric: {metrics}"
