from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import yaml

from src.gse245601_cnv_method_comparison import LABEL_MALIGNANT, LABEL_NONMALIGNANT
from src.gse245601_infercnv_score_metric_diagnostics import (
    GROUP_DISAGREEMENT,
    GROUP_GOOD_CONTROL,
    GROUP_SIGNAL_QUALITY,
    ScoreMetricConfig,
    build_group_comparison_summary,
    build_tumor10_diagnostic_table,
    compute_per_cell_failure_category,
    load_score_metric_table,
    tumor_group_for,
)

REPO_ROOT = Path(__file__).parent.parent
CONFIG_PATH = REPO_ROOT / "config" / "config.yaml"


def _real_config() -> dict:
    with open(CONFIG_PATH) as f:
        return yaml.safe_load(f)


def _real_cfg() -> ScoreMetricConfig:
    return ScoreMetricConfig.from_config(_real_config())


def _synthetic_metric_df() -> pd.DataFrame:
    rows = []
    for patient, condition, base_score in [
        ("Tumor_01", "Control", 0.003),
        ("Tumor_01", "Tamoxifen", 0.003),
        ("Tumor_02", "Control", 0.02),
        ("Tumor_02", "Tamoxifen", 0.02),
        ("Tumor_10", "Control", 0.005),
        ("Tumor_10", "Tamoxifen", 0.005),
    ]:
        for i in range(10):
            rows.append(
                {
                    "cell_id": f"{patient}_{condition}_c{i}",
                    "sample_id": f"{patient}_{condition}",
                    "patient": patient,
                    "condition": condition,
                    "cnv_score": base_score + 0.0001 * i,
                    "max_abs_deviation": 0.2,
                    "p95_abs_deviation": 0.1,
                    "p99_abs_deviation": 0.15,
                    "fraction_genes_dev_gt_0.05": 0.1,
                    "fraction_genes_dev_gt_0.1": 0.05,
                    "fraction_genes_dev_gt_0.15": 0.01,
                    "n_chromosomes_affected_0.05": 2,
                    "n_chromosomes_affected_0.1": 0,
                    "n_chromosomes_affected_0.15": 0,
                }
            )
    return pd.DataFrame(rows)


class TestTumorGroupFor:
    def test_classifies_each_configured_group(self):
        cfg = _real_cfg()
        assert tumor_group_for("Tumor_01", cfg) == GROUP_DISAGREEMENT
        assert tumor_group_for("Tumor_04", cfg) == GROUP_DISAGREEMENT
        assert tumor_group_for("Tumor_08", cfg) == GROUP_DISAGREEMENT
        assert tumor_group_for("Tumor_02", cfg) == GROUP_GOOD_CONTROL
        assert tumor_group_for("Tumor_03", cfg) == GROUP_GOOD_CONTROL
        assert tumor_group_for("Tumor_10", cfg) == GROUP_SIGNAL_QUALITY

    def test_raises_on_unrecognized_patient(self):
        cfg = _real_cfg()
        with pytest.raises(ValueError, match="not in any configured tumor group"):
            tumor_group_for("Tumor_99", cfg)


class TestLoadScoreMetricTable:
    def test_raises_when_recomputed_score_disagrees(self, tmp_path):
        df = pd.DataFrame(
            {
                "cell_id": ["c1"],
                "sample_id": ["S1"],
                "patient": ["S"],
                "condition": ["Control"],
                "primary_malignancy_label": [LABEL_NONMALIGNANT],
                "cnv_correlation_to_seed": [0.3],
                "cnv_score": [0.01],
                "cnv_score_recomputed": [0.05],  # deliberately mismatched
                "max_abs_deviation": [0.2],
                "p95_abs_deviation": [0.1],
                "p99_abs_deviation": [0.15],
            }
        )
        path = tmp_path / "bad.tsv"
        df.to_csv(path, sep="\t", index=False)
        with pytest.raises(ValueError, match="disagree"):
            load_score_metric_table(path)

    def test_loads_real_table_if_present(self):
        cfg = _real_cfg()
        if not cfg.diagnostics_tsv.exists():
            pytest.skip("score metric table not extracted in this checkout")
        df = load_score_metric_table(cfg.diagnostics_tsv)
        assert len(df) > 0
        assert set(cfg.selected_samples) == set(df["sample_id"].unique())


class TestBuildGroupComparisonSummary:
    def test_one_row_per_patient_condition_with_correct_group_labels(self):
        cfg = _real_cfg()
        metric_df = _synthetic_metric_df()
        out = build_group_comparison_summary(metric_df, cfg)
        assert len(out) == 6
        assert set(out["tumor_group"]) == {GROUP_DISAGREEMENT, GROUP_GOOD_CONTROL, GROUP_SIGNAL_QUALITY}
        tumor01 = out.loc[(out["patient"] == "Tumor_01") & (out["condition"] == "Control")].iloc[0]
        assert tumor01["tumor_group"] == GROUP_DISAGREEMENT
        assert tumor01["median_cnv_score"] == pytest.approx(0.003 + 0.0001 * 4.5, abs=1e-6)

    def test_ordered_disagreement_then_good_control_then_signal_quality(self):
        cfg = _real_cfg()
        metric_df = _synthetic_metric_df()
        out = build_group_comparison_summary(metric_df, cfg)
        groups_in_order = out["tumor_group"].tolist()
        first_good_idx = groups_in_order.index(GROUP_GOOD_CONTROL)
        first_signal_idx = groups_in_order.index(GROUP_SIGNAL_QUALITY)
        assert all(g == GROUP_DISAGREEMENT for g in groups_in_order[:first_good_idx])
        assert first_signal_idx > first_good_idx


class TestComputePerCellFailureCategory:
    def test_real_data_categories_are_valid_and_complete(self):
        cfg = _real_cfg()
        if not cfg.infercnv_labels_tsv.exists():
            pytest.skip("frozen label table not present in this checkout")
        out = compute_per_cell_failure_category(cfg)
        valid = {"passes_both_malignant", "fails_cnv_score_only", "fails_correlation_only", "fails_both"}
        assert set(out["failure_category"].unique()) <= valid
        assert out["failure_category"].isna().sum() == 0


class TestBuildTumor10DiagnosticTable:
    def test_real_data_only_contains_tumor10_and_matches_expected_counts(self):
        cfg = _real_cfg()
        if not cfg.diagnostics_tsv.exists() or not cfg.infercnv_labels_tsv.exists():
            pytest.skip("frozen tables not present in this checkout")
        metric_df = load_score_metric_table(cfg.diagnostics_tsv)
        out = build_tumor10_diagnostic_table(metric_df, cfg)
        assert set(out["patient"].unique()) == {"Tumor_10"}
        by_sample_cat = out.groupby(["sample_id", "failure_category"]).size()
        # cross-check against the independently-written Point-1 by-sample table
        by_sample_path = REPO_ROOT / "results" / "tables" / "cnv_method_comparison_by_sample.tsv"
        if by_sample_path.exists():
            by_sample = pd.read_csv(by_sample_path, sep="\t").set_index("sample_id")
            for sample_id in ("Tumor_10_Control", "Tumor_10_Tamoxifen"):
                expected_malignant = int(by_sample.loc[sample_id, "infercnv_malignant_count"])
                got_malignant = int(by_sample_cat.get((sample_id, "passes_both_malignant"), 0))
                assert got_malignant == expected_malignant
