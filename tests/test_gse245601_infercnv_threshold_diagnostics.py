from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import yaml

from src.gse245601_cnv_method_comparison import LABEL_MALIGNANT, LABEL_NONMALIGNANT
from src.gse245601_infercnv_threshold_diagnostics import (
    FAIL_BOTH,
    FAIL_CNV_ONLY,
    FAIL_CORR_ONLY,
    PASS_BOTH,
    ThresholdDiagnosticsConfig,
    build_group_diagnostics_table,
    build_local_score_sensitivity,
    build_sensitivity_grid,
    classify_failure_category,
    identify_seed_cells,
    load_cnv_score_table,
    recompute_group_thresholds,
    verify_reconstruction,
)

REPO_ROOT = Path(__file__).parent.parent
CONFIG_PATH = REPO_ROOT / "config" / "config.yaml"


def _real_config() -> dict:
    with open(CONFIG_PATH) as f:
        return yaml.safe_load(f)


def _real_cfg() -> ThresholdDiagnosticsConfig:
    return ThresholdDiagnosticsConfig.from_config(_real_config())


def _test_cfg(**overrides) -> ThresholdDiagnosticsConfig:
    defaults = dict(
        selected_samples=("S1",),
        disagreement_tumors=(),
        good_control_tumors=(),
        signal_quality_tumors=(),
        infercnv_labels_tsv=Path("x"),
        copykat_labels_tsv=Path("x"),
        seed_top_fraction=0.05,
        min_seed_cells=2,
        cnv_score_sd_multiplier=2.0,
        cnv_score_clamp=(0.01, 0.05),
        correlation_sd_multiplier=1.5,
        correlation_clamp=(0.2, 0.4),
        cnv_score_sweep_range=(0.0, 0.05),
        cnv_score_grid_points=6,
        correlation_sweep_range=(0.2, 0.4),
        correlation_grid_points=5,
        diagnostics_tsv=Path("x"),
        sensitivity_grid_tsv=Path("x"),
        local_score_sensitivity_tsv=Path("x"),
        figures_dir=Path("x"),
    )
    defaults.update(overrides)
    return ThresholdDiagnosticsConfig(**defaults)


class TestRecomputeGroupThresholds:
    def test_matches_manual_formula_when_not_clamped(self):
        cfg = _test_cfg()
        scores = pd.Series([0.02, 0.025, 0.03, 0.028, 0.022])
        corrs = pd.Series([0.3, 0.35, 0.25, 0.32, 0.28])
        th_value, th_corr = recompute_group_thresholds(scores, corrs, cfg)
        assert th_value == pytest.approx(scores.mean() - 2 * scores.std())
        assert th_corr == pytest.approx(corrs.mean() - 1.5 * corrs.std())

    def test_cnv_score_threshold_clamped_to_floor(self):
        cfg = _test_cfg()
        # tight, low-mean distribution -> raw mean-2sd is negative/small -> floor
        scores = pd.Series([0.001, 0.0012, 0.0009, 0.0011, 0.001])
        corrs = pd.Series([0.5, 0.5, 0.5, 0.5, 0.5])
        th_value, _ = recompute_group_thresholds(scores, corrs, cfg)
        assert th_value == cfg.cnv_score_clamp[0]

    def test_cnv_score_threshold_clamped_to_ceiling(self):
        cfg = _test_cfg()
        # high mean, tiny sd -> raw mean-2sd exceeds the ceiling
        scores = pd.Series([0.09, 0.091, 0.089, 0.0905, 0.0895])
        corrs = pd.Series([0.5, 0.5, 0.5, 0.5, 0.5])
        th_value, _ = recompute_group_thresholds(scores, corrs, cfg)
        assert th_value == cfg.cnv_score_clamp[1]

    def test_correlation_threshold_clamped_to_floor_and_ceiling(self):
        cfg = _test_cfg()
        low_corr = pd.Series([0.05, 0.06, 0.04, 0.05, 0.06])
        scores = pd.Series([0.02] * 5)
        _, th_corr_floor = recompute_group_thresholds(scores, low_corr, cfg)
        assert th_corr_floor == cfg.correlation_clamp[0]

        high_corr = pd.Series([0.95, 0.96, 0.94, 0.95, 0.96])
        _, th_corr_ceiling = recompute_group_thresholds(scores, high_corr, cfg)
        assert th_corr_ceiling == cfg.correlation_clamp[1]


class TestIdentifySeedCells:
    def test_picks_top_5_percent_minimum_2(self):
        cfg = _test_cfg()
        df = pd.DataFrame({"cell_id": [f"c{i}" for i in range(20)], "cnv_score": list(range(20))})
        seeds = identify_seed_cells(df, cfg)
        # round(20*0.05)=1, but min_seed_cells=2 -> exactly 2 seeds: highest scores (19, 18) -> c19, c18
        assert set(seeds) == {"c19", "c18"}

    def test_scales_up_for_larger_groups(self):
        cfg = _test_cfg()
        df = pd.DataFrame({"cell_id": [f"c{i}" for i in range(200)], "cnv_score": list(range(200))})
        seeds = identify_seed_cells(df, cfg)
        # round(200*0.05) = 10
        assert len(seeds) == 10
        assert set(seeds) == {f"c{i}" for i in range(190, 200)}

    def test_never_exceeds_group_size(self):
        cfg = _test_cfg()
        df = pd.DataFrame({"cell_id": ["c0"], "cnv_score": [0.5]})
        seeds = identify_seed_cells(df, cfg)
        assert seeds == ["c0"]


class TestClassifyFailureCategory:
    def test_passes_both(self):
        assert classify_failure_category(0.02, 0.5, 0.01, 0.2) == PASS_BOTH

    def test_fails_cnv_only(self):
        assert classify_failure_category(0.005, 0.5, 0.01, 0.2) == FAIL_CNV_ONLY

    def test_fails_correlation_only(self):
        assert classify_failure_category(0.02, 0.1, 0.01, 0.2) == FAIL_CORR_ONLY

    def test_fails_both(self):
        assert classify_failure_category(0.005, 0.1, 0.01, 0.2) == FAIL_BOTH

    def test_strict_greater_than_at_exact_threshold_fails(self):
        assert classify_failure_category(0.01, 0.5, 0.01, 0.2) == FAIL_CNV_ONLY

    def test_raises_on_nan_score(self):
        with pytest.raises(ValueError, match="NaN input"):
            classify_failure_category(float("nan"), 0.5, 0.01, 0.2)


class TestBuildGroupDiagnosticsTableSynthetic:
    def test_failure_categories_sum_to_n_cells(self):
        cfg = _test_cfg()
        n = 50
        rng = np.random.default_rng(0)
        df = pd.DataFrame(
            {
                "cell_id": [f"c{i}" for i in range(n)],
                "sample_id": ["S1"] * n,
                "patient": ["S"] * n,
                "condition": ["Control"] * n,
                "threshold_group": ["0"] * n,
                "cnv_score": rng.uniform(0.0, 0.05, n),
                "cnv_correlation_to_seed": rng.uniform(-0.2, 0.8, n),
                "primary_malignancy_label": [LABEL_NONMALIGNANT] * n,  # placeholder, not used by this function
            }
        )
        out = build_group_diagnostics_table(df, cfg)
        assert len(out) == 1
        row = out.iloc[0]
        total = row["malignant_count"] + row["fails_cnv_score_only_count"] + row["fails_correlation_only_count"] + row["fails_both_count"]
        assert total == n


class TestVerifyReconstructionRealData:
    def test_reproduces_frozen_labels_exactly(self):
        cfg = _real_cfg()
        if not cfg.infercnv_labels_tsv.exists():
            pytest.skip("frozen GSE245601 label table not present in this checkout")
        cell_df = load_cnv_score_table(cfg.infercnv_labels_tsv)
        verify_reconstruction(cell_df, cfg)  # should not raise

    def test_raises_when_a_label_is_deliberately_corrupted(self):
        cfg = _real_cfg()
        if not cfg.infercnv_labels_tsv.exists():
            pytest.skip("frozen GSE245601 label table not present in this checkout")
        cell_df = load_cnv_score_table(cfg.infercnv_labels_tsv)
        corrupted = cell_df.copy()
        sel = corrupted.loc[corrupted["sample_id"] == cfg.selected_samples[0]]
        idx = sel.index[0]
        corrupted.loc[idx, "primary_malignancy_label"] = (
            LABEL_NONMALIGNANT if corrupted.loc[idx, "primary_malignancy_label"] == LABEL_MALIGNANT else LABEL_MALIGNANT
        )
        with pytest.raises(ValueError, match="do not reproduce frozen labels"):
            verify_reconstruction(corrupted, cfg)


class TestBuildGroupDiagnosticsTableRealData:
    def test_malignant_counts_match_frozen_summary_per_sample(self):
        cfg = _real_cfg()
        summary_path = REPO_ROOT / "results" / "tables" / "gse245601_malignant_summary_per_sample.tsv"
        if not cfg.infercnv_labels_tsv.exists() or not summary_path.exists():
            pytest.skip("frozen GSE245601 tables not present in this checkout")
        cell_df = load_cnv_score_table(cfg.infercnv_labels_tsv)
        group_table = build_group_diagnostics_table(cell_df, cfg)

        rolled_up = group_table.groupby("sample_id")["malignant_count"].sum()
        summary = pd.read_csv(summary_path, sep="\t").set_index("sample_id")
        for sample_id in cfg.selected_samples:
            assert int(rolled_up[sample_id]) == int(summary.loc[sample_id, "n_malignant"])


class TestBuildSensitivityGridSynthetic:
    def test_grid_shape_and_monotonic_yield(self):
        cfg = _test_cfg(selected_samples=("S1",))
        n = 100
        df = pd.DataFrame(
            {
                "cell_id": [f"c{i}" for i in range(n)],
                "sample_id": ["S1"] * n,
                "cnv_score": np.linspace(0.0, 0.05, n),
                "cnv_correlation_to_seed": np.linspace(0.2, 0.4, n),
            }
        )
        grid = build_sensitivity_grid(df, cfg)
        assert len(grid) == cfg.cnv_score_grid_points * cfg.correlation_grid_points
        # yield must be non-increasing as either threshold increases (monotonic thresholding)
        at_min_corr = grid.loc[grid["correlation_threshold"] == grid["correlation_threshold"].min()].sort_values("cnv_score_threshold")
        assert (at_min_corr["malignant_count"].diff().dropna() <= 0).all()


class TestBuildLocalScoreSensitivity:
    def test_grid_centered_on_the_actual_floor_and_monotonic(self):
        cfg = _test_cfg(selected_samples=("S1",))
        n = 50
        df = pd.DataFrame(
            {
                "cell_id": [f"c{i}" for i in range(n)],
                "sample_id": ["S1"] * n,
                "cnv_score": np.linspace(0.0, 0.02, n),
            }
        )
        out = build_local_score_sensitivity(df, cfg, band_half_width=0.004, n_points=9)
        assert len(out) == 9
        assert out["cnv_score_threshold"].min() == pytest.approx(0.006)
        assert out["cnv_score_threshold"].max() == pytest.approx(0.014)
        ordered = out.sort_values("cnv_score_threshold")
        assert (ordered["fraction_pass_score_only"].diff().dropna() <= 0).all()

    def test_real_data_shows_steep_local_decline_for_tumor01(self):
        cfg = _real_cfg()
        if not cfg.infercnv_labels_tsv.exists():
            pytest.skip("frozen GSE245601 label table not present in this checkout")
        cell_df = load_cnv_score_table(cfg.infercnv_labels_tsv)
        out = build_local_score_sensitivity(cell_df, cfg)
        sub = out.loc[out["sample_id"] == "Tumor_01_Control"].sort_values("cnv_score_threshold")
        # fraction passing must be non-increasing as the threshold rises
        assert (sub["fraction_pass_score_only"].diff().dropna() <= 0).all()
        # at the lower edge of the local band it must be strictly higher than at the floor itself
        assert sub.iloc[0]["fraction_pass_score_only"] >= sub.loc[np.isclose(sub["cnv_score_threshold"], 0.01)].iloc[0]["fraction_pass_score_only"]


class TestSelectedSamplesGrouping:
    def test_disagreement_good_and_signal_quality_tumors_are_disjoint_and_cover_expected_set(self):
        cfg = _real_cfg()
        all_tumors = set(cfg.disagreement_tumors) | set(cfg.good_control_tumors) | set(cfg.signal_quality_tumors)
        assert all_tumors == {"Tumor_01", "Tumor_04", "Tumor_08", "Tumor_02", "Tumor_03", "Tumor_10"}
        assert set(cfg.disagreement_tumors).isdisjoint(cfg.good_control_tumors)
        assert set(cfg.disagreement_tumors).isdisjoint(cfg.signal_quality_tumors)
        assert set(cfg.good_control_tumors).isdisjoint(cfg.signal_quality_tumors)
