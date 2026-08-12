from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from src.cross_dataset_ranking import (
    assign_coverage_tier,
    build_global_ranking,
    build_high_signal_low_coverage,
    compute_dataset_percentiles,
    compute_within_dataset_percentile,
)

REPO_ROOT = Path(__file__).parent.parent


class TestComputeWithinDatasetPercentile:
    def test_best_gene_gets_percentile_one(self):
        effect = pd.Series([2.0, 0.5, 1.0])
        p = pd.Series([0.001, 0.5, 0.1])
        fdr = pd.Series([0.005, 0.6, 0.2])
        out = compute_within_dataset_percentile(effect, p, fdr)
        assert out.iloc[0] == pytest.approx(1.0)

    def test_worst_gene_gets_percentile_zero(self):
        effect = pd.Series([2.0, 0.5, 1.0])
        p = pd.Series([0.001, 0.5, 0.1])
        fdr = pd.Series([0.005, 0.6, 0.2])
        out = compute_within_dataset_percentile(effect, p, fdr)
        assert out.iloc[1] == pytest.approx(0.0)

    def test_untestable_gene_gets_nan_not_zero(self):
        effect = pd.Series([1.0, np.nan])
        p = pd.Series([0.01, np.nan])
        fdr = pd.Series([0.02, np.nan])
        out = compute_within_dataset_percentile(effect, p, fdr)
        assert pd.isna(out.iloc[1])
        assert out.iloc[0] == pytest.approx(1.0)

    def test_fdr_beats_larger_effect(self):
        # gene A: worse FDR but bigger effect; gene B: better FDR, smaller effect -> B ranks first
        effect = pd.Series([5.0, 0.1])
        p = pd.Series([0.2, 0.001])
        fdr = pd.Series([0.3, 0.002])
        out = compute_within_dataset_percentile(effect, p, fdr)
        assert out.iloc[1] > out.iloc[0]

    def test_single_testable_gene_gets_one(self):
        effect = pd.Series([1.0])
        p = pd.Series([0.5])
        fdr = pd.Series([0.5])
        out = compute_within_dataset_percentile(effect, p, fdr)
        assert out.iloc[0] == pytest.approx(1.0)

    def test_deterministic_across_repeated_calls(self):
        effect = pd.Series([2.0, 0.5, 1.0, -3.0])
        p = pd.Series([0.001, 0.5, 0.1, 0.2])
        fdr = pd.Series([0.005, 0.6, 0.2, 0.3])
        out1 = compute_within_dataset_percentile(effect, p, fdr)
        out2 = compute_within_dataset_percentile(effect, p, fdr)
        pd.testing.assert_series_equal(out1, out2)

    def test_exact_tie_broken_by_gene_when_supplied(self):
        # two genes with identical fdr/p/effect -- without an explicit gene tiebreak,
        # the outcome would silently depend on input row order; with `gene` supplied,
        # the alphabetically-earlier gene must always rank first regardless of row order
        effect = pd.Series([1.0, 1.0])
        p = pd.Series([0.01, 0.01])
        fdr = pd.Series([0.02, 0.02])
        gene_forward = pd.Series(["AAA", "ZZZ"])
        gene_reversed = pd.Series(["ZZZ", "AAA"])
        out_forward = compute_within_dataset_percentile(effect, p, fdr, gene_forward)
        out_reversed = compute_within_dataset_percentile(effect, p, fdr, gene_reversed)
        # AAA (index 0 in forward, index 1 in reversed) must get the same (best) percentile either way
        assert out_forward.iloc[0] == out_reversed.iloc[1]
        assert out_forward.iloc[0] > out_forward.iloc[1]


def _synthetic_wide():
    return pd.DataFrame(
        {
            "gene": ["G1", "G2", "G3", "G4"],
            "crispr_effect": [-2.0, 1.0, np.nan, 0.5], "crispr_p": [0.001, 0.1, np.nan, 0.5], "crispr_fdr": [0.01, 0.2, np.nan, 0.6], "crispr_testable": [True, True, False, True],
            "gse118713_log2fc": [0.9, np.nan, np.nan, 0.2], "gse118713_p": [0.01, np.nan, np.nan, 0.3], "gse118713_fdr": [0.03, np.nan, np.nan, 0.4], "gse118713_testable": [True, False, False, True],
            "gse245601_epi_log2fc": [0.5, 0.3, np.nan, np.nan], "gse245601_epi_p": [0.02, 0.3, np.nan, np.nan], "gse245601_epi_fdr": [0.04, 0.4, np.nan, np.nan],
            "gse245601_malignant_log2fc": [0.4, np.nan, np.nan, np.nan], "gse245601_malignant_p": [0.05, np.nan, np.nan, np.nan], "gse245601_malignant_fdr": [0.06, np.nan, np.nan, np.nan],
            "gse245601_testable": [True, True, False, False],
            "gse240112_tumor_log2fc": [0.6, np.nan, 0.2, np.nan], "gse240112_tumor_p": [0.01, np.nan, 0.3, np.nan], "gse240112_tumor_fdr": [0.02, np.nan, 0.4, np.nan],
            "gse240112_epi_log2fc": [0.5, 0.1, 0.3, np.nan], "gse240112_epi_p": [0.02, 0.4, 0.2, np.nan], "gse240112_epi_fdr": [0.03, 0.5, 0.3, np.nan],
            "gse240112_testable": [True, True, True, False],
            "gse111151_log2fc": [0.3, 0.1, np.nan, np.nan], "gse111151_p": [0.1, 0.5, np.nan, np.nan], "gse111151_fdr": [0.3, 0.6, np.nan, np.nan], "gse111151_testable": [True, True, False, False],
        }
    )


class TestComputeDatasetPercentiles:
    def test_gse245601_is_mean_of_both_tracks_when_both_available(self):
        wide = _synthetic_wide()
        out = compute_dataset_percentiles(wide).set_index("gene")
        expected = (out.loc["G1", "gse245601_track_a_percentile"] + out.loc["G1", "gse245601_track_b_percentile"]) / 2
        assert out.loc["G1", "gse245601_evidence_percentile"] == pytest.approx(expected)

    def test_gse245601_uses_single_track_when_only_one_available(self):
        wide = _synthetic_wide()
        out = compute_dataset_percentiles(wide).set_index("gene")
        assert out.loc["G2", "gse245601_evidence_percentile"] == pytest.approx(out.loc["G2", "gse245601_track_a_percentile"])
        assert bool(out.loc["G2", "gse245601_one_track_only"])

    def test_gse240112_dataset_percentile_equals_tumor_only_never_blends_epi(self):
        wide = _synthetic_wide()
        out = compute_dataset_percentiles(wide).set_index("gene")
        # G2 has an epithelial value but no tumor-cell value -> dataset percentile must be NaN, not derived from epi
        assert pd.isna(out.loc["G2", "gse240112_tumor_log2fc"])
        assert pd.isna(out.loc["G2", "gse240112_evidence_percentile"])
        assert out.loc["G1", "gse240112_evidence_percentile"] == pytest.approx(out.loc["G1", "gse240112_tumor_percentile"])


class TestAssignCoverageTier:
    def test_tiers_match_testable_counts(self):
        wide = compute_dataset_percentiles(_synthetic_wide())
        out = assign_coverage_tier(wide).set_index("gene")
        assert out.loc["G1", "coverage_tier"] == "A"
        assert out.loc["G3", "n_datasets_testable"] == 1
        assert out.loc["G3", "coverage_tier"] == "E"


class TestBuildGlobalRanking:
    def test_low_coverage_gene_excluded_from_primary_ranking(self):
        wide = assign_coverage_tier(compute_dataset_percentiles(_synthetic_wide()))
        full, ranked = build_global_ranking(wide, min_datasets_testable=3)
        assert "G3" not in set(ranked["gene"])  # only 1 dataset testable
        assert "G3" in set(full["gene"])  # never dropped from the full table

    def test_ranking_deterministic_and_unique_ranks(self):
        wide = assign_coverage_tier(compute_dataset_percentiles(_synthetic_wide()))
        full, ranked = build_global_ranking(wide, min_datasets_testable=1)
        assert ranked["global_rank"].tolist() == list(range(1, len(ranked) + 1))
        full2, ranked2 = build_global_ranking(wide, min_datasets_testable=1)
        pd.testing.assert_frame_equal(ranked, ranked2)

    def test_no_gene_dropped_from_full_table(self):
        wide = assign_coverage_tier(compute_dataset_percentiles(_synthetic_wide()))
        full, ranked = build_global_ranking(wide, min_datasets_testable=3)
        assert len(full) == 4


def _synthetic_inputs_for_low_coverage():
    # EXTREME: testable only in crispr (1/5), clearly best (lowest FDR) of the two crispr-testable genes -> percentile 1.0
    # MIDDLING: testable only in crispr (1/5), the worse of the two -> percentile 0.0, must NOT be surfaced
    # HIGHCOV: testable in all 5 datasets -> above the coverage cutoff, must NOT be surfaced regardless of its own values
    return pd.DataFrame(
        {
            "gene": ["EXTREME", "MIDDLING", "HIGHCOV"],
            "crispr_effect": [-3.0, -0.1, 1.0], "crispr_p": [0.0001, 0.9, 0.01], "crispr_fdr": [0.0005, 0.95, 0.02], "crispr_testable": [True, True, True],
            "gse118713_log2fc": [np.nan, np.nan, 0.5], "gse118713_p": [np.nan, np.nan, 0.01], "gse118713_fdr": [np.nan, np.nan, 0.02], "gse118713_testable": [False, False, True],
            "gse245601_epi_log2fc": [np.nan, np.nan, 0.5], "gse245601_epi_p": [np.nan, np.nan, 0.01], "gse245601_epi_fdr": [np.nan, np.nan, 0.02],
            "gse245601_malignant_log2fc": [np.nan, np.nan, np.nan], "gse245601_malignant_p": [np.nan, np.nan, np.nan], "gse245601_malignant_fdr": [np.nan, np.nan, np.nan],
            "gse245601_testable": [False, False, True],
            "gse240112_tumor_log2fc": [np.nan, np.nan, 0.5], "gse240112_tumor_p": [np.nan, np.nan, 0.01], "gse240112_tumor_fdr": [np.nan, np.nan, 0.02],
            "gse240112_epi_log2fc": [np.nan, np.nan, np.nan], "gse240112_epi_p": [np.nan, np.nan, np.nan], "gse240112_epi_fdr": [np.nan, np.nan, np.nan],
            "gse240112_testable": [False, False, True],
            "gse111151_log2fc": [np.nan, np.nan, 0.5], "gse111151_p": [np.nan, np.nan, 0.01], "gse111151_fdr": [np.nan, np.nan, 0.02], "gse111151_testable": [False, False, True],
        }
    )


class TestBuildHighSignalLowCoverage:
    def test_extreme_low_coverage_gene_surfaced(self):
        wide = assign_coverage_tier(compute_dataset_percentiles(_synthetic_inputs_for_low_coverage()))
        full, _ranked = build_global_ranking(wide, min_datasets_testable=3)
        out = build_high_signal_low_coverage(full, min_datasets_testable=3, extreme_percentile_threshold=0.99)
        assert "EXTREME" in set(out["gene"])
        assert out.set_index("gene").loc["EXTREME", "max_percentile_dataset"] == "crispr_evidence_percentile"

    def test_non_extreme_low_coverage_gene_not_surfaced(self):
        wide = assign_coverage_tier(compute_dataset_percentiles(_synthetic_inputs_for_low_coverage()))
        full, _ranked = build_global_ranking(wide, min_datasets_testable=3)
        out = build_high_signal_low_coverage(full, min_datasets_testable=3, extreme_percentile_threshold=0.99)
        assert "MIDDLING" not in set(out["gene"])

    def test_high_coverage_genes_never_included(self):
        wide = assign_coverage_tier(compute_dataset_percentiles(_synthetic_inputs_for_low_coverage()))
        full, _ranked = build_global_ranking(wide, min_datasets_testable=3)
        out = build_high_signal_low_coverage(full, min_datasets_testable=3, extreme_percentile_threshold=0.99)
        assert "HIGHCOV" not in set(out["gene"])  # testable in 5/5, above the coverage cutoff regardless of its own values


class TestRealData:
    def test_real_ranking_if_present(self):
        full_path = REPO_ROOT / "results" / "tables" / "cross_dataset_genomewide" / "all_genes_cross_dataset_evidence_with_ranking.tsv"
        ranked_path = REPO_ROOT / "results" / "tables" / "cross_dataset_genomewide" / "global_ranking_eligible.tsv"
        if not full_path.exists() or not ranked_path.exists():
            pytest.skip("cross-dataset ranking tables not present in this checkout")
        full = pd.read_csv(full_path, sep="\t")
        ranked = pd.read_csv(ranked_path, sep="\t")
        assert not full["gene"].duplicated().any()
        assert ranked["global_rank"].tolist() == list(range(1, len(ranked) + 1))
        assert (ranked["n_datasets_testable"] >= 3).all()

    def test_real_high_signal_low_coverage_if_present(self):
        path = REPO_ROOT / "results" / "tables" / "cross_dataset_genomewide" / "high_signal_low_coverage.tsv"
        if not path.exists():
            pytest.skip("high_signal_low_coverage table not present in this checkout")
        out = pd.read_csv(path, sep="\t")
        assert (out["n_datasets_testable"] < 3).all()
        assert (out["max_evidence_percentile"] >= 0.99).all()
        assert "USP17L29" in set(out["gene"])  # known genome-wide-confirmed CRISPR-only extreme hit
