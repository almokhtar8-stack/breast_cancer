from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from src.cross_dataset_stability import build_ranking_stability, classify_stability

REPO_ROOT = Path(__file__).parent.parent


def _synthetic_df(n_genes=30):
    rng = np.random.default_rng(0)
    genes = [f"G{i}" for i in range(n_genes)]
    df = pd.DataFrame({"gene": genes})
    for d, pct_col, fdr_col, p_col, testable_col in [
        ("crispr", "crispr_evidence_percentile", "crispr_fdr", "crispr_p", "crispr_testable"),
        ("gse118713", "gse118713_evidence_percentile", "gse118713_fdr", "gse118713_p", "gse118713_testable"),
        ("gse245601", "gse245601_evidence_percentile", "gse245601_epi_fdr", "gse245601_epi_p", "gse245601_testable"),
        ("gse240112", "gse240112_evidence_percentile", "gse240112_tumor_fdr", "gse240112_tumor_p", "gse240112_testable"),
        ("gse111151", "gse111151_evidence_percentile", "gse111151_fdr", "gse111151_p", "gse111151_testable"),
    ]:
        df[pct_col] = rng.uniform(0, 1, n_genes)
        df[fdr_col] = rng.uniform(0, 1, n_genes)
        df[p_col] = df[fdr_col] * 0.5
        df[testable_col] = True
    # G0: robust, strong signal in every dataset
    for col in ["crispr_evidence_percentile", "gse118713_evidence_percentile", "gse245601_evidence_percentile", "gse240112_evidence_percentile", "gse111151_evidence_percentile"]:
        df.loc[df["gene"] == "G0", col] = 0.99
    for col in ["crispr_fdr", "gse118713_fdr", "gse245601_epi_fdr", "gse240112_tumor_fdr", "gse111151_fdr"]:
        df.loc[df["gene"] == "G0", col] = 0.001
    # G1: strong ONLY in crispr, weak everywhere else -> dataset-dependent
    df.loc[df["gene"] == "G1", "crispr_evidence_percentile"] = 0.99
    df.loc[df["gene"] == "G1", "crispr_fdr"] = 0.001
    for col in ["gse118713_evidence_percentile", "gse245601_evidence_percentile", "gse240112_evidence_percentile", "gse111151_evidence_percentile"]:
        df.loc[df["gene"] == "G1", col] = 0.1
    for col in ["gse118713_fdr", "gse245601_epi_fdr", "gse240112_tumor_fdr", "gse111151_fdr"]:
        df.loc[df["gene"] == "G1", col] = 0.9
    return df


def _main_ranked(df):
    from src.cross_dataset_ranking import assign_coverage_tier, build_global_ranking

    with_tier = assign_coverage_tier(df)
    _full, ranked = build_global_ranking(with_tier, min_datasets_testable=3)
    return ranked


class TestRankSubsetMatchesMainHierarchy:
    def test_hierarchy_scheme_over_all_5_datasets_reproduces_main_ranking(self):
        # regression test for a real bug caught by Codex review: _rank_subset's
        # "hierarchy" scheme must be key-for-key identical to
        # cross_dataset_ranking.build_global_ranking's sort hierarchy -- verified
        # here by confirming that applying it to ALL 5 datasets (not a leave-one-out
        # subset) gives EXACTLY the same ranks as the main global ranking itself
        from src.cross_dataset_stability import DATASET_NAMES, _rank_subset

        df = _synthetic_df()
        main_ranked = _main_ranked(df)
        subset_ranks = _rank_subset(df, DATASET_NAMES, "hierarchy", min_testable=3)
        main_ranks = main_ranked.set_index("gene")["global_rank"]
        pd.testing.assert_series_equal(subset_ranks.sort_index(), main_ranks.sort_index(), check_names=False)


class TestBuildRankingStability:
    def test_all_scheme_columns_present(self):
        df = _synthetic_df()
        main_ranked = _main_ranked(df)
        out = build_ranking_stability(df, main_ranked)
        expected_cols = {"rank_main", "rank_scheme_median_percentile", "rank_scheme_mean_percentile", "rank_scheme_fdr_count_first", "rank_scheme_top10pct_count_first"}
        expected_cols |= {f"rank_without_{d}" for d in ["crispr", "gse118713", "gse245601", "gse240112", "gse111151"]}
        assert expected_cols.issubset(set(out.columns))

    def test_robust_gene_has_low_ranks_everywhere(self):
        df = _synthetic_df()
        main_ranked = _main_ranked(df)
        out = build_ranking_stability(df, main_ranked).set_index("gene")
        assert out.loc["G0", "worst_rank"] <= 5

    def test_summary_stats_computed(self):
        df = _synthetic_df()
        main_ranked = _main_ranked(df)
        out = build_ranking_stability(df, main_ranked)
        for col in ["best_rank", "worst_rank", "median_rank", "rank_iqr", "n_top20_appearances", "n_top10_appearances", "n_top5_appearances"]:
            assert col in out.columns


class TestClassifyStability:
    def test_robust_gene_classified_robust(self):
        df = _synthetic_df()
        main_ranked = _main_ranked(df)
        stability = build_ranking_stability(df, main_ranked)
        labels = classify_stability(stability)
        idx = stability.index[stability["gene"] == "G0"][0]
        assert labels.loc[idx] == "ROBUST"

    def test_dataset_dependent_gene_flagged(self):
        df = _synthetic_df()
        main_ranked = _main_ranked(df)
        stability = build_ranking_stability(df, main_ranked)
        labels = classify_stability(stability)
        idx = stability.index[stability["gene"] == "G1"][0]
        # G1 loses its rank when crispr is removed (its only strong dataset) -> should not be ROBUST
        assert labels.loc[idx] != "ROBUST"

    def test_every_gene_gets_a_label(self):
        df = _synthetic_df()
        main_ranked = _main_ranked(df)
        stability = build_ranking_stability(df, main_ranked)
        labels = classify_stability(stability)
        assert labels.notna().all()
        assert set(labels.unique()).issubset({"ROBUST", "MODERATELY_STABLE", "DATASET_DEPENDENT"})


class TestRealData:
    def test_real_stability_table_if_present(self):
        path = REPO_ROOT / "results" / "tables" / "cross_dataset_genomewide" / "ranking_stability.tsv"
        if not path.exists():
            pytest.skip("ranking stability table not present in this checkout")
        out = pd.read_csv(path, sep="\t")
        assert not out["gene"].duplicated().any()
        assert set(out["stability_label"].unique()).issubset({"ROBUST", "MODERATELY_STABLE", "DATASET_DEPENDENT"})
