from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from src.cross_dataset_anonymization_audit import build_anonymized_mapping, compare_rankings, run_anonymized_ranking

REPO_ROOT = Path(__file__).parent.parent


def _synthetic_wide(n_genes=40):
    rng = np.random.default_rng(1)
    genes = [f"GENE_{chr(65 + i % 26)}{i}" for i in range(n_genes)]
    df = pd.DataFrame({"gene": genes})
    for prefix, effect_col, p_col, fdr_col, testable_col in [
        ("crispr", "crispr_effect", "crispr_p", "crispr_fdr", "crispr_testable"),
        ("gse118713", "gse118713_log2fc", "gse118713_p", "gse118713_fdr", "gse118713_testable"),
    ]:
        df[effect_col] = rng.normal(0, 1, n_genes)
        df[p_col] = rng.uniform(0, 1, n_genes)
        df[fdr_col] = rng.uniform(0, 1, n_genes)
        df[testable_col] = True
    for prefix in ["gse245601", "gse240112"]:
        suffix = "epi" if prefix == "gse245601" else "tumor"
        other_suffix = "malignant" if prefix == "gse245601" else "epi"
        df[f"{prefix}_{suffix}_log2fc"] = rng.normal(0, 1, n_genes)
        df[f"{prefix}_{suffix}_p"] = rng.uniform(0, 1, n_genes)
        df[f"{prefix}_{suffix}_fdr"] = rng.uniform(0, 1, n_genes)
        df[f"{prefix}_{other_suffix}_log2fc"] = rng.normal(0, 1, n_genes)
        df[f"{prefix}_{other_suffix}_p"] = rng.uniform(0, 1, n_genes)
        df[f"{prefix}_{other_suffix}_fdr"] = rng.uniform(0, 1, n_genes)
        df[f"{prefix}_testable"] = True
    df["gse111151_log2fc"] = rng.normal(0, 1, n_genes)
    df["gse111151_p"] = rng.uniform(0, 1, n_genes)
    df["gse111151_fdr"] = rng.uniform(0, 1, n_genes)
    df["gse111151_testable"] = True
    return df


class TestBuildAnonymizedMapping:
    def test_mapping_covers_every_gene_exactly_once(self):
        genes = pd.Series([f"G{i}" for i in range(30)])
        out = build_anonymized_mapping(genes)
        assert len(out) == 30
        assert not out["gene"].duplicated().any()
        assert not out["anon_id"].duplicated().any()

    def test_mapping_is_not_alphabetical(self):
        genes = pd.Series([f"G{i:03d}" for i in range(30)])
        out = build_anonymized_mapping(genes).sort_values("anon_id")
        # if the mapping were alphabetical, sorting by anon_id would reproduce the sorted gene order exactly
        assert out["gene"].tolist() != sorted(genes.tolist())

    def test_deterministic_across_calls(self):
        genes = pd.Series([f"G{i}" for i in range(30)])
        out1 = build_anonymized_mapping(genes)
        out2 = build_anonymized_mapping(genes)
        pd.testing.assert_frame_equal(out1, out2)


class TestAnonymizedRankingMatchesOriginal:
    def test_ranks_match_outside_of_genuine_ties(self):
        from src.cross_dataset_ranking import assign_coverage_tier, build_global_ranking, compute_dataset_percentiles

        wide = _synthetic_wide()
        with_pct = compute_dataset_percentiles(wide)
        with_tier = assign_coverage_tier(with_pct)
        _full, original_ranked = build_global_ranking(with_tier, min_datasets_testable=3)

        mapping = build_anonymized_mapping(wide["gene"])
        anon_ranked = run_anonymized_ranking(wide, mapping, min_datasets_testable=3)
        comparison = compare_rankings(original_ranked, anon_ranked, mapping)

        # with continuous random effect/p/fdr values, genuine exact ties across the full
        # sort hierarchy are vanishingly unlikely -- expect all ranks to match exactly
        assert comparison["ranks_match"].all()

    def test_comparison_covers_every_original_gene(self):
        wide = _synthetic_wide()
        from src.cross_dataset_ranking import assign_coverage_tier, build_global_ranking, compute_dataset_percentiles

        with_pct = compute_dataset_percentiles(wide)
        with_tier = assign_coverage_tier(with_pct)
        _full, original_ranked = build_global_ranking(with_tier, min_datasets_testable=3)
        mapping = build_anonymized_mapping(wide["gene"])
        anon_ranked = run_anonymized_ranking(wide, mapping, min_datasets_testable=3)
        comparison = compare_rankings(original_ranked, anon_ranked, mapping)
        assert set(original_ranked["gene"]).issubset(set(comparison["gene"]))


class TestRealData:
    def test_real_anonymization_comparison_if_present(self):
        path = REPO_ROOT / "results" / "tables" / "cross_dataset_genomewide" / "anonymization_comparison.tsv"
        if not path.exists():
            pytest.skip("anonymization comparison table not present in this checkout")
        out = pd.read_csv(path, sep="\t")
        mismatch_fraction = 1 - out["ranks_match"].mean()
        # allow a tiny fraction of genuine-tie mismatches, but the overwhelming majority must match
        assert mismatch_fraction < 0.01
