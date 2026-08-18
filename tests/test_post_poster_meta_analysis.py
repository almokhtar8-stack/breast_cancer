"""Tests for the post-freeze exploratory meta-analysis (Task 1).

These exercise the statistics, not merely that the module runs: the pooling and
heterogeneity estimators are checked against hand-computable fixtures, and the
output tables are checked for the failure modes that matter here -- silently
dropped genes, unlabelled post-freeze status, and drift away from the frozen
inputs.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from src import post_poster_meta_analysis as meta

ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = ROOT / "results" / "post_poster" / "meta_analysis"
REFIT_DIR = ROOT / "results" / "post_poster" / "de_refit"

pytestmark = pytest.mark.skipif(
    not (OUT_DIR / "candidates_meta_analysis.tsv").exists(),
    reason="run `python -m src.post_poster_meta_analysis` first",
)


# ---------------------------------------------------------------------------
# statistics, against hand-computable fixtures
# ---------------------------------------------------------------------------
def test_fixed_effect_pooling_matches_closed_form():
    """Two studies, equal variances: the pooled estimate is the plain mean and
    the pooled SE is se/sqrt(2)."""
    eff = np.array([1.0, 2.0])
    var = np.array([0.25, 0.25])
    res = meta.pool(eff, var, tau2=0.0)
    assert res.pooled_effect == pytest.approx(1.5)
    assert res.pooled_se == pytest.approx(np.sqrt(0.25 / 2))


def test_fixed_effect_pooling_weights_by_inverse_variance():
    """A study with a quarter of the variance gets four times the weight."""
    eff = np.array([0.0, 1.0])
    var = np.array([1.0, 0.25])
    res = meta.pool(eff, var, tau2=0.0)
    assert res.pooled_effect == pytest.approx((0 * 1 + 1 * 4) / 5)


def test_cochran_q_is_zero_for_identical_effects():
    eff = np.array([0.4, 0.4, 0.4])
    var = np.array([0.1, 0.2, 0.3])
    assert meta.cochran_q(eff, var, 0.0) == pytest.approx(0.0)
    res = meta.pool(eff, var, 0.0)
    assert res.i2 == pytest.approx(0.0)


def test_cochran_q_hand_computed():
    """Two studies, equal variance v: Q = (e1-e2)^2 / (2v)."""
    eff = np.array([0.0, 1.0])
    var = np.array([0.5, 0.5])
    assert meta.cochran_q(eff, var, 0.0) == pytest.approx(1.0 / (2 * 0.5))


def test_dersimonian_laird_is_zero_when_q_below_its_df():
    eff = np.array([0.30, 0.31, 0.29])
    var = np.array([0.5, 0.5, 0.5])
    assert meta.tau2_dersimonian_laird(eff, var) == 0.0


def test_dersimonian_laird_matches_closed_form_for_equal_variances():
    """With k studies of common variance v, C = sum(w) - sum(w^2)/sum(w)
    reduces to (k-1)/v, so tau2 = (Q - (k-1)) * v / (k-1)."""
    eff = np.array([0.0, 1.0, 2.0])
    v = 0.2
    var = np.full(3, v)
    q = meta.cochran_q(eff, var, 0.0)
    expected = (q - 2) * v / 2
    assert meta.tau2_dersimonian_laird(eff, var) == pytest.approx(expected)


def test_paule_mandel_solves_q_equals_df():
    eff = np.array([0.0, 1.0, 2.5])
    var = np.array([0.2, 0.3, 0.25])
    t2 = meta.tau2_paule_mandel(eff, var)
    assert t2 > 0
    assert meta.cochran_q(eff, var, t2) == pytest.approx(len(eff) - 1, rel=1e-8)


def test_random_effects_interval_is_wider_than_fixed_effect_under_heterogeneity():
    eff = np.array([0.0, 1.0, 2.0])
    var = np.array([0.05, 0.05, 0.05])
    t2 = meta.tau2_dersimonian_laird(eff, var)
    assert t2 > 0
    fe = meta.pool(eff, var, 0.0)
    re = meta.pool(eff, var, t2)
    assert re.pooled_se > fe.pooled_se


def test_stouffer_equal_weights_matches_closed_form():
    """Two studies, same p, same direction: z = 2*z_i/sqrt(2) = sqrt(2)*z_i."""
    from scipy import stats
    p = np.array([0.05, 0.05])
    signs = np.array([1.0, 1.0])
    z, _ = meta.stouffer(p, signs, np.ones(2))
    assert z == pytest.approx(np.sqrt(2) * stats.norm.isf(0.025))


def test_stouffer_cancels_opposing_directions():
    z, p = meta.stouffer(np.array([0.01, 0.01]), np.array([1.0, -1.0]), np.ones(2))
    assert z == pytest.approx(0.0)
    assert p == pytest.approx(1.0)


def test_benjamini_hochberg_matches_a_hand_computed_case():
    p = np.array([0.01, 0.02, 0.03, 0.04])
    got = meta.benjamini_hochberg(p)
    # raw p*n/rank = .04, .04, .04, .04, then monotone from the right
    assert got == pytest.approx([0.04, 0.04, 0.04, 0.04])


def test_benjamini_hochberg_carries_nan_through():
    got = meta.benjamini_hochberg(np.array([0.01, np.nan, 0.02]))
    assert np.isnan(got[1])
    assert np.isfinite(got[0]) and np.isfinite(got[2])


# ---------------------------------------------------------------------------
# configuration and inputs
# ---------------------------------------------------------------------------
def test_thirteen_candidates_come_from_config_not_a_literal():
    genes = meta.thirteen_candidates()
    assert len(genes) == 13
    assert "USP34" in genes and "KDM1A" in genes


def test_gse245601_is_not_among_the_pooled_datasets():
    """GSE245601 is acute 12 h and must never enter resistance evidence."""
    keys = {d.key for d in meta.DATASETS}
    assert not any("245601" in k for k in keys)
    for spec in meta.ARMS.values():
        assert not any("245601" in d for d in spec["datasets"])


def test_gse240112_confounding_is_declared_in_the_dataset_spec():
    spec = meta.DATASETS_BY_KEY["gse240112"]
    assert spec.group_confounded_with_biobank is True
    assert "biobank" in spec.caveat


def test_gse118713_is_declared_pseudoreplicated_with_one_biological_unit():
    spec = meta.DATASETS_BY_KEY["gse118713"]
    assert spec.pseudoreplicated is True
    assert spec.n_biological_units == 1


def test_edger_only_arm_excludes_the_pseudoreplicated_dataset():
    assert "gse118713" not in meta.ARMS["edger_only"]["datasets"]
    assert meta.PRIMARY_ARM == "edger_only"


def test_only_the_information_based_se_may_weight_the_pooling():
    """f_calibrated is |log2FC|/sqrt(F): computed from the very statistic being
    summarised, so its inverse-variance weight would be outcome-dependent. It
    must stay a diagnostic."""
    assert meta.PRIMARY_SE_VARIANT == "wald"
    assert meta.DIAGNOSTIC_SE_VARIANT == "f_calibrated"
    doc = meta.__doc__
    assert "outcome-dependent" in doc
    assert "never as the primary weight" in doc


def test_module_declares_the_small_k_inference_limit():
    """Plug-in random effects at k=2-3 is anti-conservative; the module must say
    so rather than presenting the p-values as calibrated."""
    doc = meta.__doc__
    assert "anti-conservative" in doc
    assert "DESCRIPTIVE" in doc


def test_gse118713_model_and_contrast_sample_counts_are_distinguished():
    """It fits three groups jointly but the resistance contrast is 3 vs 3."""
    spec = meta.DATASETS_BY_KEY["gse118713"]
    assert spec.n_samples_in_model == 9
    assert spec.n_samples_in_contrast == 6


def test_design_effect_arm_inflates_only_gse118713_and_by_sqrt_three():
    infl = meta.ARMS["all3_de3"]["se_inflation"]
    assert set(infl) == {"gse118713"}
    assert infl["gse118713"] == pytest.approx(np.sqrt(3.0))


# ---------------------------------------------------------------------------
# outputs
# ---------------------------------------------------------------------------
@pytest.fixture(scope="module")
def candidates() -> pd.DataFrame:
    return pd.read_csv(OUT_DIR / "candidates_meta_analysis.tsv", sep="\t")


def test_every_candidate_appears_in_every_arm_and_variant(candidates):
    """A gene that cannot be pooled must still appear, with a reason -- it must
    never simply vanish from the table."""
    genes = set(meta.thirteen_candidates())
    for arm in meta.ARMS:
        for variant in meta.SE_VARIANTS:
            sub = candidates[(candidates["arm"] == arm)
                             & (candidates["se_variant"] == variant)]
            assert set(sub["gene_symbol"]) == genes, f"{arm}/{variant} lost genes"


def test_unpooled_rows_always_state_why(candidates):
    unpooled = candidates[candidates["pooled_effect"].isna()]
    assert len(unpooled) > 0, "expected USP17L29 to be unpoolable"
    assert (unpooled["exclusion_reason"].str.len() > 0).all()
    assert set(unpooled["gene_symbol"]) == {"USP17L29"}


def test_pooled_rows_used_every_dataset_in_their_arm(candidates):
    pooled = candidates[candidates["pooled_effect"].notna()]
    assert (pooled["n_datasets_used"] == pooled["n_datasets_in_arm"]).all()


def test_every_output_table_is_labelled_post_freeze():
    for path in sorted(OUT_DIR.glob("*.tsv")):
        df = pd.read_csv(path, sep="\t", nrows=5)
        assert meta.POST_FREEZE_LABEL in df.columns, f"{path.name} not labelled"


def test_confidence_interval_brackets_the_pooled_estimate(candidates):
    pooled = candidates[candidates["pooled_effect"].notna()]
    assert (pooled["ci_low"] <= pooled["pooled_effect"]).all()
    assert (pooled["pooled_effect"] <= pooled["ci_high"]).all()


def test_i2_is_within_range_and_zero_when_tau2_is_zero(candidates):
    pooled = candidates[candidates["pooled_effect"].notna()]
    assert pooled["i2"].between(0, 100).all()
    no_tau = pooled[pooled["tau2_dl"] == 0]
    assert (no_tau["i2"] == 0).all()


def test_fdr_is_never_smaller_than_its_p_value(candidates):
    pooled = candidates[candidates["p_value"].notna()]
    assert (pooled["fdr"] >= pooled["p_value"] - 1e-12).all()


def test_pooled_effect_lies_within_the_range_of_contributing_effects():
    """Inverse-variance pooling is a weighted mean, so it cannot fall outside
    the datasets' own estimates."""
    forest = pd.read_csv(OUT_DIR / "forest_plot_input_candidates.tsv", sep="\t")
    cand = pd.read_csv(OUT_DIR / "candidates_meta_analysis.tsv", sep="\t")
    pooled = cand[(cand["arm"] == "edger_only") & (cand["pooled_effect"].notna())
                  & (cand["se_variant"] == meta.PRIMARY_SE_VARIANT)]
    keys = set(meta.ARMS["edger_only"]["datasets"])
    for row in pooled.itertuples():
        contrib = forest[(forest["gene_symbol"] == row.gene_symbol)
                         & (forest["dataset_key"].isin(keys))]["log2fc"]
        assert contrib.min() - 1e-9 <= row.pooled_effect <= contrib.max() + 1e-9


def test_forest_table_has_one_row_per_candidate_and_dataset():
    forest = pd.read_csv(OUT_DIR / "forest_plot_input_candidates.tsv", sep="\t")
    genes = set(meta.thirteen_candidates())
    # USP17L29 is in no dataset, so 12 genes x 3 datasets.
    assert len(forest) == 36
    assert set(forest["gene_symbol"]) <= genes
    assert forest.groupby("gene_symbol")["dataset_key"].nunique().eq(3).all()


def test_vote_count_column_is_read_from_the_frozen_table_not_recomputed():
    comp = pd.read_csv(OUT_DIR / "vote_count_vs_pooled_effect.tsv", sep="\t")
    frozen = pd.read_csv(
        ROOT / "results/tables/evidence_freeze/final_candidate_evidence.tsv", sep="\t")
    merged = comp.merge(frozen[["gene", "resistance_fdr05_count"]],
                        left_on="gene_symbol", right_on="gene", suffixes=("", "_frozen"))
    assert len(merged) > 0
    assert (merged["resistance_fdr05_count"]
            == merged["resistance_fdr05_count_frozen"]).all()


def test_vote_count_comparison_uses_the_primary_arm_and_primary_se():
    comp = pd.read_csv(OUT_DIR / "vote_count_vs_pooled_effect.tsv", sep="\t")
    assert set(comp["primary_arm"]) == {meta.PRIMARY_ARM}
    assert set(comp["primary_se_variant"]) == {meta.PRIMARY_SE_VARIANT}
    # ranked on signed evidence strength, and the column name says so
    assert "pooled_z_rank" in comp.columns
    ranked = comp[comp["z"].notna()].sort_values("pooled_z_rank")
    assert list(ranked["z"]) == sorted(ranked["z"], reverse=True)


def test_symbol_collapse_audit_accounts_for_every_row():
    audit = pd.read_csv(OUT_DIR / "gene_symbol_collapse_audit.tsv", sep="\t")
    assert set(audit["dataset_key"]) == {d.key for d in meta.DATASETS}
    assert (audit["rows_in_refit_table"]
            == audit["rows_lost_no_gene_symbol"] + audit["rows_after_symbol_filter"]).all()
    assert (audit["rows_after_symbol_filter"]
            == audit["rows_lost_duplicate_symbol_collapse"]
            + audit["rows_out_unique_symbols"]).all()
    for row in audit.itertuples():
        raw = pd.read_csv(REFIT_DIR / f"{row.dataset_key}_model_stats.tsv.gz",
                          sep="\t", usecols=["gene_symbol"])
        assert row.rows_in_refit_table == len(raw)
        assert row.rows_out_unique_symbols == raw["gene_symbol"].nunique()


def test_refit_validation_records_row_accounting():
    val = pd.read_csv(REFIT_DIR / "refit_validation.tsv", sep="\t")
    for col in ("n_genes_before_filter", "n_genes_removed_by_filter",
                "n_genes_after_filter", "n_genes_joined_to_frozen"):
        assert col in val.columns, f"{col} missing from refit_validation.tsv"
    assert (val["n_genes_before_filter"]
            == val["n_genes_removed_by_filter"] + val["n_genes_after_filter"]).all()
    assert (val["n_genes_joined_to_frozen"] == val["n_genes"]).all()


def test_se_variant_diagnostic_records_the_wald_versus_quasi_f_gap():
    diag = pd.read_csv(OUT_DIR / "se_variant_diagnostic.tsv", sep="\t")
    assert set(diag["dataset_key"]) == {d.key for d in meta.DATASETS}
    limma = diag[diag["dataset_key"] == "gse118713"].iloc[0]
    # limma's quasi_f is the squared moderated t, so the two SEs coincide.
    assert limma["t_wald_sq_over_quasi_f_median"] == pytest.approx(1.0, rel=1e-6)
    edger = diag[diag["engine"] == "edgeR_glmQLF"]
    assert (edger["t_wald_sq_over_quasi_f_median"] < 1.0).all()


def test_refit_reproduced_the_frozen_tables():
    """The whole analysis is void if the refits drifted from frozen values."""
    val = pd.read_csv(REFIT_DIR / "refit_validation.tsv", sep="\t")
    assert len(val) == 3
    assert (val["max_abs_log2fc_diff_vs_frozen"] < 1e-8).all()
    assert (val["max_abs_p_diff_vs_frozen"] < 1e-10).all()
