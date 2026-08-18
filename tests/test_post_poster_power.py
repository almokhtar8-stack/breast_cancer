"""Tests for the post-freeze exploratory minimum-detectable-effect analysis
(Task 2).

The central checks are that the power function is right (verified against an
independent closed form and against monotonicity properties that any correct
power function must satisfy), that no observed p-value is ever turned into a
power value, and that a gene missing from a dataset produces an explicit
``not_tested`` row rather than a silent null result.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy import stats

from src import post_poster_power as power
from src import post_poster_meta_analysis as meta

ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = ROOT / "results" / "post_poster" / "power"

pytestmark = pytest.mark.skipif(
    not (OUT_DIR / "candidate_minimum_detectable_effects.tsv").exists(),
    reason="run `python -m src.post_poster_power` first",
)


# ---------------------------------------------------------------------------
# the power function itself
# ---------------------------------------------------------------------------
def test_power_at_zero_effect_equals_alpha():
    """With no true effect the rejection rate is exactly the significance
    level -- the single most basic property of a correct power function."""
    for alpha in (0.05, 0.01, 0.003846):
        assert power.power_at(0.0, se=0.2, df_total=10.0, alpha=alpha) == pytest.approx(alpha)


def test_power_increases_with_effect_size():
    p = [power.power_at(d, se=0.2, df_total=10.0, alpha=0.05)
         for d in (0.1, 0.3, 0.5, 1.0)]
    assert p == sorted(p)
    assert p[-1] > 0.99


def test_power_decreases_as_the_standard_error_grows():
    p = [power.power_at(0.5, se=s, df_total=10.0, alpha=0.05) for s in (0.1, 0.2, 0.4)]
    assert p[0] > p[1] > p[2]


def test_power_decreases_as_alpha_tightens():
    loose = power.power_at(0.5, se=0.2, df_total=10.0, alpha=0.05)
    tight = power.power_at(0.5, se=0.2, df_total=10.0, alpha=0.05 / 13)
    assert tight < loose


def test_power_matches_the_large_df_normal_closed_form():
    """As df -> infinity the noncentral-F power tends to the two-sided normal
    power, which is computable independently."""
    se, delta, alpha = 0.25, 0.6, 0.05
    got = power.power_at(delta, se=se, df_total=1e7, alpha=alpha)
    z = stats.norm.isf(alpha / 2)
    lam = delta / se
    expected = stats.norm.sf(z - lam) + stats.norm.cdf(-z - lam)
    assert got == pytest.approx(expected, rel=1e-4)


def test_minimum_detectable_effect_is_the_inverse_of_the_power_function():
    """The hand-checkable round trip: power at the returned MDE is exactly the
    target."""
    for se, df, alpha in [(0.2, 10.0, 0.05), (0.15, 11.3, 0.05),
                          (0.35, 9.16, 0.05 / 13), (0.5, 6.0, 0.01)]:
        mde = power.minimum_detectable_effect(se, df, alpha)
        assert np.isfinite(mde) and mde > 0
        assert power.power_at(mde, se, df, alpha) == pytest.approx(power.TARGET_POWER,
                                                                   abs=1e-8)


def test_minimum_detectable_effect_matches_the_large_df_closed_form():
    """At large df, MDE ~ (z_{1-a/2} + z_{power}) * se, ignoring the negligible
    opposite-tail term."""
    se, alpha = 0.2, 0.05
    mde = power.minimum_detectable_effect(se, 1e7, alpha)
    expected = (stats.norm.isf(alpha / 2) + stats.norm.isf(1 - power.TARGET_POWER)) * se
    assert mde == pytest.approx(expected, rel=1e-3)


def test_minimum_detectable_effect_scales_linearly_with_the_standard_error():
    a = power.minimum_detectable_effect(0.2, 10.0, 0.05)
    b = power.minimum_detectable_effect(0.4, 10.0, 0.05)
    assert b == pytest.approx(2 * a, rel=1e-8)


def test_minimum_detectable_effect_is_nan_for_unusable_inputs():
    assert np.isnan(power.minimum_detectable_effect(np.nan, 10.0, 0.05))
    assert np.isnan(power.minimum_detectable_effect(0.0, 10.0, 0.05))
    assert np.isnan(power.minimum_detectable_effect(0.2, 0.0, 0.05))


# ---------------------------------------------------------------------------
# alpha handling
# ---------------------------------------------------------------------------
def test_observed_bh_threshold_is_the_largest_rejected_p_value():
    p = np.array([0.001, 0.01, 0.02, 0.5])
    fdr = np.array([0.004, 0.02, 0.04, 0.5])
    assert power.observed_bh_threshold(p, fdr, level=0.05) == pytest.approx(0.02)


def test_observed_bh_threshold_is_nan_when_nothing_is_rejected():
    p = np.array([0.4, 0.6])
    fdr = np.array([0.8, 0.9])
    assert np.isnan(power.observed_bh_threshold(p, fdr, level=0.05))


def test_bonferroni_alpha_is_derived_from_the_candidate_count():
    assert power.N_CANDIDATES == len(meta.thirteen_candidates())
    alphas = pd.read_csv(OUT_DIR / "alpha_definitions.tsv", sep="\t")
    bonf = alphas[alphas["alpha_name"] == "bonferroni_13"]["alpha"].unique()
    assert bonf == pytest.approx([0.05 / 13])


def test_observed_bh_alpha_is_never_described_as_fdr_power():
    """Guards against the outcome-dependent alpha being relabelled as
    prospective FDR power."""
    alphas = pd.read_csv(OUT_DIR / "alpha_definitions.tsv", sep="\t")
    text = " ".join(alphas["alpha_definition"].unique()).lower()
    assert "descriptive only" in text
    assert "not a prospective design alpha" in text


# ---------------------------------------------------------------------------
# the interpretation column
# ---------------------------------------------------------------------------
def test_interpretation_never_claims_a_genuine_negative():
    """Failure to reject is not evidence of absence, so no label may assert
    one."""
    labels = {power.interpret(b, e, t)
              for b in (0.0, 1.0, np.nan) for e in (0.0, 1.0, np.nan)
              for t in (True, False)}
    joined = " ".join(labels).lower()
    for banned in ("genuine_negative", "informative_negative", "no_effect",
                   "absence", "proven"):
        assert banned not in joined


def test_interpretation_marks_an_untested_gene_as_not_tested():
    assert power.interpret(np.nan, np.nan, tested=False) == "not_tested"
    assert power.interpret(1.0, 1.0, tested=False) == "not_tested"


def test_interpretation_marks_an_effect_below_the_mde_as_uninformative():
    assert power.interpret(1.0, 0.0, tested=True).startswith("null_uninformative")


def test_interpretation_distinguishes_the_external_benchmark_case():
    assert power.interpret(0.0, 0.0, tested=True) == "sensitive_to_observed_magnitude"
    assert power.interpret(0.0, 1.0, tested=True) == (
        "sensitive_to_observed_magnitude_but_not_to_external_benchmark")


# ---------------------------------------------------------------------------
# outputs
# ---------------------------------------------------------------------------
@pytest.fixture(scope="module")
def table() -> pd.DataFrame:
    return pd.read_csv(OUT_DIR / "candidate_minimum_detectable_effects.tsv", sep="\t")


def test_table_covers_every_candidate_dataset_and_alpha(table):
    genes = meta.thirteen_candidates()
    expected = len(genes) * len(meta.DATASETS) * len(power.ALPHA_SPECS)
    assert len(table) == expected
    assert set(table["gene_symbol"]) == set(genes)


def test_mde_uses_only_the_information_based_standard_error(table):
    """An MDE built on f_calibrated would contain the observed effect and so
    would be partly outcome-derived. It must be refused, not merely avoided."""
    assert power.MDE_SE_VARIANT == meta.PRIMARY_SE_VARIANT == "wald"
    assert set(table["se_variant"]) == {"wald"}
    long = meta.build_long_table([d.key for d in meta.DATASETS])
    with pytest.raises(ValueError, match="outcome-dependent"):
        power.compute_power_table(long, meta.thirteen_candidates(),
                                  se_variant=meta.DIAGNOSTIC_SE_VARIANT)


def test_power_matches_an_independent_monte_carlo_estimate():
    """Validates the finite-df noncentral-F against simulation rather than
    against another call to the same distribution."""
    rng = np.random.default_rng(20260818)
    se, df, alpha, delta = 0.2, 9.0, 0.05, 0.35
    crit = stats.f.isf(alpha, dfn=1, dfd=df)
    n = 400_000
    # A statistic distributed as F(1, df, ncp) is (Z + sqrt(ncp))^2 / (chi2_df/df).
    ncp = (delta / se) ** 2
    num = (rng.standard_normal(n) + np.sqrt(ncp)) ** 2
    den = rng.chisquare(df, size=n) / df
    simulated = float(np.mean(num / den > crit))
    assert power.power_at(delta, se, df, alpha) == pytest.approx(simulated, abs=3e-3)


def test_not_tested_reason_distinguishes_filtering_from_absent_annotation():
    """USP17L29 is present before filtering in two datasets and genuinely
    unannotated in the third; the reason column must not conflate those."""
    assert "removed by filterByExpr" in power.not_tested_reason("USP17L29", "gse111151")
    assert "removed by filterByExpr" in power.not_tested_reason("USP17L29", "gse118713")
    assert "never measured" in power.not_tested_reason("USP17L29", "gse240112")


def test_untested_genes_have_an_explicit_reason_and_no_mde(table):
    untested = table[table["power_interpretation"] == "not_tested"]
    assert len(untested) > 0
    assert untested["mde80"].isna().all()
    assert (untested["not_tested_reason"].str.len() > 0).all()
    assert set(untested["gene_symbol"]) == {"USP17L29"}


def test_mde_is_positive_wherever_the_gene_was_tested(table):
    tested = table[table["tested_in_dataset"] & table["mde80"].notna()]
    assert len(tested) > 0
    assert (tested["mde80"] > 0).all()


def test_below_mde_flag_agrees_with_the_numbers_it_summarises(table):
    tested = table[table["observed_abs_log2fc_below_mde80"].notna()]
    recomputed = (tested["observed_abs_log2fc"] < tested["mde80"]).astype(float)
    assert (recomputed.to_numpy()
            == tested["observed_abs_log2fc_below_mde80"].to_numpy()).all()


def test_tighter_alpha_never_lowers_the_minimum_detectable_effect(table):
    """Bonferroni is stricter than nominal 0.05, so its MDE must be at least as
    large for every gene and dataset."""
    wide = table[table["alpha_name"] == "nominal_0.05"]
    tight = table[table["alpha_name"] == "bonferroni_13"]
    merged = wide.merge(tight, on=["gene_symbol", "dataset_key"], suffixes=("_w", "_t"))
    both = merged[merged["mde80_w"].notna() & merged["mde80_t"].notna()]
    assert len(both) > 0
    assert (both["mde80_t"] >= both["mde80_w"] - 1e-12).all()


def test_mde_does_not_depend_on_the_observed_effect(table):
    """The observed-power fallacy would make MDE a function of the observed
    p-value. It must instead depend only on SE, df and alpha."""
    sub = table[table["tested_in_dataset"] & table["mde80"].notna()].copy()
    recomputed = [
        power.minimum_detectable_effect(r.se, r.df_total, r.alpha)
        for r in sub.itertuples()
    ]
    assert np.allclose(recomputed, sub["mde80"].to_numpy(), rtol=1e-9)


def test_every_output_table_is_labelled_post_freeze():
    for path in sorted(OUT_DIR.glob("*.tsv")):
        df = pd.read_csv(path, sep="\t", nrows=5)
        assert meta.POST_FREEZE_LABEL in df.columns, f"{path.name} not labelled"


def test_summary_counts_match_the_detail_table(table):
    summary = pd.read_csv(OUT_DIR / "dataset_sensitivity_summary.tsv", sep="\t")
    detail = table
    for row in summary.itertuples():
        sub = detail[(detail["dataset_key"] == row.dataset_key)
                     & (detail["alpha_name"] == row.alpha_name)]
        assert int(sub["tested_in_dataset"].sum()) == row.n_candidates_tested
        assert int((~sub["tested_in_dataset"]).sum()) == row.n_candidates_not_tested
        assert (int((sub["observed_abs_log2fc_below_mde80"] == 1).sum())
                == row.n_observed_below_mde80)
