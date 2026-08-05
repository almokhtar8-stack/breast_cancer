import logging
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from src.gate1_checks import (
    BRANCH_CLASSIFIER,
    BRANCH_CONTINUOUS,
    BRANCH_NONE,
    _rank_by_effect_size,
    decide_gate1,
    direction_sanity_extraction,
    essentiality_contamination_check,
    load_and_validate_labels,
    load_ceg2_reference,
)


def _make_labels_df(n_genes: int, n_significant: int, fdr_threshold: float = 0.1) -> pd.DataFrame:
    """Build a synthetic, well-formed labels table with a chosen number of
    genes below fdr_threshold, everything else deterministic and finite."""
    genes = [f"GENE_{i:05d}" for i in range(n_genes)]
    fdr = np.full(n_genes, fdr_threshold + 0.5)
    fdr[:n_significant] = fdr_threshold / 2
    rng = np.random.default_rng(0)
    return pd.DataFrame(
        {
            "gene": genes,
            "n_guides": np.full(n_genes, 4),
            "effect_size": rng.normal(size=n_genes),
            "se": np.full(n_genes, 0.3),
            "p_value": np.full(n_genes, 0.5),
            "fdr": fdr,
        }
    )


# --- Gate 1 branch boundaries -----------------------------------------


@pytest.mark.parametrize(
    "n_significant, expected_branch",
    [
        (9, BRANCH_NONE),
        (10, BRANCH_CONTINUOUS),
        (29, BRANCH_CONTINUOUS),
        (30, BRANCH_CLASSIFIER),
    ],
)
def test_branch_boundaries(n_significant, expected_branch):
    df = _make_labels_df(n_genes=100, n_significant=n_significant)
    decision = decide_gate1(df, fdr_threshold=0.1, classifier_min=30, continuous_min=10)
    assert decision["n_passing"] == n_significant
    assert decision["branch_decision"] == expected_branch


# --- Proof that the significant-gene count is calculated, not hardcoded --


@pytest.mark.parametrize("n_significant", [0, 5, 15, 22, 28, 33, 47])
def test_significant_count_is_calculated_not_hardcoded(n_significant):
    df = _make_labels_df(n_genes=100, n_significant=n_significant)
    decision = decide_gate1(df, fdr_threshold=0.1, classifier_min=30, continuous_min=10)
    assert decision["n_passing"] == n_significant


def test_decide_gate1_recomputes_from_data_each_call():
    """Calling decide_gate1 on two different frames must not return the
    same cached/hardcoded count."""
    df_a = _make_labels_df(n_genes=100, n_significant=28)
    df_b = _make_labels_df(n_genes=100, n_significant=3)
    decision_a = decide_gate1(df_a, fdr_threshold=0.1, classifier_min=30, continuous_min=10)
    decision_b = decide_gate1(df_b, fdr_threshold=0.1, classifier_min=30, continuous_min=10)
    assert decision_a["n_passing"] == 28
    assert decision_b["n_passing"] == 3
    assert decision_a["n_passing"] != decision_b["n_passing"]


# --- load_and_validate_labels: structural checks ------------------------


def test_load_and_validate_labels_rejects_duplicate_gene(tmp_path: Path):
    df = _make_labels_df(n_genes=10, n_significant=2)
    df.loc[1, "gene"] = df.loc[0, "gene"]
    path = tmp_path / "labels.parquet"
    df.to_parquet(path, index=False)
    with pytest.raises(ValueError, match="duplicate"):
        load_and_validate_labels(path, n_fitted_genes_expected=10)


def test_load_and_validate_labels_rejects_null_gene(tmp_path: Path):
    df = _make_labels_df(n_genes=10, n_significant=2)
    df.loc[0, "gene"] = None
    path = tmp_path / "labels.parquet"
    df.to_parquet(path, index=False)
    with pytest.raises(ValueError, match="null"):
        load_and_validate_labels(path, n_fitted_genes_expected=10)


def test_load_and_validate_labels_rejects_missing_column(tmp_path: Path):
    df = _make_labels_df(n_genes=10, n_significant=2).drop(columns=["se"])
    path = tmp_path / "labels.parquet"
    df.to_parquet(path, index=False)
    with pytest.raises(ValueError, match="missing required columns"):
        load_and_validate_labels(path, n_fitted_genes_expected=10)


def test_load_and_validate_labels_rejects_non_finite_values(tmp_path: Path):
    df = _make_labels_df(n_genes=10, n_significant=2)
    df.loc[0, "effect_size"] = float("nan")
    path = tmp_path / "labels.parquet"
    df.to_parquet(path, index=False)
    with pytest.raises(ValueError, match="non-finite"):
        load_and_validate_labels(path, n_fitted_genes_expected=10)


def test_load_and_validate_labels_rejects_wrong_row_count(tmp_path: Path):
    df = _make_labels_df(n_genes=10, n_significant=2)
    path = tmp_path / "labels.parquet"
    df.to_parquet(path, index=False)
    with pytest.raises(ValueError, match="19103"):
        load_and_validate_labels(path, n_fitted_genes_expected=19103)


def test_load_and_validate_labels_accepts_well_formed_table(tmp_path: Path):
    df = _make_labels_df(n_genes=10, n_significant=2)
    path = tmp_path / "labels.parquet"
    df.to_parquet(path, index=False)
    loaded = load_and_validate_labels(path, n_fitted_genes_expected=10)
    assert len(loaded) == 10


# --- Deterministic ranking on ties --------------------------------------


def test_rank_by_effect_size_breaks_ties_by_gene_ascending():
    df = pd.DataFrame(
        {
            "gene": ["ZETA", "ALPHA", "MU", "BETA"],
            "effect_size": [-1.0, -1.0, -1.0, 2.0],
        }
    )
    ranked = _rank_by_effect_size(df)
    assert list(ranked["gene"]) == ["ALPHA", "MU", "ZETA", "BETA"]


def test_rank_by_effect_size_is_deterministic_across_calls():
    df = pd.DataFrame(
        {
            "gene": [f"G{i}" for i in range(20)],
            "effect_size": [0.0] * 20,
        }
    )
    first = list(_rank_by_effect_size(df)["gene"])
    second = list(_rank_by_effect_size(df.sample(frac=1, random_state=1))["gene"])
    assert first == second


# --- Essentiality-contamination check: universe and expected overlap ---


def test_essentiality_check_uses_full_fitted_gene_universe(caplog):
    n_genes = 1000
    df = _make_labels_df(n_genes=n_genes, n_significant=5)
    ceg2 = set(df["gene"].iloc[:50]) | {"NOT_IN_UNIVERSE_1", "NOT_IN_UNIVERSE_2"}
    with caplog.at_level(logging.INFO):
        summary = essentiality_contamination_check(df, ceg2, tested_set_size=50)
    assert summary["fitted_gene_universe_size"] == n_genes
    assert summary["ceg2_genes_in_universe"] == 50
    assert summary["ceg2_symbols_absent_from_universe"] == 2


def test_essentiality_check_expected_overlap_formula():
    n_genes = 200
    df = _make_labels_df(n_genes=n_genes, n_significant=5)
    ceg2 = set(df["gene"].iloc[:20])
    summary = essentiality_contamination_check(df, ceg2, tested_set_size=50)
    expected = 50 * (20 / 200)
    assert summary["expected_ceg2_overlap"] == pytest.approx(expected)
    assert summary["tested_set_size"] == 50


def test_essentiality_check_observed_overlap_matches_construction():
    """Force a known number of CEG2 genes into the bottom-50 by effect_size
    and check the observed overlap is counted correctly.

    Ground truth here is derived independently of ``_rank_by_effect_size``
    (the function under test elsewhere): effect_size is assigned in plain
    ascending order by gene index, so gene ``GENE_00000..00049`` are the
    bottom 50 by construction, not by calling the ranking helper.
    """
    n_genes = 300
    genes = [f"GENE_{i:05d}" for i in range(n_genes)]
    df = pd.DataFrame(
        {
            "gene": genes,
            "n_guides": 4,
            "effect_size": np.arange(n_genes, dtype=float),
            "se": 0.3,
            "p_value": 0.5,
            "fdr": 0.9,
        }
    )
    bottom_50_genes = genes[:50]
    ceg2 = set(bottom_50_genes[:7]) | set(genes[100:110])
    summary = essentiality_contamination_check(df, ceg2, tested_set_size=50)
    assert summary["observed_ceg2_overlap"] == 7
    assert summary["ceg2_genes_in_universe"] == 17


def test_essentiality_check_hypergeometric_p_value_matches_hand_calculation():
    """Small, fully worked-out hypergeometric case: universe=10, CEG2=3,
    tested_set=4, observed_overlap=2. By hand: P(X>=2) = [C(3,2)C(7,2) +
    C(3,3)C(7,1)] / C(10,4) = (63 + 7) / 210 = 1/3."""
    genes = [f"G{i}" for i in range(10)]
    df = pd.DataFrame(
        {
            "gene": genes,
            "n_guides": 4,
            "effect_size": np.arange(10, dtype=float),  # ascending: G0 most negative
            "se": 0.3,
            "p_value": 0.5,
            "fdr": 0.9,
        }
    )
    # Bottom-4 by effect_size (tested set) are G0-G3. CEG2 = {G0, G1, G5}:
    # 2 of the 3 CEG2 genes (G0, G1) fall in the tested set.
    ceg2 = {"G0", "G1", "G5"}
    summary = essentiality_contamination_check(df, ceg2, tested_set_size=4)
    assert summary["observed_ceg2_overlap"] == 2
    assert summary["ceg2_genes_in_universe"] == 3
    assert summary["hypergeometric_p_value_one_sided"] == pytest.approx(1 / 3, rel=1e-9)


def test_essentiality_check_rejects_tested_set_larger_than_universe():
    df = _make_labels_df(n_genes=10, n_significant=1)
    with pytest.raises(ValueError, match="smaller than requested"):
        essentiality_contamination_check(df, set(), tested_set_size=50)


def test_essentiality_check_does_not_log_gene_names(caplog):
    n_genes = 500
    df = _make_labels_df(n_genes=n_genes, n_significant=5)
    ranked = _rank_by_effect_size(df)
    bottom_genes = list(ranked["gene"].iloc[:50])
    ceg2 = set(bottom_genes[:10])
    with caplog.at_level(logging.DEBUG):
        essentiality_contamination_check(df, ceg2, tested_set_size=50)
    for gene in bottom_genes:
        assert gene not in caplog.text


def test_essentiality_check_never_names_individual_genes_worst_case(caplog):
    """No individual gene identity may appear in logs or outputs, even in
    the worst case where every CEG2-overlap gene lands in the tested
    (bottom-50) set. Uses placeholder identifiers -- this property does
    not depend on which specific genes are involved, so the project's
    real blind-control gene symbols are never used in this test."""
    n_genes = 300
    marker_genes = ("PLACEHOLDER_MARKER_A", "PLACEHOLDER_MARKER_B")
    genes = [f"GENE_{i:05d}" for i in range(n_genes - 2)] + list(marker_genes)
    effect_size = np.concatenate([np.linspace(1.0, 5.0, n_genes - 2), [-10.0, -9.0]])
    df = pd.DataFrame(
        {
            "gene": genes,
            "n_guides": 4,
            "effect_size": effect_size,
            "se": 0.3,
            "p_value": 0.5,
            "fdr": 0.9,
        }
    )
    ceg2 = set(marker_genes)
    with caplog.at_level(logging.DEBUG):
        summary = essentiality_contamination_check(df, ceg2, tested_set_size=50)
    assert summary["observed_ceg2_overlap"] == 2
    for gene in marker_genes:
        assert gene not in caplog.text
    assert all(gene not in str(v) for v in summary.values() for gene in marker_genes)


# --- Direction-sanity extraction ----------------------------------------


def test_direction_sanity_missing_gene_raises():
    df = _make_labels_df(n_genes=10, n_significant=1)
    with pytest.raises(ValueError, match="TP53"):
        direction_sanity_extraction(df, ["GENE_00000", "TP53"])


def test_direction_sanity_observed_direction_labels():
    df = pd.DataFrame(
        {
            "gene": ["NEG", "POS", "ZERO"],
            "n_guides": [4, 4, 4],
            "effect_size": [-0.5, 0.5, 0.0],
            "se": [0.1, 0.1, 0.1],
            "p_value": [0.01, 0.01, 0.9],
            "fdr": [0.05, 0.05, 0.95],
        }
    )
    result = direction_sanity_extraction(df, ["NEG", "POS", "ZERO"]).set_index("gene")
    assert result.loc["NEG", "observed_direction"] == "depleted_under_4OHT"
    assert result.loc["POS", "observed_direction"] == "enriched_under_4OHT"
    assert result.loc["ZERO", "observed_direction"] == "no_change"
    assert result.loc["NEG", "sign_interpretation"] == "sensitizes (knockout needed for tamoxifen tolerance)"
    assert result.loc["POS", "sign_interpretation"] == (
        "growth_advantage (knockout confers advantage under tamoxifen)"
    )
    assert result.loc["ZERO", "sign_interpretation"] == "no_directional_call"


def test_direction_sanity_preserves_requested_order():
    df = _make_labels_df(n_genes=10, n_significant=1)
    requested = ["GENE_00007", "GENE_00001", "GENE_00003"]
    result = direction_sanity_extraction(df, requested)
    assert list(result["gene"]) == requested


def test_direction_sanity_t_moderated_is_effect_over_se():
    df = pd.DataFrame(
        {
            "gene": ["G"],
            "n_guides": [4],
            "effect_size": [0.6],
            "se": [0.3],
            "p_value": [0.05],
            "fdr": [0.1],
        }
    )
    result = direction_sanity_extraction(df, ["G"])
    assert result.loc[0, "t_moderated"] == pytest.approx(2.0)


def test_project_direction_sanity_gene_list_excludes_blind_genes():
    """Guard against a future config edit accidentally adding a blind
    control gene to the direction-sanity extraction list.

    This only inspects the static, preregistered gene list in
    config/config.yaml -- never any screen result -- so naming the two
    blind controls here (per PREANALYSIS.md §5) does not inspect or
    report their behavior in the data.
    """
    import yaml

    preregistered_blind_controls = ("RCOR1", "KDM1A")

    with open(Path(__file__).parent.parent / "config" / "config.yaml") as f:
        config = yaml.safe_load(f)
    requested_genes = config["gate1"]["direction_sanity"]["genes"]
    for blind_gene in preregistered_blind_controls:
        assert blind_gene not in requested_genes


# --- CEG2 reference loading: no silent symbol loss ----------------------


def test_load_ceg2_reference_accepts_well_formed_file(tmp_path: Path):
    path = tmp_path / "ceg2.tsv"
    path.write_text("gene_symbol\thgnc_id\nAARS\tHGNC:20\nABCE1\tHGNC:69\n")
    symbols = load_ceg2_reference(path)
    assert symbols == {"AARS", "ABCE1"}


def test_load_ceg2_reference_rejects_duplicate_symbols(tmp_path: Path):
    path = tmp_path / "ceg2.tsv"
    path.write_text("gene_symbol\thgnc_id\nAARS\tHGNC:20\nAARS\tHGNC:20\n")
    with pytest.raises(ValueError, match="duplicate"):
        load_ceg2_reference(path)


def test_load_ceg2_reference_rejects_blank_symbol(tmp_path: Path):
    path = tmp_path / "ceg2.tsv"
    path.write_text("gene_symbol\thgnc_id\n\tHGNC:20\nABCE1\tHGNC:69\n")
    with pytest.raises(ValueError, match="null/blank"):
        load_ceg2_reference(path)


def test_load_ceg2_reference_rejects_missing_column(tmp_path: Path):
    path = tmp_path / "ceg2.tsv"
    path.write_text("symbol\thgnc_id\nAARS\tHGNC:20\n")
    with pytest.raises(ValueError, match="gene_symbol"):
        load_ceg2_reference(path)
