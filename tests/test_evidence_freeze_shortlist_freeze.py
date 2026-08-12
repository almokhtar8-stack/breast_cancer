import pandas as pd
import pytest

from src.evidence_freeze_shortlist_freeze import _crispr_band, _sample_robustness_band, build_freeze, determine_eligibility, rank_eligible


def _pool(rows: list[dict]) -> pd.DataFrame:
    defaults = {
        "crispr_direction": "sensitising_KO", "crispr_fdr": 0.05, "crispr_effect": -1.0,
        "resistance_fdr05_count": 0, "resistance_direction_consistency": "insufficient",
        "human_tumor_support": "not_significant", "ranking_stability": "DATASET_DEPENDENT",
        "resistance_pattern_3": "NA | NA | NA", "full_rna_pattern_4": "NA | NA | NA || NA",
        "gse111151_cell_line_consistency": "not_computed", "main_strength": "", "main_weakness": "",
    }
    out = []
    for r in rows:
        row = dict(defaults)
        row.update(r)
        out.append(row)
    return pd.DataFrame(out)


class TestCrisprBand:
    def test_bands_by_fdr(self):
        assert _crispr_band(0.005) == "VERY_STRONG"
        assert _crispr_band(0.03) == "STRONG"
        assert _crispr_band(0.08) == "MODERATE"
        assert _crispr_band(0.2) == "WEAK"
        assert _crispr_band(0.5) is None
        assert _crispr_band(float("nan")) is None


class TestSampleRobustnessBand:
    def test_high_medium_low_unknown(self):
        assert _sample_robustness_band("4/4") == "HIGH"
        assert _sample_robustness_band("3/4") == "MEDIUM"
        assert _sample_robustness_band("1/4") == "LOW"
        assert _sample_robustness_band("not_computed") == "UNKNOWN"
        assert _sample_robustness_band(None) == "UNKNOWN"


class TestRankEligibleCriterionOrder:
    """Regression tests for the Phase 20 Codex review finding: the sort
    previously checked human-tumor evidence BEFORE resistance direction
    consistency, and omitted sample robustness entirely, diverging from
    the declared hierarchy (CRISPR -> resistance -> consistency -> human
    -> sample robustness -> stability -> gene)."""

    def test_consistency_beats_human_evidence_when_crispr_and_resistance_tie(self):
        df = _pool(
            [
                {"gene": "CONSISTENT_NO_HUMAN", "crispr_fdr": 0.04, "resistance_fdr05_count": 1, "resistance_direction_consistency": "all_up", "human_tumor_support": "not_significant"},
                {"gene": "INCONSISTENT_WITH_HUMAN", "crispr_fdr": 0.04, "resistance_fdr05_count": 1, "resistance_direction_consistency": "majority_up", "human_tumor_support": "significant"},
            ]
        )
        eligible = determine_eligibility(df)
        ranked = rank_eligible(eligible.loc[eligible["eligible_for_freeze"]])
        assert ranked["gene"].tolist() == ["CONSISTENT_NO_HUMAN", "INCONSISTENT_WITH_HUMAN"]

    def test_sample_robustness_breaks_a_tie_after_human_evidence(self):
        df = _pool(
            [
                {"gene": "ROBUST_SAMPLE", "crispr_fdr": 0.04, "resistance_fdr05_count": 1, "resistance_direction_consistency": "all_up", "human_tumor_support": "significant", "gse111151_cell_line_consistency": "4/4"},
                {"gene": "WEAK_SAMPLE", "crispr_fdr": 0.04, "resistance_fdr05_count": 1, "resistance_direction_consistency": "all_up", "human_tumor_support": "significant", "gse111151_cell_line_consistency": "1/4"},
            ]
        )
        eligible = determine_eligibility(df)
        ranked = rank_eligible(eligible.loc[eligible["eligible_for_freeze"]])
        assert ranked["gene"].tolist() == ["ROBUST_SAMPLE", "WEAK_SAMPLE"]


class TestEligibilityGate:
    def test_tolerance_direction_never_eligible_even_with_strong_everything(self):
        df = _pool([{"gene": "TOL_STRONG", "crispr_direction": "tolerance_associated_KO", "crispr_fdr": 0.001, "resistance_fdr05_count": 3, "resistance_direction_consistency": "all_up", "human_tumor_support": "significant"}])
        audited = determine_eligibility(df)
        assert not audited.loc[0, "eligible_for_freeze"]
        assert "tolerance" in audited.loc[0, "ineligibility_reason"].lower()

    def test_meaningless_sign_with_no_real_crispr_evidence_excluded(self):
        """A gene with a negative CRISPR effect but FDR~0.8 (pure noise)
        must not qualify merely because the sign happens to be negative --
        regression test for the same class of bug caught during candidate
        adjudication (the GREB1 case)."""
        df = _pool([{"gene": "NOISY", "crispr_direction": "sensitising_KO", "crispr_fdr": 0.82, "resistance_fdr05_count": 2, "resistance_direction_consistency": "all_up", "human_tumor_support": "significant"}])
        audited = determine_eligibility(df)
        assert not audited.loc[0, "eligible_for_freeze"]
        assert "no real evidence" in audited.loc[0, "ineligibility_reason"].lower()

    def test_sensitising_real_crispr_but_zero_rna_support_excluded_from_list_a(self):
        df = _pool([{"gene": "FUNCONLY", "crispr_direction": "sensitising_KO", "crispr_fdr": 0.001, "resistance_fdr05_count": 0, "resistance_direction_consistency": "insufficient", "human_tumor_support": "not_significant"}])
        audited = determine_eligibility(df)
        assert not audited.loc[0, "eligible_for_freeze"]
        assert "functional-only" in audited.loc[0, "ineligibility_reason"].lower()

    def test_sensitising_real_crispr_with_resistance_support_is_eligible(self):
        df = _pool([{"gene": "GOOD", "crispr_direction": "sensitising_KO", "crispr_fdr": 0.04, "resistance_fdr05_count": 1, "resistance_direction_consistency": "all_up", "human_tumor_support": "not_significant"}])
        audited = determine_eligibility(df)
        assert audited.loc[0, "eligible_for_freeze"]

    def test_sensitising_real_crispr_with_only_human_support_is_eligible(self):
        df = _pool([{"gene": "GOOD_HUMAN", "crispr_direction": "sensitising_KO", "crispr_fdr": 0.04, "resistance_fdr05_count": 0, "resistance_direction_consistency": "insufficient", "human_tumor_support": "significant"}])
        audited = determine_eligibility(df)
        assert audited.loc[0, "eligible_for_freeze"]


class TestFreezeReproduction:
    def test_reproduces_real_frozen_four_genes(self):
        try:
            full_table = pd.read_csv("results/tables/evidence_freeze/final_candidate_evidence.tsv", sep="\t")
        except FileNotFoundError:
            pytest.skip("evidence-freeze full table not generated in this environment")
        frozen, _ = build_freeze(full_table)
        # USP34 ranks ahead of VEZF1: both tie on CRISPR-strength and resistance-evidence
        # bands, but USP34's resistance direction is fully concordant (all_up) while
        # VEZF1's is only majority_up -- consistency is checked before human evidence
        assert frozen["gene"].tolist() == ["USP34", "VEZF1", "EML5", "CITED2"]

    def test_freeze_is_deterministic(self):
        try:
            full_table = pd.read_csv("results/tables/evidence_freeze/final_candidate_evidence.tsv", sep="\t")
        except FileNotFoundError:
            pytest.skip("evidence-freeze full table not generated in this environment")
        f1, _ = build_freeze(full_table)
        f2, _ = build_freeze(full_table)
        pd.testing.assert_frame_equal(f1, f2)

    def test_all_frozen_genes_pass_sensitising_gate(self):
        try:
            full_table = pd.read_csv("results/tables/evidence_freeze/final_candidate_evidence.tsv", sep="\t")
        except FileNotFoundError:
            pytest.skip("evidence-freeze full table not generated in this environment")
        frozen, _ = build_freeze(full_table)
        assert (frozen["crispr_direction"] == "sensitising_KO").all()

    def test_freeze_size_within_2_to_5_and_not_manufactured(self):
        try:
            full_table = pd.read_csv("results/tables/evidence_freeze/final_candidate_evidence.tsv", sep="\t")
        except FileNotFoundError:
            pytest.skip("evidence-freeze full table not generated in this environment")
        frozen, audited = build_freeze(full_table, max_n=5)
        n_eligible = int(audited["eligible_for_freeze"].sum())
        assert 2 <= len(frozen) <= 5
        assert len(frozen) == min(n_eligible, 5)  # never padded beyond the genuinely eligible count
