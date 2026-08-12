import pandas as pd

from src.candidate_adjudication_axes import (
    ARCHETYPE_ORDER,
    assign_archetype,
    classify_axis_a_functional,
    classify_axis_b_resistance,
    classify_axis_c_human,
)


def _row(**kwargs):
    base = {
        "crispr_fdr": float("nan"), "crispr_evidence_percentile": float("nan"),
        "resistance_direction_consensus": None, "resistance_fdr05_count": float("nan"), "resistance_median_percentile": float("nan"),
        "gse245601_epi_fdr": float("nan"), "gse245601_malignant_fdr": float("nan"), "gse240112_tumor_fdr": float("nan"),
        "gse245601_evidence_percentile": float("nan"), "gse240112_evidence_percentile": float("nan"),
    }
    base.update(kwargs)
    return pd.Series(base)


class TestAxisA:
    def test_very_strong_below_001(self):
        assert classify_axis_a_functional(_row(crispr_fdr=0.005)) == "VERY_STRONG"

    def test_strong_below_005(self):
        assert classify_axis_a_functional(_row(crispr_fdr=0.03)) == "STRONG"

    def test_moderate_below_010(self):
        assert classify_axis_a_functional(_row(crispr_fdr=0.08)) == "MODERATE"

    def test_weak_needs_high_percentile(self):
        assert classify_axis_a_functional(_row(crispr_fdr=0.15, crispr_evidence_percentile=0.95)) == "WEAK"
        assert classify_axis_a_functional(_row(crispr_fdr=0.15, crispr_evidence_percentile=0.5)) == "NO_EVIDENCE"

    def test_untestable_is_no_evidence(self):
        assert classify_axis_a_functional(_row()) == "NO_EVIDENCE"


class TestAxisB:
    def test_very_strong_needs_two_fdr_and_full_consensus(self):
        assert classify_axis_b_resistance(_row(resistance_direction_consensus="all_up", resistance_fdr05_count=2)) == "VERY_STRONG"

    def test_strong_one_fdr_and_majority_consensus(self):
        assert classify_axis_b_resistance(_row(resistance_direction_consensus="majority_down", resistance_fdr05_count=1)) == "STRONG"

    def test_mixed_is_discordant_even_with_high_percentile(self):
        assert classify_axis_b_resistance(_row(resistance_direction_consensus="mixed", resistance_median_percentile=0.99)) == "DISCORDANT"

    def test_insufficient_is_no_evidence(self):
        assert classify_axis_b_resistance(_row(resistance_direction_consensus="insufficient")) == "NO_EVIDENCE"


class TestAxisC:
    def test_very_strong_needs_both_acute_and_recurrence(self):
        assert classify_axis_c_human(_row(gse245601_epi_fdr=0.01, gse240112_tumor_fdr=0.01)) == "VERY_STRONG"

    def test_strong_needs_only_one(self):
        assert classify_axis_c_human(_row(gse245601_epi_fdr=0.01, gse240112_tumor_fdr=0.5)) == "STRONG"
        assert classify_axis_c_human(_row(gse245601_epi_fdr=0.5, gse240112_tumor_fdr=0.01)) == "STRONG"

    def test_no_significance_and_low_percentile_is_no_evidence(self):
        assert classify_axis_c_human(_row(gse245601_epi_fdr=0.5, gse240112_tumor_fdr=0.5, gse245601_evidence_percentile=0.3, gse240112_evidence_percentile=0.2)) == "NO_EVIDENCE"


class TestArchetypeAssignment:
    def test_convergence_requires_both_functional_and_resistance_strong(self):
        row = pd.Series({"axis_a_functional": "STRONG", "axis_b_resistance": "STRONG", "axis_c_human": "NO_EVIDENCE"})
        assert assign_archetype(row) == "A_FUNCTIONAL_RESISTANCE_CONVERGENCE"

    def test_functional_only_when_resistance_and_human_absent(self):
        row = pd.Series({"axis_a_functional": "STRONG", "axis_b_resistance": "NO_EVIDENCE", "axis_c_human": "NO_EVIDENCE"})
        assert assign_archetype(row) == "C_FUNCTIONAL_ONLY"

    def test_resistance_biomarker_when_functional_absent(self):
        row = pd.Series({"axis_a_functional": "NO_EVIDENCE", "axis_b_resistance": "VERY_STRONG", "axis_c_human": "NO_EVIDENCE"})
        assert assign_archetype(row) == "D_RESISTANCE_BIOMARKER_PATHWAY"

    def test_low_evidence_fallback(self):
        row = pd.Series({"axis_a_functional": "WEAK", "axis_b_resistance": "WEAK", "axis_c_human": "WEAK", "human_evidence_sources": "none_FDR<0.05"})
        assert assign_archetype(row) == "H_LOW_INSUFFICIENT_EVIDENCE"

    def test_every_archetype_name_in_declared_order(self):
        for a in ARCHETYPE_ORDER:
            assert a[0] in "ABCDEFGH"

    def test_context_dependent_is_reachable_discordant_plus_functional(self):
        """Regression test for the Phase 34 Codex review finding: with the
        original precedence order (G checked after B/C/E/F), a discordant
        -resistance gene with strong functional evidence was always caught
        by C_FUNCTIONAL_ONLY first, making G_CONTEXT_DEPENDENT dead code.
        G must now be checked immediately after A."""
        row = pd.Series({"axis_a_functional": "STRONG", "axis_b_resistance": "DISCORDANT", "axis_c_human": "NO_EVIDENCE", "human_evidence_sources": "none_FDR<0.05"})
        assert assign_archetype(row) == "G_CONTEXT_DEPENDENT"

    def test_context_dependent_is_reachable_discordant_plus_human(self):
        row = pd.Series({"axis_a_functional": "NO_EVIDENCE", "axis_b_resistance": "DISCORDANT", "axis_c_human": "STRONG", "human_evidence_sources": "GSE240112_recurrence"})
        assert assign_archetype(row) == "G_CONTEXT_DEPENDENT"

    def test_discordant_with_no_other_evidence_is_low_evidence_not_context_dependent(self):
        row = pd.Series({"axis_a_functional": "NO_EVIDENCE", "axis_b_resistance": "DISCORDANT", "axis_c_human": "NO_EVIDENCE", "human_evidence_sources": "none_FDR<0.05"})
        assert assign_archetype(row) == "H_LOW_INSUFFICIENT_EVIDENCE"
