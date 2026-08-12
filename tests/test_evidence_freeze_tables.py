from pathlib import Path

import pandas as pd
import pytest

from src.evidence_freeze_tables import annotate_freeze_columns, build_compact_summary, build_multimodal7_five_layer

REPO_ROOT = Path(__file__).parent.parent
SEVEN = ["USP34", "VEZF1", "CUX1", "DPP9", "LZTR1", "SOX2", "TFAP2C"]


@pytest.fixture(scope="module")
def full_table():
    path = REPO_ROOT / "results" / "tables" / "evidence_freeze" / "final_candidate_evidence.tsv"
    if not path.exists():
        pytest.skip("evidence-freeze full table not generated in this environment")
    return pd.read_csv(path, sep="\t")


class TestFullTableStructure:
    def test_no_duplicate_genes(self, full_table):
        assert not full_table["gene"].duplicated().any()

    def test_contains_all_seven_multimodal_strong_genes(self, full_table):
        assert set(SEVEN).issubset(set(full_table["gene"]))

    def test_contains_all_four_shortlisted_genes(self, full_table):
        assert {"USP34", "VEZF1", "EML5", "CITED2"}.issubset(set(full_table["gene"]))

    def test_only_sensitising_genes_marked_therapeutic_eligible(self, full_table):
        eligible = full_table.loc[full_table["crispr_direction_supports_inhibition_strategy"]]
        assert (eligible["crispr_direction"] == "sensitising_KO").all()
        tolerance = full_table.loc[full_table["crispr_direction"] == "tolerance_associated_KO"]
        assert not tolerance["crispr_direction_supports_inhibition_strategy"].any()

    def test_freeze_shortlisted_is_subset_of_eligible_for_freeze(self, full_table):
        assert "eligible_for_freeze" in full_table.columns, "final_candidate_evidence.tsv must be the post-freeze table (run finalize_tables_with_freeze)"
        shortlisted = full_table.loc[full_table["freeze_shortlisted"]]
        assert shortlisted["eligible_for_freeze"].all()

    def test_freeze_shortlist_rank_is_1_to_4_with_no_gaps(self, full_table):
        ranks = sorted(full_table.loc[full_table["freeze_shortlisted"], "freeze_shortlist_rank"].tolist())
        assert ranks == list(range(1, len(ranks) + 1))

    def test_full_rna_pattern_always_has_divider(self, full_table):
        assert full_table["full_rna_pattern_4"].str.contains(r"\|\|").all()

    def test_resistance_pattern_never_mentions_gse245601(self, full_table):
        assert not full_table["resistance_pattern_3"].str.contains("245601", na=False).any()


class TestMultimodal7FiveLayer:
    def test_exactly_seven_rows(self, full_table):
        out = build_multimodal7_five_layer(full_table, SEVEN)
        assert len(out) == 7
        assert set(out["Gene"]) == set(SEVEN)

    def test_raises_if_not_exactly_seven(self, full_table):
        with pytest.raises(ValueError):
            build_multimodal7_five_layer(full_table, SEVEN + ["NOT_A_REAL_GENE"])

    def test_usp34_and_vezf1_eligible_others_not(self, full_table):
        out = build_multimodal7_five_layer(full_table, SEVEN).set_index("Gene")
        assert out.loc["USP34", "TherapeuticInhibitionEligible"]
        assert out.loc["VEZF1", "TherapeuticInhibitionEligible"]
        for g in ["CUX1", "DPP9", "LZTR1", "SOX2", "TFAP2C"]:
            assert not out.loc[g, "TherapeuticInhibitionEligible"]


class TestAnnotateFreezeColumns:
    """`annotate_freeze_columns` is the ONLY place freeze_shortlisted is
    set -- regression coverage for the Phase 20 Codex review finding that
    figures/summary/source-verification were previously reading a
    different, pre-freeze shortlist list than the one this function
    produces."""

    def test_only_freeze_manifest_genes_get_shortlisted_true(self):
        full_table = pd.DataFrame({"gene": ["A", "B", "C"], "global_rank": [1, 2, 3]})
        freeze_manifest = pd.DataFrame({"gene": ["B"], "freeze_rank": [1]})
        eligibility_audit = pd.DataFrame({"gene": ["A", "B", "C"], "eligible_for_freeze": [False, True, False], "ineligibility_reason": ["x", "", "y"]})
        out = annotate_freeze_columns(full_table, freeze_manifest, eligibility_audit)
        assert out.loc[out["gene"] == "B", "freeze_shortlisted"].iloc[0]
        assert not out.loc[out["gene"] == "A", "freeze_shortlisted"].iloc[0]
        assert not out.loc[out["gene"] == "C", "freeze_shortlisted"].iloc[0]
        assert out.loc[out["gene"] == "B", "freeze_shortlist_rank"].iloc[0] == 1

    def test_shortlisted_genes_sorted_first(self):
        full_table = pd.DataFrame({"gene": ["A", "B", "C"], "global_rank": [1, 2, 3]})
        freeze_manifest = pd.DataFrame({"gene": ["C"], "freeze_rank": [1]})
        eligibility_audit = pd.DataFrame({"gene": ["A", "B", "C"], "eligible_for_freeze": [False, False, True], "ineligibility_reason": ["", "", ""]})
        out = annotate_freeze_columns(full_table, freeze_manifest, eligibility_audit)
        assert out["gene"].tolist()[0] == "C"


class TestCompactSummary:
    def test_freeze_status_labels_are_consistent_with_eligibility(self, full_table):
        summary = build_compact_summary(full_table)
        frozen = summary.loc[summary["FreezeStatus"] == "FROZEN_THERAPEUTIC_SHORTLIST"]
        assert (frozen["CRISPR_Direction"] == "sensitising_KO").all()
        not_eligible = summary.loc[summary["FreezeStatus"] == "not_eligible_wrong_crispr_direction"]
        assert (not_eligible["CRISPR_Direction"] == "tolerance_associated_KO").all()
