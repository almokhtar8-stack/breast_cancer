"""Targeted tests for the literature/mechanism review phase."""

from __future__ import annotations

import subprocess
from pathlib import Path

import pandas as pd
import pytest

TABLES = Path("results/tables/literature_mechanism")
REPORTS = Path("results/reports/literature_mechanism")
FIGURES = Path("results/figures/literature_mechanism")
FROZEN_CANDIDATES = {"USP34", "VEZF1", "EML5", "CITED2"}


class TestFrozenShortlistUntouched:
    def test_all_four_candidates_present_in_claim_table(self):
        df = pd.read_csv(TABLES / "four_candidate_claim_evidence.tsv", sep="\t")
        assert set(df["candidate"]) == FROZEN_CANDIDATES

    def test_evidence_freeze_and_systems_network_untouched(self):
        result = subprocess.run(
            ["git", "status", "--porcelain", "results/tables/evidence_freeze/", "docs/THERAPEUTIC_SHORTLIST_FREEZE.md", "results/networks/systems_network/", "results/tables/systems_network/four_candidate_network_audit.tsv"],
            capture_output=True,
            text=True,
            cwd=Path(__file__).resolve().parents[1],
        )
        assert result.stdout.strip() == "", f"frozen files show as changed: {result.stdout}"


class TestClaimEvidenceTable:
    def test_evidence_levels_are_valid(self):
        df = pd.read_csv(TABLES / "four_candidate_claim_evidence.tsv", sep="\t")
        assert df["evidence_level"].isin([1, 2, 3, 4, 5]).all()

    def test_claim_ids_unique(self):
        df = pd.read_csv(TABLES / "four_candidate_claim_evidence.tsv", sep="\t")
        assert df["claim_id"].is_unique

    def test_review_or_primary_column_is_fully_populated_and_reviews_never_masquerade_as_primary(self):
        df = pd.read_csv(TABLES / "four_candidate_claim_evidence.tsv", sep="\t")
        assert df["review_or_primary"].notna().all()
        # every claim actually used as evidence in this phase is a primary
        # source or an explicitly-labeled database resource -- no row may
        # claim "review" status while also being cited as if it ran the
        # experiment (review_or_primary containing "review" would require
        # an accompanying orientation-only disclaimer, checked here by
        # requiring such rows to carry evidence_level 5)
        reviews = df.loc[df["review_or_primary"].str.contains("review", case=False, na=False)]
        assert (reviews["evidence_level"] == 5).all()

    def test_pathway_context_claims_never_counted_as_candidate_specific(self):
        df = pd.read_csv(TABLES / "four_candidate_claim_evidence.tsv", sep="\t")
        pathway_context = df.loc[df["claim_id"].str.contains("PATHWAY-CONTEXT")]
        assert len(pathway_context) >= 1
        assert (~pathway_context["is_candidate_specific"]).all()

    def test_eml5_has_no_claim_supporting_project_hypothesis(self):
        # EML5's literature is either null-result or associative-only in a
        # different tissue -- none of it should be marked as supporting a
        # breast-cancer resistance mechanism for EML5
        df = pd.read_csv(TABLES / "four_candidate_claim_evidence.tsv", sep="\t")
        eml5 = df.loc[df["candidate"] == "EML5"]
        assert not eml5["supports_project_hypothesis"].any()

    def test_at_least_one_contradiction_flagged_for_cited2_and_usp34(self):
        df = pd.read_csv(TABLES / "four_candidate_claim_evidence.tsv", sep="\t")
        for candidate in ["CITED2", "USP34"]:
            sub = df.loc[df["candidate"] == candidate]
            assert sub["contradicts_project_hypothesis"].any(), f"{candidate} should have >=1 documented contradiction"

    def test_no_claim_both_supports_and_contradicts(self):
        df = pd.read_csv(TABLES / "four_candidate_claim_evidence.tsv", sep="\t")
        assert not (df["supports_project_hypothesis"] & df["contradicts_project_hypothesis"]).any()

    def test_tamoxifen_specific_requires_breast_cancer_specific(self):
        # a claim cannot be tamoxifen-specific without being in a breast
        # cancer system (tamoxifen is a breast-cancer/ER+ drug in this
        # project's context)
        df = pd.read_csv(TABLES / "four_candidate_claim_evidence.tsv", sep="\t")
        tam = df.loc[df["tamoxifen_specific"]]
        assert tam["breast_cancer_specific"].all()


class TestVerifiedReferences:
    def test_every_reference_has_a_pmid_or_doi_or_is_explicitly_a_database_resource(self):
        df = pd.read_csv(REPORTS / "verified_references.tsv", sep="\t")
        has_id = df["PMID"].notna() | df["DOI"].notna()
        assert has_id.all(), "every named reference must carry a PMID or DOI"

    def test_verification_status_is_disclosed_for_every_reference(self):
        df = pd.read_csv(REPORTS / "verified_references.tsv", sep="\t")
        assert df["verification_status"].notna().all()
        assert df["verification_status"].isin(
            [
                "full text fetched and quote-verified",
                "bibliographic identity verified (PMID/DOI/journal/year); narrative content from secondary aggregation, primary full text paywalled/blocked",
            ]
        ).all()

    def test_no_duplicate_references(self):
        df = pd.read_csv(REPORTS / "verified_references.tsv", sep="\t")
        assert not df.duplicated(subset=["paper_title", "PMID"]).any()


class TestLiteratureComparisonTable:
    def test_all_four_candidates_present(self):
        df = pd.read_csv(TABLES / "four_candidate_literature_comparison.tsv", sep="\t")
        assert set(df["candidate"]) == FROZEN_CANDIDATES

    def test_usp34_has_zero_direct_tamoxifen_papers_after_excluding_pathway_context(self):
        df = pd.read_csv(TABLES / "four_candidate_literature_comparison.tsv", sep="\t").set_index("candidate")
        assert df.loc["USP34", "n_direct_tamoxifen_papers"] == 0
        assert "PATHWAY-CONTEXT" in df.loc["USP34", "pathway_level_context_excluded_from_counts"]

    def test_cited2_has_the_most_breast_cancer_mechanistic_papers(self):
        df = pd.read_csv(TABLES / "four_candidate_literature_comparison.tsv", sep="\t").set_index("candidate")
        assert df.loc["CITED2", "n_breast_cancer_mechanistic_papers"] == df["n_breast_cancer_mechanistic_papers"].max()

    def test_eml5_has_zero_papers_in_every_evidence_category(self):
        df = pd.read_csv(TABLES / "four_candidate_literature_comparison.tsv", sep="\t").set_index("candidate")
        count_cols = [c for c in df.columns if c.startswith("n_")]
        assert (df.loc["EML5", count_cols] == 0).all()


class TestFigureAndReportExist:
    def test_figure_exists(self):
        assert (FIGURES / "01_candidate_mechanism_evidence_map.png").exists()

    def test_report_exists_and_mentions_all_candidates(self):
        path = REPORTS / "four_candidate_mechanism_review.md"
        assert path.exists()
        text = path.read_text()
        for candidate in FROZEN_CANDIDATES:
            assert candidate in text

    def test_report_documents_no_upstream_rerun(self):
        text = " ".join((REPORTS / "four_candidate_mechanism_review.md").read_text().lower().split())
        assert "no upstream" in text


if __name__ == "__main__":
    import sys

    sys.exit(pytest.main([__file__, "-v"]))
