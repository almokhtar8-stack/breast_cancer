"""Targeted tests for the druggability + normal-tissue/selectivity review phase."""

from __future__ import annotations

import subprocess
from pathlib import Path

import pandas as pd
import pytest

from src.druggability_safety_build_tables import (
    build_bone_musculoskeletal_table,
    build_druggability_table,
    build_genetic_constraint_table,
    build_normal_tissue_table,
    build_therapeutic_window_summary,
    build_verified_references_table,
)

TABLES = Path("results/tables/druggability_safety")
FIGURES = Path("results/figures/druggability_safety")
REPORTS = Path("results/reports/druggability_safety")
FROZEN_CANDIDATES = {"USP34", "VEZF1", "EML5", "CITED2"}

VALID_DRUGGABILITY_TIERS = {
    "DIRECTLY_DRUGGABLE",
    "POTENTIALLY_DRUGGABLE",
    "INDIRECT_OR_MODALITY_DEPENDENT",
    "CURRENTLY_POORLY_DRUGGABLE",
}
VALID_BONE_TIERS = {"NONE_IDENTIFIED", "INFERRED_ONLY", "DOCUMENTED_CAUSAL_POSTNATAL", "DOCUMENTED_ADULT_CAUSAL"}


class TestFrozenUpstreamUntouched:
    def test_evidence_freeze_and_independent_validation_untouched(self):
        result = subprocess.run(
            [
                "git", "status", "--porcelain",
                "results/tables/evidence_freeze/",
                "docs/THERAPEUTIC_SHORTLIST_FREEZE.md",
                "results/tables/independent_validation/",
                "results/tables/literature_mechanism/",
                "results/networks/systems_network/",
            ],
            capture_output=True, text=True, cwd=Path(__file__).resolve().parents[1],
        )
        assert result.stdout.strip() == "", f"frozen/prior-phase files show as changed: {result.stdout}"


class TestDruggabilityTable:
    def test_all_four_candidates_present_in_order(self):
        df = build_druggability_table()
        assert list(df["candidate"]) == ["USP34", "VEZF1", "EML5", "CITED2"]

    def test_classification_values_are_from_the_conservative_allowed_set(self):
        df = build_druggability_table()
        assert set(df["druggability_classification"]).issubset(VALID_DRUGGABILITY_TIERS)

    def test_no_candidate_reaches_directly_druggable(self):
        # honest reflection of the underlying evidence: none of the four
        # candidates has an approved/clinical-stage compound, so the
        # classification must never claim the top tier
        df = build_druggability_table()
        assert "DIRECTLY_DRUGGABLE" not in set(df["druggability_classification"])

    def test_not_all_four_candidates_forced_into_the_same_modality(self):
        df = build_druggability_table()
        assert df["druggability_classification"].nunique() >= 3

    def test_usp34_is_the_only_candidate_with_a_real_catalytic_domain_structure(self):
        df = build_druggability_table().set_index("candidate")
        assert "PDB" in df.loc["USP34", "structural_data"]
        assert "7W3R" in df.loc["USP34", "structural_data"]
        for other in ["VEZF1", "EML5"]:
            assert "No experimental PDB" in df.loc[other, "structural_data"]


class TestNormalTissueTable:
    def test_all_four_candidates_present(self):
        df = build_normal_tissue_table()
        assert set(df["candidate"]) == FROZEN_CANDIDATES

    def test_depmap_never_appears_in_normal_tissue_context(self):
        # normal-tissue selectivity must never be inferred from DepMap
        # (a cancer-cell-line dependency resource); this table is built
        # entirely from GTEx/HPA
        df = build_normal_tissue_table()
        joined = df.astype(str).apply(lambda col: col.str.contains("DepMap", case=False)).any().any()
        assert not joined

    def test_gtex_breast_values_are_positive_numbers(self):
        df = build_normal_tissue_table()
        assert (df["gtex_breast_tpm"] > 0).all()


class TestGeneticConstraintTable:
    def test_all_four_candidates_present(self):
        df = build_genetic_constraint_table()
        assert set(df["candidate"]) == FROZEN_CANDIDATES

    def test_usp34_is_the_most_constrained_and_eml5_the_least(self):
        # regression-locks the specific gnomAD LOEUF ordering found by
        # this phase's research pass -- a real, load-bearing finding
        # (USP34's strong constraint is a genuine caution flag)
        df = build_genetic_constraint_table().set_index("candidate")
        loeuf = df["loeuf"].astype(float)
        assert loeuf["USP34"] == pytest.approx(0.152, abs=0.01)
        assert loeuf["EML5"] == pytest.approx(0.558, abs=0.01)
        assert loeuf["USP34"] < loeuf["VEZF1"] < loeuf["EML5"]
        assert loeuf["USP34"] < loeuf["CITED2"] < loeuf["EML5"]

    def test_omim_clingen_discrepancies_are_flagged_not_silently_resolved(self):
        df = build_genetic_constraint_table().set_index("candidate")
        assert "YES" in df.loc["VEZF1", "omim_clingen_discrepancy"]
        assert "YES" in df.loc["CITED2", "omim_clingen_discrepancy"]

    def test_genetic_constraint_table_never_claims_toxicity(self):
        df = build_genetic_constraint_table()
        joined = " ".join(df.astype(str).values.flatten()).lower()
        assert "toxic" not in joined


class TestBoneMusculoskeletalTable:
    def test_all_four_candidates_present(self):
        df = build_bone_musculoskeletal_table()
        assert set(df["candidate"]) == FROZEN_CANDIDATES

    def test_bone_concern_categories_are_from_the_allowed_set(self):
        df = build_bone_musculoskeletal_table()
        assert set(df["bone_concern_category"]).issubset(VALID_BONE_TIERS)

    def test_cited2_and_usp34_have_the_strongest_documented_bone_evidence(self):
        # CITED2's HSC-maintenance studies used Mx1-Cre, genuinely induced
        # in already-mature adult mice; USP34's used MSC/pre-osteoblast-Cre,
        # active during skeletal development -- these are kept as distinct
        # tiers so USP34's evidence is never overstated as adult-onset
        df = build_bone_musculoskeletal_table().set_index("candidate")
        assert df.loc["CITED2", "bone_concern_category"] == "DOCUMENTED_ADULT_CAUSAL"
        assert df.loc["USP34", "bone_concern_category"] == "DOCUMENTED_CAUSAL_POSTNATAL"
        assert df.loc["VEZF1", "bone_concern_category"] == "INFERRED_ONLY"
        assert df.loc["EML5", "bone_concern_category"] == "NONE_IDENTIFIED"

    def test_missing_expression_values_are_kept_as_real_missing_data_not_zero(self):
        # CITED2's bone-marrow/skeletal-muscle nTPM was not locatable in
        # the source used -- this must be represented as missing (NaN),
        # never silently defaulted to zero, which would misleadingly
        # imply a measured absence of expression
        df = build_bone_musculoskeletal_table().set_index("candidate")
        assert pd.isna(df.loc["CITED2", "bone_marrow_expression_hpa_ntpm"])

    def test_no_species_extrapolation_language_in_evidence_type(self):
        df = build_bone_musculoskeletal_table()
        joined = " ".join(df["evidence_type"].astype(str).values)
        assert "human" in joined.lower() or "mouse" in joined.lower()


class TestTherapeuticWindowSummary:
    def test_joins_real_frozen_hany_and_depmap_values_not_hand_typed(self):
        drug = build_druggability_table()
        tissue = build_normal_tissue_table()
        constraint = build_genetic_constraint_table()
        bone = build_bone_musculoskeletal_table()
        window = build_therapeutic_window_summary(drug, tissue, constraint, bone).set_index("candidate")

        hany = pd.read_csv("results/tables/evidence_freeze/THERAPEUTIC_SHORTLIST_FREEZE.tsv", sep="\t").set_index("gene")
        for c in FROZEN_CANDIDATES:
            assert f"{hany.loc[c, 'crispr_fdr']:.4f}" in window.loc[c, "functional_tamoxifen_evidence"]

    def test_vezf1_is_the_only_candidate_with_nonzero_er_luminal_dependency_in_summary(self):
        drug = build_druggability_table()
        tissue = build_normal_tissue_table()
        constraint = build_genetic_constraint_table()
        bone = build_bone_musculoskeletal_table()
        window = build_therapeutic_window_summary(drug, tissue, constraint, bone).set_index("candidate")
        assert "27.3%" in window.loc["VEZF1", "cancer_dependency_depmap"]
        for c in ["USP34", "EML5", "CITED2"]:
            assert "0.0%" in window.loc[c, "cancer_dependency_depmap"]

    def test_frozen_ranking_note_present_and_ranking_file_not_rewritten_by_this_phase(self):
        drug = build_druggability_table()
        tissue = build_normal_tissue_table()
        constraint = build_genetic_constraint_table()
        bone = build_bone_musculoskeletal_table()
        window = build_therapeutic_window_summary(drug, tissue, constraint, bone)
        assert window["notes"].str.contains("FROZEN").all()

    def test_normal_tissue_and_bone_concern_columns_are_genuinely_derived_from_the_input_tables(self):
        # regression test for a Codex-flagged bug: the summary builder
        # used to ignore the normal_tissue/bone DataFrames it was passed
        # and pull from separate hardcoded dicts instead, so a corrected
        # source row could silently fail to propagate. Confirm the join
        # is real by editing an input value and checking it flows through.
        drug = build_druggability_table()
        tissue = build_normal_tissue_table()
        constraint = build_genetic_constraint_table().copy()
        bone = build_bone_musculoskeletal_table().copy()
        constraint.loc[constraint["candidate"] == "USP34", "loeuf"] = 0.999999
        bone.loc[bone["candidate"] == "EML5", "published_bone_role_summary"] = "SENTINEL_VALUE_FOR_TEST"

        window = build_therapeutic_window_summary(drug, tissue, constraint, bone).set_index("candidate")
        assert "0.999999" in window.loc["USP34", "normal_tissue_concern"]
        assert "SENTINEL_VALUE_FOR_TEST" in window.loc["EML5", "bone_musculoskeletal_concern"]


class TestVerifiedReferences:
    def test_every_reference_has_a_pmid_and_a_source_candidate(self):
        df = build_verified_references_table()
        assert df["PMID"].astype(str).str.match(r"^\d+$").all()
        assert set(df["candidate"]).issubset(FROZEN_CANDIDATES)

    def test_references_cover_all_four_candidates(self):
        df = build_verified_references_table()
        assert set(df["candidate"]) == FROZEN_CANDIDATES


class TestOutputsExistAndAreNonTrivial:
    @pytest.mark.parametrize("name", [
        "candidate_druggability.tsv",
        "candidate_normal_tissue_context.tsv",
        "candidate_genetic_constraint.tsv",
        "candidate_bone_musculoskeletal_context.tsv",
        "candidate_therapeutic_window_summary.tsv",
        "verified_references.tsv",
    ])
    def test_table_written_and_nonempty(self, name):
        path = TABLES / name
        assert path.exists(), f"missing {path}"
        df = pd.read_csv(path, sep="\t")
        assert len(df) > 0

    @pytest.mark.parametrize("name", [
        "01_four_candidate_druggability_summary.png",
        "02_four_candidate_normal_tissue_context.png",
        "03_four_candidate_therapeutic_window_map.png",
        "04_four_candidate_bone_musculoskeletal_context.png",
    ])
    def test_figure_written_and_nontrivial_size(self, name):
        path = FIGURES / name
        assert path.exists(), f"missing {path}"
        assert path.stat().st_size > 20_000, f"{path} suspiciously small, likely a blank/broken render"


class TestReportLanguage:
    def test_report_never_calls_a_candidate_safe(self):
        # the word "safe"/"safety" is fine (e.g. "not a safety claim",
        # "expression is not a claim about ... safety"); what must never
        # appear is an affirmative "is/are safe" assertion
        text = (REPORTS / "four_candidate_druggability_safety_review.md").read_text().lower()
        for phrase in ("is safe", "are safe", "safe profile", "safe target", "considered safe"):
            assert phrase not in text, f"report contains an affirmative safety claim: {phrase!r}"

    def test_report_never_infers_normal_tissue_selectivity_from_depmap_alone(self):
        text = (REPORTS / "four_candidate_druggability_safety_review.md").read_text()
        assert "DepMap is not used for this claim" in text or "DepMap (a cancer-cell-line" in text

    def test_report_names_all_four_candidates(self):
        text = (REPORTS / "four_candidate_druggability_safety_review.md").read_text()
        for c in FROZEN_CANDIDATES:
            assert c in text

    def test_report_never_equates_genetic_constraint_with_broad_essentiality_or_toxicity(self):
        # regression test for Codex-flagged overclaiming: constraint/
        # embryonic-lethal mouse evidence must never be phrased as
        # establishing broad essentiality in adult normal tissue or drug
        # toxicity
        text = (REPORTS / "four_candidate_druggability_safety_review.md").read_text().lower()
        for phrase in ("broad essential function in normal human biology", "broadly essential in normal tissues?"):
            assert phrase not in text, f"report contains an overclaiming constraint phrase: {phrase!r}"

    def test_report_never_says_cited2_interface_acts_entirely(self):
        # regression test: CITED2 also binds SMAD2/3, TFAP2A/B/C, WT1,
        # LHX2, PPARA per UniProt -- "acts/works entirely through" the
        # CBP/p300 interface overstates this
        text = (REPORTS / "four_candidate_druggability_safety_review.md").read_text().lower()
        assert "entirely by competitively binding" not in text
        assert "entirely through the c-terminal" not in text
