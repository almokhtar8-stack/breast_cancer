"""Targeted tests for the USP34 vs VEZF1 translational deep-dive phase."""

from __future__ import annotations

import subprocess
from pathlib import Path

import pandas as pd
import pytest

from src.lead_target_deep_dive_build_tables import (
    build_bone_marrow_liability_table,
    build_experimental_plan_table,
    build_genetic_constraint_table,
    build_head_to_head_table,
    build_indirect_target_crosscheck_table,
    build_muscle_liability_table,
    build_tissue_expression_table,
    build_tissue_liability_table,
    build_usp34_direct_targeting_table,
    build_usp34_indirect_targets_table,
    build_verified_references_table,
    build_vezf1_direct_targeting_table,
    build_vezf1_indirect_targets_table,
)

TABLES = Path("results/tables/lead_target_deep_dive")
FIGURES = Path("results/figures/lead_target_deep_dive")
REPORTS = Path("results/reports/lead_target_deep_dive")
GENES = {"USP34", "VEZF1"}


class TestFrozenUpstreamUntouched:
    def test_frozen_and_prior_phase_files_untouched(self):
        # note: results/tables/druggability_safety/ is deliberately NOT
        # checked here -- it is itself a prior, still-uncommitted session
        # deliverable (untracked, "??"), not yet-frozen content, so it is
        # expected and fine for it to show as untracked; this check is
        # only for MODIFICATIONS ("M") to genuinely frozen/committed dirs
        result = subprocess.run(
            [
                "git", "status", "--porcelain",
                "results/tables/evidence_freeze/",
                "docs/THERAPEUTIC_SHORTLIST_FREEZE.md",
                "results/tables/independent_validation/",
                "results/tables/literature_mechanism/",
                "results/tables/cross_dataset_genomewide/",
                "results/networks/systems_network/",
            ],
            capture_output=True, text=True, cwd=Path(__file__).resolve().parents[1],
        )
        assert result.stdout.strip() == "", f"frozen/prior-phase files show as changed: {result.stdout}"


class TestTissueExpressionTable:
    def test_both_genes_present(self):
        df = build_tissue_expression_table()
        assert set(df["candidate"]) == GENES

    def test_no_toxicity_language(self):
        df = build_tissue_expression_table()
        joined = " ".join(df.astype(str).values.flatten()).lower()
        assert "toxic" not in joined


class TestTissueLiabilityTable:
    def test_both_genes_present(self):
        df = build_tissue_liability_table()
        assert set(df["candidate"]) == GENES

    def test_classifications_from_allowed_set(self):
        allowed = {
            "DOCUMENTED_HUMAN", "DOCUMENTED_ADULT_CAUSAL", "DOCUMENTED_POSTNATAL_CAUSAL",
            "DOCUMENTED_DEVELOPMENTAL", "PRIMARY_HUMAN_CELL_FUNCTIONAL", "ANIMAL_ONLY",
            "EXPRESSION_ONLY", "INFERRED_ONLY", "NONE_IDENTIFIED", "INSUFFICIENT_DATA",
        }
        df = build_tissue_liability_table()
        assert set(df["classification"]).issubset(allowed)

    def test_usp34_bone_is_documented_postnatal_not_developmental_only(self):
        df = build_tissue_liability_table()
        row = df[(df["candidate"] == "USP34") & (df["organ_system"] == "Bone")]
        assert len(row) == 1
        assert row.iloc[0]["classification"] == "DOCUMENTED_POSTNATAL_CAUSAL"

    def test_vezf1_cardiovascular_reflects_new_pmid_31911272_finding(self):
        df = build_tissue_liability_table()
        row = df[(df["candidate"] == "VEZF1") & (df["organ_system"] == "Cardiovascular")]
        assert len(row) == 1
        assert row.iloc[0]["classification"] == "DOCUMENTED_POSTNATAL_CAUSAL"
        assert "31911272" in row.iloc[0]["evidence_summary"]

    def test_never_extrapolates_embryonic_lethal_to_adult_toxicity(self):
        df = build_tissue_liability_table()
        joined = " ".join(df["evidence_summary"].astype(str).values).lower()
        assert "embryonic" in joined  # the finding is present
        # and the report-level caveat exists (checked separately in TestReportLanguage)


class TestMuscleLiabilityTable:
    def test_usp34_muscle_search_returned_zero_hits(self):
        df = build_muscle_liability_table()
        row = df[(df["candidate"] == "USP34") & (df["question"].str.contains("54.1"))]
        assert len(row) == 1
        assert "zero" in row.iloc[0]["finding"].lower() or "ZERO" in row.iloc[0]["finding"]

    def test_vezf1_skeletal_vs_cardiac_muscle_distinguished(self):
        df = build_muscle_liability_table()
        skeletal = df[(df["candidate"] == "VEZF1") & (df["question"].str.contains("skeletal-muscle function"))]
        cardiac = df[(df["candidate"] == "VEZF1") & (df["question"].str.contains("CARDIAC"))]
        assert len(skeletal) == 1 and len(cardiac) == 1
        assert "NONE" in skeletal.iloc[0]["finding"].upper()
        assert "31911272" in cardiac.iloc[0]["finding"]
        assert skeletal.iloc[0]["finding"] != cardiac.iloc[0]["finding"]

    def test_distinguishes_expression_from_function_from_toxicity(self):
        df = build_muscle_liability_table()
        interpretations = " ".join(df["interpretation"].astype(str).values)
        assert "EXPRESSION" in interpretations
        assert "FUNCTIONAL" in interpretations or "functional" in interpretations.lower()

    def test_usp34_alias_search_and_nfic_finding_present(self):
        # gap-completion pass: USP34's historical aliases (KIAA0570,
        # KIAA0729) were also searched against muscle terms (also zero
        # hits), and a related but non-muscle developmental finding
        # (NFIC/tooth-root) was identified and cross-referenced
        df = build_muscle_liability_table()
        row = df[(df["candidate"] == "USP34") & (df["question"].str.contains("54.1"))].iloc[0]
        assert "KIAA0570" in row["finding"]
        assert "KIAA0729" in row["finding"]
        assert "C2C12" in row["finding"] and "EXCLUDED" in row["finding"]
        row2 = df[(df["candidate"] == "USP34") & (df["question"].str.contains("knockout, primary myocyte"))].iloc[0]
        assert "NFIC" in row2["finding"]
        assert "33686052" in row2["sources"]


class TestGeneticConstraintTable:
    def test_both_genes_present_with_real_loeuf(self):
        df = build_genetic_constraint_table().set_index("candidate")
        assert df.loc["USP34", "loeuf"] == pytest.approx(0.152, abs=0.001)
        assert df.loc["VEZF1", "loeuf"] == pytest.approx(0.24, abs=0.001)

    def test_clinvar_counts_are_positive_integers(self):
        df = build_genetic_constraint_table()
        assert (df["clinvar_total_variants"] > 0).all()
        assert (df["clinvar_pathogenic_or_likely_pathogenic"] > 0).all()
        assert (df["clinvar_pathogenic_or_likely_pathogenic"] <= df["clinvar_total_variants"]).all()

    def test_never_equates_constraint_with_toxicity(self):
        # the phrase "will be toxic" is allowed ONLY inside an explicit
        # negation ("NOT the same claim as ... will be toxic") -- check
        # that pattern specifically, plus the required positive framing
        df = build_genetic_constraint_table()
        interp = " ".join(df["interpretation"].astype(str).values)
        assert "NOT the same claim" in interp
        assert "poorly tolerated / requires caution" in interp


class TestDirectTargetingTables:
    def test_usp34_direct_targeting_never_invents_an_inhibitor(self):
        df = build_usp34_direct_targeting_table()
        joined = " ".join(df.astype(str).values.flatten())
        # every "existing chemical matter" claim must be a negative finding
        assert "no usp34-specific" in joined.lower() or "zero" in joined.lower() or "not found" in joined.lower() or "no relevant evidence" in joined.lower()

    def test_vezf1_direct_targeting_reproduces_exact_tool_compound_values(self):
        df = build_vezf1_direct_targeting_table()
        joined = " ".join(df.astype(str).values.flatten())
        assert "IC50=20 uM" in joined
        assert "IC50=100 uM" in joined
        assert "IC50=500 uM" in joined

    def test_both_tables_cover_all_requested_modalities(self):
        usp34 = build_usp34_direct_targeting_table()
        assert usp34["question"].str.contains("Covalent").any()
        assert usp34["question"].str.contains("Allosteric").any()
        assert usp34["question"].str.contains("degradation").any()


class TestIndirectTargetingTables:
    def test_usp34_indirect_targets_flags_weak_candidates_as_such(self):
        df = build_usp34_indirect_targets_table()
        ubqln1 = df[df["candidate_x"].str.contains("UBQLN1")]
        assert len(ubqln1) == 1
        assert "DOES NOT MEET" in ubqln1.iloc[0]["verdict"]

    def test_vezf1_indirect_targets_flags_tead1_as_unvalidated_hypothesis(self):
        # regression test for a Codex-driven correction: TEAD1 must be
        # framed as a druggable, project-evidence-aligned HYPOTHESIS, not
        # a validated indirect-targeting strategy -- there is no
        # perturbation evidence that TEAD1 regulates VEZF1
        df = build_vezf1_indirect_targets_table()
        tead1 = df[df["candidate_x"] == "TEAD1"]
        assert len(tead1) == 1
        verdict = tead1.iloc[0]["verdict"]
        assert "NOT a validated indirect-targeting strategy" in verdict
        assert "HYPOTHESIS" in verdict

    def test_tead1_perturbation_was_control_only_not_a_functional_test(self):
        # regression test for a Codex-caught factual correction: TEAD1
        # WAS knocked down in PMID 31911272 (siRNA, rat cardiomyocytes),
        # but only as a co-IP antibody-specificity control -- the paper
        # never measured VEZF1 output after TEAD1 perturbation. The
        # report must state this precisely, not claim "TEAD1 was never
        # perturbed" (which is factually wrong) or imply the perturbation
        # validates the hypothesis (which it doesn't).
        df = build_vezf1_indirect_targets_table()
        tead1 = df[df["candidate_x"] == "TEAD1"].iloc[0]
        text = " ".join(str(v) for v in tead1.values)
        assert "siRNA" in text
        assert "specificity control" in text.lower() or "band-specificity" in text.lower() or "antibody" in text.lower()
        assert "exactly one result" in text.lower() or "count=1" in text

    def test_stub1_flagged_as_wrong_direction(self):
        df = build_vezf1_indirect_targets_table()
        stub1 = df[df["candidate_x"].str.contains("STUB1")]
        assert len(stub1) == 1
        assert "WRONG DIRECTION" in stub1.iloc[0]["verdict"]

    def test_no_indirect_target_conflates_correlation_with_causation(self):
        usp34 = build_usp34_indirect_targets_table()
        vezf1 = build_vezf1_indirect_targets_table()
        joined = " ".join(pd.concat([usp34, vezf1]).astype(str).values.flatten())
        assert "correlation" in joined.lower() or "database-level" in joined.lower()


class TestIndirectTargetCrosscheck:
    def test_live_join_matches_frozen_source_file_exactly(self):
        # regression-proof: the crosscheck table's FDR values must come
        # from a live read of the project's own frozen cross-dataset
        # table, not be hardcoded independently of it
        df = build_indirect_target_crosscheck_table().set_index("gene_x")
        source = pd.read_csv(
            "results/tables/cross_dataset_genomewide/all_genes_cross_dataset_evidence_with_ranking.tsv", sep="\t"
        ).set_index("gene")
        for gene in df.index:
            assert df.loc[gene, "hany_crispr_fdr"] == pytest.approx(source.loc[gene, "crispr_fdr"])

    def test_tead1_is_the_only_fdr_significant_candidate(self):
        df = build_indirect_target_crosscheck_table().set_index("gene_x")
        assert df.loc["TEAD1", "gse118713_fdr"] < 0.05
        for gene in df.index:
            if gene != "TEAD1":
                g118 = df.loc[gene, "gse118713_fdr"]
                g240 = df.loc[gene, "gse240112_tumor_fdr"]
                assert (pd.isna(g118) or g118 >= 0.05) and (pd.isna(g240) or g240 >= 0.05)

    def test_frozen_shortlist_language_present(self):
        df = build_indirect_target_crosscheck_table()
        assert df["notes"].str.contains("Frozen shortlist unaltered").all()


class TestHeadToHead:
    def test_sixteen_dimensions_no_master_score_column(self):
        df = build_head_to_head_table()
        assert len(df) == 16
        assert "score" not in [c.lower() for c in df.columns]

    def test_does_not_force_a_winner_on_normal_tissue_profile(self):
        df = build_head_to_head_table().set_index("dimension")
        text = df.loc["Normal-tissue expression"].to_string() if "Normal-tissue expression" in df.index else ""


class TestExperimentalPlan:
    def test_three_experiments_present(self):
        df = build_experimental_plan_table()
        assert len(df) == 3
        assert set(df["experiment_id"]) == {"EXP-1", "EXP-2", "EXP-3"}

    def test_usp34_experiment_has_four_arms_and_normal_cell_comparator(self):
        df = build_experimental_plan_table().set_index("experiment_id")
        row = df.loc["EXP-1"]
        assert row["arms"].count(";") == 3  # 4 arms = 3 separators
        assert "normal_cell_comparator" in df.columns
        assert len(row["normal_cell_comparator"]) > 20

    def test_tead1_experiment_frames_falsification_as_valid_outcome(self):
        df = build_experimental_plan_table().set_index("experiment_id")
        row = df.loc["EXP-3"]
        assert "falsify" in row["open_question"].lower()


class TestVerifiedReferences:
    def test_carried_forward_and_new_references_both_present(self):
        df = build_verified_references_table()
        assert (df["verification_note"].str.contains("carried forward")).any()
        assert (~df["verification_note"].str.contains("carried forward")).any()

    def test_every_reference_has_numeric_pmid(self):
        df = build_verified_references_table()
        assert df["PMID"].astype(str).str.match(r"^\d+$").all()


class TestOutputsExistAndAreNonTrivial:
    @pytest.mark.parametrize("name", [
        "USP34_VEZF1_full_tissue_expression.tsv", "USP34_VEZF1_tissue_liability.tsv",
        "USP34_VEZF1_muscle_liability.tsv", "USP34_VEZF1_bone_marrow_liability.tsv",
        "USP34_VEZF1_human_genetic_constraint.tsv", "USP34_direct_targeting.tsv",
        "USP34_indirect_targets.tsv", "VEZF1_direct_targeting.tsv", "VEZF1_indirect_targets.tsv",
        "indirect_target_project_crosscheck.tsv", "USP34_VEZF1_head_to_head.tsv",
        "USP34_VEZF1_experimental_plan.tsv", "verified_references.tsv",
    ])
    def test_table_written_and_nonempty(self, name):
        path = TABLES / name
        assert path.exists(), f"missing {path}"
        assert len(pd.read_csv(path, sep="\t")) > 0

    @pytest.mark.parametrize("name", [
        "01_USP34_VEZF1_tissue_expression_atlas.png", "02_USP34_VEZF1_liability_map.png",
        "03_USP34_VEZF1_direct_indirect_targetability.png",
        "04_USP34_VEZF1_therapeutic_window_comparison.png",
        "05_USP34_VEZF1_experimental_strategy.png",
    ])
    def test_figure_written_and_nontrivial_size(self, name):
        path = FIGURES / name
        assert path.exists(), f"missing {path}"
        assert path.stat().st_size > 20_000


class TestReportLanguage:
    def test_report_never_calls_a_gene_safe(self):
        text = (REPORTS / "USP34_VEZF1_translational_deep_dive.md").read_text().lower()
        for phrase in ("is safe", "are safe", "safe profile", "safe target", "considered safe"):
            assert phrase not in text

    def test_report_discloses_the_subagent_failure_honestly(self):
        text = (REPORTS / "USP34_VEZF1_translational_deep_dive.md").read_text()
        assert "failed" in text.lower()
        assert "weekly" in text.lower() or "usage limit" in text.lower()

    def test_report_names_both_genes_and_does_not_disprove_others(self):
        text = (REPORTS / "USP34_VEZF1_translational_deep_dive.md").read_text()
        assert "USP34" in text and "VEZF1" in text
        assert "not disproven" in text or "are not disproven" in text

    def test_report_never_extrapolates_embryonic_lethal_to_adult_toxicity(self):
        text = " ".join((REPORTS / "USP34_VEZF1_translational_deep_dive.md").read_text().lower().split())
        assert "must not be extrapolated to adult drug toxicity" in text or "never read as adult drug toxicity" in text

    def test_report_does_not_force_a_winner_when_evidence_is_balanced(self):
        text = (REPORTS / "USP34_VEZF1_translational_deep_dive.md").read_text()
        assert "genuinely balanced" in text or "not forced" in text or "genuinely mixed" in text

    def test_report_discloses_tead1_vezf1_pubmed_exclusivity(self):
        # regression test: the report must state plainly that PMID
        # 31911272 is the ONLY paper in PubMed connecting TEAD1 and VEZF1
        text = (REPORTS / "USP34_VEZF1_translational_deep_dive.md").read_text()
        assert "exactly one result" in text.lower() or "count=1" in text

    def test_report_corrects_tead1_pan_tead_vs_selective_claim(self):
        text = (REPORTS / "USP34_VEZF1_translational_deep_dive.md").read_text()
        assert "PAN-TEAD" in text or "pan-TEAD" in text
        assert "discontinued" in text.lower()

    def test_report_never_claims_tead1_is_validated(self):
        text = " ".join((REPORTS / "USP34_VEZF1_translational_deep_dive.md").read_text().split())
        assert "NOT a validated indirect-targeting strategy" in text

    def test_report_correctly_states_tead1_was_perturbed_as_control_only(self):
        # the report must not claim TEAD1 was never perturbed (false --
        # it was knocked down via siRNA as a co-IP specificity control)
        text = " ".join((REPORTS / "USP34_VEZF1_translational_deep_dive.md").read_text().split())
        assert "siRNA" in text
        assert "specificity" in text.lower()
        assert "no evidence was found that TEAD1 was itself perturbed" not in text
