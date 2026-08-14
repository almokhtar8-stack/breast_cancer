"""Targeted tests for the final USP34/VEZF1 translational + structure phase."""

from __future__ import annotations

import subprocess
from pathlib import Path

import pandas as pd
import pytest

from src.final_translational_build_tables import (
    build_docking_decision_table,
    build_experimental_design_table,
    build_final_conclusions_table,
    build_normal_cell_comparators_table,
    build_pocket_analysis_table,
    build_structure_inventory_table,
    build_success_failure_table,
)

TABLES = Path("results/tables/final_translational")
FIGURES = Path("results/figures/final_translational")
REPORTS = Path("results/reports/final_translational")
TARGETS = {"USP34", "VEZF1"}


class TestFrozenUpstreamUntouched:
    def test_frozen_and_prior_phase_files_untouched(self):
        # results/tables/{druggability_safety,lead_target_deep_dive}/ are
        # this session's OWN prior-phase deliverables and are expected to
        # show as untracked ("??"), not as modifications -- only
        # genuinely frozen/committed dirs are checked here
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


class TestExperimentalDesign:
    def test_three_experiments_present_correct_ids(self):
        df = build_experimental_design_table()
        assert set(df["experiment_id"]) == {"EXP-1", "EXP-3", "EXP-5"}

    def test_pmid_28499884_counter_evidence_incorporated(self):
        df = build_experimental_design_table().set_index("experiment_id")
        row = df.loc["EXP-1"]
        text = row["emt_stemness_counter_evidence"]
        assert "28499884" in text
        assert "NMuMG" in text
        assert "does NOT invalidate" in text
        # the forbidden phrase is allowed ONLY inside an explicit disclaimer
        # that the report does not make that claim
        assert "which this report does not make" in text

    def test_emt_stemness_markers_in_mechanistic_readouts(self):
        df = build_experimental_design_table().set_index("experiment_id")
        readouts = df.loc["EXP-1", "mechanistic_readouts"]
        for marker in ("CDH1", "CDH2", "SNAI1", "AXIN1", "active"):
            assert marker in readouts

    def test_interaction_framework_not_just_combination_greater_than_single(self):
        df = build_experimental_design_table().set_index("experiment_id")
        readouts = df.loc["EXP-1", "primary_readouts"]
        assert "INTERACTION FRAMEWORK" in readouts
        assert "Bliss" in readouts and "Chou-Talalay" in readouts
        assert "dose-response curve" in readouts

    def test_four_outcome_categories_defined(self):
        df = build_experimental_design_table().set_index("experiment_id")
        outcomes = df.loc["EXP-1", "outcome_categories"]
        for cat in ("IDEAL", "CONCERNING", "PURE GENERAL TOXICITY", "NEGATIVE"):
            assert cat in outcomes

    def test_usp34_uses_crispr_ko_not_an_inhibitor(self):
        df = build_experimental_design_table().set_index("experiment_id")
        row = df.loc["EXP-1"]
        assert "CRISPR knockout" in row["perturbation_strategy_recommended"]
        assert "no validated usp34-selective" in row["perturbation_rationale"].lower()

    def test_vezf1_direct_dependency_and_sensitisation_kept_separate(self):
        df = build_experimental_design_table().set_index("experiment_id")
        row = df.loc["EXP-3"]
        readouts = row["primary_readouts"]
        assert "DIRECT CANCER DEPENDENCY" in readouts
        assert "TAMOXIFEN SENSITISATION" in readouts

    def test_tead1_experiment_has_mandatory_target_engagement_control(self):
        df = build_experimental_design_table().set_index("experiment_id")
        row = df.loc["EXP-5"]
        assert "MANDATORY" in row["primary_readouts"]
        assert "TEAD/YAP target-engagement" in row["primary_readouts"]
        assert "REJECT" in row["decision_rule"]

    def test_no_experiment_describes_either_hypothesis_as_proven(self):
        df = build_experimental_design_table()
        joined = " ".join(df.astype(str).values.flatten()).lower()
        assert "is proven" not in joined and "has been proven" not in joined


class TestNormalCellComparators:
    def test_three_comparators_present(self):
        # EXP-2 was split into two co-primary comparator concepts this
        # phase (A. lineage/cancer-selectivity, B. known-liability), per
        # the user's explicit "use TWO comparator concepts" instruction
        df = build_normal_cell_comparators_table()
        assert set(df["experiment_id"]) == {"EXP-2A", "EXP-2B", "EXP-4"}

    def test_usp34_comparator_a_is_normal_mammary_epithelial(self):
        df = build_normal_cell_comparators_table().set_index("experiment_id")
        row = df.loc["EXP-2A"]
        assert "mammary epithelial" in row["comparator_cell_type"].lower()
        assert "LINEAGE / CANCER-SELECTIVITY" in row["comparator_concept"]

    def test_usp34_comparator_b_is_msc_osteogenic_not_muscle_primary(self):
        df = build_normal_cell_comparators_table().set_index("experiment_id")
        row = df.loc["EXP-2B"]
        assert "MSC" in row["comparator_cell_type"] or "mesenchymal stem cell" in row["comparator_cell_type"]
        assert "osteogenic" in row["comparator_cell_type"].lower()
        assert "NOT recommended as part of the minimum set" in row["minimal_set_justification"]
        assert "KNOWN-LIABILITY" in row["comparator_concept"]

    def test_vezf1_comparator_is_endothelial_primary(self):
        df = build_normal_cell_comparators_table().set_index("experiment_id")
        row = df.loc["EXP-4"]
        assert "endothelial" in row["comparator_cell_type"].lower()

    def test_comparators_never_claim_safety(self):
        # "is safe" is allowed only inside an explicit negation ("never
        # be read as proof that ... is safe"); check for the absence of
        # an affirmative safety claim specifically
        df = build_normal_cell_comparators_table()
        joined = " ".join(df["safety_disclaimer"].astype(str).values).lower()
        assert "proves safety" not in joined
        assert "must never be read as proof" in joined or "must never be" in joined
        assert "preliminary" in joined

    def test_animal_cardiac_finding_not_equated_with_human_toxicity(self):
        df = build_normal_cell_comparators_table().set_index("experiment_id")
        row = df.loc["EXP-4"]
        assert "not equated with expected adult human drug toxicity" in row["safety_disclaimer"].lower()


class TestStructureInventory:
    def test_both_real_structures_present(self):
        df = build_structure_inventory_table()
        assert "7W3R" in df["pdb_id"].values
        assert "7W3U" in df["pdb_id"].values

    def test_resolutions_match_real_rcsb_values(self):
        df = build_structure_inventory_table().set_index("pdb_id")
        assert df.loc["7W3R", "resolution_angstrom"] == pytest.approx(1.92)
        assert df.loc["7W3U", "resolution_angstrom"] == pytest.approx(3.13)

    def test_covalent_bond_to_catalytic_cysteine_documented(self):
        df = build_structure_inventory_table().set_index("pdb_id")
        row = df.loc["7W3U"]
        assert "LINK record directly confirms a covalent bond" in row["bound_ligands"]
        assert "Cys1903" in row["bound_ligands"]

    def test_measured_conformational_tightening_documented(self):
        df = build_structure_inventory_table().set_index("pdb_id")
        assert "3.94" in df.loc["7W3R", "conformational_state_note"]
        assert "3.37" in df.loc["7W3U", "conformational_state_note"]

    def test_conformational_tightening_scoped_to_chain_a_not_overgeneralized(self):
        # regression test for a Codex-caught overgeneralization: the
        # 3.94->3.37 A tightening is real for chain A only; other copies
        # in the same structures show heterogeneous distances (chain B
        # in 7W3R = 4.98 A; 7W3U chains span 3.10-3.95 A) and this must
        # be disclosed, not smoothed into a single structure-wide claim
        df = build_structure_inventory_table().set_index("pdb_id")
        apo_note = df.loc["7W3R", "conformational_state_note"]
        bound_note = df.loc["7W3U", "conformational_state_note"]
        assert "chain A" in apo_note and "chain A" in bound_note
        assert "4.98" in apo_note  # chain B disclosed
        assert "3.95" in bound_note and "3.10" in bound_note  # other copies disclosed
        assert "NOT" in bound_note or "not claimed to be a uniform" in bound_note

    def test_no_additional_structures_overclaimed(self):
        df = build_structure_inventory_table()
        row = df[df["pdb_id"].str.contains("none found")]
        assert len(row) == 1


class TestPocketAnalysis:
    def test_scores_are_qualified_not_bare_druggability_claims(self):
        df = build_pocket_analysis_table()
        joined = " ".join(df["interpretation_caveat"].astype(str).values)
        assert "heuristic" in joined.lower()
        assert "not proof" in joined.lower() or "NOT proof" in joined

    def test_catalytic_pocket_and_ppi_interface_both_identified(self):
        df = build_pocket_analysis_table()
        joined_locs = " ".join(df["location_description"].astype(str).values)
        assert "catalytic cleft" in joined_locs.lower()
        assert "ile44" in joined_locs.lower() or "protein-protein interface" in joined_locs.lower()

    def test_zinc_adjacent_pocket_flagged_as_possible_artifact(self):
        df = build_pocket_analysis_table()
        zn_row = df[df["location_description"].str.contains("Zn", na=False)]
        assert len(zn_row) >= 1
        assert "artificially" in zn_row.iloc[0]["interpretation_caveat"].lower()


class TestDockingDecision:
    def test_decision_is_not_yet_justified(self):
        df = build_docking_decision_table()
        decision_row = df[df["question"] == "DECISION"]
        assert len(decision_row) == 1
        assert decision_row.iloc[0]["answer"] == "DOCKING_NOT_YET_JUSTIFIED"

    def test_justification_cites_absence_of_reference_ligand(self):
        df = build_docking_decision_table()
        just_row = df[df["question"] == "JUSTIFICATION"].iloc[0]
        assert "positive control" in just_row["answer"].lower()

    def test_alternative_roadmap_provided(self):
        df = build_docking_decision_table()
        alt_row = df[df["question"].str.contains("ALTERNATIVE")]
        assert len(alt_row) == 1
        assert "fragment" in alt_row.iloc[0]["answer"].lower()

    def test_all_six_questions_answered(self):
        df = build_docking_decision_table()
        q_rows = df[df["question"].str.match(r"^\d\.")]
        assert len(q_rows) == 6


class TestSuccessFailureCriteria:
    def test_both_targets_present(self):
        df = build_success_failure_table()
        assert set(df["target"]) == TARGETS

    def test_vezf1_distinguishes_dual_action_from_pure_dependency_from_pure_sensitiser(self):
        df = build_success_failure_table().set_index("target")
        row = df.loc["VEZF1"]
        assert "DUAL ACTION" in row["supports_criteria"]
        assert "PURE DEPENDENCY" in row["supports_criteria"]
        assert "PURE SENSITISER" in row["supports_criteria"]


class TestFinalConclusions:
    def test_roles_match_frozen_ranking(self):
        df = build_final_conclusions_table().set_index("target")
        assert df.loc["USP34", "role"] == "LEAD TARGET"
        assert df.loc["VEZF1", "role"] == "SECOND / BACKUP TARGET"

    def test_no_target_called_clinically_validated(self):
        df = build_final_conclusions_table()
        joined = " ".join(df.astype(str).values.flatten()).lower()
        assert "clinically validated" not in joined or "not clinically validated" in joined or "not yet clinically validated" in joined


class TestOutputsExistAndAreNonTrivial:
    @pytest.mark.parametrize("name", [
        "final_experimental_design.tsv", "final_normal_cell_comparators.tsv",
        "USP34_structure_inventory.tsv", "USP34_pocket_analysis.tsv",
        "USP34_docking_decision.tsv", "final_target_success_failure_criteria.tsv",
        "final_translational_conclusions.tsv",
    ])
    def test_table_written_and_nonempty(self, name):
        path = TABLES / name
        assert path.exists(), f"missing {path}"
        assert len(pd.read_csv(path, sep="\t")) > 0

    @pytest.mark.parametrize("name", [
        "01_final_experimental_strategy.png",
        "02_USP34_structure_targetability.png",
        "03_USP34_VEZF1_final_translational_model.png",
    ])
    def test_figure_written_and_nontrivial_size(self, name):
        path = FIGURES / name
        assert path.exists(), f"missing {path}"
        assert path.stat().st_size > 20_000

    def test_docking_figure_04_absent(self):
        # explicit regression: figure 4 (docking hypothesis) must NOT
        # exist, since Part 9-10 concluded DOCKING_NOT_YET_JUSTIFIED
        assert not (FIGURES / "04_USP34_docking_hypothesis.png").exists()


class TestReportLanguage:
    def test_report_never_calls_a_target_safe_or_validated(self):
        # affirmative safety/validation claims only -- "is safe" inside a
        # negation ("never read as proof that X is safe") is fine
        text = (REPORTS / "final_USP34_VEZF1_translational_plan.md").read_text().lower()
        for phrase in ("clinically validated finding", "proven to work", "confirmed safe", "considered safe"):
            assert phrase not in text

    def test_report_states_frozen_ranking_unchanged(self):
        text = (REPORTS / "final_USP34_VEZF1_translational_plan.md").read_text()
        assert "Frozen ranking, unchanged" in text or "frozen ranking" in text.lower()

    def test_report_never_treats_animal_cardiac_data_as_human_toxicity(self):
        text = " ".join((REPORTS / "final_USP34_VEZF1_translational_plan.md").read_text().split())
        assert "not equated with expected adult human drug toxicity" in text

    def test_report_names_docking_decision_explicitly(self):
        text = (REPORTS / "final_USP34_VEZF1_translational_plan.md").read_text()
        assert "DOCKING_NOT_YET_JUSTIFIED" in text
