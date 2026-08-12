from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import yaml

from src.gse245601_candidate_deepdive_data import GENES, load_per_cell_table
from src.gse245601_candidate_deepdive_direction import build_patient_direction_summary
from src.gse245601_candidate_deepdive_pseudobulk import (
    build_group_pseudobulk,
    build_malignancy_condition_patient_summary,
    build_patient_all_epithelial_pseudobulk,
    build_patient_malignant_pseudobulk,
)

REPO_ROOT = Path(__file__).parent.parent
FROZEN_GENES = {"USP34", "VEZF1", "EML5", "CITED2"}


@pytest.fixture(scope="module")
def config():
    with open(REPO_ROOT / "config" / "config.yaml") as f:
        return yaml.safe_load(f)


@pytest.fixture(scope="module")
def per_cell(config):
    path = Path(config["gse245601_candidate_deepdive"]["output"]["per_cell_tsv"])
    if not path.exists():
        pytest.skip("deep-dive per-cell table not generated in this environment")
    return load_per_cell_table(config)


class TestExactFourGenes:
    def test_config_genes_are_exactly_the_frozen_four(self):
        assert set(GENES) == FROZEN_GENES

    def test_per_cell_table_has_columns_for_exactly_these_four(self, per_cell):
        for gene in GENES:
            assert f"{gene}_raw_count" in per_cell.columns
            assert f"{gene}_log_norm" in per_cell.columns


class TestFrozenShortlistUntouched:
    def test_evidence_freeze_files_not_modified(self):
        import subprocess

        result = subprocess.run(["git", "-C", str(REPO_ROOT), "status", "--porcelain", "results/tables/evidence_freeze/", "results/figures/evidence_freeze/", "docs/THERAPEUTIC_SHORTLIST_FREEZE.md"], capture_output=True, text=True)
        for line in result.stdout.splitlines():
            status = line[:2]
            assert status in ("??", ""), f"frozen evidence-freeze path unexpectedly modified: {line}"

    def test_frozen_shortlist_freeze_manifest_still_lists_the_four_genes(self):
        path = REPO_ROOT / "results" / "tables" / "evidence_freeze" / "THERAPEUTIC_SHORTLIST_FREEZE.tsv"
        if not path.exists():
            pytest.skip("evidence freeze not present in this environment")
        df = pd.read_csv(path, sep="\t")
        assert set(df["gene"]) == FROZEN_GENES


class TestPatientAndConditionMapping:
    def test_ten_patients_tumor01_to_tumor10(self, per_cell):
        expected = {f"Tumor_{i:02d}" for i in range(1, 11)}
        assert set(per_cell["patient"]) == expected

    def test_condition_is_control_or_tamoxifen_only(self, per_cell):
        assert set(per_cell["condition"]) == {"Control", "Tamoxifen"}

    def test_sample_id_matches_patient_and_condition(self, per_cell):
        mismatches = per_cell.loc[per_cell["sample_id"] != per_cell["patient"] + "_" + per_cell["condition"]]
        assert len(mismatches) == 0


class TestMalignancyLabelsUnchanged:
    def test_malignancy_status_is_frozen_two_values(self, per_cell):
        assert set(per_cell["malignancy_status"]) == {"malignant", "non-malignant epithelial"}

    def test_malignancy_counts_match_frozen_table_exactly(self, config, per_cell):
        frozen = pd.read_csv(config["gse245601_candidate_deepdive"]["inputs"]["cell_level_summary_tsv"], sep="\t")
        assert per_cell["malignancy_status"].value_counts().to_dict() == frozen["malignancy_status"].value_counts().to_dict()


class TestNoCellPairingAcrossConditions:
    """Individual cells must never be treated as matched pairs across
    Control/Tamoxifen -- only patient-level pseudobulk pairing is valid."""

    def test_cell_ids_are_condition_specific_not_shared(self, per_cell):
        control_ids = set(per_cell.loc[per_cell["condition"] == "Control", "cell_id"])
        tam_ids = set(per_cell.loc[per_cell["condition"] == "Tamoxifen", "cell_id"])
        assert control_ids.isdisjoint(tam_ids), "a cell_id appears in both conditions -- would imply the same cell was measured twice"

    def test_pseudobulk_pairing_is_by_patient_not_by_cell(self, per_cell):
        pb = build_patient_all_epithelial_pseudobulk(per_cell, ["USP34"])
        assert set(pb.columns) >= {"patient", "gene"}
        assert "cell_id" not in pb.columns


class TestPatientLevelPseudobulk:
    def test_pseudobulk_library_size_matches_frozen_track_a_total(self, config, per_cell):
        from src.gse245601_candidate_deepdive_data import load_track_a

        _, ta_meta = load_track_a(config)
        pb = build_group_pseudobulk(per_cell, ["USP34"], ["patient", "condition"])
        pb["sample_id"] = pb["patient"] + "_" + pb["condition"]
        merged = pb.merge(ta_meta[["sample_id", "total_library_size"]], on="sample_id")
        assert (merged["library_size"] == merged["total_library_size"]).all()

    def test_direction_is_from_exact_normalized_values(self):
        pb = pd.DataFrame(
            {
                "gene": ["G1", "G1"], "patient": ["P1", "P2"], "condition": ["Control", "Control"],
                "n_cells": [10, 10], "raw_umi_gene": [5, 5], "library_size": [1000, 1000], "normalized_expression": [1.0, 1.0],
            }
        )
        # minimal smoke check that pivoting/direction logic doesn't crash on a tiny frame
        from src.gse245601_candidate_deepdive_pseudobulk import _add_direction

        out = _add_direction(pb, ["G1"])
        assert "patient" in out.columns

    def test_missing_arm_is_not_comparable_not_equal(self):
        """Regression test for the Phase 25 Codex review (second pass):
        np.select's default branch fires for NaN deltas too, so a missing
        arm (normalized_expression=NaN, e.g. zero cells sampled) was
        previously mislabeled 'equal' instead of 'not_comparable'."""
        from src.gse245601_candidate_deepdive_pseudobulk import _add_direction

        pb = pd.DataFrame(
            {
                "gene": ["G1", "G1"], "patient": ["P1", "P1"], "condition": ["Control", "Tamoxifen"],
                "n_cells": [0, 20], "raw_umi_gene": [0, 5], "library_size": [0, 1000], "normalized_expression": [np.nan, 1.0],
            }
        )
        out = _add_direction(pb, ["G1"])
        assert out.loc[0, "direction_control_to_tam"] == "not_comparable"
        assert pd.isna(out.loc[0, "log_fold_change_descriptive"])

    def test_zero_delta_is_genuinely_equal(self):
        from src.gse245601_candidate_deepdive_pseudobulk import _add_direction

        pb = pd.DataFrame(
            {
                "gene": ["G1", "G1"], "patient": ["P1", "P1"], "condition": ["Control", "Tamoxifen"],
                "n_cells": [20, 20], "raw_umi_gene": [5, 5], "library_size": [1000, 1000], "normalized_expression": [1.0, 1.0],
            }
        )
        out = _add_direction(pb, ["G1"])
        assert out.loc[0, "direction_control_to_tam"] == "equal"


class TestMalignantEligibilityAndLowCellFlagging:
    def test_low_cell_count_warning_flags_below_50(self, per_cell):
        pb = build_patient_malignant_pseudobulk(per_cell, ["USP34"], min_cells=50)
        flagged = pb.loc[pb.get("low_cell_count_warning_Control", pd.Series(dtype=bool)) == True]  # noqa: E712
        for _, row in flagged.iterrows():
            assert row["n_cells_Control"] < 50

    def test_warning_is_both_directions_not_just_flagged_rows(self, per_cell):
        """Regression test for the Phase 25 Codex review finding: the
        warning must equal (n_cells < min_cells) exactly, in both
        directions -- not just "flagged implies low count" but also "low
        count implies flagged", and this must hold even for a patient
        x condition combination with ZERO malignant cells."""
        pb = build_patient_malignant_pseudobulk(per_cell, ["USP34"], min_cells=50)
        for arm in ("Control", "Tamoxifen"):
            n_col, w_col = f"n_cells_{arm}", f"low_cell_count_warning_{arm}"
            assert (pb[w_col] == (pb[n_col] < 50)).all()

    def test_zero_malignant_cell_arm_is_explicit_not_missing(self, per_cell):
        """Tumor_04 has zero malignant Control cells -- this must appear
        as n_cells_Control=0 with low_cell_count_warning_Control=True,
        never as a silently-absent row or a NaN that could be confused
        with a join/pivot bug."""
        pb = build_patient_malignant_pseudobulk(per_cell, ["USP34"], min_cells=50).set_index("patient")
        assert "Tumor_04" in pb.index
        row = pb.loc["Tumor_04"]
        assert row["n_cells_Control"] == 0
        assert bool(row["low_cell_count_warning_Control"]) is True

    def test_reliable_track_b_patients_are_02_03_07(self, config):
        assert config["gse245601_candidate_deepdive"]["track_b_eligible_patients"] == ["Tumor_02", "Tumor_03", "Tumor_07"]


class TestPrevalenceAndIntensity:
    def test_prevalence_uses_raw_count_greater_than_zero(self, per_cell):
        summary = build_malignancy_condition_patient_summary(per_cell.head(500), ["USP34"], include_pooled_all_epithelial=False)
        for _, row in summary.iterrows():
            grp = per_cell.head(500)
            grp = grp.loc[(grp["patient"] == row["patient"]) & (grp["malignancy_status"] == row["malignancy_status"]) & (grp["condition"] == row["condition"])]
            if len(grp) == 0:
                continue
            expected_frac = (grp["USP34_raw_count"] > 0).mean()
            assert np.isclose(row["fraction_expressing"], expected_frac)

    def test_positive_cell_intensity_excludes_zero_count_cells(self):
        df = pd.DataFrame({"patient": ["P1"] * 4, "malignancy_status": ["malignant"] * 4, "condition": ["Control"] * 4, "nCount_RNA": [100, 100, 100, 100], "G_raw_count": [0, 0, 5, 10], "G_log_norm": [0.0, 0.0, 1.0, 2.0]})
        out = build_malignancy_condition_patient_summary(df, ["G"], include_pooled_all_epithelial=False)
        row = out.iloc[0]
        assert np.isclose(row["mean_normalized_positive_cells_only"], 1.5)  # mean of [1.0, 2.0], zero-count cells excluded
        assert np.isclose(row["fraction_expressing"], 0.5)


class TestClusterAnalysis:
    def test_no_reclustering_seurat_clusters_reused(self, per_cell):
        from src.gse245601_candidate_deepdive_clusters import determine_cluster_support

        support = determine_cluster_support(per_cell)
        assert set(support["seurat_clusters"].astype(str)) <= set(per_cell["seurat_clusters"].astype(str))

    def test_support_rule_requires_at_least_3_tumors(self, per_cell):
        from src.gse245601_candidate_deepdive_clusters import MIN_TUMORS_SUPPORTED, determine_cluster_support

        support = determine_cluster_support(per_cell)
        supported = support.loc[support["sufficiently_represented"]]
        assert (supported["n_tumors_supported"] >= MIN_TUMORS_SUPPORTED).all()


class TestResponseStateNotPresentedAsTime:
    def test_no_pseudotime_or_rebound_language_in_docs(self):
        doc_path = REPO_ROOT / "docs" / "GSE245601_CANDIDATE_DEEPDIVE.md"
        if not doc_path.exists():
            pytest.skip("doc not generated")
        text = doc_path.read_text().lower()
        for forbidden in ["pseudotime", "rebound", "goes back up over time", "temporal progression"]:
            assert forbidden not in text

    def test_phase17_response_state_module_genuinely_absent(self):
        """Phase 17 (response-state exploration) was explicitly OMITTED
        (documented in docs/GSE245601_CANDIDATE_DEEPDIVE.md), not silently
        done wrong. This asserts absence of any of its expected artifacts,
        so an accidental future addition without the required candidate
        -gene exclusion would be caught here rather than pass silently."""
        src_dir = REPO_ROOT / "src"
        forbidden_module_names = ["response_state", "addmodulescore", "module_score"]
        for py_file in src_dir.glob("gse245601_candidate_deepdive_*.py"):
            name = py_file.stem.lower()
            for forbidden in forbidden_module_names:
                assert forbidden not in name, f"unexpected response-state module found: {py_file.name}"

        tables_dir = REPO_ROOT / "results" / "tables" / "gse245601_candidate_deepdive"
        if tables_dir.exists():
            forbidden_files = {"response_state_score.tsv", "usp34_response_state.tsv", "all_candidates_response_state.tsv"}
            existing = {p.name for p in tables_dir.glob("*.tsv")}
            assert existing.isdisjoint(forbidden_files), f"response-state output(s) found despite Phase 17 being documented as omitted: {existing & forbidden_files}"

        figures_dir = REPO_ROOT / "results" / "figures" / "gse245601_candidate_deepdive"
        if figures_dir.exists():
            forbidden_figures = {"30_USP34_response_state.png", "31_all_candidates_response_state.png"}
            existing_figs = {p.name for p in figures_dir.glob("*.png")}
            assert existing_figs.isdisjoint(forbidden_figures), f"response-state figure(s) found despite Phase 17 being documented as omitted: {existing_figs & forbidden_figures}"

    def test_if_response_state_module_exists_it_must_exclude_the_four_candidates(self):
        """Structural guard against circularity: if Phase 17 is ever
        implemented, any config-declared response-state signature gene
        list must not include USP34/VEZF1/EML5/CITED2 (a gene cannot be
        correlated against a signature score that already contains it).
        Currently a no-op (no such config key exists) -- it activates
        automatically the moment one is added."""
        with open(REPO_ROOT / "config" / "config.yaml") as f:
            config = yaml.safe_load(f)
        cfg = config.get("gse245601_candidate_deepdive", {})
        signature_genes = cfg.get("response_state_signature_genes")
        if signature_genes is None:
            pytest.skip("no response-state signature configured (Phase 17 correctly not implemented)")
        assert set(signature_genes).isdisjoint(FROZEN_GENES)


class TestPatientDirectionSummary:
    def test_no_arbitrary_flat_band(self):
        pb = pd.DataFrame(
            {
                "gene": ["G1", "G1"], "patient": ["P1", "P2"],
                "log_fold_change_descriptive": [0.5, -0.5],
                "direction_control_to_tam": ["increase", "decrease"],
            }
        )
        out = build_patient_direction_summary(pb)
        row = out.iloc[0]
        assert row["n_patients_increase"] == 1
        assert row["n_patients_decrease"] == 1
        assert row["n_patients_equal"] == 0


class TestPreviousFrozenOutputsUnchanged:
    def test_candidate_adjudication_and_cross_dataset_dirs_not_modified(self):
        import subprocess

        result = subprocess.run(
            ["git", "-C", str(REPO_ROOT), "status", "--porcelain",
             "results/tables/candidate_adjudication/", "results/tables/cross_dataset_genomewide/",
             "results/tables/gse111151/", "results/tables/gse240112_pseudobulk/", "results/tables/gse245601_pseudobulk/"],
            capture_output=True, text=True,
        )
        for line in result.stdout.splitlines():
            status = line[:2]
            assert status in ("??", ""), f"frozen upstream path unexpectedly modified: {line}"
