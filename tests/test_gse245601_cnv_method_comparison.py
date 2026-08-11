from pathlib import Path

import pandas as pd
import pytest

from src.gse245601_cnv_method_comparison import (
    LABEL_MALIGNANT,
    LABEL_NONMALIGNANT,
    LABEL_NOT_DEFINED,
    build_by_sample_table,
    build_by_tumor_table,
    load_copykat_labels,
    load_infercnv_labels,
    load_infercnv_sample_summary,
    match_cells,
    verify_overall_agreement,
)

REPO_ROOT = Path(__file__).parent.parent
CONFIG_PATH = REPO_ROOT / "config" / "config.yaml"


def _infercnv_row(cell_id, sample_id, patient, condition, label):
    return {
        "cell_id": cell_id,
        "sample_id": sample_id,
        "patient": patient,
        "condition": condition,
        "seurat_cluster": "0",
        "threshold_group": "0",
        "cnv_score": 0.01,
        "cnv_correlation_to_seed": 0.5,
        "primary_malignancy_label": label,
    }


def _copykat_row(cell_id, sample_id, label):
    return {"cell_id": cell_id, "sensitivity_malignancy_label": label, "sample_id": sample_id}


def _toy_tables():
    """4 cells in Tumor_01_Control, 4 in Tumor_01_Tamoxifen -- a small,
    hand-worked 2x2 breakdown including one not_defined cell, exercised
    against the real merge/aggregation code rather than re-deriving the
    arithmetic in the test."""
    infercnv_rows = [
        _infercnv_row("c1", "Tumor_01_Control", "Tumor_01", "Control", LABEL_MALIGNANT),
        _infercnv_row("c2", "Tumor_01_Control", "Tumor_01", "Control", LABEL_MALIGNANT),
        _infercnv_row("c3", "Tumor_01_Control", "Tumor_01", "Control", LABEL_NONMALIGNANT),
        _infercnv_row("c4", "Tumor_01_Control", "Tumor_01", "Control", LABEL_NONMALIGNANT),
        _infercnv_row("c5", "Tumor_01_Tamoxifen", "Tumor_01", "Tamoxifen", LABEL_MALIGNANT),
        _infercnv_row("c6", "Tumor_01_Tamoxifen", "Tumor_01", "Tamoxifen", LABEL_NONMALIGNANT),
        _infercnv_row("c7", "Tumor_01_Tamoxifen", "Tumor_01", "Tamoxifen", LABEL_NONMALIGNANT),
        _infercnv_row("c8", "Tumor_01_Tamoxifen", "Tumor_01", "Tamoxifen", LABEL_MALIGNANT),
    ]
    copykat_rows = [
        _copykat_row("c1", "Tumor_01_Control", LABEL_MALIGNANT),  # both malignant
        _copykat_row("c2", "Tumor_01_Control", LABEL_NOT_DEFINED),  # excluded from n_compared
        _copykat_row("c3", "Tumor_01_Control", LABEL_NONMALIGNANT),  # both nonmalignant
        _copykat_row("c4", "Tumor_01_Control", LABEL_MALIGNANT),  # infercnv non, copykat aneuploid
        _copykat_row("c5", "Tumor_01_Tamoxifen", LABEL_NONMALIGNANT),  # infercnv malignant, copykat diploid
        _copykat_row("c6", "Tumor_01_Tamoxifen", LABEL_NONMALIGNANT),  # both nonmalignant
        _copykat_row("c7", "Tumor_01_Tamoxifen", LABEL_NONMALIGNANT),  # both nonmalignant
        _copykat_row("c8", "Tumor_01_Tamoxifen", LABEL_MALIGNANT),  # both malignant
    ]
    infercnv_df = pd.DataFrame(infercnv_rows)
    copykat_df = pd.DataFrame(copykat_rows)
    summary_df = pd.DataFrame(
        [
            {"sample_id": "Tumor_01_Control", "n_epithelial_cells": 4, "status": "ok"},
            {"sample_id": "Tumor_01_Tamoxifen", "n_epithelial_cells": 4, "status": "ok"},
        ]
    )
    return infercnv_df, copykat_df, summary_df


class TestLoaders:
    def test_load_infercnv_labels_rejects_unknown_label(self, tmp_path):
        df = pd.DataFrame([_infercnv_row("c1", "Tumor_01_Control", "Tumor_01", "Control", "ambiguous")])
        path = tmp_path / "infercnv.tsv"
        df.to_csv(path, sep="\t", index=False)
        with pytest.raises(ValueError, match="unexpected primary_malignancy_label"):
            load_infercnv_labels(path)

    def test_load_infercnv_labels_rejects_duplicate_cell_id(self, tmp_path):
        df = pd.DataFrame(
            [
                _infercnv_row("c1", "Tumor_01_Control", "Tumor_01", "Control", LABEL_MALIGNANT),
                _infercnv_row("c1", "Tumor_01_Control", "Tumor_01", "Control", LABEL_NONMALIGNANT),
            ]
        )
        path = tmp_path / "infercnv.tsv"
        df.to_csv(path, sep="\t", index=False)
        with pytest.raises(ValueError, match="duplicate cell_id"):
            load_infercnv_labels(path)

    def test_load_copykat_labels_rejects_unknown_label(self, tmp_path):
        df = pd.DataFrame([_copykat_row("c1", "Tumor_01_Control", "aneuploid")])
        path = tmp_path / "copykat.tsv"
        df.to_csv(path, sep="\t", index=False)
        with pytest.raises(ValueError, match="unexpected sensitivity_malignancy_label"):
            load_copykat_labels(path)

    def test_load_real_frozen_tables_if_present(self):
        """If the real frozen tables exist on disk, they must satisfy the
        same loader validation used on synthetic data (schema/label-set
        checks only -- no comparison of values)."""
        cfg_inputs = {
            "infercnv": REPO_ROOT / "results" / "tables" / "gse245601_malignant_cell_labels.tsv",
            "copykat": REPO_ROOT / "results" / "tables" / "gse245601_copykat_sensitivity_labels.tsv",
            "summary": REPO_ROOT / "results" / "tables" / "gse245601_malignant_summary_per_sample.tsv",
        }
        if not all(p.exists() for p in cfg_inputs.values()):
            pytest.skip("frozen GSE245601 result tables not present in this checkout")
        infercnv_df = load_infercnv_labels(cfg_inputs["infercnv"])
        copykat_df = load_copykat_labels(cfg_inputs["copykat"])
        summary_df = load_infercnv_sample_summary(cfg_inputs["summary"])
        assert len(infercnv_df) > 0
        assert len(copykat_df) > 0
        assert len(summary_df) > 0


class TestMatchCells:
    def test_matches_by_sample_and_cell_id(self):
        infercnv_df, copykat_df, _ = _toy_tables()
        merged = match_cells(infercnv_df, copykat_df)
        assert len(merged) == 8
        assert set(merged["cell_id"]) == set(infercnv_df["cell_id"])

    def test_raises_on_infercnv_cell_missing_from_copykat(self):
        infercnv_df, copykat_df, _ = _toy_tables()
        copykat_df = copykat_df.loc[copykat_df["cell_id"] != "c1"]
        with pytest.raises(ValueError, match="cell_id sets do not match"):
            match_cells(infercnv_df, copykat_df)

    def test_raises_on_sample_id_mismatch_for_same_cell_id(self):
        infercnv_df, copykat_df, _ = _toy_tables()
        copykat_df = copykat_df.copy()
        copykat_df.loc[copykat_df["cell_id"] == "c1", "sample_id"] = "Tumor_02_Control"
        with pytest.raises(ValueError, match="different sample_id"):
            match_cells(infercnv_df, copykat_df)


class TestBuildBySampleTable:
    def test_2x2_breakdown_and_denominators(self):
        infercnv_df, copykat_df, summary_df = _toy_tables()
        merged = match_cells(infercnv_df, copykat_df)
        out = build_by_sample_table(merged, summary_df)
        assert len(out) == 2

        control = out.loc[out["sample_id"] == "Tumor_01_Control"].iloc[0]
        assert control["total_epithelial_cells"] == 4
        assert control["n_compared"] == 3  # c2 excluded (not_defined)
        assert control["copykat_not_defined_count"] == 1
        assert control["both_malignant_aneuploid"] == 1  # c1
        assert control["both_nonmalignant_diploid"] == 1  # c3
        assert control["infercnv_malignant_copykat_diploid"] == 0
        assert control["infercnv_nonmalignant_copykat_aneuploid"] == 1  # c4
        assert control["infercnv_malignant_count"] == 1
        assert control["copykat_aneuploid_count"] == 2
        assert control["agreement_count"] == 2
        assert control["agreement_pct"] == pytest.approx(200.0 / 3.0)

        tam = out.loc[out["sample_id"] == "Tumor_01_Tamoxifen"].iloc[0]
        assert tam["n_compared"] == 4
        assert tam["copykat_not_defined_count"] == 0
        assert tam["both_malignant_aneuploid"] == 1  # c8
        assert tam["both_nonmalignant_diploid"] == 2  # c6, c7
        assert tam["infercnv_malignant_copykat_diploid"] == 1  # c5
        assert tam["infercnv_nonmalignant_copykat_aneuploid"] == 0
        assert tam["agreement_count"] == 3
        assert tam["agreement_pct"] == pytest.approx(75.0)

    def test_2x2_counts_always_sum_to_n_compared(self):
        infercnv_df, copykat_df, summary_df = _toy_tables()
        merged = match_cells(infercnv_df, copykat_df)
        out = build_by_sample_table(merged, summary_df)
        summed = (
            out["both_malignant_aneuploid"]
            + out["both_nonmalignant_diploid"]
            + out["infercnv_malignant_copykat_diploid"]
            + out["infercnv_nonmalignant_copykat_aneuploid"]
        )
        assert (summed == out["n_compared"]).all()

    def test_raises_if_total_epithelial_does_not_match_matched_cells(self):
        infercnv_df, copykat_df, summary_df = _toy_tables()
        summary_df = summary_df.copy()
        summary_df.loc[summary_df["sample_id"] == "Tumor_01_Control", "n_epithelial_cells"] = 999
        merged = match_cells(infercnv_df, copykat_df)
        with pytest.raises(ValueError, match="total_epithelial_cells"):
            build_by_sample_table(merged, summary_df)


class TestBuildByTumorTable:
    def test_pairs_control_and_tamoxifen_with_correct_diff_sign(self):
        infercnv_df, copykat_df, summary_df = _toy_tables()
        merged = match_cells(infercnv_df, copykat_df)
        by_sample = build_by_sample_table(merged, summary_df)
        by_tumor = build_by_tumor_table(by_sample)

        assert len(by_tumor) == 1
        row = by_tumor.iloc[0]
        assert row["tumor"] == "Tumor_01"
        control_pct = by_sample.loc[by_sample["sample_id"] == "Tumor_01_Control", "infercnv_malignant_pct"].iloc[0]
        tam_pct = by_sample.loc[by_sample["sample_id"] == "Tumor_01_Tamoxifen", "infercnv_malignant_pct"].iloc[0]
        assert row["infercnv_malignant_pct_diff_tam_minus_control"] == pytest.approx(tam_pct - control_pct)

    def test_raises_if_a_tumor_is_missing_a_condition(self):
        infercnv_df, copykat_df, summary_df = _toy_tables()
        merged = match_cells(infercnv_df, copykat_df)
        by_sample = build_by_sample_table(merged, summary_df)
        by_sample = by_sample.loc[by_sample["condition"] != "Tamoxifen"]
        with pytest.raises(ValueError, match="not every tumor has both"):
            build_by_tumor_table(by_sample)


class TestVerifyOverallAgreement:
    def test_passes_when_within_tolerance(self):
        infercnv_df, copykat_df, summary_df = _toy_tables()
        merged = match_cells(infercnv_df, copykat_df)
        by_sample = build_by_sample_table(merged, summary_df)
        overall = 100.0 * by_sample["agreement_count"].sum() / by_sample["n_compared"].sum()
        result = verify_overall_agreement(by_sample, expected_pct=overall)
        assert result == pytest.approx(overall)

    def test_raises_when_outside_tolerance(self):
        infercnv_df, copykat_df, summary_df = _toy_tables()
        merged = match_cells(infercnv_df, copykat_df)
        by_sample = build_by_sample_table(merged, summary_df)
        with pytest.raises(ValueError, match="does not match"):
            verify_overall_agreement(by_sample, expected_pct=1.0, tolerance_pct=0.01)


class TestRealDataReproducesReportedConcordance:
    def test_pooled_agreement_reproduces_existing_concordance_table(self):
        """Independent, from-scratch recomputation: reads the same two raw
        per-cell label tables the original concordance script read, using
        this module's own matching/aggregation code, and checks the pooled
        agreement against the existing, separately-computed
        gse245601_malignancy_concordance.tsv (sum(n_concordant)/sum(n_compared))."""
        infercnv_path = REPO_ROOT / "results" / "tables" / "gse245601_malignant_cell_labels.tsv"
        copykat_path = REPO_ROOT / "results" / "tables" / "gse245601_copykat_sensitivity_labels.tsv"
        summary_path = REPO_ROOT / "results" / "tables" / "gse245601_malignant_summary_per_sample.tsv"
        concordance_path = REPO_ROOT / "results" / "tables" / "gse245601_malignancy_concordance.tsv"
        if not all(p.exists() for p in (infercnv_path, copykat_path, summary_path, concordance_path)):
            pytest.skip("frozen GSE245601 result tables not present in this checkout")

        infercnv_df = load_infercnv_labels(infercnv_path)
        copykat_df = load_copykat_labels(copykat_path)
        summary_df = load_infercnv_sample_summary(summary_path)
        merged = match_cells(infercnv_df, copykat_df)
        by_sample = build_by_sample_table(merged, summary_df)

        concordance = pd.read_csv(concordance_path, sep="\t")
        concordance_ok = concordance.loc[concordance["status"] == "ok"]
        expected_pct = 100.0 * concordance_ok["n_concordant"].sum() / concordance_ok["n_compared"].sum()

        result = verify_overall_agreement(by_sample, expected_pct=expected_pct, tolerance_pct=0.001)
        assert result == pytest.approx(expected_pct, abs=0.001)
        assert int(by_sample["n_compared"].sum()) == int(concordance_ok["n_compared"].sum())
        assert int(by_sample["agreement_count"].sum()) == int(concordance_ok["n_concordant"].sum())
