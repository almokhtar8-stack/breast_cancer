from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from src.gse240112_candidate_extraction import (
    benjamini_hochberg,
    build_candidate_table,
    build_paics_row,
    build_sample_level_table,
    load_de_table,
)

REPO_ROOT = Path(__file__).parent.parent


def _synthetic_de():
    return pd.DataFrame(
        {
            "gene": ["USP34", "CTDNEP1", "VEZF1"],
            "log2fc": [0.4, -0.1, 1.1],
            "avg_log_cpm": [7.0, 5.8, 5.5],
            "p_value": [0.12, 0.70, 0.004],
            "fdr": [0.23, 0.80, 0.02],
        }
    )


def _synthetic_counts_meta():
    genes = ["USP34", "CTDNEP1", "VEZF1"]
    counts = pd.DataFrame(
        {
            "PT1": [100, 200, 50],
            "PT2": [110, 210, 55],
            "PT3": [90, 190, 45],
            "RT1": [150, 180, 200],
            "RT2": [160, 175, 210],
            "RT3": [140, 185, 190],
        },
        index=genes,
    )
    counts.index.name = "gene"
    metadata = pd.DataFrame(
        {"sample_id": ["PT1", "PT2", "PT3", "RT1", "RT2", "RT3"], "group": ["PT", "PT", "PT", "RT", "RT", "RT"]}
    )
    return counts, metadata


class TestBenjaminiHochberg:
    def test_matches_known_example(self):
        p = pd.Series([0.01, 0.04, 0.03, 0.005])
        out = benjamini_hochberg(p)
        assert out.iloc[3] == pytest.approx(0.02)  # smallest p, rank1 -> 0.005*4/1=0.02
        assert (out <= 1.0).all() and (out >= 0.0).all()


class TestLoadDeTable:
    def test_raises_on_missing_columns(self, tmp_path):
        p = tmp_path / "de.tsv"
        pd.DataFrame({"gene": ["A"], "log2fc": [1.0]}).to_csv(p, sep="\t", index=False)
        with pytest.raises(ValueError, match="missing required columns"):
            load_de_table(p)


class TestBuildCandidateTable:
    def test_untested_gene_marked_with_reason_not_dropped(self):
        de = _synthetic_de()
        counts, metadata = _synthetic_counts_meta()
        detection_audit = pd.DataFrame({"gene": ["USP17L29"], "present_in_feature_space": [False]})
        out = build_candidate_table(de, counts, metadata, ["USP34", "USP17L29"], detection_audit).set_index("gene")
        assert not out.loc["USP17L29", "tested"]
        assert "absent from tumor-cell feature space" in out.loc["USP17L29", "reason_not_tested"]
        assert pd.isna(out.loc["USP17L29", "candidate_set_bh_fdr"])

    def test_candidate_set_bh_excludes_untested(self):
        de = _synthetic_de()
        counts, metadata = _synthetic_counts_meta()
        detection_audit = pd.DataFrame({"gene": [], "present_in_feature_space": []})
        out = build_candidate_table(de, counts, metadata, ["USP34", "CTDNEP1", "VEZF1"], detection_audit).set_index("gene")
        # 3-gene family: VEZF1 p=0.004 rank1 -> 0.004*3/1=0.012; USP34 p=0.12 rank2 -> 0.12*3/2=0.18; CTDNEP1 p=0.70 rank3 -> 0.70
        assert out.loc["VEZF1", "candidate_set_bh_fdr"] == pytest.approx(0.012)
        assert out.loc["USP34", "candidate_set_bh_fdr"] == pytest.approx(0.18)

    def test_direction_matches_log2fc_sign(self):
        de = _synthetic_de()
        counts, metadata = _synthetic_counts_meta()
        detection_audit = pd.DataFrame({"gene": [], "present_in_feature_space": []})
        out = build_candidate_table(de, counts, metadata, ["USP34", "CTDNEP1"], detection_audit).set_index("gene")
        assert out.loc["USP34", "direction"] == "up_in_RT"
        assert out.loc["CTDNEP1", "direction"] == "down_in_RT"

    def test_separation_flag_true_when_no_overlap(self):
        de = _synthetic_de()
        counts, metadata = _synthetic_counts_meta()
        detection_audit = pd.DataFrame({"gene": [], "present_in_feature_space": []})
        out = build_candidate_table(de, counts, metadata, ["VEZF1"], detection_audit).set_index("gene")
        # VEZF1 raw counts: PT=[50,55,45], RT=[200,210,190] -- fully separated, higher in RT
        assert out.loc["VEZF1", "all_rt_above_pt_range"]
        assert not out.loc["VEZF1", "all_rt_below_pt_range"]


class TestSampleLevelTable:
    def test_one_row_per_gene_per_sample(self):
        counts, metadata = _synthetic_counts_meta()
        out = build_sample_level_table(counts, metadata, ["USP34", "VEZF1"])
        assert len(out) == 2 * 6
        assert set(out["sample_id"]) == {"PT1", "PT2", "PT3", "RT1", "RT2", "RT3"}

    def test_gene_absent_from_counts_skipped(self):
        counts, metadata = _synthetic_counts_meta()
        out = build_sample_level_table(counts, metadata, ["USP34", "USP17L29"])
        assert "USP17L29" not in set(out["gene"])
        assert "USP34" in set(out["gene"])


class TestPaicsRow:
    def test_paics_labeled_as_benchmark(self):
        de = _synthetic_de().rename(columns={"gene": "gene"})
        de.loc[len(de)] = ["PAICS", 0.13, 7.2, 0.69, 0.79]
        counts, metadata = _synthetic_counts_meta()
        counts.loc["PAICS"] = [500, 510, 490, 600, 610, 590]
        out = build_paics_row(de, counts, metadata, "PAICS")
        assert out.iloc[0]["benchmark_label"] == "published_benchmark_not_in_13_candidate_bh_family"
        assert out.iloc[0]["tested"]


class TestRealData:
    def test_real_candidate_table_if_present(self):
        table_path = REPO_ROOT / "results" / "tables" / "gse240112" / "candidate_table.tsv"
        if not table_path.exists():
            pytest.skip("GSE240112 candidate table not present in this checkout")
        out = pd.read_csv(table_path, sep="\t")
        assert len(out) == 13
        assert set(out["gene"]) == {
            "USP34", "CTDNEP1", "EIF4ENIF1", "HMGB1", "KDM1A", "PET117", "TADA2B", "VEZF1", "ICK", "SUPT4H1", "TLK2", "TSR3", "USP17L29",
        }
        usp17l29 = out.set_index("gene").loc["USP17L29"]
        assert not usp17l29["tested"]
        # BH FDR is applied only across tested candidates -- exactly the testable count
        n_tested = int(out["tested"].sum())
        tested_fdrs = out.loc[out["tested"], "candidate_set_bh_fdr"]
        assert tested_fdrs.notna().all()
        assert n_tested <= 13
