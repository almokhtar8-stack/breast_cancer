from pathlib import Path

import pandas as pd
import pytest

from src.gse240112_candidate_detection_audit import build_detection_audit

REPO_ROOT = Path(__file__).parent.parent


def _write_tables(tmp_path):
    meta = pd.DataFrame(
        {
            "cell_id": ["c1", "c2", "c3", "c4"],
            "orig.ident": ["PT1", "PT1", "RT1", "RT1"],
            "group": ["PT", "PT", "RT", "RT"],
        }
    )
    raw = pd.DataFrame(
        {
            "cell_id": ["c1", "c2", "c3", "c4"],
            "USP34": [0, 2, 3, 0],
            "PAICS": [1, 1, 0, 0],
        }
    )
    meta_path = tmp_path / "meta.tsv"
    raw_path = tmp_path / "raw.tsv"
    meta.to_csv(meta_path, sep="\t", index=False)
    raw.to_csv(raw_path, sep="\t", index=False)
    return raw_path, meta_path


class TestBuildDetectionAudit:
    def test_present_gene_detection_counts(self, tmp_path):
        raw_path, meta_path = _write_tables(tmp_path)
        out = build_detection_audit(raw_path, meta_path, ["USP34"], "PAICS").set_index("gene")
        assert out.loc["USP34", "present_in_feature_space"]
        assert out.loc["USP34", "n_cells_detected"] == 2
        assert out.loc["USP34", "total_counts"] == 5
        assert out.loc["USP34", "pct_tumor_cells_detected"] == pytest.approx(50.0)
        assert out.loc["USP34", "n_samples_with_detection"] == 2
        assert "PT1=1" in out.loc["USP34", "per_sample_cells_detected"]
        assert "RT1=1" in out.loc["USP34", "per_sample_cells_detected"]

    def test_absent_gene_marked_with_reason(self, tmp_path):
        raw_path, meta_path = _write_tables(tmp_path)
        out = build_detection_audit(raw_path, meta_path, ["USP34", "USP17L29"], "PAICS").set_index("gene")
        assert not out.loc["USP17L29", "present_in_feature_space"]
        assert out.loc["USP17L29", "n_cells_detected"] == 0
        assert out.loc["USP17L29", "reason"] != ""

    def test_paics_flagged_as_benchmark_not_dropped(self, tmp_path):
        raw_path, meta_path = _write_tables(tmp_path)
        out = build_detection_audit(raw_path, meta_path, ["USP34"], "PAICS").set_index("gene")
        assert out.loc["PAICS", "is_paics_benchmark"]
        assert not out.loc["USP34", "is_paics_benchmark"]

    def test_raises_on_cell_id_mismatch(self, tmp_path):
        meta = pd.DataFrame({"cell_id": ["a", "b"], "orig.ident": ["PT1", "PT1"], "group": ["PT", "PT"]})
        raw = pd.DataFrame({"cell_id": ["a", "z"], "USP34": [1, 2]})
        meta_path = tmp_path / "meta.tsv"
        raw_path = tmp_path / "raw.tsv"
        meta.to_csv(meta_path, sep="\t", index=False)
        raw.to_csv(raw_path, sep="\t", index=False)
        with pytest.raises(ValueError, match="cell_id sets differ"):
            build_detection_audit(raw_path, meta_path, ["USP34"], "PAICS")


class TestRealData:
    def test_real_data_if_present(self):
        raw_path = REPO_ROOT / "data" / "processed" / "gse240112" / "tt_cancer_candidate_raw_counts.tsv"
        meta_path = REPO_ROOT / "data" / "processed" / "gse240112" / "tt_cancer_metadata.tsv"
        if not raw_path.exists() or not meta_path.exists():
            pytest.skip("GSE240112 extracted candidate tables not present in this checkout")
        candidates = ["USP34", "CTDNEP1", "EIF4ENIF1", "HMGB1", "KDM1A", "PET117", "TADA2B", "VEZF1", "ICK", "SUPT4H1", "TLK2", "TSR3", "USP17L29"]
        out = build_detection_audit(raw_path, meta_path, candidates, "PAICS").set_index("gene")
        assert len(out) == 14
        assert not out.loc["USP17L29", "present_in_feature_space"]
        assert out.loc["USP34", "present_in_feature_space"]
        assert out.loc["USP34", "n_cells_detected"] > 0
