from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from src.gse240112_candidate_visualization import build_dotplot_data, build_pseudobulk_log2cpm

REPO_ROOT = Path(__file__).parent.parent


class TestBuildDotplotData:
    def test_mean_and_pct_detected_by_group(self, tmp_path):
        raw = pd.DataFrame({"cell_id": ["a", "b", "c", "d"], "USP34": [0, 2, 4, 0]})
        meta = pd.DataFrame({"cell_id": ["a", "b", "c", "d"], "group": ["PT", "PT", "RT", "RT"]})
        raw_path = tmp_path / "raw.tsv"
        meta_path = tmp_path / "meta.tsv"
        raw.to_csv(raw_path, sep="\t", index=False)
        meta.to_csv(meta_path, sep="\t", index=False)
        out = build_dotplot_data(raw_path, meta_path, ["USP34"]).set_index(["gene", "group"])
        assert out.loc[("USP34", "PT"), "mean_raw_count"] == pytest.approx(1.0)
        assert out.loc[("USP34", "PT"), "pct_detected"] == pytest.approx(50.0)
        assert out.loc[("USP34", "RT"), "mean_raw_count"] == pytest.approx(2.0)
        assert out.loc[("USP34", "RT"), "pct_detected"] == pytest.approx(50.0)

    def test_untested_gene_skipped_not_invented(self, tmp_path):
        raw = pd.DataFrame({"cell_id": ["a", "b"], "USP34": [1, 2]})
        meta = pd.DataFrame({"cell_id": ["a", "b"], "group": ["PT", "RT"]})
        raw_path = tmp_path / "raw.tsv"
        meta_path = tmp_path / "meta.tsv"
        raw.to_csv(raw_path, sep="\t", index=False)
        meta.to_csv(meta_path, sep="\t", index=False)
        out = build_dotplot_data(raw_path, meta_path, ["USP34", "USP17L29"])
        assert "USP17L29" not in set(out["gene"])


class TestBuildPseudobulkLog2Cpm:
    def test_log2cpm_computed_per_sample(self, tmp_path):
        counts = pd.DataFrame(
            {"gene": ["USP34", "CTDNEP1"], "PT1": [10, 90], "PT2": [10, 90], "PT3": [10, 90], "RT1": [10, 90], "RT2": [10, 90], "RT3": [10, 90]}
        )
        counts_path = tmp_path / "counts.tsv"
        counts.to_csv(counts_path, sep="\t", index=False)
        out = build_pseudobulk_log2cpm(counts_path, ["USP34", "CTDNEP1"])
        assert list(out.columns) == ["PT1", "PT2", "PT3", "RT1", "RT2", "RT3"]
        assert out.loc["USP34", "PT1"] == pytest.approx(np.log2(10 / 100 * 1e6 + 1))

    def test_missing_gene_excluded(self, tmp_path):
        counts = pd.DataFrame({"gene": ["USP34"], "PT1": [10], "PT2": [10], "PT3": [10], "RT1": [10], "RT2": [10], "RT3": [10]})
        counts_path = tmp_path / "counts.tsv"
        counts.to_csv(counts_path, sep="\t", index=False)
        out = build_pseudobulk_log2cpm(counts_path, ["USP34", "USP17L29"])
        assert list(out.index) == ["USP34"]


class TestRealData:
    def test_real_dotplot_data_if_present(self):
        raw_path = REPO_ROOT / "data" / "processed" / "gse240112" / "tt_cancer_candidate_raw_counts.tsv"
        meta_path = REPO_ROOT / "data" / "processed" / "gse240112" / "tt_cancer_metadata.tsv"
        if not raw_path.exists() or not meta_path.exists():
            pytest.skip("GSE240112 extracted candidate tables not present in this checkout")
        candidates = ["USP34", "CTDNEP1", "USP17L29"]
        out = build_dotplot_data(raw_path, meta_path, candidates)
        assert set(out["gene"]) == {"USP34", "CTDNEP1"}
        assert set(out["group"]) == {"PT", "RT"}
