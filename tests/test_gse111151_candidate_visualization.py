from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from src.gse111151_candidate_visualization import build_candidate_log2cpm_matrix

REPO_ROOT = Path(__file__).parent.parent


class TestBuildCandidateLog2CpmMatrix:
    def test_computes_log2cpm_correctly(self, tmp_path):
        counts = pd.DataFrame({"gene_id": ["ENSG1", "ENSG2"], "gene_name": ["G1", "G2"], "s1": [10, 90], "s2": [50, 50]})
        counts_path = tmp_path / "counts.tsv"
        counts.to_csv(counts_path, sep="\t", index=False)
        out = build_candidate_log2cpm_matrix(counts_path, ["USP34"], {"USP34": "ENSG1"})
        assert list(out.index) == ["USP34"]
        assert out.loc["USP34", "s1"] == pytest.approx(np.log2(10 / 100 * 1e6 + 1))

    def test_missing_gene_excluded_not_invented(self, tmp_path):
        counts = pd.DataFrame({"gene_id": ["ENSG1"], "gene_name": ["G1"], "s1": [10], "s2": [50]})
        counts_path = tmp_path / "counts.tsv"
        counts.to_csv(counts_path, sep="\t", index=False)
        out = build_candidate_log2cpm_matrix(counts_path, ["USP34", "USP17L29"], {"USP34": "ENSG1", "USP17L29": "ENSG_MISSING"})
        assert list(out.index) == ["USP34"]


class TestRealData:
    def test_real_matrix_if_present(self):
        counts_path = REPO_ROOT / "results" / "tables" / "gse111151" / "counts_matrix.tsv.gz"
        if not counts_path.exists():
            pytest.skip("GSE111151 count matrix not present in this checkout")
        import yaml

        cfg = yaml.safe_load(open(REPO_ROOT / "config" / "config.yaml"))["gse111151"]
        out = build_candidate_log2cpm_matrix(counts_path, cfg["candidates"]["thirteen"], cfg["candidate_ensembl_ids"])
        # this function only checks presence in the raw counts matrix, not edgeR's filterByExpr
        # statistical criterion -- USP17L29 has a valid row (all-zero counts), so it IS included here
        assert "USP17L29" in out.index
        assert (out.loc["USP17L29"] == 0).all()
        assert len(out) == 13
