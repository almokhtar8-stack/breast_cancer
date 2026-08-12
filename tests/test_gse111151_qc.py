from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from src.gse111151_qc import (
    build_qc_summary,
    compute_log2cpm,
    compute_pca,
    compute_sample_correlations,
    load_counts,
    load_tmm_norm_factors,
    select_top_variable_genes,
)

REPO_ROOT = Path(__file__).parent.parent

SAMPLES = ["MCF-7", "MCF-7_Tam1", "T-47D", "T-47D_Tam1", "T-47D_Tam2", "ZR-75-1", "ZR-75-1_Tam1", "ZR-75-1_Tam2", "BT-474", "BT-474_Tam1", "BT-474_Tam2"]
CELL_LINES = ["MCF-7", "MCF-7", "T-47D", "T-47D", "T-47D", "ZR-75-1", "ZR-75-1", "ZR-75-1", "BT-474", "BT-474", "BT-474"]
STATUS = ["parental", "resistant", "parental", "resistant", "resistant", "parental", "resistant", "resistant", "parental", "resistant", "resistant"]


def _synthetic_counts():
    genes = [f"ENSG{i:05d}" for i in range(20)]
    rng = np.random.default_rng(0)
    data = {"gene_name": [f"GENE{i}" for i in range(20)]}
    for s in SAMPLES:
        data[s] = rng.integers(10, 1000, size=len(genes))
    counts = pd.DataFrame(data, index=genes)
    counts.index.name = "gene_id"
    metadata = pd.DataFrame(
        {
            "sample_id": SAMPLES,
            "cell_line": CELL_LINES,
            "resistance_status": STATUS,
            "library_size": counts[SAMPLES].sum(axis=0).values,
            "n_detected_genes": (counts[SAMPLES] > 0).sum(axis=0).values,
        }
    )
    return counts.drop(columns="gene_name"), metadata


class TestLoadCounts:
    def test_raises_on_missing_sample(self, tmp_path):
        counts_path = tmp_path / "counts.tsv"
        meta_path = tmp_path / "meta.tsv"
        pd.DataFrame({"gene_id": ["g1"], "gene_name": ["G1"], "MCF-7": [5]}).to_csv(counts_path, sep="\t", index=False)
        pd.DataFrame({"sample_id": ["MCF-7", "T-47D"], "cell_line": ["MCF-7", "T-47D"], "resistance_status": ["parental", "parental"]}).to_csv(meta_path, sep="\t", index=False)
        with pytest.raises(ValueError, match="not present"):
            load_counts(counts_path, meta_path)


class TestComputeLog2Cpm:
    def test_normalizes_by_library_size(self):
        counts = pd.DataFrame({"s1": [10, 90], "s2": [50, 50]}, index=["g1", "g2"])
        out = compute_log2cpm(counts, ["s1", "s2"])
        assert out.loc["g1", "s1"] == pytest.approx(np.log2(10 / 100 * 1e6 + 1))

    def test_uses_effective_lib_size_when_provided(self):
        counts = pd.DataFrame({"s1": [10, 90], "s2": [50, 50]}, index=["g1", "g2"])
        effective = pd.Series({"s1": 200.0, "s2": 100.0})
        out = compute_log2cpm(counts, ["s1", "s2"], effective_lib_sizes=effective)
        # s1 effective lib size (200) differs from its raw sum (100) -- TMM-adjusted result must differ from naive
        naive = compute_log2cpm(counts, ["s1", "s2"])
        assert out.loc["g1", "s1"] != pytest.approx(naive.loc["g1", "s1"])
        assert out.loc["g1", "s1"] == pytest.approx(np.log2(10 / 200 * 1e6 + 1))


class TestLoadTmmNormFactors:
    def test_indexed_by_sample_id(self, tmp_path):
        df = pd.DataFrame({"sample_id": ["s1", "s2"], "library_size": [100, 200], "norm_factor": [1.2, 0.9], "effective_library_size": [120.0, 180.0]})
        path = tmp_path / "tmm.tsv"
        df.to_csv(path, sep="\t", index=False)
        out = load_tmm_norm_factors(path)
        assert out.loc["s1"] == pytest.approx(120.0)
        assert out.loc["s2"] == pytest.approx(180.0)

    def test_real_file_if_present(self):
        path = REPO_ROOT / "results" / "tables" / "gse111151" / "tmm_norm_factors.tsv"
        if not path.exists():
            pytest.skip("GSE111151 TMM norm factors not present in this checkout")
        out = load_tmm_norm_factors(path)
        assert len(out) == 11
        assert (out > 0).all()


class TestQcSummaryAndCorrelation:
    def test_qc_summary_has_all_11_samples(self):
        counts, metadata = _synthetic_counts()
        out = build_qc_summary(metadata)
        assert len(out) == 11
        assert set(out["sample_id"]) == set(SAMPLES)

    def test_correlation_diagonal_is_one(self):
        counts, metadata = _synthetic_counts()
        corr = compute_sample_correlations(counts, metadata["sample_id"].tolist())
        diag = corr.loc[corr["sample_id_1"] == corr["sample_id_2"]]
        assert diag["correlation"].to_numpy() == pytest.approx(1.0)


class TestPca:
    def test_pca_variance_explained_sums_to_one(self):
        counts, metadata = _synthetic_counts()
        top = select_top_variable_genes(counts, metadata["sample_id"].tolist(), 20)
        pca_df = compute_pca(top, metadata, metadata["sample_id"].tolist())
        per_pc_var = pca_df.drop_duplicates("pc")["variance_explained_fraction"]
        assert per_pc_var.sum() == pytest.approx(1.0, abs=1e-6)
        assert len(pca_df["sample_id"].unique()) == 11


class TestRealData:
    def test_real_data_if_present(self):
        counts_path = REPO_ROOT / "results" / "tables" / "gse111151" / "counts_matrix.tsv.gz"
        meta_path = REPO_ROOT / "results" / "tables" / "gse111151" / "sample_metadata.tsv"
        if not counts_path.exists() or not meta_path.exists():
            pytest.skip("GSE111151 count matrix not present in this checkout")
        counts, metadata, gene_names = load_counts(counts_path, meta_path)
        assert len(metadata) == 11
        assert counts.shape[0] == 60619
        assert set(metadata["cell_line"]) == {"MCF-7", "T-47D", "ZR-75-1", "BT-474"}
        assert metadata["resistance_status"].value_counts().to_dict() == {"resistant": 7, "parental": 4}
