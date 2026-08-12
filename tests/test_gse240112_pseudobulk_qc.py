from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from src.gse240112_pseudobulk_qc import (
    build_qc_summary,
    compute_log2cpm,
    compute_pca,
    compute_sample_correlations,
    load_pseudobulk,
    select_top_variable_genes,
)

REPO_ROOT = Path(__file__).parent.parent


def _synthetic_counts():
    genes = [f"g{i}" for i in range(20)]
    rng = np.random.default_rng(0)
    data = {}
    for s, group in [("PT1", "PT"), ("PT2", "PT"), ("PT3", "PT"), ("RT1", "RT"), ("RT2", "RT"), ("RT3", "RT")]:
        data[s] = rng.integers(10, 1000, size=len(genes))
    counts = pd.DataFrame(data, index=genes)
    counts.index.name = "gene"
    metadata = pd.DataFrame(
        {
            "sample_id": ["PT1", "PT2", "PT3", "RT1", "RT2", "RT3"],
            "group": ["PT", "PT", "PT", "RT", "RT", "RT"],
            "n_contributing_cells": [1029, 1975, 1442, 2721, 2597, 178],
            "total_library_size": counts.sum(axis=0).values,
            "n_detected_genes": (counts > 0).sum(axis=0).values,
        }
    )
    return counts, metadata


class TestLoadPseudobulk:
    def test_raises_on_missing_sample(self, tmp_path):
        counts_path = tmp_path / "counts.tsv"
        meta_path = tmp_path / "meta.tsv"
        pd.DataFrame({"gene": ["g1"], "PT1": [5]}).to_csv(counts_path, sep="\t", index=False)
        pd.DataFrame({"sample_id": ["PT1", "RT1"], "group": ["PT", "RT"]}).to_csv(meta_path, sep="\t", index=False)
        with pytest.raises(ValueError, match="not present"):
            load_pseudobulk(counts_path, meta_path)


class TestComputeLog2Cpm:
    def test_normalizes_by_library_size(self):
        counts = pd.DataFrame({"s1": [10, 90], "s2": [50, 50]}, index=["g1", "g2"])
        out = compute_log2cpm(counts, ["s1", "s2"])
        assert out.loc["g1", "s1"] == pytest.approx(np.log2(10 / 100 * 1e6 + 1))
        assert out.loc["g1", "s2"] == pytest.approx(np.log2(50 / 100 * 1e6 + 1))


class TestQcSummaryAndCorrelation:
    def test_qc_summary_sorted_and_complete(self):
        counts, metadata = _synthetic_counts()
        out = build_qc_summary(metadata)
        assert len(out) == 6
        assert list(out["sample_id"]) == sorted(out["sample_id"], key=lambda s: (out.set_index("sample_id").loc[s, "group"], s))

    def test_correlation_diagonal_is_one(self):
        counts, metadata = _synthetic_counts()
        corr = compute_sample_correlations(counts, metadata["sample_id"].tolist())
        diag = corr.loc[corr["sample_id_1"] == corr["sample_id_2"]]
        assert diag["correlation"].to_numpy() == pytest.approx(1.0)


class TestPca:
    def test_top_variable_genes_subset(self):
        counts, metadata = _synthetic_counts()
        top = select_top_variable_genes(counts, metadata["sample_id"].tolist(), 5)
        assert len(top) == 5
        assert set(top.index).issubset(set(counts.index))

    def test_pca_variance_explained_sums_to_one(self):
        counts, metadata = _synthetic_counts()
        top = select_top_variable_genes(counts, metadata["sample_id"].tolist(), 20)
        pca_df = compute_pca(top, metadata, metadata["sample_id"].tolist())
        per_pc_var = pca_df.drop_duplicates("pc")["variance_explained_fraction"]
        assert per_pc_var.sum() == pytest.approx(1.0, abs=1e-6)
        assert len(pca_df["sample_id"].unique()) == 6


class TestRealData:
    def test_real_data_if_present(self):
        cfg_counts = REPO_ROOT / "results" / "tables" / "gse240112_pseudobulk" / "tumor_cell_counts.tsv.gz"
        cfg_meta = REPO_ROOT / "results" / "tables" / "gse240112_pseudobulk" / "tumor_cell_metadata.tsv"
        if not cfg_counts.exists() or not cfg_meta.exists():
            pytest.skip("GSE240112 tumor-cell pseudobulk outputs not present in this checkout")
        counts, metadata = load_pseudobulk(cfg_counts, cfg_meta)
        assert len(metadata) == 6
        assert set(metadata["group"]) == {"PT", "RT"}
        assert metadata["group"].value_counts().to_dict() == {"PT": 3, "RT": 3}
        assert counts.shape[0] == 27161
