from pathlib import Path

import pandas as pd
import pytest

from src.gse240112_usp34_visualization import build_per_sample_summary, load_usp34_cell_data

REPO_ROOT = Path(__file__).parent.parent


def _write_tables(tmp_path):
    meta = pd.DataFrame(
        {
            "cell_id": [f"c{i}" for i in range(6)],
            "orig.ident": ["PT1", "PT1", "PT2", "RT1", "RT1", "RT2"],
            "group": ["PT", "PT", "PT", "RT", "RT", "RT"],
            "umap_1": range(6),
            "umap_2": range(6),
        }
    )
    expr = pd.DataFrame({"cell_id": [f"c{i}" for i in range(6)], "USP34": [0.0, 2.0, 0.0, 3.0, 0.0, 4.0], "CTDNEP1": [1.0] * 6})
    meta_path = tmp_path / "meta.tsv"
    expr_path = tmp_path / "expr.tsv"
    meta.to_csv(meta_path, sep="\t", index=False)
    expr.to_csv(expr_path, sep="\t", index=False)
    return meta_path, expr_path


class TestLoadUsp34CellData:
    def test_joins_correctly(self, tmp_path):
        meta_path, expr_path = _write_tables(tmp_path)
        df = load_usp34_cell_data(meta_path, expr_path)
        assert len(df) == 6
        assert df.loc[df["cell_id"] == "c1", "log_norm_expression"].iloc[0] == 2.0

    def test_raises_if_gene_missing(self, tmp_path):
        meta_path, expr_path = _write_tables(tmp_path)
        with pytest.raises(ValueError, match="not present"):
            load_usp34_cell_data(meta_path, expr_path, gene="NOTAGENE")

    def test_raises_on_cell_id_mismatch(self, tmp_path):
        meta = pd.DataFrame({"cell_id": ["a", "b"], "orig.ident": ["PT1", "PT1"], "group": ["PT", "PT"], "umap_1": [0, 1], "umap_2": [0, 1]})
        expr = pd.DataFrame({"cell_id": ["a", "z"], "USP34": [1.0, 2.0]})
        meta_path = tmp_path / "meta.tsv"
        expr_path = tmp_path / "expr.tsv"
        meta.to_csv(meta_path, sep="\t", index=False)
        expr.to_csv(expr_path, sep="\t", index=False)
        with pytest.raises(ValueError, match="cell_id sets differ"):
            load_usp34_cell_data(meta_path, expr_path)


class TestPerSampleSummary:
    def test_counts_and_detection_by_sample(self, tmp_path):
        meta_path, expr_path = _write_tables(tmp_path)
        df = load_usp34_cell_data(meta_path, expr_path)
        summary = build_per_sample_summary(df).set_index("sample_id")
        assert summary.loc["PT1", "n_cells"] == 2
        assert summary.loc["PT1", "pct_expressing"] == pytest.approx(50.0)
        assert summary.loc["PT2", "n_cells"] == 1
        assert summary.loc["PT2", "pct_expressing"] == pytest.approx(0.0)
        assert summary.loc["RT2", "pct_expressing"] == pytest.approx(100.0)

    def test_missing_sample_reported_as_nan_not_dropped(self, tmp_path):
        meta_path, expr_path = _write_tables(tmp_path)
        df = load_usp34_cell_data(meta_path, expr_path)
        summary = build_per_sample_summary(df)
        assert set(summary["sample_id"]) == {"PT1", "PT2", "PT3", "RT1", "RT2", "RT3"}
        pt3 = summary.set_index("sample_id").loc["PT3"]
        assert pt3["n_cells"] == 0
        assert pd.isna(pt3["pct_expressing"])


class TestRealData:
    def test_real_data_if_present(self):
        meta_path = REPO_ROOT / "data" / "processed" / "gse240112" / "tt_cancer_metadata.tsv"
        expr_path = REPO_ROOT / "data" / "processed" / "gse240112" / "tt_cancer_candidate_lognorm.tsv"
        if not meta_path.exists() or not expr_path.exists():
            pytest.skip("GSE240112 extracted tumor-cell tables not present in this checkout")
        df = load_usp34_cell_data(meta_path, expr_path)
        assert len(df) == 9942
        summary = build_per_sample_summary(df)
        assert summary["n_cells"].sum() == 9942
