from pathlib import Path

import pandas as pd
import pytest

from src.gse245601_usp34_expression import MALIGNANT, NONMALIGNANT, build_group_summary, load_usp34_data

REPO_ROOT = Path(__file__).parent.parent


def _synthetic_df() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "cell_id": [f"c{i}" for i in range(8)],
            "condition": ["Control", "Control", "Control", "Control", "Tamoxifen", "Tamoxifen", "Tamoxifen", "Tamoxifen"],
            "malignancy_status": [MALIGNANT, MALIGNANT, NONMALIGNANT, NONMALIGNANT, MALIGNANT, NONMALIGNANT, NONMALIGNANT, NONMALIGNANT],
            "umap_1": range(8),
            "umap_2": range(8),
            "usp34_log_norm_expression": [0.0, 2.0, 0.0, 1.0, 3.0, 0.0, 0.0, 4.0],
        }
    )


class TestBuildGroupSummary:
    def test_six_groups_with_correct_counts_and_pct_expressing(self):
        df = _synthetic_df()
        out = build_group_summary(df).set_index(["population", "condition"])

        assert out.loc[("all_epithelial", "Control"), "n_cells"] == 4
        assert out.loc[("all_epithelial", "Tamoxifen"), "n_cells"] == 4
        assert out.loc[("malignant", "Control"), "n_cells"] == 2
        assert out.loc[("malignant", "Tamoxifen"), "n_cells"] == 1
        assert out.loc[("nonmalignant", "Control"), "n_cells"] == 2
        assert out.loc[("nonmalignant", "Tamoxifen"), "n_cells"] == 3

        # malignant Control: values [0.0, 2.0] -> 1/2 expressing
        assert out.loc[("malignant", "Control"), "pct_expressing"] == pytest.approx(50.0)
        assert out.loc[("malignant", "Control"), "mean_expression"] == pytest.approx(1.0)
        assert out.loc[("malignant", "Control"), "median_expression"] == pytest.approx(1.0)

        # nonmalignant Tamoxifen: values [0.0, 0.0, 4.0] -> 1/3 expressing
        assert out.loc[("nonmalignant", "Tamoxifen"), "pct_expressing"] == pytest.approx(100.0 / 3.0)

    def test_all_six_requested_groups_present(self):
        df = _synthetic_df()
        out = build_group_summary(df)
        expected = {
            ("all_epithelial", "Control"), ("all_epithelial", "Tamoxifen"),
            ("malignant", "Control"), ("malignant", "Tamoxifen"),
            ("nonmalignant", "Control"), ("nonmalignant", "Tamoxifen"),
        }
        assert set(zip(out["population"], out["condition"])) == expected


class TestLoadUsp34Data:
    def test_raises_on_cell_id_mismatch(self, tmp_path):
        meta = pd.DataFrame({"cell_id": ["a", "b"], "condition": ["Control", "Tamoxifen"], "malignancy_status": [MALIGNANT, NONMALIGNANT], "umap_1": [0, 1], "umap_2": [0, 1]})
        expr = pd.DataFrame({"cell_id": ["a", "c"], "usp34_log_norm_expression": [1.0, 2.0]})
        meta_path = tmp_path / "meta.tsv"
        expr_path = tmp_path / "expr.tsv"
        meta.to_csv(meta_path, sep="\t", index=False)
        expr.to_csv(expr_path, sep="\t", index=False)
        with pytest.raises(ValueError, match="cell_id sets differ"):
            load_usp34_data(meta_path, expr_path)

    def test_raises_on_unexpected_malignancy_status(self, tmp_path):
        meta = pd.DataFrame({"cell_id": ["a"], "condition": ["Control"], "malignancy_status": ["ambiguous"], "umap_1": [0], "umap_2": [0]})
        expr = pd.DataFrame({"cell_id": ["a"], "usp34_log_norm_expression": [1.0]})
        meta_path = tmp_path / "meta.tsv"
        expr_path = tmp_path / "expr.tsv"
        meta.to_csv(meta_path, sep="\t", index=False)
        expr.to_csv(expr_path, sep="\t", index=False)
        with pytest.raises(ValueError, match="unexpected malignancy_status"):
            load_usp34_data(meta_path, expr_path)

    def test_real_data_if_present(self):
        cfg_meta = REPO_ROOT / "results" / "tables" / "gse245601_pseudobulk" / "malignant_vs_nonmalignant_cell_level_summary.tsv"
        cfg_expr = REPO_ROOT / "results" / "tables" / "gse245601_usp34_expression" / "usp34_per_cell_expression.tsv"
        if not cfg_meta.exists() or not cfg_expr.exists():
            pytest.skip("USP34 expression outputs not present in this checkout")
        df = load_usp34_data(cfg_meta, cfg_expr)
        assert len(df) == 29175
        assert (df["usp34_log_norm_expression"] >= 0).all()
