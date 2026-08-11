from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import yaml

from src.gse245601_cnv_heatmap_comparison import (
    AGREE_NORMAL_LIKE,
    AGREE_TUMOR_LIKE,
    COPYKAT_NOT_DEFINED_STATUS,
    DISAGREE,
    CnvHeatmapComparisonConfig,
    build_annotation_table,
    build_cell_order,
    build_contact_sheet,
    build_sample_figure,
    build_summary_table,
    compute_agreement_status,
    load_compact_matrices,
    run_cnv_heatmap_comparison,
    verify_matched_columns,
)
from src.gse245601_cnv_method_comparison import LABEL_MALIGNANT, LABEL_NONMALIGNANT, LABEL_NOT_DEFINED

REPO_ROOT = Path(__file__).parent.parent
CONFIG_PATH = REPO_ROOT / "config" / "config.yaml"

EXPECTED_SELECTED_SAMPLES = (
    "Tumor_01_Control",
    "Tumor_01_Tamoxifen",
    "Tumor_02_Control",
    "Tumor_02_Tamoxifen",
    "Tumor_03_Control",
    "Tumor_03_Tamoxifen",
    "Tumor_04_Control",
    "Tumor_04_Tamoxifen",
    "Tumor_08_Control",
    "Tumor_08_Tamoxifen",
    "Tumor_10_Control",
    "Tumor_10_Tamoxifen",
)


def _real_config() -> dict:
    with open(CONFIG_PATH) as f:
        return yaml.safe_load(f)


def _real_cfg() -> CnvHeatmapComparisonConfig:
    return CnvHeatmapComparisonConfig.from_config(_real_config())


class TestExpectedSelectedSampleSet:
    def test_selected_samples_match_exactly_in_order(self):
        cfg = _real_cfg()
        assert cfg.selected_samples == EXPECTED_SELECTED_SAMPLES

    def test_control_and_tamoxifen_adjacent_for_each_tumor(self):
        cfg = _real_cfg()
        for i in range(0, len(cfg.selected_samples), 2):
            control, tamoxifen = cfg.selected_samples[i], cfg.selected_samples[i + 1]
            control_patient, control_condition = control.rsplit("_", 1)
            tam_patient, tam_condition = tamoxifen.rsplit("_", 1)
            assert control_condition == "Control"
            assert tam_condition == "Tamoxifen"
            assert control_patient == tam_patient


class TestComputeAgreementStatus:
    def test_both_malignant_is_agree_tumor_like(self):
        assert compute_agreement_status(LABEL_MALIGNANT, LABEL_MALIGNANT) == AGREE_TUMOR_LIKE

    def test_both_nonmalignant_is_agree_normal_like(self):
        assert compute_agreement_status(LABEL_NONMALIGNANT, LABEL_NONMALIGNANT) == AGREE_NORMAL_LIKE

    def test_infercnv_malignant_copykat_diploid_is_disagree(self):
        assert compute_agreement_status(LABEL_MALIGNANT, LABEL_NONMALIGNANT) == DISAGREE

    def test_infercnv_nonmalignant_copykat_aneuploid_is_disagree(self):
        assert compute_agreement_status(LABEL_NONMALIGNANT, LABEL_MALIGNANT) == DISAGREE

    def test_copykat_not_defined_overrides_infercnv_malignant(self):
        assert compute_agreement_status(LABEL_MALIGNANT, LABEL_NOT_DEFINED) == COPYKAT_NOT_DEFINED_STATUS

    def test_copykat_not_defined_overrides_infercnv_nonmalignant(self):
        assert compute_agreement_status(LABEL_NONMALIGNANT, LABEL_NOT_DEFINED) == COPYKAT_NOT_DEFINED_STATUS

    def test_raises_on_unexpected_infercnv_label(self):
        with pytest.raises(ValueError, match="unexpected infercnv_label"):
            compute_agreement_status("ambiguous", LABEL_MALIGNANT)

    def test_raises_on_unexpected_copykat_label(self):
        with pytest.raises(ValueError, match="unexpected copykat_label"):
            compute_agreement_status(LABEL_MALIGNANT, "ambiguous")


class TestBuildAnnotationTable:
    def _tables(self):
        infercnv_df = pd.DataFrame(
            {
                "cell_id": ["c1", "c2", "c3"],
                "sample_id": ["S1", "S1", "S1"],
                "primary_malignancy_label": [LABEL_MALIGNANT, LABEL_NONMALIGNANT, LABEL_MALIGNANT],
                "cnv_score": [0.05, 0.01, 0.03],
            }
        )
        copykat_df = pd.DataFrame(
            {
                "cell_id": ["c1", "c2", "c3"],
                "sample_id": ["S1", "S1", "S1"],
                "sensitivity_malignancy_label": [LABEL_MALIGNANT, LABEL_NONMALIGNANT, LABEL_NOT_DEFINED],
            }
        )
        return infercnv_df, copykat_df

    def test_agreement_status_computed_per_cell(self):
        infercnv_df, copykat_df = self._tables()
        out = build_annotation_table("S1", infercnv_df, copykat_df)
        out = out.set_index("cell_id")
        assert out.loc["c1", "agreement_status"] == AGREE_TUMOR_LIKE
        assert out.loc["c2", "agreement_status"] == AGREE_NORMAL_LIKE
        assert out.loc["c3", "agreement_status"] == COPYKAT_NOT_DEFINED_STATUS

    def test_raises_on_mismatched_cell_sets(self):
        infercnv_df, copykat_df = self._tables()
        copykat_df = copykat_df.loc[copykat_df["cell_id"] != "c1"]
        with pytest.raises(ValueError, match="cell_id sets do not match"):
            build_annotation_table("S1", infercnv_df, copykat_df)


class TestBuildCellOrder:
    def test_malignant_first_then_by_cnv_score_descending(self):
        annotation_df = pd.DataFrame(
            {
                "cell_id": ["a", "b", "c", "d"],
                "infercnv_label": [LABEL_NONMALIGNANT, LABEL_MALIGNANT, LABEL_MALIGNANT, LABEL_NONMALIGNANT],
                "cnv_score": [0.5, 0.1, 0.9, 0.2],
            }
        )
        order = build_cell_order(annotation_df)
        # malignant block (c=0.9, b=0.1) desc, then non-malignant block (a=0.5, d=0.2) desc
        assert order == ["c", "b", "a", "d"]

    def test_raises_on_unexpected_label(self):
        annotation_df = pd.DataFrame({"cell_id": ["a"], "infercnv_label": ["ambiguous"], "cnv_score": [0.1]})
        with pytest.raises(ValueError, match="unexpected infercnv_label"):
            build_cell_order(annotation_df)


class TestVerifyMatchedColumns:
    def test_passes_when_columns_match(self):
        ic = pd.DataFrame({"c1": [1.0], "c2": [1.0]}, index=["gene1"])
        ck = pd.DataFrame({"c2": [0.0], "c1": [0.0]}, index=["gene1"])
        verify_matched_columns(ic, ck, ["c1", "c2"])  # should not raise

    def test_raises_when_infercnv_missing_a_cell(self):
        ic = pd.DataFrame({"c1": [1.0]}, index=["gene1"])
        ck = pd.DataFrame({"c1": [0.0], "c2": [0.0]}, index=["gene1"])
        with pytest.raises(ValueError, match="InferCNV matrix columns"):
            verify_matched_columns(ic, ck, ["c1", "c2"])

    def test_raises_when_copykat_has_an_extra_cell(self):
        ic = pd.DataFrame({"c1": [1.0], "c2": [1.0]}, index=["gene1"])
        ck = pd.DataFrame({"c1": [0.0], "c2": [0.0], "c3": [0.0]}, index=["gene1"])
        with pytest.raises(ValueError, match="CopyKAT matrix columns"):
            verify_matched_columns(ic, ck, ["c1", "c2"])

    def test_raises_on_duplicate_column_in_infercnv_matrix(self):
        ic = pd.DataFrame([[1.0, 1.0]], columns=["c1", "c1"], index=["gene1"])
        ck = pd.DataFrame({"c1": [0.0]}, index=["gene1"])
        with pytest.raises(ValueError, match="InferCNV matrix has duplicate"):
            verify_matched_columns(ic, ck, ["c1"])

    def test_raises_on_duplicate_expected_cell_id(self):
        ic = pd.DataFrame({"c1": [1.0]}, index=["gene1"])
        ck = pd.DataFrame({"c1": [0.0]}, index=["gene1"])
        with pytest.raises(ValueError, match="expected_cell_ids contains duplicate"):
            verify_matched_columns(ic, ck, ["c1", "c1"])


class TestLoadCompactMatricesValidation:
    def _write_sample(self, tmp_path, ic_genes, ck_genes, gene_order_genes):
        sample_dir = tmp_path / "S1"
        sample_dir.mkdir()
        pd.DataFrame({"gene": ic_genes, "c1": [1.0] * len(ic_genes)}).to_csv(
            sample_dir / "infercnv_matrix.tsv.gz", sep="\t", index=False
        )
        pd.DataFrame({"gene": ck_genes, "c1": [0.0] * len(ck_genes)}).to_csv(
            sample_dir / "copykat_matrix.tsv.gz", sep="\t", index=False
        )
        pd.DataFrame({"gene": gene_order_genes, "chr": [1] * len(gene_order_genes), "start": range(len(gene_order_genes))}).to_csv(
            sample_dir / "gene_order.tsv", sep="\t", index=False
        )
        return tmp_path

    def test_raises_on_duplicate_gene_row_in_infercnv_matrix(self, tmp_path):
        working_dir = self._write_sample(tmp_path, ["g1", "g1"], ["g1", "g2"], ["g1", "g2"])
        cfg = CnvHeatmapComparisonConfig(
            selected_samples=("S1",),
            infercnv_labels_tsv=Path("x"),
            copykat_labels_tsv=Path("x"),
            method_comparison_by_sample_tsv=Path("x"),
            working_dir=working_dir,
            extraction_r_script=Path("x"),
            micromamba_binary="x",
            micromamba_env="x",
            figures_dir=Path("x"),
            tables_dir=Path("x"),
            summary_tsv=Path("x"),
            gene_exclusion_tsv=Path("x"),
            contact_sheet_png=Path("x"),
        )
        with pytest.raises(ValueError, match="duplicate gene rows"):
            load_compact_matrices("S1", cfg)

    def test_raises_on_mismatched_gene_sets(self, tmp_path):
        working_dir = self._write_sample(tmp_path, ["g1", "g2"], ["g1", "g3"], ["g1", "g2"])
        cfg = CnvHeatmapComparisonConfig(
            selected_samples=("S1",),
            infercnv_labels_tsv=Path("x"),
            copykat_labels_tsv=Path("x"),
            method_comparison_by_sample_tsv=Path("x"),
            working_dir=working_dir,
            extraction_r_script=Path("x"),
            micromamba_binary="x",
            micromamba_env="x",
            figures_dir=Path("x"),
            tables_dir=Path("x"),
            summary_tsv=Path("x"),
            gene_exclusion_tsv=Path("x"),
            contact_sheet_png=Path("x"),
        )
        with pytest.raises(ValueError, match="gene sets do not match"):
            load_compact_matrices("S1", cfg)


class TestBuildSampleFigureOutputsCreated:
    def test_writes_a_nonempty_png(self, tmp_path):
        genes = [f"g{i}" for i in range(20)]
        cells = [f"c{i}" for i in range(6)]
        rng = np.random.default_rng(0)
        ic_mat = pd.DataFrame(1.0 + 0.05 * rng.standard_normal((len(genes), len(cells))), index=genes, columns=cells)
        ck_mat = pd.DataFrame(0.05 * rng.standard_normal((len(genes), len(cells))), index=genes, columns=cells)
        gene_order = pd.DataFrame({"gene": genes, "chr": [1] * 10 + [2] * 10, "start": list(range(20))})
        annotation_df = pd.DataFrame(
            {
                "cell_id": cells,
                "infercnv_label": [LABEL_MALIGNANT, LABEL_NONMALIGNANT] * 3,
                "copykat_label": [LABEL_MALIGNANT, LABEL_NONMALIGNANT] * 3,
                "cnv_score": [0.1, 0.2, 0.3, 0.4, 0.5, 0.6],
                "agreement_status": [AGREE_TUMOR_LIKE, AGREE_NORMAL_LIKE] * 3,
            }
        )
        out_path = tmp_path / "fig.png"
        build_sample_figure("SyntheticSample", ic_mat, ck_mat, gene_order, annotation_df, out_path)
        assert out_path.exists()
        assert out_path.stat().st_size > 1000


class TestBuildContactSheet:
    def test_raises_on_odd_sample_count(self, tmp_path):
        with pytest.raises(ValueError, match="even number"):
            build_contact_sheet(("Tumor_01_Control",), tmp_path, tmp_path / "sheet.png")

    def test_raises_on_non_adjacent_pair(self, tmp_path):
        with pytest.raises(ValueError, match="adjacent"):
            build_contact_sheet(("Tumor_01_Control", "Tumor_02_Tamoxifen"), tmp_path, tmp_path / "sheet.png")

    def test_raises_on_repeated_tumor(self, tmp_path):
        samples = ("Tumor_01_Control", "Tumor_01_Tamoxifen", "Tumor_01_Control", "Tumor_01_Tamoxifen")
        with pytest.raises(ValueError, match="more than one pair"):
            build_contact_sheet(samples, tmp_path, tmp_path / "sheet.png")

    def test_single_pair_does_not_crash_on_axes_indexing(self, tmp_path):
        """Regression test for a squeeze=False fix: plt.subplots(1, 2)
        returns a 1-D axes array by default, which would break
        axes[row, col] indexing for exactly one Control/Tamoxifen pair."""
        import matplotlib.pyplot as plt

        for sample_id in ("Tumor_01_Control", "Tumor_01_Tamoxifen"):
            fig, ax = plt.subplots(figsize=(1, 1))
            ax.plot([0, 1])
            fig.savefig(tmp_path / f"{sample_id}_heatmap_comparison.png")
            plt.close(fig)

        build_contact_sheet(("Tumor_01_Control", "Tumor_01_Tamoxifen"), tmp_path, tmp_path / "sheet.png")
        assert (tmp_path / "sheet.png").exists()


class TestRealDataSmokeTest:
    def test_tumor04_end_to_end_matches_point1_counts(self):
        """Runs the real (non-R) part of the pipeline on the smallest
        sample pair and cross-checks this module's independently-computed
        agreement counts against the already-verified
        cnv_method_comparison_by_sample.tsv (Point 1 output) for the same
        two samples -- two separately-written computations over the same
        frozen label tables, which must agree exactly."""
        cfg = _real_cfg()
        working_dir_exists = (cfg.working_dir / "Tumor_04_Control").exists()
        if not working_dir_exists:
            pytest.skip("compact matrices not extracted in this checkout (R extraction step not run)")

        result = run_cnv_heatmap_comparison(run_r_extraction_step=False, samples=("Tumor_04_Control", "Tumor_04_Tamoxifen"))
        for sample_id, path in result["figure_paths"].items():
            assert Path(path).exists()
            assert Path(path).stat().st_size > 1000

        infercnv_df = pd.read_csv(cfg.infercnv_labels_tsv, sep="\t")
        copykat_df = pd.read_csv(cfg.copykat_labels_tsv, sep="\t")
        by_sample = pd.read_csv(cfg.method_comparison_by_sample_tsv, sep="\t").set_index("sample_id")

        for sample_id in ("Tumor_04_Control", "Tumor_04_Tamoxifen"):
            annotation_df = build_annotation_table(sample_id, infercnv_df, copykat_df)
            counts = annotation_df["agreement_status"].value_counts()
            expected = by_sample.loc[sample_id]
            assert int(counts.get(AGREE_TUMOR_LIKE, 0)) == int(expected["both_malignant_aneuploid"])
            assert int(counts.get(AGREE_NORMAL_LIKE, 0)) == int(expected["both_nonmalignant_diploid"])
            assert int(counts.get(COPYKAT_NOT_DEFINED_STATUS, 0)) == int(expected["copykat_not_defined_count"])

    def test_summary_table_matches_point1_fractions(self):
        cfg = _real_cfg()
        if not cfg.gene_exclusion_tsv.exists() or not cfg.method_comparison_by_sample_tsv.exists():
            pytest.skip("gene exclusion report or Point-1 by-sample table not present in this checkout")
        summary = build_summary_table(cfg)
        assert list(summary["sample"]) == list(cfg.selected_samples)

        by_sample = pd.read_csv(cfg.method_comparison_by_sample_tsv, sep="\t").set_index("sample_id")
        for row in summary.itertuples(index=False):
            expected = by_sample.loc[row.sample]
            assert row.n_matched_cells == int(expected["total_epithelial_cells"])
            assert row.infercnv_malignant_fraction == pytest.approx(expected["infercnv_malignant_pct"] / 100.0)
            assert row.copykat_aneuploid_fraction == pytest.approx(expected["copykat_aneuploid_pct"] / 100.0)
            assert row.agreement_fraction == pytest.approx(expected["agreement_pct"] / 100.0)

    def test_load_compact_matrices_real_data_if_present(self):
        cfg = _real_cfg()
        if not (cfg.working_dir / "Tumor_04_Control").exists():
            pytest.skip("compact matrices not extracted in this checkout")
        ic, ck, gene_order = load_compact_matrices("Tumor_04_Control", cfg)
        assert ic.shape[0] == len(gene_order)
        assert ck.shape[0] == len(gene_order)
        assert set(ic.columns) == set(ck.columns)
        assert not ic.isna().any().any()
        assert not ck.isna().any().any()
