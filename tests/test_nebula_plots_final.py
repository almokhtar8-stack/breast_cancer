import ast
from pathlib import Path

import pandas as pd
import pytest
import yaml
from PIL import Image

from src.candidate_evidence_summary import EVIDENCE_CLASS_PRIMARY
from src.nebula_plots_final import (
    NebulaPlotsFinalConfig,
    build_fig4_heatmap_inputs,
    check_fig4_layout,
    check_fig5_layout,
    contiguous_evidence_blocks,
    fig4_block_boundaries,
    fig4_col_geometry,
    fig4_row_geometry,
    fig4_sample_group_boundaries,
    load_fig1_input,
    load_fig1_paics_inset,
    load_fig2_input,
    load_fig3_input,
    load_fig4_evidence_input,
    load_filtered_tpm,
    load_fig5_expression_input,
    load_fig5_summary_input,
    load_sample_metadata,
    plot_fig1,
    plot_fig2,
    plot_fig3,
    plot_fig4,
    plot_fig5,
    run_nebula_plots_final,
)

REPO_ROOT = Path(__file__).parent.parent
CONFIG_PATH = REPO_ROOT / "config" / "config.yaml"
FINAL_SRC = REPO_ROOT / "src" / "nebula_plots_final.py"

# Frozen expected pixel dimensions of the visually-approved figure set
# (matching v6 Figures 1-3 and v8 Figures 4-5 exactly -- confirmed
# byte-identical at approval time).
EXPECTED_DIMENSIONS = {
    "fig1_crispr_landscape.png": (2967, 2833),
    "fig2_pca.png": (2102, 1613),
    "fig3_volcano.png": (2354, 1918),
    "fig4_candidate_expression_heatmap.png": (4638, 3328),
    "fig5_usp34_expression.png": (2669, 1664),
}


def _real_cfg() -> NebulaPlotsFinalConfig:
    with open(CONFIG_PATH) as f:
        config = yaml.safe_load(f)
    return NebulaPlotsFinalConfig.from_config(config)


# --- Independence: no dependency on experimental v2-v8 modules --------------


class TestModuleIndependence:
    def test_no_v2_through_v8_imports(self):
        with open(FINAL_SRC) as f:
            tree = ast.parse(f.read())
        imported_modules = []
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                imported_modules.extend(n.name for n in node.names)
            elif isinstance(node, ast.ImportFrom):
                if node.module:
                    imported_modules.append(node.module)
        forbidden = [m for m in imported_modules if m and "nebula_plots_v" in m]
        assert forbidden == [], f"nebula_plots_final.py must not import experimental modules, found: {forbidden}"

    def test_config_has_no_experimental_sections(self):
        with open(CONFIG_PATH) as f:
            config = yaml.safe_load(f)
        for version in ("nebula_plots_v2", "nebula_plots_v3", "nebula_plots_v4", "nebula_plots_v5", "nebula_plots_v6", "nebula_plots_v7", "nebula_plots_v8"):
            assert version not in config, f"experimental config section {version!r} should have been removed"
        assert "nebula_plots_final" in config
        assert "nebula_plots" in config, "the original committed v1 section must remain untouched"


# --- Fig4 shared row/column geometry (approved v8 alignment fix) ------------


class TestFigure4SharedGeometry:
    def test_row_geometry_matches_imshow_convention(self):
        geom = fig4_row_geometry(13)
        assert list(geom["centers"]) == [float(i) for i in range(13)]
        assert geom["ylim"] == (12.5, -0.5)

    def test_block_boundaries_from_fixed_block_sizes(self):
        genes = (
            ["USP34"]
            + ["CTDNEP1", "EIF4ENIF1", "HMGB1", "KDM1A", "PET117", "TADA2B", "VEZF1"]
            + ["ICK", "SUPT4H1", "TLK2", "TSR3"]
            + ["USP17L29"]
        )
        classes = (
            [EVIDENCE_CLASS_PRIMARY]
            + ["SECONDARY_CONTEXT_SUPPORT"] * 7
            + ["NO_SIGNIFICANT_RNA_SUPPORT"] * 4
            + ["RNA_UNAVAILABLE"]
        )
        evidence_class = pd.Series(classes, index=genes)
        blocks = contiguous_evidence_blocks(evidence_class)
        assert [end - start for start, end, _cls in blocks] == [1, 7, 4, 1]
        geom = fig4_row_geometry(13)
        assert fig4_block_boundaries(blocks, geom["edges"]) == [0.5, 7.5, 11.5]

    def test_column_geometry_and_group_boundaries(self):
        geom = fig4_col_geometry(9)
        assert geom["xlim"] == (-0.5, 8.5)
        assert fig4_sample_group_boundaries(9, geom["edges"]) == [2.5, 5.5]


# --- Fig4/Fig5 rendered layout: no collisions, aligned boundaries -----------


class TestFigure4And5RenderedLayout:
    @staticmethod
    @pytest.fixture(scope="class")
    def rendered_fig4():
        cfg = _real_cfg()
        evidence_df = load_fig4_evidence_input(cfg)
        filtered_tpm_df = load_filtered_tpm(cfg)
        sample_meta_df = load_sample_metadata(cfg)
        fig5_summary_df = load_fig5_summary_input(cfg)
        z_score_df, available_mask, evidence_class = build_fig4_heatmap_inputs(evidence_df, filtered_tpm_df, sample_meta_df)
        fig, artifacts = plot_fig4(z_score_df, available_mask, evidence_class, sample_meta_df, fig5_summary_df, cfg, cfg.expected_primary_gene)
        return fig, artifacts

    @staticmethod
    @pytest.fixture(scope="class")
    def rendered_fig5():
        cfg = _real_cfg()
        summary_df = load_fig5_summary_input(cfg)
        expression_df = load_fig5_expression_input(cfg)
        fig, artifacts = plot_fig5(summary_df, expression_df, cfg, cfg.expected_primary_gene)
        return fig, artifacts

    def test_fig4_no_layout_violations(self, rendered_fig4):
        fig, artifacts = rendered_fig4
        assert check_fig4_layout(fig, artifacts) == []

    def test_fig4_row_axes_share_identical_ylim(self, rendered_fig4):
        _fig, artifacts = rendered_fig4
        ylims = list(artifacts["row_ylims"].values())
        assert all(v == ylims[0] for v in ylims)

    def test_fig4_separator_positions_identical_across_regions(self, rendered_fig4):
        _fig, artifacts = rendered_fig4
        positions = list(artifacts["separator_positions"].values())
        assert all(p == positions[0] for p in positions)
        assert positions[0] == [0.5, 7.5, 11.5]

    def test_fig4_group_strip_matches_heatmap_columns(self, rendered_fig4):
        _fig, artifacts = rendered_fig4
        assert artifacts["col_boundaries"]["group_strip"] == artifacts["col_boundaries"]["heat"]

    def test_fig5_no_layout_violations(self, rendered_fig5):
        fig, artifacts = rendered_fig5
        assert check_fig5_layout(fig, artifacts) == []

    def test_fig5_comparison_annotations_both_present_and_distinct(self, rendered_fig5):
        _fig, artifacts = rendered_fig5
        mcf7_text = artifacts["mcf7_tamr_bracket_text"].get_text()
        fasr_text = artifacts["tamr_fasr_bracket_text"].get_text()
        assert "TAMR vs MCF7" in mcf7_text
        assert "TAMR vs FASR" in fasr_text
        assert "n.s." in fasr_text


# --- Scientific invariants (re-verified against the frozen, committed inputs)-


class TestScientificInvariants:
    def test_fig1_hits_and_paics_separation(self):
        cfg = _real_cfg()
        fig1_df = load_fig1_input(cfg)
        paics_df = load_fig1_paics_inset(cfg)
        assert len(fig1_df) == 28
        assert (fig1_df["crispr_effect_size"] < 0).sum() == 13
        assert (fig1_df["crispr_effect_size"] > 0).sum() == 15
        assert "PAICS" not in set(fig1_df["gene_symbol"])
        assert paics_df.iloc[0]["gene_symbol"] == "PAICS"
        assert paics_df.iloc[0]["benchmark_label"] == "published_benchmark_not_gate1_hit"

    def test_fig2_pca_sample_and_variance_invariants(self):
        cfg = _real_cfg()
        fig2_df = load_fig2_input(cfg)
        assert fig2_df["sample_id"].nunique() == 9
        assert (fig2_df.drop_duplicates("sample_id")["group"].value_counts() == 3).all()
        assert round(float(fig2_df["pc1_variance_explained_pct"].iloc[0]), 1) == 59.8
        assert round(float(fig2_df["pc2_variance_explained_pct"].iloc[0]), 1) == 28.8

    def test_fig3_usp34_frozen_limma_values(self):
        cfg = _real_cfg()
        fig3_df = load_fig3_input(cfg)
        usp34_row = fig3_df.loc[fig3_df["gene_symbol"] == "USP34"].iloc[0]
        assert round(float(usp34_row["log2fc"]), 2) == 0.59
        assert f"{usp34_row['fdr']:.3g}" == "0.00731"

    def test_fig4_thirteen_sensitising_twelve_available(self):
        cfg = _real_cfg()
        evidence_df = load_fig4_evidence_input(cfg)
        filtered_tpm_df = load_filtered_tpm(cfg)
        sample_meta_df = load_sample_metadata(cfg)
        z_score_df, available_mask, evidence_class = build_fig4_heatmap_inputs(evidence_df, filtered_tpm_df, sample_meta_df)

        assert len(evidence_df) == 13
        assert int(available_mask.sum()) == 12
        assert not available_mask["USP17L29"]
        assert list(z_score_df.columns) == [
            "MCF7_Rep1", "MCF7_Rep2", "MCF7_Rep3",
            "TAMR_Rep1", "TAMR_Rep2", "TAMR_Rep3",
            "FASR_Rep1", "FASR_Rep2", "FASR_Rep3",
        ]
        blocks = contiguous_evidence_blocks(evidence_class)
        assert [end - start for start, end, _cls in blocks] == [1, 7, 4, 1]
        primary_rows = evidence_class[evidence_class == EVIDENCE_CLASS_PRIMARY]
        assert list(primary_rows.index) == ["USP34"]

    def test_fig5_frozen_summary_values(self):
        cfg = _real_cfg()
        summary_df = load_fig5_summary_input(cfg)
        row = summary_df.iloc[0]
        assert round(float(row["crispr_effect_size"]), 2) == -1.39
        assert f"{row['crispr_fdr']:.3g}" == "0.0417"
        assert round(float(row["tamr_vs_mcf7_log2fc"]), 2) == 0.59
        assert f"{row['tamr_vs_mcf7_fdr']:.3g}" == "0.00731"
        assert round(float(row["tamr_vs_fasr_log2fc"]), 2) == 0.02
        assert f"{row['tamr_vs_fasr_fdr']:.3g}" == "0.912"


# --- Real-repository smoke test: dimensions, filenames, manifest ------------


class TestRunAgainstRealConfig:
    def test_final_figures_written_with_expected_filenames_and_dimensions(self):
        result = run_nebula_plots_final()
        cfg = _real_cfg()

        expected_files = [
            cfg.fig1_png, cfg.fig1_pdf, cfg.fig2_png, cfg.fig2_pdf, cfg.fig3_png, cfg.fig3_pdf,
            cfg.fig4_png, cfg.fig4_pdf, cfg.fig5_png, cfg.fig5_pdf,
        ]
        for path in expected_files:
            assert path.exists(), f"{path} was not written"
            assert path.stat().st_size > 0

        for filename, expected_dims in EXPECTED_DIMENSIONS.items():
            path = cfg.output_dir / filename
            with Image.open(path) as image:
                assert image.size == expected_dims, f"{filename}: expected {expected_dims}, got {image.size}"

        assert cfg.contact_sheet_png.exists()
        assert cfg.manifest_tsv.exists()

        manifest_df = result["manifest"]
        assert len(manifest_df) == 5
        assert set(manifest_df["approved_source_version"]) == {"v6", "v8"}
        fig123 = manifest_df[manifest_df["figure"].isin(["fig1_crispr_landscape", "fig2_pca", "fig3_volcano"])]
        assert (fig123["approved_source_version"] == "v6").all()
        fig45 = manifest_df[manifest_df["figure"].isin(["fig4_candidate_expression_heatmap", "fig5_usp34_expression"])]
        assert (fig45["approved_source_version"] == "v8").all()

    def test_manifest_has_required_columns(self):
        result = run_nebula_plots_final()
        manifest_df = result["manifest"]
        required_columns = {"figure", "filename", "approved_source_version", "width_px", "height_px", "sha256", "scientific_source", "purpose"}
        assert required_columns.issubset(set(manifest_df.columns))

    def test_no_transparent_output_files(self):
        cfg = _real_cfg()
        run_nebula_plots_final()
        for path in [cfg.fig1_png, cfg.fig2_png, cfg.fig3_png, cfg.fig4_png, cfg.fig5_png]:
            image = Image.open(path)
            if image.mode in ("RGBA", "LA"):
                alpha = image.getchannel("A")
                assert alpha.getextrema() == (255, 255), f"{path} has transparent pixels"
