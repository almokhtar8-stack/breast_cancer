import copy
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import yaml

from src.candidate_evidence_summary import (
    BENCHMARK_GENE_SYMBOL,
    BENCHMARK_LABEL,
    DIRECTION_SENSITISING,
    DIRECTION_TOLERANCE,
    EVIDENCE_CLASS_NO_SIGNIFICANT_RNA,
    EVIDENCE_CLASS_PRIMARY,
    EVIDENCE_CLASS_RNA_UNAVAILABLE,
    EVIDENCE_CLASS_SECONDARY,
)
from src.nebula_plots import (
    NebulaPlotsConfig,
    _validate_direction_matches_sign,
    build_fig1_inputs,
    build_fig2_inputs,
    build_fig3_inputs,
    build_fig4_inputs,
    build_fig5_expression_input,
    build_fig5_summary_input,
    load_evidence_summary,
    load_paics_benchmark,
    load_pca_coordinates,
    load_sample_metadata,
    load_shortlist,
    plot_fig3_volcano,
    run_nebula_plots,
    run_qc_checks,
)

REPO_ROOT = Path(__file__).parent.parent
CONFIG_PATH = REPO_ROOT / "config" / "config.yaml"


def _evidence_row(
    gene_symbol,
    crispr_effect_size,
    evidence_class,
    crispr_fdr=0.01,
    mcf7_log2fc=0.5,
    mcf7_fdr=0.01,
    fasr_log2fc=0.5,
    fasr_fdr=0.01,
    de_na_reason=np.nan,
):
    direction = DIRECTION_SENSITISING if crispr_effect_size < 0 else DIRECTION_TOLERANCE
    return {
        "gene_symbol": gene_symbol,
        "crispr_effect_size": crispr_effect_size,
        "crispr_fdr": crispr_fdr,
        "crispr_direction": direction,
        "crispr_direction_description": "desc",
        "tamr_vs_mcf7_log2fc": mcf7_log2fc if not (isinstance(de_na_reason, str)) else np.nan,
        "tamr_vs_mcf7_fdr": mcf7_fdr if not (isinstance(de_na_reason, str)) else np.nan,
        "tamr_vs_fasr_log2fc": fasr_log2fc if not (isinstance(de_na_reason, str)) else np.nan,
        "tamr_vs_fasr_fdr": fasr_fdr if not (isinstance(de_na_reason, str)) else np.nan,
        "mcf7_baseline_log2_tpm_plus1": 5.0,
        "gse118713_de_na_reason": de_na_reason,
        "evidence_class": evidence_class,
        "evidence_class_label": "label",
        "primary_rna_direction": "not_applicable",
        "rna_pattern_description": "pattern",
    }


def _cfg(tmp_path, **overrides):
    defaults = dict(
        evidence_summary_tsv=tmp_path / "es.tsv",
        shortlist_tsv=tmp_path / "shortlist.tsv",
        sensitisation_candidates_tsv=tmp_path / "sens.tsv",
        tolerance_hits_tsv=tmp_path / "tol.tsv",
        paics_benchmark_tsv=tmp_path / "paics.tsv",
        pca_coordinates_tsv=tmp_path / "pca.tsv",
        differential_expression_tsv_gz=tmp_path / "de2.tsv.gz",
        filtered_gene_tpm_tsv_gz=tmp_path / "tpm.tsv.gz",
        frozen_filtered_gene_tpm_sha256="deadbeef",
        sample_metadata_tsv=tmp_path / "meta.tsv",
        rna_significance_fdr=0.05,
        expected_n_hits=28,
        expected_n_sensitising=13,
        expected_n_tolerance=15,
        expected_n_samples=9,
        expected_primary_gene="USP34",
        poster_background="#0B0B12",
        panel_background="#F8FAFC",
        palette={
            "nebula_purple": "#6E44FF",
            "soft_lilac": "#B799FF",
            "cosmic_magenta": "#D946EF",
            "rose_pink": "#F472B6",
            "neutral_grey": "#94A3B8",
        },
        group_colors={"MCF7": "#B799FF", "TAMR": "#6E44FF", "FASR": "#F472B6"},
        direction_colors={DIRECTION_SENSITISING: "#6E44FF", DIRECTION_TOLERANCE: "#F472B6"},
        output_dir=tmp_path / "figs",
        plot_input_dir=tmp_path / "plot_inputs",
        fig1_png=tmp_path / "f1.png",
        fig1_png_transparent=tmp_path / "f1t.png",
        fig1_pdf=tmp_path / "f1.pdf",
        fig1_input_tsv=tmp_path / "f1.tsv",
        fig1_paics_inset_tsv=tmp_path / "f1p.tsv",
        fig2_png=tmp_path / "f2.png",
        fig2_png_transparent=tmp_path / "f2t.png",
        fig2_pdf=tmp_path / "f2.pdf",
        fig2_input_tsv=tmp_path / "f2.tsv",
        fig3_png=tmp_path / "f3.png",
        fig3_png_transparent=tmp_path / "f3t.png",
        fig3_pdf=tmp_path / "f3.pdf",
        fig3_input_tsv=tmp_path / "f3.tsv",
        fig4_png=tmp_path / "f4.png",
        fig4_png_transparent=tmp_path / "f4t.png",
        fig4_pdf=tmp_path / "f4.pdf",
        fig4_input_tsv=tmp_path / "f4.tsv",
        fig5_png=tmp_path / "f5.png",
        fig5_png_transparent=tmp_path / "f5t.png",
        fig5_pdf=tmp_path / "f5.pdf",
        fig5_summary_input_tsv=tmp_path / "f5s.tsv",
        fig5_expression_input_tsv=tmp_path / "f5e.tsv",
    )
    defaults.update(overrides)
    return NebulaPlotsConfig(**defaults)


def _make_evidence_df():
    rows = []
    rows.append(_evidence_row("USP34", -1.4, EVIDENCE_CLASS_PRIMARY, mcf7_fdr=0.007))
    for i, gene in enumerate(["CTDNEP1", "EIF4ENIF1", "HMGB1", "KDM1A", "PET117", "TADA2B", "VEZF1"]):
        rows.append(_evidence_row(gene, -1.5 - i * 0.1, EVIDENCE_CLASS_SECONDARY, mcf7_fdr=0.5, fasr_fdr=0.01))
    for i, gene in enumerate(["ICK", "SUPT4H1", "TLK2", "TSR3"]):
        rows.append(_evidence_row(gene, -1.2 - i * 0.1, EVIDENCE_CLASS_NO_SIGNIFICANT_RNA, mcf7_fdr=0.5, fasr_fdr=0.5))
    rows.append(_evidence_row("USP17L29", -2.1, EVIDENCE_CLASS_RNA_UNAVAILABLE, de_na_reason="filtered_out_no_de_fit"))
    for i, gene in enumerate(
        ["PTEN", "CSK", "NF1", "KMT2C", "SOX2", "LZTR1", "DPP9", "KDM6A", "RASA2", "CUX1", "TFAP2C", "PTPN3", "ILK", "SPRED2", "TSC1"]
    ):
        rows.append(_evidence_row(gene, 1.0 + i * 0.1, "not_applicable_not_a_sensitising_knockout"))
    return pd.DataFrame(rows)


# --- build_fig1_inputs -------------------------------------------------


class TestBuildFig1Inputs:
    def test_sorted_by_effect_size_ascending_no_score(self):
        evidence_df = _make_evidence_df()
        paics_df = pd.DataFrame(
            [{"gene_symbol": "PAICS", "crispr_effect_size": -0.3, "crispr_fdr": 0.85, "crispr_direction": DIRECTION_SENSITISING, "benchmark_label": BENCHMARK_LABEL}]
        )
        fig1_df, paics_inset_df = build_fig1_inputs(evidence_df, paics_df)
        assert len(fig1_df) == 28
        assert list(fig1_df["crispr_effect_size"]) == sorted(fig1_df["crispr_effect_size"])
        assert "PAICS" not in set(fig1_df["gene_symbol"])
        assert len(paics_inset_df) == 1
        assert paics_inset_df.iloc[0]["gene_symbol"] == "PAICS"
        forbidden = {"score", "composite_score", "rank", "weighted_score"}
        assert forbidden.isdisjoint(set(fig1_df.columns))


# --- build_fig2_inputs ---------------------------------------------------


class TestBuildFig2Inputs:
    def test_pivots_and_carries_variance(self):
        pca_df = pd.DataFrame(
            [
                {"sample_id": "MCF7_Rep1", "group": "MCF7", "pc": "PC1", "coordinate": 1.0, "variance_explained_fraction": 0.598},
                {"sample_id": "MCF7_Rep1", "group": "MCF7", "pc": "PC2", "coordinate": 2.0, "variance_explained_fraction": 0.288},
                {"sample_id": "TAMR_Rep1", "group": "TAMR", "pc": "PC1", "coordinate": -1.0, "variance_explained_fraction": 0.598},
                {"sample_id": "TAMR_Rep1", "group": "TAMR", "pc": "PC2", "coordinate": -2.0, "variance_explained_fraction": 0.288},
            ]
        )
        fig2_df = build_fig2_inputs(pca_df)
        assert set(fig2_df["sample_id"]) == {"MCF7_Rep1", "TAMR_Rep1"}
        assert fig2_df["pc1_variance_explained_pct"].iloc[0] == pytest.approx(59.8, abs=0.05)
        assert fig2_df["pc2_variance_explained_pct"].iloc[0] == pytest.approx(28.8, abs=0.05)


# --- build_fig3_inputs -----------------------------------------------------


class TestBuildFig3Inputs:
    def test_left_join_flags_gate1_hits_without_dropping_background(self):
        background = pd.DataFrame(
            [
                {"gene_id": "E1", "gene_symbol": "USP34", "log2fc": 0.59, "fdr": 0.007, "contrast": "TAMR_vs_MCF7"},
                {"gene_id": "E2", "gene_symbol": "RANDOMGENE", "log2fc": 0.1, "fdr": 0.9, "contrast": "TAMR_vs_MCF7"},
            ]
        )
        evidence_df = pd.DataFrame(
            [_evidence_row("USP34", -1.4, EVIDENCE_CLASS_PRIMARY, mcf7_fdr=0.007)]
        )
        fig3_df = build_fig3_inputs(background, evidence_df)
        assert len(fig3_df) == 2
        indexed = fig3_df.set_index("gene_symbol")
        assert indexed.loc["USP34", "is_gate1_hit"]
        assert not indexed.loc["RANDOMGENE", "is_gate1_hit"]

    # --- fix 3: FDR must be finite and in (0, 1] before -log10 transform --

    def _background(self, fdr):
        return pd.DataFrame([{"gene_id": "E1", "gene_symbol": "SOME_GENE", "log2fc": 0.1, "fdr": fdr, "contrast": "TAMR_vs_MCF7"}])

    def _empty_evidence(self):
        return pd.DataFrame(columns=["gene_symbol", "crispr_direction", "evidence_class"])

    def test_raises_on_non_finite_fdr(self):
        with pytest.raises(ValueError, match="non-finite FDR"):
            build_fig3_inputs(self._background(float("nan")), self._empty_evidence())

    def test_raises_on_zero_fdr(self):
        with pytest.raises(ValueError, match=r"outside \(0, 1\]"):
            build_fig3_inputs(self._background(0.0), self._empty_evidence())

    def test_raises_on_negative_fdr(self):
        with pytest.raises(ValueError, match=r"outside \(0, 1\]"):
            build_fig3_inputs(self._background(-0.01), self._empty_evidence())

    def test_raises_on_fdr_above_one(self):
        with pytest.raises(ValueError, match=r"outside \(0, 1\]"):
            build_fig3_inputs(self._background(1.5), self._empty_evidence())

    def test_accepts_fdr_exactly_one(self):
        fig3_df = build_fig3_inputs(self._background(1.0), self._empty_evidence())
        assert len(fig3_df) == 1


# --- fix 6: secondary-context genes cannot be misread as significant in --
# --- this (TAMR_vs_MCF7) contrast -----------------------------------------


class TestFig3SecondaryContextWording:
    def _fig3_df(self):
        background = pd.DataFrame(
            [
                {"gene_id": "E1", "gene_symbol": "USP34", "log2fc": 0.59, "fdr": 0.007, "contrast": "TAMR_vs_MCF7"},
                {"gene_id": "E2", "gene_symbol": "CTDNEP1", "log2fc": -0.3, "fdr": 0.5, "contrast": "TAMR_vs_MCF7"},
                {"gene_id": "E3", "gene_symbol": "BACKGROUND_GENE", "log2fc": 0.1, "fdr": 0.9, "contrast": "TAMR_vs_MCF7"},
            ]
        )
        evidence_df = pd.DataFrame(
            [
                _evidence_row("USP34", -1.4, EVIDENCE_CLASS_PRIMARY, mcf7_fdr=0.007),
                _evidence_row("CTDNEP1", -1.0, EVIDENCE_CLASS_SECONDARY, mcf7_fdr=0.5, fasr_fdr=0.01),
            ]
        )
        return build_fig3_inputs(background, evidence_df)

    def test_secondary_gene_label_carries_a_marker_and_a_caption_explains_it(self, tmp_path):
        fig3_df = self._fig3_df()
        cfg = _cfg(tmp_path)
        fig = plot_fig3_volcano(fig3_df, cfg, "USP34", ["CTDNEP1"])
        ax = fig.axes[0]

        annotated_texts = [a.get_text() for a in ax.texts]
        # The secondary-context gene's on-plot label must carry a marker
        # (not a bare gene symbol) distinguishing it from a gene significant
        # in *this* TAMR_vs_MCF7 contrast.
        assert any("CTDNEP1" in t and t != "CTDNEP1" for t in annotated_texts)
        # A caption using wording equivalent to "Secondary-context CRISPR
        # candidates" must be present and must explain the marker.
        assert any("Secondary-context CRISPR candidates" in t for t in annotated_texts)
        assert any("not necessarily significant in this TAMR_vs_MCF7" in t for t in annotated_texts)

    def test_primary_gene_label_has_no_secondary_marker(self, tmp_path):
        fig3_df = self._fig3_df()
        cfg = _cfg(tmp_path)
        fig = plot_fig3_volcano(fig3_df, cfg, "USP34", ["CTDNEP1"])
        ax = fig.axes[0]
        annotated_texts = [a.get_text() for a in ax.texts]
        assert "USP34" in annotated_texts
        assert "USP34†" not in annotated_texts

    def test_no_caption_when_there_are_no_secondary_genes(self, tmp_path):
        fig3_df = self._fig3_df()
        cfg = _cfg(tmp_path)
        fig = plot_fig3_volcano(fig3_df, cfg, "USP34", [])
        ax = fig.axes[0]
        annotated_texts = [a.get_text() for a in ax.texts]
        assert not any("Secondary-context CRISPR candidates" in t for t in annotated_texts)


# --- build_fig4_inputs: deterministic evidence-class then alphabetical ----


class TestBuildFig4Inputs:
    def test_ordering_is_class_then_alphabetical(self):
        sensitising_df = pd.DataFrame(
            [
                _evidence_row("ZZZ_SECONDARY", -1.0, EVIDENCE_CLASS_SECONDARY, fasr_fdr=0.01, mcf7_fdr=0.5),
                _evidence_row("AAA_SECONDARY", -1.0, EVIDENCE_CLASS_SECONDARY, fasr_fdr=0.01, mcf7_fdr=0.5),
                _evidence_row("USP34", -1.4, EVIDENCE_CLASS_PRIMARY, mcf7_fdr=0.007),
                _evidence_row("NO_RNA_GENE", -1.0, EVIDENCE_CLASS_NO_SIGNIFICANT_RNA, mcf7_fdr=0.5, fasr_fdr=0.5),
            ]
        )
        fig4_df = build_fig4_inputs(sensitising_df)
        assert list(fig4_df["gene_symbol"]) == ["USP34", "AAA_SECONDARY", "ZZZ_SECONDARY", "NO_RNA_GENE"]


# --- build_fig5_summary_input / build_fig5_expression_input ---------------


class TestBuildFig5Inputs:
    def test_summary_input_pulls_exact_row(self):
        evidence_df = _make_evidence_df()
        summary_df = build_fig5_summary_input(evidence_df, "USP34")
        assert len(summary_df) == 1
        assert summary_df.iloc[0]["crispr_effect_size"] == -1.4

    def test_summary_input_raises_if_gene_missing(self):
        evidence_df = _make_evidence_df()
        with pytest.raises(ValueError, match="expected exactly one"):
            build_fig5_summary_input(evidence_df, "NOT_A_GENE")

    def test_expression_input_uses_actual_tpm_not_inferred_from_log2fc(self):
        filtered_tpm_df = pd.DataFrame(
            [{"gene_id": "E1", "gene_symbol": "USP34", "symbol_mapping_status": "resolved",
              "MCF7_Rep1": 28.97, "MCF7_Rep2": 36.58, "MCF7_Rep3": 35.5,
              "TAMR_Rep1": 48.56, "TAMR_Rep2": 51.06, "TAMR_Rep3": 53.39,
              "FASR_Rep1": 49.08, "FASR_Rep2": 49.55, "FASR_Rep3": 52.23}]
        )
        sample_meta_df = pd.DataFrame(
            [
                {"sample_id": "MCF7_Rep1", "group": "MCF7", "replicate": 1},
                {"sample_id": "MCF7_Rep2", "group": "MCF7", "replicate": 2},
                {"sample_id": "MCF7_Rep3", "group": "MCF7", "replicate": 3},
                {"sample_id": "TAMR_Rep1", "group": "TAMR", "replicate": 1},
                {"sample_id": "TAMR_Rep2", "group": "TAMR", "replicate": 2},
                {"sample_id": "TAMR_Rep3", "group": "TAMR", "replicate": 3},
                {"sample_id": "FASR_Rep1", "group": "FASR", "replicate": 1},
                {"sample_id": "FASR_Rep2", "group": "FASR", "replicate": 2},
                {"sample_id": "FASR_Rep3", "group": "FASR", "replicate": 3},
            ]
        )
        expr_df = build_fig5_expression_input(filtered_tpm_df, sample_meta_df, "USP34")
        assert len(expr_df) == 9
        row = expr_df.loc[expr_df["sample_id"] == "MCF7_Rep1"].iloc[0]
        assert row["tpm"] == pytest.approx(28.97)
        assert row["log2_tpm_plus1"] == pytest.approx(np.log2(29.97))

    def test_expression_input_raises_if_sample_column_missing(self):
        filtered_tpm_df = pd.DataFrame(
            [{"gene_id": "E1", "gene_symbol": "USP34", "MCF7_Rep1": 10.0}]
        )
        sample_meta_df = pd.DataFrame([{"sample_id": "MISSING_SAMPLE", "group": "MCF7", "replicate": 1}])
        with pytest.raises(ValueError, match="not found"):
            build_fig5_expression_input(filtered_tpm_df, sample_meta_df, "USP34")

    # --- fix 4: TPM must be finite and non-negative before log2(TPM+1) ----

    def test_raises_on_non_finite_tpm(self):
        filtered_tpm_df = pd.DataFrame([{"gene_id": "E1", "gene_symbol": "USP34", "MCF7_Rep1": float("nan")}])
        sample_meta_df = pd.DataFrame([{"sample_id": "MCF7_Rep1", "group": "MCF7", "replicate": 1}])
        with pytest.raises(ValueError, match="not finite"):
            build_fig5_expression_input(filtered_tpm_df, sample_meta_df, "USP34")

    def test_raises_on_negative_tpm(self):
        filtered_tpm_df = pd.DataFrame([{"gene_id": "E1", "gene_symbol": "USP34", "MCF7_Rep1": -1.0}])
        sample_meta_df = pd.DataFrame([{"sample_id": "MCF7_Rep1", "group": "MCF7", "replicate": 1}])
        with pytest.raises(ValueError, match="negative"):
            build_fig5_expression_input(filtered_tpm_df, sample_meta_df, "USP34")

    def test_accepts_zero_tpm(self):
        filtered_tpm_df = pd.DataFrame([{"gene_id": "E1", "gene_symbol": "USP34", "MCF7_Rep1": 0.0}])
        sample_meta_df = pd.DataFrame([{"sample_id": "MCF7_Rep1", "group": "MCF7", "replicate": 1}])
        expr_df = build_fig5_expression_input(filtered_tpm_df, sample_meta_df, "USP34")
        assert expr_df.iloc[0]["tpm"] == 0.0
        assert expr_df.iloc[0]["log2_tpm_plus1"] == pytest.approx(0.0)


# --- run_qc_checks ----------------------------------------------------------


class TestRunQcChecks:
    def _pieces(self, tmp_path):
        evidence_df = _make_evidence_df()
        sensitising_df = evidence_df.loc[evidence_df["crispr_direction"] == DIRECTION_SENSITISING].reset_index(drop=True)
        tolerance_df = evidence_df.loc[evidence_df["crispr_direction"] == DIRECTION_TOLERANCE].reset_index(drop=True)
        shortlist_df = pd.DataFrame([{"gene_symbol": "USP34"}])
        paics_df = pd.DataFrame([{"gene_symbol": BENCHMARK_GENE_SYMBOL, "benchmark_label": BENCHMARK_LABEL}])
        sample_ids = [f"MCF7_Rep{i}" for i in range(1, 4)] + [f"TAMR_Rep{i}" for i in range(1, 4)] + [f"FASR_Rep{i}" for i in range(1, 4)]
        groups = ["MCF7"] * 3 + ["TAMR"] * 3 + ["FASR"] * 3
        pca_df = pd.DataFrame({"sample_id": sample_ids, "group": groups})
        sample_meta_df = pd.DataFrame({"sample_id": sample_ids, "group": groups})
        cfg = _cfg(tmp_path)
        return evidence_df, shortlist_df, sensitising_df, tolerance_df, paics_df, pca_df, sample_meta_df, cfg

    def test_passes_on_valid_data(self, tmp_path):
        run_qc_checks(*self._pieces(tmp_path))  # must not raise

    def test_raises_if_paics_in_evidence_df(self, tmp_path):
        # Replace one gene (keeping row count at 28) with PAICS, so this
        # exercises the PAICS-specific guard rather than the row-count check.
        evidence_df, shortlist_df, sensitising_df, tolerance_df, paics_df, pca_df, sample_meta_df, cfg = self._pieces(tmp_path)
        evidence_df = evidence_df.copy()
        evidence_df.loc[evidence_df["gene_symbol"] == "TSC1", "gene_symbol"] = BENCHMARK_GENE_SYMBOL
        with pytest.raises(ValueError, match="must never be counted"):
            run_qc_checks(evidence_df, shortlist_df, sensitising_df, tolerance_df, paics_df, pca_df, sample_meta_df, cfg)

    def test_raises_if_primary_is_not_exactly_usp34(self, tmp_path):
        evidence_df, shortlist_df, sensitising_df, tolerance_df, paics_df, pca_df, sample_meta_df, cfg = self._pieces(tmp_path)
        shortlist_df = pd.DataFrame([{"gene_symbol": "SOME_OTHER_GENE"}])
        with pytest.raises(ValueError, match="shortlist is not exactly"):
            run_qc_checks(evidence_df, shortlist_df, sensitising_df, tolerance_df, paics_df, pca_df, sample_meta_df, cfg)

    def test_raises_if_wrong_sample_count(self, tmp_path):
        evidence_df, shortlist_df, sensitising_df, tolerance_df, paics_df, pca_df, sample_meta_df, cfg = self._pieces(tmp_path)
        pca_df = pd.DataFrame({"sample_id": [f"S{i}" for i in range(8)]})
        with pytest.raises(ValueError, match="expected 9 PCA samples"):
            run_qc_checks(evidence_df, shortlist_df, sensitising_df, tolerance_df, paics_df, pca_df, sample_meta_df, cfg)

    def test_raises_if_wrong_sensitising_count(self, tmp_path):
        evidence_df, shortlist_df, sensitising_df, tolerance_df, paics_df, pca_df, sample_meta_df, cfg = self._pieces(tmp_path)
        sensitising_df = sensitising_df.iloc[:-1]
        with pytest.raises(ValueError, match="sensitising candidates table"):
            run_qc_checks(evidence_df, shortlist_df, sensitising_df, tolerance_df, paics_df, pca_df, sample_meta_df, cfg)

    # --- fix 1: PAICS must never appear in any candidate-facing table -----

    def test_raises_if_paics_in_sensitising_df(self, tmp_path):
        evidence_df, shortlist_df, sensitising_df, tolerance_df, paics_df, pca_df, sample_meta_df, cfg = self._pieces(tmp_path)
        sensitising_df = sensitising_df.copy()
        sensitising_df.loc[sensitising_df["gene_symbol"] == "ICK", "gene_symbol"] = BENCHMARK_GENE_SYMBOL
        with pytest.raises(ValueError, match="must never appear in sensitising_df"):
            run_qc_checks(evidence_df, shortlist_df, sensitising_df, tolerance_df, paics_df, pca_df, sample_meta_df, cfg)

    def test_raises_if_paics_in_tolerance_df(self, tmp_path):
        evidence_df, shortlist_df, sensitising_df, tolerance_df, paics_df, pca_df, sample_meta_df, cfg = self._pieces(tmp_path)
        tolerance_df = tolerance_df.copy()
        tolerance_df.loc[tolerance_df["gene_symbol"] == "PTEN", "gene_symbol"] = BENCHMARK_GENE_SYMBOL
        with pytest.raises(ValueError, match="must never appear in tolerance_df"):
            run_qc_checks(evidence_df, shortlist_df, sensitising_df, tolerance_df, paics_df, pca_df, sample_meta_df, cfg)

    def test_raises_if_paics_in_shortlist_df(self, tmp_path):
        evidence_df, shortlist_df, sensitising_df, tolerance_df, paics_df, pca_df, sample_meta_df, cfg = self._pieces(tmp_path)
        shortlist_df = pd.DataFrame([{"gene_symbol": BENCHMARK_GENE_SYMBOL}])
        with pytest.raises(ValueError, match="must never appear in shortlist_df"):
            run_qc_checks(evidence_df, shortlist_df, sensitising_df, tolerance_df, paics_df, pca_df, sample_meta_df, cfg)

    # --- fix 2: crispr_direction must agree with crispr_effect_size sign --

    def test_raises_if_evidence_df_direction_contradicts_sign(self, tmp_path):
        evidence_df, shortlist_df, sensitising_df, tolerance_df, paics_df, pca_df, sample_meta_df, cfg = self._pieces(tmp_path)
        evidence_df = evidence_df.copy()
        evidence_df.loc[evidence_df["gene_symbol"] == "USP34", "crispr_direction"] = DIRECTION_TOLERANCE
        with pytest.raises(ValueError, match="evidence_df: crispr_direction does not match"):
            run_qc_checks(evidence_df, shortlist_df, sensitising_df, tolerance_df, paics_df, pca_df, sample_meta_df, cfg)

    def test_raises_if_sensitising_df_direction_contradicts_sign(self, tmp_path):
        evidence_df, shortlist_df, sensitising_df, tolerance_df, paics_df, pca_df, sample_meta_df, cfg = self._pieces(tmp_path)
        sensitising_df = sensitising_df.copy()
        sensitising_df.loc[sensitising_df["gene_symbol"] == "USP34", "crispr_direction"] = DIRECTION_TOLERANCE
        with pytest.raises(ValueError, match="sensitising_df: crispr_direction does not match"):
            run_qc_checks(evidence_df, shortlist_df, sensitising_df, tolerance_df, paics_df, pca_df, sample_meta_df, cfg)

    # --- fix 5: PCA group balance and PCA/metadata sample-set agreement ---

    def test_raises_if_pca_group_imbalanced(self, tmp_path):
        evidence_df, shortlist_df, sensitising_df, tolerance_df, paics_df, pca_df, sample_meta_df, cfg = self._pieces(tmp_path)
        pca_df = pca_df.copy()
        pca_df.loc[pca_df["sample_id"] == "FASR_Rep3", "group"] = "TAMR"
        with pytest.raises(ValueError, match="expected 3 PCA samples for group"):
            run_qc_checks(evidence_df, shortlist_df, sensitising_df, tolerance_df, paics_df, pca_df, sample_meta_df, cfg)

    def test_raises_if_pca_and_metadata_sample_sets_differ(self, tmp_path):
        evidence_df, shortlist_df, sensitising_df, tolerance_df, paics_df, pca_df, sample_meta_df, cfg = self._pieces(tmp_path)
        sample_meta_df = sample_meta_df.copy()
        sample_meta_df.loc[sample_meta_df["sample_id"] == "FASR_Rep3", "sample_id"] = "FASR_Rep4"
        with pytest.raises(ValueError, match="PCA sample set does not match sample metadata"):
            run_qc_checks(evidence_df, shortlist_df, sensitising_df, tolerance_df, paics_df, pca_df, sample_meta_df, cfg)


# --- fix 1 (round 3): crispr_effect_size must be finite and non-zero ------


class TestValidateDirectionMatchesSignEdgeCases:
    def _df(self, effect_size, direction):
        return pd.DataFrame([{"gene_symbol": "G1", "crispr_effect_size": effect_size, "crispr_direction": direction}])

    def test_raises_on_nan_effect_size(self):
        with pytest.raises(ValueError, match="crispr_effect_size is not finite"):
            _validate_direction_matches_sign(self._df(float("nan"), DIRECTION_SENSITISING), "ctx")

    def test_raises_on_positive_infinity(self):
        with pytest.raises(ValueError, match="crispr_effect_size is not finite"):
            _validate_direction_matches_sign(self._df(float("inf"), DIRECTION_TOLERANCE), "ctx")

    def test_raises_on_negative_infinity(self):
        with pytest.raises(ValueError, match="crispr_effect_size is not finite"):
            _validate_direction_matches_sign(self._df(float("-inf"), DIRECTION_SENSITISING), "ctx")

    def test_raises_on_exact_zero(self):
        with pytest.raises(ValueError, match="crispr_effect_size is exactly zero"):
            _validate_direction_matches_sign(self._df(0.0, DIRECTION_SENSITISING), "ctx")

    def test_accepts_correctly_signed_finite_nonzero_values(self):
        _validate_direction_matches_sign(self._df(-1.4, DIRECTION_SENSITISING), "ctx")  # must not raise
        _validate_direction_matches_sign(self._df(1.4, DIRECTION_TOLERANCE), "ctx")  # must not raise


# --- fix 2 (round 3): PCA group label must agree between PC1 and PC2 ------


class TestLoadPcaCoordinatesGroupConsistency:
    def _rows(self, n=9):
        rows = []
        groups = (["MCF7"] * 3 + ["TAMR"] * 3 + ["FASR"] * 3)[:n]
        for i in range(n):
            for pc, coord in (("PC1", float(i)), ("PC2", float(-i))):
                rows.append(
                    {"sample_id": f"S{i}", "group": groups[i], "pc": pc, "coordinate": coord, "variance_explained_fraction": 0.5}
                )
        return rows

    def test_raises_if_pc1_and_pc2_group_disagree_for_a_sample(self, tmp_path):
        cfg = _cfg(tmp_path, expected_n_samples=9)
        rows = self._rows()
        # Corrupt S0's PC2 row to report a different group than its PC1 row.
        for r in rows:
            if r["sample_id"] == "S0" and r["pc"] == "PC2":
                r["group"] = "TAMR"
        pd.DataFrame(rows).to_csv(cfg.pca_coordinates_tsv, sep="\t", index=False)
        with pytest.raises(ValueError, match="PCA group label mismatch"):
            load_pca_coordinates(cfg)

    def test_passes_when_pc1_and_pc2_groups_agree(self, tmp_path):
        cfg = _cfg(tmp_path, expected_n_samples=9)
        pd.DataFrame(self._rows()).to_csv(cfg.pca_coordinates_tsv, sep="\t", index=False)
        df = load_pca_coordinates(cfg)  # must not raise
        assert df["sample_id"].nunique() == 9


# --- file-based loaders: structural validation ------------------------------


class TestLoaders:
    def test_load_evidence_summary_raises_on_wrong_row_count(self, tmp_path):
        cfg = _cfg(tmp_path)
        pd.DataFrame([_evidence_row("A", -1.0, EVIDENCE_CLASS_PRIMARY)]).to_csv(cfg.evidence_summary_tsv, sep="\t", index=False)
        with pytest.raises(ValueError, match="expected 28"):
            load_evidence_summary(cfg)

    def test_load_evidence_summary_raises_if_paics_present(self, tmp_path):
        cfg = _cfg(tmp_path, expected_n_hits=2)
        df = pd.DataFrame(
            [_evidence_row("A", -1.0, EVIDENCE_CLASS_PRIMARY), _evidence_row(BENCHMARK_GENE_SYMBOL, -1.0, EVIDENCE_CLASS_PRIMARY)]
        )
        df.to_csv(cfg.evidence_summary_tsv, sep="\t", index=False)
        with pytest.raises(ValueError, match="must never appear"):
            load_evidence_summary(cfg)

    def test_load_shortlist_raises_if_not_exactly_primary_gene(self, tmp_path):
        cfg = _cfg(tmp_path)
        pd.DataFrame([{"gene_symbol": "WRONG_GENE"}]).to_csv(cfg.shortlist_tsv, sep="\t", index=False)
        with pytest.raises(ValueError, match="primary shortlist"):
            load_shortlist(cfg)

    def test_load_paics_benchmark_raises_if_mislabelled(self, tmp_path):
        cfg = _cfg(tmp_path)
        pd.DataFrame([{"gene_symbol": BENCHMARK_GENE_SYMBOL, "benchmark_label": "wrong_label"}]).to_csv(
            cfg.paics_benchmark_tsv, sep="\t", index=False
        )
        with pytest.raises(ValueError, match="missing expected label"):
            load_paics_benchmark(cfg)

    def test_load_pca_coordinates_raises_on_wrong_sample_count(self, tmp_path):
        cfg = _cfg(tmp_path)
        pd.DataFrame({"sample_id": ["A", "B"], "group": ["MCF7", "MCF7"], "pc": ["PC1", "PC1"], "coordinate": [1.0, 2.0], "variance_explained_fraction": [0.5, 0.5]}).to_csv(
            cfg.pca_coordinates_tsv, sep="\t", index=False
        )
        with pytest.raises(ValueError, match="expected 9 PCA samples"):
            load_pca_coordinates(cfg)

    def _valid_pca_rows(self, n=9):
        rows = []
        for i in range(n):
            for pc, coord in (("PC1", float(i)), ("PC2", float(-i))):
                rows.append(
                    {"sample_id": f"S{i}", "group": "MCF7", "pc": pc, "coordinate": coord, "variance_explained_fraction": 0.5}
                )
        return rows

    def test_load_pca_coordinates_raises_on_duplicated_pc_row(self, tmp_path):
        cfg = _cfg(tmp_path, expected_n_samples=9)
        rows = self._valid_pca_rows()
        rows.append({"sample_id": "S0", "group": "MCF7", "pc": "PC1", "coordinate": 99.0, "variance_explained_fraction": 0.5})
        pd.DataFrame(rows).to_csv(cfg.pca_coordinates_tsv, sep="\t", index=False)
        with pytest.raises(ValueError, match="expected exactly one PC1 coordinate per sample"):
            load_pca_coordinates(cfg)

    def test_load_pca_coordinates_raises_on_missing_pc_row(self, tmp_path):
        cfg = _cfg(tmp_path, expected_n_samples=9)
        rows = [r for r in self._valid_pca_rows() if not (r["sample_id"] == "S0" and r["pc"] == "PC2")]
        pd.DataFrame(rows).to_csv(cfg.pca_coordinates_tsv, sep="\t", index=False)
        with pytest.raises(ValueError, match="expected exactly one PC2 coordinate per sample"):
            load_pca_coordinates(cfg)

    def test_load_pca_coordinates_raises_on_non_finite_coordinate(self, tmp_path):
        cfg = _cfg(tmp_path, expected_n_samples=9)
        rows = self._valid_pca_rows()
        rows[0]["coordinate"] = float("nan")
        pd.DataFrame(rows).to_csv(cfg.pca_coordinates_tsv, sep="\t", index=False)
        with pytest.raises(ValueError, match="non-finite values"):
            load_pca_coordinates(cfg)

    def test_load_sample_metadata_raises_on_wrong_sample_count(self, tmp_path):
        cfg = _cfg(tmp_path)
        pd.DataFrame({"sample_id": ["A"], "group": ["MCF7"], "replicate": [1]}).to_csv(cfg.sample_metadata_tsv, sep="\t", index=False)
        with pytest.raises(ValueError, match="expected 9 samples"):
            load_sample_metadata(cfg)


# --- real-repository smoke test: full pipeline, isolated output paths ------


class TestRunAgainstRealConfig:
    def test_real_pipeline_produces_expected_figures_and_qc(self, tmp_path):
        with open(CONFIG_PATH) as f:
            config = yaml.safe_load(f)
        config = copy.deepcopy(config)
        config["nebula_plots"]["output_dir"] = str(tmp_path / "figs")
        config["nebula_plots"]["plot_input_dir"] = str(tmp_path / "plot_inputs")
        figs = config["nebula_plots"]["figures"]
        for fig_key, paths in figs.items():
            for key, original in paths.items():
                paths[key] = str(tmp_path / f"{fig_key}_{key}{Path(original).suffix}")

        tmp_config_path = tmp_path / "config.yaml"
        with open(tmp_config_path, "w") as f:
            yaml.safe_dump(config, f)

        result = run_nebula_plots(str(tmp_config_path))

        assert len(result["evidence_summary"]) == 28
        assert "PAICS" not in set(result["evidence_summary"]["gene_symbol"])
        assert list(result["shortlist"]["gene_symbol"]) == ["USP34"]
        assert len(result["sensitising"]) == 13
        assert len(result["tolerance"]) == 15
        assert len(result["fig4_input"]) == 13
        assert result["fig2_input"]["sample_id"].nunique() == 9
        assert result["paics_benchmark"].iloc[0]["benchmark_label"] == BENCHMARK_LABEL

        for fig_key, paths in figs.items():
            for key, path_str in paths.items():
                path = Path(path_str)
                assert path.exists(), f"{fig_key}.{key} was not written"
                assert path.stat().st_size > 0
