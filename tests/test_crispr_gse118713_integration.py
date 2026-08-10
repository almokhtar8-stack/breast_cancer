from pathlib import Path

import pandas as pd
import pytest

from src.crispr_gse118713_integration import (
    DE_NA_BLINDED,
    DE_NA_FILTERED_OUT,
    DE_NA_NOT_MAPPED,
    MAPPING_AMBIGUOUS,
    MAPPING_UNIQUE_FILTERED,
    MAPPING_UNIQUE_FILTERED_OUT,
    MAPPING_UNMATCHED,
    MASTER_TABLE_COLUMNS,
    IntegrationConfig,
    build_master_table,
    build_symbol_index,
    compute_qc_summary,
    load_crispr_hits,
    run_integration,
)


def _annotation_row(gene_id, gene_symbol, status, mcf7=(1.0, 1.0, 1.0), tamr=(1.0, 1.0, 1.0), fasr=(1.0, 1.0, 1.0)):
    return {
        "gene_id": gene_id,
        "gene_symbol": gene_symbol,
        "symbol_mapping_status": status,
        "MCF7_Rep1": mcf7[0],
        "MCF7_Rep2": mcf7[1],
        "MCF7_Rep3": mcf7[2],
        "TAMR_Rep1": tamr[0],
        "TAMR_Rep2": tamr[1],
        "TAMR_Rep3": tamr[2],
        "FASR_Rep1": fasr[0],
        "FASR_Rep2": fasr[1],
        "FASR_Rep3": fasr[2],
    }


def _hit(gene, effect_size=1.0, se=0.2, p_value=0.001, fdr=0.01, n_guides=4):
    return {"gene": gene, "n_guides": n_guides, "effect_size": effect_size, "se": se, "p_value": p_value, "fdr": fdr}


def _specificity_row(gene_id, mcf7_fc=0.5, mcf7_fdr=0.2, fasr_fc=0.3, fasr_fdr=0.4):
    return {
        "gene_id": gene_id,
        "tamr_vs_mcf7_log2fc": mcf7_fc,
        "tamr_vs_mcf7_fdr": mcf7_fdr,
        "tamr_vs_fasr_log2fc": fasr_fc,
        "tamr_vs_fasr_fdr": fasr_fdr,
    }


# --- build_symbol_index: no silent alias substitution --------------------


class TestBuildSymbolIndex:
    def test_indexes_only_resolved_rows(self):
        annotation = pd.DataFrame(
            [
                _annotation_row("G1", "SYM1", "resolved"),
                _annotation_row("G2", "SYM2", "missing"),
                _annotation_row("G3", "SYM3", "ambiguous"),
            ]
        )
        index = build_symbol_index(annotation)
        assert index == {"SYM1": ["G1"]}

    def test_symbol_shared_by_two_gene_ids_is_kept_as_two_candidates(self):
        annotation = pd.DataFrame(
            [
                _annotation_row("G1", "SHARED", "resolved"),
                _annotation_row("G2", "SHARED", "resolved"),
            ]
        )
        index = build_symbol_index(annotation)
        assert sorted(index["SHARED"]) == ["G1", "G2"]


# --- build_master_table: mapping logic and row preservation --------------


class TestBuildMasterTable:
    def test_unique_and_filter_passing_gene_gets_full_de_row(self):
        hits = pd.DataFrame([_hit("SYM1")])
        annotation = pd.DataFrame([_annotation_row("G1", "SYM1", "resolved", mcf7=(1.0, 3.0, 7.0))])
        filtered_ids = {"G1"}
        specificity = pd.DataFrame([_specificity_row("G1", mcf7_fc=0.7, mcf7_fdr=0.01, fasr_fc=-0.2, fasr_fdr=0.3)])

        out = build_master_table(hits, annotation, filtered_ids, specificity)
        assert len(out) == 1
        row = out.iloc[0]
        assert row["mapping_status"] == MAPPING_UNIQUE_FILTERED
        assert row["gse118713_gene_id"] == "G1"
        assert bool(row["passed_gse118713_expression_filter"]) is True
        assert row["tamr_vs_mcf7_log2fc"] == pytest.approx(0.7)
        assert row["tamr_vs_fasr_fdr"] == pytest.approx(0.3)
        assert pd.isna(row["gse118713_de_na_reason"])
        # mean(log2(x+1)) over TPM 1, 3, 7
        import numpy as np

        expected = float(np.mean(np.log2(np.array([1.0, 3.0, 7.0]) + 1.0)))
        assert row["mcf7_baseline_log2_tpm_plus1"] == pytest.approx(expected)

    def test_ambiguous_symbol_excluded_from_rna_fields_but_row_kept(self):
        hits = pd.DataFrame([_hit("SHARED")])
        annotation = pd.DataFrame(
            [
                _annotation_row("G1", "SHARED", "resolved"),
                _annotation_row("G2", "SHARED", "resolved"),
            ]
        )
        out = build_master_table(hits, annotation, filtered_ids={"G1", "G2"}, specificity_df=pd.DataFrame(columns=["gene_id", "tamr_vs_mcf7_log2fc", "tamr_vs_mcf7_fdr", "tamr_vs_fasr_log2fc", "tamr_vs_fasr_fdr"]))
        assert len(out) == 1
        row = out.iloc[0]
        assert row["mapping_status"] == MAPPING_AMBIGUOUS
        assert pd.isna(row["gse118713_gene_id"])
        assert pd.isna(row["tamr_vs_mcf7_log2fc"])
        assert row["gse118713_de_na_reason"] == DE_NA_NOT_MAPPED

    def test_unmatched_symbol_row_kept_with_na_rna_fields(self):
        hits = pd.DataFrame([_hit("NOWHERE")])
        annotation = pd.DataFrame([_annotation_row("G1", "OTHER", "resolved")])
        out = build_master_table(hits, annotation, filtered_ids={"G1"}, specificity_df=pd.DataFrame(columns=["gene_id", "tamr_vs_mcf7_log2fc", "tamr_vs_mcf7_fdr", "tamr_vs_fasr_log2fc", "tamr_vs_fasr_fdr"]))
        assert len(out) == 1
        row = out.iloc[0]
        assert row["mapping_status"] == MAPPING_UNMATCHED
        assert pd.isna(row["gse118713_gene_id"])
        assert row["gse118713_de_na_reason"] == DE_NA_NOT_MAPPED

    def test_unique_but_filtered_out_gene_has_baseline_but_no_de(self):
        hits = pd.DataFrame([_hit("SYM1")])
        annotation = pd.DataFrame([_annotation_row("G1", "SYM1", "resolved", mcf7=(0.0, 0.0, 0.0))])
        out = build_master_table(
            hits,
            annotation,
            filtered_ids=set(),  # G1 did not pass the expression filter
            specificity_df=pd.DataFrame(columns=["gene_id", "tamr_vs_mcf7_log2fc", "tamr_vs_mcf7_fdr", "tamr_vs_fasr_log2fc", "tamr_vs_fasr_fdr"]),
        )
        row = out.iloc[0]
        assert row["mapping_status"] == MAPPING_UNIQUE_FILTERED_OUT
        assert bool(row["passed_gse118713_expression_filter"]) is False
        assert row["mcf7_baseline_log2_tpm_plus1"] == pytest.approx(0.0)
        assert pd.isna(row["tamr_vs_mcf7_log2fc"])
        assert row["gse118713_de_na_reason"] == DE_NA_FILTERED_OUT

    def test_mapped_and_filter_passing_but_absent_from_specificity_is_blinded_gap(self):
        # Mirrors RCOR1/KDM1A: uniquely mapped, passed the expression
        # filter, but absent from the specificity table because the limma
        # pipeline redacted them before writing any DE output.
        hits = pd.DataFrame([_hit("SYM1")])
        annotation = pd.DataFrame([_annotation_row("G1", "SYM1", "resolved")])
        out = build_master_table(
            hits,
            annotation,
            filtered_ids={"G1"},
            specificity_df=pd.DataFrame(columns=["gene_id", "tamr_vs_mcf7_log2fc", "tamr_vs_mcf7_fdr", "tamr_vs_fasr_log2fc", "tamr_vs_fasr_fdr"]),
        )
        row = out.iloc[0]
        assert row["mapping_status"] == MAPPING_UNIQUE_FILTERED
        assert pd.isna(row["tamr_vs_mcf7_log2fc"])
        assert row["gse118713_de_na_reason"] == DE_NA_BLINDED

    def test_no_hit_is_ever_dropped(self):
        hits = pd.DataFrame([_hit("MATCHED"), _hit("AMBIG"), _hit("MISSING")])
        annotation = pd.DataFrame(
            [
                _annotation_row("G1", "MATCHED", "resolved"),
                _annotation_row("G2", "AMBIG", "resolved"),
                _annotation_row("G3", "AMBIG", "resolved"),
            ]
        )
        out = build_master_table(
            hits,
            annotation,
            filtered_ids={"G1", "G2", "G3"},
            specificity_df=pd.DataFrame(columns=["gene_id", "tamr_vs_mcf7_log2fc", "tamr_vs_mcf7_fdr", "tamr_vs_fasr_log2fc", "tamr_vs_fasr_fdr"]),
        )
        assert len(out) == 3
        assert set(out["gene_symbol"]) == {"MATCHED", "AMBIG", "MISSING"}

    def test_required_columns_present(self):
        hits = pd.DataFrame([_hit("SYM1")])
        annotation = pd.DataFrame([_annotation_row("G1", "SYM1", "resolved")])
        out = build_master_table(
            hits,
            annotation,
            filtered_ids={"G1"},
            specificity_df=pd.DataFrame(columns=["gene_id", "tamr_vs_mcf7_log2fc", "tamr_vs_mcf7_fdr", "tamr_vs_fasr_log2fc", "tamr_vs_fasr_fdr"]),
        )
        assert list(out.columns) == list(MASTER_TABLE_COLUMNS)

    def test_gene_symbol_and_mapped_gene_id_are_unique(self):
        hits = pd.DataFrame([_hit("SYM1"), _hit("SYM2")])
        annotation = pd.DataFrame(
            [
                _annotation_row("G1", "SYM1", "resolved"),
                _annotation_row("G2", "SYM2", "resolved"),
            ]
        )
        out = build_master_table(
            hits,
            annotation,
            filtered_ids={"G1", "G2"},
            specificity_df=pd.DataFrame(columns=["gene_id", "tamr_vs_mcf7_log2fc", "tamr_vs_mcf7_fdr", "tamr_vs_fasr_log2fc", "tamr_vs_fasr_fdr"]),
        )
        assert out["gene_symbol"].is_unique
        mapped_ids = out["gse118713_gene_id"].dropna()
        assert mapped_ids.is_unique


# --- compute_qc_summary: aggregate counts only ----------------------------


class TestComputeQcSummary:
    def test_counts_match_construction(self):
        hits = pd.DataFrame([_hit(g) for g in ("A", "B", "C", "D")])
        annotation = pd.DataFrame(
            [
                _annotation_row("G1", "A", "resolved"),
                _annotation_row("G2", "B", "resolved"),
                _annotation_row("G3", "C", "resolved"),
                _annotation_row("G4", "D", "resolved"),
            ]
        )
        specificity = pd.DataFrame(
            [
                _specificity_row("G1", mcf7_fc=1.0, mcf7_fdr=0.01, fasr_fc=1.0, fasr_fdr=0.01),  # sig both
                _specificity_row("G2", mcf7_fc=-1.0, mcf7_fdr=0.2, fasr_fc=-1.0, fasr_fdr=0.02),  # sig fasr only
            ]
        )
        # G3 filtered out; G4 unmatched intentionally omitted from annotation lookup by symbol change
        out = build_master_table(hits, annotation, filtered_ids={"G1", "G2"}, specificity_df=specificity)
        summary = compute_qc_summary(out)
        assert summary["n_crispr_hits_total"] == 4
        assert summary["n_mapped_unique_and_filtered"] == 2
        assert summary["n_mapped_unique_but_filtered_out"] == 2
        assert summary["n_mapping_ambiguous"] == 0
        assert summary["n_mapping_unmatched"] == 0
        assert summary["n_tamr_vs_mcf7_fdr_lt_0_05"] == 1
        assert summary["n_tamr_vs_fasr_fdr_lt_0_05"] == 2
        assert summary["n_significant_in_both_comparisons"] == 1

    def test_handles_all_na_gracefully(self):
        hits = pd.DataFrame([_hit("A")])
        annotation = pd.DataFrame([_annotation_row("G1", "A", "resolved")])
        out = build_master_table(
            hits,
            annotation,
            filtered_ids=set(),
            specificity_df=pd.DataFrame(columns=["gene_id", "tamr_vs_mcf7_log2fc", "tamr_vs_mcf7_fdr", "tamr_vs_fasr_log2fc", "tamr_vs_fasr_fdr"]),
        )
        summary = compute_qc_summary(out)
        assert summary["n_tamr_vs_mcf7_fdr_lt_0_05"] == 0
        assert summary["n_significant_in_both_comparisons"] == 0


# --- load_crispr_hits: cross-check against the committed Gate-1 decision --


class TestLoadCrisprHits:
    def test_mismatched_recorded_decision_raises(self, tmp_path: Path):
        labels_path = tmp_path / "labels.parquet"
        pd.DataFrame(
            {
                "gene": ["A", "B", "C"],
                "n_guides": [4, 4, 4],
                "effect_size": [1.0, -1.0, 0.5],
                "se": [0.1, 0.1, 0.1],
                "p_value": [0.001, 0.001, 0.5],
                "fdr": [0.01, 0.01, 0.5],
            }
        ).to_parquet(labels_path, index=False)
        decision_path = tmp_path / "gate1_decision.tsv"
        pd.DataFrame([{"n_passing": 99}]).to_csv(decision_path, sep="\t", index=False)

        cfg = IntegrationConfig(
            labels_path=labels_path,
            fdr_threshold=0.1,
            n_fitted_genes_expected=3,
            gate1_decision_tsv=decision_path,
            expected_n_hits=99,
            gene_tpm_parquet=tmp_path / "unused.parquet",
            frozen_gene_tpm_sha256="unused",
            filtered_gene_tpm_tsv=tmp_path / "unused.tsv.gz",
            specificity_tsv_gz=tmp_path / "unused.tsv.gz",
            master_table_tsv=tmp_path / "out.tsv",
            master_table_csv=tmp_path / "out.csv",
            master_table_parquet=tmp_path / "out.parquet",
            qc_summary_tsv=tmp_path / "qc.tsv",
            qc_summary_md=tmp_path / "qc.md",
        )
        with pytest.raises(ValueError, match="does not match"):
            load_crispr_hits(cfg)


# --- real-repository smoke test: exact deliverable shape ------------------


class TestRunIntegrationAgainstRealConfig:
    def test_real_pipeline_produces_exactly_28_rows_with_required_columns(self):
        df = run_integration("config/config.yaml")
        assert len(df) == 28
        assert list(df.columns) == list(MASTER_TABLE_COLUMNS)
        assert df["gene_symbol"].is_unique
        mapped_ids = df["gse118713_gene_id"].dropna()
        assert mapped_ids.is_unique
        assert df["mapping_status"].isin(
            [MAPPING_UNIQUE_FILTERED, MAPPING_UNIQUE_FILTERED_OUT, MAPPING_AMBIGUOUS, MAPPING_UNMATCHED]
        ).all()
