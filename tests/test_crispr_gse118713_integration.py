import copy
import hashlib
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import yaml

from src.gate1_checks import decide_gate1
from src.crispr_gse118713_integration import (
    DE_NA_FILTERED_OUT,
    DE_NA_HISTORICALLY_BLINDED,
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
    load_gse118713_annotation,
    load_gse118713_specificity,
    run_integration,
)

REPO_ROOT = Path(__file__).parent.parent
CONFIG_PATH = REPO_ROOT / "config" / "config.yaml"
EMPTY_SPECIFICITY = pd.DataFrame(
    columns=["gene_id", "tamr_vs_mcf7_log2fc", "tamr_vs_mcf7_fdr", "tamr_vs_fasr_log2fc", "tamr_vs_fasr_fdr"]
)


def _sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _sha256_file(path: Path) -> str:
    return _sha256_bytes(path.read_bytes())


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


def _make_labels(genes, fdrs, effect_sizes=None):
    n = len(genes)
    effect_sizes = effect_sizes if effect_sizes is not None else [1.0] * n
    return pd.DataFrame(
        {
            "gene": genes,
            "n_guides": [4] * n,
            "effect_size": effect_sizes,
            "se": [0.2] * n,
            "p_value": [0.01] * n,
            "fdr": fdrs,
        }
    )


def _integration_cfg(tmp_path: Path, **overrides) -> IntegrationConfig:
    defaults = dict(
        labels_path=tmp_path / "unused_labels.parquet",
        frozen_labels_sha256="unused",
        fdr_threshold=0.1,
        n_fitted_genes_expected=3,
        classifier_min=30,
        continuous_min=10,
        gate1_decision_tsv=tmp_path / "unused_decision.tsv",
        expected_n_hits=99,
        gene_tpm_parquet=tmp_path / "unused.parquet",
        frozen_gene_tpm_sha256="unused",
        filtered_gene_tpm_tsv=tmp_path / "unused.tsv.gz",
        frozen_filtered_gene_tpm_sha256="unused",
        specificity_tsv_gz=tmp_path / "unused.tsv.gz",
        frozen_specificity_sha256="unused",
        blinded_gene_ids=(),
        master_table_tsv=tmp_path / "out.tsv",
        master_table_csv=tmp_path / "out.csv",
        master_table_parquet=tmp_path / "out.parquet",
        qc_summary_tsv=tmp_path / "qc.tsv",
        qc_summary_md=tmp_path / "qc.md",
    )
    defaults.update(overrides)
    return IntegrationConfig(**defaults)


# --- build_symbol_index: no silent alias substitution, dedup safety ------


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

    def test_symbol_shared_by_two_distinct_gene_ids_is_kept_as_two_candidates(self):
        annotation = pd.DataFrame(
            [
                _annotation_row("G1", "SHARED", "resolved"),
                _annotation_row("G2", "SHARED", "resolved"),
            ]
        )
        index = build_symbol_index(annotation)
        assert sorted(index["SHARED"]) == ["G1", "G2"]

    def test_duplicate_identical_gene_id_symbol_rows_do_not_create_false_ambiguity(self):
        # Same (gene_id, symbol) pair appears twice -- e.g. a duplicated
        # annotation row -- must collapse to ONE candidate, not two.
        annotation = pd.DataFrame(
            [
                _annotation_row("G1", "SYM1", "resolved"),
                _annotation_row("G1", "SYM1", "resolved"),
            ]
        )
        index = build_symbol_index(annotation)
        assert index == {"SYM1": ["G1"]}


# --- build_master_table: mapping logic, row preservation, blind handling -


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
        out = build_master_table(hits, annotation, filtered_ids={"G1", "G2"}, specificity_df=EMPTY_SPECIFICITY)
        assert len(out) == 1
        row = out.iloc[0]
        assert row["mapping_status"] == MAPPING_AMBIGUOUS
        assert pd.isna(row["gse118713_gene_id"])
        assert pd.isna(row["tamr_vs_mcf7_log2fc"])
        assert row["gse118713_de_na_reason"] == DE_NA_NOT_MAPPED

    def test_unmatched_symbol_row_kept_with_na_rna_fields(self):
        hits = pd.DataFrame([_hit("NOWHERE")])
        annotation = pd.DataFrame([_annotation_row("G1", "OTHER", "resolved")])
        out = build_master_table(hits, annotation, filtered_ids={"G1"}, specificity_df=EMPTY_SPECIFICITY)
        assert len(out) == 1
        row = out.iloc[0]
        assert row["mapping_status"] == MAPPING_UNMATCHED
        assert pd.isna(row["gse118713_gene_id"])
        assert row["gse118713_de_na_reason"] == DE_NA_NOT_MAPPED

    def test_unique_but_filtered_out_gene_has_baseline_but_no_de(self):
        hits = pd.DataFrame([_hit("SYM1")])
        annotation = pd.DataFrame([_annotation_row("G1", "SYM1", "resolved", mcf7=(0.0, 0.0, 0.0))])
        out = build_master_table(hits, annotation, filtered_ids=set(), specificity_df=EMPTY_SPECIFICITY)
        row = out.iloc[0]
        assert row["mapping_status"] == MAPPING_UNIQUE_FILTERED_OUT
        assert bool(row["passed_gse118713_expression_filter"]) is False
        assert row["mcf7_baseline_log2_tpm_plus1"] == pytest.approx(0.0)
        assert pd.isna(row["tamr_vs_mcf7_log2fc"])
        assert row["gse118713_de_na_reason"] == DE_NA_FILTERED_OUT

    def test_unexplained_missing_specificity_row_raises(self):
        # Uniquely mapped, passed the filter, but missing from the
        # specificity table AND not a configured historical blind ID --
        # this must raise, never be silently assumed to be blinding.
        hits = pd.DataFrame([_hit("SYM1")])
        annotation = pd.DataFrame([_annotation_row("G1", "SYM1", "resolved")])
        with pytest.raises(ValueError, match="unexplained gap"):
            build_master_table(
                hits, annotation, filtered_ids={"G1"}, specificity_df=EMPTY_SPECIFICITY, blinded_gene_ids=()
            )

    def test_only_configured_blind_id_gets_historical_reason(self):
        hits = pd.DataFrame([_hit("SYM1"), _hit("SYM2")])
        annotation = pd.DataFrame(
            [
                _annotation_row("G1", "SYM1", "resolved"),  # configured blind id
                _annotation_row("G2", "SYM2", "resolved"),  # NOT a configured blind id
            ]
        )
        # Neither gene is in the specificity table.
        with pytest.raises(ValueError, match="unexplained gap"):
            build_master_table(
                hits, annotation, filtered_ids={"G1", "G2"}, specificity_df=EMPTY_SPECIFICITY, blinded_gene_ids=("G1",)
            )

    def test_configured_blind_id_alone_gets_historical_reason_not_an_error(self):
        hits = pd.DataFrame([_hit("SYM1")])
        annotation = pd.DataFrame([_annotation_row("G1", "SYM1", "resolved")])
        out = build_master_table(
            hits, annotation, filtered_ids={"G1"}, specificity_df=EMPTY_SPECIFICITY, blinded_gene_ids=("G1",)
        )
        row = out.iloc[0]
        assert row["mapping_status"] == MAPPING_UNIQUE_FILTERED
        assert pd.isna(row["tamr_vs_mcf7_log2fc"])
        assert row["gse118713_de_na_reason"] == DE_NA_HISTORICALLY_BLINDED

    def test_no_hit_is_ever_dropped(self):
        hits = pd.DataFrame([_hit("MATCHED"), _hit("AMBIG"), _hit("MISSING")])
        annotation = pd.DataFrame(
            [
                _annotation_row("G1", "MATCHED", "resolved"),
                _annotation_row("G2", "AMBIG", "resolved"),
                _annotation_row("G3", "AMBIG", "resolved"),
            ]
        )
        specificity = pd.DataFrame([_specificity_row("G1")])  # only MATCHED (G1) needs a DE row
        out = build_master_table(hits, annotation, filtered_ids={"G1", "G2", "G3"}, specificity_df=specificity)
        assert len(out) == 3
        assert set(out["gene_symbol"]) == {"MATCHED", "AMBIG", "MISSING"}

    def test_required_columns_present(self):
        hits = pd.DataFrame([_hit("SYM1")])
        annotation = pd.DataFrame([_annotation_row("G1", "SYM1", "resolved")])
        specificity = pd.DataFrame([_specificity_row("G1")])
        out = build_master_table(hits, annotation, filtered_ids={"G1"}, specificity_df=specificity)
        assert list(out.columns) == list(MASTER_TABLE_COLUMNS)

    def test_gene_symbol_and_mapped_gene_id_are_unique(self):
        hits = pd.DataFrame([_hit("SYM1"), _hit("SYM2")])
        annotation = pd.DataFrame(
            [
                _annotation_row("G1", "SYM1", "resolved"),
                _annotation_row("G2", "SYM2", "resolved"),
            ]
        )
        specificity = pd.DataFrame([_specificity_row("G1"), _specificity_row("G2")])
        out = build_master_table(hits, annotation, filtered_ids={"G1", "G2"}, specificity_df=specificity)
        assert out["gene_symbol"].is_unique
        assert out["gse118713_gene_id"].dropna().is_unique


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
                _specificity_row("G1", mcf7_fc=1.0, mcf7_fdr=0.01, fasr_fc=1.0, fasr_fdr=0.01),
                _specificity_row("G2", mcf7_fc=-1.0, mcf7_fdr=0.2, fasr_fc=-1.0, fasr_fdr=0.02),
            ]
        )
        out = build_master_table(hits, annotation, filtered_ids={"G1", "G2"}, specificity_df=specificity)
        summary = compute_qc_summary(out)
        assert summary["n_crispr_hits_total"] == 4
        assert summary["n_mapped_unique_and_filtered"] == 2
        assert summary["n_mapped_unique_but_filtered_out"] == 2
        assert summary["n_mapping_ambiguous"] == 0
        assert summary["n_mapping_unmatched"] == 0
        assert summary["n_de_available"] == 2
        assert summary["n_tamr_vs_mcf7_fdr_lt_0_05"] == 1
        assert summary["n_tamr_vs_fasr_fdr_lt_0_05"] == 2
        assert summary["n_significant_in_both_comparisons"] == 1

    def test_handles_all_na_gracefully(self):
        hits = pd.DataFrame([_hit("A")])
        annotation = pd.DataFrame([_annotation_row("G1", "A", "resolved")])
        out = build_master_table(hits, annotation, filtered_ids=set(), specificity_df=EMPTY_SPECIFICITY)
        summary = compute_qc_summary(out)
        assert summary["n_tamr_vs_mcf7_fdr_lt_0_05"] == 0
        assert summary["n_significant_in_both_comparisons"] == 0
        assert summary["n_de_available"] == 0


# --- load_crispr_hits: full Gate-1 decision provenance cross-check --------


class TestLoadCrisprHits:
    def _write_labels(self, path, genes, fdrs):
        labels = _make_labels(genes, fdrs)
        labels.to_parquet(path, index=False)
        return labels

    def _write_decision(self, path, decision, labels_path, labels_sha256):
        row = dict(decision)
        row["labels_input_path"] = str(labels_path)
        row["labels_file_sha256"] = labels_sha256
        pd.DataFrame([row]).to_csv(path, sep="\t", index=False)

    def test_matching_decision_record_succeeds(self, tmp_path):
        labels_path = tmp_path / "labels.parquet"
        genes = ["A", "B", "C"]
        fdrs = [0.01, 0.01, 0.5]
        self._write_labels(labels_path, genes, fdrs)
        labels_sha256 = _sha256_file(labels_path)

        decision_path = tmp_path / "decision.tsv"
        decision = decide_gate1(_make_labels(genes, fdrs), fdr_threshold=0.1, classifier_min=30, continuous_min=10)
        self._write_decision(decision_path, decision, labels_path, labels_sha256)

        cfg = _integration_cfg(
            tmp_path,
            labels_path=labels_path,
            frozen_labels_sha256=labels_sha256,
            n_fitted_genes_expected=3,
            gate1_decision_tsv=decision_path,
            expected_n_hits=2,
        )
        hits = load_crispr_hits(cfg)
        assert len(hits) == 2

    def test_fdr_threshold_mismatch_raises(self, tmp_path):
        labels_path = tmp_path / "labels.parquet"
        genes, fdrs = ["A", "B", "C"], [0.01, 0.01, 0.5]
        self._write_labels(labels_path, genes, fdrs)
        labels_sha256 = _sha256_file(labels_path)

        decision_path = tmp_path / "decision.tsv"
        decision = decide_gate1(_make_labels(genes, fdrs), fdr_threshold=0.1, classifier_min=30, continuous_min=10)
        decision["fdr_threshold"] = 0.2  # tampered
        self._write_decision(decision_path, decision, labels_path, labels_sha256)

        cfg = _integration_cfg(
            tmp_path,
            labels_path=labels_path,
            frozen_labels_sha256=labels_sha256,
            n_fitted_genes_expected=3,
            gate1_decision_tsv=decision_path,
            expected_n_hits=2,
        )
        with pytest.raises(ValueError, match="does not match"):
            load_crispr_hits(cfg)

    def test_total_fitted_genes_mismatch_raises(self, tmp_path):
        labels_path = tmp_path / "labels.parquet"
        genes, fdrs = ["A", "B", "C"], [0.01, 0.01, 0.5]
        self._write_labels(labels_path, genes, fdrs)
        labels_sha256 = _sha256_file(labels_path)

        decision_path = tmp_path / "decision.tsv"
        decision = decide_gate1(_make_labels(genes, fdrs), fdr_threshold=0.1, classifier_min=30, continuous_min=10)
        decision["total_fitted_genes"] = 999  # tampered
        self._write_decision(decision_path, decision, labels_path, labels_sha256)

        cfg = _integration_cfg(
            tmp_path,
            labels_path=labels_path,
            frozen_labels_sha256=labels_sha256,
            n_fitted_genes_expected=3,
            gate1_decision_tsv=decision_path,
            expected_n_hits=2,
        )
        with pytest.raises(ValueError, match="does not match"):
            load_crispr_hits(cfg)

    def test_labels_checksum_mismatch_raises(self, tmp_path):
        labels_path = tmp_path / "labels.parquet"
        genes, fdrs = ["A", "B", "C"], [0.01, 0.01, 0.5]
        self._write_labels(labels_path, genes, fdrs)
        real_sha256 = _sha256_file(labels_path)

        decision_path = tmp_path / "decision.tsv"
        decision = decide_gate1(_make_labels(genes, fdrs), fdr_threshold=0.1, classifier_min=30, continuous_min=10)
        self._write_decision(decision_path, decision, labels_path, real_sha256)

        cfg = _integration_cfg(
            tmp_path,
            labels_path=labels_path,
            frozen_labels_sha256="0" * 64,  # tampered / stale config checksum
            n_fitted_genes_expected=3,
            gate1_decision_tsv=decision_path,
            expected_n_hits=2,
        )
        with pytest.raises(ValueError, match="checksum mismatch"):
            load_crispr_hits(cfg)

    def test_expected_n_hits_mismatch_raises(self, tmp_path):
        labels_path = tmp_path / "labels.parquet"
        genes, fdrs = ["A", "B", "C"], [0.01, 0.01, 0.5]
        self._write_labels(labels_path, genes, fdrs)
        labels_sha256 = _sha256_file(labels_path)

        decision_path = tmp_path / "decision.tsv"
        decision = decide_gate1(_make_labels(genes, fdrs), fdr_threshold=0.1, classifier_min=30, continuous_min=10)
        self._write_decision(decision_path, decision, labels_path, labels_sha256)

        cfg = _integration_cfg(
            tmp_path,
            labels_path=labels_path,
            frozen_labels_sha256=labels_sha256,
            n_fitted_genes_expected=3,
            gate1_decision_tsv=decision_path,
            expected_n_hits=99,  # does not match the real n_passing=2
        )
        with pytest.raises(ValueError, match="expected_n_hits"):
            load_crispr_hits(cfg)


# --- load_gse118713_annotation: TPM/id validation --------------------------


class TestLoadGse118713Annotation:
    def _write(self, tmp_path, rows):
        path = tmp_path / "annotation.parquet"
        pd.DataFrame(rows).to_parquet(path, index=False)
        return path

    def test_duplicate_gene_id_raises(self, tmp_path):
        path = self._write(tmp_path, [_annotation_row("G1", "SYM1", "resolved"), _annotation_row("G1", "SYM1", "resolved")])
        cfg = _integration_cfg(tmp_path, gene_tpm_parquet=path, frozen_gene_tpm_sha256=_sha256_file(path))
        with pytest.raises(ValueError, match="duplicate gene_id"):
            load_gse118713_annotation(cfg)

    def test_negative_tpm_raises(self, tmp_path):
        path = self._write(tmp_path, [_annotation_row("G1", "SYM1", "resolved", mcf7=(-1.0, 1.0, 1.0))])
        cfg = _integration_cfg(tmp_path, gene_tpm_parquet=path, frozen_gene_tpm_sha256=_sha256_file(path))
        with pytest.raises(ValueError, match="negative TPM"):
            load_gse118713_annotation(cfg)

    def test_nonfinite_tpm_raises(self, tmp_path):
        path = self._write(tmp_path, [_annotation_row("G1", "SYM1", "resolved", mcf7=(float("nan"), 1.0, 1.0))])
        cfg = _integration_cfg(tmp_path, gene_tpm_parquet=path, frozen_gene_tpm_sha256=_sha256_file(path))
        with pytest.raises(ValueError, match="non-finite TPM"):
            load_gse118713_annotation(cfg)

    def test_checksum_mismatch_raises(self, tmp_path):
        path = self._write(tmp_path, [_annotation_row("G1", "SYM1", "resolved")])
        cfg = _integration_cfg(tmp_path, gene_tpm_parquet=path, frozen_gene_tpm_sha256="0" * 64)
        with pytest.raises(ValueError, match="checksum mismatch"):
            load_gse118713_annotation(cfg)

    def test_valid_annotation_loads(self, tmp_path):
        path = self._write(tmp_path, [_annotation_row("G1", "SYM1", "resolved")])
        cfg = _integration_cfg(tmp_path, gene_tpm_parquet=path, frozen_gene_tpm_sha256=_sha256_file(path))
        df = load_gse118713_annotation(cfg)
        assert len(df) == 1


# --- load_gse118713_specificity: DE validation -----------------------------


class TestLoadGse118713Specificity:
    def _write(self, tmp_path, rows):
        path = tmp_path / "specificity.tsv.gz"
        pd.DataFrame(rows).to_csv(path, sep="\t", index=False, compression="gzip")
        return path

    def test_duplicate_gene_id_raises(self, tmp_path):
        path = self._write(tmp_path, [_specificity_row("G1"), _specificity_row("G1")])
        cfg = _integration_cfg(tmp_path, specificity_tsv_gz=path, frozen_specificity_sha256=_sha256_file(path))
        with pytest.raises(ValueError, match="duplicate gene_id"):
            load_gse118713_specificity(cfg)

    def test_nonfinite_logfc_raises(self, tmp_path):
        path = self._write(tmp_path, [_specificity_row("G1", mcf7_fc=float("inf"))])
        cfg = _integration_cfg(tmp_path, specificity_tsv_gz=path, frozen_specificity_sha256=_sha256_file(path))
        with pytest.raises(ValueError, match="non-finite"):
            load_gse118713_specificity(cfg)

    def test_fdr_above_one_raises(self, tmp_path):
        path = self._write(tmp_path, [_specificity_row("G1", mcf7_fdr=1.5)])
        cfg = _integration_cfg(tmp_path, specificity_tsv_gz=path, frozen_specificity_sha256=_sha256_file(path))
        with pytest.raises(ValueError, match=r"\[0, 1\]"):
            load_gse118713_specificity(cfg)

    def test_fdr_below_zero_raises(self, tmp_path):
        path = self._write(tmp_path, [_specificity_row("G1", fasr_fdr=-0.1)])
        cfg = _integration_cfg(tmp_path, specificity_tsv_gz=path, frozen_specificity_sha256=_sha256_file(path))
        with pytest.raises(ValueError, match=r"\[0, 1\]"):
            load_gse118713_specificity(cfg)

    def test_checksum_mismatch_raises(self, tmp_path):
        path = self._write(tmp_path, [_specificity_row("G1")])
        cfg = _integration_cfg(tmp_path, specificity_tsv_gz=path, frozen_specificity_sha256="0" * 64)
        with pytest.raises(ValueError, match="checksum mismatch"):
            load_gse118713_specificity(cfg)

    def test_valid_specificity_loads(self, tmp_path):
        path = self._write(tmp_path, [_specificity_row("G1")])
        cfg = _integration_cfg(tmp_path, specificity_tsv_gz=path, frozen_specificity_sha256=_sha256_file(path))
        df = load_gse118713_specificity(cfg)
        assert len(df) == 1


# --- real-repository integration: exact deliverable shape, isolated I/O --


class TestRunIntegrationAgainstRealConfig:
    def _isolated_config_path(self, tmp_path: Path) -> Path:
        """Real frozen inputs, but outputs redirected to tmp_path so a
        normal test run never rewrites tracked repository files."""
        with open(CONFIG_PATH) as f:
            config = yaml.safe_load(f)
        config = copy.deepcopy(config)
        config["crispr_gse118713_integration"]["output"] = {
            "master_table_tsv": str(tmp_path / "master.tsv"),
            "master_table_csv": str(tmp_path / "master.csv"),
            "master_table_parquet": str(tmp_path / "master.parquet"),
            "qc_summary_tsv": str(tmp_path / "qc.tsv"),
            "qc_summary_md": str(tmp_path / "qc.md"),
        }
        tmp_config_path = tmp_path / "config.yaml"
        with open(tmp_config_path, "w") as f:
            yaml.safe_dump(config, f)
        return tmp_config_path

    def test_real_pipeline_produces_exactly_28_rows_with_required_columns(self, tmp_path):
        tmp_config_path = self._isolated_config_path(tmp_path)
        df = run_integration(str(tmp_config_path))

        assert len(df) == 28
        assert list(df.columns) == list(MASTER_TABLE_COLUMNS)
        assert df["gene_symbol"].is_unique
        assert df["gse118713_gene_id"].dropna().is_unique
        assert df["mapping_status"].isin(
            [MAPPING_UNIQUE_FILTERED, MAPPING_UNIQUE_FILTERED_OUT, MAPPING_AMBIGUOUS, MAPPING_UNMATCHED]
        ).all()
        # No unexplained/historical-blind gaps: the post-unblinding
        # specificity source has full coverage of every filtered gene.
        assert (df["gse118713_de_na_reason"] != DE_NA_HISTORICALLY_BLINDED).all()
        # Outputs went to the isolated tmp_path config, not tracked files.
        assert (tmp_path / "master.tsv").exists()

    def test_tsv_csv_parquet_outputs_are_content_equivalent(self, tmp_path):
        tmp_config_path = self._isolated_config_path(tmp_path)
        run_integration(str(tmp_config_path))

        tsv_df = pd.read_csv(tmp_path / "master.tsv", sep="\t")
        csv_df = pd.read_csv(tmp_path / "master.csv")
        parquet_df = pd.read_parquet(tmp_path / "master.parquet")

        # Compare as strings for the text formats (TSV/CSV round-trip NA
        # representations identically to each other) and separately check
        # parquet preserves the same shape and gene set.
        pd.testing.assert_frame_equal(tsv_df, csv_df, check_dtype=False)
        assert list(parquet_df["gene_symbol"]) == list(tsv_df["gene_symbol"])
        assert len(parquet_df) == len(tsv_df) == len(csv_df)

    def test_generated_qc_is_reproducible_from_the_master_table(self, tmp_path):
        tmp_config_path = self._isolated_config_path(tmp_path)
        df = run_integration(str(tmp_config_path))

        recomputed = compute_qc_summary(df)
        written = pd.read_csv(tmp_path / "qc.tsv", sep="\t").iloc[0].to_dict()
        for key, value in recomputed.items():
            assert written[key] == value, f"QC field {key!r} not reproducible from the master table"
