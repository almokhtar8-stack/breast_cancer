from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from src.gse118713_expression_filter import (
    FilterConfig,
    filter_expression,
    load_frozen_matrix,
    load_sample_metadata,
    validate_frozen_matrix,
    verify_checksum,
    write_filtered_matrix,
    write_filtering_summary,
)

SAMPLE_IDS = [
    "MCF7_Rep1", "MCF7_Rep2", "MCF7_Rep3",
    "TAMR_Rep1", "TAMR_Rep2", "TAMR_Rep3",
    "FASR_Rep1", "FASR_Rep2", "FASR_Rep3",
]
GROUPS = ["MCF7"] * 3 + ["TAMR"] * 3 + ["FASR"] * 3


def _meta_df() -> pd.DataFrame:
    return pd.DataFrame({"sample_id": SAMPLE_IDS, "group": GROUPS, "replicate": [1, 2, 3] * 3})


def _cfg(tmp_path: Path, gene_tpm_parquet: Path, sha256: str) -> FilterConfig:
    return FilterConfig(
        gene_tpm_parquet=gene_tpm_parquet,
        sample_metadata_tsv=tmp_path / "meta.tsv",
        expected_sha256=sha256,
        expected_n_genes=3,
        expected_n_samples=9,
        expected_groups=("MCF7", "TAMR", "FASR"),
        expected_replicates_per_group=3,
        min_tpm=1.0,
        min_samples=3,
        filtered_gene_tpm_tsv=tmp_path / "filtered.tsv.gz",
        filtering_summary_tsv=tmp_path / "filtering.tsv",
    )


def _gene_df() -> pd.DataFrame:
    data = {"gene_id": ["G1", "G2", "G3"], "gene_symbol": ["A", "B", "C"], "symbol_mapping_status": ["resolved"] * 3}
    # G1: detected (>=1) in all 9 samples -> retained
    # G2: detected in exactly 2 samples -> removed (below min_samples=3)
    # G3: detected in exactly 3 samples -> retained (boundary)
    for i, sample_id in enumerate(SAMPLE_IDS):
        data[sample_id] = [5.0, 2.0 if i < 2 else 0.0, 1.0 if i < 3 else 0.5]
    return pd.DataFrame(data)


class TestFilterExpression:
    def test_tpm_filtering_rule_keeps_ge_min_samples(self):
        df = _gene_df()
        filtered, record = filter_expression(df, SAMPLE_IDS, min_tpm=1.0, min_samples=3)
        assert set(filtered["gene_id"]) == {"G1", "G3"}
        assert record["genes_in"] == 3
        assert record["genes_retained"] == 2
        assert record["genes_removed"] == 1

    def test_no_silent_row_loss_accounting(self):
        df = _gene_df()
        _, record = filter_expression(df, SAMPLE_IDS, min_tpm=1.0, min_samples=3)
        assert record["genes_in"] == record["genes_retained"] + record["genes_removed"]

    def test_boundary_exactly_min_samples_is_retained(self):
        df = _gene_df()
        filtered, _ = filter_expression(df, SAMPLE_IDS, min_tpm=1.0, min_samples=3)
        assert "G3" in set(filtered["gene_id"])  # detected in exactly 3 samples

    def test_per_sample_detection_counts_recorded(self):
        df = _gene_df()
        _, record = filter_expression(df, SAMPLE_IDS, min_tpm=1.0, min_samples=3)
        # G1 detected (>=1) in all 9 samples; G2 in 2; G3 in 3 -> total per sample varies by column
        for sample_id in SAMPLE_IDS:
            assert f"{sample_id}_genes_detected" in record

    def test_threshold_recorded_not_hardcoded_downstream(self):
        df = _gene_df()
        _, record_strict = filter_expression(df, SAMPLE_IDS, min_tpm=2.0, min_samples=3)
        _, record_loose = filter_expression(df, SAMPLE_IDS, min_tpm=0.1, min_samples=1)
        assert record_strict["genes_retained"] != record_loose["genes_retained"]
        assert record_strict["min_tpm_threshold"] == 2.0
        assert record_loose["min_tpm_threshold"] == 0.1


class TestValidateFrozenMatrix:
    def _base_cfg(self, tmp_path):
        return _cfg(tmp_path, tmp_path / "x.parquet", sha256="unused")

    def test_rejects_duplicate_gene_id(self, tmp_path):
        cfg = self._base_cfg(tmp_path)
        df = _gene_df()
        df.loc[1, "gene_id"] = "G1"  # duplicate
        with pytest.raises(ValueError, match="duplicate gene_id"):
            validate_frozen_matrix(df, cfg, SAMPLE_IDS)

    def test_rejects_wrong_gene_count(self, tmp_path):
        cfg = _cfg(tmp_path, tmp_path / "x.parquet", sha256="unused")
        object.__setattr__(cfg, "expected_n_genes", 999)
        df = _gene_df()
        with pytest.raises(ValueError, match="expected exactly 999 gene rows"):
            validate_frozen_matrix(df, cfg, SAMPLE_IDS)

    def test_rejects_missing_sample_column(self, tmp_path):
        cfg = self._base_cfg(tmp_path)
        df = _gene_df().drop(columns=["FASR_Rep3"])
        with pytest.raises(ValueError, match="missing expected sample columns"):
            validate_frozen_matrix(df, cfg, SAMPLE_IDS)

    def test_rejects_leaked_mean_column(self, tmp_path):
        cfg = self._base_cfg(tmp_path)
        df = _gene_df()
        df["MCF7.mean.TPM"] = 1.0
        with pytest.raises(ValueError, match="possible mean leakage"):
            validate_frozen_matrix(df, cfg, SAMPLE_IDS)

    def test_rejects_negative_tpm(self, tmp_path):
        cfg = self._base_cfg(tmp_path)
        df = _gene_df()
        df.loc[0, "MCF7_Rep1"] = -1.0
        with pytest.raises(ValueError, match="negative TPM"):
            validate_frozen_matrix(df, cfg, SAMPLE_IDS)

    def test_rejects_non_finite_tpm(self, tmp_path):
        cfg = self._base_cfg(tmp_path)
        df = _gene_df()
        df.loc[0, "MCF7_Rep1"] = np.inf
        with pytest.raises(ValueError, match="non-finite TPM"):
            validate_frozen_matrix(df, cfg, SAMPLE_IDS)

    def test_rejects_null_gene_id(self, tmp_path):
        cfg = self._base_cfg(tmp_path)
        df = _gene_df()
        df.loc[0, "gene_id"] = None
        with pytest.raises(ValueError, match="null or blank gene_id"):
            validate_frozen_matrix(df, cfg, SAMPLE_IDS)


class TestChecksum:
    def test_verify_checksum_passes_on_match(self, tmp_path):
        p = tmp_path / "matrix.parquet"
        p.write_bytes(b"hello world")
        import hashlib

        expected = hashlib.sha256(b"hello world").hexdigest()
        cfg = _cfg(tmp_path, p, sha256=expected)
        assert verify_checksum(cfg) == expected

    def test_verify_checksum_fails_on_mismatch(self, tmp_path):
        p = tmp_path / "matrix.parquet"
        p.write_bytes(b"hello world")
        cfg = _cfg(tmp_path, p, sha256="0" * 64)
        with pytest.raises(ValueError, match="checksum mismatch"):
            verify_checksum(cfg)


class TestSampleMetadata:
    def test_rejects_duplicated_sample_id(self, tmp_path):
        meta = _meta_df()
        meta.loc[1, "sample_id"] = meta.loc[0, "sample_id"]
        meta_path = tmp_path / "meta.tsv"
        meta.to_csv(meta_path, sep="\t", index=False)
        cfg = _cfg(tmp_path, tmp_path / "x.parquet", sha256="unused")
        object.__setattr__(cfg, "sample_metadata_tsv", meta_path)
        with pytest.raises(ValueError, match="duplicated sample_id"):
            load_sample_metadata(cfg)

    def test_rejects_wrong_replicate_count(self, tmp_path):
        meta = _meta_df()
        meta.loc[0, "group"] = "TAMR"  # now MCF7 has 2, TAMR has 4
        meta_path = tmp_path / "meta.tsv"
        meta.to_csv(meta_path, sep="\t", index=False)
        cfg = _cfg(tmp_path, tmp_path / "x.parquet", sha256="unused")
        object.__setattr__(cfg, "sample_metadata_tsv", meta_path)
        with pytest.raises(ValueError, match="replicates per group"):
            load_sample_metadata(cfg)

    def test_rejects_missing_sample_count(self, tmp_path):
        meta = _meta_df().iloc[:-1]
        meta_path = tmp_path / "meta.tsv"
        meta.to_csv(meta_path, sep="\t", index=False)
        cfg = _cfg(tmp_path, tmp_path / "x.parquet", sha256="unused")
        object.__setattr__(cfg, "sample_metadata_tsv", meta_path)
        with pytest.raises(ValueError, match="expected exactly 9 samples"):
            load_sample_metadata(cfg)

    def test_accepts_valid_metadata(self, tmp_path):
        meta = _meta_df()
        meta_path = tmp_path / "meta.tsv"
        meta.to_csv(meta_path, sep="\t", index=False)
        cfg = _cfg(tmp_path, tmp_path / "x.parquet", sha256="unused")
        object.__setattr__(cfg, "sample_metadata_tsv", meta_path)
        loaded = load_sample_metadata(cfg)
        assert len(loaded) == 9


class TestWriters:
    def test_write_filtered_matrix_and_summary(self, tmp_path):
        df = _gene_df()
        filtered, record = filter_expression(df, SAMPLE_IDS, min_tpm=1.0, min_samples=3)
        cfg = _cfg(tmp_path, tmp_path / "x.parquet", sha256="unused")
        write_filtered_matrix(filtered, cfg)
        write_filtering_summary(record, cfg)
        assert cfg.filtered_gene_tpm_tsv.exists()
        assert cfg.filtering_summary_tsv.exists()
        reloaded = pd.read_csv(cfg.filtered_gene_tpm_tsv, sep="\t")
        assert set(reloaded["gene_id"]) == {"G1", "G3"}

    def test_filtered_matrix_gzip_is_byte_identical_across_writes(self, tmp_path):
        # Frozen outputs must be byte-reproducible, not merely numerically
        # equal -- a gzip header embedding the write timestamp would defeat
        # SHA256-based freeze verification even though the content is
        # unchanged.
        df = _gene_df()
        filtered, _ = filter_expression(df, SAMPLE_IDS, min_tpm=1.0, min_samples=3)
        cfg1 = _cfg(tmp_path / "run1", tmp_path / "run1" / "x.parquet", sha256="unused")
        cfg2 = _cfg(tmp_path / "run2", tmp_path / "run2" / "x.parquet", sha256="unused")
        import time

        write_filtered_matrix(filtered, cfg1)
        time.sleep(1.1)  # force a different wall-clock second between writes
        write_filtered_matrix(filtered, cfg2)

        import hashlib

        h1 = hashlib.sha256(cfg1.filtered_gene_tpm_tsv.read_bytes()).hexdigest()
        h2 = hashlib.sha256(cfg2.filtered_gene_tpm_tsv.read_bytes()).hexdigest()
        assert h1 == h2
