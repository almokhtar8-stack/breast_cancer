from collections import Counter
from dataclasses import replace
from pathlib import Path

import pandas as pd
import pytest
import yaml

from src.gse118713_prep import (
    SYMBOL_AMBIGUOUS,
    SYMBOL_MISSING,
    Gse118713Config,
    aggregate_to_gene_level,
    load_transcript_tpm,
    strip_ensembl_version,
)

CONFIG_PATH = Path(__file__).parent.parent / "config" / "config.yaml"

SAMPLE_COLUMNS = (
    "Rep1.MCF7.TPM",
    "Rep2.MCF7.TPM",
    "Rep3.MCF7.TPM",
    "Rep1.TAMR.TPM",
    "Rep2.TAMR.TPM",
    "Rep3.TAMR.TPM",
    "Rep1.FASR.TPM",
    "Rep2.FASR.TPM",
    "Rep3.FASR.TPM",
)
SAMPLE_IDS = (
    "MCF7_Rep1",
    "MCF7_Rep2",
    "MCF7_Rep3",
    "TAMR_Rep1",
    "TAMR_Rep2",
    "TAMR_Rep3",
    "FASR_Rep1",
    "FASR_Rep2",
    "FASR_Rep3",
)
GROUPS = ("MCF7",) * 3 + ("TAMR",) * 3 + ("FASR",) * 3
REPLICATES = (1, 2, 3) * 3
MEAN_COLUMNS = ("MCF7.mean.TPM", "TAMR.mean.TPM", "FASR.mean.TPM")


def _test_cfg(tmp_path) -> Gse118713Config:
    return Gse118713Config(
        source_path=tmp_path / "unused.tsv.gz",
        gene_tpm_parquet=tmp_path / "gene_tpm.parquet",
        sample_metadata_tsv=tmp_path / "meta.tsv",
        qc_tsv=tmp_path / "qc.tsv",
        sample_columns=SAMPLE_COLUMNS,
        sample_ids=SAMPLE_IDS,
        groups=GROUPS,
        replicates=REPLICATES,
        mean_columns=MEAN_COLUMNS,
    )


def _row(transcript_id, gene_id, symbol, values, mean_values=(0.0, 0.0, 0.0)):
    row = {"transcript.id": transcript_id, "gene.id": gene_id, "gene.symbol": symbol}
    row.update(dict(zip(SAMPLE_COLUMNS, values)))
    row.update(dict(zip(MEAN_COLUMNS, mean_values)))
    return row


def _synthetic_df() -> pd.DataFrame:
    rows = [
        _row("ENST001.1", "ENSG0001.3", "GENEA", [10, 10, 10, 1, 1, 1, 2, 2, 2]),
        _row("ENST002.1", "ENSG0001.3", "GENEA", [5, 5, 5, 1, 1, 1, 2, 2, 2]),
        _row("ENST010.1", "ENSG0005.1", "GENEA", [3, 3, 3, 0, 0, 0, 0, 0, 0]),
        _row("ENST003.1", "ENSG0002.1", None, [4, 4, 4, 0, 0, 0, 0, 0, 0]),
        _row("ENST004.1", "ENSG0003.2", "GENEC1", [1, 1, 1, 1, 1, 1, 1, 1, 1]),
        _row("ENST005.1", "ENSG0003.2", "GENEC2", [1, 1, 1, 1, 1, 1, 1, 1, 1]),
        _row("ENST006", "ENSG0004", "GENED", [7, 7, 7, 7, 7, 7, 7, 7, 7]),
    ]
    return pd.DataFrame(rows)


def test_transcript_tpm_sums_to_gene_tpm(tmp_path):
    cfg = _test_cfg(tmp_path)
    gene_df, _ = aggregate_to_gene_level(_synthetic_df(), cfg)
    gene_a = gene_df.set_index("gene_id").loc["ENSG0001"]
    assert gene_a["MCF7_Rep1"] == pytest.approx(15.0)  # 10 + 5, two transcripts of ENSG0001
    assert gene_a["TAMR_Rep1"] == pytest.approx(2.0)  # 1 + 1
    assert gene_a["FASR_Rep1"] == pytest.approx(4.0)  # 2 + 2


def test_total_tpm_conserved_after_aggregation(tmp_path):
    cfg = _test_cfg(tmp_path)
    df = _synthetic_df()
    gene_df, qc = aggregate_to_gene_level(df, cfg)
    for source_col, sample_id in zip(SAMPLE_COLUMNS, SAMPLE_IDS):
        assert gene_df[sample_id].sum() == pytest.approx(df[source_col].sum())
    assert qc.tpm_conserved_all_samples is True


def test_ensembl_version_suffix_removed():
    ids = pd.Series(["ENSG00000004059.8", "ENSG00000003056", "ENSG0001.12"])
    assert list(strip_ensembl_version(ids)) == ["ENSG00000004059", "ENSG00000003056", "ENSG0001"]


def test_duplicated_symbol_does_not_collapse_genes(tmp_path):
    cfg = _test_cfg(tmp_path)
    gene_df, _ = aggregate_to_gene_level(_synthetic_df(), cfg)
    genes_with_symbol_a = set(gene_df.loc[gene_df["gene_symbol"] == "GENEA", "gene_id"])
    assert genes_with_symbol_a == {"ENSG0001", "ENSG0005"}
    assert len(gene_df) == gene_df["gene_id"].nunique()


def test_missing_symbol_flagged(tmp_path):
    cfg = _test_cfg(tmp_path)
    gene_df, qc = aggregate_to_gene_level(_synthetic_df(), cfg)
    row = gene_df.set_index("gene_id").loc["ENSG0002"]
    assert row["symbol_mapping_status"] == SYMBOL_MISSING
    assert pd.isna(row["gene_symbol"])
    assert qc.missing_symbol_count == 1


def test_ambiguous_symbol_flagged_not_resolved(tmp_path):
    cfg = _test_cfg(tmp_path)
    gene_df, qc = aggregate_to_gene_level(_synthetic_df(), cfg)
    row = gene_df.set_index("gene_id").loc["ENSG0003"]
    assert row["symbol_mapping_status"] == SYMBOL_AMBIGUOUS
    assert pd.isna(row["gene_symbol"])
    assert qc.ambiguous_symbol_count == 1


def test_rejects_blank_gene_id(tmp_path):
    cfg = _test_cfg(tmp_path)
    df = _synthetic_df()
    df.loc[0, "gene.id"] = ""
    with pytest.raises(ValueError, match="blank gene.id"):
        aggregate_to_gene_level(df, cfg)


def test_rejects_null_gene_id(tmp_path):
    cfg = _test_cfg(tmp_path)
    df = _synthetic_df()
    df.loc[0, "gene.id"] = None
    with pytest.raises(ValueError, match="null gene.id"):
        aggregate_to_gene_level(df, cfg)


def test_mean_columns_excluded_from_output(tmp_path):
    df = _synthetic_df()
    source_path = tmp_path / "source.tsv.gz"
    df.to_csv(source_path, sep="\t", index=False)
    cfg = replace(_test_cfg(tmp_path), source_path=source_path)

    loaded = load_transcript_tpm(cfg)
    assert set(MEAN_COLUMNS).issubset(loaded.columns)  # present in the source file

    gene_df, _ = aggregate_to_gene_level(loaded, cfg)
    assert not set(MEAN_COLUMNS) & set(gene_df.columns)  # never carried into gene-level output
    assert set(SAMPLE_IDS).issubset(gene_df.columns)


def test_load_transcript_tpm_rejects_negative_values(tmp_path):
    df = _synthetic_df()
    df.loc[0, "Rep1.MCF7.TPM"] = -1.0
    source_path = tmp_path / "source.tsv.gz"
    df.to_csv(source_path, sep="\t", index=False)
    cfg = replace(_test_cfg(tmp_path), source_path=source_path)
    with pytest.raises(ValueError, match="negative"):
        load_transcript_tpm(cfg)


def test_real_config_has_exact_nine_sample_metadata():
    with open(CONFIG_PATH) as f:
        config = yaml.safe_load(f)
    cfg = Gse118713Config.from_config(config)

    assert len(cfg.sample_ids) == 9
    assert len(set(cfg.sample_ids)) == 9
    assert Counter(cfg.groups) == {"MCF7": 3, "TAMR": 3, "FASR": 3}
    for group in ("MCF7", "TAMR", "FASR"):
        reps = sorted(r for r, g in zip(cfg.replicates, cfg.groups) if g == group)
        assert reps == [1, 2, 3]
    assert set(cfg.mean_columns) == {"MCF7.mean.TPM", "TAMR.mean.TPM", "FASR.mean.TPM"}
    assert set(cfg.mean_columns).isdisjoint(set(cfg.sample_columns))


def test_no_hardcoded_absolute_path_in_module_source():
    text = (Path(__file__).parent.parent / "src" / "gse118713_prep.py").read_text()
    assert "/ibex" not in text
