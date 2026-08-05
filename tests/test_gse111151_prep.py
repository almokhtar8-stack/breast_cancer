from pathlib import Path

import pandas as pd
import pytest
import yaml

from src.gse111151_prep import (
    Gse111151Config,
    SampleSpec,
    build_log2cpm_matrix,
    read_one_sample,
)

CONFIG_PATH = Path(__file__).parent.parent / "config" / "config.yaml"
COLUMNS = ["EnsEMBL_GenID", "gene_name", "counts", "CPM (log2)", "CPM_batch (log2)"]


def _write_sample(path, rows):
    pd.DataFrame(rows, columns=COLUMNS).to_csv(path, sep="\t", index=False)


def _base_rows():
    return [
        ["ENSG0001", "GENEA", 10, 1.0, 1.1],
        ["ENSG0002", "GENEB", 0, -5.0, -5.1],
        ["ENSG0003", "GENEB", 5, 0.5, 0.6],  # duplicate symbol GENEB on a distinct id
    ]


def _cfg_for(tmp_path, filenames, sample_ids) -> Gse111151Config:
    samples = tuple(
        SampleSpec(
            gsm=f"GSM{i}",
            filename=filenames[i],
            sample_id=sample_ids[i],
            parental_line=f"LINE{i}",
            status="parental",
            derivative_id=None,
            paired_parental_sample_id=None,
        )
        for i in range(len(filenames))
    )
    return Gse111151Config(
        sample_dir=tmp_path,
        log2cpm_parquet=tmp_path / "out.parquet",
        sample_metadata_tsv=tmp_path / "meta.tsv",
        qc_tsv=tmp_path / "qc.tsv",
        samples=samples,
    )


def test_log2cpm_used_batch_and_counts_excluded(tmp_path):
    _write_sample(tmp_path / "a.txt.gz", _base_rows())
    rows_b = [[r[0], r[1], r[2] + 1, r[3] + 10, r[4] + 10] for r in _base_rows()]
    _write_sample(tmp_path / "b.txt.gz", rows_b)
    cfg = _cfg_for(tmp_path, ["a.txt.gz", "b.txt.gz"], ["A", "B"])
    tables = {"A": read_one_sample(tmp_path / "a.txt.gz", "A"), "B": read_one_sample(tmp_path / "b.txt.gz", "B")}

    matrix, qc = build_log2cpm_matrix(tables, cfg)

    assert qc.log2cpm_column_used is True
    assert qc.log2cpm_batch_column_used is False
    assert qc.counts_column_used is False

    indexed = matrix.set_index("gene_id")
    assert indexed.loc["ENSG0001", "A"] == pytest.approx(1.0)  # CPM (log2), not counts=10 or batch=1.1
    assert indexed.loc["ENSG0001", "B"] == pytest.approx(11.0)


def test_raw_counts_not_used_as_prepared_matrix(tmp_path):
    _write_sample(tmp_path / "a.txt.gz", _base_rows())
    cfg = _cfg_for(tmp_path, ["a.txt.gz"], ["A"])
    tables = {"A": read_one_sample(tmp_path / "a.txt.gz", "A")}

    matrix, _ = build_log2cpm_matrix(tables, cfg)

    assert "counts" not in matrix.columns
    value = matrix.set_index("gene_id").loc["ENSG0001", "A"]
    assert value != 10
    assert value == pytest.approx(1.0)


def test_duplicate_symbol_not_collapsed(tmp_path):
    _write_sample(tmp_path / "a.txt.gz", _base_rows())
    cfg = _cfg_for(tmp_path, ["a.txt.gz"], ["A"])
    tables = {"A": read_one_sample(tmp_path / "a.txt.gz", "A")}

    matrix, qc = build_log2cpm_matrix(tables, cfg)

    assert len(matrix) == 3  # ENSG0002 and ENSG0003 both named GENEB, kept as separate rows
    assert qc.duplicate_symbol_count == 2


def test_no_silent_row_loss(tmp_path):
    _write_sample(tmp_path / "a.txt.gz", _base_rows())
    cfg = _cfg_for(tmp_path, ["a.txt.gz"], ["A"])
    tables = {"A": read_one_sample(tmp_path / "a.txt.gz", "A")}

    _, qc = build_log2cpm_matrix(tables, cfg)

    assert qc.rows_in == qc.rows_out == len(_base_rows())


def test_rejects_duplicate_gene_id(tmp_path):
    rows = _base_rows() + [["ENSG0001", "GENEA", 1, 0.0, 0.0]]
    _write_sample(tmp_path / "a.txt.gz", rows)
    with pytest.raises(ValueError, match="duplicate"):
        read_one_sample(tmp_path / "a.txt.gz", "A")


def test_rejects_blank_gene_id(tmp_path):
    rows = _base_rows()
    rows[0][0] = ""
    _write_sample(tmp_path / "a.txt.gz", rows)
    with pytest.raises(ValueError, match="blank or null"):
        read_one_sample(tmp_path / "a.txt.gz", "A")


def test_rejects_non_integer_counts(tmp_path):
    rows = _base_rows()
    rows[0][2] = 10.5
    _write_sample(tmp_path / "a.txt.gz", rows)
    with pytest.raises(ValueError, match="non-integer counts"):
        read_one_sample(tmp_path / "a.txt.gz", "A")


def test_rejects_mismatched_gene_ids_across_files(tmp_path):
    _write_sample(tmp_path / "a.txt.gz", _base_rows())
    rows_b = _base_rows()
    rows_b[-1][0] = "ENSG9999"
    _write_sample(tmp_path / "b.txt.gz", rows_b)
    cfg = _cfg_for(tmp_path, ["a.txt.gz", "b.txt.gz"], ["A", "B"])
    tables = {"A": read_one_sample(tmp_path / "a.txt.gz", "A"), "B": read_one_sample(tmp_path / "b.txt.gz", "B")}
    with pytest.raises(ValueError, match="gene ID set mismatch"):
        build_log2cpm_matrix(tables, cfg)


def test_changed_gene_order_is_detected_but_aligned_by_identifier(tmp_path):
    _write_sample(tmp_path / "a.txt.gz", _base_rows())
    _write_sample(tmp_path / "b.txt.gz", list(reversed(_base_rows())))
    cfg = _cfg_for(tmp_path, ["a.txt.gz", "b.txt.gz"], ["A", "B"])
    tables = {"A": read_one_sample(tmp_path / "a.txt.gz", "A"), "B": read_one_sample(tmp_path / "b.txt.gz", "B")}

    matrix, qc = build_log2cpm_matrix(tables, cfg)

    assert qc.row_order_identical_across_files is False
    indexed = matrix.set_index("gene_id")
    # correct despite the row-order difference, because alignment is by gene_id, not position
    assert indexed.loc["ENSG0002", "A"] == pytest.approx(-5.0)
    assert indexed.loc["ENSG0002", "B"] == pytest.approx(-5.0)
    assert indexed.loc["ENSG0001", "A"] == pytest.approx(1.0)
    assert indexed.loc["ENSG0001", "B"] == pytest.approx(1.0)


def test_real_config_has_exact_eleven_samples_and_seven_pairs():
    with open(CONFIG_PATH) as f:
        config = yaml.safe_load(f)
    cfg = Gse111151Config.from_config(config)

    assert len(cfg.samples) == 11
    parental = [s for s in cfg.samples if s.status == "parental"]
    resistant = [s for s in cfg.samples if s.status == "resistant"]
    assert len(parental) == 4
    assert len(resistant) == 7

    parental_by_line = {s.parental_line: s.sample_id for s in parental}
    pairs = {(r.sample_id, r.paired_parental_sample_id) for r in resistant}
    assert len(pairs) == 7
    for r in resistant:
        assert parental_by_line[r.parental_line] == r.paired_parental_sample_id


def test_no_hardcoded_absolute_path_in_module_source():
    text = (Path(__file__).parent.parent / "src" / "gse111151_prep.py").read_text()
    assert "/ibex" not in text
