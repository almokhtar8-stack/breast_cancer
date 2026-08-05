from pathlib import Path

import pandas as pd
import pytest
import yaml

from src.labels import (
    ARM_E2,
    ARM_E2_4OHT,
    CONTROL_T0,
    REPLICATES,
    REQUIRED_ARMS,
    _validate_header,
    _validate_values,
    add_log2_cpm,
    compute_effect_table,
    filter_low_representation_genes,
    fit_gene_treatment_effects,
    load_raw_counts,
)

CONFIG_PATH = Path(__file__).parent.parent / "config" / "config.yaml"


def _raw_data_path() -> Path:
    with open(CONFIG_PATH) as f:
        config = yaml.safe_load(f)
    return Path(config["data"]["raw"]["hany_data_s1"])


def test_raw_input_row_count():
    raw = load_raw_counts(_raw_data_path())
    assert len(raw) == 76_441


def test_target_conditions_present():
    raw = load_raw_counts(_raw_data_path())
    for arm in (ARM_E2, ARM_E2_4OHT):
        for rep in REPLICATES:
            assert f"{arm}__{rep}" in raw.columns


def test_sign_convention_negative_for_depleted_gene():
    """A gene deliberately more depleted under 4-OHT must get a negative
    effect_size, per PREANALYSIS.md §2's sign convention."""
    rows = [
        # gene, sgRNA, T0 R1, T0 R2, E2 R1, E2 R2, 4-OHT R1, 4-OHT R2
        ("DEPLETED", "g1", 120, 130, 200, 210, 40, 45),
        ("DEPLETED", "g2", 110, 100, 190, 205, 35, 42),
        ("DEPLETED", "g3", 140, 125, 220, 195, 50, 38),
        ("ENRICHED", "g4", 120, 130, 40, 45, 200, 210),
        ("ENRICHED", "g5", 110, 100, 35, 42, 190, 205),
        ("ENRICHED", "g6", 140, 125, 50, 38, 220, 195),
    ]
    columns = [
        "Gene",
        "sgRNA",
        f"{CONTROL_T0}__R1",
        f"{CONTROL_T0}__R2",
        f"{ARM_E2}__R1",
        f"{ARM_E2}__R2",
        f"{ARM_E2_4OHT}__R1",
        f"{ARM_E2_4OHT}__R2",
    ]
    df = pd.DataFrame(rows, columns=columns)

    sample_cols = [f"{arm}__{rep}" for arm in (ARM_E2, ARM_E2_4OHT) for rep in REPLICATES]
    normalised = add_log2_cpm(df, sample_cols)

    from src.labels import _to_long

    long_df = _to_long(normalised)
    gene_fits = fit_gene_treatment_effects(long_df)
    effect_table = compute_effect_table(gene_fits).set_index("gene")

    assert effect_table.loc["DEPLETED", "effect_size"] < 0
    assert effect_table.loc["ENRICHED", "effect_size"] > 0


def _well_formed_columns() -> list[str]:
    return ["sgRNA", "Gene"] + [f"{arm}__{rep}" for arm in REQUIRED_ARMS for rep in REPLICATES]


def test_validate_header_accepts_well_formed_columns():
    _validate_header(_well_formed_columns())  # must not raise


def test_validate_header_rejects_missing_required_column():
    columns = [c for c in _well_formed_columns() if c != f"{ARM_E2_4OHT}__R1"]
    with pytest.raises(ValueError, match="missing"):
        _validate_header(columns)


def test_validate_header_rejects_duplicate_column():
    columns = _well_formed_columns() + [f"{ARM_E2}__R1"]
    with pytest.raises(ValueError, match="duplicate"):
        _validate_header(columns)


def test_validate_values_rejects_duplicate_sgrna():
    sample_cols = [f"{arm}__{rep}" for arm in REQUIRED_ARMS for rep in REPLICATES]
    df = pd.DataFrame(
        [["g1", "GENE_A"] + [10.0] * len(sample_cols), ["g1", "GENE_B"] + [10.0] * len(sample_cols)],
        columns=["sgRNA", "Gene"] + sample_cols,
    )
    with pytest.raises(ValueError, match="duplicate sgRNA"):
        _validate_values(df, sample_cols)


def test_validate_values_rejects_non_finite_counts():
    sample_cols = [f"{arm}__{rep}" for arm in REQUIRED_ARMS for rep in REPLICATES]
    df = pd.DataFrame(
        [["g1", "GENE_A"] + [float("nan")] * len(sample_cols)],
        columns=["sgRNA", "Gene"] + sample_cols,
    )
    with pytest.raises(ValueError, match="non-numeric or non-finite"):
        _validate_values(df, sample_cols)


def test_low_representation_filter_counts_distinct_guides_not_rows():
    """A gene with a duplicated guide row must still count as under-represented."""
    sample_cols = [f"{arm}__{rep}" for arm in REQUIRED_ARMS for rep in REPLICATES]
    rows = [
        ["g1", "TWOGUIDES"] + [10.0] * len(sample_cols),
        ["g1", "TWOGUIDES"] + [10.0] * len(sample_cols),  # duplicate row, same guide
        ["g2", "TWOGUIDES"] + [10.0] * len(sample_cols),
    ]
    df = pd.DataFrame(rows, columns=["sgRNA", "Gene"] + sample_cols)

    kept, step = filter_low_representation_genes(df, min_guides=3)

    assert kept.empty
    assert step.genes_out == 0
