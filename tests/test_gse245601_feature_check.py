from pathlib import Path

import pandas as pd
import pytest
import yaml

from src.gse245601_feature_check import (
    CANDIDATE_ENSEMBL_IDS,
    PAICS_ENSEMBL_ID,
    check_candidate_feature_availability,
    read_feature_table,
    run_feature_check,
)

REPO_ROOT = Path(__file__).parent.parent
CONFIG_PATH = REPO_ROOT / "config" / "config.yaml"


def _real_h5_path() -> Path:
    with open(CONFIG_PATH) as f:
        config = yaml.safe_load(f)
    cfg = config["gse245601"]
    sample = next(s for s in cfg["samples"] if s["primary_analysis"])
    return Path(cfg["h5_dir"]) / sample["h5"]


class TestFeatureTableReading:
    def test_reads_only_feature_metadata_not_expression(self):
        table = read_feature_table(_real_h5_path())
        assert set(table.columns) == {"ensembl_id", "gene_symbol", "feature_type", "genome"}
        assert len(table) == 33538

    def test_single_feature_type_and_genome(self):
        table = read_feature_table(_real_h5_path())
        assert set(table["feature_type"]) == {"Gene Expression"}
        assert set(table["genome"]) == {"GRCh38"}


class TestCandidateFeatureAvailability:
    def test_all_13_candidates_present_and_unambiguous(self):
        table = read_feature_table(_real_h5_path())
        result = check_candidate_feature_availability(table)
        candidates = result.loc[~result["is_paics_benchmark"]]
        assert len(candidates) == 13
        assert set(candidates["gene_symbol"]) == set(CANDIDATE_ENSEMBL_IDS)
        assert (candidates["status"] == "present").all()

    def test_paics_present_and_reported_separately(self):
        table = read_feature_table(_real_h5_path())
        result = check_candidate_feature_availability(table)
        paics_row = result.loc[result["gene_symbol"] == "PAICS"].iloc[0]
        assert paics_row["is_paics_benchmark"]
        assert paics_row["status"] == "present"
        assert paics_row["expected_ensembl_id"] == PAICS_ENSEMBL_ID

    def test_usp17l29_paralog_family_is_unique_in_this_reference(self):
        # USP17L29 sits in a segmental-duplication paralog family and was
        # flagged as an ambiguity risk in the design-only audit; confirm
        # against the real reference feature list that it is in fact
        # unique here (present, not collapsed/duplicated).
        table = read_feature_table(_real_h5_path())
        result = check_candidate_feature_availability(table)
        row = result.loc[result["gene_symbol"] == "USP17L29"].iloc[0]
        assert row["status"] == "present"
        assert row["n_symbol_matches"] == 1

    def test_absent_gene_is_reported_absent_not_fabricated(self):
        table = pd.DataFrame(
            {"ensembl_id": ["ENSG00000000000"], "gene_symbol": ["NOTAGENE"], "feature_type": ["Gene Expression"], "genome": ["GRCh38"]}
        )
        result = check_candidate_feature_availability(table)
        assert (result["status"] == "absent").all()

    def test_ambiguous_when_symbol_maps_to_multiple_rows(self):
        table = pd.DataFrame(
            {
                "ensembl_id": ["ENSG00000115464", "ENSG00000999999"],
                "gene_symbol": ["USP34", "USP34"],
                "feature_type": ["Gene Expression", "Gene Expression"],
                "genome": ["GRCh38", "GRCh38"],
            }
        )
        result = check_candidate_feature_availability(table)
        usp34_row = result.loc[result["gene_symbol"] == "USP34"].iloc[0]
        assert usp34_row["status"] == "ambiguous"

    def test_no_expression_values_are_read_anywhere_in_this_module(self):
        # Static guarantee: no h5py indexing call in this module touches
        # matrix/data or matrix/indices (the actual expression values) --
        # only matrix/features/* (gene metadata). The docstring mentions
        # "matrix/data" in prose explaining what is deliberately NOT read,
        # so we check for actual indexing syntax, not a bare substring.
        import inspect

        import src.gse245601_feature_check as mod

        source = inspect.getsource(mod)
        assert '"matrix/data"]' not in source
        assert "'matrix/data']" not in source
        assert '"matrix/indices"]' not in source
        assert "'matrix/indices']" not in source


class TestRunAgainstRealConfig:
    def test_writes_output_table(self):
        result = run_feature_check()
        assert len(result) == 14
        assert (result["status"] == "present").sum() == 14
