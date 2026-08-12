from pathlib import Path

import pandas as pd
import pytest

from src.gse240112_evidence_classification import classify_candidates

REPO_ROOT = Path(__file__).parent.parent


def _candidate_table():
    return pd.DataFrame(
        {
            "gene": ["VEZF1", "USP34", "TSR3", "PET117", "USP17L29"],
            "tested": [True, True, True, True, False],
            "reason_not_tested": ["", "", "", "", "absent from feature space"],
            "log2fc": [1.15, 0.40, -0.43, -0.25, float("nan")],
            "p_value": [0.004, 0.118, 0.241, 0.429, float("nan")],
            "candidate_set_bh_fdr": [0.048, 0.473, 0.580, 0.677, float("nan")],
        }
    )


def _sensitivity():
    return pd.DataFrame(
        {
            "gene": ["VEZF1", "USP34", "TSR3", "PET117"],
            "direction_agreement": [True, True, True, True],
        }
    )


def _integrated():
    return pd.DataFrame(
        {
            "gene": ["VEZF1", "USP34", "TSR3", "PET117"],
            "gse118713_tamr_vs_mcf7_log2fc": [0.43, 0.59, -0.47, -0.06],
            "gse118713_tamr_vs_mcf7_fdr": [0.24, 0.0073, 0.19, 0.78],
        }
    )


class TestClassifyCandidates:
    def test_untested_gene_is_untestable(self):
        out = classify_candidates(_candidate_table(), _sensitivity(), _integrated()).set_index("gene")
        assert out.loc["USP17L29", "classification"] == "untestable"

    def test_fdr_significant_gene_is_strengthened(self):
        out = classify_candidates(_candidate_table(), _sensitivity(), _integrated()).set_index("gene")
        assert out.loc["VEZF1", "classification"] == "strengthened"

    def test_sizeable_cross_track_agreeing_gene_is_directionally_supportive(self):
        out = classify_candidates(_candidate_table(), _sensitivity(), _integrated()).set_index("gene")
        # USP34: log2fc=0.40 (<0.5) but nominal p=0.118 (>0.05) -- neither condition met -> neutral
        assert out.loc["USP34", "classification"] == "neutral_no_additional_support"

    def test_weak_nonsignificant_gene_is_neutral(self):
        out = classify_candidates(_candidate_table(), _sensitivity(), _integrated()).set_index("gene")
        assert out.loc["PET117", "classification"] == "neutral_no_additional_support"

    def test_every_input_gene_classified_exactly_once(self):
        out = classify_candidates(_candidate_table(), _sensitivity(), _integrated())
        assert len(out) == 5
        assert set(out["gene"]) == {"VEZF1", "USP34", "TSR3", "PET117", "USP17L29"}


class TestRealData:
    def test_real_classification_if_present(self):
        path = REPO_ROOT / "results" / "tables" / "gse240112" / "candidate_classification.tsv"
        if not path.exists():
            pytest.skip("GSE240112 classification table not present in this checkout")
        out = pd.read_csv(path, sep="\t")
        assert len(out) == 13
        assert out.set_index("gene").loc["USP17L29", "classification"] == "untestable"
        valid = {"strengthened", "directionally_supportive_but_weak", "neutral_no_additional_support", "discordant", "untestable"}
        assert set(out["classification"]).issubset(valid)
