from pathlib import Path

import pandas as pd
import pytest

from src.gse111151_evidence_classification import classify_candidates, compute_cell_line_consistency

REPO_ROOT = Path(__file__).parent.parent


def _candidate_table():
    return pd.DataFrame(
        {
            "gene": ["STRONG", "WEAK", "OPPOSITE", "FLAT", "MISSING"],
            "tested": [True, True, True, True, False],
            "reason_not_tested": ["", "", "", "", "absent"],
            "log2fc": [0.6, 0.3, -0.5, 0.05, float("nan")],
            "p_value": [0.01, 0.20, 0.03, 0.90, float("nan")],
            "candidate_set_bh_fdr": [0.03, 0.5, 0.5, 0.9, float("nan")],
        }
    )


def _sample_level():
    rows = []
    # STRONG: FDR-significant, consistency irrelevant
    # WEAK: 3/4 cell lines consistent with positive overall direction
    for gene, deltas in [("STRONG", [1, 1, 1, 1]), ("WEAK", [1, 1, 1, -1]), ("OPPOSITE", [-1, -1, -1, 1]), ("FLAT", [1, -1, 1, -1])]:
        for i, cl in enumerate(["MCF-7", "T-47D", "ZR-75-1", "BT-474"]):
            rows.append({"gene": gene, "sample_id": f"{cl}", "cell_line": cl, "resistance_status": "parental", "log2cpm": 5.0})
            rows.append({"gene": gene, "sample_id": f"{cl}_Tam1", "cell_line": cl, "resistance_status": "resistant", "log2cpm": 5.0 + deltas[i]})
    return pd.DataFrame(rows)


def _gse118713_bulk():
    return pd.DataFrame(
        {
            "gene_symbol": ["STRONG", "WEAK", "OPPOSITE", "FLAT"],
            "tamr_vs_mcf7_log2fc": [0.5, 0.4, 0.5, 0.1],
            "tamr_vs_mcf7_fdr": [0.01, 0.3, 0.01, 0.5],
        }
    )


class TestComputeCellLineConsistency:
    def test_counts_consistent_cell_lines(self):
        sl = _sample_level()
        n_consistent, n_with_both = compute_cell_line_consistency(sl, "WEAK", overall_log2fc=0.3)
        assert n_with_both == 4
        assert n_consistent == 3


class TestClassifyCandidates:
    def test_untested_gene_is_untestable(self):
        out = classify_candidates(_candidate_table(), _sample_level(), _gse118713_bulk()).set_index("gene")
        assert out.loc["MISSING", "classification"] == "untestable"

    def test_fdr_significant_gene_is_independently_supported(self):
        out = classify_candidates(_candidate_table(), _sample_level(), _gse118713_bulk()).set_index("gene")
        assert out.loc["STRONG", "classification"] == "independently_supported"

    def test_consistent_cell_lines_gives_directionally_supportive(self):
        out = classify_candidates(_candidate_table(), _sample_level(), _gse118713_bulk()).set_index("gene")
        assert out.loc["WEAK", "classification"] == "directionally_supportive_but_weak"

    def test_high_p_value_with_consistency_alone_is_not_supportive(self):
        # cell-line consistency without a companion nominal p<0.3 must not qualify -- a 3/4 or 4/4
        # split happens too often by chance alone (regression test for the pre-Codex-review rule fix)
        ct = _candidate_table().copy()
        ct.loc[ct["gene"] == "WEAK", "p_value"] = 0.95
        ct.loc[ct["gene"] == "WEAK", "candidate_set_bh_fdr"] = 0.95
        out = classify_candidates(ct, _sample_level(), _gse118713_bulk()).set_index("gene")
        assert out.loc["WEAK", "classification"] == "neutral_no_additional_support"

    def test_opposing_significant_bulk_direction_is_discordant(self):
        out = classify_candidates(_candidate_table(), _sample_level(), _gse118713_bulk()).set_index("gene")
        assert out.loc["OPPOSITE", "classification"] == "discordant"

    def test_flat_inconsistent_gene_is_neutral(self):
        out = classify_candidates(_candidate_table(), _sample_level(), _gse118713_bulk()).set_index("gene")
        assert out.loc["FLAT", "classification"] == "neutral_no_additional_support"

    def test_missing_integration_table_does_not_crash(self):
        out = classify_candidates(_candidate_table(), _sample_level(), pd.DataFrame())
        assert len(out) == 5


class TestRealData:
    def test_real_classification_if_present(self):
        path = REPO_ROOT / "results" / "tables" / "gse111151" / "candidate_classification.tsv"
        if not path.exists():
            pytest.skip("GSE111151 classification table not present in this checkout")
        out = pd.read_csv(path, sep="\t")
        assert len(out) == 13
        assert out.set_index("gene").loc["USP17L29", "classification"] == "untestable"
        valid = {"independently_supported", "directionally_supportive_but_weak", "neutral_no_additional_support", "discordant", "untestable"}
        assert set(out["classification"]).issubset(valid)
