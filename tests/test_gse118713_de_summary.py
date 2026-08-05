import numpy as np
import pandas as pd
import pytest

from src.gse118713_de_summary import build_de_summary

CONTRASTS = ("TAMR_vs_MCF7", "FASR_vs_MCF7", "TAMR_vs_FASR")


def _row(gene_id, contrast, log2fc, fdr):
    return {"gene_id": gene_id, "gene_symbol": gene_id, "log2fc": log2fc, "fdr": fdr, "contrast": contrast}


class TestBuildDeSummary:
    def test_genes_tested_matches_row_count_per_contrast(self):
        rows = []
        for contrast in CONTRASTS:
            for i in range(5):
                rows.append(_row(f"G{i}", contrast, 0.1 * i, 0.5))
        de_df = pd.DataFrame(rows)
        summary = build_de_summary(de_df, CONTRASTS)
        assert (summary["genes_tested"] == 5).all()

    def test_significance_counts_match_manual_threshold(self):
        rows = [
            _row("G1", "TAMR_vs_MCF7", 2.0, 0.01),   # significant positive
            _row("G2", "TAMR_vs_MCF7", -3.0, 0.02),  # significant negative
            _row("G3", "TAMR_vs_MCF7", 1.0, 0.5),    # not significant
        ]
        for contrast in ("FASR_vs_MCF7", "TAMR_vs_FASR"):
            rows.append(_row("G1", contrast, 1.0, 0.5))
        de_df = pd.DataFrame(rows)
        summary = build_de_summary(de_df, CONTRASTS).set_index("contrast")
        row = summary.loc["TAMR_vs_MCF7"]
        assert row["n_significant_fdr_lt_0_05"] == 2
        assert row["n_significant_positive"] == 1
        assert row["n_significant_negative"] == 1

    def test_median_abs_effect_and_min_fdr_computed_from_data(self):
        rows = [
            _row("G1", "TAMR_vs_MCF7", 1.0, 0.3),
            _row("G2", "TAMR_vs_MCF7", -3.0, 0.01),
            _row("G3", "TAMR_vs_MCF7", 5.0, 0.4),
        ]
        for contrast in ("FASR_vs_MCF7", "TAMR_vs_FASR"):
            rows.append(_row("G1", contrast, 1.0, 0.5))
        de_df = pd.DataFrame(rows)
        summary = build_de_summary(de_df, CONTRASTS).set_index("contrast")
        row = summary.loc["TAMR_vs_MCF7"]
        assert row["median_abs_log2fc"] == pytest.approx(np.median([1.0, 3.0, 5.0]))
        assert row["min_fdr"] == pytest.approx(0.01)

    def test_no_hardcoded_results_output_tracks_input(self):
        rows_a = [_row(f"G{i}", c, 0.1, 0.01) for i in range(10) for c in CONTRASTS]
        rows_b = [_row(f"G{i}", c, 0.1, 0.9) for i in range(10) for c in CONTRASTS]
        summary_a = build_de_summary(pd.DataFrame(rows_a), CONTRASTS)
        summary_b = build_de_summary(pd.DataFrame(rows_b), CONTRASTS)
        assert (summary_a["n_significant_fdr_lt_0_05"] == 10).all()
        assert (summary_b["n_significant_fdr_lt_0_05"] == 0).all()

    def test_missing_contrast_raises(self):
        rows = [_row("G1", "TAMR_vs_MCF7", 1.0, 0.01)]
        with pytest.raises(ValueError, match="missing contrasts"):
            build_de_summary(pd.DataFrame(rows), CONTRASTS)
