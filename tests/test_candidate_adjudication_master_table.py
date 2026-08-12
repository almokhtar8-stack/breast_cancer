import inspect

import numpy as np
import pandas as pd

import src.candidate_adjudication_leaders as leaders_module
import src.candidate_adjudication_near_misses as near_misses_module
from src.candidate_adjudication_master_table import _direction_label, _rank_within_dataset

RESISTANCE_DATASET_COLUMN_FRAGMENTS = ["gse118713", "gse240112_tumor", "gse111151"]


class TestRankWithinDataset:
    def test_lowest_fdr_gets_rank_one(self):
        effect = pd.Series([1.0, -2.0, 0.5])
        p = pd.Series([0.01, 0.001, 0.5])
        fdr = pd.Series([0.05, 0.01, 0.6])
        gene = pd.Series(["A", "B", "C"])
        ranks = _rank_within_dataset(effect, p, fdr, gene)
        assert ranks.iloc[1] == 1  # B has the lowest FDR (0.01) so it gets rank 1

    def test_untestable_gene_gets_nan_rank(self):
        effect = pd.Series([1.0, np.nan])
        p = pd.Series([0.01, np.nan])
        fdr = pd.Series([0.05, np.nan])
        gene = pd.Series(["A", "B"])
        ranks = _rank_within_dataset(effect, p, fdr, gene)
        assert pd.isna(ranks.iloc[1])

    def test_ties_broken_by_gene_symbol_ascending(self):
        effect = pd.Series([1.0, 1.0])
        p = pd.Series([0.01, 0.01])
        fdr = pd.Series([0.05, 0.05])
        gene = pd.Series(["ZGENE", "AGENE"])
        ranks = _rank_within_dataset(effect, p, fdr, gene)
        assert ranks.iloc[1] == 1  # AGENE (alphabetically first) wins the tie
        assert ranks.iloc[0] == 2


class TestDirectionLabel:
    def test_positive_is_up(self):
        assert _direction_label(1.5) == "up"

    def test_negative_is_down(self):
        assert _direction_label(-0.3) == "down"

    def test_nan_is_not_testable(self):
        assert _direction_label(float("nan")) == "not_testable"

    def test_zero_is_unchanged(self):
        assert _direction_label(0.0) == "unchanged"


class TestResistanceDatasetScopeDiscipline:
    def test_resistance_leaderboard_ranking_and_fdr_counts_come_only_from_the_frozen_resistance_consensus_input(self):
        """`build_top_resistance_genes` is allowed to *display* GSE245601/
        GSE240112-epithelial columns as extra context (e.g. for the
        human-support label), but the resistance ranking/FDR-count columns
        it reports must be passed through unchanged from the already
        -frozen resistance-consensus table (built with only GSE118713,
        GSE240112 tumor-cell, and GSE111151), never recomputed from
        GSE245601 or the epithelial track."""
        resistance = pd.DataFrame(
            {
                "gene": ["G1", "G2"], "resistance_datasets_testable": [3, 3], "resistance_up_count": [3, 0],
                "resistance_down_count": [0, 3], "resistance_fdr05_count": [3, 2], "resistance_nominal_p05_count": [3, 3],
                "resistance_top10pct_count": [3, 2], "resistance_top20pct_count": [3, 2], "resistance_median_percentile": [0.99, 0.95],
                "resistance_direction_consensus": ["all_up", "all_down"],
            }
        )
        wide = pd.DataFrame(
            {
                "gene": ["G1", "G2"], "gse118713_log2fc": [1.0, -1.0], "gse118713_fdr": [0.01, 0.02],
                "gse240112_tumor_log2fc": [1.0, -1.0], "gse240112_tumor_fdr": [0.01, 0.02],
                "gse111151_log2fc": [1.0, -1.0], "gse111151_fdr": [0.01, 0.5],
                "crispr_effect": [0.1, -0.1], "crispr_fdr": [0.9, 0.9], "crispr_direction": ["tolerance_associated_KO", "sensitising_KO"],
                # deliberately extreme/misleading GSE245601 values -- must NOT change resistance_fdr05_count or consensus
                "gse245601_epi_fdr": [1e-30, 1e-30], "gse240112_epi_log2fc": [999.0, -999.0],
            }
        )
        ranked = pd.DataFrame({"gene": ["G1", "G2"], "global_rank": [1, 2]})
        categories = pd.DataFrame({"gene": ["G1", "G2"], "evidence_category": ["RNA_RESISTANCE_CONSENSUS", "RNA_RESISTANCE_CONSENSUS"]})

        out = leaders_module.build_top_resistance_genes(resistance, wide, ranked, categories, top_n=2)
        assert out["resistance_fdr05_count"].tolist() == [3, 2]
        assert out["resistance_direction_consensus"].tolist() == ["all_up", "all_down"]

    def test_near_misses_module_does_not_use_gse245601_for_resistance_reasons(self):
        src_text = inspect.getsource(near_misses_module.find_broad_near_misses)
        for line in src_text.splitlines():
            if "resistance" in line.lower():
                assert "gse245601" not in line.lower(), f"resistance-labeled line references GSE245601: {line}"
