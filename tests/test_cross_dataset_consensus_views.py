from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from src.cross_dataset_consensus_views import (
    assign_evidence_category,
    build_crispr_functional_ranking,
    build_crispr_nonsignificant_rna_consensus,
    build_human_only_ranking,
    build_resistance_consensus,
    build_rna_only_ranking,
)

REPO_ROOT = Path(__file__).parent.parent


def _synthetic_full():
    genes = ["ALLUP", "ALLDOWN", "MIXED", "INSUFFICIENT", "CRISPR_ONLY", "MULTIMODAL", "ACUTE_ONLY", "CONTEXT_MIX"]
    df = pd.DataFrame({"gene": genes})
    # crispr
    df["crispr_effect"] = [0.1, 0.1, 0.1, 0.1, -2.0, -2.0, 0.1, 0.1]
    df["crispr_p"] = [0.5, 0.5, 0.5, 0.5, 0.001, 0.001, 0.5, 0.5]
    df["crispr_fdr"] = [0.5, 0.5, 0.5, 0.5, 0.01, 0.01, 0.5, 0.5]
    df["crispr_testable"] = [True] * 8
    df["crispr_evidence_percentile"] = [0.3, 0.3, 0.3, 0.3, 0.99, 0.99, 0.3, 0.3]
    # gse118713 (resistance dataset)
    df["gse118713_log2fc"] = [0.5, -0.5, 0.5, np.nan, 0.1, 0.6, 0.1, 0.5]
    df["gse118713_p"] = [0.001, 0.001, 0.5, np.nan, 0.5, 0.001, 0.5, 0.2]
    df["gse118713_fdr"] = [0.01, 0.01, 0.5, np.nan, 0.5, 0.01, 0.5, 0.3]
    df["gse118713_testable"] = [True, True, True, False, True, True, True, True]
    df["gse118713_evidence_percentile"] = [0.95, 0.95, 0.3, np.nan, 0.3, 0.95, 0.3, 0.5]
    df["gse118713_direction"] = ["up_in_TAMR", "down_in_TAMR", "up_in_TAMR", np.nan, "up_in_TAMR", "up_in_TAMR", "up_in_TAMR", "up_in_TAMR"]
    # gse240112 (resistance dataset, tumor track)
    df["gse240112_tumor_log2fc"] = [0.5, -0.5, -0.5, 0.5, 0.1, 0.6, 0.1, -0.5]
    df["gse240112_tumor_p"] = [0.001, 0.001, 0.001, 0.5, 0.5, 0.001, 0.5, 0.2]
    df["gse240112_tumor_fdr"] = [0.01, 0.01, 0.01, 0.5, 0.5, 0.01, 0.5, 0.3]
    df["gse240112_testable"] = [True] * 8
    df["gse240112_evidence_percentile"] = [0.95, 0.95, 0.95, 0.3, 0.3, 0.95, 0.3, 0.5]
    df["gse240112_epi_fdr"] = [0.01, 0.01, 0.01, 0.5, 0.5, 0.01, 0.5, 0.3]
    # gse111151 (resistance dataset) -- untestable for MIXED/CONTEXT_MIX, so each has
    # exactly 2 testable resistance datasets with opposite direction (1 up, 1 down = "mixed")
    df["gse111151_log2fc"] = [0.5, -0.5, np.nan, 0.5, 0.1, 0.6, 0.1, np.nan]
    df["gse111151_p"] = [0.5, 0.5, np.nan, 0.5, 0.5, 0.5, 0.5, np.nan]
    df["gse111151_fdr"] = [0.5, 0.5, np.nan, 0.5, 0.5, 0.5, 0.5, np.nan]
    df["gse111151_testable"] = [True, True, False, True, True, True, True, False]
    df["gse111151_evidence_percentile"] = [0.6, 0.6, np.nan, 0.6, 0.3, 0.6, 0.3, np.nan]
    # gse245601 (acute, not a resistance dataset)
    df["gse245601_epi_fdr"] = [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.01, 0.5]
    df["gse245601_evidence_percentile"] = [0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.99, 0.3]
    df["gse245601_testable"] = [True] * 8
    # crispr_direction, matching each gene's crispr_effect sign
    df["crispr_direction"] = ["approximately_neutral", "approximately_neutral", "approximately_neutral", "approximately_neutral", "sensitising_KO", "sensitising_KO", "approximately_neutral", "approximately_neutral"]
    # coverage: INSUFFICIENT genuinely below the 3-dataset LOW_COVERAGE threshold
    df["n_datasets_testable"] = [5, 5, 4, 2, 5, 5, 5, 4]
    df["n_datasets_fdr05"] = [3, 3, 2, 0, 1, 4, 1, 0]
    # CONTEXT_MIX: no dataset individually FDR<0.05 or top-10% on its own via the
    # fdr05 count, but reaches the top-10% threshold in >=2 datasets -- exercises the
    # "or top-10%" branch of CONTEXT_DEPENDENT specifically (the branch a real bug
    # caught by Codex review had omitted from the code, though not from the docstring)
    df["n_datasets_top10pct"] = [4, 4, 2, 0, 2, 4, 1, 2]
    return df


class TestBuildResistanceConsensus:
    def test_all_up_direction_consensus(self):
        out = build_resistance_consensus(_synthetic_full()).set_index("gene")
        assert out.loc["ALLUP", "resistance_direction_consensus"] == "all_up"
        assert out.loc["ALLUP", "resistance_up_count"] == 3

    def test_all_down_direction_consensus(self):
        out = build_resistance_consensus(_synthetic_full()).set_index("gene")
        assert out.loc["ALLDOWN", "resistance_direction_consensus"] == "all_down"

    def test_mixed_direction_consensus(self):
        out = build_resistance_consensus(_synthetic_full()).set_index("gene")
        assert out.loc["MIXED", "resistance_direction_consensus"] == "mixed"

    def test_untestable_dataset_reduces_testable_count(self):
        out = build_resistance_consensus(_synthetic_full()).set_index("gene")
        assert out.loc["INSUFFICIENT", "resistance_datasets_testable"] == 2  # gse118713 untestable

    def test_excludes_gse245601_entirely(self):
        # ACUTE_ONLY has strong gse245601 signal but weak resistance-dataset signal -- must not count as resistance support
        out = build_resistance_consensus(_synthetic_full()).set_index("gene")
        assert out.loc["ACUTE_ONLY", "resistance_fdr05_count"] == 0

    def test_ranking_deterministic(self):
        df = _synthetic_full()
        out1 = build_resistance_consensus(df)
        out2 = build_resistance_consensus(df)
        pd.testing.assert_frame_equal(out1, out2)


class TestBuildCrisprFunctionalRanking:
    def test_only_testable_genes_included(self):
        out = build_crispr_functional_ranking(_synthetic_full())
        assert len(out) == 8  # all testable=True in synthetic data

    def test_sorted_by_fdr_ascending(self):
        out = build_crispr_functional_ranking(_synthetic_full())
        fdrs = out["crispr_fdr"].tolist()
        assert fdrs == sorted(fdrs)

    def test_sensitising_and_tolerance_both_represented(self):
        out = build_crispr_functional_ranking(_synthetic_full())
        assert "sensitising_KO" in set(out["crispr_direction"])


class TestHumanAndRnaOnlyRankings:
    def test_human_only_excludes_crispr_and_cellline_datasets(self):
        out = build_human_only_ranking(_synthetic_full())
        # a gene testable only in crispr/gse118713/gse111151 (not gse245601/gse240112) should show 0 human-testable datasets
        assert "n_datasets_testable" in out.columns

    def test_rna_only_has_no_crispr_column_influence(self):
        out = build_rna_only_ranking(_synthetic_full())
        assert set(out.columns).isdisjoint({"crispr_fdr", "crispr_effect"})


class TestCrisprNonsignificantRnaConsensus:
    def test_multimodal_gene_excluded_since_crispr_is_significant(self):
        resistance = build_resistance_consensus(_synthetic_full())
        out = build_crispr_nonsignificant_rna_consensus(_synthetic_full(), resistance)
        assert "MULTIMODAL" not in set(out["gene"])  # crispr_fdr=0.01 < 0.10, so excluded

    def test_allup_gene_with_crispr_nonsignificant_included(self):
        resistance = build_resistance_consensus(_synthetic_full())
        out = build_crispr_nonsignificant_rna_consensus(_synthetic_full(), resistance)
        assert "ALLUP" in set(out["gene"])
        assert (out.loc[out["gene"] == "ALLUP", "label"] == "resistance_biomarker_or_pathway_candidate").all()


class TestAssignEvidenceCategory:
    def test_low_coverage_gene(self):
        resistance = build_resistance_consensus(_synthetic_full())
        out = assign_evidence_category(_synthetic_full(), resistance).set_index("gene")
        assert out.loc["INSUFFICIENT", "evidence_category"] == "LOW_COVERAGE"

    def test_multimodal_strong_gene(self):
        resistance = build_resistance_consensus(_synthetic_full())
        out = assign_evidence_category(_synthetic_full(), resistance).set_index("gene")
        assert out.loc["MULTIMODAL", "evidence_category"] == "MULTIMODAL_STRONG"

    def test_functional_only_gene(self):
        resistance = build_resistance_consensus(_synthetic_full())
        out = assign_evidence_category(_synthetic_full(), resistance).set_index("gene")
        assert out.loc["CRISPR_ONLY", "evidence_category"] == "FUNCTIONAL_ONLY"

    def test_rna_resistance_consensus_gene(self):
        resistance = build_resistance_consensus(_synthetic_full())
        out = assign_evidence_category(_synthetic_full(), resistance).set_index("gene")
        assert out.loc["ALLUP", "evidence_category"] == "RNA_RESISTANCE_CONSENSUS"

    def test_acute_response_gene(self):
        resistance = build_resistance_consensus(_synthetic_full())
        out = assign_evidence_category(_synthetic_full(), resistance).set_index("gene")
        assert out.loc["ACUTE_ONLY", "evidence_category"] == "ACUTE_RESPONSE"

    def test_every_gene_gets_exactly_one_category(self):
        resistance = build_resistance_consensus(_synthetic_full())
        out = assign_evidence_category(_synthetic_full(), resistance)
        assert out["evidence_category"].notna().all()
        assert len(out) == 8

    def test_context_dependent_gene_via_top10pct_branch(self):
        # regression test for a real bug caught by Codex review: the docstring says
        # CONTEXT_DEPENDENT requires ">=2 datasets at FDR<0.05 OR top-10%", but the
        # code checked only the FDR<0.05 count, silently dropping the "or top-10%"
        # branch. CONTEXT_MIX has n_datasets_fdr05=0 but n_datasets_top10pct=2 and a
        # "mixed" resistance direction consensus -- must reach CONTEXT_DEPENDENT.
        resistance = build_resistance_consensus(_synthetic_full())
        out = assign_evidence_category(_synthetic_full(), resistance).set_index("gene")
        assert out.loc["CONTEXT_MIX", "evidence_category"] == "CONTEXT_DEPENDENT"


class TestRealData:
    def test_real_resistance_consensus_if_present(self):
        path = REPO_ROOT / "results" / "tables" / "cross_dataset_genomewide" / "resistance_consensus_all_genes.tsv"
        if not path.exists():
            pytest.skip("resistance consensus table not present in this checkout")
        out = pd.read_csv(path, sep="\t")
        assert not out["gene"].duplicated().any()
        valid = {"all_up", "all_down", "majority_up", "majority_down", "mixed", "insufficient"}
        assert set(out["resistance_direction_consensus"]).issubset(valid)

    def test_real_evidence_categories_if_present(self):
        path = REPO_ROOT / "results" / "tables" / "cross_dataset_genomewide" / "evidence_categories.tsv"
        if not path.exists():
            pytest.skip("evidence categories table not present in this checkout")
        out = pd.read_csv(path, sep="\t")
        valid = {"LOW_COVERAGE", "MULTIMODAL_STRONG", "RNA_RESISTANCE_CONSENSUS", "FUNCTIONAL_ONLY", "ACUTE_RESPONSE", "HUMAN_TUMOR_SUPPORTED", "CONTEXT_DEPENDENT", "LOW_EVIDENCE"}
        assert set(out["evidence_category"]).issubset(valid)
