from pathlib import Path

import pandas as pd
import pytest

from src.evidence_freeze_visualization import plot_therapeutic_shortlist_head_to_head

REPO_ROOT = Path(__file__).parent.parent


def _minimal_df(genes_and_directions: dict[str, str]) -> pd.DataFrame:
    rows = []
    for gene, direction in genes_and_directions.items():
        rows.append(
            {
                "gene": gene, "crispr_direction": direction, "crispr_fdr": 0.04,
                "resistance_fdr05_count": 1, "human_tumor_support": "significant",
                "resistance_direction_consistency": "all_up", "ranking_stability": "DATASET_DEPENDENT",
            }
        )
    return pd.DataFrame(rows)


class TestTherapeuticHeadToHeadGuardsDirection:
    def test_raises_if_a_tolerance_gene_is_included(self, tmp_path):
        df = _minimal_df({"GOOD": "sensitising_KO", "BAD": "tolerance_associated_KO"})
        with pytest.raises(ValueError):
            plot_therapeutic_shortlist_head_to_head(df, ["GOOD", "BAD"], tmp_path / "out.png")

    def test_succeeds_when_all_sensitising(self, tmp_path):
        df = _minimal_df({"GOOD1": "sensitising_KO", "GOOD2": "sensitising_KO"})
        out_path = tmp_path / "out.png"
        plot_therapeutic_shortlist_head_to_head(df, ["GOOD1", "GOOD2"], out_path)
        assert out_path.exists()


class TestRealFiguresExist:
    def test_four_final_review_figures_present(self):
        final_review = REPO_ROOT / "results" / "figures" / "evidence_freeze" / "final_review"
        if not final_review.exists():
            pytest.skip("evidence-freeze figures not generated in this environment")
        expected = {"01_four_rna_direction_matrix.png", "02_five_layer_support_matrix.png", "03_therapeutic_shortlist_head_to_head.png", "04_frozen_shortlist_summary.png"}
        actual = {p.name for p in final_review.glob("*.png")}
        assert expected.issubset(actual)
