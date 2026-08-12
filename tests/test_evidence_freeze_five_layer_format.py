import pandas as pd

from src.evidence_freeze_five_layer_format import (
    FULL_RNA_DATASET_ORDER,
    RESISTANCE_DATASET_ORDER,
    acute_direction,
    arrow_for,
    full_rna_pattern_4,
    resistance_direction_consistency,
    resistance_fdr05_count,
    resistance_pattern_3,
)


def _row(**kwargs):
    base = {
        "gse118713_log2fc": float("nan"), "gse118713_fdr": float("nan"), "gse118713_p": float("nan"),
        "gse240112_log2fc": float("nan"), "gse240112_fdr": float("nan"), "gse240112_p": float("nan"),
        "gse111151_log2fc": float("nan"), "gse111151_fdr": float("nan"), "gse111151_p": float("nan"),
        "gse245601_epi_log2fc": float("nan"), "gse245601_epi_fdr": float("nan"),
    }
    base.update(kwargs)
    return pd.Series(base)


class TestCanonicalOrder:
    def test_resistance_order_is_gse118713_gse240112_gse111151(self):
        assert RESISTANCE_DATASET_ORDER == ["gse118713", "gse240112", "gse111151"]

    def test_full_order_appends_gse245601_last(self):
        assert FULL_RNA_DATASET_ORDER == ["gse118713", "gse240112", "gse111151", "gse245601"]


class TestArrowFor:
    def test_positive_significant_gets_up_star(self):
        assert arrow_for(1.5, 0.01) == "↑*"

    def test_positive_nonsignificant_gets_up_no_star(self):
        assert arrow_for(1.5, 0.5) == "↑"

    def test_negative_gets_down(self):
        assert arrow_for(-1.5, 0.01) == "↓*"

    def test_nan_gets_na_not_zero_or_flat(self):
        assert arrow_for(float("nan"), float("nan")) == "NA"

    def test_zero_gets_flat_arrow(self):
        assert arrow_for(0.0, 0.9) == "→"


class TestResistancePatternExcludesGse245601:
    def test_pattern_has_exactly_three_arrows(self):
        row = _row(gse118713_log2fc=1, gse118713_fdr=0.01, gse240112_log2fc=1, gse240112_fdr=0.01, gse111151_log2fc=1, gse111151_fdr=0.01)
        pattern = resistance_pattern_3(row)
        assert pattern.count("|") == 2  # 3 arrows joined by 2 pipes
        assert "gse245601" not in pattern.lower()

    def test_gse245601_extreme_value_does_not_leak_into_resistance_pattern(self):
        row = _row(gse118713_log2fc=1, gse118713_fdr=0.01, gse240112_log2fc=1, gse240112_fdr=0.01, gse111151_log2fc=-1, gse111151_fdr=0.01, gse245601_epi_log2fc=999, gse245601_epi_fdr=1e-30)
        pattern = resistance_pattern_3(row)
        assert pattern == "↑* | ↑* | ↓*"


class TestFullPatternHasFourSlotsWithDivider:
    def test_full_pattern_contains_divider_and_four_arrows(self):
        row = _row(gse118713_log2fc=1, gse118713_fdr=0.01, gse240112_log2fc=1, gse240112_fdr=0.5, gse111151_log2fc=-1, gse111151_fdr=0.5, gse245601_epi_log2fc=-1, gse245601_epi_fdr=0.01)
        full = full_rna_pattern_4(row)
        assert full == "↑* | ↑ | ↓ || ↓*"
        assert "||" in full

    def test_acute_direction_matches_last_slot_of_full_pattern(self):
        row = _row(gse245601_epi_log2fc=2, gse245601_epi_fdr=0.001)
        assert full_rna_pattern_4(row).split("||")[1].strip() == acute_direction(row)


class TestResistanceFdr05CountExcludesGse245601:
    def test_count_ignores_gse245601_even_if_significant(self):
        row = _row(gse118713_fdr=0.01, gse240112_fdr=0.5, gse111151_fdr=0.5, gse245601_epi_fdr=1e-20)
        assert resistance_fdr05_count(row) == 1


class TestResistanceDirectionConsistency:
    def test_all_up(self):
        row = _row(gse118713_log2fc=1, gse240112_log2fc=2, gse111151_log2fc=0.5)
        assert resistance_direction_consistency(row) == "all_up"

    def test_majority_up_when_two_of_three_agree(self):
        row = _row(gse118713_log2fc=1, gse240112_log2fc=-2, gse111151_log2fc=0.5)
        assert resistance_direction_consistency(row) == "majority_up"

    def test_mixed_when_exactly_tied(self):
        """A tie is only possible when one of the three is untestable,
        leaving an even 1-vs-1 split among the remaining two."""
        row = _row(gse118713_log2fc=1, gse240112_log2fc=-1, gse111151_log2fc=float("nan"))
        assert resistance_direction_consistency(row) == "mixed"

    def test_gse245601_direction_never_affects_this_calculation(self):
        row_a = _row(gse118713_log2fc=1, gse240112_log2fc=1, gse111151_log2fc=1, gse245601_epi_log2fc=-999)
        row_b = _row(gse118713_log2fc=1, gse240112_log2fc=1, gse111151_log2fc=1, gse245601_epi_log2fc=999)
        assert resistance_direction_consistency(row_a) == resistance_direction_consistency(row_b) == "all_up"

    def test_no_testable_resistance_dataset_is_insufficient(self):
        row = _row()
        assert resistance_direction_consistency(row) == "insufficient"
