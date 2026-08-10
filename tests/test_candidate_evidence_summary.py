import copy
import math
from pathlib import Path

import pandas as pd
import pytest
import yaml

import src.candidate_evidence_summary as candidate_evidence_summary
from src.candidate_evidence_summary import (
    BENCHMARK_GENE_SYMBOL,
    BENCHMARK_LABEL,
    DIRECTION_INDETERMINATE,
    DIRECTION_SENSITISING,
    DIRECTION_TOLERANCE,
    EVIDENCE_CLASS_NO_SIGNIFICANT_RNA,
    EVIDENCE_CLASS_NOT_APPLICABLE,
    EVIDENCE_CLASS_PRIMARY,
    EVIDENCE_CLASS_RNA_UNAVAILABLE,
    EVIDENCE_CLASS_SECONDARY,
    RNA_DIRECTION_ELEVATED,
    RNA_DIRECTION_NOT_APPLICABLE,
    RNA_DIRECTION_REDUCED,
    EvidenceSummaryConfig,
    build_evidence_summary,
    build_paics_benchmark,
    build_plot_inputs,
    classify_crispr_direction,
    classify_evidence_class,
    classify_rna_direction,
    describe_rna_pattern,
    load_master_table,
    run_candidate_evidence_summary,
    select_primary_shortlist,
    select_secondary_candidates,
    split_by_direction,
)

REPO_ROOT = Path(__file__).parent.parent
CONFIG_PATH = REPO_ROOT / "config" / "config.yaml"
RNA_FDR = 0.05


def _row(
    gene_symbol="G1",
    crispr_effect_size=-1.0,
    crispr_fdr=0.01,
    mcf7_log2fc=1.0,
    mcf7_fdr=0.01,
    fasr_log2fc=1.0,
    fasr_fdr=0.01,
    de_na_reason=pd.NA,
):
    return pd.Series(
        {
            "gene_symbol": gene_symbol,
            "crispr_effect_size": crispr_effect_size,
            "crispr_fdr": crispr_fdr,
            "mapping_status": "unique_filtered_out" if not pd.isna(de_na_reason) else "unique_filtered",
            "tamr_vs_mcf7_log2fc": mcf7_log2fc,
            "tamr_vs_mcf7_fdr": mcf7_fdr,
            "tamr_vs_fasr_log2fc": fasr_log2fc,
            "tamr_vs_fasr_fdr": fasr_fdr,
            "mcf7_baseline_log2_tpm_plus1": 5.0,
            "gse118713_de_na_reason": de_na_reason,
        }
    )


def _row_with_direction(**kwargs):
    row = _row(**kwargs)
    row["crispr_direction"] = classify_crispr_direction(row["crispr_effect_size"])
    return row


# --- CRISPR direction: verified sign convention ---------------------------


class TestClassifyCrisprDirection:
    def test_negative_is_sensitising(self):
        assert classify_crispr_direction(-1.5) == DIRECTION_SENSITISING

    def test_positive_is_tolerance_associated(self):
        assert classify_crispr_direction(2.0) == DIRECTION_TOLERANCE

    def test_zero_is_indeterminate(self):
        assert classify_crispr_direction(0.0) == DIRECTION_INDETERMINATE


# --- evidence-class hierarchy: TAMR_vs_MCF7 outranks TAMR_vs_FASR ---------


class TestClassifyEvidenceClass:
    def test_hierarchy_is_deterministic_same_inputs_same_output(self):
        row = _row_with_direction(mcf7_fdr=0.01, fasr_fdr=0.01)
        results = {classify_evidence_class(row, RNA_FDR) for _ in range(5)}
        assert results == {EVIDENCE_CLASS_PRIMARY}

    def test_significant_in_both_is_primary_not_secondary(self):
        # TAMR_vs_MCF7 has priority: significance in both contrasts still
        # resolves to PRIMARY, never SECONDARY.
        row = _row_with_direction(mcf7_fdr=0.01, fasr_fdr=0.01)
        assert classify_evidence_class(row, RNA_FDR) == EVIDENCE_CLASS_PRIMARY

    def test_significant_in_mcf7_only_is_primary(self):
        row = _row_with_direction(mcf7_fdr=0.01, fasr_fdr=0.5)
        assert classify_evidence_class(row, RNA_FDR) == EVIDENCE_CLASS_PRIMARY

    def test_significant_in_fasr_only_is_secondary(self):
        row = _row_with_direction(mcf7_fdr=0.5, fasr_fdr=0.01)
        assert classify_evidence_class(row, RNA_FDR) == EVIDENCE_CLASS_SECONDARY

    def test_significant_in_neither_is_no_significant_rna_support(self):
        row = _row_with_direction(mcf7_fdr=0.5, fasr_fdr=0.5)
        assert classify_evidence_class(row, RNA_FDR) == EVIDENCE_CLASS_NO_SIGNIFICANT_RNA

    def test_de_unavailable_is_rna_unavailable_regardless_of_fdr_values(self):
        row = _row_with_direction(mcf7_fdr=0.01, fasr_fdr=0.01, de_na_reason="filtered_out_no_de_fit")
        assert classify_evidence_class(row, RNA_FDR) == EVIDENCE_CLASS_RNA_UNAVAILABLE

    def test_boundary_fdr_exactly_at_threshold_is_not_significant(self):
        row = _row_with_direction(mcf7_fdr=RNA_FDR, fasr_fdr=RNA_FDR)
        assert classify_evidence_class(row, RNA_FDR) == EVIDENCE_CLASS_NO_SIGNIFICANT_RNA

    def test_positive_effect_is_not_applicable_even_if_both_significant(self):
        row = _row_with_direction(crispr_effect_size=2.0, mcf7_fdr=0.01, fasr_fdr=0.01)
        assert classify_evidence_class(row, RNA_FDR) == EVIDENCE_CLASS_NOT_APPLICABLE

    def test_zero_effect_is_not_applicable(self):
        row = _row_with_direction(crispr_effect_size=0.0, mcf7_fdr=0.01, fasr_fdr=0.01)
        assert classify_evidence_class(row, RNA_FDR) == EVIDENCE_CLASS_NOT_APPLICABLE


# --- RNA direction annotation ----------------------------------------------


class TestClassifyRnaDirection:
    def test_positive_log2fc_is_elevated(self):
        assert classify_rna_direction(0.5) == RNA_DIRECTION_ELEVATED

    def test_negative_log2fc_is_reduced(self):
        assert classify_rna_direction(-0.5) == RNA_DIRECTION_REDUCED


# --- RNA pattern description: neutral, no causal/specificity language -----


class TestDescribeRnaPattern:
    def test_de_unavailable(self):
        row = _row(de_na_reason="filtered_out_no_de_fit")
        desc = describe_rna_pattern(row, RNA_FDR)
        assert "No GSE118713" in desc
        assert "filter" in desc.lower()

    def test_both_significant_mentions_both_contrasts(self):
        row = _row(mcf7_fdr=0.01, mcf7_log2fc=1.0, fasr_fdr=0.01, fasr_log2fc=-1.0)
        desc = describe_rna_pattern(row, RNA_FDR)
        assert "both MCF7 and FASR" in desc
        assert "up vs MCF7" in desc
        assert "down vs FASR" in desc

    def test_never_calls_fasr_significance_tamoxifen_specific(self):
        row = _row(mcf7_fdr=0.5, fasr_fdr=0.01)
        desc = describe_rna_pattern(row, RNA_FDR)
        assert "specific" not in desc.lower()
        assert "tamoxifen" not in desc.lower()

    def test_never_treats_fasr_nonsignificance_as_specificity_proof(self):
        row = _row(mcf7_fdr=0.01, fasr_fdr=0.5)
        desc = describe_rna_pattern(row, RNA_FDR)
        assert "specific" not in desc.lower()

    def test_no_causal_language(self):
        for row in (
            _row(mcf7_fdr=0.01, fasr_fdr=0.01),
            _row(mcf7_fdr=0.01, fasr_fdr=0.5),
            _row(mcf7_fdr=0.5, fasr_fdr=0.5),
            _row(de_na_reason="filtered_out_no_de_fit"),
        ):
            desc = describe_rna_pattern(row, RNA_FDR).lower()
            for causal_word in ("causes", "drives", "because", "due to", "restores sensitivity", "reverses resistance"):
                assert causal_word not in desc


# --- build_evidence_summary: no drops, no scoring, hierarchy-aware --------


class TestBuildEvidenceSummary:
    def test_every_row_preserved_and_classified(self):
        master_df = pd.DataFrame(
            [
                _row("PRIMARY_GENE", crispr_effect_size=-1.0, mcf7_fdr=0.01, fasr_fdr=0.01),
                _row("SECONDARY_GENE", crispr_effect_size=-1.0, mcf7_fdr=0.5, fasr_fdr=0.01),
                _row("UNAVAILABLE_GENE", crispr_effect_size=-1.0, de_na_reason="filtered_out_no_de_fit"),
                _row("TOLERANCE_GENE", crispr_effect_size=2.0, mcf7_fdr=0.01, fasr_fdr=0.01),
            ]
        )
        out = build_evidence_summary(master_df, RNA_FDR)
        assert len(out) == 4
        indexed = out.set_index("gene_symbol")
        assert indexed.loc["PRIMARY_GENE", "evidence_class"] == EVIDENCE_CLASS_PRIMARY
        assert indexed.loc["PRIMARY_GENE", "primary_rna_direction"] == RNA_DIRECTION_ELEVATED
        assert indexed.loc["SECONDARY_GENE", "evidence_class"] == EVIDENCE_CLASS_SECONDARY
        assert indexed.loc["SECONDARY_GENE", "primary_rna_direction"] == RNA_DIRECTION_NOT_APPLICABLE
        assert indexed.loc["UNAVAILABLE_GENE", "evidence_class"] == EVIDENCE_CLASS_RNA_UNAVAILABLE
        assert indexed.loc["TOLERANCE_GENE", "crispr_direction"] == DIRECTION_TOLERANCE
        assert indexed.loc["TOLERANCE_GENE", "evidence_class"] == EVIDENCE_CLASS_NOT_APPLICABLE

    def test_no_score_or_rank_column_present(self):
        master_df = pd.DataFrame([_row("A")])
        out = build_evidence_summary(master_df, RNA_FDR)
        forbidden = {"score", "composite_score", "rank", "weighted_score"}
        assert forbidden.isdisjoint(set(out.columns))


# --- split_by_direction: two clearly separate biological groups -----------


class TestSplitByDirection:
    def test_splits_into_sensitising_and_tolerance(self):
        master_df = pd.DataFrame(
            [
                _row("SENS", crispr_effect_size=-1.0),
                _row("TOL", crispr_effect_size=2.0),
            ]
        )
        evidence_df = build_evidence_summary(master_df, RNA_FDR)
        sensitising_df, tolerance_df = split_by_direction(evidence_df)
        assert list(sensitising_df["gene_symbol"]) == ["SENS"]
        assert list(tolerance_df["gene_symbol"]) == ["TOL"]

    def test_tolerance_group_carries_not_inhibition_candidate_reason(self):
        master_df = pd.DataFrame([_row("TOL", crispr_effect_size=2.0)])
        evidence_df = build_evidence_summary(master_df, RNA_FDR)
        _, tolerance_df = split_by_direction(evidence_df)
        reason = tolerance_df.iloc[0]["not_inhibition_candidate_reason"]
        assert "not a candidate for inhibition" in reason
        assert "tolerance" in reason.lower()

    def test_positive_effect_genes_never_leak_into_sensitising_group(self):
        master_df = pd.DataFrame(
            [
                _row("TOL1", crispr_effect_size=2.0, mcf7_fdr=0.01, fasr_fdr=0.01),
                _row("TOL2", crispr_effect_size=1.5),
                _row("SENS", crispr_effect_size=-1.0),
            ]
        )
        evidence_df = build_evidence_summary(master_df, RNA_FDR)
        sensitising_df, _ = split_by_direction(evidence_df)
        assert set(sensitising_df["gene_symbol"]) == {"SENS"}

    def test_indeterminate_gene_excluded_from_both_groups(self):
        master_df = pd.DataFrame(
            [
                _row("ZERO", crispr_effect_size=0.0),
                _row("SENS", crispr_effect_size=-1.0),
            ]
        )
        evidence_df = build_evidence_summary(master_df, RNA_FDR)
        sensitising_df, tolerance_df = split_by_direction(evidence_df)
        assert "ZERO" not in set(sensitising_df["gene_symbol"])
        assert "ZERO" not in set(tolerance_df["gene_symbol"])

    def test_sensitising_group_ordered_by_evidence_class_priority(self):
        master_df = pd.DataFrame(
            [
                _row("SECONDARY_GENE", crispr_effect_size=-1.0, mcf7_fdr=0.5, fasr_fdr=0.01),
                _row("PRIMARY_GENE", crispr_effect_size=-1.0, mcf7_fdr=0.01, fasr_fdr=0.5),
                _row("NO_RNA_GENE", crispr_effect_size=-1.0, mcf7_fdr=0.5, fasr_fdr=0.5),
                _row("UNAVAILABLE_GENE", crispr_effect_size=-1.0, de_na_reason="filtered_out_no_de_fit"),
            ]
        )
        evidence_df = build_evidence_summary(master_df, RNA_FDR)
        sensitising_df, _ = split_by_direction(evidence_df)
        assert list(sensitising_df["gene_symbol"]) == [
            "PRIMARY_GENE",
            "SECONDARY_GENE",
            "NO_RNA_GENE",
            "UNAVAILABLE_GENE",
        ]


# --- select_primary_shortlist: PRIMARY_RESISTANCE_SUPPORT only, no fallback


class TestSelectPrimaryShortlist:
    def _sensitising_df(self, rows):
        master_df = pd.DataFrame(rows)
        evidence_df = build_evidence_summary(master_df, RNA_FDR)
        sensitising_df, _ = split_by_direction(evidence_df)
        return sensitising_df

    def test_only_primary_class_selected(self):
        sensitising_df = self._sensitising_df(
            [
                _row("PRIMARY_GENE", crispr_effect_size=-1.0, mcf7_fdr=0.01, fasr_fdr=0.5),
                _row("SECONDARY_GENE", crispr_effect_size=-1.0, mcf7_fdr=0.5, fasr_fdr=0.01),
                _row("NO_RNA_GENE", crispr_effect_size=-1.0, mcf7_fdr=0.5, fasr_fdr=0.5),
            ]
        )
        shortlist = select_primary_shortlist(sensitising_df)
        assert list(shortlist["gene_symbol"]) == ["PRIMARY_GENE"]

    def test_no_fallback_to_secondary_when_primary_empty(self):
        sensitising_df = self._sensitising_df(
            [
                _row("SECONDARY_GENE", crispr_effect_size=-1.0, mcf7_fdr=0.5, fasr_fdr=0.01),
                _row("NO_RNA_GENE", crispr_effect_size=-1.0, mcf7_fdr=0.5, fasr_fdr=0.5),
            ]
        )
        shortlist = select_primary_shortlist(sensitising_df)
        assert len(shortlist) == 0

    def test_single_gene_not_padded(self):
        sensitising_df = self._sensitising_df(
            [_row("ONLY_PRIMARY", crispr_effect_size=-1.0, mcf7_fdr=0.01, fasr_fdr=0.5)]
        )
        shortlist = select_primary_shortlist(sensitising_df)
        assert len(shortlist) == 1
        assert shortlist.iloc[0]["gene_symbol"] == "ONLY_PRIMARY"

    def test_reason_uses_required_wording(self):
        sensitising_df = self._sensitising_df(
            [_row("A_GENE", crispr_effect_size=-1.0, crispr_fdr=0.02, mcf7_fdr=0.01, fasr_fdr=0.03)]
        )
        shortlist = select_primary_shortlist(sensitising_df)
        reason = shortlist.iloc[0]["shortlist_reason"]
        assert "resistance-associated expression support" in reason
        assert "not a proven resistance mechanism" in reason
        assert "restored tamoxifen" in reason


class TestSelectSecondaryCandidates:
    def test_only_secondary_class_selected(self):
        master_df = pd.DataFrame(
            [
                _row("PRIMARY_GENE", crispr_effect_size=-1.0, mcf7_fdr=0.01, fasr_fdr=0.5),
                _row("SECONDARY_GENE", crispr_effect_size=-1.0, mcf7_fdr=0.5, fasr_fdr=0.01),
                _row("NO_RNA_GENE", crispr_effect_size=-1.0, mcf7_fdr=0.5, fasr_fdr=0.5),
            ]
        )
        evidence_df = build_evidence_summary(master_df, RNA_FDR)
        sensitising_df, _ = split_by_direction(evidence_df)
        secondary_df = select_secondary_candidates(sensitising_df)
        assert list(secondary_df["gene_symbol"]) == ["SECONDARY_GENE"]

    def test_sorted_alphabetically_not_by_fdr(self):
        master_df = pd.DataFrame(
            [
                _row("ZZZ", crispr_effect_size=-1.0, crispr_fdr=0.001, mcf7_fdr=0.5, fasr_fdr=0.01),
                _row("AAA", crispr_effect_size=-1.0, crispr_fdr=0.09, mcf7_fdr=0.5, fasr_fdr=0.01),
            ]
        )
        evidence_df = build_evidence_summary(master_df, RNA_FDR)
        sensitising_df, _ = split_by_direction(evidence_df)
        secondary_df = select_secondary_candidates(sensitising_df)
        assert list(secondary_df["gene_symbol"]) == ["AAA", "ZZZ"]

    def test_note_mentions_no_significant_mcf7_and_contextual_fasr(self):
        master_df = pd.DataFrame([_row("A", crispr_effect_size=-1.0, mcf7_fdr=0.5, fasr_fdr=0.01)])
        evidence_df = build_evidence_summary(master_df, RNA_FDR)
        sensitising_df, _ = split_by_direction(evidence_df)
        secondary_df = select_secondary_candidates(sensitising_df)
        note = secondary_df.iloc[0]["secondary_context_note"]
        assert "no significant TAMR_vs_MCF7" in note
        assert "Contextual" in note


# --- PAICS benchmark: strictly separate from candidate discovery ----------


class TestBuildPaicsBenchmark:
    def _labels_df(self, effect_size=-0.296841, fdr=0.852723, n_guides=4):
        return pd.DataFrame(
            [
                {"gene": "OTHER", "n_guides": 4, "effect_size": -1.0, "se": 0.1, "p_value": 0.01, "fdr": 0.01},
                {
                    "gene": BENCHMARK_GENE_SYMBOL,
                    "n_guides": n_guides,
                    "effect_size": effect_size,
                    "se": 0.3,
                    "p_value": 0.34,
                    "fdr": fdr,
                },
            ]
        )

    def _de_df(self):
        return pd.DataFrame(
            [
                {
                    "gene_id": "ENSG00000128050",
                    "gene_symbol": BENCHMARK_GENE_SYMBOL,
                    "log2fc": -0.323554,
                    "se": 0.1,
                    "moderated_t": -2.5,
                    "p_value": 0.028,
                    "fdr": 0.0794,
                    "ave_expr": 7.4,
                    "contrast": "TAMR_vs_MCF7",
                    "direction": "down",
                },
                {
                    "gene_id": "ENSG00000128050",
                    "gene_symbol": BENCHMARK_GENE_SYMBOL,
                    "log2fc": 0.714911,
                    "se": 0.1,
                    "moderated_t": 5.6,
                    "p_value": 0.0002,
                    "fdr": 0.000831,
                    "ave_expr": 7.4,
                    "contrast": "TAMR_vs_FASR",
                    "direction": "up",
                },
            ]
        )

    def _master_df_without_paics(self):
        return pd.DataFrame([_row("SOME_HIT", crispr_effect_size=-1.0)])

    def test_paics_benchmark_row_built_correctly(self):
        paics_df = build_paics_benchmark(self._labels_df(), self._de_df(), self._master_df_without_paics())
        assert len(paics_df) == 1
        row = paics_df.iloc[0]
        assert row["gene_symbol"] == BENCHMARK_GENE_SYMBOL
        assert row["crispr_direction"] == DIRECTION_SENSITISING
        assert row["benchmark_label"] == BENCHMARK_LABEL
        assert row["tamr_vs_mcf7_log2fc"] == pytest.approx(-0.323554)
        assert row["tamr_vs_fasr_log2fc"] == pytest.approx(0.714911)

    def test_paics_cannot_enter_the_28_hit_shortlist(self):
        # If PAICS is ever found in the frozen master table, the benchmark
        # function must refuse to run rather than silently merging it in.
        master_df_with_paics = pd.DataFrame(
            [_row("SOME_HIT", crispr_effect_size=-1.0), _row(BENCHMARK_GENE_SYMBOL, crispr_effect_size=-1.0)]
        )
        with pytest.raises(ValueError, match="strictly separate"):
            build_paics_benchmark(self._labels_df(), self._de_df(), master_df_with_paics)

    def test_paics_not_present_in_real_master_table(self):
        # Sanity check against the real, frozen 28-hit table.
        with open(CONFIG_PATH) as f:
            config = yaml.safe_load(f)
        cfg = EvidenceSummaryConfig.from_config(config)
        master_df = load_master_table(cfg)
        assert BENCHMARK_GENE_SYMBOL not in set(master_df["gene_symbol"])

    def test_non_finite_crispr_effect_size_raises(self):
        with pytest.raises(ValueError, match="crispr_effect_size is not finite"):
            build_paics_benchmark(
                self._labels_df(effect_size=math.nan), self._de_df(), self._master_df_without_paics()
            )

    def test_non_finite_crispr_fdr_raises(self):
        with pytest.raises(ValueError, match="crispr_fdr is not finite"):
            build_paics_benchmark(self._labels_df(fdr=math.inf), self._de_df(), self._master_df_without_paics())

    def test_crispr_fdr_out_of_range_raises(self):
        with pytest.raises(ValueError, match=r"crispr_fdr is outside \[0, 1\]"):
            build_paics_benchmark(self._labels_df(fdr=1.2), self._de_df(), self._master_df_without_paics())

    def test_non_finite_de_value_raises(self):
        de_df = self._de_df()
        de_df.loc[de_df["contrast"] == "TAMR_vs_MCF7", "log2fc"] = math.nan
        with pytest.raises(ValueError, match="tamr_vs_mcf7_log2fc is not finite"):
            build_paics_benchmark(self._labels_df(), de_df, self._master_df_without_paics())

    def test_de_fdr_out_of_range_raises(self):
        de_df = self._de_df()
        de_df.loc[de_df["contrast"] == "TAMR_vs_FASR", "fdr"] = -0.01
        with pytest.raises(ValueError, match=r"tamr_vs_fasr_fdr is outside \[0, 1\]"):
            build_paics_benchmark(self._labels_df(), de_df, self._master_df_without_paics())


# --- wording: pharmacological inhibition is not asserted to equal knockout


class TestCautiousInhibitionWording:
    def test_module_docstring_does_not_overclaim_inhibition_mimics_knockout(self):
        docstring = " ".join(candidate_evidence_summary.__doc__.split())
        assert "mimics its knockout" not in docstring
        assert "phenocopy" in docstring
        assert "requires experimental validation" in docstring


# --- build_plot_inputs: plain tidy tables, no figures ----------------------


class TestBuildPlotInputs:
    def _pipeline(self, rows):
        master_df = pd.DataFrame(rows)
        evidence_df = build_evidence_summary(master_df, RNA_FDR)
        sensitising_df, tolerance_df = split_by_direction(evidence_df)
        shortlist_df = select_primary_shortlist(sensitising_df)
        return evidence_df, shortlist_df

    def test_contrasts_table_is_long_format_two_rows_per_gene(self):
        evidence_df, shortlist_df = self._pipeline(
            [_row("A", crispr_effect_size=-1.0), _row("B", crispr_effect_size=2.0)]
        )
        plot_inputs = build_plot_inputs(evidence_df, shortlist_df)
        contrasts = plot_inputs["gse118713_contrasts"]
        assert len(contrasts) == 4  # 2 genes x 2 contrasts
        assert set(contrasts["contrast"]) == {"TAMR_vs_MCF7", "TAMR_vs_FASR"}
        assert "crispr_direction" in contrasts.columns
        assert "evidence_class" in contrasts.columns

    def test_shortlist_matrix_contains_only_shortlisted_genes(self):
        evidence_df, shortlist_df = self._pipeline(
            [
                _row("A", crispr_effect_size=-1.0, mcf7_fdr=0.01, fasr_fdr=0.01),
                _row("B", crispr_effect_size=-1.0, mcf7_fdr=0.5, fasr_fdr=0.5),
            ]
        )
        plot_inputs = build_plot_inputs(evidence_df, shortlist_df)
        assert list(plot_inputs["shortlist_matrix"]["gene_symbol"]) == ["A"]

    def test_crispr_effects_table_covers_all_genes_regardless_of_direction(self):
        evidence_df, shortlist_df = self._pipeline(
            [_row("A", crispr_effect_size=-1.0), _row("B", crispr_effect_size=2.0)]
        )
        plot_inputs = build_plot_inputs(evidence_df, shortlist_df)
        assert set(plot_inputs["crispr_effects"]["gene_symbol"]) == {"A", "B"}

    def test_plot_inputs_do_not_duplicate_sensitising_or_tolerance_keys(self):
        # Codex flagged the previous **plot_inputs spread as producing
        # duplicate "sensitising"/"tolerance" keys against the caller's
        # own dict entries. build_plot_inputs must not emit those keys.
        evidence_df, shortlist_df = self._pipeline(
            [_row("A", crispr_effect_size=-1.0), _row("B", crispr_effect_size=2.0)]
        )
        plot_inputs = build_plot_inputs(evidence_df, shortlist_df)
        assert "sensitising" not in plot_inputs
        assert "tolerance" not in plot_inputs


# --- load_master_table: structural validation ------------------------------


class TestLoadMasterTable:
    def _cfg(self, tmp_path, path, expected_n_hits):
        return EvidenceSummaryConfig(
            master_table_tsv=path,
            expected_n_hits=expected_n_hits,
            rna_significance_fdr=RNA_FDR,
            labels_parquet=tmp_path / "labels.parquet",
            gse118713_de_tsv_gz=tmp_path / "de.tsv.gz",
            evidence_summary_tsv=tmp_path / "es.tsv",
            sensitisation_candidates_tsv=tmp_path / "sens.tsv",
            tolerance_hits_tsv=tmp_path / "tol.tsv",
            shortlist_tsv=tmp_path / "sl.tsv",
            secondary_candidates_tsv=tmp_path / "sec.tsv",
            paics_benchmark_tsv=tmp_path / "paics.tsv",
            plot_crispr_effects_tsv=tmp_path / "p1.tsv",
            plot_gse118713_contrasts_tsv=tmp_path / "p2.tsv",
            plot_shortlist_matrix_tsv=tmp_path / "p3.tsv",
        )

    def test_row_count_mismatch_raises(self, tmp_path):
        path = tmp_path / "master.tsv"
        pd.DataFrame([_row("A")]).to_csv(path, sep="\t", index=False)
        cfg = self._cfg(tmp_path, path, expected_n_hits=28)
        with pytest.raises(ValueError, match="expected exactly 28"):
            load_master_table(cfg)

    def test_missing_column_raises(self, tmp_path):
        path = tmp_path / "master.tsv"
        pd.DataFrame([{"gene_symbol": "A"}]).to_csv(path, sep="\t", index=False)
        cfg = self._cfg(tmp_path, path, expected_n_hits=1)
        with pytest.raises(ValueError, match="missing required columns"):
            load_master_table(cfg)

    def test_duplicate_gene_symbol_raises(self, tmp_path):
        path = tmp_path / "master.tsv"
        pd.DataFrame([_row("A"), _row("A")]).to_csv(path, sep="\t", index=False)
        cfg = self._cfg(tmp_path, path, expected_n_hits=2)
        with pytest.raises(ValueError, match="duplicate gene_symbol"):
            load_master_table(cfg)

    def test_non_finite_crispr_effect_size_raises(self, tmp_path):
        path = tmp_path / "master.tsv"
        pd.DataFrame([_row("A", crispr_effect_size=math.nan)]).to_csv(path, sep="\t", index=False)
        cfg = self._cfg(tmp_path, path, expected_n_hits=1)
        with pytest.raises(ValueError, match="non-finite crispr_effect_size"):
            load_master_table(cfg)

    def test_non_finite_crispr_fdr_raises(self, tmp_path):
        path = tmp_path / "master.tsv"
        pd.DataFrame([_row("A", crispr_fdr=math.inf)]).to_csv(path, sep="\t", index=False)
        cfg = self._cfg(tmp_path, path, expected_n_hits=1)
        with pytest.raises(ValueError, match="non-finite crispr_fdr"):
            load_master_table(cfg)

    def test_crispr_fdr_out_of_range_raises(self, tmp_path):
        path = tmp_path / "master.tsv"
        pd.DataFrame([_row("A", crispr_fdr=1.5)]).to_csv(path, sep="\t", index=False)
        cfg = self._cfg(tmp_path, path, expected_n_hits=1)
        with pytest.raises(ValueError, match=r"crispr_fdr values outside \[0, 1\]"):
            load_master_table(cfg)

    def test_de_available_non_finite_rna_value_raises(self, tmp_path):
        path = tmp_path / "master.tsv"
        pd.DataFrame([_row("A", mcf7_log2fc=math.nan)]).to_csv(path, sep="\t", index=False)
        cfg = self._cfg(tmp_path, path, expected_n_hits=1)
        with pytest.raises(ValueError, match="non-finite tamr_vs_mcf7_log2fc"):
            load_master_table(cfg)

    def test_de_available_rna_fdr_out_of_range_raises(self, tmp_path):
        path = tmp_path / "master.tsv"
        pd.DataFrame([_row("A", fasr_fdr=-0.1)]).to_csv(path, sep="\t", index=False)
        cfg = self._cfg(tmp_path, path, expected_n_hits=1)
        with pytest.raises(ValueError, match=r"tamr_vs_fasr_fdr values outside \[0, 1\]"):
            load_master_table(cfg)

    def test_rna_unavailable_row_is_not_checked_for_finite_rna_values(self, tmp_path):
        # DE-unavailable rows carry NaN RNA fields by design (see
        # crispr_gse118713_integration) -- this must not be flagged.
        path = tmp_path / "master.tsv"
        pd.DataFrame(
            [_row("A", mcf7_log2fc=math.nan, mcf7_fdr=math.nan, de_na_reason="filtered_out_no_de_fit")]
        ).to_csv(path, sep="\t", index=False)
        cfg = self._cfg(tmp_path, path, expected_n_hits=1)
        assert len(load_master_table(cfg)) == 1


# --- real-repository smoke test: isolated output paths ---------------------


class TestRunAgainstRealConfig:
    def test_real_master_table_produces_expected_hierarchy_and_shortlist(self, tmp_path):
        with open(CONFIG_PATH) as f:
            config = yaml.safe_load(f)
        config = copy.deepcopy(config)
        config["candidate_evidence_summary"]["output"] = {
            "evidence_summary_tsv": str(tmp_path / "es.tsv"),
            "sensitisation_candidates_tsv": str(tmp_path / "sens.tsv"),
            "tolerance_hits_tsv": str(tmp_path / "tol.tsv"),
            "shortlist_tsv": str(tmp_path / "sl.tsv"),
            "secondary_candidates_tsv": str(tmp_path / "sec.tsv"),
            "paics_benchmark_tsv": str(tmp_path / "paics.tsv"),
            "plot_crispr_effects_tsv": str(tmp_path / "p1.tsv"),
            "plot_gse118713_contrasts_tsv": str(tmp_path / "p2.tsv"),
            "plot_shortlist_matrix_tsv": str(tmp_path / "p3.tsv"),
        }
        tmp_config_path = tmp_path / "config.yaml"
        with open(tmp_config_path, "w") as f:
            yaml.safe_dump(config, f)

        result = run_candidate_evidence_summary(str(tmp_config_path))
        evidence_df = result["evidence_summary"]

        assert len(evidence_df) == 28
        direction_counts = evidence_df["crispr_direction"].value_counts().to_dict()
        assert direction_counts[DIRECTION_SENSITISING] == 13
        assert direction_counts[DIRECTION_TOLERANCE] == 15

        sensitising_only = evidence_df.loc[evidence_df["crispr_direction"] == DIRECTION_SENSITISING]
        class_counts = sensitising_only["evidence_class"].value_counts().to_dict()
        assert class_counts[EVIDENCE_CLASS_PRIMARY] == 1
        assert class_counts[EVIDENCE_CLASS_SECONDARY] == 7
        assert class_counts[EVIDENCE_CLASS_NO_SIGNIFICANT_RNA] == 4
        assert class_counts[EVIDENCE_CLASS_RNA_UNAVAILABLE] == 1

        assert list(result["shortlist"]["gene_symbol"]) == ["USP34"]
        assert len(result["secondary"]) == 7

        paics_df = result["paics_benchmark"]
        assert len(paics_df) == 1
        assert paics_df.iloc[0]["gene_symbol"] == BENCHMARK_GENE_SYMBOL
        assert paics_df.iloc[0]["benchmark_label"] == BENCHMARK_LABEL
        assert "PAICS" not in set(result["shortlist"]["gene_symbol"])
        assert "PAICS" not in set(evidence_df["gene_symbol"])

        assert (tmp_path / "es.tsv").exists()
        assert (tmp_path / "sec.tsv").exists()
        assert (tmp_path / "paics.tsv").exists()
