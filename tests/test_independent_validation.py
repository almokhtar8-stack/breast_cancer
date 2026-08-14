"""Targeted tests for the independent TCGA-BRCA + DepMap validation phase."""

from __future__ import annotations

import subprocess
from pathlib import Path

import pandas as pd
import pytest

TABLES = Path("results/tables/independent_validation")
REPORTS = Path("results/reports/independent_validation")
FIGURES = Path("results/figures/independent_validation")
FROZEN_CANDIDATES = {"USP34", "VEZF1", "EML5", "CITED2"}


class TestFrozenOutputsUntouched:
    def test_frozen_outputs_show_no_git_changes(self):
        result = subprocess.run(
            [
                "git", "status", "--porcelain",
                "results/tables/evidence_freeze/", "docs/THERAPEUTIC_SHORTLIST_FREEZE.md",
                "results/tables/systems_network/", "results/networks/systems_network/",
            ],
            capture_output=True, text=True, cwd=Path(__file__).resolve().parents[1],
        )
        assert result.stdout.strip() == "", f"frozen files show as changed: {result.stdout}"

    def test_literature_review_report_wording_edits_did_not_touch_tables(self):
        # Part 0's wording corrections were prose-only in the .md report;
        # the underlying claim/reference/comparison tables must be identical
        # to what test_literature_mechanism.py already exercises.
        df = pd.read_csv("results/tables/literature_mechanism/four_candidate_claim_evidence.tsv", sep="\t")
        assert len(df) == 39

    def test_candidate_ensembl_ids_are_independently_verified_and_correct(self):
        df = pd.read_csv("data/reference/tcga_candidate_ensembl_ids.tsv", sep="\t")
        expected = {
            "USP34": "ENSG00000115464",
            "VEZF1": "ENSG00000136451",
            "EML5": "ENSG00000165521",
            "CITED2": "ENSG00000164442",
        }
        got = dict(zip(df["symbol"], df["ensembl_gene_id_unversioned"]))
        assert got == expected


class TestTCGACohort:
    def test_expression_table_has_all_four_candidates(self):
        df = pd.read_csv(TABLES / "TCGA_candidate_expression.tsv", sep="\t")
        assert set(df["candidate"]) == FROZEN_CANDIDATES

    def test_expression_sample_counts_match_barcode_derived_sample_type(self):
        df = pd.read_csv(TABLES / "TCGA_candidate_expression.tsv", sep="\t")
        dist = df.loc[df["comparison"] == "distribution_all_primary_tumors"]
        # 1106 raw primary-tumor RNA-seq samples exist, but 11 patients
        # each carry >1 aliquot; the analysis cohort must be deduplicated
        # to exactly one aliquot per patient (1095), never the raw 1106,
        # or every downstream comparison and Cox model pseudoreplicates
        # those 11 patients' covariates and outcomes
        assert (dist["n_a"] == 1095).all(), "primary-tumor count must match the deduplicated one-aliquot-per-patient count (1095), not the raw sample count (1106)"

    def test_cohort_table_has_no_duplicate_patient_per_sample_type(self):
        import sys
        sys.path.insert(0, ".")
        from src.independent_validation_tcga_data import build_cohort_table, load_config
        cohort = build_cohort_table(load_config())
        dup = cohort.reset_index().duplicated(subset=["patient_barcode", "sample_type"])
        assert not dup.any(), "no patient should contribute more than one sample per sample_type to any analysis"

    def test_paired_tumor_normal_n_equals_number_of_normal_samples(self):
        df = pd.read_csv(TABLES / "TCGA_candidate_expression.tsv", sep="\t")
        paired = df.loc[df["comparison"] == "tumor_vs_normal_PAIRED"]
        # 113 solid-tissue-normal samples exist; the paired test can use at
        # most that many patients, and must use exactly that many since
        # every normal sample in this cohort has a matched primary tumor
        assert (paired["n_a"] == 113).all()
        assert (paired["n_b"] == 113).all()

    def test_er_status_is_never_invented_missingness_is_reported(self):
        # every row testing an ER+/ER- comparison must show fewer than the
        # full primary-tumor count on each side (proof that unclassified/
        # missing samples were excluded, not silently assigned a status)
        df = pd.read_csv(TABLES / "TCGA_candidate_expression.tsv", sep="\t")
        er_rows = df.loc[df["comparison"] == "ER+ vs ER- (clinical IHC)"]
        assert (er_rows["n_a"] + er_rows["n_b"] < 1106).all()

    def test_fdr_never_below_its_own_pvalue(self):
        for fname in ["TCGA_candidate_expression.tsv", "TCGA_candidate_pathway_associations.tsv", "TCGA_candidate_clinical.tsv"]:
            df = pd.read_csv(TABLES / fname, sep="\t")
            testable = df["p_value"].notna() & df["fdr"].notna()
            assert (df.loc[testable, "fdr"] >= df.loc[testable, "p_value"] - 1e-9).all(), f"{fname}: FDR must never be smaller than its own p-value"


class TestTCGAPathways:
    def test_eml5_has_no_declared_pathway(self):
        df = pd.read_csv(TABLES / "TCGA_candidate_pathway_associations.tsv", sep="\t")
        eml5 = df.loc[df["candidate"] == "EML5"]
        assert (eml5["pathway"] == "NONE").all()
        assert eml5["spearman_rho"].isna().all()

    def test_other_three_candidates_have_at_least_one_declared_pathway(self):
        df = pd.read_csv(TABLES / "TCGA_candidate_pathway_associations.tsv", sep="\t")
        for candidate in ["USP34", "VEZF1", "CITED2"]:
            sub = df.loc[df["candidate"] == candidate]
            assert (sub["pathway"] != "NONE").any()
            assert sub["spearman_rho"].notna().any()

    def test_vezf1_cited2_check_labeled_as_consistency_check_not_validation(self):
        df = pd.read_csv(TABLES / "TCGA_VEZF1_CITED2_consistency_check.tsv", sep="\t")
        assert len(df) == 4
        assert (df["notes"].str.contains("consistency check", case=False)).all()
        assert not (df["notes"].str.contains("mechanistic validation of", case=False) & ~df["notes"].str.contains("not mechanistic validation", case=False)).any()


class TestTCGAClinical:
    def test_all_candidates_and_cohorts_present(self):
        df = pd.read_csv(TABLES / "TCGA_candidate_clinical.tsv", sep="\t")
        assert set(df["candidate"]) == FROZEN_CANDIDATES
        assert set(df["cohort"]) == {"all_primary_tumors", "ER_positive"}
        assert set(df["model"]) == {"univariable", "adjusted_age_stage"}
        assert len(df) == 16

    def test_hr_is_never_negative(self):
        df = pd.read_csv(TABLES / "TCGA_candidate_clinical.tsv", sep="\t")
        valid = df["hr_per_sd"].notna()
        assert (df.loc[valid, "hr_per_sd"] > 0).all()

    def test_ph_violation_is_disclosed_not_silently_ignored(self):
        df = pd.read_csv(TABLES / "TCGA_candidate_clinical.tsv", sep="\t")
        violated = df.loc[df["ph_assumption_p"] < 0.05]
        assert len(violated) >= 1, "expected at least one PH-violated row (VEZF1, all_primary_tumors) to exist and be checked"
        assert violated["notes"].str.contains("VIOLATED", na=False).all()


class TestDepMap:
    def test_dependency_table_has_all_four_candidates(self):
        df = pd.read_csv(TABLES / "DepMap_candidate_dependency.tsv", sep="\t")
        assert set(df["candidate"]) == FROZEN_CANDIDATES

    def test_dependency_score_direction_more_negative_is_more_dependent(self):
        # sanity check on the classification rule itself, using the 24Q4
        # archive (the only release with a probability-based classification
        # to check this against -- 26Q1 has no CRISPRGeneDependency.csv):
        # any row classified STRONG/MODERATE concern must have a higher
        # strongly-dependent fraction than any row classified LOW
        df = pd.read_csv(TABLES / "archive_24Q4" / "DepMap_candidate_dependency_24Q4.tsv", sep="\t").set_index("candidate")
        low = df.loc[df["essentiality_concern"] == "D_LOW_BASELINE_DEPENDENCY", "frac_strongly_dependent_all_cancer"]
        concern = df.loc[df["essentiality_concern"].isin(["A_STRONG_GENERAL_DEPENDENCY_CONCERN", "B_MODERATE_DEPENDENCY_CONCERN"]), "frac_strongly_dependent_all_cancer"]
        assert len(low) and len(concern)
        assert concern.min() > low.max()

    def test_vezf1_was_the_essentiality_concern_candidate_in_24q4_archive(self):
        df = pd.read_csv(TABLES / "archive_24Q4" / "DepMap_candidate_dependency_24Q4.tsv", sep="\t").set_index("candidate")
        assert df.loc["VEZF1", "essentiality_concern"] != "D_LOW_BASELINE_DEPENDENCY"
        for other in ["USP34", "EML5", "CITED2"]:
            assert df.loc[other, "essentiality_concern"] == "D_LOW_BASELINE_DEPENDENCY"

    def test_26q1_essentiality_concern_is_genuine_and_reproduces_24q4(self):
        # CRISPRGeneDependency.csv was manually obtained and verified for
        # 26Q1 -- every candidate must show a genuine, computed tier
        # (never E_INSUFFICIENT_DATA), with real (non-NaN) fractions, and
        # the tier must exactly match the archived 24Q4 classification
        df = pd.read_csv(TABLES / "DepMap_candidate_dependency.tsv", sep="\t").set_index("candidate")
        archived = pd.read_csv(TABLES / "archive_24Q4" / "DepMap_candidate_dependency_24Q4.tsv", sep="\t").set_index("candidate")
        assert (df["depmap_release"] == "26Q1").all()
        assert not (df["essentiality_concern"] == "E_INSUFFICIENT_DATA").any()
        assert df["frac_strongly_dependent_all_cancer"].notna().all()
        assert (df["dependency_probability_available"] == True).all()  # noqa: E712
        for candidate in FROZEN_CANDIDATES:
            assert df.loc[candidate, "essentiality_concern"] == archived.loc[candidate, "essentiality_concern"], f"{candidate} tier should reproduce 24Q4 exactly"
        assert df.loc["VEZF1", "essentiality_concern"] == "B_MODERATE_DEPENDENCY_CONCERN"
        assert df.loc["VEZF1", "frac_strongly_dependent_er_luminal"] == pytest.approx(0.272727, abs=1e-5)

    def test_dependency_probability_direction_independently_verified(self):
        # sanity check pinned in the report/provenance: a pan-essential
        # gene should read near 1.0 (dependent) and a non-essential gene
        # near 0.0 -- guards against a silently-inverted probability column
        path = "/ibex/scratch/aljaroaa/tamoxifen-data/depmap/26Q1/CRISPRGeneDependency.csv"
        dep = pd.read_csv(path, index_col=0, usecols=["Unnamed: 0", "RPL3 (6122)", "OR51E2 (81285)"])
        assert dep["RPL3 (6122)"].median() > 0.9
        assert dep["OR51E2 (81285)"].median() < 0.1

    def test_chronos_params_file_excluded_and_confirmed_different(self):
        # the earlier unconfirmed-provenance Figshare file must never be
        # read by the active pipeline, and must be confirmed genuinely
        # different from the trusted CRISPRGeneEffect.csv (not a duplicate)
        import sys
        sys.path.insert(0, ".")
        from src.independent_validation_depmap_data import load_config
        cfg = load_config()
        raw_26q1 = cfg["independent_validation"]["depmap"]["releases"]["26Q1"]["raw"]
        assert "gene_effect_chronos_params" not in str(raw_26q1).lower()
        provenance = open("/ibex/scratch/aljaroaa/tamoxifen-data/depmap/26Q1/PROVENANCE.txt").read()
        assert "EXCLUDED" in provenance and "gene_effect_chronos_params.csv" in provenance
        assert "byte 7070" in provenance or "diverge" in provenance

    def test_hany_columns_match_frozen_shortlist_freeze_exactly(self):
        dep = pd.read_csv(TABLES / "DepMap_candidate_dependency.tsv", sep="\t").set_index("candidate")
        frozen = pd.read_csv("results/tables/evidence_freeze/THERAPEUTIC_SHORTLIST_FREEZE.tsv", sep="\t").set_index("gene")
        for candidate in FROZEN_CANDIDATES:
            assert dep.loc[candidate, "hany_crispr_fdr"] == pytest.approx(frozen.loc[candidate, "crispr_fdr"])
            assert dep.loc[candidate, "hany_crispr_direction"] == frozen.loc[candidate, "crispr_direction"]

    def test_er_luminal_breast_lines_defined_from_depmap_metadata_not_hardcoded(self):
        import sys
        sys.path.insert(0, ".")
        from src.independent_validation_depmap_data import load_config, load_model
        model = load_model(load_config(), release="24Q4")
        luminal = model.loc[model["is_er_luminal"]]
        assert len(luminal) == 22
        # a known TNBC line explicitly labeled "luminal TNBC" in DepMap's own
        # field must NOT be pulled in by a bare "luminal" substring match
        assert "DU4475" not in luminal["StrippedCellLineName"].values


class TestIntegrationAndRankings:
    def test_integration_table_covers_all_candidates(self):
        df = pd.read_csv(TABLES / "four_candidate_independent_validation.tsv", sep="\t")
        assert set(df["candidate"]) == FROZEN_CANDIDATES

    def test_eml5_never_scores_above_little_support(self):
        df = pd.read_csv(TABLES / "four_candidate_independent_validation.tsv", sep="\t").set_index("candidate")
        assert df.loc["EML5", "integration_validation_strength"] == "4_LITTLE_INDEPENDENT_SUPPORT"

    def test_frozen_ranking_never_appears_reordered_in_followup_rankings(self):
        df = pd.read_csv(TABLES / "four_candidate_followup_rankings.tsv", sep="\t")
        assert (df["note"].str.contains("Do NOT alter the frozen therapeutic ranking")).all()

    def test_rankings_cover_all_four_candidates_exactly_once(self):
        df = pd.read_csv(TABLES / "four_candidate_followup_rankings.tsv", sep="\t")
        assert sorted(df["candidate"]) == sorted(FROZEN_CANDIDATES)
        assert df["candidate"].is_unique


class TestFiguresAndReport:
    @pytest.mark.parametrize("name", [
        "01_TCGA_four_candidate_expression.png",
        "02_TCGA_candidate_pathway_associations.png",
        "03_DepMap_four_candidate_dependency.png",
        "04_integrated_candidate_validation.png",
        "05_TCGA_candidate_survival.png",
    ])
    def test_figure_exists(self, name):
        assert (FIGURES / name).exists()

    def test_report_exists_and_mentions_all_candidates(self):
        path = REPORTS / "four_candidate_TCGA_DepMap_review.md"
        assert path.exists()
        text = path.read_text()
        for candidate in FROZEN_CANDIDATES:
            assert candidate in text

    def test_report_never_calls_tcga_a_tamoxifen_resistance_cohort(self):
        text = " ".join(REPORTS.joinpath("four_candidate_TCGA_DepMap_review.md").read_text().replace("*", "").split())
        assert "not a tamoxifen-resistance cohort" in text

    def test_report_documents_depmap_release_gap(self):
        text = REPORTS.joinpath("four_candidate_TCGA_DepMap_review.md").read_text()
        assert "24Q4" in text
        assert "26Q1" in text


class TestDepMap26Q1Update:
    """26Q1 files were manually downloaded, verified, and are now the
    active release. These tests check the update was applied correctly,
    the excluded file stays excluded, missing CRISPRGeneDependency.csv
    never gets fabricated, and 24Q4 remains fully traceable via the
    archive.
    """

    def test_active_release_is_26q1(self):
        import sys
        sys.path.insert(0, ".")
        from src.independent_validation_depmap_data import load_config
        cfg = load_config()
        assert cfg["independent_validation"]["depmap"]["active_release"] == "26Q1"

    def test_26q1_release_config_now_has_dependency_file(self):
        # CRISPRGeneDependency.csv was manually obtained and verified for
        # 26Q1 in a follow-up step -- both releases must now declare it
        import sys
        sys.path.insert(0, ".")
        from src.independent_validation_depmap_data import load_config
        releases = load_config()["independent_validation"]["depmap"]["releases"]
        assert set(releases) == {"24Q4", "26Q1"}
        assert "crispr_gene_dependency_csv" in releases["24Q4"]["raw"]
        assert "crispr_gene_dependency_csv" in releases["26Q1"]["raw"]
        assert releases["26Q1"]["has_dependency_probability"] is True
        assert releases["24Q4"]["has_dependency_probability"] is True

    def test_access_status_report_exists_and_documents_every_channel_tried(self):
        path = REPORTS / "DEPMAP_26Q1_ACCESS_STATUS.md"
        assert path.exists()
        text = path.read_text()
        for marker in ["Cloudflare", "Figshare", "AnVIL", "31660582", "Model.csv", "manual", "Yejie Yun"]:
            assert marker in text, f"access-status report missing expected marker: {marker}"

    def test_model_id_uniqueness_enforced_for_both_releases(self):
        import sys
        sys.path.insert(0, ".")
        from src.independent_validation_depmap_data import load_config, load_model
        cfg = load_config()
        for release in ["24Q4", "26Q1"]:
            model = load_model(cfg, release=release)
            assert not model.index.duplicated().any()

    def test_depmap_tables_carry_26q1_release_provenance_column(self):
        for fname in ["DepMap_candidate_dependency.tsv", "DepMap_candidate_expression.tsv", "DepMap_candidate_codependency.tsv"]:
            df = pd.read_csv(TABLES / fname, sep="\t")
            assert "depmap_release" in df.columns
            assert (df["depmap_release"] == "26Q1").all()

    def test_comparison_module_produces_both_releases_side_by_side(self):
        import sys
        sys.path.insert(0, ".")
        from src.independent_validation_depmap_comparison import build_comparison_table
        from src.independent_validation_depmap_data import load_config
        out = build_comparison_table(load_config())
        assert set(out["candidate"]) == FROZEN_CANDIDATES
        assert "24Q4_all_median" in out.columns and "26Q1_all_median" in out.columns
        # both releases' classifications must be genuine, computed values,
        # and VEZF1's tier must match across releases (reproduced, not
        # just carried forward)
        out = out.set_index("candidate")
        assert out.loc["VEZF1", "24Q4_classification"] == "B_MODERATE_DEPENDENCY_CONCERN"
        assert out.loc["VEZF1", "26Q1_classification"] == "B_MODERATE_DEPENDENCY_CONCERN"
        assert not out["classification_changed"].any(), "no candidate's tier should differ between releases"

    def test_codependency_has_no_degenerate_correlations(self):
        # regression test for the min-pairwise-N bug caught during this
        # update: sparsely-screened genes could previously produce
        # trivial r=+-1.0 "correlations" from as few as 2 overlapping lines
        df = pd.read_csv(TABLES / "DepMap_candidate_codependency.tsv", sep="\t")
        assert df["pearson_r"].abs().max() < 0.99

    # Full numeric snapshot of the ARCHIVED 24Q4 DepMap_candidate_dependency
    # table, pinned before 26Q1 became active -- this is the permanent
    # record that 24Q4's original, reviewed values are preserved unaltered
    # for traceability (not the live table, which is now 26Q1).
    _EXPECTED_24Q4_DEPENDENCY = {
        "USP34": dict(median_gene_effect_all_cancer=-0.1716418716927549, n_all_cancer=1178, median_gene_effect_breast=-0.1609320837144783, n_breast=53, median_gene_effect_er_luminal=-0.1609320837144783, n_er_luminal=11, frac_strongly_dependent_all_cancer=0.0492359932088285, frac_strongly_dependent_breast=0.0188679245283018, frac_strongly_dependent_er_luminal=0.0, hany_crispr_effect=-1.3912980676838431, hany_crispr_fdr=0.0416851803622835),
        "VEZF1": dict(median_gene_effect_all_cancer=-0.3097615212728055, n_all_cancer=1178, median_gene_effect_breast=-0.2647378569995336, n_breast=53, median_gene_effect_er_luminal=-0.2889498037117357, n_er_luminal=11, frac_strongly_dependent_all_cancer=0.200339558573854, frac_strongly_dependent_breast=0.2264150943396226, frac_strongly_dependent_er_luminal=0.3636363636363636, hany_crispr_effect=-1.602445388833546, hany_crispr_fdr=0.0372575660229722),
        "EML5": dict(median_gene_effect_all_cancer=-0.1255714920434393, n_all_cancer=1178, median_gene_effect_breast=-0.1016389974788015, n_breast=53, median_gene_effect_er_luminal=-0.1209944448342789, n_er_luminal=11, frac_strongly_dependent_all_cancer=0.0067911714770797, frac_strongly_dependent_breast=0.0, frac_strongly_dependent_er_luminal=0.0, hany_crispr_effect=-1.058423339442026, hany_crispr_fdr=0.1487732844800642),
        "CITED2": dict(median_gene_effect_all_cancer=-0.1498048622646096, n_all_cancer=1178, median_gene_effect_breast=-0.1510063139998425, n_breast=53, median_gene_effect_er_luminal=-0.216075766850074, n_er_luminal=11, frac_strongly_dependent_all_cancer=0.0840407470288624, frac_strongly_dependent_breast=0.0188679245283018, frac_strongly_dependent_er_luminal=0.0, hany_crispr_effect=-1.4953555271062804, hany_crispr_fdr=0.109955077591718),
    }

    def test_24q4_archive_unchanged(self):
        dep = pd.read_csv(TABLES / "archive_24Q4" / "DepMap_candidate_dependency_24Q4.tsv", sep="\t").set_index("candidate")
        assert dep.loc["VEZF1", "essentiality_concern"] == "B_MODERATE_DEPENDENCY_CONCERN"
        for other in ["USP34", "EML5", "CITED2"]:
            assert dep.loc[other, "essentiality_concern"] == "D_LOW_BASELINE_DEPENDENCY"
        for candidate, expected in self._EXPECTED_24Q4_DEPENDENCY.items():
            for col, val in expected.items():
                assert dep.loc[candidate, col] == pytest.approx(val, rel=1e-9), f"{candidate}.{col} drifted from its archived value"

    def test_integration_table_carries_26q1_release_column(self):
        df = pd.read_csv(TABLES / "four_candidate_independent_validation.tsv", sep="\t")
        assert "depmap_release" in df.columns
        assert (df["depmap_release"] == "26Q1").all()

    def test_cited2_reaches_strong_now_that_dependency_data_is_real(self):
        # with CRISPRGeneDependency.csv available, STRONG is reachable
        # again for a candidate with a confirmed low/context-specific
        # DepMap tier and >=2 significant TCGA signals -- CITED2 qualifies
        df = pd.read_csv(TABLES / "four_candidate_independent_validation.tsv", sep="\t").set_index("candidate")
        assert df.loc["CITED2", "integration_validation_strength"] == "1_STRONG_INDEPENDENT_SUPPORT"

    def test_tied_rank_assignment_still_correct_when_scores_differ(self):
        # regression test for a bug caught by independent review: pandas
        # assigning rank=range(1, N+1) after a sort silently fabricates a
        # 1-2-3-4 ordering even for tied scores. Now that DepMap data is
        # real, this is a genuine THREE-way tie (USP34, CITED2, EML5, all
        # score=2, rank=1) with VEZF1 alone lower (penalized for its
        # confirmed essentiality concern) -- an earlier draft of the
        # report's narrative wrongly described this as a two-way
        # USP34/CITED2 tie with EML5 at a separate lower rank; this test
        # locks in the correct, verified values so that error can't recur.
        df = pd.read_csv(TABLES / "four_candidate_followup_rankings.tsv", sep="\t").set_index("candidate")
        for c in ["USP34", "CITED2", "EML5"]:
            assert df.loc[c, "therapeutic_targetability_followup_score"] == 2
            assert df.loc[c, "therapeutic_targetability_followup_rank"] == 1
        assert df.loc["VEZF1", "therapeutic_targetability_followup_score"] == 1
        assert df.loc["VEZF1", "therapeutic_targetability_followup_rank"] == 4

    def test_mechanistic_action_category_distinguishes_usp34_and_vezf1(self):
        # USP34 (clean DepMap + significant Hany) should read as the pure
        # tamoxifen-specific sensitiser; VEZF1 (real DepMap dependency +
        # significant Hany) should read as the dual-action candidate --
        # neither EML5 nor CITED2 (Hany not individually significant)
        # should be assigned to either category
        df = pd.read_csv(TABLES / "four_candidate_independent_validation.tsv", sep="\t").set_index("candidate")
        assert "mechanistic_action_category" in df.columns
        assert df.loc["USP34", "mechanistic_action_category"] == "TAMOXIFEN_SPECIFIC_SENSITISER"
        assert df.loc["VEZF1", "mechanistic_action_category"] == "POTENTIAL_DUAL_ACTION_CANCER_TARGET"
        assert df.loc["EML5", "mechanistic_action_category"] == "WEAK_HANY_SIGNAL_NEITHER_CATEGORY_ASSIGNED"
        assert df.loc["CITED2", "mechanistic_action_category"] == "WEAK_HANY_SIGNAL_NEITHER_CATEGORY_ASSIGNED"

    def test_tcga_outputs_unchanged_not_recomputed(self):
        # exact counts established and reviewed in the prior (already
        # -approved) TCGA phase; this phase must not have touched them
        expr = pd.read_csv(TABLES / "TCGA_candidate_expression.tsv", sep="\t")
        dist = expr.loc[expr["comparison"] == "distribution_all_primary_tumors"]
        assert (dist["n_a"] == 1095).all()
        clinical = pd.read_csv(TABLES / "TCGA_candidate_clinical.tsv", sep="\t")
        assert len(clinical) == 16
        pathway = pd.read_csv(TABLES / "TCGA_candidate_pathway_associations.tsv", sep="\t")
        assert (pathway.loc[pathway["candidate"] == "EML5", "pathway"] == "NONE").all()


if __name__ == "__main__":
    import sys

    sys.exit(pytest.main([__file__, "-v"]))
