"""Targeted tests for the final USP34/VEZF1 GDSC pharmacogenomics phase."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

TABLES = Path("results/tables/final_pharmacogenomics")
REPORTS = Path("results/reports/final_pharmacogenomics")
FIGURES = Path("results/figures/final_pharmacogenomics")
GENES = ["USP34", "VEZF1"]


class TestConfigAndProvenance:
    def test_gdsc_section_declares_release_and_direction_source(self):
        from src.final_pharmacogenomics_gdsc_data import load_config
        cfg = load_config()
        g = cfg["final_pharmacogenomics"]["gdsc"]
        assert g["release"] == "8.5"
        assert "GDSC_Fitted_Data_Description.pdf" in g["response_metric_definitions"]["source_of_definitions"]

    def test_response_metric_directions_documented_lower_is_more_sensitive(self):
        from src.final_pharmacogenomics_gdsc_data import load_config
        g = load_config()["final_pharmacogenomics"]["gdsc"]["response_metric_definitions"]
        assert "LOWER = more sensitive" in g["LN_IC50"]
        assert "LOWER = more sensitive" in g["AUC"]

    def test_expression_source_is_reused_depmap_not_a_new_dataset(self):
        from src.final_pharmacogenomics_gdsc_data import load_config
        g = load_config()["final_pharmacogenomics"]["gdsc"]
        assert "26Q1" in g["expression_source"]
        assert "reused" in g["expression_source"] or "already" in g["expression_source"]

    def test_provenance_table_has_no_missing_fields(self):
        df = pd.read_csv(TABLES / "GDSC_data_provenance.tsv", sep="\t")
        assert df["value"].notna().all()

    def test_raw_gdsc_data_is_outside_git(self):
        import subprocess
        result = subprocess.run(
            ["git", "ls-files", "--error-unmatch"],
            cwd="/ibex/scratch/aljaroaa/tamoxifen-data/gdsc",
            capture_output=True, text=True,
        )
        assert result.returncode != 0, "raw GDSC files must not be tracked by git"


class TestCellLineJoin:
    def test_join_is_exact_id_based_and_logs_losses(self, caplog):
        import logging
        from src.final_pharmacogenomics_gdsc_data import build_breast_expression_joined, load_config
        caplog.set_level(logging.INFO)
        df = build_breast_expression_joined(load_config())
        assert df["is_breast"].all()
        assert df["SANGER_MODEL_ID"].nunique() <= 51
        assert any("matched to DepMap 26Q1" in r.message for r in caplog.records)

    def test_er_luminal_subset_is_smaller_than_full_breast_set(self):
        from src.final_pharmacogenomics_gdsc_data import build_breast_expression_joined, load_config
        df = build_breast_expression_joined(load_config())
        n_all = df["SANGER_MODEL_ID"].nunique()
        n_luminal = df.loc[df["is_er_luminal"], "SANGER_MODEL_ID"].nunique()
        assert 0 < n_luminal < n_all


class TestCompoundAvailability:
    def test_tamoxifen_and_fulvestrant_present_4oht_and_endoxifen_absent(self):
        df = pd.read_csv(TABLES / "GDSC_compound_availability.tsv", sep="\t")
        avail = df.groupby("compound")["present_in_gdsc"].any()
        assert avail["Tamoxifen"] == True  # noqa: E712
        assert avail["Fulvestrant"] == True  # noqa: E712
        assert avail["4-hydroxytamoxifen / 4-OHT"] == False  # noqa: E712
        assert avail["Endoxifen"] == False  # noqa: E712

    def test_fulvestrant_has_multiple_drug_ids_not_collapsed(self):
        df = pd.read_csv(TABLES / "GDSC_compound_availability.tsv", sep="\t")
        fulv = df[df["compound"] == "Fulvestrant"]
        assert fulv["drug_id"].nunique() >= 2


class TestAssociationTables:
    @pytest.mark.parametrize("gene", GENES)
    def test_every_row_meets_minimum_n(self, gene):
        df = pd.read_csv(TABLES / f"{gene}_GDSC_drug_associations.tsv", sep="\t")
        assert (df["n"] >= 15).all()

    @pytest.mark.parametrize("gene", GENES)
    def test_fdr_never_below_its_own_pvalue(self, gene):
        df = pd.read_csv(TABLES / f"{gene}_GDSC_drug_associations.tsv", sep="\t")
        assert (df["fdr"] >= df["p_value"] - 1e-9).all()

    @pytest.mark.parametrize("gene", GENES)
    def test_fdr_computed_separately_per_dataset_and_metric(self, gene):
        # regression guard: pooling GDSC1/GDSC2 or LN_IC50/AUC into one FDR
        # family would be wrong, since they are materially different
        # screening campaigns/metrics -- check each (dataset, metric) family
        # independently reproduces a fresh BH correction of its own p-values
        from statsmodels.stats.multitest import multipletests
        df = pd.read_csv(TABLES / f"{gene}_GDSC_drug_associations.tsv", sep="\t")
        for (ds, met), sub in df.groupby(["dataset", "metric"]):
            recompute = multipletests(sub["p_value"].to_numpy(), method="fdr_bh")[1]
            assert sub["fdr"].to_numpy() == pytest.approx(recompute, abs=1e-9)

    def test_no_pseudoreplicated_drug_id_dataset_cell_line_groups(self):
        # regression test for the caught pseudoreplication trap: a compound
        # (e.g. AZD7762) can carry more than one DRUG_ID even within one
        # GDSC release, so grouping by DRUG_NAME alone would inflate N.
        # The underlying joined table must have exactly one row per
        # (SANGER_MODEL_ID, DRUG_ID, DATASET).
        from src.final_pharmacogenomics_gdsc_data import build_breast_expression_joined, load_config
        df = build_breast_expression_joined(load_config())
        dupe_counts = df.groupby(["SANGER_MODEL_ID", "DRUG_ID", "DATASET"]).size()
        assert dupe_counts.max() == 1

    def test_azd7762_has_two_distinct_drug_ids_within_gdsc1(self):
        # confirms the real data property that motivated the DRUG_ID (not
        # DRUG_NAME) grouping rule above
        df = pd.read_csv(TABLES / "USP34_GDSC_drug_associations.tsv", sep="\t")
        azd = df[(df["drug_name"] == "AZD7762") & (df["dataset"] == "GDSC1")]
        assert azd["drug_id"].nunique() >= 2

    def test_azd7762_only_one_of_its_two_gdsc1_drug_ids_is_fdr_significant(self):
        # regression test: an earlier draft's report prose implied BOTH
        # GDSC1 re-screening batches of AZD7762 were FDR-significant --
        # only DRUG_ID 1402 is; DRUG_ID 1022 is directionally consistent
        # but not significant. Caught by Codex audit, fixed in the report.
        df = pd.read_csv(TABLES / "USP34_GDSC_drug_associations.tsv", sep="\t")
        azd = df[(df["drug_name"] == "AZD7762") & (df["dataset"] == "GDSC1")]
        sig_ids = set(azd.loc[azd["fdr"] < 0.05, "drug_id"])
        nonsig_ids = set(azd.loc[azd["fdr"] >= 0.05, "drug_id"])
        assert sig_ids == {1402}
        assert 1022 in nonsig_ids

    def test_association_grouping_is_drug_id_dataset_only_not_metadata_columns(self):
        # regression test for a real bug caught by Codex audit: grouping by
        # (DRUG_ID, DRUG_NAME, PUTATIVE_TARGET, PATHWAY_NAME, DATASET) let a
        # handful of NaN-metadata rows silently drop out of the FDR family
        # for a given DRUG_ID (fewer tests than the true 1,278). Grouping by
        # (DRUG_ID, DATASET) only recovers those rows.
        u = pd.read_csv(TABLES / "USP34_GDSC_drug_associations.tsv", sep="\t")
        assert len(u) == 1278

    def test_usp34_has_significant_hits_vezf1_has_none(self):
        u = pd.read_csv(TABLES / "USP34_GDSC_drug_associations.tsv", sep="\t")
        v = pd.read_csv(TABLES / "VEZF1_GDSC_drug_associations.tsv", sep="\t")
        assert (u["fdr"] < 0.05).sum() == 9
        assert (v["fdr"] < 0.05).sum() == 0

    def test_all_usp34_significant_hits_are_higher_expression_more_sensitive(self):
        u = pd.read_csv(TABLES / "USP34_GDSC_drug_associations.tsv", sep="\t")
        sig = u[u["fdr"] < 0.05]
        assert (sig["spearman_rho"] < 0).all()

    def test_tamoxifen_not_significant_for_either_gene(self):
        for gene in GENES:
            df = pd.read_csv(TABLES / f"{gene}_GDSC_drug_associations.tsv", sep="\t")
            tam = df[df["drug_name"] == "Tamoxifen"]
            assert len(tam) > 0
            assert (tam["fdr"] >= 0.05).all()

    def test_no_tead_hippo_compound_exists_in_gdsc(self):
        u = pd.read_csv(TABLES / "USP34_GDSC_drug_associations.tsv", sep="\t")
        assert not u["drug_name"].str.contains("TEAD|Hippo|Verteporfin", case=False, na=False).any()
        assert not u["pathway"].astype(str).str.contains("Hippo", case=False, na=False).any()


class TestTopAssociationsAndSubsets:
    def test_top_associations_only_contains_fdr_significant_or_top_effect_tiers(self):
        df = pd.read_csv(TABLES / "GDSC_top_associations.tsv", sep="\t")
        assert set(df["tier"]) <= {"FDR_SIGNIFICANT", "TOP_EFFECT_SIZE_NOT_NECESSARILY_SIGNIFICANT"}

    def test_er_luminal_subset_always_marked_exploratory(self):
        df = pd.read_csv(TABLES / "GDSC_ER_luminal_subset.tsv", sep="\t")
        assert (df["exploratory_only"] == True).all()  # noqa: E712
        assert df["note"].str.contains("EXPLORATORY|too small", case=False, regex=True).all()

    def test_er_luminal_subset_uses_exact_drug_id_azd7762(self):
        # AZD7762 has 2 DRUG_IDs within GDSC1 (see test above); the ER+/
        # luminal subset row for the significant one (1402) must reflect
        # the correct, non-pseudoreplicated n (<= number of ER+/luminal
        # lines in the panel, never doubled by pooling both DRUG_IDs)
        df = pd.read_csv(TABLES / "GDSC_ER_luminal_subset.tsv", sep="\t")
        azd = df[(df["drug_name"] == "AZD7762") & (df["drug_id"] == 1402)]
        assert len(azd) > 0
        assert (azd["n_er_luminal_lines"] <= azd["n_total_er_luminal_lines_in_panel"]).all()

    def test_ldn193189_does_not_replicate_in_er_luminal_subset(self):
        df = pd.read_csv(TABLES / "GDSC_ER_luminal_subset.tsv", sep="\t")
        ldn = df[(df["drug_name"] == "LDN-193189") & (df["metric"] == "LN_IC50")]
        assert len(ldn) == 1
        assert ldn.iloc[0]["p_value"] > 0.05


class TestCrosscheckAndClassification:
    def test_crosscheck_does_not_alter_frozen_shortlist_file(self):
        import subprocess
        result = subprocess.run(
            ["git", "status", "--porcelain", "results/tables/cross_dataset_genomewide/", "results/tables/evidence_freeze/"],
            capture_output=True, text=True, cwd=Path(__file__).resolve().parents[1],
        )
        assert result.stdout.strip() == ""

    def test_crosscheck_table_has_direction_columns_for_non_hany_datasets(self):
        # regression test: an earlier version of this table had no direction
        # info for gse118713/gse240112/gse111151/gse245601, which meant the
        # report could not honestly assess cross-dataset convergence
        df = pd.read_csv(TABLES / "GDSC_project_crosscheck.tsv", sep="\t")
        for col in ["gse118713_direction", "gse240112_tumor_direction", "gse111151_direction", "gse245601_epi_direction"]:
            assert col in df.columns

    def test_crosscheck_shows_multiple_genes_significant_in_gse118713_not_just_chek1(self):
        # regression test: an earlier draft's report claimed CHEK1 was "the
        # only convergent signal" -- CHEK2, FGFR2, FGFR4, JAK1, and NAMPT
        # also reach FDR<0.05 in the frozen GSE118713 dataset. Caught by
        # Codex audit; report rewritten to list all of them honestly.
        df = pd.read_csv(TABLES / "GDSC_project_crosscheck.tsv", sep="\t").set_index("drug_target_gene")
        sig_genes = set(df.index[df["gse118713_fdr"] < 0.05])
        assert {"CHEK2", "FGFR2", "FGFR4", "JAK1", "NAMPT"} <= sig_genes

    def test_fgfr2_and_fgfr4_significant_and_same_direction_in_both_resistance_datasets(self):
        df = pd.read_csv(TABLES / "GDSC_project_crosscheck.tsv", sep="\t").set_index("drug_target_gene")
        for gene in ["FGFR2", "FGFR4"]:
            row = df.loc[gene]
            assert row["gse118713_fdr"] < 0.05 and row["gse240112_tumor_fdr"] < 0.05
            assert row["gse118713_direction"] == "down_in_TAMR"
            assert row["gse240112_tumor_direction"] == "down"

    def test_chek1_crosscheck_shows_gse118713_nominal_hit_not_hany(self):
        # regression test: an earlier draft of this phase mislabeled
        # CHEK1's real GSE118713 nominal hit (FDR=0.0497) as a Hany CRISPR
        # hit -- CHEK1's actual Hany CRISPR FDR is 0.812 (not significant).
        # Caught by this test before the report was finalized; fixed in
        # both the report and final_pharmacogenomics_interpretation_data.py.
        df = pd.read_csv(TABLES / "GDSC_project_crosscheck.tsv", sep="\t").set_index("drug_target_gene")
        assert df.loc["CHEK1", "hany_crispr_fdr"] == pytest.approx(0.812, abs=1e-2)
        assert df.loc["CHEK1", "gse118713_fdr"] == pytest.approx(0.0497, abs=1e-3)

    def test_indirect_targeting_classification_valid_categories_only(self):
        from src.final_pharmacogenomics_interpretation_data import INDIRECT_TARGETING_CATEGORIES
        df = pd.read_csv(TABLES / "GDSC_indirect_targeting_classification.tsv", sep="\t")
        assert set(df["classification"]) <= INDIRECT_TARGETING_CATEGORIES

    def test_no_hit_classified_as_direct_target_or_upstream_regulator(self):
        df = pd.read_csv(TABLES / "GDSC_indirect_targeting_classification.tsv", sep="\t")
        assert not df["classification"].isin(["DIRECT_TARGET", "KNOWN_UPSTREAM_REGULATOR", "KNOWN_REQUIRED_PARTNER"]).any()

    def test_vezf1_tead_row_stays_pharmacogenomic_association_only(self):
        df = pd.read_csv(TABLES / "GDSC_indirect_targeting_classification.tsv", sep="\t")
        tead_row = df[df["drug_name"].str.contains("TEAD", case=False, na=False)]
        assert len(tead_row) == 1
        assert tead_row.iloc[0]["classification"] == "PHARMACOGENOMIC_ASSOCIATION_ONLY"


class TestFinalInterpretation:
    def test_covers_both_genes(self):
        df = pd.read_csv(TABLES / "GDSC_final_interpretation.tsv", sep="\t")
        assert set(df["candidate"]) == {"USP34", "VEZF1"}

    def test_gdsc_does_not_alter_usp34_lead_or_vezf1_backup(self):
        df = pd.read_csv(TABLES / "GDSC_final_interpretation.tsv", sep="\t")
        usp34_lead_row = df[(df["candidate"] == "USP34") & (df["question"].str.contains("LEAD"))]
        vezf1_backup_row = df[(df["candidate"] == "VEZF1") & (df["question"].str.contains("BACKUP"))]
        assert len(usp34_lead_row) == 1 and usp34_lead_row.iloc[0]["answer"].startswith("No")
        assert len(vezf1_backup_row) == 1 and vezf1_backup_row.iloc[0]["answer"].startswith("No")


class TestFiguresAndReport:
    @pytest.mark.parametrize("name", [
        "01_USP34_GDSC_drug_response.png",
        "02_VEZF1_GDSC_drug_response.png",
        "03_USP34_VEZF1_pharmacogenomic_summary.png",
    ])
    def test_figure_exists(self, name):
        assert (FIGURES / name).exists()

    def test_report_exists_and_mentions_both_candidates(self):
        path = REPORTS / "USP34_VEZF1_GDSC_review.md"
        assert path.exists()
        text = path.read_text()
        assert "USP34" in text and "VEZF1" in text

    def test_report_never_states_causal_language(self):
        text = " ".join(REPORTS.joinpath("USP34_VEZF1_GDSC_review.md").read_text().replace("*", "").split())
        for forbidden in ["USP34 causes", "VEZF1 causes", "USP34 inhibits", "should be combined with tamoxifen"]:
            assert forbidden not in text

    def test_report_uses_allowed_association_wording(self):
        text = REPORTS.joinpath("USP34_VEZF1_GDSC_review.md").read_text()
        assert "associated with" in text or "association" in text

    def test_report_states_zero_hippo_tead_compounds_honestly(self):
        text = REPORTS.joinpath("USP34_VEZF1_GDSC_review.md").read_text()
        assert "zero" in text.lower() and "TEAD" in text

    def test_report_documents_gdsc_pseudoreplication_fix(self):
        text = REPORTS.joinpath("USP34_VEZF1_GDSC_review.md").read_text()
        assert "DRUG_ID" in text and "pseudoreplicat" in text.lower()

    def test_report_correctly_attributes_chek1_echo_to_gse118713_not_hany(self):
        text = REPORTS.joinpath("USP34_VEZF1_GDSC_review.md").read_text()
        assert "GSE118713" in text
        assert "0.812" in text, "report must disclose CHEK1's real (non-significant) Hany CRISPR FDR alongside the GSE118713 echo"

    def test_report_does_not_claim_chek1_is_the_only_convergent_signal(self):
        text = " ".join(REPORTS.joinpath("USP34_VEZF1_GDSC_review.md").read_text().replace("*", "").split())
        assert "only convergent signal" not in text
        assert "FGFR2" in text and "FGFR4" in text


if __name__ == "__main__":
    import sys

    sys.exit(pytest.main([__file__, "-v"]))
