"""Tests for the GSE245601 pseudobulk pipeline (gse245601_PREANALYSIS.md
section 13): design/eligibility audit, pseudobulk construction integrity,
QC computations, candidate extraction (BH scope, PAICS exclusion),
malignant-vs-nonmalignant context, and integration-table join integrity.

Real-data checks skip (not fail) if the pipeline outputs are not present
in a given checkout, matching this project's established pattern.
"""

from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import yaml

from src.gse245601_candidate_extraction import benjamini_hochberg, build_candidate_table, compute_patient_direction
from src.gse245601_candidate_integration import build_integrated_table, load_crispr_bulk_evidence
from src.gse245601_malignant_vs_nonmalignant import build_malignant_vs_nonmalignant_candidate_table
from src.gse245601_pseudobulk_qc import compute_log2cpm, compute_sample_correlations

REPO_ROOT = Path(__file__).parent.parent
CONFIG_PATH = REPO_ROOT / "config" / "config.yaml"


def _config() -> dict:
    with open(CONFIG_PATH) as f:
        return yaml.safe_load(f)


CANDIDATES_13 = [
    "USP34", "CTDNEP1", "EIF4ENIF1", "HMGB1", "KDM1A", "PET117", "TADA2B",
    "VEZF1", "TLK2", "ICK", "SUPT4H1", "TSR3", "USP17L29",
]


# --- Phase 2 audit table: re-verify every check the R script already ran ---


class TestDesignEligibilityAudit:
    def test_all_audit_checks_passed(self):
        path = REPO_ROOT / "results" / "tables" / "gse245601_pseudobulk" / "design_eligibility_audit.tsv"
        if not path.exists():
            pytest.skip("design/eligibility audit table not present in this checkout")
        audit = pd.read_csv(path, sep="\t")
        failed = audit.loc[~audit["pass"]]
        assert len(failed) == 0, f"failed audit checks:\n{failed}"

    def test_candidate_mapping_checks_cover_all_13_plus_paics(self):
        path = REPO_ROOT / "results" / "tables" / "gse245601_pseudobulk" / "design_eligibility_audit.tsv"
        if not path.exists():
            pytest.skip("design/eligibility audit table not present in this checkout")
        audit = pd.read_csv(path, sep="\t")
        mapping_checks = audit.loc[audit["check"].str.startswith("candidate_gene_unique_in_feature_space__"), "check"]
        genes_checked = {c.replace("candidate_gene_unique_in_feature_space__", "") for c in mapping_checks}
        assert genes_checked == set(CANDIDATES_13) | {"PAICS"}


# --- Phase 3 pseudobulk construction integrity ---


class TestPseudobulkConstructionIntegrity:
    def test_track_a_patient_treatment_parsing_and_one_each(self):
        path = REPO_ROOT / "results" / "tables" / "gse245601_pseudobulk" / "track_a_epithelial_metadata.tsv"
        if not path.exists():
            pytest.skip("Track A metadata not present in this checkout")
        meta = pd.read_csv(path, sep="\t")
        assert len(meta) == 20
        assert set(meta["patient"]) == {f"Tumor_{i:02d}" for i in range(1, 11)}
        for patient, sub in meta.groupby("patient"):
            assert sorted(sub["condition"]) == ["Control", "Tamoxifen"]
            assert len(sub) == 2
        # sample_id must be exactly "{patient}_{condition}"
        assert (meta["sample_id"] == meta["patient"] + "_" + meta["condition"]).all()

    def test_track_a_no_duplicate_cell_assignment(self):
        """A duplicated cell would inflate n_contributing_cells beyond the
        independently-known total epithelial cell count (29,175, already
        cross-verified against the frozen malignant summary table)."""
        path = REPO_ROOT / "results" / "tables" / "gse245601_pseudobulk" / "track_a_epithelial_metadata.tsv"
        if not path.exists():
            pytest.skip("Track A metadata not present in this checkout")
        meta = pd.read_csv(path, sep="\t")
        assert meta["n_contributing_cells"].sum() == 29175

    def test_track_b_eligibility_is_exactly_tumor_02_03_07(self):
        path = REPO_ROOT / "results" / "tables" / "gse245601_pseudobulk" / "track_b_malignant_metadata.tsv"
        if not path.exists():
            pytest.skip("Track B metadata not present in this checkout")
        meta = pd.read_csv(path, sep="\t")
        assert set(meta["patient"]) == {"Tumor_02", "Tumor_03", "Tumor_07"}
        assert len(meta) == 6
        for patient, sub in meta.groupby("patient"):
            assert sorted(sub["condition"]) == ["Control", "Tamoxifen"]

    def test_track_b_no_duplicate_cell_assignment(self):
        path = REPO_ROOT / "results" / "tables" / "gse245601_pseudobulk" / "track_b_malignant_metadata.tsv"
        if not path.exists():
            pytest.skip("Track B metadata not present in this checkout")
        meta = pd.read_csv(path, sep="\t")
        assert meta["n_contributing_cells"].sum() == 1401  # 317 (Tumor_02) + 858 (Tumor_03) + 226 (Tumor_07)

    def test_pseudobulk_column_sums_equal_metadata_library_size(self):
        """Pseudobulk sums equal the contributing single-cell raw counts:
        each pseudobulk sample's own count-matrix column sum must exactly
        equal its recorded total_library_size."""
        counts_path = REPO_ROOT / "results" / "tables" / "gse245601_pseudobulk" / "track_a_epithelial_counts.tsv.gz"
        meta_path = REPO_ROOT / "results" / "tables" / "gse245601_pseudobulk" / "track_a_epithelial_metadata.tsv"
        if not counts_path.exists():
            pytest.skip("Track A counts not present in this checkout")
        counts = pd.read_csv(counts_path, sep="\t").set_index("gene")
        meta = pd.read_csv(meta_path, sep="\t").set_index("sample_id")
        for sample_id in counts.columns:
            assert int(counts[sample_id].sum()) == int(meta.loc[sample_id, "total_library_size"])

    def test_malignant_cells_are_a_subset_of_epithelial_per_sample(self):
        """Every Track B (malignant) sample's contributing-cell count must
        not exceed the corresponding Track A (all-epithelial) sample's
        count for the same sample_id."""
        a_path = REPO_ROOT / "results" / "tables" / "gse245601_pseudobulk" / "track_a_epithelial_metadata.tsv"
        b_path = REPO_ROOT / "results" / "tables" / "gse245601_pseudobulk" / "track_b_malignant_metadata.tsv"
        if not a_path.exists() or not b_path.exists():
            pytest.skip("pseudobulk metadata not present in this checkout")
        a = pd.read_csv(a_path, sep="\t").set_index("sample_id")
        b = pd.read_csv(b_path, sep="\t").set_index("sample_id")
        for sample_id in b.index:
            assert b.loc[sample_id, "n_contributing_cells"] <= a.loc[sample_id, "n_contributing_cells"]


# --- Phase 4 QC computations ---


class TestPseudobulkQcComputations:
    def test_compute_log2cpm_sums_to_expected_scale(self):
        counts = pd.DataFrame({"s1": [10, 20, 70], "s2": [0, 50, 50]}, index=["g1", "g2", "g3"])
        log2cpm = compute_log2cpm(counts, ["s1", "s2"])
        cpm_s1 = (2 ** log2cpm["s1"]) - 1
        assert cpm_s1.sum() == pytest.approx(1e6, rel=1e-6)

    def test_sample_correlation_matrix_is_symmetric_with_unit_diagonal(self):
        rng = np.random.default_rng(0)
        counts = pd.DataFrame(rng.integers(1, 1000, size=(50, 4)), index=[f"g{i}" for i in range(50)], columns=["a", "b", "c", "d"])
        corr_long = compute_sample_correlations(counts, ["a", "b", "c", "d"])
        corr = corr_long.pivot(index="sample_id_1", columns="sample_id_2", values="correlation")
        assert np.allclose(np.diag(corr.to_numpy()), 1.0)
        assert np.allclose(corr.to_numpy(), corr.to_numpy().T)


# --- Phase 6 candidate extraction: BH scope, PAICS exclusion ---


class TestBenjaminiHochberg:
    def test_matches_known_reference_values(self):
        # classic textbook example: p = [0.01, 0.02, 0.03, 0.04, 0.5]
        p = pd.Series([0.01, 0.02, 0.03, 0.04, 0.5])
        bh = benjamini_hochberg(p)
        expected = pd.Series([0.05, 0.05, 0.05, 0.05, 0.5])
        assert np.allclose(bh.to_numpy(), expected.to_numpy())

    def test_monotonic_nondecreasing_after_sorting_by_pvalue(self):
        rng = np.random.default_rng(1)
        p = pd.Series(rng.uniform(0, 1, 30))
        bh = benjamini_hochberg(p)
        order = p.sort_values().index
        assert (bh.loc[order].diff().dropna() >= -1e-12).all()


class TestBuildCandidateTableBhScope:
    def _de_df(self, genes, p_values):
        return pd.DataFrame({"gene": genes, "log2fc": [0.1] * len(genes), "avg_log_cpm": [5.0] * len(genes), "p_value": p_values, "fdr": p_values})

    def _pseudobulk(self, genes, patients=("Tumor_01", "Tumor_02")):
        sample_ids = [f"{p}_{c}" for p in patients for c in ("Control", "Tamoxifen")]
        rng = np.random.default_rng(2)
        counts = pd.DataFrame(rng.integers(10, 1000, size=(len(genes), len(sample_ids))), index=genes, columns=sample_ids)
        metadata = pd.DataFrame(
            {"sample_id": sample_ids, "patient": [s.rsplit("_", 1)[0] for s in sample_ids], "condition": [s.rsplit("_", 1)[1] for s in sample_ids]}
        )
        return counts, metadata

    def test_bh_family_excludes_untested_candidate(self):
        tested_genes = ["g1", "g2", "g3"]  # only 3 of the "13" are actually testable here
        de_df = self._de_df(tested_genes, [0.01, 0.02, 0.5])
        counts, metadata = self._pseudobulk(tested_genes)
        candidate_genes = tested_genes + ["g4_untested"]
        out = build_candidate_table(de_df, counts, metadata, candidate_genes, "test_track")

        assert set(out.loc[out["tested"], "gene"]) == set(tested_genes)
        assert out.loc[out["gene"] == "g4_untested", "tested"].iloc[0] == False  # noqa: E712
        assert np.isnan(out.loc[out["gene"] == "g4_untested", "candidate_set_bh_fdr"].iloc[0])
        # BH over exactly the 3 tested p-values must match a from-scratch computation
        manual_bh = benjamini_hochberg(pd.Series([0.01, 0.02, 0.5]))
        computed = out.loc[out["gene"].isin(tested_genes)].set_index("gene")["candidate_set_bh_fdr"]
        for gene, p_idx in zip(tested_genes, range(3)):
            assert computed[gene] == pytest.approx(manual_bh.iloc[p_idx])

    def test_paics_never_appears_in_candidate_table(self):
        tested_genes = ["g1", "g2"]
        de_df = self._de_df(tested_genes + ["PAICS"], [0.01, 0.02, 0.001])
        counts, metadata = self._pseudobulk(tested_genes + ["PAICS"])
        out = build_candidate_table(de_df, counts, metadata, tested_genes, "test_track")
        assert "PAICS" not in set(out["gene"])


# --- Phase 8: malignant-vs-nonmalignant paired-by-patient test ---


class TestMalignantVsNonmalignantCandidateTable:
    def test_uses_patient_as_unit_not_per_cell(self):
        genes = ["g1", "PAICS"]
        patients = ["Tumor_02", "Tumor_03", "Tumor_07"]
        sample_ids = [f"{p}_{status}" for p in patients for status in ("malignant", "nonmalignant")]
        rng = np.random.default_rng(3)
        counts = pd.DataFrame(rng.integers(10, 1000, size=(len(genes), len(sample_ids))), index=genes, columns=sample_ids)
        metadata = pd.DataFrame(
            {
                "sample_id": sample_ids,
                "patient": [s.rsplit("_", 1)[0] for s in sample_ids],
                "malignancy_status": [s.rsplit("_", 1)[1] for s in sample_ids],
            }
        )
        out = build_malignant_vs_nonmalignant_candidate_table(counts, metadata, ["g1"], "PAICS")
        assert out.loc[out["gene"] == "g1", "n_patients"].iloc[0] == 3
        assert bool(out.loc[out["gene"] == "PAICS", "is_paics_benchmark"].iloc[0])
        # PAICS must be excluded from the candidate-set BH column
        assert np.isnan(out.loc[out["gene"] == "PAICS", "candidate_set_bh_fdr"].iloc[0])

    def test_real_data_five_eligible_patients_if_present(self):
        path = REPO_ROOT / "results" / "tables" / "gse245601_pseudobulk" / "malignant_vs_nonmalignant_metadata.tsv"
        if not path.exists():
            pytest.skip("malignant-vs-nonmalignant pseudobulk metadata not present in this checkout")
        meta = pd.read_csv(path, sep="\t")
        assert set(meta["patient"]) == {"Tumor_02", "Tumor_03", "Tumor_07", "Tumor_09", "Tumor_10"}
        assert len(meta) == 10


# --- Phase 9: integration table join integrity ---


class TestIntegratedTableJoinIntegrity:
    def _crispr_bulk(self, genes):
        return pd.DataFrame(
            {
                "gene_symbol": genes,
                "crispr_effect_size": [-1.0] * len(genes),
                "crispr_fdr": [0.01] * len(genes),
                "crispr_direction": ["sensitising_knockout"] * len(genes),
                "tamr_vs_mcf7_log2fc": [0.2] * len(genes),
                "tamr_vs_mcf7_fdr": [0.1] * len(genes),
                "evidence_class": ["PRIMARY_RESISTANCE_SUPPORT"] * len(genes),
            }
        )

    def _sc_track(self, genes, track_label):
        return pd.DataFrame(
            {
                "gene": genes,
                "track": [track_label] * len(genes),
                "tested": [True] * len(genes),
                "log2fc": [0.1] * len(genes),
                "p_value": [0.2] * len(genes),
                "candidate_set_bh_fdr": [0.5] * len(genes),
                "direction": ["up"] * len(genes),
                "n_patients_up": [3] * len(genes),
                "n_patients_down": [1] * len(genes),
            }
        )

    def test_join_preserves_all_13_candidates_exactly_once(self):
        genes = CANDIDATES_13
        crispr_bulk = self._crispr_bulk(genes)
        track_a = self._sc_track(genes, "track_a_epithelial")
        track_b = self._sc_track(genes, "track_b_malignant")
        malignant_context = pd.DataFrame(
            {
                "gene": genes + ["PAICS"],
                "is_paics_benchmark": [False] * len(genes) + [True],
                "tested": [True] * (len(genes) + 1),
                "mean_delta_malignant_minus_nonmalignant": [0.1] * (len(genes) + 1),
                "p_value": [0.3] * (len(genes) + 1),
                "candidate_set_bh_fdr": [0.6] * len(genes) + [np.nan],
                "direction": ["up_in_malignant"] * (len(genes) + 1),
                "n_patients": [5] * (len(genes) + 1),
            }
        )
        out = build_integrated_table(crispr_bulk, track_a, track_b, malignant_context, genes, track_b_n_pairs=3)
        assert len(out) == 13
        assert out["gene_symbol"].tolist() == genes  # exact order preserved
        assert not out["gene_symbol"].duplicated().any()
        assert "PAICS" not in set(out["gene_symbol"])
        assert (out["sc_track_b_exploratory_n3"]).all()
        assert (out["sc_track_b_n_patient_pairs"] == 3).all()

    def test_real_crispr_bulk_evidence_has_all_13_candidates(self):
        evidence_path = REPO_ROOT / "results" / "tables" / "candidate_evidence_summary.tsv"
        if not evidence_path.exists():
            pytest.skip("frozen candidate evidence summary not present in this checkout")
        loaded = load_crispr_bulk_evidence(evidence_path, CANDIDATES_13)
        assert len(loaded) == 13
        assert set(loaded["gene_symbol"]) == set(CANDIDATES_13)

    def test_real_integrated_table_if_present(self):
        path = REPO_ROOT / "results" / "tables" / "gse245601_candidate_integration" / "integrated_crispr_bulk_singlecell_candidates.tsv"
        if not path.exists():
            pytest.skip("integrated candidate table not present in this checkout")
        out = pd.read_csv(path, sep="\t")
        assert len(out) == 13
        assert set(out["gene_symbol"]) == set(CANDIDATES_13)
        assert not out["gene_symbol"].duplicated().any()


# --- patient-direction consistency helper ---


class TestComputePatientDirection:
    def test_direction_signs_match_manual_calculation(self):
        # a stable "background" gene is required so library-size (CPM)
        # normalization is meaningful -- with only the gene of interest in
        # the matrix, CPM is trivially 1e6 in every sample regardless of
        # its raw count.
        genes = ["g1", "background"]
        sample_ids = ["Tumor_01_Control", "Tumor_01_Tamoxifen", "Tumor_02_Control", "Tumor_02_Tamoxifen"]
        counts = pd.DataFrame(
            {
                "Tumor_01_Control": [10, 1000],
                "Tumor_01_Tamoxifen": [100, 1000],
                "Tumor_02_Control": [100, 1000],
                "Tumor_02_Tamoxifen": [10, 1000],
            },
            index=genes,
        )
        metadata = pd.DataFrame(
            {"sample_id": sample_ids, "patient": ["Tumor_01", "Tumor_01", "Tumor_02", "Tumor_02"], "condition": ["Control", "Tamoxifen", "Control", "Tamoxifen"]}
        )
        out = compute_patient_direction(counts, metadata, "g1")
        assert out.set_index("patient").loc["Tumor_01", "direction"] == "up"
        assert out.set_index("patient").loc["Tumor_02", "direction"] == "down"
