"""Targeted tests for the systems-biology / pathway / network phase."""

from __future__ import annotations

import subprocess
from pathlib import Path

import pandas as pd
import pytest

TABLES = Path("results/tables/systems_network")
NETWORKS = Path("results/networks/systems_network")
FROZEN_CANDIDATES = {"USP34", "VEZF1", "EML5", "CITED2"}


class TestFrozenShortlistUntouched:
    def test_frozen_four_genes_exact(self):
        from src.systems_network_ranking import DATASET_LABELS  # noqa: F401
        from src.systems_network_candidate_pathways import CANDIDATES

        assert set(CANDIDATES) == FROZEN_CANDIDATES

    def test_evidence_freeze_files_untouched_by_git_status(self):
        result = subprocess.run(["git", "status", "--porcelain", "results/tables/evidence_freeze/", "docs/THERAPEUTIC_SHORTLIST_FREEZE.md"], capture_output=True, text=True, cwd=Path(__file__).resolve().parents[1])
        assert result.stdout.strip() == "", f"frozen evidence-freeze files show as changed: {result.stdout}"

    def test_gse245601_deepdive_files_untouched(self):
        result = subprocess.run(["git", "status", "--porcelain", "results/tables/gse245601_candidate_deepdive/", "results/figures/gse245601_candidate_deepdive/"], capture_output=True, text=True, cwd=Path(__file__).resolve().parents[1])
        assert result.stdout.strip() == "", f"frozen deep-dive outputs show as changed: {result.stdout}"


class TestDatasetSeparation:
    def test_resistance_consensus_excludes_gse245601(self):
        from src.systems_network_consensus import RESISTANCE_DATASETS

        assert "gse245601" not in RESISTANCE_DATASETS
        assert set(RESISTANCE_DATASETS) == {"gse118713", "gse240112", "gse111151"}

    def test_multimodal_table_shows_acute_as_separate_column(self):
        df = pd.read_csv(TABLES / "multimodal_pathway_convergence.tsv", sep="\t")
        assert "gse245601_NES" in df.columns
        assert "resistance_median_NES" in df.columns
        assert df["gse245601_NES"] is not df["resistance_median_NES"]

    def test_acute_response_category_never_requires_resistance_signal(self):
        df = pd.read_csv(TABLES / "multimodal_pathway_convergence.tsv", sep="\t")
        acute_rows = df.loc[df["convergence_category"] == "ACUTE_RESPONSE"]
        # ACUTE_RESPONSE requires the ABSENCE of a reproducible (>=2 dataset)
        # resistance signal -- not_tested_in_resistance/LOW_EVIDENCE/
        # SINGLE_DATASET are all valid non-reproducible states; only
        # STRONG_CONSENSUS/DIRECTIONAL_CONSENSUS (reproducible) must never
        # co-occur with ACUTE_RESPONSE.
        assert not acute_rows["resistance_consensus_class"].isin(["STRONG_CONSENSUS", "DIRECTIONAL_CONSENSUS"]).any()


class TestRankedGeneVectors:
    def test_full_gene_universe_not_significant_only(self):
        for name in ["gse118713", "gse240112", "gse111151", "gse245601"]:
            df = pd.read_csv(TABLES / f"{name}_ranked_genes.tsv", sep="\t")
            assert len(df) > 10000, f"{name} ranked list looks truncated to a significant-only subset ({len(df)} genes)"
            assert (df["fdr"] >= 0.05).any(), f"{name} ranked list contains only significant genes"

    def test_no_duplicate_gene_symbols_in_any_ranking(self):
        for name in ["gse118713", "gse240112", "gse111151", "gse245601", "crispr"]:
            df = pd.read_csv(TABLES / f"{name}_ranked_genes.tsv", sep="\t")
            assert df["gene"].is_unique, f"{name} ranked list has duplicate gene symbols"

    def test_gse118713_uses_true_model_statistic(self):
        df = pd.read_csv(TABLES / "gse118713_ranked_genes.tsv", sep="\t")
        assert (df["ranking_stat_method"] == "moderated_t").all()

    def test_fallback_datasets_document_their_statistic_choice(self):
        for name in ["gse240112", "gse111151", "gse245601"]:
            df = pd.read_csv(TABLES / f"{name}_ranked_genes.tsv", sep="\t")
            assert (df["ranking_stat_method"] == "sign(log2fc)*-log10(p_value)").all()


class TestCrisprSignSemantics:
    def test_sensitising_is_negative_effect(self):
        df = pd.read_csv(TABLES / "recurrent_leading_edge_genes.tsv", sep="\t")
        sens = df.loc[df["crispr_direction"] == "sensitising_KO"]
        assert (sens["crispr_effect"] < 0).all()
        assert (sens["crispr_fdr"] < 0.05).all()

    def test_tolerance_is_positive_effect(self):
        df = pd.read_csv(TABLES / "recurrent_leading_edge_genes.tsv", sep="\t")
        tol = df.loc[df["crispr_direction"] == "tolerance_associated_KO"]
        assert (tol["crispr_effect"] > 0).all()
        assert (tol["crispr_fdr"] < 0.05).all()

    def test_nonsignificant_never_labeled_sensitising_or_tolerance_ko(self):
        df = pd.read_csv(TABLES / "pathway_crispr_overlay.tsv", sep="\t")
        ns = df.loc[df["crispr_direction"].str.startswith("nonsignificant", na=False)]
        assert (ns["crispr_fdr"] >= 0.05).all()


class TestCodexReviewRegressions:
    """Regression tests for the 5 poster-blocking issues found by the
    Phase-29 Codex independent review (2026-08-13), so each fix stays fixed."""

    def test_directional_consensus_requires_nominal_significance_in_two_datasets(self):
        df = pd.read_csv(TABLES / "resistance_pathway_consensus.tsv", sep="\t")
        directional = df.loc[df["consensus_category"] == "DIRECTIONAL_CONSENSUS"]
        assert (directional["datasets_nominal"] >= 2).all()

    def test_single_dataset_reachable_even_when_all_three_tested(self):
        df = pd.read_csv(TABLES / "resistance_pathway_consensus.tsv", sep="\t")
        single = df.loc[df["consensus_category"] == "SINGLE_DATASET"]
        assert len(single) > 0
        assert (single["datasets_FDR05"] == 1).all()

    def test_candidate_leading_edge_requires_that_datasets_own_enrichment(self):
        membership = pd.read_csv(TABLES / "candidate_pathway_membership.tsv", sep="\t")
        gsea = {d: pd.read_csv(TABLES / f"gsea_{d}.tsv", sep="\t").set_index(["collection", "pathway"]) for d in ["gse118713", "gse240112", "gse111151", "gse245601"]}
        rows = membership.loc[membership["candidate_is_leading_edge_datasets"].fillna("") != ""]
        checked = 0
        for _, row in rows.iterrows():
            for dataset in row["candidate_is_leading_edge_datasets"].split(","):
                key = (row["collection"], row["pathway"])
                if key in gsea[dataset].index:
                    assert gsea[dataset].loc[key, "nom_pvalue"] < 0.05, f"{row['candidate']} {key} {dataset} not actually enriched"
                    checked += 1
        assert checked > 0

    def test_multimodal_mixed_takes_precedence_over_functional_only(self):
        df = pd.read_csv(TABLES / "multimodal_pathway_convergence.tsv", sep="\t")
        mixed_resistance_and_crispr = df.loc[(df["resistance_consensus_class"] == "MIXED") & (df["crispr_fdr"] < 0.05)]
        assert len(mixed_resistance_and_crispr) > 0
        assert (mixed_resistance_and_crispr["convergence_category"] == "MIXED").all()

    def test_string_physical_ppi_uses_string_own_physical_network(self):
        edges = pd.read_csv(NETWORKS / "edges.tsv", sep="\t")
        physical = edges.loc[edges["interaction_type"] == "physical_PPI"]
        assert len(physical) > 0
        assert physical["evidence_notes"].str.contains("network_type=physical membership=True").all()

    def test_max_same_pathway_dataset_count_never_exceeds_aggregate_count(self):
        # Codex review, second pass: max_same_pathway_dataset_count (strict,
        # same-pathway recurrence) must be <= leading_edge_dataset_count
        # (loose, any-qualifying-pathway aggregate) for every gene, and the
        # two must genuinely differ for at least one gene (otherwise the
        # distinction is vacuous).
        df = pd.read_csv(TABLES / "recurrent_leading_edge_genes.tsv", sep="\t")
        assert (df["max_same_pathway_dataset_count"] <= df["leading_edge_dataset_count"]).all()
        assert (df["max_same_pathway_dataset_count"] < df["leading_edge_dataset_count"]).any()

    def test_node_universe_category_b_ranked_by_strict_recurrence(self):
        # the node-universe top-40 selection must use the strict metric, not
        # the loose aggregate that let low-same-pathway-recurrence genes
        # (e.g. the previously-included RET/EGFR/VEGFA/DKK1) outrank
        # genuinely 3-dataset-same-pathway-recurrent genes
        le = pd.read_csv(TABLES / "recurrent_leading_edge_genes.tsv", sep="\t").set_index("gene")
        nodes = pd.read_csv(NETWORKS / "node_universe.tsv", sep="\t")
        b_genes = nodes.loc[nodes["node_category"].str.contains("resistance_leading_edge"), "gene"]
        b_genes = [g for g in b_genes if g in le.index]
        assert len(b_genes) > 0
        min_strict = le.loc[b_genes, "max_same_pathway_dataset_count"].min()
        # every gene in category B must be at least as strictly-recurrent as
        # the weakest gene actually selected (sanity: selection is sorted)
        excluded_stronger = le.loc[~le.index.isin(b_genes)].loc[lambda d: d["max_same_pathway_dataset_count"] > min_strict]
        # any excluded gene strictly stronger than the weakest included one
        # must have been excluded on a tie-break criterion, not have simply
        # been skipped -- i.e. none should exist once B is truly sorted desc
        assert len(excluded_stronger) == 0

    def test_string_functional_association_not_in_physical_network(self):
        edges = pd.read_csv(NETWORKS / "edges.tsv", sep="\t")
        functional = edges.loc[edges["interaction_type"] == "functional_association"]
        assert len(functional) > 0
        assert functional["evidence_notes"].str.contains("network_type=physical membership=False").all()


class TestPathwayDirection:
    def test_strong_consensus_requires_two_significant_same_direction(self):
        df = pd.read_csv(TABLES / "resistance_pathway_consensus.tsv", sep="\t")
        strong = df.loc[df["consensus_category"] == "STRONG_CONSENSUS"]
        assert (strong["datasets_FDR05"] >= 2).all()
        assert (strong["direction_consistency"] != "mixed").all()

    def test_mixed_category_has_conflicting_direction(self):
        df = pd.read_csv(TABLES / "resistance_pathway_consensus.tsv", sep="\t")
        mixed = df.loc[df["consensus_category"] == "MIXED"]
        assert (mixed["direction_consistency"] == "mixed").all()


class TestNetworkNodes:
    def test_no_duplicate_nodes(self):
        df = pd.read_csv(NETWORKS / "nodes.tsv", sep="\t")
        assert df["gene"].is_unique

    def test_node_universe_within_target_range(self):
        df = pd.read_csv(NETWORKS / "nodes.tsv", sep="\t")
        assert 50 <= len(df) <= 300

    def test_all_frozen_candidates_are_nodes(self):
        df = pd.read_csv(NETWORKS / "nodes.tsv", sep="\t")
        assert FROZEN_CANDIDATES <= set(df["gene"])

    def test_node_universe_rule_is_deterministic(self):
        from src.systems_network_node_universe import run_node_universe

        first = run_node_universe().sort_values("gene").reset_index(drop=True)
        second = run_node_universe().sort_values("gene").reset_index(drop=True)
        pd.testing.assert_frame_equal(first, second)


class TestNetworkEdges:
    def test_every_edge_has_source_and_type(self):
        df = pd.read_csv(NETWORKS / "edges.tsv", sep="\t")
        assert df["database_source"].notna().all()
        assert df["interaction_type"].notna().all()
        assert (df["database_source"] != "").all()

    def test_no_self_loops(self):
        df = pd.read_csv(NETWORKS / "edges.tsv", sep="\t")
        assert (df["source_gene"] != df["target_gene"]).all()

    def test_edge_endpoints_are_valid_nodes(self):
        nodes = set(pd.read_csv(NETWORKS / "nodes.tsv", sep="\t")["gene"])
        edges = pd.read_csv(NETWORKS / "edges.tsv", sep="\t")
        assert set(edges["source_gene"]) <= nodes
        assert set(edges["target_gene"]) <= nodes

    def test_correlation_never_labeled_physical_ppi(self):
        df = pd.read_csv(NETWORKS / "edges.tsv", sep="\t")
        physical = df.loc[df["interaction_type"] == "physical_PPI"]
        assert (physical["database_source"] == "STRING").all()
        notes = physical["evidence_notes"].fillna("")
        for note in notes:
            assert "escore" in note and "dscore" in note, "physical_PPI edge missing experimental/database evidence citation"

    def test_pathway_co_membership_edges_carry_a_pathway_name(self):
        df = pd.read_csv(NETWORKS / "edges.tsv", sep="\t")
        pw = df.loc[df["interaction_type"] == "pathway_co_membership"]
        assert (pw["pathway"] != "").all()

    def test_valid_interaction_types_only(self):
        df = pd.read_csv(NETWORKS / "edges.tsv", sep="\t")
        allowed = {"physical_PPI", "functional_association", "regulatory", "pathway_co_membership"}
        assert set(df["interaction_type"]) <= allowed


class TestCytoscapeExports:
    def test_cytoscape_files_present_and_valid(self):
        cyto = NETWORKS / "cytoscape"
        for name in ["network_nodes.tsv", "network_edges.tsv"] + [f"{c}_{k}.tsv" for c in ["USP34", "VEZF1", "EML5", "CITED2"] for k in ("nodes", "edges")]:
            path = cyto / name
            assert path.exists(), f"missing {path}"
            df = pd.read_csv(path, sep="\t")
            assert len(df.columns) > 0

    def test_cytoscape_edge_files_have_source_target_columns(self):
        cyto = NETWORKS / "cytoscape"
        df = pd.read_csv(cyto / "network_edges.tsv", sep="\t")
        assert {"source_gene", "target_gene"} <= set(df.columns)


class TestCandidatePathwayMapping:
    def test_no_pathway_included_without_evidence(self):
        df = pd.read_csv(TABLES / "candidate_pathway_membership.tsv", sep="\t")
        has_member = df["candidate_is_member"]
        has_le = df["candidate_is_leading_edge_datasets"].fillna("") != ""
        has_interactor = df["interactor_genes_in_pathway"].fillna("") != ""
        assert (has_member | has_le | has_interactor).all()

    def test_eml5_has_no_pathway_evidence(self):
        df = pd.read_csv(TABLES / "candidate_pathway_membership.tsv", sep="\t")
        assert (df["candidate"] != "EML5").all() or len(df.loc[df["candidate"] == "EML5"]) == 0


class TestNoGSE245601InResistanceCrossContamination:
    def test_leading_edge_genes_source_only_resistance_datasets(self):
        df = pd.read_csv(TABLES / "recurrent_leading_edge_genes.tsv", sep="\t")
        all_datasets = set()
        for val in df["leading_edge_datasets"].dropna():
            all_datasets.update(val.split(","))
        assert all_datasets <= {"gse118713", "gse240112", "gse111151"}


if __name__ == "__main__":
    import sys

    sys.exit(pytest.main([__file__, "-v"]))
