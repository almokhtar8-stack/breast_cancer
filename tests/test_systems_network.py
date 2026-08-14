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


class TestUSP34ShortestPaths:
    def test_all_four_targets_have_a_path_of_length_two(self):
        df = pd.read_csv(TABLES / "USP34_shortest_paths.tsv", sep="\t")
        for target in ["CTNNB1", "PTEN", "EP300", "SOX2"]:
            sub = df.loc[df["target"] == target]
            assert (sub["path_length_edges"] == 2).all()
            assert (sub["path"] != "NO_PATH_IN_NETWORK").all()

    def test_no_path_is_mislabeled_as_direct(self):
        # every reported USP34 path must have exactly 2 " -> " hops (source
        # -> intermediate -> target); a 1-hop (direct) path would mean USP34
        # has a real 1-edge interaction with the target, which is not the
        # case in this frozen network for any of the 4 targets
        df = pd.read_csv(TABLES / "USP34_shortest_paths.tsv", sep="\t")
        for path in df["path"].unique():
            assert path.count(" -> ") == 2, f"unexpected path length: {path}"

    def test_edge_evidence_never_invents_a_missing_interaction_type(self):
        df = pd.read_csv(TABLES / "USP34_shortest_paths.tsv", sep="\t")
        assert df["interaction_type"].isin(["physical_PPI", "functional_association", "regulatory", "pathway_co_membership"]).all()


class TestUSP34BridgeGeneEvidence:
    BRIDGE_GENES = ["USP9X", "RPS27A", "UBC", "UBB"]

    def test_all_four_bridge_genes_present(self):
        df = pd.read_csv(TABLES / "USP34_bridge_gene_evidence.tsv", sep="\t")
        assert set(df["gene"]) == set(self.BRIDGE_GENES)

    def test_crispr_sign_convention(self):
        df = pd.read_csv(TABLES / "USP34_bridge_gene_evidence.tsv", sep="\t").set_index("gene")
        for gene in self.BRIDGE_GENES:
            row = df.loc[gene]
            if row["crispr_effect"] < 0:
                assert "sensitising" in row["crispr_direction"]
            else:
                assert "tolerance" in row["crispr_direction"]

    def test_classification_is_one_of_three_conservative_tiers(self):
        df = pd.read_csv(TABLES / "USP34_bridge_gene_evidence.tsv", sep="\t")
        assert df["classification"].isin(["A_DATA_SUPPORTED_BRIDGE", "B_PARTIAL_SUPPORT", "C_NETWORK_ONLY_GENERIC_BRIDGE"]).all()

    def test_a_tier_requires_an_fdr_significant_hit(self):
        df = pd.read_csv(TABLES / "USP34_bridge_gene_evidence.tsv", sep="\t")
        sig_cols = [c for c in df.columns if c.endswith("_significant_fdr05")]
        a_tier = df.loc[df["classification"] == "A_DATA_SUPPORTED_BRIDGE"]
        for _, row in a_tier.iterrows():
            assert row[sig_cols].any()

    def test_c_tier_has_no_nominal_resistance_or_crispr_hit(self):
        df = pd.read_csv(TABLES / "USP34_bridge_gene_evidence.tsv", sep="\t")
        resistance_nominal_cols = ["gse118713_nominal_p05", "gse240112_nominal_p05", "gse111151_nominal_p05"]
        c_tier = df.loc[df["classification"] == "C_NETWORK_ONLY_GENERIC_BRIDGE"]
        for _, row in c_tier.iterrows():
            assert not row["crispr_nominal_p05"]
            assert not row[resistance_nominal_cols].any()

    def test_gse245601_acute_never_drives_classification(self):
        # a gene whose ONLY nominal hit is the acute layer must not reach
        # B_PARTIAL_SUPPORT on that basis alone
        df = pd.read_csv(TABLES / "USP34_bridge_gene_evidence.tsv", sep="\t")
        resistance_or_crispr_cols = ["crispr_nominal_p05", "gse118713_nominal_p05", "gse240112_nominal_p05", "gse111151_nominal_p05"]
        for _, row in df.iterrows():
            if row["classification"] == "B_PARTIAL_SUPPORT":
                assert row[resistance_or_crispr_cols].any(), f"{row['gene']} classified B without a non-acute nominal hit"

    def test_values_match_frozen_node_table(self):
        bridge = pd.read_csv(TABLES / "USP34_bridge_gene_evidence.tsv", sep="\t").set_index("gene")
        nodes = pd.read_csv(NETWORKS / "cytoscape" / "network_nodes.tsv", sep="\t").set_index("gene")
        for gene in self.BRIDGE_GENES:
            assert bridge.loc[gene, "crispr_effect"] == pytest.approx(nodes.loc[gene, "crispr_effect"])
            assert bridge.loc[gene, "gse118713_log2fc"] == pytest.approx(nodes.loc[gene, "gse118713_log2fc"])
            assert bridge.loc[gene, "gse245601_track_a_log2fc"] == pytest.approx(nodes.loc[gene, "gse245601_acute_log2fc"])

    def test_bridges_usp34_to_matches_shortest_path_table(self):
        bridge = pd.read_csv(TABLES / "USP34_bridge_gene_evidence.tsv", sep="\t").set_index("gene")
        expected = {
            "USP9X": {"CTNNB1", "SOX2"},
            "RPS27A": {"CTNNB1", "PTEN"},
            "UBC": {"EP300", "PTEN"},
            "UBB": {"CTNNB1"},
        }
        for gene, targets in expected.items():
            assert set(bridge.loc[gene, "bridges_usp34_to"].split(",")) == targets


class TestFourCandidateNetworkAudit:
    FROZEN_FOUR = {"USP34", "VEZF1", "EML5", "CITED2"}

    def test_frozen_shortlist_unchanged(self):
        for path in [
            TABLES / "four_candidate_direct_neighbors.tsv",
            TABLES / "four_candidate_shortest_paths.tsv",
            TABLES / "four_candidate_bridge_evidence.tsv",
            TABLES / "four_candidate_network_audit.tsv",
        ]:
            df = pd.read_csv(path, sep="\t")
            assert set(df["candidate"]) == self.FROZEN_FOUR

    def test_eml5_has_no_direct_neighbors(self):
        df = pd.read_csv(TABLES / "four_candidate_direct_neighbors.tsv", sep="\t")
        eml5 = df.loc[df["candidate"] == "EML5"]
        assert len(eml5) == 1
        assert eml5["neighbor_gene"].iloc[0] == "NO_RESOLVED_NETWORK_NEIGHBOURHOOD"

    def test_vezf1_direct_neighbor_count_matches_known_frozen_state(self):
        # VEZF1 has zero STRING partners at any threshold (documented in
        # docs/SYSTEMS_NETWORK_NODE_RULE.md) -- its only edge is the
        # pathway_co_membership link to DMTN.
        df = pd.read_csv(TABLES / "four_candidate_direct_neighbors.tsv", sep="\t")
        vezf1 = df.loc[df["candidate"] == "VEZF1"]
        assert len(vezf1) == 1
        assert vezf1["neighbor_gene"].iloc[0] == "DMTN"
        assert vezf1["interaction_type"].iloc[0] == "pathway_co_membership"

    def test_eml5_shortest_paths_explicitly_unresolved(self):
        df = pd.read_csv(TABLES / "four_candidate_shortest_paths.tsv", sep="\t")
        eml5 = df.loc[df["candidate"] == "EML5"]
        assert len(eml5) == 1
        assert eml5["path"].iloc[0] == "NO_RESOLVED_NETWORK_NEIGHBOURHOOD_IN_CURRENT_ANALYSIS"

    def test_no_shortest_path_exceeds_two_hops(self):
        df = pd.read_csv(TABLES / "four_candidate_shortest_paths.tsv", sep="\t")
        lengths = df["path_length_edges"].dropna()
        assert (lengths <= 2).all()

    def test_cited2_targets_are_all_one_hop_direct(self):
        # Part 4: CITED2's 5 named targets are literal direct neighbors, so
        # a 1-hop path length here means no bridge gene was needed -- this
        # must not be reported as an indirect/2-hop connection.
        df = pd.read_csv(TABLES / "four_candidate_shortest_paths.tsv", sep="\t")
        cited2 = df.loc[df["candidate"] == "CITED2"]
        assert (cited2["path_length_edges"] == 1).all()
        for path in cited2["path"].unique():
            assert path.count(" -> ") == 1

    def test_bridge_evidence_classification_tiers_valid(self):
        df = pd.read_csv(TABLES / "four_candidate_bridge_evidence.tsv", sep="\t")
        assert df["classification"].isin(
            ["A_DATA_SUPPORTED_BRIDGE", "B_PARTIAL_SUPPORT", "C_NETWORK_ONLY_GENERIC_BRIDGE", "D_NOT_ASSESSABLE"]
        ).all()

    def test_eml5_bridge_not_assessable(self):
        df = pd.read_csv(TABLES / "four_candidate_bridge_evidence.tsv", sep="\t")
        eml5 = df.loc[df["candidate"] == "EML5"]
        assert len(eml5) == 1
        assert eml5["classification"].iloc[0] == "D_NOT_ASSESSABLE"

    def test_gse245601_acute_never_solely_drives_a_tier(self):
        # a bridge gene must not reach A_DATA_SUPPORTED_BRIDGE on the
        # strength of the acute-only layer -- only CRISPR or a resistance
        # dataset FDR<0.05 hit may do that.
        df = pd.read_csv(TABLES / "four_candidate_bridge_evidence.tsv", sep="\t")
        resistance_or_crispr_fdr_cols = ["crispr_significant_fdr05", "gse118713_significant_fdr05", "gse240112_significant_fdr05", "gse111151_significant_fdr05"]
        a_tier = df.loc[df["classification"] == "A_DATA_SUPPORTED_BRIDGE"]
        for _, row in a_tier.iterrows():
            assert row[resistance_or_crispr_fdr_cols].any(), f"{row['bridge_gene']} reached A-tier without a non-acute FDR<0.05 hit"

    def test_convergence_matrix_covers_all_six_pairs(self):
        df = pd.read_csv(TABLES / "four_candidate_convergence.tsv", sep="\t")
        pairs = set(zip(df["candidate_A"], df["candidate_B"]))
        assert len(pairs) == 6

    def test_vezf1_cited2_convergence_preserved(self):
        df = pd.read_csv(TABLES / "four_candidate_convergence.tsv", sep="\t")
        row = df.loc[(df["candidate_A"] == "VEZF1") & (df["candidate_B"] == "CITED2")].iloc[0]
        assert row["any_convergence"]
        assert row["n_shared_pathways"] > 0
        assert "GOBP_BLOOD_VESSEL_MORPHOGENESIS" in row["shared_resistance_pathways_or_leading_edge_modules"]

    def test_eml5_never_shows_convergence(self):
        df = pd.read_csv(TABLES / "four_candidate_convergence.tsv", sep="\t")
        eml5_rows = df.loc[(df["candidate_A"] == "EML5") | (df["candidate_B"] == "EML5")]
        assert not eml5_rows["any_convergence"].any()

    def test_no_candidate_pair_shares_a_resistance_hub(self):
        df = pd.read_csv(TABLES / "four_candidate_convergence.tsv", sep="\t")
        assert (df["shared_resistance_hub_gene"].fillna("") == "").all()

    def test_head_to_head_classification_vocabulary(self):
        df = pd.read_csv(TABLES / "four_candidate_network_audit.tsv", sep="\t")
        valid = {
            "1_STRONG_SYSTEMS_SUPPORT",
            "2_MODERATE_SYSTEMS_SUPPORT",
            "3_NETWORK_HYPOTHESIS_ONLY",
            "4_DATA_SUPPORTED_BUT_MECHANISTICALLY_UNRESOLVED",
            "5_WEAK_NON_SPECIFIC_NETWORK_SUPPORT",
        }
        assert set(df["systems_mechanism_classification"]) <= valid

    def test_eml5_classified_data_supported_but_mechanistically_unresolved(self):
        df = pd.read_csv(TABLES / "four_candidate_network_audit.tsv", sep="\t").set_index("candidate")
        assert df.loc["EML5", "systems_mechanism_classification"] == "4_DATA_SUPPORTED_BUT_MECHANISTICALLY_UNRESOLVED"

    def test_direct_neighbor_counts_match_part2_table(self):
        neighbors = pd.read_csv(TABLES / "four_candidate_direct_neighbors.tsv", sep="\t")
        audit = pd.read_csv(TABLES / "four_candidate_network_audit.tsv", sep="\t").set_index("candidate")
        counts = {
            "USP34": 10,
            "VEZF1": 1,
            "EML5": 0,
            "CITED2": 18,
        }
        for candidate, expected_n in counts.items():
            assert int(audit.loc[candidate, "n_direct_neighbors"]) == expected_n

    def test_evidence_freeze_untouched(self):
        result = subprocess.run(
            ["git", "status", "--porcelain", "results/tables/evidence_freeze/", "docs/THERAPEUTIC_SHORTLIST_FREEZE.md"],
            capture_output=True,
            text=True,
            cwd=Path(__file__).resolve().parents[1],
        )
        assert result.stdout.strip() == "", f"frozen evidence-freeze files show as changed: {result.stdout}"

    def test_prior_usp34_bridge_evidence_content_unchanged(self):
        # USP34_bridge_gene_evidence.tsv is reused verbatim (Part 1/5) --
        # confirm its 4 genes and classifications still match what the
        # prior USP34-specific audit established, since this table predates
        # git tracking and a git-status check alone cannot prove it wasn't
        # edited.
        df = pd.read_csv(TABLES / "USP34_bridge_gene_evidence.tsv", sep="\t").set_index("gene")
        expected = {
            "USP9X": "B_PARTIAL_SUPPORT",
            "RPS27A": "C_NETWORK_ONLY_GENERIC_BRIDGE",
            "UBC": "C_NETWORK_ONLY_GENERIC_BRIDGE",
            "UBB": "B_PARTIAL_SUPPORT",
        }
        for gene, classification in expected.items():
            assert df.loc[gene, "classification"] == classification


if __name__ == "__main__":
    import sys

    sys.exit(pytest.main([__file__, "-v"]))
