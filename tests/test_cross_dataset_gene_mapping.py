from pathlib import Path

import pandas as pd
import pytest

from src.cross_dataset_gene_mapping import (
    build_full_gene_universe,
    build_gene_mapping_audit,
    flag_ambiguous_symbols,
    load_crispr_screen,
    load_gse111151,
    load_gse118713,
    load_gse245601_track,
    load_gse240112_track,
    resolve_ensgr_duplicates,
)

REPO_ROOT = Path(__file__).parent.parent


class TestResolveEnsgrDuplicates:
    def test_identical_stats_ensgr_pair_resolved(self):
        df = pd.DataFrame(
            {
                "ensembl_id": ["ENSG00000002586", "ENSGR0000002586"],
                "symbol": ["CD99", "CD99"],
                "effect": [-0.56, -0.56],
                "p_value": [0.002, 0.002],
                "fdr": [0.012, 0.012],
            }
        )
        out, resolved = resolve_ensgr_duplicates(df, "ensembl_id", "symbol")
        assert len(out) == 1
        assert out["ensembl_id"].iloc[0] == "ENSG00000002586"
        assert resolved == ["CD99"]

    def test_different_stats_not_resolved(self):
        df = pd.DataFrame(
            {
                "ensembl_id": ["ENSG00000002586", "ENSGR0000002586"],
                "symbol": ["CD99", "CD99"],
                "effect": [-0.56, 0.10],
                "p_value": [0.002, 0.5],
                "fdr": [0.012, 0.7],
            }
        )
        out, resolved = resolve_ensgr_duplicates(df, "ensembl_id", "symbol")
        assert len(out) == 2
        assert resolved == []

    def test_no_ensgr_rows_unchanged(self):
        df = pd.DataFrame({"ensembl_id": ["ENSG1", "ENSG2"], "symbol": ["A", "B"], "effect": [1.0, 2.0], "p_value": [0.1, 0.2], "fdr": [0.1, 0.2]})
        out, resolved = resolve_ensgr_duplicates(df, "ensembl_id", "symbol")
        assert len(out) == 2
        assert resolved == []


class TestFlagAmbiguousSymbols:
    def test_duplicate_symbol_excluded_and_flagged(self):
        df = pd.DataFrame({"symbol": ["A", "A", "B"], "effect": [1.0, 2.0, 3.0]})
        clean, ambiguous = flag_ambiguous_symbols(df, "symbol")
        assert list(clean["symbol"]) == ["B"]
        assert ambiguous == ["A"]

    def test_no_duplicates_nothing_removed(self):
        df = pd.DataFrame({"symbol": ["A", "B"], "effect": [1.0, 2.0]})
        clean, ambiguous = flag_ambiguous_symbols(df, "symbol")
        assert len(clean) == 2
        assert ambiguous == []


class TestBuildFullGeneUniverse:
    def test_union_not_intersection(self):
        datasets = {
            "d1": pd.DataFrame({"symbol": ["A", "B"]}),
            "d2": pd.DataFrame({"symbol": ["B", "C"]}),
        }
        out = build_full_gene_universe(datasets).set_index("gene")
        assert set(out.index) == {"A", "B", "C"}
        assert out.loc["A", "n_datasets_present"] == 1
        assert out.loc["B", "n_datasets_present"] == 2
        assert out.loc["C", "n_datasets_present"] == 1

    def test_not_seeded_by_any_named_gene_list(self):
        # the universe must include a gene even if it appears in only one obscure dataset
        # and is absent from every other -- no allowlist gates membership
        datasets = {"d1": pd.DataFrame({"symbol": ["OBSCURE1"]}), "d2": pd.DataFrame({"symbol": []})}
        out = build_full_gene_universe(datasets)
        assert "OBSCURE1" in set(out["gene"])

    def test_flags_present_and_testable_consistent(self):
        datasets = {"d1": pd.DataFrame({"symbol": ["A"]})}
        out = build_full_gene_universe(datasets).set_index("gene")
        assert bool(out.loc["A", "d1_present"])
        assert bool(out.loc["A", "d1_testable"])


class TestBuildGeneMappingAudit:
    def test_all_categories_represented(self):
        out = build_gene_mapping_audit(["ENSGR1"], ["AMB1"], ["ENSGR2"], ["AMB2"])
        assert len(out) == 4
        assert set(out["issue"]) == {"pseudoautosomal_ENSGR_duplicate", "duplicate_symbol_multiple_distinct_ensembl_ids"}
        assert out.set_index("symbol").loc["AMB1", "excluded_from_testable"]
        assert not out.set_index("symbol").loc["ENSGR1", "excluded_from_testable"]


class TestRealDataLoaders:
    def test_crispr_loader_if_present(self):
        path = REPO_ROOT / "data" / "processed" / "labels.parquet"
        if not path.exists():
            pytest.skip("frozen CRISPR labels not present in this checkout")
        out = load_crispr_screen(path)
        assert len(out) == 19103
        assert not out["symbol"].duplicated().any()
        assert out[["effect", "p_value", "fdr"]].notna().all().all()

    def test_gse118713_loader_resolves_known_ensgr_pair_if_present(self):
        path = REPO_ROOT / "results" / "tables" / "gse118713_differential_expression_unredacted.tsv.gz"
        if not path.exists():
            pytest.skip("frozen GSE118713 unredacted DE table not present in this checkout")
        clean, ensgr_resolved, ambiguous = load_gse118713(path)
        assert "CD99" in ensgr_resolved
        assert not clean["symbol"].duplicated().any()
        assert len(ambiguous) > 0

    def test_gse118713_loader_uses_only_primary_contrast_not_secondary(self):
        # TAMR_vs_FASR and FASR_vs_MCF7 must never leak into this dataset's one
        # independent contribution -- confirmed by row count matching the
        # primary-contrast-only subset, not the full 3-contrasts-stacked file
        # (ambiguous-symbol rows are excluded from `clean`, and some ambiguous
        # symbols occur more than twice, so row counts are compared as
        # bounds/ratios rather than an exact clean+ambiguous reconstruction)
        path = REPO_ROOT / "results" / "tables" / "gse118713_differential_expression_unredacted.tsv.gz"
        if not path.exists():
            pytest.skip("frozen GSE118713 unredacted DE table not present in this checkout")
        full_file = pd.read_csv(path, sep="\t")
        assert set(full_file["contrast"].unique()) == {"TAMR_vs_MCF7", "TAMR_vs_FASR", "FASR_vs_MCF7"}
        primary_only_count = (full_file["contrast"] == "TAMR_vs_MCF7").sum()
        clean, ensgr_resolved, ambiguous = load_gse118713(path)
        # clean must be smaller than the full 3-contrast file (proves filtering happened)
        # but close to (not equal to, due to ENSGR/ambiguous exclusions) the single-contrast count
        assert len(clean) < len(full_file)
        assert primary_only_count - len(clean) < 200  # only ENSGR (15) + ambiguous-symbol rows (a few dozen) excluded, not thousands
        assert len(ensgr_resolved) == 15
        assert len(ambiguous) == 56

    def test_gse111151_loader_if_present(self):
        path = REPO_ROOT / "results" / "tables" / "gse111151" / "genomewide_de.tsv.gz"
        if not path.exists():
            pytest.skip("frozen GSE111151 DE table not present in this checkout")
        clean, ensgr_resolved, ambiguous = load_gse111151(path)
        assert not clean["symbol"].duplicated().any()

    def test_gse245601_and_gse240112_loaders_if_present(self):
        base = REPO_ROOT / "results" / "tables"
        a = base / "gse245601_pseudobulk" / "track_a_genomewide_de.tsv.gz"
        t = base / "gse240112_pseudobulk" / "tumor_cell_genomewide_de.tsv.gz"
        if not a.exists() or not t.exists():
            pytest.skip("frozen GSE245601/GSE240112 DE tables not present in this checkout")
        out_a = load_gse245601_track(a)
        out_t = load_gse240112_track(t)
        assert not out_a["symbol"].duplicated().any()
        assert not out_t["symbol"].duplicated().any()


class TestRealPipeline:
    def test_universe_and_audit_if_present(self):
        universe_path = REPO_ROOT / "results" / "tables" / "cross_dataset_genomewide" / "full_gene_universe.tsv"
        audit_path = REPO_ROOT / "results" / "tables" / "cross_dataset_genomewide" / "gene_mapping_audit.tsv"
        if not universe_path.exists() or not audit_path.exists():
            pytest.skip("cross-dataset gene universe not present in this checkout")
        universe = pd.read_csv(universe_path, sep="\t")
        assert not universe["gene"].duplicated().any()
        # universe must be far larger than any previously-discussed candidate list (13 or 28 genes)
        assert len(universe) > 1000
        audit = pd.read_csv(audit_path, sep="\t")
        assert len(audit) > 0
