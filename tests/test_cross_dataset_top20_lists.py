from pathlib import Path

import pandas as pd
import pytest

from src.cross_dataset_top20_lists import (
    build_top20_crispr_direction,
    build_top20_global,
    build_top20_multimodal,
    build_top20_resistance_consensus,
)

REPO_ROOT = Path(__file__).parent.parent


def _ranked():
    genes = [f"G{i}" for i in range(1, 26)]
    return pd.DataFrame(
        {
            "gene": genes,
            "global_rank": range(1, 26),
            "coverage_tier": ["A"] * 25,
            "crispr_evidence_percentile": [0.9] * 25,
            "gse118713_evidence_percentile": [0.8] * 25,
            "gse245601_evidence_percentile": [0.7] * 25,
            "gse240112_evidence_percentile": [0.6] * 25,
            "gse111151_evidence_percentile": [0.5] * 25,
            "n_datasets_fdr05": [2] * 25,
            "n_datasets_top10pct": [1] * 25,
            "median_evidence_percentile": [0.7] * 25,
            "n_datasets_testable": [5] * 25,
        }
    )


def _resistance_consensus(genes):
    return pd.DataFrame({"gene": genes, "resistance_direction_consensus": ["all_up"] * len(genes)})


def _categories(genes):
    cats = ["MULTIMODAL_STRONG" if i < 3 else "LOW_EVIDENCE" for i in range(len(genes))]
    return pd.DataFrame({"gene": genes, "evidence_category": cats})


class TestBuildTop20Global:
    def test_returns_exactly_20_when_available(self):
        ranked = _ranked()
        resistance = _resistance_consensus(ranked["gene"])
        categories = _categories(ranked["gene"])
        wide = pd.DataFrame({"gene": ranked["gene"], "gse245601_one_track_only": [False] * 25, "gse240112_outlier_fragility": [False] * 25})
        out = build_top20_global(ranked, resistance, categories, wide)
        assert len(out) == 20
        assert out["rank"].tolist() == list(range(1, 21))

    def test_caveat_flags_borderline_coverage(self):
        ranked = _ranked()
        ranked.loc[0, "n_datasets_testable"] = 3
        resistance = _resistance_consensus(ranked["gene"])
        categories = _categories(ranked["gene"])
        wide = pd.DataFrame({"gene": ranked["gene"], "gse245601_one_track_only": [False] * 25, "gse240112_outlier_fragility": [False] * 25})
        out = build_top20_global(ranked, resistance, categories, wide)
        assert "borderline coverage" in out.iloc[0]["main_caveat"]

    def test_no_caveat_is_explicit_none_not_blank(self):
        ranked = _ranked()
        resistance = _resistance_consensus(ranked["gene"])
        categories = _categories(ranked["gene"])
        wide = pd.DataFrame({"gene": ranked["gene"], "gse245601_one_track_only": [False] * 25, "gse240112_outlier_fragility": [False] * 25})
        out = build_top20_global(ranked, resistance, categories, wide)
        assert (out["main_caveat"] == "none").all()


class TestBuildTop20Multimodal:
    def test_only_multimodal_strong_genes_included(self):
        ranked = _ranked()
        categories = _categories(ranked["gene"])
        out = build_top20_multimodal(ranked, categories)
        assert len(out) == 3  # only 3 genes are MULTIMODAL_STRONG in the fixture
        assert set(out["gene"]) == {"G1", "G2", "G3"}


class TestBuildTop20ResistanceConsensus:
    def test_head_20(self):
        rc = pd.DataFrame({"gene": [f"G{i}" for i in range(30)], "resistance_fdr05_count": range(30)})
        out = build_top20_resistance_consensus(rc)
        assert len(out) == 20


class TestBuildTop20CrisprDirection:
    def test_filters_by_direction(self):
        cf = pd.DataFrame({"gene": ["A", "B", "C"], "crispr_direction": ["sensitising_KO", "tolerance_associated_KO", "sensitising_KO"], "crispr_fdr": [0.01, 0.01, 0.02]})
        out = build_top20_crispr_direction(cf, "sensitising_KO")
        assert set(out["gene"]) == {"A", "C"}


class TestRealData:
    def test_real_top20_global_if_present(self):
        path = REPO_ROOT / "results" / "tables" / "cross_dataset_genomewide" / "top20_global.tsv"
        if not path.exists():
            pytest.skip("top20_global table not present in this checkout")
        out = pd.read_csv(path, sep="\t")
        assert len(out) == 20
        assert not out["gene"].duplicated().any()
        assert out["rank"].tolist() == list(range(1, 21))

    def test_real_top20_rna_only_unique_genes_if_present(self):
        path = REPO_ROOT / "results" / "tables" / "cross_dataset_genomewide" / "top20_rna_only.tsv"
        if not path.exists():
            pytest.skip("top20_rna_only table not present in this checkout")
        out = pd.read_csv(path, sep="\t")
        assert len(out) == 20
        assert not out["gene"].duplicated().any()

    def test_real_all_other_top20_lists_have_unique_genes_if_present(self):
        names = [
            "top20_multimodal.tsv", "top20_resistance_consensus.tsv", "top20_crispr_sensitising.tsv",
            "top20_crispr_tolerance.tsv", "top20_human_tumor.tsv", "top20_crispr_nonsignificant_rna.tsv",
        ]
        base = REPO_ROOT / "results" / "tables" / "cross_dataset_genomewide"
        found_any = False
        for name in names:
            path = base / name
            if not path.exists():
                continue
            found_any = True
            out = pd.read_csv(path, sep="\t")
            assert not out["gene"].duplicated().any(), f"{name} has duplicate genes"
        if not found_any:
            pytest.skip("no additional top-20 tables present in this checkout")
