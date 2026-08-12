from pathlib import Path

import pandas as pd
import pytest

from src.cross_dataset_surprise_discovery import build_fallen_genes, build_surprise_discovery

REPO_ROOT = Path(__file__).parent.parent


def _ranked():
    genes = [f"G{i}" for i in range(1, 26)]
    return pd.DataFrame({"gene": genes, "global_rank": range(1, 26)})


class TestBuildSurpriseDiscovery:
    def test_flags_genes_outside_old28(self):
        ranked = _ranked()
        old28 = {"G1", "G2"}
        out = build_surprise_discovery(ranked, old28, top_n=20)
        assert len(out) == 20
        assert out.set_index("gene").loc["G1", "is_surprise"] == False  # noqa: E712
        assert out.set_index("gene").loc["G3", "is_surprise"] == True  # noqa: E712

    def test_previously_known_flag_is_inverse_of_surprise(self):
        ranked = _ranked()
        old28 = {"G1"}
        out = build_surprise_discovery(ranked, old28, top_n=10)
        assert (out["is_surprise"] != out["previously_in_28_crispr_hits"]).all()


class TestBuildFallenGenes:
    def test_gene_ranked_far_down_flagged(self):
        ranked = pd.DataFrame({"gene": ["A", "B"], "global_rank": [5, 5000]})
        full = pd.DataFrame({"gene": ["A", "B"], "n_datasets_testable": [5, 5]})
        out = build_fallen_genes(ranked, {"A", "B"}, full, rank_threshold=1000)
        assert list(out["gene"]) == ["B"]
        assert out.iloc[0]["status"] == "ranked_but_far_down"

    def test_gene_below_coverage_threshold_flagged(self):
        ranked = pd.DataFrame({"gene": ["A"], "global_rank": [5]})
        full = pd.DataFrame({"gene": ["A", "C"], "n_datasets_testable": [5, 1]})
        out = build_fallen_genes(ranked, {"A", "C"}, full, rank_threshold=1000)
        assert list(out["gene"]) == ["C"]
        assert "below_coverage_threshold_1_of_5" in out.iloc[0]["status"]

    def test_gene_not_in_universe_at_all(self):
        ranked = pd.DataFrame({"gene": ["A"], "global_rank": [5]})
        full = pd.DataFrame({"gene": ["A"], "n_datasets_testable": [5]})
        out = build_fallen_genes(ranked, {"A", "MISSING"}, full, rank_threshold=1000)
        assert list(out["gene"]) == ["MISSING"]
        assert out.iloc[0]["status"] == "not_in_universe"

    def test_well_ranked_gene_not_flagged(self):
        ranked = pd.DataFrame({"gene": ["A"], "global_rank": [5]})
        full = pd.DataFrame({"gene": ["A"], "n_datasets_testable": [5]})
        out = build_fallen_genes(ranked, {"A"}, full, rank_threshold=1000)
        assert len(out) == 0


class TestRealData:
    def test_real_surprise_discovery_if_present(self):
        path = REPO_ROOT / "results" / "tables" / "cross_dataset_genomewide" / "surprise_discovery_top20.tsv"
        if not path.exists():
            pytest.skip("surprise discovery table not present in this checkout")
        out = pd.read_csv(path, sep="\t")
        assert len(out) == 20
        # per the analysis, the overwhelming majority of the unbiased top20 were not in the old CRISPR-first candidate set
        assert out["is_surprise"].sum() >= 15
