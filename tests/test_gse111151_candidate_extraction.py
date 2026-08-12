from pathlib import Path

import pandas as pd
import pytest

from src.gse111151_candidate_extraction import (
    benjamini_hochberg,
    build_candidate_table,
    build_paics_row,
    build_sample_level_table,
    load_de_table,
)

REPO_ROOT = Path(__file__).parent.parent

ENSEMBL_IDS = {"USP34": "ENSG00000001", "CTDNEP1": "ENSG00000002", "VEZF1": "ENSG00000003", "PAICS": "ENSG00000004"}


def _synthetic_de():
    return pd.DataFrame(
        {
            "gene_id": ["ENSG00000001", "ENSG00000002", "ENSG00000003"],
            "gene_name": ["USP34", "CTDNEP1", "VEZF1"],
            "log2fc": [0.4, -0.1, 1.1],
            "avg_log_cpm": [7.0, 5.8, 5.5],
            "p_value": [0.12, 0.70, 0.004],
            "fdr": [0.23, 0.80, 0.02],
        }
    )


def _synthetic_counts_meta():
    genes = ["ENSG00000001", "ENSG00000002", "ENSG00000003"]
    counts = pd.DataFrame(
        {
            "MCF-7": [100, 200, 50], "MCF-7_Tam1": [150, 190, 200],
            "T-47D": [110, 210, 55], "T-47D_Tam1": [160, 195, 210], "T-47D_Tam2": [140, 185, 190],
        },
        index=genes,
    )
    counts.index.name = "gene_id"
    metadata = pd.DataFrame(
        {
            "sample_id": ["MCF-7", "MCF-7_Tam1", "T-47D", "T-47D_Tam1", "T-47D_Tam2"],
            "cell_line": ["MCF-7", "MCF-7", "T-47D", "T-47D", "T-47D"],
            "resistance_status": ["parental", "resistant", "parental", "resistant", "resistant"],
        }
    )
    return counts, metadata


class TestBenjaminiHochberg:
    def test_matches_known_example(self):
        p = pd.Series([0.01, 0.04, 0.03, 0.005])
        out = benjamini_hochberg(p)
        assert out.iloc[3] == pytest.approx(0.02)


class TestLoadDeTable:
    def test_raises_on_missing_columns(self, tmp_path):
        p = tmp_path / "de.tsv"
        pd.DataFrame({"gene_id": ["A"], "log2fc": [1.0]}).to_csv(p, sep="\t", index=False)
        with pytest.raises(ValueError, match="missing required columns"):
            load_de_table(p)


class TestBuildCandidateTable:
    def test_untested_gene_marked_with_reason_not_dropped(self):
        de = _synthetic_de()
        counts, metadata = _synthetic_counts_meta()
        eids = dict(ENSEMBL_IDS, USP17L29="ENSG00000099")
        out = build_candidate_table(de, counts, metadata, ["USP34", "USP17L29"], eids).set_index("gene")
        assert not out.loc["USP17L29", "tested"]
        assert "gene_id absent" in out.loc["USP17L29", "reason_not_tested"]
        assert pd.isna(out.loc["USP17L29", "candidate_set_bh_fdr"])

    def test_candidate_set_bh_excludes_untested(self):
        de = _synthetic_de()
        counts, metadata = _synthetic_counts_meta()
        out = build_candidate_table(de, counts, metadata, ["USP34", "CTDNEP1", "VEZF1"], ENSEMBL_IDS).set_index("gene")
        assert out.loc["VEZF1", "candidate_set_bh_fdr"] == pytest.approx(0.004 * 3 / 1)
        assert out.loc["USP34", "candidate_set_bh_fdr"] == pytest.approx(0.12 * 3 / 2)

    def test_direction_matches_log2fc_sign(self):
        de = _synthetic_de()
        counts, metadata = _synthetic_counts_meta()
        out = build_candidate_table(de, counts, metadata, ["USP34", "CTDNEP1"], ENSEMBL_IDS).set_index("gene")
        assert out.loc["USP34", "direction"] == "up_in_resistant"
        assert out.loc["CTDNEP1", "direction"] == "down_in_resistant"

    def test_lookup_uses_ensembl_id_not_symbol(self):
        # gene_name column in DE table is deliberately irrelevant to lookup -- only gene_id matters
        de = _synthetic_de()
        de["gene_name"] = "SOMEOTHERSYMBOL"
        counts, metadata = _synthetic_counts_meta()
        out = build_candidate_table(de, counts, metadata, ["USP34"], ENSEMBL_IDS).set_index("gene")
        assert out.loc["USP34", "tested"]
        assert out.loc["USP34", "log2fc"] == pytest.approx(0.4)


class TestSampleLevelTable:
    def test_one_row_per_gene_per_sample(self):
        counts, metadata = _synthetic_counts_meta()
        out = build_sample_level_table(counts, metadata, ["USP34", "VEZF1"], ENSEMBL_IDS)
        assert len(out) == 2 * 5
        assert set(out["sample_id"]) == set(metadata["sample_id"])

    def test_gene_absent_from_counts_skipped(self):
        counts, metadata = _synthetic_counts_meta()
        eids = dict(ENSEMBL_IDS, USP17L29="ENSG00000099")
        out = build_sample_level_table(counts, metadata, ["USP34", "USP17L29"], eids)
        assert "USP17L29" not in set(out["gene"])
        assert "USP34" in set(out["gene"])


class TestPaicsRow:
    def test_paics_labeled_as_benchmark(self):
        de = _synthetic_de()
        de.loc[len(de)] = ["ENSG00000004", "PAICS", 0.13, 7.2, 0.69, 0.79]
        counts, metadata = _synthetic_counts_meta()
        counts.loc["ENSG00000004"] = [500, 600, 510, 590, 610]
        out = build_paics_row(de, counts, metadata, "PAICS", "ENSG00000004")
        assert out.iloc[0]["benchmark_label"] == "published_benchmark_not_in_13_candidate_bh_family"
        assert out.iloc[0]["tested"]


class TestRealData:
    def test_real_candidate_table_if_present(self):
        table_path = REPO_ROOT / "results" / "tables" / "gse111151" / "candidate_table.tsv"
        if not table_path.exists():
            pytest.skip("GSE111151 candidate table not present in this checkout")
        out = pd.read_csv(table_path, sep="\t")
        assert len(out) == 13
        assert set(out["gene"]) == {
            "USP34", "CTDNEP1", "EIF4ENIF1", "HMGB1", "KDM1A", "PET117", "TADA2B", "VEZF1", "ICK", "SUPT4H1", "TLK2", "TSR3", "USP17L29",
        }
        tested_fdrs = out.loc[out["tested"], "candidate_set_bh_fdr"]
        assert tested_fdrs.notna().all()
