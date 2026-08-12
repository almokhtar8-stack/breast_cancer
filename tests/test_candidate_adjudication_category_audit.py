from pathlib import Path

import pandas as pd

from src.candidate_adjudication_category_audit import find_near_misses, reconstruct_multimodal_strong

REPO_ROOT = Path(__file__).parent.parent
EXPECTED_SEVEN = {"USP34", "VEZF1", "CUX1", "DPP9", "LZTR1", "SOX2", "TFAP2C"}


def _synthetic_wide(genes, crispr_fdr, n_fdr05_overall):
    return pd.DataFrame(
        {
            "gene": genes,
            "coverage_tier": ["A"] * len(genes),
            "crispr_fdr": crispr_fdr,
            "n_datasets_fdr05": n_fdr05_overall,
        }
    )


def _synthetic_resistance(genes, fdr05_count, consensus):
    return pd.DataFrame({"gene": genes, "resistance_fdr05_count": fdr05_count, "resistance_direction_consensus": consensus})


def _synthetic_wide_raw(rows: list[dict]) -> pd.DataFrame:
    """Builds a wide-table-shaped fixture with the raw per-dataset
    testable/FDR columns `reconstruct_multimodal_strong` now recomputes
    from directly (Phase 34 Codex-review fix -- it no longer trusts the
    frozen table's own precomputed n_datasets_fdr05/resistance_fdr05_count
    aggregate columns). Each row dict may set: gene, coverage_tier,
    crispr_fdr, gse118713_fdr, gse240112_tumor_fdr, gse111151_fdr,
    gse245601_epi_fdr, and any `*_testable` flag (defaults: all testable)."""
    cols = ["crispr", "gse118713", "gse245601", "gse240112", "gse111151"]
    out_rows = []
    for r in rows:
        row = {
            "gene": r["gene"], "coverage_tier": r.get("coverage_tier", "A"),
            "crispr_fdr": r.get("crispr_fdr"), "gse118713_fdr": r.get("gse118713_fdr", 0.9),
            "gse245601_epi_fdr": r.get("gse245601_epi_fdr", 0.9), "gse240112_tumor_fdr": r.get("gse240112_tumor_fdr", 0.9),
            "gse111151_fdr": r.get("gse111151_fdr", 0.9),
        }
        for c in cols:
            row[f"{c}_testable"] = r.get(f"{c}_testable", True)
        out_rows.append(row)
    return pd.DataFrame(out_rows)


class TestReconstructMultimodalStrong:
    def test_reproduces_real_frozen_seven_genes_exactly(self):
        cdx = REPO_ROOT / "results" / "tables" / "cross_dataset_genomewide"
        if not (cdx / "all_genes_cross_dataset_evidence_with_ranking.tsv").exists():
            return
        wide = pd.read_csv(cdx / "all_genes_cross_dataset_evidence_with_ranking.tsv", sep="\t")
        resistance = pd.read_csv(cdx / "resistance_consensus_all_genes.tsv", sep="\t")
        reconstructed = reconstruct_multimodal_strong(wide, resistance)
        assert set(reconstructed["gene"]) == EXPECTED_SEVEN

    def test_two_resistance_datasets_fdr05_qualifies(self):
        wide = _synthetic_wide_raw([{"gene": "G1", "crispr_fdr": 0.05, "gse118713_fdr": 0.01, "gse240112_tumor_fdr": 0.01, "gse111151_fdr": 0.5}])
        resistance = _synthetic_resistance(["G1"], [0], ["all_up"])  # deliberately wrong precomputed count -- must be ignored
        out = reconstruct_multimodal_strong(wide, resistance)
        assert list(out["gene"]) == ["G1"]

    def test_one_resistance_dataset_plus_overall_two_qualifies(self):
        wide = _synthetic_wide_raw([{"gene": "G1", "crispr_fdr": 0.03, "gse118713_fdr": 0.01, "gse240112_tumor_fdr": 0.5, "gse111151_fdr": 0.5}])
        resistance = _synthetic_resistance(["G1"], [0], ["all_up"])
        out = reconstruct_multimodal_strong(wide, resistance)
        assert list(out["gene"]) == ["G1"]

    def test_one_resistance_dataset_alone_does_not_qualify(self):
        """CRISPR FDR=0.08 is under the 0.10 category gate but NOT under
        0.05, so it contributes nothing to the overall FDR<0.05 count;
        with only 1 resistance dataset significant and nothing else,
        neither cond_a (>=2 resistance) nor cond_b (>=1 resistance AND
        >=2 overall) is satisfied."""
        wide = _synthetic_wide_raw([{"gene": "G1", "crispr_fdr": 0.08, "gse118713_fdr": 0.01, "gse240112_tumor_fdr": 0.5, "gse111151_fdr": 0.5}])
        resistance = _synthetic_resistance(["G1"], [0], ["majority_up"])
        out = reconstruct_multimodal_strong(wide, resistance)
        assert len(out) == 0

    def test_crispr_fdr_above_010_never_qualifies_regardless_of_resistance(self):
        wide = _synthetic_wide_raw([{"gene": "G1", "crispr_fdr": 0.15, "gse118713_fdr": 0.01, "gse240112_tumor_fdr": 0.01, "gse111151_fdr": 0.01}])
        resistance = _synthetic_resistance(["G1"], [0], ["all_up"])
        out = reconstruct_multimodal_strong(wide, resistance)
        assert len(out) == 0

    def test_low_coverage_gene_never_qualifies_even_with_strong_fdrs(self):
        """A gene testable in only 2/5 datasets must be excluded by the
        LOW_COVERAGE gate even if both its testable FDRs are extremely
        significant (recomputed directly from the raw *_testable flags,
        not trusted from a precomputed coverage_tier/n_datasets_testable)."""
        wide = _synthetic_wide_raw(
            [{"gene": "G1", "crispr_fdr": 0.001, "gse118713_fdr": 0.001, "gse240112_tumor_fdr": 0.5, "gse111151_fdr": 0.5,
              "gse245601_testable": False, "gse240112_testable": False, "gse111151_testable": False}]
        )
        resistance = _synthetic_resistance(["G1"], [0], ["all_up"])
        out = reconstruct_multimodal_strong(wide, resistance)
        assert len(out) == 0

    def test_no_gene_identity_special_cased(self):
        """The reconstruction must depend only on numeric columns -- run it
        on a gene named 'USP34' with weak evidence and a gene named
        'FAKE999' with strong evidence, and confirm only the evidence
        -strong one qualifies, regardless of name."""
        wide = _synthetic_wide_raw(
            [
                {"gene": "USP34", "crispr_fdr": 0.5, "gse118713_fdr": 0.9, "gse240112_tumor_fdr": 0.9, "gse111151_fdr": 0.9},
                {"gene": "FAKE999", "crispr_fdr": 0.02, "gse118713_fdr": 0.01, "gse240112_tumor_fdr": 0.01, "gse111151_fdr": 0.9},
            ]
        )
        resistance = _synthetic_resistance(["USP34", "FAKE999"], [0, 0], ["mixed", "all_down"])
        out = reconstruct_multimodal_strong(wide, resistance)
        assert list(out["gene"]) == ["FAKE999"]


class TestFindNearMisses:
    def test_near_miss_boundary_detected(self):
        wide = _synthetic_wide(["G1"], [0.078], [1])
        resistance = _synthetic_resistance(["G1"], [1], ["majority_up"])
        out = find_near_misses(wide, resistance, multimodal_genes=set())
        assert "G1" in set(out["gene"])
        assert out.loc[out["gene"] == "G1", "missing_condition"].iloc[0].startswith("has_CRISPR_FDR<0.10")

    def test_far_from_threshold_gene_not_a_near_miss(self):
        wide = _synthetic_wide(["G1"], [0.9], [0])
        resistance = _synthetic_resistance(["G1"], [0], ["insufficient"])
        out = find_near_misses(wide, resistance, multimodal_genes=set())
        assert "G1" not in set(out["gene"])
