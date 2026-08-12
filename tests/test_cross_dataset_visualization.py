from pathlib import Path

import pandas as pd
import pytest

from src.cross_dataset_visualization import (
    plot_all_gene_evidence_heatmap,
    plot_coverage_vs_evidence,
    plot_evidence_category_counts,
    plot_surprise_candidates,
)

REPO_ROOT = Path(__file__).parent.parent


def _ranked(n=30):
    return pd.DataFrame(
        {
            "gene": [f"G{i}" for i in range(n)],
            "crispr_evidence_percentile": [0.9 - i * 0.01 for i in range(n)],
            "gse118713_evidence_percentile": [0.8] * n,
            "gse245601_evidence_percentile": [0.7] * n,
            "gse240112_evidence_percentile": [0.6] * n,
            "gse111151_evidence_percentile": [0.5] * n,
            "median_evidence_percentile": [0.7 - i * 0.01 for i in range(n)],
        }
    )


class TestPlotSurpriseCandidates:
    def test_excludes_old28_genes(self, tmp_path):
        ranked = _ranked()
        old28 = {"G0", "G1", "G2"}
        out_path = tmp_path / "surprise.png"
        plot_surprise_candidates(ranked, old28, out_path, n=5)
        assert out_path.exists()

    def test_all_old28_still_produces_valid_plot(self, tmp_path):
        ranked = _ranked(n=5)
        old28 = set(ranked["gene"])
        out_path = tmp_path / "surprise_empty.png"
        # every gene excluded -> surprises is empty; must not crash
        plot_surprise_candidates(ranked, old28, out_path, n=5)
        assert out_path.exists()


class TestPlotAllGeneEvidenceHeatmap:
    def test_produces_file(self, tmp_path):
        ranked = _ranked()
        out_path = tmp_path / "heatmap.png"
        plot_all_gene_evidence_heatmap(ranked, out_path, n_genes=10)
        assert out_path.exists()


class TestPlotEvidenceCategoryCounts:
    def test_produces_file(self, tmp_path):
        categories = pd.DataFrame({"gene": ["A", "B", "C"], "evidence_category": ["MULTIMODAL_STRONG", "LOW_EVIDENCE", "LOW_EVIDENCE"]})
        out_path = tmp_path / "categories.png"
        plot_evidence_category_counts(categories, out_path)
        assert out_path.exists()


class TestPlotCoverageVsEvidence:
    def test_produces_file(self, tmp_path):
        full = pd.DataFrame({"coverage_tier": ["A", "B", "C", "D", "E"] * 4, "equal_dataset_mean_percentile": [0.5] * 20})
        out_path = tmp_path / "coverage.png"
        plot_coverage_vs_evidence(full, out_path)
        assert out_path.exists()


class TestRealFigures:
    def test_final_review_figures_present_if_generated(self):
        final_review = REPO_ROOT / "results" / "figures" / "cross_dataset_genomewide" / "final_review"
        if not final_review.exists() or not any(final_review.iterdir()):
            pytest.skip("cross-dataset final_review figures not present in this checkout")
        pngs = list(final_review.glob("*.png"))
        assert len(pngs) == 12
