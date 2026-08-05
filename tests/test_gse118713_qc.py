from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from src.gse118713_qc import (
    compute_pca,
    compute_sample_correlations,
    compute_sample_summary,
    plot_pca,
    plot_sample_correlation,
    select_top_variable_genes,
    write_pca_coordinates,
    write_sample_correlations,
)

SAMPLE_IDS = [
    "MCF7_Rep1", "MCF7_Rep2", "MCF7_Rep3",
    "TAMR_Rep1", "TAMR_Rep2", "TAMR_Rep3",
    "FASR_Rep1", "FASR_Rep2", "FASR_Rep3",
]
GROUPS = ["MCF7"] * 3 + ["TAMR"] * 3 + ["FASR"] * 3


def _meta_df() -> pd.DataFrame:
    return pd.DataFrame({"sample_id": SAMPLE_IDS, "group": GROUPS})


def _flat_gene_df(n=10, seed=0) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    data = {"gene_id": [f"G{i}" for i in range(n)], "gene_symbol": [f"S{i}" for i in range(n)]}
    for sample_id in SAMPLE_IDS:
        data[sample_id] = rng.uniform(1, 50, n)
    return pd.DataFrame(data)


class TestSampleSummary:
    def test_all_nine_samples_present(self):
        df = _flat_gene_df()
        meta = _meta_df()
        summary = compute_sample_summary(df, meta, SAMPLE_IDS)
        assert list(summary["sample_id"]) == SAMPLE_IDS

    def test_no_automatic_sample_removal_even_with_degenerate_sample(self):
        df = _flat_gene_df()
        df["MCF7_Rep1"] = 0.0  # a fully-zero sample -- must still appear in QC output
        meta = _meta_df()
        summary = compute_sample_summary(df, meta, SAMPLE_IDS)
        assert "MCF7_Rep1" in set(summary["sample_id"])
        assert len(summary) == 9

    def test_summary_values_match_manual_computation(self):
        df = _flat_gene_df()
        meta = _meta_df()
        summary = compute_sample_summary(df, meta, SAMPLE_IDS)
        row = summary[summary["sample_id"] == "MCF7_Rep1"].iloc[0]
        values = df["MCF7_Rep1"].to_numpy()
        assert row["total_tpm"] == pytest.approx(values.sum())
        assert row["n_genes_tpm_ge_1"] == int((values >= 1).sum())
        assert row["median_tpm"] == pytest.approx(np.median(values))


class TestSelectTopVariableGenes:
    def test_selects_requested_count(self):
        df = _flat_gene_df(n=20)
        top = select_top_variable_genes(df, SAMPLE_IDS, n_genes=5)
        assert len(top) == 5

    def test_deterministic_tie_break_by_gene_id(self):
        # Two genes with identical values across all samples -> tied variance (zero).
        data = {"gene_id": ["G2", "G1"], "gene_symbol": ["s2", "s1"]}
        for sample_id in SAMPLE_IDS:
            data[sample_id] = [3.0, 3.0]
        df = pd.DataFrame(data)
        top = select_top_variable_genes(df, SAMPLE_IDS, n_genes=1)
        assert list(top["gene_id"]) == ["G1"]  # lower gene_id wins the tie

    def test_repeated_calls_are_identical(self):
        df = _flat_gene_df(n=30)
        top1 = select_top_variable_genes(df, SAMPLE_IDS, n_genes=10)
        top2 = select_top_variable_genes(df, SAMPLE_IDS, n_genes=10)
        pd.testing.assert_frame_equal(top1, top2)

    def test_highest_variance_gene_included(self):
        low_vals = [10.0] * 9
        high_vals = [1.0, 100.0, 1.0, 100.0, 1.0, 100.0, 1.0, 100.0, 1.0]
        df = pd.DataFrame({"gene_id": ["LOW", "HIGH"], "gene_symbol": ["l", "h"]})
        for i, sample_id in enumerate(SAMPLE_IDS):
            df[sample_id] = [low_vals[i], high_vals[i]]
        top = select_top_variable_genes(df, SAMPLE_IDS, n_genes=1)
        assert list(top["gene_id"]) == ["HIGH"]


class TestSampleCorrelations:
    def test_self_correlation_is_one(self):
        df = _flat_gene_df()
        corr = compute_sample_correlations(df, SAMPLE_IDS)
        for sample_id in SAMPLE_IDS:
            r = corr[(corr["sample_id_1"] == sample_id) & (corr["sample_id_2"] == sample_id)]["spearman_r"].iloc[0]
            assert r == pytest.approx(1.0)

    def test_symmetric(self):
        df = _flat_gene_df()
        corr = compute_sample_correlations(df, SAMPLE_IDS)
        a, b = SAMPLE_IDS[0], SAMPLE_IDS[1]
        r_ab = corr[(corr["sample_id_1"] == a) & (corr["sample_id_2"] == b)]["spearman_r"].iloc[0]
        r_ba = corr[(corr["sample_id_1"] == b) & (corr["sample_id_2"] == a)]["spearman_r"].iloc[0]
        assert r_ab == pytest.approx(r_ba)


class TestPca:
    def test_variance_explained_sums_to_one(self):
        df = _flat_gene_df(n=15)
        meta = _meta_df()
        pca_df = compute_pca(df, meta, SAMPLE_IDS)
        variances = pca_df.drop_duplicates("pc")["variance_explained_fraction"]
        assert variances.sum() == pytest.approx(1.0, abs=1e-6)

    def test_all_samples_present(self):
        df = _flat_gene_df(n=15)
        meta = _meta_df()
        pca_df = compute_pca(df, meta, SAMPLE_IDS)
        assert set(pca_df["sample_id"]) == set(SAMPLE_IDS)

    def test_groups_separate_along_pc1_when_group_shift_dominates(self):
        rng = np.random.default_rng(1)
        n = 50
        data = {"gene_id": [f"G{i}" for i in range(n)], "gene_symbol": [f"S{i}" for i in range(n)]}
        base = rng.uniform(5, 50, n)
        for sample_id, group in zip(SAMPLE_IDS, GROUPS):
            shift = {"MCF7": 0, "TAMR": 20, "FASR": -10}[group]
            data[sample_id] = np.clip(base + shift, 0.01, None)
        df = pd.DataFrame(data)
        meta = _meta_df()
        pca_df = compute_pca(df, meta, SAMPLE_IDS)
        wide = pca_df.pivot(index=["sample_id", "group"], columns="pc", values="coordinate").reset_index()
        pc1_by_group = wide.groupby("group")["PC1"].mean()
        # TAMR and FASR should be well separated on PC1 given the dominant group shift.
        assert abs(pc1_by_group["TAMR"] - pc1_by_group["FASR"]) > 1.0


class TestFiguresFromSavedTables:
    def test_plot_pca_reads_only_from_file(self, tmp_path: Path):
        df = _flat_gene_df(n=15)
        meta = _meta_df()
        pca_df = compute_pca(df, meta, SAMPLE_IDS)
        pca_tsv = tmp_path / "pca.tsv"
        pca_df.to_csv(pca_tsv, sep="\t", index=False)

        out_pdf = tmp_path / "pca.pdf"
        plot_pca(pca_tsv, out_pdf)  # only a path is passed in, not the DataFrame
        assert out_pdf.exists()
        assert out_pdf.stat().st_size > 0

    def test_plot_sample_correlation_reads_only_from_file(self, tmp_path: Path):
        df = _flat_gene_df()
        corr_df = compute_sample_correlations(df, SAMPLE_IDS)
        corr_tsv = tmp_path / "corr.tsv"
        corr_df.to_csv(corr_tsv, sep="\t", index=False)

        out_pdf = tmp_path / "corr.pdf"
        plot_sample_correlation(corr_tsv, out_pdf)
        assert out_pdf.exists()
        assert out_pdf.stat().st_size > 0

    def test_figure_regenerable_from_table_alone_after_process_restart(self, tmp_path: Path):
        # Simulates a fresh process: write the table, then only ever touch the
        # path (never an in-memory object) to produce the figure.
        df = _flat_gene_df()
        corr_df = compute_sample_correlations(df, SAMPLE_IDS)
        cfg_like_path = tmp_path / "gse118713_sample_correlations.tsv"
        write_sample_correlations(corr_df, type("Cfg", (), {"sample_correlations_tsv": cfg_like_path})())
        del corr_df
        out_pdf = tmp_path / "corr.pdf"
        plot_sample_correlation(cfg_like_path, out_pdf)
        assert out_pdf.exists()


class TestFigureDeterminism:
    # Frozen figures must be byte-reproducible, not merely visually equal --
    # matplotlib's PDF backend embeds a CreationDate by default, which would
    # defeat SHA256-based freeze verification even though the plotted
    # content is unchanged.
    def test_pca_pdf_byte_identical_across_writes(self, tmp_path: Path):
        import hashlib
        import time

        df = _flat_gene_df(n=15)
        meta = _meta_df()
        pca_df = compute_pca(df, meta, SAMPLE_IDS)
        pca_tsv = tmp_path / "pca.tsv"
        pca_df.to_csv(pca_tsv, sep="\t", index=False)

        out1 = tmp_path / "pca1.pdf"
        out2 = tmp_path / "pca2.pdf"
        plot_pca(pca_tsv, out1)
        time.sleep(1.1)  # force a different wall-clock second between writes
        plot_pca(pca_tsv, out2)

        h1 = hashlib.sha256(out1.read_bytes()).hexdigest()
        h2 = hashlib.sha256(out2.read_bytes()).hexdigest()
        assert h1 == h2

    def test_correlation_pdf_byte_identical_across_writes(self, tmp_path: Path):
        import hashlib
        import time

        df = _flat_gene_df()
        corr_df = compute_sample_correlations(df, SAMPLE_IDS)
        corr_tsv = tmp_path / "corr.tsv"
        corr_df.to_csv(corr_tsv, sep="\t", index=False)

        out1 = tmp_path / "corr1.pdf"
        out2 = tmp_path / "corr2.pdf"
        plot_sample_correlation(corr_tsv, out1)
        time.sleep(1.1)  # force a different wall-clock second between writes
        plot_sample_correlation(corr_tsv, out2)

        h1 = hashlib.sha256(out1.read_bytes()).hexdigest()
        h2 = hashlib.sha256(out2.read_bytes()).hexdigest()
        assert h1 == h2
