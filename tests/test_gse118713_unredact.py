import pandas as pd
import pytest

from src.gse118713_unredact import NUMERIC_ATOL, UnredactionConfig, compare_redacted_vs_unredacted


def _de_row(gene_id, gene_symbol, contrast, log2fc=1.0, se=0.1, p_value=0.01, fdr=0.02, ave_expr=5.0):
    return {
        "gene_id": gene_id,
        "gene_symbol": gene_symbol,
        "log2fc": log2fc,
        "se": se,
        "moderated_t": log2fc / se,
        "p_value": p_value,
        "fdr": fdr,
        "ave_expr": ave_expr,
        "contrast": contrast,
        "direction": "up" if log2fc > 0 else "down",
    }


CONTRASTS = ("TAMR_vs_MCF7", "FASR_vs_MCF7", "TAMR_vs_FASR")


def _write_de_table(path, gene_specs):
    """gene_specs: list of (gene_id, gene_symbol, kwargs) present in every contrast."""
    rows = []
    for gene_id, gene_symbol, kwargs in gene_specs:
        for contrast in CONTRASTS:
            rows.append(_de_row(gene_id, gene_symbol, contrast, **kwargs))
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False, compression="gzip")


def _cfg(tmp_path, old_path, new_path, blinded_gene_ids=("BLIND1", "BLIND2")):
    return UnredactionConfig(
        filtered_gene_tpm_tsv=tmp_path / "unused_filtered.tsv.gz",
        sample_metadata_tsv=tmp_path / "unused_meta.tsv",
        limma_script=tmp_path / "unused.R",
        old_de_tsv_gz=old_path,
        new_de_tsv_gz=new_path,
        new_redaction_record_tsv=tmp_path / "unused_record.tsv",
        new_specificity_tsv_gz=tmp_path / "unused_specificity.tsv.gz",
        comparison_tsv=tmp_path / "comparison.tsv",
        blinded_gene_ids=blinded_gene_ids,
    )


class TestCompareRedactedVsUnredacted:
    def test_identical_shared_rows_report_zero_diffs_and_correct_added_genes(self, tmp_path):
        old_path = tmp_path / "old.tsv.gz"
        new_path = tmp_path / "new.tsv.gz"
        shared = [("G1", "SYM1", {}), ("G2", "SYM2", {"log2fc": -0.5})]
        _write_de_table(old_path, shared)
        _write_de_table(new_path, shared + [("BLIND1", "RCOR1", {}), ("BLIND2", "KDM1A", {})])

        cfg = _cfg(tmp_path, old_path, new_path)
        out = compare_redacted_vs_unredacted(cfg)

        assert len(out) == 3  # one row per contrast
        for _, row in out.iterrows():
            assert row["n_old_genes"] == 2
            assert row["n_new_genes"] == 4
            assert row["n_shared_genes"] == 2
            assert row["n_added_genes"] == 2
            assert set(row["added_gene_ids"].split(",")) == {"BLIND1", "BLIND2"}
            assert row["n_fdr_changed_for_shared_genes"] == 0
            assert row["max_abs_diff_log2fc"] == pytest.approx(0.0)
            assert row["max_abs_diff_fdr"] == pytest.approx(0.0)
        assert cfg.comparison_tsv.exists()

    def test_tiny_floating_point_noise_within_tolerance_passes(self, tmp_path):
        old_path = tmp_path / "old.tsv.gz"
        new_path = tmp_path / "new.tsv.gz"
        _write_de_table(old_path, [("G1", "SYM1", {"fdr": 0.02})])
        _write_de_table(
            new_path,
            [
                ("G1", "SYM1", {"fdr": 0.02 + NUMERIC_ATOL / 10}),
                ("BLIND1", "RCOR1", {}),
                ("BLIND2", "KDM1A", {}),
            ],
        )
        cfg = _cfg(tmp_path, old_path, new_path)
        out = compare_redacted_vs_unredacted(cfg)
        assert (out["n_fdr_changed_for_shared_genes"] == 0).all()

    def test_a_real_fdr_change_beyond_tolerance_raises(self, tmp_path):
        old_path = tmp_path / "old.tsv.gz"
        new_path = tmp_path / "new.tsv.gz"
        _write_de_table(old_path, [("G1", "SYM1", {"fdr": 0.02})])
        _write_de_table(
            new_path,
            [
                ("G1", "SYM1", {"fdr": 0.05}),  # real change, not noise
                ("BLIND1", "RCOR1", {}),
                ("BLIND2", "KDM1A", {}),
            ],
        )
        cfg = _cfg(tmp_path, old_path, new_path)
        with pytest.raises(ValueError, match="reproducibility tolerance"):
            compare_redacted_vs_unredacted(cfg)

    def test_missing_previously_reported_gene_raises(self, tmp_path):
        old_path = tmp_path / "old.tsv.gz"
        new_path = tmp_path / "new.tsv.gz"
        _write_de_table(old_path, [("G1", "SYM1", {}), ("G2", "SYM2", {})])
        # G2 silently dropped in the "new" table -- must never happen
        _write_de_table(new_path, [("G1", "SYM1", {}), ("BLIND1", "RCOR1", {}), ("BLIND2", "KDM1A", {})])
        cfg = _cfg(tmp_path, old_path, new_path)
        with pytest.raises(ValueError, match="missing from"):
            compare_redacted_vs_unredacted(cfg)

    def test_unexpected_extra_gene_raises(self, tmp_path):
        old_path = tmp_path / "old.tsv.gz"
        new_path = tmp_path / "new.tsv.gz"
        _write_de_table(old_path, [("G1", "SYM1", {})])
        # An extra gene beyond the two configured blind IDs appears -- not allowed
        _write_de_table(
            new_path,
            [("G1", "SYM1", {}), ("BLIND1", "RCOR1", {}), ("BLIND2", "KDM1A", {}), ("EXTRA", "MYSTERY", {})],
        )
        cfg = _cfg(tmp_path, old_path, new_path)
        with pytest.raises(ValueError, match="expected exactly"):
            compare_redacted_vs_unredacted(cfg)

    def test_partial_unblinding_raises(self, tmp_path):
        old_path = tmp_path / "old.tsv.gz"
        new_path = tmp_path / "new.tsv.gz"
        _write_de_table(old_path, [("G1", "SYM1", {})])
        # only one of the two configured blind genes appears
        _write_de_table(new_path, [("G1", "SYM1", {}), ("BLIND1", "RCOR1", {})])
        cfg = _cfg(tmp_path, old_path, new_path)
        with pytest.raises(ValueError, match="expected exactly"):
            compare_redacted_vs_unredacted(cfg)

    def test_gene_symbol_mismatch_for_shared_gene_raises(self, tmp_path):
        old_path = tmp_path / "old.tsv.gz"
        new_path = tmp_path / "new.tsv.gz"
        _write_de_table(old_path, [("G1", "SYM1", {})])
        _write_de_table(new_path, [("G1", "DIFFERENT_SYMBOL", {}), ("BLIND1", "RCOR1", {}), ("BLIND2", "KDM1A", {})])
        cfg = _cfg(tmp_path, old_path, new_path)
        with pytest.raises(ValueError, match="gene_symbol"):
            compare_redacted_vs_unredacted(cfg)


class TestRealUnredactionOutputs:
    """Light checks against the committed unredaction outputs, if present.

    Skipped rather than failed when the outputs don't exist, since
    generating them requires an R/limma environment this suite does not
    assume is always available.
    """

    def test_committed_comparison_shows_exactly_rcor1_and_kdm1a_added(self):
        import yaml
        from pathlib import Path

        config_path = Path("config/config.yaml")
        if not config_path.exists():
            pytest.skip("no config.yaml in this environment")
        with open(config_path) as f:
            config = yaml.safe_load(f)
        comparison_path = Path(config["gse118713_phase2b"]["unredaction"]["comparison_tsv"])
        if not comparison_path.exists():
            pytest.skip("unredaction has not been run in this environment")

        comparison = pd.read_csv(comparison_path, sep="\t")
        blinded_ids = set(config["gse118713_phase2b"]["blinding"]["blinded_gene_ids"])
        assert len(comparison) == 3  # one row per contrast
        for _, row in comparison.iterrows():
            assert set(row["added_gene_ids"].split(",")) == blinded_ids
            assert row["n_fdr_changed_for_shared_genes"] == 0
            assert row["n_added_genes"] == 2
