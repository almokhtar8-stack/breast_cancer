import hashlib
import time

import pandas as pd
import pytest

from src.gse118713_tamr_specificity import (
    SpecificityConfig,
    build_specificity_table,
    load_de_table,
    write_specificity_table,
)


def _de_row(gene_id, gene_symbol, contrast, log2fc, fdr, p_value=0.01):
    return {
        "gene_id": gene_id,
        "gene_symbol": gene_symbol,
        "log2fc": log2fc,
        "se": 0.1,
        "moderated_t": log2fc / 0.1,
        "p_value": p_value,
        "fdr": fdr,
        "ave_expr": 5.0,
        "contrast": contrast,
        "direction": "up" if log2fc > 0 else "down",
    }


class TestBuildSpecificityTable:
    def test_joins_all_three_contrasts_by_gene(self):
        rows = []
        for contrast in ("TAMR_vs_MCF7", "FASR_vs_MCF7", "TAMR_vs_FASR"):
            rows.append(_de_row("G1", "SYM1", contrast, 1.0, 0.01))
        de_df = pd.DataFrame(rows)
        out = build_specificity_table(de_df)
        assert len(out) == 1
        assert out.loc[0, "tamr_vs_mcf7_log2fc"] == 1.0
        assert out.loc[0, "fasr_vs_mcf7_log2fc"] == 1.0
        assert out.loc[0, "tamr_vs_fasr_log2fc"] == 1.0

    def test_direct_tamr_vs_fasr_is_the_specificity_score(self):
        rows = [
            _de_row("G1", "SYM1", "TAMR_vs_MCF7", 2.0, 0.2),
            _de_row("G1", "SYM1", "FASR_vs_MCF7", 0.05, 0.9),  # nonsignificant, near zero
            _de_row("G1", "SYM1", "TAMR_vs_FASR", 1.95, 0.001),  # highly significant, large
        ]
        de_df = pd.DataFrame(rows)
        out = build_specificity_table(de_df)
        assert out.loc[0, "tamr_vs_fasr_log2fc"] == pytest.approx(1.95)
        assert out.loc[0, "tamr_vs_fasr_fdr"] == pytest.approx(0.001)

    def test_no_false_specificity_from_fasr_nonsignificance(self):
        # FASR_vs_MCF7 nonsignificant (fdr=0.9) but TAMR_vs_FASR also
        # nonsignificant (fdr=0.8) -- a naive "FASR is NS so TAMR is
        # specific" rule would wrongly call this gene specific. The table
        # must report the real (nonsignificant) TAMR_vs_FASR value, not a
        # derived specificity label based on FASR alone.
        rows = [
            _de_row("G1", "SYM1", "TAMR_vs_MCF7", 0.5, 0.85),
            _de_row("G1", "SYM1", "FASR_vs_MCF7", 0.02, 0.90),
            _de_row("G1", "SYM1", "TAMR_vs_FASR", 0.48, 0.80),
        ]
        de_df = pd.DataFrame(rows)
        out = build_specificity_table(de_df)
        assert out.loc[0, "tamr_vs_fasr_fdr"] == pytest.approx(0.80)
        # No column encodes a specificity verdict derived only from FASR.
        forbidden_cols = {"is_specific", "specificity_label", "fasr_nonsignificant"}
        assert forbidden_cols.isdisjoint(set(out.columns))

    def test_direction_concordance_flag(self):
        rows_concordant = [
            _de_row("G1", "SYM1", "TAMR_vs_MCF7", 2.0, 0.01),
            _de_row("G1", "SYM1", "FASR_vs_MCF7", -1.0, 0.01),
            _de_row("G1", "SYM1", "TAMR_vs_FASR", 3.0, 0.01),
        ]
        rows_discordant = [
            _de_row("G2", "SYM2", "TAMR_vs_MCF7", 2.0, 0.01),
            _de_row("G2", "SYM2", "FASR_vs_MCF7", -1.0, 0.01),
            _de_row("G2", "SYM2", "TAMR_vs_FASR", -3.0, 0.01),
        ]
        de_df = pd.DataFrame(rows_concordant + rows_discordant)
        out = build_specificity_table(de_df).set_index("gene_id")
        assert bool(out.loc["G1", "same_direction_tamr_vs_mcf7_and_tamr_vs_fasr"]) is True
        assert bool(out.loc["G2", "same_direction_tamr_vs_mcf7_and_tamr_vs_fasr"]) is False

    def test_mismatched_gene_sets_across_contrasts_raises(self):
        rows = [
            _de_row("G1", "SYM1", "TAMR_vs_MCF7", 1.0, 0.01),
            _de_row("G1", "SYM1", "FASR_vs_MCF7", 1.0, 0.01),
            _de_row("G2", "SYM2", "TAMR_vs_FASR", 1.0, 0.01),  # different gene entirely
        ]
        de_df = pd.DataFrame(rows)
        with pytest.raises(ValueError, match="gene set does not match"):
            build_specificity_table(de_df)

    def test_duplicate_gene_within_contrast_raises(self):
        rows = [
            _de_row("G1", "SYM1", "TAMR_vs_MCF7", 1.0, 0.01),
            _de_row("G1", "SYM1", "TAMR_vs_MCF7", 2.0, 0.02),  # duplicate
            _de_row("G1", "SYM1", "FASR_vs_MCF7", 1.0, 0.01),
            _de_row("G1", "SYM1", "TAMR_vs_FASR", 1.0, 0.01),
        ]
        de_df = pd.DataFrame(rows)
        with pytest.raises(ValueError, match="duplicate gene_id"):
            build_specificity_table(de_df)

    def test_output_sorted_by_tamr_vs_fasr_fdr(self):
        rows = []
        for gene, fdr in (("G1", 0.5), ("G2", 0.01), ("G3", 0.2)):
            rows.append(_de_row(gene, gene, "TAMR_vs_MCF7", 1.0, 0.1))
            rows.append(_de_row(gene, gene, "FASR_vs_MCF7", 1.0, 0.1))
            rows.append(_de_row(gene, gene, "TAMR_vs_FASR", 1.0, fdr))
        de_df = pd.DataFrame(rows)
        out = build_specificity_table(de_df)
        assert list(out["gene_id"]) == ["G2", "G3", "G1"]


class TestWriteSpecificityTableDeterminism:
    def test_gzip_output_byte_identical_across_writes(self, tmp_path):
        # Frozen outputs must be byte-reproducible, not merely numerically
        # equal -- a gzip header embedding the write timestamp would defeat
        # SHA256-based freeze verification even though the content is
        # unchanged.
        rows = []
        for contrast in ("TAMR_vs_MCF7", "FASR_vs_MCF7", "TAMR_vs_FASR"):
            rows.append(_de_row("G1", "SYM1", contrast, 1.0, 0.01))
        specificity_df = build_specificity_table(pd.DataFrame(rows))

        cfg1 = SpecificityConfig(
            differential_expression_tsv_gz=tmp_path / "unused.tsv.gz",
            output_tsv_gz=tmp_path / "run1" / "out.tsv.gz",
        )
        cfg2 = SpecificityConfig(
            differential_expression_tsv_gz=tmp_path / "unused.tsv.gz",
            output_tsv_gz=tmp_path / "run2" / "out.tsv.gz",
        )
        write_specificity_table(specificity_df, cfg1)
        time.sleep(1.1)  # force a different wall-clock second between writes
        write_specificity_table(specificity_df, cfg2)

        h1 = hashlib.sha256(cfg1.output_tsv_gz.read_bytes()).hexdigest()
        h2 = hashlib.sha256(cfg2.output_tsv_gz.read_bytes()).hexdigest()
        assert h1 == h2


class TestLoadDeTable:
    def test_rejects_missing_contrast(self, tmp_path):
        rows = [_de_row("G1", "SYM1", "TAMR_vs_MCF7", 1.0, 0.01)]
        path = tmp_path / "de.tsv.gz"
        pd.DataFrame(rows).to_csv(path, sep="\t", index=False, compression="gzip")
        with pytest.raises(ValueError, match="missing required contrasts"):
            load_de_table(path)

    def test_rejects_missing_columns(self, tmp_path):
        path = tmp_path / "de.tsv.gz"
        pd.DataFrame({"gene_id": ["G1"], "contrast": ["TAMR_vs_MCF7"]}).to_csv(
            path, sep="\t", index=False, compression="gzip"
        )
        with pytest.raises(ValueError, match="missing required columns"):
            load_de_table(path)
