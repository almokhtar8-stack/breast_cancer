"""Integration tests for scripts/analysis/gse118713_limma.R.

Runs the real R script (bc environment: R 4.5.3, limma 3.66.0, statmod
1.5.2) as a subprocess against small synthetic fixtures, and checks the
statistical properties PREANALYSIS.md's Phase 2B amendment requires:
correct unpaired three-group design, correct contrast signs, the
TAMR_vs_FASR-equals-coefficient-difference identity, per-contrast BH
correction, deterministic sorting, and clear failure reporting for bad
input. No DESeq2/edgeR usage is confirmed by a static source check.
"""

from __future__ import annotations

import gzip
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from statsmodels.stats.multitest import multipletests

REPO_ROOT = Path(__file__).parent.parent
SCRIPT = REPO_ROOT / "scripts" / "analysis" / "gse118713_limma.R"
LIB = REPO_ROOT / "scripts" / "analysis" / "gse118713_limma_lib.R"

SAMPLE_IDS = [
    "MCF7_Rep1", "MCF7_Rep2", "MCF7_Rep3",
    "TAMR_Rep1", "TAMR_Rep2", "TAMR_Rep3",
    "FASR_Rep1", "FASR_Rep2", "FASR_Rep3",
]
GROUPS = ["MCF7"] * 3 + ["TAMR"] * 3 + ["FASR"] * 3

pytestmark = pytest.mark.skipif(shutil.which("Rscript") is None, reason="Rscript not available")


def _write_meta(path: Path, sample_ids=SAMPLE_IDS, groups=GROUPS, replicates=None) -> None:
    replicates = replicates or [1, 2, 3] * 3
    pd.DataFrame({"sample_id": sample_ids, "group": groups, "replicate": replicates}).to_csv(
        path, sep="\t", index=False
    )


def _write_expr(path: Path, values: dict[str, list[float]], gene_ids=None, gene_symbols=None) -> None:
    n = len(next(iter(values.values())))
    gene_ids = gene_ids or [f"G{i}" for i in range(n)]
    gene_symbols = gene_symbols or [f"S{i}" for i in range(n)]
    df = pd.DataFrame({"gene_id": gene_ids, "gene_symbol": gene_symbols, **values})
    df.to_csv(path, sep="\t", index=False, compression="gzip")


def _synthetic_expression(n_genes=40, seed=0) -> dict[str, list[float]]:
    rng = np.random.default_rng(seed)
    base = rng.uniform(5, 200, n_genes)
    values: dict[str, list[float]] = {sid: [] for sid in SAMPLE_IDS}
    tamr_shift = rng.normal(3, 1, n_genes)
    fasr_shift = rng.normal(-2, 1, n_genes)
    for sample_id, group in zip(SAMPLE_IDS, GROUPS):
        noise = rng.normal(0, 0.3, n_genes)
        shift = {"MCF7": 0, "TAMR": tamr_shift, "FASR": fasr_shift}[group]
        values[sample_id] = list(np.clip(base + shift + noise, 0.01, None))
    return values


def _run(expr_path: Path, meta_path: Path, out_path: Path):
    # The R CLI requires the blinding arguments on every invocation (see
    # TestBlindedGeneRedaction) -- generic statistics-only tests pass an
    # explicit empty blinded-gene list rather than relying on any implicit
    # "no blinding" default.
    record_path = out_path.parent / f"{out_path.name}.redaction_record.tsv"
    return subprocess.run(
        ["Rscript", str(SCRIPT), str(expr_path), str(meta_path), str(out_path), "", str(record_path)],
        capture_output=True,
        text=True,
    )


def _read_output(out_path: Path) -> pd.DataFrame:
    with gzip.open(out_path, "rt") as f:
        return pd.read_csv(f, sep="\t")


class TestValidRun:
    def test_runs_successfully_and_produces_all_contrasts(self, tmp_path):
        expr_path = tmp_path / "expr.tsv.gz"
        meta_path = tmp_path / "meta.tsv"
        out_path = tmp_path / "out.tsv.gz"
        _write_expr(expr_path, _synthetic_expression())
        _write_meta(meta_path)

        result = _run(expr_path, meta_path, out_path)
        assert result.returncode == 0, result.stderr
        de = _read_output(out_path)
        assert set(de["contrast"].unique()) == {"TAMR_vs_MCF7", "FASR_vs_MCF7", "TAMR_vs_FASR"}
        assert set(de["gene_id"].unique()) == {f"G{i}" for i in range(40)}
        # every filtered gene must appear in every contrast -- no silent gene loss
        for contrast in de["contrast"].unique():
            assert len(de[de["contrast"] == contrast]) == 40

    def test_required_columns_present(self, tmp_path):
        expr_path = tmp_path / "expr.tsv.gz"
        meta_path = tmp_path / "meta.tsv"
        out_path = tmp_path / "out.tsv.gz"
        _write_expr(expr_path, _synthetic_expression())
        _write_meta(meta_path)
        _run(expr_path, meta_path, out_path)
        de = _read_output(out_path)
        expected_cols = {
            "gene_id", "gene_symbol", "log2fc", "se", "moderated_t",
            "p_value", "fdr", "ave_expr", "contrast", "direction",
        }
        assert expected_cols.issubset(set(de.columns))


class TestContrastSignsAndAlgebra:
    def test_contrast_signs_match_known_direction(self, tmp_path):
        # One gene, deterministically higher in TAMR, lower in FASR, vs MCF7.
        expr_path = tmp_path / "expr.tsv.gz"
        meta_path = tmp_path / "meta.tsv"
        out_path = tmp_path / "out.tsv.gz"
        values = {}
        for sample_id, group in zip(SAMPLE_IDS, GROUPS):
            level = {"MCF7": 10.0, "TAMR": 40.0, "FASR": 2.0}[group]
            values.setdefault(sample_id, []).append(level)
        _write_expr(expr_path, values, gene_ids=["G0"], gene_symbols=["S0"])
        _write_meta(meta_path)
        _run(expr_path, meta_path, out_path)
        de = _read_output(out_path).set_index("contrast")

        assert de.loc["TAMR_vs_MCF7", "log2fc"] > 0  # TAMR > MCF7
        assert de.loc["FASR_vs_MCF7", "log2fc"] < 0  # FASR < MCF7
        assert de.loc["TAMR_vs_FASR", "log2fc"] > 0  # TAMR > FASR

    def test_tamr_vs_fasr_equals_difference_of_coefficients(self, tmp_path):
        expr_path = tmp_path / "expr.tsv.gz"
        meta_path = tmp_path / "meta.tsv"
        out_path = tmp_path / "out.tsv.gz"
        _write_expr(expr_path, _synthetic_expression(n_genes=25, seed=2))
        _write_meta(meta_path)
        _run(expr_path, meta_path, out_path)
        de = _read_output(out_path)

        wide = de.pivot(index="gene_id", columns="contrast", values="log2fc")
        implied = wide["TAMR_vs_MCF7"] - wide["FASR_vs_MCF7"]
        assert np.allclose(implied.to_numpy(), wide["TAMR_vs_FASR"].to_numpy(), atol=1e-9)


class TestUnpairedThreeGroupDesign:
    def test_design_matrix_is_one_hot_by_group_not_replicate(self, tmp_path):
        r_code = f"""
        source("{LIB}")
        meta <- data.frame(
          sample_id = c({", ".join(f'"{s}"' for s in SAMPLE_IDS)}),
          group = c({", ".join(f'"{g}"' for g in GROUPS)})
        )
        design <- build_design(meta)
        cat(paste(colnames(design), collapse=","), "\\n")
        cat(paste(colSums(design), collapse=","), "\\n")
        cat(ncol(design), "\\n")
        """
        result = subprocess.run(["Rscript", "-e", r_code], capture_output=True, text=True)
        assert result.returncode == 0, result.stderr
        lines = result.stdout.strip().splitlines()
        columns = [c.strip() for c in lines[0].split(",")]
        col_sums = [float(x) for x in lines[1].split(",")]
        n_cols = int(lines[2])

        assert columns == ["MCF7", "TAMR", "FASR"]
        assert n_cols == 3  # no separate replicate/blocking term
        assert col_sums == [3.0, 3.0, 3.0]


class TestBhCorrectionSeparatePerContrast:
    def test_fdr_matches_independent_bh_per_contrast(self, tmp_path):
        expr_path = tmp_path / "expr.tsv.gz"
        meta_path = tmp_path / "meta.tsv"
        out_path = tmp_path / "out.tsv.gz"
        _write_expr(expr_path, _synthetic_expression(n_genes=60, seed=3))
        _write_meta(meta_path)
        _run(expr_path, meta_path, out_path)
        de = _read_output(out_path)

        for contrast, sub in de.groupby("contrast"):
            recomputed = multipletests(sub["p_value"].to_numpy(), method="fdr_bh")[1]
            assert np.allclose(recomputed, sub["fdr"].to_numpy(), atol=1e-9)


class TestDeterministicSorting:
    def test_sort_order_contrast_fdr_effect_gene_id(self, tmp_path):
        expr_path = tmp_path / "expr.tsv.gz"
        meta_path = tmp_path / "meta.tsv"
        out_path = tmp_path / "out.tsv.gz"
        _write_expr(expr_path, _synthetic_expression(n_genes=50, seed=4))
        _write_meta(meta_path)
        _run(expr_path, meta_path, out_path)
        de = _read_output(out_path)

        expected_contrast_order = ["TAMR_vs_MCF7", "FASR_vs_MCF7", "TAMR_vs_FASR"]
        observed_contrast_order = list(dict.fromkeys(de["contrast"]))
        assert observed_contrast_order == expected_contrast_order

        for _, sub in de.groupby("contrast", sort=False):
            fdrs = sub["fdr"].to_numpy()
            assert np.all(np.diff(fdrs) >= -1e-12)  # non-decreasing FDR within a contrast

    def test_repeated_runs_produce_identical_output(self, tmp_path):
        expr_path = tmp_path / "expr.tsv.gz"
        meta_path = tmp_path / "meta.tsv"
        out1 = tmp_path / "out1.tsv.gz"
        out2 = tmp_path / "out2.tsv.gz"
        _write_expr(expr_path, _synthetic_expression(n_genes=30, seed=5))
        _write_meta(meta_path)
        _run(expr_path, meta_path, out1)
        _run(expr_path, meta_path, out2)
        pd.testing.assert_frame_equal(_read_output(out1), _read_output(out2))

    def test_repeated_runs_produce_byte_identical_gzip(self, tmp_path):
        # Frozen outputs must be byte-reproducible, not merely numerically
        # equal -- a gzip header embedding the write timestamp would defeat
        # SHA256-based freeze verification even though the content is
        # unchanged. R's gzfile() writes mtime=0 by default; this locks
        # that behaviour in.
        import hashlib
        import time

        expr_path = tmp_path / "expr.tsv.gz"
        meta_path = tmp_path / "meta.tsv"
        out1 = tmp_path / "out1.tsv.gz"
        out2 = tmp_path / "out2.tsv.gz"
        _write_expr(expr_path, _synthetic_expression(n_genes=15, seed=6))
        _write_meta(meta_path)
        _run(expr_path, meta_path, out1)
        time.sleep(1.1)  # force a different wall-clock second between writes
        _run(expr_path, meta_path, out2)
        h1 = hashlib.sha256(out1.read_bytes()).hexdigest()
        h2 = hashlib.sha256(out2.read_bytes()).hexdigest()
        assert h1 == h2


class TestRejectsBadInput:
    def test_rejects_duplicated_sample_id(self, tmp_path):
        expr_path = tmp_path / "expr.tsv.gz"
        meta_path = tmp_path / "meta.tsv"
        out_path = tmp_path / "out.tsv.gz"
        _write_expr(expr_path, _synthetic_expression(n_genes=10))
        bad_ids = list(SAMPLE_IDS)
        bad_ids[1] = bad_ids[0]  # duplicate
        _write_meta(meta_path, sample_ids=bad_ids)
        result = _run(expr_path, meta_path, out_path)
        assert result.returncode != 0
        assert "duplicated sample_id" in result.stderr
        assert not out_path.exists()

    def test_rejects_duplicated_gene_id(self, tmp_path):
        expr_path = tmp_path / "expr.tsv.gz"
        meta_path = tmp_path / "meta.tsv"
        out_path = tmp_path / "out.tsv.gz"
        _write_expr(expr_path, _synthetic_expression(n_genes=10), gene_ids=["G0"] * 10)
        _write_meta(meta_path)
        result = _run(expr_path, meta_path, out_path)
        assert result.returncode != 0
        assert "duplicated gene_id" in result.stderr

    def test_rejects_negative_tpm(self, tmp_path):
        expr_path = tmp_path / "expr.tsv.gz"
        meta_path = tmp_path / "meta.tsv"
        out_path = tmp_path / "out.tsv.gz"
        values = _synthetic_expression(n_genes=5)
        values[SAMPLE_IDS[0]][0] = -5.0
        _write_expr(expr_path, values)
        _write_meta(meta_path)
        result = _run(expr_path, meta_path, out_path)
        assert result.returncode != 0
        assert "negative TPM" in result.stderr

    def test_rejects_non_finite_tpm(self, tmp_path):
        expr_path = tmp_path / "expr.tsv.gz"
        meta_path = tmp_path / "meta.tsv"
        out_path = tmp_path / "out.tsv.gz"
        values = _synthetic_expression(n_genes=5)
        values[SAMPLE_IDS[0]][0] = float("inf")
        _write_expr(expr_path, values)
        _write_meta(meta_path)
        result = _run(expr_path, meta_path, out_path)
        assert result.returncode != 0

    def test_rejects_wrong_replicate_count_per_group(self, tmp_path):
        expr_path = tmp_path / "expr.tsv.gz"
        meta_path = tmp_path / "meta.tsv"
        out_path = tmp_path / "out.tsv.gz"
        _write_expr(expr_path, _synthetic_expression(n_genes=10))
        bad_groups = list(GROUPS)
        bad_groups[0] = "TAMR"  # MCF7 now has 2, TAMR has 4
        _write_meta(meta_path, groups=bad_groups)
        result = _run(expr_path, meta_path, out_path)
        assert result.returncode != 0
        assert "replicates per group" in result.stderr

    def test_rejects_missing_sample_column_in_expression(self, tmp_path):
        expr_path = tmp_path / "expr.tsv.gz"
        meta_path = tmp_path / "meta.tsv"
        out_path = tmp_path / "out.tsv.gz"
        values = _synthetic_expression(n_genes=5)
        del values["FASR_Rep3"]
        _write_expr(expr_path, values)
        _write_meta(meta_path)
        result = _run(expr_path, meta_path, out_path)
        assert result.returncode != 0
        assert "missing sample columns" in result.stderr

    def test_wrong_number_of_arguments_reported_clearly(self, tmp_path):
        result = subprocess.run(["Rscript", str(SCRIPT), "one_arg"], capture_output=True, text=True)
        assert result.returncode != 0
        assert "usage" in result.stderr.lower() or "ERROR" in result.stderr


class TestNoCountBasedMethods:
    def test_no_deseq2_or_edger_reference_in_source(self):
        for path in (SCRIPT, LIB):
            text = path.read_text()
            for forbidden in ("DESeq2", "DESeqDataSet", "edgeR", "DGEList", "calcNormFactors"):
                assert forbidden not in text, f"{forbidden} referenced in {path}"


class TestNoHardcodedResults:
    def test_different_inputs_yield_different_outputs(self, tmp_path):
        expr_path_a = tmp_path / "a.tsv.gz"
        expr_path_b = tmp_path / "b.tsv.gz"
        meta_path = tmp_path / "meta.tsv"
        out_a = tmp_path / "out_a.tsv.gz"
        out_b = tmp_path / "out_b.tsv.gz"
        _write_expr(expr_path_a, _synthetic_expression(n_genes=20, seed=10))
        _write_expr(expr_path_b, _synthetic_expression(n_genes=20, seed=11))
        _write_meta(meta_path)
        _run(expr_path_a, meta_path, out_a)
        _run(expr_path_b, meta_path, out_b)
        de_a = _read_output(out_a).sort_values(["contrast", "gene_id"]).reset_index(drop=True)
        de_b = _read_output(out_b).sort_values(["contrast", "gene_id"]).reset_index(drop=True)
        assert not np.allclose(de_a["log2fc"].to_numpy(), de_b["log2fc"].to_numpy())


def _run_with_redaction(expr_path, meta_path, out_path, blinded_csv, record_path):
    return subprocess.run(
        ["Rscript", str(SCRIPT), str(expr_path), str(meta_path), str(out_path), blinded_csv, str(record_path)],
        capture_output=True,
        text=True,
    )


class TestBlindedGeneRedaction:
    """Guard tests for the preregistered-blind-control redaction path.

    PREANALYSIS.md S5 / CLAUDE.md: RCOR1 and KDM1A must be fit exactly like
    every other gene (so BH correction and every other gene's statistics
    are unaffected) but withheld from every generated result table until
    the model is frozen. See redact_blinded_genes() in
    gse118713_limma_lib.R.
    """

    def test_redacted_genes_absent_from_output(self, tmp_path):
        expr_path = tmp_path / "expr.tsv.gz"
        meta_path = tmp_path / "meta.tsv"
        out_path = tmp_path / "out.tsv.gz"
        record_path = tmp_path / "record.tsv"
        _write_expr(expr_path, _synthetic_expression(n_genes=10, seed=20))
        _write_meta(meta_path)
        result = _run_with_redaction(expr_path, meta_path, out_path, "G3,G7", record_path)
        assert result.returncode == 0, result.stderr
        de = _read_output(out_path)
        assert "G3" not in set(de["gene_id"])
        assert "G7" not in set(de["gene_id"])
        assert set(de["gene_id"]) == {f"G{i}" for i in range(10)} - {"G3", "G7"}

    def test_non_redacted_genes_unaffected_by_redaction(self, tmp_path):
        # The whole point of post-hoc redaction: fit + BH correction happen
        # on the complete gene set first, so every other gene's logFC, SE,
        # t, p-value, FDR, and AveExpr must be bit-for-bit identical to an
        # unredacted run.
        expr_path = tmp_path / "expr.tsv.gz"
        meta_path = tmp_path / "meta.tsv"
        out_unredacted = tmp_path / "out_unredacted.tsv.gz"
        out_redacted = tmp_path / "out_redacted.tsv.gz"
        record_path = tmp_path / "record.tsv"
        _write_expr(expr_path, _synthetic_expression(n_genes=15, seed=21))
        _write_meta(meta_path)
        _run(expr_path, meta_path, out_unredacted)
        _run_with_redaction(expr_path, meta_path, out_redacted, "G2,G9", record_path)

        unredacted = _read_output(out_unredacted)
        redacted = _read_output(out_redacted)
        merged = unredacted.merge(redacted, on=["gene_id", "contrast"], suffixes=("_full", "_redacted"))
        assert len(merged) == len(redacted)  # every reported gene matched a full-run row
        for col in ["log2fc", "se", "moderated_t", "p_value", "fdr", "ave_expr"]:
            assert (merged[f"{col}_full"] == merged[f"{col}_redacted"]).all(), col

    def test_redaction_record_counts(self, tmp_path):
        expr_path = tmp_path / "expr.tsv.gz"
        meta_path = tmp_path / "meta.tsv"
        out_path = tmp_path / "out.tsv.gz"
        record_path = tmp_path / "record.tsv"
        _write_expr(expr_path, _synthetic_expression(n_genes=12, seed=22))
        _write_meta(meta_path)
        _run_with_redaction(expr_path, meta_path, out_path, "G0,G11", record_path)
        record = pd.read_csv(record_path, sep="\t")
        assert record.loc[0, "genes_fitted"] == 12
        assert record.loc[0, "genes_withheld"] == 2
        assert record.loc[0, "genes_reported"] == 10

    def test_empty_blinded_list_redacts_nothing_but_still_writes_record(self, tmp_path):
        expr_path = tmp_path / "expr.tsv.gz"
        meta_path = tmp_path / "meta.tsv"
        out_path = tmp_path / "out.tsv.gz"
        record_path = tmp_path / "record.tsv"
        _write_expr(expr_path, _synthetic_expression(n_genes=8, seed=23))
        _write_meta(meta_path)
        result = _run_with_redaction(expr_path, meta_path, out_path, "", record_path)
        assert result.returncode == 0, result.stderr
        de = _read_output(out_path)
        assert set(de["gene_id"]) == {f"G{i}" for i in range(8)}
        record = pd.read_csv(record_path, sep="\t")
        assert record.loc[0, "genes_withheld"] == 0
        assert record.loc[0, "genes_fitted"] == record.loc[0, "genes_reported"] == 8

    def test_redaction_never_prints_gene_identities_to_stdout(self, tmp_path):
        expr_path = tmp_path / "expr.tsv.gz"
        meta_path = tmp_path / "meta.tsv"
        out_path = tmp_path / "out.tsv.gz"
        record_path = tmp_path / "record.tsv"
        _write_expr(expr_path, _synthetic_expression(n_genes=10, seed=24), gene_ids=[f"BLIND{i}" for i in range(10)])
        _write_meta(meta_path)
        result = _run_with_redaction(expr_path, meta_path, out_path, "BLIND2,BLIND5", record_path)
        assert result.returncode == 0, result.stderr
        assert "BLIND2" not in result.stdout
        assert "BLIND5" not in result.stdout

    def test_wrong_argument_count_rejected_clearly(self, tmp_path):
        # 4 args (missing either the blinded-list or the record path) must
        # not silently fall back to the 3-arg unredacted path.
        expr_path = tmp_path / "expr.tsv.gz"
        meta_path = tmp_path / "meta.tsv"
        out_path = tmp_path / "out.tsv.gz"
        _write_expr(expr_path, _synthetic_expression(n_genes=5))
        _write_meta(meta_path)
        result = subprocess.run(
            ["Rscript", str(SCRIPT), str(expr_path), str(meta_path), str(out_path), "G0"],
            capture_output=True,
            text=True,
        )
        assert result.returncode != 0
        assert "usage" in result.stderr.lower()

    def test_real_configured_blinded_genes_absent_from_committed_output_files(self):
        # If the real Phase 2B outputs exist on disk (i.e. the pipeline has
        # been run against the frozen matrix), the two preregistered blind
        # genes must never appear in them -- this is the actual guard
        # against the leakage the redaction path exists to prevent.
        import yaml

        with open(REPO_ROOT / "config" / "config.yaml") as f:
            config = yaml.safe_load(f)
        blinded_gene_ids = config["gse118713_phase2b"]["blinding"]["blinded_gene_ids"]
        assert len(blinded_gene_ids) == 2

        rel_paths = (
            config["gse118713_phase2b"]["limma"]["differential_expression_tsv_gz"],
            config["gse118713_phase2b"]["specificity"]["output_tsv_gz"],
        )
        missing = [p for p in rel_paths if not (REPO_ROOT / p).exists()]
        if missing:
            # Skip visibly rather than silently checking only whichever
            # subset happens to exist -- a partial checkout must not look
            # like a full pass of this guard.
            pytest.skip(
                f"real Phase 2B output file(s) not present in this checkout: {missing} -- "
                "run `python -m src.gse118713_phase2b` to exercise this guard"
            )
        for rel_path in rel_paths:
            df = pd.read_csv(REPO_ROOT / rel_path, sep="\t")
            for blind_id in blinded_gene_ids:
                assert blind_id not in set(df["gene_id"]), f"blinded gene present in {rel_path}"
