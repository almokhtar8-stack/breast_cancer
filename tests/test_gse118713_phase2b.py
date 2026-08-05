import hashlib
import shutil
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import yaml

from src.gse118713_phase2b import run_limma_script, run_phase2b

REPO_ROOT = Path(__file__).parent.parent
SAMPLE_IDS = [
    "MCF7_Rep1", "MCF7_Rep2", "MCF7_Rep3",
    "TAMR_Rep1", "TAMR_Rep2", "TAMR_Rep3",
    "FASR_Rep1", "FASR_Rep2", "FASR_Rep3",
]
GROUPS = ["MCF7"] * 3 + ["TAMR"] * 3 + ["FASR"] * 3

pytestmark = pytest.mark.skipif(shutil.which("Rscript") is None, reason="Rscript not available")


class TestRunLimmaScript:
    def test_reports_failure_clearly(self, tmp_path):
        bad_script = tmp_path / "always_fails.R"
        bad_script.write_text('message("ERROR: synthetic failure"); quit(status = 1, save = "no")\n')
        with pytest.raises(RuntimeError, match="synthetic failure"):
            run_limma_script(
                expression_tsv_gz=tmp_path / "unused.tsv.gz",
                metadata_tsv=tmp_path / "unused_meta.tsv",
                output_tsv_gz=tmp_path / "out.tsv.gz",
                script_path=bad_script,
                blinded_gene_ids=[],
                redaction_record_tsv=tmp_path / "record.tsv",
            )

    def test_succeeds_silently_on_valid_script(self, tmp_path):
        ok_script = tmp_path / "always_succeeds.R"
        ok_script.write_text('cat("ok\\n")\n')
        run_limma_script(
            expression_tsv_gz=tmp_path / "unused.tsv.gz",
            metadata_tsv=tmp_path / "unused_meta.tsv",
            output_tsv_gz=tmp_path / "out.tsv.gz",
            script_path=ok_script,
            blinded_gene_ids=[],
            redaction_record_tsv=tmp_path / "record.tsv",
        )  # must not raise

    def test_blinding_args_are_required(self, tmp_path):
        # No implicit "no blinding" default: every caller must make an
        # explicit choice, even if that choice is an empty list.
        ok_script = tmp_path / "always_succeeds.R"
        ok_script.write_text('cat("ok\\n")\n')
        with pytest.raises(TypeError):
            run_limma_script(
                expression_tsv_gz=tmp_path / "unused.tsv.gz",
                metadata_tsv=tmp_path / "unused_meta.tsv",
                output_tsv_gz=tmp_path / "out.tsv.gz",
                script_path=ok_script,
            )


class TestFullPipelineOnSyntheticFixture:
    def _build_config(self, tmp_path: Path, n_genes: int = 30, blinded_gene_ids: list[str] | None = None) -> Path:
        blinded_gene_ids = blinded_gene_ids or []
        rng = np.random.default_rng(42)
        base = rng.uniform(1, 200, n_genes)
        gene_ids = [f"ENSG{i:05d}" for i in range(n_genes)]

        data = {
            "gene_id": gene_ids,
            "gene_symbol": [f"SYM{i}" for i in range(n_genes)],
            "symbol_mapping_status": ["resolved"] * n_genes,
        }
        for sample_id, group in zip(SAMPLE_IDS, GROUPS):
            shift = {"MCF7": 0, "TAMR": 4, "FASR": -2}[group]
            data[sample_id] = np.clip(base + shift + rng.normal(0, 0.5, n_genes), 0.01, None)
        # A handful of genes below the detection floor -- they must be
        # filtered out and not reach the DE step.
        for sample_id in SAMPLE_IDS:
            data[sample_id][-3:] = 0.0

        gene_df = pd.DataFrame(data)
        parquet_path = tmp_path / "gene_tpm.parquet"
        gene_df.to_parquet(parquet_path, index=False)
        checksum = hashlib.sha256(parquet_path.read_bytes()).hexdigest()

        meta_df = pd.DataFrame({"sample_id": SAMPLE_IDS, "column": SAMPLE_IDS, "group": GROUPS, "replicate": [1, 2, 3] * 3})
        meta_path = tmp_path / "sample_metadata.tsv"
        meta_df.to_csv(meta_path, sep="\t", index=False)

        config = {
            "gse118713": {
                "output": {
                    "gene_tpm_parquet": str(parquet_path),
                    "sample_metadata_tsv": str(meta_path),
                    "qc_tsv": str(tmp_path / "prep_qc.tsv"),
                },
                "frozen_gene_tpm_sha256": checksum,
            },
            "gse118713_phase2b": {
                "expected_n_genes": n_genes,
                "expected_n_samples": 9,
                "expected_groups": ["MCF7", "TAMR", "FASR"],
                "expected_replicates_per_group": 3,
                "filtering": {
                    "min_tpm": 1.0,
                    "min_samples": 3,
                    "filtered_gene_tpm_tsv": str(tmp_path / "filtered.tsv.gz"),
                    "filtering_summary_tsv": str(tmp_path / "filtering_summary.tsv"),
                },
                "qc": {
                    "n_pca_genes": 10,
                    "sample_qc_tsv": str(tmp_path / "sample_qc.tsv"),
                    "sample_correlations_tsv": str(tmp_path / "sample_correlations.tsv"),
                    "pca_coordinates_tsv": str(tmp_path / "pca_coordinates.tsv"),
                    "pca_figure_pdf": str(tmp_path / "pca.pdf"),
                    "correlation_figure_pdf": str(tmp_path / "correlation.pdf"),
                },
                "limma": {
                    "script": str(REPO_ROOT / "scripts" / "analysis" / "gse118713_limma.R"),
                    "contrasts": ["TAMR_vs_MCF7", "FASR_vs_MCF7", "TAMR_vs_FASR"],
                    "differential_expression_tsv_gz": str(tmp_path / "de.tsv.gz"),
                },
                "specificity": {"output_tsv_gz": str(tmp_path / "specificity.tsv.gz")},
                "summary": {"output_tsv": str(tmp_path / "de_summary.tsv")},
                "blinding": {
                    "blinded_gene_ids": blinded_gene_ids,
                    "redaction_record_tsv": str(tmp_path / "blinding_redaction.tsv"),
                },
            },
        }
        config_path = tmp_path / "config.yaml"
        config_path.write_text(yaml.dump(config))
        return config_path

    def test_end_to_end_writes_every_output(self, tmp_path):
        # Two of the 27 filtered synthetic genes stand in for the real
        # preregistered blind controls, to prove the end-to-end wiring
        # (config -> orchestrator -> R script) actually redacts them from
        # every committed table, not just that the R script can do it in
        # isolation.
        blinded_gene_ids = ["ENSG00000", "ENSG00005"]
        config_path = self._build_config(tmp_path, blinded_gene_ids=blinded_gene_ids)
        outputs = run_phase2b(config_path)

        assert (tmp_path / "filtering_summary.tsv").exists()
        assert (tmp_path / "sample_qc.tsv").exists()
        assert (tmp_path / "sample_correlations.tsv").exists()
        assert (tmp_path / "pca_coordinates.tsv").exists()
        assert (tmp_path / "pca.pdf").exists()
        assert (tmp_path / "correlation.pdf").exists()
        assert (tmp_path / "de.tsv.gz").exists()
        assert (tmp_path / "specificity.tsv.gz").exists()
        assert (tmp_path / "de_summary.tsv").exists()
        assert (tmp_path / "blinding_redaction.tsv").exists()

        filtering_summary = pd.read_csv(tmp_path / "filtering_summary.tsv", sep="\t")
        assert filtering_summary.loc[0, "genes_in"] == 30
        assert filtering_summary.loc[0, "genes_removed"] == 3  # the three zeroed-out genes

        redaction = pd.read_csv(tmp_path / "blinding_redaction.tsv", sep="\t")
        assert redaction.loc[0, "genes_fitted"] == 27
        assert redaction.loc[0, "genes_withheld"] == 2
        assert redaction.loc[0, "genes_reported"] == 25

        de = pd.read_csv(tmp_path / "de.tsv.gz", sep="\t")
        specificity = pd.read_csv(tmp_path / "specificity.tsv.gz", sep="\t")
        for blind_id in blinded_gene_ids:
            assert blind_id not in set(de["gene_id"])
            assert blind_id not in set(specificity["gene_id"])

        de_summary = outputs["de_summary"]
        assert set(de_summary["contrast"]) == {"TAMR_vs_MCF7", "FASR_vs_MCF7", "TAMR_vs_FASR"}
        assert (de_summary["genes_tested"] == 25).all()

    def test_checksum_mismatch_stops_pipeline_before_any_output(self, tmp_path):
        config_path = self._build_config(tmp_path)
        config = yaml.safe_load(config_path.read_text())
        config["gse118713"]["frozen_gene_tpm_sha256"] = "0" * 64
        config_path.write_text(yaml.dump(config))

        with pytest.raises(ValueError, match="checksum mismatch"):
            run_phase2b(config_path)
        assert not (tmp_path / "de.tsv.gz").exists()

    def test_blinded_gene_ids_never_appear_in_python_logs(self, tmp_path, caplog):
        # The R subprocess necessarily receives the blinded gene IDs as an
        # argv entry, but the Python orchestrator must never write that
        # entry into application logs.
        blinded_gene_ids = ["ENSG00000", "ENSG00005"]  # matches the synthetic fixture's gene_id pattern
        config_path = self._build_config(tmp_path, blinded_gene_ids=blinded_gene_ids)
        with caplog.at_level("DEBUG"):
            run_phase2b(config_path)
        log_text = "\n".join(record.getMessage() for record in caplog.records)
        for blind_id in blinded_gene_ids:
            assert blind_id not in log_text

    def test_unmatched_blinded_gene_id_fails_loudly_not_silently(self, tmp_path):
        # A typo'd or already-filtered-out configured ID must never result
        # in a quietly smaller (or absent) redaction.
        config_path = self._build_config(tmp_path, blinded_gene_ids=["ENSG00000", "NOT_A_REAL_GENE_ID"])
        with pytest.raises(RuntimeError, match="count mismatch"):
            run_phase2b(config_path)

    def test_duplicate_blinded_gene_id_fails_loudly(self, tmp_path):
        config_path = self._build_config(tmp_path, blinded_gene_ids=["ENSG00000", "ENSG00000"])
        with pytest.raises(RuntimeError, match="duplicate"):
            run_phase2b(config_path)
