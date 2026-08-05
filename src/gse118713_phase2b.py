"""Phase 2B orchestrator: QC, filtering, limma DE, specificity and summary
for the frozen GSE118713 gene-level TPM matrix.

Runs, in order: frozen-matrix validation and expression filtering
(``src.gse118713_expression_filter``), sample QC (``src.gse118713_qc``),
limma differential expression via ``scripts/analysis/gse118713_limma.R``
(invoked as a subprocess -- limma/statmod have no Python binding used
here), the TAMR-specificity join (``src.gse118713_tamr_specificity``), and
the per-contrast summary (``src.gse118713_de_summary``).

See PREANALYSIS.md's 2026-08-05 Phase 2B statistical-plan amendment for
the preregistered rules this pipeline implements.
"""

from __future__ import annotations

import logging
import subprocess
import sys
from pathlib import Path

import yaml

from src.gse118713_de_summary import run_de_summary
from src.gse118713_expression_filter import FilterConfig, load_frozen_matrix, load_sample_metadata, run_expression_filtering
from src.gse118713_qc import run_sample_qc
from src.gse118713_tamr_specificity import run_tamr_specificity

logger = logging.getLogger(__name__)


def _load_config(config_path: str | Path) -> dict:
    with open(config_path) as f:
        return yaml.safe_load(f)


def run_limma_script(
    expression_tsv_gz: str | Path,
    metadata_tsv: str | Path,
    output_tsv_gz: str | Path,
    script_path: str | Path,
    blinded_gene_ids: list[str],
    redaction_record_tsv: str | Path,
) -> None:
    """Invoke the limma R script and fail loudly (with stderr) on any error.

    ``blinded_gene_ids`` and ``redaction_record_tsv`` are required (not
    optional) so no caller can produce a reportable result without making
    an explicit blinding decision -- pass an empty list to explicitly
    redact nothing. Those genes are fit exactly like every other gene (BH
    correction runs on the complete gene set) and are withheld from
    ``output_tsv_gz`` only at the reporting stage -- see
    ``redact_blinded_genes()`` in ``gse118713_limma_lib.R``.

    Only counts are ever logged here. The blinded gene IDs are passed to
    the subprocess via argv but deliberately excluded from the logged
    command line so they never reach application logs.
    """
    Path(output_tsv_gz).parent.mkdir(parents=True, exist_ok=True)
    Path(redaction_record_tsv).parent.mkdir(parents=True, exist_ok=True)
    blinded_gene_ids_csv = ",".join(blinded_gene_ids)
    cmd = [
        "Rscript",
        str(script_path),
        str(expression_tsv_gz),
        str(metadata_tsv),
        str(output_tsv_gz),
        blinded_gene_ids_csv,
        str(redaction_record_tsv),
    ]
    loggable_cmd = list(cmd)
    loggable_cmd[5] = f"<redacted: {len(blinded_gene_ids)} gene id(s)>"
    logger.info("run_limma_script: %s", " ".join(loggable_cmd))
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.stdout:
        logger.info("run_limma_script stdout: %s", result.stdout.strip())
    if result.returncode != 0:
        logger.error("run_limma_script stderr: %s", result.stderr.strip())
        raise RuntimeError(
            f"gse118713_limma.R failed with exit code {result.returncode}:\n{result.stderr}"
        )
    logger.info("run_limma_script: completed successfully")


def run_phase2b(config_path: str | Path = "config/config.yaml") -> dict[str, object]:
    config = _load_config(config_path)
    filter_cfg = FilterConfig.from_config(config)
    phase2b = config["gse118713_phase2b"]

    meta = load_sample_metadata(filter_cfg)
    sample_ids = list(meta["sample_id"])
    df = load_frozen_matrix(filter_cfg, meta)

    filtered_df, filtering_summary = run_expression_filtering(config_path)

    qc_outputs = run_sample_qc(df, filtered_df, meta, sample_ids, config_path)

    blinding = phase2b["blinding"]
    run_limma_script(
        expression_tsv_gz=filter_cfg.filtered_gene_tpm_tsv,
        metadata_tsv=filter_cfg.sample_metadata_tsv,
        output_tsv_gz=phase2b["limma"]["differential_expression_tsv_gz"],
        script_path=phase2b["limma"]["script"],
        blinded_gene_ids=list(blinding["blinded_gene_ids"]),
        redaction_record_tsv=blinding["redaction_record_tsv"],
    )

    specificity_df = run_tamr_specificity(config_path)
    summary_df = run_de_summary(config_path)

    return {
        "filtering_summary": filtering_summary,
        **qc_outputs,
        "specificity": specificity_df,
        "de_summary": summary_df,
    }


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    try:
        run_phase2b()
    except Exception as exc:  # noqa: BLE001 -- top-level CLI failure reporting
        logger.error("run_phase2b failed: %s", exc)
        sys.exit(1)
