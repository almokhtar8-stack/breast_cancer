"""Controlled unredaction of the GSE118713 differential-expression blind.

Source: PREANALYSIS.md's 2026-08-10 provenance-correction amendment
authorises releasing RCOR1/KDM1A's withheld GSE118713 rows using the exact
frozen analysis (``scripts/analysis/gse118713_limma.R`` /
``gse118713_limma_lib.R``), reading the same frozen inputs already used
for the redacted run (config ``gse118713_phase2b.filtering.
filtered_gene_tpm_tsv``, ``gse118713.output.sample_metadata_tsv``), with
only the export-stage redaction step disabled (``blinded_gene_ids`` passed
as an empty list). No modelling, filtering, or statistical choice is
changed from the frozen Phase 2B specification; only the export step
differs.

This module never overwrites the original redacted outputs
(``results/tables/gse118713_differential_expression.tsv.gz``,
``results/tables/gse118713_tamr_specificity.tsv.gz``) -- it writes a
separate, versioned "_unredacted" derivative (config
``gse118713_phase2b.unredaction``) and produces a row-for-row provenance
comparison between the two, confirming every previously-exported gene's
statistics are unchanged and exactly the configured blinded gene IDs were
added.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from src.gse118713_phase2b import run_limma_script
from src.gse118713_tamr_specificity import (
    REQUIRED_CONTRASTS,
    SpecificityConfig,
    build_specificity_table,
    load_de_table,
    write_specificity_table,
)

logger = logging.getLogger(__name__)

COMPARISON_FIELDS: tuple[str, ...] = ("log2fc", "se", "moderated_t", "p_value", "fdr", "ave_expr")

# Two independent Rscript invocations of the identical fit can differ by
# BLAS/LAPACK floating-point non-associativity (observed magnitude ~1e-13
# here), not by anything redaction-related -- the fit and BH correction are
# untouched by blinded_gene_ids. NUMERIC_ATOL is set two orders of
# magnitude above that observed noise floor and many orders below any
# scientifically meaningful difference, so it distinguishes reproducibility
# jitter from an actual discrepancy rather than masking one.
NUMERIC_ATOL = 1e-8


def _load_config(config_path: str | Path) -> dict:
    with open(config_path) as f:
        return yaml.safe_load(f)


@dataclass(frozen=True)
class UnredactionConfig:
    """Resolved, config-driven paths. No hardcoded paths."""

    filtered_gene_tpm_tsv: Path
    sample_metadata_tsv: Path
    limma_script: Path
    old_de_tsv_gz: Path
    new_de_tsv_gz: Path
    new_redaction_record_tsv: Path
    new_specificity_tsv_gz: Path
    comparison_tsv: Path
    blinded_gene_ids: tuple[str, ...]

    @classmethod
    def from_config(cls, config: dict) -> "UnredactionConfig":
        gse = config["gse118713"]
        phase2b = config["gse118713_phase2b"]
        unred = phase2b["unredaction"]
        return cls(
            filtered_gene_tpm_tsv=Path(phase2b["filtering"]["filtered_gene_tpm_tsv"]),
            sample_metadata_tsv=Path(gse["output"]["sample_metadata_tsv"]),
            limma_script=Path(phase2b["limma"]["script"]),
            old_de_tsv_gz=Path(phase2b["limma"]["differential_expression_tsv_gz"]),
            new_de_tsv_gz=Path(unred["differential_expression_unredacted_tsv_gz"]),
            new_redaction_record_tsv=Path(unred["redaction_record_unredacted_tsv"]),
            new_specificity_tsv_gz=Path(unred["specificity_unredacted_tsv_gz"]),
            comparison_tsv=Path(unred["comparison_tsv"]),
            blinded_gene_ids=tuple(phase2b["blinding"]["blinded_gene_ids"]),
        )


def run_unredacted_limma(cfg: UnredactionConfig) -> None:
    """Rerun the frozen limma script on the exact same inputs, redaction disabled.

    Same filtered TPM matrix, same sample metadata, same script -- only
    ``blinded_gene_ids`` differs (empty instead of the configured RCOR1/
    KDM1A pair), which disables ``redact_blinded_genes()`` entirely
    (``gse118713_limma_lib.R``'s zero-length-list branch returns the
    complete fitted table unchanged). No fitting, contrast, or BH-universe
    choice is touched.
    """
    run_limma_script(
        expression_tsv_gz=cfg.filtered_gene_tpm_tsv,
        metadata_tsv=cfg.sample_metadata_tsv,
        output_tsv_gz=cfg.new_de_tsv_gz,
        script_path=cfg.limma_script,
        blinded_gene_ids=[],
        redaction_record_tsv=cfg.new_redaction_record_tsv,
    )
    logger.info("run_unredacted_limma: wrote unredacted DE table to %s", cfg.new_de_tsv_gz)


def build_unredacted_specificity(cfg: UnredactionConfig) -> pd.DataFrame:
    de_df = load_de_table(cfg.new_de_tsv_gz)
    specificity_df = build_specificity_table(de_df)
    write_specificity_table(
        specificity_df,
        SpecificityConfig(
            differential_expression_tsv_gz=cfg.new_de_tsv_gz,
            output_tsv_gz=cfg.new_specificity_tsv_gz,
        ),
    )
    logger.info("build_unredacted_specificity: wrote %d genes to %s", len(specificity_df), cfg.new_specificity_tsv_gz)
    return specificity_df


def compare_redacted_vs_unredacted(cfg: UnredactionConfig) -> pd.DataFrame:
    """Row-for-row provenance comparison, per contrast.

    Raises if any previously-reported gene is missing from the unredacted
    table, or if the unredacted table's added gene IDs are anything other
    than exactly the configured ``blinded_gene_ids``. Reports, per
    contrast, the shared-row max absolute difference for every numeric DE
    field and the count of shared-row FDR changes -- all expected to be
    zero, since BH correction was already finalized on the complete gene
    set before either export.
    """
    old_df = load_de_table(cfg.old_de_tsv_gz)
    new_df = load_de_table(cfg.new_de_tsv_gz)
    expected_added = set(cfg.blinded_gene_ids)

    records: list[dict[str, object]] = []
    for contrast in REQUIRED_CONTRASTS:
        old_c = old_df.loc[old_df["contrast"] == contrast].set_index("gene_id")
        new_c = new_df.loc[new_df["contrast"] == contrast].set_index("gene_id")

        old_ids = set(old_c.index)
        new_ids = set(new_c.index)
        added_ids = new_ids - old_ids
        removed_ids = old_ids - new_ids
        shared_ids = sorted(old_ids & new_ids)

        if removed_ids:
            raise ValueError(
                f"{contrast}: {len(removed_ids)} previously-reported gene(s) missing from "
                f"unredacted table: {sorted(removed_ids)}"
            )
        if added_ids != expected_added:
            raise ValueError(
                f"{contrast}: unredacted table added {sorted(added_ids)}, expected exactly "
                f"the configured blinded gene ids {sorted(expected_added)}"
            )

        symbol_mismatches = int((old_c.loc[shared_ids, "gene_symbol"].to_numpy() != new_c.loc[shared_ids, "gene_symbol"].to_numpy()).sum())
        if symbol_mismatches:
            raise ValueError(f"{contrast}: {symbol_mismatches} shared gene(s) have a different gene_symbol between old and new tables")

        shared_old = old_c.loc[shared_ids, list(COMPARISON_FIELDS)].to_numpy(dtype=float)
        shared_new = new_c.loc[shared_ids, list(COMPARISON_FIELDS)].to_numpy(dtype=float)
        abs_diff = np.abs(shared_old - shared_new)
        max_abs_diff = {field: float(np.max(abs_diff[:, i])) for i, field in enumerate(COMPARISON_FIELDS)}
        exceeds_tolerance = {field: v for field, v in max_abs_diff.items() if v > NUMERIC_ATOL}
        if exceeds_tolerance:
            raise ValueError(
                f"{contrast}: shared-gene values differ by more than the {NUMERIC_ATOL} "
                f"reproducibility tolerance: {exceeds_tolerance}"
            )

        old_fdr = old_c.loc[shared_ids, "fdr"].to_numpy(dtype=float)
        new_fdr = new_c.loc[shared_ids, "fdr"].to_numpy(dtype=float)
        n_fdr_changed = int(np.sum(np.abs(old_fdr - new_fdr) > NUMERIC_ATOL))

        record: dict[str, object] = {
            "contrast": contrast,
            "n_old_genes": len(old_ids),
            "n_new_genes": len(new_ids),
            "n_shared_genes": len(shared_ids),
            "n_added_genes": len(added_ids),
            "added_gene_ids": ",".join(sorted(added_ids)),
            "n_fdr_changed_for_shared_genes": n_fdr_changed,
        }
        for field, value in max_abs_diff.items():
            record[f"max_abs_diff_{field}"] = value
        records.append(record)
        logger.info("compare_redacted_vs_unredacted[%s]: %s", contrast, record)

    comparison_df = pd.DataFrame(records)
    cfg.comparison_tsv.parent.mkdir(parents=True, exist_ok=True)
    comparison_df.to_csv(cfg.comparison_tsv, sep="\t", index=False)
    logger.info("compare_redacted_vs_unredacted: wrote comparison to %s", cfg.comparison_tsv)
    return comparison_df


def run_unredaction(config_path: str | Path = "config/config.yaml") -> pd.DataFrame:
    """Run the full controlled-unredaction pipeline and write all outputs."""
    config = _load_config(config_path)
    cfg = UnredactionConfig.from_config(config)

    run_unredacted_limma(cfg)
    build_unredacted_specificity(cfg)
    return compare_redacted_vs_unredacted(cfg)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run_unredaction()
