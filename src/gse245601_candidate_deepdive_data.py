"""GSE245601 candidate deep-dive: shared data-loading utilities. Joins
the R-extracted per-cell candidate expression
(`results/tables/gse245601_candidate_deepdive/candidate_per_cell_expression.tsv`,
built by `scripts/analysis/gse245601_18_extract_candidate_deepdive.R`)
with the already-frozen per-cell metadata table
(`results/tables/gse245601_pseudobulk/malignant_vs_nonmalignant_cell_level_summary.tsv`)
by `cell_id`. Every downstream table/figure in this deep-dive reads
through this module, so there is exactly one source of per-cell truth.

Cells are descriptive units only. No function in this module or any
consumer of it may treat individual cells as independent biological
replicates for an inferential (p-value-bearing) claim -- the tumor/
patient remains the unit for inference (see
`src.gse245601_candidate_deepdive_patient` for the inferential layer).
"""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

logger = logging.getLogger(__name__)

GENES = ["USP34", "VEZF1", "EML5", "CITED2"]


def _load_config(config_path: str | Path) -> dict:
    with open(config_path) as f:
        return yaml.safe_load(f)


def load_per_cell_table(config: dict) -> pd.DataFrame:
    cfg = config["gse245601_candidate_deepdive"]
    per_cell = pd.read_csv(cfg["output"]["per_cell_tsv"], sep="\t")
    metadata = pd.read_csv(cfg["inputs"]["cell_level_summary_tsv"], sep="\t")

    merged = metadata.merge(per_cell, on="cell_id", how="inner", validate="one_to_one")
    if len(merged) != len(metadata) or len(merged) != len(per_cell):
        raise ValueError(f"cell join lost or duplicated rows: metadata={len(metadata)}, per_cell={len(per_cell)}, merged={len(merged)}")
    logger.info("load_per_cell_table: %d epithelial cells, %d patients, malignancy=%s", len(merged), merged["patient"].nunique(), dict(merged["malignancy_status"].value_counts()))
    return merged


def load_track_a(config: dict) -> tuple[pd.DataFrame, pd.DataFrame]:
    """All-epithelial pseudobulk raw counts + metadata (Track A)."""
    cfg = config["gse245601_candidate_deepdive"]["inputs"]
    counts = pd.read_csv(cfg["track_a_counts_tsv"], sep="\t").set_index("gene")
    metadata = pd.read_csv(cfg["track_a_metadata_tsv"], sep="\t")
    return counts, metadata


def load_track_b(config: dict) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Strict-malignant pseudobulk raw counts + metadata (Track B)."""
    cfg = config["gse245601_candidate_deepdive"]["inputs"]
    counts = pd.read_csv(cfg["track_b_counts_tsv"], sep="\t").set_index("gene")
    metadata = pd.read_csv(cfg["track_b_metadata_tsv"], sep="\t")
    return counts, metadata


def compute_log2cpm(counts: pd.DataFrame, sample_ids: list[str]) -> pd.DataFrame:
    """Same convention as `src.gse245601_pseudobulk_qc.compute_log2cpm`:
    library size = each pseudobulk sample's own column sum, log2(CPM+1).
    Reimplemented locally (not imported) so this deep-dive module has no
    hidden coupling to a module that could change independently -- values
    are cross-checked against the frozen output in Phase 3."""
    lib_sizes = counts[sample_ids].sum(axis=0)
    cpm = counts[sample_ids].div(lib_sizes, axis=1) * 1e6
    return np.log2(cpm + 1)
