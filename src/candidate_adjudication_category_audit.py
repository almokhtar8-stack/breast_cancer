"""Candidate adjudication Phase 1: independently reconstruct the
MULTIMODAL_STRONG evidence-category rule from
``src.cross_dataset_consensus_views.assign_evidence_category`` and verify,
programmatically (never by hard-coding the seven names), which genes
satisfy it. Also surfaces genes that came close but did not qualify, so
the adjudication phase does not start from an assumed-correct answer.

Data source: `results/tables/cross_dataset_genomewide/all_genes_cross_dataset_evidence_with_ranking.tsv`
and `resistance_consensus_all_genes.tsv` (both already frozen, committed
by the prior cross-dataset genome-wide integration phase; read-only here).
"""

from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd
import yaml

logger = logging.getLogger(__name__)

CRISPR_FDR_SIGNIFICANT = 0.10  # matches src.cross_dataset_consensus_views.CRISPR_FDR_SIGNIFICANT


def _load_config(config_path: str | Path) -> dict:
    with open(config_path) as f:
        return yaml.safe_load(f)


RESISTANCE_RAW_FDR_COLS = ["gse118713_fdr", "gse240112_tumor_fdr", "gse111151_fdr"]
ALL_FIVE_RAW_FDR_COLS = ["crispr_fdr", "gse118713_fdr", "gse245601_epi_fdr", "gse240112_tumor_fdr", "gse111151_fdr"]
ALL_FIVE_TESTABLE_COLS = ["crispr_testable", "gse118713_testable", "gse245601_testable", "gse240112_testable", "gse111151_testable"]


def reconstruct_multimodal_strong(wide: pd.DataFrame, resistance_consensus: pd.DataFrame) -> pd.DataFrame:
    """Independent re-derivation of the MULTIMODAL_STRONG condition,
    written from scratch against the rule's own docstring rather than by
    importing ``assign_evidence_category``. Deliberately does NOT reuse
    the frozen wide table's own precomputed ``n_datasets_fdr05`` or the
    resistance table's own precomputed ``resistance_fdr05_count`` -- both
    are recomputed here directly from the five raw per-dataset FDR
    columns (the three resistance-dataset FDR columns are a subset of
    those five, used a second time for the resistance-specific count; and,
    for the LOW_COVERAGE gate, the five raw ``*_testable`` columns), so
    that a bug in the original aggregation
    step could not silently propagate into and be hidden by this
    'verification'. Returns one row per gene with the boolean
    sub-conditions exposed (not just the final label), so it is possible
    to see *why* a gene did or did not qualify."""
    merged = wide.merge(resistance_consensus[["gene", "resistance_direction_consensus"]], on="gene", how="left")

    n_testable_reconstructed = merged[ALL_FIVE_TESTABLE_COLS].sum(axis=1)
    low_coverage = n_testable_reconstructed < 3

    crispr_strong = (merged["crispr_fdr"] < CRISPR_FDR_SIGNIFICANT).fillna(False)
    resistance_fdr05_count_reconstructed = sum((merged[c] < 0.05).fillna(False).astype(int) for c in RESISTANCE_RAW_FDR_COLS)
    n_fdr05_overall_reconstructed = sum((merged[c] < 0.05).fillna(False).astype(int) for c in ALL_FIVE_RAW_FDR_COLS)

    cond_a = resistance_fdr05_count_reconstructed >= 2
    cond_b = (resistance_fdr05_count_reconstructed >= 1) & (n_fdr05_overall_reconstructed >= 2)
    qualifies = (~low_coverage) & crispr_strong & (cond_a | cond_b)

    out = merged[["gene", "coverage_tier", "crispr_fdr", "resistance_direction_consensus"]].copy()
    out["n_datasets_testable_reconstructed"] = n_testable_reconstructed
    out["low_coverage_reconstructed"] = low_coverage
    out["n_datasets_fdr05_reconstructed"] = n_fdr05_overall_reconstructed
    out["resistance_fdr05_count_reconstructed"] = resistance_fdr05_count_reconstructed
    out["crispr_strong_fdr_lt_010"] = crispr_strong
    out["resistance_fdr05_ge2"] = cond_a
    out["resistance_fdr05_ge1_and_overall_fdr05_ge2"] = cond_b
    out["is_multimodal_strong"] = qualifies
    logger.info("reconstruct_multimodal_strong: %d genes qualify (raw-column recomputation, not reusing the frozen aggregate columns)", int(qualifies.sum()))
    return out.loc[qualifies].sort_values("gene").reset_index(drop=True)


def find_near_misses(wide: pd.DataFrame, resistance_consensus: pd.DataFrame, multimodal_genes: set[str]) -> pd.DataFrame:
    """Genes that satisfy exactly one of the two required conditions
    (CRISPR FDR<0.10, OR the resistance-FDR pattern) but not both, ranked
    by how close they are to qualifying. 'Distance' is reported as the
    specific missing condition, never a synthesized numeric score."""
    merged = wide.merge(resistance_consensus[["gene", "resistance_fdr05_count", "resistance_direction_consensus"]], on="gene", how="left")
    merged = merged.loc[~merged["gene"].isin(multimodal_genes)].copy()

    crispr_strong = (merged["crispr_fdr"] < CRISPR_FDR_SIGNIFICANT).fillna(False)
    resistance_fdr05_count = merged["resistance_fdr05_count"].fillna(0)
    n_fdr05_overall = merged["n_datasets_fdr05"].fillna(0)
    resistance_pattern_met = (resistance_fdr05_count >= 2) | ((resistance_fdr05_count >= 1) & (n_fdr05_overall >= 2))

    has_crispr_only = crispr_strong & ~resistance_pattern_met & (resistance_fdr05_count >= 1)
    has_resistance_only = (~crispr_strong) & resistance_pattern_met & (merged["crispr_fdr"] < 0.20).fillna(False)

    near = merged.loc[has_crispr_only | has_resistance_only].copy()
    near["missing_condition"] = "unknown"
    near.loc[has_crispr_only.reindex(near.index, fill_value=False), "missing_condition"] = "has_CRISPR_FDR<0.10_but_only_1_resistance_dataset_FDR<0.05_and_<2_overall"
    near.loc[has_resistance_only.reindex(near.index, fill_value=False), "missing_condition"] = "has_resistance_FDR_pattern_but_CRISPR_FDR_in_[0.10,0.20)"
    out = near[["gene", "coverage_tier", "crispr_fdr", "n_datasets_fdr05", "resistance_fdr05_count", "resistance_direction_consensus", "missing_condition"]].copy()
    out = out.sort_values(by=["missing_condition", "crispr_fdr", "gene"]).reset_index(drop=True)
    logger.info("find_near_misses: %d genes near but not in MULTIMODAL_STRONG", len(out))
    return out


def run_category_audit(config_path: str | Path = "config/config.yaml") -> dict[str, pd.DataFrame]:
    config = _load_config(config_path)
    cdx_out = config["cross_dataset_genomewide"]["output"]
    tables_dir = Path(cdx_out["wide_matrix_tsv"]).parent

    wide = pd.read_csv(tables_dir / "all_genes_cross_dataset_evidence_with_ranking.tsv", sep="\t")
    resistance = pd.read_csv(cdx_out["resistance_consensus_tsv"], sep="\t")
    categories = pd.read_csv(tables_dir / "evidence_categories.tsv", sep="\t")

    reconstructed = reconstruct_multimodal_strong(wide, resistance)
    official = set(categories.loc[categories["evidence_category"] == "MULTIMODAL_STRONG", "gene"])
    mismatch = official.symmetric_difference(set(reconstructed["gene"]))
    if mismatch:
        raise ValueError(f"reconstructed MULTIMODAL_STRONG set does not match the frozen evidence_categories.tsv output: {sorted(mismatch)}")

    near_misses = find_near_misses(wide, resistance, official)

    adj_cfg = config["candidate_adjudication"]["output"]
    out_dir = Path(adj_cfg["tables_dir"])
    out_dir.mkdir(parents=True, exist_ok=True)
    reconstructed.to_csv(out_dir / "multimodal_strong_reconstructed.tsv", sep="\t", index=False)
    near_misses.to_csv(out_dir / "multimodal_strong_near_misses_raw.tsv", sep="\t", index=False)

    logger.info("run_category_audit: %d MULTIMODAL_STRONG genes confirmed to match the frozen category assignment exactly", len(reconstructed))
    return {"multimodal_strong": reconstructed, "near_misses": near_misses}


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run_category_audit()
