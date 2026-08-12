"""Cross-dataset genome-wide integration, Phases 6-9: within-dataset
evidence percentiles, multi-track collapsing to one vote per independent
dataset, coverage tiers, and the transparent global equal-dataset
ranking.

No arbitrary weighted score is computed anywhere in this module (Phase
25). Every summary statistic is either a percentile (rank-derived,
0=weakest to 1=strongest within that one dataset) or a plain count/
median/mean of those percentiles, each explicitly labeled by name.

Percentile method (Phase 6, fixed and deterministic): within a dataset's
testable genes only, sort by FDR ascending, then nominal p ascending,
then |effect| descending; rank 1 = strongest. Percentile = 1 - (rank-1)/(N-1)
for N>1 testable genes (rank 1 -> 1.0, rank N -> 0.0); a single testable
gene gets percentile 1.0 (there is no weaker gene to rank below).

Track collapsing (Phase 7, fixed rule -- not a choice made per gene):
GSE245601's dataset-level percentile is the arithmetic mean of Track A
and Track B's own within-track percentiles (whichever are available);
GSE240112's dataset-level percentile is the tumor-cell track's percentile
alone (the all-epithelial track is sensitivity-only, never blended into
the vote). GSE118713 uses only the TAMR_vs_MCF7 contrast (already
enforced upstream in `cross_dataset_gene_mapping.load_gse118713`).

Data source: `results/tables/cross_dataset_genomewide/all_genes_cross_dataset_evidence.tsv`
(built by `src.cross_dataset_evidence_tables`).
"""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

logger = logging.getLogger(__name__)

COVERAGE_TIER_LABELS = {5: "A", 4: "B", 3: "C", 2: "D", 1: "E", 0: "none"}
DATASET_NAMES = ["crispr", "gse118713", "gse245601", "gse240112", "gse111151"]


def _load_config(config_path: str | Path) -> dict:
    with open(config_path) as f:
        return yaml.safe_load(f)


def compute_within_dataset_percentile(effect: pd.Series, p_value: pd.Series, fdr: pd.Series, gene: pd.Series | None = None) -> pd.Series:
    """Deterministic rank -> percentile for one dataset/track. Untestable
    genes (all three inputs NaN) get NaN, never 0 -- 0 would falsely
    imply "tested and weakest" rather than "not tested here". Genes tied
    on the full (fdr, p, |effect|) key are broken by gene symbol ascending
    if ``gene`` is supplied -- an explicit tie-break rather than relying
    on the caller's row order (which happens to already be alphabetical
    upstream, but that was an implicit invariant, not a guarantee this
    function itself enforced; made explicit after the Phase 27 Codex
    review flagged it as a robustness gap)."""
    testable = fdr.notna() | p_value.notna()
    n = int(testable.sum())
    out = pd.Series(np.nan, index=effect.index, dtype=float)
    if n == 0:
        return out
    if n == 1:
        out.loc[testable] = 1.0
        return out

    sub = pd.DataFrame({"fdr": fdr[testable], "p": p_value[testable], "abs_effect": effect[testable].abs()})
    sub["_neg_abs_effect"] = -sub["abs_effect"]
    by = ["fdr", "p", "_neg_abs_effect"]
    ascending = [True, True, True]
    if gene is not None:
        sub["_gene"] = gene[testable]
        by.append("_gene")
        ascending.append(True)
    ordered = sub.sort_values(by=by, ascending=ascending, na_position="last", kind="mergesort")
    rank = pd.Series(np.arange(1, n + 1), index=ordered.index)
    percentile = 1.0 - (rank - 1) / (n - 1)
    out.loc[percentile.index] = percentile.values
    return out


def compute_dataset_percentiles(wide: pd.DataFrame) -> pd.DataFrame:
    out = wide.copy()

    gene = out["gene"]
    out["crispr_evidence_percentile"] = compute_within_dataset_percentile(out["crispr_effect"], out["crispr_p"], out["crispr_fdr"], gene)
    out["gse118713_evidence_percentile"] = compute_within_dataset_percentile(out["gse118713_log2fc"], out["gse118713_p"], out["gse118713_fdr"], gene)

    out["gse245601_track_a_percentile"] = compute_within_dataset_percentile(out["gse245601_epi_log2fc"], out["gse245601_epi_p"], out["gse245601_epi_fdr"], gene)
    out["gse245601_track_b_percentile"] = compute_within_dataset_percentile(out["gse245601_malignant_log2fc"], out["gse245601_malignant_p"], out["gse245601_malignant_fdr"], gene)
    both = out[["gse245601_track_a_percentile", "gse245601_track_b_percentile"]]
    out["gse245601_evidence_percentile"] = both.mean(axis=1, skipna=True)
    out.loc[both.isna().all(axis=1), "gse245601_evidence_percentile"] = np.nan
    out["gse245601_one_track_only"] = both.notna().sum(axis=1) == 1

    out["gse240112_tumor_percentile"] = compute_within_dataset_percentile(out["gse240112_tumor_log2fc"], out["gse240112_tumor_p"], out["gse240112_tumor_fdr"], gene)
    out["gse240112_epi_percentile"] = compute_within_dataset_percentile(out["gse240112_epi_log2fc"], out["gse240112_epi_p"], out["gse240112_epi_fdr"], gene)
    out["gse240112_evidence_percentile"] = out["gse240112_tumor_percentile"]  # tumor-cell only -- Phase 7 explicit rule, epithelial never blended in

    out["gse111151_evidence_percentile"] = compute_within_dataset_percentile(out["gse111151_log2fc"], out["gse111151_p"], out["gse111151_fdr"], gene)

    logger.info("compute_dataset_percentiles: computed for %d genes", len(out))
    return out


def assign_coverage_tier(wide: pd.DataFrame) -> pd.DataFrame:
    out = wide.copy()
    testable_cols = [f"{d}_testable" for d in DATASET_NAMES]
    out["n_datasets_testable"] = out[testable_cols].sum(axis=1).astype(int)
    out["coverage_tier"] = out["n_datasets_testable"].map(COVERAGE_TIER_LABELS)
    logger.info("assign_coverage_tier: %s", out["coverage_tier"].value_counts().to_dict())
    return out


# representative percentile/fdr column used as "the dataset's one vote" for
# significance counting and top-10%/20% counting -- same track selection as
# compute_dataset_percentiles (Track A for GSE245601, tumor-cell for GSE240112)
_DATASET_PERCENTILE_COL = {
    "crispr": "crispr_evidence_percentile", "gse118713": "gse118713_evidence_percentile",
    "gse245601": "gse245601_evidence_percentile", "gse240112": "gse240112_evidence_percentile",
    "gse111151": "gse111151_evidence_percentile",
}
_DATASET_FDR_COL = {
    "crispr": "crispr_fdr", "gse118713": "gse118713_fdr", "gse245601": "gse245601_epi_fdr",
    "gse240112": "gse240112_tumor_fdr", "gse111151": "gse111151_fdr",
}
_DATASET_P_COL = {
    "crispr": "crispr_p", "gse118713": "gse118713_p", "gse245601": "gse245601_epi_p",
    "gse240112": "gse240112_tumor_p", "gse111151": "gse111151_p",
}


def build_global_ranking(wide_with_percentiles: pd.DataFrame, min_datasets_testable: int = 3) -> pd.DataFrame:
    """Transparent, reconstructable global equal-dataset ranking. Sort
    hierarchy (Phase 9, fixed -- not a weighted formula):
    1. coverage_tier (A>B>C>D>E)
    2. n_datasets_fdr05 (desc)
    3. n_datasets_top10pct (desc)
    4. n_datasets_top20pct (desc)
    5. median_evidence_percentile (desc)
    6. equal_dataset_mean_percentile (desc)
    7. gene symbol (asc) -- deterministic tie-break only
    Genes with fewer than ``min_datasets_testable`` testable datasets are
    excluded from this primary ranking (they remain in the full wide
    table and in high_signal_low_coverage.tsv, never silently discarded
    from the project)."""
    df = wide_with_percentiles.copy()
    percentile_cols = list(_DATASET_PERCENTILE_COL.values())

    df["n_datasets_fdr05"] = sum((df[_DATASET_FDR_COL[d]] < 0.05).fillna(False).astype(int) for d in DATASET_NAMES)
    df["n_datasets_nominal_p05"] = sum((df[_DATASET_P_COL[d]] < 0.05).fillna(False).astype(int) for d in DATASET_NAMES)
    df["n_datasets_top10pct"] = (df[percentile_cols] >= 0.90).sum(axis=1)
    df["n_datasets_top20pct"] = (df[percentile_cols] >= 0.80).sum(axis=1)
    df["median_evidence_percentile"] = df[percentile_cols].median(axis=1, skipna=True)
    df["equal_dataset_mean_percentile"] = df[percentile_cols].mean(axis=1, skipna=True)
    df["min_evidence_percentile"] = df[percentile_cols].min(axis=1, skipna=True)
    df["max_evidence_percentile"] = df[percentile_cols].max(axis=1, skipna=True)

    tier_order = {"A": 5, "B": 4, "C": 3, "D": 2, "E": 1, "none": 0}
    df["_tier_sort"] = df["coverage_tier"].map(tier_order)

    eligible = df.loc[df["n_datasets_testable"] >= min_datasets_testable].copy()
    eligible = eligible.sort_values(
        by=["_tier_sort", "n_datasets_fdr05", "n_datasets_top10pct", "n_datasets_top20pct", "median_evidence_percentile", "equal_dataset_mean_percentile", "gene"],
        ascending=[False, False, False, False, False, False, True],
        kind="mergesort",
    ).reset_index(drop=True)
    eligible["global_rank"] = np.arange(1, len(eligible) + 1)

    df = df.drop(columns=["_tier_sort"])
    logger.info("build_global_ranking: %d genes eligible (>=%d datasets testable) of %d total", len(eligible), min_datasets_testable, len(df))
    return df, eligible


def build_high_signal_low_coverage(df: pd.DataFrame, min_datasets_testable: int, extreme_percentile_threshold: float = 0.99) -> pd.DataFrame:
    """Phase 8: genes below the primary ranking's coverage requirement
    (fewer than ``min_datasets_testable`` datasets testable) that
    nonetheless show an extreme (>=``extreme_percentile_threshold``)
    within-dataset percentile in at least one dataset -- never silently
    hidden just because they lack cross-dataset coverage."""
    percentile_cols = list(_DATASET_PERCENTILE_COL.values())
    low_cov = df.loc[df["n_datasets_testable"] < min_datasets_testable].copy()
    low_cov["max_evidence_percentile"] = low_cov[percentile_cols].max(axis=1, skipna=True)
    has_any = low_cov[percentile_cols].notna().any(axis=1)
    low_cov.loc[has_any, "max_percentile_dataset"] = low_cov.loc[has_any, percentile_cols].idxmax(axis=1)

    out = low_cov.loc[low_cov["max_evidence_percentile"] >= extreme_percentile_threshold].copy()
    out = out.sort_values(by=["max_evidence_percentile", "gene"], ascending=[False, True], kind="mergesort").reset_index(drop=True)
    cols = ["gene", "n_datasets_testable", "coverage_tier", "max_evidence_percentile", "max_percentile_dataset"] + percentile_cols
    logger.info("build_high_signal_low_coverage: %d genes (extreme signal in <%d testable datasets)", len(out), min_datasets_testable)
    return out[cols]


def run_ranking(config_path: str | Path = "config/config.yaml") -> dict[str, pd.DataFrame]:
    config = _load_config(config_path)
    cfg = config["cross_dataset_genomewide"]
    out = cfg["output"]
    min_testable = cfg["min_datasets_testable_for_primary_ranking"]

    wide = pd.read_csv(out["wide_matrix_tsv"], sep="\t")
    with_pct = compute_dataset_percentiles(wide)
    with_tier = assign_coverage_tier(with_pct)
    full, ranked = build_global_ranking(with_tier, min_testable)

    full_path = Path(out["wide_matrix_tsv"]).parent / "all_genes_cross_dataset_evidence_with_ranking.tsv"
    full.to_csv(full_path, sep="\t", index=False)
    ranked_path = Path(out["wide_matrix_tsv"]).parent / "global_ranking_eligible.tsv"
    ranked.to_csv(ranked_path, sep="\t", index=False)

    high_signal_low_coverage = build_high_signal_low_coverage(full, min_testable)
    high_signal_low_coverage.to_csv(out["high_signal_low_coverage_tsv"], sep="\t", index=False)
    logger.info("wrote %s, %s, and %s", full_path, ranked_path, out["high_signal_low_coverage_tsv"])

    return {"full": full, "ranked": ranked, "high_signal_low_coverage": high_signal_low_coverage}


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run_ranking()
