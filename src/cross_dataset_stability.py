"""Cross-dataset genome-wide integration, Phases 18-19: ranking-stability
sensitivity analysis under several reasonable equal-treatment schemes,
and an explicit leave-one-dataset-out check for the Top-20 genes.

No scheme is treated as "the" ranking -- the main global ranking
(`src.cross_dataset_ranking.build_global_ranking`) remains the primary
answer; this module asks whether that answer would change under
reasonable alternative equal-treatment rules, and whether it survives
losing any single independent dataset.

Data source: `results/tables/cross_dataset_genomewide/all_genes_cross_dataset_evidence_with_ranking.tsv`.
"""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

logger = logging.getLogger(__name__)

DATASET_NAMES = ["crispr", "gse118713", "gse245601", "gse240112", "gse111151"]
PCT_COL = {
    "crispr": "crispr_evidence_percentile", "gse118713": "gse118713_evidence_percentile",
    "gse245601": "gse245601_evidence_percentile", "gse240112": "gse240112_evidence_percentile", "gse111151": "gse111151_evidence_percentile",
}
FDR_COL = {
    "crispr": "crispr_fdr", "gse118713": "gse118713_fdr", "gse245601": "gse245601_epi_fdr",
    "gse240112": "gse240112_tumor_fdr", "gse111151": "gse111151_fdr",
}
TESTABLE_COL = {d: f"{d}_testable" for d in DATASET_NAMES}


def _load_config(config_path: str | Path) -> dict:
    with open(config_path) as f:
        return yaml.safe_load(f)


def _rank_subset(df: pd.DataFrame, dataset_names: list[str], primary: str, min_testable: int = 3) -> pd.Series:
    """Rank genes using only ``dataset_names``' percentile/FDR columns.
    The ``"hierarchy"`` scheme (used for all 5 leave-one-dataset-out
    variants) is deliberately identical, key-for-key, to
    ``cross_dataset_ranking.build_global_ranking``'s main sort hierarchy
    -- (1) n_testable-within-subset (the leave-one-out analogue of
    coverage tier, since tier is itself derived from testable count),
    (2) n_fdr05, (3) n_top10pct, (4) n_top20pct, (5) median percentile,
    (6) mean percentile, (7) gene ascending. An earlier version of this
    function omitted (1) and (4), which the Phase 27 Codex review found
    made the leave-one-out reruns not genuinely comparable to the main
    ranking -- fixed here, not just documented as a known gap. Returns a
    rank Series (1=best) indexed by gene, NaN for genes below
    ``min_testable`` among the given dataset subset."""
    pct_cols = [PCT_COL[d] for d in dataset_names]
    n_testable = df[[TESTABLE_COL[d] for d in dataset_names]].sum(axis=1)
    n_fdr05 = sum((df[FDR_COL[d]] < 0.05).fillna(False).astype(int) for d in dataset_names)
    n_top10 = (df[pct_cols] >= 0.90).sum(axis=1)
    n_top20 = (df[pct_cols] >= 0.80).sum(axis=1)
    median_pct = df[pct_cols].median(axis=1, skipna=True)
    mean_pct = df[pct_cols].mean(axis=1, skipna=True)

    work = pd.DataFrame(
        {"gene": df["gene"], "n_testable": n_testable, "n_fdr05": n_fdr05, "n_top10": n_top10, "n_top20": n_top20, "median_pct": median_pct, "mean_pct": mean_pct}
    )
    eligible = work.loc[work["n_testable"] >= min(min_testable, len(dataset_names))].copy()

    sort_specs = {
        "hierarchy": (["n_testable", "n_fdr05", "n_top10", "n_top20", "median_pct", "mean_pct", "gene"], [False, False, False, False, False, False, True]),
        "median_percentile": (["median_pct", "mean_pct", "gene"], [False, False, True]),
        "mean_percentile": (["mean_pct", "median_pct", "gene"], [False, False, True]),
        "fdr_count_first": (["n_fdr05", "median_pct", "gene"], [False, False, True]),
        "top10pct_count_first": (["n_top10", "median_pct", "gene"], [False, False, True]),
    }
    by, ascending = sort_specs[primary]
    eligible = eligible.sort_values(by=by, ascending=ascending, kind="mergesort").reset_index(drop=True)
    eligible["rank"] = np.arange(1, len(eligible) + 1)
    return eligible.set_index("gene")["rank"]


def build_ranking_stability(df: pd.DataFrame, main_ranked: pd.DataFrame) -> pd.DataFrame:
    """One row per gene present in the main global ranking: rank under
    each of the 4 alternative equal-treatment schemes (median-percentile,
    mean-percentile, FDR-count-first, top10pct-count-first) plus each of
    the 5 leave-one-dataset-out variants, plus summary stats (best/worst/
    median rank, IQR, Top-20/10/5 appearance counts)."""
    schemes: dict[str, pd.Series] = {
        "rank_main": main_ranked.set_index("gene")["global_rank"],
        "rank_scheme_median_percentile": _rank_subset(df, DATASET_NAMES, "median_percentile"),
        "rank_scheme_mean_percentile": _rank_subset(df, DATASET_NAMES, "mean_percentile"),
        "rank_scheme_fdr_count_first": _rank_subset(df, DATASET_NAMES, "fdr_count_first"),
        "rank_scheme_top10pct_count_first": _rank_subset(df, DATASET_NAMES, "top10pct_count_first"),
    }
    for leave_out in DATASET_NAMES:
        remaining = [d for d in DATASET_NAMES if d != leave_out]
        schemes[f"rank_without_{leave_out}"] = _rank_subset(df, remaining, "hierarchy")

    genes = main_ranked["gene"]
    table = pd.DataFrame({"gene": genes})
    for col_name, ranks in schemes.items():
        table[col_name] = table["gene"].map(ranks)

    rank_cols = list(schemes.keys())
    all_ranks = table[rank_cols]
    table["best_rank"] = all_ranks.min(axis=1, skipna=True)
    table["worst_rank"] = all_ranks.max(axis=1, skipna=True)
    table["median_rank"] = all_ranks.median(axis=1, skipna=True)
    table["rank_iqr"] = all_ranks.quantile(0.75, axis=1) - all_ranks.quantile(0.25, axis=1)
    table["n_top20_appearances"] = (all_ranks <= 20).sum(axis=1)
    table["n_top10_appearances"] = (all_ranks <= 10).sum(axis=1)
    table["n_top5_appearances"] = (all_ranks <= 5).sum(axis=1)

    logger.info("build_ranking_stability: %d genes across %d ranking variants", len(table), len(rank_cols))
    return table


def classify_stability(stability: pd.DataFrame) -> pd.Series:
    """ROBUST: stays in the global Top 20 under all 5 leave-one-out
    variants. MODERATELY_STABLE: stays in the Top 20 under >=3 of 5.
    DATASET_DEPENDENT: drops out of the Top 20 under >=3 of 5 (i.e. its
    apparent strength depends heavily on one particular dataset)."""
    leave_one_out_cols = [f"rank_without_{d}" for d in DATASET_NAMES]
    n_top20_loo = (stability[leave_one_out_cols] <= 20).sum(axis=1)
    out = pd.Series("DATASET_DEPENDENT", index=stability.index)
    out.loc[n_top20_loo >= 3] = "MODERATELY_STABLE"
    out.loc[n_top20_loo == 5] = "ROBUST"
    return out


def run_stability(config_path: str | Path = "config/config.yaml") -> pd.DataFrame:
    config = _load_config(config_path)
    cfg = config["cross_dataset_genomewide"]
    out = cfg["output"]
    tables_dir = Path(out["wide_matrix_tsv"]).parent

    df = pd.read_csv(tables_dir / "all_genes_cross_dataset_evidence_with_ranking.tsv", sep="\t")
    main_ranked = pd.read_csv(tables_dir / "global_ranking_eligible.tsv", sep="\t")

    stability = build_ranking_stability(df, main_ranked)
    stability["stability_label"] = classify_stability(stability)

    stability.to_csv(out["ranking_stability_tsv"], sep="\t", index=False)
    logger.info("wrote %s", out["ranking_stability_tsv"])
    return stability


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run_stability()
