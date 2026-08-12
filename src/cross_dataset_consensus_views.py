"""Cross-dataset genome-wide integration, Phases 10-15: alternative,
equally-valid views of the same wide evidence matrix -- resistance-state
RNA consensus (established/chronic resistance datasets only, no CRISPR
requirement), full-genome functional CRISPR ranking, human-only and
RNA-only views, CRISPR-independent discovery, and a single, deterministic,
programmatic evidence-category assignment (unifying Phase 11's functional
-vs-RNA cross-reference and Phase 15's multimodal categories into one
precedence-ordered decision tree, documented below -- not two competing,
overlapping category systems).

None of these views is presented as superior to another (task's Phase 12
instruction, generalized to every view here) -- each answers a different
biological question from the same underlying evidence.

Data source: `results/tables/cross_dataset_genomewide/all_genes_cross_dataset_evidence_with_ranking.tsv`
(built by `src.cross_dataset_ranking`).
"""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

logger = logging.getLogger(__name__)

# Resistance-state datasets: established/chronic acquired-resistance or
# recurrence context only. GSE245601 (acute 12h ex vivo) is explicitly
# excluded from this consensus per the task's Phase 10 instruction.
RESISTANCE_DATASETS = ["gse118713", "gse240112", "gse111151"]
RESISTANCE_LOG2FC_COL = {"gse118713": "gse118713_log2fc", "gse240112": "gse240112_tumor_log2fc", "gse111151": "gse111151_log2fc"}
RESISTANCE_FDR_COL = {"gse118713": "gse118713_fdr", "gse240112": "gse240112_tumor_fdr", "gse111151": "gse111151_fdr"}
RESISTANCE_P_COL = {"gse118713": "gse118713_p", "gse240112": "gse240112_tumor_p", "gse111151": "gse111151_p"}
RESISTANCE_PCT_COL = {"gse118713": "gse118713_evidence_percentile", "gse240112": "gse240112_evidence_percentile", "gse111151": "gse111151_evidence_percentile"}

HUMAN_DATASETS = ["gse245601", "gse240112"]
RNA_DATASETS = ["gse118713", "gse245601", "gse240112", "gse111151"]

CRISPR_FDR_SIGNIFICANT = 0.10  # matches this project's established Gate-1 CRISPR threshold (config.yaml gate1.fdr_threshold), not a new invented cutoff


def _load_config(config_path: str | Path) -> dict:
    with open(config_path) as f:
        return yaml.safe_load(f)


def _direction_consensus_row(directions: list[str]) -> str:
    clean = [d for d in directions if d in ("up", "down")]
    if len(clean) == 0:
        return "insufficient"
    n_up = clean.count("up")
    n_down = clean.count("down")
    if n_up == len(clean):
        return "all_up"
    if n_down == len(clean):
        return "all_down"
    if n_up > n_down:
        return "majority_up"
    if n_down > n_up:
        return "majority_down"
    return "mixed"


def build_resistance_consensus(df: pd.DataFrame) -> pd.DataFrame:
    """Phase 10: uses ONLY gse118713/gse240112/gse111151 (established
    resistance/recurrence context); never requires CRISPR significance;
    never includes GSE245601 (acute response)."""
    out = df[["gene"]].copy()

    directions = pd.DataFrame(index=df.index)
    for d in RESISTANCE_DATASETS:
        log2fc = df[RESISTANCE_LOG2FC_COL[d]]
        directions[d] = np.select([log2fc > 0, log2fc < 0], ["up", "down"], default=None)

    # derived from the representative log2fc columns (not the df's own *_testable
    # flags) so "testable" here stays exactly consistent with the direction/fdr
    # values actually used in this consensus, including GSE240112's tumor-cell-
    # only convention
    out["resistance_datasets_testable"] = directions.notna().sum(axis=1)
    out["resistance_up_count"] = (directions == "up").sum(axis=1)
    out["resistance_down_count"] = (directions == "down").sum(axis=1)
    out["resistance_fdr05_count"] = sum((df[RESISTANCE_FDR_COL[d]] < 0.05).fillna(False).astype(int) for d in RESISTANCE_DATASETS)
    out["resistance_nominal_p05_count"] = sum((df[RESISTANCE_P_COL[d]] < 0.05).fillna(False).astype(int) for d in RESISTANCE_DATASETS)
    out["resistance_top10pct_count"] = sum((df[RESISTANCE_PCT_COL[d]] >= 0.90).fillna(False).astype(int) for d in RESISTANCE_DATASETS)
    out["resistance_top20pct_count"] = sum((df[RESISTANCE_PCT_COL[d]] >= 0.80).fillna(False).astype(int) for d in RESISTANCE_DATASETS)
    out["resistance_median_percentile"] = df[[RESISTANCE_PCT_COL[d] for d in RESISTANCE_DATASETS]].median(axis=1, skipna=True)

    out["resistance_direction_consensus"] = directions.apply(lambda row: _direction_consensus_row(row.tolist()), axis=1)

    out = out.sort_values(
        by=["resistance_fdr05_count", "resistance_direction_consensus", "resistance_top10pct_count", "resistance_top20pct_count", "resistance_median_percentile", "gene"],
        key=lambda col: col.map({"all_up": 4, "all_down": 4, "majority_up": 2, "majority_down": 2, "mixed": 0, "insufficient": -1}) if col.name == "resistance_direction_consensus" else col,
        ascending=[False, False, False, False, False, True],
        kind="mergesort",
    ).reset_index(drop=True)
    logger.info("build_resistance_consensus: %d genes", len(out))
    return out


def build_crispr_functional_ranking(df: pd.DataFrame) -> pd.DataFrame:
    """Phase 11: full CRISPR genome (all 19,103 fitted genes), not only
    prior FDR<0.10 hits, ranked within each direction separately
    (sensitising KO and tolerance-associated KO are not comparable on
    one scale -- opposite biological meanings)."""
    testable = df.loc[df["crispr_testable"]].copy()
    testable["crispr_fdr_significant"] = testable["crispr_fdr"] < CRISPR_FDR_SIGNIFICANT
    out = testable[["gene", "crispr_effect", "crispr_p", "crispr_fdr", "crispr_direction", "crispr_fdr_significant", "crispr_evidence_percentile"]].copy()
    out = out.sort_values(by=["crispr_fdr", "crispr_p", "gene"], ascending=[True, True, True], kind="mergesort").reset_index(drop=True)
    logger.info("build_crispr_functional_ranking: %d testable genes (%d sensitising, %d tolerance)", len(out), (out["crispr_direction"] == "sensitising_KO").sum(), (out["crispr_direction"] == "tolerance_associated_KO").sum())
    return out


def _build_subset_ranking(df: pd.DataFrame, dataset_names: list[str], min_datasets_testable: int = 1) -> pd.DataFrame:
    """Shared helper for Phase 12 (human-only) and Phase 13 (RNA-only):
    same equal-dataset-percentile-mean logic as the global ranking, but
    restricted to the given subset of datasets, with no re-weighting."""
    testable_cols = [f"{d}_testable" for d in dataset_names]
    pct_cols = {
        "crispr": "crispr_evidence_percentile", "gse118713": "gse118713_evidence_percentile",
        "gse245601": "gse245601_evidence_percentile", "gse240112": "gse240112_evidence_percentile", "gse111151": "gse111151_evidence_percentile",
    }
    fdr_cols = {
        "crispr": "crispr_fdr", "gse118713": "gse118713_fdr", "gse245601": "gse245601_epi_fdr",
        "gse240112": "gse240112_tumor_fdr", "gse111151": "gse111151_fdr",
    }
    used_pct_cols = [pct_cols[d] for d in dataset_names]

    out = df[["gene"]].copy()
    out["n_datasets_testable"] = df[testable_cols].sum(axis=1).astype(int)
    out["n_datasets_fdr05"] = sum((df[fdr_cols[d]] < 0.05).fillna(False).astype(int) for d in dataset_names)
    out["median_evidence_percentile"] = df[used_pct_cols].median(axis=1, skipna=True)
    out["equal_dataset_mean_percentile"] = df[used_pct_cols].mean(axis=1, skipna=True)

    eligible = out.loc[out["n_datasets_testable"] >= min_datasets_testable].copy()
    eligible = eligible.sort_values(
        by=["n_datasets_fdr05", "median_evidence_percentile", "equal_dataset_mean_percentile", "gene"],
        ascending=[False, False, False, True], kind="mergesort",
    ).reset_index(drop=True)
    eligible["rank"] = np.arange(1, len(eligible) + 1)
    return eligible


def build_human_only_ranking(df: pd.DataFrame) -> pd.DataFrame:
    """Phase 12: GSE245601 + GSE240112 only -- both human tumor tissue.
    CRISPR (MCF7-derived screen), GSE118713 and GSE111151 (both cell
    lines) are excluded. A robustness perspective, not claimed superior."""
    return _build_subset_ranking(df, HUMAN_DATASETS, min_datasets_testable=1)


def build_rna_only_ranking(df: pd.DataFrame) -> pd.DataFrame:
    """Phase 13: all four RNA datasets, CRISPR excluded entirely."""
    return _build_subset_ranking(df, RNA_DATASETS, min_datasets_testable=2)


def build_crispr_nonsignificant_rna_consensus(df: pd.DataFrame, resistance_consensus: pd.DataFrame) -> pd.DataFrame:
    """Phase 14: genes NOT CRISPR-significant (FDR>=0.10, or untestable
    in CRISPR) but with strong repeated resistance-RNA evidence. Labeled
    `resistance_biomarker_or_pathway_candidates`, never dismissed."""
    merged = df.merge(resistance_consensus, on="gene", how="left")
    crispr_nonsig = (merged["crispr_fdr"] >= CRISPR_FDR_SIGNIFICANT) | merged["crispr_fdr"].isna()

    pattern_a = merged["resistance_fdr05_count"] >= 2
    pattern_b = (merged["resistance_datasets_testable"] == 3) & merged["resistance_direction_consensus"].isin(["all_up", "all_down"]) & (merged["resistance_median_percentile"] >= 0.80)
    pattern_c = (merged["gse240112_evidence_percentile"] >= 0.90).fillna(False) & (merged["gse111151_evidence_percentile"] >= 0.90).fillna(False)
    pattern_d = (merged["resistance_median_percentile"] >= 0.90).fillna(False) & (merged["crispr_effect"].abs() < 1.0).fillna(True)

    out = merged.loc[crispr_nonsig & (pattern_a | pattern_b | pattern_c | pattern_d)].copy()
    out["pattern_A_ge2_fdr05"] = pattern_a.loc[out.index]
    out["pattern_B_all3_same_direction_high_percentile"] = pattern_b.loc[out.index]
    out["pattern_C_human_recurrence_and_independent_cellline"] = pattern_c.loc[out.index]
    out["pattern_D_strong_rna_consensus_weak_crispr"] = pattern_d.loc[out.index]
    out["label"] = "resistance_biomarker_or_pathway_candidate"
    out = out.sort_values(by=["resistance_fdr05_count", "resistance_median_percentile", "gene"], ascending=[False, False, True], kind="mergesort").reset_index(drop=True)
    logger.info("build_crispr_nonsignificant_rna_consensus: %d genes", len(out))
    return out


def assign_evidence_category(df: pd.DataFrame, resistance_consensus: pd.DataFrame) -> pd.DataFrame:
    """Phase 15 (unifying Phase 11's functional-vs-RNA cross-reference):
    ONE deterministic, precedence-ordered decision tree -- evaluated
    top-to-bottom, first match wins, every gene gets exactly one
    category:

    1. LOW_COVERAGE: fewer than 3 datasets testable.
    2. MULTIMODAL_STRONG: CRISPR FDR<0.10 AND (>=2 resistance datasets
       FDR<0.05, OR >=1 resistance dataset FDR<0.05 together with >=2
       datasets overall at FDR<0.05).
    3. RNA_RESISTANCE_CONSENSUS: CRISPR not significant, but resistance
       direction consensus is all/majority AND >=1 resistance dataset
       FDR<0.05.
    4. FUNCTIONAL_ONLY: CRISPR FDR<0.10, no resistance-dataset FDR<0.05,
       and no acute (GSE245601) FDR<0.05 either.
    5. ACUTE_RESPONSE: GSE245601 (Track A) FDR<0.05 or top-10%, but no
       resistance-dataset FDR<0.05.
    6. HUMAN_TUMOR_SUPPORTED: GSE245601 and/or GSE240112 FDR<0.05 or
       top-10%, not already captured above.
    7. CONTEXT_DEPENDENT: >=2 datasets at FDR<0.05 or top-10% overall,
       but resistance direction consensus is "mixed".
    8. LOW_EVIDENCE: everything else with >=3 datasets testable.
    """
    merged = df.merge(resistance_consensus[["gene", "resistance_fdr05_count", "resistance_direction_consensus"]], on="gene", how="left")

    crispr_strong = (merged["crispr_fdr"] < CRISPR_FDR_SIGNIFICANT).fillna(False)
    acute_strong = (merged["gse245601_epi_fdr"] < 0.05).fillna(False) | (merged["gse245601_evidence_percentile"] >= 0.90).fillna(False)
    acute_fdr_strict = (merged["gse245601_epi_fdr"] < 0.05).fillna(False)
    resistance_fdr05 = merged["resistance_fdr05_count"].fillna(0) >= 1
    resistance_fdr05_ge2 = merged["resistance_fdr05_count"].fillna(0) >= 2
    n_fdr05_overall = merged["n_datasets_fdr05"].fillna(0)
    n_top10pct_overall = merged["n_datasets_top10pct"].fillna(0) if "n_datasets_top10pct" in merged.columns else pd.Series(0, index=merged.index)
    human_strong = ((merged["gse245601_epi_fdr"] < 0.05) | (merged["gse245601_evidence_percentile"] >= 0.90) | (merged["gse240112_tumor_fdr"] < 0.05) | (merged["gse240112_evidence_percentile"] >= 0.90)).fillna(False)
    direction_mixed = merged["resistance_direction_consensus"] == "mixed"
    direction_consensus_strong = merged["resistance_direction_consensus"].isin(["all_up", "all_down", "majority_up", "majority_down"])
    low_coverage = merged["n_datasets_testable"] < 3

    category = pd.Series("LOW_EVIDENCE", index=merged.index, dtype=object)
    category.loc[low_coverage] = "LOW_COVERAGE"
    remaining = ~low_coverage

    m = remaining & crispr_strong & (resistance_fdr05_ge2 | (resistance_fdr05 & (n_fdr05_overall >= 2)))
    category.loc[m] = "MULTIMODAL_STRONG"
    remaining &= ~m

    m = remaining & (~crispr_strong) & direction_consensus_strong & resistance_fdr05
    category.loc[m] = "RNA_RESISTANCE_CONSENSUS"
    remaining &= ~m

    m = remaining & crispr_strong & (~resistance_fdr05) & (~acute_fdr_strict)
    category.loc[m] = "FUNCTIONAL_ONLY"
    remaining &= ~m

    m = remaining & acute_strong & (~resistance_fdr05)
    category.loc[m] = "ACUTE_RESPONSE"
    remaining &= ~m

    m = remaining & human_strong
    category.loc[m] = "HUMAN_TUMOR_SUPPORTED"
    remaining &= ~m

    m = remaining & ((n_fdr05_overall >= 2) | (n_top10pct_overall >= 2)) & direction_mixed
    category.loc[m] = "CONTEXT_DEPENDENT"
    remaining &= ~m

    out = df[["gene"]].copy()
    out["evidence_category"] = category
    logger.info("assign_evidence_category: %s", category.value_counts().to_dict())
    return out


def run_consensus_views(config_path: str | Path = "config/config.yaml") -> dict[str, pd.DataFrame]:
    config = _load_config(config_path)
    cfg = config["cross_dataset_genomewide"]
    out = cfg["output"]

    df = pd.read_csv(Path(out["wide_matrix_tsv"]).parent / "all_genes_cross_dataset_evidence_with_ranking.tsv", sep="\t")

    resistance = build_resistance_consensus(df)
    crispr_functional = build_crispr_functional_ranking(df)
    human_only = build_human_only_ranking(df)
    rna_only = build_rna_only_ranking(df)
    crispr_ns_rna = build_crispr_nonsignificant_rna_consensus(df, resistance)
    categories = assign_evidence_category(df, resistance)

    resistance.to_csv(out["resistance_consensus_tsv"], sep="\t", index=False)
    crispr_functional.to_csv(out["crispr_functional_all_genes_tsv"], sep="\t", index=False)
    human_only.to_csv(out["human_only_tsv"], sep="\t", index=False)
    rna_only.to_csv(out["rna_only_tsv"], sep="\t", index=False)
    crispr_ns_rna.to_csv(out["crispr_nonsignificant_rna_consensus_tsv"], sep="\t", index=False)
    categories_path = Path(out["wide_matrix_tsv"]).parent / "evidence_categories.tsv"
    categories.to_csv(categories_path, sep="\t", index=False)
    logger.info("run_consensus_views: wrote 6 tables")

    return {
        "resistance_consensus": resistance, "crispr_functional": crispr_functional, "human_only": human_only,
        "rna_only": rna_only, "crispr_nonsignificant_rna_consensus": crispr_ns_rna, "evidence_categories": categories,
    }


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run_consensus_views()
