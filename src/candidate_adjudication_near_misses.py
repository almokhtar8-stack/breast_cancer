"""Candidate adjudication Phase 4: genes that plausibly belong near
MULTIMODAL_STRONG but do not carry that label, found by widening the
Phase 1 rule along the same axes the task names explicitly (near-FDR
CRISPR support, near-FDR resistance support, strong evidence trapped in a
low-coverage gene, or a resistance-percentile leader with a
below-threshold CRISPR p-value) -- never by inventing a new rescue rule
or a numeric composite score. Each near-miss reason is a named boolean
condition, reported explicitly, not hidden inside a score.

Data source: `results/tables/cross_dataset_genomewide/all_genes_cross_dataset_evidence_with_ranking.tsv`,
`resistance_consensus_all_genes.tsv`, `evidence_categories.tsv`,
`high_signal_low_coverage.tsv`, `global_ranking_eligible.tsv` (all frozen,
read-only).
"""

from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd
import yaml

logger = logging.getLogger(__name__)

CRISPR_FDR_SIGNIFICANT = 0.10


def _load_config(config_path: str | Path) -> dict:
    with open(config_path) as f:
        return yaml.safe_load(f)


def find_broad_near_misses(wide: pd.DataFrame, resistance: pd.DataFrame, categories: pd.DataFrame, ranked: pd.DataFrame, high_signal_low_coverage: pd.DataFrame, multimodal_genes: set[str], top_n: int = 20) -> pd.DataFrame:
    merged = wide.merge(resistance, on="gene", how="left").merge(categories, on="gene", how="left").merge(ranked[["gene", "global_rank"]], on="gene", how="left")
    merged = merged.loc[~merged["gene"].isin(multimodal_genes)].copy()

    crispr_fdr = merged["crispr_fdr"]
    resistance_fdr05_count = merged["resistance_fdr05_count"].fillna(0)
    resistance_top10pct_count = merged["resistance_top10pct_count"].fillna(0) if "resistance_top10pct_count" in merged.columns else pd.Series(0, index=merged.index)
    resistance_median_pct = merged["resistance_median_percentile"] if "resistance_median_percentile" in merged.columns else pd.Series(pd.NA, index=merged.index)

    # every gene reaching this point already fails MULTIMODAL_STRONG (the
    # caller pre-filters `multimodal_genes` out), so no extra "would
    # already qualify" exclusion is needed here
    reason_a = (crispr_fdr < 0.15) & (resistance_fdr05_count >= 1)
    reason_b = (crispr_fdr < CRISPR_FDR_SIGNIFICANT) & (resistance_top10pct_count >= 2) & (resistance_fdr05_count < 1)
    reason_c = (merged["evidence_category"] == "FUNCTIONAL_ONLY") & (resistance_median_pct.fillna(0) >= 0.90)

    merged["near_miss_reason_a_crispr_borderline"] = reason_a
    merged["near_miss_reason_b_resistance_top10pct_not_fdr"] = reason_b
    merged["near_miss_reason_c_functional_only_high_resistance_percentile"] = reason_c

    any_reason = reason_a | reason_b | reason_c
    candidates = merged.loc[any_reason].copy()

    def why(row: pd.Series) -> str:
        reasons = []
        if row["near_miss_reason_a_crispr_borderline"]:
            if row["crispr_fdr"] < CRISPR_FDR_SIGNIFICANT:
                reasons.append(f"CRISPR FDR={row['crispr_fdr']:.3f} (<0.10) with exactly 1 resistance dataset FDR<0.05, but overall FDR<0.05 dataset count is only {int(row['n_datasets_fdr05'])} (needs >=2)")
            else:
                reasons.append(f"CRISPR FDR={row['crispr_fdr']:.3f} just above the 0.10 gate, with >=1 resistance dataset FDR<0.05")
        if row["near_miss_reason_b_resistance_top10pct_not_fdr"]:
            reasons.append(">=2 resistance datasets in their own top-10%% but none reaches FDR<0.05")
        if row["near_miss_reason_c_functional_only_high_resistance_percentile"]:
            reasons.append("classified FUNCTIONAL_ONLY, but resistance median percentile >=0.90")
        return "; ".join(reasons)

    candidates["why_not_MULTIMODAL_STRONG"] = candidates.apply(why, axis=1)
    candidates["crispr_evidence"] = candidates.apply(lambda r: f"fdr={r['crispr_fdr']:.4f}, direction={r['crispr_direction']}" if pd.notna(r["crispr_fdr"]) else "not_testable", axis=1)
    candidates["resistance_rna_support"] = candidates.apply(
        lambda r: f"{int(r['resistance_fdr05_count']) if pd.notna(r['resistance_fdr05_count']) else 0}/3 FDR<0.05, consensus={r.get('resistance_direction_consensus', 'NA')}", axis=1
    )
    candidates["human_evidence"] = candidates.apply(
        lambda r: f"gse245601_epi_fdr={r['gse245601_epi_fdr']:.3f}, gse240112_tumor_fdr={r['gse240112_tumor_fdr']:.3f}" if pd.notna(r["gse245601_epi_fdr"]) and pd.notna(r["gse240112_tumor_fdr"]) else "partial", axis=1
    )
    candidates["distance_to_category_rule"] = candidates.apply(
        lambda r: "1 condition short (strict boundary miss)" if (CRISPR_FDR_SIGNIFICANT <= r["crispr_fdr"] < 0.15 and r["resistance_fdr05_count"] >= 1) else "2 or more axes short", axis=1
    )

    out = candidates[["gene", "why_not_MULTIMODAL_STRONG", "crispr_evidence", "resistance_rna_support", "human_evidence", "global_rank", "distance_to_category_rule"]]
    out = out.sort_values(by=["global_rank"], na_position="last").head(top_n).reset_index(drop=True)
    logger.info("find_broad_near_misses: %d near-miss genes surfaced (of %d candidates meeting >=1 named reason)", len(out), int(any_reason.sum()))
    return out


def check_low_coverage_multimodal_pattern(high_signal_low_coverage: pd.DataFrame) -> pd.DataFrame:
    """Tier D/E genes (excluded from MULTIMODAL_STRONG purely by the
    LOW_COVERAGE gate, never by evidence strength) that nonetheless show
    an extreme signal in BOTH a functional (CRISPR) and a resistance-RNA
    column among their few testable datasets."""
    crispr_col = "crispr_evidence_percentile"
    resistance_cols = [c for c in ["gse118713_evidence_percentile", "gse240112_evidence_percentile", "gse111151_evidence_percentile"] if c in high_signal_low_coverage.columns]
    has_crispr = high_signal_low_coverage[crispr_col].notna() & (high_signal_low_coverage[crispr_col] >= 0.90)
    has_resistance = high_signal_low_coverage[resistance_cols].notna().any(axis=1) & (high_signal_low_coverage[resistance_cols].max(axis=1) >= 0.90)
    out = high_signal_low_coverage.loc[has_crispr & has_resistance].copy()
    logger.info("check_low_coverage_multimodal_pattern: %d tier D/E genes show both functional and resistance-RNA extreme signal", len(out))
    return out


def run_near_misses(config_path: str | Path = "config/config.yaml") -> dict[str, pd.DataFrame]:
    config = _load_config(config_path)
    cdx_out = config["cross_dataset_genomewide"]["output"]
    tables_dir = Path(cdx_out["wide_matrix_tsv"]).parent

    wide = pd.read_csv(tables_dir / "all_genes_cross_dataset_evidence_with_ranking.tsv", sep="\t")
    resistance = pd.read_csv(cdx_out["resistance_consensus_tsv"], sep="\t")
    categories = pd.read_csv(tables_dir / "evidence_categories.tsv", sep="\t")
    ranked = pd.read_csv(tables_dir / "global_ranking_eligible.tsv", sep="\t")
    high_signal_low_coverage = pd.read_csv(cdx_out["high_signal_low_coverage_tsv"], sep="\t")
    multimodal_genes = set(config["candidate_adjudication"]["multimodal7"]["genes"])

    near_misses = find_broad_near_misses(wide, resistance, categories, ranked, high_signal_low_coverage, multimodal_genes)
    low_coverage_pattern = check_low_coverage_multimodal_pattern(high_signal_low_coverage)

    out_dir = Path(config["candidate_adjudication"]["output"]["tables_dir"])
    out_dir.mkdir(parents=True, exist_ok=True)
    near_misses.to_csv(out_dir / "multimodal_near_misses.tsv", sep="\t", index=False)
    low_coverage_pattern.to_csv(out_dir / "multimodal_low_coverage_pattern_genes.tsv", sep="\t", index=False)
    return {"near_misses": near_misses, "low_coverage_pattern": low_coverage_pattern}


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run_near_misses()
