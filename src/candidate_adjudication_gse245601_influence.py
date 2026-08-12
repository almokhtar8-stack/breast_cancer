"""Candidate adjudication Phase 18: quantifies why removing GSE245601
disrupts the global Top 20 more than removing any other single dataset
(already observed in the frozen `ranking_stability.tsv`), and checks
whether an alternative GSE245601 track-collapsing scheme (Track A only,
Track B only, instead of the frozen mean-of-both) would change the Top 20
materially. This is a sensitivity audit only -- it does not alter the
primary ranking (Phase 18's explicit instruction).

Data source: `results/tables/cross_dataset_genomewide/all_genes_cross_dataset_evidence_with_ranking.tsv`,
`ranking_stability.tsv`, `top20_global.tsv` (frozen, read-only). Reuses
`src.cross_dataset_ranking.assign_coverage_tier` / `build_global_ranking`
(the same frozen hierarchy) on modified copies of the GSE245601 percentile
/FDR columns -- never a new ranking rule.
"""

from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd
import yaml

from src.cross_dataset_ranking import assign_coverage_tier, build_global_ranking

logger = logging.getLogger(__name__)


def _load_config(config_path: str | Path) -> dict:
    with open(config_path) as f:
        return yaml.safe_load(f)


def compute_top20_overlap_by_leave_one_out(stability: pd.DataFrame) -> pd.DataFrame:
    main_top20 = set(stability.loc[stability["rank_main"] <= 20, "gene"])
    rows = []
    for col, label in [
        ("rank_without_crispr", "crispr"), ("rank_without_gse118713", "gse118713"),
        ("rank_without_gse245601", "gse245601"), ("rank_without_gse240112", "gse240112"),
        ("rank_without_gse111151", "gse111151"),
    ]:
        alt_top20 = set(stability.loc[stability[col] <= 20, "gene"])
        rows.append({"dataset_removed": label, "top20_overlap_with_main": len(main_top20 & alt_top20), "top20_size": 20})
    out = pd.DataFrame(rows).sort_values("top20_overlap_with_main")
    logger.info("compute_top20_overlap_by_leave_one_out: %s", dict(zip(out["dataset_removed"], out["top20_overlap_with_main"])))
    return out


def compute_coverage_dropout_risk(wide: pd.DataFrame) -> pd.DataFrame:
    """For each dataset, how many genes sit at exactly 3/5 testable
    datasets AND count that dataset among the three -- these fall below
    the coverage-tier eligibility floor if that dataset is removed. Only
    relevant to genes near the tier boundary, not to the Top 20 itself
    (which is 5/5-testable Tier A throughout, so removing any one dataset
    only demotes those genes to Tier B, never below the eligibility floor)."""
    dataset_names = ["crispr", "gse118713", "gse245601", "gse240112", "gse111151"]
    testable_cols = {d: f"{d}_testable" for d in dataset_names}
    n_testable = wide[list(testable_cols.values())].sum(axis=1)
    exactly3 = wide.loc[n_testable == 3]
    rows = []
    for d, col in testable_cols.items():
        rows.append({"dataset": d, "n_testable_total": int(wide[col].sum()), "n_genes_at_risk_of_coverage_dropout_if_removed": int(exactly3[col].sum())})
    return pd.DataFrame(rows).sort_values("n_genes_at_risk_of_coverage_dropout_if_removed", ascending=False).reset_index(drop=True)


def compute_median_percentile_shift(wide: pd.DataFrame, top20_genes: list[str]) -> pd.DataFrame:
    """For the actual Top-20 genes, how much would the 5-value median/mean
    percentile shift if GSE245601's percentile were dropped from the set
    -- the tie-break statistic GSE245601 can move even when it contributes
    zero FDR<0.05 votes."""
    cols = ["crispr_evidence_percentile", "gse118713_evidence_percentile", "gse245601_evidence_percentile", "gse240112_evidence_percentile", "gse111151_evidence_percentile"]
    sub = wide.loc[wide["gene"].isin(top20_genes), ["gene"] + cols].set_index("gene")
    rows = []
    for gene, vals in sub.iterrows():
        full_median = vals.dropna().median()
        without = vals.drop("gse245601_evidence_percentile").dropna()
        new_median = without.median()
        rows.append({"gene": gene, "gse245601_percentile": vals["gse245601_evidence_percentile"], "median_percentile_with_gse245601": full_median, "median_percentile_without_gse245601": new_median, "shift": new_median - full_median})
    out = pd.DataFrame(rows).sort_values("shift", key=abs, ascending=False).reset_index(drop=True)
    return out


def build_track_scheme_comparison(wide: pd.DataFrame, main_top20_genes: list[str]) -> pd.DataFrame:
    """Rebuilds the global ranking three times: (1) frozen mean-of-Track-A/B
    (baseline, reproduces top20_global.tsv exactly), (2) Track A only,
    (3) Track B only -- for BOTH the percentile used in the tie-break AND
    the FDR column used for the significance-count keys, since Track A is
    the frozen "representative" for that count. Only the GSE245601 columns
    are substituted; every other dataset's columns and the sort hierarchy
    itself are untouched."""
    variants = {}

    baseline = wide.copy()
    cov = assign_coverage_tier(baseline)
    _, ranked = build_global_ranking(cov)
    variants["frozen_mean_of_A_and_B"] = set(ranked.head(20)["gene"])

    track_a_only = wide.copy()
    track_a_only["gse245601_evidence_percentile"] = track_a_only["gse245601_track_a_percentile"]
    cov = assign_coverage_tier(track_a_only)
    _, ranked = build_global_ranking(cov)
    variants["track_a_only"] = set(ranked.head(20)["gene"])

    track_b_only = wide.copy()
    track_b_only["gse245601_evidence_percentile"] = track_b_only["gse245601_track_b_percentile"]
    track_b_only["gse245601_epi_fdr"] = track_b_only["gse245601_malignant_fdr"]  # so the FDR-count key also reflects Track B, not Track A
    cov = assign_coverage_tier(track_b_only)
    _, ranked = build_global_ranking(cov)
    variants["track_b_only"] = set(ranked.head(20)["gene"])

    rows = []
    main_set = set(main_top20_genes)
    for name, genes in variants.items():
        rows.append({"scheme": name, "top20_overlap_with_frozen_main": len(main_set & genes), "genes": ",".join(sorted(genes))})
    out = pd.DataFrame(rows)
    logger.info("build_track_scheme_comparison: %s", {r["scheme"]: r["top20_overlap_with_frozen_main"] for _, r in out.iterrows()})
    return out


def run_gse245601_influence(config_path: str | Path = "config/config.yaml") -> dict[str, pd.DataFrame]:
    config = _load_config(config_path)
    cdx_out = config["cross_dataset_genomewide"]["output"]
    tables_dir = Path(cdx_out["wide_matrix_tsv"]).parent

    wide = pd.read_csv(tables_dir / "all_genes_cross_dataset_evidence_with_ranking.tsv", sep="\t")
    stability = pd.read_csv(cdx_out["ranking_stability_tsv"], sep="\t")
    top20 = pd.read_csv(tables_dir / "top20_global.tsv", sep="\t")["gene"].tolist()

    overlap = compute_top20_overlap_by_leave_one_out(stability)
    coverage_risk = compute_coverage_dropout_risk(wide)
    median_shift = compute_median_percentile_shift(wide, top20)
    track_comparison = build_track_scheme_comparison(wide, top20)

    out_dir = Path(config["candidate_adjudication"]["output"]["tables_dir"])
    out_dir.mkdir(parents=True, exist_ok=True)
    overlap.to_csv(out_dir / "gse245601_influence_top20_overlap.tsv", sep="\t", index=False)
    coverage_risk.to_csv(out_dir / "gse245601_influence_coverage_risk.tsv", sep="\t", index=False)
    median_shift.to_csv(out_dir / "gse245601_influence_median_shift.tsv", sep="\t", index=False)
    track_comparison.to_csv(out_dir / "gse245601_influence_track_scheme_comparison.tsv", sep="\t", index=False)
    return {"overlap": overlap, "coverage_risk": coverage_risk, "median_shift": median_shift, "track_comparison": track_comparison}


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run_gse245601_influence()
