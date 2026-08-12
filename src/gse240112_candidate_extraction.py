"""Extracts the 13 frozen candidates (and, separately, the PAICS
benchmark) from the genome-wide GSE240112 tumor-cell pseudobulk edgeR
results (``scripts/analysis/gse240112_03_pseudobulk_edger.R``), per
docs/GSE240112_PREANALYSIS.md sections D/E/H.

BH multiple-testing correction is applied to exactly the candidate genes
present in the filtered DE table (genes removed by edgeR's
``filterByExpr`` are excluded from the BH family and reported as
untestable with an explicit reason, never silently dropped). PAICS is a
published benchmark and is never included in that BH family. USP17L29 is
absent from the tumor-cell feature space entirely (confirmed genuinely
undetected in the raw sequencing data across all 6 GSE240112 PT/RT
samples -- see gse240112_candidate_detection_audit.py) and is reported as
untestable for the same "never silently dropped" reason.

Sample-level direction uses simple per-sample library-size log2(CPM+1)
from the raw pseudobulk counts (not edgeR's internal TMM-adjusted
values, which are not exported per sample) -- a simpler, supplementary
question ("how many of the 3 PT/3 RT samples moved in the same
direction") than the edgeR group coefficient.

Data source: GSE240112 (Fang et al., Genome Medicine 2024, PMID 39558215),
tumor-cell pseudobulk RT vs PT, version as downloaded 2026-08-12.
"""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from src.gse240112_pseudobulk_qc import compute_log2cpm, load_pseudobulk

logger = logging.getLogger(__name__)

REQUIRED_DE_COLUMNS = ("gene", "log2fc", "avg_log_cpm", "p_value", "fdr")


def _load_config(config_path: str | Path) -> dict:
    with open(config_path) as f:
        return yaml.safe_load(f)


def load_de_table(path: str | Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    missing = [c for c in REQUIRED_DE_COLUMNS if c not in df.columns]
    if missing:
        raise ValueError(f"DE table missing required columns: {missing}")
    logger.info("load_de_table: %d genes tested (%s)", len(df), path)
    return df


def benjamini_hochberg(p_values: pd.Series) -> pd.Series:
    """Standard BH step-up FDR, applied to exactly the p-values passed in."""
    p = p_values.to_numpy(dtype=float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order]
    bh = ranked * n / (np.arange(1, n + 1))
    bh = np.minimum.accumulate(bh[::-1])[::-1]
    bh = np.clip(bh, 0, 1)
    out = np.empty(n)
    out[order] = bh
    return pd.Series(out, index=p_values.index)


def compute_sample_direction(counts: pd.DataFrame, metadata: pd.DataFrame, gene: str) -> pd.DataFrame:
    """Per-sample log2(CPM+1) for one gene, one row per PT/RT sample."""
    sample_ids = metadata["sample_id"].tolist()
    log2cpm = compute_log2cpm(counts, sample_ids)
    if gene not in log2cpm.index:
        return pd.DataFrame(columns=["sample_id", "group", "log2cpm"])
    gene_vals = log2cpm.loc[gene]
    meta_indexed = metadata.set_index("sample_id")
    rows = [{"sample_id": s, "group": meta_indexed.loc[s, "group"], "log2cpm": float(gene_vals[s])} for s in sample_ids]
    return pd.DataFrame(rows)


def build_candidate_table(
    de_df: pd.DataFrame, counts: pd.DataFrame, metadata: pd.DataFrame, candidate_genes: list[str], detection_audit: pd.DataFrame
) -> pd.DataFrame:
    """One row per candidate gene requested. A gene absent from the
    tumor-cell feature space entirely (per the detection audit) or
    removed by edgeR's filterByExpr is reported as ``tested=False`` with
    an explicit reason, excluded from the candidate-set BH family, never
    silently dropped."""
    rows = []
    feature_absent = set(detection_audit.loc[~detection_audit["present_in_feature_space"], "gene"])

    tested = de_df.loc[de_df["gene"].isin(candidate_genes)].copy()
    tested_genes = set(tested["gene"])
    untested_genes = [g for g in candidate_genes if g not in tested_genes]

    if len(tested) > 0:
        tested["candidate_set_bh_fdr"] = benjamini_hochberg(tested["p_value"])
    else:
        tested["candidate_set_bh_fdr"] = pd.Series(dtype=float)

    for gene in candidate_genes:
        if gene in untested_genes:
            if gene in feature_absent:
                reason = "gene absent from tumor-cell feature space (confirmed genuinely undetected in raw sequencing data, not a symbol-mapping error -- see candidate_detection_audit.tsv)"
            else:
                reason = "filtered out by edgeR filterByExpr (too lowly expressed in tumor-cell pseudobulk to test)"
            rows.append(
                {
                    "gene": gene,
                    "tested": False,
                    "reason_not_tested": reason,
                    "log2fc": np.nan,
                    "avg_log_cpm": np.nan,
                    "p_value": np.nan,
                    "genomewide_fdr": np.nan,
                    "candidate_set_bh_fdr": np.nan,
                    "direction": "not_tested",
                    "pt_log2cpm_mean": np.nan,
                    "rt_log2cpm_mean": np.nan,
                    "all_rt_above_pt_range": False,
                    "all_rt_below_pt_range": False,
                }
            )
            continue

        row = tested.loc[tested["gene"] == gene].iloc[0]
        sample_dir = compute_sample_direction(counts, metadata, gene)
        pt_vals = sample_dir.loc[sample_dir["group"] == "PT", "log2cpm"]
        rt_vals = sample_dir.loc[sample_dir["group"] == "RT", "log2cpm"]
        rows.append(
            {
                "gene": gene,
                "tested": True,
                "reason_not_tested": "",
                "log2fc": float(row["log2fc"]),
                "avg_log_cpm": float(row["avg_log_cpm"]),
                "p_value": float(row["p_value"]),
                "genomewide_fdr": float(row["fdr"]),
                "candidate_set_bh_fdr": float(row["candidate_set_bh_fdr"]),
                "direction": "up_in_RT" if row["log2fc"] > 0 else ("down_in_RT" if row["log2fc"] < 0 else "unchanged"),
                "pt_log2cpm_mean": float(pt_vals.mean()),
                "rt_log2cpm_mean": float(rt_vals.mean()),
                # simple non-parametric separation check: does every RT sample sit entirely above (or below) the full PT range -- a single outlier sample cannot produce this on its own
                "all_rt_above_pt_range": bool(rt_vals.min() > pt_vals.max()),
                "all_rt_below_pt_range": bool(rt_vals.max() < pt_vals.min()),
            }
        )

    out = pd.DataFrame(rows)
    logger.info(
        "build_candidate_table: %d/%d candidates tested, %d untested",
        int(out["tested"].sum()),
        len(candidate_genes),
        len(untested_genes),
    )
    return out


def build_sample_level_table(counts: pd.DataFrame, metadata: pd.DataFrame, candidate_genes: list[str]) -> pd.DataFrame:
    """Long-format per-sample log2(CPM+1) for every candidate gene present
    in the pseudobulk count matrix -- one row per (gene, sample). Genes
    absent from the count matrix are skipped (already reported as
    untested in build_candidate_table)."""
    rows = []
    for gene in candidate_genes:
        if gene not in counts.index:
            continue
        sample_dir = compute_sample_direction(counts, metadata, gene)
        sample_dir["gene"] = gene
        rows.append(sample_dir)
    if not rows:
        return pd.DataFrame(columns=["gene", "sample_id", "group", "log2cpm"])
    out = pd.concat(rows, ignore_index=True)[["gene", "sample_id", "group", "log2cpm"]]
    logger.info("build_sample_level_table: %d gene x sample rows", len(out))
    return out


def build_paics_row(de_df: pd.DataFrame, counts: pd.DataFrame, metadata: pd.DataFrame, paics_gene: str) -> pd.DataFrame:
    tested = de_df.loc[de_df["gene"] == paics_gene]
    if len(tested) == 0:
        return pd.DataFrame(
            [
                {
                    "gene": paics_gene,
                    "benchmark_label": "published_benchmark_not_in_13_candidate_bh_family",
                    "tested": False,
                    "log2fc": np.nan,
                    "avg_log_cpm": np.nan,
                    "p_value": np.nan,
                    "genomewide_fdr": np.nan,
                    "direction": "not_tested",
                    "pt_log2cpm_mean": np.nan,
                    "rt_log2cpm_mean": np.nan,
                }
            ]
        )
    row = tested.iloc[0]
    sample_dir = compute_sample_direction(counts, metadata, paics_gene)
    pt_vals = sample_dir.loc[sample_dir["group"] == "PT", "log2cpm"]
    rt_vals = sample_dir.loc[sample_dir["group"] == "RT", "log2cpm"]
    return pd.DataFrame(
        [
            {
                "gene": paics_gene,
                "benchmark_label": "published_benchmark_not_in_13_candidate_bh_family",
                "tested": True,
                "log2fc": float(row["log2fc"]),
                "avg_log_cpm": float(row["avg_log_cpm"]),
                "p_value": float(row["p_value"]),
                "genomewide_fdr": float(row["fdr"]),
                "direction": "up_in_RT" if row["log2fc"] > 0 else ("down_in_RT" if row["log2fc"] < 0 else "unchanged"),
                "pt_log2cpm_mean": float(pt_vals.mean()),
                "rt_log2cpm_mean": float(rt_vals.mean()),
            }
        ]
    )


def run_candidate_extraction(config_path: str | Path = "config/config.yaml") -> dict[str, pd.DataFrame]:
    config = _load_config(config_path)
    cfg = config["gse240112"]
    candidates_13 = cfg["candidates"]["thirteen"]
    paics_gene = cfg["candidates"]["paics"]

    de_df = load_de_table(cfg["output"]["de"]["genomewide_tsv"])
    counts, metadata = load_pseudobulk(cfg["output"]["tumor_cell"]["counts_tsv"], cfg["output"]["tumor_cell"]["metadata_tsv"])
    detection_audit = pd.read_csv(cfg["output"]["candidate_detection_audit_tsv"], sep="\t")

    candidate_table = build_candidate_table(de_df, counts, metadata, candidates_13, detection_audit)
    sample_level = build_sample_level_table(counts, metadata, candidates_13 + [paics_gene])
    paics_row = build_paics_row(de_df, counts, metadata, paics_gene)

    Path(cfg["output"]["candidate_table_tsv"]).parent.mkdir(parents=True, exist_ok=True)
    candidate_table.to_csv(cfg["output"]["candidate_table_tsv"], sep="\t", index=False)
    paics_row.to_csv(cfg["output"]["paics_tsv"], sep="\t", index=False)
    sample_level_path = Path(cfg["output"]["candidate_table_tsv"]).parent / "candidate_sample_level_log2cpm.tsv"
    sample_level.to_csv(sample_level_path, sep="\t", index=False)
    logger.info("run_candidate_extraction: wrote candidate table, PAICS row, and sample-level table")

    return {"candidate_table": candidate_table, "paics": paics_row, "sample_level": sample_level}


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run_candidate_extraction()
