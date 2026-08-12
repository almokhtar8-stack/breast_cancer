"""Extracts the 13 frozen candidates (and, separately, the PAICS
benchmark) from the genome-wide GSE111151 cell-line-blocked edgeR results
(``scripts/analysis/gse111151_02_edger.R``), per
docs/GSE111151_PREANALYSIS.md sections D/H/I.

BH multiple-testing correction is applied to exactly the candidate genes
present in the filtered DE table (genes removed by edgeR's
``filterByExpr`` are excluded from the BH family and reported as
untestable with an explicit reason, never silently dropped). PAICS is a
published benchmark and is never included in that BH family. Candidate
lookup uses Ensembl gene ID (unique in this dataset), not gene symbol
(308 duplicated symbols exist among the 60,619 genes -- see
docs/GSE111151_DATA_AUDIT.md section 4).

Sample-level values use TMM-adjusted log2(CPM+1) (effective library size
= raw library size x edgeR's per-sample norm.factor, exported by
``scripts/analysis/gse111151_02_edger.R``) -- **not** a naive raw-
library-size-only CPM. TMM factors range 0.82-1.30 across these 11
samples; using naive CPM was found to produce a spurious apparent
per-cell-line direction for several candidates (e.g. it made USP34 look
consistently *down* in every cell line even though the properly
TMM-normalized model estimates it *up*) before this was caught and
corrected -- see ``src.gse111151_qc.load_tmm_norm_factors`` and
docs/GSE111151_ANALYSIS_REPORT.md limitations.

Data source: GSE111151 (Hultsch et al., BMC Cancer 2018, PMID 30143015),
cell-line-blocked resistant-vs-parental DE, version as downloaded
2026-08-12.
"""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from src.gse111151_qc import compute_log2cpm, load_counts, load_tmm_norm_factors

logger = logging.getLogger(__name__)

REQUIRED_DE_COLUMNS = ("gene_id", "gene_name", "log2fc", "avg_log_cpm", "p_value", "fdr")


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


def compute_sample_values(counts: pd.DataFrame, metadata: pd.DataFrame, gene_id: str, effective_lib_sizes: pd.Series | None = None) -> pd.DataFrame:
    """Per-sample TMM-adjusted log2(CPM+1) for one gene (by Ensembl ID),
    one row per sample, with cell_line and resistance_status attached."""
    sample_ids = metadata["sample_id"].tolist()
    log2cpm = compute_log2cpm(counts, sample_ids, effective_lib_sizes)
    if gene_id not in log2cpm.index:
        return pd.DataFrame(columns=["sample_id", "cell_line", "resistance_status", "log2cpm"])
    gene_vals = log2cpm.loc[gene_id]
    meta_indexed = metadata.set_index("sample_id")
    rows = [
        {
            "sample_id": s,
            "cell_line": meta_indexed.loc[s, "cell_line"],
            "resistance_status": meta_indexed.loc[s, "resistance_status"],
            "log2cpm": float(gene_vals[s]),
        }
        for s in sample_ids
    ]
    return pd.DataFrame(rows)


def build_candidate_table(
    de_df: pd.DataFrame, counts: pd.DataFrame, metadata: pd.DataFrame, candidate_genes: list[str], candidate_ensembl_ids: dict[str, str]
) -> pd.DataFrame:
    """One row per candidate gene symbol requested, looked up by its
    frozen Ensembl ID. A gene removed by filterByExpr (or, in principle,
    absent from the count matrix -- ruled out for all 14 genes in the
    Phase 1 data audit) is reported as ``tested=False`` with an explicit
    reason, excluded from the candidate-set BH family, never silently
    dropped."""
    rows = []
    de_indexed = de_df.set_index("gene_id")

    gene_id_to_symbol = {candidate_ensembl_ids[g]: g for g in candidate_genes}
    tested_ids = [gid for gid in gene_id_to_symbol if gid in de_indexed.index]
    untested_symbols = [g for g in candidate_genes if candidate_ensembl_ids[g] not in tested_ids]

    if tested_ids:
        tested_p = de_indexed.loc[tested_ids, "p_value"]
        candidate_bh = benjamini_hochberg(tested_p)
    else:
        candidate_bh = pd.Series(dtype=float)

    for gene in candidate_genes:
        gene_id = candidate_ensembl_ids[gene]
        if gene in untested_symbols:
            reason = "gene_id absent from count matrix" if gene_id not in counts.index else "filtered out by edgeR filterByExpr (too lowly expressed to test)"
            rows.append(
                {
                    "gene": gene,
                    "gene_id": gene_id,
                    "tested": False,
                    "reason_not_tested": reason,
                    "log2fc": np.nan,
                    "avg_log_cpm": np.nan,
                    "p_value": np.nan,
                    "genomewide_fdr": np.nan,
                    "candidate_set_bh_fdr": np.nan,
                    "direction": "not_tested",
                }
            )
            continue

        row = de_indexed.loc[gene_id]
        rows.append(
            {
                "gene": gene,
                "gene_id": gene_id,
                "tested": True,
                "reason_not_tested": "",
                "log2fc": float(row["log2fc"]),
                "avg_log_cpm": float(row["avg_log_cpm"]),
                "p_value": float(row["p_value"]),
                "genomewide_fdr": float(row["fdr"]),
                "candidate_set_bh_fdr": float(candidate_bh[gene_id]),
                "direction": "up_in_resistant" if row["log2fc"] > 0 else ("down_in_resistant" if row["log2fc"] < 0 else "unchanged"),
            }
        )

    out = pd.DataFrame(rows)
    logger.info("build_candidate_table: %d/%d candidates tested, %d untested", int(out["tested"].sum()), len(candidate_genes), len(untested_symbols))
    return out


def build_sample_level_table(
    counts: pd.DataFrame, metadata: pd.DataFrame, candidate_genes: list[str], candidate_ensembl_ids: dict[str, str], effective_lib_sizes: pd.Series | None = None
) -> pd.DataFrame:
    rows = []
    for gene in candidate_genes:
        gene_id = candidate_ensembl_ids[gene]
        if gene_id not in counts.index:
            continue
        sample_vals = compute_sample_values(counts, metadata, gene_id, effective_lib_sizes)
        sample_vals["gene"] = gene
        rows.append(sample_vals)
    if not rows:
        return pd.DataFrame(columns=["gene", "sample_id", "cell_line", "resistance_status", "log2cpm"])
    out = pd.concat(rows, ignore_index=True)[["gene", "sample_id", "cell_line", "resistance_status", "log2cpm"]]
    logger.info("build_sample_level_table: %d gene x sample rows", len(out))
    return out


def build_paics_row(de_df: pd.DataFrame, counts: pd.DataFrame, metadata: pd.DataFrame, paics_gene: str, paics_gene_id: str) -> pd.DataFrame:
    de_indexed = de_df.set_index("gene_id")
    if paics_gene_id not in de_indexed.index:
        return pd.DataFrame(
            [{"gene": paics_gene, "gene_id": paics_gene_id, "benchmark_label": "published_benchmark_not_in_13_candidate_bh_family", "tested": False, "log2fc": np.nan, "avg_log_cpm": np.nan, "p_value": np.nan, "genomewide_fdr": np.nan, "direction": "not_tested"}]
        )
    row = de_indexed.loc[paics_gene_id]
    return pd.DataFrame(
        [
            {
                "gene": paics_gene,
                "gene_id": paics_gene_id,
                "benchmark_label": "published_benchmark_not_in_13_candidate_bh_family",
                "tested": True,
                "log2fc": float(row["log2fc"]),
                "avg_log_cpm": float(row["avg_log_cpm"]),
                "p_value": float(row["p_value"]),
                "genomewide_fdr": float(row["fdr"]),
                "direction": "up_in_resistant" if row["log2fc"] > 0 else ("down_in_resistant" if row["log2fc"] < 0 else "unchanged"),
            }
        ]
    )


def run_candidate_extraction(config_path: str | Path = "config/config.yaml") -> dict[str, pd.DataFrame]:
    config = _load_config(config_path)
    cfg = config["gse111151"]
    candidates_13 = cfg["candidates"]["thirteen"]
    paics_gene = cfg["candidates"]["paics"]
    ensembl_ids = cfg["candidate_ensembl_ids"]

    de_df = load_de_table(cfg["output"]["de"]["genomewide_tsv"])
    counts, metadata, _gene_names = load_counts(cfg["output"]["counts_tsv"], cfg["output"]["metadata_tsv"])
    effective_lib_sizes = load_tmm_norm_factors(cfg["output"]["de"]["tmm_norm_factors_tsv"])

    candidate_table = build_candidate_table(de_df, counts, metadata, candidates_13, ensembl_ids)
    sample_level = build_sample_level_table(counts, metadata, candidates_13 + [paics_gene], ensembl_ids, effective_lib_sizes)
    paics_row = build_paics_row(de_df, counts, metadata, paics_gene, ensembl_ids[paics_gene])

    Path(cfg["output"]["candidate_table_tsv"]).parent.mkdir(parents=True, exist_ok=True)
    candidate_table.to_csv(cfg["output"]["candidate_table_tsv"], sep="\t", index=False)
    paics_row.to_csv(cfg["output"]["paics_tsv"], sep="\t", index=False)
    sample_level.to_csv(cfg["output"]["sample_level_tsv"], sep="\t", index=False)
    logger.info("run_candidate_extraction: wrote candidate table, PAICS row, and sample-level table")

    return {"candidate_table": candidate_table, "paics": paics_row, "sample_level": sample_level}


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run_candidate_extraction()
