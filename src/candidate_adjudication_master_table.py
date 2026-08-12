"""Candidate adjudication Phases 2-3: the exact, non-summarized evidence
table for the seven MULTIMODAL_STRONG genes, and an independent
value-by-value provenance trace back to each dataset's own original
result file. Every number here is read from an already-frozen output
file (the genome-wide integration's wide matrix/ranking/stability tables,
or a dataset's own raw DE table for provenance and for values -- like
GSE118713's secondary TAMR_vs_FASR contrast, or GSE111151's per-cell-line
sample consistency -- that the wide matrix does not carry).

Data sources: `results/tables/cross_dataset_genomewide/*` (frozen);
`results/tables/gse118713_differential_expression_unredacted.tsv.gz`;
`results/tables/gse111151/genomewide_de.tsv.gz` and its per-sample counts
matrix; `results/tables/gse240112_pseudobulk/*`;
`results/tables/gse245601_pseudobulk/*`. All read-only here.
"""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from src.gse111151_qc import compute_log2cpm, load_counts, load_tmm_norm_factors

logger = logging.getLogger(__name__)


def _load_config(config_path: str | Path) -> dict:
    with open(config_path) as f:
        return yaml.safe_load(f)


def _rank_within_dataset(effect: pd.Series, p_value: pd.Series, fdr: pd.Series, gene: pd.Series) -> pd.Series:
    """Same deterministic sort as `src.cross_dataset_ranking.compute_within_dataset_percentile`
    (fdr asc, p asc, |effect| desc, gene asc), reported as an integer rank
    (1 = strongest) instead of a percentile -- the two are the same
    information, rank is more directly interpretable in a per-gene table."""
    testable = fdr.notna() | p_value.notna()
    n = int(testable.sum())
    out = pd.Series(np.nan, index=effect.index, dtype=float)
    if n == 0:
        return out
    sub = pd.DataFrame({"fdr": fdr[testable], "p": p_value[testable], "_neg_abs_effect": -effect[testable].abs(), "_gene": gene[testable]})
    ordered = sub.sort_values(by=["fdr", "p", "_neg_abs_effect", "_gene"], ascending=True, na_position="last", kind="mergesort")
    rank = pd.Series(np.arange(1, n + 1), index=ordered.index)
    out.loc[rank.index] = rank.values
    return out


def _direction_label(log2fc: float) -> str:
    if pd.isna(log2fc):
        return "not_testable"
    return "up" if log2fc > 0 else ("down" if log2fc < 0 else "unchanged")


def load_gse111151_cell_line_consistency(genes: list[str], config: dict) -> pd.DataFrame:
    """Per gene: how many of the 4 cell lines show the same direction
    (resistant vs its own parental) as the overall cell-line-blocked
    model estimate, using TMM-adjusted log2(CPM+1) at the sample level
    (never naive library-size CPM -- see docs/GSE111151_ANALYSIS_REPORT.md
    limitations for why that distinction matters)."""
    g111 = config["gse111151"]
    counts, metadata, _gene_names = load_counts(g111["output"]["counts_tsv"], g111["output"]["metadata_tsv"])
    effective_lib = load_tmm_norm_factors(g111["output"]["de"]["tmm_norm_factors_tsv"])

    de = pd.read_csv(g111["output"]["de"]["genomewide_tsv"], sep="\t")
    log2cpm = compute_log2cpm(counts, metadata["sample_id"].tolist(), effective_lib)

    rows = []
    for gene in genes:
        de_row = de.loc[de["gene_name"] == gene]
        if len(de_row) == 0:
            rows.append({"gene": gene, "gse111151_consistent_cell_lines_n": np.nan, "gse111151_consistent_cell_lines_total": np.nan})
            continue
        gene_id = de_row["gene_id"].iloc[0]
        overall_log2fc = float(de_row["log2fc"].iloc[0])
        if gene_id not in log2cpm.index or pd.isna(overall_log2fc):
            rows.append({"gene": gene, "gse111151_consistent_cell_lines_n": np.nan, "gse111151_consistent_cell_lines_total": np.nan})
            continue
        vals = log2cpm.loc[gene_id]
        meta_idx = metadata.set_index("sample_id")
        cell_lines = sorted(metadata["cell_line"].unique())
        n_consistent, n_total = 0, 0
        overall_dir = np.sign(overall_log2fc)
        for cl in cell_lines:
            cl_samples = meta_idx.loc[meta_idx["cell_line"] == cl]
            parental = cl_samples.loc[cl_samples["resistance_status"] == "parental"].index
            resistant = cl_samples.loc[cl_samples["resistance_status"] == "resistant"].index
            if len(parental) == 0 or len(resistant) == 0:
                continue
            parental_mean = vals.loc[parental].mean()
            resistant_mean = vals.loc[resistant].mean()
            delta = resistant_mean - parental_mean
            n_total += 1
            if np.sign(delta) == overall_dir and delta != 0:
                n_consistent += 1
        rows.append({"gene": gene, "gse111151_consistent_cell_lines_n": n_consistent, "gse111151_consistent_cell_lines_total": n_total})
    return pd.DataFrame(rows)


def load_gse118713_secondary_contrast(genes: list[str], config: dict) -> pd.DataFrame:
    """TAMR_vs_FASR: secondary context only, never a second independent
    vote -- reported here for completeness, kept in separately-labeled
    columns."""
    path = config["cross_dataset_genomewide"]["inputs"]["gse118713_de_tsv"]
    df = pd.read_csv(path, sep="\t")
    sub = df.loc[(df["contrast"] == "TAMR_vs_FASR") & (df["gene_symbol"].isin(genes))]
    sub = sub.rename(columns={"gene_symbol": "gene", "log2fc": "gse118713_secondary_tamr_vs_fasr_log2fc", "p_value": "gse118713_secondary_tamr_vs_fasr_p", "fdr": "gse118713_secondary_tamr_vs_fasr_fdr"})
    return sub[["gene", "gse118713_secondary_tamr_vs_fasr_log2fc", "gse118713_secondary_tamr_vs_fasr_p", "gse118713_secondary_tamr_vs_fasr_fdr"]].drop_duplicates(subset="gene")


def build_multimodal7_master_table(genes: list[str], config: dict) -> pd.DataFrame:
    cdx_out = config["cross_dataset_genomewide"]["output"]
    tables_dir = Path(cdx_out["wide_matrix_tsv"]).parent

    wide = pd.read_csv(tables_dir / "all_genes_cross_dataset_evidence_with_ranking.tsv", sep="\t")
    ranked = pd.read_csv(tables_dir / "global_ranking_eligible.tsv", sep="\t")
    stability = pd.read_csv(cdx_out["ranking_stability_tsv"], sep="\t")
    categories = pd.read_csv(tables_dir / "evidence_categories.tsv", sep="\t")
    resistance = pd.read_csv(cdx_out["resistance_consensus_tsv"], sep="\t")
    crispr_functional = pd.read_csv(cdx_out["crispr_functional_all_genes_tsv"], sep="\t")

    # per-dataset within-dataset ranks (same key sort as the frozen percentile function),
    # computed on the full wide table (correct denominator) and attached back by gene --
    # never by positional .loc, which would silently misalign once `sub`'s index has been
    # reset by an intervening merge.
    wide = wide.copy()
    wide["gse118713_rank"] = _rank_within_dataset(wide["gse118713_log2fc"], wide["gse118713_p"], wide["gse118713_fdr"], wide["gene"])
    wide["gse240112_rank"] = _rank_within_dataset(wide["gse240112_tumor_log2fc"], wide["gse240112_tumor_p"], wide["gse240112_tumor_fdr"], wide["gene"])
    wide["gse111151_rank"] = _rank_within_dataset(wide["gse111151_log2fc"], wide["gse111151_p"], wide["gse111151_fdr"], wide["gene"])

    sub = wide.loc[wide["gene"].isin(genes)].copy()

    # CRISPR full-screen rank (over all 19,103 testable genes, not the 15,255-gene ranking table)
    crispr_functional = crispr_functional.reset_index(drop=True)
    crispr_functional["crispr_rank_full_screen"] = crispr_functional.index + 1
    sub = sub.merge(crispr_functional[["gene", "crispr_rank_full_screen"]], on="gene", how="left")

    sub = sub.merge(ranked[["gene", "global_rank"]], on="gene", how="left")
    sub = sub.merge(stability[["gene", "stability_label", "n_top20_appearances"]].rename(columns={"stability_label": "global_stability_class", "n_top20_appearances": "global_top20_leave_one_out_appearances"}), on="gene", how="left")
    sub = sub.merge(categories, on="gene", how="left")
    sub = sub.merge(
        resistance[["gene", "resistance_datasets_testable", "resistance_up_count", "resistance_down_count", "resistance_fdr05_count", "resistance_nominal_p05_count", "resistance_direction_consensus"]],
        on="gene", how="left",
    )
    sub = sub.merge(load_gse118713_secondary_contrast(genes, config), on="gene", how="left")
    sub = sub.merge(load_gse111151_cell_line_consistency(genes, config), on="gene", how="left")

    sub["crispr_direction_class"] = sub["crispr_direction"]
    sub["crispr_fdr_lt_005"] = sub["crispr_fdr"] < 0.05
    sub["crispr_fdr_lt_010"] = sub["crispr_fdr"] < 0.10
    sub["crispr_nominal_p_lt_005"] = sub["crispr_p"] < 0.05

    sub["gse118713_direction"] = sub["gse118713_log2fc"].map(_direction_label)
    sub["gse240112_direction"] = sub["gse240112_tumor_log2fc"].map(_direction_label)
    sub["gse111151_direction"] = sub["gse111151_log2fc"].map(_direction_label)
    sub["gse245601_epi_direction"] = sub["gse245601_epi_log2fc"].map(_direction_label)
    sub["gse245601_malignant_direction"] = sub["gse245601_malignant_log2fc"].map(_direction_label)

    def resistance_pattern(row: pd.Series) -> str:
        parts = []
        for d in ("gse118713_direction", "gse240112_direction", "gse111151_direction"):
            v = row[d]
            parts.append(v.upper() if v in ("up", "down") else "NA")
        return "_".join(parts)

    sub["resistance_direction_pattern"] = sub.apply(resistance_pattern, axis=1)

    def human_summary(row: pd.Series) -> tuple[int, str, str]:
        testable = int(row["gse245601_testable"]) + int(row["gse240112_testable"])
        parts = []
        for label, dirn, fdr in (
            ("GSE245601_epi", row["gse245601_epi_direction"], row["gse245601_epi_fdr"]),
            ("GSE240112_tumor", row["gse240112_direction"], row["gse240112_tumor_fdr"]),
        ):
            if dirn == "not_testable":
                parts.append(f"{label}:not_testable")
            else:
                sig = "FDR<0.05" if pd.notna(fdr) and fdr < 0.05 else "not_significant"
                parts.append(f"{label}:{dirn}({sig})")
        pattern = "_".join(p.split(":")[1].split("(")[0].upper() if "not_testable" not in p else "NA" for p in parts)
        return testable, "; ".join(parts), pattern

    human_rows = sub.apply(lambda r: pd.Series(human_summary(r), index=["human_tumor_datasets_testable", "human_tumor_support_summary", "human_tumor_direction_pattern"]), axis=1)
    sub = pd.concat([sub, human_rows], axis=1)

    id_cols = [
        "gene", "global_rank", "coverage_tier", "n_datasets_testable", "n_datasets_fdr05", "n_datasets_nominal_p05",
        "n_datasets_top10pct", "n_datasets_top20pct", "equal_dataset_mean_percentile", "median_evidence_percentile",
        "global_stability_class", "global_top20_leave_one_out_appearances", "evidence_category",
    ]
    crispr_cols = [
        "crispr_effect", "crispr_p", "crispr_fdr", "crispr_rank_full_screen", "crispr_evidence_percentile",
        "crispr_direction_class", "crispr_fdr_lt_005", "crispr_fdr_lt_010", "crispr_nominal_p_lt_005",
    ]
    gse118713_cols = [
        "gse118713_log2fc", "gse118713_p", "gse118713_fdr", "gse118713_rank", "gse118713_evidence_percentile",
        "gse118713_direction", "gse118713_testable", "gse118713_secondary_tamr_vs_fasr_log2fc",
        "gse118713_secondary_tamr_vs_fasr_p", "gse118713_secondary_tamr_vs_fasr_fdr",
    ]
    gse245601_cols = [
        "gse245601_epi_log2fc", "gse245601_epi_p", "gse245601_epi_fdr", "gse245601_evidence_percentile",
        "gse245601_malignant_log2fc", "gse245601_malignant_p", "gse245601_malignant_fdr",
        "gse245601_track_direction_agreement", "gse245601_testable", "gse245601_one_track_only",
    ]
    gse240112_cols = [
        "gse240112_tumor_log2fc", "gse240112_tumor_p", "gse240112_tumor_fdr", "gse240112_rank",
        "gse240112_evidence_percentile", "gse240112_direction", "gse240112_testable",
        "gse240112_epi_log2fc", "gse240112_track_direction_agreement", "gse240112_outlier_fragility",
    ]
    gse111151_cols = [
        "gse111151_log2fc", "gse111151_p", "gse111151_fdr", "gse111151_rank", "gse111151_evidence_percentile",
        "gse111151_direction", "gse111151_testable", "gse111151_consistent_cell_lines_n", "gse111151_consistent_cell_lines_total",
    ]
    resistance_cols = [
        "resistance_datasets_testable", "resistance_fdr05_count", "resistance_nominal_p05_count",
        "resistance_up_count", "resistance_down_count", "resistance_direction_pattern", "resistance_direction_consensus",
    ]
    human_cols = ["human_tumor_datasets_testable", "human_tumor_support_summary", "human_tumor_direction_pattern"]

    all_cols = id_cols + crispr_cols + gse118713_cols + gse245601_cols + gse240112_cols + gse111151_cols + resistance_cols + human_cols
    out = sub[all_cols].sort_values("gene").reset_index(drop=True)
    logger.info("build_multimodal7_master_table: %d genes x %d columns", len(out), len(out.columns))
    return out


def run_master_table(config_path: str | Path = "config/config.yaml") -> pd.DataFrame:
    config = _load_config(config_path)
    genes = config["candidate_adjudication"]["multimodal7"]["genes"]
    out_dir = Path(config["candidate_adjudication"]["output"]["tables_dir"])
    out_dir.mkdir(parents=True, exist_ok=True)
    table = build_multimodal7_master_table(genes, config)
    table.to_csv(out_dir / "multimodal7_exact_evidence.tsv", sep="\t", index=False)
    return table


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run_master_table()
