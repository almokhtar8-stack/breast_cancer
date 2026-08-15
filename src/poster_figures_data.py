"""Data loaders for the poster figure set (science-first rebuild).

This module performs NO new analysis, NO new candidate discovery, and NO
re-ranking. Every plotted number is read directly from an already-frozen
project source table or an already-validated per-sample/per-cell-line data
matrix:

- data/processed/labels.parquet (genome-wide Hany CRISPR labels, 19,103 genes)
- data/processed/gse118713_gene_tpm.parquet (per-sample TPM matrix)
- results/tables/gse118713_differential_expression.tsv.gz (validated limma DE)
- results/tables/gse240112/candidate_sample_level_log2cpm.tsv (per-tumor pseudobulk)
- results/tables/gse240112/candidate_table.tsv (validated edgeR DE)
- results/tables/evidence_freeze/THERAPEUTIC_SHORTLIST_FREEZE.tsv
- results/tables/cross_dataset_genomewide/all_genes_cross_dataset_evidence_with_ranking.tsv
- results/tables/independent_validation/DepMap_candidate_dependency.tsv (26Q1)
- src/independent_validation_depmap_data.py::load_gene_effect (real per-line
  Chronos values, 26Q1 -- reused unchanged, not recomputed)
- src/final_pharmacogenomics_gdsc_data.py::build_breast_expression_joined
  (real per-cell-line GDSC x DepMap-expression join, reused unchanged)

The two USP34 PDB structures (7W3R, 7W3U) used for Figure 5 were already
established as frozen structural evidence in the final_translational phase;
this module only re-loads the same public coordinate files for rendering.
"""

from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd
import yaml

logger = logging.getLogger(__name__)

CANDIDATES = ["USP34", "VEZF1", "EML5", "CITED2"]
LEAD, BACKUP = "USP34", "VEZF1"

FREEZE_TABLE = Path("results/tables/evidence_freeze/THERAPEUTIC_SHORTLIST_FREEZE.tsv")
CROSS_DATASET_TABLE = Path("results/tables/cross_dataset_genomewide/all_genes_cross_dataset_evidence_with_ranking.tsv")
DEPMAP_TABLE = Path("results/tables/independent_validation/DepMap_candidate_dependency.tsv")
LABELS_PARQUET = Path("data/processed/labels.parquet")
GSE118713_TPM = Path("data/processed/gse118713_gene_tpm.parquet")
GSE118713_DE = Path("results/tables/gse118713_differential_expression.tsv.gz")
GSE240112_SAMPLE_LOG2CPM = Path("results/tables/gse240112/candidate_sample_level_log2cpm.tsv")
GSE240112_CANDIDATE_TABLE = Path("results/tables/gse240112/candidate_table.tsv")


def load_config(config_path: str | Path = "config/config.yaml") -> dict:
    with open(config_path) as f:
        return yaml.safe_load(f)


def load_shortlist_freeze() -> pd.DataFrame:
    """Frozen 4-gene ranking with real Hany CRISPR effect/FDR -- the single
    source of truth for the numbers pinned in tests/test_poster_figures.py.
    """
    df = pd.read_csv(FREEZE_TABLE, sep="\t")
    assert set(df["gene"]) == set(CANDIDATES)
    return df.set_index("gene").loc[CANDIDATES].reset_index()


def load_genomewide_crispr() -> pd.DataFrame:
    """Real genome-wide Hany CRISPR per-gene labels (19,103 genes fitted).
    Same file the Gate-1 decision and every downstream phase was built
    from -- read here only for the discovery-plot visualization, not
    recomputed."""
    df = pd.read_parquet(LABELS_PARQUET)
    logger.info("load_genomewide_crispr: %d genes", len(df))
    return df


def load_gse118713_usp34_samples() -> pd.DataFrame:
    """Real per-sample (3 replicates x 3 conditions) USP34 TPM values from
    the frozen GSE118713 expression matrix, long format."""
    df = pd.read_parquet(GSE118713_TPM)
    row = df[df["gene_symbol"] == "USP34"]
    assert len(row) == 1, f"expected exactly one USP34 row, got {len(row)}"
    sample_cols = [c for c in df.columns if c.startswith(("MCF7_", "TAMR_", "FASR_"))]
    long = row[sample_cols].T.reset_index()
    long.columns = ["sample", "tpm"]
    long["condition"] = long["sample"].str.split("_").str[0]
    logger.info("load_gse118713_usp34_samples: %d samples across %d conditions", len(long), long["condition"].nunique())
    return long


def load_gse118713_usp34_stats() -> pd.DataFrame:
    """Validated limma differential-expression statistics (log2FC/p/FDR)
    for USP34 across all three GSE118713 contrasts."""
    df = pd.read_csv(GSE118713_DE, sep="\t")
    out = df[df["gene_symbol"] == "USP34"][["contrast", "log2fc", "p_value", "fdr", "direction"]]
    assert len(out) == 3
    return out.reset_index(drop=True)


def load_gse240112_vezf1_samples() -> pd.DataFrame:
    """Real per-tumor (pseudobulk) VEZF1 log2CPM values, primary (PT) vs
    recurrent (RT), from the frozen GSE240112 candidate expression table."""
    df = pd.read_csv(GSE240112_SAMPLE_LOG2CPM, sep="\t")
    out = df[df["gene"] == "VEZF1"].reset_index(drop=True)
    assert len(out) == 6, f"expected 6 tumor samples (3 PT + 3 RT), got {len(out)}"
    return out


def load_gse240112_vezf1_stats() -> pd.Series:
    """Validated edgeR differential-expression statistics for VEZF1 in the
    GSE240112 primary-vs-recurrent tumor track."""
    df = pd.read_csv(GSE240112_CANDIDATE_TABLE, sep="\t")
    row = df[df["gene"] == "VEZF1"].iloc[0]
    return row


def load_cross_dataset_evidence() -> pd.DataFrame:
    """Per-dataset (GSE118713, GSE240112 tumor track, GSE111151, GSE245601
    epithelial track) log2FC/p/FDR for the 4 frozen candidates only."""
    df = pd.read_csv(CROSS_DATASET_TABLE, sep="\t").set_index("gene")
    sub = df.loc[CANDIDATES].reset_index()
    logger.info("load_cross_dataset_evidence: %d candidates x %d columns", len(sub), sub.shape[1])
    return sub


def load_depmap_dependency() -> pd.DataFrame:
    """Current (26Q1) DepMap candidate dependency summary table -- asserts
    the active release explicitly (never the archived 24Q4 table)."""
    df = pd.read_csv(DEPMAP_TABLE, sep="\t")
    assert (df["depmap_release"] == "26Q1").all(), "poster figures must use the CURRENT (26Q1) DepMap release, never the archived 24Q4 table"
    assert set(df["candidate"]) == set(CANDIDATES)
    return df.set_index("candidate").loc[CANDIDATES].reset_index()


def load_depmap_gene_effect_er_luminal() -> pd.DataFrame:
    """Real per-cell-line Chronos gene-effect values (26Q1) for the 4
    candidates, restricted to the n=11 ER+/luminal screened lines --
    reuses the already-validated independent_validation loader unchanged."""
    from src.independent_validation_depmap_data import load_gene_effect, load_model

    cfg = load_config()
    model = load_model(cfg, "26Q1")
    eff = load_gene_effect(cfg, "26Q1", CANDIDATES)
    luminal = model.loc[model["is_er_luminal"]]
    joined = eff.join(luminal[["is_er_luminal"]], how="inner")
    assert len(joined) == 11, f"expected 11 ER+/luminal screened lines, got {len(joined)}"
    joined = joined[CANDIDATES].reset_index().rename(columns={"index": "cell_line"})
    long = joined.melt(id_vars="cell_line", var_name="gene", value_name="chronos_effect")
    logger.info("load_depmap_gene_effect_er_luminal: %d lines x %d genes", len(joined), len(CANDIDATES))
    return long


def load_gdsc_usp34_azd7762() -> pd.DataFrame:
    """Real per-cell-line USP34 expression vs AZD7762 (DRUG_ID 1402, GDSC1)
    response, reusing the already-validated final_pharmacogenomics join
    unchanged (the strongest, FDR-significant USP34 GDSC hit)."""
    from src.final_pharmacogenomics_gdsc_data import build_breast_expression_joined

    cfg = load_config()
    df = build_breast_expression_joined(cfg)
    azd = df[(df["DRUG_ID"] == 1402) & (df["DATASET"] == "GDSC1")].copy()
    out = azd[["SANGER_MODEL_ID", "USP34", "LN_IC50", "AUC", "is_er_luminal"]].reset_index(drop=True)
    logger.info("load_gdsc_usp34_azd7762: %d breast cell lines", len(out))
    return out


def load_gdsc_usp34_azd7762_stats() -> pd.Series:
    """Real, already-computed Spearman rho/FDR for USP34 vs AZD7762
    (DRUG_ID 1402, GDSC1, LN_IC50) from the frozen final_pharmacogenomics
    association table -- read here, never recomputed or hand-typed."""
    df = pd.read_csv("results/tables/final_pharmacogenomics/USP34_GDSC_drug_associations.tsv", sep="\t")
    row = df[(df["drug_id"] == 1402) & (df["dataset"] == "GDSC1") & (df["metric"] == "LN_IC50")]
    assert len(row) == 1
    return row.iloc[0]


def usp34_structure_paths(cfg: dict | None = None) -> dict[str, Path]:
    cfg = cfg or load_config()
    d = Path(cfg["data"]["raw"][cfg["poster"]["structures"]["raw_dir_key"]])
    return {pdb_id: d / f"{pdb_id}.pdb" for pdb_id in cfg["poster"]["structures"]["pdb_ids"]}


DATASET_LABELS = {
    "gse118713": "GSE118713 (TAMR/FASR)",
    "gse240112_tumor": "GSE240112 (primary vs recurrent)",
    "gse111151": "GSE111151 (independent resistant)",
    "gse245601_epi": "GSE245601 (acute 12h)",
}


def build_forest_table() -> pd.DataFrame:
    """Tidy log2FC/FDR table for the Figure-2 cross-dataset consistency
    panel: one row per (gene, dataset)."""
    cross = load_cross_dataset_evidence()
    rows = []
    for _, r in cross.iterrows():
        gene = r["gene"]
        for key, value_col, fdr_col in [
            ("gse118713", "gse118713_log2fc", "gse118713_fdr"),
            ("gse240112_tumor", "gse240112_tumor_log2fc", "gse240112_tumor_fdr"),
            ("gse111151", "gse111151_log2fc", "gse111151_fdr"),
            ("gse245601_epi", "gse245601_epi_log2fc", "gse245601_epi_fdr"),
        ]:
            rows.append(dict(
                gene=gene, dataset_key=key, dataset_label=DATASET_LABELS[key],
                log2fc=r[value_col], fdr=r[fdr_col],
            ))
    out = pd.DataFrame(rows)
    logger.info("build_forest_table: %d genes x %d datasets = %d rows", len(CANDIDATES), len(DATASET_LABELS), len(out))
    return out
