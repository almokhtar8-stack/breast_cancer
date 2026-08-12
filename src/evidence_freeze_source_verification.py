"""Evidence freeze Phase 19: independently re-verify every displayed
value for the four frozen therapeutic-shortlist genes against each
dataset's own original frozen result file -- a fresh `pd.read_csv`/
`pd.read_parquet` of the original file for each comparison, deliberately
not reusing `src.cross_dataset_gene_mapping` or
`src.evidence_freeze_tables` code, so a shared bug cannot hide itself
from this check (same discipline as
`src.candidate_adjudication_provenance`, extended here to explicitly
include the GSE245601 malignant track and all three datasets' p-values).

Data source: `data/processed/labels.parquet`;
`results/tables/gse118713_differential_expression_unredacted.tsv.gz`;
`results/tables/gse240112_pseudobulk/tumor_cell_genomewide_de.tsv.gz`;
`results/tables/gse111151/genomewide_de.tsv.gz`;
`results/tables/gse245601_pseudobulk/track_{a,b}_genomewide_de.tsv.gz`
(all frozen, read-only).
"""

from __future__ import annotations

import logging
import math
from pathlib import Path

import pandas as pd
import yaml

logger = logging.getLogger(__name__)


def _load_config(config_path: str | Path) -> dict:
    with open(config_path) as f:
        return yaml.safe_load(f)


def _values_match(a: float, b: float) -> bool:
    if pd.isna(a) and pd.isna(b):
        return True
    if pd.isna(a) or pd.isna(b):
        return False
    return math.isclose(float(a), float(b), rel_tol=1e-9, abs_tol=1e-12)


def build_source_value_verification(genes: list[str], full_table: pd.DataFrame, config: dict) -> pd.DataFrame:
    cdx_in = config["cross_dataset_genomewide"]["inputs"]
    ft = full_table.set_index("gene")

    crispr = pd.read_parquet(cdx_in["crispr_labels_parquet"]).set_index("gene")
    gse118713 = pd.read_csv(cdx_in["gse118713_de_tsv"], sep="\t")
    gse118713 = gse118713.loc[gse118713["contrast"] == "TAMR_vs_MCF7"].set_index("gene_symbol")
    gse240112 = pd.read_csv(cdx_in["gse240112_tumor_cell_tsv"], sep="\t").set_index("gene")
    gse111151 = pd.read_csv(cdx_in["gse111151_de_tsv"], sep="\t").set_index("gene_name")
    gse245601_a = pd.read_csv(cdx_in["gse245601_track_a_tsv"], sep="\t").set_index("gene")
    gse245601_b = pd.read_csv(cdx_in["gse245601_track_b_tsv"], sep="\t").set_index("gene")

    rows = []

    def add(gene: str, dataset: str, metric: str, freeze_value, source_file: str, source_value):
        rows.append({"gene": gene, "dataset": dataset, "metric": metric, "freeze_value": freeze_value, "source_file": source_file, "source_value": source_value, "match": _values_match(freeze_value, source_value)})

    for gene in genes:
        f = ft.loc[gene]

        c = crispr.loc[gene]
        add(gene, "crispr", "effect", f["crispr_effect"], "data/processed/labels.parquet", c["effect_size"])
        add(gene, "crispr", "p_value", f["crispr_p"], "data/processed/labels.parquet", c["p_value"])
        add(gene, "crispr", "fdr", f["crispr_fdr"], "data/processed/labels.parquet", c["fdr"])

        g = gse118713.loc[gene]
        add(gene, "gse118713", "log2fc", f["gse118713_log2fc"], "gse118713_differential_expression_unredacted.tsv.gz[TAMR_vs_MCF7]", g["log2fc"])
        add(gene, "gse118713", "p_value", f["gse118713_p"], "gse118713_differential_expression_unredacted.tsv.gz[TAMR_vs_MCF7]", g["p_value"])
        add(gene, "gse118713", "fdr", f["gse118713_fdr"], "gse118713_differential_expression_unredacted.tsv.gz[TAMR_vs_MCF7]", g["fdr"])

        t = gse240112.loc[gene]
        add(gene, "gse240112_tumor", "log2fc", f["gse240112_log2fc"], "gse240112_pseudobulk/tumor_cell_genomewide_de.tsv.gz", t["log2fc"])
        add(gene, "gse240112_tumor", "p_value", f["gse240112_p"], "gse240112_pseudobulk/tumor_cell_genomewide_de.tsv.gz", t["p_value"])
        add(gene, "gse240112_tumor", "fdr", f["gse240112_fdr"], "gse240112_pseudobulk/tumor_cell_genomewide_de.tsv.gz", t["fdr"])

        r = gse111151.loc[gene]
        if isinstance(r, pd.DataFrame):
            r = r.iloc[0]
        add(gene, "gse111151", "log2fc", f["gse111151_log2fc"], "gse111151/genomewide_de.tsv.gz", r["log2fc"])
        add(gene, "gse111151", "p_value", f["gse111151_p"], "gse111151/genomewide_de.tsv.gz", r["p_value"])
        add(gene, "gse111151", "fdr", f["gse111151_fdr"], "gse111151/genomewide_de.tsv.gz", r["fdr"])

        if gene in gse245601_a.index:
            a = gse245601_a.loc[gene]
            add(gene, "gse245601_track_a", "log2fc", f["gse245601_epi_log2fc"], "gse245601_pseudobulk/track_a_genomewide_de.tsv.gz", a["log2fc"])
            add(gene, "gse245601_track_a", "p_value", f["gse245601_epi_p"], "gse245601_pseudobulk/track_a_genomewide_de.tsv.gz", a["p_value"])
            add(gene, "gse245601_track_a", "fdr", f["gse245601_epi_fdr"], "gse245601_pseudobulk/track_a_genomewide_de.tsv.gz", a["fdr"])
        else:
            # gene filtered out (not testable) in this track -- the freeze table must also show NaN, never a fabricated value
            add(gene, "gse245601_track_a", "log2fc", f["gse245601_epi_log2fc"], "gse245601_pseudobulk/track_a_genomewide_de.tsv.gz (gene not in filterByExpr-tested set)", float("nan"))

        if gene in gse245601_b.index:
            b = gse245601_b.loc[gene]
            add(gene, "gse245601_track_b", "log2fc", f["gse245601_malignant_log2fc"], "gse245601_pseudobulk/track_b_genomewide_de.tsv.gz", b["log2fc"])
            add(gene, "gse245601_track_b", "p_value", f["gse245601_malignant_p"], "gse245601_pseudobulk/track_b_genomewide_de.tsv.gz", b["p_value"])
            add(gene, "gse245601_track_b", "fdr", f["gse245601_malignant_fdr"], "gse245601_pseudobulk/track_b_genomewide_de.tsv.gz", b["fdr"])
        else:
            add(gene, "gse245601_track_b", "log2fc", f["gse245601_malignant_log2fc"], "gse245601_pseudobulk/track_b_genomewide_de.tsv.gz (gene not in filterByExpr-tested set)", float("nan"))

    out = pd.DataFrame(rows)
    n_mismatch = int((~out["match"]).sum())
    logger.info("build_source_value_verification: %d comparisons across %d genes, %d mismatches", len(out), len(genes), n_mismatch)
    return out


def run_source_verification(config_path: str | Path = "config/config.yaml") -> pd.DataFrame:
    config = _load_config(config_path)
    out_dir = Path(config["evidence_freeze"]["output"]["tables_dir"])
    full_table = pd.read_csv(out_dir / "final_candidate_evidence.tsv", sep="\t")
    # THERAPEUTIC_SHORTLIST_FREEZE.tsv is the single source of truth for shortlist
    # membership (see src.evidence_freeze_visualization.run_visualization for the
    # same discipline) -- never re-derived from final_candidate_evidence.tsv's own
    # freeze_shortlisted column in case the two tables were regenerated out of order
    freeze_manifest = pd.read_csv(out_dir / "THERAPEUTIC_SHORTLIST_FREEZE.tsv", sep="\t")
    genes = freeze_manifest["gene"].tolist()

    verification = build_source_value_verification(genes, full_table, config)
    verification.to_csv(out_dir / "source_value_verification.tsv", sep="\t", index=False)
    if not verification["match"].all():
        bad = verification.loc[~verification["match"]]
        raise ValueError(f"source-value mismatch(es) found -- freeze must not proceed until investigated:\n{bad}")
    return verification


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run_source_verification()
