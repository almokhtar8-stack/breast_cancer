"""Candidate adjudication Phase 3: independently re-load every dataset's
own original frozen result file (not the integrated wide matrix) for the
seven MULTIMODAL_STRONG genes, and compare every value the wide matrix
reports against the value found in that original file. This is a
provenance trace, not a re-run of the join/percentile pipeline -- each
comparison here uses a fresh, independent `pd.read_csv` of the original
file and a plain dictionary/row lookup, deliberately not reusing
`src.cross_dataset_gene_mapping` or `src.cross_dataset_evidence_tables`
code, so a bug shared between the builder and the checker cannot hide
itself.

Data sources: the same five original per-dataset files referenced in
`docs/CROSS_DATASET_GENOMEWIDE_DATA_AUDIT.md` (`data/processed/labels.parquet`;
`results/tables/gse118713_differential_expression_unredacted.tsv.gz`;
`results/tables/gse245601_pseudobulk/track_{a,b}_genomewide_de.tsv.gz`;
`results/tables/gse240112_pseudobulk/tumor_cell_genomewide_de.tsv.gz`;
`results/tables/gse111151/genomewide_de.tsv.gz`); all frozen, read-only.
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


def build_provenance_table(genes: list[str], config: dict) -> pd.DataFrame:
    cdx_in = config["cross_dataset_genomewide"]["inputs"]
    cdx_out = config["cross_dataset_genomewide"]["output"]
    wide = pd.read_csv(Path(cdx_out["wide_matrix_tsv"]).parent / "all_genes_cross_dataset_evidence_with_ranking.tsv", sep="\t")
    wide = wide.loc[wide["gene"].isin(genes)].set_index("gene")

    crispr = pd.read_parquet(cdx_in["crispr_labels_parquet"]).set_index("gene")
    gse118713 = pd.read_csv(cdx_in["gse118713_de_tsv"], sep="\t")
    gse118713 = gse118713.loc[gse118713["contrast"] == "TAMR_vs_MCF7"].set_index("gene_symbol")
    gse245601_a = pd.read_csv(cdx_in["gse245601_track_a_tsv"], sep="\t").set_index("gene")
    gse245601_b = pd.read_csv(cdx_in["gse245601_track_b_tsv"], sep="\t").set_index("gene")
    gse240112_tumor = pd.read_csv(cdx_in["gse240112_tumor_cell_tsv"], sep="\t").set_index("gene")
    gse111151 = pd.read_csv(cdx_in["gse111151_de_tsv"], sep="\t").set_index("gene_name")

    rows = []

    def add(gene: str, dataset: str, metric: str, integrated_value, source_file: str, source_value):
        matches = _values_match(integrated_value, source_value)
        rows.append({"gene": gene, "dataset": dataset, "metric": metric, "integrated_value": integrated_value, "source_file": source_file, "source_value": source_value, "matches": matches})

    for gene in genes:
        w = wide.loc[gene]

        if gene in crispr.index:
            c = crispr.loc[gene]
            add(gene, "crispr", "effect", w["crispr_effect"], "data/processed/labels.parquet", c["effect_size"])
            add(gene, "crispr", "p_value", w["crispr_p"], "data/processed/labels.parquet", c["p_value"])
            add(gene, "crispr", "fdr", w["crispr_fdr"], "data/processed/labels.parquet", c["fdr"])

        if gene in gse118713.index:
            g = gse118713.loc[gene]
            add(gene, "gse118713", "log2fc", w["gse118713_log2fc"], "gse118713_differential_expression_unredacted.tsv.gz[TAMR_vs_MCF7]", g["log2fc"])
            add(gene, "gse118713", "p_value", w["gse118713_p"], "gse118713_differential_expression_unredacted.tsv.gz[TAMR_vs_MCF7]", g["p_value"])
            add(gene, "gse118713", "fdr", w["gse118713_fdr"], "gse118713_differential_expression_unredacted.tsv.gz[TAMR_vs_MCF7]", g["fdr"])

        if gene in gse245601_a.index:
            a = gse245601_a.loc[gene]
            add(gene, "gse245601_track_a", "log2fc", w["gse245601_epi_log2fc"], "gse245601_pseudobulk/track_a_genomewide_de.tsv.gz", a["log2fc"])
            add(gene, "gse245601_track_a", "p_value", w["gse245601_epi_p"], "gse245601_pseudobulk/track_a_genomewide_de.tsv.gz", a["p_value"])
            add(gene, "gse245601_track_a", "fdr", w["gse245601_epi_fdr"], "gse245601_pseudobulk/track_a_genomewide_de.tsv.gz", a["fdr"])
        if gene in gse245601_b.index:
            b = gse245601_b.loc[gene]
            add(gene, "gse245601_track_b", "log2fc", w["gse245601_malignant_log2fc"], "gse245601_pseudobulk/track_b_genomewide_de.tsv.gz", b["log2fc"])
            add(gene, "gse245601_track_b", "p_value", w["gse245601_malignant_p"], "gse245601_pseudobulk/track_b_genomewide_de.tsv.gz", b["p_value"])
            add(gene, "gse245601_track_b", "fdr", w["gse245601_malignant_fdr"], "gse245601_pseudobulk/track_b_genomewide_de.tsv.gz", b["fdr"])

        if gene in gse240112_tumor.index:
            t = gse240112_tumor.loc[gene]
            add(gene, "gse240112_tumor", "log2fc", w["gse240112_tumor_log2fc"], "gse240112_pseudobulk/tumor_cell_genomewide_de.tsv.gz", t["log2fc"])
            add(gene, "gse240112_tumor", "p_value", w["gse240112_tumor_p"], "gse240112_pseudobulk/tumor_cell_genomewide_de.tsv.gz", t["p_value"])
            add(gene, "gse240112_tumor", "fdr", w["gse240112_tumor_fdr"], "gse240112_pseudobulk/tumor_cell_genomewide_de.tsv.gz", t["fdr"])

        if gene in gse111151.index:
            de_row = gse111151.loc[gene]
            if isinstance(de_row, pd.DataFrame):
                de_row = de_row.iloc[0]
            add(gene, "gse111151", "log2fc", w["gse111151_log2fc"], "gse111151/genomewide_de.tsv.gz", de_row["log2fc"])
            add(gene, "gse111151", "p_value", w["gse111151_p"], "gse111151/genomewide_de.tsv.gz", de_row["p_value"])
            add(gene, "gse111151", "fdr", w["gse111151_fdr"], "gse111151/genomewide_de.tsv.gz", de_row["fdr"])

    out = pd.DataFrame(rows)
    logger.info("build_provenance_table: %d value comparisons, %d mismatches", len(out), int((~out["matches"]).sum()))
    return out


def run_provenance(config_path: str | Path = "config/config.yaml") -> pd.DataFrame:
    config = _load_config(config_path)
    genes = config["candidate_adjudication"]["multimodal7"]["genes"]
    out_dir = Path(config["candidate_adjudication"]["output"]["tables_dir"])
    out_dir.mkdir(parents=True, exist_ok=True)
    table = build_provenance_table(genes, config)
    table.to_csv(out_dir / "multimodal7_value_provenance.tsv", sep="\t", index=False)
    if not table["matches"].all():
        bad = table.loc[~table["matches"]]
        raise ValueError(f"provenance mismatch(es) found, integration join/transformation must be investigated:\n{bad}")
    return table


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run_provenance()
