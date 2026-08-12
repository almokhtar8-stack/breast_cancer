"""GSE240112 Phase 14 secondary sensitivity analysis: compares the
primary tumor-cell (malignant-only) pseudobulk RT-vs-PT result against an
all-epithelial-cell pseudobulk RT-vs-PT result (reconstructed from the
raw GEO Cell Ranger matrices via
``scripts/analysis/gse240112_04_cellranger_epithelial.R``, since no
author-provided all-cell-type PT/RT object exists -- see
docs/GSE240112_DATA_AUDIT.md section 4). This reduces dependence on the
Phase 7 tumor-cell definition for the headline RT-vs-PT direction; it is
never used to retroactively change the Phase 7 labels.

Data source: GSE240112 (Fang et al., Genome Medicine 2024, PMID 39558215),
tumor-cell and all-epithelial pseudobulk, version as downloaded
2026-08-12.
"""

from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd
import yaml

logger = logging.getLogger(__name__)


def _load_config(config_path: str | Path) -> dict:
    with open(config_path) as f:
        return yaml.safe_load(f)


def build_sensitivity_comparison(
    tumor_candidate_table: pd.DataFrame, epithelial_de: pd.DataFrame, candidate_genes: list[str]
) -> pd.DataFrame:
    """One row per candidate: tumor-cell log2FC/FDR (from the frozen
    primary candidate table) alongside the all-epithelial log2FC/nominal
    p-value (genome-wide FDR; the epithelial track is not a preregistered
    inferential family and gets no separate candidate-set BH), plus a
    simple direction-agreement flag. A gene untested in either track is
    reported with NA, not dropped."""
    tumor_indexed = tumor_candidate_table.set_index("gene")
    epi_indexed = epithelial_de.set_index("gene")

    rows = []
    for gene in candidate_genes:
        tumor_tested = gene in tumor_indexed.index and bool(tumor_indexed.loc[gene, "tested"])
        epi_tested = gene in epi_indexed.index

        tumor_log2fc = float(tumor_indexed.loc[gene, "log2fc"]) if tumor_tested else float("nan")
        tumor_fdr = float(tumor_indexed.loc[gene, "candidate_set_bh_fdr"]) if tumor_tested else float("nan")
        epi_log2fc = float(epi_indexed.loc[gene, "log2fc"]) if epi_tested else float("nan")
        epi_fdr = float(epi_indexed.loc[gene, "fdr"]) if epi_tested else float("nan")

        if tumor_tested and epi_tested:
            direction_agreement = (tumor_log2fc > 0) == (epi_log2fc > 0)
        else:
            direction_agreement = None

        rows.append(
            {
                "gene": gene,
                "tumor_cell_tested": tumor_tested,
                "tumor_cell_log2fc": tumor_log2fc,
                "tumor_cell_candidate_bh_fdr": tumor_fdr,
                "all_epithelial_tested": epi_tested,
                "all_epithelial_log2fc": epi_log2fc,
                "all_epithelial_genomewide_fdr": epi_fdr,
                "direction_agreement": direction_agreement,
            }
        )
    out = pd.DataFrame(rows)
    n_agree = out["direction_agreement"].sum() if out["direction_agreement"].notna().any() else 0
    n_comparable = out["direction_agreement"].notna().sum()
    logger.info("build_sensitivity_comparison: %d/%d comparable candidates agree in direction", n_agree, n_comparable)
    return out


def run_epithelial_sensitivity(config_path: str | Path = "config/config.yaml") -> pd.DataFrame:
    config = _load_config(config_path)
    cfg = config["gse240112"]
    candidates = cfg["candidates"]["thirteen"]

    tumor_candidate_table = pd.read_csv(cfg["output"]["candidate_table_tsv"], sep="\t")
    epithelial_de = pd.read_csv(cfg["output"]["epithelial"]["genomewide_de_tsv"], sep="\t")

    out = build_sensitivity_comparison(tumor_candidate_table, epithelial_de, candidates)
    out_path = Path(cfg["output"]["epithelial"]["sensitivity_comparison_tsv"])
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(out_path, sep="\t", index=False)
    logger.info("wrote %s", out_path)
    return out


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run_epithelial_sensitivity()
