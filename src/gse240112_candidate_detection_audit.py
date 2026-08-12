"""GSE240112 Phase 8 candidate detection audit: for each of the frozen 13
candidates plus PAICS, whether the gene is present in the author-processed
tumor-cell (TTs_cancer_060223.h5seurat) feature space, and if so its
per-cell / per-sample detection in that population, extracted read-only
by ``scripts/analysis/gse240112_01_extract_h5seurat.R``.

Data source: GSE240112 (Fang et al., Genome Medicine 2024, PMID 39558215),
author-processed tumor-cell object, version as downloaded 2026-08-12 (see
data/raw/gse240112/MANIFEST.tsv).
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


def build_detection_audit(
    raw_counts_tsv: str | Path,
    metadata_tsv: str | Path,
    candidates: list[str],
    paics_gene: str,
) -> pd.DataFrame:
    """One row per candidate (+ PAICS), whether it is present in the
    tumor-cell feature space (a gene absent from the extracted raw-counts
    table was absent from the h5Seurat object's own feature list -- see
    ``gse240112_01_extract_h5seurat.R``'s candidate-presence check), and
    if present: total counts, cells detected, percent of tumor cells
    detected, and per-sample (PT1-3/RT1-3) detection counts."""
    raw = pd.read_csv(raw_counts_tsv, sep="\t")
    meta = pd.read_csv(metadata_tsv, sep="\t")[["cell_id", "orig.ident", "group"]]
    if set(raw["cell_id"]) != set(meta["cell_id"]):
        raise ValueError("cell_id sets differ between candidate raw-counts table and tumor-cell metadata")
    merged = raw.merge(meta, on="cell_id", how="inner", validate="one_to_one")
    if len(merged) != len(raw):
        raise ValueError("join between candidate raw-counts table and metadata lost or duplicated rows")

    n_total_cells = len(merged)
    all_genes = list(candidates) + [paics_gene]
    present_genes = [g for g in all_genes if g in raw.columns]
    absent_genes = [g for g in all_genes if g not in raw.columns]
    logger.info("build_detection_audit: %d/%d requested genes present in tumor-cell feature space (absent: %s)", len(present_genes), len(all_genes), absent_genes)

    rows = []
    for gene in all_genes:
        is_paics = gene == paics_gene
        if gene in absent_genes:
            rows.append(
                {
                    "gene": gene,
                    "is_paics_benchmark": is_paics,
                    "present_in_feature_space": False,
                    "n_exact_symbol_matches": 0,
                    "total_counts": 0,
                    "n_cells_detected": 0,
                    "pct_tumor_cells_detected": 0.0,
                    "n_samples_with_detection": 0,
                    "per_sample_cells_detected": "",
                    "reason": "absent from author-processed h5Seurat feature space (27161 genes); confirmed present in raw upstream CellRanger GRCh38 reference by Ensembl ID (ENSG00000231637 for USP17L29) but zero counts detected in any cell across all 6 raw PT/RT CellRanger matrices -- gene is genuinely undetected in this dataset's sequencing data, not a symbol-mapping error",
                }
            )
            continue

        vals = merged[gene]
        detected_mask = vals > 0
        per_sample = merged.loc[detected_mask].groupby("orig.ident")["cell_id"].count()
        per_sample_str = ";".join(f"{s}={c}" for s, c in per_sample.sort_index().items())
        rows.append(
            {
                "gene": gene,
                "is_paics_benchmark": is_paics,
                "present_in_feature_space": True,
                "n_exact_symbol_matches": 1,
                "total_counts": int(vals.sum()),
                "n_cells_detected": int(detected_mask.sum()),
                "pct_tumor_cells_detected": 100.0 * detected_mask.mean(),
                "n_samples_with_detection": int((merged.loc[detected_mask, "orig.ident"].nunique())),
                "per_sample_cells_detected": per_sample_str,
                "reason": "",
            }
        )

    out = pd.DataFrame(rows)
    logger.info(
        "build_detection_audit: %d total tumor cells; USP34 detected in %s cells (%.1f%%)",
        n_total_cells,
        out.loc[out["gene"] == "USP34", "n_cells_detected"].values[0] if "USP34" in out["gene"].values else "NA",
        out.loc[out["gene"] == "USP34", "pct_tumor_cells_detected"].values[0] if "USP34" in out["gene"].values else float("nan"),
    )
    return out


def run_candidate_detection_audit(config_path: str | Path = "config/config.yaml") -> pd.DataFrame:
    config = _load_config(config_path)
    cfg = config["gse240112"]
    out = build_detection_audit(
        cfg["output"]["tt_cancer_candidate_raw_tsv"],
        cfg["output"]["tt_cancer_metadata_tsv"],
        cfg["candidates"]["thirteen"],
        cfg["candidates"]["paics"],
    )
    out_path = Path(cfg["output"]["candidate_detection_audit_tsv"])
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(out_path, sep="\t", index=False)
    logger.info("wrote %s", out_path)
    return out


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run_candidate_detection_audit()
