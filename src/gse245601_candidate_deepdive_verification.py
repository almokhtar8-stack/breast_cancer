"""GSE245601 candidate deep-dive Phases 2-3: verify the four candidate
genes' presence/detection rates in the epithelial expression matrix, and
reproduce the already-frozen Track A/Track B acute effect (log2FC, p,
FDR) directly from the pseudobulk source files before any new
interpretation is attempted.

Data source: `results/tables/gse245601_candidate_deepdive/candidate_per_cell_expression.tsv`
(this phase's own R extraction), `results/tables/evidence_freeze/final_candidate_evidence.tsv`
(frozen, the reference to reproduce), `results/tables/gse245601_pseudobulk/track_{a,b}_genomewide_de.tsv.gz`
(frozen, read-only).
"""

from __future__ import annotations

import logging
import math
from pathlib import Path

import pandas as pd
import yaml

from src.gse245601_candidate_deepdive_data import GENES, load_per_cell_table

logger = logging.getLogger(__name__)


def _load_config(config_path: str | Path) -> dict:
    with open(config_path) as f:
        return yaml.safe_load(f)


def build_gene_presence_table(per_cell: pd.DataFrame, genes: list[str]) -> pd.DataFrame:
    rows = []
    malignant = per_cell.loc[per_cell["malignancy_status"] == "malignant"]
    nonmalignant = per_cell.loc[per_cell["malignancy_status"] == "non-malignant epithelial"]
    for gene in genes:
        col = f"{gene}_raw_count"
        n_all = len(per_cell)
        n_expr_all = int((per_cell[col] > 0).sum())
        n_expr_mal = int((malignant[col] > 0).sum())
        n_expr_nonmal = int((nonmalignant[col] > 0).sum())
        rows.append(
            {
                "gene": gene, "feature_count_in_object": 1, "total_epithelial_cells": n_all,
                "n_cells_detectable": n_expr_all, "pct_expressing_all_epithelial": 100 * n_expr_all / n_all,
                "n_malignant_cells": len(malignant), "n_expressing_malignant": n_expr_mal,
                "pct_expressing_malignant": 100 * n_expr_mal / len(malignant) if len(malignant) else float("nan"),
                "n_nonmalignant_cells": len(nonmalignant), "n_expressing_nonmalignant": n_expr_nonmal,
                "pct_expressing_nonmalignant": 100 * n_expr_nonmal / len(nonmalignant) if len(nonmalignant) else float("nan"),
            }
        )
    out = pd.DataFrame(rows)
    logger.info("build_gene_presence_table: %s", dict(zip(out["gene"], out["pct_expressing_all_epithelial"].round(1))))
    return out


def _values_close(a: float, b: float, rel_tol: float = 1e-6, abs_tol: float = 1e-9) -> bool:
    if pd.isna(a) and pd.isna(b):
        return True
    if pd.isna(a) or pd.isna(b):
        return False
    return math.isclose(float(a), float(b), rel_tol=rel_tol, abs_tol=abs_tol)


def build_frozen_effect_verification(genes: list[str], config: dict) -> pd.DataFrame:
    """Independently reloads Track A/B genome-wide DE straight from the
    frozen `.tsv.gz` files (not from the evidence-freeze table's own
    already-joined columns) and compares against
    `final_candidate_evidence.tsv`, which was itself independently
    provenance-checked against these same source files during the
    evidence-freeze phase -- a second, independent cross-check here."""
    cfg = config["gse245601_candidate_deepdive"]["inputs"]
    track_a_de = pd.read_csv(cfg["track_a_genomewide_de_tsv"], sep="\t").set_index("gene")
    track_b_de = pd.read_csv(cfg["track_b_genomewide_de_tsv"], sep="\t").set_index("gene")
    frozen = pd.read_csv("results/tables/evidence_freeze/final_candidate_evidence.tsv", sep="\t").set_index("gene")

    rows = []
    for gene in genes:
        f = frozen.loc[gene]
        if gene in track_a_de.index:
            a = track_a_de.loc[gene]
            rows.append({"gene": gene, "track": "A_all_epithelial", "metric": "log2fc", "frozen_value": f["gse245601_epi_log2fc"], "source_value": a["log2fc"], "match": _values_close(f["gse245601_epi_log2fc"], a["log2fc"])})
            rows.append({"gene": gene, "track": "A_all_epithelial", "metric": "p_value", "frozen_value": f["gse245601_epi_p"], "source_value": a["p_value"], "match": _values_close(f["gse245601_epi_p"], a["p_value"])})
            rows.append({"gene": gene, "track": "A_all_epithelial", "metric": "fdr", "frozen_value": f["gse245601_epi_fdr"], "source_value": a["fdr"], "match": _values_close(f["gse245601_epi_fdr"], a["fdr"])})
        else:
            rows.append({"gene": gene, "track": "A_all_epithelial", "metric": "not_tested", "frozen_value": f["gse245601_epi_log2fc"], "source_value": float("nan"), "match": pd.isna(f["gse245601_epi_log2fc"])})

        if gene in track_b_de.index:
            b = track_b_de.loc[gene]
            rows.append({"gene": gene, "track": "B_strict_malignant", "metric": "log2fc", "frozen_value": f["gse245601_malignant_log2fc"], "source_value": b["log2fc"], "match": _values_close(f["gse245601_malignant_log2fc"], b["log2fc"])})
            rows.append({"gene": gene, "track": "B_strict_malignant", "metric": "p_value", "frozen_value": f["gse245601_malignant_p"], "source_value": b["p_value"], "match": _values_close(f["gse245601_malignant_p"], b["p_value"])})
            rows.append({"gene": gene, "track": "B_strict_malignant", "metric": "fdr", "frozen_value": f["gse245601_malignant_fdr"], "source_value": b["fdr"], "match": _values_close(f["gse245601_malignant_fdr"], b["fdr"])})
        else:
            rows.append({"gene": gene, "track": "B_strict_malignant", "metric": "not_tested", "frozen_value": f["gse245601_malignant_log2fc"], "source_value": float("nan"), "match": pd.isna(f["gse245601_malignant_log2fc"])})

    out = pd.DataFrame(rows)
    n_mismatch = int((~out["match"]).sum())
    logger.info("build_frozen_effect_verification: %d comparisons, %d mismatches", len(out), n_mismatch)
    return out


def run_verification(config_path: str | Path = "config/config.yaml") -> dict[str, pd.DataFrame]:
    config = _load_config(config_path)
    out_dir = Path(config["gse245601_candidate_deepdive"]["output"]["tables_dir"])
    out_dir.mkdir(parents=True, exist_ok=True)

    per_cell = load_per_cell_table(config)
    presence = build_gene_presence_table(per_cell, GENES)
    presence.to_csv(out_dir / "gene_presence_verification.tsv", sep="\t", index=False)

    reproduction = build_frozen_effect_verification(GENES, config)
    reproduction.to_csv(out_dir / "frozen_effect_verification.tsv", sep="\t", index=False)
    if not reproduction["match"].all():
        bad = reproduction.loc[~reproduction["match"]]
        raise ValueError(f"frozen GSE245601 effect reproduction FAILED -- do not continue interpretation:\n{bad}")

    return {"presence": presence, "reproduction": reproduction}


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run_verification()
