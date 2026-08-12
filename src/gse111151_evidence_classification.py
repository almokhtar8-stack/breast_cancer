"""GSE111151 Phase 9: conservative per-candidate independent-validation
classification. A nominal p<0.05 alone is never treated as "independently
supported" (explicit task instruction); the fixed, documented rule
considers candidate-set FDR, effect magnitude, and sample-level (cell
-line) consistency together:

- untestable: gene not tested in the cell-line-blocked edgeR model.
- independently_supported: candidate-set BH FDR < 0.05.
- discordant: nominal p_value < 0.05 AND direction opposes a
  *significant* GSE118713 bulk acquired-resistance direction (a real,
  non-noise-level disagreement with the most directly comparable prior
  layer -- not just "different sign from a nonsignificant result").
- directionally_supportive_but_weak: not FDR-significant, but either
  (a) nominal p_value < 0.05 alone, or (b) nominal p_value < 0.3 AND at
  least 3 of the (up to) 4 cell-line blocks show the resistant
  subline(s) shifted in the same direction as the overall model
  coefficient (TMM-adjusted per-sample values, not eyeballed). Cell-line
  consistency is deliberately never used on its own with no p-value
  ceiling: with only 4 blocks, getting 3/4 or 4/4 agreement by chance
  under a true null is not rare (binomial P(>=3 of 4)~31% for a fair
  coin), so requiring a companion nominal signal avoids overclaiming a
  noise pattern as "supportive" -- an earlier draft of this rule used
  consistency alone with no p-value ceiling and let genes with p>0.9
  qualify, caught and corrected before the Phase 14 Codex review.
- neutral_no_additional_support: everything else tested.

Data source: GSE111151 (Hultsch et al., BMC Cancer 2018, PMID 30143015),
candidate table + sample-level table, version as computed 2026-08-12.
"""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

logger = logging.getLogger(__name__)

CELL_LINE_ORDER = ["MCF-7", "T-47D", "ZR-75-1", "BT-474"]


def _load_config(config_path: str | Path) -> dict:
    with open(config_path) as f:
        return yaml.safe_load(f)


def compute_cell_line_consistency(sample_level: pd.DataFrame, gene: str, overall_log2fc: float) -> tuple[int, int]:
    """(n_cell_lines_consistent, n_cell_lines_with_both_arms) for one
    gene: for each of the 4 cell lines with both a parental and >=1
    resistant sample, does the resistant-minus-parental delta (mean of
    resistant sublines) have the same sign as the overall model
    coefficient?"""
    sub = sample_level.loc[sample_level["gene"] == gene]
    overall_positive = overall_log2fc > 0
    n_consistent = 0
    n_with_both = 0
    for cell_line in CELL_LINE_ORDER:
        cl_sub = sub.loc[sub["cell_line"] == cell_line]
        parental = cl_sub.loc[cl_sub["resistance_status"] == "parental", "log2cpm"]
        resistant = cl_sub.loc[cl_sub["resistance_status"] == "resistant", "log2cpm"]
        if len(parental) == 0 or len(resistant) == 0:
            continue
        n_with_both += 1
        delta = resistant.mean() - parental.iloc[0]
        if (delta > 0) == overall_positive:
            n_consistent += 1
    return n_consistent, n_with_both


def classify_candidates(candidate_table: pd.DataFrame, sample_level: pd.DataFrame, gse118713_bulk: pd.DataFrame) -> pd.DataFrame:
    """``gse118713_bulk`` is the already-frozen
    ``results/tables/candidate_evidence_summary.tsv`` (columns
    ``gene_symbol``, ``tamr_vs_mcf7_log2fc``, ``tamr_vs_mcf7_fdr``), read
    directly -- not this project's own not-yet-built 5-layer integration
    table, to avoid a circular dependency (the integration table itself
    depends on this classification's output)."""
    bulk_indexed = gse118713_bulk.set_index("gene_symbol") if len(gse118713_bulk) > 0 else pd.DataFrame()

    rows = []
    for _, row in candidate_table.iterrows():
        gene = row["gene"]
        if not bool(row["tested"]):
            rows.append({"gene": gene, "classification": "untestable", "justification": row["reason_not_tested"], "n_cell_lines_consistent": np.nan, "n_cell_lines_with_both_arms": np.nan})
            continue

        fdr = float(row["candidate_set_bh_fdr"])
        p = float(row["p_value"])
        log2fc = float(row["log2fc"])
        n_consistent, n_with_both = compute_cell_line_consistency(sample_level, gene, log2fc)

        if fdr < 0.05:
            rows.append(
                {
                    "gene": gene,
                    "classification": "independently_supported",
                    "justification": f"candidate-set BH FDR={fdr:.3f} < 0.05",
                    "n_cell_lines_consistent": n_consistent,
                    "n_cell_lines_with_both_arms": n_with_both,
                }
            )
            continue

        bulk_fdr = bulk_indexed.loc[gene, "tamr_vs_mcf7_fdr"] if len(bulk_indexed) > 0 and gene in bulk_indexed.index else np.nan
        bulk_log2fc = bulk_indexed.loc[gene, "tamr_vs_mcf7_log2fc"] if len(bulk_indexed) > 0 and gene in bulk_indexed.index else np.nan
        bulk_significant = pd.notna(bulk_fdr) and bulk_fdr < 0.05
        opposes_significant_bulk = bulk_significant and pd.notna(bulk_log2fc) and (np.sign(log2fc) != np.sign(bulk_log2fc))

        if p < 0.05 and opposes_significant_bulk:
            rows.append(
                {
                    "gene": gene,
                    "classification": "discordant",
                    "justification": f"nominal p={p:.3f} < 0.05 but direction opposes the significant GSE118713 bulk result (FDR={bulk_fdr:.3f})",
                    "n_cell_lines_consistent": n_consistent,
                    "n_cell_lines_with_both_arms": n_with_both,
                }
            )
            continue

        if p < 0.05 or (p < 0.3 and n_consistent >= 3):
            rows.append(
                {
                    "gene": gene,
                    "classification": "directionally_supportive_but_weak",
                    "justification": f"not candidate-set FDR-significant (FDR={fdr:.3f}) but nominal p={p:.3f}, {n_consistent}/{n_with_both} cell lines consistent",
                    "n_cell_lines_consistent": n_consistent,
                    "n_cell_lines_with_both_arms": n_with_both,
                }
            )
            continue

        rows.append(
            {
                "gene": gene,
                "classification": "neutral_no_additional_support",
                "justification": f"FDR={fdr:.3f}, nominal p={p:.3f}, only {n_consistent}/{n_with_both} cell lines consistent",
                "n_cell_lines_consistent": n_consistent,
                "n_cell_lines_with_both_arms": n_with_both,
            }
        )

    out = pd.DataFrame(rows)
    logger.info("classify_candidates: %s", out["classification"].value_counts().to_dict())
    return out


def run_evidence_classification(config_path: str | Path = "config/config.yaml") -> pd.DataFrame:
    config = _load_config(config_path)
    cfg = config["gse111151"]
    ces_cfg = config["candidate_evidence_summary"]

    candidate_table = pd.read_csv(cfg["output"]["candidate_table_tsv"], sep="\t")
    sample_level = pd.read_csv(cfg["output"]["sample_level_tsv"], sep="\t")
    gse118713_bulk = pd.read_csv(ces_cfg["output"]["evidence_summary_tsv"], sep="\t")

    out = classify_candidates(candidate_table, sample_level, gse118713_bulk)
    out_path = Path(cfg["output"]["classification_tsv"])
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(out_path, sep="\t", index=False)
    logger.info("wrote %s", out_path)
    return out


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run_evidence_classification()
