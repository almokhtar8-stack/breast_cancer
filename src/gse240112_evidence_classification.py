"""GSE240112 Phase 20: conservative per-candidate evidence classification.
Fixed, documented rule applied uniformly (not hand-picked per gene) to
avoid any post-hoc judgment-call ambiguity:

- untestable: gene not tested in the tumor-cell pseudobulk DE.
- strengthened: candidate-set BH FDR < 0.05.
- directionally_supportive_but_weak: not FDR-significant, but nominal
  p_value < 0.05 OR (|log2fc| >= 0.5 AND the all-epithelial sensitivity
  track agrees in direction) -- a sizeable, cross-track-consistent effect
  even without reaching FDR significance in this small (n=3 vs 3) design.
- neutral_no_additional_support: everything else that was tested.

"Discordant" is reserved for a tested gene whose nominal p_value < 0.05
AND whose direction opposes its own GSE118713 bulk direction (when that
bulk result was itself significant) -- i.e. a real, non-noise-level
disagreement, not a difference between two nonsignificant results.

Data source: GSE240112 (Fang et al., Genome Medicine 2024, PMID 39558215),
tumor-cell candidate table + integrated 4-layer evidence table, version as
downloaded/computed 2026-08-12.
"""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

logger = logging.getLogger(__name__)


def _load_config(config_path: str | Path) -> dict:
    with open(config_path) as f:
        return yaml.safe_load(f)


def classify_candidates(candidate_table: pd.DataFrame, sensitivity: pd.DataFrame, integrated: pd.DataFrame) -> pd.DataFrame:
    sens_indexed = sensitivity.set_index("gene")
    integ_indexed = integrated.set_index("gene")

    rows = []
    for _, row in candidate_table.iterrows():
        gene = row["gene"]
        if not bool(row["tested"]):
            rows.append({"gene": gene, "classification": "untestable", "justification": row["reason_not_tested"]})
            continue

        fdr = float(row["candidate_set_bh_fdr"])
        p = float(row["p_value"])
        log2fc = float(row["log2fc"])

        if fdr < 0.05:
            rows.append({"gene": gene, "classification": "strengthened", "justification": f"candidate-set BH FDR={fdr:.3f} < 0.05"})
            continue

        epi_agree = bool(sens_indexed.loc[gene, "direction_agreement"]) if gene in sens_indexed.index and pd.notna(sens_indexed.loc[gene, "direction_agreement"]) else False
        sizeable_and_consistent = abs(log2fc) >= 0.5 and epi_agree

        bulk_fdr = integ_indexed.loc[gene, "gse118713_tamr_vs_mcf7_fdr"] if gene in integ_indexed.index else np.nan
        bulk_log2fc = integ_indexed.loc[gene, "gse118713_tamr_vs_mcf7_log2fc"] if gene in integ_indexed.index else np.nan
        bulk_significant = pd.notna(bulk_fdr) and bulk_fdr < 0.05
        opposes_significant_bulk = bulk_significant and pd.notna(bulk_log2fc) and (np.sign(log2fc) != np.sign(bulk_log2fc))

        if p < 0.05 and opposes_significant_bulk:
            rows.append(
                {
                    "gene": gene,
                    "classification": "discordant",
                    "justification": f"nominal p={p:.3f} < 0.05 but direction opposes the significant GSE118713 bulk result (FDR={bulk_fdr:.3f})",
                }
            )
            continue

        if p < 0.05 or sizeable_and_consistent:
            rows.append(
                {
                    "gene": gene,
                    "classification": "directionally_supportive_but_weak",
                    "justification": f"not candidate-set FDR-significant (FDR={fdr:.3f}) but nominal p={p:.3f}, |log2fc|={abs(log2fc):.3f}, epithelial-track agreement={epi_agree}",
                }
            )
            continue

        rows.append({"gene": gene, "classification": "neutral_no_additional_support", "justification": f"FDR={fdr:.3f}, nominal p={p:.3f}, no sizeable cross-track-consistent effect"})

    out = pd.DataFrame(rows)
    logger.info("classify_candidates: %s", out["classification"].value_counts().to_dict())
    return out


def run_evidence_classification(config_path: str | Path = "config/config.yaml") -> pd.DataFrame:
    config = _load_config(config_path)
    cfg = config["gse240112"]

    candidate_table = pd.read_csv(cfg["output"]["candidate_table_tsv"], sep="\t")
    sensitivity = pd.read_csv(cfg["output"]["epithelial"]["sensitivity_comparison_tsv"], sep="\t")
    integrated = pd.read_csv(cfg["output"]["integration"]["integrated_tsv"], sep="\t")

    out = classify_candidates(candidate_table, sensitivity, integrated)
    out_path = Path(cfg["output"]["candidate_table_tsv"]).parent / "candidate_classification.tsv"
    out.to_csv(out_path, sep="\t", index=False)
    logger.info("wrote %s", out_path)
    return out


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run_evidence_classification()
