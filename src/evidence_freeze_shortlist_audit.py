"""Evidence freeze Phase 2: independently re-verify the provisional
multimodal therapeutic shortlist that the candidate-adjudication phase
actually produced (`results/tables/candidate_adjudication/shortlist_A_multimodal_therapeutic.tsv`),
by re-loading each member's exact values from
`final_candidate_decision_table.tsv` -- never by trusting memory of the
prior session's report. A gene is only `eligible_for_inhibition_strategy`
if its CRISPR direction is `sensitising_KO`; a `tolerance_associated_KO`
gene must never be marked eligible even if it previously appeared on a
shortlist by some other route.

Data source: `results/tables/candidate_adjudication/shortlist_A_multimodal_therapeutic.tsv`,
`final_candidate_decision_table.tsv` (both frozen, committed in `fdd1a44`,
read-only here).
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


def build_shortlist_audit(prior_shortlist: pd.DataFrame, decision_table: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for _, srow in prior_shortlist.iterrows():
        gene = srow["gene"]
        drow = decision_table.loc[decision_table["gene"] == gene]
        if len(drow) == 0:
            raise ValueError(f"shortlisted gene {gene} not found in final_candidate_decision_table.tsv -- cannot verify")
        d = drow.iloc[0]
        eligible = bool(d["crispr_sensitising"])
        rows.append(
            {
                "gene": gene,
                "shortlisted_before_freeze": True,
                "crispr_effect": d["crispr_effect"],
                "crispr_fdr": d["crispr_fdr"],
                "crispr_direction": "sensitising_KO" if d["crispr_sensitising"] else "tolerance_associated_KO",
                "eligible_for_inhibition_strategy": eligible,
                "reason_for_previous_shortlisting": (
                    "one of the two MULTIMODAL_STRONG sensitising genes" if gene in ("USP34", "VEZF1")
                    else "near-miss multimodal gene (sensitising direction + resistance/human RNA support, CRISPR FDR just above the MULTIMODAL_STRONG 0.10 gate)"
                ),
                "resistance_support": f"{int((decision_table.loc[decision_table['gene'] == gene, ['gse118713_fdr', 'gse240112_fdr', 'gse111151_fdr']] < 0.05).sum(axis=1).iloc[0])}/3 FDR<0.05, consensus={d['resistance_consensus']}",
                "human_support": d["human_tumor_summary"],
                "main_limitation": d["main_weakness"],
                "freeze_decision": "CONFIRMED" if eligible else "REJECTED_wrong_crispr_direction",
            }
        )
    out = pd.DataFrame(rows)
    n_wrong_direction = int((~out["eligible_for_inhibition_strategy"]).sum())
    if n_wrong_direction:
        logger.warning("build_shortlist_audit: %d previously-shortlisted gene(s) FAIL the sensitising-KO eligibility gate -- see freeze_decision column", n_wrong_direction)
    logger.info("build_shortlist_audit: %d genes audited, %d eligible for the inhibition strategy", len(out), int(out["eligible_for_inhibition_strategy"].sum()))
    return out


def run_shortlist_audit(config_path: str | Path = "config/config.yaml") -> pd.DataFrame:
    config = _load_config(config_path)
    adj_tables = Path(config["candidate_adjudication"]["output"]["tables_dir"])
    out_dir = Path(config["evidence_freeze"]["output"]["tables_dir"])
    out_dir.mkdir(parents=True, exist_ok=True)

    prior_shortlist = pd.read_csv(adj_tables / "shortlist_A_multimodal_therapeutic.tsv", sep="\t")
    decision_table = pd.read_csv(adj_tables / "final_candidate_decision_table.tsv", sep="\t")

    audit = build_shortlist_audit(prior_shortlist, decision_table)
    audit.to_csv(out_dir / "therapeutic_shortlist_audit.tsv", sep="\t", index=False)
    return audit


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run_shortlist_audit()
