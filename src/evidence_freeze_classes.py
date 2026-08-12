"""Evidence freeze Phase 13: four separate, non-exclusive frozen
contextual gene lists, so no interesting gene is forced into a role its
evidence doesn't support. Reuses the candidate-adjudication phase's own
already-built, already-Codex-reviewed leaderboards -- never a new sort.

A. THERAPEUTIC INHIBITION SHORTLIST -- sensitising KO required (this
   phase's frozen output, `THERAPEUTIC_SHORTLIST_FREEZE.tsv`).
B. RESISTANCE BIOMARKER/PATHWAY LEADERS -- `shortlist_B_resistance_biomarker.tsv`,
   CRISPR sensitising direction not required.
C. FUNCTIONAL SENSITISATION LEADERS -- `shortlist_C_functional_sensitisation.tsv`,
   RNA support not required.
D. HUMAN-TUMOR LEADERS -- `shortlist_D_human_tumor.tsv`.

Data source: `results/tables/candidate_adjudication/shortlist_{B,C,D}_*.tsv`,
this phase's own `THERAPEUTIC_SHORTLIST_FREEZE.tsv` (read-only).
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


def build_frozen_candidate_classes(freeze_manifest: pd.DataFrame, list_b: pd.DataFrame, list_c: pd.DataFrame, list_d: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for _, r in freeze_manifest.iterrows():
        rows.append({"candidate_class": "A_THERAPEUTIC_INHIBITION_SHORTLIST", "rank_within_class": r["freeze_rank"], "gene": r["gene"], "class_defining_evidence": r["freeze_reason"]})
    for _, r in list_b.iterrows():
        rows.append({"candidate_class": "B_RESISTANCE_BIOMARKER_PATHWAY_LEADERS", "rank_within_class": r["shortlist_rank"], "gene": r["gene"], "class_defining_evidence": f"resistance consensus={r['resistance_direction_consensus']}, {int(r['resistance_fdr05_count'])}/3 FDR<0.05, CRISPR sensitising not required"})
    for _, r in list_c.iterrows():
        rows.append({"candidate_class": "C_FUNCTIONAL_SENSITISATION_LEADERS", "rank_within_class": r["shortlist_rank"], "gene": r["gene"], "class_defining_evidence": f"CRISPR FDR={r['crispr_fdr']:.3g}, sensitising_KO, RNA support not required"})
    for _, r in list_d.iterrows():
        rows.append({"candidate_class": "D_HUMAN_TUMOR_LEADERS", "rank_within_class": r["shortlist_rank"], "gene": r["gene"], "class_defining_evidence": f"human_only_rank={int(r['human_only_rank'])}, {int(r['n_datasets_fdr05'])}/2 human datasets FDR<0.05"})
    out = pd.DataFrame(rows)
    logger.info("build_frozen_candidate_classes: %s", out["candidate_class"].value_counts().to_dict())
    return out


def run_classes(config_path: str | Path = "config/config.yaml") -> pd.DataFrame:
    config = _load_config(config_path)
    adj_tables = Path(config["candidate_adjudication"]["output"]["tables_dir"])
    ef_tables = Path(config["evidence_freeze"]["output"]["tables_dir"])

    freeze_manifest = pd.read_csv(ef_tables / "THERAPEUTIC_SHORTLIST_FREEZE.tsv", sep="\t")
    list_b = pd.read_csv(adj_tables / "shortlist_B_resistance_biomarker.tsv", sep="\t")
    list_c = pd.read_csv(adj_tables / "shortlist_C_functional_sensitisation.tsv", sep="\t")
    list_d = pd.read_csv(adj_tables / "shortlist_D_human_tumor.tsv", sep="\t")

    out = build_frozen_candidate_classes(freeze_manifest, list_b, list_c, list_d)
    out.to_csv(ef_tables / "frozen_candidate_classes.tsv", sep="\t", index=False)
    return out


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run_classes()
