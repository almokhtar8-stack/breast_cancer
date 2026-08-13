"""Systems-network phase 5: resistance pathway consensus.

Integrates ONLY the three chronic-resistance/recurrence datasets
(GSE118713, GSE240112, GSE111151). GSE245601 (acute 12h response) is
deliberately excluded from this consensus calculation -- it is a
biologically distinct context (see docs/SYSTEMS_NETWORK_INPUT_AUDIT.md /
CLAUDE.md dataset semantics) and is shown only alongside the consensus,
never folded into it.
"""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

logger = logging.getLogger(__name__)

RESISTANCE_DATASETS = ["gse118713", "gse240112", "gse111151"]


def _load_config(config_path: str | Path) -> dict:
    with open(config_path) as f:
        return yaml.safe_load(f)


def _classify(row: pd.Series) -> str:
    """Codex review fix: the original rule granted DIRECTIONAL_CONSENSUS to
    any same-sign pathway tested in >=2 datasets with NO significance
    requirement at all (1,972 pathways, 1,432 of them with zero datasets at
    FDR<0.05) -- same-sign NES alone is not evidence of reproducibility.
    DIRECTIONAL_CONSENSUS now additionally requires nominal p<0.05 in >=2
    datasets. The original rule also checked `datasets_tested>=2` before
    `datasets_FDR05==1`, making SINGLE_DATASET unreachable for any pathway
    tested in more than one dataset even if only one was ever significant;
    SINGLE_DATASET is now checked as a fallback whenever the (now stricter)
    DIRECTIONAL_CONSENSUS bar is not met, regardless of datasets_tested."""
    if row["direction_consistency"] == "mixed":
        return "MIXED"
    if row["datasets_FDR05"] >= 2:
        return "STRONG_CONSENSUS"
    if row["datasets_tested"] >= 2 and row["datasets_nominal"] >= 2:
        return "DIRECTIONAL_CONSENSUS"
    if row["datasets_FDR05"] == 1:
        return "SINGLE_DATASET"
    return "LOW_EVIDENCE"


def build_resistance_pathway_consensus(gsea_tables: dict[str, pd.DataFrame]) -> pd.DataFrame:
    frames = []
    for label in RESISTANCE_DATASETS:
        df = gsea_tables[label].copy()
        df["dataset"] = label
        frames.append(df)
    combined = pd.concat(frames, ignore_index=True)

    rows = []
    for (collection, pathway), grp in combined.groupby(["collection", "pathway"], observed=True):
        n_tested = len(grp)
        n_fdr05 = int((grp["fdr"] < 0.05).sum())
        n_nominal = int((grp["nom_pvalue"] < 0.05).sum())
        n_pos = int((grp["NES"] > 0).sum())
        n_neg = int((grp["NES"] < 0).sum())
        if n_pos > 0 and n_neg > 0:
            direction_consistency = "mixed"
        elif n_pos > 0:
            direction_consistency = "all_positive"
        else:
            direction_consistency = "all_negative"

        lead_sets = [set(g.split(";")) for g in grp["leading_edge_genes"].dropna()]
        leading_edge_overlap = len(set.intersection(*lead_sets)) if len(lead_sets) >= 2 else np.nan

        rows.append(
            {
                "collection": collection,
                "pathway": pathway,
                "datasets_tested": n_tested,
                "datasets_tested_names": ",".join(sorted(grp["dataset"])),
                "datasets_FDR05": n_fdr05,
                "datasets_nominal": n_nominal,
                "NES_positive_count": n_pos,
                "NES_negative_count": n_neg,
                "direction_consistency": direction_consistency,
                "median_NES": float(grp["NES"].median()),
                "mean_NES": float(grp["NES"].mean()),
                "leading_edge_overlap": leading_edge_overlap,
            }
        )
    out = pd.DataFrame(rows)
    out["consensus_category"] = out.apply(_classify, axis=1)
    out = out.sort_values(["consensus_category", "datasets_FDR05", "median_NES"], ascending=[True, False, False])

    counts = out["consensus_category"].value_counts().to_dict()
    logger.info("build_resistance_pathway_consensus: %d pathway x collection rows, categories=%s", len(out), counts)
    return out


def run_consensus(config_path: str | Path = "config/config.yaml") -> pd.DataFrame:
    config = _load_config(config_path)
    cfg = config["systems_network"]
    tables_dir = Path(cfg["output"]["tables_dir"])

    gsea_tables = {label: pd.read_csv(tables_dir / f"gsea_{label}.tsv", sep="\t") for label in RESISTANCE_DATASETS}
    consensus = build_resistance_pathway_consensus(gsea_tables)
    consensus.to_csv(tables_dir / "resistance_pathway_consensus.tsv", sep="\t", index=False)
    return consensus


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run_consensus()
