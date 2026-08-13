"""Systems-network phase 12: multimodal pathway convergence.

One row per pathway (collection, pathway) across the full union tested in
any of the five GSEA runs (three resistance RNA datasets, GSE245601 acute,
CRISPR). Combines: RNA resistance consensus category (Phase 5, GSE245601
NEVER included in that consensus), acute GSE245601 response (shown
separately), CRISPR pathway-level functional enrichment (Phase 11, same
preranked-GSEA methodology as the RNA datasets so NES is directly
comparable), and candidate pathway association (Phase 7).

Category assignment (fixed precedence, declared before inspecting results):
  1. MULTIMODAL_PATHWAY: reproducible RNA resistance signal
     (STRONG_CONSENSUS or DIRECTIONAL_CONSENSUS) AND CRISPR FDR<0.05
  2. RESISTANCE_ONLY: reproducible RNA resistance signal, CRISPR FDR>=0.05
  3. MIXED: RNA resistance direction itself is MIXED (Phase 5) and at
     least one of CRISPR/acute is also significant
  4. FUNCTIONAL_ONLY: CRISPR FDR<0.05, no reproducible AND no mixed RNA
     resistance signal (i.e. the pathway has no resistance-dataset signal
     at all, not merely a contradictory one)
  5. ACUTE_RESPONSE: no RNA resistance signal, no CRISPR signal, but
     GSE245601 acute FDR<0.05
  6. NONE: none of the above (excluded from the "top" summaries, kept in
     the full table for transparency)

Codex review fix: an earlier version checked FUNCTIONAL_ONLY (old #3)
before MIXED (old #5), so every one of the 134 pathways with a genuinely
contradictory (MIXED) RNA direction AND CRISPR FDR<0.05 was mislabeled
FUNCTIONAL_ONLY -- indistinguishable from a pathway with no RNA signal at
all. MIXED is now checked first among the two, matching this docstring's
intent (and the module's own name for the category).
"""

from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd
import yaml

logger = logging.getLogger(__name__)

CANDIDATES = ["USP34", "VEZF1", "EML5", "CITED2"]


def _load_config(config_path: str | Path) -> dict:
    with open(config_path) as f:
        return yaml.safe_load(f)


def _classify(row: pd.Series) -> str:
    has_resistance = row["resistance_consensus_class"] in ("STRONG_CONSENSUS", "DIRECTIONAL_CONSENSUS")
    has_crispr = row["crispr_fdr"] < 0.05 if pd.notna(row["crispr_fdr"]) else False
    has_acute = row["gse245601_fdr"] < 0.05 if pd.notna(row["gse245601_fdr"]) else False
    is_mixed_resistance = row["resistance_consensus_class"] == "MIXED"

    if has_resistance and has_crispr:
        return "MULTIMODAL_PATHWAY"
    if has_resistance:
        return "RESISTANCE_ONLY"
    if is_mixed_resistance and (has_crispr or has_acute):
        return "MIXED"
    if has_crispr:
        return "FUNCTIONAL_ONLY"
    if has_acute:
        return "ACUTE_RESPONSE"
    return "NONE"


def build_multimodal_convergence(
    consensus: pd.DataFrame,
    gsea_crispr: pd.DataFrame,
    gsea_acute: pd.DataFrame,
    candidate_membership: pd.DataFrame,
) -> pd.DataFrame:
    consensus_idx = consensus.set_index(["collection", "pathway"])
    crispr_idx = gsea_crispr.set_index(["collection", "pathway"])
    acute_idx = gsea_acute.set_index(["collection", "pathway"])

    all_keys = set(consensus_idx.index) | set(crispr_idx.index) | set(acute_idx.index)
    logger.info("build_multimodal_convergence: %d unique (collection,pathway) keys across all sources", len(all_keys))

    candidate_by_pathway: dict[tuple[str, str], set[str]] = {}
    connected = candidate_membership.loc[
        candidate_membership["candidate_is_member"] | (candidate_membership["candidate_is_leading_edge_datasets"].fillna("") != "")
    ]
    for _, row in connected.iterrows():
        candidate_by_pathway.setdefault((row["collection"], row["pathway"]), set()).add(row["candidate"])

    rows = []
    for key in all_keys:
        collection, pathway = key
        rc = consensus_idx.loc[key] if key in consensus_idx.index else None
        cr = crispr_idx.loc[key] if key in crispr_idx.index else None
        ac = acute_idx.loc[key] if key in acute_idx.index else None

        row = {
            "collection": collection,
            "pathway": pathway,
            "resistance_consensus_class": rc["consensus_category"] if rc is not None else "not_tested_in_resistance",
            "resistance_median_NES": float(rc["median_NES"]) if rc is not None else float("nan"),
            "resistance_datasets_FDR05": int(rc["datasets_FDR05"]) if rc is not None else 0,
            "gse245601_NES": float(ac["NES"]) if ac is not None else float("nan"),
            "gse245601_fdr": float(ac["fdr"]) if ac is not None else float("nan"),
            "crispr_NES": float(cr["NES"]) if cr is not None else float("nan"),
            "crispr_fdr": float(cr["fdr"]) if cr is not None else float("nan"),
            "candidates_connected": ",".join(sorted(candidate_by_pathway.get(key, set()))),
        }
        rows.append(row)

    out = pd.DataFrame(rows)
    out["convergence_category"] = out.apply(_classify, axis=1)
    out = out.sort_values(["convergence_category", "resistance_datasets_FDR05"], ascending=[True, False])

    counts = out["convergence_category"].value_counts().to_dict()
    logger.info("build_multimodal_convergence: category counts=%s", counts)
    return out


def run_multimodal(config_path: str | Path = "config/config.yaml") -> pd.DataFrame:
    config = _load_config(config_path)
    cfg = config["systems_network"]
    tables_dir = Path(cfg["output"]["tables_dir"])

    consensus = pd.read_csv(tables_dir / "resistance_pathway_consensus.tsv", sep="\t")
    gsea_crispr = pd.read_csv(tables_dir / "gsea_crispr.tsv", sep="\t")
    gsea_acute = pd.read_csv(tables_dir / "gsea_gse245601.tsv", sep="\t")
    candidate_membership = pd.read_csv(tables_dir / "candidate_pathway_membership.tsv", sep="\t")

    out = build_multimodal_convergence(consensus, gsea_crispr, gsea_acute, candidate_membership)
    out.to_csv(tables_dir / "multimodal_pathway_convergence.tsv", sep="\t", index=False)
    return out


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run_multimodal()
