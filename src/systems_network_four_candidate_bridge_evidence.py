"""Comparative systems-network audit, Part 5: evidence audit of every
UNIQUE bridge/connector gene identified for the four frozen candidates in
Parts 2 and 4 -- reusing the completed USP34 audit
(results/tables/systems_network/USP34_bridge_gene_evidence.tsv) verbatim
and extending the same method to CITED2's and VEZF1's connector genes.

"Bridge/connector gene" here means: the gene(s) mediating a candidate's
network story, whether that gene sits on a 2-hop shortest path (USP34's
USP9X/RPS27A/UBC/UBB) or is itself the candidate's direct (1-hop) STRING
partner that carries the candidate into resistance-relevant pathway biology
(CITED2's EP300/CREBBP/TFAP2C/HIF1A/TP53; VEZF1's sole connection, DMTN).
EML5 has no network neighbors at all (Part 2) and therefore has no bridge
gene to audit -- reported as NOT_ASSESSABLE, not fabricated.

Reads only already-frozen outputs -- does not rerun any upstream
ranking/GSEA/consensus phase:
  - results/tables/systems_network/USP34_bridge_gene_evidence.tsv (reused
    verbatim for USP34's 4 genes)
  - results/networks/systems_network/cytoscape/network_nodes.tsv (frozen
    five-layer per-gene display, same source used for USP34)
  - results/tables/systems_network/{gse118713,gse240112,gse111151,
    gse245601}_ranked_genes.tsv (frozen Phase 2 per-dataset p-values)
  - data/processed/labels.parquet (frozen CRISPR label source)

CRISPR sign convention: negative effect_size = sensitising_KO, positive =
tolerance_associated_KO. GSE245601 is acute (12h) response only and is
reported as a fifth layer for completeness but never counted toward the
classification below (same rule as the USP34 audit).

Classification (identical rule to the USP34 audit, applied uniformly):
  A. DATA_SUPPORTED_BRIDGE      -- FDR<0.05 in CRISPR and/or a resistance
                                    dataset.
  B. PARTIAL_SUPPORT            -- no FDR<0.05 hit, but >=1 nominal
                                    (p<0.05) hit in CRISPR or a resistance
                                    dataset (a GSE245601-only nominal hit
                                    does not qualify).
  C. NETWORK_ONLY_GENERIC_BRIDGE -- no FDR<0.05 and no qualifying nominal
                                    hit anywhere in CRISPR/resistance data.
  D. NOT_ASSESSABLE             -- gene absent/unmeasured, or (EML5) no
                                    bridge gene exists to assess.
"""

from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd

logger = logging.getLogger(__name__)

RESISTANCE_DATASETS = ["gse118713", "gse240112", "gse111151"]
NOMINAL_P = 0.05
FDR_SIG = 0.05

NODES_PATH = Path("results/networks/systems_network/cytoscape/network_nodes.tsv")
RANKED_PATHS = {
    "gse118713": Path("results/tables/systems_network/gse118713_ranked_genes.tsv"),
    "gse240112": Path("results/tables/systems_network/gse240112_ranked_genes.tsv"),
    "gse111151": Path("results/tables/systems_network/gse111151_ranked_genes.tsv"),
    "gse245601": Path("results/tables/systems_network/gse245601_ranked_genes.tsv"),
}
USP34_BRIDGE_EVIDENCE_PATH = Path("results/tables/systems_network/USP34_bridge_gene_evidence.tsv")
OUT_TABLE = Path("results/tables/systems_network/four_candidate_bridge_evidence.tsv")

# (candidate, bridge_gene, what the bridge connects the candidate to,
#  geometric relationship -- "1-hop_direct_partner" vs "2-hop_intermediate")
NEW_BRIDGES = [
    ("CITED2", "EP300", "REACTOME_SIGNALING_BY_RECEPTOR_TYROSINE_KINASES (STRONG_CONSENSUS leading-edge, gse240112)", "1-hop_direct_partner"),
    ("CITED2", "CREBBP", "no STRONG_CONSENSUS Hallmark/Reactome leading-edge pathway found for CREBBP in any resistance dataset", "1-hop_direct_partner"),
    ("CITED2", "TFAP2C", "HALLMARK_ESTROGEN_RESPONSE_EARLY/LATE (STRONG_CONSENSUS leading-edge, gse118713+gse111151)", "1-hop_direct_partner"),
    ("CITED2", "HIF1A", "REACTOME_SIGNALING_BY_RECEPTOR_TYROSINE_KINASES + HALLMARK_INFLAMMATORY_RESPONSE (STRONG_CONSENSUS leading-edge, gse111151/gse240112)", "1-hop_direct_partner"),
    ("CITED2", "TP53", "HALLMARK_E2F_TARGETS + REACTOME_CELL_CYCLE_MITOTIC/CHECKPOINTS (STRONG_CONSENSUS leading-edge, all 3 resistance datasets)", "1-hop_direct_partner"),
    ("VEZF1", "DMTN", "HALLMARK_HEME_METABOLISM (STRONG_CONSENSUS; both VEZF1 and DMTN are independently leading-edge)", "1-hop_direct_partner_pathway_co_membership"),
]


def _direction(log2fc: float) -> str:
    return "up" if log2fc > 0 else ("down" if log2fc < 0 else "flat")


def build_new_bridge_rows(nodes: pd.DataFrame, ranked: dict[str, pd.DataFrame]) -> pd.DataFrame:
    nodes_idx = nodes.set_index("gene")
    crispr_labels = pd.read_parquet("data/processed/labels.parquet").set_index("gene")

    out_rows = []
    for candidate, gene, bridges_to, relationship in NEW_BRIDGES:
        n = nodes_idx.loc[gene]
        crispr_row = crispr_labels.loc[gene]
        crispr_p = float(crispr_row["p_value"])
        crispr_fdr = float(crispr_row["fdr"])
        crispr_effect = float(crispr_row["effect_size"])
        crispr_sig_fdr = crispr_fdr < FDR_SIG
        crispr_sig_nom = crispr_p < NOMINAL_P

        rec = {
            "candidate": candidate,
            "bridge_gene": gene,
            "bridge_relationship": relationship,
            "bridges_candidate_to": bridges_to,
            "crispr_effect": crispr_effect,
            "crispr_p": crispr_p,
            "crispr_fdr": crispr_fdr,
            "crispr_direction": n["crispr_direction"],
            "crispr_significant_fdr05": crispr_sig_fdr,
            "crispr_nominal_p05": crispr_sig_nom,
        }

        any_fdr_sig = crispr_sig_fdr
        any_nominal_resistance_sig = False
        directions = []

        for dataset in RESISTANCE_DATASETS:
            r = ranked[dataset].set_index("gene").loc[gene]
            log2fc = float(r["log2fc"])
            p = float(r["p_value"])
            fdr = float(r["fdr"])
            sig_fdr = fdr < FDR_SIG
            sig_nom = p < NOMINAL_P
            direction = _direction(log2fc)
            rec[f"{dataset}_log2fc"] = log2fc
            rec[f"{dataset}_p"] = p
            rec[f"{dataset}_fdr"] = fdr
            rec[f"{dataset}_direction"] = direction
            rec[f"{dataset}_significant_fdr05"] = sig_fdr
            rec[f"{dataset}_nominal_p05"] = sig_nom
            directions.append(direction)
            any_fdr_sig = any_fdr_sig or sig_fdr
            any_nominal_resistance_sig = any_nominal_resistance_sig or sig_nom

        acute = ranked["gse245601"].set_index("gene").loc[gene]
        acute_log2fc = float(acute["log2fc"])
        acute_p = float(acute["p_value"])
        acute_fdr = float(acute["fdr"])
        rec["gse245601_track_a_log2fc"] = acute_log2fc
        rec["gse245601_track_a_p"] = acute_p
        rec["gse245601_track_a_fdr"] = acute_fdr
        rec["gse245601_track_a_direction"] = _direction(acute_log2fc)
        rec["gse245601_track_a_significant_fdr05"] = acute_fdr < FDR_SIG
        rec["gse245601_track_a_nominal_p05"] = acute_p < NOMINAL_P
        rec["gse245601_acute_note"] = "ACUTE 12h response -- not resistance evidence, not counted in classification"

        n_up = directions.count("up")
        n_down = directions.count("down")
        rec["resistance_direction_consistency"] = "consistent" if (n_up == 0 or n_down == 0) else "mixed"
        rec["resistance_pattern_frozen"] = n["resistance_pattern"]

        if any_fdr_sig:
            classification = "A_DATA_SUPPORTED_BRIDGE"
        elif crispr_sig_nom or any_nominal_resistance_sig:
            classification = "B_PARTIAL_SUPPORT"
        else:
            classification = "C_NETWORK_ONLY_GENERIC_BRIDGE"
        rec["classification"] = classification

        out_rows.append(rec)

    return pd.DataFrame(out_rows)


def build_four_candidate_bridge_evidence(usp34_bridge_evidence: pd.DataFrame, nodes: pd.DataFrame, ranked: dict[str, pd.DataFrame]) -> pd.DataFrame:
    usp34 = usp34_bridge_evidence.copy()
    usp34.insert(0, "candidate", "USP34")
    usp34 = usp34.rename(columns={"gene": "bridge_gene", "bridges_usp34_to": "bridges_candidate_to"})
    usp34.insert(2, "bridge_relationship", "2-hop_intermediate")
    # column order must match new_rows; reindex after concat instead of here

    new_rows = build_new_bridge_rows(nodes, ranked)

    eml5_row = {
        "candidate": "EML5",
        "bridge_gene": "NONE",
        "bridge_relationship": "not_applicable",
        "bridges_candidate_to": "EML5 has zero network neighbors (Part 2); no bridge gene exists to audit",
        "classification": "D_NOT_ASSESSABLE",
    }

    combined = pd.concat([usp34, new_rows, pd.DataFrame([eml5_row])], ignore_index=True)

    lead_cols = ["candidate", "bridge_gene", "bridge_relationship", "bridges_candidate_to", "classification"]
    other_cols = [c for c in combined.columns if c not in lead_cols]
    combined = combined[lead_cols + other_cols]
    return combined


def run(
    usp34_bridge_evidence_path: Path = USP34_BRIDGE_EVIDENCE_PATH,
    nodes_path: Path = NODES_PATH,
    ranked_paths: dict[str, Path] = None,
    out_table: Path = OUT_TABLE,
) -> pd.DataFrame:
    ranked_paths = ranked_paths or RANKED_PATHS
    usp34_bridge_evidence = pd.read_csv(usp34_bridge_evidence_path, sep="\t")
    nodes = pd.read_csv(nodes_path, sep="\t")
    ranked = {k: pd.read_csv(v, sep="\t") for k, v in ranked_paths.items()}

    out = build_four_candidate_bridge_evidence(usp34_bridge_evidence, nodes, ranked)
    out_table.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(out_table, sep="\t", index=False)
    logger.info(
        "wrote %s (%d rows; classification counts: %s)",
        out_table,
        len(out),
        out["classification"].value_counts().to_dict(),
    )
    return out


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run()
