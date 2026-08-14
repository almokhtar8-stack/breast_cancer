"""Comparative systems-network audit, Part 8: final head-to-head comparison
of the four frozen candidates.

Every numeric/factual column is pulled programmatically from already-frozen
or already-computed-in-this-audit tables (never hand-typed):
  - results/networks/systems_network/cytoscape/network_nodes.tsv (frozen
    five-layer per-candidate evidence)
  - results/tables/systems_network/four_candidate_direct_neighbors.tsv
    (Part 2, this audit)
  - results/tables/systems_network/four_candidate_bridge_evidence.tsv
    (Part 5, this audit)
  - results/tables/systems_network/four_candidate_convergence.tsv (Part 7,
    this audit)
  - results/networks/systems_network/network_node_metrics.tsv (degree, for
    the generic-hub concern column)

The four qualitative synthesis columns (mechanistic_specificity,
generic_hub_concern, candidate_convergence summary text, and the final
systems-mechanism classification) require cross-referencing all of the
above at once and are not reducible to a single groupby -- they are set
explicitly per candidate below, with the specific frozen numbers that
justify each call documented inline (same treatment as any other
judgment-call table already produced in this phase, e.g.
resistance_pattern in the node table). No candidate ranking/order implies a
change to the frozen therapeutic shortlist.
"""

from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd

logger = logging.getLogger(__name__)

CANDIDATES = ["USP34", "VEZF1", "EML5", "CITED2"]

NODES_PATH = Path("results/networks/systems_network/cytoscape/network_nodes.tsv")
NEIGHBORS_PATH = Path("results/tables/systems_network/four_candidate_direct_neighbors.tsv")
BRIDGE_EVIDENCE_PATH = Path("results/tables/systems_network/four_candidate_bridge_evidence.tsv")
CONVERGENCE_PATH = Path("results/tables/systems_network/four_candidate_convergence.tsv")
METRICS_PATH = Path("results/networks/systems_network/network_node_metrics.tsv")
OUT_TABLE = Path("results/tables/systems_network/four_candidate_network_audit.tsv")

# Qualitative synthesis calls -- see module docstring. Each rationale cites
# the specific frozen/computed number it rests on so it can be checked
# against the source tables rather than taken on faith.
SYNTHESIS = {
    "USP34": {
        "strongest_pathway_module": "GOBP_CANONICAL_WNT_SIGNALING_PATHWAY / GOBP_WNT_SIGNALING_PATHWAY (STRONG_CONSENSUS; USP34 itself is leading-edge in gse118713+gse111151)",
        "mechanistic_specificity": "MIXED",
        "mechanistic_specificity_rationale": "Direct curated WNT pathway membership is specific and candidate-driven (USP34 itself is leading-edge, not just an interactor). But USP34's only network route to CTNNB1/PTEN/EP300/SOX2 runs through RPS27A/UBB/UBC -- all high-degree (15/17/20), high-betweenness generic ubiquitin hubs per network_node_metrics.tsv -- with no independent CRISPR/RNA support (Part 5: C tier). USP9X is lower-degree (6) but only B-tier (one contradicted nominal hit).",
        "generic_hub_concern": "HIGH (network-bridge layer only; direct pathway layer is not hub-driven)",
        "systems_mechanism_classification": "2_MODERATE_SYSTEMS_SUPPORT",
        "classification_rationale": "Own-candidate CRISPR FDR=0.042 (sensitising_KO) + own-candidate RNA FDR=0.0073 (gse118713, up) + direct STRONG_CONSENSUS leading-edge WNT pathway membership are real, multi-layer, candidate-level support. But the broader 2-hop network story to CTNNB1/PTEN/EP300/SOX2 is not independently validated (best bridge is B-tier) and is generic-hub-mediated -- not enough to call STRONG.",
    },
    "VEZF1": {
        "strongest_pathway_module": "GOBP_BLOOD_VESSEL_MORPHOGENESIS / HALLMARK_HEME_METABOLISM (both STRONG_CONSENSUS; VEZF1 itself is leading-edge in gse240112)",
        "mechanistic_specificity": "SPECIFIC / INTERPRETABLE",
        "mechanistic_specificity_rationale": "VEZF1's only network connection (DMTN, degree=1) is not a generic hub by any measure, and DMTN is itself very strongly independently resistance-associated (gse118713 FDR=0.00024, gse240112 FDR=0.00021, both up). Every piece of VEZF1's story is candidate/pathway-specific; there is no hub-mediated component to discount.",
        "generic_hub_concern": "NONE",
        "systems_mechanism_classification": "2_MODERATE_SYSTEMS_SUPPORT",
        "classification_rationale": "Own-candidate CRISPR FDR=0.037 (sensitising_KO) + own-candidate RNA FDR=0.0195 (gse240112, up) + 2 direct STRONG_CONSENSUS leading-edge pathways + its sole network connector is itself A-tier data-supported. Held at MODERATE rather than STRONG only because network breadth is minimal (1 edge total) and own-candidate RNA significance is limited to 1 of 3 resistance datasets -- the evidence that exists is high-quality but thin in volume.",
    },
    "EML5": {
        "strongest_pathway_module": "NO RESOLVED CONNECTION -- zero pathway rows, zero STRING partners even down to STRING's lowest confidence band (docs/SYSTEMS_NETWORK_NODE_RULE.md)",
        "mechanistic_specificity": "NO RESOLVED NETWORK",
        "mechanistic_specificity_rationale": "EML5 has zero edges in the frozen network and zero rows in candidate_pathway_membership.tsv. There is nothing to assess for hub-driven vs specific connectivity because there is no connectivity.",
        "generic_hub_concern": "NOT APPLICABLE (no network)",
        "systems_mechanism_classification": "4_DATA_SUPPORTED_BUT_MECHANISTICALLY_UNRESOLVED",
        "classification_rationale": "EML5 has the strongest own-candidate RESISTANCE RNA signal of all four candidates (gse118713 FDR=0.000129, gse240112 FDR=0.040, both up -- 2 of 3 resistance datasets significant, consistent direction) but zero systems/network/pathway footprint. This is a textbook data-supported-but-mechanistically-unresolved case: the expression evidence is real, but this analysis layer has nothing further to say about mechanism.",
    },
    "CITED2": {
        "strongest_pathway_module": "HALLMARK_UV_RESPONSE_DN (direct curated membership + CITED2 itself leading-edge in gse240112 + STRONG_CONSENSUS RNA + CRISPR pathway-level FDR=0.0088 -- CITED2's only genuinely MULTIMODAL_PATHWAY connection)",
        "mechanistic_specificity": "MIXED",
        "mechanistic_specificity_rationale": "TFAP2C (degree=3) and HIF1A (degree=17) are low/moderate-degree, independently data-supported (A-tier) connectors -- specific signal. But 2 of CITED2's 5 audited connectors, EP300 (degree=28) and TP53 (degree=34), are respectively the #3 and #1 highest-degree genes in the entire 119-node network (network_node_metrics.tsv) and are both flagged is_bridge_gene_high_betweenness=True -- their prominence in CITED2's network picture partly reflects generic network centrality, not necessarily CITED2-specific biology, even though both do independently reach FDR<0.05 in a resistance dataset.",
        "generic_hub_concern": "MODERATE (2 of 5 audited connectors are the network's top-2 degree hubs)",
        "systems_mechanism_classification": "2_MODERATE_SYSTEMS_SUPPORT",
        "classification_rationale": "Own-candidate RNA FDR=0.0087 (gse240112, up); own-candidate CRISPR not significant (FDR=0.110). Broadest direct network neighborhood of the four (18 neighbors) and 4 of 5 audited connectors reach A-tier independent support (Part 5), including a dual CRISPR+RNA-significant, low-degree, specific one (TFAP2C, CRISPR FDR=0.048, gse118713 FDR=0.001). Held at MODERATE rather than STRONG because the two most connectivity-heavy routes into headline pathways (E2F_TARGETS, ESTROGEN_DEPENDENT_GENE_EXPRESSION) run through the network's biggest generic hubs (EP300, TP53) and only one pathway (UV_RESPONSE_DN) is a genuinely direct, multimodal, candidate-driven connection.",
    },
}


def build_head_to_head(
    nodes: pd.DataFrame,
    neighbors: pd.DataFrame,
    bridge_evidence: pd.DataFrame,
    convergence: pd.DataFrame,
) -> pd.DataFrame:
    nodes_idx = nodes.set_index("gene")

    rows = []
    for candidate in CANDIDATES:
        n = nodes_idx.loc[candidate]
        cand_neighbors = neighbors.loc[neighbors["candidate"] == candidate]
        n_direct = 0 if (len(cand_neighbors) == 1 and cand_neighbors["neighbor_gene"].iloc[0] == "NO_RESOLVED_NETWORK_NEIGHBOURHOOD") else len(cand_neighbors)

        strongest_edges = cand_neighbors.sort_values("confidence", ascending=False).head(3)
        strongest_interactions = "; ".join(f"{r['neighbor_gene']} ({r['interaction_type']}, {r['confidence']:.3f})" for _, r in strongest_edges.iterrows() if pd.notna(r["confidence"])) or "none"

        cand_bridges = bridge_evidence.loc[bridge_evidence["candidate"] == candidate]
        cand_bridges_valid = cand_bridges.loc[cand_bridges["bridge_gene"] != "NONE"]
        tier_rank = {"A_DATA_SUPPORTED_BRIDGE": 0, "B_PARTIAL_SUPPORT": 1, "C_NETWORK_ONLY_GENERIC_BRIDGE": 2, "D_NOT_ASSESSABLE": 3}
        if len(cand_bridges_valid) == 0:
            strongest_bridge, bridge_supported = "none (no network neighbors)", "not_applicable"
        else:
            best = cand_bridges_valid.assign(_rank=cand_bridges_valid["classification"].map(tier_rank)).sort_values("_rank").iloc[0]
            strongest_bridge = f"{best['bridge_gene']} ({best['classification']})"
            bridge_supported = {"A_DATA_SUPPORTED_BRIDGE": "yes", "B_PARTIAL_SUPPORT": "partial", "C_NETWORK_ONLY_GENERIC_BRIDGE": "no"}[best["classification"]]

        conv = convergence.loc[(convergence["candidate_A"] == candidate) | (convergence["candidate_B"] == candidate)]
        conv_partners = conv.loc[conv["any_convergence"], ["candidate_A", "candidate_B"]].values.flatten()
        conv_partners = sorted(set(conv_partners) - {candidate})
        convergence_summary = f"converges with {', '.join(conv_partners)}" if conv_partners else "no convergence with any other candidate"

        rec = {
            "candidate": candidate,
            "frozen_crispr_direction": n["crispr_direction"],
            "frozen_crispr_fdr": n["crispr_fdr"],
            "resistance_rna_datasets_fdr05": sum(n[f"{d}_fdr"] < 0.05 for d in ["gse118713", "gse240112", "gse111151"]),
            "resistance_rna_direction_summary": ",".join(f"{d}:{'up' if n[f'{d}_log2fc'] > 0 else 'down'}{'*' if n[f'{d}_fdr'] < 0.05 else ''}" for d in ["gse118713", "gse240112", "gse111151"]),
            "acute_gse245601_fdr": n["gse245601_acute_fdr"],
            "acute_gse245601_significant": bool(n["gse245601_acute_fdr"] < 0.05),
            "n_direct_neighbors": n_direct,
            "strongest_direct_interactions": strongest_interactions,
            "strongest_pathway_module": SYNTHESIS[candidate]["strongest_pathway_module"],
            "strongest_bridge": strongest_bridge,
            "bridge_independently_supported": bridge_supported,
            "generic_hub_concern": SYNTHESIS[candidate]["generic_hub_concern"],
            "mechanistic_specificity": SYNTHESIS[candidate]["mechanistic_specificity"],
            "mechanistic_specificity_rationale": SYNTHESIS[candidate]["mechanistic_specificity_rationale"],
            "candidate_convergence": convergence_summary,
            "systems_mechanism_classification": SYNTHESIS[candidate]["systems_mechanism_classification"],
            "classification_rationale": SYNTHESIS[candidate]["classification_rationale"],
        }
        rows.append(rec)

    return pd.DataFrame(rows)


def run(
    nodes_path: Path = NODES_PATH,
    neighbors_path: Path = NEIGHBORS_PATH,
    bridge_evidence_path: Path = BRIDGE_EVIDENCE_PATH,
    convergence_path: Path = CONVERGENCE_PATH,
    out_table: Path = OUT_TABLE,
) -> pd.DataFrame:
    nodes = pd.read_csv(nodes_path, sep="\t")
    neighbors = pd.read_csv(neighbors_path, sep="\t")
    bridge_evidence = pd.read_csv(bridge_evidence_path, sep="\t")
    convergence = pd.read_csv(convergence_path, sep="\t")

    out = build_head_to_head(nodes, neighbors, bridge_evidence, convergence)
    out_table.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(out_table, sep="\t", index=False)
    logger.info("wrote %s (%d rows)", out_table, len(out))
    return out


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run()
