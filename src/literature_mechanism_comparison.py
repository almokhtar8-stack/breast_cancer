"""Literature/mechanism review, Part 9: head-to-head literature comparison
across the four frozen candidates.

Reads only results/tables/literature_mechanism/four_candidate_claim_evidence.tsv
(built by src/literature_mechanism_build_tables.py from the verified
literature search). Boolean/count columns are computed directly from that
table; the four qualitative synthesis columns (overall_literature_depth,
fit_with_our_data, mechanism_confidence, recommended_next_step) require
judgment across all of a candidate's claims at once and are set explicitly
below, with the rationale documented inline, in the same manner as the
prior systems-network audit's synthesis columns.
"""

from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd

logger = logging.getLogger(__name__)

CANDIDATES = ["USP34", "VEZF1", "EML5", "CITED2"]

CLAIMS_PATH = Path("results/tables/literature_mechanism/four_candidate_claim_evidence.tsv")
OUT_TABLE = Path("results/tables/literature_mechanism/four_candidate_literature_comparison.tsv")

SYNTHESIS = {
    "USP34": {
        "candidate_partner_evidence": "USP9X-SOX2-CTNNB1 not tested together; USP34-SOX2 (laryngeal SCC, PMID 32783291) and USP34-c-Myc (HCC, PMID 40260316) are independently documented partner relationships in non-breast cancers, offering only weak, indirect precedent for the project's STRING-hop hypothesis.",
        "overall_literature_depth": "MODERATE (8 claims; one well-characterized primary mechanism (Axin/Wnt), several other pathways documented, but pleiotropic/scattered across cancer types, none breast+endocrine)",
        "fit_with_our_data": "PARTIAL -- the Axin/Wnt mechanism (positive Wnt regulation) is directionally consistent with our own USP34-sensitising_KO + WNT-pathway-membership finding, but was never tested in breast tissue and is directly contradicted in direction by a separate mammary-epithelium paper (PMID 28499884).",
        "mechanism_confidence": "LOW-MODERATE",
        "recommended_next_step": "Literature review only for now: test whether project's own RNA data shows AXIN1/NKD1/TNFRSF19 (Lui 2011 Wnt targets) moving concordantly with USP34 in resistant samples, before further mechanism claims.",
    },
    "VEZF1": {
        "candidate_partner_evidence": "VEZF1 DIRECTLY represses CITED2 (ChIP-verified, PMID 29794136) -- a real, literature-documented regulatory link between two of the project's own four candidates, discovered independently of the project's network analysis (which found no direct VEZF1-CITED2 interaction). DMTN: no relationship found (explicit negative result).",
        "overall_literature_depth": "MODERATE-STRONG for vascular/developmental biology (11 claims, well-characterized founding gene); NONE for breast/endocrine biology specifically",
        "fit_with_our_data": "WEAK-PARTIAL -- VEZF1's angiogenic gene program is real and well-documented, but no paper connects it to tumor angiogenesis, hypoxia adaptation, or any drug resistance; the only 2 cancer papers (both HCC) act through unrelated, non-vascular mechanisms (PAQR4, TNS1/O-GlcNAcylation), undercutting the 'vascular program drives resistance' framing our own pathway enrichment suggests.",
        "mechanism_confidence": "LOW for a resistance-specific mechanism; MODERATE for VEZF1 being a genuine, non-generic angiogenic transcription factor in general.",
        "recommended_next_step": "The VEZF1-CITED2 literature link (independent of, and orthogonal to, our own network finding of no direct interaction) is worth flagging for the CITED2 track too; otherwise treat the vascular-resistance narrative as unresolved pending project-specific target-gene concordance checks (VEGFR2, endothelin-1, TIMP3/MMP2 in our own resistant-sample RNA data).",
    },
    "EML5": {
        "candidate_partner_evidence": "None -- EML5 has no documented protein partner from any primary experimental study; only weak, unvalidated STRING/BioGRID computational associations (MAD2L1, IFT140, TTC8, FBXW11).",
        "overall_literature_depth": "VERY THIN (5 claims, only 2 EML5-specific primary papers exist in all of PubMed/PMC, neither in cancer)",
        "fit_with_our_data": "UNRESOLVED -- literature independently corroborates that EML5 is an understudied gene; it does not explain the project's own strong, reproducible resistance-RNA signal.",
        "mechanism_confidence": "NONE",
        "recommended_next_step": "Not ready for mechanism-driven follow-up; if pursued at all, basic functional characterization (e.g. does EML5 even function as a canonical EML-family microtubule regulator) would need to precede any resistance-specific hypothesis.",
    },
    "CITED2": {
        "candidate_partner_evidence": "CITED2-EP300 (direct physical interaction, PMID 40313859, non-breast); CITED2-TFAP2C and CITED2-CREBBP/p300 (direct physical interaction via AP-2 dimerization domain, PMID 12586840, non-breast); CITED2-HIF1A (structural competition mechanism, PMID 12778114, non-breast); CITED2-FOXO3 (FOXO3a induces CITED2 transcription, demonstrated IN BREAST CANCER CELLS, PMID 18158893); CITED2-TP53 (CITED2 suppresses p53 accumulation, demonstrated IN BREAST CANCER, PMID 27627783). Every named partner in the project's network has at least one literature-documented mechanistic relationship with CITED2, though most are non-breast.",
        "overall_literature_depth": "STRONGEST OF THE FOUR (15 claims, including 1 LEVEL-1 direct tamoxifen paper, several breast-cancer mechanistic papers, and a genuine multi-tissue literature base)",
        "fit_with_our_data": "PARTIAL AND GENUINELY MIXED -- multiple breast-cancer papers support a pro-tumor/pro-resistance role (chemoresistance via p53 suppression, ER-activity enhancement reducing 4-OHT sensitivity, metastasis via IKKa/VEGFA), but the single most directly relevant paper (tamoxifen-resistance-selected cell model) is internally contradicted by its own clinical outcome data, and two non-breast cancer papers (colon, liver) show CITED2 acting as a tumor SUPPRESSOR in the opposite direction -- this is real literature depth with real, unresolved tension, not one-directional support.",
        "mechanism_confidence": "MODERATE, with an explicit caveat that the single directly-relevant (tamoxifen) study is self-contradictory and CITED2's role is documented as tissue-context-dependent (oncogenic in breast/prostate, tumor-suppressive in colon/liver).",
        "recommended_next_step": "Strongest literature case for continued mechanistic attention among the four, but the CITED2-tamoxifen-resistance direction question (does high CITED2 predict WORSE or BETTER tamoxifen response?) needs to be resolved against the project's own data before further claims are made.",
    },
}


def build_comparison_table(claims: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for candidate in CANDIDATES:
        c_all = claims.loc[claims["candidate"] == candidate]
        # exclude pathway-context claims (about a different gene) from every
        # rollup that would otherwise misrepresent them as candidate-specific
        c = c_all.loc[c_all["is_candidate_specific"]]
        c_named = c.loc[~c["paper_title"].isin(["(no paper found -- explicit negative result)", "(database resource, not a paper)", "(database annotation, not a paper)", "(database aggregation, not a paper)"])]

        direct_tamoxifen = c.loc[c["tamoxifen_specific"]]
        endocrine_resistance = c.loc[c["endocrine_resistance_specific"]]
        er_positive = c.loc[c["ER_positive_specific"]]
        breast_mechanistic = c.loc[c["breast_cancer_specific"] & (c["evidence_level"] <= 2)]
        pathway_evidence = c.loc[c["evidence_level"] == 3]
        other_cancers_only = c_named.loc[(c_named["breast_cancer_specific"] == False) & (c_named["evidence_level"].isin([2, 3]))]
        contradictory = c.loc[c["contradicts_project_hypothesis"]]

        pathway_context = c_all.loc[~c_all["is_candidate_specific"]]
        pathway_context_note = (
            "; ".join(f"{r.claim_id}: pathway-level evidence about a DIFFERENT gene ({r.paper_title}), not {candidate} itself -- excluded from all counts above" for r in pathway_context.itertuples())
            if len(pathway_context)
            else "none"
        )

        rows.append(
            {
                "candidate": candidate,
                "pathway_level_context_excluded_from_counts": pathway_context_note,
                "n_direct_tamoxifen_papers": len(direct_tamoxifen),
                "direct_tamoxifen_evidence": "; ".join(direct_tamoxifen["claim_id"]) or "none",
                "n_endocrine_resistance_papers": len(endocrine_resistance),
                "endocrine_resistance_evidence": "; ".join(endocrine_resistance["claim_id"]) or "none",
                "n_ER_positive_papers": len(er_positive),
                "ER_positive_evidence": "; ".join(er_positive["claim_id"]) or "none",
                "n_breast_cancer_mechanistic_papers": len(breast_mechanistic),
                "breast_cancer_mechanistic_evidence": "; ".join(breast_mechanistic["claim_id"]) or "none",
                "n_relevant_pathway_papers": len(pathway_evidence),
                "relevant_pathway_evidence": "; ".join(pathway_evidence["claim_id"]) or "none",
                "candidate_partner_evidence": SYNTHESIS[candidate]["candidate_partner_evidence"],
                "n_other_cancers_only_papers": len(other_cancers_only),
                "evidence_from_other_cancers_only": "; ".join(other_cancers_only["claim_id"]) or "none",
                "n_contradictory_papers": len(contradictory),
                "contradictory_evidence": "; ".join(contradictory["claim_id"]) or "none found (absence of contradiction is not proof)",
                "overall_literature_depth": SYNTHESIS[candidate]["overall_literature_depth"],
                "fit_with_our_data": SYNTHESIS[candidate]["fit_with_our_data"],
                "mechanism_confidence": SYNTHESIS[candidate]["mechanism_confidence"],
                "recommended_next_step": SYNTHESIS[candidate]["recommended_next_step"],
            }
        )
    return pd.DataFrame(rows)


def run(claims_path: Path = CLAIMS_PATH, out_table: Path = OUT_TABLE) -> pd.DataFrame:
    claims = pd.read_csv(claims_path, sep="\t")
    out = build_comparison_table(claims)
    out_table.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(out_table, sep="\t", index=False)
    logger.info("wrote %s (%d rows)", out_table, len(out))
    return out


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run()
