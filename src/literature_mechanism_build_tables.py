"""Literature/mechanism review phase: builds the claim-evidence table and
supporting outputs from a manually curated, citation-verified literature
survey of the four frozen therapeutic candidates (USP34, VEZF1, EML5,
CITED2).

Data source and provenance: every row below was gathered by four
independent literature-search passes (one per candidate) against PubMed,
PMC, NCBI eutils, Europe PMC, UniProt, and Human Protein Atlas, each paper
verified by directly fetching its PubMed/PMC/journal/database record
(title, authors, journal, year, PMID, DOI cross-checked against the
fetched record, not recalled from memory). Rows explicitly flagged
"narrative content not independently fetched" in `limitations` had their
bibliographic identifiers (PMID/DOI/title/authors/journal/year) confirmed
via NCBI eutils or Europe PMC/Semantic Scholar structured API records, but
the mechanistic description was drawn from secondary aggregation because
the publisher's full-text page returned a paywall/cookie-wall error on
direct fetch -- this is disclosed per-row, never silently smoothed over.
No PMID, DOI, author, journal, or finding in this module was invented.

This module performs no CRISPR/RNA-seq/network computation of any kind and
reads no upstream project data file -- it is pure literature curation, the
appropriate "source data" for a literature-review phase (analogous to how
config/config.yaml or PREANALYSIS.md hold curated, not computed, values
elsewhere in this project).
"""

from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd

logger = logging.getLogger(__name__)

OUT_TABLES_DIR = Path("results/tables/literature_mechanism")
OUT_REPORTS_DIR = Path("results/reports/literature_mechanism")

CLAIM_COLUMNS = [
    "candidate",
    "claim_id",
    "mechanistic_claim",
    "evidence_level",
    "direct_or_indirect",
    "breast_cancer_specific",
    "ER_positive_specific",
    "tamoxifen_specific",
    "endocrine_resistance_specific",
    "experimental_perturbation",
    "model_system",
    "direction_of_effect",
    "supports_project_hypothesis",
    "contradicts_project_hypothesis",
    "PMID",
    "DOI",
    "year",
    "journal",
    "paper_title",
    "exact_evidence_summary",
    "limitations",
    "review_or_primary",
    "fully_verified_primary_text",
    "is_candidate_specific",
]

# ---------------------------------------------------------------------------
# CITED2
# ---------------------------------------------------------------------------
CITED2_CLAIMS = [
    dict(
        claim_id="CITED2-01",
        mechanistic_claim="CITED2 upregulation was selected for during acquisition of tamoxifen resistance in a breast cancer cell model (retroviral insertional mutagenesis), but high CITED2 mRNA in patients associates with LONGER metastasis-free survival and GREATER clinical benefit from tamoxifen -- an internal contradiction the authors could not resolve.",
        evidence_level=1,
        direct_or_indirect="direct",
        breast_cancer_specific=True, ER_positive_specific=True, tamoxifen_specific=True, endocrine_resistance_specific=True,
        experimental_perturbation="retroviral insertional mutagenesis screen selecting resistant ZR-75-1 clones; correlative expression/outcome analysis in 2 independent patient cohorts",
        model_system="ZR-75-1 ER+ breast cancer cell line (69 tamoxifen-resistant derivatives); 620 lymph-node-negative primary tumors; 296 advanced-BC patients on 1st-line tamoxifen",
        direction_of_effect="contradictory within the paper: cell-model selection favors CITED2 up in resistance; clinical outcome data show CITED2 up = better MFS (HR=0.71,p=0.017) and better tamoxifen benefit (OR=1.91-2.20, p<=0.017)",
        supports_project_hypothesis=False, contradicts_project_hypothesis=False,  # genuinely mixed -- see summary
        PMID="19904269", DOI="10.1038/sj.bjc.6605423", year=2009, journal="British Journal of Cancer",
        paper_title="CITED2 and NCOR2 in anti-oestrogen resistance and progression of breast cancer",
        exact_evidence_summary="Cell lines with viral integration at the CITED2 locus showed strongly increased CITED2 mRNA. Multivariate MFS HR=0.71 (p=0.017, favorable). Tamoxifen-benefit OR=1.91 multivariate (p=0.017), OR=2.20 univariate (p=0.001), both favoring high CITED2. Authors: transfection experiments to directly confirm causality were 'not yet successful.'",
        limitations="Correlative clinical data cannot establish causation; cell-selection and clinical-outcome results point in opposite directions, unreconciled by authors.",
        review_or_primary="primary", fully_verified_primary_text=True,
    ),
    dict(
        claim_id="CITED2-02",
        mechanistic_claim="CITED2 overexpression increases ligand-independent ER transcriptional activity and reduces 4-hydroxytamoxifen growth-inhibitory sensitivity in ER+ breast cancer cell lines.",
        evidence_level=2, direct_or_indirect="direct",
        breast_cancer_specific=True, ER_positive_specific=True, tamoxifen_specific=True, endocrine_resistance_specific=False,
        experimental_perturbation="stable retroviral CITED2 overexpression vs empty vector; ERE-luciferase reporter; qRT-PCR (TFF1/PGR); ChIP; MTS proliferation +/- 1uM 4-OHT",
        model_system="MCF-7, T47D, CAMA-1 (ER+ human breast cancer lines)",
        direction_of_effect="CITED2 overexpression -> increased ER activity, estrogen-independent growth, reduced 4-OHT sensitivity",
        supports_project_hypothesis=True, contradicts_project_hypothesis=False,
        PMID="23811274", DOI="10.1016/j.bbrc.2013.06.063", year=2013, journal="Biochemical and Biophysical Research Communications",
        paper_title="CITED2 Modulates Estrogen Receptor Transcriptional Activity in Breast Cancer Cells",
        exact_evidence_summary="'reduced sensitivity to the growth inhibitory effects' of 4-OHT in CITED2-overexpressing cells; 2-4 fold increase in ER transcriptional activity without estrogen; CITED2 mRNA elevated in 3/8 ER+ tumors from patients with <5-year survival vs 0/11 with >5-year survival.",
        limitations="Small clinical subset (n=19), no formal statistics on the survival association; overexpression-only design (no knockdown arm); short-term drug-sensitivity assay, not an acquired-resistance model.",
        review_or_primary="primary", fully_verified_primary_text=True,
    ),
    dict(
        claim_id="CITED2-03",
        mechanistic_claim="CITED2 is an independent adverse prognostic factor in breast carcinoma and drives chemoresistance to epirubicin/5-FU by suppressing p53 accumulation; CITED2 is inversely correlated with ER status in this cohort.",
        evidence_level=2, direct_or_indirect="direct",
        breast_cancer_specific=True, ER_positive_specific=False, tamoxifen_specific=False, endocrine_resistance_specific=False,
        experimental_perturbation="CITED2 plasmid overexpression and siRNA knockdown; epirubicin/5-FU/paclitaxel viability assays; p53 Western blot after 5-FU",
        model_system="109-patient (+56 validation) IHC cohort; MCF-7 (ER+) and SKBR-3 (ER-/HER2+) cell lines",
        direction_of_effect="CITED2 up -> more proliferation, more chemoresistance to epirubicin/5-FU (not paclitaxel), less p53 accumulation; inversely correlated with ER status (p=0.033)",
        supports_project_hypothesis=True, contradicts_project_hypothesis=False,
        PMID="27627783", DOI="10.1111/cas.13081", year=2016, journal="Cancer Science",
        paper_title="CITED2 in breast carcinoma as a potent prognostic predictor associated with proliferation, migration and chemoresistance",
        exact_evidence_summary="CITED2 'significantly associated with increased incidence of recurrence and breast cancer-specific death,' independently prognostic for DFS (p=0.0036). 'CITED2 causes chemoresistance of breast carcinoma through inhibition of p53-dependent apoptosis.' Wound-healing assay showed CITED2 overexpression DECREASED migration -- authors flag this as paradoxical.",
        limitations="High-CITED2 tumors skew ER-negative/HER2+ in this cohort, limiting direct ER+/tamoxifen relevance; chemoresistance shown is cytotoxic-chemotherapy, not endocrine therapy; migration result internally paradoxical.",
        review_or_primary="primary", fully_verified_primary_text=True,
    ),
    dict(
        claim_id="CITED2-04",
        mechanistic_claim="CITED2 knockdown reduces breast cancer cell migration, invasion, and bone/brain metastasis via IKKa/NF-kB signaling, without affecting proliferation.",
        evidence_level=2, direct_or_indirect="direct",
        breast_cancer_specific=True, ER_positive_specific=False, tamoxifen_specific=False, endocrine_resistance_specific=False,
        experimental_perturbation="lentiviral shRNA CITED2 knockdown (>75%); transwell migration/invasion; intracardiac xenograft metastasis; ChIP on IKKa promoter; IKKa rescue",
        model_system="MDA-MB-231, MDA-MB-468 (ER-negative basal lines); athymic nude mice",
        direction_of_effect="CITED2 knockdown -> reduced migration/invasion/bone metastasis, no proliferation change; IKKa restoration rescues invasiveness",
        supports_project_hypothesis=True, contradicts_project_hypothesis=False,
        PMID="27216153", DOI="10.1158/1541-7786.MCR-16-0081", year=2016, journal="Molecular Cancer Research",
        paper_title="CITED2 Modulates Breast Cancer Metastatic Ability Through Effects on IKKa",
        exact_evidence_summary="shCITED2 cells: 'significantly reduced metastasis to bone' and 'similar growth rate to that of scramble cells.' ChIP: CITED2 localizes to the IKKa promoter.",
        limitations="ER-negative basal lines only -- cannot be extrapolated to ER+/luminal or endocrine-therapy biology without further evidence.",
        review_or_primary="primary", fully_verified_primary_text=True,
    ),
    dict(
        claim_id="CITED2-05",
        mechanistic_claim="CITED2 knockdown reduces primary breast tumor growth and vessel formation by reducing TGF-beta-induced VEGFA expression, in a Smad4-dependent (not universal) manner.",
        evidence_level=2, direct_or_indirect="direct",
        breast_cancer_specific=True, ER_positive_specific=False, tamoxifen_specific=False, endocrine_resistance_specific=False,
        experimental_perturbation="lentiviral shRNA CITED2 knockdown (>60%); orthotopic mammary fat pad xenograft; TGF-beta stimulation; ChIP on VEGFA promoter",
        model_system="MDA-MB-231 (Smad4-functional) and MDA-MB-468 (Smad4-null), both ER-negative; nude mice",
        direction_of_effect="CITED2 knockdown -> reduced tumor volume, vessel area, VEGFA121/165/189 -- effect only in the Smad4-functional line",
        supports_project_hypothesis=True, contradicts_project_hypothesis=False,
        PMID="28008154", DOI="10.18632/oncotarget.14048", year=2016, journal="Oncotarget",
        paper_title="Down-regulation of CITED2 attenuates breast tumor growth, vessel formation and TGF-beta-induced expression of VEGFA",
        exact_evidence_summary="Effect entirely dependent on functional Smad4 -- no vascularization/VEGFA change in the Smad4-null line, showing this is context-dependent, not universal.",
        limitations="ER-negative only; no endocrine-therapy angle tested; context-dependent (Smad4 status) effect.",
        review_or_primary="primary", fully_verified_primary_text=True,
    ),
    dict(
        claim_id="CITED2-06",
        mechanistic_claim="CITED2 physically bridges TFAP2 (AP-2) family transcription factors to the p300/CBP CH1 domain, binding the conserved AP-2 dimerization domain to enable coactivator recruitment.",
        evidence_level=3, direct_or_indirect="direct",
        breast_cancer_specific=False, ER_positive_specific=False, tamoxifen_specific=False, endocrine_resistance_specific=False,
        experimental_perturbation="co-immunoprecipitation, GST pulldown, transfection/reporter assays",
        model_system="U2-OS (human osteosarcoma) cells and in vitro biochemistry",
        direction_of_effect="CITED2 binds TFAP2A/B/C dimerization domain and bridges to p300/CBP CH1, enabling p300/CBP coactivation of TFAP2A-driven transcription",
        supports_project_hypothesis=True, contradicts_project_hypothesis=False,
        PMID="12586840", DOI="10.1074/jbc.M208144200", year=2003, journal="Journal of Biological Chemistry",
        paper_title="Physical and Functional Interactions among AP-2 Transcription Factors, p300/CREB-binding Protein, and CITED2",
        exact_evidence_summary="p300, TFAP2A, and endogenous CITED2 co-immunoprecipitate from transfected U2-OS cells, 'indicating that they can interact physically in vivo'; CITED2 interacts with the AP-2 dimerization domain.",
        limitations="Non-breast, non-cancer developmental biology system; full text paywalled -- mechanistic detail from secondary aggregation, bibliographic identity cross-verified via Semantic Scholar API.",
        review_or_primary="primary", fully_verified_primary_text=False,
    ),
    dict(
        claim_id="CITED2-07",
        mechanistic_claim="CITED2's transactivation domain structurally displaces HIF-1a from an overlapping binding site on the p300/CBP CH1 domain -- the structural basis for CITED2 as a HIF-1 negative-feedback regulator.",
        evidence_level=3, direct_or_indirect="direct",
        breast_cancer_specific=False, ER_positive_specific=False, tamoxifen_specific=False, endocrine_resistance_specific=False,
        experimental_perturbation="solution NMR structure determination of the CITED2-TAD/p300-CH1 complex; competition binding assays",
        model_system="purified recombinant protein domains, in vitro biophysics",
        direction_of_effect="CITED2 TAD folds on p300 CH1 and directly displaces HIF-1a C-TAD",
        supports_project_hypothesis=True, contradicts_project_hypothesis=False,
        PMID="12778114", DOI="10.1038/nsb936", year=2003, journal="Nature Structural Biology",
        paper_title="Structural basis for negative regulation of hypoxia-inducible factor-1alpha by CITED2",
        exact_evidence_summary="'the CITED2 TAD disrupts a complex of the HIF-1a C-TAD and the CH1 domain of p300/CBP by binding CH1 with high affinity... a conserved LPXL sequence motif in CITED2 and HIF-1a interacts with an overlapping binding site on CH1.'",
        limitations="Pure structural/biochemical study, no cellular or cancer context whatsoever.",
        review_or_primary="primary", fully_verified_primary_text=True,
    ),
    dict(
        claim_id="CITED2-08",
        mechanistic_claim="Under hypoxic stress, FOXO3a is activated and transcriptionally induces CITED2, which suppresses HIF1-induced pro-apoptotic target genes (NIX, RTP801), promoting survival -- demonstrated in breast cancer cells as well as fibroblasts.",
        evidence_level=2, direct_or_indirect="direct",
        breast_cancer_specific=True, ER_positive_specific=True, tamoxifen_specific=False, endocrine_resistance_specific=False,
        experimental_perturbation="FOXO3a knockdown/activation; CITED2 transcriptional induction measured; hypoxia exposure; apoptosis (NIX/RTP801) readouts",
        model_system="fibroblasts and breast cancer cells",
        direction_of_effect="FOXO3a activation -> CITED2 induction -> reduced HIF1-driven apoptosis (pro-survival)",
        supports_project_hypothesis=True, contradicts_project_hypothesis=False,
        PMID="18158893", DOI="10.1016/j.molcel.2007.10.035", year=2007, journal="Molecular Cell",
        paper_title="FOXO3a is activated in response to hypoxic stress and inhibits HIF1-induced apoptosis via regulation of CITED2",
        exact_evidence_summary="Verbatim abstract: 'In fibroblasts as well as in breast cancer cells, FOXO3a inhibits HIF1-induced apoptosis by stimulating the transcription of CITED2, which results in reduced expression of the proapoptotic HIF1 target genes NIX and RTP801.'",
        limitations="Study perturbs FOXO3a, not CITED2, directly; specific breast cell-line identity and quantitative details from secondary aggregation, abstract text directly verified via NCBI eutils.",
        review_or_primary="primary", fully_verified_primary_text=False,
    ),
    dict(
        claim_id="CITED2-09",
        mechanistic_claim="CITED2 recruits p300 to promote MYC-driven transactivation of E2F3 and suppression of p21CIP1, acting as a proliferation/quiescence switch; a CITED2/MYC/E2F3/p21 signature predicts poor lung cancer prognosis.",
        evidence_level=3, direct_or_indirect="direct",
        breast_cancer_specific=False, ER_positive_specific=False, tamoxifen_specific=False, endocrine_resistance_specific=False,
        experimental_perturbation="lentiviral CITED2 knockdown/overexpression; MYC/EGFR/p21 co-manipulation; subcutaneous+orthotopic xenografts; doxycycline-inducible knockdown",
        model_system="A549, H1975, CL1-0/CL1-5, H1299 lung cancer lines; nude mice; 72 NSCLC patient samples",
        direction_of_effect="CITED2 up -> E2F3 up, p21 down, G1/S progression, tumor growth; knockdown -> G1/S arrest, tumor shrinkage, increased mouse survival (p=0.0004); CITED2-positive patients had worse OS",
        supports_project_hypothesis=True, contradicts_project_hypothesis=False,
        PMID="22814619", DOI="10.1038/cdd.2012.91", year=2012, journal="Cell Death & Differentiation",
        paper_title="CITED2 functions as a molecular switch of cytokine-induced proliferation and quiescence",
        exact_evidence_summary="'CITED2 interacts with MYC's C-terminus and recruits p300 to promote MYC-p300-mediated transactivation of E2F3'; knockdown 'significantly decreased the growth of orthotopically implanted lung tumors and increased the overall survival rate of recipient mice (P=0.0004).'",
        limitations="Entirely a lung cancer (NSCLC) system -- no breast cancer line tested; the only substantive CITED2-E2F link found in any cancer type.",
        review_or_primary="primary", fully_verified_primary_text=True,
    ),
    dict(
        claim_id="CITED2-10",
        mechanistic_claim="CITED2 acts as a molecular chaperone guiding PRMT5 and p300 to nucleolin, activating a nucleolin-AKT translational axis that drives prostate cancer EMT and metastasis.",
        evidence_level=3, direct_or_indirect="direct",
        breast_cancer_specific=False, ER_positive_specific=False, tamoxifen_specific=False, endocrine_resistance_specific=False,
        experimental_perturbation="siRNA/shRNA knockdown and overexpression of CITED2; orthotopic xenograft",
        model_system="PC3, DU145, VCaP, LNCaP, C42B, 22RV1 (prostate cancer lines); HEK293T",
        direction_of_effect="CITED2 overexpression -> increased AKT translation (via nucleolin), increased EMT/migration, enhanced xenograft metastasis; correlates with poor survival and Gleason score",
        supports_project_hypothesis=True, contradicts_project_hypothesis=False,
        PMID="30291252", DOI="10.1038/s41467-018-06606-2", year=2018, journal="Nature Communications",
        paper_title="Aberrant expression of CITED2 promotes prostate cancer metastasis by activating the nucleolin-AKT pathway",
        exact_evidence_summary="'CITED2 acts as a molecular chaperone to guide PRMT5 and p300 to nucleolin, thereby activating nucleolin'; xenograft metastasis 'significantly enhanced by CITED2 overexpression.'",
        limitations="Prostate cancer, not breast -- related-mechanism context only.",
        review_or_primary="primary", fully_verified_primary_text=True,
    ),
    dict(
        claim_id="CITED2-11",
        mechanistic_claim="In an in vivo genome-wide CRISPR screen, CITED2 was identified as a driver of prostate cancer bone metastasis, validated with an organ-on-a-chip bone-tropism invasion platform.",
        evidence_level=3, direct_or_indirect="direct",
        breast_cancer_specific=False, ER_positive_specific=False, tamoxifen_specific=False, endocrine_resistance_specific=False,
        experimental_perturbation="genome-wide in vivo CRISPR activation/inhibition screening; organ-on-a-chip functional validation",
        model_system="prostate cancer cell lines/mouse models",
        direction_of_effect="CITED2 identified as a driver (tolerance/gain-of-function direction) of bone-tropic metastasis",
        supports_project_hypothesis=True, contradicts_project_hypothesis=False,
        PMID="38454137", DOI="10.1038/s41388-024-02995-5", year=2024, journal="Oncogene",
        paper_title="In vivo genome-wide CRISPR screening identifies CITED2 as a driver of prostate cancer bone metastasis",
        exact_evidence_summary="Bibliographically confirmed (title/authors Arriaga et al./journal/year/PMID/DOI/pages 1303-1315) across two independent sources; a Correction was later issued (PMID 38627522), content of correction not reviewed.",
        limitations="Full text not independently fetched -- bibliographic-level confirmation only; prostate, not breast; methodologically parallel to our own CRISPR-screen approach, which is an analogy, not evidence transfer.",
        review_or_primary="primary", fully_verified_primary_text=False,
    ),
    dict(
        claim_id="CITED2-12",
        mechanistic_claim="CITED2 physically binds EP300 and this interaction is required for downstream HSPA6 expression supporting proliferation and survival.",
        evidence_level=3, direct_or_indirect="direct",
        breast_cancer_specific=False, ER_positive_specific=False, tamoxifen_specific=False, endocrine_resistance_specific=False,
        experimental_perturbation="siRNA CITED2 knockdown; co-IP; immunofluorescence; RNA-seq; HSPA6 rescue",
        model_system="immortalized human spermatogonial stem cell line",
        direction_of_effect="CITED2 knockdown -> reduced EP300 protein, reduced HSPA6, reduced proliferation, increased apoptosis; HSPA6 re-expression rescues",
        supports_project_hypothesis=True, contradicts_project_hypothesis=False,
        PMID="40313859", DOI="10.1155/sci/2362489", year=2025, journal="Stem Cells International",
        paper_title="CITED2 Binding to EP300 Regulates Human Spermatogonial Stem Cell Proliferation and Survival Through HSPA6",
        exact_evidence_summary="Co-IP 'confirmed a specific interaction between CITED2 and EP300'; ~75% of CITED2-positive cells co-expressed EP300; knockdown reduced EP300 protein and HSPA6, rescued by HSPA6 restoration.",
        limitations="Normal (non-cancer) human stem cell biology, not breast tissue -- direct physical-interaction evidence for CITED2-EP300 only.",
        review_or_primary="primary", fully_verified_primary_text=True,
    ),
    dict(
        claim_id="CITED2-13",
        mechanistic_claim="CONTRADICTION: in colon cancer cells, CITED2 knockdown INCREASES invasiveness (via induced MMP-13) -- the opposite direction from CITED2's pro-invasive role in breast/prostate cancer.",
        evidence_level=3, direct_or_indirect="direct",
        breast_cancer_specific=False, ER_positive_specific=False, tamoxifen_specific=False, endocrine_resistance_specific=False,
        experimental_perturbation="CITED2 knockdown; gene expression profiling; butyrate (HDAC inhibitor) treatment; ectopic CITED2 re-expression",
        model_system="RKO colon cancer cell line",
        direction_of_effect="CITED2 knockdown -> INCREASED invasiveness, MMP-13 induction; ectopic CITED2 arrested growth; butyrate induces CITED2, suppresses MMP-13",
        supports_project_hypothesis=False, contradicts_project_hypothesis=True,
        PMID="18054336", DOI="10.1016/j.febslet.2007.11.072", year=2007, journal="FEBS Letters",
        paper_title="A role for CITED2, a CBP/p300 interacting protein, in colon cancer cell invasion",
        exact_evidence_summary="Bibliographic record confirmed via Semantic Scholar API (authors Bai L, Merchant J). CITED2 knockdown increases invasiveness -- anti-invasive/growth-suppressive role, opposite to breast/prostate literature.",
        limitations="Colon, not breast; full primary text paywalled, mechanistic detail from secondary aggregation only.",
        review_or_primary="primary", fully_verified_primary_text=False,
    ),
    dict(
        claim_id="CITED2-14",
        mechanistic_claim="CONTRADICTION: in hepatocellular carcinoma, CITED2 is a direct PPARgamma-induced tumor-suppressor effector -- knockdown increases HCC cell viability/clonogenicity, overexpression suppresses HCC growth -- opposite the oncogenic direction seen in breast cancer.",
        evidence_level=3, direct_or_indirect="direct",
        breast_cancer_specific=False, ER_positive_specific=False, tamoxifen_specific=False, endocrine_resistance_specific=False,
        experimental_perturbation="CITED2 knockdown and overexpression; PPARgamma agonist (rosiglitazone); ChIP-PCR; viability/clonogenicity/cell-cycle assays",
        model_system="LO2 (immortalized hepatocyte) and Hep3B (HCC) cell lines",
        direction_of_effect="CITED2 knockdown -> INCREASED viability/clonogenicity/G1-S transition; overexpression -> SUPPRESSED HCC growth (tumor-suppressive)",
        supports_project_hypothesis=False, contradicts_project_hypothesis=True,
        PMID="23212831", DOI="10.1002/cncr.27865", year=2013, journal="Cancer",
        paper_title="CITED2 is a novel direct effector of peroxisome proliferator-activated receptor gamma in suppressing hepatocellular carcinoma cell growth",
        exact_evidence_summary="Bibliographic record confirmed via Semantic Scholar API. CITED2 identified as 'the most prominent PPARgamma-bound target' among 114 candidates by ChIP-PCR; functions as growth suppressor in this system.",
        limitations="Liver, not breast; full primary text paywalled, mechanistic detail from secondary aggregation only; directly opposite direction to the breast/prostate oncogenic pattern.",
        review_or_primary="primary", fully_verified_primary_text=False,
    ),
    dict(
        claim_id="CITED2-15",
        mechanistic_claim="Acute loss of Cited2 impairs Nanog expression and rapidly decreases self-renewal/survival of mouse embryonic stem cells.",
        evidence_level=3, direct_or_indirect="direct",
        breast_cancer_specific=False, ER_positive_specific=False, tamoxifen_specific=False, endocrine_resistance_specific=False,
        experimental_perturbation="conditional Cre-loxP knockout; RNAi knockdown; ChIP on Nanog promoter",
        model_system="mouse embryonic stem cells",
        direction_of_effect="Cited2 loss -> decreased Nanog, impaired self-renewal, spontaneous differentiation within 48h",
        supports_project_hypothesis=False, contradicts_project_hypothesis=False,
        PMID="25377420", DOI="10.1002/stem.1889", year=2015, journal="Stem Cells",
        paper_title="Acute Loss of Cited2 Impairs Nanog Expression and Decreases Self-Renewal of Mouse Embryonic Stem Cells",
        exact_evidence_summary="'significant enrichment of the Nanog proximal promoter' by CITED2 ChIP; Cited2-depleted cells show 'impaired self-renewal capacity.'",
        limitations="Normal mouse ESC biology, not cancer or breast tissue -- background only, no cancer-stem-cell-specific CITED2 breast paper was found.",
        review_or_primary="primary", fully_verified_primary_text=True,
    ),
]

# ---------------------------------------------------------------------------
# VEZF1
# ---------------------------------------------------------------------------
VEZF1_CLAIMS = [
    dict(
        claim_id="VEZF1-01",
        mechanistic_claim="Vezf1 knockout is embryonic lethal with vascular remodeling defects/hemorrhage; heterozygotes show lymphatic hypervascularization/edema (haploinsufficiency) -- VEZF1 is essential for developmental blood-vessel and lymphatic morphogenesis.",
        evidence_level=3, direct_or_indirect="direct",
        breast_cancer_specific=False, ER_positive_specific=False, tamoxifen_specific=False, endocrine_resistance_specific=False,
        experimental_perturbation="targeted gene inactivation (Vezf1-null mice), homozygous and heterozygous, analyzed E8.5-birth",
        model_system="mouse, in vivo embryo",
        direction_of_effect="Vezf1 loss -> vascular remodeling defects, loss of vascular integrity, hemorrhage (homozygous); lymphatic hypervascularization/edema (heterozygous)",
        supports_project_hypothesis=True, contradicts_project_hypothesis=False,
        PMID="15882861", DOI="10.1016/j.ydbio.2005.04.003", year=2005, journal="Developmental Biology",
        paper_title="Dosage-dependent requirement for mouse Vezf1 in vascular system development",
        exact_evidence_summary="'Homozygous mutant embryos display vascular remodeling defects and loss of vascular integrity leading to localized hemorrhaging.' Heterozygotes: 'lymphatic hypervascularization associated with hemorrhaging and edema in the jugular region.'",
        limitations="Purely developmental mouse embryology, no adult/cancer/vascular-pathology context.",
        review_or_primary="primary", fully_verified_primary_text=True,
    ),
    dict(
        claim_id="VEZF1-02",
        mechanistic_claim="Founding paper: VEZF1 is a zinc-finger transcription factor whose expression is restricted to vascular endothelial cells and their precursors during mouse development.",
        evidence_level=3, direct_or_indirect="direct",
        breast_cancer_specific=False, ER_positive_specific=False, tamoxifen_specific=False, endocrine_resistance_specific=False,
        experimental_perturbation="none - gene identification/expression cloning, in situ hybridization",
        model_system="mouse embryo, retroviral entrapment vector screen",
        direction_of_effect="not applicable (discovery/expression paper)",
        supports_project_hypothesis=False, contradicts_project_hypothesis=False,
        PMID="9986727", DOI="10.1006/dbio.1998.9144", year=1999, journal="Developmental Biology",
        paper_title="Vezf1: a Zn finger transcription factor restricted to endothelial cells and their precursors",
        exact_evidence_summary="Gene-trap screen identified a gene expressed exclusively in vascular endothelial cells/precursors in yolk-sac blood islands, hence 'vascular endothelial zinc finger 1.'",
        limitations="Full text paywalled (ScienceDirect 403); bibliographic details verified via NCBI PubMed esummary, narrative content from secondary search-snippet only.",
        review_or_primary="primary", fully_verified_primary_text=False,
    ),
    dict(
        claim_id="VEZF1-03",
        mechanistic_claim="Vezf1 knockdown reduces endothelial cell proliferation, migration, network formation and in vivo angiogenesis, acting in part through downstream target gene stathmin/OP18.",
        evidence_level=3, direct_or_indirect="direct",
        breast_cancer_specific=False, ER_positive_specific=False, tamoxifen_specific=False, endocrine_resistance_specific=False,
        experimental_perturbation="antisense oligodeoxynucleotide knockdown of Vezf1 and, separately, stathmin/OP18",
        model_system="cultured endothelial cells + in vivo angiogenesis assay (non-cancer)",
        direction_of_effect="Vezf1 knockdown -> decreased EC proliferation/migration/tube formation, increased apoptosis, decreased in vivo angiogenesis (VEZF1 is pro-angiogenic)",
        supports_project_hypothesis=True, contradicts_project_hypothesis=False,
        PMID="15031128", DOI="", year=2004, journal="Arteriosclerosis, Thrombosis, and Vascular Biology",
        paper_title="Vascular endothelial zinc finger 1 is involved in the regulation of angiogenesis: possible contribution of stathmin/OP18 as a downstream target gene",
        exact_evidence_summary="AS-ODN knockdown 'significantly inhibited the proliferation, migration, and network formation of cultured ECs as well as angiogenesis in vivo.' Knockdown of stathmin/OP18 phenocopied Vezf1 knockdown.",
        limitations="Non-cancer endothelial system; DOI not independently verified.",
        review_or_primary="primary", fully_verified_primary_text=False,
    ),
    dict(
        claim_id="VEZF1-04",
        mechanistic_claim="Vezf1/DB1 is an endothelial-cell-specific transcription factor that directly regulates the endothelin-1 gene promoter -- one of its earliest identified direct target genes.",
        evidence_level=3, direct_or_indirect="direct",
        breast_cancer_specific=False, ER_positive_specific=False, tamoxifen_specific=False, endocrine_resistance_specific=False,
        experimental_perturbation="promoter-reporter assays, EMSA (exact perturbation not independently fetched)",
        model_system="endothelial cell lines",
        direction_of_effect="VEZF1 activates/regulates endothelin-1 promoter",
        supports_project_hypothesis=False, contradicts_project_hypothesis=False,
        PMID="11504723", DOI="10.1074/jbc.M105166200", year=2001, journal="Journal of Biological Chemistry",
        paper_title="Vezf1/DB1 is an endothelial cell-specific transcription factor that regulates expression of the endothelin-1 promoter",
        exact_evidence_summary="Bibliographic metadata confirmed via PubMed esummary; content drawn from secondary search-engine summaries, not independently fetched full text.",
        limitations="Full text not fetched directly -- mechanistic detail beyond title UNVERIFIED.",
        review_or_primary="primary", fully_verified_primary_text=False,
    ),
    dict(
        claim_id="VEZF1-05",
        mechanistic_claim="VEZF1 is required for normal genomic DNA methylation via direct transcriptional activation of DNA methyltransferase Dnmt3b -- an epigenetic regulatory role, not a heme-metabolism role.",
        evidence_level=3, direct_or_indirect="direct",
        breast_cancer_specific=False, ER_positive_specific=False, tamoxifen_specific=False, endocrine_resistance_specific=False,
        experimental_perturbation="Vezf1-/- mouse ES cells vs wild-type; rescue by Vezf1 re-expression; ChIP for Vezf1 at Dnmt3b locus",
        model_system="mouse embryonic stem cells",
        direction_of_effect="Vezf1 loss -> ~4-fold reduction in Dnmt3b1 transcript -> widespread loss of DNA methylation; re-expression rescues",
        supports_project_hypothesis=False, contradicts_project_hypothesis=False,
        PMID="18676812", DOI="10.1101/gad.1658408", year=2008, journal="Genes & Development",
        paper_title="Vezf1 regulates genomic DNA methylation through its effects on expression of DNA methyltransferase Dnmt3b",
        exact_evidence_summary="'Widespread loss of DNA methylation in Vezf1-/- cells' with a 'fourfold reduction in the Dnmt3b1 transcript.' ChIP identified a functional Vezf1-binding site in Dnmt3b intron 3.",
        limitations="ES cell system; no erythroid/heme-pathway relevance despite VEZF1's HALLMARK_HEME_METABOLISM gene-set co-membership in the project's data -- this paper is about global epigenetic regulation, not erythropoiesis/heme biosynthesis.",
        review_or_primary="primary", fully_verified_primary_text=True,
    ),
    dict(
        claim_id="VEZF1-06",
        mechanistic_claim="VEZF1 directly represses the antiangiogenic factor CITED2 in endothelial cells; loss of VEZF1 causes an endothelial differentiation/survival defect rescued by knocking down CITED2.",
        evidence_level=3, direct_or_indirect="direct",
        breast_cancer_specific=False, ER_positive_specific=False, tamoxifen_specific=False, endocrine_resistance_specific=False,
        experimental_perturbation="Vezf1-/- mouse ES cells differentiated to endothelial cells; ChIP for Vezf1 at Cited2 promoter; shRNA Cited2 knockdown rescue in Vezf1-/- background",
        model_system="mouse embryonic stem cell-derived endothelial cells",
        direction_of_effect="VEZF1 loss -> Cited2 up 4-5 fold, >80% cell death during EC differentiation, reduced tube formation; Cited2 knockdown rescues tube formation",
        supports_project_hypothesis=True, contradicts_project_hypothesis=False,
        PMID="29794136", DOI="10.1074/jbc.RA118.002911", year=2018, journal="Journal of Biological Chemistry",
        paper_title="The transcription factor Vezf1 represses the expression of the antiangiogenic factor Cited2 in endothelial cells",
        exact_evidence_summary="'Vezf1-/- ESCs have significantly increased expression of...Cited2,' 4-5 fold higher than WT; 'over 80% Vezf1-/- cells died' during EC differentiation; 'repression of Cited2...largely rescued the defective tube formation.' ChIP showed direct Vezf1 binding at the Cited2 promoter.",
        limitations="ES-cell-derived EC differentiation model, not adult/tumor vasculature or cancer -- BUT this is an independently documented literature link between two of the project's own four frozen candidates (VEZF1 directly represses CITED2), discovered by the literature search, not assumed by the project's own network analysis (which found no direct interaction between them).",
        review_or_primary="primary", fully_verified_primary_text=True,
    ),
    dict(
        claim_id="VEZF1-07",
        mechanistic_claim="Nuclear RhoB-GTP regulates lineage-specific gene expression in blood vs lymphatic endothelium by controlling VEZF1-mediated transcription; VEZF1 directly regulates VEGFR2, NRP1 (blood-vessel genes) and TIMP3, MMP2, Vasohibin1 (lymphatic genes).",
        evidence_level=3, direct_or_indirect="direct",
        breast_cancer_specific=False, ER_positive_specific=False, tamoxifen_specific=False, endocrine_resistance_specific=False,
        experimental_perturbation="RhoB-null mice; VEZF1-heterozygous mice; RhoB siRNA/overexpression/dominant-negative in primary human blood/lymphatic ECs; oxygen-induced retinopathy and dermal wounding models; small-molecule VEZF1-DNA-binding inhibitor",
        model_system="mouse (retinopathy, wound healing) + primary human blood/lymphatic endothelial cells",
        direction_of_effect="RhoB loss decreases pathological retinal angiogenesis but increases lymphangiogenesis after injury, via altered VEZF1 target-gene expression",
        supports_project_hypothesis=False, contradicts_project_hypothesis=False,
        PMID="", DOI="10.1038/ncomms3824", year=2013, journal="Nature Communications",
        paper_title="RhoB controls coordination of adult angiogenesis and lymphangiogenesis following injury by regulating VEZF1-mediated transcription",
        exact_evidence_summary="'Nuclear RhoB-GTP controls expression of distinct gene sets in each endothelial lineage by regulating VEZF1-mediated transcription.' Direct VEZF1 target genes: VEGFR2, NRP1 (blood); TIMP3, MMP2, Vasohibin1 (lymphatic). Tumor relevance mentioned only in discussion ('RhoB can alter tumour formation... but yet promoting tumour angiogenesis'), not tested experimentally.",
        limitations="Closest literature comes to a tumor-angiogenesis link for VEZF1, but only a discussion-level mention -- actual experiments are in ischemic retinopathy and wound healing, non-cancer contexts.",
        review_or_primary="primary", fully_verified_primary_text=True,
    ),
    dict(
        claim_id="VEZF1-08",
        mechanistic_claim="In hepatocellular carcinoma, VEZF1 is an oncoprotein: knockdown suppresses, overexpression promotes, proliferation and metastasis, acting by transcriptionally activating PAQR4; VEZF1 stability is controlled by STUB1-mediated ubiquitination.",
        evidence_level=3, direct_or_indirect="direct",
        breast_cancer_specific=False, ER_positive_specific=False, tamoxifen_specific=False, endocrine_resistance_specific=False,
        experimental_perturbation="VEZF1 silencing and overexpression in HCC cell lines; STUB1 manipulation",
        model_system="human HCC cell lines and clinical HCC tissue",
        direction_of_effect="VEZF1 knockdown -> suppressed proliferation/metastasis; overexpression -> promoted (via PAQR4 transactivation)",
        supports_project_hypothesis=True, contradicts_project_hypothesis=False,
        PMID="36241701", DOI="10.1038/s41417-022-00540-8", year=2023, journal="Cancer Gene Therapy",
        paper_title="VEZF1, destabilized by STUB1, affects cellular growth and metastasis of hepatocellular carcinoma by transcriptionally regulating PAQR4",
        exact_evidence_summary="'VEZF1 has been recognized as an oncoprotein in certain types of cancer.' Silencing suppressed, overexpression promoted, HCC proliferation/metastasis; 'VEZF1 transcriptionally activated...PAQR4.' STUB1 ubiquitinates/destabilizes VEZF1.",
        limitations="Liver cancer only, not breast; no endocrine/estrogen/tamoxifen context; mechanism is entirely non-vascular (PAQR4, not angiogenesis).",
        review_or_primary="primary", fully_verified_primary_text=True,
    ),
    dict(
        claim_id="VEZF1-09",
        mechanistic_claim="In HCC, GFAT1 drives O-GlcNAcylation of VEZF1, stabilizing VEZF1 protein and increasing VEZF1-driven transcription of TNS1, promoting tumor growth/metastasis; a peptide blocking this modification suppresses HCC in vivo.",
        evidence_level=3, direct_or_indirect="direct",
        breast_cancer_specific=False, ER_positive_specific=False, tamoxifen_specific=False, endocrine_resistance_specific=False,
        experimental_perturbation="GFAT1 WT vs catalytically-dead overexpression; GFAT1 knockdown; OGT/OGA inhibitors; non-glycosylatable VEZF1 mutant; cell-penetrating peptide (CPPtat-V1); xenograft/orthotopic/metastasis mouse models",
        model_system="human HCC cell lines (SNU449, Huh7, PLC/PRF/5, HepG2, MHCC-97H) + HCC mouse models",
        direction_of_effect="GFAT1 up -> VEZF1 O-GlcNAcylation up -> VEZF1 ubiquitination down/stability up -> TNS1 up -> HCC proliferation/invasion/metastasis up; CPPtat-V1 reverses and suppresses HCC in vivo",
        supports_project_hypothesis=True, contradicts_project_hypothesis=False,
        PMID="40858565", DOI="10.1038/s41419-025-07975-5", year=2025, journal="Cell Death & Disease",
        paper_title="GFAT1 promotes the progression of hepatocellular carcinoma via enhancing the O-GlcNAcylation of VEZF1",
        exact_evidence_summary="ChIP-Seq/RNA-Seq identified TNS1 as a VEZF1 target; non-glycosylatable VEZF1 mutant showed 'shortened half-life and increased ubiquitination.' CPPtat-V1 'significantly suppressed HCC progression' in vivo.",
        limitations="Liver cancer, not breast; the 'stress' link is a general nutrient/hexosamine metabolic-flux-sensing pathway, not drug/therapy stress or hypoxia; no endocrine relevance.",
        review_or_primary="primary", fully_verified_primary_text=True,
    ),
    dict(
        claim_id="VEZF1-10-DMTN-NULL",
        mechanistic_claim="No documented physical, transcriptional-regulatory, or functional relationship between VEZF1 and DMTN (dematin) was found anywhere in the accessible literature -- searched exhaustively and confirmed absent.",
        evidence_level=5, direct_or_indirect="not applicable",
        breast_cancer_specific=False, ER_positive_specific=False, tamoxifen_specific=False, endocrine_resistance_specific=False,
        experimental_perturbation="none - exhaustive search, no paper exists",
        model_system="not applicable",
        direction_of_effect="not applicable",
        supports_project_hypothesis=False, contradicts_project_hypothesis=False,
        PMID="", DOI="", year=pd.NA, journal="", paper_title="(no paper found -- explicit negative result)",
        exact_evidence_summary="Searches for 'VEZF1 DMTN' / 'VEZF1 dematin' returned only unrelated dematin papers (erythrocyte membrane, malaria) and unrelated VEZF1 papers (Dnmt3b), zero cross-mentions. The closest VEZF1-blood-lineage paper (Das et al. 2023, ETV2-VEZF1 hematoendothelial interaction, DOI 10.3389/fcell.2023.1109648) was directly checked and contains no mention of erythroid genes, heme metabolism, or DMTN/dematin -- only endothelial/early hematopoietic markers.",
        limitations="Absence of evidence in searchable literature is not proof of absence in unpublished/paywalled data, but no relationship is documented in any accessible source checked. Project's VEZF1-DMTN link is a pathway co-membership (gene-set) inference only, not a validated biological relationship.",
        review_or_primary="not applicable", fully_verified_primary_text=True,
    ),
    dict(
        claim_id="VEZF1-11-BREAST-NULL",
        mechanistic_claim="Human Protein Atlas pan-cancer (TCGA) survival analysis finds VEZF1 a significant prognostic marker (p<0.001) in liver, lung, and renal cancers, but breast cancer does NOT meet that significance threshold -- no bulk-tumor prognostic signal for VEZF1 in breast cancer at the population level.",
        evidence_level=4, direct_or_indirect="indirect",
        breast_cancer_specific=True, ER_positive_specific=False, tamoxifen_specific=False, endocrine_resistance_specific=False,
        experimental_perturbation="none - correlative only (TCGA/GEO expression vs survival)",
        model_system="human tumor cohorts (TCGA/GEO), bulk expression database",
        direction_of_effect="unfavorable in LIHC/LUAD/KICH; favorable in KIRC; no significant association for breast cancer",
        supports_project_hypothesis=False, contradicts_project_hypothesis=False,
        PMID="", DOI="", year=pd.NA, journal="Human Protein Atlas (database resource)",
        paper_title="(database resource, not a paper)",
        exact_evidence_summary="'VEZF1 is a prognostic marker in:' LIHC, LUAD, KICH, KIRC (all p<0.001); breast cancer included in the general expression panel but not listed among prognostic cancers.",
        limitations="Bulk pan-cancer resource, not breast-specific or subtype-stratified; does not test ER+ subgroup, tamoxifen-treated cohorts, or resistant/recurrent disease specifically -- a null result here does not rule out a role restricted to a resistant/recurrent subpopulation, which is what the project's own GSE240112 finding addresses.",
        review_or_primary="not applicable (database)", fully_verified_primary_text=True,
    ),
]

# ---------------------------------------------------------------------------
# USP34
# ---------------------------------------------------------------------------
USP34_CLAIMS = [
    dict(
        claim_id="USP34-01",
        mechanistic_claim="USP34 deubiquitinates AXIN1, opposing tankyrase/RNF146-dependent PARsylation-driven ubiquitination, promoting AXIN1 nuclear accumulation, which is REQUIRED FOR (positively regulates, not inhibits) beta-catenin-mediated transcription.",
        evidence_level=3, direct_or_indirect="direct",
        breast_cancer_specific=False, ER_positive_specific=False, tamoxifen_specific=False, endocrine_resistance_specific=False,
        experimental_perturbation="siRNA/shRNA USP34 knockdown; catalytically-dead USP34 (C1903S) mutant; co-AP + LC-MS/MS of AXIN complexes; K48 ubiquitin-chain cleavage/deubiquitination assays; cycloheximide-chase for AXIN1 stability; IF for AXIN1 nuclear localization; TOPFlash/TCF luciferase reporter; qPCR for NKD1/TNFRSF19",
        model_system="HEK293T; RKO/SW480/HCT116 (human colorectal carcinoma); no breast cell line",
        direction_of_effect="USP34 knockdown -> AXIN1 degraded -> DECREASED Wnt3A-induced TOPFlash reporter and NKD1/TNFRSF19 target-gene expression; also inhibited constitutive beta-catenin signaling in RKO cells. USP34 = positive regulator of canonical Wnt/beta-catenin transcriptional output.",
        supports_project_hypothesis=True, contradicts_project_hypothesis=False,
        PMID="21383061", DOI="10.1128/mcb.01094-10", year=2011, journal="Molecular and Cellular Biology",
        paper_title="The ubiquitin-specific protease USP34 regulates axin stability and Wnt/beta-catenin signaling",
        exact_evidence_summary="'Analysis of purified axin-containing protein complexes by liquid chromatography-tandem mass spectrometry revealed the presence of the ubiquitin protease USP34... USP34 functions downstream of the beta-catenin destruction complex to control the stability of axin and opposes its tankyrase-dependent ubiquitination... interfering with USP34 function by RNA interference leads to the degradation of axin and to the inhibition of beta-catenin-mediated transcription.' Discussion: USP34 'is required at a step subsequent to beta-catenin stabilization' and 'positively modulate[s] Wnt signaling.'",
        limitations="No breast tissue/cell line; no in vivo/tumor model in this paper; discovered via unbiased AXIN interactome MS, not a targeted breast-cancer hypothesis test. NOTE: multiple secondary sources (a 2020 review, PMC7311976) mischaracterize this mechanism as Wnt-INHIBITORY (likely conflating USP34 with USP7); the primary source, verified 3 independent ways (abstract, TOPFlash data, target-gene qPCR), is unambiguous that USP34 loss DECREASES Wnt output.",
        review_or_primary="primary", fully_verified_primary_text=True,
    ),
    dict(
        claim_id="USP34-02",
        mechanistic_claim="USP34 inhibition in normal (non-transformed) mammary epithelial cells induces EMT and stemness, accompanied by INCREASED active (non-phospho) beta-catenin -- the opposite direction from USP34-01.",
        evidence_level=2, direct_or_indirect="direct",
        breast_cancer_specific=False, ER_positive_specific=False, tamoxifen_specific=False, endocrine_resistance_specific=False,
        experimental_perturbation="USP34 knockdown (siRNA/shRNA) in NMuMG mouse mammary epithelial cells; transplantation of knockdown cells into mouse mammary fat pads in vivo",
        model_system="NMuMG (normal murine mammary gland epithelial line) + in vivo mammary fat pad transplantation",
        direction_of_effect="USP34 knockdown -> increased N-cadherin, phospho-Smad3, Snail, active-beta-catenin; decreased E-cadherin; increased mammosphere formation and Nanog/Oct4/Sox2; knockdown cells regenerated ductal structures in vivo",
        supports_project_hypothesis=False, contradicts_project_hypothesis=True,
        PMID="28499884", DOI="10.1016/j.cellsig.2017.05.009", year=2017, journal="Cellular Signalling",
        paper_title="Inhibition of ubiquitin-specific protease 34 (USP34) induces epithelial-mesenchymal transition and promotes stemness in mammary epithelial cells",
        exact_evidence_summary="'Inhibition of USP34 in NMuMG cells induced EMT, as evidenced by upregulation of N-cadherin, phospho-Smad3, Snail and active-beta-catenin' and promoted stemness; knockdown cells transplanted in vivo 'reconstituted gland tissue with ductal development.'",
        limitations="Non-transformed mammary epithelial line, not human, not cancer, not ER+, not tamoxifen-tested. 'Active beta-catenin' increase co-occurs with increased phospho-Smad3 (TGF-beta/Smad activation), so the beta-catenin change may reflect EMT/TGF-beta crosstalk rather than the Axin-dependent mechanism in USP34-01 -- an unresolved directional contradiction between the two most relevant primary papers, not smoothed over here.",
        review_or_primary="primary", fully_verified_primary_text=True,
    ),
    dict(
        claim_id="USP34-03",
        mechanistic_claim="USP34 stabilizes eIF3m via deubiquitination; eIF3m upregulates MTCH2, maintaining mitochondrial function and TNBC proliferation.",
        evidence_level=2, direct_or_indirect="direct",
        breast_cancer_specific=True, ER_positive_specific=False, tamoxifen_specific=False, endocrine_resistance_specific=False,
        experimental_perturbation="USP34 silencing in TNBC cell lines; eIF3m and MTCH2 overexpression rescue",
        model_system="human TNBC cell lines",
        direction_of_effect="USP34 silencing -> mitochondrial dysfunction -> decreased TNBC proliferation; USP34 described as 'abnormally overexpressed' in TNBC",
        supports_project_hypothesis=True, contradicts_project_hypothesis=False,
        PMID="42023842", DOI="10.1080/01478885.2026.2648740", year=2026, journal="Journal of Histotechnology",
        paper_title="USP34 modulates mitochondrial function in triple-negative breast cancer cells through the eIf3m/MTCH2 axis",
        exact_evidence_summary="'USP34 has been predicted to be abnormally overexpressed in TNBC... silencing USP34 inhibits cell proliferation through mitochondrial dysfunction... USP34 stabilizes eIF3m protein via deubiquitination, eIF3m promotes MTCH2 expression by binding its 5'UTR region.'",
        limitations="Full text not accessible (abstract only, via structured database record); mechanism entirely unrelated to Wnt/Axin/beta-catenin; TNBC is ER-negative by definition -- no bearing on tamoxifen/ER+ biology.",
        review_or_primary="primary", fully_verified_primary_text=False,
    ),
    dict(
        claim_id="USP34-04",
        mechanistic_claim="USP34 physically interacts with and deubiquitinates/stabilizes SOX2, promoting laryngeal squamous cell carcinoma survival and cisplatin resistance.",
        evidence_level=3, direct_or_indirect="direct",
        breast_cancer_specific=False, ER_positive_specific=False, tamoxifen_specific=False, endocrine_resistance_specific=False,
        experimental_perturbation="USP34 knockdown in LSCC cells and cisplatin-resistant LSCC sublines; SOX2 overexpression rescue",
        model_system="human laryngeal squamous cell carcinoma (LSCC) tumor tissue and cell lines",
        direction_of_effect="USP34 knockdown -> reduced SOX2 protein (via increased polyubiquitination) -> decreased LSCC growth; SOX2 overexpression reverses; USP34 knockdown in cisplatin-resistant cells restored cisplatin sensitivity",
        supports_project_hypothesis=True, contradicts_project_hypothesis=False,
        PMID="32783291", DOI="10.1002/kjm2.12285", year=2020, journal="The Kaohsiung Journal of Medical Sciences",
        paper_title="The deubiquitinase USP34 stabilizes SOX2 and induces cell survival and drug resistance in laryngeal squamous cell carcinoma",
        exact_evidence_summary="'USP34 and SOX2 were elevated in LSCC tumor tissues' with positive correlation; USP34 'stabilized SOX2' through reduced polyubiquitination; knockdown 'could enhance the drug sensitivity of cisplatin in the resistant cells.'",
        limitations="Not breast cancer, not endocrine therapy (cisplatin only); this SOX2 link is independent of, and unrelated to, the Wnt/Axin mechanism (USP34-01). Offers non-breast precedent for a USP34-SOX2 relationship, relevant context for (but not validation of) the project's own USP34->SOX2 STRING-hop hypothesis.",
        review_or_primary="primary", fully_verified_primary_text=False,
    ),
    dict(
        claim_id="USP34-05",
        mechanistic_claim="USP34 deubiquitinates and stabilizes c-Myc, promoting hepatocellular carcinoma proliferation, migration, invasion, and glycolysis.",
        evidence_level=3, direct_or_indirect="direct",
        breast_cancer_specific=False, ER_positive_specific=False, tamoxifen_specific=False, endocrine_resistance_specific=False,
        experimental_perturbation="siRNA USP34 knockdown in HCC lines; c-Myc overexpression rescue",
        model_system="human HCC cell lines",
        direction_of_effect="USP34 knockdown -> c-Myc protein decreased (ubiquitination increased, degradation accelerated) -> proliferation/migration/invasion/glycolysis all decreased; c-Myc overexpression reverses",
        supports_project_hypothesis=True, contradicts_project_hypothesis=False,
        PMID="40260316", DOI="10.5152/tjg.2025.24335", year=2025, journal="The Turkish Journal of Gastroenterology",
        paper_title="The Knockdown of USP34 Inhibits the Progression of Hepatocellular Carcinoma by Accelerating c-Myc Degradation",
        exact_evidence_summary="'Interference with USP34 resulted in increased ubiquitination levels of c-Myc'; knockdown 'did not affect the mRNA level of c-Myc, but notably decreased the protein level'; silencing 'suppressed the proliferation, migration, and invasion of HCC cells.'",
        limitations="Not breast/endocrine; c-Myc is a canonical Wnt target gene but no Wnt-pathway component was directly tested in this paper.",
        review_or_primary="primary", fully_verified_primary_text=True,
    ),
    dict(
        claim_id="USP34-06",
        mechanistic_claim="USP34 deubiquitinates and stabilizes Smad1 and RUNX2, required for BMP2-driven osteogenic differentiation and bone formation -- independent of Wnt, demonstrating USP34's substrate breadth.",
        evidence_level=3, direct_or_indirect="direct",
        breast_cancer_specific=False, ER_positive_specific=False, tamoxifen_specific=False, endocrine_resistance_specific=False,
        experimental_perturbation="USP34 depletion in human MSCs; conditional Usp34 knockout mice (Prx1-Cre/Osx-Cre); Smurf1 co-depletion rescue",
        model_system="human mesenchymal stem cells; mouse conditional knockout, in vivo bone phenotyping",
        direction_of_effect="USP34 depletion -> inhibited osteogenic differentiation; conditional KO mice -> low bone mass, blunted BMP2 response",
        supports_project_hypothesis=False, contradicts_project_hypothesis=False,
        PMID="30181118", DOI="10.15252/embj.201899398", year=2018, journal="The EMBO Journal",
        paper_title="Ubiquitin-specific protease USP34 controls osteogenic differentiation and bone formation by regulating BMP2 signaling",
        exact_evidence_summary="'USP34 stabilizes both Smad1 and RUNX2... depletion of Smurf1 restores the osteogenic potential of Usp34-deficient MSCs.' Conditional knockout 'leads to low bone mass in mice' and 'blunts BMP2-induced responses.'",
        limitations="Bone biology, not cancer; demonstrates USP34 is a pleiotropic, context-dependent DUB with multiple validated substrates beyond Axin, not a Wnt-dedicated enzyme.",
        review_or_primary="primary", fully_verified_primary_text=True,
    ),
    dict(
        claim_id="USP34-07",
        mechanistic_claim="An unbiased DUB siRNA screen identified USP34 (alongside CYLD) as a regulator of thrombin-induced GPCR-p38 MAPK inflammatory signaling in endothelial cells.",
        evidence_level=4, direct_or_indirect="direct",
        breast_cancer_specific=False, ER_positive_specific=False, tamoxifen_specific=False, endocrine_resistance_specific=False,
        experimental_perturbation="siRNA screen of 96 deubiquitinating enzymes; USP34 knockdown individually validated",
        model_system="human endothelial cells, thrombin (PAR1 GPCR agonist) stimulation",
        direction_of_effect="USP34 knockdown -> decreased thrombin-induced IL-6/inflammatory cytokine production, no effect on barrier permeability",
        supports_project_hypothesis=False, contradicts_project_hypothesis=False,
        PMID="37865315", DOI="10.1016/j.jbc.2023.105370", year=2023, journal="Journal of Biological Chemistry",
        paper_title="An siRNA library screen identifies CYLD and USP34 as deubiquitinases that regulate GPCR-p38 MAPK signaling and distinct inflammatory responses",
        exact_evidence_summary="USP34 depletion 'decreased inflammatory cytokine production without affecting barrier function.'",
        limitations="Endothelial biology, unrelated to cancer/Wnt/breast; further underscores USP34's pleiotropy (now 6 distinct documented pathways: Axin/Wnt, Smad1-RUNX2/BMP, eIF3m/MTCH2, c-Myc, SOX2, GPCR-p38).",
        review_or_primary="primary", fully_verified_primary_text=True,
    ),
    dict(
        claim_id="USP34-08-PATHWAY-CONTEXT",
        mechanistic_claim="PATHWAY-LEVEL CONTEXT, NOT USP34-SPECIFIC: loss of CXXC4 (a Wnt/beta-catenin pathway antagonist) activates canonical Wnt/beta-catenin signaling and directly causes tamoxifen resistance in ER+ breast cancer cells.",
        evidence_level=2, direct_or_indirect="direct",
        breast_cancer_specific=True, ER_positive_specific=True, tamoxifen_specific=True, endocrine_resistance_specific=True,
        experimental_perturbation="CXXC4 shRNA knockdown / lentiviral overexpression",
        model_system="MCF-7 and BT474 (ER+ breast cancer) plus tamoxifen-resistant sublines MCF-7/TMR and BT474/TMR",
        direction_of_effect="CXXC4 knockdown -> Wnt/beta-catenin activation (GSK-3beta phosphorylation, cyclin D1/c-Myc up) -> accelerated proliferation and tamoxifen INSENSITIVITY; CXXC4 overexpression -> increased tamoxifen sensitivity",
        supports_project_hypothesis=True, contradicts_project_hypothesis=False,
        PMID="31911277", DOI="10.1016/j.tranon.2019.12.005", year=2020, journal="Translational Oncology",
        paper_title="Downregulation of CXXC Finger Protein 4 Leads to a Tamoxifen-resistant Phenotype in Breast Cancer Cells Through Activation of the Wnt/beta-catenin Pathway",
        exact_evidence_summary="'CXXC4 knockdown accelerates cell proliferation' and 'renders breast cancer cells insensitive to tamoxifen,' via Wnt/beta-catenin activation (GSK-3beta phosphorylation, cyclin D1, c-Myc upregulation). USP34 does not appear anywhere in this paper (confirmed via direct full-text fetch).",
        limitations="This is CXXC4 evidence, not USP34 evidence -- cited only to establish that canonical Wnt/beta-catenin activation is an independently documented, CAUSAL driver of tamoxifen resistance in ER+ breast cancer cell lines, which is the pathway our project's own USP34 finding sits inside, not a claim about USP34 itself.",
        review_or_primary="primary", fully_verified_primary_text=True,
    ),
]

# ---------------------------------------------------------------------------
# EML5
# ---------------------------------------------------------------------------
EML5_CLAIMS = [
    dict(
        claim_id="EML5-01",
        mechanistic_claim="EML5 protein/mRNA was elevated in neurons/glia of resected human epileptic temporal neocortex vs control brain tissue -- the only human-tissue EML5-specific disease-expression study found.",
        evidence_level=4, direct_or_indirect="direct",
        breast_cancer_specific=False, ER_positive_specific=False, tamoxifen_specific=False, endocrine_resistance_specific=False,
        experimental_perturbation="none - Western blot/IHC/IF/qPCR expression comparison, no EML5 perturbation",
        model_system="human surgical brain tissue (36 epilepsy patients vs 8 controls)",
        direction_of_effect="EML5 UP in diseased (epileptic) tissue (protein OD 1.80 vs 1.19, p<0.05; mRNA >5-fold by FQ-PCR)",
        supports_project_hypothesis=False, contradicts_project_hypothesis=False,
        PMID="26730336", DOI="", year=2015, journal="Iranian Journal of Basic Medical Sciences",
        paper_title="Echinoderm microtubule-associated protein-like protein 5 in anterior temporal neocortex of patients with intractable epilepsy",
        exact_evidence_summary="Western blot OD 1.8030+/-0.1335 (epilepsy) vs 1.1852+/-0.2253 (control), P<0.05; mRNA >5-fold higher in epilepsy by FQ-PCR.",
        limitations="Small control group (n=8); no functional perturbation; disease context is epilepsy/CNS, not cancer; only generically consistent with 'EML5 can be up in a pathological state,' not specific to resistance biology.",
        review_or_primary="primary", fully_verified_primary_text=True,
    ),
    dict(
        claim_id="EML5-02",
        mechanistic_claim="Cloning/characterization of rat Eml5: WD40+HELP-domain protein homologous to sea-urchin EMAP, expressed in hippocampus/cerebellum/olfactory bulb; a microtubule-regulatory role is INFERRED by homology to other EML paralogs, not directly demonstrated for Eml5 in this paper.",
        evidence_level=4, direct_or_indirect="indirect (homology-based inference for the functional claim; direct for expression pattern)",
        breast_cancer_specific=False, ER_positive_specific=False, tamoxifen_specific=False, endocrine_resistance_specific=False,
        experimental_perturbation="none",
        model_system="rat, in situ hybridization/immunocytochemistry",
        direction_of_effect="not applicable (developmental/tissue expression characterization)",
        supports_project_hypothesis=False, contradicts_project_hypothesis=False,
        PMID="15225882", DOI="10.1016/j.gene.2004.04.012", year=2004, journal="Gene",
        paper_title="Eml5, a novel WD40 domain protein expressed in rat brain",
        exact_evidence_summary="'Eml5 contains 11 putative WD40 domains and 3 hydrophobic stretches...HELP domains, which have been suggested to be involved in microtubule binding... it is likely that Eml5 plays a role in the regulation of cytoskeletal rearrangements' (explicitly hedged, no direct microtubule-binding assay performed on Eml5 itself).",
        limitations="Rat, not human; no direct microtubule co-sedimentation/binding assay on Eml5 in this paper; pre-genomic-era paper.",
        review_or_primary="primary", fully_verified_primary_text=True,
    ),
    dict(
        claim_id="EML5-03-ANNOTATION",
        mechanistic_claim="UniProt (Q05BV3) FUNCTION annotation: 'May modify the assembly dynamics of microtubules' -- evidence code ECO:0000250 (inferred by sequence similarity to a characterized homolog), NOT experimentally demonstrated for human EML5 itself.",
        evidence_level=5, direct_or_indirect="indirect (homology/computational inference)",
        breast_cancer_specific=False, ER_positive_specific=False, tamoxifen_specific=False, endocrine_resistance_specific=False,
        experimental_perturbation="none",
        model_system="not applicable (database annotation)",
        direction_of_effect="not applicable",
        supports_project_hypothesis=False, contradicts_project_hypothesis=False,
        PMID="", DOI="", year=pd.NA, journal="UniProtKB/Swiss-Prot Q05BV3",
        paper_title="(database annotation, not a paper)",
        exact_evidence_summary="'FUNCTION: May modify the assembly dynamics of microtubules...{ECO:0000250}.' 29 WD40 repeats + 3 HELP domains (Pfam PF03451), structurally consistent with the EML family but functionally uncharacterized by direct experiment on human EML5.",
        limitations="Curated database inference, not a cited primary experimental paper on EML5's microtubule activity -- NO SUPPORT at the experimental level.",
        review_or_primary="not applicable (database)", fully_verified_primary_text=True,
    ),
    dict(
        claim_id="EML5-04-INTERACTOME",
        mechanistic_claim="STRING v12 and BioGRID/HPA show only weak, low-medium-confidence, mostly text-mining-driven computational associations for EML5 (e.g. MAD2L1, IFT140, TTC8, FBXW11) -- none independently validated by any primary paper found; none exceed STRING's high-confidence threshold.",
        evidence_level=5, direct_or_indirect="indirect (computational/database aggregation)",
        breast_cancer_specific=False, ER_positive_specific=False, tamoxifen_specific=False, endocrine_resistance_specific=False,
        experimental_perturbation="none",
        model_system="not applicable",
        direction_of_effect="not applicable",
        supports_project_hypothesis=False, contradicts_project_hypothesis=False,
        PMID="", DOI="", year=pd.NA, journal="STRING v12 / BioGRID / BioGRID ORCS / Human Protein Atlas (database queries)",
        paper_title="(database aggregation, not a paper)",
        exact_evidence_summary="STRING API (required_score=150): MAD2L1 combined score 0.415 (experimental subscore 0.31); IFT140 0.425 (0.312); TTC8 0.52 -- all below STRING's 0.7 high-confidence threshold, dominated by text-mining. BioGRID ORCS: EML5 flagged in 8 of ~1,400 curated genome-wide CRISPR screens, but the specific screens could not be extracted (JS-rendered interface) -- unverified lead, not citable.",
        limitations="Raw, unfiltered, pan-genomic database aggregation; no dedicated hypothesis-driven study of EML5 exists to validate any of these associations.",
        review_or_primary="not applicable (database)", fully_verified_primary_text=True,
    ),
    dict(
        claim_id="EML5-05-BREAST-NULL",
        mechanistic_claim="Human Protein Atlas: EML5 shows no clear cancer prognostic signal across 34 cohorts generally; the ONE significant association is glioblastoma (higher EML5 = BETTER survival, p<0.001); no breast-cancer-specific prognostic association is reported.",
        evidence_level=4, direct_or_indirect="direct (EML5 expression itself, in the GBM cohort)",
        breast_cancer_specific=True, ER_positive_specific=False, tamoxifen_specific=False, endocrine_resistance_specific=False,
        experimental_perturbation="none - correlative only",
        model_system="TCGA pan-cancer cohorts",
        direction_of_effect="EML5 UP associated with BETTER survival in GBM (opposite valence to the project's finding of EML5 UP associated with resistance/recurrence in breast cancer -- different tissue context, not directly comparable)",
        supports_project_hypothesis=False, contradicts_project_hypothesis=False,
        PMID="", DOI="", year=pd.NA, journal="Human Protein Atlas (database resource)",
        paper_title="(database resource, not a paper)",
        exact_evidence_summary="HPA overall reliability score for the EML5 annotation explicitly marked 'Uncertain.' No breast-cancer prognostic/expression signal reported.",
        limitations="Database resource, not a peer-reviewed dedicated EML5 paper; breast cancer specifically shows NO signal in this resource, which is itself an informative negative result.",
        review_or_primary="not applicable (database)", fully_verified_primary_text=True,
    ),
]

ALL_CLAIMS = CITED2_CLAIMS + VEZF1_CLAIMS + USP34_CLAIMS + EML5_CLAIMS
for c in ALL_CLAIMS:
    c["candidate"] = c["claim_id"].split("-")[0]
    # "PATHWAY-CONTEXT" claims are explicitly about a DIFFERENT gene (e.g.
    # CXXC4, not USP34) cited only to show the pathway the candidate sits
    # inside is independently causal for tamoxifen resistance -- these must
    # never count toward "candidate has direct tamoxifen/ER+/endocrine
    # evidence" rollups, or they would misrepresent evidence about another
    # gene as evidence about the candidate itself.
    c["is_candidate_specific"] = "PATHWAY-CONTEXT" not in c["claim_id"]


def build_claim_evidence_table() -> pd.DataFrame:
    df = pd.DataFrame(ALL_CLAIMS)
    return df[CLAIM_COLUMNS]


def build_literature_evidence_table(claims: pd.DataFrame) -> pd.DataFrame:
    """Per-candidate, per-evidence-level rollup: how many claims, at what
    level, support vs contradict the project's working hypothesis for that
    candidate."""
    rows = []
    for candidate, grp in claims.groupby("candidate"):
        for level in sorted(grp["evidence_level"].unique()):
            sub = grp.loc[grp["evidence_level"] == level]
            rows.append(
                {
                    "candidate": candidate,
                    "evidence_level": level,
                    "n_claims": len(sub),
                    "n_supporting": int(sub["supports_project_hypothesis"].sum()),
                    "n_contradicting": int(sub["contradicts_project_hypothesis"].sum()),
                    "claim_ids": ";".join(sub["claim_id"]),
                    "breast_cancer_specific_count": int(sub["breast_cancer_specific"].sum()),
                    "tamoxifen_specific_count": int(sub["tamoxifen_specific"].sum()),
                }
            )
    return pd.DataFrame(rows).sort_values(["candidate", "evidence_level"])


def build_verified_references_table(claims: pd.DataFrame) -> pd.DataFrame:
    refs = claims.loc[
        claims["paper_title"] != "(no paper found -- explicit negative result)",
        ["candidate", "paper_title", "PMID", "DOI", "year", "journal", "review_or_primary", "fully_verified_primary_text", "claim_id"],
    ].copy()
    refs = refs.loc[~refs["paper_title"].str.startswith("(database", na=False)]
    refs = refs.drop_duplicates(subset=["paper_title", "PMID", "DOI"])
    refs = refs.rename(columns={"claim_id": "first_claim_id_citing_this_reference"})
    refs["verification_status"] = refs["fully_verified_primary_text"].map(
        {True: "full text fetched and quote-verified", False: "bibliographic identity verified (PMID/DOI/journal/year); narrative content from secondary aggregation, primary full text paywalled/blocked"}
    )
    return refs.sort_values(["candidate", "year"])


def run(out_tables_dir: Path = OUT_TABLES_DIR, out_reports_dir: Path = OUT_REPORTS_DIR) -> dict[str, pd.DataFrame]:
    out_tables_dir.mkdir(parents=True, exist_ok=True)
    out_reports_dir.mkdir(parents=True, exist_ok=True)

    claims = build_claim_evidence_table()
    claims.to_csv(out_tables_dir / "four_candidate_claim_evidence.tsv", sep="\t", index=False)
    logger.info("wrote four_candidate_claim_evidence.tsv (%d rows)", len(claims))

    lit_evidence = build_literature_evidence_table(claims)
    lit_evidence.to_csv(out_tables_dir / "four_candidate_literature_evidence.tsv", sep="\t", index=False)
    logger.info("wrote four_candidate_literature_evidence.tsv (%d rows)", len(lit_evidence))

    refs = build_verified_references_table(claims)
    refs.to_csv(out_reports_dir / "verified_references.tsv", sep="\t", index=False)
    logger.info("wrote verified_references.tsv (%d unique references)", len(refs))

    return {"claims": claims, "literature_evidence": lit_evidence, "references": refs}


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run()
