"""Curated interpretation content for the final USP34/VEZF1 GDSC
pharmacogenomics phase -- synthesizes (never re-derives) the real,
already-computed statistical results in
results/tables/final_pharmacogenomics/{USP34,VEZF1}_GDSC_drug_associations.tsv
and GDSC_top_associations.tsv. Every specific number quoted here was
computed by src/final_pharmacogenomics_build_tables.py from real GDSC
Release 8.5 + DepMap 26Q1 data, never hand-invented.

Absolute causality rule (enforced throughout): a GDSC expression-response
correlation is reported only as "X expression is associated with drug Y
response in GDSC breast cancer lines" -- never "X causes sensitivity",
never "Y inhibits X" (unless Y is a direct, annotated target of X, which
is not the case for any finding below), never "Y should be combined with
tamoxifen."
"""

from __future__ import annotations

INDIRECT_TARGETING_CATEGORIES = {
    "DIRECT_TARGET", "KNOWN_UPSTREAM_REGULATOR", "KNOWN_REQUIRED_PARTNER",
    "PATHWAY_LEVEL_CONNECTION", "PHARMACOGENOMIC_ASSOCIATION_ONLY", "NO_KNOWN_CONNECTION",
}

# Curated classification for the specific drugs/targets that reached
# FDR<0.05 or were otherwise flagged as top hits in the real analysis --
# judgment applied to REAL, already-computed findings, never used to
# select which findings to report.
DRUG_INDIRECT_TARGETING_ROWS = [
    dict(gene="USP34", drug_name="AZD7762", target="CHEK1, CHEK2", fdr=0.0078,
         classification="PATHWAY_LEVEL_CONNECTION",
         rationale="No direct biochemical or genetic evidence that USP34 regulates CHK1/CHK2 or vice versa was identified anywhere in this project's prior research passes (druggability_safety, lead_target_deep_dive). The connection is a cell-cycle/DNA-damage-response pathway-level co-occurrence. CHEK1 itself is directionally consistent but NOT FDR-significant in this project's own frozen Hany CRISPR screen (sensitising_KO direction, FDR=0.812 -- not significant); it IS nominally FDR-significant in this project's frozen GSE118713 transcriptomic tamoxifen-resistance dataset (FDR=0.0497, down_in_TAMR direction) -- a real, independently-generated echo worth noting, but this remains a correlational pharmacogenomic finding, not a demonstrated mechanistic link between USP34 and CHK1/CHK2."),
    dict(gene="USP34", drug_name="LDN-193189", target="BMP", fdr=0.0439,
         classification="PATHWAY_LEVEL_CONNECTION",
         rationale="Biologically coherent with USP34's own independently-established BMP2/Smad1/RUNX2 osteogenesis mechanism (Guo et al. 2018, PMID 30181118, from this project's druggability_safety phase) -- but that mechanism was demonstrated in bone/mesenchymal biology, not breast cancer, and this GDSC correlation does not itself demonstrate that USP34 regulates BMP signaling in breast cancer cells specifically. Also did not replicate in the small (N=14) ER+/luminal subset (p=0.37), so this signal should not be over-weighted."),
    dict(gene="USP34", drug_name="FGFR_3831 / AZD1332 / AZD6738 / JAK1_3715 / Sphingosine Kinase 1 Inhibitor II / Ara-G", target="FGFR1-4 / NTRK1-3 / ATR / JAK1 / SPHK1 / anti-metabolite", fdr="0.024-0.045",
         classification="PHARMACOGENOMIC_ASSOCIATION_ONLY",
         rationale="No known direct, upstream-regulatory, or pathway-level mechanistic connection to USP34 biology (Wnt/AXIN1, bone/osteogenesis, or CNS/developmental) was identified for these targets. Reported as real, FDR-significant correlations with no current mechanistic explanation -- consistent with a real, if currently unexplained, biological association, not promoted beyond that. Several of these targets' own genes (FGFR2, FGFR4, JAK1) independently show FDR-significant differential expression in this project's own frozen GSE118713/GSE240112 resistance datasets -- see GDSC_project_crosscheck.tsv -- a correlational echo across independent datasets, not a demonstrated mechanistic link to USP34."),
    dict(gene="VEZF1", drug_name="(none reach FDR<0.05)", target="n/a", fdr="n/a (best FDR=0.132, Paclitaxel)",
         classification="NO_KNOWN_CONNECTION",
         rationale="No VEZF1 drug association survives FDR correction in the full breast-line analysis. The single best nominal hit (Paclitaxel, microtubule/mitosis, p=0.0005, FDR=0.132) does not reach significance and has no known mechanistic link to VEZF1's documented vascular/transcription-factor biology."),
    dict(gene="VEZF1", drug_name="TEAD/Hippo pathway compounds", target="n/a", fdr="n/a -- NO TEAD/YAP/Hippo-pathway compound exists in GDSC at all",
         classification="PHARMACOGENOMIC_ASSOCIATION_ONLY",
         rationale="GDSC Release 8.5's compound list was searched directly (TARGET_PATHWAY containing 'Hippo', TARGET containing 'TEAD/YAP/TAZ/WWTR1', DRUG_NAME containing 'TEAD'/'verteporfin') and returns ZERO compounds -- there is no TEAD/Hippo drug in GDSC to test the VEZF1-TEAD1 hypothesis pharmacogenomically at all. This is an honest absence, not a negative finding about the hypothesis itself. GDSC's WNT-signaling compounds (including tankyrase inhibitors XAV939/MN-64/WIKI4/AZ6102, mechanistically adjacent to USP34's own AXIN1 axis) were checked for both genes and show no FDR-significant association for either USP34 or VEZF1 (best USP34 WNT-pathway p=0.0076 for TWS119/GSK3, FDR=0.138, not significant)."),
]

FINAL_INTERPRETATION_ROWS = [
    dict(candidate="USP34", question="Does USP34 expression associate with any existing drug response?",
         answer="Yes. 9 of 1,278 drug-tests (breast lines, N>=15, both metrics, both screens) reach FDR<0.05, all in GDSC1, all with the SAME direction (higher USP34 expression -> more sensitive, i.e. lower LN_IC50/AUC)."),
    dict(candidate="USP34", question="Any FDR-significant hits?",
         answer="Yes, 9: AZD7762 (CHEK1/CHEK2, FDR=0.0078-0.0084, both LN_IC50 and AUC), FGFR_3831 (FGFR1-4, FDR=0.0248), Sphingosine Kinase 1 Inhibitor II (FDR=0.0248), JAK1_3715 (FDR=0.0248), AZD1332 (NTRK1-3, FDR=0.0248), AZD6738 (ATR, FDR=0.0248), LDN-193189 (BMP, FDR=0.0439), Ara-G (anti-metabolite, FDR=0.0458). Full detail: USP34_GDSC_drug_associations.tsv / GDSC_top_associations.tsv."),
    dict(candidate="USP34", question="Any consistent cross-screen signal?",
         answer="Partial. AZD7762 (CHEK1/CHEK2) is the strongest example: FDR-significant in GDSC1 under DRUG_ID 1402 only (the other GDSC1 re-screening batch, DRUG_ID 1022, is directionally consistent but NOT FDR-significant: FDR=0.174-0.350) AND shows a nominally significant, SAME-direction signal in GDSC2 (rho=-0.372, p=0.011) that does not survive its own within-screen FDR correction. This is the closest thing to cross-screen consistency found; no other hit drug reaches even nominal significance in the other screen."),
    dict(candidate="USP34", question="Most biologically coherent pathway?",
         answer="Two candidates, both PATHWAY_LEVEL_CONNECTION (not DIRECT_TARGET): (1) CHEK1/CHEK2 (AZD7762) -- CHEK1 itself is directionally consistent (sensitising_KO) but not FDR-significant in this project's own frozen Hany CRISPR screen (FDR=0.812); it IS nominally FDR-significant in this project's frozen GSE118713 transcriptomic tamoxifen-resistance dataset (FDR=0.0497, down_in_TAMR direction), a real cross-dataset echo. (2) BMP (LDN-193189) -- directly coherent with USP34's own established BMP2/Smad1/RUNX2 osteogenesis mechanism from the druggability_safety phase, though this signal did not replicate in the small ER+/luminal subset."),
    dict(candidate="USP34", question="Any realistic hypothesis worth future testing?",
         answer="Yes, as a hypothesis, not a validated combination: whether USP34 perturbation sensitises resistant ER+ cells to CHK1/CHK2 inhibition (a DNA-damage-response combination hypothesis) is worth a dedicated follow-up given the cross-dataset (GDSC correlation + GSE118713 nominal down-in-TAMR signal for CHEK1) convergence on the CHEK1/CHEK2 axis -- but this is a NEW hypothesis about a CHK1/CHK2 co-target, not evidence bearing on the tamoxifen-sensitisation hypothesis itself, and is not incorporated into EXP-1's core design in this phase."),
    dict(candidate="USP34", question="Does GDSC alter USP34 = LEAD?",
         answer="No. Nothing in this analysis strengthens or weakens the frozen USP34=LEAD conclusion -- tamoxifen and fulvestrant themselves show no significant USP34 association in GDSC breast lines (consistent with, not contradictory to, the Hany finding, since GDSC tests baseline/unperturbed expression-vs-response correlation across a cross-sectional cell-line panel, a different question from Hany's CRISPR-perturbation-plus-tamoxifen design). The new CHEK1/CHEK2 finding is an interesting adjacent hypothesis, not a change to the frozen lead conclusion."),
    dict(candidate="VEZF1", question="Does VEZF1 expression associate with any drug response?",
         answer="No FDR-significant association with any of the 1,278 drug-tests performed (breast lines, N>=15). Best nominal signal: Paclitaxel (microtubule/mitosis), p=0.0005, FDR=0.132 -- does not survive correction."),
    dict(candidate="VEZF1", question="Any FDR-significant hits?", answer="None."),
    dict(candidate="VEZF1", question="Any consistent cross-screen signal?", answer="No -- with zero FDR-significant hits in either screen, there is nothing to check for cross-screen consistency."),
    dict(candidate="VEZF1", question="Any Hippo/TEAD relevance?",
         answer="Cannot be tested: GDSC Release 8.5 contains ZERO TEAD/YAP/Hippo-pathway compounds of any kind (confirmed by direct compound-list search). This is an absence of the necessary chemical matter in GDSC, not a negative pharmacogenomic finding about the VEZF1-TEAD1 hypothesis."),
    dict(candidate="VEZF1", question="Any credible indirect pharmacological opportunity?",
         answer="None identified in this GDSC analysis. No drug-response association reaches significance, and the one mechanistically-relevant drug class this project has flagged for VEZF1 (TEAD inhibitors) is entirely absent from GDSC."),
    dict(candidate="VEZF1", question="Does GDSC materially increase VEZF1 translational feasibility?",
         answer="No. This is a negative result, reported as such rather than reframed -- consistent with VEZF1's already-documented poor direct/indirect druggability from the prior lead_target_deep_dive phase."),
    dict(candidate="VEZF1", question="Does GDSC alter VEZF1 = BACKUP?",
         answer="No. The frozen VEZF1=SECOND/BACKUP conclusion is unchanged; this analysis adds a further negative data point (no pharmacogenomic drug-response signal in GDSC) consistent with, not contradicting, VEZF1's already-established poor current druggability."),
]
