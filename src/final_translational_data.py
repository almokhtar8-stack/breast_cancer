"""Curated source data for the final USP34/VEZF1 translational + structure
phase. USP34 is the frozen lead target, VEZF1 the frozen second/backup --
neither ranking is altered here. No frozen upstream evidence (discovery,
TCGA, DepMap, druggability_safety, lead_target_deep_dive) is touched or
recomputed; this phase only designs forward-looking experiments and adds
one genuinely new analysis -- a real, locally-run structural/pocket
analysis of the two experimental USP34 structures.

Structural/pocket provenance: PDB files 7W3R and 7W3U were downloaded
directly from RCSB (files.rcsb.org) this session. All resolution, chain,
residue-range, missing-residue, HETATM/ligand, and LINK (covalent bond)
facts below were extracted directly by parsing those downloaded
coordinate files (not recalled from memory or re-typed from the earlier
PMID 35588869 citation). Catalytic-residue distances were computed
directly from the atomic coordinates (Python, Euclidean distance, no
external geometry library). Pocket detection was performed with fpocket
4.2.3 (conda-forge), installed and run locally this session on both
structures -- this is a real, reproducible, open-source local tool run
per the user's explicit instruction to prefer such tools "where
available." No pocket is called "druggable" without the fpocket
druggability-score qualification attached; a score is a geometric/
physicochemical heuristic, not proof of an achievable drug.
"""

from __future__ import annotations

TARGETS = ["USP34", "VEZF1"]

# ---------------------------------------------------------------------------
# 1. LOCKED FINAL BIOLOGICAL QUESTIONS
# ---------------------------------------------------------------------------
LOCKED_QUESTIONS = {
    "USP34": (
        "Does reducing USP34 restore or enhance tamoxifen response in "
        "acquired tamoxifen-resistant ER+ breast cancer cells without "
        "simply causing nonspecific baseline toxicity? -- NOT described as "
        "proven; this is the question the experiment in this report is "
        "designed to answer, not an established finding."
    ),
    "VEZF1": (
        "Does reducing VEZF1 (A) impair resistant ER+ cancer-cell fitness "
        "on its own, AND/OR (B) further restore/enhance tamoxifen "
        "sensitivity? -- NOT described as proven; both sub-questions are "
        "explicitly separated because the DepMap 26Q1 baseline-dependency "
        "signal (27.3%) and the Hany tamoxifen-sensitisation signal "
        "(FDR=0.037) are two distinct, independently-generated pieces of "
        "evidence that this experiment is designed to test jointly, not "
        "an already-confirmed dual mechanism."
    ),
}

# ---------------------------------------------------------------------------
# 2/4/6. EXPERIMENTAL DESIGN (EXP-1 USP34 main, EXP-3 VEZF1 main, EXP-5 TEAD1)
# ---------------------------------------------------------------------------
EXPERIMENTAL_DESIGN_ROWS = [
    dict(
        experiment_id="EXP-1", target="USP34",
        model_system="Acquired tamoxifen-resistant ER+ breast cancer cell line (matching this project's own TAMR-model framing, e.g. an MCF-7- or T47D-derived long-term-4-OHT-selected resistant line) WITH a matched parental ER+ line as comparator -- both arms of the resistant-vs-parental comparison are required to distinguish a resistance-specific effect from a generic ER+ cell effect.",
        arm_1="Control (vehicle)", arm_2="Tamoxifen / 4-OHT (single agent, resistant-line-appropriate dose)",
        arm_3="USP34 perturbation (see perturbation_strategy)", arm_4="USP34 perturbation + tamoxifen/4-OHT",
        perturbation_strategy_recommended="CRISPR knockout (KO) via lentiviral sgRNA/Cas9 for the FIRST validation pass, NOT CRISPRi, siRNA, or an inhibitor.",
        perturbation_rationale=(
            "No validated USP34-selective chemical probe or inhibitor currently exists (confirmed across three independent research passes in this project) -- "
            "so a pharmacological arm is not available for a first validation and must not be simulated by a proxy. Between the remaining genetic options: "
            "siRNA/shRNA give only partial, transient, and often inconsistent knockdown for a large (3546-aa) multi-domain protein, which would confound "
            "interpretation of a negative result (is it truly no phenotype, or insufficient knockdown?). CRISPRi (promoter-region dCas9-KRAB) reduces transcription "
            "but a preexisting protein pool of a stable, abundant enzyme like a DUB can persist for the duration of a short assay, again risking a false-negative "
            "reading. CRISPR KO gives a clean, complete, stable loss-of-function -- the most INFORMATIVE first-pass answer to 'does USP34 loss have this "
            "phenotype at all', even though it will not by itself distinguish complete genetic loss from the partial pharmacological inhibition a future drug "
            "would achieve -- that distinction is explicitly deferred to a later, dose-titratable, inducible-degron or analog-sensitive-allele follow-up once a "
            "positive CRISPR-KO signal justifies the added complexity. This is stated plainly as a limitation of the first-pass design, not a proof of drug-like relevance. "
            "COMPLETE KO IS EXPLICITLY NOT EQUIVALENT TO PARTIAL PHARMACOLOGICAL INHIBITION -- a later titratable CRISPRi/siRNA/chemical follow-up would be required "
            "to model what a partial inhibitor would actually do."
        ),
        emt_stemness_counter_evidence=(
            "PMID 28499884 (Wang et al., Cell Signal 2017;36:230-239, 'Inhibition of ubiquitin-specific protease 34 (USP34) induces epithelial-mesenchymal transition "
            "and promotes stemness in mammary epithelial cells' -- verified directly via NCBI eutils this session, exact abstract text checked, not summarized from "
            "memory). MODEL: NMuMG cells (normal, non-transformed mouse mammary epithelial cell line -- NOT a cancer line, NOT ER+, NOT a tamoxifen-resistance model); "
            "correlative expression comparison also used MDA-MB-231 and 4T1 (both ER-negative, mesenchymal-like/triple-negative lines); in vivo mammary-fat-pad "
            "transplantation in mice. PERTURBATION: USP34 knockdown (KD, not KO). FINDINGS (verbatim from the abstract): USP34 KD in NMuMG cells induced EMT, "
            "evidenced by upregulation of N-cadherin, phospho-Smad3, Snail, and active-beta-catenin, and downregulation of Axin1 and E-cadherin; USP34 KD also "
            "produced invasive behavior and enhanced mammosphere-forming ability with upregulated Nanog/Oct4/Sox2; USP34-KD cells transplanted into cleared mammary "
            "fat pads reconstituted the mammary gland with ductal-tree development within 3 months. INTERPRETATION (explicitly, per the user's required framing): "
            "this does NOT invalidate the USP34 tamoxifen-sensitisation hypothesis. It is counter-evidence / a potential context-specific liability, because this "
            "study is not equivalent to acquired-tamoxifen-resistant-ER+-breast-cancer + USP34-perturbation + tamoxifen -- it used a normal (non-transformed, "
            "non-cancer, non-ER+) mammary epithelial line, not a resistant cancer model, and never combined USP34 loss with tamoxifen. Prior mammary-epithelial work "
            "suggests USP34 loss can promote EMT/stem-like features in SOME contexts, motivating explicit EMT/stemness monitoring during tamoxifen-resensitisation "
            "experiments -- it does not support the stronger claim that 'USP34 inhibition promotes aggressive breast cancer,' which this report does not make."
        ),
        primary_readouts=(
            "Viability (dose-response over 5-7 days), apoptosis (Annexin V/cleaved caspase-3), clonogenic survival (10-14 day colony assay). "
            "INTERACTION FRAMEWORK (required, not optional -- sensitisation is explicitly NOT defined merely as 'combination > either single arm'): "
            "generate a FULL tamoxifen/4-OHT dose-response curve under (i) non-targeting-control perturbation and (ii) USP34-KO, and compare the two curves directly "
            "(EC50/IC50 shift, not just an endpoint viability comparison at one dose); apply a predefined interaction/additivity model (Bliss independence and/or the "
            "Chou-Talalay combination index) to the full dose-response matrix, not a single-dose 'combination vs single arm' comparison, before any sensitisation claim is made."
        ),
        mechanistic_readouts=(
            "USP34 protein loss (western blot, confirming KO); AXIN1 ubiquitination status and active (non-phospho) beta-catenin levels (the one experimentally-validated "
            "USP34 substrate axis, PMID 21383061) -- included ONLY as a mechanistic readout tied to a real prior finding, not assumed relevant to tamoxifen resistance "
            "without this data; ESR1 target-gene qPCR panel (TFF1, PGR, GREB1) to test whether USP34 loss perturbs ER transcriptional output specifically or acts "
            "through a parallel pathway. EMT/STEMNESS MONITORING (added this phase, directly motivated by PMID 28499884, using the SAME marker panel that paper "
            "reported): E-cadherin/CDH1, N-cadherin/CDH2, SNAI1/Snail, AXIN1, active (non-phospho) beta-catenin -- note AXIN1 and active-beta-catenin are read out "
            "for BOTH the Wnt-mechanism question above AND the EMT-liability question, since PMID 28499884 implicates the same axis. Optional functional assays if "
            "feasible: migration/invasion (Transwell/Boyden chamber) and mammosphere formation (matching PMID 28499884's own stemness readout)."
        ),
        outcome_categories=(
            "IDEAL: USP34 perturbation increases tamoxifen response WITHOUT meaningful EMT/stemness induction (EMT/stemness markers and mammosphere formation "
            "unchanged or minimally changed vs non-targeting control). "
            "CONCERNING: tamoxifen sensitisation occurs BUT substantial EMT/stemness/invasion also increases (E-cadherin down, N-cadherin/Snail/active-beta-catenin "
            "up, enhanced mammosphere formation) -- this would not disprove the sensitisation finding but would flag a real, monitorable liability requiring further "
            "characterization before any combination strategy is pursued. "
            "PURE GENERAL TOXICITY: USP34 perturbation strongly kills cells regardless of tamoxifen (arm 3 approx. equals arm 4, both far below arm 1) -- weakens "
            "the sensitisation hypothesis regardless of the EMT read. "
            "NEGATIVE: no meaningful tamoxifen sensitisation is observed (arm 4 not meaningfully different from arm 2) -- the core hypothesis is not supported in "
            "this model, independent of the EMT question."
        ),
        replicates_and_controls="Non-targeting sgRNA control (matched lentiviral delivery); >=3 independent biological replicates per arm; parental-vs-resistant comparison run in parallel with identical reagents/timing; rescue experiment (re-expression of sgRNA-resistant USP34 cDNA, ideally including a catalytic-dead C1903S point mutant to test whether any phenotype requires catalytic activity specifically) recommended if an initial phenotype is observed, to rule out an off-target CRISPR effect.",
    ),
    dict(
        experiment_id="EXP-3", target="VEZF1",
        model_system="Same acquired tamoxifen-resistant ER+ line (+ matched parental line) as EXP-1, for direct cross-target comparability.",
        arm_1="Control (vehicle)", arm_2="Tamoxifen / 4-OHT (single agent)",
        arm_3="VEZF1 suppression (see perturbation_strategy)", arm_4="VEZF1 suppression + tamoxifen/4-OHT",
        perturbation_strategy_recommended="CRISPRi (dCas9-KRAB targeted to the VEZF1 promoter/TSS) as the primary first-pass modality, with inducible CRISPR KO as a secondary/confirmatory option if CRISPRi signal is ambiguous.",
        perturbation_rationale=(
            "VEZF1 is a transcription factor with no validated direct pharmacology (confirmed across this project's prior phases: no experimental structure, only "
            "one weak 2018 tool-compound series). Because there is no drug to eventually match the perturbation to, the priority for a TF is INTERPRETABLE "
            "functional evidence rather than mimicking a future drug's mechanism. CRISPRi gives clean, tunable, reversible transcriptional silencing without the "
            "confound of a truncated/non-functional protein product that CRISPR KO of a DNA-binding protein can sometimes produce (a partial DBD-truncated "
            "product could retain dominant-negative or gain-of-function activity, complicating interpretation for a TF specifically, unlike for an enzyme like "
            "USP34 where a null allele is comparatively easy to interpret). siRNA/shRNA are secondary options if CRISPRi reagents are not readily available; "
            "they are workable for a TF but give a less complete and less durable knockdown than CRISPRi over the multi-day assay windows viability/clonogenic "
            "assays require."
        ),
        primary_readouts="Viability, apoptosis, clonogenic survival, tamoxifen dose-response shift -- identical panel to EXP-1 for direct comparability. Explicitly analyze arm 3 vs arm 1 SEPARATELY from arm 4 vs (arm 2 + arm 3 effects) to distinguish DIRECT CANCER DEPENDENCY (does VEZF1 suppression alone impair resistant cells, consistent with the real 27.3% DepMap ER+/luminal baseline-dependency signal) from TAMOXIFEN SENSITISATION (does the combination outperform what arms 2 and 3 alone would predict, consistent with the real Hany FDR=0.037 signal) -- these are two distinct questions and must not be collapsed into one 'VEZF1 matters' readout.",
        mechanistic_readouts="ER transcriptional response (TFF1/PGR/GREB1 qPCR, as in EXP-1, to compare mechanistic overlap or divergence between the two targets); a VEZF1 target-program readout IS defensible and should be included -- VEZF1's own documented target genes (VEGFR2/KDR, TIMP3, MMP2) from this project's prior literature-mechanism phase, read out by qPCR, to confirm CRISPRi/knockdown is functionally silencing the VEZF1 transcriptional program, not merely reducing VEZF1 mRNA without functional consequence.",
        replicates_and_controls="Non-targeting sgRNA/scrambled siRNA control; >=3 independent biological replicates; parental-vs-resistant comparison in parallel; rescue via sgRNA-resistant VEZF1 cDNA re-expression recommended if a phenotype is observed.",
    ),
    dict(
        experiment_id="EXP-5", target="TEAD1 (VEZF1 indirect-target hypothesis test -- NOT a therapeutic experiment)",
        model_system="Same ER+ tamoxifen-resistant line, used here purely for TARGET VALIDATION, not efficacy testing.",
        arm_1="Control (vehicle)", arm_2="TEAD inhibitor (pan-TEAD auto-palmitoylation-pocket compound, e.g. the VT3989 chemotype class -- no TEAD1-selective compound remains in active clinical development, so a pan-TEAD tool compound is the only available chemical route)",
        arm_3="VEZF1 suppression (same CRISPRi reagent as EXP-3, as a reference/positive-control arm for what a real VEZF1 reduction looks like on the readouts below)", arm_4="TEAD inhibitor + tamoxifen (ONLY if arm 2 shows a genuine VEZF1 effect in the readouts below -- otherwise do not proceed to this arm)",
        perturbation_strategy_recommended="Pharmacological pan-TEAD inhibitor (arm 2) as the primary perturbation, since the entire point of this experiment is to test whether an available/developable DRUG-LIKE intervention on TEAD affects VEZF1 -- a genetic TEAD1 knockdown arm could be added as a secondary, TEAD1-specific (vs pan-TEAD) mechanistic follow-up if the pharmacological arm is positive, to determine whether the effect (if any) is TEAD1-specific or a pan-TEAD/YAP-axis effect.",
        primary_readouts="This experiment's readouts ARE the primary question, not a downstream efficacy measure: VEZF1 RNA (qPCR), VEZF1 protein (western blot), and the VEZF1-dependent transcriptional program (VEGFR2/TIMP3/MMP2 panel, same as EXP-3's mechanistic readout, for direct comparability). A MANDATORY positive-control readout is also required: a canonical TEAD/YAP target-engagement gene panel (e.g. CTGF/CCN2, CYR61/CCN1) to confirm the TEAD inhibitor is pharmacologically active in this cell line at all -- without this control, a negative VEZF1 result would be uninterpretable (target not engaged, vs. target engaged but VEZF1 unaffected).",
        mechanistic_readouts="n/a -- this experiment's readouts are themselves the mechanistic test.",
        replicates_and_controls=">=3 independent biological replicates; the VEZF1-suppression arm (arm 3) serves as an internal positive-control reference for what a genuine VEZF1-program reduction looks like on the same readout panel, so a TEAD-inhibitor result can be compared against a known-positive benchmark rather than judged in isolation.",
        decision_rule="IF no VEZF1 RNA/protein/target-program change occurs despite confirmed TEAD/YAP target engagement: REJECT TEAD1 as an indirect VEZF1-targeting strategy, explicitly and permanently for this specific hypothesis (do not retest with a different pan-TEAD compound without new biological rationale). IF VEZF1 function decreases: proceed to arm 4 (TEAD inhibitor +/- tamoxifen) as a follow-up efficacy experiment. TEAD1 must NOT be incorporated into the main USP34/VEZF1 therapeutic model (EXP-1/EXP-3 framing) before this experiment produces a positive result.",
    ),
]

# ---------------------------------------------------------------------------
# 3/5. NORMAL-CELL COMPARATORS (EXP-2 USP34, EXP-4 VEZF1)
# ---------------------------------------------------------------------------
NORMAL_CELL_COMPARATOR_ROWS = [
    dict(
        experiment_id="EXP-2A", target="USP34",
        comparator_concept="A. LINEAGE / CANCER-SELECTIVITY COMPARATOR",
        comparator_cell_type="Normal human mammary epithelial cells (e.g. primary human mammary epithelial cells, HMECs, if experimentally feasible; MCF10A as a widely-used non-transformed line alternative), given the SAME USP34 perturbation (CRISPR KO) as EXP-1.",
        rationale=(
            "This is the direct tissue-of-origin comparator for the cancer model, and its purpose is different from EXP-2B's: it asks whether COMPARABLE "
            "USP34 perturbation disproportionately affects malignant ER+ cells relative to normal mammary cells (a cancer-selectivity/therapeutic-index "
            "question), not whether USP34 is required for a specific normal-tissue differentiation program. It is also the most directly relevant normal "
            "cell type given PMID 28499884's EMT/stemness counter-evidence, which was generated in a mammary epithelial model (NMuMG, mouse) -- this human "
            "comparator directly extends that context and lets the same EMT/stemness marker panel from EXP-1 (E-cadherin/CDH1, N-cadherin/CDH2, SNAI1, "
            "AXIN1, active beta-catenin, mammosphere formation) be checked in a NORMAL mammary background, not only in the resistant-cancer background."
        ),
        readouts="Viability/proliferation under USP34 perturbation relative to the ER+ resistant/parental lines in EXP-1 (a direct selectivity ratio, not an isolated readout); the SAME EMT/stemness marker panel as EXP-1 (E-cadherin/CDH1, N-cadherin/CDH2, SNAI1, AXIN1, active beta-catenin, optional mammosphere formation) to test whether PMID 28499884's normal-mammary-epithelial EMT/stemness finding reproduces in this human system under the SAME perturbation used for the cancer experiment.",
        minimal_set_justification="Promoted to a co-primary comparator (not a secondary/optional add-on) in this revision, specifically because of the PMID 28499884 counter-evidence -- a normal mammary epithelial readout is now directly informative for BOTH the lineage-selectivity question AND the EMT/stemness-liability question, not merely a generic epithelial control.",
        safety_disclaimer="This is a PRELIMINARY THERAPEUTIC-WINDOW / selectivity experiment, not a comprehensive safety study -- a favorable selectivity ratio here must never be described as proof that USP34 inhibition is safe in normal mammary tissue in vivo.",
    ),
    dict(
        experiment_id="EXP-2B", target="USP34",
        comparator_concept="B. KNOWN-LIABILITY COMPARATOR",
        comparator_cell_type="Primary human bone-marrow-derived mesenchymal stem cells (MSCs) undergoing osteogenic differentiation induction, given the SAME USP34 perturbation (CRISPR KO) as EXP-1.",
        rationale=(
            "Of USP34's candidate liability signals (high skeletal-muscle expression with NO demonstrated functional muscle liability identified in any "
            "search performed across this project's phases; and real, causal, replicated bone/osteogenic-mesenchymal biology across THREE independent findings "
            "-- BMP2/osteogenesis via Smad1/RUNX2 stabilization, tooth-root/NFIC stabilization, and a 2026 human haploinsufficiency syndrome with limb "
            "anomalies), only the bone/osteogenic axis has ANY demonstrated functional requirement for USP34 in a normal cell type. Skeletal muscle does NOT "
            "need to be the primary comparator: expression is high, but no functional skeletal-muscle liability was identified in this project's repeated "
            "searches. This MSC assay directly extends the exact published assay (PMID 30181118) that established USP34's bone-formation role, giving a "
            "literature-anchored positive-control expectation to compare against."
        ),
        readouts="Viability; osteogenic differentiation (RUNX2 expression, alkaline phosphatase activity, Alizarin Red mineralization at a standard functional-differentiation endpoint, e.g. day 14-21); relative sensitivity -- does the level of USP34 loss needed to produce a measurable cancer-cell effect in EXP-1 also measurably impair normal osteogenic differentiation at the same perturbation level, or is there a detectable window between the two?",
        minimal_set_justification="MSC/osteogenic remains the single highest-priority, evidence-anchored KNOWN-LIABILITY comparator. A skeletal-muscle myoblast/myotube comparator is explicitly NOT recommended as part of the minimum set, given the absence of any functional muscle evidence to test.",
        safety_disclaimer="This is a PRELIMINARY THERAPEUTIC-WINDOW experiment, designed to detect the most plausible already-identified liability signal -- it is explicitly NOT a comprehensive safety study and a negative/favorable result here must never be described as proof that USP34 inhibition is safe.",
    ),
    dict(
        experiment_id="EXP-4", target="VEZF1",
        comparator_cell_type="Primary human vascular endothelial cells (e.g. HUVECs or a more disease-relevant microvascular endothelial line) as the primary comparator, with iPSC-derived cardiomyocytes as a secondary comparator if resources allow.",
        rationale=(
            "VEZF1's own namesake biology and its best-established, most-replicated developmental role is vascular/endothelial (founding papers, embryonic "
            "knockout phenotype, target genes VEGFR2/TIMP3/MMP2) -- endothelial cells are the single most directly literature-anchored normal-cell comparator "
            "available. The newer, real cardiac-muscle finding (PMID 31911272, zebrafish knockdown impairing postnatal cardiac contractile function via "
            "Myh7/TEAD-1) is genuinely important and is why iPSC-cardiomyocytes are recommended as a secondary comparator rather than omitted -- but it is a "
            "single non-human (zebrafish), partial-knockdown finding, and extending it to a human primary/iPSC cardiomyocyte system is the appropriate "
            "translational next step, not an assumption that the zebrafish phenotype will reproduce. Normal mammary epithelial cells (MCF10A) are a third, "
            "lower-priority option for the same generic-epithelial-effect control rationale as in EXP-2, but are not prioritized above the two tissue-specific "
            "options given VEZF1's much stronger vascular/cardiac biological anchoring versus generic epithelial biology."
        ),
        readouts="Endothelial: viability/proliferation, tube-formation/angiogenic-capacity assay (directly extending the assay used in this project's own prior literature review of VEZF1's endothelial biology), VEGFR2/TIMP3/MMP2 target-gene expression under VEZF1 perturbation. Cardiomyocyte (if included): viability, and where feasible a contractile-function readout (e.g. video-based contractility/Ca2+-imaging platform) directly analogous to the zebrafish beta-adrenergic-response readout in PMID 31911272, plus Myh7/beta-MHC expression.",
        minimal_set_justification="Endothelial cells alone are the recommended true minimum (highest literature anchoring, most tractable primary-cell system); iPSC-cardiomyocytes are recommended as an add-on specifically because they test the newest and most functionally causal liability finding in this whole project (a partial-knockdown, postnatal, contractile-function phenotype) rather than being included reflexively.",
        safety_disclaimer="Animal (zebrafish) cardiac knockdown findings are NOT equated with expected adult human drug toxicity anywhere in this design -- this experiment exists specifically to generate the missing human-cell data point, not to confirm an assumption. This is a preliminary therapeutic-window experiment, not a safety claim.",
    ),
]

# ---------------------------------------------------------------------------
# 7. USP34 STRUCTURAL INVENTORY -- extracted directly from downloaded PDB
# coordinate files this session (RCSB files.rcsb.org/download/{ID}.pdb),
# not retyped from a paper's claims.
# ---------------------------------------------------------------------------
STRUCTURE_INVENTORY_ROWS = [
    dict(
        pdb_id="7W3R", title="USP34 catalytic domain (apo)", method="X-RAY DIFFRACTION", resolution_angstrom=1.92,
        chains="A, B (2 copies of USP34 catalytic domain in the asymmetric unit)",
        construct_residue_range="1892-2261 (chain A, as modeled; full catalytic-domain construct per RCSB header)",
        missing_disordered_residues="105 residues listed in REMARK 465 (not modeled/disordered) within this range -- a substantial fraction, consistent with a flexible, only partially-ordered catalytic-domain construct in the unliganded state",
        state="APO (no ubiquitin or ubiquitin-mimetic bound)",
        catalytic_residues_resolved="Cys1903 and His2164 both directly confirmed present/resolved (verified by direct coordinate parsing this session)",
        bound_ligands="ZN (structural zinc ion, chain-A-coordinated by Cys2018/His2020/Cys2062/Cys2065 per LINK records -- a distinct site from the catalytic Cys1903/His2164 dyad, NOT previously highlighted in this project's prior structural summaries); ordered waters",
        conformational_state_note="Measured directly from coordinates this session, chain A: Cys1903(SG)-His2164(ND1) distance = 3.94 A. [Codex-review correction: chain B of this same apo structure measures 4.98 A at the identical atom pair -- a real, disclosed ~1 A difference between the two crystallographically independent copies in this structure. The chain-A value is reported as the headline figure because it is the copy compared against 7W3U chain A below, not because it is representative of the whole structure; both values are given so the real copy-to-copy heterogeneity is not hidden.]",
        source="downloaded this session from files.rcsb.org/download/7W3R.pdb; original structure paper PMID 35588869",
    ),
    dict(
        pdb_id="7W3U", title="USP34 catalytic domain in complex with a ubiquitin-propargylamide (UbPA) activity-based probe", method="X-RAY DIFFRACTION", resolution_angstrom=3.13,
        chains="A, B, C (3 copies of USP34 catalytic domain); D, E, F (3 copies of the covalently-linked ubiquitin probe)",
        construct_residue_range="1892-2268 (chain A, as modeled); ubiquitin probe chains span residues 1-75 of ubiquitin's canonical 76 -- the native C-terminal Gly76 position is where the propargylamide probe chemistry is covalently attached, so 'residues 1-75 + probe adduct' rather than an unmodified 'full-length ubiquitin' is the technically precise description",
        missing_disordered_residues="141 residues listed in REMARK 465 within the USP34 chains -- more disorder than the apo structure, plausibly reflecting the lower (3.13 A) resolution rather than a genuine increase in flexibility upon binding; this caveat is stated explicitly rather than over-interpreted",
        state="COVALENT PROBE-BOUND (ubiquitin C-terminus covalently trapped at the catalytic cysteine)",
        catalytic_residues_resolved="Cys1903 and His2164 both directly confirmed present/resolved",
        bound_ligands="ZN (same structural zinc site as 7W3R, Cys2018/His2020/Cys2062/Cys2065, present in all three USP34 chains); AYE (prop-2-en-1-amine / allylamine, the propargylamide-probe-derived warhead remnant) -- LINK record directly confirms a covalent bond, SG of Cys1903 to C2 of AYE, bond distance 1.59-2.48 A across the three copies in the asymmetric unit -- this is DIRECT crystallographic proof that Cys1903 is the catalytic nucleophile and is covalently reactive/modifiable, not an inference from sequence alone",
        conformational_state_note="Measured directly from coordinates this session, chain A: Cys1903(SG)-His2164(ND1) distance = 3.37 A, a measurable ~0.6 A tightening relative to the SAME chain (A) of the apo structure (3.94 A) -- a real, independently-computed numerical corroboration of the literature's qualitative 'apo-inactive, realigns on ubiquitin engagement' description (PMID 35588869), not merely a repeated citation of that claim. [Codex-review correction, full transparency across all 3 copies in this asymmetric unit: chain B = 3.95 A, chain C = 3.10 A -- i.e. two of the three probe-bound copies (A and C) are tighter than the apo chain-A value, while chain B is essentially unchanged from apo. The chain-A-to-chain-A comparison is a real, valid single-copy observation; it is NOT claimed to be a uniform, structure-wide effect, and the raw heterogeneity is reported here rather than smoothed over.]",
        source="downloaded this session from files.rcsb.org/download/7W3U.pdb; original structure paper PMID 35588869",
    ),
    dict(
        pdb_id="(search performed, none found)", title="Additional experimentally-solved USP34 structures beyond 7W3R/7W3U", method="n/a", resolution_angstrom=None,
        chains="n/a", construct_residue_range="n/a",
        missing_disordered_residues="n/a", state="n/a", catalytic_residues_resolved="n/a",
        bound_ligands="n/a", conformational_state_note="No additional experimentally-solved USP34 structure was identified in this session's RCSB search (repeat of the prior phases' search, reconfirmed) -- 7W3R and 7W3U remain the only two, both catalytic-domain-only (covering ~12% of the 3546-aa full-length protein). No AlphaFold model exists for the full-length protein (API returns 404, reconfirmed in an earlier phase).",
        source="RCSB Search API, this session",
    ),
]

# ---------------------------------------------------------------------------
# 8. USP34 POCKET ANALYSIS -- real fpocket 4.2.3 output, run locally this
# session on both downloaded structures. Druggability scores are
# fpocket's own 0-1 heuristic (not a probability of clinical success);
# every score is reported with this qualification attached, per instruction.
# ---------------------------------------------------------------------------
POCKET_ANALYSIS_ROWS = [
    dict(
        structure="7W3R (apo)", pocket_label="Pocket 6 (top-ranked by fpocket druggability score among 54 total pockets detected)",
        fpocket_druggability_score=0.845, volume_A3=1471.4, num_alpha_spheres=136,
        key_lining_residues="Cys1903, His2164 (both catalytic residues), Gln1976 (the residue found this session to be spatially closest, 3.50 A, to His2164 among all Asp/Asn/Glu/Gln side chains checked -- a geometric observation from real coordinates, NOT a literature-confirmed catalytic-triad assignment, reported with that explicit caveat), plus Glu1917/Gln1920/Lys1926/Ser1994-Thr2000/Arg2082/Asp2161/Tyr2165-2166/Asp2188/Asp2252",
        location_description="This IS the catalytic cleft / extended ubiquitin-binding groove -- fpocket's single highest-scoring pocket in the apo structure directly encompasses the catalytic dyad.",
        interpretation_caveat="A 0.845 fpocket druggability score is a favorable GEOMETRIC/PHYSICOCHEMICAL heuristic, not proof a drug-like molecule can be found -- qualified explicitly. The very large volume (1471 cubic Angstrom, versus a typical compact drug-like pocket of ~200-400 cubic Angstrom) is characteristic of an extended protein-substrate (ubiquitin-binding) groove, not a small-molecule-sized cavity; a real inhibitor would engage only a sub-region of this groove, and defining which sub-region is druggable would require either a fragment screen or a co-crystallized reference ligand, neither of which exists for USP34.",
    ),
    dict(
        structure="7W3U (ubiquitin-bound)", pocket_label="Pocket 1 (top-ranked among 94 total pockets detected)",
        fpocket_druggability_score=0.873, volume_A3=None, num_alpha_spheres=103,
        key_lining_residues="USP34 chain B: Thr1984, Ile2009, Asn2011, Glu2026/2027, Phe2028, Tyr2029, Thr2030, Val2031, Arg2032, Gln2034, Ser2044, Glu2047, Val2048, Lys2073, Ala2075, Met2089. Ubiquitin chain E: Lys6, Thr7, Gly10, Ile44, Gly47, Thr66, His68, Leu69, Val70 -- includes the ubiquitin Ile44 hydrophobic patch, a well-characterized, precedented ubiquitin-recognition hotspot targeted by other ubiquitin-pathway PPI-disruptor programs.",
        location_description="This is the USP34-ubiquitin PROTEIN-PROTEIN INTERFACE, centered on ubiquitin's Ile44 patch, NOT the catalytic Cys1903/His2164 dyad itself (which sits in a separate, lower-scoring sub-pocket in this bound structure, consistent with the catalytic cysteine being covalently occupied by the probe in this particular crystal).",
        interpretation_caveat="The Ile44-patch-adjacent location is a genuinely precedented class of druggable/disruptable PPI interface in the ubiquitin system generally (qualification: precedent exists for the CLASS of interface, not for a validated USP34 compound specifically). This pocket is partly defined by the bound ubiquitin probe itself, so its exact shape/druggability in the ABSENCE of a bound ubiquitin (i.e. in a free, substrate-accessible enzyme) is not directly observed by this structure and should not be assumed identical.",
    ),
    dict(
        structure="7W3U (ubiquitin-bound)", pocket_label="Pocket 2 (second-ranked, chain A)",
        fpocket_druggability_score=0.636, volume_A3=790.9, num_alpha_spheres=113,
        key_lining_residues="Residues 2009-2089 (chain A) -- this range directly brackets the structural zinc-binding module identified this session (Cys2018/His2020/Cys2062/Cys2065)",
        location_description="Adjacent to/overlapping the structural Zn-binding module, distinct from the catalytic dyad and distinct from the Ile44-patch interface above.",
        interpretation_caveat="A pocket detected near a structural metal-coordination site can sometimes be scored artificially favorably by geometric pocket-detection algorithms because the metal-coordination geometry itself creates a cavity-like signature -- this is flagged explicitly as a possible confound, not asserted as a genuine allosteric drug pocket. Disrupting a structural (fold-stabilizing) zinc site is also mechanistically more likely to cause global unfolding/destabilization than a clean allosteric modulation, a materially different (and less attractive) drug mechanism that would need explicit confirmation before pursuit.",
    ),
    dict(
        structure="7W3R (apo)", pocket_label="Small sub-pocket fragment overlapping Cys1903 (fpocket 'Pocket 18' of 54 -- reported for full transparency since it also directly contains the catalytic cysteine)",
        fpocket_druggability_score=0.000, volume_A3=372.5, num_alpha_spheres=25,
        key_lining_residues="Includes Cys1903 among others (25 alpha spheres total)",
        location_description="A smaller alpha-sphere cluster that fpocket separately enumerates within the same general catalytic-cleft region already captured, more comprehensively, by the top-ranked Pocket 6 above.",
        interpretation_caveat="Included for full transparency: fpocket's automatic pocket-fragmentation can split one biological cleft into multiple listed 'pockets' of very different scores; this fragment's near-zero score should NOT be read as contradicting Pocket 6's high score for the same general region -- Pocket 6 (which also contains Cys1903, more completely) is the more informative call, and this entry exists only so the raw fpocket output is not selectively reported.",
    ),
]

# ---------------------------------------------------------------------------
# 9. DOCKING DECISION
# ---------------------------------------------------------------------------
DOCKING_QUESTIONS_ROWS = [
    dict(question="1. Does USP34 contain a conventional active-site pocket?", answer="Partially. The catalytic Cys1903/His2164 dyad sits within fpocket's top-ranked, high-druggability-scored (0.845) pocket in the apo structure -- but that pocket is large (1471 cubic Angstrom) and groove-shaped (an extended ubiquitin-binding surface), not a conventional compact small-molecule active site. A sub-pocket within this groove may be conventionally druggable, but which sub-pocket is not established by geometry alone."),
    dict(question="2. Could catalytic-cysteine covalent inhibition be feasible?", answer="Plausible, with real direct structural support. The LINK record in 7W3U crystallographically proves Cys1903 is a reactive, covalently-modifiable nucleophile (it is covalently bonded to the probe's allylamine warhead in the solved structure) -- this is the strongest, most direct evidence identified for any targeting modality on either USP34 or VEZF1 in this project. No optimized covalent INHIBITOR (as opposed to the activity-based research PROBE used to solve this structure) exists yet."),
    dict(question="3. Does USP34 have unique neighboring residues that might support selectivity?", answer="Not established this session. A proper selectivity assessment requires structural alignment against other USP-family catalytic domains (USP7, USP30, USP1) to identify non-conserved pocket-lining residues -- this specific comparative-alignment analysis was not performed this session (disclosed scope limitation) and would be a defensible next computational step."),
    dict(question="4. Are there plausible allosteric sites?", answer="One candidate identified, with an important caveat. The zinc-module-adjacent pocket (7W3U Pocket 2, score 0.636) is spatially distinct from the catalytic site, but disrupting a structural zinc-coordination site is mechanistically more likely to destabilize the fold than to allosterically modulate activity -- this is a candidate for further investigation, not a confirmed allosteric drug site. The real, literature-documented apo-inactive/ubiquitin-bound-active conformational transition (independently re-confirmed by direct distance measurement this session in one crystallographic copy, chain A: 3.94 A to 3.37 A -- though other copies in the same structures show more variable distances, see USP34_structure_inventory.tsv, so this is real single-copy evidence, not a uniform structure-wide effect) is the more mechanistically-grounded basis for an allosteric/conformational-trap strategy, though no specific allosteric pocket stabilizing the inactive state was structurally characterized this session."),
    dict(question="5. Is the ubiquitin-binding interface itself targetable?", answer="Plausibly yes. fpocket's single highest-scoring pocket in the ubiquitin-bound structure (0.873) is centered on the USP34-ubiquitin interface at ubiquitin's Ile44 patch -- a precedented class of druggable/disruptable PPI interface in the ubiquitin system generally. This would be a PPI-disruptor (e.g. stapled-peptide or fragment-based) strategy, distinct from and complementary to the catalytic-site covalent strategy above."),
    dict(question="6. Would a small molecule likely need to stabilize an inactive conformation?", answer="Plausibly relevant given the real, measured apo-to-bound conformational tightening in one crystallographic copy (chain A: 3.94 A to 3.37 A at the Cys1903-His2164 dyad) -- an allosteric ligand that locks the apo-inactive conformation (analogous conceptually, not structurally validated, to USP14's IU1) is a mechanistically coherent hypothesis, but the effect is not uniform across all copies in these structures (other chains range 3.10-4.98 A) and no specific pocket was identified this session that is structurally confirmed to control this transition."),
]

DOCKING_DECISION_ROW = dict(
    decision="DOCKING_NOT_YET_JUSTIFIED",
    justification=(
        "Per the four required conditions (a plausible pocket exists; an appropriate structure/conformation is available; a scientifically defensible ligand set "
        "exists; docking would answer a real question), condition (a) and (b) are met (real pockets and both apo/bound structures exist, established above), but "
        "(c) is NOT met: no USP34-validated ligand, chemical probe, fragment-screen hit, or even a closely-homologous USP-family co-crystallized small molecule "
        "was identified anywhere in this project's three research passes to serve as a positive control for validating a docking protocol's pose-prediction or "
        "scoring accuracy on this specific pocket -- docking without any way to check whether the method reproduces a KNOWN correct answer on a related system "
        "risks generating numerically specific-looking but uncalibrated results, which is exactly the 'misleading at this stage' scenario this phase was asked "
        "to watch for. Additionally, the top-scoring catalytic-adjacent pocket is unusually large/groove-shaped (1471 cubic Angstrom) with no established "
        "sub-pocket boundary, so even a well-calibrated docking run would face an ill-defined search space. No molecular-docking software (AutoDock Vina or "
        "equivalent) was installed or run this session -- only the pocket-DETECTION tool fpocket was used, which is a distinct, more defensible task (geometric "
        "cavity identification) than pose/affinity prediction (docking)."
    ),
    alternative_recommended="Structure-based roadmap (see report Part 9/10) instead of blind or arbitrary docking: (1) a genuine fragment-based screening campaign against the real 7W3R/7W3U structures, using the two independently-supported chemical starting points identified this session -- covalent-cysteine-directed fragments (justified by the crystallographically-proven Cys1903 reactivity) and Ile44-patch-directed PPI fragments (justified by the real fpocket-scored interface pocket); (2) a structural alignment against USP7/USP30/USP1 to assess pocket selectivity BEFORE any compound design; (3) revisit the DOCKING_JUSTIFIED decision once either a fragment-screen hit or a homologous-DUB reference ligand exists to calibrate a docking protocol against.",
)

# ---------------------------------------------------------------------------
# 12. SUCCESS / FAILURE CRITERIA
# ---------------------------------------------------------------------------
SUCCESS_FAILURE_ROWS = [
    dict(target="USP34", hypothesis_type="Tamoxifen-specific sensitiser (single hypothesis type -- USP34's DepMap baseline dependency is 0.0%, so a direct-cancer-dependency hypothesis is not separately proposed)",
         supports_criteria="USP34 perturbation alone causes only a MODEST effect on resistant-cell viability/clonogenic survival; tamoxifen alone is weak in the resistant line (confirming resistance); the COMBINATION produces a clearly larger effect than either single arm, preferably visible as a leftward shift in the tamoxifen dose-response curve and/or a synergy score (Bliss/CI) indicating more-than-additive interaction; the EXP-2 normal-cell (MSC/osteogenic) comparator is substantially LESS affected by the same perturbation level.",
         weakens_criteria="USP34 perturbation kills cells equally with or without tamoxifen (no combination benefit -- suggests a generic dependency or toxicity, not a sensitisation mechanism); the combination shows no additional benefit over the better single arm; the EXP-2 normal-cell comparator shows strong functional impairment (e.g. blunted osteogenic differentiation) at a similar perturbation level, indicating a narrow or absent therapeutic window.",
         negative_criteria="Neither a direct effect nor a combination effect is observed in the resistant line at all -- the hypothesis as stated is not supported by this experiment."),
    dict(target="VEZF1", hypothesis_type="Dual-action target (both sub-hypotheses A and B tested jointly, per the locked question in Part 1)",
         supports_criteria="DUAL ACTION supported if: VEZF1 suppression impairs resistant cancer cells ALONE (consistent with real 27.3% DepMap dependency) AND the combination further improves tamoxifen response beyond that direct effect. PURE DEPENDENCY (not dual-action, but still a real finding) if: VEZF1 suppression kills cells alone and tamoxifen adds no further meaningful benefit on top of that. PURE SENSITISER (the alternative single-mechanism outcome) if: little-to-no baseline single-agent effect, but a large combination-specific effect appears.",
         weakens_criteria="A large single-agent VEZF1-suppression effect that is NOT further enhanced by tamoxifen at all would argue against the 'AND' (dual-action) framing specifically, even though it would still support a pure-dependency finding; strong functional impairment of the EXP-4 normal-cell comparators (endothelial tube-formation, cardiomyocyte contractility) at a similar perturbation level would narrow the apparent therapeutic window.",
         negative_criteria="NEGATIVE if neither a direct single-agent effect NOR a combination-enhancement effect occurs in the resistant line -- both the DepMap-based dependency hypothesis and the Hany-based sensitisation hypothesis would be unsupported by this specific cell-based follow-up, notwithstanding the real underlying CRISPR/DepMap screen signals that motivated the experiment."),
]

# ---------------------------------------------------------------------------
# 13. POSTER-READY FINAL TRANSLATIONAL CONCLUSIONS
# ---------------------------------------------------------------------------
FINAL_CONCLUSIONS_ROWS = [
    dict(target="USP34", role="LEAD TARGET",
         evidence_chain="Functional CRISPR sensitisation (Hany FDR=0.042, frozen) + low baseline ER+/luminal cancer-cell dependency (DepMap 26Q1 = 0.0%) + real catalytic targetability (crystallographically-confirmed reactive Cys1903, PDB 7W3R/7W3U) = LEAD COMBINATION-TARGET HYPOTHESIS",
         poster_line="USP34: functional CRISPR sensitisation + low baseline ER+ dependency + real catalytic targetability -> lead combination-target hypothesis, not yet clinically validated.",
         caveats_carried_forward="No validated USP34-selective inhibitor exists yet (first-pass validation must be genetic, per EXP-1); real bone/osteogenic liability evidence exists and is the priority normal-cell comparator (EXP-2); strong human genetic constraint (LOEUF=0.152, most constrained candidate examined in this project) is a caution, not a toxicity prediction."),
    dict(target="VEZF1", role="SECOND / BACKUP TARGET",
         evidence_chain="Strong functional CRISPR sensitisation (Hany FDR=0.037, nominally the strongest of the two, frozen) + real baseline ER+/luminal cancer-cell dependency (DepMap 26Q1 = 27.3%) = DUAL-ACTION BIOLOGICAL HYPOTHESIS, limited by poor direct druggability",
         poster_line="VEZF1: strong CRISPR sensitisation + baseline ER+ dependency -> dual-action biological hypothesis, limited by poor direct druggability; TEAD1 remains an unvalidated indirect-targeting hypothesis pending EXP-5.",
         caveats_carried_forward="Cardiovascular/developmental liability evidence (postnatal zebrafish cardiac-contractility finding, PMID 31911272) is stronger and more directly causal than any skeletal-muscle or bone finding for VEZF1; TEAD1 is explicitly NOT a validated indirect-targeting strategy -- EXP-5 exists specifically to test, and could reject, that hypothesis before it is used for anything therapeutic."),
]

