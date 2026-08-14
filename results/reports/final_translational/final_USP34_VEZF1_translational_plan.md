# Final USP34 / VEZF1 translational + structure plan

**Frozen ranking, unchanged by this phase**: USP34 = LEAD TARGET, VEZF1 =
SECOND/BACKUP TARGET. No frozen upstream evidence (discovery, TCGA,
DepMap, druggability/safety, lead-target deep-dive) is modified,
recomputed, or reopened here. No additional expression datasets were
added. This phase is forward-looking experimental design, plus one
genuinely new analysis: a real, locally-run structural/pocket analysis of
the two experimental USP34 structures (see Part 7-8 for exact provenance
-- PDB files downloaded from RCSB this session, pocket detection run
locally with fpocket 4.2.3, catalytic-residue distances computed directly
from the downloaded coordinates).

---

## 1. Locked final biological questions

**USP34**: Does reducing USP34 restore or enhance tamoxifen response in
acquired tamoxifen-resistant ER+ breast cancer cells **without simply
causing nonspecific baseline toxicity**? Not described as proven anywhere
in this report -- this is the question EXP-1 is designed to answer.

**VEZF1**: Does reducing VEZF1 **(A)** impair resistant ER+ cancer-cell
fitness on its own, **AND/OR (B)** further restore/enhance tamoxifen
sensitivity? Both sub-questions are tested jointly and explicitly kept
separate in EXP-3's analysis plan (Part 4) -- neither is described as
proven.

Full table: `results/tables/final_translational/final_experimental_design.tsv`.

---

## 2. USP34 experiment (EXP-1)

**Model system**: acquired tamoxifen-resistant ER+ breast cancer cells,
with a matched parental ER+ line run in parallel (required to distinguish
a resistance-specific effect from a generic ER+ cell effect).

**Four arms**: (1) control; (2) tamoxifen/4-OHT; (3) USP34 perturbation;
(4) USP34 perturbation + tamoxifen/4-OHT.

**Recommended first-pass perturbation: CRISPR knockout (KO)**, not
CRISPRi, siRNA/shRNA, or an inhibitor (none exists). Rationale: no
validated USP34-selective chemical probe or inhibitor currently exists
(reconfirmed across three independent research passes in this project),
so a pharmacological arm is not available and must not be simulated by a
proxy claim. siRNA/shRNA give only partial, often inconsistent knockdown
of a large (3546-aa) protein, risking an uninterpretable false-negative.
CRISPRi reduces transcription, but a pre-existing protein pool of a
stable enzyme can persist through a short assay window. CRISPR KO gives
the cleanest, most informative first-pass answer to "does USP34 loss
produce this phenotype at all" -- explicitly **not** the same question as
"what would partial pharmacological inhibition do," which is deferred to
a later dose-titratable follow-up (e.g. an inducible degron or
analog-sensitive allele) only once a positive KO signal justifies it.

**USP34 EMT/stemness counter-evidence -- incorporated this revision
(PMID 28499884)**: Wang et al., *Cell Signal* 2017;36:230-239, "Inhibition
of ubiquitin-specific protease 34 (USP34) induces epithelial-mesenchymal
transition and promotes stemness in mammary epithelial cells" (title/
abstract independently re-fetched and verified this session, not
summarized from memory). **Model**: NMuMG cells -- a normal,
non-transformed MOUSE mammary epithelial line, **not** a cancer line, not
ER+, not a tamoxifen-resistance model (correlative expression comparisons
also used MDA-MB-231/4T1, both ER-negative mesenchymal-like lines; plus
an in vivo mammary-fat-pad transplantation). **Perturbation**: USP34
knockdown (KD, not KO). **Findings, verbatim from the abstract**: USP34
KD induced EMT -- upregulation of N-cadherin, phospho-Smad3, Snail, and
active-beta-catenin; downregulation of Axin1 and E-cadherin -- plus
invasive behavior and enhanced mammosphere formation with upregulated
Nanog/Oct4/Sox2; transplanted USP34-KD cells still reconstituted a normal
ductal tree in vivo within 3 months. **This does NOT invalidate the
USP34 tamoxifen-sensitisation hypothesis.** It is counter-evidence / a
potential context-specific liability, because this study is not
equivalent to acquired-tamoxifen-resistant-ER+-breast-cancer +
USP34-perturbation + tamoxifen -- it used a normal, non-cancer, non-ER+
mammary line and never combined USP34 loss with tamoxifen. Prior
mammary-epithelial work suggests USP34 loss can promote EMT/stem-like
features in some contexts, motivating explicit EMT/stemness monitoring
during tamoxifen-resensitisation experiments. This report does **not**
write "USP34 inhibition promotes aggressive breast cancer."

**Primary readouts**: viability (5-7 day dose-response), apoptosis
(Annexin V/cleaved caspase-3), clonogenic survival (10-14 day colony
assay). **Interaction framework, now required rather than optional**:
sensitisation is explicitly **not** defined merely as "combination >
either single arm" -- a full tamoxifen/4-OHT dose-response curve is
generated under both non-targeting-control and USP34-KO backgrounds and
compared directly (EC50/IC50 shift), with a predefined Bliss-independence
and/or Chou-Talalay combination-index analysis applied to the full
dose-response matrix before any sensitisation claim is made.

**Mechanistic readouts**: USP34 protein loss (western blot, confirms
KO); AXIN1 ubiquitination status and active beta-catenin (the one
experimentally-validated USP34 substrate axis, PMID 21383061) -- included
because it is real prior evidence, not assumed relevant to resistance
without this readout; ESR1 target-gene qPCR (TFF1, PGR, GREB1) to test
whether USP34 loss perturbs ER transcriptional output directly.
**EMT/stemness monitoring, added this revision, directly motivated by
PMID 28499884 and using the same marker panel that paper reported**:
E-cadherin/CDH1, N-cadherin/CDH2, SNAI1/Snail, AXIN1, active beta-catenin
(the latter two doing double duty for the Wnt-mechanism question above
and the EMT-liability question, since PMID 28499884 implicates the same
axis); optional migration/invasion and mammosphere-formation assays.

**Outcome categories, explicitly defined**: **IDEAL** -- sensitisation
without meaningful EMT/stemness induction. **CONCERNING** -- sensitisation
occurs but EMT/stemness/invasion also increase substantially (a real,
monitorable liability, not a disproof of sensitisation). **PURE GENERAL
TOXICITY** -- USP34 perturbation kills regardless of tamoxifen. **NEGATIVE**
-- no meaningful sensitisation, independent of the EMT question.

**Controls/replicates**: non-targeting sgRNA control; >=3 biological
replicates; parental-vs-resistant run in parallel; rescue via
sgRNA-resistant USP34 cDNA (ideally including a catalytic-dead C1903S
point mutant to test catalytic-activity dependence) recommended if an
initial phenotype is observed.

---

## 3. USP34 normal-cell comparators (EXP-2A + EXP-2B)

Two comparator concepts are used, per the explicit instruction that a
single comparator conflates two different questions.

**EXP-2A -- lineage/cancer-selectivity comparator**: normal human mammary
epithelial cells (primary HMECs if feasible, MCF10A as an alternative),
given the same USP34 perturbation as EXP-1. Purpose: does comparable
USP34 perturbation disproportionately affect malignant ER+ cells relative
to normal mammary cells? **Promoted to a co-primary comparator this
revision** (not a secondary add-on), specifically because it is now
directly informative for the PMID 28499884 EMT/stemness question too --
the same EMT/stemness marker panel from EXP-1 is read out here, in a
normal human mammary background, directly extending PMID 28499884's own
(mouse) mammary-epithelial context.

**EXP-2B -- known-liability comparator: primary human bone-marrow-derived
MSCs undergoing osteogenic differentiation induction** (same USP34
CRISPR-KO perturbation, for direct comparability) -- not a panel of every
possible normal cell type, and explicitly **not** a skeletal-muscle
myoblast/myotube comparator.

**Why MSC/osteogenic, not muscle, for EXP-2B**: of USP34's candidate
liability axes, only bone/osteogenic mesenchymal differentiation has ANY
demonstrated functional requirement for USP34 in a normal cell type,
replicated across three independent findings in this project (BMP2/
Smad1/RUNX2-mediated osteogenesis; NFIC-mediated tooth-root formation;
a 2026 human haploinsufficiency syndrome with limb anomalies). Skeletal
muscle has high USP34 RNA expression (54.1 nTPM) but **zero** functional
evidence in any of three independent, alias-inclusive PubMed searches
across this project's phases -- it does not need to be the primary
comparator. MSC/osteogenic directly extends the exact published assay
(PMID 30181118) that established this real liability, giving a
literature-anchored expectation to compare against.

**Readouts**: viability/proliferation under USP34 perturbation;
osteogenic differentiation markers (RUNX2 expression, alkaline
phosphatase activity, Alizarin Red mineralization at day 14-21); relative
sensitivity -- does the perturbation level needed to see a cancer-cell
effect in EXP-1 also measurably impair normal osteogenic differentiation,
or is there a detectable window?

**Explicit framing**: this is a **preliminary therapeutic-window
experiment**, not a safety study, and a favorable result must never be
read as proof that USP34 inhibition is safe. A secondary MCF10A (normal
mammary epithelial) comparator is a reasonable add-on if resources allow,
controlling for a generic pan-epithelial effect, but is not the priority.

---

## 4. VEZF1 experiment (EXP-3)

Same four-arm structure as EXP-1: (1) control; (2) tamoxifen/4-OHT; (3)
VEZF1 suppression; (4) VEZF1 suppression + tamoxifen/4-OHT, in the same
resistant + parental model system.

**Recommended perturbation: CRISPRi** (dCas9-KRAB at the VEZF1 TSS) as
primary, with inducible CRISPR KO as a secondary/confirmatory option.
Rationale: VEZF1 has no validated direct pharmacology at all, so unlike
USP34 there is no future drug mechanism to match the perturbation to --
the priority is interpretable functional evidence. CRISPRi gives clean,
tunable, reversible silencing without the specific confound a
DNA-binding-domain-truncated CRISPR-KO allele can create for a
transcription factor (a partial product could retain dominant-negative
or gain-of-function activity, complicating interpretation in a way that
is less of a concern for an enzyme like USP34). siRNA/shRNA are
acceptable secondary options.

**Explicit separation (required by Part 1's locked question)**:
- **DIRECT CANCER DEPENDENCY**: arm 3 vs arm 1 -- does VEZF1 suppression
  alone impair resistant cells, consistent with the real 27.3% DepMap
  26Q1 ER+/luminal baseline-dependency signal?
- **TAMOXIFEN SENSITISATION**: arm 4 vs what arms 2 and 3 alone would
  predict -- does the combination outperform expectation, consistent
  with the real Hany FDR=0.037 signal?
These are analyzed as two distinct questions and never collapsed into a
single "VEZF1 matters" readout.

**Readouts**: same viability/apoptosis/clonogenic/dose-response panel as
EXP-1, plus ER transcriptional response (TFF1/PGR/GREB1, for mechanistic
comparability with USP34), plus a VEZF1 target-program readout -- **this
is defensible and included**: VEZF1's own documented target genes
(VEGFR2/KDR, TIMP3, MMP2, from this project's own literature-mechanism
phase), confirming CRISPRi is functionally silencing the VEZF1 program,
not merely reducing mRNA without functional consequence.

**Controls/replicates**: non-targeting sgRNA/scrambled siRNA control;
>=3 biological replicates; parental-vs-resistant in parallel; rescue via
sgRNA-resistant VEZF1 cDNA recommended if a phenotype is observed.

---

## 5. VEZF1 normal-cell comparator (EXP-4)

**Recommended minimum comparator: primary human vascular endothelial
cells** (e.g. HUVECs or a microvascular line) as the priority choice,
**with iPSC-derived cardiomyocytes as a secondary comparator** if
resources allow -- not a larger panel.

**Why endothelial + cardiomyocyte**: VEZF1's founding, most-replicated
biology is vascular/endothelial (its namesake, embryonic knockout
phenotype, documented target genes VEGFR2/TIMP3/MMP2) -- the single most
literature-anchored comparator available. The newer, real cardiac-muscle
finding (PMID 31911272: zebrafish knockdown causally impairing postnatal
cardiac contractile function via Myh7/TEAD-1) is genuinely important and
is why cardiomyocytes are recommended as a secondary addition rather than
omitted -- but it is a single, partial-knockdown, non-human finding, and
testing it in a human iPSC-cardiomyocyte system is the appropriate
translational next step, not an assumption that the zebrafish phenotype
will reproduce. **Animal (zebrafish) cardiac knockdown findings are not
equated with expected adult human drug toxicity anywhere in this
design.**

**Readouts**: endothelial -- viability, tube-formation/angiogenic-capacity
assay, VEGFR2/TIMP3/MMP2 expression under VEZF1 perturbation.
Cardiomyocyte (if included) -- viability, contractile-function readout
(video-based contractility or Ca2+ imaging, directly analogous to PMID
31911272's beta-adrenergic-response readout), Myh7/beta-MHC expression.

**Explicit framing**: preliminary therapeutic-window experiment, not a
safety claim.

---

## 6. TEAD1 hypothesis test (EXP-5) -- hypothesis only, not a therapeutic arm

TEAD1 is **not** established as a VEZF1 indirect target (frozen
conclusion from the prior phase, unchanged). This experiment exists
purely to test, and could reject, the hypothesis.

**Arms**: (1) control; (2) pan-TEAD inhibitor (a VT3989-chemotype
compound -- no TEAD1-selective compound remains in active clinical
development, so a pan-TEAD tool compound is the only available chemical
route); (3) VEZF1 suppression (same CRISPRi reagent as EXP-3, included as
an internal positive-control reference for what a genuine VEZF1-program
reduction looks like on the readout panel below); (4) TEAD inhibitor +
tamoxifen -- **only if arm 2 shows a genuine VEZF1 effect**.

**Readouts**: VEZF1 RNA (qPCR), VEZF1 protein (western blot), the
VEZF1-dependent target program (VEGFR2/TIMP3/MMP2, same panel as EXP-3
for direct comparability). **Mandatory control**: a canonical TEAD/YAP
target-engagement readout (e.g. CTGF/CCN2, CYR61/CCN1) to confirm the
inhibitor is pharmacologically active in this cell line at all -- without
this, a negative VEZF1 result would be uninterpretable (target not
engaged vs. target engaged but VEZF1 unaffected).

**Decision rule**: if no VEZF1 RNA/protein/target-program change occurs
despite confirmed TEAD/YAP engagement -- **reject TEAD1 as an indirect
VEZF1-targeting strategy**, explicitly and permanently for this specific
hypothesis. If VEZF1 function decreases -- proceed to arm 4. **TEAD1 must
not be incorporated into the main USP34/VEZF1 therapeutic model before
this experiment produces a positive result.**

Full comparator table: `results/tables/final_translational/final_normal_cell_comparators.tsv`.

---

## 7. USP34 structural druggability

Full table: `results/tables/final_translational/USP34_structure_inventory.tsv`.
Figure: `02_USP34_structure_targetability.png`.

**Provenance**: both structures were freshly downloaded this session
directly from `files.rcsb.org/download/{7W3R,7W3U}.pdb` and parsed
directly (resolution, chains, residue ranges, missing residues, HETATM
records, and LINK/covalent-bond records extracted from the raw coordinate
files, not retyped from the earlier PMID 35588869 citation).

**7W3R -- USP34 catalytic domain, apo.** X-ray, **1.92 A**. Chains A, B
(construct residues 1892-2261 as modeled). **105 residues** listed as
missing/disordered (REMARK 465) within that range -- substantial, typical
of a flexible unliganded catalytic domain. Catalytic Cys1903 and His2164
both directly confirmed resolved. Bound ligand: a structural **Zn2+ ion**,
coordinated by Cys2018/His2020/Cys2062/Cys2065 (LINK records) -- **a
distinct site from the catalytic dyad, not previously highlighted in
this project's earlier structural summaries** and identified fresh this
session.

**7W3U -- USP34 catalytic domain + ubiquitin-propargylamide (UbPA)
probe.** X-ray, **3.13 A**. Chains A, B, C (USP34, residues 1892-2268)
and D, E, F (covalently-linked ubiquitin, residues 1-75 of ubiquitin's
canonical 76 -- the native C-terminal Gly76 position is where the probe
chemistry is covalently attached, so "residues 1-75 + probe adduct" is
the technically precise description, not an unmodified full-length
ubiquitin). 141 residues missing (more than the apo structure, plausibly reflecting
the lower resolution rather than genuinely greater flexibility -- stated
as a caveat, not over-interpreted). Same structural Zn2+ site present in
all three USP34 chains. **Direct crystallographic proof of covalent
catalytic-cysteine engagement**: LINK records show Cys1903(SG) covalently
bonded to the probe's allylamine warhead (residue AYE, C2 atom) at
1.59-2.48 A across the three copies in the asymmetric unit -- this is
observed fact from the deposited structure, not an inference from
sequence homology.

**Measured conformational transition** (computed directly from
coordinates this session, not merely cited from the literature, and
independently re-verified by Codex from a fresh download): in **chain A**,
the Cys1903(SG)-His2164(ND1) distance is **3.94 A in the apo structure and
3.37 A in the probe-bound structure** -- a real, independently-computed
~0.6 A tightening in this specific copy that numerically corroborates the
literature's qualitative "apo-inactive, realigns on ubiquitin engagement"
description (PMID 35588869). **This is a single-copy observation, not a
uniform structure-wide effect, and is reported with that caveat**: the
apo structure's chain B measures 4.98 A at the same atom pair (looser than
chain A), and the probe-bound structure's three copies range from 3.10 A
(chain C, tighter than chain A) to 3.95 A (chain B, essentially unchanged
from apo chain A) -- real crystallographic heterogeneity that is disclosed
here rather than smoothed into a single clean number.

**No additional experimentally-solved USP34 structure** was identified
in a repeated RCSB search this session; 7W3R/7W3U remain the only two,
both covering only the catalytic domain (~12% of the full-length
protein). No AlphaFold model exists (reconfirmed).

---

## 8. USP34 pocket / ligandability analysis

Full table: `results/tables/final_translational/USP34_pocket_analysis.tsv`.
Figure: `02_USP34_structure_targetability.png`.

**Method**: fpocket 4.2.3 (conda-forge), installed and run locally this
session on both downloaded structures -- a real, reproducible, open-source
tool, per the explicit instruction to prefer such tools where available.
No pocket is described as "druggable" without its fpocket score
attached; a score is a geometric/physicochemical heuristic, not proof a
drug-like molecule can be found.

**7W3R (apo), 54 total pockets detected. Top-ranked: druggability
0.845**, volume 1471 A^3, 136 alpha spheres. **This pocket directly
contains Cys1903 and His2164** (the catalytic dyad) plus Gln1976 -- the
residue found this session to sit closest (3.50 A) to His2164 among all
Asp/Asn/Glu/Gln side chains checked by direct 3D distance calculation
(a geometric observation from real coordinates, **not** a
literature-confirmed catalytic-triad assignment, reported with that
explicit caveat). The very large volume (roughly 4-7x a typical compact
drug-like pocket) is characteristic of an extended
protein-substrate/ubiquitin-binding groove, not a small-molecule-sized
cavity -- a real inhibitor would engage only a sub-region, and defining
that sub-region would require a fragment screen or reference ligand,
neither of which exists.

**7W3U (ubiquitin-bound), 94 total pockets. Top-ranked: druggability
0.873**, at the **USP34-ubiquitin interface**, centered on ubiquitin's
**Ile44 hydrophobic patch** -- a precedented class of druggable/
disruptable interface in the ubiquitin system generally (precedent for
the class, not a validated USP34-specific compound). **Second-ranked:
druggability 0.636**, adjacent to the newly-identified structural Zn
module (residues 2009-2089) -- flagged explicitly as a possible
scoring artifact of metal-coordination geometry, and mechanistically more
likely to destabilize the fold than to allosterically modulate activity
if disrupted; not promoted as a confirmed allosteric site.

**Answers to the six specific questions**:
1. *Conventional active-site pocket?* Partially -- real, high-scoring, but
   large/groove-shaped, not compact.
2. *Covalent cysteine inhibition feasible?* Plausible, with the strongest
   direct evidence in this project: Cys1903 is crystallographically
   proven reactive.
3. *Selectivity-supporting neighboring residues?* Not established this
   session -- would require structural alignment against USP7/USP30/USP1
   (a defensible next step, not performed here).
4. *Plausible allosteric sites?* One candidate (Zn-module-adjacent
   pocket), with the destabilization caveat above; the apo-to-bound
   conformational transition is the more mechanistically-grounded basis
   for an allosteric strategy, though no specific pocket controlling it
   was structurally characterized.
5. *Ubiquitin-binding interface targetable?* Plausibly yes -- real,
   high-scoring, Ile44-patch-centered, precedented PPI-interface class.
6. *Would a small molecule need to stabilize an inactive conformation?*
   Plausibly relevant, given the real measured 3.94->3.37 A tightening in
   one crystallographic copy (chain A) -- though the other copies in these
   structures show more variable distances (3.10-4.98 A), so this is not
   a uniform effect -- and not yet structurally confirmed to a specific pocket.

---

## 9-10. Docking decision

Full table: `results/tables/final_translational/USP34_docking_decision.tsv`.

**DECISION: DOCKING_NOT_YET_JUSTIFIED.**

Of the four required conditions, (a) a plausible pocket exists and (b) an
appropriate structure/conformation is available are both met (Part 8).
**Condition (c) is not met**: no USP34-validated ligand, chemical probe,
fragment-screen hit, or closely-homologous USP-family co-crystallized
small molecule was identified anywhere in this project's three research
passes to serve as a positive control for validating a docking protocol's
pose-prediction or scoring accuracy on this specific pocket. Docking
without any way to check whether the method reproduces a known-correct
answer on a related system risks generating numerically specific-looking
but uncalibrated results -- exactly the "misleading at this stage"
scenario this phase was asked to watch for. The top-scoring
catalytic-adjacent pocket is also unusually large/groove-shaped with no
established sub-pocket boundary, so even a well-calibrated run would face
an ill-defined search space. No docking software (AutoDock Vina or
equivalent) was installed or run this session -- only the pocket
**detection** tool fpocket was used, a distinct and more defensible task
than pose/affinity prediction.

**Structure-based roadmap instead (per the user's own fallback
instruction)**:
1. A genuine fragment-based screening campaign against the real
   7W3R/7W3U structures, using the two independently-supported chemical
   starting points identified this session: covalent-cysteine-directed
   fragments (justified by the crystallographically-proven Cys1903
   reactivity) and Ile44-patch-directed PPI fragments (justified by the
   real fpocket-scored interface pocket).
2. Structural alignment against USP7/USP30/USP1 to assess pocket
   selectivity **before** any compound design.
3. Revisit this DOCKING_JUSTIFIED/NOT_YET_JUSTIFIED decision once either
   a fragment-screen hit or a homologous-DUB reference ligand exists to
   calibrate a docking protocol against.

No figure 04 (docking hypothesis) was produced, consistent with the
user's own instruction that it is conditional on a JUSTIFIED decision.

---

## 11. Final experimental priority

Figure: `01_final_experimental_strategy.png`.

EXP-1 (USP34 perturbation +/- tamoxifen) -> EXP-2A/EXP-2B (USP34
normal-cell comparators) -> EXP-3 (VEZF1 perturbation +/- tamoxifen) ->
EXP-4 (VEZF1 normal-cell comparator) -> EXP-5 (TEAD inhibition -> does
VEZF1 change?).

**If only one experiment can be done: EXP-1.** USP34 is the frozen lead
target; it has the cleaner single-mechanism hypothesis (tamoxifen-specific
sensitiser, not a dual mechanism requiring two readouts to interpret);
and it is the only one of the two genes with a real structural handle for
future medicinal chemistry, so a positive EXP-1 result would immediately
motivate the (currently unfunded) structure-based work in Parts 7-10.

---

## 12. Success / failure criteria

Full table: `results/tables/final_translational/final_target_success_failure_criteria.tsv`.

**USP34 supports sensitisation** if: perturbation alone causes a modest
effect; tamoxifen alone is weak in the resistant line; the combination is
clearly larger than either alone (ideally a visible tamoxifen dose-response
shift and/or synergy score); the EXP-2A/EXP-2B normal-cell comparators are
substantially less affected, AND EMT/stemness markers are not
substantially induced (see IDEAL vs CONCERNING in Part 2). **Weakens**
if: perturbation kills equally with or without tamoxifen; no combination
benefit; strong normal-cell functional impairment at a similar
perturbation level.

**VEZF1 -- dual action** supported if: suppression impairs resistant
cells alone AND the combination further improves tamoxifen response.
**Pure dependency** (still a real, useful finding, just not "dual") if:
suppression kills cells alone and tamoxifen adds no further benefit.
**Pure sensitiser** if: little baseline effect but a large
combination-specific effect. **Negative** if neither effect occurs.

---

## 13. Poster-ready final model

Figure: `03_USP34_VEZF1_final_translational_model.png`.

**USP34**: functional CRISPR sensitisation (Hany FDR=0.042) + low
baseline ER+ cancer dependency (DepMap 26Q1=0.0%) + real catalytic
targetability (PDB 7W3R/7W3U, Cys1903 crystallographically confirmed
reactive) = **lead combination-target hypothesis, not clinically
validated.**

**VEZF1**: strong CRISPR sensitisation (Hany FDR=0.037) + baseline ER+/
luminal cancer dependency (DepMap 26Q1=27.3%) = **dual-action biological
hypothesis, limited by poor direct druggability; TEAD1 remains an
unvalidated indirect-targeting hypothesis pending EXP-5.**

---

## Outputs

**Report**: `results/reports/final_translational/final_USP34_VEZF1_translational_plan.md`.

**Tables** (`results/tables/final_translational/`): `final_experimental_design.tsv`,
`final_normal_cell_comparators.tsv`, `USP34_structure_inventory.tsv`,
`USP34_pocket_analysis.tsv`, `USP34_docking_decision.tsv`,
`final_target_success_failure_criteria.tsv`, `final_translational_conclusions.tsv`.

**Figures** (`results/figures/final_translational/`): `01_final_experimental_strategy.png`,
`02_USP34_structure_targetability.png`, `03_USP34_VEZF1_final_translational_model.png`.
No figure 04 (docking hypothesis) -- Part 9-10 concluded DOCKING_NOT_YET_JUSTIFIED.

**Code**: `src/final_translational_data.py` (curated experimental designs
+ real structural/pocket findings, source-cited), `src/final_translational_build_tables.py`
(deterministic table builder), `src/final_translational_visualization.py`
(figure builder).
