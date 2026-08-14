# Druggability + normal-tissue / selectivity review

USP34, VEZF1, EML5, CITED2 -- frozen therapeutic shortlist and frozen
ranking (USP34 > VEZF1 > EML5 > CITED2 in the original poster-deliverable
ranking; USP34/CITED2/EML5 tied at raw follow-up rank 1, VEZF1 alone at
rank 4 in the independent-validation follow-up ranking), **both
unchanged by this phase**. No CRISPR/RNA-seq/TCGA/DepMap computation was
rerun; no frozen output was modified. This phase asks five new questions
that the prior phases (project function, TCGA human relevance, DepMap
baseline dependency) did not: which candidates are realistically
targetable, what intervention type is plausible, whether targeting them
could plausibly affect normal tissue (with bone/musculoskeletal biology
explicitly included), which candidates have the most favorable
cancer-vs-normal rationale, and which 1-2 deserve experimental follow-up.

**Method and provenance.** Every specific fact below (UniProt domain call,
PDB/AlphaFold accession, ChEMBL/Pharos/Open Targets query result, GTEx/HPA
expression value, gnomAD constraint number, OMIM/ClinGen classification,
MGI/IMPC phenotype, PubMed finding) was gathered by four independent,
source-verified research passes (one per candidate) against UniProt, RCSB
PDB, AlphaFold DB, ChEMBL, Pharos, Open Targets, GTEx Portal, Human
Protein Atlas, gnomAD, OMIM, ClinGen, MGI/IMPC, and PubMed/NCBI eutils,
each confirmed by directly fetching the source record, never recalled
from memory. Items the searches could not verify are recorded literally
as "NOT FOUND" rather than inferred. Full per-candidate detail with exact
citations is in `results/tables/druggability_safety/` (five tables plus
`verified_references.tsv`); this report summarizes and interprets those
tables, never retyping a number the tables don't already contain.

**Hard constraints honored throughout:** no candidate is called "safe";
DepMap (a cancer-cell-line dependency resource) is never used to infer
normal-tissue selectivity; high expression is never equated with
essentiality; genetic constraint is never equated with drug toxicity; no
claim rests on a single database without cross-checking where a second
source exists.

---

## Part 1 -- frozen biological interpretation carried forward, unchanged

| Candidate | Hany effect | Hany FDR | DepMap 26Q1 ER+/luminal strongly dependent | DepMap tier | Prior interpretation |
|---|---|---|---|---|---|
| USP34 | -1.391298 | 0.041685 | 0.0% | D_LOW | potential tamoxifen-specific sensitiser profile |
| VEZF1 | -1.602445 | 0.037258 | 27.3% | B_MODERATE | potential dual-action hypothesis (baseline dependency + tamoxifen sensitisation) |
| EML5 | -1.058423 | 0.148773 | 0.0% | D_LOW | strongest low-baseline-dependency profile; weakest mechanistic anchor |
| CITED2 | -1.495356 | 0.109955 | 0.0% | D_LOW | strongest independent mechanistic/literature support among non-Hany-significant candidates |

Source: `results/tables/evidence_freeze/THERAPEUTIC_SHORTLIST_FREEZE.tsv`,
`results/tables/independent_validation/DepMap_candidate_dependency.tsv`
(both frozen, read-only, unmodified by this phase).

---

## Part 2 -- druggability review

Full detail: `results/tables/druggability_safety/candidate_druggability.tsv`.
Figure: `results/figures/druggability_safety/01_four_candidate_druggability_summary.png`.

### USP34 -- POTENTIALLY_DRUGGABLE
A cysteine-protease deubiquitinase (DUB, peptidase family C19; UniProt
Q70CQ2) with a defined catalytic dyad (Cys1903/His2164) and **real
experimental structures** of the catalytic domain: PDB 7W3R (apo, 1.92 A)
and 7W3U (+ubiquitin-propargylamide probe, 3.13 A), both from Xu et al.
2022 (PMID 35588869). The apo enzyme is conformationally *inactive*
(misaligned catalytic histidine) until ubiquitin engagement -- a feature
that could in principle support an allosteric/conformational-trap
inhibitor strategy, analogous to how USP14's IU1 works, though this has
not been tested for USP34. No AlphaFold model exists (API returns 404),
an empirical gap, not a demonstrated absence of foldable structure. **No
USP34-specific chemical probe, inhibitor, or degrader exists in any
database checked** -- ChEMBL returns zero true USP34 entries (one text
hit is a mislabeled USP35 record); a Progenra SBIR screening effort (IPF
indication) discloses no compound structures. Related DUBs (USP30, USP7)
have preclinical/early-clinical tool compounds, offering class-level, not
USP34-specific, precedent. USP34 is the only one of the four candidates
with a real catalytic-domain structure and defined two-state mechanism --
the strongest structure-based drug-discovery starting point of the four,
even though no chemical matter exists yet.

### VEZF1 -- CURRENTLY_POORLY_DRUGGABLE
A C2H2 zinc-finger transcription factor (UniProt Q14119, six tandem zinc
fingers) with **no catalytic activity and no experimental structure of
any kind** (RCSB search returns zero hits, confirmed twice). The AlphaFold
model that exists is low-confidence (mean pLDDT 59.4; nearly half the
model "very low confidence"). One VEZF1-specific tool-compound series
exists (He et al. 2018, PMID 29970794): three compounds identified by
docking against a homology model of the DNA-binding interface, blocking
VEZF1-DNA binding at IC50 = 20-500 uM -- unoptimized, low-potency hits,
not chemical probes. Pharos classifies VEZF1 as target-development-level
"Tbio" (known biology, no qualifying bioactive ligand). Plausible paths
forward are modality-dependent rather than classical inhibition:
degradation/PROTAC of the transcription factor (no VEZF1 example exists,
but STUB1-mediated natural degradation shows the protein is
degradation-controllable), antisense/siRNA (a natural miR-191-VEZF1
regulatory relationship is documented, PMID 31064890), or further
medicinal-chemistry optimization of the He et al. 2018 series.

### EML5 -- CURRENTLY_POORLY_DRUGGABLE
A non-catalytic, WD40/HELP-domain microtubule-associated scaffold protein
(UniProt Q05BV3, 29 annotated WD-repeat features). Pharos classifies EML5
as target-development-level "Tdark" -- its lowest tier, essentially no
functional annotation. **Zero ChEMBL entries, zero PDB structures, no
validated protein interactor, no inhibitor/probe/degrader of any kind
found.** An AlphaFold model does exist with comparatively high confidence
(pLDDT 87.0), so the fold itself is well-predicted -- the limitation is
total absence of pharmacology and validated targetable interfaces, not
structural intractability per se. Family-level context only: EML5's
paralog EML4 is the fusion partner in EML4-ALK (lung cancer), but the
drugged moiety there is ALK's kinase domain, not any EML-family region --
this offers no direct precedent for drugging EML5 itself.

### CITED2 -- INDIRECT_OR_MODALITY_DEPENDENT
An intrinsically disordered transcriptional co-regulator (UniProt Q99967,
270 aa) with **no DNA-binding domain and no catalytic activity** -- it
works primarily by competitively binding the CBP/p300 CH1/TAZ1 surface
(the same surface HIF-1alpha uses); UniProt also lists CITED2 interactions
with SMAD2/3, TFAP2A/B/C, WT1, LHX2 and PPARA, none independently
co-structured. Unlike VEZF1 and EML5, the CBP/p300 interface specifically
**is structurally characterized**: native-complex NMR
structures PDB 1P4Q (PMID 12778114) and 1R8U (PMID 14594809) show a
coupled folding-and-binding interaction over an extended, non-pocket-like
groove. Open Targets' formal small-molecule tractability assessment
returns FALSE across every bucket checked (approved drug, clinical
candidate, structure-with-ligand, high/med-quality pocket, druggable
family) -- consistent with this being a difficult, extended
protein-protein interface rather than a classical active site. Zero
ChEMBL entries, zero registered chemical probes, zero clinical
candidates. The most plausible path is a PPI-interface disruptor
(stapled/macrocyclic peptide mimicking the CITED2 activation-domain
element, for which a real structural template exists) or targeted
degradation -- CITED2 is classified INDIRECT_OR_MODALITY_DEPENDENT rather
than CURRENTLY_POORLY_DRUGGABLE specifically because that structural
template exists, unlike VEZF1/EML5.

**No candidate reaches DIRECTLY_DRUGGABLE.** The four require materially
different intervention strategies -- they are not forced into one
modality.

---

## Part 3 -- normal-tissue / selectivity review

Full detail: `results/tables/druggability_safety/candidate_normal_tissue_context.tsv`.
Figure: `results/figures/druggability_safety/02_four_candidate_normal_tissue_context.png`.

All four candidates show **broad, non-tissue-restricted expression** in
GTEx v8 / HPA, with one partial exception (EML5, see below). This is
expression breadth, not a claim about essentiality or safety in any
tissue.

- **USP34**: expressed in all 54 GTEx v8 tissues (no near-zero tissue);
  HPA calls it "low tissue specificity." Breast 18.6 TPM, blood 4.3 TPM
  (comparative low point), liver 4.1 TPM (lowest GTEx tissue), heart 5.4
  TPM, kidney 6.7 TPM, brain 6.6-27.6 TPM across regions, GI tract
  12.5-24.0 TPM, reproductive tissues 25-29 TPM.
- **VEZF1**: broad, low-specificity expression, <9-fold range across all
  50 GTEx tissues sampled. Breast 30.5 TPM, blood 9.8 TPM, liver 6.7 TPM
  (lowest GTEx tissue), heart 14.1 TPM, kidney 14.1 TPM. HPA notes
  cell-type enrichment in endothelial cells specifically (consistent with
  VEZF1's namesake vascular biology) without being tissue-exclusive.
- **EML5**: the most tissue-restricted of the four -- HPA calls it
  "tissue enhanced (ovary, retina)." Breast 0.27 TPM, blood 0.02 TPM,
  liver 0.03 TPM, heart 0.02 TPM -- near the floor of detection in most
  organs, with a real mRNA/protein discordance noted (ovary mRNA highest,
  but protein "not detected" by HPA immunohistochemistry).
- **CITED2**: broadly and, in absolute terms, often highly expressed
  (492 TPM in fibroblasts, 361 TPM in ovary; breast 133.7 TPM is itself
  among the higher values in the entire GTEx panel). Heart (30.0 TPM) and
  whole blood (26.3 TPM) are comparatively lower but clearly non-zero.
  Not detected in blood plasma (not a circulating protein).

### Bone / musculoskeletal system (explicitly required, not omitted for GTEx's sampling gap)

GTEx does not sample bone directly; HPA's bulk-tissue panel includes bone
marrow and skeletal muscle but **no osteoblast, osteocyte, osteoclast,
chondrocyte, or skeletal mesenchymal-progenitor cell-type entry exists in
HPA's single-cell atlas for any of the four candidates** -- a genuine
data gap, disclosed rather than papered over. Full detail:
`results/tables/druggability_safety/candidate_bone_musculoskeletal_context.tsv`.
Figure: `results/figures/druggability_safety/04_four_candidate_bone_musculoskeletal_context.png`
(expression and function shown as two separate panels, on purpose).

**USP34** -- bone marrow 35.5 nTPM, skeletal muscle 54.1 nTPM (HPA,
3rd-highest of 50 HPA tissues). **Direct, causal published role**: Guo et
al. 2018 (PMID 30181118, *EMBO J*) -- USP34 stabilizes Smad1/RUNX2
(opposing SMURF1-mediated degradation); depletion in human MSCs in vitro
inhibits osteogenic differentiation; conditional (MSC- or
pre-osteoblast-Cre) knockout mice show **low bone mass and impaired
BMP2-driven bone regeneration in vivo**. Separately, Wigoda et al. 2026
(PMID 42315110) reports human heterozygous de novo USP34 loss-of-function
causes a neurodevelopmental syndrome that includes **distal limb
anomalies** -- human, germline, developmental evidence, not an
adult-onset or pharmacological-inhibition finding. Note also that the
mouse knockout itself used MSC/pre-osteoblast-Cre, which is active during
skeletal development/differentiation, not a deletion induced in
already-mature adult bone -- so this is best described as postnatal
-developmental causal evidence, not adult-onset causal evidence (contrast
with CITED2's genuinely adult-induced Mx1-Cre studies below). Evidence
type: DOCUMENTED_CAUSAL_POSTNATAL (mouse, postnatal developmental
conditional knockout) + DOCUMENTED_DEVELOPMENTAL (human, germline).

**VEZF1** -- bone marrow 53.5 nTPM (HPA), but immunohistochemistry is
discordant between two antibodies in the same hematopoietic-cell
population (one "not detected," one "medium") and mass spectrometry found
"not detected" in all 3 replicates -- a genuine cross-method
discrepancy, not a clean positive protein-level finding. **No direct
bone/skeletal paper exists.** The closest relevant finding (Das et al.
2023, PMID 36923254) shows *Vezf1*-knockout mouse embryoid bodies have
reduced hematoendothelial marker genes -- developmental,
lineage-specification evidence upstream of adult bone-marrow-niche
biology, which this review explicitly declines to over-read as direct
bone-marrow-niche evidence. The founding knockout paper (PMID 15882861)
reports embryonic-lethal vascular defects, with **no skeletal phenotype
category** flagged in MGI's curation. Evidence type: INFERRED_ONLY.

**EML5** -- bone marrow 0.3 nTPM, skeletal muscle not detected (0 nTPM).
**No published bone/skeletal role of any kind** -- a PubMed search for
"EML5 AND (bone OR osteoblast OR osteoclast OR skeletal)" returns zero
results, out of only 17 total EML5 papers in all of PubMed. Mouse: MGI
records zero curated phenotypes despite six alleles existing; IMPC has
not yet tested this gene (0/24 systems). Evidence type: NONE_IDENTIFIED.

**CITED2** -- no single bone-marrow/skeletal-muscle bulk-tissue nTPM value
was locatable in the HPA extract used (reported honestly as missing data
in the table/figure rather than omitted or estimated). Despite that gap,
CITED2 has **the most extensive, mechanistically direct, multi-model
bone/marrow literature of the four candidates**, spanning four
independent axes: (1) **hematopoietic stem cell/bone-marrow
maintenance** -- adult conditional (Mx1-Cre) knockout causes HSC loss,
loss of quiescence, and multilineage bone-marrow failure in mice (PMID
19951693, 22308296, 34715054), with human confirmatory data in primary
CD34+ cells and AML xenografts (PMID 25184385); (2) **osteoclast
differentiation** -- CITED2 is described as "the molecular switch
triggering terminal differentiation of osteoclasts," with in vivo
lineage-specific deletion in osteoclast precursors causing failure to
commit to the osteoclast fate (PMID 33288951, *Nat Metab* 2020); (3)
**fracture healing** -- identified as a negative regulator of fracture
healing via MMP suppression in a rat model (PMID 19607804); (4)
**cartilage mechanotransduction** -- mediates MMP-1/MMP-13 repression
under physiological joint loading in human chondrocytes and rat in vivo
models (PMID 12960175, 20826544). Also implicated in
adipogenic-vs-osteogenic lineage commitment of human MSCs (PMID
33691475). Constitutive knockout is separately embryonic lethal
(cardiac/neural-tube, PMID 11694877/12149478/15750185) without a
distinctly reported skeletal phenotype in those specific papers (absence
of report, not a documented negative finding). Evidence type:
DOCUMENTED_ADULT_CAUSAL, the strongest of the four.

**This is a genuinely important, non-obvious finding of this phase**:
CITED2 -- which had the lowest DepMap baseline dependency and the
weakest Hany significance of the four -- carries the most substantiated,
causal, adult/postnatal bone- and bone-marrow-relevant biology in the
literature. This is a real caution specific to the hematopoietic/skeletal
axis, separate from and not implied by its otherwise-favorable
low-cancer-dependency profile.

---

## Part 4 -- genetic constraint (kept strictly separate from expression, essentiality, dependency, and toxicity)

Full detail: `results/tables/druggability_safety/candidate_genetic_constraint.tsv`.

| Candidate | gnomAD LOEUF | pLI | Human phenotype (OMIM) | ClinGen dosage-sensitivity curation | Mouse KO |
|---|---|---|---|---|---|
| USP34 | **0.152** (most constrained of the four) | 1.0 | Gene-only entry; ClinGen "Limited" congenital-heart-disease link; 2026 haploinsufficiency-neurodevelopmental paper (see Part 3) | No formal 0-3 score on record | Conditional KO: reduced bone formation (see Part 3) |
| VEZF1 | 0.24 | 1.0 | Single-family dilated-cardiomyopathy report | **DISPUTED**: ClinGen formally rates "No Known Disease Relationship" (2026-03-04), contradicting the OMIM listing | Homozygous: embryonic lethal (vascular); heterozygous: ~20%-penetrant lymphatic hypervascularization |
| EML5 | 0.558 (least constrained of the four) | ~0 | None (gene-only OMIM entry) | "Awaiting Review" (not yet curated) | Zero curated phenotypes from six alleles |
| CITED2 | 0.17 (point estimate; 90% CI upper bound 0.54, wide due to the gene's small size) | 0.94 | VSD2 (#614431), ASD8 (#614433) -- rare hypomorphic in-frame variants | **DISPUTED**: ClinGen formally rates "No Evidence for Haploinsufficiency" (2023-11-29), citing incomplete penetrance and variants found in both patients and controls | Constitutive KO: embryonic/perinatal lethal (cardiac, neural-tube, adrenal) |

Two genuine OMIM-vs-ClinGen discrepancies were identified (VEZF1
cardiomyopathy, CITED2 congenital heart defects) and are reported as
unresolved tensions, not silently resolved in either direction: OMIM
reflects individual published gene-disease associations, while ClinGen
reflects a later, systematic, expert-panel evidence review -- the two are
not always in agreement, and both are reported.

**USP34 is the most strongly constrained against loss-of-function
variation in the human population of the four candidates (LOEUF=0.152,
pLI=1.0) -- well within gnomAD's "constrained" range**, and this
constraint co-occurs with the only postnatal-causal mouse bone-formation
phenotype and a confirmed human germline haploinsufficiency syndrome
among the four. **Genetic constraint reflects purifying selection
against germline loss-of-function during development and is not
predictive of adult pharmacological-inhibition toxicity** -- USP34's
strong constraint does not establish that a small-molecule inhibitor
would be unsafe in adults, but it is a real, independent flag of strong
purifying selection against germline loss-of-function -- evidence of an
important developmental/reproductive-fitness role, not itself proof of
broad essentiality in adult normal tissue -- and should be weighed
alongside (not instead of) the bone-formation and expression findings
above.

**EML5 is the least constrained of the four (LOEUF=0.558, pLI~0)**,
consistent with its absence of any human disease link -- but this should
be read together with EML5 being the single least-studied and
least-annotated candidate (Pharos "Tdark," only 17 total PubMed papers),
so low constraint may partly reflect genuine tolerance of
loss-of-function and may partly reflect that EML5 has simply not been
extensively assessed. Both are stated; neither is asserted as the
explanation.

---

## Part 5 -- cancer-vs-normal therapeutic-window interpretation, per candidate

Figure: `results/figures/druggability_safety/03_four_candidate_therapeutic_window_map.png`
(real DepMap cancer-dependency data on one axis; an explicitly-labeled,
exploratory qualitative normal-tissue+bone concern index on the other --
never collapsed into a single hidden score).

**USP34.** Low ER+/luminal baseline cancer dependency (0.0%) plus strong
Hany sensitisation (FDR=0.042) does support a relatively
tamoxifen-specific combination hypothesis on the cancer-dependency axis
alone. But normal-tissue biology meaningfully tempers this: USP34 is
broadly expressed across nearly all normal tissues, is the most strongly
LoF-constrained gene of the four, has a documented causal role in mouse
postnatal bone formation, and a documented human developmental
haploinsufficiency syndrome with limb anomalies. USP34/WNT biology
therefore **does create a plausible bone-remodelling/skeletal concern**
for systemic inhibition -- not because DepMap suggests it (it does not;
DepMap is not used for this claim), but because of the independently
sourced USP34-BMP2-osteogenesis literature and the human genetics
findings above.

**VEZF1.** The 27.3% baseline ER+/luminal dependency is a real cancer
vulnerability signal, distinct from and additional to the tamoxifen
sensitisation signal -- this is the strongest basis among the four
candidates for a genuine dual-action-target hypothesis (see Part 6). Is
there also evidence of an important normal-tissue developmental role for
VEZF1? Yes, though this falls short of demonstrating adult normal-tissue
essentiality: VEZF1 is similarly strongly constrained (LOEUF=0.24,
pLI=1.0) and homozygous loss is embryonic lethal in mice (vascular) --
both consistent with an important developmental role, neither a direct
measurement of adult-tissue essentiality. There is **no direct
bone/skeletal paper** and the one human disease association
(cardiomyopathy) is formally disputed by ClinGen. VEZF1's vascular/
developmental biology creates a **plausible but only inferred**
bone-marrow-vasculature concern (via the hematoendothelial-lineage
literature, PMID 36923254) -- weaker and less direct than USP34's or
CITED2's bone evidence, and explicitly not extrapolated further than the
source supports.

**CITED2.** Its strong transcriptional/co-regulatory biology (an
intrinsically disordered protein working through an extended
protein-protein interface) does make direct small-molecule inhibition
difficult -- Open Targets' formal tractability assessment agrees (all
buckets false). Indirect modulation -- a PPI-interface disruptor
exploiting the real 1P4Q/1R8U structural template, or targeted
degradation -- is more plausible than classical inhibition. Important
normal developmental/stress-response functions do constrain systemic
inhibition: constitutive knockout is embryonic lethal (cardiac/neural
tube) in mice, and -- more relevant to an adult patient receiving a
systemic drug -- CITED2 has **documented, causal, postnatal roles in
hematopoietic stem cell maintenance and osteoclast differentiation** (see
Part 3). CITED2 is, unambiguously, involved in bone/marrow/cartilage
biology: this is the single most load-bearing normal-tissue finding of
this entire review.

**EML5.** Its clean DepMap profile (0.0% baseline dependency) does not
translate into a realistic targetability advantage -- the dominant
limitation is a near-total absence of mechanism, structure, and
drug-binding opportunity (Pharos "Tdark," zero ChEMBL entries, zero PDB
structures, zero validated interactors), not normal-tissue safety data
supporting it. There is no credible bone/musculoskeletal biology of any
kind for EML5 (Part 3) -- this changes its interpretation only by
removing a normal-tissue concern that was never established one way or
the other for the other three candidates as clearly as it is absent
here, not by adding a positive finding.

---

## Part 6 -- "better than sensitisation" signals

Every candidate was assessed against five categories: (A) pure tamoxifen
sensitiser, (B) direct ER+ cancer vulnerability, (C) dual-action target
(direct cancer-cell impairment + tamoxifen sensitisation), (D)
context-specific synthetic-lethal-type target, (E) mechanistically
interesting but poorly druggable.

**VEZF1 is the clearest candidate in category C.** It combines a real,
non-trivial baseline ER+/luminal cancer-cell dependency (27.3% strongly
dependent, DepMap 26Q1, B_MODERATE tier) with a significant Hany
tamoxifen-sensitisation signal (FDR=0.037) -- these are two additive,
mechanistically distinct signals, not one signal counted twice. This
supports the project's existing "potential dual-action hypothesis"
framing and, on the evidence gathered in this phase, is not undercut by
any normal-tissue finding strong enough to override it (VEZF1's
normal-tissue evidence is real but comparatively the weakest-documented
bone-specific concern of the two candidates with any positive bone
signal at all -- see Part 3). VEZF1 remains category E as well
(mechanistically interesting, poorly druggable) simultaneously -- these
categories are not mutually exclusive, and VEZF1's biological interest
(category C) is not matched by its current druggability (Part 2).

No other candidate clears category B or C on the evidence gathered here:
USP34 and EML5/CITED2 all show 0.0% baseline ER+/luminal dependency by
DepMap (category A/E territory, not B/C). CITED2's strongest signal in
this phase is its bone/marrow mechanistic biology, which is a normal
-tissue finding, not a cancer-vulnerability finding, and should not be
mistaken for one.

---

## Part 7 -- separate evidence dimensions (no fake master score)

Full detail (all ten dimensions, per candidate):
`results/tables/druggability_safety/candidate_therapeutic_window_summary.tsv`.
This table is explicitly a join of independent dimensions -- functional
tamoxifen evidence, cancer dependency, mechanism strength, human-tumor
support, direct druggability, alternative-modality feasibility,
normal-tissue concern, bone/musculoskeletal concern, genetic-constraint
concern, and novelty -- kept as separate columns throughout. The frozen
therapeutic ranking, the frozen mechanistic follow-up ranking (three-way
tie: USP34/CITED2/EML5 rank 1, VEZF1 rank 4), this phase's druggability
assessment, and this phase's normal-tissue assessment remain four
distinct concepts and are never collapsed into one number in any table in
this project. Figure 03's y-axis "concern index" is the one place a sum
appears anywhere in this phase, and it is explicitly labeled exploratory,
restricted to two qualitative tiers, and never presented as a ranking
criterion.

---

## Part 8 -- final candidate recommendation

**A. Best candidate for "tamoxifen-specific sensitiser": USP34.** Lowest
baseline cancer dependency (0.0%) combined with the more favorable of the
two significant Hany signals in terms of a clean D_LOW DepMap profile
(VEZF1's Hany FDR is nominally lower but comes with a real baseline
dependency that argues against a *purely* sensitiser mechanism).

**B. Best candidate for "potential dual-action cancer target": VEZF1.**
The only candidate combining a real baseline ER+/luminal dependency with
a significant Hany signal (Part 6).

**C. Best candidate for "strongest mechanism / biological rationale":
CITED2**, on two fronts simultaneously: the deepest breast-cancer/ER
literature base of the four (per the prior literature-mechanism phase)
*and*, newly in this phase, the most extensive, causal, multi-model
bone/marrow biology -- though this second front is a normal-tissue
caution, not a cancer-mechanism strength, and the two should not be
conflated when citing this candidate's "strongest mechanism."

**D. Best candidate for "lowest apparent baseline cancer dependency":**
a three-way tie -- USP34, EML5, and CITED2 all show 0.0% ER+/luminal
strongly-dependent fraction in DepMap 26Q1.

**E. Best candidate for "most realistically druggable": USP34** -- the
only candidate with a real catalytic-domain structure, a defined
enzymatic mechanism, and a documented two-state (inactive/active)
conformational feature that could support future inhibitor design, even
though no USP34-specific chemical matter exists today.

**F. Candidate with "lowest apparent normal-tissue concern" (by available
data): EML5** -- comparatively tissue-restricted expression, the least
genetically constrained of the four, no established human disease
phenotype, and zero bone/skeletal literature of any kind. This is
explicitly **not** a safety claim: EML5 is also the least-studied
candidate by a wide margin (Pharos "Tdark," 17 total PubMed papers), so
part of its "lowest apparent concern" reflects data sparsity rather than
a demonstrated benign profile. Absence of evidence is disclosed as such,
not read as evidence of absence.

**G. Candidate with "highest potential bone/musculoskeletal concern":
CITED2**, ahead of USP34. CITED2 has the most extensive causal,
multi-model (human + mouse + rat) evidence across four independent
bone/marrow/cartilage axes (Part 3), and -- unlike USP34's finding below
-- its HSC-maintenance studies used a conditional-deletion system induced
in already-mature adult mice (Mx1-Cre), not a developmental promoter.
USP34 is second, with a directly causal mouse conditional-knockout
bone-formation phenotype (using a developmental MSC/pre-osteoblast-Cre
promoter, active during skeletal differentiation rather than in
already-formed adult bone) plus a human developmental limb-anomaly
finding. Both are materially higher-concern than VEZF1 (inferred only)
and EML5 (no evidence identified).

**H. Best overall 1-2 candidates for experimental follow-up: USP34 and
VEZF1.** USP34 offers the cleanest tamoxifen-specific mechanistic
hypothesis and the only real structural handle for future drug
discovery, tempered by a genuine, newly-identified bone/genetic
-constraint caution that should inform any in vivo dosing or duration
design (e.g., monitoring bone parameters) rather than halt exploration.
VEZF1 offers the most therapeutically interesting biological hypothesis
(dual-action target, Part 6) despite currently poor druggability, making
it the stronger near-term target for mechanistic/functional follow-up
work (e.g., confirming the dual-action hypothesis with orthogonal assays)
rather than immediate medicinal chemistry. CITED2, despite its strong
mechanistic case, is de-prioritized for *this specific* follow-up
recommendation because of its newly-identified, uniquely extensive
preclinical (mouse/rat/human-ex-vivo) evidence of a causal bone-marrow/
skeletal role -- a genuine caution for further development and a real
trade-off worth weighing, not an established human toxicity finding and
not a re-ranking of the frozen shortlist, which this phase does not
alter.

---

## Quality-control checklist

1. Frozen four-candidate shortlist and both frozen rankings unchanged --
   confirmed by direct read of `THERAPEUTIC_SHORTLIST_FREEZE.tsv` and
   `four_candidate_followup_rankings.tsv`, never rewritten.
2. No candidate called "safe" anywhere in this report or its tables --
   confirmed by direct text search of this file (zero occurrences of the
   word "safe" outside this checklist item and the explicit prohibition
   restated at the top).
3. DepMap never used as normal-tissue evidence -- every DepMap reference
   in Parts 3/4/6 is explicitly scoped to cancer-cell-line dependency; the
   genetic-constraint and bone tables use gnomAD/OMIM/ClinGen/MGI/IMPC/
   GTEx/HPA/PubMed only.
4. High expression never equated with essentiality (Figure 02's own
   subtitle states this explicitly); genetic constraint never equated
   with drug toxicity (Part 4's explicit caveat paragraph).
5. Every PMID/accession/numeric value in this report traces to
   `results/tables/druggability_safety/*.tsv`, which in turn trace to
   `verified_references.tsv` and the source URLs recorded in each table's
   `sources` column -- nothing in this report was hand-typed independent
   of those tables.
6. Bone/musculoskeletal system explicitly assessed for all four
   candidates despite GTEx's sampling gap, using HPA/PubMed/MGI/IMPC as
   the user required; missing data (e.g., CITED2's bone-marrow bulk
   nTPM) is shown as an explicit gap in both the table and Figure 04, not
   silently omitted or estimated.
7. Mouse/developmental evidence never extrapolated to expected human
   adult toxicity (explicit caveats on USP34's and CITED2's embryonic
   -lethal/developmental findings throughout Part 3).
8. Evidence dimensions (druggability, normal-tissue, bone,
   genetic-constraint, mechanism, cancer-dependency) kept in separate
   table columns throughout; the one numeric sum in this phase (Figure
   03's y-axis) is explicitly labeled exploratory and is not used as a
   ranking criterion anywhere in Part 8.
9. Two genuine OMIM-vs-ClinGen discrepancies (VEZF1, CITED2) reported as
   unresolved tensions, not silently resolved.
10. All four candidates classified into different druggability tiers
    where the evidence actually differs (POTENTIALLY_DRUGGABLE,
    CURRENTLY_POORLY_DRUGGABLE x2, INDIRECT_OR_MODALITY_DEPENDENT) -- not
    forced into one intervention type.

---

## Outputs

**Tables** (`results/tables/druggability_safety/`): `candidate_druggability.tsv`,
`candidate_normal_tissue_context.tsv`, `candidate_genetic_constraint.tsv`,
`candidate_bone_musculoskeletal_context.tsv`,
`candidate_therapeutic_window_summary.tsv`, `verified_references.tsv`.

**Figures** (`results/figures/druggability_safety/`):
`01_four_candidate_druggability_summary.png`,
`02_four_candidate_normal_tissue_context.png`,
`03_four_candidate_therapeutic_window_map.png`,
`04_four_candidate_bone_musculoskeletal_context.png`.

**Code**: `src/druggability_safety_data.py` (curated, source-cited dossier
data), `src/druggability_safety_build_tables.py` (deterministic table
builder, joins this phase's curated data with frozen upstream files),
`src/druggability_safety_visualization.py` (figure builder).
