# USP34 vs VEZF1 -- full translational deep-dive

USP34 and VEZF1 are the two lead candidates selected for deeper
follow-up. **CITED2 and EML5 are not disproven** -- they remain frozen
shortlist members, simply not the subject of this phase. The frozen
four-candidate therapeutic ranking, frozen TCGA results, and frozen
DepMap 26Q1 results are all unchanged and unrerun.

## A note on how this phase was actually conducted

This phase was originally scoped as four parallel background research
passes (USP34 tissue/muscle/genetics, VEZF1 tissue/muscle/genetics, USP34
direct/indirect targeting, VEZF1 direct/indirect targeting). **All four
failed mid-run** because the session hit its weekly API usage limit
before any of them produced results -- no data from those failed runs is
used anywhere in this report or its tables. Rather than fabricate the
missing coverage, the work was instead completed by (a) carrying forward
every already-verified fact from the prior `druggability_safety` phase
(re-cited, never re-derived or altered) and (b) a smaller, still fully
source-verified, direct research pass in this same session -- live STRING
interactome queries, targeted PubMed/NCBI-eutils searches, ClinVar
counts, and re-verification fetches -- prioritized toward this phase's
two highest-value new questions: **systematic indirect-target discovery**
and a **muscle-specific literature search** for both genes. Coverage is
real but narrower than the full 18-section brief requested (for example,
per-cell-type single-cell atlas queries beyond HPA, and an exhaustive
per-organ PubMed sweep, were not completed for every organ system). Every
place this narrower scope applies is flagged explicitly in the tables
(`data_status`, `sources`, or `evidence_summary` columns) rather than
silently gapfilled.

**Gap-completion and validation pass (this revision)**: the user asked
for a focused follow-up targeting five specific gaps left by the pass
above: (1) cell-type-resolved tissue-atlas coverage beyond HPA's
bulk/immune panel, (2) a more systematic USP34 muscle search including
historical gene aliases, (3) precise characterization of the VEZF1
cardiac finding, (4) critical validation of whether TEAD1 is a genuinely
*validated* VEZF1 indirect target or only a *hypothesis*, and (5) exact
TEAD1 pharmacology (TEAD1-specific vs pan-TEAD, true clinical status).
This revision addresses all five directly, and **materially downgrades
the TEAD1-VEZF1 claim** from how the prior revision framed it -- see Part
9. No frozen upstream evidence was touched; no analysis already
adequately completed was redone.

---

## 0. Preflight consistency check

Checked whether any druggability_safety report/table/figure asserts the
raw gnomAD LOEUF ordering "USP34 > VEZF1 > CITED2 > EML5" (which would
incorrectly rank VEZF1 as more constrained than CITED2 -- the real values
are USP34=0.152, CITED2=0.17, VEZF1=0.24, EML5=0.558, so CITED2 is
actually more constrained than VEZF1). **No such claim exists anywhere**:
the genetic-constraint table lists candidates in fixed
USP34/VEZF1/EML5/CITED2 candidate order (not sorted by LOEUF), and the
report's prose only ever asserts the two unambiguously-true endpoint
statements ("USP34 is the most strongly constrained," "EML5 is the least
constrained"), never a full four-way ordering that would misplace VEZF1
relative to CITED2. **No correction was required.**

---

## 1. Tissue expression highlights

Full table: `results/tables/lead_target_deep_dive/USP34_VEZF1_full_tissue_expression.tsv`.
Figure: `01_USP34_VEZF1_tissue_expression_atlas.png`.

Both genes are broadly, similarly expressed across nearly all normal
tissues (both "low tissue specificity" per HPA) -- neither is
tissue-restricted. VEZF1 runs numerically higher than USP34 in most of
the tissues compared (breast, blood, liver, heart, skeletal muscle, bone
marrow), though both are well within the range of broadly-expressed
genes generally. **Expression breadth alone is not treated as a liability
claim anywhere in this report.**

**Gap-completion additions**: HPA's single-cell-type ENRICHMENT calls
(distinct from, and more informative than, the bulk numbers above) were
newly queried this pass via HPA's JSON API for both genes. **USP34** is
enriched in Sertoli cells (786.9 nCPM, its single highest cell type of
any kind), adrenal cortex cells, and testicular spermatocytes/
spermatogonia; its brain snRNA-seq data shows the highest single region
as choroid plexus (68.1 nTPM). **VEZF1** is enriched in adrenal cortex
cells, skin and splenic endothelial cells (consistent with, though not
exclusive to, its namesake vascular biology), and testicular
spermatogonia; its brain snRNA-seq cluster is "white matter - signal
transduction," consistent with its corpus-callosum enrichment. **Neither
gene shows an enrichment call in any muscle, cardiac, neuronal-subtype,
or bone-lineage cell type** -- this null result is itself informative and
consistent with the expression-only/no-cell-type-specificity picture
built up elsewhere in this report.

A genuine data gap, still disclosed rather than papered over: a true
per-cell-type numeric TABLE (not just the enrichment/highest-value
summary calls above) could not be retrieved this pass -- HPA's
interactive single-cell visualization is JavaScript-rendered and not
exposed via the static JSON/HTML endpoints available to this pass's
tooling. Bone-lineage cell types (osteoblasts/osteocytes/osteoclasts/
skeletal-MSCs) remain absent from HPA specifically; this pass confirmed
that at least three dedicated human bone/bone-marrow single-cell atlases
exist in the literature (a pediatric bone-marrow multiomic atlas, a
673,750-cell bone-marrow reference atlas with annotated osteoclasts, and
a human primary-osteoblast single-cell dissection), but none could be
queried for a per-gene USP34/VEZF1 value with this pass's tooling -- they
are primary-paper supplementary datasets, not simple gene-lookup portals.
This is an explicit, disclosed limitation, not a claim that the data
doesn't exist anywhere.

---

## 2. Full normal-tissue liability analysis

Full table: `results/tables/lead_target_deep_dive/USP34_VEZF1_tissue_liability.tsv`.
Figure: `02_USP34_VEZF1_liability_map.png`.

**USP34's strongest liabilities**: Bone (`DOCUMENTED_CAUSAL_POSTNATAL` --
human MSC knockdown + mouse conditional-KO both impair osteogenic
differentiation/bone formation) and Neurological/CNS
(`DOCUMENTED_DEVELOPMENTAL` -- a 2026 human paper reporting heterozygous
de novo loss-of-function causes a neurodevelopmental syndrome). Skeletal
muscle and marrow are `EXPRESSION_ONLY` -- expressed, but zero functional
evidence found. Cardiovascular is `INSUFFICIENT_DATA` (not deeply searched
this pass).

**VEZF1's strongest liability**: Cardiovascular
(`DOCUMENTED_POSTNATAL_CAUSAL` -- see Part 6 below, this is the single
most important new finding of this phase). Vascular/endothelial
(embryonic) is `DOCUMENTED_DEVELOPMENTAL` (embryonic-lethal homozygous
mouse knockout) -- this embryonic-lethal finding must not be extrapolated
to adult drug toxicity; it is developmental-mouse evidence, never read
here as adult drug toxicity. Bone/marrow is `INFERRED_ONLY` -- no direct
paper. Skeletal muscle is `NONE_IDENTIFIED`.

Every row in the liability table explicitly separates: is the gene merely
expressed there (A)? human genetic evidence of reduced-dosage effect (B)?
primary-human-cell perturbation evidence (C)? adult-conditional vs
developmental-only animal evidence (D/E)? reversible or
developmental/irreversible (G)? heterozygous or homozygous (H)? relevant
to partial drug inhibition or only complete germline knockout (I)? -- per
the user's exact question list.

---

## 3. Special muscle liability deep-dive

Full table: `results/tables/lead_target_deep_dive/USP34_VEZF1_muscle_liability.tsv`.

**USP34 -- VERDICT D: no relevant evidence identified in the searches
performed** (of the four options A-D the user specified: A documented
functional liability, B suggestive functional evidence, C expression-only
concern, D no relevant evidence identified). HPA reports high
skeletal-muscle expression (54.1 nTPM bulk; "myonuclei" 766.5 nCPM,
2nd-highest of 154 HPA single-cell types, carried forward from the prior
phase). A live PubMed search this pass for USP34 combined with skeletal
muscle / myocyte / myoblast / myotube / satellite cell / muscle
differentiation / muscle regeneration / sarcomere / muscle weakness /
myopathy / muscular dystrophy / exercise / muscle atrophy returned **zero
hits** (NCBI esearch, confirmed live 2026-08-14). **Gap-completion**: the
search was repeated substituting USP34's two historical gene aliases,
**KIAA0570 and KIAA0729** (confirmed as genuine USP34 aliases via
Sigma-Aldrich's antibody catalog cross-reference), against the same
muscle terms -- again **zero hits**. This is now a twice-confirmed,
reproducible null result. Separately, this pass identified a real, related
(but non-muscle) developmental finding: USP34 deubiquitinates and
stabilizes NFIC to drive tooth-root morphogenesis, in a mouse conditional
(Sp7-Cre) knockout AND human dental pulp cells (PMID 33686052, *Int J
Oral Sci* 2021) -- a THIRD documented USP34 postnatal
mesenchymal-differentiation role (alongside bone/BMP2 and the human
limb-anomaly syndrome), reinforcing that USP34's genuine developmental
theme is osteogenic/odontogenic mesenchymal differentiation, never
muscle. One AI-summarized web-search result claiming a C2C12-myoblast
differential-expression finding for Usp34 could **not** be independently
re-confirmed by a direct follow-up search and is explicitly excluded
here, not cited as evidence. **Conclusion: HIGH EXPRESSION with NO
FUNCTIONAL REQUIREMENT evidence and NO ACTUAL TOXICITY evidence.** These
three are explicitly distinguished, per the user's instruction, and must
not be conflated.

**VEZF1 -- skeletal muscle: VERDICT D (no relevant evidence identified);
cardiac muscle: VERDICT A (documented functional liability), fully
characterized this pass.** The skeletal-muscle search returned 8 raw
hits; each was individually reviewed. None describe a genuine
**skeletal**-muscle function -- one is an unrelated "Numb" gene paper, one
is a rare tumor gene-fusion case report, the rest are cardiac/vascular
papers. **Gap-completion, full characterization of the cardiac finding**
(PMID 31911272, *Vezf1 regulates cardiac structure and contractile
function*, EBioMedicine 2020, re-verified via three independent routes --
PubMed, the EBioMedicine listing, and an independent web search, all
agreeing): **Organism** -- zebrafish (knockdown, causal) + primary
cardiomyocytes (mechanistic) + human myocardium/mouse hearts
(expression-association only, correlative, never perturbed). **Context**
-- explicitly a POSTNATAL model of compensatory cardiac growth/stress
response, not embryonic patterning. **Perturbation** -- morpholino
knockdown (a partial loss-of-function, not a genetic null). **Phenotype**
-- reduced cardiac growth and impaired ventricular CONTRACTILE response
to beta-adrenergic stimuli (a functional, not merely structural,
readout); Ca2+-transient kinetics were explicitly unaffected, ruling out
a generic excitation-contraction defect. **Mechanism** -- Myh7/beta-MHC,
a core cardiac SARCOMERIC gene, via an MCAT promoter element, with TEAD-1
reported as a binding partner at that site (see Part 9 for how this is
scoped for indirect-targeting purposes -- it is weaker evidence than the
prior version of this report implied). **Scope** -- this is cardiac
muscle/cardiomyocyte biology specifically, distinct from VEZF1's separate
embryonic-lethal VASCULAR knockout finding (PMID 15882861, below) and
from skeletal muscle (VEZF1 has none). **Relevance and limitations for
adult pharmacologic VEZF1 inhibition**: a PARTIAL knockdown causally
impaired a POSTNATAL functional readout, which is more directly relevant
to a partial pharmacological inhibitor than an embryonic-lethal knockout
would be -- but this is zebrafish, not human, and must not be read as an
expected human drug-toxicity finding. It is the single most
drug-relevant liability signal identified for either gene in this
report.

---

## 4. Human genetic liability

Full table: `results/tables/lead_target_deep_dive/USP34_VEZF1_human_genetic_constraint.tsv`.
Figure: `04_USP34_VEZF1_therapeutic_window_comparison.png` (bottom row).

| | USP34 | VEZF1 |
|---|---|---|
| gnomAD LOEUF | 0.152 (more constrained) | 0.24 |
| pLI | 1.0 | 1.0 |
| ClinVar total variants | 620 | 86 |
| ClinVar pathogenic/likely-pathogenic | 29 | 17 |
| OMIM/ClinGen | Gene-only entry; ClinGen "Limited" heart-disease link; disputed phenotype-map status | Disease link (dilated cardiomyopathy) disputed by ClinGen ("No Known Disease Relationship," 2026-03-04) |
| Heterozygous human phenotype | Yes -- 2026 neurodevelopmental syndrome (PMID 42315110) | Disputed single-family cardiomyopathy report |
| Homozygous phenotype | no relevant evidence identified | Mouse: embryonic lethal (vascular) |

Both genes are strongly constrained against germline loss-of-function.
**This is interpreted, per the user's explicit instruction, as "human
genetics suggests reduced gene dosage may be poorly tolerated / requires
caution" -- never as "constrained = the drug will be toxic."** ClinVar
totals (620 vs 86) were queried live this pass but not individually
reviewed variant-by-variant within this pass's time budget; only the
pathogenic/likely-pathogenic filtered counts were confirmed.

---

## 5. Normal-cell functional essentiality beyond DepMap

No relevant evidence identified in the searches performed for either
gene -- no published CRISPR screen, Perturb-seq, or functional-genomics
dataset in non-cancer human cells (primary cells, iPSC-derived cells,
organoids, normal mammary epithelial cells, endothelial cells, muscle
cells, hematopoietic cells) reporting on USP34 or VEZF1 was located this
pass or in the prior phase. This is reported as an honest gap, not
inferred one way or the other.

---

## 6. USP34 direct targeting

Full table: `results/tables/lead_target_deep_dive/USP34_direct_targeting.tsv`.
Figure: `03_USP34_VEZF1_direct_indirect_targetability.png` (left panel).

USP34 is the only one of the two genes with real experimental structures
(PDB 7W3R apo, 7W3U +ubiquitin-propargylamide probe, both PMID 35588869)
and a defined two-state catalytic mechanism (apo-inactive, ubiquitin
-bound-active; catalytic dyad Cys1903/His2164). No additional structures
beyond these two were found this pass. Verdicts on the six requested
modalities:

- **(A) Direct competitive inhibition**: plausible in principle (real
  catalytic-cleft structures exist), but zero USP34-specific chemical
  matter exists in ChEMBL/BindingDB/PubChem.
- **(B) Covalent inhibition**: plausible in principle -- Cys1903 is a
  real, structurally characterized reactive cysteine (the DUB's own
  catalytic nucleophile) -- but no covalent lead compound exists; whether
  any other ligandable cysteine exists was not established.
- **(C) Allosteric inhibition**: plausible in principle, given the real
  documented apo-inactive/ubiquitin-bound-active conformational
  transition, but untested.
- **(D/E) Targeted degradation / PROTAC**: speculative, contingent on a
  binder that does not yet exist.
- **(F) Molecular glue**: speculative, no evidence found.

**Verdict: no USP34 inhibitor is invented here.** USP34 has the stronger
structural handle of the two genes in this phase but zero validated
chemical matter of any kind.

---

## 7. USP34 indirect targeting

Full table: `results/tables/lead_target_deep_dive/USP34_indirect_targets.tsv`.

STRING's top USP34 interactors beyond generic ubiquitin-pathway proteins
(UBC/UBA52/RPS27A, present in nearly every DUB's interactome and not
USP34-specific) are **AXIN1** (0.651, already-known substrate) and
**UBQLN1** (0.618, database/computational only). A targeted PubMed search
for USP34 combined with "required for" / "necessary for" / CRISPR /
synthetic lethal returned only 2 hits, neither a genuine new regulator.

- **AXIN1**: the causal direction runs the wrong way -- USP34
  deubiquitinates AXIN1 (USP34 acts ON AXIN1), not the reverse -- so this
  is not a viable indirect-inhibition strategy for USP34 itself.
- **UBQLN1**: correlation/database-level only (STRING text-mining-heavy
  score); the primary USP34-AXIN1 mechanism paper does not mention it at
  all (checked directly); also shows a tolerance-associated (wrong)
  direction in this project's own frozen CRISPR data.

**Gap-completion: a deeper, wider STRING query (required_score=400,
limit=30) and an expanded PubMed sweep** (USP34 + HSP90/chaperone/
kinase/phosphorylation/"protein stability," 14 hits, all describing USP34
acting on downstream substrates -- PIN1, SOX2, FOXC1-via-BPTF, PAR1 --
never an upstream regulator of USP34 itself) surfaced a richer
interactome, of which two are worth naming explicitly and neither is
promoted as a real lead:

- **CSNK1A1** (casein kinase 1 alpha): the highest experimental-evidence
  STRING sub-score of any USP34 interactor found (escore=0.607, versus
  0.329 for AXIN1 and 0.554 for UBQLN1) and independently druggable via a
  clinically-validated mechanism (CRBN-mediated CK1-alpha degradation is
  how lenalidomide works in del(5q) MDS) -- but **zero literature**
  connects it functionally to USP34. Surfaced as a screening hypothesis
  only.
- **BPTF**: PMID 38686761 (glioma, U251 cells) reports BPTF, FOXC1, and
  USP34 in one pathway (BPTF knockdown -> less FOXC1 -> mediated by
  USP34-dependent deubiquitination), with real genetic perturbation
  evidence for the BPTF-FOXC1 link -- but whether BPTF is genuinely
  *upstream* of USP34 (which would make it a real indirect target) or
  USP34 simply acts downstream of a BPTF-driven process is **not resolved
  by the available abstract**, and the context is glioma, not breast.
  Ambiguous, not promoted.
- (Tankyrase/RNF146, mechanistically adjacent to the AXIN1 axis, was also
  checked and ruled out -- it stabilizes AXIN1 in the SAME direction as
  USP34, not an opposing one, so it cannot substitute for or reduce
  USP34's effect.)

**USP34's indirect-targeting landscape remains genuinely thin after this
deeper search** -- reported honestly rather than padded. No candidate
identified in either pass meets the bar of genuine perturbation evidence
+ druggable + relevant. CSNK1A1 and BPTF are the most structurally
plausible leads for a future dedicated look, but neither currently has
confirmed functional evidence linking it to USP34.

---

## 8. VEZF1 direct targeting

Full table: `results/tables/lead_target_deep_dive/VEZF1_direct_targeting.tsv`.
Figure: `03_USP34_VEZF1_direct_indirect_targetability.png` (left panel).

Re-verified this pass: the He et al. 2018 (PMID 29970794) tool-compound
figures are exactly as previously recorded (T4 IC50=20uM, T6=100uM,
NSC1012=500uM, homology-model-guided docking hits blocking VEZF1-DNA
binding) -- no follow-up medicinal chemistry or co-crystal structure was
found citing this paper. No covalent-targeting evidence exists. IMiD
-class molecular glues are validated for degrading OTHER zinc-finger TFs
(class-level precedent only, not VEZF1-specific). No VEZF1 PROTAC exists
(contingent on a binder that isn't there). miR-191 naturally targets
VEZF1 mRNA (PMID 31064890) -- proof the transcript is a druggable RNA
node, not a demonstrated therapeutic ASO. **New this pass**: two real
physical-interaction leads were identified -- p68RacGAP (PMID 14966113,
co-IP-based) and TEAD-1 (PMID 31911272) -- neither yet exploited as a
PPI-disruptor strategy; see Part 9 for a critical, downgraded assessment
of the TEAD-1 claim specifically.

**Verdict: direct small-molecule targeting of VEZF1 remains poorly
feasible near-term.** The most realistic near-term paths are RNA-targeting
(building on the miR-191 precedent) or the indirect route below.

---

## 9. VEZF1 indirect targeting -- TEAD1 critical validation (rewritten this pass)

Full table: `results/tables/lead_target_deep_dive/VEZF1_indirect_targets.tsv`
and `indirect_target_project_crosscheck.tsv`.
Figure: `03_USP34_VEZF1_direct_indirect_targetability.png` (right panel).

**The prior version of this report overstated TEAD1 as "the strongest
indirect-targeting lead" in a way that could be read as implying it was
validated. It is not. This section replaces that framing with a
critical, evidence-type-by-evidence-type assessment, per the user's exact
request.**

**Is there evidence that TEAD1 perturbation changes VEZF1
expression/activity? NO -- corrected and precisely characterized below.**
Classifying the available evidence against the six-category scheme (A
direct biochemical interaction, B co-occupancy, C transcriptional
regulation, D genetic perturbation, E drug perturbation, F correlation
only): PMID 31911272's own wording -- "TEAD-1 is a binding partner of
Vezf1" at a shared MCAT element in the Myh7 promoter -- places this claim
at (A)/(B), a direct-interaction or co-occupancy report. **Correction
after full-text verification** (an initial pass of this report incorrectly
stated TEAD1 was never perturbed at all; Codex review located the full
text via PMC6948172 and found this was wrong in a specific way): TEAD1
**was** experimentally perturbed in the paper -- neonatal/adult rat
ventricular cardiomyocytes were transfected with two distinct TEAD1
siRNAs, and TEAD1 knockdown was used in Figure 6F **solely to confirm the
specificity of the ~40-kDa TEAD1 band in a co-immunoprecipitation
experiment** (an antibody-validation control). The paper does **not**
report measuring VEZF1 expression, VEZF1 activity, or the Myh7-promoter
response after TEAD1 knockdown -- category (D) genetic perturbation of
TEAD1 did occur, but was never used to test whether TEAD1 is required for
VEZF1's regulatory effect. Category (E) drug perturbation was not
attempted. There is still no reported experiment measuring a
TEAD1-perturbation-to-VEZF1-output relationship in either direction.

**Is this evidence specific to breast cancer, ER+ biology, or mammary
epithelium? NO.** A dedicated PubMed search this pass for
`(TEAD1 AND VEZF1)` returns **exactly one result in the entirety of
PubMed: PMID 31911272 itself** (NCBI esearch, live query, confirmed
2026-08-14, count=1). No other paper, in any species, tissue, or cancer
context, has ever studied TEAD1 and VEZF1 together. Separately, TEAD1
does have its own independent breast/ER-relevant literature (39 PubMed
hits for TEAD1 + breast/ER/tamoxifen), including two genuinely relevant
background papers -- PMID 42050304 (TEAD1 binds the ERalpha promoter, in
an *osteogenesis* context) and PMID 40500464 (estrogen-receptor
activation remodels TEAD1 expression, in a *hepatic-steatosis* context)
-- but **neither involves VEZF1**, so they establish only that TEAD1 has
*some* independent relevance to estrogen-receptor biology in other
tissues, not that the VEZF1-TEAD1 relationship itself is
breast-cancer-relevant.

**TEAD1 pharmacology, corrected (the prior version overstated this)**:
the actively-progressing clinical compound, **VT3989**, is a **PAN-TEAD**
auto-palmitoylation-pocket inhibitor -- it blocks TEAD1-4 collectively,
not TEAD1 selectively -- and is in real Phase 1/2 trials with FDA Fast
Track designation for mesothelioma/NF2-mutated solid tumors. A
genuinely **TEAD1-selective** compound, **IK-930**, did exist but its
clinical program was **discontinued** after a review showed only modest
clinical activity. No compound with confirmed TEAD1-only selectivity
remains in active clinical development as of the searches performed, and
no breast-cancer or ER+ clinical evidence for either compound was
identified. TEAD1 does have real ChEMBL target entries confirming
chemical tractability in principle (CHEMBL3334416 direct TEAD1;
CHEMBL3430909 YAP1/TEAD1 PPI; CHEMBL6066051 cereblon/TEAD1 PPI, an active
molecular-glue-degrader chemotype; CHEMBL6195504 VHL/TEAD1 PPI).

**What remains independently true and unchanged**: TEAD1 is
independently significant in this project's own frozen resistance-RNA
datasets -- GSE118713 FDR=0.0151, GSE240112 tumor FDR=0.0464 (both
down-in-resistance, direction agreement=True); see Part 10 for the full
cross-check including GSE111151 and GSE245601, which are not
significant. This RNA-level signal is real and independent of the
VEZF1-TEAD1 biochemical claim -- it does **not**, by itself, confirm that
TEAD1 regulates VEZF1.

**Verdict (answers the user's exact question): TEAD1 is currently a
druggable (pan-TEAD, not TEAD1-specific), project-evidence-aligned
HYPOTHESIS worth testing -- it is NOT a validated indirect-targeting
strategy.** It remains the most interesting available VEZF1
indirect-targeting lead only because every other candidate examined is
weaker still (see below), not because this lead is strong in absolute
terms.

Two other candidates were checked and are weaker: **p68RacGAP**
(possibly, but not confirmed, related to ARHGAP22) is a real physical
interactor but without a perturbation-of-X-reduces-VEZF1-output
experiment, and has no known inhibitor. **STUB1** is the wrong direction
entirely -- it degrades VEZF1, so inhibiting STUB1 would be predicted to
*increase*, not decrease, VEZF1.

**None of this adds TEAD1 to the frozen therapeutic shortlist.** It
documents that an independently-identified indirect-target hypothesis
happens to align with part of this project's own pre-existing evidence,
which is exactly the kind of finding Part 10 below exists to surface --
not a claim that the hypothesis is confirmed.

---

## 10. Cross-check against project evidence

Full table: `results/tables/lead_target_deep_dive/indirect_target_project_crosscheck.tsv`
(live-joined against the project's own frozen
`cross_dataset_genomewide` table -- not hand-typed, and re-derived
automatically if that frozen table is ever regenerated).

| Gene X | Parent | Hany FDR | GSE118713 FDR | GSE240112 tumor FDR | GSE245601 epi FDR | GSE111151 FDR |
|---|---|---|---|---|---|---|
| AXIN1 | USP34 | 0.688 | 0.312 | 0.190 | not queried | not queried |
| UBQLN1 | USP34 | 0.994 | 0.406 | 0.999 | not queried | not queried |
| TEAD1 | VEZF1 | 0.801 | **0.015** | **0.046** | 0.358 (nominal p=0.050) | 0.967 (no signal) |
| ARHGAP22 | VEZF1 | 0.791 | not testable | 0.244 | 0.757 | not queried |
| STUB1 | VEZF1 | 0.993 | 0.222 | 0.289 | not queried | not queried |

**Gap-completion**: TEAD1's full evidence profile across all 5 testable
resistance-RNA datasets is now shown (previously only 2 of 5 were
reported). TEAD1 reaches FDR<0.05 in 2 datasets, nominal p<0.05 without
surviving FDR correction in a 3rd (GSE245601 epithelial), and shows no
signal at all in GSE111151. No genome-wide TCGA or DepMap table exists in
this project (confirmed by directory search) -- both are scoped to the
four frozen candidates only, so TEAD1 cannot be looked up in either. TEAD1
remains the only indirect-target candidate reaching FDR<0.05 in any of
this project's own independent resistance datasets. The frozen shortlist
is unaltered by this table -- it is a lookup, not a re-ranking exercise,
and this RNA-level alignment does **not** by itself validate the
TEAD1-VEZF1 biochemical hypothesis (see Part 9).

---

## 11. Therapeutic window -- head-to-head

Full table: `results/tables/lead_target_deep_dive/USP34_VEZF1_head_to_head.tsv`
(16 separate dimensions, no master score). Figure:
`04_USP34_VEZF1_therapeutic_window_comparison.png` (four real-data
dimensions plotted; qualitative dimensions live in Figure 02, not
repeated as a fake composite).

Selected trade-offs (full 16-dimension detail in the table):
USP34 has the stronger **direct druggability** and the cleaner
**tamoxifen-specific** profile (0% baseline dependency), but carries a
**documented bone liability** and a **documented CNS/developmental human
finding**. VEZF1 has the stronger **cancer-dependency** and **dual-action**
signal, a **documented cardiovascular liability** (new this pass), and
essentially no direct druggability -- but the only candidate with a
druggable, project-evidence-aligned (though **unvalidated**) **indirect
-targeting hypothesis** (TEAD1; see Part 9 for why this is a hypothesis,
not a validated strategy).

---

## 12. Proposed experiments

Full table: `results/tables/lead_target_deep_dive/USP34_VEZF1_experimental_plan.tsv`.
Figure: `05_USP34_VEZF1_experimental_strategy.png`.

**EXP-1 (USP34)**: 4 arms (control / tamoxifen / USP34 inhibition
[genetic -- CRISPRi/shRNA/catalytic-dead C1903S, since no chemical probe
exists] / USP34 inhibition + tamoxifen) in an acquired tamoxifen
-resistant ER+ model, reading viability, apoptosis, clonogenic survival,
ER signaling, and a mechanistic AXIN1/RUNX2/Smad1 readout, with Bliss
synergy analysis. **Normal-cell comparator**: human primary MSCs under
osteogenic induction (direct extension of PMID 30181118), directly
testing whether the level of USP34 inhibition needed for a cancer-cell
effect also measurably impairs normal osteogenic differentiation -- the
central therapeutic-window question for USP34.

**EXP-2 (VEZF1, direct)**: the same 4-arm structure with VEZF1
inhibition/CRISPRi, explicitly separating (A) direct cancer-cell effect,
(B) additional tamoxifen sensitisation, and (C) synergy -- testing the
dual-action hypothesis directly. **Normal-cell comparator**: primary
human cardiomyocytes or a zebrafish larval assay directly extending PMID
31911272's own model, given the new cardiac-muscle finding.

**EXP-3 (TEAD1, VEZF1 indirect-target validation)**: TEAD1 inhibitor
(existing clinical-stage chemotype) +/- tamoxifen. **Primary readout**:
does TEAD1 inhibition reduce VEZF1 expression or VEZF1 target-gene
activity (VEGFR2, TIMP3/MMP2, CITED2 repression) in ER+ breast cancer
cells? This is the single most important unconfirmed link in the whole
phase -- the VEZF1-TEAD1 relationship was demonstrated in cardiac/
zebrafish tissue, not breast cancer cells, so this experiment tests
transferability and is written to be genuinely informative whether it
confirms or falsifies the hypothesis.

---

## 13. Final lead selection

**A. Best tamoxifen-specific combination target: USP34** (clean D_LOW
DepMap profile + FDR=0.042 sensitisation).
**B. Strongest direct cancer vulnerability: VEZF1** (27.3%
ER+/luminal DepMap dependency vs USP34's 0%).
**C. Strongest dual-action hypothesis: VEZF1** (only candidate combining
real baseline dependency + tamoxifen sensitisation).
**D. Easiest to drug directly: USP34** (real catalytic-domain structures;
VEZF1 has none).
**E. Strongest indirect-targeting opportunity: VEZF1**, via TEAD1 --
druggable now (as a pan-TEAD compound; no TEAD1-selective compound
remains in active clinical development), and independently significant in
2 of 5 of this project's own resistance-RNA datasets, though the
VEZF1-TEAD1 link itself rests on a single cardiac/zebrafish paper with no
perturbation evidence and is explicitly a hypothesis, not a validated
mechanism (Part 9). USP34's indirect-targeting landscape remains thinner
still (CSNK1A1/BPTF are unconfirmed screening leads at best), so VEZF1
wins this dimension by default, not because TEAD1 is strong in absolute
terms.
**F. Better normal-tissue profile: genuinely mixed, not clearly
better for either.** USP34 carries a documented postnatal bone liability
and a documented human CNS/developmental finding; VEZF1 carries a newly
-identified documented postnatal cardiovascular liability. Neither is
"cleaner" once the new cardiac finding is weighed against USP34's bone
finding -- this is stated plainly rather than forcing a false winner.
**G. Greater muscle concern: neither has a skeletal-muscle liability**
(USP34 is expression-only with zero function; VEZF1 has none identified
for skeletal muscle). If the question is muscle broadly including
cardiac, **VEZF1** has the documented finding (PMID 31911272); USP34 has
none.
**H. Greater bone/marrow concern: USP34** (documented postnatal-causal
bone-formation impairment; VEZF1 is inferred-only for bone).
**I. Greater cardiovascular/developmental concern: VEZF1** (documented
postnatal cardiac-contractile-function impairment, new this pass, plus
embryonic-lethal vascular knockout); USP34's comparable axis (CNS/
developmental) is real but in a different organ system, so this is not a
direct apples-to-apples "worse" comparison -- both have a real
developmental-axis finding, in different systems.
**J. Best current therapeutic-window hypothesis:** genuinely balanced,
not forced. USP34 offers a cleaner cancer-selectivity story but real,
documented normal-tissue causal evidence (bone) with no direct
inhibitor. VEZF1 offers a stronger cancer-biology hypothesis (dual
-action) and a genuinely exciting near-term-testable indirect route
(TEAD1) but a newly-identified real cardiovascular liability and no
direct druggability of its own.

**Gap-completion reassessment (the five specific questions asked this
pass):**
1. *Is USP34 still the lead?* **Yes, unchanged.** Nothing in this pass's
   new findings (the NFIC/tooth paper, the alias-confirmed muscle null
   result, the thinner-than-ever indirect-targeting landscape) weakens
   USP34's structural/druggability advantage or strengthens a competing
   case for VEZF1 as direct lead.
2. *Is VEZF1 still the backup?* **Yes, unchanged**, but specifically via
   the TEAD1 route rather than direct targeting, exactly as before.
3. *Does TEAD1 materially strengthen VEZF1's translational feasibility?*
   **Only modestly, and less than the prior version of this report
   implied.** It gives VEZF1 a concrete, near-term, testable next
   experiment (EXP-3) that USP34 lacks -- that much is still true and
   valuable. But it does not strengthen VEZF1's underlying biology or
   safety case, because the TEAD1 link itself is unvalidated.
4. *Is TEAD1 currently a validated indirect strategy or only a hypothesis
   worth testing?* **A hypothesis worth testing. Not validated.** (Part 9)
5. *Does either gene now have a clearly better therapeutic-window
   hypothesis?* **No.** The gap-completion pass added detail (the full
   characterization of VEZF1's cardiac finding, USP34's third
   developmental finding) but did not resolve the genuine trade-off
   described in J below -- if anything, the more precise characterization
   of both liabilities (VEZF1's partial-knockdown/postnatal cardiac
   phenotype being unusually drug-relevant; USP34's replicated
   mesenchymal-differentiation theme across three independent findings)
   makes both liabilities more concrete without making either clearly
   worse than the other.

**K. LEAD TARGET: USP34**, on translational-feasibility grounds --
it is the only one of the two with an actual structural handle for
future medicinal chemistry, its liabilities (bone, CNS/developmental)
are at least well-characterized enough to design a normal-cell comparator
experiment against directly (EXP-1), and its tamoxifen-specific,
low-baseline-dependency profile is the cleanest starting hypothesis of
the two.

**L. SECOND/BACKUP TARGET: VEZF1**, specifically via the TEAD1 indirect
route (EXP-3) rather than direct VEZF1 targeting, which remains poorly
feasible. VEZF1's dual-action biology is the more interesting cancer
-biology hypothesis of the two, and TEAD1 gives it a concrete, near-term
-testable path that USP34 currently lacks (USP34's only path is de novo
structure-based drug discovery, a multi-year program). **This is a real,
disclosed trade-off, not a forced ranking**: if EXP-3 confirms the
VEZF1-TEAD1 link in breast cancer cells, VEZF1/TEAD1 could become the
faster-moving program of the two despite USP34's stronger current
structural position.

---

## Outputs

**Report**: `results/reports/lead_target_deep_dive/USP34_VEZF1_translational_deep_dive.md`.

**Tables** (`results/tables/lead_target_deep_dive/`):
`USP34_VEZF1_full_tissue_expression.tsv`, `USP34_VEZF1_tissue_liability.tsv`,
`USP34_VEZF1_muscle_liability.tsv`, `USP34_VEZF1_bone_marrow_liability.tsv`,
`USP34_VEZF1_human_genetic_constraint.tsv`, `USP34_direct_targeting.tsv`,
`USP34_indirect_targets.tsv`, `VEZF1_direct_targeting.tsv`,
`VEZF1_indirect_targets.tsv`, `indirect_target_project_crosscheck.tsv`,
`USP34_VEZF1_head_to_head.tsv`, `USP34_VEZF1_experimental_plan.tsv`,
`verified_references.tsv`.

**Figures** (`results/figures/lead_target_deep_dive/`):
`01_USP34_VEZF1_tissue_expression_atlas.png`,
`02_USP34_VEZF1_liability_map.png`,
`03_USP34_VEZF1_direct_indirect_targetability.png`,
`04_USP34_VEZF1_therapeutic_window_comparison.png`,
`05_USP34_VEZF1_experimental_strategy.png`.

**Code**: `src/lead_target_deep_dive_data.py` (curated, source-cited data),
`src/lead_target_deep_dive_build_tables.py` (deterministic table builder,
live-joins the project's own frozen evidence for the cross-check table),
`src/lead_target_deep_dive_visualization.py` (figure builder).
