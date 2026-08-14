# Focused external literature + mechanism review

USP34, VEZF1, EML5, CITED2 -- frozen therapeutic shortlist, unchanged. This
phase asks: *do published experimental studies support, contradict, or
fail to explain the biological mechanisms suggested by our own data?* No
upstream CRISPR/RNA-seq/scRNA-seq/pathway/network analysis was rerun; no
frozen output was modified; no TCGA/DepMap/druggability/docking/structural
work was done. All findings below come from four independent,
citation-verified literature searches (PubMed/PMC/NCBI eutils/Europe
PMC/UniProt/Human Protein Atlas), one per candidate, each paper confirmed
by directly fetching its record where accessible. Every PMID/DOI cited is
real and traceable to `results/tables/literature_mechanism/verified_references.tsv`;
none were invented. Full per-claim detail is in
`results/tables/literature_mechanism/four_candidate_claim_evidence.tsv`
(39 rows, claim_id references throughout this report).

Three evidence classes are kept strictly separate throughout, as required:
**OUR DATA** (this project's CRISPR/RNA-seq/network results), **LITERATURE**
(independently published experiments), **INFERENCE** (this report's own
synthesis/judgment, always flagged).

---

## PART 1-4 SUMMARY -- see per-claim table for full detail

**CITED2** (15 claims, claim_ids CITED2-01 to CITED2-15): the richest
literature base of the four. One LEVEL-1 direct human-tamoxifen paper
(van Agthoven et al. 2009, PMID 19904269) that is *internally
self-contradictory* (resistance-cell-model selection favors CITED2 up;
clinical outcome data show CITED2 up = better MFS and better tamoxifen
benefit). Independent breast-cancer mechanistic support exists for
CITED2-ER activity (PMID 23811274), CITED2-p53/chemoresistance (PMID
27627783), CITED2-metastasis via IKKa (PMID 27216153) and VEGFA/angiogenesis
(PMID 28008154), and CITED2-hypoxia/FOXO3 survival specifically "in breast
cancer cells" (PMID 18158893). Physical-interaction mechanisms with
EP300/CREBBP/TFAP2C/HIF1A are documented (PMIDs 12586840, 12778114,
40313859) but only in non-breast systems. Two clear CONTRADICTIONS exist
in other cancer types: CITED2 acts as an invasion-*suppressor* in colon
cancer (PMID 18054336) and a growth-*suppressor* in hepatocellular
carcinoma (PMID 23212831) -- opposite to its breast/prostate pattern.
Sub-question M (UV-response/general stress biology) returned **LEVEL 5 --
NO SUPPORT FOUND**: no paper links CITED2 to UV or genotoxic stress
response at all; this is a real literature gap relative to the project's
own HALLMARK_UV_RESPONSE_DN finding.

**VEZF1** (11 claims, VEZF1-01 to VEZF1-11): genuinely strong, well-verified
developmental/vascular-biology literature (founding paper PMID 9986727;
knockout phenotype PMID 15882861; angiogenic knockdown PMID 15031128;
target genes VEGFR2/endothelin-1/TIMP3/MMP2 PMIDs 11504723, ncomms3824).
**No relevant breast-cancer, ER+, or tamoxifen/endocrine-resistance
literature was identified in our search.** DMTN check: no independent
literature support for a functional VEZF1-DMTN relationship (physical,
regulatory, or functional) was identified in our search -- recorded as
an explicit negative result (claim VEZF1-10-DMTN-NULL). A
genuinely important literature discovery, independent of the project's own
network analysis: **VEZF1 directly represses CITED2** (ChIP-verified, PMID
29794136, claim VEZF1-06) -- a real biological link between two of the
project's four candidates that the project's own STRING-based network never
found (STRING has zero VEZF1 partners at any threshold). VEZF1's only two
cancer papers (both HCC) act through entirely non-vascular mechanisms
(PAQR4 transactivation PMID 36241701; TNS1/O-GlcNAcylation PMID 40858565).

**USP34** (8 claims, USP34-01 to USP34-08): USP34 experimentally regulates
canonical WNT/Axin/beta-catenin biology, but the phenotypic direction is
context-dependent. In HEK293T and colorectal cancer lines, USP34
deubiquitinates AXIN1, and USP34 knockdown *decreases* Wnt3A-induced
Wnt/beta-catenin transcriptional output (Lui et al. 2011, PMID 21383061,
verified 3 independent ways against the primary text, correcting a
mischaracterization found in a secondary review). This was demonstrated
only in HEK293T and colorectal cancer lines, **never in breast tissue**.
A separate mammary-epithelium paper (PMID 28499884) found USP34 loss
*increases* active beta-catenin/EMT/stemness -- the opposite direction
from the Axin/Wnt mechanism paper, a genuine, unresolved directional
contradiction that is not smoothed over here. USP34 has documented roles
in at least 5 other, unrelated
pathways (BMP2/osteogenesis, TNBC mitochondrial function, HCC/c-Myc,
laryngeal SCC/SOX2 and cisplatin resistance, GPCR-p38 inflammation) --
USP34 is a pleiotropic DUB, not a Wnt-dedicated enzyme. No paper tests
USP34 in ER+ breast cancer or tamoxifen. Pathway-level (non-USP34) context:
CXXC4 loss causally drives tamoxifen resistance via Wnt/beta-catenin
activation directly in MCF-7/BT474 tamoxifen-resistant sublines (PMID
31911277) -- independent confirmation that the pathway USP34 sits inside
is a real, causal driver of tamoxifen resistance in breast cancer, even
though USP34 itself has not been tested there.

**EML5** (5 claims, EML5-01 to EML5-05): after an exhaustive search across
PubMed/PMC/UniProt/NCBI Gene/Human Protein Atlas/STRING/BioGRID, only two
EML5-specific primary papers exist in the entire literature -- a 2004 rat
brain cloning paper (PMID 15225882) and a 2015 human epilepsy-tissue
expression paper (PMID 26730336) -- **neither in cancer, neither with any
functional perturbation, neither touching breast tissue or endocrine
signaling.** All "microtubule function" claims for EML5 trace to homology
inference (UniProt ECO:0000250), never direct experimental demonstration,
even in EML5's own dedicated papers. Human Protein Atlas shows **no**
breast-cancer prognostic signal for EML5 (only cancer-type hit: favorable
prognosis in glioblastoma -- opposite tissue, opposite direction from the
project's own finding). **Verdict: LITERATURE-MECHANISM UNRESOLVED.**

---

## PART 5 -- CONTRADICTION SEARCH (dedicated summary)

| Candidate | Contradiction | Claim ID | Severity |
|---|---|---|---|
| CITED2 | Clinical tamoxifen-benefit/MFS data show HIGH CITED2 = BETTER outcome, opposite to the resistance-cell-model-selection direction in the same paper | CITED2-01 | High -- directly on-topic |
| CITED2 | CITED2 acts as invasion-suppressor in colon cancer (knockdown increases invasiveness) | CITED2-13 | Moderate -- different tissue |
| CITED2 | CITED2 acts as growth-suppressor in hepatocellular carcinoma (PPARgamma effector) | CITED2-14 | Moderate -- different tissue |
| USP34 | Mammary-epithelium USP34 loss *increases* active beta-catenin/EMT/stemness, opposite direction to the Axin/Wnt mechanism paper | USP34-02 | High -- directly relevant tissue (mammary), opposite direction |
| VEZF1 | None found | -- | Absence of contradiction is not proof -- VEZF1 also has almost no cancer literature at all to contradict |
| EML5 | None found (no literature exists to contradict) | -- | Not applicable |

---

## PART 6 -- Claim-evidence table

`results/tables/literature_mechanism/four_candidate_claim_evidence.tsv` --
39 rows, every mechanistic claim in this report is traceable to a
`claim_id` in that table.

---

## PART 7 -- OUR DATA vs EXTERNAL LITERATURE vs INTEGRATED HYPOTHESIS

### USP34

**OUR DATA**: CRISPR sensitising_KO, FDR=0.042. RNA up in GSE118713
(TAMR vs MCF7), FDR=0.0073 (1 of 3 resistance datasets). Direct curated
membership + leading-edge status in 4 STRONG_CONSENSUS/DIRECTIONAL_CONSENSUS
GO:BP WNT-signaling gene sets. 2-hop STRING paths to CTNNB1/PTEN/EP300/SOX2
via ubiquitin-hub genes, none independently supported in our own data
(prior systems-network audit).

**EXTERNAL LITERATURE**: A real, specific, primary-verified USP34-AXIN1-Wnt
mechanism exists (PMID 21383061) -- but only in HEK293T/colorectal cells,
and contradicted in direction by a mammary-epithelium paper (PMID 28499884).
No breast-cancer or endocrine-resistance USP34 paper exists. USP34 is
independently documented as functionally pleiotropic (5+ unrelated
substrate pathways).

**INTEGRATED HYPOTHESIS** (inference only): our project's direct WNT
pathway-membership finding for USP34 is consistent with the literature
finding that USP34 experimentally regulates canonical WNT/Axin/beta-catenin
biology -- but that mechanism has never been tested in breast tissue, the
phenotypic direction is context-dependent (directionally contradicted in
the one mammary paper that does exist), and there is independent (non-USP34)
literature proof that Wnt/beta-catenin activation causally drives tamoxifen
resistance in breast cancer cells (PMID 31911277, CXXC4). The most
defensible combined statement: *USP34 sits inside a pathway with
independently demonstrated causal relevance to tamoxifen resistance, and
USP34 itself has a real, literature-documented (non-breast) mechanism for
engaging that pathway, though its phenotypic direction is context-dependent
-- but no study has yet connected USP34 specifically to that causal
relevance.*

### VEZF1

**OUR DATA**: CRISPR sensitising_KO, FDR=0.037. RNA up in GSE240112, FDR=0.0195
(1 of 3 resistance datasets). Direct curated membership + leading-edge
status in 2 STRONG_CONSENSUS pathways (GOBP_BLOOD_VESSEL_MORPHOGENESIS,
HALLMARK_HEME_METABOLISM). Only network connector: DMTN (pathway
co-membership, not STRING), itself very strongly independently
RNA-resistance-significant (prior audit).

**EXTERNAL LITERATURE**: VEZF1 is a real, well-characterized vascular
developmental transcription factor with documented direct target genes
(VEGFR2, endothelin-1, TIMP3/MMP2, Dnmt3b, and -- newly discovered by this
search -- CITED2, PMID 29794136). **No** breast-cancer, ER+, or endocrine
literature exists. **No** documented VEZF1-DMTN relationship exists at
all (explicit negative finding). VEZF1's only cancer literature (HCC, 2
papers) works through non-vascular mechanisms entirely unrelated to its
namesake angiogenic function.

**INTEGRATED HYPOTHESIS** (inference only): the project's "VEZF1 -> vascular
program" pathway finding rests on real, non-generic biology (VEZF1 truly
is a specific, well-validated angiogenesis transcription factor, not
pathway-enrichment noise) -- but literature does not independently connect
that biology to tumor angiogenesis, hypoxia adaptation, or any drug
resistance mechanism in any system. No independent literature support for
a functional VEZF1-DMTN relationship was identified; the connection is
currently supported only by the project network/pathway construction.
The one substantive new fact this search adds --
VEZF1 directly represses CITED2 -- is worth carrying forward as a literature
fact independent of, and not previously captured by, the project's own
network (which found no VEZF1-CITED2 edge at all).

### CITED2

**OUR DATA**: CRISPR non-significant (FDR=0.110, nominal sensitising
direction). RNA up in GSE240112, FDR=0.0087 (1 of 3 resistance datasets).
Broadest direct network neighborhood of the four candidates (18 STRING
neighbors). Single strongest own-data pathway finding: direct curated
membership + leading-edge + genuinely multimodal (RNA+CRISPR) status in
HALLMARK_UV_RESPONSE_DN. Indirect (via 1-hop partners) connections to
estrogen-response/E2F/G2M/P53 Hallmark sets.

**EXTERNAL LITERATURE**: The richest and most directly relevant literature
base of the four candidates, including one LEVEL-1 human-tamoxifen paper
that is internally contradictory, several independent breast-cancer
mechanistic papers broadly supporting a pro-tumor/pro-resistance role
(chemoresistance via p53 suppression, reduced 4-OHT sensitivity via
enhanced ER activity, metastasis via IKKa/VEGFA), physical-interaction
mechanisms with EP300/CREBBP/TFAP2C/HIF1A (non-breast systems), and two
clear tissue-context-dependent contradictions (CITED2 as tumor suppressor
in colon and liver cancer). No UV-response/stress-biology literature link
was found at all.

**INTEGRATED HYPOTHESIS** (inference only): CITED2 is genuinely the
best-supported candidate in external literature for a general pro-tumor,
pro-chemoresistance role in breast cancer -- but the literature's ONE
direct test of the tamoxifen-resistance question specifically produced
contradictory clinical vs. mechanistic results, meaning "CITED2 drives
tamoxifen resistance" is not a settled literature fact, it is a genuinely
open, actively contested question that this project's own multimodal
HALLMARK_UV_RESPONSE_DN finding does not resolve either way (that pathway
has no independent literature support for CITED2 at all). The
partner-mechanism literature (EP300/CREBBP/TFAP2C/HIF1A/FOXO3/TP53) is
real but entirely non-breast, so the proposed chain CITED2->coactivator->
estrogen/stress programs->resistance is *plausible by assembly of separate
papers*, not demonstrated end-to-end by any single study.

### EML5

**OUR DATA**: Strongest own-candidate resistance-RNA signal of all four
candidates (2 of 3 datasets significant, one at FDR=0.000129, both up).
CRISPR non-significant. Zero network/pathway footprint of any kind.

**EXTERNAL LITERATURE**: Essentially nothing. Two EML5-specific primary
papers exist in all of PubMed/PMC (rat brain characterization, human
epilepsy tissue), neither touching cancer, breast tissue, or endocrine
biology. No breast-cancer prognostic signal in Human Protein Atlas.

**INTEGRATED HYPOTHESIS**: none can be responsibly proposed. This
candidate remains **DATA-SUPPORTED BUT MECHANISTICALLY UNRESOLVED** by
both our own systems-network analysis AND independent external literature
-- the two analyses independently agree on this, which strengthens
confidence that EML5's resistance-RNA signal is real but currently
unexplained, not that the unresolved status reflects a gap in this
project's own methods.

---

## PART 8 -- Mechanistic model per candidate

**USP34: C. PLAUSIBLE CROSS-SYSTEM MECHANISM**
> USP34 -> AXIN1 deubiquitination/stabilization (PMID 21383061, colorectal/HEK293T)
> -> Wnt/beta-catenin transcriptional output (phenotypic direction
>    context-dependent, not a fixed activate/inhibit relationship)
> -> possible endocrine phenotype (never tested in breast; directionally
>    contradicted by a mammary-epithelium paper, PMID 28499884)

**VEZF1: D. NETWORK/LITERATURE HYPOTHESIS ONLY**
> VEZF1 -> vascular/developmental angiogenic program (real, well-documented,
>   non-cancer biology)
> -> possible tumor-state adaptation (no literature connects this to
>    resistance in any system; no independent literature support for a
>    functional VEZF1-DMTN relationship was identified -- the DMTN
>    connector is currently supported only by the project network/pathway
>    construction)

**CITED2: B. STRONG BREAST-CANCER MECHANISTIC SUPPORT, INDIRECT/CONTESTED FOR ENDOCRINE RESISTANCE**
> CITED2 -> EP300/CREBBP and/or TFAP2C, FOXO3, TP53 (documented partners,
>   mostly non-breast for the physical interactions, breast-specific for
>   the FOXO3-hypoxia-survival and p53-chemoresistance links)
> -> ER transcriptional activity / hypoxia survival / p53 suppression
>   (breast-cancer-demonstrated)
> -> endocrine-resistance phenotype (CONTESTED -- the one direct test found
>    contradictory clinical vs. mechanistic evidence; not a settled finding)

**EML5: E. MECHANISTICALLY UNRESOLVED**
> EML5 -> unresolved

---

## PART 9 -- Head-to-head literature assessment

Full table: `results/tables/literature_mechanism/four_candidate_literature_comparison.tsv`.

| | USP34 | VEZF1 | CITED2 | EML5 |
|---|---|---|---|---|
| Direct tamoxifen evidence | 0 papers | 0 papers | **2 papers** (contested direction) | 0 papers |
| Endocrine-resistance evidence | 0 | 0 | 1 | 0 |
| ER+ evidence | 0 | 0 | **3** | 0 |
| Breast-cancer mechanistic evidence | 1 | 0 | **6** | 0 |
| Relevant pathway/mechanism evidence (other systems) | 4 | 9 | 6 | 0 |
| Candidate-partner evidence | USP34-SOX2, USP34-c-Myc (both non-breast) | **VEZF1 directly represses CITED2** (literature fact, not in our network) | Every named STRING partner (EP300/CREBBP/TFAP2C/HIF1A/FOXO3/TP53) has a documented mechanistic relationship | none |
| Evidence from other cancers only | 5 | 9 | 9 | 0 |
| Contradictory evidence | 1 (mammary EMT paper, opposite direction) | 0 found | **2** (colon, liver -- opposite tissue-context direction) | 0 |
| Overall literature depth | MODERATE | MODERATE-STRONG (vascular biology only) | **STRONGEST** | VERY THIN |
| Fit with our own data | PARTIAL | WEAK-PARTIAL | PARTIAL, genuinely mixed | UNRESOLVED |
| Mechanism confidence | LOW-MODERATE | LOW (resistance-specific); MODERATE (general vascular biology) | MODERATE (with contested-direction caveat) | NONE |

**LITERATURE-BASED mechanistic follow-up order** (INFERENCE; does not
alter the frozen therapeutic ranking, which remains USP34 > VEZF1 > EML5 >
CITED2 per the earlier freeze):

1. **CITED2** -- richest literature base, real breast-cancer mechanisms,
   but the direction question (does high CITED2 predict better or worse
   tamoxifen response?) must be resolved before further investment.
2. **USP34** -- the single most mechanistically specific, primary-verified
   biochemical pathway (Axin/Wnt) of the four, undermined only by never
   having been tested in breast tissue and one directional contradiction.
3. **VEZF1** -- real, specific, non-generic biology, but no relevant papers
   were identified in our search bridging that biology to any resistance
   phenotype; the newly found VEZF1-CITED2 literature link is worth
   carrying into the CITED2 track.
4. **EML5** -- not ready for mechanism-driven literature follow-up; would
   need basic functional characterization first, independent of resistance
   biology.

---

## PART 10 -- Should we continue to druggability?

| Candidate | Answer | Rationale (INFERENCE) |
|---|---|---|
| USP34 | **NOT YET** | Real biochemical mechanism exists but untested in breast/endocrine context and directionally contradicted in the one mammary paper found; premature to move to targetability review. |
| VEZF1 | **NOT YET** | Real, specific biology, but no relevant papers were identified in our search connecting it to any resistance mechanism in any system; the pathway-enrichment signal currently has no external validation. |
| CITED2 | **YES, BUT WITH CAUTION** | Richest, most direct evidence base of the four, including breast-cancer and ER+-specific mechanistic support -- but the single most relevant study is self-contradictory and CITED2 is documented as tissue-context-dependent (oncogenic in breast/prostate, tumor-suppressive in colon/liver), so any druggability rationale must explicitly account for the contested direction rather than assume it. |
| EML5 | **NO CURRENT MECHANISTIC BASIS** | Both our own network analysis and an exhaustive independent literature search agree: there is currently nothing to target rationally. |

This is not the druggability analysis itself -- no compounds/inhibitors
were searched.

---

## FINAL REPORT

**1. Is there direct published evidence connecting CITED2 to tamoxifen or
endocrine resistance?**
LITERATURE: Yes, but it is internally contradictory. van Agthoven et al.
2009 (PMID 19904269, Br J Cancer) directly tested this in a tamoxifen-
resistant ZR-75-1 cell model and in two independent clinical cohorts:
CITED2 upregulation was selected for during resistance acquisition in the
cell model, but in patients, HIGH CITED2 predicted LONGER metastasis-free
survival (HR=0.71, p=0.017) and GREATER clinical benefit from tamoxifen
(OR=1.91-2.20, p<=0.017) -- the opposite of a simple "CITED2 drives
resistance" story. Separately, Lau et al. 2013 (PMID 23811274) showed
CITED2 overexpression DOES reduce 4-OHT growth-inhibitory sensitivity in
MCF-7/T47D/CAMA-1 cells in vitro. INFERENCE: the literature answer is
"yes, tested directly, but the direction of effect is unresolved," not a
clean "yes, CITED2 drives resistance."

**2. What is the strongest experimentally demonstrated CITED2 mechanism
relevant to our findings?**
LITERATURE: CITED2 overexpression increases ligand-independent ER
transcriptional activity (2-4 fold) and reduces 4-OHT sensitivity in three
ER+ breast cancer lines (Lau et al. 2013, PMID 23811274) -- the most
directly relevant, breast-specific, ER+-specific mechanistic result found.

**3. Is VEZF1 genuinely implicated in breast/endocrine resistance, or is
our result novel/indirect?**
LITERATURE: No relevant breast, ER+, or endocrine-resistance VEZF1 papers
were identified in our search. OUR DATA independently and reproducibly
implicates VEZF1 (CRISPR + RNA + 2 STRONG_CONSENSUS pathways). INFERENCE:
this project's VEZF1 finding is genuinely novel relative to the literature,
not a confirmation of prior work, and the literature offers no external
validation of a resistance mechanism, only validation that VEZF1's
angiogenic biology itself is real.

**4. Does DMTN strengthen the VEZF1 mechanism?**
LITERATURE: No. An exhaustive search found no documented relationship
(physical, regulatory, or functional) between VEZF1 and DMTN.
INFERENCE: no independent literature support for a functional VEZF1-DMTN
relationship was identified; the connection is currently supported only
by the project network/pathway construction, not by independently
validated biology.

**5. Is USP34-WNT experimentally established independently of our
network?**
LITERATURE: Yes, specifically and rigorously: Lui et al. 2011 (PMID
21383061, Mol Cell Biol) showed USP34 deubiquitinates AXIN1, opposing
tankyrase-dependent degradation, and that USP34 knockdown DECREASES
Wnt3A-induced TOPFlash reporter activity and Wnt target genes (NKD1,
TNFRSF19) in HEK293T and colorectal cancer cell lines. USP34 experimentally
regulates canonical WNT/Axin/beta-catenin biology, but (per Q6 below) the
phenotypic direction is context-dependent, not a fixed activate/inhibit
relationship.

**6. Has USP34-WNT been connected specifically to breast cancer or
endocrine resistance?**
LITERATURE: No. No relevant papers were identified testing the
USP34-AXIN1-Wnt mechanism in any breast cell line or endocrine-therapy
context. A separate, non-Wnt-pathway
USP34 study in normal mammary epithelium (Oh et al. 2017, PMID 28499884)
found USP34 loss INCREASES active beta-catenin and promotes EMT/stemness
-- opposite in direction to the Axin/Wnt mechanism, an unresolved
contradiction. Pathway-level (non-USP34) evidence independently shows Wnt/
beta-catenin activation causally drives tamoxifen resistance in ER+ breast
cancer cells (CXXC4 loss, PMID 31911277).

**7. Is there a credible EML5 mechanism, or is EML5 still unresolved?**
LITERATURE + INFERENCE: Still unresolved, and the literature search
independently corroborates our own network analysis's conclusion, rather
than resolving it. Only two EML5-specific primary papers exist in all
searchable literature (PMID 15225882, PMID 26730336), neither in cancer.
Verdict stands: LITERATURE-MECHANISM UNRESOLVED.

**8. Which candidate has the strongest external mechanistic validation?**
LITERATURE: CITED2 -- by volume and directness (15 verified claims,
including breast-cancer- and ER+-specific mechanistic papers), though with
an important contested-direction caveat on the single most relevant study
(PMID 19904269).

**9. Which candidate is potentially the most novel?**
INFERENCE: VEZF1. Its network/pathway signal in this project is real and
statistically robust (OUR DATA), but literature offers zero precedent for
a resistance role of any kind -- if VEZF1 turns out to be real, it would
be a genuinely new finding, not a confirmation of known biology.

**10. Which candidates are ready for a druggability review?**
INFERENCE: CITED2 only, and with an explicit caution flag about the
contested direction of effect (see Part 10 table). USP34 and VEZF1: not
yet. EML5: no current mechanistic basis at all.

**11. Which mechanistic claims from our network analysis became stronger
after literature review?**
INFERENCE:
- USP34's direct WNT-pathway-membership finding is now backed by a real,
  specific, primary-verified biochemical mechanism (Axin deubiquitination,
  PMID 21383061) -- stronger, though still untested in breast tissue.
- CITED2's connections to its network partners (EP300, CREBBP, TFAP2C,
  HIF1A, FOXO3, TP53) are now all independently literature-documented as
  real physical/regulatory relationships (claims CITED2-06, 07, 08, 12) --
  stronger, though mostly non-breast.
- VEZF1's vascular/developmental pathway membership is now confirmed as
  real, specific, non-generic biology by extensive independent literature
  (11 verified papers) -- stronger as a statement about VEZF1's general
  function, though not as a resistance mechanism specifically.
- A previously unknown fact emerged that neither the project's own network
  analysis nor the original hypothesis anticipated: **VEZF1 directly
  represses CITED2** (PMID 29794136) -- a literature-only finding that adds
  a new, real biological link between two of the four frozen candidates.

**12. Which claims became weaker or should be removed?**
INFERENCE:
- USP34's 2-hop network-bridge story to CTNNB1/PTEN/EP300/SOX2 remains
  weak (unchanged by literature, since no paper tests it), and is now
  additionally complicated by a directional contradiction in the one
  breast/mammary USP34 paper that does exist (PMID 28499884) -- this
  should be weighted down, not up.
- VEZF1's "vascular program -> resistance" framing is weaker after
  literature review, not stronger: the underlying biology is real, but
  literature supplies no bridge at all from angiogenesis/vascular biology
  to drug resistance, hypoxia adaptation, or any tumor phenotype in any
  system -- the resistance-relevance step of the chain is now explicitly
  unsupported rather than merely untested.
- CITED2's implicit "CITED2 up = worse tamoxifen response" framing (which
  would be a natural but unstated assumption from the project's own
  resistance-associated-up-in-GSE240112 finding) should be weakened/flagged:
  the one paper that directly tested clinical tamoxifen benefit vs. CITED2
  found the OPPOSITE association. This does not remove CITED2 as a
  candidate of interest, but it does mean "CITED2 drives resistance" should
  not be stated as if it were literature-confirmed.
- EML5's own-data resistance-RNA signal is unchanged in strength by
  literature (there is no literature to weaken or strengthen it), but its
  candidacy for prioritized systems/mechanism follow-up is correctly
  ranked last, now with independent literature confirmation (not just our
  own null network result) that there is nothing external to build on yet.

---

**Stopping here per instructions.** Waiting for review before any commit.
No frozen output (systems-network, evidence-freeze, or otherwise) was
modified -- confirmed by `git status` below.
