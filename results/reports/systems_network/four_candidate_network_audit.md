# Four-candidate comparative systems-network audit

USP34, VEZF1, EML5, CITED2 -- frozen therapeutic shortlist, unchanged. This
is a mechanistic/network interpretation layer over already-frozen
systems-network outputs. No CRISPR, bulk RNA, scRNA, pathway-enrichment, or
candidate-selection analysis was rerun; no frozen evidence/network file was
modified (confirmed by `git status` at the end of this report). Docking,
protein structure, toxicity, normal-tissue analysis, and final drug
selection remain explicitly out of scope.

Throughout, three evidence classes are kept visually and textually
separate:
- **OUR DATA** -- this project's own CRISPR screen and RNA-seq/scRNA-seq
  resistance/acute datasets.
- **NETWORK DATABASE EVIDENCE** -- STRING interaction scores, TRRUST
  regulatory edges, MSigDB gene-set membership. Curated by third parties;
  never causal evidence about tamoxifen resistance on its own.
- **INFERENCE/HYPOTHESIS** -- this analysis's own synthesis/judgment calls
  (classification tiers, "specific vs generic" calls, follow-up priority).
  Always flagged as such.

---

## Part 1 -- USP34 (reused from the prior USP34-specific audit)

Source: `results/tables/systems_network/USP34_shortest_paths.tsv`,
`USP34_bridge_gene_evidence.tsv`,
`results/reports/systems_network/USP34_bridge_gene_evidence.md`,
figures `07_USP34_shortest_paths.png`, `08_USP34_bridge_gene_evidence.png`
(all untouched by this audit).

**OUR DATA**: USP34 itself reaches CRISPR FDR=0.042 (`sensitising_KO`) and
RNA FDR=0.0073 in GSE118713 (up in TAMR vs MCF7; 1 of 3 resistance
datasets significant, GSE240112/GSE111151 non-significant but same
direction). GSE245601 acute: not significant (FDR=0.90).

**NETWORK DATABASE EVIDENCE**: USP34 reaches CTNNB1, PTEN, EP300 and SOX2
only via 2-hop STRING paths (never a direct 1-hop interaction), almost
entirely through ubiquitin/deubiquitination machinery -- RPS27A, UBB, UBC
(physical_PPI, STRING score 0.95-0.98) and USP9X (functional_association,
0.79). Re-audited here: **none of the four bridge genes reach FDR<0.05 in
our own CRISPR or resistance data.** USP9X and UBB reach only
`B_PARTIAL_SUPPORT` (one contradicted or non-robust nominal hit); RPS27A
and UBC are `C_NETWORK_ONLY_GENERIC_BRIDGE`.

**INFERENCE**: the 2-hop bridge network to CTNNB1/PTEN/EP300/SOX2 remains
an unvalidated network hypothesis, not a mechanism -- preserved unchanged
from the prior audit.

---

## Part 2 -- direct (1-hop) network neighborhood per candidate

Table: `results/tables/systems_network/four_candidate_direct_neighbors.tsv`
(30 rows). Source: `results/networks/systems_network/cytoscape/network_edges.tsv`
(STRING required_score>=0.7 / TRRUST / pathway co-membership -- the exact
thresholds already used throughout the systems-network phase; none lowered
here).

| Candidate | n direct neighbors | Breakdown |
|---|---|---|
| **USP34** | 10 | 4 physical_PPI (RPS27A 0.961, UBA52 0.968, UBB 0.954, UBC 0.982) + 6 functional_association (USP9X 0.793, MKLN1 0.790, TRIP12 0.723, ASH1L 0.822, SMG1 0.722, H7C0V5_HUMAN 0.740 -- an unreviewed UniProt fragment identifier, reported as-is) |
| **VEZF1** | 1 | DMTN only, via **pathway_co_membership** (not STRING) -- VEZF1 has zero STRING partners even down to STRING's lowest confidence band (score>=0.15), confirmed in the prior 30-phase audit and unchanged here |
| **EML5** | 0 | **NO RESOLVED NETWORK NEIGHBOURHOOD IN CURRENT ANALYSIS** -- zero STRING partners at any threshold, zero TRRUST edges, zero pathway co-membership edges |
| **CITED2** | 18 | 7 physical_PPI (EP300 0.999, CREBBP 0.995, TFAP2A 0.989, TFAP2C 0.945, LHX2 0.895, TFAP2B 0.880, HIF1A 0.915) + 6 functional_association (EPAS1 0.790, TP53 0.723, FOXO3 0.933, WT1 0.711, FOXO1 0.816, PCGF2 0.741) + 5 pathway_co_membership (COL1A1, ERBB2, PPARG, SMAD3, TGFBR2) |

Node classes of neighbors span candidate_partner, leading_edge,
resistance_gene, and "multiple" (e.g. TFAP2C is simultaneously a CITED2
STRING partner AND a genome-wide CRISPR FDR<0.05 hit -- see Part 5). No
neighbor is a CRISPR sensitiser for VEZF1 or EML5 (there are none to
classify).

---

## Part 3 -- candidate-pathway/module connections

Source: `results/tables/systems_network/candidate_pathway_membership.tsv`,
`resistance_pathway_consensus.tsv`, `multimodal_pathway_convergence.tsv`
(all frozen; read-only). Categories: **A** direct curated membership, **B**
direct network connection (candidate's own 1-hop STRING/TRRUST partner is a
curated pathway member), **C** indirect/2-hop connection (a candidate's
2-hop bridge target is itself a pathway member), **D** leading-edge/module
association (candidate itself drives the enrichment signal), **E** no
resolved connection.

### USP34
- **A + D**: `GOBP_CANONICAL_WNT_SIGNALING_PATHWAY`, `GOBP_WNT_SIGNALING_PATHWAY`,
  `GOBP_REGULATION_OF_CANONICAL_WNT_SIGNALING_PATHWAY`,
  `GOBP_REGULATION_OF_WNT_SIGNALING_PATHWAY` -- all STRONG_CONSENSUS, USP34
  itself is a curated member AND significant leading-edge gene in
  gse118713+gse111151. `REACTOME_TCF_DEPENDENT_SIGNALING_IN_RESPONSE_TO_WNT`
  -- DIRECTIONAL_CONSENSUS, USP34 leading-edge in gse118713.
- **C** (new finding from this audit -- USP34's 2-hop targets are
  themselves leading-edge genes for STRONG_CONSENSUS pathways USP34 has no
  direct route to): via CTNNB1 -> `HALLMARK_WNT_BETA_CATENIN_SIGNALING`,
  `HALLMARK_TGF_BETA_SIGNALING`, `REACTOME_SIGNALING_BY_RECEPTOR_TYROSINE_KINASES`;
  via PTEN -> `HALLMARK_UV_RESPONSE_DN`, `HALLMARK_APICAL_JUNCTION`; via
  EP300 -> `REACTOME_SIGNALING_BY_RECEPTOR_TYROSINE_KINASES`. SOX2 (the
  fourth 2-hop target) maps to no STRONG_CONSENSUS Hallmark/Reactome
  pathway leading edge.
- **B**: indirect-via-interactor rows exist but are dominated by large,
  generic ubiquitin/proteasome pathways (RPS27A/UBA52/UBB/UBC appear in
  hundreds of Reactome terms) -- reported in the full table but not
  singled out here as they add no new interpretable biology.

### VEZF1
- **A + D**: `GOBP_BLOOD_VESSEL_MORPHOGENESIS` (STRONG_CONSENSUS, sig in
  gse118713+gse240112, all_positive, median NES=1.63; VEZF1 leading-edge in
  gse240112) and `HALLMARK_HEME_METABOLISM` (STRONG_CONSENSUS, sig in 2 of
  3 resistance datasets, all_positive, median NES=1.47; VEZF1 leading-edge
  in gse240112). VEZF1 is a direct curated member of all 24 pathway rows in
  its table (no interactor-mediated rows exist -- VEZF1 has no STRING
  partners, Part 2).
- **C**: none found -- VEZF1's only network edge (DMTN) does not open a new
  pathway route beyond HALLMARK_HEME_METABOLISM (DMTN is itself leading-edge
  for that same pathway, reinforcing rather than extending it).
- **E**: no resolved connection to any WNT pathway -- not forced.

### CITED2
- **A + D, genuinely MULTIMODAL**: `HALLMARK_UV_RESPONSE_DN` -- CITED2 is a
  direct curated member AND leading-edge (gse240112) AND the pathway
  reaches CRISPR-level significance (FDR=0.0088, same positive direction as
  RNA). This is CITED2's only pathway that is simultaneously direct,
  candidate-driven, and multimodal (RNA+CRISPR) -- the single strongest
  pathway-level finding for any candidate in this audit.
- **A + D**: `GOBP_BLOOD_VESSEL_MORPHOGENESIS` -- CITED2 is ALSO an
  independent direct curated member and leading-edge gene (gse240112) of
  this pathway, the same one VEZF1 independently drives (Part 7
  convergence). `HALLMARK_GLYCOLYSIS` (DIRECTIONAL_CONSENSUS),
  `HALLMARK_HYPOXIA` (SINGLE_DATASET).
- **B**: `HALLMARK_E2F_TARGETS` (via TP53, STRONG_CONSENSUS, MULTIMODAL_PATHWAY
  in the frozen table but candidate-tagged only for direct connections, not
  interactor-mediated), `HALLMARK_ESTROGEN_RESPONSE_EARLY`/`_LATE` (via
  TFAP2C, STRONG_CONSENSUS), `HALLMARK_G2M_CHECKPOINT` (via HIF1A,
  STRONG_CONSENSUS), `HALLMARK_P53_PATHWAY` (via FOXO3/TP53,
  STRONG_CONSENSUS), `REACTOME_ESTROGEN_DEPENDENT_GENE_EXPRESSION` (via
  CREBBP/EP300, STRONG_CONSENSUS, MULTIMODAL_PATHWAY). **These are the
  pathways named in the task as CITED2's expected story -- confirmed
  present, but every one of them is indirect (via a 1-hop partner), not a
  direct CITED2 membership.** This distinction matters: CITED2 itself is
  not a curated member of HALLMARK_E2F_TARGETS or the estrogen-response
  sets.

### EML5
- **E, exclusively**: zero rows in `candidate_pathway_membership.tsv`. No
  pathway membership, no leading-edge status, no interactor-mediated route
  (EML5 has no STRING partners, Part 2). Nothing invented.

---

## Part 4 -- shortest paths to candidate-specific resistance nodes

Table: `results/tables/systems_network/four_candidate_shortest_paths.tsv`
(21 rows). Targets are candidate-specific, not forced (see script docstring
for the exact derivation rule per candidate). No path exceeds 2 hops.

- **USP34** -> CTNNB1, PTEN, EP300, SOX2: all **2-hop** (reused verbatim
  from `USP34_shortest_paths.tsv`; 3/2/2/2 equally-short paths
  respectively, all via ubiquitin-family/USP9X intermediates).
- **CITED2** -> EP300, CREBBP, TFAP2C, HIF1A, TP53: all **1-hop, direct**
  (4 physical_PPI + 1 functional_association). These are literally CITED2's
  Part-2 direct neighbors -- no bridge gene needed, no weak/remote flag.
- **VEZF1** -> DMTN: **1-hop, direct** (pathway_co_membership). This is the
  *only* node any path from VEZF1 can reach -- VEZF1's connected component
  in the frozen network is exactly `{VEZF1, DMTN}` (confirmed via
  `networkx.node_connected_component`). No other vascular/developmental
  gene is wired to VEZF1 in this network, even though several are pathway
  co-members (Part 3) -- co-membership only becomes a pairwise edge when a
  pathway's node-universe leading edge is <=10 genes
  (`MAX_PATHWAY_GENES_FOR_EDGES`, an existing, documented cap, not altered
  here), and VEZF1's larger pathways exceed that cap.
- **EML5**: **NO_RESOLVED_NETWORK_NEIGHBOURHOOD_IN_CURRENT_ANALYSIS** --
  zero edges, no path of any length exists, none manufactured.

---

## Part 5 -- bridge/connector gene evidence audit

Table: `results/tables/systems_network/four_candidate_bridge_evidence.tsv`
(11 rows: 4 reused USP34 + 5 new CITED2 + 1 new VEZF1 + 1 EML5 placeholder).
Same conservative rule as the original USP34 audit: **A** requires
FDR<0.05 in CRISPR and/or a resistance dataset; **B** requires only a
nominal (p<0.05) hit; **C** requires neither; **D** not assessable.
GSE245601 acute is reported per gene but never drives classification.

| Candidate | Bridge gene | Relationship | Classification | Key evidence |
|---|---|---|---|---|
| USP34 | USP9X | 2-hop intermediate | B_PARTIAL_SUPPORT | one contradicted nominal RNA hit; CRISPR null |
| USP34 | RPS27A | 2-hop intermediate | C_NETWORK_ONLY | flat everywhere |
| USP34 | UBC | 2-hop intermediate | C_NETWORK_ONLY | flat everywhere |
| USP34 | UBB | 2-hop intermediate | B_PARTIAL_SUPPORT | nominal CRISPR (p=0.015) + nominal acute only |
| CITED2 | EP300 | 1-hop direct partner | **A_DATA_SUPPORTED** | gse240112 FDR=0.042, up |
| CITED2 | CREBBP | 1-hop direct partner | C_NETWORK_ONLY | flat everywhere; no STRONG_CONSENSUS leading-edge pathway either |
| CITED2 | TFAP2C | 1-hop direct partner | **A_DATA_SUPPORTED** | **CRISPR FDR=0.048 (tolerance_associated_KO) + gse118713 FDR=0.001 (down)** -- dual-layer, the single best-supported connector in this audit |
| CITED2 | HIF1A | 1-hop direct partner | **A_DATA_SUPPORTED** | gse118713 FDR=0.037 + gse240112 FDR=0.029, both up (2 of 3 datasets) |
| CITED2 | TP53 | 1-hop direct partner | **A_DATA_SUPPORTED** | gse118713 FDR=0.0047 + gse240112 FDR=0.026, both down (2 of 3 datasets) |
| VEZF1 | DMTN | 1-hop direct partner (pathway_co_membership) | **A_DATA_SUPPORTED** | gse118713 FDR=0.0002 + gse240112 FDR=0.0002, both strongly up |
| EML5 | none | not applicable | D_NOT_ASSESSABLE | zero network neighbors |

**Notable asymmetry**: CITED2 and VEZF1's connectors are far more
independently data-supported than USP34's. Four of five CITED2 connectors
and VEZF1's only connector reach A-tier; none of USP34's four reach A-tier.
This is a genuine, unexpected finding, not something assumed going in.

---

## Part 6 -- generic hub check

Source: `results/networks/systems_network/network_node_metrics.tsv`
(degree/betweenness on the 119-node frozen network).

| Gene | Degree | High-betweenness bridge flag | Role |
|---|---|---|---|
| CTNNB1 | 38 (#1 in network) | True | USP34's 2-hop target |
| TP53 | 34 (#2) | True | CITED2's 1-hop partner (A-tier) |
| EP300 | 28 (#3) | True | CITED2's 1-hop partner (A-tier) |
| UBC | 20 | True | USP34's 1-hop bridge (C-tier) |
| UBB | 17 | True | USP34's 1-hop bridge (B-tier) |
| HIF1A | 17 | False | CITED2's 1-hop partner (A-tier) |
| RPS27A | 15 | True | USP34's 1-hop bridge (C-tier) |
| CREBBP | 16 | False | CITED2's 1-hop partner (C-tier) |
| USP9X | 6 | False | USP34's 1-hop bridge (B-tier) |
| TFAP2C | 3 | False | CITED2's 1-hop partner (A-tier, dual CRISPR+RNA) |
| DMTN | 1 | False | VEZF1's only connector (A-tier) |

**Per-candidate call** (INFERENCE, rationale cites the numbers above):

| Candidate | Call | Why |
|---|---|---|
| USP34 | **MIXED** | Direct WNT pathway membership is candidate-specific, not hub-driven. But its only network route to CTNNB1/PTEN/EP300/SOX2 runs through RPS27A/UBB/UBC -- all high-degree, high-betweenness generic hubs, none independently supported (Part 5). |
| VEZF1 | **SPECIFIC / INTERPRETABLE** | Its one connector (DMTN, degree=1) is the opposite of a generic hub, and is independently the most strongly resistance-associated bridge gene in this entire audit. |
| CITED2 | **MIXED** | TFAP2C (degree=3) and HIF1A (degree=17) are specific and independently supported. But EP300 (#3 highest-degree) and TP53 (#1 highest-degree) -- both A-tier -- are also the network's two biggest generic hubs; their prominence in CITED2's picture is partly structural, not necessarily CITED2-specific, even though both independently reach FDR<0.05. |
| EML5 | **NO RESOLVED NETWORK** | Nothing to assess. |

---

## Part 7 -- candidate convergence

Table: `results/tables/systems_network/four_candidate_convergence.tsv` (6
pairs). Checked: direct interaction, shared direct neighbor, shared bridge
gene, shared resistance pathway/leading-edge module, shared transcriptional
regulator, shared resistance hub.

| Pair | Direct interaction | Shared neighbors | Shared bridges | Shared pathways | Any convergence |
|---|---|---|---|---|---|
| USP34-VEZF1 | No | 0 | 0 | 0 | **No** |
| USP34-EML5 | No | 0 | 0 | 0 | **No** |
| USP34-CITED2 | No | 0 | 0 | 1 (weak, generic: `GOBP_POSITIVE_REGULATION_OF_SIGNALING`, not in resistance consensus) | Weak |
| VEZF1-EML5 | No | 0 | 0 | 0 | **No** |
| VEZF1-CITED2 | No | 0 | 0 | **11** (incl. `GOBP_BLOOD_VESSEL_MORPHOGENESIS`, STRONG_CONSENSUS, both independently direct-member+leading-edge) | **Yes -- real** |
| EML5-CITED2 | No | 0 | 0 | 0 | **No** |

**Shared transcriptional regulator**: none found -- no candidate has a
TRRUST regulatory (TF-target) edge in the frozen network at all (checked
directly against `network_edges.tsv`; reported honestly rather than
stretched from looser evidence). **Shared resistance hub**: none -- zero
genes in the entire 119-node network connect directly to >=2 candidates
(`n_candidates_directly_connected` max = 1), corroborating the earlier
Louvain-community finding that all four candidates land in different
network modules.

**The VEZF1-CITED2 convergence on `GOBP_BLOOD_VESSEL_MORPHOGENESIS` is
real and non-generic**: both candidates are independently curated members
AND independently significant leading-edge genes (each via their own
gse240112 expression signal) of the same STRONG_CONSENSUS pathway -- not
merely two genes co-occurring in a broad catch-all GO term. USP34-CITED2's
one shared term, by contrast, is a large, generic signaling category not
itself in the resistance consensus table and is not upgraded to
"convergence" here.

---

## Part 8 -- final head-to-head

Table: `results/tables/systems_network/four_candidate_network_audit.tsv`.

| | USP34 | VEZF1 | EML5 | CITED2 |
|---|---|---|---|---|
| Frozen CRISPR | sensitising_KO, FDR=0.042 | sensitising_KO, FDR=0.037 | nonsig sensitising, FDR=0.149 | nonsig sensitising, FDR=0.110 |
| Resistance RNA (datasets FDR<0.05 / 3) | 1 (gse118713, up) | 1 (gse240112, up) | **2** (gse118713+gse240112, both up) | 1 (gse240112, up) |
| Acute GSE245601 | not significant | not significant | not significant | not significant |
| Direct neighbors | 10 | 1 | 0 | 18 |
| Strongest direct interaction(s) | UBC/UBA52/RPS27A/UBB (physical_PPI, 0.95-0.98) | DMTN (pathway_co_membership only) | none | EP300/CREBBP/TFAP2A (physical_PPI, 0.99+) |
| Strongest pathway/module | GOBP_CANONICAL/REGULATION_OF_WNT_SIGNALING (direct+leading-edge, STRONG_CONSENSUS) | GOBP_BLOOD_VESSEL_MORPHOGENESIS + HALLMARK_HEME_METABOLISM (direct+leading-edge, STRONG_CONSENSUS) | none (E) | HALLMARK_UV_RESPONSE_DN (direct+leading-edge+**MULTIMODAL**) |
| Strongest bridge | USP9X (B) | DMTN (**A**) | none | TFAP2C (**A**, dual CRISPR+RNA) |
| Bridge independently supported? | partial | **yes** | N/A | **yes** |
| Generic-hub concern | HIGH (bridge layer only) | NONE | N/A | MODERATE (2 of 5) |
| Mechanistic specificity | MIXED | SPECIFIC/INTERPRETABLE | NO RESOLVED NETWORK | MIXED |
| Candidate convergence | weak (CITED2) | **real (CITED2)** | none | **real (VEZF1)** |
| Systems-mechanism classification | MODERATE SYSTEMS SUPPORT | MODERATE SYSTEMS SUPPORT | **DATA-SUPPORTED BUT MECHANISTICALLY UNRESOLVED** | MODERATE SYSTEMS SUPPORT |

No candidate reached STRONG SYSTEMS SUPPORT under this audit's (documented,
conservative) bar: that would require own-candidate multi-dataset RNA
significance AND own-candidate CRISPR significance AND a broad,
independently-validated, low-generic-hub-risk network story simultaneously
-- none of the four combine all of that at once. This does **not** change
the frozen therapeutic ranking.

---

## MECHANISTIC FOLLOW-UP PRIORITY (INFERENCE -- separate from the frozen shortlist)

1. **CITED2** -- broadest network reach (18 neighbors), one genuinely
   direct+multimodal pathway (HALLMARK_UV_RESPONSE_DN), 4 of 5 connectors
   independently data-supported including a dual-layer one (TFAP2C), and a
   real, non-generic convergence with VEZF1. The richest overall systems
   story, with one clearly flagged caveat (EP300/TP53 hub risk).
2. **VEZF1** -- smallest network footprint of the three connected
   candidates, but the *cleanest*: its only connector is independently very
   strongly supported and carries zero generic-hub risk, and its
   convergence with CITED2 is real. Small, but not confounded.
3. **USP34** -- real, direct, candidate-driven WNT pathway membership (a
   genuine finding, distinct from its network-bridge story) is worth
   keeping, but the widely-cited "USP34 talks to CTNNB1/PTEN/EP300/SOX2"
   *network* narrative specifically should not be advanced without
   independent validation of USP9X or another mechanism -- three of its
   four bridges are unsupported generic ubiquitin hubs.
4. **EML5** -- deserves attention, but of a different kind: not a
   network/mechanism follow-up (there is no network story to chase), but a
   basic literature/gene-function review given its unusually strong,
   reproducible RNA resistance signal (2 of 3 datasets, one at
   FDR=0.000129) despite zero systems-network footprint.

---

## FINAL REPORT

**1. Strongest mechanistic/network story for USP34?**
INFERENCE, built on OUR DATA + NETWORK DATABASE EVIDENCE: USP34 is a direct
curated member and significant leading-edge gene (candidate itself, not a
proxy) of canonical WNT-signaling GO:BP terms, reproducibly significant in
2 of 3 resistance datasets (STRONG_CONSENSUS), plus its own CRISPR
FDR=0.042. That direct pathway-membership finding is real. Its separate,
more visually striking 2-hop network story into CTNNB1/PTEN/EP300/SOX2 is
NETWORK DATABASE EVIDENCE only -- unsupported by our own data in 3 of 4
bridge genes and should be treated as a hypothesis, not a mechanism.

**2. Strongest story for VEZF1?**
OUR DATA + NETWORK DATABASE EVIDENCE, small but clean: VEZF1 is a direct,
leading-edge member of `GOBP_BLOOD_VESSEL_MORPHOGENESIS` and
`HALLMARK_HEME_METABOLISM` (both STRONG_CONSENSUS), backed by its own
CRISPR (FDR=0.037) and RNA (FDR=0.0195) significance, and its sole network
connector (DMTN) is itself one of the most strongly resistance-associated
genes in the whole audit (FDR<0.001 in 2 datasets) with zero generic-hub
confound.

**3. Strongest story for CITED2?**
OUR DATA + NETWORK DATABASE EVIDENCE: `HALLMARK_UV_RESPONSE_DN` is
CITED2's only pathway that is simultaneously direct membership, candidate
leading-edge, AND multimodal (RNA STRONG_CONSENSUS + CRISPR
FDR=0.0088) -- the single cleanest pathway-level finding of any candidate.
Layered on top: 4 of 5 direct network connectors independently reach
FDR<0.05 in our data, most notably TFAP2C (CRISPR FDR=0.048 AND RNA
FDR=0.001).

**4. Any meaningful mechanistic story for EML5, or still unresolved?**
Still unresolved, and should be reported as such. OUR DATA gives EML5 the
strongest own-candidate resistance-RNA signal of all four candidates (2 of
3 datasets significant, one at FDR=0.000129), but NETWORK DATABASE EVIDENCE
contributes nothing: zero STRING partners at any threshold, zero pathway
membership. This analysis layer has no mechanism to offer for EML5.

**5. Most specific (vs generic) network support?**
**VEZF1.** Its only network evidence is non-hub, non-generic, and
independently well-supported (Part 6).

**6. Strongest experimentally supported bridge, overall?**
**TFAP2C (CITED2's connector).** It is the only bridge/connector gene in
the whole audit with BOTH CRISPR (FDR=0.048) and RNA (FDR=0.001)
significance -- OUR DATA, two independent layers, same qualitative
direction of resistance-relevance. DMTN (VEZF1's connector) is close behind
with the single strongest RNA signal (FDR<0.001 in 2 datasets) but no
CRISPR signal.

**7. Relies most heavily on generic network hubs?**
**USP34.** Three of its four bridge/connector genes (RPS27A, UBB, UBC) are
high-degree, high-betweenness generic ubiquitin hubs with zero independent
data support (Part 5/6).

**8. Where do the four candidates genuinely converge?**
Only **VEZF1 <-> CITED2**, on `GOBP_BLOOD_VESSEL_MORPHOGENESIS` -- both
independently direct curated members and independently significant
leading-edge genes of the same STRONG_CONSENSUS pathway (11 shared GO:BP
terms overall). No other pair shares a direct interaction, a direct
neighbor, a bridge gene, a transcriptional regulator, or a resistance hub
(Part 7). USP34-CITED2 share only one broad, generic, non-consensus GO
term -- not treated as convergence.

**9. Which 1-2 candidates deserve deep literature/mechanism review next?**
**CITED2 and VEZF1** (INFERENCE, see Mechanistic Follow-Up Priority above)
-- the two candidates whose network stories are both richly supported by
OUR DATA and not primarily artifacts of generic network structure.

**10. Which candidate should NOT be overinterpreted despite looking
interesting in network visualizations?**
**USP34.** Its subnetwork is the most visually dense of the four (24 nodes,
186 edges in the candidate-subnetwork export) and superficially the most
"connected-looking," but that density is driven almost entirely by generic
ubiquitin/proteasome hub genes with no independent CRISPR or RNA support.
Secondary caution: within CITED2's otherwise well-supported picture,
specifically the EP300 and TP53 routes (not CITED2 as a whole) carry the
same generic-hub risk and should not be read as CITED2-specific biology on
their own.

---

## Outputs

Tables: `four_candidate_network_audit.tsv`, `four_candidate_direct_neighbors.tsv`,
`four_candidate_shortest_paths.tsv`, `four_candidate_bridge_evidence.tsv`,
`four_candidate_convergence.tsv` (all under `results/tables/systems_network/`).

Figures: `results/figures/systems_network/final_review/09_four_candidate_network_comparison.png`,
`.../10_four_candidate_mechanism_map.png`.

Reused, untouched: `USP34_shortest_paths.tsv`, `USP34_bridge_gene_evidence.tsv`,
`USP34_bridge_gene_evidence.md`, figures `07_*.png`, `08_*.png`.

**Stopping here per instructions.** Not starting druggability, structural
biology, or docking. No file under `results/tables/evidence_freeze/`,
`docs/THERAPEUTIC_SHORTLIST_FREEZE.md`, or the prior systems-network
outputs was modified.
