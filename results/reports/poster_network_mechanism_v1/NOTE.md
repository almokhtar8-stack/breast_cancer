# Network/Mechanism Figure v1 -- Short Data Note

**Purpose:** local gene/protein mechanism context for the four focus
candidates (KDM1A, TLK2, USP34, VEZF1) -- which nearby genes/proteins each
candidate connects to, how strong that connection is, and where a
plausible indirect targeting route exists. This is deliberately NOT the
pathway-biology figure (`poster_pathway_v1`): no EMT/estrogen-response/WNT
program-level scores are shown here, only local gene-level mechanism.

## 1. Sources used (all read unmodified, no new network build, no new GSEA)

- `results/tables/systems_network/four_candidate_direct_neighbors.tsv` and
  `four_candidate_shortest_paths.tsv` -- the project's existing, frozen
  STRING-derived local network. This network was built for an **earlier**
  four-candidate therapeutic shortlist -- USP34, VEZF1, EML5, CITED2
  (`docs/SYSTEMS_NETWORK_NODE_RULE.md`, Category A) -- not the current
  poster four (KDM1A, TLK2, USP34, VEZF1). USP34 and VEZF1 carry over
  directly; EML5 and CITED2 are not part of this figure at all.
- `results/tables/systems_network/gsea_crispr.tsv` -- the frozen GSEA run
  on the genome-wide CRISPR ranking, read only to check leading-edge
  *pathway membership* for KDM1A and TLK2.
- `src/post_audit_sensitivity_data.load_significant_sensitising_hits()` --
  the same already-tested loader used by the CRISPR discovery figure, for
  each candidate's rank/effect/FDR among the 13 significant sensitising
  hits.

## 2. Why KDM1A and TLK2 have no network panel

KDM1A and TLK2 were never part of the systems-network build: KDM1A was
held as a blind positive control at the time that network was built
(CLAUDE.md; blinding retired 2026-08-10), and TLK2 was never on that
earlier four-candidate shortlist at all. Checked directly against every
relevant frozen table for this figure:
`four_candidate_direct_neighbors.tsv`, `four_candidate_shortest_paths.tsv`,
`candidate_pathway_membership.tsv`, `recurrent_leading_edge_genes.tsv`,
`pathway_crispr_overlay.tsv`, and `data/reference/interactions/string_candidate_partners.tsv`
-- neither gene appears in any of them. No PPI/network edge is fabricated
for either gene; this is stated explicitly on their panels.

## 3. What each subpanel shows

- **USP34** (richest evidence): 3 frozen direct STRING partners
  (RPS27A, UBC, USP9X; solid edges, line width scaled to STRING
  confidence) plus 2 frozen rank-1 shortest-path bridge targets --
  CTNNB1 (via RPS27A) and SOX2 (via USP9X), dashed edges -- the
  project's own frozen route into WNT/stemness biology. RPS27A and UBC
  are themselves annotated in the frozen network tables as WNT-pathway
  members (`REACTOME_SIGNALING_BY_WNT`, `GOBP_WNT_SIGNALING_PATHWAY`,
  among others).
- **VEZF1** (real but minimal): its single frozen network edge, DMTN
  (pathway co-membership, `A_DATA_SUPPORTED_BRIDGE` tier in
  `four_candidate_bridge_evidence.tsv`), shared leading-edge pathways
  HALLMARK_HEME_METABOLISM / BLOOD_VESSEL_MORPHOGENESIS, and DMTN's own
  frozen human recurrence association (up in gse118713 and gse240112,
  both FDR<0.001). No further edges exist in the frozen tables for VEZF1
  -- this is not a display choice, it is the entire frozen network.
- **KDM1A** and **TLK2** (concept panels, no network): candidate node
  connected by dashed edges to rounded boxes, a visually distinct node
  type meaning "frozen GSEA leading-edge pathway membership," never a
  gene-gene interaction. KDM1A: leading-edge in
  REACTOME_CHROMATIN_MODIFYING_ENZYMES and
  REACTOME_ESTROGEN_DEPENDENT_GENE_EXPRESSION (CRISPR-ranked GSEA). TLK2:
  leading-edge in GOBP_REGULATION_OF_CHROMATIN_ORGANIZATION only -- the
  same pathway KDM1A is also leading-edge in, noted on the panel.

## 4. Pathway selection for KDM1A/TLK2 was curated, not a blind FDR cut

A first pass took the lowest-FDR leading-edge pathways per gene
automatically; for KDM1A this surfaced `REACTOME_SARS_COV_INFECTIONS`
(FDR=0.0027) ahead of the estrogen-signaling pathway (FDR=0.0114) --
real, frozen, but mechanistically uninformative for this project (a
generic transcriptional-machinery overlap, not a breast-cancer-relevant
program). Pathways shown were hand-selected from each gene's full frozen
leading-edge list for mechanistic interpretability, the same curation
approach already used for Panel C of `poster_pathway_v1`. Every NES/FDR
value itself is still read from the frozen table, never hand-typed;
nothing was added that isn't already a real leading-edge membership.

## 5. Encoding

- Circle = gene node. Rounded box = pathway-membership node (never a
  gene-gene edge).
- Solid edge = direct, project-supported STRING partner. Dashed edge =
  indirect (2-hop shortest-path bridge, or pathway-membership link).
- Line width on USP34's direct edges scales with STRING confidence score.
- Candidate identity colors are the frozen Okabe-Ito set already used
  throughout the poster (`poster_final_data.FOCUS_COLORS`, reused via
  `poster_exploration_v2_data.FOCUS_COLORS`); neighbor/partner and
  pathway-node colors are fixed neutral tones, deliberately not drawn
  from that identity palette so they never compete with candidate
  identity.
- Each panel's subtitle (rank, FDR) reuses the same 13-hit sensitising
  table and FDR<0.10 gate as the CRISPR discovery figure.

## 6. What was excluded

- No DepMap dependency, structural-tractability, or drug-association data
  -- out of scope for this figure (later, separate figures per the task
  ordering: network/mechanism, then pathway, then DepMap/structure).
- No candidate-candidate network edges are drawn (the frozen
  `four_candidate_convergence.tsv`/`candidate_candidate_connections.tsv`
  show no direct interaction or shared bridge between USP34 and VEZF1;
  this is not omitted evidence, there simply is none).
- EML5 and CITED2, present in the frozen network build, are not part of
  this figure -- they are not among the current four poster focus genes.

## 7. Confirmation no frozen result changed

No STRING edge, shortest path, GSEA NES/FDR, or CRISPR effect/FDR value
was recomputed, re-thresholded, or re-ranked. All values are read
directly from already-frozen tables at render time.
