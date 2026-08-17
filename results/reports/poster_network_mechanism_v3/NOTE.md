# Network/Mechanism Figure v3 -- VISUAL REFINEMENT ONLY of v2

**This is a rendering-only change. No science changed.** v3 calls
`poster_network_mechanism_v2.build_network()` and
`poster_network_mechanism_v2.network_stats()` directly, unmodified. It
does not rerun STRING, does not change the interaction source, confidence
threshold, node set, edge set, filtering rule, or pathway annotation
logic. Every fact in `results/reports/poster_network_mechanism_v2/
NOTE.md` still applies unchanged; this note only documents what changed
about the *picture*.

## What is identical to v2 (verified by test)

- Exact same 47 nodes, 147 edges, 3 connected components.
- Exact same node identities, kinds (candidate / Level-1 partner /
  Level-2 bridge), and `pathways` (program-membership) annotations.
- Exact same edges, scores, and `interaction_type`
  (physical_PPI / functional_association) labels.
- Exact same STRING source data
  (`data/reference/interactions/string_v2_level{1,2}_{functional,
  physical}.tsv`), unchanged, not re-downloaded.
- KDM1A and USP34 remain in the same (large) connected component; TLK2
  remains its own separate component; VEZF1 remains a degree-0 isolated
  singleton.

## What changed (rendering only)

1. **Layout.** v2 ran one global `kamada_kawai_layout` over the whole
   (disconnected) graph, which has no principled basis for placing
   separate connected components relative to each other and scattered
   TLK2 and VEZF1 far from the main cluster with large empty gaps. v3
   instead lays out each connected component SEPARATELY with its own
   real `kamada_kawai_layout` (topology within each component is
   unchanged, genuine graph layout), then composes the three components
   compactly: the KDM1A/USP34 component centered and dominant, the TLK2
   component scaled down and placed just above it, and the VEZF1
   singleton placed nearby on the right with a clear but modest gap.
   Component placement (translation/scale of each whole component) is a
   layout/composition choice; no node was individually hand-positioned
   within its component, and no edge was added, removed, or rerouted.
2. **Canvas.** Tighter figure aspect ratio matched to the actual layout's
   data extent, smaller outer margins, and a shorter title/subtitle block
   -- the network now occupies most of the canvas instead of floating in
   a mostly-empty page.
3. **Node scale.** Candidate nodes are substantially larger (base marker
   area increased ~2.5x over v2) and use bigger, bolder labels (fontsize
   18 vs. v2's 12.5). Hub/bridge sizing still comes from the same
   `degree`/`betweenness` values already computed by
   `build_network()` -- only the visual scale constants changed, not the
   underlying metric or which nodes count as hubs.
4. **Labels.** Same labeling rule as v2 (candidates + every Level-1
   partner + Level-2 nodes that are a canonical biology marker gene or
   have degree>=6), with `adjustText` collision avoidance again used so
   no labels overlap. High-degree/canonical-gene labels render larger and
   bold; peripheral labels stay smaller. No gene that was labeled in v2
   became unlabeled in v3, and vice versa.
5. **Edges.** Rendered more quietly than v2 (thinner lines, lower alpha)
   so nodes read as the dominant visual element; physical_PPI edges are
   still drawn a touch stronger than functional_association edges, using
   the same real `interaction_type` field from v2 -- no edge is
   highlighted as a "mechanistic path."
6. **Pathway halos.** Same wedge-per-program encoding as v2, same
   `PATHWAY_COLORS`, rendered with a slightly thinner ring and partial
   transparency so they read as secondary to candidate/hub color and
   size.
7. **Legend.** Consolidated v2's two separate legend blocks (network
   encoding; program/halo key) into one combined legend in a single
   corner.
8. **Title/subtitle.** Shortened to "Molecular networks reveal distinct
   candidate neighborhoods" / "High-confidence STRING associations
   (score >= 0.7); post-freeze exploratory analysis." Full methods detail
   (species, endpoint, exact query rules) remains documented only in the
   v2 NOTE.md, not restated in the figure.
9. **Micro-annotations.** Three short italic captions ("KDM1A-USP34:
   connected component", "TLK2: separate chromatin neighborhood",
   "VEZF1: no high-confidence STRING partners") placed in the whitespace
   near each component -- a direct restatement of the graph's own
   already-computed connected-component structure, not a new claim.

## Scientific interpretation (unchanged from v2)

Every edge remains an undirected STRING functional-association/physical-
interaction edge. **No edge implies activation, inhibition, or a causal
relationship.** This figure is exploratory and hypothesis-generating
only; frozen candidate rankings, CRISPR/RNA/pathway results, and the
underlying network construction are all unchanged from v2.
