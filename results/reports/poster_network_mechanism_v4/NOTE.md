# Network/Mechanism Figure v4 -- Five-Panel Presentation Rebuild of v2

**Post-freeze exploratory analysis. Presentation rebuild only -- no
science changed.**

## 1. Source

The existing v2 STRING graph, obtained by importing and calling
`poster_network_mechanism_v2.build_network()` unmodified: 47 nodes, 147
edges, 3 connected components, built from the already-downloaded local
STRING tables (`data/reference/interactions/string_v2_level{1,2}_
{functional,physical}.tsv`; species 9606, `interaction_partners`
endpoint, required_score >= 0.7, same rule for all four candidates).

## 2. No STRING requery

No network/API call of any kind occurs in `src/poster_network_
mechanism_v4.py` -- it imports no HTTP library and reads no URL. Only
the v2 module's already-tested loader chain is reused (verified by
test).

## 3. Local-subgraph derivation rule (identical for all four candidates)

For candidate X, the local panel shows the v2-graph subgraph induced on:

    {X}
  ∪ {direct (Level-1) neighbors of X in the v2 graph}
  ∪ {v2 nodes of kind "level2_bridge" adjacent to >= 1 of those
     neighbors}

with ALL v2 edges among that node set. No gene list is hand-typed; the
graph decides membership deterministically. Level-1 partners of OTHER
candidates (e.g. DNMT1, a KDM1A partner) do not enter another
candidate's panel unless the rule itself admits them, and no candidate
appears inside another candidate's panel.

## 4. Panel contents (as derived at build time)

| Panel | Candidate | Nodes | Edges |
|---|---|---|---|
| A | KDM1A | 29 | 105 |
| B | TLK2  | 6  | 13  |
| C | USP34 | 13 | 25  |
| D | VEZF1 | 1  | 0   |

The size differences are real and informative -- panels were not padded
or equalized. Panel D honestly shows VEZF1 alone: at the frozen
score>=0.7 STRING threshold it has zero partners. The old VEZF1--DMTN
relation was MSigDB pathway co-membership, not a STRING interaction,
and is therefore NOT drawn here.

## 5. Shortest-path calculation (Panel E)

Computed at render time with `networkx.shortest_path` /
`all_shortest_paths` on the full v2 graph -- nothing hand-typed. Result:

- **KDM1A <-> USP34**: shortest path `KDM1A -- DNMT1 -- UBC -- USP34`,
  3 edges. There are 4 equally short 3-edge routes, all via DNMT1
  (second hop UBC, UBB, UBA52, or RPS27A); the panel displays the
  networkx-returned path and states the count of equally short routes
  rather than hiding it.
- **All pairs involving TLK2 or VEZF1**: no path. TLK2 sits in its own
  separate connected component; VEZF1 is an isolated singleton at this
  threshold. Both facts are stated on the panel; no path was fabricated.

## 6. Undirected network -- no direction implied

Every edge is an undirected STRING functional-association /
physical-interaction edge. No arrows are drawn anywhere (all connectors
are plain lines), and nothing in the figure states or implies
activation, inhibition, or causation.

## 7. Rendering choices

- Each candidate panel gets its OWN Kamada-Kawai layout of its own
  subgraph (edge weight = 1 - STRING score), recentered on the
  candidate -- not a crop of a global layout.
- Node sizing: candidate largest; Level-1 partners medium; Level-2
  bridge nodes smaller; modest local-degree scaling. In the dense KDM1A
  panel, non-candidate markers are scaled down by a density factor so
  no node is hidden (no node is removed).
- Labels: candidate names bold white inside the candidate node; other
  labels seeded on the far side of their node from the candidate and
  then resolved with adjustText collision avoidance (thin leader lines
  where a label moved). In Panel A only, Level-2 nodes with local
  degree <= 1 are left unlabeled for legibility -- the nodes themselves
  remain drawn.
- Pathway halos from v2 are OMITTED from this figure for network
  readability (pathway analysis is the next poster section); the
  underlying `pathways` node annotations in the v2 data are untouched.
- Thin gray panel borders; one three-item legend (candidate / direct
  partner / Level-2 bridge).

## 8. Frozen-science integrity

No CRISPR/RNA/pathway value, candidate ranking, or network-construction
rule was recomputed or altered. v1/v2/v3 outputs are not overwritten.
The frozen therapeutic shortlist and science-freeze tag are unchanged
(verified by test).
