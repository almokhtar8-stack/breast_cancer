# Network/Mechanism Figure v2 -- POST-FREEZE EXPLORATORY Data Note

**This is a new, post-freeze exploratory mechanism analysis.** It does not
revise, recompute, or supersede any frozen CRISPR/RNA/pathway result or
candidate ranking. It replaces the v1 figure, which reused an older,
asymmetric STRING network built for a different four-gene shortlist
(USP34, VEZF1, EML5, CITED2) and therefore could say nothing real about
KDM1A or TLK2. v2 builds one new, consistent network for the CURRENT four
poster candidates -- KDM1A, TLK2, USP34, VEZF1 -- using identical rules
for all four.

## 1. Interaction source

STRING (string-db.org) REST API, `interaction_partners` endpoint,
`species=9606` (human), `required_score=700` (STRING's own "high
confidence" band, 0-1000 scale -- 0.7 on the 0-1 scale used in the output
tables). Fetched 2026-08-17 by
`scripts/download_string_network_v2_four_focus.py`, which writes local
TSVs; the analysis module (`src/poster_network_mechanism_v2.py`) reads
only those local files, never the network, at build/render time.

## 2. Exact network query rules (identical for all four candidates)

- **Level 1**: `interaction_partners` queried for all four candidate
  genes in one batched call, separately for `network_type=functional`
  (primary edge/score source) and `network_type=physical` (used only to
  label which functional edges are also direct physical interactions --
  never to add edges the functional query didn't return). This endpoint
  returns each gene's partners independently, so a hub-heavy candidate
  cannot crowd out a sparse one (unlike STRING's `network` endpoint with
  `add_nodes`, which was tried first and let KDM1A's very well-studied
  neighborhood consume the entire node budget, leaving USP34/VEZF1
  effectively empty -- rejected for that reason).
- **Level-1 display cap**: top 12 partners per candidate by STRING
  combined score (deterministic alphabetical tie-break), same K for all
  four. This only ever removes a candidate's own low-ranking partners; it
  never adds partners to a sparse candidate. Uncapped Level-1 partner
  counts were KDM1A=151, USP34=10, TLK2=5, VEZF1=0 -- USP34 and TLK2 show
  every partner they have (below the cap); only KDM1A's list is
  meaningfully truncated.
- **Level 2**: `interaction_partners` queried again (functional +
  physical) for the pooled, uncapped Level-1 partner set (166 genes). A
  resulting edge is displayed only if its target is (a) already a node in
  the graph -- a genuine bridge/convergence hit -- or (b) a member of a
  pre-specified canonical marker-gene list for the biological programs
  named in the task (estrogen/endocrine signaling, EMT, WNT/beta-catenin,
  cell cycle): `ESR1, ESR2, GREB1, FOXA1, PGR, CTNNB1, APC, AXIN1, GSK3B,
  TCF7L2, CDH1, VIM, SNAI1, SNAI2, ZEB1, TWIST1, RB1, E2F1, CCND1, CCNE1,
  CDK4, CDK6`. This list was fixed before the Level-2 query was filtered
  against it -- nothing was added after seeing which genes would result.
  Source nodes for Level-2 edges are restricted to the capped Level-1
  display set, so every bridge shown originates from a node that is
  itself displayed.
- One low-information node, `H7C0V5_HUMAN` (an uncharacterized STRING/
  TrEMBL entry, not an interpretable gene symbol), was excluded from
  display at both levels.

## 3. Network is undirected

All edges are undirected STRING association/interaction edges. **A
STRING edge does not imply activation, inhibition, or any specific
regulatory direction** -- it means the two proteins/genes have
evidence of interacting or being functionally associated, nothing more.
No direction, sign, or mechanism-of-action is drawn or implied anywhere
in this figure.

## 4. Final network size

**47 nodes, 147 edges, 3 connected components.**

| Candidate | Level-1 partners (capped/uncapped) | Component |
|---|---|---|
| KDM1A | 12 / 151 | large component (n=40), joint with USP34 |
| USP34 | 9 / 10 | same large component (n=40) |
| TLK2 | 5 / 5 | its own small component (n=6) |
| VEZF1 | 0 / 0 | isolated singleton |

## 5. Centrality / structure (networkx, unweighted)

- **Top hubs by degree**: HDAC1 (23), DNMT1 (17), HDAC2 (16), SNAI1 (16),
  EZH2 (14), RCOR1 (13), CTBP1 (13), KDM1A (12) -- all in KDM1A's
  CoREST/NuRD-complex neighborhood, reflecting real, densely
  interconnected chromatin-complex biology, not a construction artifact.
- **Top bridges by betweenness centrality**: DNMT1 (0.188), HDAC1
  (0.146), USP34 (0.11), CTNNB1 (0.108), SNAI1 (0.093), AR (0.064),
  RPS27A (0.064), CTBP1 (0.056).
- **KDM1A <-> USP34 shortest path** (the only candidate pair connected at
  all): `KDM1A -- DNMT1 -- UBC -- USP34` (3 edges; DNMT1 is a Level-1
  KDM1A partner, UBC is a Level-1 USP34 partner, and the DNMT1-UBC edge
  is a real Level-2 bridge). All other candidate pairs have no path in
  this network (KDM1A-TLK2, KDM1A-VEZF1, TLK2-USP34, TLK2-VEZF1,
  USP34-VEZF1: no connection).
- **Shared direct (Level-1) neighbors between any two candidates**: none.
  Any convergence between candidates in this network is indirect (via
  Level-2 bridges), not through a shared direct partner.

## 6. Does USP34 connect toward CTNNB1/WNT?

Yes, at Level 2: USP34's own Level-1 partners RPS27A, UBB, UBC, and
USP9X each have a real, frozen-at-query-time STRING edge to CTNNB1
(scores 0.735-0.880), and RPS27A/USP9X also connect to AXIN1 (WNT
destruction-complex component, scores 0.743-0.768). These are all
ubiquitin/proteasome-machinery genes acting as generic-hub bridges, not
a USP34-specific mechanism -- shown as such (thin Level-2 edges through
a shared intermediate, not a direct USP34-CTNNB1 edge, because none
exists at this threshold).

## 7. Does VEZF1 connect to anything beyond DMTN under this rule?

**No.** Under the same required_score>=0.7 STRING query used for the
other three candidates, VEZF1 returns **zero** interaction partners --
not even DMTN, which appeared in the older v1 figure via a
*pathway-co-membership* edge (MSigDB Hallmark, not STRING) rather than a
protein-interaction edge. This network is STRING-only, so that edge type
does not carry over. VEZF1 is a genuine, honestly-reported isolated
singleton in this analysis -- not an omission or a display choice.

## 8. Do the four candidates converge?

**Partially.** KDM1A and USP34 fall into the same connected component (a
real, if generic-hub-mediated, 3-edge bridge) and both Level-1
neighborhoods independently reach WNT/EMT/estrogen-signaling programs
(see halos in the figure and Section 6). TLK2 forms its own small,
separate, chromatin-assembly-only cluster (ASF1A, ASF1B, CABIN1, SRSF1,
TLK1) with no bridge to the other three at this threshold. VEZF1 is
fully isolated. This is reported as found -- convergence was not forced,
and the absence of a TLK2/VEZF1 bridge is itself a real, stated result.

## 9. Pathway/program annotation

Node "program membership" halos are read from the project's already-
frozen, locally-cached MSigDB GMT files (`data/reference/genesets/
hallmark.gmt`, `reactome.gmt`), unmodified -- no new enrichment test:
`HALLMARK_ESTROGEN_RESPONSE_EARLY/LATE` (Estrogen response),
`HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION` (EMT),
`HALLMARK_WNT_BETA_CATENIN_SIGNALING` (WNT/beta-catenin),
`HALLMARK_E2F_TARGETS` (E2F/cell cycle),
`REACTOME_CHROMATIN_MODIFYING_ENZYMES` (Chromatin regulation). A node's
halo is a direct gene-set-membership annotation, never a drawn
gene-to-pathway edge/box.

## 10. Druggability

The project's only druggability table (`results/tables/
druggability_safety/candidate_druggability.tsv`) covers the OLD
candidate shortlist (USP34, VEZF1, EML5, CITED2) at the protein level
only -- it has no entries for any Level-1/Level-2 partner or bridge node
appearing in this network. No druggable-node marker is shown, since
there is no reliable project-supported data to base one on for any node
other than the two candidates it already covers (USP34, VEZF1), which
are already distinctly sized/colored as candidates. This is stated
honestly rather than adding a marker with no real backing.

## 11. Visual encoding

- Node size: candidate > Level-1 partner > Level-2 bridge, with a
  modest additional boost for degree/betweenness (hubs/bridges read
  slightly larger within their tier).
- Node color: candidates use the frozen Okabe-Ito identity colors already
  used throughout the poster; Level-1/Level-2 nodes use fixed neutral
  blue-grey tones (never a candidate identity color).
- Edge: solid line for every real STRING edge (physical_PPI drawn
  thicker/more opaque than functional_association); a Level-2
  relationship is drawn through its real intermediate node, never as a
  dashed shortcut between non-adjacent nodes.
- Layout: `networkx.kamada_kawai_layout`, edge weight = `1 - score` (so
  higher-confidence interactions pull nodes closer) -- a genuine
  graph-layout algorithm, no node was manually positioned.
- Labels: placed via `adjustText`'s collision-avoidance solver against
  every other label and every node position; a thin leader line is drawn
  only where a label had to move off its default position.

## 12. Frozen-science integrity

No CRISPR effect/FDR, RNA differential-expression value, GSEA NES/FDR,
or candidate ranking was recomputed, re-thresholded, or altered. This
figure adds a new, clearly-labeled exploratory data source
(`data/reference/interactions/string_v2_level{1,2}_{functional,
physical}.tsv`) and reads it alongside already-frozen pathway gene sets;
it does not touch any tracked scientific result, and the frozen
candidate ranking is unchanged.
