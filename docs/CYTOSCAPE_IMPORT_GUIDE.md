# Cytoscape import guide

Files: `results/networks/systems_network/cytoscape/`
(`network_nodes.tsv` + `network_edges.tsv` = the full focused network;
`USP34_nodes.tsv`/`USP34_edges.tsv` etc = the four per-candidate
subnetworks, Phase 17).

## Import steps

1. **File -> Import -> Network from File...**, select the edges file
   (e.g. `network_edges.tsv`).
2. In the import dialog, set:
   - **Source Interaction** column = `source_gene`
   - **Target Interaction** column = `target_gene`
   - Mark `interaction_type` as an **Edge Attribute** (not source/target).
   - Mark `database_source`, `confidence`, `pathway`, `evidence_notes` as
     Edge Attributes too.
3. **File -> Import -> Table from File...**, select the matching nodes
   file (e.g. `network_nodes.tsv`).
   - Import target: **Node Table**.
   - Key column for network = `gene` (this matches the node names created
     from the edge import).
4. Repeat for any of the four per-candidate subnetwork pairs if a focused
   view is wanted instead of the full network.

## Suggested visual mappings

- **Node fill color** -> `node_class` (discrete mapping): `candidate` =
  dark red, `crispr_sensitiser` = red, `crispr_tolerance_associated` =
  blue, `leading_edge` / `resistance_gene` = gray-blue, `candidate_partner`
  = light gray, `multiple` = orange.
- **Node border width** -> `candidate_connection_count` (continuous,
  thicker border = connected to more of the four candidates).
- **Node size** -> `number_consensus_pathways` (continuous).
- **Edge line style** -> `interaction_type` (discrete): `physical_PPI` =
  solid, `regulatory` = dashed (with arrow, since TRRUST edges are
  directional TF->target), `pathway_co_membership` = dotted,
  `functional_association` = solid, thin.
- **Edge width** -> `confidence` (continuous; STRING edges only, `pathway_
  co_membership`/`regulatory` edges have no numeric confidence and should
  render at the default width).
- **Layout**: `yFiles Organic` or `Prefuse Force Directed` both work well
  on this network's size (~120 nodes, ~1,150 edges for the full network;
  each per-candidate subnetwork is far smaller and better read with
  `Circular` or `Attributes Layout` centered on the candidate node).

## Notes

- `confidence` is populated only for `STRING`-sourced edges (0-1 combined
  score); it is blank for `TRRUST` and `pathway_co_membership` edges by
  design (those sources do not have a comparable numeric confidence
  score) -- do not treat a blank confidence as zero.
- A gene named `H7C0V5_HUMAN` appears once, as a USP34 STRING partner --
  this is STRING's own identifier for an unreviewed/fragment UniProt
  entry, not a standard HGNC symbol; it is reported as-is rather than
  silently dropped or renamed.
- `EML5_edges.tsv` has a header row but zero data rows (EML5 has no
  curated network connections at all, in either functional or physical
  STRING networks, TRRUST, or pathway co-membership -- a real finding, not
  a build error; see `results/networks/systems_network/EML5_nodes.tsv`,
  which does have the single EML5 node). Import **File -> Import -> Table
  from File...** on `EML5_nodes.tsv` directly as a new Node Table (rather
  than starting from the empty edge file) to create the isolated node in
  Cytoscape.
