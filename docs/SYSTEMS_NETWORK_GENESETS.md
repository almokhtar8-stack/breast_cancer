# Systems-network phase: gene-set collections

Phase 3. Downloaded once via `scripts/download_genesets.sh`; all GSEA
modules in this phase read the local files below only.

| Collection | Local file | Resource | Version/date | Gene-set count | Source |
|---|---|---|---|---|---|
| Hallmark | `data/reference/genesets/hallmark.gmt` | MSigDB H | release 2024.1.Hs (2024-08-09) | 50 | `data.broadinstitute.org/gsea-msigdb/msigdb/release/2024.1.Hs/h.all.v2024.1.Hs.symbols.gmt` |
| Reactome | `data/reference/genesets/reactome.gmt` | MSigDB C2:CP:REACTOME | release 2024.1.Hs (2024-08-09) | 1,736 | `.../c2.cp.reactome.v2024.1.Hs.symbols.gmt` |
| GO Biological Process | `data/reference/genesets/go_bp.gmt` | MSigDB C5:GO:BP | release 2024.1.Hs (2024-08-09) | 7,608 | `.../c5.go.bp.v2024.1.Hs.symbols.gmt` |

KEGG was not included: the current MSigDB release does not redistribute
KEGG gene sets under an open license (KEGG requires a separate commercial
license for gene-set-level redistribution), and no other KEGG source was
available in this environment. This is an optional resource per the task
spec and its absence does not block the analysis — Hallmark + Reactome +
GO:BP already provide broad, non-overlapping pathway/process coverage.

Reactome is sourced through MSigDB's curated C2:CP:REACTOME collection
(identifiers and gene membership pulled directly from Reactome by MSigDB's
curation pipeline) rather than a separate direct Reactome GMT download, so
that all three collections share one consistent gene-symbol identifier
system and one release date — avoiding a symbol-harmonization mismatch
across collections.

## Gene-set size filtering

A standard GSEA pre-filter (`min_size=15`, `max_size=500` genes per set) is
applied uniformly to all three collections in every GSEA run (Phase 4/11).
This is a fixed methodological choice made before any dataset's results were
inspected — it excludes gene sets too small for a stable enrichment-score
estimate and too large to be biologically specific, and is not tuned per
dataset or per collection.
