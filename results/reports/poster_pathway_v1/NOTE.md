# Pathway Biology Figure v1 -- Short Data Note

**Purpose:** answers "what biological programs distinguish endocrine
resistance/recurrence from acute tamoxifen response?" -- the pathway
entry point of the poster story, separate from the CRISPR discovery
figure and the transcriptomic hero heatmap.

**1. Exact frozen source tables** (read unmodified, no new enrichment
run, no re-ranking): `results/tables/systems_network/gsea_{dataset}.tsv`
for `dataset` in `gse118713, gse111151, gse240112, gse245601` -- the
project's already-frozen GSEA output (one row per Hallmark/Reactome/GO_BP
pathway per dataset). Read via
`poster_exploration_v2_data.load_pathway_trajectories()`, the same
already-tested loader used elsewhere in this project's figure bank.

**2. Exact pathways displayed:**
- Panel A: `HALLMARK_ESTROGEN_RESPONSE_EARLY`, `HALLMARK_ESTROGEN_RESPONSE_LATE`
- Panel B: `HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION`
- Panel C (supporting): `HALLMARK_TGF_BETA_SIGNALING`, `HALLMARK_APICAL_JUNCTION`,
  `HALLMARK_E2F_TARGETS` -- chosen after inspecting all Hallmark pathways
  present in every one of the four frozen tables for ones that (a) use
  directly comparable NES values and (b) reinforce, rather than dilute,
  Panels A/B's story: TGF-β signaling and apical junction (an epithelial
  integrity/junction program, the closest available proxy to "ECM" in the
  Hallmark collection) show the same resistance/recurrence-up,
  acute-down reversal as EMT; E2F targets shows the same suppressed-
  everywhere pattern as estrogen response. No ECM-named Hallmark gene set
  exists in MSigDB, so apical junction was used instead of an invented
  ECM label.

**3. What NES means:** the Normalized Enrichment Score from the frozen
GSEA run -- positive means the pathway's genes are shifted toward the
up-regulated end of that context's ranked gene list, negative toward the
down-regulated end, normalized for gene-set size. NES is comparable
across the four contexts (same normalization), unlike raw expression.

**4. Context definitions** (poster labels used on the figure):
- GSE118713 -- "Cell-line resistance model" (MCF7 vs TAMR)
- GSE111151 -- "Independent resistant sublines" (4 parental lines vs 7 TamR sublines)
- GSE240112 -- "Primary vs recurrent tumours"
- GSE245601 -- "Acute 12 h tamoxifen"

**5. GSE240112 caveat:** recurrence-ASSOCIATED, unpaired human tumour
comparison -- not an experimentally established chronic-resistance model.
Labelled "Primary vs recurrent tumours" on the figure, never "resistant."

**6. GSE245601 caveat:** acute 12h ex vivo tamoxifen exposure, not
resistance. Labelled "Acute 12 h tamoxifen" and visually set apart with a
light background band in Panels A and B -- this is the context that
breaks from the resistance/recurrence pattern in Panel B.

**7. Significance convention:** frozen FDR only, no recomputation.
Panels A and B show all 4 contexts x 3 pathways (12 points total); every
one of those 12 points has frozen FDR well under 0.05 (largest is
gse240112 EMT at FDR = 0.017), so no per-point significance marker was
added there to avoid clutter -- stated here instead, as directed. Panel
C's 3 supporting pathways are more mixed (e.g. TGF-β signaling is not
significant at FDR < 0.05 for the recurrence and acute contexts): those
points use a filled/open convention (filled = FDR < 0.05, open = FDR ≥
0.05), explained once in the panel's own legend. A small, fixed,
per-context vertical offset (never applied to the NES/x-position) is used
in Panel C only to keep two near-identical NES values (TGF-β signaling:
gse111151 = 1.46 vs gse240112 = 1.47) from fully occluding one another --
a standard beeswarm-style separation, not a data change.

**8. Confirmation no pathway result changed:** no NES, FDR, leading-edge
gene, or enrichment ranking was recomputed, re-thresholded, or
re-ordered. The figure reuses the already-frozen, already-tested
`load_pathway_trajectories()` loader unmodified.

**On the optional GSEA-curve alternate:** considered but not built. The
frozen ranked-gene tables and gene-set membership needed to reconstruct a
real running-enrichment-score curve exist (`{dataset}_ranked_genes.tsv`,
`data/reference/genesets/hallmark.gmt`, already used elsewhere via
`build_enrichment_curve()`), so it would have been technically possible
without any new analysis -- but Panels A and B already communicate the
full estrogen-suppression and EMT-reversal story clearly and at poster
distance; a curve view would only re-present the same two conclusions in
a more technical form, which the task explicitly warned against turning
into another dense technical panel.
