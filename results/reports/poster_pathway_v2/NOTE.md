# Pathway Figure v2 -- Short Data Note

**Purpose:** the pathway-analysis step of the poster story, positioned
after the candidate-network section. It answers an orthogonal question:
"do the biological programs implicated by the candidate networks
actually change in the transcriptomic models?" It does NOT claim
candidate -> pathway -> resistance causality; network + pathway
agreement supports a hypothesis only.

**1. Exact frozen source tables** (read unmodified; no new enrichment,
no re-ranking, no recomputation): `results/tables/systems_network/
gsea_{dataset}.tsv` for `gse118713, gse111151, gse240112, gse245601`,
read via the already-tested
`poster_exploration_v2_data.load_pathway_trajectories()` loader.

**2. Exact pathways included** (rows, top to bottom):
- `HALLMARK_ESTROGEN_RESPONSE_EARLY` ("Estrogen response — early")
- `HALLMARK_ESTROGEN_RESPONSE_LATE` ("Estrogen response — late")
- `HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION` ("EMT")
- `HALLMARK_WNT_BETA_CATENIN_SIGNALING` ("WNT / β-catenin")
- `HALLMARK_E2F_TARGETS` ("E2F targets")

**3. Why these pathways:** theme-first selection. The rows are the
directly comparable frozen Hallmark readouts of the biological themes
raised by the network/mechanism section (estrogen response, EMT,
WNT/beta-catenin via USP34's bridge biology, E2F/cell cycle) --
pre-specified before inspecting plot aesthetics. WNT is included despite
being weak/non-significant in two contexts BECAUSE it was a pre-specified
network theme; omitting it for weakness would be selective reporting.
`HALLMARK_G2M_CHECKPOINT` (the optional cell-cycle companion) was
inspected and OMITTED: it mirrors E2F targets' direction (negative in
all four contexts) with weaker signal and is non-significant in
gse240112 (FDR=0.71) and gse245601 (FDR=0.11) -- a redundancy decision,
not significance cherry-picking, since the pre-specified cell-cycle
theme is already covered by E2F targets.

**4. Context definitions** (poster labels; accessions in small gray):
- GSE118713 -- "Cell-line resistance model" (MCF7 vs TAMR)
- GSE111151 -- "Independent resistant sublines"
- GSE240112 -- "Primary vs recurrent tumours"
- GSE245601 -- "Acute 12 h tamoxifen"

**5. NES interpretation:** Normalized Enrichment Score from the frozen
GSEA run. Negative = program suppressed (genes shifted toward the
down-regulated end of that context's ranking); positive = enriched.
Encoded as color (blue = suppressed, terracotta = enriched, the same
diverging palette endpoints as the expression heatmap) and marker size
(|NES|); the legend spells out "Suppressed <- 0 -> Enriched" so no GSEA
background is needed.

**6. Significance convention:** frozen FDR only, threshold 0.05 (same
convention as poster_pathway_v1). Filled circle = FDR < 0.05; open
circle = FDR >= 0.05. Open circles in the figure: WNT/β-catenin in
gse118713 (FDR=0.168) and gse245601 (FDR=0.710). Every other displayed
point is significant.

**7. GSE240112 caveat:** recurrence-ASSOCIATED, unpaired human tumour
comparison -- not a tamoxifen-resistance experiment. Labeled "Primary vs
recurrent tumours", never "resistant".

**8. GSE245601 caveat:** acute 12 h tamoxifen exposure, not established
resistance. Visually separated with a light background band and its own
"Acute response" group caption; the first three columns are grouped as
"Resistance / recurrence-associated".

**9. Theme coverage:** WNT/beta-catenin has a directly comparable frozen
readout (`HALLMARK_WNT_BETA_CATENIN_SIGNALING`, present in all four
tables) and is shown. "Chromatin regulation" -- a genuine network-section
theme (KDM1A/TLK2) -- has NO single comparable Hallmark gene set across
all four frozen tables, and none was invented; that theme therefore has
no row in this figure.

**10. Confirmation:** no enrichment result was changed, recomputed,
re-thresholded, or re-ordered. The loaded frozen values matched the
previously reported directional patterns (estrogen early/late negative
in all four contexts; EMT positive in the three resistance/recurrence
contexts and negative in acute) -- no discrepancy to report.

**Panel count decision:** single panel. The optional "EMT direction"
lollipop panel was considered and NOT built: the EMT matrix row already
shows the resistance/recurrence-up vs acute-down reversal unambiguously
(three warm dots flipping to blue inside the highlighted acute band),
and a second panel would duplicate it.
