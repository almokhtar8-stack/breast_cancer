# GSE240112 cell-type / tumor-cell definitions (Phases 6-7)

**Date written:** 2026-08-12. Documents decisions already frozen in
`docs/GSE240112_PREANALYSIS.md` section G before any candidate-gene value
was inspected; this file adds the supporting composition numbers produced
once the pipelines ran.

## Phase 7: tumor-cell definition (primary population)

**Case A applies.** The `TTs_cancer_060223.h5seurat` object is, by
construction, the authors' own cancer-cell subset: `cell.annot ==
"Breast cancer cells"` for all 9,942 cells, generated (per the paper
Methods) via `Seurat::FindAllMarkers` against curated reference marker
panels (PanglaoDB, CellMarker, ProteinAtlas, literature) -- explicitly not
a CNV-based method. This is used directly as the primary tumor-cell
population with no independent re-classification, no threshold tuning,
and no comparison of alternative malignancy-calling methods, per the
user's explicit time-boxing instruction. Per-sample cell counts:

| Sample | Group | Tumor cells |
|---|---|---|
| PT1 | PT | 1,029 |
| PT2 | PT | 1,975 |
| PT3 | PT | 1,442 |
| RT1 | RT | 2,721 |
| RT2 | RT | 2,597 |
| RT3 | RT | 178 |

RT3's low cell count was flagged before any DE result was inspected
(PREANALYSIS.md section G). Pseudobulk QC (`results/tables/gse240112_pseudobulk/tumor_cell_pca_coordinates.tsv`)
shows RT3 sits cleanly on the RT side of PC1 (85.7% of variance; PT ~
-83 to -84, RT ~ +79 to +92) and correlates more strongly with RT1/RT2
(Spearman 0.92-0.93) than with any PT sample (0.87-0.88) -- it behaves as
a smaller, noisier, but genuine RT replicate, not a technical outlier.
Retained, not excluded.

## Phase 6/14: broad cell-type compartments (secondary population)

No author-provided all-cell-type PT/RT object exists publicly (the only
PT/RT-covering processed object is the cancer-cell-only `TTs_cancer`
object -- see `docs/GSE240112_DATA_AUDIT.md` section 4). Broad
compartments were therefore reconstructed from the raw GEO Cell Ranger
matrices (`scripts/analysis/gse240112_04_cellranger_epithelial.R`):
standard per-sample QC (`min.cells=3`, `min.features=200` at object
creation per the paper's own stated Seurat call; `nFeature_RNA>200`,
`percent.mt<20`; doublet removal via `scDblFinder`), merged without
cross-sample integration (Phase 6 explicitly discourages integrating away
real PT/RT biology, and no scientific justification for it was found),
standard `NormalizeData`/`FindVariableFeatures`/`ScaleData`/`RunPCA`/
`FindNeighbors`/`FindClusters`(res=0.5)/`RunUMAP`. Each cluster was
assigned to the broad compartment (epithelial / immune / stromal /
endothelial) with the highest mean log-normalized expression of a fixed,
non-tuned canonical marker panel (epithelial: EPCAM, KRT8, KRT18, KRT19;
immune: PTPRC; stromal: PDGFRB, COL1A1, DCN; endothelial: PECAM1, VWF) --
this is compartment calling only, not a malignancy classifier, and it
never redefines or replaces the Phase 7 Case A tumor-cell labels used for
the primary analysis.

Per-sample QC funnel and compartment composition: see
`results/tables/gse240112/epithelial_compartment_composition.tsv` and
`results/figures/gse240112_pseudobulk/broad_celltype_umap.png`.

**Finding (not a pipeline bug -- verified directly against raw counts):**
every one of the 15 clusters across all 6 raw PT/RT samples (22,447
QC-passed, doublet-removed cells) was called "epithelial." This was
initially suspected to be a scoring-code bug and was re-verified with a
second, independent, transparent implementation (per-cluster mean
log-normalized marker expression, replacing the first attempt's
`Seurat::AddModuleScore`) -- the result was identical. Direct inspection
of raw counts confirms this is real: the immune marker `PTPRC` (CD45) has
nonzero counts in only 17 of 22,447 cells (0.08%), and stromal
(`PDGFRB`/`COL1A1`/`DCN`) and endothelial (`PECAM1`/`VWF`) markers show
comparably negligible detection, while `EPCAM` is detected in 20,483/
22,447 cells (91.3%). The raw GEO-deposited PT/RT scRNA samples are
therefore themselves essentially epithelial-enriched (consistent with the
authors' own cancer-cell-only processed object requiring comparatively
little further curation from an already largely-epithelial input) --
whether by tissue-dissociation protocol or an unreported sorting step is
not stated in the paper Methods available to this audit.

**Consequence for Phase 14:** the "all-epithelial" sensitivity population
is compositionally close to "all QC-passed cells from these samples" --
it is still a genuinely broader, differently-curated population than the
author's marker-based `TTs_cancer` calls (Phase 7 Case A), so it still
serves Phase 14's stated purpose of reducing dependence on that specific
curation, but it is not an independent non-tumor-epithelial comparator in
the way "epithelial vs. immune/stromal" framing might suggest.
