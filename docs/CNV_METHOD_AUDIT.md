# CNV Method Audit: InferCNV vs. CopyKAT in the GSE245601 pipeline

**Status:** audit-only. No parameter, script, or result was changed while
writing this document. No CRISPR/bulk-RNA candidate gene was inspected.
No InferCNV or CopyKAT run was re-executed.

**Scope:** line-by-line comparison of
[`scripts/analysis/gse245601_05_infercnv_malignant.R`](../scripts/analysis/gse245601_05_infercnv_malignant.R)
and
[`scripts/analysis/gse245601_06_copykat_sensitivity.R`](../scripts/analysis/gse245601_06_copykat_sensitivity.R)
against the official cloned source repositories at `external_refs/inferCNV/`
and `external_refs/copykat/`.

**Audited reference versions** (git clones, read-only):

| Package | Cloned commit | Cloned `DESCRIPTION` version | Version actually installed/run in this project |
|---|---|---|---|
| infercnv | `65e6bf554600c92e7d5c121e98247d9c3f888a0d` (2025-11-14) | 1.23.0 | 1.22.0 (Bioconductor, per prior session record) |
| copykat | `01924c0999bb660064da8d3e1561615ac6efb76f` (2026-06-16) | 1.2.5 | matches the CRAN/GitHub version installed for this project |

**Caveat:** the cloned `infercnv` repo (1.23.0) is a slightly newer commit
than the 1.22.0 build actually used to produce our results. `inst/NEWS`
was checked for changes between these versions; no entry was found that
alters the default value of any parameter audited below (`cutoff`,
`min_cells_per_gene`, `cluster_by_groups`, `scale_data`, `HMM`,
`analysis_mode`, `denoise`, `window_length`, `smooth_method`,
`max_centered_threshold`). The comparison below is therefore treated as
valid for 1.22.0, but this is a documented assumption, not a verified
byte-identical match.

---

## 1. InferCNV audit

### 1.1 What our script does (from `gse245601_05_infercnv_malignant.R`)

- **Input count matrix:** raw UMI counts (`GetAssayData(sample_obj, layer = "counts")`)
  from the merged, QC'd, annotated Seurat object
  (`data/processed/gse245601/seurat_clustered/annotated.rds`), subset to
  one sample's reference + observation cells, then further subset to
  genes present in the gene-order file (`common_genes <- intersect(...)`).
- **Unit of "sample":** `merged$sample_id` is one value per **GSM**, i.e.
  one patient **at one condition** (e.g. `Tumor_01_Control` and
  `Tumor_01_Tamoxifen` are two separate values, confirmed by the 20
  distinct output directories under `data/processed/gse245601/infercnv/`
  and by `results/tables/gse245601_qc_summary.tsv`). InferCNV is run
  **20 times, once per GSM, never combined across patients or across the
  two conditions of the same patient.**
- **Reference (diploid) cells:** immune (`t_nk`, `b_plasma`, `myeloid`) +
  `endothelial` cells, from that same GSM only, identified from the
  candidate-gene-blind marker-based `broad_cell_type` annotation (Step 6).
- **Observation cells:** `broad_cell_type == "epithelial"` cells, same GSM only.
- **Gene-order file:** `data/processed/gse245601/infercnv/gene_order.tsv`,
  built by `scripts/analysis/gse245601_04_gene_order_file.R` from
  `TxDb.Hsapiens.UCSC.hg38.knownGene` + `org.Hs.eg.db`; contains 24 distinct
  chromosome values including `chrX` (836 genes) and `chrY` (86 genes) —
  see §1.3 on why these are dropped anyway.
- **Explicitly set `run()` parameters:** `cutoff=0.1`, `min_cells_per_gene=3`,
  `cluster_by_groups=TRUE`, `scale_data=FALSE`, `HMM=FALSE`,
  `analysis_mode="samples"`, `denoise=TRUE`, `no_prelim_plot=TRUE`,
  `png_res=60`, `num_threads=8`, `plot_steps=FALSE`, `resume_mode=FALSE`,
  `save_rds=FALSE`.
- **Explicitly set `CreateInfercnvObject()` arguments:** `raw_counts_matrix`,
  `annotations_file`, `gene_order_file`, `ref_group_names="reference"`
  (a single reference group; no tumor/patient-specific reference splitting).

### 1.2 Parameter table

| Parameter | Our value | Official default / documented example | Same / different | Likely effect |
|---|---|---|---|---|
| `raw_counts_matrix` | Seurat `counts` layer (raw UMI, 10x) | matrix/data.frame/file accepted; no "correct" default | n/a | Correct input type (raw counts, not normalized) — matches the package's stated requirement. |
| `gene_order_file` | TxDb-derived, chr/start/end per gene | user-supplied file, any source | n/a | Standard, documented way to build this file when one isn't separately downloaded (no canonical "official" gene-order file ships with the package for GRCh38). |
| `ref_group_names` | `"reference"` (one group: immune+endothelial) | any string(s) matching values in the annotations file | n/a | Single reference group, not per-cell-type-separated references. Simpler than possible; not wrong, but means InferCNV can't separately validate against multiple normal cell types. |
| `delim` (CreateInfercnvObject) | not set → `"\t"` | `"\t"` | same (default) | — |
| `max_cells_per_group` | not set → `NULL` | `NULL` (no subsampling) | same (default) | All reference/observation cells used, no downsampling. |
| `min_max_counts_per_cell` | not set → `c(100, +Inf)` | `c(100, +Inf)` | same (default) | Cells with <100 total counts across gene-order genes would be dropped; irrelevant in practice since our QC already requires nCount_RNA≥5000 per cell. |
| `chr_exclude` | not set → `c('chrX','chrY','chrM')` | `c('chrX','chrY','chrM')` | same (default) | **chrX/chrY genes present in our gene-order file (836+86=922 genes) are silently dropped by InferCNV itself** at `CreateInfercnvObject()` time. This is package-default behavior, not something our script configures, and is standard practice (avoids X-inactivation/sex-chromosome dosage confounds). |
| `cutoff` | `0.1` | function default `1` (Smart-seq2); **vignette explicitly states `cutoff=0.1` for 10x Genomics** | different from the coded default, **same as the documented value for our data type** | Correct choice for 10x/droplet data (our data is 10x GRCh38 Cell Ranger output); `cutoff=1` would be wrong here. |
| `min_cells_per_gene` | `3` | `3` | same (default) | — |
| `cluster_by_groups` | `TRUE` | `TRUE` (current default; changed to `TRUE` from `FALSE` in v1.15.1 per `inst/NEWS`) | same | Reference and observation cells are each clustered within their own group for heatmap/image ordering (does not itself classify malignancy). |
| `scale_data` | `FALSE` | `FALSE` | same (default) | No additional cross-cell scaling beyond InferCNV's own normalization. |
| `HMM` | `FALSE` | `FALSE` | same (default) | No hidden-Markov discrete CNV-state calling; output stays as continuous smoothed relative expression only. |
| `analysis_mode` | `"samples"` | programmatic (`match.arg`) default is `"subclusters"` (first-listed value); the parameter's own Rd documentation states "default: samples (fastest, but **subclusters is ideal**)", and the package's own worked `run.Rd` example explicitly sets `analysis_mode="samples"` | matches the value the Rd prose calls "default" and the worked example uses — but the same Rd text explicitly says the *other* mode ("subclusters") is "ideal", i.e. more thorough. So `"samples"` is a documented, explicitly-exampled choice, not an undocumented one — but it is the faster/lighter option, not the one the package authors call best. | No per-tumor subclustering is performed; each observation group (all epithelial cells of one GSM) is treated as one block for image-filtering purposes. This choice does not, by itself, invalidate the analysis, but it does mean any within-sample tumor heterogeneity is not resolved by InferCNV's own subclustering; that responsibility falls entirely to our downstream per-(GSM, Seurat-cluster) thresholding logic (§1.3). |
| `denoise` | `TRUE` | `FALSE` (coded default); **the package's own vignette "full default analysis" example sets `denoise=TRUE`** | different from the coded default, **identical to the vignette's example** | Values within `sd_amplifier` (default 1.5) SD of the reference mean are set to zero/neutral ("whitening"), reducing background noise before our CNV-score computation. |
| `window_length` | not set → `101` | `101` | same (default) | 101-gene pyramidal smoothing window along chromosome order. |
| `smooth_method` | not set → `"pyramidinal"` (first of `c('pyramidinal','runmeans','coordinates')`) | `"pyramidinal"` | same (default) | — |
| `max_centered_threshold` | not set → `3` | `3` | same (default) | Final relative-expression values are clipped to `[1-3, 1+3]` (i.e. `[-2, 4]`) before being reported. |
| `num_ref_groups` | not set → `NULL` | `NULL` | same (default) | — |
| `no_prelim_plot` | `TRUE` | `FALSE` (coded default); vignette example sets `TRUE` | different from coded default, matches vignette example | Suppresses only the preliminary (pre-final-processing) diagnostic plot; the **final** `infercnv.png` heatmap is still produced (confirmed present on disk, §4). |
| `png_res` | `60` | `300` (coded default); vignette example sets `60` | different from coded default, matches vignette example | Lower-resolution PNG (600×566 px, confirmed on disk) — a speed/disk tradeoff from the package's own example, not a scientific parameter. |
| `num_threads` | `8` | `4` | different | Performance-only; does not change the computed result. |
| `resume_mode` | `FALSE` | `TRUE` | different | Forces a fresh run rather than resuming from any leftover per-step `.rds` checkpoint; relevant because each of our 20 per-GSM output directories is used exactly once. Not a scientific parameter. |
| `save_rds` | `FALSE` | `TRUE` | different | Suppresses ~20 intermediate per-processing-step `.rds` checkpoint files. **`save_final_rds` was left at its default (`TRUE`, not overridden by our script)**, so the final object is still written to disk as `run.final.infercnv_obj` (confirmed present, 61 MB, §4) — no scientific output is lost by this choice. |
| `plot_steps` | `FALSE` | `FALSE` | same (default) | — |

**Bottom line on parameters:** not every scientific parameter is at the
package's coded default — `cutoff`, `denoise`, and `analysis_mode` are
not. But each of those differs from the coded default in a way that
matches a value the package's own vignette/documentation explicitly uses
or recommends for this exact scenario (10x data, the vignette's "full
default analysis" example), rather than an arbitrary or undocumented
value. `no_prelim_plot`, `png_res`, `num_threads`, `resume_mode`, and
`save_rds` are non-scientific I/O/performance/output settings that do not
change the CNV computation itself. No parameter was found that diverges
from both the coded default *and* every documented example/recommendation
in a way that would plausibly change the underlying computation.

### 1.3 Part A vs. Part B: what InferCNV itself produces vs. our downstream logic

This distinction is the central methodological fact of this audit.

**Part A — produced directly by `infercnv::run()`:**
InferCNV performs reference normalization, gene filtering, smoothing,
per-group clustering/ordering, denoising, and plotting; the numeric
result our script actually reads and builds on is the `infercnv-class`
object's `@expr.data` slot: a gene × cell matrix of relative expression,
per gene order-sorted along the genome, smoothed (101-gene pyramidal
window), reference-centered (values near `1.0` = no change vs. the
reference group's average), denoised (near-reference values whitened to
exactly neutral), and clipped to `max_centered_threshold`. **This is the
package output our downstream classification consumes** — it is not a
claim that InferCNV's internal processing (clustering, plotting, image
generation) is nothing more than producing this matrix, only that
`@expr.data` is the only piece of it our script uses. With
`analysis_mode="samples"` and `HMM=FALSE` (our settings), InferCNV does
**not** perform tumor subclustering, does **not** run an HMM CNV-state
caller, and — critically — **does not itself output any
malignant/non-malignant or tumor/normal cell-level classification.**
`cluster_by_groups=TRUE` controls whether cells are hierarchically
clustered/ordered *within* the pre-supplied reference/observation groups
(affecting image ordering and, per the option's documented role, the
per-group image filtering) — it is a within-group ordering/filtering
setting, not a classifier, and it does not produce or influence any
malignant/non-malignant label.

**Part B — our own downstream logic (lines 161–216 of the script,
reconstructed from the original paper's author code, not part of the
`infercnv` R package at all):**
1. `expr_mat <- run_result@expr.data[, obs_cells]` — pull out only the
   epithelial (observation) columns of InferCNV's own output.
2. `cnv_score = colMeans((expr_mat - 1)^2)` — our own per-cell summary
   statistic (mean squared deviation from neutral). Not computed or
   suggested by `infercnv` itself.
3. Per (GSM, Seurat cluster) group (with a documented small-cluster
   pooling fallback, `MIN_CLUSTER_CELLS_FOR_OWN_THRESHOLD=10`): take the
   top 5% of cells (minimum 2) by `cnv_score` as a malignant "seed", and
   compute each cell's Kendall-tau correlation to the seed's mean profile.
4. Adaptive per-group thresholds: `th_value = mean(score) - 2*SD`
   (clamped to `[0.01, 0.05]`), `th_corr = mean(corr) - 1.5*SD` (clamped
   to `[0.2, 0.4]`). A cell is labeled `"malignant"` iff both thresholds
   are exceeded, else `"non-malignant epithelial"`.

None of step 2–4 is InferCNV package code. It is this project's
reconstruction (per `docs/gse245601_PREANALYSIS.md` §9, "ng_2021_and_thresholding")
of a separate, non-package downstream classification script the original
paper's authors used on top of raw InferCNV output. **Every malignant/
non-malignant label in `results/tables/gse245601_malignant_cell_labels.tsv`
is Part B, not a native InferCNV output.**

---

## 2. CopyKAT audit

### 2.1 What our script does (from `gse245601_06_copykat_sensitivity.R`)

- **Input count matrix:** raw UMI counts (`GetAssayData(sample_obj, layer = "counts")`,
  coerced to a dense `matrix`), same source Seurat object as InferCNV.
- **Per-sample execution:** each of the 20 GSMs (patient × condition) is
  run as an independent `copykat()` call, in its own working directory
  (`setwd(sample_dir)` / `setwd(orig_wd)`), matching the package README's
  explicit guidance: *"It is suggested to run one sample at a time.
  Combining different samples would pick up batch effects."*
- **Input cell population:** the *same* reference (immune+endothelial) +
  epithelial cells that InferCNV used for that GSM (`copykat_input_cells <- c(ref_cells, epi_cells)`).
- **Known normal cells:** **not supplied.** `norm.cell.names` is left at
  its default (`""`), so CopyKAT performs its own internal
  confident-normal-cell/baseline inference from the mixed population,
  rather than being told which cells are reference — this is what makes
  it methodologically independent of InferCNV's externally-supplied
  reference, not merely a different implementation of the same idea.

### 2.2 Parameter table

| Parameter | Our value | Official default / documented example | Same / different | Likely effect |
|---|---|---|---|---|
| `rawmat` | raw UMI counts, gene symbols × cells | required, no default | n/a | Matches the documented requirement ("raw data matrix; genes in rows; cell names in columns"). |
| `id.type` | `"S"` | `"S"` | same (default) | Gene symbols (matches Cell Ranger's default feature naming). |
| `cell.line` | `"no"` | `"no"` | same (default) | Correctly declares this as mixed tumor+normal data, not a pure cell line. |
| `ngene.chr` | `5` | `5` | same (default; also the README's own worked example uses `5`) | Minimum 5 genes per chromosome required for a cell to be classified; matches documented guidance ("using at least 5 genes ... is not very stringent"). |
| `min.gene.per.cell` | not set → `200` | `200` | same (default) | — |
| `LOW.DR` | not set → `0.05` | `0.05` | same (default) | Note: the package README's *prose* states `UP.DR=0.2` as default, but the actual `copykat()` function signature codes `UP.DR=0.1` — a minor inconsistency in the **official repo's own documentation**, not in our script. Our script does not override `UP.DR`, so it uses whatever the installed package's coded default is. |
| `UP.DR` | not set → `0.1` (per function signature) | `0.1` (code) / `0.2` (README prose — inconsistent) | same as coded default | See note above. |
| `win.size` | `25` | `25` | same (default; README worked example also uses `25`) | Minimum 25 genes per segmentation bin. |
| `norm.cell.names` | not set → `""` (empty) | `""` (empty = automatic detection) | same (default) | **Intentional and central to this being an independent check** — CopyKAT infers its own normal/baseline cells rather than reusing our marker-based reference. |
| `KS.cut` | `0.1` | `0.1` (coded default); **README's own worked example also uses `0.1`**, and states "usually it works in a range of 0.05-0.15"; a separate `vignettes/copycat-vignettes.Rmd` example uses `0.2` | same as coded default and as the README's worked example; within the package's own stated usual range either way | Segmentation sensitivity: `0.1` gives more/finer segments (higher sensitivity) than `0.2`. Both are within the range the package documents as reasonable; `0.1` is the more conservative, more literally "default" choice. |
| `sam.name` | `sample_id` (e.g. `"Tumor_01_Control"`) | any string | n/a | Output file-naming prefix only, not scientific. |
| `distance` | `"euclidean"` | `"euclidean"` | same (default; README worked example also uses `"euclidean"`) | — |
| `test.emd` | not set → `"FALSE"` | `"FALSE"` | same (default) | — |
| `output.seg` | `FALSE` (logical) | `"FALSE"` (string, coded default) | same intended (off) setting — our script passes a logical `FALSE` where the function's coded default is the string `"FALSE"`; R's internal comparisons coerce these equivalently in this codepath, so this is not a scientific difference, but it is not literally the identical type as the coded default. | — |
| `plot.genes` | `FALSE` | `"TRUE"` | **different** | Suppresses CopyKAT's optional per-gene-labeled heatmap variant only; the main `_copykat_bin_by_cell_heatmap.jpeg` output is unaffected and confirmed present for all 20 samples (§4). Not a scientific parameter — output/plotting only. |
| `genome` | `"hg20"` | `"hg20"` | same (default) | hg20 = GRCh38 (package's naming convention); consistent with our GRCh38 10x reference. |
| `n.cores` | `8` | `1` | different | Performance-only; does not change the computed result. |

**Bottom line on parameters:** every scientific parameter (`id.type`,
`cell.line`, `ngene.chr`, `min.gene.per.cell`, `LOW.DR`, `UP.DR`,
`win.size`, `norm.cell.names`, `KS.cut`, `distance`, `genome`) is at the
package's own coded default, and where the package's own documentation
gives a worked example, our values match that example too. The only
deviations (`plot.genes=FALSE`, `n.cores=8`, `sam.name`) are output/
performance/naming settings with no effect on the CNV computation or the
aneuploid/diploid/not.defined classification itself.

### 2.3 How `aneuploid` / `diploid` / `not.defined` are produced and used

**Produced natively by `copykat()` itself** (confirmed in `R/copykat.R`,
e.g. lines 382–492 and 713–820): CopyKAT clusters all input cells by
their segmented CNV profiles, identifies the cluster with the lowest
overall CNV signal as the "confident diploid" cluster and the
highest-signal cluster as "confident aneuploid", checks the correlation
between their consensus profiles (`cor(conses.diploid, conses.aneuploid)`)
against empirical cutoffs of `0.4` and `0.6` to gauge separation
confidence, then assigns every remaining cell to whichever consensus
profile it is closer to by Earth Mover's/Wasserstein distance
(`transport::wasserstein1d`). Cells that fail CopyKAT's own upstream
per-cell QC (`ngene.chr`, `min.gene.per.cell`) are explicitly re-added to
the returned prediction table with the literal label `"not.defined"`
(`R/copykat.R` line 490-492: `ndef <- colnames(rawmat)[... %!in% names(com.preN)]`,
then `rep("not.defined", length(ndef))`) — in the two code branches shown
in `R/copykat.R` (the two internal execution paths this version of the
function can take), every input cell that fails to receive a
diploid/aneuploid call is explicitly re-added with the `"not.defined"`
label rather than silently dropped from the returned table. This was also
verified empirically for our actual 20 runs (§2.3 below, and §4): the
count of epithelial cells with a CopyKAT-derived label plus the count of
`not.defined` cells accounts for 100% of epithelial cells in every GSM.
This **is** CopyKAT's classification — unlike InferCNV, there is no
downstream reconstruction step for CopyKAT; `ck_result$prediction$copykat.pred`
is used directly.

**How our script uses it** (`gse245601_06_copykat_sensitivity.R`, lines
95–102): a simple relabeling —
`aneuploid → "malignant"`, `diploid → "non-malignant epithelial"`,
anything else → `"not_defined"`. Inspecting the actual raw
`copykat.pred` values written to
`data/processed/gse245601/copykat/*/*_copykat_prediction.txt` for all 20
samples confirms only three literal values ever occur in this dataset:
`aneuploid`, `diploid`, and `not.defined` (2 cells total, both in
`Tumor_09_*`) — the package's other possible low-confidence labels
(`"c1:diploid:low.conf"`, `"c2:aneuploid:low.conf"`, which the code path
in §2.3 can in principle emit when the diploid/aneuploid consensus
correlation falls in `[0.4, 0.6)`) **do not appear anywhere in this
project's actual CopyKAT output.** So while our catch-all
`ifelse(...,"not_defined")` would also silently absorb any such
low-confidence label into `"not_defined"` if one occurred, this did not
in fact happen for any of the 20 GSMs — confirmed by direct inspection of
the output files, not merely by reading the script.

---

## 3. Methodological comparison: InferCNV vs. CopyKAT

| Aspect | InferCNV (as run here) | CopyKAT (as run here) |
|---|---|---|
| **Normal/reference baseline** | Externally supplied: marker-based `broad_cell_type` immune+endothelial cells from the *same GSM*, fixed before any InferCNV computation. Never re-estimated by the algorithm itself. | Internally inferred: CopyKAT clusters *all* input cells (reference+epithelial, mixed) by CNV signal and picks its own "confident diploid" cluster; our marker-based labels are deliberately **not** given to it (`norm.cell.names=""`). |
| **Chromosome ordering** | Per-gene, from our TxDb-derived `gene_order.tsv` (chr/start/end); chrX/chrY/chrM dropped by package default. | Per ~220 kb genomic bin, from CopyKAT's own bundled hg20 gene-position annotation (`annotateGenes.hg20.R`) — the "220KB windows" bin size is stated in `copykat.Rd`'s own `@return` description (not independently re-derived in this audit). |
| **Smoothing** | Continuous 101-gene pyramidal (triangular-weighted) moving average along chromosome order (`window_length=101`, `smooth_method="pyramidinal"`), as a step separate from any segmentation (none is performed here — `HMM=FALSE`). | A KS-test-based iterative segmentation procedure, with a minimum window of `win.size=25` genes per segment; per the package's own documentation this single procedure both smooths the per-bin signal and defines discrete breakpoints — the audit did not independently verify from source whether CopyKAT has a wholly separate smoothing pass, only that segmentation and effective noise-reduction are governed by the same `win.size`/`KS.cut` parameters. |
| **Segmentation** | Not performed by default (no discrete breakpoints) unless `HMM=TRUE` (not used here — `HMM=FALSE`). Output stays continuous. | Explicit discrete segmentation via a Kolmogorov-Smirnov test (`KS.cut` controls sensitivity/number of breakpoints). |
| **CNV estimation** | Relative expression per gene per cell vs. the reference group's mean, log-residual, smoothed, denoised (`denoise=TRUE`, whitens near-neutral values), clipped to `max_centered_threshold=3`. Centered at `1.0` = no change. | Log2 ratio per genomic bin vs. the internally-estimated baseline, after segmentation. |
| **Cell clustering** | With `analysis_mode="samples"` (our setting): **no** tumor subclustering performed by InferCNV; `cluster_by_groups=TRUE` only orders cells for the heatmap dendrogram within the pre-supplied reference/observation groups. | Hierarchical clustering of cells by their segmented CNV profiles, using the specified `distance="euclidean"` metric, to separate "confident diploid" vs. "confident aneuploid" clusters. |
| **Malignant/aneuploid classification** | **Not produced by InferCNV at all** in this configuration. Entirely our own downstream reconstruction (Part B, §1.3): CNV score + Kendall-tau-to-seed correlation + adaptive per-(GSM, Seurat-cluster) thresholds. | **Produced natively by `copykat()`** as `copykat.pred` ∈ {`aneuploid`, `diploid`, `not.defined`, and in principle the two `low.conf` labels — none observed here} — a single, package-internal, unsupervised two-cluster-plus-distance decision. |
| **Is classification built into the method or added downstream?** | **Added downstream** (by us / the reconstructed paper procedure). InferCNV's own scope ends at the smoothed relative-expression matrix. | **Built into the method.** CopyKAT's `prediction` output is its final product; no further modeling is layered on top by our script beyond a label rename. |

### Why the two methods can disagree even when both are run correctly

This is a direct, structural consequence of the differences above, not
evidence that either tool was misconfigured:

1. **Different reference definitions.** InferCNV's reference is
   externally fixed (marker-defined immune+endothelial cells); CopyKAT's
   reference is whatever cluster its own internal algorithm judges to be
   the lowest-CNV-signal group among *all* input cells, including the
   epithelial population. If a sample's true tumor content is very low
   or very high, or if the epithelial and stromal/immune populations are
   not cleanly separable in CopyKAT's own bin-level CNV space, CopyKAT's
   inferred "confident diploid" cluster can end up composed differently
   from InferCNV's marker-based reference — with no way for the two
   methods to agree by construction in that scenario.
2. **One classification is native, the other is reconstructed.**
   CopyKAT's aneuploid/diploid call is the package's own tested,
   internally-consistent output. InferCNV's malignant/non-malignant call
   is *our* reconstruction of a separate, non-package script (Part B).
   Any imprecision in reconstructing that seed-selection/thresholding
   procedure (documented as best-effort in `docs/gse245601_PREANALYSIS.md`
   §9, tier **[C]** in places) is a plausible, code-supported
   contributor to disagreement that has no counterpart on the CopyKAT
   side.
3. **Different underlying signal representations.** InferCNV outputs a
   continuous, gene-level, pyramidally-smoothed relative-expression
   profile. CopyKAT outputs a discretely-segmented, bin-level log-ratio
   profile derived via a KS-test. These are not the same numerical object
   being thresholded two different ways — they are two different
   derived signals from the same raw counts.

---

## 4. Existing saved outputs — what we already have

All paths below are under `data/processed/gse245601/` (gitignored,
present locally only) unless noted.

### 4.1 InferCNV

For **all 20 GSMs** (`Tumor_01_Control` … `Tumor_10_Tamoxifen`), each
per-GSM directory `infercnv/<GSM>/` contains:

- `infercnv.png` — **the official InferCNV chromosome-scale CNV heatmap**,
  produced directly by `infercnv::run()` (confirmed: `no_prelim_plot=TRUE`
  suppresses only the *preliminary* plot; the final plot is unaffected).
  600×566 px (low-res, `png_res=60`, matching the package's own example).
  Shows reference cells (top) and observation/epithelial cells (bottom),
  each internally ordered by `cluster_by_groups=TRUE` hierarchical
  clustering, across chromosome-ordered genes.
- `infercnv.heatmap_thresholds.txt` — the color-scale breakpoints used
  for that heatmap (16 values, ~0.80 to ~1.20, centered on 1.0).
- `infercnv.observation_groupings.txt` — per-cell dendrogram
  group/annotation-group/color assignment used to draw the heatmap
  (this is InferCNV's own internal grouping for the plot, **not** our
  malignant/non-malignant label).
- `run.final.infercnv_obj` (~62 MB) — the complete final `infercnv`
  object (`save_final_rds` left at its default `TRUE`), from which
  `@expr.data` and the full gene/cell CNV matrix can be recovered without
  re-running InferCNV.
- `gene_order_used.tsv`, `annotations.tsv` — the exact per-run gene-order
  and reference/observation assignment inputs (provenance).

**What this heatmap does *not* show:** our downstream malignant vs.
non-malignant label (that classification did not exist yet when the plot
was drawn — it is computed afterward, in-memory, from `@expr.data`). The
heatmap's own color bars reflect InferCNV's internal
dendrogram/subgroup structure only.

### 4.2 CopyKAT

For **all 20 GSMs**, each per-GSM directory `copykat/<GSM>/` contains:

- `<GSM>_copykat_bin_by_cell_heatmap.jpeg` — **the official CopyKAT
  chromosome-scale (bin-level) CNV heatmap**, produced directly by
  `copykat()`. Confirmed present for all 20 GSMs.
- `<GSM>_copykat_prediction.txt` — the per-cell `aneuploid`/`diploid`/
  `not.defined` calls (this is what `ck_result$prediction` in our script
  reads; also what feeds `results/tables/gse245601_copykat_sensitivity_labels.tsv`).
- `<GSM>_copykat_raw_results_gene_by_cell.txt`,
  `<GSM>_copykat_raw_results_bin_by_cell.txt`,
  `<GSM>_copykat_final_results_bin_by_cell.txt` — the underlying CNV
  ratio matrices at gene and bin level, pre- and post-segmentation
  (each ~550–860 MB per GSM).

### 4.3 Can the same epithelial cells currently be compared side-by-side?

**Cell-identity level: yes, already.** Both scripts derive `obs_cells` /
`epi_cells` identically (`sample_obj$broad_cell_type == "epithelial"` on
the same `annotated.rds` object, same `sample_id`/GSM), and CopyKAT's
`cell.names` output is a verbatim echo of the input matrix's column
names — so the per-cell `primary_malignancy_label` (InferCNV-derived) and
`sensitivity_malignancy_label` (CopyKAT-derived) tables already share
identical `cell_id` strings for every epithelial cell in a GSM (this is
exactly what makes the existing concordance calculation in §5 possible
without any additional barcode reconciliation work).

**Visual/heatmap level: no, not yet, for two reasons:**

1. **Different plotting engines, different cell orderings, different
   color scales.** `infercnv.png` orders cells by InferCNV's own
   `cluster_by_groups` dendrogram and uses a linear 1.0-centered color
   scale with `infercnv.heatmap_thresholds.txt`'s breakpoints. The
   CopyKAT `_bin_by_cell_heatmap.jpeg` orders cells by CopyKAT's own
   hierarchical clustering and uses its own `heatmap.3.R`-based color
   scheme and gene-bin resolution. Neither plot currently marks cells by
   our downstream `primary_malignancy_label` / `sensitivity_malignancy_label`,
   and neither is restricted to only the malignant-labeled subset — both
   show every reference+observation cell that was fed to that run.
2. **No existing figure restricts either heatmap to one tumor's two
   conditions side by side**, or overlays/orders cells by our two derived
   labels for direct visual agreement/disagreement inspection.

### 4.4 What would need to be generated for a fair side-by-side comparison

(Described here for planning only — **nothing below was generated in
this audit.**)

- A **new plotting step** (not yet written) that, per selected GSM(s):
  (a) reads `run.final.infercnv_obj`'s `@expr.data` for the same
  epithelial cells used in the CopyKAT bin-level matrix, (b) reads the
  corresponding CopyKAT `*_copykat_raw_results_gene_by_cell.txt` (or
  `bin_by_cell`) matrix for those same `cell_id`s, (c) draws both as
  chromosome-ordered heatmaps with a **shared cell ordering** (e.g. both
  ordered by `primary_malignancy_label` or both by a common hierarchical
  clustering) and a **shared color scale convention**, with a color side
  bar showing `primary_malignancy_label` and `sensitivity_malignancy_label`
  side by side so agreement/disagreement is visible directly on the plot.
- This requires reading (not recomputing) `run.final.infercnv_obj` and
  the existing CopyKAT `*_results_gene_by_cell.txt` files — both already
  exist on disk (§4.1–4.2), so no re-run of either tool would be needed
  to build this comparison, only a new visualization script.

---

## 5. Tracing the ~54.1% InferCNV/CopyKAT agreement

### 5.1 Where the number comes from

Computed in `gse245601_06_copykat_sensitivity.R`, lines 104–117, per GSM,
then summed across all 20 GSMs at lines 131–133:

```r
merged_compare <- merge(
  primary_labels[primary_labels$sample_id == sample_id, c("cell_id", "primary_malignancy_label")],
  epi_pred[, c("cell_id", "sensitivity_malignancy_label")],
  by = "cell_id"
)
comparable <- merged_compare[merged_compare$sensitivity_malignancy_label != "not_defined", ]
n_concordant <- sum(comparable$primary_malignancy_label == comparable$sensitivity_malignancy_label)
n_compared <- nrow(comparable)
...
overall_concordant <- sum(concordance$n_concordant, na.rm = TRUE)
overall_compared <- sum(concordance$n_compared, na.rm = TRUE)
# OVERALL concordance: overall_concordant / overall_compared
```

Re-computed directly from `results/tables/gse245601_malignancy_concordance.tsv`
during this audit:

- `sum(n_concordant)` = 15,777
- `sum(n_compared)` = 29,173
- `15777 / 29173` = **0.5408 = 54.08%**, matching the previously reported
  "~54.1%" exactly.

**Note on the README's separate "~56%" figure:** that number is a
different, also-valid statistic — the **unweighted mean of the 20
per-sample `concordance_rate` values** (0.5645 = 56.4%), vs. the 54.1%
figure which is the **cell-count-weighted (pooled/micro-average)** rate.
Both are computed correctly from the same table; they differ only in
aggregation method (macro- vs. micro-average) and are not in conflict.

### 5.2 Denominator, exclusions, and barcode matching — verified

- **Denominator = comparable epithelial cells only**, i.e. cells that (a)
  are epithelial, (b) have an InferCNV-derived `primary_malignancy_label`
  for that GSM, and (c) do **not** have CopyKAT's `not_defined` catch-all
  label (§2.3).
- **Total epithelial cells across all 20 GSMs:** 29,175
  (`results/tables/gse245601_cell_metadata_frozen.tsv.gz`, confirmed
  during this audit).
- **Total `n_compared`:** 29,173.
- **Difference: exactly 2** — and direct inspection of all 20
  `*_copykat_prediction.txt` files during this audit found **exactly 2**
  cells labeled `not.defined` by CopyKAT among epithelial cells (1 in
  `Tumor_09_Control`, 1 in `Tumor_09_Tamoxifen`). This is an exact match,
  confirming: **no sample was dropped, no cells were lost to a
  barcode-matching failure, and the only cells excluded from the
  denominator are the ones CopyKAT itself explicitly declined to
  classify.**
- **All 20 GSMs have `status = "ok"`** in both
  `gse245601_malignant_summary_per_sample.tsv` (InferCNV) and
  `gse245601_malignancy_concordance.tsv` (CopyKAT) — no sample hit the
  `insufficient_reference`, `infercnv_run_failed`, `too_few_cells`, or
  `copykat_failed` fallback paths, so none of those exclusion rules
  contributed to the 54.1% figure.
- **Barcode/cell-ID matching:** both `primary_malignancy_label` and
  `sensitivity_malignancy_label` derive `cell_id` from the same Seurat
  object's `colnames()` for the same GSM (§4.3); cell IDs are prefixed
  with `sample_id` (e.g. `Tumor_01_Control_TGCATCCGTGCATGTT-1`), so
  cross-sample collisions are structurally impossible, and the `merge()`
  is additionally pre-restricted to one `sample_id` per iteration. No
  evidence of a matching bug was found.

### 5.3 What plausibly explains the disagreement

Three candidate explanations were checked against the code/output; only
one is supported:

1. **Barcode mismatch or an implementation bug — ruled out.** §5.2 shows
   an exact, fully-reconciled accounting of the denominator (29,175 total
   epithelial → 29,173 compared, difference exactly matches CopyKAT's own
   2 `not.defined` calls). No sample failed, no cells were silently
   dropped, no cross-sample ID collision is possible.
2. **Parameter misconfiguration — not supported.** §1.2 and §2.2 show that
   for every scientific parameter that is *not* at the package's coded
   default (InferCNV's `cutoff`, `denoise`, `analysis_mode`), the value
   used instead matches a value the package's own documentation
   explicitly uses or recommends for this data type — not an arbitrary or
   undocumented value. Neither tool was run with an unusual or
   undocumented setting.
3. **Intrinsic algorithmic differences, amplified by our reconstructed
   InferCNV threshold — supported by direct evidence.** The single most
   striking case in the data is `Tumor_01`: InferCNV's downstream
   procedure calls essentially no cells malignant in either arm (1/3,334
   epithelial cells in `Control`, 2/2,774 in `Tamoxifen` —
   `results/tables/gse245601_malignant_summary_per_sample.tsv`), while
   CopyKAT calls 3,369/3,605 input cells (~93%) `aneuploid` in `Control`
   but 0/2,968 `aneuploid` in `Tamoxifen` (all `diploid`) — confirmed by
   direct inspection of `Tumor_01_Control_copykat_prediction.txt` and
   `Tumor_01_Tamoxifen_copykat_prediction.txt` during this audit. This
   produces the sample's near-0% (`Control`, 34/3,334 concordant) and
   near-100% (`Tamoxifen`, 2,772/2,774 concordant, because both methods
   happen to agree "almost nothing is malignant" there) concordance pair
   seen in `gse245601_malignancy_concordance.tsv`. This pattern —
   near-total disagreement in one arm of a patient and near-total
   agreement in the other arm of the *same* patient, using the *same*
   pipeline code both times — rules out a global misconfiguration or a
   systematic pipeline bug as the explanation (since the same script,
   with the same parameters, produced both outcomes for the same
   patient), and is consistent with **the hypothesis that the
   disagreement tracks which baseline/classification each method happened
   to settle on for that specific sample's cell mixture** (§3, points
   1–2: CopyKAT's internally-inferred baseline can differ from InferCNV's
   externally-fixed reference, and InferCNV's own malignant call is our
   downstream reconstruction, not a package output). This is a plausible,
   code-supported explanation, not a proven one — see the next paragraph
   for what would be needed to confirm it mechanistically.

**What this audit cannot determine from the code alone:** *why* CopyKAT's
internal baseline-clustering settled on "almost everything is aneuploid"
in `Tumor_01_Control` specifically (e.g., whether its own
confident-diploid cluster was small/ambiguous in that sample) would
require inspecting CopyKAT's intermediate clustering objects for that one
sample — data that exists on disk (`*_copykat_raw_results_*.txt`) but was
not opened in this pass, consistent with the "audit-only, no new
analysis" instruction. This would be a reasonable, narrowly-scoped next
diagnostic step, not a re-run of either tool.

---

## 6. Codex independent review

**Verdict: PASS WITH NOTES.**

Codex (via `mcp__codex__codex`) was given the full text of this document,
the two full pipeline scripts, and pasted excerpts of the official
`inferCNV`/`copykat` source (function signatures, Rd docs, vignette/README
worked examples, the `not.defined`/aneuploid/diploid decision code) —
its own filesystem sandbox could not access this repository directly
(`bwrap: pivot_root: Invalid argument`, consistent with this project's
established Codex-review limitation), so the inline-evidence fallback
was used. Codex confirmed: all directly checkable script parameter
values matched the two pasted scripts; the InferCNV
native/non-classifying vs. our-downstream-reconstruction distinction was
correct; CopyKAT's native classification claim was correct; and the
54.08%/56.4%/29,175-vs-29,173/Tumor_01 numerical claims were internally
consistent with the pasted evidence.

Codex flagged the following, all of which were corrected in this document
before this section was written:

1. **`analysis_mode`** — the audit had characterized `"samples"` as
   matching "the package's own documented recommendation"; corrected to
   note the same Rd text explicitly calls the *other* mode
   ("subclusters") "ideal" — `"samples"` is documented and exampled, not
   the package's stated best option.
2. **"`@expr.data` is the entirety of what InferCNV computes"** — too
   broad; corrected to state `@expr.data` is the specific output our
   downstream script consumes, not a claim that InferCNV's internal
   processing (clustering, plotting, etc.) amounts to nothing more than
   producing that matrix.
3. **`cluster_by_groups=TRUE` description** — corrected from "only
   affects the heatmap dendrogram/color bars" to also note its documented
   role in per-group image filtering/ordering, while keeping the correct
   conclusion that it does not produce a malignant/normal label.
4. **`output.seg` table row** — corrected to note our script passes a
   logical `FALSE` where the coded default is the string `"FALSE"`; same
   intended setting, not a literally identical default.
5. **"Every input cell always appears" (CopyKAT)** — corrected to scope
   this claim to the two code branches actually shown and to this
   project's own empirically-verified 20 runs, rather than asserting it
   as a universal guarantee of the `copykat()` function.
6. **"All scientific parameters at/matching documented defaults"** — this
   summary sentence (used in both §1.2's bottom line and §5.3's
   ruled-out-explanations list) was corrected to be precise that
   InferCNV's `cutoff`, `denoise`, and `analysis_mode` are *not* at the
   coded default — they match a *documented example/recommendation*
   instead, which is a different and weaker claim.
7. **CopyKAT bin size / smoothing-vs-segmentation wording** — the "~220 kb
   bins" figure is sourced from `copykat.Rd`'s own `@return` field (now
   cited explicitly); "smoothing and segmentation are the same step" was
   softened to note this audit did not independently verify from source
   whether CopyKAT has a fully separate smoothing pass.
8. **Tumor_01 causal attribution** — reworded from stating the
   baseline-difference hypothesis as established fact to explicitly
   framing it as a plausible, code-supported hypothesis consistent with
   the evidence, not a proven mechanism (the document already correctly
   scoped this as unproven in its closing paragraph; the earlier
   statement of the same idea was inconsistently more assertive and has
   been aligned with it).

No factual error was found in the parameter *values* transcribed from the
two pipeline scripts, in the 54.08%/56.4% arithmetic, or in the
29,175-vs-29,173 cell accounting. All issues were about precision of
wording/scope, not incorrect numbers or incorrect script transcription.
