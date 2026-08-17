# Poster Story v1 -- Planning Memo

This phase throws away the "figure bank" mindset entirely and rebuilds the
visual story around ONE central heatmap, with a small number (6-8) of
genuinely strong supporting figures. No new science; every number below
traces to an already-frozen table. See `DATA_AUDIT.md` (this folder) for
the full source inventory this plan is built on.

## The story

**"From a CRISPR screen to therapeutic vulnerabilities"** told as six
beats, each with exactly one figure carrying it:

1. **Functional discovery** -- which genes emerge from the CRISPR screen,
   and where do the 4 focus genes sit among them.
2. **Cross-context corroboration** (THE HERO) -- one heatmap showing how
   the 4 focus genes behave across every transcriptomic context this
   project has: two resistance-model cell-line panels, one human
   recurrence panel, one acute human ex-vivo panel.
3. **Biological interpretation** -- pathway convergence (estrogen
   response, EMT) across the same contexts.
4. **Translational interpretation** -- structural/pharmacological
   maturity, compared honestly across all 4 genes (including VEZF1's real
   absence of a structure).
5. **Disease / clinical context** -- real human recurrence signal
   (GSE240112) and real malignant-vs-non-malignant specificity
   (GSE245601, copyKAT-derived malignancy calls), plus baseline DepMap
   dependency as an orthogonal axis.
6. **Final synthesis** -- why these 4 genes, and what distinct role each
   plays, without a scorecard.

## What the central heatmap must encode

Four real transcriptomic contexts exist in this project, each already
structured as a natural two-condition PAIR:

| Dataset | Pair | What it is |
|---|---|---|
| GSE118713 | MCF7 (baseline) vs TAMR (resistant) | resistance-model cell line |
| GSE111151 | Parental (4 backgrounds) vs Resistant (7 sublines) | resistance-model cell line |
| GSE240112 | Primary (n=3) vs Recurrent (n=3) | human recurrence, UNPAIRED |
| GSE245601 | Control vs Tamoxifen (12h) | acute human ex-vivo, patient-matched |

**The comparability problem, stated plainly:** GSE118713 is TPM, GSE111151
and GSE240112 are log2(CPM), and GSE245601's per-gene values are computed
here via a disclosed log2(CPM+1) transform (see `DATA_AUDIT.md`). These
are NOT on a shared absolute scale, and a single heatmap colorbar spanning
raw values from all four would misrepresent that as comparability. log2
fold-change (a ratio), on the other hand, IS a standard, defensible,
comparable unit across studies -- this is why meta-analysis figures in
real papers routinely show cross-study log2FC heatmaps.

**Resolution chosen after building and comparing 5 candidates (Phase 2):**
a **paired-row hybrid** design. Each dataset contributes two visually
linked rows (bracket-connected, matching the requested design): the
baseline/reference condition row is always neutral (it is the reference
point, not "zero signal" -- labeled as such), and the second row is
colored by that gene's real log2 fold-change relative to its own paired
baseline in that same dataset. This keeps the requested "two paired rows"
visual structure, keeps every color honest (log2FC, a comparable unit),
and never blends TPM against log2CPM on one colorbar. Full comparison of
all 5 candidates and the final decision are in `FINAL_FIGURE_RECOMMENDATION.md`.

## Gene-set decision

**Hero heatmap = 4 focus genes** (KDM1A, TLK2, USP34, VEZF1). Reasoning:
this project's central claim after the post-audit sensitivity analysis is
that these 4 occupy *distinct* roles, not that they are the top-13
CRISPR hits -- a 4-column heatmap is legible at poster distance and tells
that exact story. **Supporting alternate = all 13 significant sensitising
hits**, using the SAME already-frozen log2FC columns
(`gse118713_log2fc`, `gse111151_log2fc`, `gse240112_tumor_log2fc`,
`gse245601_epi_log2fc` in `all_genes_cross_dataset_evidence_with_ranking.tsv`)
-- this shows the 4 focus genes are not cherry-picked in isolation, but
this version is denser and is kept as a supporting/appendix panel, not
the hero.

## New real data layer found this phase: malignant vs. non-malignant

`results/tables/gse245601_pseudobulk/malignant_vs_nonmalignant_*` and
`results/tables/gse245601_candidate_integration/malignant_vs_nonmalignant_candidates.tsv`
contain real, already-frozen pseudobulk comparisons (5 patients:
Tumor_02/03/07/09/10) of malignant vs. non-malignant epithelial cells,
covering all 13 significant sensitising genes including all 4 focus
genes. The malignancy calls themselves come from copyKAT
(`gse245601_copykat_sensitivity_labels.tsv`, cell-level calls feeding the
`cnv_score` used in `malignant_vs_nonmalignant_cell_level_summary.tsv`).
This is real, frozen, and was never visualized in any prior poster phase.

**Decision:** use it as ONE panel inside the disease/clinical-context
figure (a real specificity check: is the signal tumor-cell-specific?),
not as a separate "copyKAT figure." A standalone copyKAT QC figure would
be methods-supplement material, not poster-grade -- exactly the kind of
figure this phase is explicitly told to avoid. The copyKAT provenance is
stated once in that figure's caption, not turned into its own panel.

## What gets dropped entirely

- **Network figures as a main-sequence panel.** Confirmed (again) this
  phase: KDM1A and TLK2 have ZERO frozen network analysis anywhere in this
  project (the systems-network phase only ever ran on the original
  4-candidate set). A network figure cannot carry the 4-focus-gene story.
  It is built once, honestly, as an explicit backup/appendix item only --
  never promoted to the main shortlist.
- **Every dashboard/scorecard/role-card design.** No colored boxes
  asserting a conclusion in text anywhere in this phase.
- **Small-multiple grids larger than ~4 panels.** If a comparison needs
  more than 4 small panels, it is reworked into one denser, more elegant
  encoding (as the hero heatmap does for what would otherwise be 4
  datasets x 4 genes x 2 conditions = 32 numbers).
- **A separate TCGA panel as a hero.** Already established (science-freeze
  phase) that TCGA associations for USP34/VEZF1 are weak/non-significant
  (FDR 0.21, 0.90) and TCGA was never run for KDM1A/TLK2. Not worth a
  panel at all in a small, high-signal set; omitted entirely rather than
  demoted to "secondary," since this set has no room for a weak panel.
- **Re-designing figures that already met the bar.** The prior phase
  (poster_exploration_v3) already produced pathway-convergence and
  structural-comparison figures independently judged strong,
  poster-grade heroes. Re-designing them from scratch here would waste
  effort and risk making them worse for no reason. This phase reuses
  their code directly (imported, not copy-pasted) into the new namespace,
  and concentrates NEW design effort on the two things that did not exist
  before: the hero heatmap, and the disease/clinical-context figure using
  the newly-surfaced malignant-vs-non-malignant layer.

## Planned figure set (6-8 main + explicit backup)

1. CRISPR discovery (reused design, new namespace)
2. **Hero: cross-context heatmap** (NEW, 4-gene paired-row hybrid)
3. Pathway convergence (reused design)
4. Structural / pharmacological comparison (reused design)
5. Disease / clinical context: recurrence + malignant-vs-non-malignant (NEW)
6. Baseline DepMap dependency (reused design)
7. Final synthesis (reused design, lightly adapted)

Backup/appendix only: 13-hit broader heatmap variant; USP34 network
neighborhood (real, honest, but explicitly not part of the main sequence).
