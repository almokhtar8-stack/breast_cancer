# Exploration-v2 Summary: Shortlist and Proposed Poster Architecture

This is a synthesis document, not a poster assembly. No poster has been
built. It records the shortlist ranking and a proposed architecture for
later review; see `FIGURE_CANDIDATE_GUIDE.md` for full per-figure detail
and `DATA_FOR_VISUALIZATION_AUDIT.md` for data provenance.

## Ranked shortlist (all 28 candidates)

**MUST CONSIDER (12) -- the strongest single candidate in each niche:**

| Figure | Why it's in the top tier |
|---|---|
| A1. Genome-wide ranked landscape | Establishes the discovery scale; inset resolves the 4 focus genes' true rank honestly |
| A2. 13-hit lollipop | Cleanest single demonstration that KDM1A/TLK2 are functionally stronger than USP34/VEZF1 |
| B1. GSE118713 sample panel | Real per-sample data, all 4 genes, immediately shows heterogeneity |
| B2. GSE111151 trajectories | The strongest "honest null" figure in the bank -- real recorded pairing, not invented |
| B3. GSE240112 tumours | Correctly unpaired, VEZF1's real evidence made visible |
| B5a. Dataset-centric overview | One image answers "why did each candidate get retained" across all 4 datasets |
| C2. Estrogen/EMT hero | The single most memorable pathway result (estrogen down everywhere; EMT flips at acute) |
| E1. UpSet evidence intersection | Shows the real evidence-combination structure without a composite score |
| E2. Quantitative evidence map | The richest single evidence-integration panel; 4 real encoded dimensions |
| E3. Historical vs. post-audit | Visualizes the audit's central finding directly (KDM1A/TLK2 excluded by the RNA gate) |
| F1. DepMap heatmap | Real per-line values, real cell-line names -- no aggregation |
| G1. Three structures + VEZF1 | The most visually striking figure in the entire bank; fixes the prior phase's USP34-only narrowness |

**STRONG (9) -- excellent, but redundant with (or an alternative to) a
MUST-CONSIDER figure in the same niche:**
A4 (hexbin density, alt. to A1), B4 (GSE245601 acute, companion to B1-B3),
B5b (gene-centric, alt. orientation to B5a), C3 (real GSEA curves, technical-depth
companion to C2), F2 (cell-line fingerprint, alt. to F1), F3 (human+DepMap
combo, consolidator), G2 (pharmacology framing, alt. to G1), G3 (pocket
close-ups, alt. or inset to G1), H3 (translational maturity ladder, a
genuinely different data-grounded axis pairing).

**MAYBE (4):** C1b (pathway small multiples, only if C1a/C2 prove too
dense), D1 and D2 (network -- real and honest, but KDM1A/TLK2 coverage gap
limits fit for a KDM1A/TLK2-forward poster), H2 (experimental schematic --
useful only if the poster wants an explicit forward-looking panel).

**SUPPLEMENTARY (2):** A3 (redesigned volcano, a clean but redundant
benchmark), F4 (TCGA, deliberately minor given weak/non-significant results).

**REJECT (0):** every built figure earned at least a supporting role; none
were decorative-without-information (see Codex review for an independent
check of this claim).

## What changed from the previous ("poster_final") attempt

The previous 6-figure set was rejected for: wrong plot choices, too much
audit/report logic baked into panels, too many scorecards/summary boxes,
captions doing work the data should do, generic heatmaps where richer
biological plots were possible, and an overall "dashboard" feel. This
exploration bank addresses each directly:

- **No scorecards or role cards anywhere.** Every figure plots real
  sample/patient/cell-line/gene-level values (Sections A, B, F, G) or a
  small number of explicitly-named, non-composite evidence dimensions
  (Sections E, H3) -- never a colored box asserting a conclusion in text.
- **Real biological units are the default, not the exception.** GSE118713
  (9 samples), GSE111151 (11 samples, real pairing), GSE240112 (6 tumours,
  correctly unpaired), GSE245601 (3 patients, correctly paired), DepMap
  (11 real-named cell lines) all appear with individual points, not
  pre-aggregated percentages.
- **Structural comparison is now genuinely multi-candidate.** KDM1A
  (6NQU) and TLK2 (5O0Y) were fetched and rendered this phase, matching
  USP34's existing frozen renders in style -- the prior phase's
  USP34-only structural narrative is gone.
- **A real GSEA implementation, not just a summary table.** C3 reconstructs
  actual running-enrichment-score curves from frozen inputs.
- **Minimal captions inside figures**; detailed caveats live in this guide
  and the audit doc, not as paragraphs baked into each PNG.

## Proposed poster architecture (data-driven, not assumed to be six panels)

Based on what actually looked strongest after building and inspecting all
28 candidates, a **7-panel architecture** is proposed (more panels than the
previous six, because the "show the data" mandate genuinely needs more
visual real estate than a summary-driven poster does):

1. **CRISPR hit landscape** -- A1 (hero) with A2 as an adjacent or
   alternate panel if space allows a second CRISPR view.
2. **Real sample/patient-level transcriptomic evidence** -- B5a as the
   main integrated panel; B2 (GSE111151) promoted alongside it if space
   allows, since it is this bank's strongest "honest null" story.
3. **Pathway dynamics** -- C2 (Estrogen/EMT hero) as the main panel; C3
   (real enrichment curves) as a technical-depth companion if the poster
   has room for it.
4. **Evidence-integration** -- E2 (quantitative evidence map) as the
   primary panel; E1 (UpSet) or E3 (historical-vs-post-audit) as a second
   evidence panel -- **E3 is the strongest single "how did the audit
   change the story" visual in the entire bank and deserves serious
   consideration for main-sequence inclusion**, not just a supporting role.
5. **DepMap / baseline dependency** -- F1 (heatmap) as the main panel.
6. **Structural comparison** -- G1 as the hero image (this bank's single
   most striking figure); G3 as an optional close-up inset.
7. **Translational closing** -- open question, see below.

**Network (Section D):** small or omitted. Both D1 and D2 are real and
honest, but neither can include KDM1A or TLK2 (no frozen network data
exists for either), which conflicts with a poster built around the 4
current focus genes. Recommend appendix/backup material only, not main
sequence, unless a future phase runs the network analysis on KDM1A/TLK2.

**TCGA (Section F, F4):** small or omitted, exactly as F4 was designed --
both USP34 and VEZF1 associations are weak/non-significant (FDR 0.21 and
0.90). Include only as a small secondary panel if a reviewer specifically
asks about human tumour expression beyond GSE240112, never as a main panel.

**GDSC:** remains supplementary, consistent with its existing
SUPPORTING/SUPPLEMENTARY status from the prior poster-final phase's
`FIGURE_BANK_REVIEW.md` -- nothing in this exploration phase changes that;
no new GDSC figure was built.

**Translational closing panel -- open question (H1's finding):** testing
whether the bank already tells the whole story without a forced synthesis
figure suggests two viable options, left open for the poster owner's
choice rather than decided here:
- *(a)* End on G1 (the structural hero image) plus a short written
  conclusion, since E3 and E2 already carry the "no universal winner"
  message quantitatively earlier in the sequence; or
- *(b)* Include H3 (translational maturity ladder) as a genuine 7th data-
  grounded panel if the poster wants an explicit closing data point rather
  than ending on structure.

This is a genuine change from assuming six panels: the strongest material
this phase produced clusters more heavily in evidence-integration (Section
E, 3 strong candidates) and transcriptomics (Section B, multiple strong
real-data candidates) than a six-panel template can comfortably hold
without returning to the density problem the reset was meant to fix.
