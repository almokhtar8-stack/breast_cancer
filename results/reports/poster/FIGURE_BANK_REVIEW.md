# Figure Bank Review

Ratings and CORE / SUPPORTING / SUPPLEMENTARY recommendations for every
candidate figure in `results/figures/poster_candidates/`. Nothing is
deleted based on this ranking -- all 15 candidates remain in the bank for
visual selection. See `SCIENTIFIC_FIGURE_AUDIT.md` for the underlying data
audit that produced this candidate set.

Ratings scale: HIGH / MEDIUM / LOW.

---

## 1a. `01a_crispr_ranked_effect`
- Scientific importance: **HIGH** -- shows the 4 candidates' genome-wide rank among 19,103 fitted genes.
- Visual impact: **HIGH** -- the elbow shape and tight cluster of 4 colored points among the sensitising tail is immediately legible.
- Uniqueness of message: **MEDIUM** -- overlaps with 1c (same underlying data, different axis choice).
- Poster readability: **HIGH** -- large fonts, minimal background clutter (rasterized).
- Redundancy: **MEDIUM** (vs 1c).
- **Recommendation: CORE** (pick this OR 1c, not both -- see note below).

## 1c. `01c_crispr_volcano_refined`
- Scientific importance: **HIGH** -- same discovery claim as 1a, classic volcano form.
- Visual impact: **HIGH** -- USP34/VEZF1 visually dominant as requested, EML5/CITED2 secondary.
- Uniqueness of message: **MEDIUM** -- same underlying data as 1a.
- Poster readability: **HIGH**.
- Redundancy: **MEDIUM** (vs 1a).
- **Recommendation: CORE** (alternative to 1a -- choose whichever reads better at poster scale; 1a communicates *how extreme* the effect is among all genes, 1c communicates *significance* more conventionally).

## 03. `03_GSE118713_resistance_landscape`
- Scientific importance: **HIGH** -- proves the resistant phenotype is genome-wide, not single-gene, and situates USP34 inside it.
- Visual impact: **HIGH** -- 3 real panels (PCA + volcano + inset), no icons.
- Uniqueness of message: **HIGH** -- PCA panel does not exist anywhere else in the project's poster figures.
- Poster readability: **MEDIUM-HIGH** -- 3 panels is a lot of information in one figure; the PCA and volcano panels can each stand alone if space is tight.
- Redundancy: **LOW**.
- **Recommendation: CORE.**

## 04. `04_GSE240112_recurrence`
- Scientific importance: **HIGH** -- the strongest human-tumor evidence for VEZF1, now shown against the full genome-wide recurrence signal.
- Visual impact: **HIGH**.
- Uniqueness of message: **HIGH**.
- Poster readability: **HIGH**.
- Redundancy: **LOW** (panel B overlaps with the old poster figure 2 panel B, but panel A is new).
- **Recommendation: CORE.**

## 05. `05_cross_dataset_candidate_effects`
- Scientific importance: **HIGH** -- the only figure that shows all 4 candidates across all 4 datasets at once, with the acute-context caveat visually enforced (square markers, red italic labels).
- Visual impact: **MEDIUM-HIGH** -- information-dense; reads well at close range, needs a large print size at poster scale.
- Uniqueness of message: **HIGH**.
- Poster readability: **MEDIUM** -- 16 rows of small text labels is the densest figure in the bank.
- Redundancy: **LOW**.
- **Recommendation: CORE.**

## 06. `06_resistance_pathway_landscape`
- Scientific importance: **HIGH** -- this is the single largest completed analysis (5,847 pathways scored) that had **no** figure anywhere in the project before this phase.
- Visual impact: **HIGH** -- a real heatmap with real NES values, FDR flagged, CRISPR context correctly separated.
- Uniqueness of message: **HIGH** -- pathway/systems-biology message exists nowhere else in the bank.
- Poster readability: **HIGH**.
- Redundancy: **NONE**.
- **Recommendation: CORE.** Strongest pathway-biology candidate in the bank.

## 07. `07_candidate_mechanism_map`
- Scientific importance: **MEDIUM** -- correctly shows sparse, mostly-absent candidate-candidate connectivity (a real, if modest, finding), and EML5's total lack of pathway footprint.
- Visual impact: **MEDIUM** -- clean and honest, but visually quiet compared to the data-dense figures.
- Uniqueness of message: **MEDIUM** -- content overlaps with figure 06 (same STRONG_CONSENSUS pathway data, different framing).
- Poster readability: **HIGH**.
- Redundancy: **MEDIUM** (vs 06).
- **Recommendation: SUPPORTING.** Useful for explaining *why* EML5 was deprioritized and how little the 4 candidates converge, but not essential if poster space is limited and figure 06 is already included.

## 08. `08_TCGA_human_validation`
- Scientific importance: **HIGH** -- TCGA was almost entirely absent from the previous poster set despite being a major completed analysis (n=1095 tumors).
- Visual impact: **HIGH** -- real forest plot with real 95% CIs, honest about which candidates are/aren't significant.
- Uniqueness of message: **HIGH**.
- Poster readability: **HIGH**.
- Redundancy: **NONE**.
- **Recommendation: CORE.** Strongest human-tumor-validation candidate in the bank.

## 09. `09_Hany_vs_DepMap_context_map`
- Scientific importance: **HIGH** -- directly answers "USP34 vs VEZF1: what's actually different?" with a real two-axis map across all 28 Gate-1 genes, not just the 4 candidates.
- Visual impact: **HIGH** -- immediately shows USP34 and VEZF1 occupying different quadrants for a genuine biological reason.
- Uniqueness of message: **HIGH**.
- Poster readability: **HIGH**.
- Redundancy: **LOW** (complements, doesn't duplicate, 09b).
- **Recommendation: CORE.** Likely the single best "USP34 vs VEZF1" figure in the entire bank.

## 09b. `09b_USP34_VEZF1_line_dependencies`
- Scientific importance: **MEDIUM-HIGH** -- shows the actual 11 per-line values and the 2 specific lines driving VEZF1's 27% figure.
- Visual impact: **HIGH** -- the slope-graph format is unusual and reads well.
- Uniqueness of message: **MEDIUM** -- narrower version of the same claim as figure 09 and the frozen `poster/03_depmap_distributions.png`.
- Poster readability: **HIGH**.
- Redundancy: **MEDIUM** (vs 09 and vs the existing `poster/03`).
- **Recommendation: SUPPORTING.** Good supplementary evidence for a "drill-down" panel if 09 is the headline DepMap figure.

## 10. `10_tissue_liability_context`
- Scientific importance: **MEDIUM** -- genuinely separates expression from documented functional liability (avoids a fabricated "toxicity score"), but the finding itself ("broad, non-specific expression for both genes") is a modest, mostly-negative result.
- Visual impact: **MEDIUM**.
- Uniqueness of message: **MEDIUM**.
- Poster readability: **MEDIUM** -- panel B's row count (10 organ systems) is dense at small size.
- Redundancy: **LOW**.
- **Recommendation: SUPPLEMENTARY.** Valuable as safety-diligence backup material; not a headline scientific claim.

## 11. `11_GDSC_USP34_pharmacogenomics`
- Scientific importance: **MEDIUM** -- association-only pharmacogenomics on a secondary axis, explicitly labeled as such.
- Visual impact: **HIGH** -- the lollipop + scatter combination is clean and the AZD7762 relationship is visually convincing.
- Uniqueness of message: **MEDIUM** -- overlaps with the frozen `poster/04_pharmacogenomics.png`.
- Poster readability: **HIGH**.
- Redundancy: **MEDIUM** (vs the existing frozen poster figure 4, which shows the same AZD7762 scatter as panel B here).
- **Recommendation: SUPPORTING.** Panel A (top associations) is new and adds value beyond the existing frozen figure; panel B duplicates it.

## 12. `12_USP34_structure_surface`
- Scientific importance: **HIGH** -- same real structural facts as before (catalytic dyad, covalent probe), now shown as an actual publication-style hero image.
- Visual impact: **HIGH** -- the surface render with the highlighted cleft patch is the most visually striking single image in the bank.
- Uniqueness of message: **HIGH** (vs the old 4-tile poster figure).
- Poster readability: **HIGH** -- large structure fills most of the frame.
- Redundancy: **MEDIUM** (vs 12b -- same two structures, different composition).
- **Recommendation: CORE.** Strongest structural candidate for a poster that wants one large, striking structural image.

## 12b. `12b_USP34_structure_comparison`
- Scientific importance: **HIGH** -- same facts as 12, framed as the apo-vs-bound comparison the original request specifically asked for.
- Visual impact: **HIGH** -- three large, matched-orientation panels; the ubiquitin probe (green) is unambiguous in panel B.
- Uniqueness of message: **MEDIUM** (vs 12).
- Poster readability: **HIGH**.
- Redundancy: **MEDIUM** (vs 12).
- **Recommendation: CORE (alternative to 12).** Pick 12 for a single striking hero image, or 12b if the apo-vs-bound structural comparison is the more important story to tell; both are real, publication-quality PyMOL renders, not "screenshots."

## 13. `13_validation_experiment`
- Scientific importance: **MEDIUM** -- this is the one schematic in the bank by design, not primary evidence.
- Visual impact: **MEDIUM-HIGH** -- the cell-glyph style (nucleus, KO mark, treatment ring) reads as biological rather than a generic flowchart.
- Uniqueness of message: **HIGH** within the bank (nothing else shows the proposed forward experiment).
- Poster readability: **HIGH**.
- Redundancy: **NONE**.
- **Recommendation: SUPPORTING.** Necessary to close the poster's narrative arc but intentionally not a main evidence panel.

---

## Summary table

| Figure | Importance | Visual impact | Uniqueness | Readability | Redundancy | Verdict |
|---|---|---|---|---|---|---|
| 1a | HIGH | HIGH | MEDIUM | HIGH | MEDIUM | CORE |
| 1c | HIGH | HIGH | MEDIUM | HIGH | MEDIUM | CORE |
| 03 | HIGH | HIGH | HIGH | MED-HIGH | LOW | CORE |
| 04 | HIGH | HIGH | HIGH | HIGH | LOW | CORE |
| 05 | HIGH | MED-HIGH | HIGH | MEDIUM | LOW | CORE |
| 06 | HIGH | HIGH | HIGH | HIGH | NONE | CORE |
| 07 | MEDIUM | MEDIUM | MEDIUM | HIGH | MEDIUM | SUPPORTING |
| 08 | HIGH | HIGH | HIGH | HIGH | NONE | CORE |
| 09 | HIGH | HIGH | HIGH | HIGH | LOW | CORE |
| 09b | MED-HIGH | HIGH | MEDIUM | HIGH | MEDIUM | SUPPORTING |
| 10 | MEDIUM | MEDIUM | MEDIUM | MEDIUM | LOW | SUPPLEMENTARY |
| 11 | MEDIUM | HIGH | MEDIUM | HIGH | MEDIUM | SUPPORTING |
| 12 | HIGH | HIGH | HIGH | HIGH | MEDIUM | CORE |
| 12b | HIGH | HIGH | MEDIUM | HIGH | MEDIUM | CORE |
| 13 | MEDIUM | MED-HIGH | HIGH | HIGH | NONE | SUPPORTING |

**1a vs 1c and 12 vs 12b are each an either/or choice, not two separate poster panels** -- both members of each pair tell the same core story from a different angle; use whichever reads better once laid out at actual poster scale.

## Which previous poster figures this bank suggests retiring

- `poster/03_depmap_distributions.png` (single-percentage-driven 4-gene box) is superseded by the combination of `09_Hany_vs_DepMap_context_map` (adds the CRISPR axis and all 28 Gate-1 genes) and `09b_USP34_VEZF1_line_dependencies` (keeps the per-line drill-down, in a cleaner slope-graph form).
- `poster/04_pharmacogenomics.png` (AZD7762 scatter alone) is superseded by `11_GDSC_USP34_pharmacogenomics` panel B (same scatter) plus panel A (new top-associations context).
- `poster/01_crispr_discovery.png` (volcano) is superseded by whichever of `01a`/`01c` is chosen -- both are stripped-down, more publication-styled versions of the same message.
- `poster/02_expression_evidence.png` panels A and B are superseded by `03_GSE118713_resistance_landscape` and `04_GSE240112_recurrence`, which add genome-wide context (PCA, full volcano) the old panels didn't have; panel C (the forest plot) is superseded by the larger, cleaner `05_cross_dataset_candidate_effects`.
- `poster/05_structure.png` (4-equal-tile PyMOL composite) is superseded by `12` or `12b`.
- `poster/06_experimental_strategy.png` is superseded by `13_validation_experiment` (biological-glyph redesign).

None of the frozen `poster/` files have been modified or deleted -- this is a recommendation for the next (not-yet-started) poster-assembly phase, to be acted on only after visual review.
