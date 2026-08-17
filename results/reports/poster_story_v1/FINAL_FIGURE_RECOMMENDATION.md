# Final Figure Recommendation -- Poster Story v1

## 1. What is the best central heatmap?

**Candidate 4: the paired-row hybrid design** (`HERO_cross_context_heatmap.png`,
4 focus genes). All 5 candidates were built and compared side by side in
`CONTACT_HEATMAP_CANDIDATES.png`:

| # | Design | Verdict |
|---|---|---|
| 1 | Raw paired rows, within-block min-max color | Honest but color is not cross-dataset comparable, and that caveat has to live in the caption every time -- a real weakness for a poster where captions get ignored. |
| 2 | Delta-only, single row per dataset | Scientifically cleanest, but drops the explicit visual "pairing" the story needs -- reads more like a generic heatmap than a before/after story. |
| 3 | Paired rows + delta arrows, within-block color | The arrows add real information but visually compete with the cell text; busier than Candidate 4 for the same information content. |
| **4** | **Paired rows, reference row neutral gray, comparison row colored by log2FC** | **Selected.** Keeps the requested visual pairing (bracket-connected, dataset blocks clearly separated), uses an honestly comparable unit (log2FC) for every color, and is the cleanest, most poster-legible of the five. |
| 5 | Same design as 4, all 13 significant sensitising hits | Real and valuable, but denser -- kept as the supporting/appendix variant, not the hero (see gene-set decision below). |

**Why Candidate 4 wins:** it is the only design that satisfies all three
requirements simultaneously -- visually paired rows with brackets, a
single honestly-comparable color unit across all four datasets, and
poster-legible density (4 columns, 8 rows, generous whitespace).
Candidates 1 and 3 are honest but weaker because their color scale
genuinely isn't comparable across blocks; Candidate 2 is honest and
clean but loses the requested paired-row story structure.

## 2. Gene-set decision (re-evaluated, not assumed)

**Hero = 4 focus genes. Supporting/appendix = 13-hit broader version.**
Confirmed correct after building both: the 4-gene version is
poster-legible at a glance and matches the "distinct roles" narrative
that is this project's actual scientific conclusion; the 13-gene version
(Candidate 5) is real and valuable evidence that the 4 focus genes are
not cherry-picked in isolation, but is too dense to be the primary visual
-- kept as `CANDIDATE5_hybrid_13hit.png` for appendix/backup use only.

## 3. Best final figure sequence (7 figures, heatmap as central anchor)

1. **`FIG1_crispr_discovery`** -- genome-wide CRISPR lollipop, establishes the discovery funnel.
2. **`HERO_cross_context_heatmap`** -- THE central anchor. Answers "how do these 4 genes behave everywhere we looked."
3. **`FIG2_structural_comparison`** -- translational interpretation: structural/pharmacological maturity, real structures for 3 genes, honest absence for VEZF1.
4. **`FIG3_pathway_convergence`** -- biological interpretation: estrogen response down everywhere, EMT divergent by context.
5. **`FIG4_disease_clinical_context`** -- NEW this phase: real human recurrence signal + real malignant-vs-non-malignant tumour-cell specificity (copyKAT-derived).
6. **`FIG5_depmap_context`** -- baseline dependency, real per-line heatmap.
7. **`FIG6_final_synthesis`** -- why these genes, framed as an interpretation framework, not a leaderboard.

See `CONTACT_MAIN_STORY.png` for all 7 in this order on one sheet.

## 4. Which previous figure concepts should be retired completely?

- **Every small-multiple grid larger than ~4 panels** (v2's B1/B3/B4/B5a/b,
  C1b -- the ~16-panel-equivalent transcriptomic dumps). The hero heatmap
  now carries that entire story (4 datasets x 4 genes x 2 conditions = 32
  numbers) in one elegant figure instead.
- **v2's UpSet plot and quantitative evidence map (E1/E2).** Real and
  informative, but report-grade at their original density; the hero
  heatmap plus `FIG6_final_synthesis` now carry the "why these genes"
  message more elegantly.
- **v3's separate B1-B4 four-panel transcriptomic figure.** Fully
  subsumed by the hero heatmap, which shows the same 4 datasets with a
  cleaner, honestly-comparable encoding.
- **Any standalone copyKAT/QC figure.** Deliberately never built (see
  `STORY_PLAN.md`) -- copyKAT's role is disclosed in `FIG4`'s caption,
  not turned into its own methods-supplement panel.
- **TCGA as any kind of panel.** Omitted entirely (weak/non-significant,
  never assessed for KDM1A/TLK2) -- not worth space in a small, high-signal set.

## 5. Which visual hero(s) are strongest?

**The hero heatmap (`HERO_cross_context_heatmap`) and
`FIG2_structural_comparison`** are the two strongest, most memorable
images in the set -- the heatmap because it is genuinely novel (nothing
like it existed in v2/v3) and carries the entire cross-context story in
one legible image; the structural comparison because it is the most
visually striking single image (real ray-traced molecular structures)
and was already independently judged the strongest hero in the prior
phase. `FIG3_pathway_convergence` is a close third.

## 6. Main panel vs. backup panel vs. omitted

**Main panel (7, in the sequence above):** FIG1, HERO, FIG2, FIG3, FIG4, FIG5, FIG6.

**Backup / appendix only (2):**
- `CANDIDATE5_hybrid_13hit.png` -- broader 13-hit heatmap, real evidence the 4 focus genes aren't cherry-picked, but too dense for the main sequence.
- `BACKUP_network_usp34.png` -- USP34's real STRING neighborhood; explicitly not main-sequence because KDM1A and TLK2 have zero frozen network data (a real project gap, not a design choice).

**Omitted entirely:** TCGA (weak/never assessed for 2 of 4 genes), any
standalone copyKAT/QC panel, any scorecard/role-card/dashboard design,
any small-multiple grid larger than 4 panels.
