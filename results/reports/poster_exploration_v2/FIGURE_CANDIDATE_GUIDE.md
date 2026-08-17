# Figure Candidate Guide (poster-exploration-v2)

28 exploratory figures across 8 sections (A-H). Every figure's source data
is documented in `DATA_FOR_VISUALIZATION_AUDIT.md`; every number is read
from a frozen table via `src/poster_exploration_v2_data.py` (see that
module's docstrings for exact sources). This guide records, per figure:
question answered, source, biological unit, encoding, interpretation,
caveat, strengths/weaknesses, likely poster role, size, and a keep
recommendation. Detailed shortlist ranking and rationale for the ~8-12
strongest candidates is in `EXPLORATION_V2_SUMMARY.md`.

---

## Section A -- CRISPR discovery

### A1. `A1_genomewide_ranked_landscape`
- Question: Where do the 4 focus genes sit among all 19,103 fitted genes?
- Source: `data/processed/labels.parquet`, all genes, ranked by effect.
- Unit: gene (n=19,103).
- Encoding: x=rank, y=effect size; gray background = non-Gate-1; darker
  gray = Gate-1; colored dots = focus genes; zoomed inset = rank 1-80.
- Interpretation: focus genes sit deep in the sensitising tail but are
  visually indistinguishable at genome scale without the inset.
- Caveat: inset x-range differs from main panel (necessarily, given scale).
- Strengths: shows the full discovery funnel honestly; inset resolves the
  "how deep in the tail" question precisely.
- Weaknesses: inset is small; requires close reading.
- Poster role: strong opener (visually establishes "genome-wide screen").
- Size: wide, landscape.
- Keep: **yes, MUST CONSIDER.**

### A2. `A2_lollipop_13_hits`
- Question: How do the 13 significant sensitising hits compare functionally?
- Source: `post_audit_sensitivity_data.load_significant_sensitising_hits()`.
- Unit: gene (n=13).
- Encoding: y=gene (sorted by effect), x=effect size, point size=FDR strength.
- Interpretation: KDM1A/TLK2 visually dominate; USP34/VEZF1 are mid-to-low.
- Caveat: point-size FDR encoding is relative, not an absolute scale bar.
- Strengths: clean, immediately legible, real Cleveland-plot convention.
- Weaknesses: none significant.
- Poster role: strong "why not just pick the top CRISPR hit" panel.
- Size: portrait, moderate.
- Keep: **yes, MUST CONSIDER.**

### A3. `A3_redesigned_volcano`
- Question: benchmark/comparator volcano (minimal, polished).
- Source: same as A1.
- Unit: gene (n=19,103).
- Encoding: classic volcano, x=effect, y=-log10(FDR).
- Interpretation: standard, familiar framing; less information-dense than A1.
- Caveat: none.
- Strengths: instantly recognizable format for a genomics audience.
- Weaknesses: redundant with A1's information content; less novel.
- Poster role: supplementary/backup only if A1 doesn't work for a reviewer.
- Size: square, moderate.
- Keep: **SUPPLEMENTARY** (redundant with stronger A1).

### A4. `A4_hexbin_density`
- Question: how far into the density tail do focus genes sit?
- Source: same as A1.
- Unit: gene (n=19,103), hexbin-aggregated.
- Encoding: hexbin density (grayscale) + focus-gene + RCOR1 overlay.
- Interpretation: focus genes and RCOR1 both sit in extremely sparse tail bins.
- Caveat: hexbin obscures exact individual-gene positions (a density view,
  not a scatter).
- Strengths: visually striking, includes the RCOR1 honesty check others don't.
- Weaknesses: less precise than A1's inset for exact rank reading.
- Poster role: strong alternative opener to A1; RCOR1 inclusion is valuable.
- Size: square.
- Keep: **STRONG.**

---

## Section B -- transcriptomics (real sample-level data)

### B1. `B1_gse118713_sample_panel`
- Question: how do all 4 focus genes actually behave across MCF7/TAMR/FASR?
- Source: `data/processed/gse118713_gene_tpm.parquet`.
- Unit: sample (n=9: 3 MCF7 + 3 TAMR + 3 FASR).
- Encoding: 2x2 small multiples, real jittered points + mean bar per gene.
- Interpretation: real heterogeneity visible; KDM1A/TLK2 show little
  separation; USP34 shows a real (if modest) TAMR/FASR shift.
- Caveat: n=3 per condition, one derivation event (stated in footnote).
- Strengths: shows the actual data, not a summary statistic.
- Weaknesses: dense at small size (4 mini-panels).
- Poster role: core transcriptomics panel.
- Size: square/portrait.
- Keep: **yes, MUST CONSIDER.**

### B2. `B2_gse111151_trajectories`
- Question: parental->resistant trajectories, real design.
- Source: `data/processed/gse111151_log2cpm.parquet` + real pairing metadata.
- Unit: sample (n=11: 4 parental + 7 resistant), real slopegraph pairing.
- Encoding: 2x2 small multiples, vertical connector = real parental->
  derivative pairing.
- Interpretation: heterogeneous, mostly-null result across genes -- an
  honest, visible null layer.
- Caveat: BT-474 is HER2-amplified, a different subtype (stated in footnote).
- Strengths: uses REAL recorded pairing, not invented; visually honest null.
- Weaknesses: none significant.
- Poster role: core transcriptomics panel -- the strongest "honest null"
  figure in the whole bank.
- Size: square/portrait.
- Keep: **yes, MUST CONSIDER.**

### B3. `B3_gse240112_tumours`
- Question: real per-tumour recurrence-associated values, all 4 genes.
- Source: `results/tables/gse240112/candidate_sample_level_log2cpm.tsv`.
- Unit: tumour (n=6: 3 PT + 3 RT), UNPAIRED.
- Encoding: 2x2 small multiples, no connecting lines (deliberately).
- Interpretation: VEZF1 (bold) shows real separation; others do not.
- Caveat: unpaired, different biobanks, n=3 vs n=3 (stated in footnote/title).
- Strengths: correctly never implies pairing; VEZF1 emphasis matches the
  frozen conclusion.
- Weaknesses: n=3 per arm is visually thin.
- Poster role: core VEZF1-evidence panel.
- Size: square/portrait.
- Keep: **yes, MUST CONSIDER.**

### B4. `B4_gse245601_acute_paired`
- Question: real patient-matched acute (12h) response, all 4 genes.
- Source: frozen `patient_malignant_pseudobulk.tsv` (USP34/VEZF1) + a
  disclosed formula-identical extension for KDM1A/TLK2 (see audit doc).
- Unit: patient (n=3 eligible), Control->Tamoxifen, genuinely paired.
- Encoding: connected lines per patient (appropriate here, unlike B3).
- Interpretation: heterogeneous per-patient direction; small n visible, not
  hidden.
- Caveat: n=3 patients only; explicitly labeled ACUTE, not resistance.
- Strengths: extends real coverage to all 4 focus genes; correct use of
  connecting lines (contrast with B3's deliberate absence of lines).
- Weaknesses: small n.
- Poster role: strong acute-context companion to B1-B3.
- Size: square/portrait.
- Keep: **STRONG.**

### B5a. `B5a_integrated_dataset_centric`
- Question: integrated overview, dataset-major ordering.
- Source: same as B1-B4, reorganized.
- Unit: mixed (samples/tumours/patients per row).
- Encoding: 4x4 grid, condensed mean+spread strips per cell.
- Interpretation: lets a viewer compare all 4 genes within one dataset at a
  glance.
- Caveat: condensed strips lose some of B1-B4's individual-point richness.
- Strengths: excellent overview; answers "why did USP34/VEZF1/KDM1A/TLK2
  diverge" in one image.
- Weaknesses: dense (16 mini-panels); best at large poster size.
- Poster role: candidate hero transcriptomics panel (may replace B1-B4
  individually if poster space is limited).
- Size: large, landscape.
- Keep: **yes, MUST CONSIDER.**

### B5b. `B5b_integrated_gene_centric`
- Question: same as B5a, gene-major ordering.
- Source/unit/encoding: identical data to B5a, transposed.
- Interpretation: easier to read "this gene's whole story" per row.
- Caveat: same as B5a.
- Strengths: complements B5a; gene-centric reading may suit a poster
  narrated candidate-by-candidate.
- Weaknesses: same density concern as B5a.
- Poster role: alternative to B5a -- pick whichever orientation suits the
  final poster's narrative flow.
- Size: large, landscape.
- Keep: **STRONG** (pick one of B5a/B5b, not necessarily both).

---

## Section C -- pathway biology

### C1a. `C1a_pathway_trajectories_lines`
- Question: how do 7 pathways trend across the 4 transcriptomic contexts?
- Source: `results/tables/systems_network/gsea_{dataset}.tsv`.
- Unit: pathway x dataset (7 pathways x 4 contexts).
- Encoding: one line per pathway, shaded column = acute context.
- Interpretation: estrogen response falls everywhere; EMT/E2F/Wnt rise in
  non-acute contexts.
- Caveat: NES scale is not directly comparable across differently-sized
  gene sets (a general GSEA caveat, not fixed here).
- Strengths: real convergence story visible in one image.
- Weaknesses: 7 lines is near the upper limit of legibility.
- Poster role: strong core pathway panel.
- Size: wide, landscape.
- Keep: **yes, MUST CONSIDER.**

### C1b. `C1b_pathway_small_multiples`
- Question: same as C1a, alternative layout.
- Source/unit: same as C1a.
- Encoding: one bar-chart row per pathway (sparkline style).
- Interpretation: same story, per-pathway rows.
- Caveat: each row's y-axis autoscales independently (not a shared scale).
- Strengths: avoids C1a's line-crossing density.
- Weaknesses: independent y-scales per row could mislead a fast reader
  into comparing magnitudes across rows; needs a clear caption if used.
- Poster role: comparator/alternative to C1a.
- Size: tall, portrait.
- Keep: **MAYBE** (only if C1a proves too dense; needs the y-scale caveat
  stated explicitly in any poster caption).

### C2. `C2_estrogen_emt_hero`
- Question: high-impact, simplified 2-pathway story.
- Source: same as C1a, restricted to Estrogen response + EMT.
- Unit: pathway x dataset.
- Encoding: 2 large panels, line plot + bar chart.
- Interpretation: the two clearest, most memorable pathway results.
- Caveat: none beyond C1a's general GSEA caveat.
- Strengths: the single most "remember this" pathway figure in the bank.
- Weaknesses: omits the 4 secondary pathways (ECM/E2F/G2M/Wnt) shown in C1a/b.
- Poster role: **strong hero-figure candidate**, possibly THE pathway panel
  if poster space allows only one.
- Size: wide, landscape.
- Keep: **yes, MUST CONSIDER.**

### C3. `C3_enrichment_curves`
- Question: what do the real GSEA running-enrichment curves look like?
- Source: real ranked-gene lists + Hallmark gene sets (reconstructed
  running-sum statistic, see audit doc).
- Unit: gene (genome-wide per dataset), reconstructed curve.
- Encoding: classic GSEA running-ES curve + rank-position tick marks.
- Interpretation: visually demonstrates enrichment shape, not just a
  summary NES/FDR number; makes the acute-vs-non-acute EMT contrast vivid.
- Caveat: curve is a visualization of an already-frozen statistic, not a
  new test (stated in footnote).
- Strengths: the only figure in the bank showing "real GSEA machinery" --
  strong for a technically sophisticated audience.
- Weaknesses: less immediately intuitive to a non-genomics audience than
  C1/C2's simpler line/bar plots.
- Poster role: excellent technical-depth panel or poster-adjacent handout.
- Size: square, moderate.
- Keep: **STRONG.**

---

## Section D -- network / systems biology

### D1. `D1_candidate_program_network`
- Question: which biological programs connect to which candidates?
- Source: `candidate_pathway_membership.tsv`, original 4 candidates only.
- Unit: candidate-pathway edge (n=4 USP34 + 2 VEZF1 + 0 EML5 + 49 CITED2,
  aggregated).
- Encoding: bipartite network, real edges only.
- Interpretation: USP34/VEZF1 have real but narrow membership; USP34-VEZF1
  share no direct interaction or pathway.
- Caveat: **KDM1A/TLK2 have no frozen network data and are absent** --
  explicitly disclosed in-figure, not a redesigned 4-focus-gene network.
- Strengths: completely honest about a real project coverage gap.
- Weaknesses: cannot show the two currently-most-important focus genes
  (KDM1A/TLK2) -- limits its poster usefulness given the current framing.
- Poster role: supporting/backup only, given the coverage gap.
- Size: wide, landscape.
- Keep: **MAYBE** (valuable as an honesty artifact, but the missing
  KDM1A/TLK2 coverage limits its fit for a KDM1A/TLK2-forward poster).

### D2. `D2_usp34_local_neighborhood`
- Question: what does USP34's real STRING neighborhood look like?
- Source: `four_candidate_direct_neighbors.tsv`, USP34 subset.
- Unit: gene-gene edge (n=10 real neighbors).
- Encoding: radial network, edge width/opacity = STRING confidence.
- Interpretation: mostly ubiquitin-system genes (UBB/UBC/UBA52/RPS27A);
  none independently reach CRISPR significance.
- Caveat: 1 of 30 total candidate-neighbor rows across the whole table
  (not this subset) does reach FDR<0.1 -- stated exactly, not rounded to
  "none."
- Strengths: visually clean, real data, correctly scoped claim.
- Weaknesses: USP34-only (same network coverage gap as D1).
- Poster role: supporting panel for a USP34-focused sub-story.
- Size: square, moderate.
- Keep: **MAYBE** (same coverage-gap caveat as D1).

**Section D overall: does network deserve poster space?** Marginal. Both
figures are honest and real, but the total absence of KDM1A/TLK2 network
coverage means neither can carry the "4 focus genes" narrative the rest of
the poster tells. Recommend supporting/appendix role only, not main sequence.

---

## Section E -- evidence intersection

### E1. `E1_upset_evidence_intersection`
- Question: which evidence-set combinations do the 13 hits actually occupy?
- Source: `build_evidence_sets_13()` (6 boolean evidence sets).
- Unit: gene (n=13), grouped into 6 observed combinations.
- Encoding: manual UpSet plot (stacked bar + dot-matrix), matplotlib-native.
- Interpretation: KDM1A alone occupies the "validated inhibitor + structure
  + low dependency" combination; no gene combines high dependency with
  strong RNA support.
- Caveat: high/low dependency kept as separate, non-overlapping sets
  (stated explicitly).
- Strengths: real UpSet-style intersection, no upsetplot dependency added;
  correctly avoids conflating advantage/liability sets.
- Weaknesses: 6-set UpSet plots reward close reading; not a glanceable panel.
- Poster role: strong for the "why these genes" evidence-integration story.
- Size: wide, landscape.
- Keep: **yes, MUST CONSIDER.**

### E2. `E2_quantitative_evidence_map`
- Question: quantitative 2D map of CRISPR strength vs. baseline dependency,
  with RNA/structure/pharmacology as marker channels.
- Source: `build_evidence_matrix()` + Table 06b.
- Unit: gene (n=13).
- Encoding: x/y position, marker shape (RNA support), ring color
  (structure), halo (validated inhibitor) -- 4 real dimensions, no score.
- Interpretation: KDM1A dominates on strength+pharmacology but has zero
  dependency and zero non-acute RNA support; TLK2 has strong dependency;
  USP34/VEZF1 are RNA-triangle-marked.
- Caveat: none beyond the general "not a composite" statement in footnote.
- Strengths: the richest single-panel evidence-integration figure in the
  bank; every visual channel is real and named.
- Weaknesses: 4 simultaneous encoding channels require a careful legend
  read.
- Poster role: **strongest evidence-integration candidate overall.**
- Size: wide, landscape.
- Keep: **yes, MUST CONSIDER.**

### E3. `E3_historical_vs_postaudit`
- Question: how did the original RNA-eligibility gate change candidate
  standing vs. a CRISPR-only reordering?
- Source: `build_selection_rule_sensitivity()` (Rule 0 vs Rule 1, already frozen).
- Unit: gene (n=13).
- Encoding: two aligned tracks, connecting line only for genes that pass
  Rule 0.
- Interpretation: KDM1A (rank 1) and TLK2 (rank 2) were excluded purely by
  the RNA gate; VEZF1/USP34 passed despite weaker CRISPR ranks.
- Caveat: Rule 0 is reproduced from the ORIGINAL frozen gate code
  unmodified (not reimplemented) -- stated in the underlying data module.
- Strengths: the single clearest "how did the audit change the story"
  figure in the entire bank.
- Weaknesses: requires the viewer to understand what "Rule 0"/"Rule 1"
  mean (a caption dependency).
- Poster role: **exceptionally strong** -- directly visualizes the audit's
  central finding.
- Size: portrait, tall.
- Keep: **yes, MUST CONSIDER.**

---

## Section F -- human tumour / DepMap

### F1. `F1_depmap_heatmap`
- Question: real per-line DepMap dependency, all 4 genes, real names.
- Source: DepMap 26Q1 raw matrix + Model.csv names.
- Unit: cell line (n=11) x gene (n=4).
- Encoding: heatmap, real values annotated, real cell-line names.
- Interpretation: TLK2 column is uniformly darker (stronger dependency)
  across nearly every line; KDM1A/USP34 columns are pale throughout.
- Caveat: none -- fully real, no aggregation.
- Strengths: shows real heterogeneity AND real identities, not an
  anonymized summary.
- Weaknesses: none significant.
- Poster role: **strong core DepMap panel.**
- Size: portrait, moderate-tall.
- Keep: **yes, MUST CONSIDER.**

### F2. `F2_cellline_fingerprint`
- Question: does each cell line have a distinct 4-gene "fingerprint"?
- Source: same as F1.
- Unit: cell line (n=11), each as one connected line across 4 genes.
- Encoding: parallel-coordinates style, real per-line trajectories.
- Interpretation: MCF-7 and ACC-3133 show the deepest TLK2 dependency;
  fingerprints diverge substantially between lines.
- Caveat: none.
- Strengths: shows heterogeneity within a gene's own distribution, a
  different story than F1's grid.
- Weaknesses: 11 overlapping lines is busy even after de-collision.
- Poster role: strong alternative/companion to F1.
- Size: wide, landscape.
- Keep: **STRONG.**

### F3. `F3_human_depmap_combo`
- Question: combined view -- human/resistance-model evidence next to
  baseline dependency.
- Source: GSE240112 (VEZF1) + GSE118713 (USP34) + DepMap (all 4).
- Unit: mixed (tumour / sample / cell line).
- Encoding: 3-panel composite.
- Interpretation: places "why VEZF1/USP34 were retained" directly next to
  "what baseline dependency looks like" for context.
- Caveat: panel A/B use different genes (VEZF1/USP34 respectively) by
  design -- stated in footnote.
- Strengths: efficient use of poster space; answers two questions at once.
- Weaknesses: 3 different plot types in one figure demands a careful legend.
- Poster role: candidate consolidated F-section panel if poster space is
  tight.
- Size: wide, landscape.
- Keep: **STRONG.**

### F4. `F4_tcga_secondary`
- Question: TCGA evidence, deliberately kept small/secondary.
- Source: `TCGA_candidate_expression.tsv`, USP34/VEZF1 only.
- Unit: patient (paired, n as in TCGA-BRCA).
- Encoding: small forest plot.
- Interpretation: both weak/non-significant (FDR 0.21, 0.90).
- Caveat: original-4-only coverage; explicitly NOT a hero figure.
- Strengths: correctly sized to its actual evidentiary weight.
- Weaknesses: not much visual information (by design).
- Poster role: supplementary/backup, exactly as intended.
- Size: small, wide.
- Keep: **SUPPLEMENTARY** (by design).

---

## Section G -- structural / pharmacological comparison

### G1. `G1_three_structures_plus_vezf1`
- Question: how do the 3 solved structures compare, with VEZF1's absence
  shown honestly?
- Source: 6NQU (KDM1A, fetched this phase), 5O0Y (TLK2, fetched this
  phase), 7W3U (USP34, already frozen).
- Unit: structure (n=3 real PDB structures + 1 documented absence).
- Encoding: matched-style PyMOL cartoon renders, real ligands colored and
  labeled precisely (inhibitor / ATP-analog / covalent probe).
- Interpretation: visually striking three-way comparison; VEZF1 shown as a
  genuine gap, not a fabricated structure.
- Caveat: KDM1A's rendered view is cropped to the compact catalytic domain
  (excludes an ~100 A Tower-domain helix that would otherwise dwarf the
  ligand-bound region) -- disclosed in the source module's comments.
- Strengths: the most visually impressive figure in the entire bank; fixes
  the previous phase's "USP34-only structure" narrowness.
- Weaknesses: PNG-only (rasterized render).
- Poster role: **hero structural panel, strongest candidate in the bank.**
- Size: wide, landscape, large.
- Keep: **yes, MUST CONSIDER.**

### G2. `G2_structure_pharmacology_maturity`
- Question: same 3 structures + compact pharmacology-maturity labels.
- Source: same structures + Table 06b facets.
- Unit: same as G1 + one categorical label per gene.
- Encoding: same renders, small text labels below (no cards, no paragraphs).
- Interpretation: same structural comparison, framed toward drug-development
  stage rather than raw structural facts.
- Caveat: none beyond G1's.
- Strengths: same visual strength as G1, different narrative emphasis.
- Weaknesses: near-duplicate of G1 -- pick one, not both.
- Poster role: alternative to G1 depending on which caption emphasis suits
  the final poster (structural facts vs. drug-development stage).
- Size: wide, landscape, large.
- Keep: **STRONG** (choose G1 or G2, not necessarily both).

### G3. `G3_pocket_closeups`
- Question: what do the three catalytic/binding pockets actually look like?
- Source: same 3 structures, zoomed on the ligand/catalytic residue.
- Unit: same 3 structures, cropped views.
- Encoding: 3 large close-ups + VEZF1 omission stated in title.
- Interpretation: shows the literal binding geometry at high detail.
- Caveat: USP34's close-up is reused from the already-frozen figure bank
  (comp_closeup) rather than newly rendered -- same real structure, just
  not re-rendered from scratch.
- Strengths: may be more visually impressive than G1/G2's whole-protein
  views at poster viewing distance.
- Weaknesses: loses the whole-domain context G1 provides.
- Poster role: strong alternative hero panel, or a companion inset to G1.
- Size: wide, landscape.
- Keep: **STRONG.**

---

## Section H -- translational / follow-up

### H1. No summary figure (a decision, not a figure)
Tested whether the bank already tells the whole story without a forced
synthesis panel. Conclusion: E3 (historical vs. post-audit) and E2
(quantitative evidence map) already carry the "no universal winner"
message quantitatively; G1 carries the structural conclusion. A poster
ending on a strong structural hero image (G1) plus a 1-2 sentence written
conclusion is a legitimate, and possibly stronger, choice than forcing a
7th synthesis figure. See `EXPLORATION_V2_SUMMARY.md` for the final
recommendation on this question.

### H2. `H2_experimental_schematic`
- Question: what would the proposed follow-up experiments look like?
- Source: no data table (this is the one legitimately schematic figure) --
  arm labels are grounded in the real pharmacology-maturity facts from
  Table 06b (only KDM1A gets a pharmacological arm).
- Unit: N/A (schematic).
- Encoding: real biological visual language (cell glyphs, perturbation
  fill, treatment ring) -- not colored rectangles or a flowchart.
- Interpretation: makes explicit that KDM1A alone currently supports a
  pharmacological combination arm; the other 3 are genetic-perturbation-only.
- Caveat: proposed, not yet executed -- clearly a forward-looking schematic.
- Strengths: concise, biologically styled, avoids "PowerPoint" aesthetic.
- Weaknesses: inherently schematic (no real data plotted).
- Poster role: reasonable closing panel if a translational-logic figure is
  wanted at all.
- Size: portrait, moderate.
- Keep: **MAYBE** (only if the poster wants an explicit forward-looking panel).

### H3. `H3_translational_maturity_landscape`
- Question: CRISPR strength vs. an ordered pharmacological-maturity ladder.
- Source: `build_evidence_matrix()` + Table 06b (ordinal categorical, not a
  composite score).
- Unit: gene (n=13).
- Encoding: y = 5-level ordered category (no structure -> clinical-stage),
  x = CRISPR strength.
- Interpretation: KDM1A is uniquely at the top-right; TLK2/USP34 share the
  "probe/analog-bound structure" tier despite very different CRISPR strength.
- Caveat: y-axis is explicitly a single ordered read of one Table-06b
  column combination, not a weighted score (stated in footnote).
- Strengths: a genuinely different axis pairing than E2/F-series (drug-
  development readiness, not baseline dependency); data-grounded, not
  schematic.
- Weaknesses: the 5-level ladder compresses real nuance (e.g., "probe-bound"
  covers two very different kinds of evidence for TLK2 vs. USP34, disclosed
  elsewhere but not on this axis).
- Poster role: strong data-grounded alternative or companion to H2.
- Size: portrait, moderate.
- Keep: **STRONG.**

---

## Summary count by keep recommendation

- **MUST CONSIDER (12):** A1, A2, B1, B2, B3, B5a, C2, E1, E2, E3, F1, G1
  (several of comparable strength within a section; final shortlist
  narrows further to ~8-12, see `EXPLORATION_V2_SUMMARY.md`)
- **STRONG (9):** A4, B4, B5b, C3, F2, F3, G2, G3, H3
- **MAYBE (3):** C1b, D1, D2, H2
- **SUPPLEMENTARY (2):** A3, F4
- **REJECT (0):** none -- every built figure earned at least a supporting role.
