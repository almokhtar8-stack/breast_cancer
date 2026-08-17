# Poster Final Figure Guide

Six figures, selected and built from the already-frozen science (the
2026-08-15 `science-freeze-2026-08-15` tag plus the 2026-08-16 wording-only
cleanup patch). No scientific result, CRISPR/RNA/TCGA/DepMap/GDSC number,
selection-rule output, or sensitivity-analysis conclusion was recomputed or
altered to build these figures. Every plotted number is read directly from
`src/poster_final_data.py`, which documents its exact frozen source table
for each value.

Files: `results/figures/poster_final/F1_crispr_discovery.*` through
`F6_final_translational_framework.*` (PNG 300dpi + PDF; SVG for the five
purely-vector figures -- F5 is PNG-only, see its section below), plus
`POSTER_FINAL_CONTACT_SHEET.png` (all 6, for layout planning) and
`RETIRED_CANDIDATES_CONTACT_SHEET.png` (identical copy of
`results/figures/poster_candidates/POSTER_FIGURE_CONTACT_SHEET.png`, the
15-figure bank these six were selected/rebuilt from).

Framing carried through every figure, per the post-audit sensitivity
analysis and the 2026-08-16 re-audit cleanup: **there is no universal
winner** among the four focus genes (KDM1A, TLK2, USP34, VEZF1). KDM1A is
the strongest functional sensitiser and the most pharmacologically mature;
TLK2 has the strongest baseline ER+/luminal dependency among the four focus
genes (not automatically an advantage); USP34 is the most structurally
tractable novel DUB lead despite not being the top CRISPR hit; VEZF1 has
the strongest human recurrence-associated signal among the four focus genes
but belongs to a harder-to-drug target class. No composite or weighted
score appears anywhere in this figure set.

---

## F1. Genome-wide CRISPR discovery

**Shows:** All 19,103 fitted genes from the Hany et al. CRISPR screen
(effect size vs. -log10 FDR), with the 28 Gate-1 hits (FDR<0.1) marked, the
four focus genes highlighted and labeled with their exact rank (by effect
size and by FDR) out of the 13 significant sensitising hits, and RCOR1 --
the *other* pre-registered blind positive control (PREANALYSIS.md Section
5) -- shown as an open marker because it was **not** recovered at the
Gate-1 threshold (FDR=0.225).

**Question answered:** Where do the four focus genes sit in the full,
honest, genome-wide discovery landscape -- and is USP34 actually the
strongest hit?

**Why it earned a place:** It is the project's starting point and its most
important honesty check. It shows explicitly, in the same panel, that
**USP34 is not the top effect-size or top-FDR hit** (rank 12/13 and 10/13),
that KDM1A is (rank 1/13 both), and that the screen's own validity check
(RCOR1) partially failed -- a real, disclosed limitation, not hidden.

**Key message when presenting:** "CRISPR identifies multiple
tamoxifen-sensitising vulnerabilities, not one. USP34 is one of several
significant hits -- not the strongest by effect size or by significance.
One of the two pre-registered blind controls (KDM1A) was strongly
recovered; the other (RCOR1) was not, at the pre-declared threshold."

**Replaces:** `poster/01_crispr_discovery.png` (4-candidate-only,
frozen-shortlist view) and figure-bank `01a_crispr_ranked_effect` /
`01c_crispr_volcano_refined` (USP34/VEZF1-only highlighting, no KDM1A/TLK2,
no blind-control honesty check).

---

## F2. Transcriptomic / pathway systems view

**Shows:** A pathway (hallmark/GO:BP) x dataset NES heatmap across all four
transcriptomic contexts used in this project, with datasets explicitly
grouped by evidence category (resistance model / recurrence-associated /
acute) via colored axis labels, and FDR<0.05 significance marked per cell.

**Question answered:** Does the biology converge more strongly at the
pathway level than at the single-gene level, and does that convergence
hold consistently across genuinely different kinds of evidence?

**Why it earned a place:** It is the strongest systems-biology figure in
the whole project. Estrogen-response programs collapse consistently
(negative NES, mostly FDR<0.05) across all four datasets including the
acute 12h context; EMT/ECM/Wnt programs rise consistently across the three
non-acute datasets but the acute-12h dataset diverges on several of them
(e.g. EMT, UV-response-down) -- a real, visible, disclosed acute-vs-chronic
distinction, not smoothed over.

**Key message when presenting:** "The four transcriptomic datasets
converge much more strongly at the pathway level than at any single-gene
level. GSE240112 is a recurrence-associated human-tumour comparison, not a
resistance model -- it is never grouped with the two cell-line resistance
datasets."

**Replaces:** `poster/02_expression_evidence.png` and figure-bank
`06_resistance_pathway_landscape` (both omit the acute-12h context and
lacked an explicit resistance-model/recurrence-associated/acute category
distinction).

---

## F3. Candidate evidence divergence / integration

**Shows:** A 2x2 grid, one real evidence dimension per panel, for the four
focus genes only: (A) CRISPR sensitisation strength; (B) a corroboration
dot-matrix across GSE118713, GSE111151, GSE240112 and TCGA (with an
explicit "n/a" for KDM1A/TLK2's TCGA cells -- TCGA was only ever run on the
original 4-candidate set, a real evidence-depth gap, not a null result);
(C) baseline ER+/luminal DepMap dependency, labeled double-edged; (D) a
structural/pharmacological tractability dot-matrix (experimental structure,
ligand/probe-bound structure, validated selective inhibitor, clinical-stage
pharmacology) read directly from Table 06b's own cited facets.

**Question answered:** Why do different candidates survive, given that no
single metric picks a clear winner?

**Why it earned a place:** This is the figure that makes the "no universal
winner" claim inspectable rather than asserted. No panel is a composite or
weighted score; each is a single, named, real evidence axis. It makes
KDM1A's near-total pharmacological lead (panel D) and USP34's comparatively
weak CRISPR rank (panel A) equally visible in the same figure.

**Key message when presenting:** "This is an evidence-integration figure,
not a ranking. Different candidates lead on different axes -- KDM1A leads
CRISPR strength and pharmacology, TLK2 leads baseline dependency, USP34 and
VEZF1 have the most direct human/RNA corroboration among the four, and
VEZF1's own dependency profile is deliberately shown on the same
double-edged terms as TLK2's."

**Replaces:** `poster/04_pharmacogenomics.png`-style head-to-head framing
and any implicit "USP34 > VEZF1 > EML5 > CITED2" ranking figure -- this
figure explicitly could not be built as a single-score ranking, by design.

---

## F4. Human (TCGA) / DepMap orthogonal validation

**Shows:** (A) TCGA-BRCA paired tumour-vs-normal log2FC + 95% CI for the
two focus genes TCGA evidence exists for (USP34, VEZF1), with KDM1A/TLK2
explicitly marked "not assessed in this project" rather than omitted or
zeroed; (B) DepMap 26Q1 per-cell-line Chronos gene-effect distributions
(11 ER+/luminal screened lines) for all four focus genes.

**Question answered:** Does orthogonal, independent evidence (human tumour
expression; cell-line essentiality) distinguish tamoxifen-specific
sensitisation from generic cancer-cell dependency?

**Why it earned a place:** It is the clearest demonstration that CRISPR
sensitisation strength (Figure F1) and baseline dependency (this figure,
panel B) are different quantities: KDM1A has the strongest sensitisation
signal but low baseline dependency, while TLK2 has the weaker CRISPR signal
of the two yet the strongest baseline dependency of the four focus genes.

**Key message when presenting:** "High baseline dependency is not
automatically an advantage -- for TLK2 (81.8% of 11 lines) and for VEZF1
(27.3%) alike, it may reflect a narrower, less tamoxifen-specific
therapeutic window rather than a purely favourable signal."

**Replaces:** `poster/03_depmap_distributions.png` (candidate-only, no
KDM1A/TLK2, no explicit double-edged framing) and figure-bank
`08_TCGA_human_validation` / `09_Hany_vs_DepMap_context_map` (this figure
combines both into one orthogonal-validation panel with the TCGA coverage
gap made explicit).

---

## F5. USP34 structure and covalent tractability

**Shows:** Three real, ray-traced PyMOL renders of the two frozen USP34 PDB
structures (7W3R apo, 7W3U covalent ubiquitin-probe complex), matched
orientation, with Cys1903/His2164 highlighted and a catalytic-cleft
close-up, plus a caption strip stating the tractability facts precisely.

**Question answered:** Given that USP34 is not the top CRISPR hit (Figure
F1), why does it remain a structurally distinctive candidate?

**Why it earned a place:** It is the most visually striking figure in the
set and documents a real, specific, experimentally-observed fact (direct
catalytic-cysteine covalent reactivity via an activity-based ubiquitin
probe) without overclaiming druggability.

**Key message when presenting:** "Direct catalytic-cysteine reactivity is
experimentally observed -- but only via a covalent activity-based ubiquitin
probe, a chemical-biology tool, not a small-molecule ligand. No validated
selective USP34 inhibitor exists, and docking was not pursued in this
project because no validated ligand set exists for calibration."

**Replaces:** `poster/05_structure.png` and figure-bank
`12_USP34_structure_surface` / `12b_USP34_structure_comparison` (this
figure combines the clearest elements of both into one three-panel
composition with a tighter, more precise caption).

**Format note:** PNG only, no PDF/SVG -- the panel content is a rasterized
ray-traced render (PyMOL `cmd.ray` + `cmd.png`), the same format
constraint already documented for the earlier `poster/05_structure.png`
and figure-bank structural figures; a vector re-export would not contain
meaningfully different information and is not attempted.

---

## F6. Final candidate logic / experimental strategy

**Shows:** (A) A real, data-grounded scatter of all 13 significant
sensitising genes -- CRISPR sensitisation strength (-log10 FDR) on the
x-axis, baseline ER+/luminal DepMap dependency (%) on the y-axis -- with
the four focus genes colored and labeled; (B) four role cards (KDM1A, TLK2,
USP34, VEZF1), each with a one-paragraph role summary and a proposed
follow-up assay, using the same role language already documented in the
post-audit sensitivity report.

**Question answered:** What is the translational landing point of this
project, and what would be tested next for each candidate?

**Why it earned a place, and why a data-grounded panel A rather than a pure
schematic:** Panel A visually proves, with real numbers, that
sensitisation strength and baseline dependency are not the same axis (the
four focus genes occupy four different quadrants) -- a stronger and more
defensible claim than a schematic could make on its own. Panel B is
schematic by necessity (proposed follow-up assays have not been run), but
its role text is not a new ranking rule -- it restates, in poster form, the
distinct-candidate-roles interpretation already written in
`results/reports/post_audit/POST_AUDIT_SENSITIVITY_REPORT.md` Section 24.

**Key message when presenting:** "This project yields a testable,
multi-candidate translational framework, not a single overclaimed winner.
Each candidate's proposed next experiment follows directly from what
distinguishes it in panel A and in Figure F3."

**Replaces:** `poster/06_experimental_strategy.png` (candidate-only,
schematic-only, no data-grounded role-compass panel) and figure-bank
`13_validation_experiment` (single-schematic, no KDM1A/TLK2).

---

## How the external audits changed the final figure strategy

**What was changed.** The original 6-figure poster set (`results/figures/
poster/`) and the earlier candidate figure bank (`results/figures/
poster_candidates/`) were built before the post-audit sensitivity analysis
and the two rounds of external re-audit. Both were built around the frozen
4-candidate shortlist (USP34 > VEZF1 > EML5 > CITED2) as if it were an
uncontested ranking. Every final figure in this set instead foregrounds
KDM1A and TLK2 -- the two strongest functional CRISPR hits that the
original RNA-eligibility gate excluded from that shortlist -- alongside
USP34 and VEZF1, and none of the six figures computes or displays a single
composite score.

**What was retired.** The implicit "USP34 is the winner" framing is retired
outright: Figure F1 states plainly, with the real rank numbers, that USP34
is neither the top effect-size nor the top-FDR CRISPR hit. The GDSC
pharmacogenomics figure (`11_GDSC_USP34_pharmacogenomics`) is retired from
the main sequence entirely, consistent with its SUPPORTING/SUPPLEMENTARY
status already documented in `results/reports/poster/
FIGURE_BANK_REVIEW.md` (updated 2026-08-16 with the GDSC1-only,
sign-reversal-caveated downgrade) -- it is not part of any of the six final
figures and no new GDSC analysis was performed.

**What is now emphasized.** (1) KDM1A's status as the strongest functional
sensitiser and most pharmacologically mature candidate is stated plainly in
F1, F3, and F6, not buried. (2) TLK2's strong baseline DepMap dependency is
never presented as a pure advantage -- F3 panel C and F4 panel B both carry
the double-edged framing explicitly, and F6 shows the same caveat applies
symmetrically to VEZF1. (3) GSE240112 is labeled "recurrence-associated"
and visually separated from the two resistance-model datasets in every
figure that uses it (F2, F3), never "chronic resistance." (4) The
KDM1A-vs-USP34 comparison in F3 is evidence-specific (probe-based catalytic
reactivity vs. established covalent small-molecule pharmacology), never a
generic "strongest covalent target" claim. (5) TCGA's real coverage gap
(only ever run on the original 4-candidate set) is shown explicitly in F3
and F4 rather than silently omitted or defaulted to zero.

**Why the final figure sequence is more defensible.** Every claim a viewer
could take away from these six figures traces to a specific frozen table
cited in `src/poster_final_data.py`, and every acknowledged weakness in the
post-audit sensitivity analysis (USP34's CRISPR rank, RCOR1's non-recovery,
TCGA's coverage gap, the acute-vs-chronic transcriptomic distinction, the
double-edged nature of high baseline dependency) is made visible in at
least one figure rather than smoothed over. No figure asserts a single
winner; Figure F6 explicitly states "no universal winner" as its headline.

---

## Ranked recommendation

**Main poster set (all 6, in this order):** F1 -> F2 -> F3 -> F4 -> F5 ->
F6. This order follows the project's own narrative arc (discovery ->
systems context -> why these candidates -> orthogonal validation ->
structural deep-dive on the most novel lead -> translational framework),
matching the poster title "From a CRISPR Screen to Therapeutic
Vulnerabilities."

**Supporting / optional (not main-sequence, available as backup slides or
a poster appendix):** figure-bank `07_candidate_mechanism_map` (why EML5
was deprioritized, useful if asked), `09b_USP34_VEZF1_line_dependencies`
(a drill-down on F4 panel B for USP34/VEZF1 specifically),
`10_tissue_liability_context` (safety-diligence backup), `13
_validation_experiment` (an alternative pure-schematic take on F6 panel B,
useful if a reviewer specifically wants a non-data-grounded flow diagram).

**Retired (superseded by the six final figures above, or explicitly
excluded per the science-freeze GDSC downgrade):** `poster/01-06_*.png`
(the original 4-candidate-only set), figure-bank `01a`, `01c`, `03`, `04`,
`05`, `06`, `08`, `09`, `12`, `12b` (each superseded by the corresponding
final figure above), and `11_GDSC_USP34_pharmacogenomics` (excluded from
the main narrative per its SUPPORTING/SUPPLEMENTARY status, kept in the
figure bank, not deleted). See `RETIRED_CANDIDATES_CONTACT_SHEET.png` for
the full 15-figure visual index this set was selected from.
