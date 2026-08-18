---
title: Figure set note — what changed per figure and why
status: post_freeze_exploratory
analysis_date: 2026-08-18
branch: poster-final (unmerged)
---

# The final figure set

Seven figures, in `results/figures/poster_final/`, each in PNG, PDF and SVG,
with `figure_manifest.tsv` (sources, plotted values, SHA-256 per file, printed
point sizes), `figure_captions.tsv`, `verification_against_frozen.tsv` and
`cvd_simulation/`.

No statistic was recomputed anywhere. Every plotted value is read at render
time from a table already committed to this repository, and **every figure
verifies its own values against the frozen source before it draws anything**
and raises rather than substituting — 87 checks in total, all passing. That
discipline is not decorative: during earlier volcano work the same gate caught
two wrong source files, one of which would have printed a materially different
false discovery rate on the poster.

---

## Figure 1 · Methods workflow — **new**

**Replaces:** nothing. The project had no methods diagram; the methods were
prose, and harder to follow than those of any poster in the reference bank.

**What it shows:** the screen and its output; the four transcriptomic datasets;
the candidate check and the pathway branch coming off the same data; and the
dependency and structural assessments as separate downstream checks. A dashed
vertical line marks where the analysis plan and the candidate list were fixed
relative to the external data.

**Why it is better:** it does the work several paragraphs of methods prose were
doing, and it makes the project's strongest methodological claim visible
instead of buried. It also marks what was added *after* the freeze in a dashed
box, so the diagram cannot be read as claiming the whole study was specified
in advance — a point the independent reviewer raised at checkpoint 1 and which
would otherwise have been a real overclaim.

---

## Figure 2 · Screen, with certainty visible

**Replaces:** frozen figure 01 (`poster_crispr_discovery_v1`), which is
untouched and stays in the repository.

**Diagnosed problem:** figure 01 sorts the 13 hits by effect size alone, so a
reader cannot see why USP34 (rank 12 of 13 by effect) and VEZF1 (rank 8) were
carried forward while USP17L29 (rank 2) and TADA2B (rank 3) were not. A
careful reader was in fact confused by exactly this.

**What changed:** effect size and certainty are now two axes, so a large effect
with weak certainty is visibly distinct from a small effect with strong
certainty. Marker **shape** carries the second, more important disclosure: the
four candidates did not come from one selection rule. Circles entered by the
frozen multimodal rule; squares were added after an external audit and have no
qualifying corroboration. EML5 and CITED2, displaced from the original frozen
shortlist by the same reinterpretation, are named in the caption.

**Why it is better:** the figure now answers the question the old one provoked.
A reader can see that KDM1A and TLK2 are the most *certain* hits, that USP34 is
not, and that two different rules produced the set.

---

## Figure 3 · Candidate corroboration

**Replaces:** frozen figure 02 (`poster_hero_heatmap_v6`), untouched.

**Diagnosed problem:** the heatmap standardises colour within each dataset, so
noise is stretched across the full blue-to-red range. A reader misread genes at
false discovery rate 0.49 as strongly changed.

**What changed:** significance is on an axis, where it cannot be misread. Four
volcano panels share one x and one y range. Filled points reach FDR 0.05,
hollow rings do not — exactly two filled points exist across all sixteen
combinations. The gate and the zero-displacement rule are inherited unchanged
from `poster_candidate_volcano_v2` (merged from the `figure2-volcano` branch,
as the brief directed, rather than rebuilt): coincident candidates are drawn as
concentric rings at their true position, and a test asserts that plotted x
equals source log2 fold change exactly for all sixteen points.

**What was added beyond the merged work:** a fourth row carrying the pooled
evidence — the random-effects meta-analysis (no candidate reaches pooled FDR
0.05; smallest 0.100) and the GSE111151 minimum-detectable-effect result. The
independent reviewer identified the absence of these as the figure set's one
real gap: without them, "the candidates do not corroborate" reads as far more
decisive than the evidence permits. The strip also states the limit of the
power argument — it explains GSE111151's null, not every null.

**Why it is better:** it makes the 2-of-16 count legible at a glance *and*
prevents that count from being over-read as a strong negative.

---

## Figure 4 · Programme-level signal

**Replaces:** frozen figure 04 (`poster_pathway_v2`), untouched.

**Diagnosed problem:** row labels were raw gene-set identifiers
(`HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION`) that a non-specialist cannot
parse, and the project's one novel observation had no more prominence than
anything else on the figure.

**What changed:** each row is labelled by the job the programme performs, with
the formal gene-set name retained small underneath for traceability. The
adhesion-and-motility row is boxed and is the visual centre. A rule separates
the three long-term settings from the acute one.

**Why it is better:** a stranger can now read the row labels, and the
dissociation — positive enrichment in three long-term settings, negative after
twelve hours — is the first thing the eye lands on.

**What it deliberately does not say:** the reviewer warned that calling this a
"mobility programme" measurement, or a temporal "dissociation", invites a
causal reading the design cannot support. The caption therefore states that
this is a gene-set enrichment score rather than a measured migration
phenotype, and that the acute dataset differs from the other three in more than
duration.

---

## Figure 5 · Network connectivity

**Replaces:** frozen figure 03 (`poster_network_mechanism_v4`), untouched.

**Diagnosed problem:** overcrowded, and its layout was not deterministic —
label placement used `adjustText`, whose solver moves label pixels between
runs, which is why figure 03's PNG is documented as not byte-reproducible.

**What changed:** the layout is seeded and each connected component is laid out
in its own slot, so an isolated node cannot push the main component into a
corner. Labels sit at fixed offsets; **no collision solver is used anywhere**,
and a test parses the module's syntax tree to prove it (prose mentioning
`adjustText` in the docstring is exempt on purpose — that docstring records why
it is *not* used). Label density is cut to the candidates, the bridge and the
ubiquitin nodes. The PNG is byte-reproducible across renders, test-enforced.

**What panel B adds:** it draws the connection as what it is. KDM1A and USP34
are three associations apart with four equally short routes — but every route
passes through DNMT1 and then through one of four genes that all encode
ubiquitin. Drawing that is what stops "connected in a network" from sounding
like evidence of a mechanism.

**Why it is better:** deterministic, readable, and it makes the weakness of the
connection the point of the figure rather than a caveat underneath it. The
figure was **not** replaced with a screenshot from any external site, which
would have lost the pinned query, the manifest, the hash and the determinism
test.

---

## Figure 6 · Baseline dependency

**Replaces:** frozen figure 05 (`poster_depmap_v2`), untouched.

**Diagnosed problem:** a scatter whose x-axis restated figure 1, carrying one
number per gene, with no numeric x scale a reader could use.

**What changed:** a bar chart of **counts, not percentages**, out of a stated
denominator of 11, sorted so the outlier is immediate. A faint full-width bar
shows the denominator behind every gene, so "0 of 11" is visibly zero-out-of-
eleven rather than an absent bar. 81.8% is simply 9 of 11, and a percentage of
eleven implies a precision the sample size does not support.

**What the scatter carried that bars cannot:** its geometry asserted that
sensitisation and baseline requirement are separate measurements. That claim is
recovered in the panel title — but in the reviewer's wording, not the original:
they are "two different measurements" that "**need not coincide**", not
"independent axes". Different counts do not establish statistical independence
and the figure does not claim it.

---

## Figure 7 · Chemical reachability

**Replaces:** frozen figure 06 (`poster_druggability_v1`), untouched.

**Diagnosed problem:** text-heavy relative to its renders.

**What changed:** a five-level evidence track a reader can scan in one pass,
with the three committed PyMOL renders given the space. The levels are kept
separate and never summed into a score, because the four candidates differ in
**kind** of evidence, not degree: an inhibitor-bound co-crystal, an
ATP-analogue-bound kinase domain, a covalent activity-based probe and a
homology-model screening hit are not four points on one axis.

**VEZF1 has no experimental structure, and none was substituted.** Its panel is
a dashed empty frame saying so. No predicted model appears anywhere, in that
slot or beside it, because the absence is itself the finding.

**What was added:** the published counter-evidence for the lead candidate —
losing USP34 has been reported to push breast cells toward a more mobile state
(PMID 28499884) — appears on the figure. The reviewer's point was that giving
USP34 visual prominence without it would make the poster look selectively
curated.

---

## Cross-cutting changes

**One palette, defined once.** `src/poster_palette.py` is the only place a
colour is written down. Tests assert that no hex literal appears outside it,
that no poster-chrome colour reaches a figure, and that nothing imports colour
from a v1 or v2 renderer. Two colours appear in both the chrome and the data
tiers by design — the violet is USP34's colour *and* the poster violet; the
warm white is the diverging ramp's midpoint *and* the chrome warm white — and
the test pins exactly those two, so a third overlap fails.

**Colour identifies the gene, never its rank or significance.** Significance is
fill: solid passes, hollow ring in the same colour fails. Identical on Figures
2, 3, 4 and 7.

**Colour-vision safety is computed, not asserted.** Deuteranopia and protanopia
are simulated with the Machado, Oliveira & Fernandes (2009) matrices applied in
linear sRGB — no colour-vision library exists in this environment and none was
added. All four genes stay distinguishable under both: the tightest pair is
USP34 (purple) against VEZF1 (light blue) at ΔE 32.3 under deuteranopia and
40.9 under protanopia, separating on **lightness** (ΔL\* 32–38) rather than
hue, exactly as anticipated. **No darkening of the purple was required.**
Every rendered PNG also has both simulations written to `cvd_simulation/`.

**Captions are poster text, not figure text.** This changed late and for a
concrete reason: with captions baked into the figures, the explanatory lines
printed at 12–13 pt at the intended placement, well under the 20 pt floor.
They are now recorded into `figure_captions.tsv` and typeset in
`POSTER_TEXT.md` at the poster's own 24 pt, and every remaining string on
every figure clears the floor — checked programmatically per figure in the
manifest (`meets_20pt_floor`), not by eye.

**Reproducibility.** `SOURCE_DATE_EPOCH` is pinned and SVG element-id salts are
fixed, so PNG and PDF bytes are identical across renders for all seven figures
— test-enforced, and a stronger guarantee than the frozen set, where figure 03
was not byte-reproducible at all.
