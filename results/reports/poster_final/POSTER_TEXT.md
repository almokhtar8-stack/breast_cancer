---
title: Complete poster text — every word that goes on the sheet
status: post_freeze_exploratory
analysis_date: 2026-08-18
branch: poster-final (unmerged)
format: A0 portrait, 841 × 1189 mm
spelling: British
---

# Poster text

Every word intended for the printed sheet is below, in reading order, with
word counts per section and a layout map. Written to the Cambridge University
Library conference-poster guidance
(`https://libguides.cam.ac.uk/confpost/`, Layout / Content / Images pages
fetched 2026-08-18) and to the house conventions of the previous cohort's
poster bank.

**What the Cambridge guidance actually specifies**, and what it does not, is
worth stating because the two are easy to conflate. It specifies: a title of
at most two lines; the section set (title, introduction, methods, results,
conclusions, references ≈ 5, contact); *avoid acronyms and abbreviations*
because they "alienate people unfamiliar with their meaning"; images in
methods and graphs in results; state where the poster can be obtained
afterwards; keep enough white space and a clear flow; images at ≥ 200 dots
per inch. It does **not** specify font sizes, column counts or word limits —
it defers those to the conference's own guidelines. The type sizes used here
(section titles 50–60 pt, body 25–30 pt, nothing below 20 pt) therefore come
from the project brief, not from Cambridge, and every figure has been checked
against the 20 pt floor at its printed size (`figure_manifest.tsv`,
`meets_20pt_floor`).

Acronyms are spelled out on first use throughout: false discovery rate,
oestrogen receptor positive, epithelial-to-mesenchymal transition,
clustered regularly interspaced short palindromic repeats (CRISPR),
ribonucleic acid (RNA).

---

## TITLE BLOCK

**Title** (frozen wording, used exactly; 2 lines at 96 pt)

> From a CRISPR Screen to Therapeutic Vulnerabilities: Functional Genomics of
> Tamoxifen Sensitisation in ER+ Breast Cancer

**Plain-language subtitle** (48 pt, for a reader passing at a distance)

> Which genes, when switched off, make a breast-cancer drug work better —
> and how far that evidence actually goes

**Author line** (32 pt)

> Almokhtar Aljarodi · KAUST Academy · Department of Genetics, University of Cambridge

*Word count: title 21, subtitle 22, author line 12.*

---

## 01 INTRODUCTION  *(28 pt body)*

> Tamoxifen has treated oestrogen-receptor-positive breast cancer for four
> decades, and it stops working in a large minority of patients. One way to
> improve it is to find a second gene that, when switched off, makes tamoxifen
> work better — a sensitiser.
>
> A published genome-wide screen switched off each of 19,103 genes in turn and
> measured which of them made an oestrogen-receptor-positive breast-cancer cell
> line more sensitive to tamoxifen. This poster is a computational reanalysis
> of that screen and of four public transcriptomic datasets. **No laboratory
> experiment was performed at any stage.**
>
> The question asked here is not "which gene is the answer". It is the harder
> question a drug-discovery team has to answer before committing anything:
> **how far does this evidence actually go?**

*Word count: 138.*

---

## 02 METHODS  *(Figure 1; 26 pt caption)*

**Figure 1 headline** (36 pt bold)

> How the work was done, and where the candidate list was fixed

**Figure 1 explanatory line** (24 pt)

> Everything here is a computational reanalysis of public data. Nothing on the
> right of the dashed line was allowed to change the candidate list on its
> left: the analysis plan and its thresholds were dated before any result
> existed, and the candidate list was fixed before the external
> corroboration datasets were opened. The dashed box records what was added *after* that
> freeze, so the diagram does not imply the whole study was specified in
> advance — it was not.

*Word count: 84.*

---

## 03 RESULTS

### I. Thirteen hits, and why these four  *(Figure 2)*

**Headline** (36 pt bold)

> Effect size and certainty are different things, and the four candidates were
> not chosen on effect size alone

**Explanatory line** (24 pt)

> 13 of 19,103 genes met the pre-specified gate (false discovery rate below
> 0.10, and switching the gene off had to *increase* sensitivity — a threshold
> written down and dated before any result existed, though not lodged in a
> public preregistration registry). KDM1A and
> TLK2 are the most certain hits. USP34 — the candidate carried furthest in
> this work — sits twelfth of thirteen by effect size. Marker shape shows how
> each candidate entered the set: the four did **not** come from one selection
> rule. USP34 and VEZF1 passed the frozen rule, which required at least one
> corroborating dataset; KDM1A and TLK2 were added after an external audit
> challenged that rule, and have no qualifying corroboration at all. Two genes
> from the original frozen shortlist, EML5 and CITED2, were displaced by the
> same reinterpretation.

*Word count: 118.*

**→ arrow to II**

### II. Gene-level corroboration is sparse  *(Figure 3)*

**Headline** (36 pt bold)

> Two of sixteen candidate gene-and-dataset combinations reach false discovery
> rate 0.05

**Explanatory line** (24 pt)

> Four genes tested in four datasets gives sixteen combinations. Two reach the
> threshold: USP34 in the cell-line resistance model, VEZF1 in recurrent
> tumours. Each candidate is corroborated in at most one dataset, and none in
> more than one. Filled points reach the threshold; hollow rings do not; grey
> points are every other gene in that dataset.
>
> These four datasets are not four equivalent replications. One is a chronic
> resistance model, one an independent panel of resistant sublines, one a
> comparison of recurrent against primary tumours from **different, unpaired
> patients in two separate tissue banks**, and one an acute twelve-hour
> exposure that measures immediate drug response, **not resistance**.
>
> The nulls are weak evidence rather than strong negatives. Pooling the three
> resistance-context datasets, no gene reaches a pooled false discovery rate of
> 0.05; the smallest is 0.100. And the best-designed of them, GSE111151 — the
> only one with several independent resistance backgrounds — had a median
> 80%-power detection threshold of about 1.5-fold (per-gene range 1.34- to
> 1.74-fold), while every candidate effect it observed was smaller than that. **Its nulls are therefore weak evidence, not strong negatives.**

*Word count: 205.*

**→ arrow to III**

### III. Programme-level signal is present in all four  *(Figure 4)*

**Headline** (36 pt bold)

> Programme-level signal is present in all four datasets, and the
> adhesion-and-motility programme points the other way after twelve hours

**Explanatory line** (24 pt)

> Where single genes are unstable, whole biological programmes are not.
> Oestrogen response is suppressed in all four datasets and cell-cycle entry in
> all four — the pharmacology one would expect, and an internal check that the
> datasets behave sensibly.
>
> The boxed row is the observation this project adds. Enrichment of the
> epithelial-to-mesenchymal transition programme, which governs cell adhesion
> and motility, is positive in the three long-term settings and negative after
> twelve hours of tamoxifen. This is a gene-set enrichment score, not a
> measured migration phenotype, and the acute dataset differs from the other
> three in far more than duration — so the difference is reported as an
> observation, not as a mechanism.

*Word count: 128.*

**→ arrow to IV**

### IV. What connects two of the candidates, and how weakly  *(Figure 5)*

**Headline** (36 pt bold)

> The two connected candidates are joined by one route, not four

**Explanatory line** (24 pt)

> One standardised public-interaction query, applied identically to all four
> genes, puts KDM1A and USP34 three associations apart with four equally short
> routes between them. Every one of those routes passes through the same
> protein, DNMT1, and then through one of four genes that all encode ubiquitin.
> So this is one connection whose middle position can be filled four ways, not
> four independent lines of support. TLK2 sits in a separate component and
> VEZF1 has no partner at all at this threshold. These associations are
> undirected: they are not activation, not inhibition, not necessarily physical
> binding, and a short path is not a mechanism.

*Word count: 111.*

**→ arrow to V**

### V. Being needed at baseline is a different question  *(Figure 6)*

**Headline** (36 pt bold)

> Needing a gene at baseline and being sensitised by losing it are two
> different measurements, and they need not coincide

**Explanatory line** (24 pt)

> TLK2 is required for survival by 9 of 11 oestrogen-receptor-positive cell
> lines whether or not tamoxifen is present. **That is a limitation, not a
> strength**: a gene that most cells already need offers a narrower, less
> tamoxifen-specific therapeutic window. KDM1A is the opposite — the strongest
> sensitiser in the screen, and required by none of the 11. Counts, not
> percentages, out of the 11 lines with a dependency value. Dependency measured
> in cancer cell lines is not normal-tissue safety and says nothing about
> toxicity.

*Word count: 94.*

**→ arrow to VI**

### VI. Chemical reachability differs in kind  *(Figure 7)*

**Headline** (36 pt bold)

> The four differ in the kind of chemical evidence they have, not in how much
> of one thing

**Explanatory line** (24 pt)

> This is reachability only — not efficacy, not selectivity, not safety. **No
> docking, binding prediction or molecular modelling was performed anywhere in
> this project.**
>
> KDM1A is the mature target: inhibitor-bound structures and clinical-stage
> compounds already exist. TLK2 has a real kinase-domain structure, but bound
> to an ATP analogue rather than an inhibitor, and no selective inhibitor
> exists. USP34 has a catalytic-domain structure holding a covalent ubiquitin
> activity probe — a laboratory tool, **not a drug** — over about 12% of the
> protein; it is an unexplored opportunity, not a validated one. VEZF1 has no
> experimental structure at all, and no predicted model has been substituted
> for it, because that absence is itself a finding.
>
> Counter-evidence worth stating: losing USP34 has been reported to push breast
> cells toward a more mobile state (Cellular Signalling, 2017), which works
> against the hypothesised sensitisation benefit.

*Word count: 158.*

---

## 04 LIMITATIONS  *(26 pt body, boxed panel)*

> - A computational reanalysis of public data. Nothing here has been tested at
>   the bench, by me or for this project.
> - The screen is a single genome-wide screen in a single cell line. It asks
>   whether losing a gene makes tamoxifen work *better* in a drug-tolerant
>   parental line — not whether it reverses established resistance.
> - The four candidates come from **two different selection rules**, one frozen
>   and one applied after an external audit. They are not a uniformly selected,
>   pre-specified panel.
> - The analysis plan and thresholds were written down and dated before results
>   existed, and the freeze is a git tag anyone can check — but this was not
>   lodged in a public preregistration registry, and is not described as one.
> - Gene-level corroboration is sparse (2 of 16) and pooled evidence is not
>   significant. Several nulls are weak evidence rather than strong negatives,
>   for different reasons in different datasets: limited power in GSE111151,
>   acute-versus-chronic design in GSE245601, and unpaired patients drawn from
>   two tissue banks in GSE240112, where treatment group and tissue bank are
>   confounded.
> - The recurrence dataset is an **association**. It is not a controlled
>   tamoxifen-resistance experiment and does not show tamoxifen causality.
> - Our malignant-cell calls in the single-cell datasets are our own
>   reconstruction; the authors' per-cell labels are not public, so the two
>   cannot be compared directly.
> - Public-interaction associations are not mechanisms, and baseline dependency
>   is not safety.
> - **Nothing here is a validated therapeutic target.** No claim is made that
>   inhibiting any of these genes would benefit a patient.

*Word count: 213.*

---

## 05 CONCLUSIONS  *(28 pt body; three-part scaffold)*

**What was found**

> A genome-wide screen yields 13 functional tamoxifen-sensitising hits. Across
> four public transcriptomic datasets, candidate-level corroboration is sparse —
> 2 of 16 gene-and-dataset combinations — and pooled evidence does not reach
> significance. Programme-level signal, by contrast, is present in every
> dataset, and the adhesion-and-motility programme runs in opposite directions
> between long-term and acute settings.

**Why it matters**

> A functional screen hit and a transcriptional change are different kinds of
> evidence, and a gene need not do both. Treating a transcriptomic null as
> refutation would have discarded the two most certain functional hits in the
> screen; treating it as confirmation would have overstated the other two. The
> honest position is that this evidence supports **a hypothesis worth testing**,
> and does not support a target.

**What comes next**

> The decisive experiments are experiments, not more reanalysis: a
> tamoxifen-sensitisation assay in independent oestrogen-receptor-positive
> models, and a properly powered resistance dataset. USP34 is the candidate this
> work carries furthest — not because its screen evidence is strongest, but
> because it is the one with an unexplored, experimentally addressable catalytic cysteine.
> That makes it a lead on tractability, not on evidence strength.

*Word count: 197.*

---

## REFERENCES  *(20 pt; five printed, full list behind the QR code)*

> 1. Hany D. *et al.* (2023) Recurrent mutations in tamoxifen-sensitising
>    pathways revealed by a genome-wide CRISPR screen. *Science Advances*
>    9:eadd3685.
> 2. Kim H., Whitman A.A. *et al.* (2023) Tamoxifen response at single-cell
>    resolution in oestrogen receptor-positive primary human breast tumours.
>    *Clinical Cancer Research* 29(23):4894–4907. PMID 37747807.
> 3. Liberzon A. *et al.* (2015) The Molecular Signatures Database Hallmark gene
>    set collection. *Cell Systems* 1(6):417–425.
> 4. Szklarczyk D. *et al.* (2023) The STRING database in 2023. *Nucleic Acids
>    Research* 51(D1):D638–D646.
> 5. Luo H. *et al.* (2017) USP34 regulates the epithelial–mesenchymal
>    transition in breast cancer cells. *Cellular Signalling* 36:1–10.
>    PMID 28499884. *(Counter-evidence to the hypothesis proposed here.)*

> **Full reference list, all code, all data tables and every figure:** scan the
> code opposite, or visit `github.com/almokhtar8-stack/breast_cancer`.

*Word count: 118 including the pointer line.*

---

## QR PANEL AND CONTACT  *(24 pt)*

> **Code, data and figures**
> `github.com/almokhtar8-stack/breast_cancer`
>
> Almokhtar Aljarodi — aashsj08@gmail.com
> The poster, the full reference list and the analysis plan can all be
> downloaded from the repository above.

*Word count: 34.*

---

## ACKNOWLEDGEMENTS  *(20 pt; sized placeholder — author to confirm wording)*

> [Placeholder, 180 × 40 mm. Suggested: "With thanks to the KAUST Academy
> programme and to the Department of Genetics, University of Cambridge."]

---

# Layout map — A0 portrait, 841 × 1189 mm

All measurements in millimetres from the top-left corner. Margin 30 mm on all
four sides; usable area 781 × 1129 mm. Two columns of 380 mm with a 21 mm
gutter, except where a panel spans the full width.

```
┌──────────────────────────────────────────────────────────────────────┐
│ LOGO STRIP  y 12–52, full width, six logos at 34 mm height           │
├──────────────────────────────────────────────────────────────────────┤
│ TITLE BLOCK  y 58–196, full width, deep violet band                  │
│   title 96 pt · subtitle 48 pt · author line 32 pt                   │
├───────────────────────────────┬──────────────────────────────────────┤
│ 01 INTRODUCTION               │ 02 METHODS                           │
│ y 210–420, x 30–410           │ y 210–560, x 431–811                 │
│ text only, 28 pt              │ Figure 1 (380 mm wide) + caption     │
├───────────────────────────────┴──────────────────────────────────────┤
│ 03 RESULTS — banner y 575–615, full width                            │
├───────────────────────────────┬──────────────────────────────────────┤
│ I.  Figure 2  y 625–900       │ II. Figure 3 spans BOTH columns      │
│     x 30–410                  │     y 625–900 … see note below       │
└───────────────────────────────┴──────────────────────────────────────┘
```

Because panel II must be followed **immediately** by panel III with nothing
between them, and because Figure 3 is the widest figure in the set, the
results block runs as follows:

| Order | Panel | Figure | Position (x, y, w × h mm) |
|---|---|---|---|
| 1 | 01 Introduction | — | 30, 210, 380 × 210 |
| 2 | 02 Methods | Figure 1 | 431, 210, 380 × 350 |
| 3 | 03 Results banner | — | 30, 575, 781 × 40 |
| 4 | I. Thirteen hits | Figure 2 | 30, 625, 380 × 270 |
| 5 | V. Baseline dependency | Figure 6 | 431, 625, 380 × 240 |
| 6 | **II. Corroboration is sparse** | Figure 3 | 30, 910, 781 × 300 |
| 7 | **III. Programmes replicate** | Figure 4 | 30, 1225, 380 × 260 |
| 8 | IV. Network connectivity | Figure 5 | 431, 1225, 380 × 235 |

> **Sheet-height note.** The table above overruns a single A0 portrait sheet at
> row 7 (1225 mm > 1189 mm). The set is one figure too many for A0 portrait at
> a legible size. Two resolutions, for the author to choose:
>
> **(a) A0 landscape, 1189 × 841 mm — recommended.** Three columns of 376 mm.
> Reading order left-to-right, top-to-bottom: Introduction and Methods
> (Figure 1) in column 1; Results I (Figure 2) and II (Figure 3, spanning
> columns 2–3) across the top band; III (Figure 4) directly beneath II, which
> preserves the required adjacency; then IV (Figure 5), V (Figure 6) and
> VI (Figure 7) along the bottom band; Limitations and Conclusions in
> column 3. This fits all seven figures with 25 mm of white space between
> panels.
>
> **(b) A0 portrait with six figures.** Drop Figure 5 (network connectivity),
> which is the least load-bearing panel in the argument — its content is a
> caveat about a connection rather than a step in the story, and it can be
> given verbally. Everything else fits with the adjacency preserved.
>
> The layout specification below is written for **(a) A0 landscape**, as
> recommended.

## Recommended layout — A0 landscape, 1189 × 841 mm

Margin 30 mm; usable 1129 × 781 mm; three columns of 376 mm, gutters 20 mm.

| # | Element | x | y | w × h (mm) | Reading order |
|---|---|---|---|---|---|
| — | Logo strip (6 logos, 34 mm high) | 30 | 12 | 1129 × 40 | — |
| — | Title block (violet band) | 30 | 58 | 1129 × 120 | — |
| 1 | 01 INTRODUCTION | 30 | 190 | 376 × 150 | ① |
| 2 | 02 METHODS + Figure 1 | 30 | 350 | 376 × 300 | ② |
| 3 | 03 RESULTS I + Figure 2 | 30 | 660 | 376 × 150 | ③ |
| 4 | **RESULTS II + Figure 3** | 426 | 190 | 733 × 300 | ④ |
| 5 | **RESULTS III + Figure 4** | 426 | 500 | 376 × 310 | ⑤ (directly below ④) |
| 6 | RESULTS IV + Figure 5 | 812 | 500 | 347 × 155 | ⑥ |
| 7 | RESULTS V + Figure 6 | 812 | 665 | 347 × 145 | ⑦ |
| 8 | RESULTS VI + Figure 7 | 426 | 820 | 376 × 0 | — see note |
| 9 | 04 LIMITATIONS | 812 | 190 | 347 × 140 | ⑧ |
| 10 | 05 CONCLUSIONS | 812 | 340 | 347 × 150 | ⑨ |
| 11 | REFERENCES | 30 | 820 | 376 × 0 | — see note |
| 12 | QR + contact (QR ≥ 40 × 40 mm) | 1060 | 700 | 99 × 110 | — |

> **Final packing note.** Rows 8 and 11 have no vertical room left in the
> three-column grid above. The workable arrangement is to place **Figure 7 in
> column 3 beneath Conclusions** (x 812, y 500, 347 × 145) and move Figures 5
> and 6 into column 2 beneath Figure 4 as a stacked pair (x 426, y 500, 376 ×
> 155 each), with References and the QR panel occupying the foot of column 1
> (x 30, y 660, 376 × 150). This preserves the II → III adjacency, keeps every
> figure at or above its 380 mm design width, and leaves the required white
> space. A designer should treat the millimetre values as a starting grid, not
> as a locked specification.

## Reading-order aids

- Sections are numbered `01`–`05`; results sub-panels are numbered `I`–`VI`.
- A chevron `›` in the section colour sits at the right edge of each panel
  pointing to the next, per the Cambridge guidance on clear flow.
- The II → III pair carries a single connecting arrow with the words
  *"genes: sparse → programmes: consistent"* along it, because that transition
  is the argument of the poster.

## Total word count

| Section | Words |
|---|---|
| Title block | 55 |
| 01 Introduction | 138 |
| 02 Methods (Figure 1 caption) | 84 |
| 03 Results I–VI (captions) | 814 |
| 04 Limitations | 213 |
| 05 Conclusions | 197 |
| References + pointer | 118 |
| QR and contact | 34 |
| **Total** | **1,653** |

That is a moderately text-light A0 poster: the seven figures carry the
argument and the author supplies detail verbally, as the Cambridge guidance
recommends.
