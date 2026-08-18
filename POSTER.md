# Conference poster

**From a CRISPR Screen to Therapeutic Vulnerabilities: Functional Genomics of
Tamoxifen Sensitisation in ER+ Breast Cancer**

*Almokhtar Aljarodi · KAUST Academy · Department of Genetics, University of Cambridge*

---

## What this work shows, in one paragraph

A published genome-wide CRISPR screen switched off each of 19,103 genes in turn
and found 13 whose loss makes tamoxifen work better in an
oestrogen-receptor-positive breast-cancer cell line. Four of those genes were
carried forward and checked against four independent public transcriptomic
datasets. **At gene level the evidence does not hold together**: only 2 of 16
gene-and-dataset combinations reach a false discovery rate of 0.05, each
candidate is corroborated in at most one dataset, and pooling the three
resistance-context datasets leaves nothing below a pooled false discovery rate
of 0.100. At programme level it is more consistent: oestrogen response and
cell-cycle entry are suppressed in all four datasets, and the
epithelial-to-mesenchymal transition programme is enriched in the three
long-term settings while running the other way after twelve hours of drug. The
four candidates differ in the *kind* of chemical evidence they have rather than
the amount. **This is a computational reanalysis of public data with no
laboratory work, and nothing here is a validated therapeutic target.**

---

## The figure set

All figures are in `results/figures/poster_final/` as PNG, PDF and SVG, with a
manifest recording sources, plotted values and SHA-256 per file.

| # | Figure | What it establishes |
|---|---|---|
| 1 | [Methods workflow](results/figures/poster_final/F1_methods_workflow.png) | How the work was done, and where the candidate list was fixed relative to the external data |
| 2 | [Screen and certainty](results/figures/poster_final/F2_screen_certainty.png) | Effect size and certainty are different things, and the four candidates were not chosen on effect size alone |
| 3 | [Candidate corroboration](results/figures/poster_final/F3_candidate_corroboration.png) | Two of sixteen gene-and-dataset combinations reach FDR 0.05 — and why the nulls are weak evidence rather than strong negatives |
| 4 | [Programme signal](results/figures/poster_final/F4_programme_signal.png) | Programme-level signal is present in all four datasets; the adhesion/motility programme points the other way after twelve hours |
| 5 | [Network connectivity](results/figures/poster_final/F5_network_connectivity.png) | The two connected candidates are joined by one route, not four |
| 6 | [Baseline dependency](results/figures/poster_final/F6_baseline_dependency.png) | Needing a gene at baseline and being sensitised by losing it are different measurements that need not coincide |
| 7 | [Chemical reachability](results/figures/poster_final/F7_reachability.png) | The four differ in the kind of chemical evidence they have, not in how much of one thing |

The **six frozen poster figures** from the earlier round are preserved
unchanged in [`poster/final_figures/`](poster/final_figures/) with their own
manifest. Nothing in this set overwrites them.

---

## Documents

| Document | What it is |
|---|---|
| [`POSTER_TEXT.md`](results/reports/poster_final/POSTER_TEXT.md) | Every word intended for the printed sheet, with word counts and a layout map |
| [`SCIENTIFIC_REVIEW.md`](results/reports/poster_final/SCIENTIFIC_REVIEW.md) | Independent validity review: every printed number re-derived from source, prohibited claims checked, disclosures checked |
| [`DEFENCE_GUIDE.md`](results/reports/poster_final/DEFENCE_GUIDE.md) | Rehearsal notes: the ninety-second opening, fifteen likely questions with answers, the three weakest points, and what to say when the answer is not known |
| [`FIGURE_SET_NOTE.md`](results/reports/poster_final/FIGURE_SET_NOTE.md) | What changed per figure and why |
| [`CODEX_REVIEWS.md`](results/reports/poster_final/CODEX_REVIEWS.md) | The three independent review checkpoints and what each changed |

---

## What it does not claim

- No candidate is a validated therapeutic target.
- No candidate is shown to cause or drive resistance.
- Nothing here says that inhibiting any of these genes would benefit a patient.
- Nothing is validated, confirmed or replicated at gene level.
- The twelve-hour dataset measures acute drug response, **not** resistance.
- The primary-versus-recurrent comparison is **unpaired** — different patients,
  two different tissue banks, with treatment group confounded by bank.
- Structural evidence indicates chemical reachability only, never efficacy.
- **No docking, binding prediction or molecular modelling was performed.**

---

## Reproducing the figures

```bash
# rebuild all seven, with manifest, verification table and colour-vision simulations
python scripts/poster/final_build_all.py

# rebuild one
python scripts/poster/final_3_corroboration.py

# the tests that guard them
python -m pytest tests/test_poster_final_palette.py tests/test_poster_final_figures.py
```

Every figure verifies its plotted values against the frozen source tables
before rendering and **fails loudly rather than substituting**. PNG and PDF
output is byte-reproducible.

---

## Assets still needed before printing

- **Logos** — none are in the repository. See [`assets/logos/README.md`](assets/logos/README.md)
  for the six required, which are missing, and the rules for using them.
- **QR code** — not supplied and not generated (no QR library is available and
  none was installed). See [`assets/qr/README.md`](assets/qr/README.md) for the
  exact specification.
