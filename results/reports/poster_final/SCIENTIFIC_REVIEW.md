---
title: Scientific validity review of the finished poster
status: post_freeze_exploratory
analysis_date: 2026-08-18
branch: poster-final (unmerged)
reviewer: this analysis, plus an independent Codex pass (checkpoint 2)
---

# Scientific validity review

A review of what the finished figures and text claim, checked against the
frozen tables. No new analysis was run for this review.

Two passes were made: an independent Codex review with direct filesystem
access (logged in `CODEX_REVIEWS.md`), and my own programmatic re-derivation
of every printed number from its source loader. Both are reported; where they
disagreed, the disagreement is recorded rather than resolved silently.

---

## 1. Numbers printed on a figure or in the text

**Result: no mismatches found.** Every number was re-derived from its frozen
source and compared to the printed value.

| Claim as printed | Source | Re-derived | Verdict |
|---|---|---|---|
| 13 hits | `load_significant_sensitising_hits()` | 13 | exact |
| 19,103 genes fitted | `load_genomewide_crispr()` | 19,103 | exact |
| USP34 twelfth of thirteen by effect size | frozen CRISPR effects, ascending | rank 12 | exact |
| 2 of 16 combinations at FDR 0.05 | four frozen DE tables | 2 (USP34/GSE118713, VEZF1/GSE240112) | exact |
| USP34 FDR 0.0073 (GSE118713) | `evidence_long.tsv` | 0.007314 | exact to 4 dp |
| VEZF1 FDR 0.0195 (GSE240112) | `evidence_long.tsv` | 0.019490 | exact to 4 dp |
| Smallest pooled FDR 0.100 | `candidates_meta_analysis.tsv`, `all3`/`wald` | 0.100257 | exact to 3 dp |
| No candidate reaches pooled FDR 0.05 | same | 0 of 13 | exact |
| GSE111151 median 80%-power threshold ≈ 1.5-fold | `dataset_sensitivity_summary.tsv` | 1.544-fold | exact to 1 dp |
| Per-gene range 1.34–1.74 fold | same | 1.34–1.74 | exact |
| TLK2 9 of 11, VEZF1 3 of 11, KDM1A 0 of 11, USP34 0 of 11 | DepMap loader | 9 / 3 / 0 / 0 | exact |
| 47 proteins, 147 associations, 3 components | pinned STRING query | 47 / 147 / 3 | exact |
| Four shortest routes, all via DNMT1 | recomputed with networkx | 4, all via DNMT1 | exact |
| USP34 catalytic domain ≈ 12% of the protein | structural audit table | "~12% of the 3546-aa protein" | exact |

Each figure additionally re-checks its own values at render time and raises
rather than plotting on a mismatch: **87 checks across the seven figures, all
passing** (`results/figures/poster_final/verification_against_frozen.tsv`).

### One number was corrected during review

Codex flagged the original wording *"could only have detected a 1.5-fold
change at 80% power"* as imprecise: 1.5-fold is the **median** per-gene
threshold, and the per-gene range is 1.34- to 1.74-fold. The text now reads
"had a median 80%-power detection threshold of about 1.5-fold (per-gene range
1.34- to 1.74-fold)". This is a real tightening, not a cosmetic one — the
original phrasing implied a single universal detection floor.

---

## 2. Prohibited claims

**Result: none found.** Every prohibited claim was searched for in the poster
text and in the figure captions, and none appears as an affirmative claim. The
poster in fact states the negation of most of them explicitly:

| Prohibited | Where the poster states the opposite |
|---|---|
| A validated target | "Nothing here is a validated therapeutic target." (Limitations) |
| Causes / drives resistance | "does not show tamoxifen causality" (Limitations); "association, not causation" (Q6 of the defence guide) |
| Patient benefit from inhibition | "No claim is made that inhibiting any of these genes would benefit a patient." (Limitations) |
| Validated / confirmed / replicated at gene level | "At the level of individual genes, this evidence does not hold together" (opening); "corroboration", never "replication", at gene level |
| Acute dataset measures resistance | "an acute twelve-hour exposure that measures immediate drug response, **not resistance**" (Results II); "not resistance" on Figure 3 itself |
| Acute dataset is single-cell in the analysis | described as pseudobulk throughout; Figure 3's panel subtitle reads "not resistance", and the analysis unit is stated as per-tumour |
| Primary-vs-recurrent is paired or matched | "**different, unpaired patients in two separate tissue banks**" (Results II, Limitations) |
| Structural evidence indicates efficacy | "This is reachability only — not efficacy, not selectivity, not safety." (Results VI) |
| Docking or modelling performed | "**No docking, binding prediction or molecular modelling was performed anywhere in this project.**" (Results VI) |

### Five near-misses were found and four were fixed

Codex identified phrasings a hurried or hostile reader could over-read. Four
were changed:

| Was | Now | Why |
|---|---|---|
| "the external **validation** datasets were opened" | "the external **corroboration** datasets were opened" | These datasets do not validate the genes; the poster's own later wording is "corroboration", and the two should not differ |
| "works against the **benefit proposed here**" | "works against the **hypothesised sensitisation benefit**" | "benefit proposed here" is the closest the poster came to a benefit claim |
| "Counter-evidence to the **mechanism** proposed here" (reference 5) | "Counter-evidence to the **hypothesis** proposed here" | The body text says the observation is *not* a mechanism; the reference annotation contradicted it |
| "an unexplored, **chemically reachable** catalytic site" | "an unexplored, **experimentally addressable** catalytic cysteine" | 7W3U demonstrates covalent engagement by a ubiquitin activity probe at a large groove-shaped surface — not demonstrated small-molecule druggability |

The fifth, "KDM1A is the mature target", was kept. It is true, it refers to
pre-existing pharmacology rather than to anything this project validated, and
the surrounding sentences make that unambiguous.

### A sixth correction: "pre-registered" was overstated

Codex challenged "pre-registered gate" as technically inaccurate: the project
has a dated analysis plan under version control and a git freeze tag, which is
**not** the same as a public preregistration in a registry. This was correct
and the term has been changed to **"pre-specified"** everywhere it appeared —
in the poster text, on Figure 1 and on Figure 2 — with an added limitation
line stating plainly that the plan was dated before results existed but was
not lodged in a public registry. This mattered: "pre-registered" is a term of
art in this field and using it loosely is the kind of thing a methodologist
notices immediately.

---

## 3. Required disclosures

**Result: all six present.**

| Required | Where | Strength |
|---|---|---|
| Computational reanalysis of public data, no laboratory work | Introduction; Figure 1 caption; Limitations; defence-guide opening | Stated four times, plainly |
| Two of sixteen combinations reach significance | Results II headline (the largest text in the results block); opening paragraph | Prominent |
| Lead ranks twelfth of thirteen by effect size **and** leads on tractability | Results I ("sits twelfth of thirteen by effect size"); Conclusions ("a lead on tractability, not on evidence strength") | Both present, in separate sections |
| One candidate required by most cell lines regardless of drug, as a limitation | Results V: "**That is a limitation, not a strength**" | Unambiguous |
| Best-designed resistance dataset returned nulls; power shows underpowered rather than negative | Results II, third paragraph; Figure 3's pooled-evidence strip | Quantified |
| Analysis plan and thresholds predate results; candidate list fixed before external data opened | Figure 1 (drawn on the diagram); Methods caption; Limitations | Prominent, and now correctly qualified |

---

## 4. Significance encoding

**Consistent across all figures.** Colour identifies the gene and never its
rank or its significance; passing a threshold is encoded by **fill** and
failing by a **hollow ring in the same colour**. This holds on Figure 2
(markers), Figure 3 (volcano points), Figure 4 (filled cell versus empty cell
marked "n.s.") and Figure 7 (evidence track). It is enforced by
`significance_marker()` in `src/poster_final_common.py` and by
`tests/test_poster_final_figures.py::test_significance_encoding_is_fill_not_colour`.

Palette isolation is likewise test-enforced: no hex literal appears outside
`src/poster_palette.py`, no poster-chrome colour reaches a figure, and no
colour is imported from a v1 or v2 renderer
(`tests/test_poster_final_palette.py`, 31 tests).

---

## 5. References

Five references are printed. Four were verified against the repository's own
provenance records; one could not be fully verified and is flagged rather than
presented as certain.

| # | Reference | Status |
|---|---|---|
| 1 | Hany D. *et al.* (2023) *Science Advances* 9:eadd3685 | **Recovered and verified.** The brief listed this as incomplete in existing documentation; it is recorded in the root `README.md` (lines 207, 237, 373) and `docs/DATA_PROVENANCE.md`, and is the source of `data/processed/labels.parquet`. **The article title as printed is a reconstruction and should be confirmed against the journal before printing** — the repository records the journal, volume, article number and the Data S1 supplement it uses, but not the full title string. |
| 2 | Kim H., Whitman A.A. *et al.* (2023) *Clinical Cancer Research* 29(23):4894–4907, PMID 37747807 | **Recovered and verified.** The brief listed this as incomplete; it is recorded verbatim in `config/config.yaml` (lines 186–189) alongside the GSE245601 provenance. Author list, journal, volume, pages and PMID all come from that record. Title as printed is likewise a reconstruction and should be confirmed. |
| 3 | Liberzon A. *et al.* (2015) *Cell Systems* 1(6):417–425 (Hallmark collection) | Standard citation for the gene-set collection used; consistent with `data/reference/genesets/hallmark.gmt`. Not independently re-verified against the journal in this pass. |
| 4 | Szklarczyk D. *et al.* (2023) *Nucleic Acids Research* 51(D1):D638–D646 (STRING) | Standard citation for the interaction database used. Not independently re-verified against the journal in this pass. |
| 5 | Luo H. *et al.* (2017) *Cellular Signalling* 36:1–10, PMID 28499884 | **PMID supplied by the brief and used as given.** The repository's structural audit cites this line of counter-evidence; the PMID, journal and year come from the brief. Author and page range should be confirmed before printing. |

**Recommendation:** references 1, 2 and 5 carry reconstructed title/author
strings. The identifiers (journal, volume, article number, PMID) are solid;
the human-readable strings around them should be checked against the journals
in five minutes before the poster goes to print. Nothing has been fabricated —
where the repository did not record a string, that is said here rather than
invented.

---

## 6. Is the narrative a fair representation?

**Yes.** Both this review and the independent Codex pass reached the same
conclusion.

A well-informed sceptic reading only this poster would come away knowing: that
it is one screen in one cell line; that the four candidates come from two
different selection rules and not one; that gene-level corroboration is 2 of
16; that pooled evidence is not significant; that the nulls are weak rather
than decisive, and *why* they are weak differs by dataset; that the acute
dataset is not resistance; that the recurrence comparison is unpaired and
confounded with tissue bank; that baseline dependency is a limitation rather
than a bonus; that a network path is not a mechanism; that structural evidence
is reachability and not efficacy; that no laboratory work was done; and that
published counter-evidence exists for the lead candidate.

The poster's own summary of its evidential position — *a hypothesis worth
testing, not a supported target* — is the correct one for this evidence.

### Where the poster is still vulnerable

Three sentences remain the most attackable, and the defence guide prepares
answers for all three:

1. **"Its nulls are therefore weak evidence, not strong negatives."** A
   statistician can argue that a null still updates belief downward even when
   power is limited, and that a post-hoc minimum-detectable-effect calculation
   does not establish that a result carries *no* negative information. The
   wording was softened from the original categorical "uninformative, not
   negative" in response to this, but the claim is still a claim.
2. **"a lead on tractability, not on evidence strength."** A covalent ubiquitin
   activity probe demonstrates catalytic-cysteine reactivity, not small-molecule
   druggability, and the pocket is a large groove rather than a compact drug
   site. See the standing disagreement below.
3. **"Programme-level signal is present in all four datasets."** Broad,
   correlated gene sets can reach significance when individual gene effects are
   unstable, so programme-level consistency is not automatically stronger
   evidence than gene-level inconsistency — merely different. The poster says
   this, but a reader skimming the headline may not absorb it.

---

## 7. Standing disagreement, recorded rather than resolved

The project brief lists as a **required** disclosure that "the lead candidate
ranks twelfth of thirteen by effect size in the screen **and leads on
tractability**". At both review checkpoints, Codex judged "leads on
tractability" to be stronger than the structural evidence supports, on the
grounds that USP34's evidence is a covalent activity **probe** in an isolated
domain, not demonstrated druggability, and that KDM1A is by a wide margin the
more chemically mature target.

I have kept the sentence, because the brief requires it and the author
reaffirmed the brief. It has been calibrated as far as it can be without
removing it: the claim is paired in the same paragraph with the statement that
USP34's screen evidence is the weakest of the four, the probe is never called
a drug or an inhibitor, "chemically reachable" was replaced with
"experimentally addressable", and Figure 7 shows KDM1A holding every level of
the evidence track that USP34 does not.

**A supervisor should be told that this specific phrase is the one an
independent reviewer pushed back on twice.** If the author is willing to drop
the brief's wording, the defensible replacement is: *"USP34 is the candidate
this work carries furthest, on the novelty of an unexplored catalytic site
rather than on the strength of its evidence."*

---

## 8. Summary

| Check | Result |
|---|---|
| Numbers match their source tables | **No discrepancies found** (15 checks here, 87 at render time) |
| Prohibited claims | **None found**; 4 near-misses reworded, 1 term of art ("pre-registered") corrected |
| Required disclosures | **All six present** |
| Significance encoding consistent | **Yes**, and test-enforced |
| References real and verifiable | 2 recovered as the brief asked; **3 carry reconstructed title/author strings flagged for a five-minute check**; none fabricated |
| Narrative fair to a sceptic | **Yes** |
| Outstanding disagreement | One, recorded in §7, for the supervisor to arbitrate |
