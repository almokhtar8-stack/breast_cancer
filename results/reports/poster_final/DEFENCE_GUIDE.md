---
title: Defence guide — rehearsal notes for presenting the poster
status: post_freeze_exploratory
analysis_date: 2026-08-18
branch: poster-final (unmerged)
---

# Defence guide

Rehearse the opening until it is fluent, skim the fifteen questions, and read
the three weak points last so they are fresh. The single most useful habit at
a poster is this: **when you do not know, say so and say what would settle
it.** A stated gap costs nothing. A guess that a specialist catches costs the
whole conversation.

---

## The ninety-second opening (verbatim)

> A published genome-wide screen switched off every gene in turn, one at a
> time, and asked which of them makes tamoxifen work better in an
> oestrogen-receptor-positive breast-cancer cell line. Thirteen genes passed a
> threshold that had been written down before anyone looked at the results.
>
> The interesting part is what happened next. I took four independent public
> transcriptomic datasets and asked whether those genes also *change* in
> resistance. Four genes across four datasets is sixteen combinations, and
> two of them reach significance. Two out of sixteen. At the level of
> individual genes, this evidence does not hold together.
>
> But at the level of biological programmes it does. Oestrogen response is
> suppressed in all four datasets, cell-cycle entry in all four — and the
> programme controlling cell adhesion and motility runs one way in the three
> long-term resistance settings and the opposite way after twelve hours of
> drug. That dissociation is the observation this project adds.
>
> So my conclusion is deliberately modest. No gene here is a validated target.
> What I have is a hypothesis with a testable shape, and one candidate —
> USP34 — that is worth carrying not because its screen evidence is strongest,
> it is twelfth of thirteen, but because it is the one with an unexplored,
> chemically reachable catalytic site. It is a lead on tractability, not on
> evidence strength.
>
> Everything here is a computational reanalysis of public data. I have not run
> a single experiment, and I would like to.

*(Roughly 240 words; about 90 seconds at a natural pace. If the listener is
clearly a specialist, drop the last paragraph and go straight to Figure 3.)*

---

## One sentence per figure

| Figure | Say this |
|---|---|
| **1 · Methods** | "Everything on the left was decided before I opened the data on the right, and the dashed box is what I added afterwards — I am not claiming the whole thing was pre-planned." |
| **2 · Screen** | "Effect size is the horizontal axis and certainty is the vertical one, and they are not the same thing — the two most certain hits are not the two biggest." |
| **3 · Corroboration** | "Sixteen combinations, two filled dots — and the grey strip underneath is why the empty ones are weak evidence rather than proof of absence." |
| **4 · Programmes** | "The programmes behave consistently across all four datasets, and the boxed row runs one way in long-term resistance and the other way after twelve hours." |
| **5 · Network** | "Two of the candidates are connected, but all four supposedly independent routes go through the same protein and then a ubiquitin gene — so it is one connection, not four." |
| **6 · Dependency** | "TLK2 is needed by nine of eleven cell lines with or without the drug, which is a narrower therapeutic window, not a better one." |
| **7 · Reachability** | "These four differ in the kind of chemical evidence they have, not the amount — and one of them has no experimental structure at all, which is itself the finding." |

---

## The fifteen most likely questions

Ordered by how likely they are to be asked.

### 1. "Has any of this been tested experimentally?"
> No. Not one bench experiment, by me or for this project. Everything is a
> reanalysis of published data. That is the single biggest limitation and it is
> on the poster. The experiments I would run first are a tamoxifen-sensitisation
> assay for USP34 and KDM1A in two or three independent
> oestrogen-receptor-positive lines.

### 2. "Is anything actually replicated?"
> Not at gene level, no — and I say so on the poster. Two of sixteen
> gene-and-dataset combinations reach significance, each candidate in at most
> one dataset, and pooling the three resistance datasets gives nothing below a
> false discovery rate of 0.100. What *is* consistent is the programme level:
> oestrogen response and cell-cycle suppression appear in all four datasets. I
> would not call that replication of the gene findings; it is a different and
> weaker claim, and I try to keep them separate.

### 3. "Why isn't KDM1A your lead? It is the strongest hit in your own screen."
> It is, on both effect and certainty — and if the question is "which gene is
> the best functional sensitiser", the answer is KDM1A. Two things move it out
> of the lead position for *this* project. It is null in every transcriptomic
> dataset I have. And it is already a mature drug target: there are
> inhibitor-bound structures and nine clinical-stage LSD1 inhibitors, so a
> reanalysis like mine adds very little to it. USP34 is the opposite — weaker
> screen evidence, but genuinely unexplored chemically. If your interest is
> the biology rather than the novelty, KDM1A is the better gene and I would say
> so.

### 4. "How were these four candidates chosen, when two genes had larger effects?"
> Honestly, by two different rules, and the poster marks which is which. USP34
> and VEZF1 came through a frozen rule that required a screen hit plus at least
> one corroborating dataset. KDM1A and TLK2 were added later, after an external
> reviewer argued that requiring transcriptional corroboration was wrong — a
> gene does not have to change its expression for its loss to sensitise cells.
> They are the two strongest screen hits and they have no corroboration at all.
> The genes with larger effects than USP34, like USP17L29 and TADA2B, failed
> the corroboration rule and were not rescued by the audit. I would not present
> these four as a single uniformly selected panel, because they are not one.

### 5. "How were the thresholds chosen, and when?"
> Before any results existed. The analysis plan is a dated document in the
> repository with the thresholds written into it — false discovery rate below
> 0.10 with a sensitising direction for the screen, 0.05 for the transcriptomic
> comparisons. There is a git tag marking the freeze, and the candidate list
> was fixed before I opened the external validation datasets. You can check
> both: the tag and the plan are in the repository behind the code opposite.

### 6. "Your recurrence dataset compares different patients. What does that actually show?"
> Association, not causation, and I would not claim more. It compares primary
> against recurrent tumours from *different, unpaired* patients drawn from two
> different tissue banks — so treatment group and tissue bank are confounded,
> and I cannot separate a recurrence effect from a bank effect. It is the only
> human-tumour signal any of my candidates has, which is why it is on the
> poster, but it is the weakest kind of evidence there.

### 7. "Why is TLK2 required by most of your cell lines? Doesn't that kill it?"
> It weakens it, and I present it as a limitation rather than a strength. TLK2
> is needed for survival by nine of eleven oestrogen-receptor-positive lines
> whether or not tamoxifen is present. A gene most cells already need offers a
> narrower, less tamoxifen-specific therapeutic window — you would be closer to
> a general cytotoxic than to a sensitiser. It is worth saying that dependency
> measured in cancer cell lines is not the same as normal-tissue toxicity; I
> have no safety data at all.

### 8. "Why do your malignant-cell calls differ from the source publication?"
> Because I reconstructed them and they did not publish theirs per cell. I
> inferred copy-number states to separate malignant from non-malignant cells,
> using my own thresholds; the authors' per-cell labels are not in the public
> record — I checked the GEO record, all seventy-six supplementary worksheets
> and the authors' repository, and the processed data sits behind controlled
> access. So the difference cannot be quantified, only acknowledged. I tried to
> quantify it and stopped when I established the labels were not available,
> rather than inventing a proxy.

### 9. "Two of sixteen — isn't that just noise?"
> It is consistent with noise, and I do not argue otherwise. What I argue is
> narrower: the nulls are weak evidence rather than strong negatives. The
> best-designed of the resistance datasets could only have detected a 1.5-fold
> change at eighty percent power, and every candidate effect it observed was
> smaller than that — so it could not have found what it was looking for. That
> does not make the candidates real. It means this design cannot settle it.

### 10. "Isn't the pathway result just a generic tamoxifen response?"
> Partly, yes, and that is the point of showing it. Oestrogen-response and
> cell-cycle suppression are exactly what you would expect if the datasets are
> behaving — I use them as an internal control. The part that is not generic is
> the adhesion-and-motility row, which is positive in the three long-term
> settings and negative after twelve hours. I would add that broad, correlated
> gene sets can reach significance when individual gene effects are unstable,
> so I do not treat programme-level consistency as stronger evidence than
> gene-level consistency — just different.

### 11. "What is the mechanism?"
> I have a hypothesis, not a mechanism. The honest position is that a network
> association is not a mechanism, and my Figure 5 is there mostly to show how
> weak the connection is: the four supposedly independent shortest routes
> between KDM1A and USP34 all pass through the same protein and then a
> ubiquitin gene. There is also published counter-evidence I should mention —
> losing USP34 has been reported to push breast cells toward a *more* mobile
> state, which runs against the direction I would want.

### 12. "Did you do any docking or modelling?"
> No. None. No docking, no pose prediction, no affinity estimation, no homology
> modelling. Figure 7 reports what experimental structures exist and what is
> bound in them, nothing more. USP34's structure holds a covalent ubiquitin
> activity probe — a laboratory tool, not a drug — over about twelve percent of
> the protein.

### 13. "Why did you use a twelve-hour dataset at all if it isn't resistance?"
> Because it answers a different and useful question — what the immediate drug
> response looks like — and because leaving it out would have hidden a result
> rather than clarified one. It is never counted as resistance evidence; it is
> excluded from the resistance consensus by a rule written down in advance. It
> is also the dataset that produces the adhesion-and-motility contrast, which
> only exists because acute and long-term are being compared.

### 14. "Your screen is one screen in one cell line."
> Correct, and that bounds everything downstream. It is MCF7-derived, a single
> genome-wide screen, and it asks whether losing a gene makes tamoxifen work
> better in a drug-tolerant parental line — not whether it reverses established
> resistance. A second screen in a different oestrogen-receptor-positive
> background is on my list of things I cannot currently do.

### 15. "So what would you actually do with this?"
> Nothing clinical, and nothing that costs a lot. I would run a sensitisation
> assay for USP34 and KDM1A in independent lines, because that is cheap and
> decisive. If USP34 replicates functionally, the covalent catalytic cysteine
> makes it a reasonable chemical-biology starting point. If it does not, I
> would say so publicly — the repository is set up so that the negative result
> is as visible as the positive one.

---

## The three weakest points, and how to disclose each before being asked

**1. No experimental validation of any kind.**
Disclose it in the opening, in your own words, before anyone asks. Something
like: *"I should say up front that none of this has been tested at the bench."*
Said first, it reads as rigour; extracted under questioning, it reads as
something you were hiding.

**2. The four candidates come from two different selection rules.**
Bring it up when you point at Figure 2, using the marker shapes as the excuse:
*"The circles and squares are two different selection rules — these four are
not one uniformly chosen panel, and I mark that rather than smoothing it
over."* A specialist who spots this unaided will assume the worst about
everything else on the poster.

**3. The lead is the weakest candidate on the screen evidence.**
Say it in the same breath as naming USP34, every time: *"twelfth of thirteen by
effect size, and it leads on tractability rather than evidence strength."*
Never let USP34 be introduced as "my top hit"; it is not, and the correction is
much more damaging if someone else makes it.

---

## When you do not know

Use a fixed formula so you are not improvising under pressure:

> "I don't know. What would settle it is ___, and I haven't done that."

Or, if you genuinely have no route to an answer:

> "That's outside what this analysis can tell you. I'd be guessing, and I'd
> rather not."

Do not fill silence with speculation about mechanism, clinical relevance, or
patient benefit. Those are exactly the areas where an overclaim is both easiest
to make and most costly in a therapeutic field.

---

## Questions you genuinely cannot answer from the current evidence

Flagged deliberately. If any of these is asked, the honest answer is a version
of "I can't answer that from what I have."

1. **Does any of this hold in a patient?** No human outcome data was analysed
   at all. There is no route from this evidence to a clinical claim.
2. **Is any candidate causal for resistance?** Nothing in the design can
   establish causation. The screen shows a functional sensitisation effect in
   one cell line; the transcriptomic work shows association at best.
3. **Would inhibiting USP34 be safe?** Unknown and un-analysed. Cell-line
   dependency is not toxicity, and no normal-tissue or safety data was used.
4. **How much do our malignant-cell calls differ from the authors'?**
   Unquantifiable — their per-cell labels are not public. This was checked
   exhaustively and stopped there rather than approximated.
5. **Which of the thirteen hits is the "real" one?** The evidence does not
   separate them. USP34 is carried furthest on tractability, which is a
   different criterion from being right.
6. **Would a bigger version of the same analysis settle it?** Probably not.
   The limiting factor is that the available resistance datasets are
   underpowered and heterogeneous in design, not that they were analysed badly.
