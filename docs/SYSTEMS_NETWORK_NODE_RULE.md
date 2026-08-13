# Systems-network phase 14: node-universe inclusion rule

**Precision correction (post Codex review):** the top-40 truncation for
categories B and D is a bounded-list-size decision made after observing
that an unbounded "any recurrence" rule returns >1,600 genes (documented in
category B below) -- it is not a fully blind pre-registration made without
seeing any aggregate size information. What IS true, and is the actual
guarantee here: (1) the ranking criteria used to pick the top 40 within
each category (dataset/pathway recurrence counts, now with an added
alphabetical final tie-breaker for determinism) were fixed before looking
at *which specific genes* would result, and (2) the same fixed N=40 is
applied identically to both categories rather than tuned per-category to
hit a target. The resulting 119-node total was not itself a target -- no
threshold here was adjusted after seeing the 119 figure. Target range:
50-300 nodes (task spec) was satisfied, not engineered.

## Categories (union, deterministic, deduplicated)

**A. Frozen candidates (fixed, n=4).**
USP34, VEZF1, EML5, CITED2 -- exactly the frozen therapeutic shortlist,
unchanged.

**B. Recurrent resistance leading-edge genes -- top 40.**
From `recurrent_leading_edge_genes.tsv` (Phase 9: genes in the GSEA leading
edge of a STRONG_CONSENSUS pathway in >=1 of GSE118713/GSE240112/GSE111151),
ranked by `(leading_edge_dataset_count desc, pathway_count desc)` and
truncated at the top 40. A raw "any recurrence" cut is not usable here: it
returns >1,600 genes because GO:BP contributes thousands of near-synonym
developmental terms sharing large overlapping gene memberships (documented
in Phase 9). A fixed top-N cut is a transparent, deterministic alternative
to an arbitrary significance threshold -- not a claim that gene #41 is
biologically unimportant, only that a bounded, interpretable network
requires a bounded gene list, and rank order is what is actually being
used to build it. **Codex review fix:** `recurrent_leading_edge_genes.tsv`
itself was corrected to only count a gene as leading-edge-in-a-dataset when
that dataset's own nominal p-value for the pathway is <0.05 (previously it
counted a gene from any dataset sharing a STRONG_CONSENSUS pathway label,
even a dataset where that specific pathway was not itself significant).

**C. Genome-wide CRISPR hits at FDR<0.05 -- all of them, no additional cut.**
11 `sensitising_KO` + 11 `tolerance_associated_KO` genes (n=22) out of
19,103 genome-wide tested. Already a stringent, pre-existing threshold
(the same FDR<0.05 bar used throughout every prior phase of this project);
no further truncation needed since the set is already small.

**D. Genes driving MULTIMODAL_PATHWAY pathways -- top 40.**
From `multimodal_pathway_convergence.tsv` (Phase 12: pathways with both a
reproducible RNA resistance signal and CRISPR pathway-level FDR<0.05; 31
such pathways as of the Codex-review-fixed Phase 5/12 classification logic
-- this count moved from an earlier draft's 77 when the MIXED-precedence
and per-dataset-significance fixes were applied, see git history for
`src/systems_network_multimodal.py`). Genes are the GSEA leading-edge
members of those pathways in the three resistance datasets (each gated by
that dataset's own nominal significance for the pathway, Phase 9's fix),
ranked by the number of distinct MULTIMODAL_PATHWAY pathways in which each
gene is leading-edge, truncated at the top 40 (same rationale as category B).

**E. Direct high-confidence STRING interaction partners of the candidates --
all of them, no truncation.**
`data/reference/interactions/string_candidate_partners.tsv`, score>=0.7
(STRING's own "high confidence" band), n=23 partner edges. VEZF1 and EML5
have zero partners at this threshold (Phase 13 finding, not a rule
artifact -- confirmed down to score>=0.15, STRING's lowest band, they still
have none above ~0.65).

## What is explicitly excluded

- No gene is added because its name "sounds relevant."
- No WGCNA/de-novo co-expression network (see CLAUDE.md statistical
  limitation note -- sample sizes are insufficient for that; see
  `docs/SYSTEMS_NETWORK_INPUT_AUDIT.md`).
- Category B/D's top-40 cut is reported as exactly that -- a bounded list,
  not an exhaustive "all significant" set.
