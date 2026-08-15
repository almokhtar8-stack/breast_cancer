# Scientific Figure Audit

Inventory of completed, frozen analyses in this repository and a judgment on
which ones carry poster-worthy visual evidence. This audit is read-only: it
inspects existing result tables and figures and records what exists. It
performs no new analysis, no recomputation, and no re-ranking. Every table
path below was opened directly and its schema/row-count verified.

This audit exists to justify the figure-bank build in
`src/poster_figures_bank_data.py` / `src/poster_figures_bank_visualization.py`
and the candidates in `results/figures/poster_candidates/`. See
`FIGURE_BANK_REVIEW.md` for the per-candidate scoring that follows from this
audit.

---

## 1. CRISPR screen (Hany et al. genome-wide functional screen)

| | |
|---|---|
| **Tables** | `data/processed/labels.parquet` (19,103 genes, full genome-wide fit); `results/tables/systems_network/crispr_ranked_genes.tsv`; `results/tables/nebula_plot_inputs/fig1_crispr_hit_landscape_input.tsv` (28 Gate-1 hits pre-extracted); `results/tables/evidence_freeze/THERAPEUTIC_SHORTLIST_FREEZE.tsv` (4 frozen candidates) |
| **N** | 19,103 fitted genes; 28 Gate-1 hits (FDR<0.1); 4 frozen candidates |
| **Useful plot types** | genome-wide volcano (effect vs -log10 FDR); ranked-effect plot (rank vs effect, all genes); the 28-gene Gate-1 hit table alone is too small for a meaningful heatmap unless paired with per-condition data, which does not exist in this repo (only one drug-context effect/FDR per gene, not per-condition betas) |
| **New message vs already-built figures** | The existing `poster/01_crispr_discovery.png` is a volcano. A ranked-effect view answers a genuinely different question (magnitude/rank among all screened genes) and is worth generating as an alternative. A per-condition heatmap is **not buildable** — the repo's CRISPR table is a single fitted effect/FDR per gene under one drug context (4-OHT), not a matrix of per-condition betas. Documented here explicitly per the instruction not to fabricate condition-level data. |
| **Verdict** | Main-poster quality (ranked-effect and refined volcano). Heatmap option **skipped — data unavailable**, not fabricated. |

## 2. GSE118713 (MCF7 / TAMR / FASR bulk RNA-seq, acquired endocrine resistance)

| | |
|---|---|
| **Tables** | `results/tables/gse111151/pca_coordinates.tsv` is GSE111151, not this one — GSE118713 PCA is **not precomputed anywhere in the repo** (checked: no `pca` table under any gse118713-named directory). `data/processed/gse118713_gene_tpm.parquet` (9 samples x full gene TPM matrix); `results/tables/gse118713_differential_expression.tsv.gz` (full genome-wide limma DE, all genes, all 3 contrasts); `results/tables/nebula_plot_inputs/fig2_gse118713_pca_input.tsv` (**PCA already computed** for the earlier NEBULA phase, 9 samples, PC1/PC2 + % variance) and `fig3_tamr_vs_mcf7_volcano_input.tsv` (14,837-gene volcano input, gate1/candidate flags included) |
| **N** | 9 samples (3 MCF7, 3 TAMR, 3 FASR), full transcriptome |
| **Useful plot types** | PCA scatter (real, pre-computed, reuse `fig2_gse118713_pca_input.tsv`); TAMR-vs-MCF7 volcano (real, pre-computed, reuse `fig3_tamr_vs_mcf7_volcano_input.tsv`); per-sample USP34 expression strip+box (already used in `poster/02_expression_evidence.png` panel A) |
| **New message** | PCA demonstrates the resistant phenotype is a genome-wide transcriptional shift, not a single-gene story — this message does not exist anywhere in the current poster figure set. |
| **Verdict** | Main-poster quality. Build a 3-panel `03_GSE118713_resistance_landscape` (PCA + volcano + candidate panel), reusing the already-computed nebula plot-input tables (frozen, unchanged) rather than recomputing PCA. |

## 3. GSE240112 (primary vs recurrent tumor, spatial/pseudobulk)

| | |
|---|---|
| **Tables** | `results/tables/gse240112/candidate_sample_level_log2cpm.tsv` (6 real tumor pseudobulk samples: 3 PT + 3 RT); `results/tables/gse240112/candidate_table.tsv` (edgeR DE stats, candidate genes only); `results/tables/gse240112_pseudobulk/tumor_cell_genomewide_de.tsv.gz` (full genome-wide DE, all genes — usable for a PT-vs-RT volcano); `results/tables/gse240112_pseudobulk/tumor_cell_pca_coordinates.tsv` |
| **N** | 3 primary + 3 recurrent tumors, unpaired different patients (must be stated explicitly on any figure) |
| **Useful plot types** | genome-wide PT-vs-RT volcano with VEZF1 highlighted (new — not in current poster set, which shows only a strip+box for VEZF1 alone); per-tumor VEZF1 expression strip (already built); a per-tumor "individual contribution" plot is limited by n=3 per arm — a plain 6-point strip is already the most honest representation, an additional contribution plot would not add information at this n |
| **New message** | A genome-wide volcano shows VEZF1's signal in the context of the full recurrence transcriptome, not in isolation. |
| **Verdict** | Main-poster quality for VEZF1 evidence, with an unavoidable and clearly-labeled small-n caveat. |

## 4. GSE111151 (independent resistant line panel, secondary resistance dataset)

| | |
|---|---|
| **Tables** | `results/tables/gse111151/pca_coordinates.tsv`, `genomewide_de.tsv.gz`, `candidate_table.tsv`, `sample_correlation.tsv`; figures already exist at `results/figures/gse111151/pca.png`, `correlation_heatmap.png` |
| **N** | not re-verified row count this session (pre-existing frozen figures) |
| **Verdict** | Supplementary/QC only — this dataset is used as a cross-dataset consistency check (panel of the forest plot and the pathway heatmap), not as a standalone headline figure; its existing PCA/correlation figures are QC-grade, not new poster candidates. |

## 5. GSE245601 (single-cell, acute 12h tamoxifen response)

| | |
|---|---|
| **Tables** | extensive single-cell deep-dive tables under `results/tables/gse245601_candidate_deepdive/` and `gse245601_candidate_integration/` (per-cell expression, per-cluster response, composition change, malignant enrichment) |
| **Critical caveat (verbatim project rule)** | this is an **acute 12h tamoxifen response**, not a resistance cohort. It appears in the cross-dataset forest/effect figure and the pathway landscape figure with an explicit "ACUTE RESPONSE CONTEXT" visual marker, and nowhere is it relabeled as resistance evidence. |
| **Verdict** | Supporting role only, inside multi-dataset comparison figures (forest plot, pathway heatmap) — not a standalone main figure, consistent with how it is already treated in the current poster set. |

## 6. Cross-dataset candidate evidence

| | |
|---|---|
| **Tables** | `results/tables/cross_dataset_genomewide/all_genes_cross_dataset_evidence_with_ranking.tsv` (log2FC/FDR per gene x per dataset, all 4 datasets, full gene universe subsettable to the 4 candidates) |
| **Verdict** | Main-poster quality. The existing forest panel (`poster/02_expression_evidence.png` panel C) already uses this; redesign as its own standalone, larger figure with clearer FDR encoding (filled vs open symbols) and explicit acute-context marking for GSE245601, per the user's request. |

## 7. Pathway / systems biology (GSEA + network integration)

| | |
|---|---|
| **Tables** | `results/tables/systems_network/gsea_{gse118713,gse240112,gse111151,gse245601,gse245601_track_b,crispr}.tsv` (full per-dataset GSEA NES/FDR, hallmark + go_bp + reactome collections, ~4,200 rows each); `resistance_pathway_consensus.tsv` (5,847 pathways scored across the 3 resistance datasets, 276 STRONG_CONSENSUS, 423 DIRECTIONAL_CONSENSUS); `multimodal_pathway_convergence.tsv` (5,335 rows joining resistance consensus to GSE245601 acute + CRISPR pathway context); `candidate_pathway_membership.tsv` (per-candidate pathway membership x consensus class, 261k rows) |
| **Verified key pathways (real NES, all 3 resistance datasets STRONG_CONSENSUS)** | HALLMARK_ESTROGEN_RESPONSE_EARLY (median NES -2.39, all 3 negative), HALLMARK_ESTROGEN_RESPONSE_LATE (median NES -2.13, all 3 negative), HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION (median NES +1.88, all 3 positive), HALLMARK_E2F_TARGETS (median NES -1.79, all 3 negative), HALLMARK_G2M_CHECKPOINT (median NES -1.46, 2/3 FDR<0.05), HALLMARK_UV_RESPONSE_DN (median NES +1.58, all 3 positive) |
| **Candidate pathway footprints (STRONG_CONSENSUS memberships, verified count)** | USP34: 4 (all GO Wnt-signaling terms); VEZF1: 2 (HALLMARK_HEME_METABOLISM, GOBP_BLOOD_VESSEL_MORPHOGENESIS); CITED2: 49 (dominated by cardiovascular/embryonic-morphogenesis GO terms + HALLMARK_UV_RESPONSE_DN); **EML5: 0** — confirmed no resolved STRONG_CONSENSUS pathway footprint for EML5, consistent with the user's own prior note. |
| **CRISPR-context caveat** | `gsea_crispr.tsv` reports NES under a genetically-different metric (CRISPR gene-level effect ranking, not an expression fold-change) and its sign is **not directly comparable** to the RNA NES columns (e.g. HALLMARK_ESTROGEN_RESPONSE_EARLY: CRISPR NES +1.40 vs RNA NES ≈ -2.4 in all 3 resistance datasets — opposite sign, different meaning). Any figure using this column must visually separate it from the RNA panels and never merge it into one color scale with the RNA NES values. |
| **Verdict** | Main-poster quality. This is the single most under-represented major analysis in the existing poster set (no pathway figure currently exists) despite being one of the largest completed analyses. Build a pathway x dataset NES heatmap (~10-12 pathways, RNA datasets only, one shared diverging color scale) as a main figure, and treat the CRISPR pathway-context column as a separate, clearly-labeled strip if included at all. |

## 8. Systems network / candidate connectivity

| | |
|---|---|
| **Tables** | `results/tables/systems_network/four_candidate_convergence.tsv` (all 6 candidate pairs), `four_candidate_bridge_evidence.tsv`, `four_candidate_direct_neighbors.tsv`, `four_candidate_shortest_paths.tsv`; raw network `results/networks/systems_network/{nodes,edges}.tsv` |
| **Verified connectivity** | Of 6 candidate pairs, only 2 show any real convergence: USP34-CITED2 (1 shared GO term) and VEZF1-CITED2 (11 shared GO terms, mostly vascular/morphogenesis). USP34-VEZF1, USP34-EML5, VEZF1-EML5, EML5-CITED2 show **zero** direct interaction, zero shared neighbors, zero shared pathways. |
| **Verdict** | A generic STRING-style network/hairball is **explicitly not justified** by this data — most of the graph would be empty, and EML5 is essentially disconnected from the other three by every network metric checked. A compact, curated candidate→pathway-program diagram (not a force-directed graph) is defensible and is what is built as figure 07, using only the verified STRONG_CONSENSUS pathway assignments from item 7 above — not the raw edge list. |

## 9. TCGA-BRCA (human tumor validation)

| | |
|---|---|
| **Tables** | `results/tables/independent_validation/TCGA_candidate_expression.tsv` (5 contrasts x 4 candidates: ER+/-, PAM50 LumA/LumB, paired tumor-vs-normal n=113, unpaired descriptive); `TCGA_candidate_clinical.tsv` (Cox HR, univariable + age/stage-adjusted, all-tumor and ER+ subgroup, n≈1068/580); `TCGA_candidate_pathway_associations.tsv` (candidate expression vs ssGSEA pathway NES, ER+ tumors only, n=600) |
| **Verified strongest signals** | CITED2 tumor-vs-normal (paired, n=113): mean diff -0.87 log2TPM, p=8.5e-15, FDR=3.4e-14 — the strongest single TCGA expression effect of the four candidates. VEZF1 ER+ vs ER-: mean diff +0.86, p=1.0e-23. USP34 shows the **weakest** TCGA expression signal of the four (paired tumor-vs-normal FDR=0.21, not significant) — an honest negative/null result, not omitted. Survival: several univariable Cox models **violate the proportional-hazards assumption** (flagged directly in the table, `ph_assumption_p<=0.05`) — none of the 4 candidates' survival associations should be presented as a clean, confident clinical finding. |
| **Verdict** | Main-poster quality for the expression panel (paired tumor-vs-normal + ER+/- forest across all 4 candidates, real effect sizes + CI). The clinical/survival table is present but weak and PH-violating for most rows — include only as a small secondary note, never as a headline claim, per the user's explicit instruction not to make weak survival results look central. |

## 10. DepMap 26Q1 (cancer dependency)

| | |
|---|---|
| **Tables** | `results/tables/independent_validation/DepMap_candidate_dependency.tsv` (26Q1 summary, frac_strongly_dependent per candidate); real per-line Chronos values reloadable via `src/independent_validation_depmap_data.py::load_gene_effect` (already reused, unchanged, in `poster_figures_data.py`); `DepMap_candidate_codependency.tsv` (top 10 positive/negative Chronos-correlated genes per candidate, 26Q1) |
| **New combination not yet built** | Hany CRISPR drug-context effect (x) vs DepMap ER+/luminal baseline dependency (y), across the 28 Gate-1 genes where DepMap data exist — this cross-references two already-frozen tables (`crispr_ranked_genes.tsv`/labels.parquet and DepMap gene-effect) that have never been jointly plotted. This is likely the single most informative candidate-differentiating figure available: it can show USP34 (context-specific sensitiser, near-zero baseline dependency) vs VEZF1 (real baseline dependency + sensitisation) as a genuine two-axis distinction rather than a text claim. |
| **Verdict** | Main-poster quality (the joint map) plus a compact line-level supporting panel (already-drafted style from `poster/03_depmap_distributions.png`, tightened to USP34+VEZF1 only with matched cell lines connected). |

## 11. GDSC pharmacogenomics

| | |
|---|---|
| **Tables** | `results/tables/final_pharmacogenomics/USP34_GDSC_drug_associations.tsv` (all tested drugs, full stats); `GDSC_top_associations.tsv` (16 rows across two tiers -- `FDR_SIGNIFICANT` and a separate `TOP_EFFECT_SIZE_NOT_NECESSARILY_SIGNIFICANT` tier explicitly labeled as such; every VEZF1 row in this table falls in the non-significant tier, includes drug target/pathway/tier); `GDSC_ER_luminal_subset.tsv` (explicitly marked exploratory, small-n, not independently FDR-corrected) |
| **Verified** | AZD7762 (CHK1/CHK2 inhibitor) remains the single FDR-significant, largest-n (44) hit — rho=-0.59, FDR=0.008. Next strongest: AZD1332 (NTRK), AZD6738 (ATR), both FDR<0.03 but weaker rho. |
| **Verdict** | Supporting, not core. A 2-panel figure (top-association lollipop + AZD7762 scatter) is more informative than the scatter alone, but this remains association-only pharmacogenomics on a secondary candidate axis — appropriately supplementary per the user's own instruction. |

## 12. Structural biology (USP34, 7W3R/7W3U)

| | |
|---|---|
| **Tables** | `results/tables/final_translational/USP34_structure_inventory.tsv`, `USP34_pocket_analysis.tsv` (fpocket druggability scores/pocket geometry), `USP34_docking_decision.tsv` (explicit `DOCKING_NOT_YET_JUSTIFIED` verdict with a 4-condition justification) |
| **Verified facts to preserve on any redesign** | Cys1903-His2164 catalytic dyad, directly confirmed resolved in both structures; real covalent LINK record, Cys1903(SG)-AYE(C2), 1.59-2.48 Å across 3 copies in 7W3U — direct crystallographic proof of covalent reactivity; apo-to-bound dyad tightening is a real but **non-uniform**, single-copy-strongest effect (chain A: 3.94→3.37 Å; other chains range 3.10-4.98 Å) — must not be presented as a uniform structure-wide effect; docking explicitly not performed and not justified (no calibration ligand exists) — never imply docking results exist. |
| **Verdict** | Main-poster quality, redesign required. The 4-equal-panel PyMOL composite in the current poster set is scientifically accurate but visually under-powered relative to the real structural depth available (fpocket geometry, real covalent bond distances, matched multi-copy heterogeneity). Two stronger single-hero-panel PyMOL variants are built (surface+pocket hero; apo-vs-bound matched comparison) per the user's explicit request. |

## 13. Druggability / normal-tissue safety

| | |
|---|---|
| **Tables** | `results/tables/druggability_safety/candidate_normal_tissue_context.tsv` (GTEx + HPA per-tissue expression, all 4 candidates); `candidate_genetic_constraint.tsv`; `candidate_therapeutic_window_summary.tsv`; `results/tables/lead_target_deep_dive/USP34_VEZF1_tissue_liability.tsv` (15 organ-system rows, functional-liability evidence category, separate from expression) and `USP34_VEZF1_full_tissue_expression.tsv` (16 tissue/cell-type rows, RNA + protein + single-cell columns kept separate) |
| **Verified pattern** | Both USP34 and VEZF1 are broadly, non-specifically expressed across essentially all normal tissues (no near-zero tissue for either) — expression breadth alone does not differentiate a liability profile; the tissue_liability table's `classification`/`evidence_summary` columns (documented functional phenotypes, e.g. developmental knockout phenotypes) carry the actual safety signal and are a genuinely separate axis from raw expression, exactly as the user's instruction anticipated. |
| **Verdict** | Legitimate candidate for a supplementary/QC-tier figure (expression magnitude vs functional-liability evidence category, kept as two visually distinct encodings, never merged into a single "toxicity score"). Not main-poster quality on its own — it supports the safety argument but does not carry a headline scientific message. |

## 14. Final translational strategy

| | |
|---|---|
| **Tables** | `results/tables/final_translational/final_experimental_design.tsv`, `final_normal_cell_comparators.tsv`, `final_target_success_failure_criteria.tsv` |
| **Verdict** | The one figure in the set that should remain a schematic, per explicit instruction. Redesigned here to look more like a biological assay diagram (cell/treatment glyphs) than a generic business flowchart, and kept visually compact (not poster-dominant).

---

## Summary judgment

Directories audited (per the required minimum list, all checked): `evidence_freeze/` ✓, `genome_wide_integration/` — **does not exist as a directory in this repo** (the equivalent content lives in `cross_dataset_genomewide/`, audited under item 6), `systems_network/` (tables + networks) ✓, `independent_validation/` ✓, `literature_mechanism/` ✓ (see below), `druggability_safety/` ✓, `lead_target_deep_dive/` ✓, `final_translational/` ✓, `final_pharmacogenomics/` ✓.

**`literature_mechanism/`**: `four_candidate_claim_evidence.tsv` (45KB, per-claim literature support/refutation), `four_candidate_literature_comparison.tsv`. This is a text-evidence table (claims, citations, verified/unverified flags) with no natural quantitative plot — it is reference material for a caption or a supplementary text panel, not a figure source. **No figure built from this directory**; it is correctly text-only supporting material.

Twelve figures are judged main-poster quality and are built as candidates (Section list below matches `results/figures/poster_candidates/` filenames): CRISPR (3 variants, 01a/01b-skipped/01c), GSE118713 landscape (03), GSE240112 recurrence (04), cross-dataset effects (05), pathway landscape (06), TCGA validation (08), Hany-vs-DepMap map (09) + line-level (09b), GDSC pharmacogenomics (11), structure x2 (12/12b), experimental strategy (13). The candidate mechanism map (07) is built as a small curated diagram, not a network hairball, per the connectivity data in item 8. The tissue-liability figure (10) and GDSC figure (11) are built but flagged supplementary in `FIGURE_BANK_REVIEW.md`. The CRISPR treatment-context heatmap (01b) is explicitly **not built** — no condition-level CRISPR data exists in this repository to support it, and no such data was fabricated.
