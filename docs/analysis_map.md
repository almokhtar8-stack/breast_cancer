# Analysis map

A conceptual index for navigating `src/` and `scripts/`. For each poster layer:
what goes in, which module implements it, what canonical output it produces, and
which tests cover it. This exists so a newcomer does not have to read ~100
modules to find the six things that matter for the poster.

For the narrative phase-by-phase history see
[`PROJECT_WORKFLOW.md`](PROJECT_WORKFLOW.md); for full module inventory see
[`CODE_MAP.md`](CODE_MAP.md).

---

## 01 · CRISPR discovery — *frozen*

| | |
|---|---|
| **Input** | `data/processed/labels.parquet` (Hany et al. genome-scale screen; 19,103 fitted genes; E2+4-OHT vs E2) |
| **Data loader** | `src/post_audit_sensitivity_data.py` → `load_genomewide_crispr()`, `load_significant_sensitising_hits()` (FDR < 0.10, effect < 0 → 13 hits) |
| **Figure module** | `src/poster_crispr_discovery_v1.py` |
| **Wrapper** | `scripts/poster/01_crispr_discovery.py` |
| **Canonical output** | `results/figures/poster_crispr_discovery_v1/CRISPR_discovery_main.*` → `poster/final_figures/01_crispr_discovery.*` |
| **Tests** | `tests/test_poster_crispr_discovery_v1.py`, `tests/test_post_audit_sensitivity.py` |

## 02 · Candidate expression — *derived from frozen*

| | |
|---|---|
| **Input** | Processed expression for GSE118713, GSE111151, GSE240112, GSE245601 (29 biological rows total) |
| **Data loader** | `src/poster_exploration_v2_data.py`; row assembly in `src/poster_hero_heatmap_v4.py` (`_build_gse*_rows`) |
| **Figure module** | `src/poster_hero_heatmap_v6.py` (presentation layer over v5 → v4) |
| **Wrapper** | `scripts/poster/02_candidate_expression.py` |
| **Canonical output** | `results/figures/poster_hero_heatmap_v6/HERO_sample_level_heatmap_v6.*` → `poster/final_figures/02_candidate_expression.*` |
| **Tests** | `tests/test_poster_hero_heatmap_v4.py`, `..._v5.py`, `..._v6.py` |

Colour is a gene-wise z-score computed **within each dataset block**; the figure
never implies cross-dataset absolute expression. GSE245601 rows are per-tumour
pseudobulk (3 patients × Control/Tamoxifen), not individual cells.

## 03 · Molecular networks — **post-freeze exploratory**

| | |
|---|---|
| **Input** | `data/reference/interactions/string_v2_level{1,2}_{functional,physical}.tsv` (STRING, species 9606, required_score ≥ 0.7) |
| **Download script** | `scripts/download_string_network_v2_four_focus.py` (run manually; never called at analysis time) |
| **Graph build** | `src/poster_network_mechanism_v2.py` → `build_network()` (47 nodes, 147 edges, 3 components) |
| **Figure module** | `src/poster_network_mechanism_v4.py` (five-panel presentation of the v2 graph) |
| **Wrapper** | `scripts/poster/03_molecular_networks.py` |
| **Canonical output** | `results/figures/poster_network_mechanism_v4/NETWORK_mechanism_v4.*` → `poster/final_figures/03_molecular_networks.*` |
| **Cytoscape export** | `src/cytoscape_v4_export.py` → `results/tables/cytoscape_v4_network_{nodes,edges}.tsv`, `cytoscape_v4_candidate_shortest_paths.tsv` |
| **Tests** | `tests/test_poster_network_mechanism_v2.py`, `..._v3.py`, `..._v4.py`, `tests/test_cytoscape_v4_export.py` |

Superseded development versions `v1`/`v2`/`v3` are retained for provenance;
`v4` is canonical. STRING edges are undirected functional associations.

## 04 · Pathway remodeling — *derived from frozen*

| | |
|---|---|
| **Input** | `results/tables/systems_network/gsea_{gse118713,gse111151,gse240112,gse245601}.tsv` (frozen GSEA output) |
| **Data loader** | `src/poster_exploration_v2_data.py` → `load_pathway_trajectories()` |
| **Figure module** | `src/poster_pathway_v2.py` (v1 retained for provenance) |
| **Wrapper** | `scripts/poster/04_pathway_remodeling.py` |
| **Canonical output** | `results/figures/poster_pathway_v2/PATHWAY_v2.*` → `poster/final_figures/04_pathway_remodeling.*` |
| **Tests** | `tests/test_poster_pathway_v1.py`, `tests/test_poster_pathway_v2.py` |

Pathways were selected theme-first from the network layer (estrogen response,
EMT, WNT/β-catenin, E2F) before inspecting significance; WNT is retained despite
being weak/mixed precisely because it was pre-specified.

## 05 · Baseline dependency — *derived from frozen*

| | |
|---|---|
| **Input** | DepMap **26Q1** `CRISPRGeneEffect.csv`, `CRISPRGeneDependency.csv`, `Model.csv` (external, see `data/README.md`) + frozen CRISPR effects |
| **Frozen rules** | ER+/luminal subset via `src/independent_validation_depmap_data.py` → `load_model()`; strong-dependency threshold `config.independent_validation.depmap.strong_dependency_probability_threshold` = 0.5 |
| **Cell-line extract** | `src/poster_depmap_v1.py` → `load_cellline_table()` (cached to `results/tables/poster_depmap_v1/`, 11 evaluable lines) |
| **Figure module** | `src/poster_depmap_v2.py` (v1 heatmap retained for provenance/detail) |
| **Wrapper** | `scripts/poster/05_depmap_dependency.py` |
| **Canonical output** | `results/figures/poster_depmap_v2/DEPMAP_v2.*` → `poster/final_figures/05_depmap_dependency.*` |
| **Tests** | `tests/test_poster_depmap_v1.py`, `tests/test_poster_depmap_v2.py`, `tests/test_independent_validation.py` |

## 06 · Structural tractability — *derived from frozen*

| | |
|---|---|
| **Input** | `results/tables/post_audit_sensitivity/06b_structural_tractability_audit.tsv` (audited evidence) + local PDB files 6NQU, 5O0Y, 7W3U |
| **Download scripts** | `scripts/download/download_usp34_structures.py`, `scripts/download/download_kdm1a_tlk2_structures.py` |
| **Structure renders** | `scripts/render_druggability_structures.py` (PyMOL; static renders only — no docking, no pose prediction) |
| **Figure module** | `src/poster_druggability_v1.py` |
| **Wrapper** | `scripts/poster/06_structural_tractability.py` |
| **Canonical output** | `results/figures/poster_druggability_v1/DRUGGABILITY_v1.*` → `poster/final_figures/06_structural_tractability.*` |
| **Tests** | `tests/test_poster_druggability_v1.py`, `tests/test_post_audit_sensitivity.py` |

PDB IDs displayed by the figure are recovered from the audit table's own text at
build time, so no structure ID is hand-asserted. VEZF1 has no experimental
structure and is deliberately rendered without one.

---

## Supporting / exploratory layers not in the poster

| Layer | Modules | Why not in the poster |
|---|---|---|
| GDSC pharmacogenomics | `src/final_pharmacogenomics_*.py` | Exploratory; USP34 signals were GDSC1-only with no GDSC2 replication |
| TCGA-BRCA | `src/independent_validation_tcga_*.py` | Supporting/orthogonal only; weak and incomplete for the current focus genes |
| Single-cell CNV method audit | `src/gse245601_*`, `scripts/analysis/*.R` | Method-audit provenance for GSE245601 malignant-compartment inference |
| Earlier poster iterations | `src/poster_exploration_v*`, `poster_story_v1_*`, `poster_final_*`, `poster_hero_heatmap_v1..v5`, `poster_network_mechanism_v1..v3`, `poster_pathway_v1`, `poster_depmap_v1` | Superseded development versions, retained for provenance |
