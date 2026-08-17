# DepMap Figure v1 -- Short Data Note

**Scientific question:** are the four focus genes already required for
baseline growth/survival of ER+/luminal breast-cancer cells, or is their
tamoxifen-sensitising CRISPR effect relatively distinct from baseline
dependency?

**1. DepMap release / source:** DepMap **26Q1**, the project's frozen
`active_release` (`config/config.yaml`, `independent_validation.depmap`).
Files: `CRISPRGeneEffect.csv` (Chronos gene effect),
`CRISPRGeneDependency.csv` (dependency probability), `Model.csv`
(cell-line metadata), from the release directory recorded in config
(provenance in that directory's `PROVENANCE.txt`). No new release, no
re-download, no re-run with different parameters.

**2. Frozen code/rules reused unmodified:**
`independent_validation_depmap_data.load_model()` for the ER+/luminal
rule and `raw_dir()` for release paths; the strong-dependency threshold
read from config; `post_audit_sensitivity_data.load_significant_sensitising_hits()`
for CRISPR ranks. The per-line matrix reproduces the identical rule used
by `post_audit_sensitivity_data.load_depmap_summary_for_genes()`, and its
per-gene medians/percentages match that frozen summary and
`results/tables/post_audit_sensitivity/06_top_candidate_translational_comparison.tsv`
exactly.

**3. ER+/luminal inclusion rule (frozen, not redefined):** breast lineage
(`OncotreeLineage == "Breast"`) AND DepMap's own curated
`ModelSubtypeFeatures` field contains `"ER+"` or the `"ER,"` shorthand.
"luminal" alone is not sufficient. 26Q1 contains 2154 models, 96 breast,
**22 annotated ER+/luminal**; of those, **11 have CRISPR data** and are
evaluable. n = 11 as expected -- no discrepancy.

**4. The 11 evaluable cell lines** (with DepMap subtype field):

| Cell line | ModelID | ModelSubtypeFeatures |
|---|---|---|
| CAMA1 | ACH-000783 | luminal ER+, PR+, HER2+ |
| EFM19 | ACH-000330 | luminal ER+, PR+, HER2+ |
| HCC1428 | ACH-000352 | luminal ER, PR+ |
| KPL1 | ACH-000028 | ER+ |
| MCF7 | ACH-000019 | ER+ |
| MDAMB361 | ACH-000934 | luminal ER+, HER2+ |
| MDAMB415 | ACH-000876 | luminal ER+, HER2+ |
| MFM223 | ACH-001819 | luminal ER, PR+ |
| T47D | ACH-000147 | luminal ER, PR+ |
| UACC3133 | ACH-001683 | ER+, HER2+ |
| ZR751 | ACH-000097 | luminal ER, PR+ |

**5. Chronos interpretation:** Chronos gene-effect score near 0 = little
or no baseline dependency; more negative = stronger baseline dependency;
approximately -1 is the conventional strong/common-dependency reference
level. Encoded here as a sequential palette (light = near 0, dark =
more negative) deliberately distinct from the expression heatmap's
diverging up/down semantics.

**6. Strong-dependency threshold (frozen, recovered unambiguously):**
`config.independent_validation.depmap.strong_dependency_probability_threshold`
= **0.5**, applied as **dependency probability > 0.5**. This is a
DepMap dependency-PROBABILITY cutoff, *not* a Chronos cutoff -- the
figure marks qualifying cells with a small ring (○) and says so in the
legend rather than implying a Chronos threshold.

**7. Missing-data rule:** each matrix is subset to annotated ER+/luminal
models and `dropna()` handled per matrix independently (identical to the
frozen loader). A model must have CRISPR gene-effect AND dependency data
to be evaluable; the 11 annotated ER+/luminal models without CRISPR data
are excluded and that loss is logged (22 in -> 11 evaluable, 11 lost).

**8. Candidate-level results** (all computed at build time, none
hand-typed):

| Gene | Median Chronos (ER+/luminal) | Lines with strong baseline dependency |
|---|---|---|
| KDM1A | -0.137 | 0 / 11 (0.0%) |
| TLK2 | -0.808 | 9 / 11 (81.8%) |
| USP34 | -0.063 | 0 / 11 (0.0%) |
| VEZF1 | -0.199 | 3 / 11 (27.3%) |

These matched the previously reported frozen values exactly -- no
discrepancy to report. **Scope note:** TLK2 is the strongest baseline
dependency **among the four current focus genes only**. The wider
13-gene significant-sensitising universe contains stronger
baseline-dependency genes (e.g. SUPT4H1); no claim is made about that
universe here.

**9. Baseline dependency is NOT tamoxifen-specific sensitisation.** The
project's CRISPR screen asks whether knockout makes 4-OHT/tamoxifen work
better (drug-context). DepMap asks whether knockout already impairs
fitness in unconditioned baseline culture. These are different
experiments; DepMap does not validate, replicate, or confirm the
tamoxifen-sensitisation result. Neither direction is automatically
"better": strong sensitisation with low baseline dependency may indicate
a more treatment-context-specific phenotype, while strong sensitisation
with high baseline dependency may reflect a broader cancer vulnerability.
Baseline dependency is not called "toxicity" anywhere.

**10. Cancer-cell dependency is NOT normal-tissue safety.** DepMap
contains cancer cell lines only. Low dependency in this subset says
nothing about normal-tissue tolerance, and the figure states this
explicitly.

**11. No frozen result changed.** No candidate ranking, CRISPR value,
DepMap release, subtype definition, or threshold was altered. The one new
file under `results/tables/` is
`poster_depmap_v1/depmap_26Q1_er_luminal_cellline_matrix.tsv`, a cached
reshaped extract of the frozen DepMap matrices (the source CSVs are
~440 MB each); it is a derived cache, not a new scientific result, and
deleting it triggers an identical rebuild.
