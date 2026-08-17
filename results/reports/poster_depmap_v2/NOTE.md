# DepMap Figure v2 -- Short Data Note

**Purpose:** a communication rebuild of `poster_depmap_v1`. v1's 11 x 4
Chronos heatmap is scientifically complete but makes the viewer decode
Chronos, colour intensity, dependency probability, ring markers, 11
cell-line rows and a second bar plot before reaching the conclusion. v2
answers the question directly with one dot per candidate:

> Does tamoxifen sensitisation occur in genes that cancer cells already
> depend on at baseline?

**Same frozen data as v1.** No DepMap re-run, no new release, no new
threshold, no changed candidate ranking. v2 imports v1's loaders
(`poster_depmap_v1.load_cellline_table` / `dependency_summary`) and the
frozen CRISPR loader directly.

**1. Frozen CRISPR source:**
`post_audit_sensitivity_data.load_significant_sensitising_hits()` --
the same 13-hit sensitising table (FDR < 0.10, effect < 0) used by the
CRISPR discovery figure. Column read: `effect_size`.

**2. Frozen DepMap source:** DepMap **26Q1** (the project's `active_release`
in `config/config.yaml`): `CRISPRGeneEffect.csv`,
`CRISPRGeneDependency.csv`, `Model.csv`, accessed through the frozen
ER+/luminal rule in `independent_validation_depmap_data.load_model()`.

**3. X-axis transformation:**
`sensitisation_strength = -1 x crispr_effect`.
A more negative frozen CRISPR effect means a stronger sensitising
knockout, so negating it makes "further right = stronger tamoxifen
sensitisation" read intuitively. **This is a display transform only** --
the frozen `effect_size` values are never modified, and the x-axis shows
no numeric ticks (only a Weaker -> Stronger direction cue) so no viewer
has to reason about negative effect sizes. No "strong sensitiser" cutoff
was invented, so no vertical threshold line is drawn.

**4. Raw frozen CRISPR values and plotted x coordinates:**

| Gene | Frozen CRISPR effect | Plotted x (= -effect) | CRISPR rank |
|---|---|---|---|
| KDM1A | -2.167336 | 2.167336 | 1 / 13 |
| TLK2 | -1.848198 | 1.848198 | 4 / 13 |
| VEZF1 | -1.602445 | 1.602445 | 8 / 13 |
| USP34 | -1.391298 | 1.391298 | 12 / 13 |

**5. The exact 11-cell-line ER+/luminal subset** (frozen rule: breast
lineage AND DepMap's curated `ModelSubtypeFeatures` contains "ER+" or
"ER,"; 22 annotated, 11 CRISPR-evaluable): CAMA1, EFM19, HCC1428, KPL1,
MCF7, MDAMB361, MDAMB415, MFM223, T47D, UACC3133, ZR751.

**6. Strong-dependency criterion:** DepMap dependency probability >
**0.5**, read from
`config.independent_validation.depmap.strong_dependency_probability_threshold`.
This is a dependency-probability cutoff, not a Chronos cutoff.

**7. Y-axis values (n/N and %, computed at build time, none hand-typed):**

| Gene | Lines with strong baseline dependency | Plotted y |
|---|---|---|
| KDM1A | 0 / 11 | 0.0% |
| TLK2 | 9 / 11 | 81.8% |
| USP34 | 0 / 11 | 0.0% |
| VEZF1 | 3 / 11 | 27.3% |

These matched the previously reported frozen values exactly -- no
discrepancy to report.

**8. Supporting detail -- median Chronos (ER+/luminal), deliberately kept
out of the figure:** KDM1A -0.137, TLK2 -0.808, USP34 -0.063,
VEZF1 -0.199. Chronos near 0 = little/no baseline dependency; more
negative = stronger; about -1 is the conventional strong-dependency
reference level. The per-cell-line Chronos matrix, ring markers and
colorbar remain available in `poster_depmap_v1`.

**9. Baseline dependency is NOT normal-tissue safety.** DepMap contains
cancer cell lines only; a low percentage here says nothing about
normal-tissue tolerance, and it does not prove tamoxifen specificity.
Baseline dependency is also not "toxicity". Neither axis direction is
"good" or "bad", and no region of the plot is labelled safe, specific,
superior, or a therapeutic window.

**Correct interpretation:** KDM1A combines strong tamoxifen sensitisation
with little baseline dependency in this cancer-cell-line subset; TLK2
combines strong sensitisation with substantial baseline dependency; USP34
shows weaker sensitisation with little baseline dependency; VEZF1 is
intermediate on both. This separates treatment-associated sensitisation
from general cancer-cell dependency -- one evidence layer, not a ranking.
**Scope note:** TLK2 is the strongest baseline dependency among these
four focus genes only; the wider 13-hit universe contains stronger
baseline-dependency genes.

**10. No frozen result changed.** No candidate ranking, CRISPR value,
DepMap release, subtype definition, or threshold was altered, and v1's
outputs were not overwritten.
