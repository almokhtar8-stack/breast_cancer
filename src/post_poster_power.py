"""post_freeze_exploratory -- minimum-detectable-effect analysis for the three
chronic resistance / recurrence datasets.

WHY THIS EXISTS
---------------
``docs/FINAL_PUBLIC_REPO_AUDIT.md`` and ``poster/README.md`` both record that
GSE111151 candidate differential expression was "largely null". "We found
nothing" and "we could not have found anything at this sample size" are
different scientific statements, and the repository could not distinguish them.
This module quantifies, per gene and per dataset, the smallest effect the study
could have detected -- using each gene's own fitted dispersion, the real library
sizes and the real design, never textbook assumptions.

WHAT THIS IS NOT
----------------
This is NOT observed (post-hoc) power. Observed power -- power evaluated at the
observed effect estimate -- is a monotone transformation of the p-value and
carries no information beyond it. No observed p-value is ever converted into a
power value here. What is computed instead is the minimum detectable effect
(MDE) at 80% power for a prespecified alpha: a property of the design and the
fitted nuisance parameters, evaluated independently of how large the observed
effect happened to be.

Two limits are stated rather than papered over:

1. The MDE is CONDITIONAL. It holds the gene's standard error fixed at its
   fitted value. A real workflow-level power calculation would resimulate the
   whole pipeline -- filtering, TMM, dispersion estimation, empirical-Bayes
   moderation and BH -- under an assumed alternative, and the standard error
   itself moves with the mean and the effect. These numbers are a conditional
   approximation, tight enough to separate "underpowered" from "well powered"
   and not tight enough to quote to two decimal places.
2. "The study was powered for an effect this size and did not find one" is NOT
   evidence of absence. At 80% power one in five true effects is missed, and
   the observed estimate is itself noisy. The interpretation column is worded
   accordingly, and equivalence testing against a prespecified smallest effect
   of interest -- which this project never prespecified -- would be the
   analysis that could actually support absence.

DATA SOURCE
-----------
``results/post_poster/de_refit/*_model_stats.tsv.gz`` from
``src/post_poster_de_refit.R`` (refits verified to reproduce the frozen tables
to machine precision). No network access; deterministic.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Sequence

import numpy as np
import pandas as pd
from scipy import optimize, stats

from src.post_poster_meta_analysis import (
    DATASETS,
    POST_FREEZE_LABEL,
    PRIMARY_SE_VARIANT,
    ROOT,
    build_long_table,
    thirteen_candidates,
)

logger = logging.getLogger(__name__)

OUT_DIR = ROOT / "results" / "post_poster" / "power"

TARGET_POWER = 0.80
# Derived, never asserted: the candidate family size comes from config via the
# meta-analysis module, so a change there cannot silently desynchronise the
# Bonferroni alpha used here.
N_CANDIDATES = len(thirteen_candidates())

# The minimum detectable effect is a property of the DESIGN. It may therefore
# only be computed from an information-based standard error. `f_calibrated`
# (= |log2FC| / sqrt(F)) contains the observed effect by construction, so an MDE
# built on it would be partly outcome-derived -- exactly the circularity this
# module exists to avoid. It is deliberately not offered here, even though the
# meta-analysis carries it as a pooling sensitivity.
MDE_SE_VARIANT = PRIMARY_SE_VARIANT

# Prespecified significance levels. The first two are design quantities fixed
# before looking at any p-value. The third is the realised BH rejection
# threshold, reported DESCRIPTIVELY only -- it is outcome-dependent (it is a
# function of the observed p-value distribution, including the gene being
# assessed), so it is deliberately NOT labelled "power at FDR 0.05".
ALPHA_SPECS: dict[str, str] = {
    "nominal_0.05": "fixed nominal alpha = 0.05, per-gene type-I error",
    "bonferroni_13": f"Bonferroni across the {N_CANDIDATES}-candidate family, 0.05/{N_CANDIDATES}",
    "observed_bh_threshold": ("largest nominal p-value that attained BH FDR <= 0.05 in "
                              "this dataset; DESCRIPTIVE ONLY -- outcome-dependent, not a "
                              "prospective design alpha, and not equivalent to FDR control"),
}


def power_at(delta: float, se: float, df_total: float, alpha: float) -> float:
    """Power of the two-sided test of H0: effect = 0 against a true effect
    ``delta``, for a statistic referred to F(1, df_total).

    Uses the noncentral F with noncentrality (delta/se)^2, which is the
    Wald/moderated-t formulation shared by edgeR's quasi-F (df1 = 1 for a
    single-coefficient contrast) and limma's moderated t.
    """
    if not (np.isfinite(se) and se > 0 and np.isfinite(df_total) and df_total > 0):
        return np.nan
    crit = stats.f.isf(alpha, dfn=1, dfd=df_total)
    ncp = (delta / se) ** 2
    # scipy's noncentral F is wrong at exactly nc = 0 (scipy 1.17 returns
    # alpha - 1 rather than alpha); it is correct for any nc > 0. At nc = 0 the
    # distribution IS the central F, so dispatch there rather than nudging the
    # noncentrality, which would silently bias the result.
    if ncp == 0.0:
        return float(stats.f.sf(crit, dfn=1, dfd=df_total))
    return float(stats.ncf.sf(crit, dfn=1, dfd=df_total, nc=ncp))


def minimum_detectable_effect(se: float, df_total: float, alpha: float,
                              target_power: float = TARGET_POWER,
                              max_multiple: float = 500.0) -> float:
    """Smallest |log2 fold change| detectable with ``target_power``.

    Solved numerically by inverting :func:`power_at` on the noncentrality
    parameter. Returns NaN where the inputs cannot support a test.
    """
    if not (np.isfinite(se) and se > 0 and np.isfinite(df_total) and df_total > 0):
        return np.nan

    def f(delta: float) -> float:
        return power_at(delta, se, df_total, alpha) - target_power

    lo, hi = 0.0, se * 3.0
    tries = 0
    while f(hi) < 0:
        hi *= 2.0
        tries += 1
        if hi > se * max_multiple:
            logger.warning("MDE search did not bracket target power (se=%.4g, df=%.4g, "
                           "alpha=%.4g) -- returning NaN", se, df_total, alpha)
            return np.nan
    if tries == 0 and f(lo) > 0:
        return 0.0
    return float(optimize.brentq(f, lo, hi, xtol=1e-10, rtol=1e-12))


def observed_bh_threshold(p_values: np.ndarray, fdr: np.ndarray,
                          level: float = 0.05) -> float:
    """The largest nominal p-value that attained BH FDR <= ``level``.

    NaN when the dataset rejected nothing at that level, in which case no MDE is
    reported for this alpha rather than a made-up fallback being substituted.
    """
    ok = np.isfinite(p_values) & np.isfinite(fdr) & (fdr <= level)
    return float(np.max(p_values[ok])) if ok.any() else np.nan


def alpha_values(long: pd.DataFrame) -> pd.DataFrame:
    """Resolve each alpha spec to a number per dataset."""
    rows = []
    for spec in DATASETS:
        sub = long[long["dataset_key"] == spec.key]
        bh = observed_bh_threshold(sub["p_value_frozen_method"].to_numpy(float),
                                   sub["fdr_frozen_method"].to_numpy(float))
        n_sig = int((sub["fdr_frozen_method"] <= 0.05).sum())
        for name, value in [("nominal_0.05", 0.05),
                            ("bonferroni_13", 0.05 / N_CANDIDATES),
                            ("observed_bh_threshold", bh)]:
            rows.append({"dataset_key": spec.key, "dataset": spec.label,
                         "alpha_name": name, "alpha": value,
                         "alpha_definition": ALPHA_SPECS[name],
                         "n_genes_tested": len(sub),
                         "n_genes_bh_fdr_le_0.05": n_sig,
                         POST_FREEZE_LABEL: True})
    out = pd.DataFrame(rows)
    for r in out.itertuples():
        if not np.isfinite(r.alpha):
            logger.warning("%s: alpha %s is undefined (dataset rejected nothing at "
                           "FDR 0.05); MDE will be NaN for that alpha",
                           r.dataset_key, r.alpha_name)
    return out


def external_reference_effects(long: pd.DataFrame, genes: Sequence[str]) -> pd.DataFrame:
    """For each (gene, dataset), the largest |log2FC| that gene shows in the
    OTHER datasets, and where it came from.

    This is an explicitly external, post-hoc benchmark -- not a truth, and not a
    prespecified effect of interest.
    """
    sub = long[long["gene_symbol"].isin(genes)][
        ["gene_symbol", "dataset_key", "dataset", "log2fc"]].copy()
    rows = []
    for gene, g in sub.groupby("gene_symbol", sort=False):
        for key in g["dataset_key"]:
            other = g[(g["dataset_key"] != key) & g["log2fc"].notna()]
            if other.empty:
                rows.append({"gene_symbol": gene, "dataset_key": key,
                             "external_reference_abs_log2fc": np.nan,
                             "external_reference_dataset": ""})
                continue
            idx = other["log2fc"].abs().idxmax()
            rows.append({"gene_symbol": gene, "dataset_key": key,
                         "external_reference_abs_log2fc": float(abs(other.loc[idx, "log2fc"])),
                         "external_reference_dataset": str(other.loc[idx, "dataset"])})
    return pd.DataFrame(rows)


def not_tested_reason(gene: str, dataset_key: str,
                      refit_dir: Path | None = None) -> str:
    """Distinguish "removed by the expression filter" from "never present".

    Both look like a missing row downstream, but they mean different things: a
    filtered gene WAS measured and fell below the detection threshold, whereas
    an absent symbol was never in the annotation at all. Conflating them would
    let "never tested" masquerade as "tested and undetectable".
    """
    refit_dir = refit_dir or (ROOT / "results" / "post_poster" / "de_refit")
    pre_filter = {
        # The counts matrices as delivered, before filterByExpr.
        "gse111151": ROOT / "results/tables/gse111151/counts_matrix.tsv.gz",
        "gse240112": ROOT / "results/tables/gse240112_pseudobulk/tumor_cell_counts.tsv.gz",
        # GSE118713 starts from an already-filtered TPM matrix; the unfiltered
        # gene universe is the frozen gene-level TPM table.
        "gse118713": ROOT / "data/processed/gse118713_gene_tpm.parquet",
    }
    path = pre_filter.get(dataset_key)
    if path is None or not path.exists():
        return "gene absent from this dataset's fitted table; pre-filter universe not checked"
    if path.suffix == ".parquet":
        universe = set(pd.read_parquet(path).get("gene_symbol", pd.Series(dtype=str)))
    else:
        cols = pd.read_csv(path, sep="\t", nrows=0).columns
        name_col = "gene_name" if "gene_name" in cols else "gene"
        universe = set(pd.read_csv(path, sep="\t", usecols=[name_col])[name_col])
    if gene in universe:
        return "present before filtering; removed by filterByExpr / expression filter"
    return "not present in this dataset's gene annotation at all; never measured"


def interpret(below_mde: float, mde_exceeds_external: float, tested: bool) -> str:
    """The single per-gene interpretation column.

    Deliberately avoids "genuine negative": failure to reject when nominally
    able to detect a hypothetical effect is not evidence of absence.
    """
    if not tested:
        return "not_tested"
    if not np.isfinite(below_mde):
        return "not_estimable"
    if bool(below_mde):
        return "null_uninformative_observed_effect_below_mde80"
    if np.isfinite(mde_exceeds_external) and bool(mde_exceeds_external):
        return "sensitive_to_observed_magnitude_but_not_to_external_benchmark"
    return "sensitive_to_observed_magnitude"


def compute_power_table(long: pd.DataFrame, genes: Sequence[str],
                        se_variant: str = MDE_SE_VARIANT) -> pd.DataFrame:
    """Per (gene, dataset, alpha) minimum detectable effect and interpretation.

    Genes absent from a dataset are emitted as explicit ``not_tested`` rows, so
    a filtered-out gene never silently becomes a null result.
    """
    if se_variant != MDE_SE_VARIANT:
        raise ValueError(
            f"minimum detectable effects may only be computed from the "
            f"information-based SE {MDE_SE_VARIANT!r}; {se_variant!r} is "
            "outcome-dependent and would make the MDE a function of the "
            "observed effect")
    alphas = alpha_values(long)
    ext = external_reference_effects(long, genes)
    se_col = f"se_{se_variant}"

    measured = long[long["gene_symbol"].isin(genes)].copy()
    rows = []
    for spec in DATASETS:
        a_ds = alphas[alphas["dataset_key"] == spec.key]
        present = measured[measured["dataset_key"] == spec.key].set_index("gene_symbol")
        for gene in genes:
            tested = gene in present.index
            g = present.loc[gene] if tested else None
            e = ext[(ext["gene_symbol"] == gene) & (ext["dataset_key"] == spec.key)]
            ext_eff = float(e["external_reference_abs_log2fc"].iloc[0]) if len(e) else np.nan
            ext_ds = str(e["external_reference_dataset"].iloc[0]) if len(e) else ""

            for a in a_ds.itertuples():
                base = {
                    "gene_symbol": gene, "dataset": spec.label,
                    "dataset_key": spec.key, "engine": spec.engine,
                    "se_variant": se_variant,
                    "alpha_name": a.alpha_name, "alpha": a.alpha,
                    "alpha_definition": a.alpha_definition,
                    "target_power": TARGET_POWER,
                    "tested_in_dataset": tested,
                    "n_biological_units": spec.n_biological_units,
                    "pseudoreplicated": spec.pseudoreplicated,
                    "group_confounded_with_biobank": spec.group_confounded_with_biobank,
                    "dataset_caveat": spec.caveat,
                    "external_reference_abs_log2fc": ext_eff,
                    "external_reference_dataset": ext_ds,
                    POST_FREEZE_LABEL: True,
                }
                if not tested:
                    rows.append({**base, "observed_log2fc": np.nan,
                                 "observed_abs_log2fc": np.nan, "se": np.nan,
                                 "df_total": np.nan, "dispersion_proxy": np.nan,
                                 "observed_p_value": np.nan, "observed_fdr": np.nan,
                                 "mde80": np.nan,
                                 "observed_abs_log2fc_below_mde80": np.nan,
                                 "mde80_exceeds_external_reference_effect": np.nan,
                                 "power_interpretation": "not_tested",
                                 "not_tested_reason": not_tested_reason(gene, spec.key)})
                    continue

                se = float(g[se_col])
                df_total = float(g["df_total"])
                mde = minimum_detectable_effect(se, df_total, float(a.alpha))
                obs = float(g["log2fc"])
                below = float(abs(obs) < mde) if np.isfinite(mde) else np.nan
                exceeds = (float(mde > ext_eff)
                           if np.isfinite(mde) and np.isfinite(ext_eff) else np.nan)
                rows.append({
                    **base, "observed_log2fc": obs, "observed_abs_log2fc": abs(obs),
                    "se": se, "df_total": df_total,
                    "dispersion_proxy": float(g["avg_log_cpm"]),
                    "observed_p_value": float(g["p_value_frozen_method"]),
                    "observed_fdr": float(g["fdr_frozen_method"]),
                    "mde80": mde,
                    "observed_abs_log2fc_below_mde80": below,
                    "mde80_exceeds_external_reference_effect": exceeds,
                    "power_interpretation": interpret(below, exceeds, True),
                    "not_tested_reason": "",
                })

    out = pd.DataFrame(rows)
    logger.info("power table: %d rows (%d genes x %d datasets x %d alphas); "
                "%d not_tested, %d not_estimable",
                len(out), len(genes), len(DATASETS), len(ALPHA_SPECS),
                int((out["power_interpretation"] == "not_tested").sum()),
                int((out["power_interpretation"] == "not_estimable").sum()))
    return out.sort_values(["gene_symbol", "dataset_key", "alpha_name"]).reset_index(drop=True)


def dataset_sensitivity_summary(power: pd.DataFrame) -> pd.DataFrame:
    """Per (dataset, alpha): the median MDE over the candidates, and how many
    candidates the dataset could say anything at all about."""
    rows = []
    for (ds, key, alpha_name), g in power.groupby(
            ["dataset", "dataset_key", "alpha_name"], sort=False):
        tested = g[g["tested_in_dataset"]]
        rows.append({
            "dataset": ds, "dataset_key": key, "alpha_name": alpha_name,
            "alpha": float(g["alpha"].iloc[0]),
            "n_candidates_tested": int(len(tested)),
            "n_candidates_not_tested": int((~g["tested_in_dataset"]).sum()),
            "median_mde80": float(tested["mde80"].median()),
            "min_mde80": float(tested["mde80"].min()) if len(tested) else np.nan,
            "max_mde80": float(tested["mde80"].max()) if len(tested) else np.nan,
            "n_observed_below_mde80": int(
                (tested["observed_abs_log2fc_below_mde80"] == 1).sum()),
            "n_sensitive_to_observed_magnitude": int(
                tested["power_interpretation"].str.startswith("sensitive").sum()),
            POST_FREEZE_LABEL: True,
        })
    return pd.DataFrame(rows).sort_values(["alpha_name", "dataset_key"]).reset_index(drop=True)


def main(out_dir: Path = OUT_DIR) -> None:
    """Write the power tables. Deterministic; no network access."""
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    out_dir.mkdir(parents=True, exist_ok=True)

    genes = thirteen_candidates()
    long = build_long_table([d.key for d in DATASETS])

    alphas = alpha_values(long)
    alphas.to_csv(out_dir / "alpha_definitions.tsv", sep="\t", index=False)
    logger.info("wrote alpha_definitions.tsv")

    power = compute_power_table(long, genes)
    power.to_csv(out_dir / "candidate_minimum_detectable_effects.tsv", sep="\t", index=False)
    logger.info("wrote candidate_minimum_detectable_effects.tsv (%d rows, SE = %s)",
                len(power), MDE_SE_VARIANT)

    summary = dataset_sensitivity_summary(power)
    summary.to_csv(out_dir / "dataset_sensitivity_summary.tsv", sep="\t", index=False)
    logger.info("wrote dataset_sensitivity_summary.tsv (%d rows)", len(summary))


if __name__ == "__main__":
    main()
