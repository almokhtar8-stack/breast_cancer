"""post_freeze_exploratory -- random-effects meta-analysis of the per-gene
resistance effect across the three chronic resistance / recurrence datasets.

WHY THIS EXISTS
---------------
The frozen evidence table ranks partly on ``resistance_fdr05_count`` (column 36
of ``results/tables/evidence_freeze/final_candidate_evidence.tsv``), a count of
how many datasets cleared FDR 0.05. Vote counting discards effect size and
treats FDR 0.049 and 0.051 as categorically different. This module replaces the
count with a pooled estimate. It does NOT produce a new shortlist, a new gene
order, or a new lead -- it describes the same evidence differently, and any
disagreement with the frozen ranking is reported as a discrepancy for a human
to decide on.

DATA SOURCES (all local; no network access)
-------------------------------------------
``results/post_poster/de_refit/{gse118713,gse111151,gse240112}_model_stats.tsv.gz``
written by ``src/post_poster_de_refit.R``, which refits each frozen model with
identical parameters and is verified to reproduce the committed frozen log2 fold
changes and p-values to machine precision.

  GSE118713  limma 3.66.0 on log2(TPM + 1), TAMR vs MCF7, n = 3 vs 3
             replicates of ONE cell-line lineage (pseudoreplicated).
             Read from the UNREDACTED frozen table: KDM1A and RCOR1 were
             withheld from the published redacted table by the preregistered
             blinding, retired 2026-08-10.
  GSE111151  edgeR glmQLF, ~cell_line + resistance_status, 11 samples,
             4 cell lines, resistant vs parental.
  GSE240112  edgeR glmQLF, ~group, tumour-cell pseudobulk, RT vs PT, 3 vs 3,
             UNPAIRED and perfectly confounded with biobank.

GSE245601 is deliberately excluded: it measures acute 12 h response, and
``docs/THERAPEUTIC_SHORTLIST_FREEZE.md`` gives the reasoning for why it is never
summed into resistance evidence. That decision is not revisited here.

THE TWO THINGS A READER MUST KNOW BEFORE USING THESE NUMBERS
-------------------------------------------------------------
1. The three datasets do not estimate the same quantity. GSE118713's effect is
   a difference of mean log2(TPM + 1) -- attenuated toward zero in an
   expression-level-dependent way by the +1 offset, and normalised by TPM
   rather than TMM. The other two are negative-binomial GLM log2 fold changes
   of counts. Pooling across that boundary is reported, but the ``edger_only``
   arm (same engine, more comparable estimand) is the primary numeric result.
2. edgeR's quasi-likelihood F test is deviance-based, not a Wald test, so no
   standard error reproduces it exactly. The pooling therefore uses ``wald``:
   the model-based standard error from the fit's own information matrix,
   scaled by the posterior QL variance. A second quantity, ``f_calibrated``
   (= |log2FC| / sqrt(F)), is carried alongside it -- but as a DIAGNOSTIC ON
   THE STATISTIC-EQUIVALENT SCALE, never as the primary weight. See below.

WHY ``wald`` IS PRIMARY AND ``f_calibrated`` IS ONLY A DIAGNOSTIC
-----------------------------------------------------------------
Because glmQLFTest computes F = LR / s2.post exactly (verified: max absolute
difference 0 over all genes), the posterior QL variance cancels out of the
ratio t_wald^2 / F, which is therefore a pure Wald-versus-likelihood-ratio
comparison. On simulated well-behaved counts the reconstruction is essentially
exact (median ratio 0.9993); dropping the s2.post scaling degrades it to 0.965;
and on the real data ``fit$dispersion`` -- the dispersion the fit actually used
-- gives a ratio closer to 1 than the trended, tagwise or common dispersion
does. The Wald covariance is therefore correctly reconstructed. On these small,
heavily overdispersed datasets the Wald statistic nonetheless runs
systematically below the likelihood-ratio statistic (median t_wald^2 / F = 0.724
in GSE111151, 0.896 in GSE240112), the familiar Wald/LR divergence for log-link
count models at small n.

That divergence is a reason to report ``f_calibrated``, and NOT a reason to pool
on it. |log2FC| / sqrt(F) is computed from the very coefficient and test
statistic being summarised, so using it as an inverse-variance weight makes the
weight outcome-dependent: it is an algebraic rescaling that reproduces F, not a
curvature/information estimate of Var(beta_hat), and it is undefined at exactly
zero effect or zero F -- which would preferentially drop the most null-looking
rows. Inverse-variance pooling requires an information-based SE, so ``wald`` is
primary. ``f_calibrated`` is retained as a labelled sensitivity that shows how
much the Wald/LR gap could matter, and the two are identical for GSE118713 by
construction (limma's quasi_f is the squared moderated t).

A LIMIT OF THE POOLED P-VALUES
------------------------------
Pooling uses a plug-in random-effects estimator: tau^2 is estimated, then the
pooled effect is referred to a normal critical value as though tau^2 were known.
At k = 2-3 studies this ignores the (large) uncertainty in tau^2 and is
anti-conservative. Every pooled p-value and interval here is therefore
DESCRIPTIVE, not a calibrated test. Hartung-Knapp is not applied because it
behaves pathologically at k = 2; the honest response is to read these intervals
as indicative and never to treat a value near 0.05 as a finding.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
import pandas as pd
import yaml
from scipy import optimize, stats

logger = logging.getLogger(__name__)

ROOT = Path(__file__).resolve().parents[1]
CONFIG_PATH = ROOT / "config" / "config.yaml"
REFIT_DIR = ROOT / "results" / "post_poster" / "de_refit"
OUT_DIR = ROOT / "results" / "post_poster" / "meta_analysis"

POST_FREEZE_LABEL = "post_freeze_exploratory"


@dataclass(frozen=True)
class DatasetSpec:
    """One resistance/recurrence dataset entering the meta-analysis."""

    key: str
    label: str
    engine: str
    # Number of independent biological units contributing to the contrast.
    # GSE118713 has three replicates but only ONE cell-line lineage per arm,
    # so its independent-unit count is 1, not 3.
    n_biological_units: int
    # Samples in the fitted model, which is not always the samples in the
    # contrast: GSE118713 fits all three groups jointly but the resistance
    # contrast itself is TAMR (3) vs MCF7 (3).
    n_samples_in_model: int
    n_samples_in_contrast: int
    pseudoreplicated: bool
    group_confounded_with_biobank: bool
    caveat: str


DATASETS: tuple[DatasetSpec, ...] = (
    DatasetSpec(
        key="gse118713", label="GSE118713", engine="limma_eBayes_log2TPM1",
        n_biological_units=1, n_samples_in_model=9, n_samples_in_contrast=6,
        pseudoreplicated=True,
        group_confounded_with_biobank=False,
        caveat=("3 replicates per arm of a single MCF7/TAMR lineage; replicates are not "
                "independent biological units, so this dataset's p-values and SEs are "
                "anticonservative by an unidentifiable factor"),
    ),
    DatasetSpec(
        key="gse111151", label="GSE111151", engine="edgeR_glmQLF",
        n_biological_units=4, n_samples_in_model=11, n_samples_in_contrast=11,
        pseudoreplicated=False,
        group_confounded_with_biobank=False,
        caveat="4 cell lines, each resistant subline compared within its own cell line",
    ),
    DatasetSpec(
        key="gse240112", label="GSE240112", engine="edgeR_glmQLF",
        n_biological_units=6, n_samples_in_model=6, n_samples_in_contrast=6,
        pseudoreplicated=False,
        group_confounded_with_biobank=True,
        caveat=("group is perfectly confounded with biobank at n = 3 vs 3; no weighting "
                "scheme separates recurrence from biobank, so this dataset's effect is "
                "not attributable to recurrence alone"),
    ),
)

DATASETS_BY_KEY = {d.key: d for d in DATASETS}

# Analysis arms. Each names the datasets pooled and how GSE118713 is treated.
ARMS: dict[str, dict] = {
    "all3": dict(
        datasets=("gse118713", "gse111151", "gse240112"),
        se_inflation={},
        description="all three datasets pooled as-is; mixes limma-on-log2(TPM+1) with "
                    "edgeR-on-counts, and gives GSE118713 its face-value SE",
    ),
    "edger_only": dict(
        datasets=("gse111151", "gse240112"),
        se_inflation={},
        description="GSE111151 + GSE240112 only; same DE engine and same NB log2 "
                    "fold-change estimand, and simultaneously the 'without the "
                    "pseudoreplicated dataset' arm",
    ),
    "all3_de3": dict(
        datasets=("gse118713", "gse111151", "gse240112"),
        se_inflation={"gse118713": float(np.sqrt(3.0))},
        description="all three, with GSE118713's SE multiplied by sqrt(3) -- the cluster "
                    "design effect at ICC = 1, i.e. treating its 3 replicates as 1 unit. "
                    "A HEURISTIC STRESS TEST, not a corrected analysis: variance "
                    "inflation cannot create the missing experimental units",
    ),
}

SE_VARIANTS = ("wald", "f_calibrated")
# See the module docstring. `wald` is the only information-based standard error
# here and is therefore the only one that may weight an inverse-variance
# pooling; `f_calibrated` is outcome-dependent by construction and is carried as
# a diagnostic on the statistic-equivalent scale.
PRIMARY_SE_VARIANT = "wald"
DIAGNOSTIC_SE_VARIANT = "f_calibrated"
PRIMARY_ARM = "edger_only"


# ---------------------------------------------------------------------------
# loading
# ---------------------------------------------------------------------------
def thirteen_candidates(config_path: Path = CONFIG_PATH) -> list[str]:
    """The 13 sensitising CRISPR candidates (FDR < 0.10, effect < 0), read from
    config rather than hardcoded."""
    config = yaml.safe_load(config_path.read_text())
    genes = list(config["gse240112"]["candidates"]["thirteen"])
    if len(genes) != 13:
        raise ValueError(f"expected 13 candidates in config, got {len(genes)}")
    return genes


def load_dataset_stats(key: str, refit_dir: Path = REFIT_DIR) -> pd.DataFrame:
    """Load one refit model-stats table and collapse to one row per gene symbol.

    Duplicate symbols (different Ensembl ids mapping to the same symbol) are
    collapsed by keeping the most highly expressed row, breaking ties on the
    gene id so the result is deterministic. Rows in, rows out and rows lost are
    logged; nothing is dropped silently.
    """
    path = refit_dir / f"{key}_model_stats.tsv.gz"
    df = pd.read_csv(path, sep="\t")
    n_in = len(df)

    df = df.dropna(subset=["gene_symbol"])
    logger.info("%s: dropped %d rows with no gene symbol", key, n_in - len(df))
    n_symbol = len(df)

    df = df.sort_values(["gene_symbol", "avg_log_cpm", "gene"],
                        ascending=[True, False, True])
    collapsed = df.drop_duplicates(subset="gene_symbol", keep="first").copy()
    logger.info(
        "%s: rows in = %d, after symbol filter = %d, after duplicate-symbol collapse = %d "
        "(lost %d duplicate-symbol rows)",
        key, n_in, n_symbol, len(collapsed), n_symbol - len(collapsed),
    )

    spec = DATASETS_BY_KEY[key]
    collapsed["dataset"] = spec.label
    collapsed["dataset_key"] = key
    collapsed["engine"] = spec.engine
    collapsed["n_biological_units"] = spec.n_biological_units
    collapsed["pseudoreplicated"] = spec.pseudoreplicated
    collapsed["group_confounded_with_biobank"] = spec.group_confounded_with_biobank
    collapsed["dataset_caveat"] = spec.caveat

    # f_calibrated SE: |log2FC| / sqrt(F) reproduces the reported test statistic
    # by construction. Undefined where F is 0 or the effect is exactly 0.
    with np.errstate(divide="ignore", invalid="ignore"):
        f_cal = np.abs(collapsed["log2fc"]) / np.sqrt(collapsed["quasi_f"])
    f_cal = f_cal.replace([np.inf, -np.inf], np.nan)
    collapsed["se_wald"] = collapsed["se_log2fc_wald"]
    collapsed["se_f_calibrated"] = f_cal
    return collapsed


def symbol_collapse_audit(keys: Sequence[str],
                          refit_dir: Path = REFIT_DIR) -> pd.DataFrame:
    """Persist the rows-in / rows-out / rows-lost of the duplicate-symbol
    collapse, so the counts survive in a committed table and not only in a log
    line that nobody re-reads."""
    rows = []
    for key in keys:
        raw = pd.read_csv(refit_dir / f"{key}_model_stats.tsv.gz", sep="\t")
        with_symbol = raw.dropna(subset=["gene_symbol"])
        collapsed = load_dataset_stats(key, refit_dir)
        spec = DATASETS_BY_KEY[key]
        rows.append({
            "dataset": spec.label, "dataset_key": key,
            "rows_in_refit_table": len(raw),
            "rows_lost_no_gene_symbol": len(raw) - len(with_symbol),
            "rows_after_symbol_filter": len(with_symbol),
            "rows_lost_duplicate_symbol_collapse": len(with_symbol) - len(collapsed),
            "rows_out_unique_symbols": len(collapsed),
            "collapse_rule": ("keep the row with the highest avg_log_cpm; ties broken on "
                              "gene id ascending, so the result is deterministic"),
            "n_genes_no_wald_se": int((~np.isfinite(collapsed["se_wald"])).sum()),
            "n_genes_no_f_calibrated_se": int(
                (~np.isfinite(collapsed["se_f_calibrated"])).sum()),
            POST_FREEZE_LABEL: True,
        })
    return pd.DataFrame(rows)


def build_long_table(keys: Sequence[str], refit_dir: Path = REFIT_DIR) -> pd.DataFrame:
    """One row per (gene_symbol, dataset) for the requested datasets."""
    frames = [load_dataset_stats(k, refit_dir) for k in keys]
    cols = ["gene_symbol", "gene", "dataset", "dataset_key", "engine", "log2fc",
            "se_wald", "se_f_calibrated", "df_total", "quasi_f",
            "p_value_frozen_method", "fdr_frozen_method", "avg_log_cpm",
            "n_biological_units", "pseudoreplicated",
            "group_confounded_with_biobank", "dataset_caveat"]
    return pd.concat([f[cols] for f in frames], ignore_index=True)


# ---------------------------------------------------------------------------
# heterogeneity estimators
# ---------------------------------------------------------------------------
def cochran_q(effects: np.ndarray, variances: np.ndarray, tau2: float = 0.0) -> float:
    """Cochran's Q at a given between-study variance (Q at tau2 = 0 is the
    conventional heterogeneity statistic)."""
    w = 1.0 / (variances + tau2)
    mu = float(np.sum(w * effects) / np.sum(w))
    return float(np.sum(w * (effects - mu) ** 2))


def tau2_dersimonian_laird(effects: np.ndarray, variances: np.ndarray) -> float:
    """DerSimonian-Laird moment estimator of the between-study variance."""
    k = len(effects)
    if k < 2:
        return 0.0
    w = 1.0 / variances
    q = cochran_q(effects, variances, 0.0)
    c = float(np.sum(w) - np.sum(w ** 2) / np.sum(w))
    if c <= 0:
        return 0.0
    return max(0.0, (q - (k - 1)) / c)


def tau2_paule_mandel(effects: np.ndarray, variances: np.ndarray,
                      upper: float = 1e4) -> float:
    """Paule-Mandel estimator: the tau2 at which Q(tau2) equals its df.

    Preferred over DerSimonian-Laird as a sensitivity, though with k = 2 or 3
    studies neither estimator is meaningfully identifiable.
    """
    k = len(effects)
    if k < 2:
        return 0.0
    target = float(k - 1)
    if cochran_q(effects, variances, 0.0) <= target:
        return 0.0

    def f(t2: float) -> float:
        return cochran_q(effects, variances, t2) - target

    if f(upper) > 0:
        return upper
    return float(optimize.brentq(f, 0.0, upper, xtol=1e-12, rtol=1e-10))


@dataclass(frozen=True)
class PooledResult:
    k: int
    pooled_effect: float
    pooled_se: float
    ci_low: float
    ci_high: float
    z: float
    p_value: float
    tau2: float
    q: float
    q_df: int
    q_p: float
    i2: float


def pool(effects: np.ndarray, variances: np.ndarray, tau2: float) -> PooledResult:
    """Inverse-variance pooling at a fixed between-study variance.

    tau2 = 0 gives the fixed-effect result; a fitted tau2 gives the
    random-effects result.
    """
    k = len(effects)
    w = 1.0 / (variances + tau2)
    est = float(np.sum(w * effects) / np.sum(w))
    se = float(np.sqrt(1.0 / np.sum(w)))
    z = est / se if se > 0 else np.nan
    p = float(2.0 * stats.norm.sf(abs(z))) if np.isfinite(z) else np.nan
    q = cochran_q(effects, variances, 0.0)
    q_df = k - 1
    q_p = float(stats.chi2.sf(q, q_df)) if q_df > 0 else np.nan
    i2 = float(max(0.0, (q - q_df) / q) * 100.0) if q > 0 and q_df > 0 else 0.0
    return PooledResult(
        k=k, pooled_effect=est, pooled_se=se,
        ci_low=est - 1.959963984540054 * se, ci_high=est + 1.959963984540054 * se,
        z=z, p_value=p, tau2=tau2, q=q, q_df=q_df, q_p=q_p, i2=i2,
    )


def stouffer(p_values: np.ndarray, signs: np.ndarray,
             weights: np.ndarray) -> tuple[float, float]:
    """Signed weighted Stouffer combination.

    Returns (z, two-sided p). Tests a DIRECTIONAL null -- "is there consistent
    signed evidence" -- which is a different estimand from a common effect size,
    and it inherits every invalidity of its input p-values (a p-value from a
    pseudoreplicated or fully confounded design does not become valid by being
    combined).
    """
    p_values = np.clip(p_values, np.finfo(float).tiny, 1.0)
    z_i = signs * stats.norm.isf(p_values / 2.0)
    z = float(np.sum(weights * z_i) / np.sqrt(np.sum(weights ** 2)))
    return z, float(2.0 * stats.norm.sf(abs(z)))


def benjamini_hochberg(p: np.ndarray) -> np.ndarray:
    """BH FDR. NaN p-values are carried through as NaN and excluded from the
    ranking rather than being silently treated as 1."""
    out = np.full(len(p), np.nan)
    ok = np.isfinite(p)
    if not ok.any():
        return out
    pv = p[ok]
    n = len(pv)
    order = np.argsort(pv)
    ranked = pv[order]
    adj = ranked * n / np.arange(1, n + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    res = np.empty(n)
    res[order] = np.clip(adj, 0, 1)
    out[ok] = res
    return out


# ---------------------------------------------------------------------------
# the meta-analysis itself
# ---------------------------------------------------------------------------
def meta_analyse(long: pd.DataFrame, arm: str, se_variant: str,
                 genes: Iterable[str] | None = None) -> pd.DataFrame:
    """Random-effects meta-analysis for one arm and one SE definition.

    Only genes measured in EVERY dataset of the arm are pooled; genes measured
    in a subset are still returned, with ``pooled_effect`` NaN and
    ``exclusion_reason`` recording why, so nothing disappears silently.
    """
    spec = ARMS[arm]
    keys = list(spec["datasets"])
    inflation = spec["se_inflation"]

    sub = long[long["dataset_key"].isin(keys)].copy()
    se_col = f"se_{se_variant}"
    sub["se_used"] = sub[se_col] * sub["dataset_key"].map(inflation).fillna(1.0)

    if genes is not None:
        genes = list(genes)
        sub = sub[sub["gene_symbol"].isin(genes)]
        universe = genes
    else:
        universe = sorted(sub["gene_symbol"].unique())

    # Group once rather than filtering per gene: a per-gene boolean mask over
    # the long table is quadratic and genome-wide scope makes that intractable.
    groups = {name: frame for name, frame in sub.groupby("gene_symbol", sort=False)}
    empty = sub.iloc[0:0]

    rows = []
    for gene in universe:
        g = groups.get(gene, empty)
        present = set(g["dataset_key"])
        usable = g[np.isfinite(g["se_used"]) & (g["se_used"] > 0)
                   & np.isfinite(g["log2fc"])]
        missing = [k for k in keys if k not in present]
        unusable = [k for k in present if k not in set(usable["dataset_key"])]

        base = {
            "gene_symbol": gene,
            "arm": arm,
            "se_variant": se_variant,
            "n_datasets_in_arm": len(keys),
            "n_datasets_with_gene": len(present),
            "n_datasets_used": len(usable),
            "datasets_used": ",".join(sorted(usable["dataset_key"])),
            "datasets_missing_gene": ",".join(sorted(missing)),
            "datasets_unusable_se": ",".join(sorted(unusable)),
            "any_dataset_pseudoreplicated": bool(usable["pseudoreplicated"].any())
            if len(usable) else False,
            "any_group_confounded_with_biobank": bool(
                usable["group_confounded_with_biobank"].any()) if len(usable) else False,
            POST_FREEZE_LABEL: True,
        }

        if len(usable) < len(keys):
            reason = []
            if missing:
                reason.append(f"not tested in {','.join(sorted(missing))}")
            if unusable:
                reason.append(f"no usable {se_variant} SE in {','.join(sorted(unusable))}")
            rows.append({**base, "exclusion_reason": "; ".join(reason),
                         "pooled_effect": np.nan, "pooled_se": np.nan,
                         "ci_low": np.nan, "ci_high": np.nan, "z": np.nan,
                         "p_value": np.nan, "tau2_dl": np.nan, "tau2_pm": np.nan,
                         "q": np.nan, "q_df": np.nan, "q_p": np.nan, "i2": np.nan,
                         "pooled_effect_fixed": np.nan, "pooled_effect_pm": np.nan,
                         "p_value_pm": np.nan, "direction_consistent": np.nan,
                         "stouffer_z_equal": np.nan, "stouffer_p_equal": np.nan,
                         "stouffer_z_units": np.nan, "stouffer_p_units": np.nan})
            continue

        eff = usable["log2fc"].to_numpy(float)
        var = usable["se_used"].to_numpy(float) ** 2
        t2_dl = tau2_dersimonian_laird(eff, var)
        t2_pm = tau2_paule_mandel(eff, var)
        re_dl = pool(eff, var, t2_dl)
        re_pm = pool(eff, var, t2_pm)
        fe = pool(eff, var, 0.0)

        signs = np.sign(eff)
        pvals = usable["p_value_frozen_method"].to_numpy(float)
        units = usable["n_biological_units"].to_numpy(float)
        z_eq, p_eq = stouffer(pvals, signs, np.ones(len(eff)))
        z_un, p_un = stouffer(pvals, signs, np.sqrt(units))

        rows.append({
            **base, "exclusion_reason": "",
            "pooled_effect": re_dl.pooled_effect, "pooled_se": re_dl.pooled_se,
            "ci_low": re_dl.ci_low, "ci_high": re_dl.ci_high,
            "z": re_dl.z, "p_value": re_dl.p_value,
            "tau2_dl": t2_dl, "tau2_pm": t2_pm,
            "q": re_dl.q, "q_df": re_dl.q_df, "q_p": re_dl.q_p, "i2": re_dl.i2,
            "pooled_effect_fixed": fe.pooled_effect,
            "pooled_effect_pm": re_pm.pooled_effect, "p_value_pm": re_pm.p_value,
            "direction_consistent": bool(np.all(signs == signs[0])),
            "stouffer_z_equal": z_eq, "stouffer_p_equal": p_eq,
            "stouffer_z_units": z_un, "stouffer_p_units": p_un,
        })

    out = pd.DataFrame(rows)
    for src, dst in [("p_value", "fdr"), ("p_value_pm", "fdr_pm"),
                     ("stouffer_p_equal", "stouffer_fdr_equal"),
                     ("stouffer_p_units", "stouffer_fdr_units")]:
        out[dst] = benjamini_hochberg(out[src].to_numpy(float))
    logger.info(
        "arm=%s se=%s: %d genes considered, %d pooled, %d not pooled",
        arm, se_variant, len(out), int(out["pooled_effect"].notna().sum()),
        int(out["pooled_effect"].isna().sum()),
    )
    return out.sort_values("gene_symbol").reset_index(drop=True)


def forest_table(long: pd.DataFrame, genes: Sequence[str]) -> pd.DataFrame:
    """Forest-plot-ready per-(gene, dataset) table: the point estimate and 95%
    interval each dataset contributes, under both SE definitions."""
    sub = long[long["gene_symbol"].isin(genes)].copy()
    for variant in SE_VARIANTS:
        se = sub[f"se_{variant}"]
        sub[f"ci_low_{variant}"] = sub["log2fc"] - 1.959963984540054 * se
        sub[f"ci_high_{variant}"] = sub["log2fc"] + 1.959963984540054 * se
    sub[POST_FREEZE_LABEL] = True
    cols = ["gene_symbol", "dataset", "dataset_key", "engine", "log2fc",
            "se_wald", "ci_low_wald", "ci_high_wald",
            "se_f_calibrated", "ci_low_f_calibrated", "ci_high_f_calibrated",
            "quasi_f", "df_total", "p_value_frozen_method", "fdr_frozen_method",
            "avg_log_cpm", "n_biological_units", "pseudoreplicated",
            "group_confounded_with_biobank", "dataset_caveat", POST_FREEZE_LABEL]
    return sub[cols].sort_values(["gene_symbol", "dataset"]).reset_index(drop=True)


def se_variant_diagnostic(long: pd.DataFrame) -> pd.DataFrame:
    """Per dataset, how far the Wald statistic sits from the reported quasi-F.

    This is the evidence for choosing ``f_calibrated`` as the primary SE, and
    it belongs in a table rather than in prose so it can be re-checked.
    """
    rows = []
    for spec in DATASETS:
        g = long[long["dataset_key"] == spec.key].copy()
        t_wald = g["log2fc"] / g["se_wald"]
        ok = np.isfinite(t_wald) & np.isfinite(g["quasi_f"]) & (g["quasi_f"] > 0)
        ratio = (t_wald[ok] ** 2) / g.loc[ok, "quasi_f"]
        rows.append({
            "dataset": spec.label, "dataset_key": spec.key, "engine": spec.engine,
            "n_genes": len(g), "n_genes_with_ratio": int(ok.sum()),
            "n_genes_no_wald_se": int((~np.isfinite(g["se_wald"])).sum()),
            "t_wald_sq_over_quasi_f_median": float(ratio.median()),
            "t_wald_sq_over_quasi_f_q25": float(ratio.quantile(0.25)),
            "t_wald_sq_over_quasi_f_q75": float(ratio.quantile(0.75)),
            "se_ratio_wald_over_f_calibrated_median": float(
                (g["se_wald"] / g["se_f_calibrated"]).median(skipna=True)),
            "interpretation": ("wald and f_calibrated are identical by construction "
                               "(quasi_f is the squared moderated t)"
                               if spec.engine.startswith("limma")
                               else "wald statistic runs below the deviance-based quasi-F; "
                                    "wald SE is the conservative variant"),
            POST_FREEZE_LABEL: True,
        })
    return pd.DataFrame(rows)


def vote_count_comparison(candidate_meta: pd.DataFrame,
                          evidence_path: Path | None = None) -> pd.DataFrame:
    """Rank the candidates under vote counting versus under the pooled effect.

    ``resistance_fdr05_count`` is read from the frozen evidence table and is
    never recomputed or modified here.
    """
    evidence_path = evidence_path or (
        ROOT / "results/tables/evidence_freeze/final_candidate_evidence.tsv")
    ev = pd.read_csv(evidence_path, sep="\t")
    keep = ["gene", "resistance_fdr05_count", "resistance_nominal_p05_count",
            "crispr_effect", "crispr_fdr", "global_rank", "resistance_rank"]
    ev = ev[[c for c in keep if c in ev.columns]]

    merged = candidate_meta.merge(ev, left_on="gene_symbol", right_on="gene", how="left")
    n_lost = int(merged["gene"].isna().sum())
    logger.info("vote-count comparison: %d candidates, %d with no frozen evidence row",
                len(merged), n_lost)

    # Vote-count rank: more datasets at FDR < 0.05 is "better" (rank 1 = best).
    merged["vote_count_rank"] = merged["resistance_fdr05_count"].rank(
        ascending=False, method="min")
    # Ranked on the pooled z, i.e. signed evidence strength -- NOT on effect
    # magnitude. The column name says so, because the vote count it replaces was
    # also an evidence-strength summary rather than an effect-size one.
    merged["pooled_z_rank"] = merged["z"].rank(ascending=False, method="min")
    merged["rank_shift"] = merged["vote_count_rank"] - merged["pooled_z_rank"]
    merged[POST_FREEZE_LABEL] = True
    return merged.sort_values("pooled_z_rank").reset_index(drop=True)


def main(out_dir: Path = OUT_DIR) -> None:
    """Write every meta-analysis table. Deterministic; no network access."""
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    out_dir.mkdir(parents=True, exist_ok=True)

    candidates = thirteen_candidates()
    all_keys = [d.key for d in DATASETS]
    long = build_long_table(all_keys)
    logger.info("long table: %d (gene, dataset) rows", len(long))

    audit = symbol_collapse_audit(all_keys)
    audit.to_csv(out_dir / "gene_symbol_collapse_audit.tsv", sep="\t", index=False)
    logger.info("wrote gene_symbol_collapse_audit.tsv (%d rows)", len(audit))

    forest = forest_table(long, candidates)
    forest.to_csv(out_dir / "forest_plot_input_candidates.tsv", sep="\t", index=False)
    logger.info("wrote forest_plot_input_candidates.tsv (%d rows)", len(forest))

    candidate_frames, genomewide_frames = [], []
    for arm in ARMS:
        for variant in SE_VARIANTS:
            cand = meta_analyse(long, arm, variant, genes=candidates)
            candidate_frames.append(cand)
            gw = meta_analyse(long, arm, variant, genes=None)
            genomewide_frames.append(gw)
            gw.to_csv(out_dir / f"genomewide_{arm}_{variant}.tsv.gz",
                      sep="\t", index=False)
            logger.info("wrote genomewide_%s_%s.tsv.gz (%d genes)", arm, variant, len(gw))

    cand_all = pd.concat(candidate_frames, ignore_index=True)
    cand_all.to_csv(out_dir / "candidates_meta_analysis.tsv", sep="\t", index=False)
    logger.info("wrote candidates_meta_analysis.tsv (%d rows)", len(cand_all))

    diagnostic = se_variant_diagnostic(long)
    diagnostic.to_csv(out_dir / "se_variant_diagnostic.tsv", sep="\t", index=False)
    logger.info("wrote se_variant_diagnostic.tsv (%d rows)", len(diagnostic))

    primary = cand_all[(cand_all["arm"] == PRIMARY_ARM)
                       & (cand_all["se_variant"] == PRIMARY_SE_VARIANT)]
    comparison = vote_count_comparison(primary)
    comparison["primary_arm"] = PRIMARY_ARM
    comparison["primary_se_variant"] = PRIMARY_SE_VARIANT
    comparison.to_csv(out_dir / "vote_count_vs_pooled_effect.tsv", sep="\t", index=False)
    logger.info("wrote vote_count_vs_pooled_effect.tsv (%d rows, arm=%s se=%s)",
                len(comparison), PRIMARY_ARM, PRIMARY_SE_VARIANT)

    arms = pd.DataFrame([
        {"arm": a, "datasets": ",".join(s["datasets"]),
         "se_inflation": ";".join(f"{k}={v:.6f}" for k, v in s["se_inflation"].items()),
         "description": s["description"], POST_FREEZE_LABEL: True}
        for a, s in ARMS.items()
    ])
    arms.to_csv(out_dir / "arm_definitions.tsv", sep="\t", index=False)
    logger.info("wrote arm_definitions.tsv")


if __name__ == "__main__":
    main()
