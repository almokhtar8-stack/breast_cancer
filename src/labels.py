"""Guide- and gene-level tamoxifen-response labels from the Hany screen.

Source: normalised sgRNA counts from Hany 2023 Data S1, median-ratio
normalised by MAGeCK-VISPR. DOI: 10.1126/sciadv.add3685.

Reads Data S1's guide-level MAGeCK-VISPR-normalised sgRNA count matrix,
restricts it to the E2-alone and E2 + 4-OHT arms (PREANALYSIS.md §1-§3),
filters out unstable guides and genes, and estimates a gene-by-treatment
interaction term (PREANALYSIS.md §2).

Normalisation: none is applied here beyond what Data S1 already carries.
The counts are already median-ratio normalised across samples by the
original authors' MAGeCK-VISPR pipeline, so a second total-count rescaling
(e.g. counts-per-million) would be redundant and cannot be defended — it
would rescale an already between-sample-comparable matrix a second time,
on top of whatever scaling MAGeCK-VISPR already chose. The modelled
response is log2(normalised_count + 1) (``LOG2_PSEUDOCOUNT``): a
pseudocount of 1 keeps zero and near-zero counts finite on the log scale
while being negligible relative to the observed count magnitudes (median
~170 across the modelled arms in this file).

Interaction model: Data S1's values are non-integer (median-ratio
normalised by the original authors, not raw integer read counts), which
makes a discrete-count negative-binomial GLM a poor fit. This module
therefore uses the moderated-linear-model option named in PREANALYSIS.md
§2: for each gene, a guide-blocked linear model is fit on log2-normalised
counts (``log2norm ~ guide + treatment``), so the treatment coefficient
reflects within-guide 4-OHT-vs-E2 change with between-guide baseline
differences removed, rather than a simple difference of raw means. Per-gene
residual variances are then shrunk toward a common prior via an
empirical-Bayes moment estimator (Smyth 2004; the method behind
``limma::eBayes``), borrowing strength across genes before computing
moderated t-statistics and p-values.

Sign convention (PREANALYSIS.md §2): ``effect_size`` is the treatment
coefficient for E2 + 4-OHT relative to E2 alone, on the log2-normalised-
count scale. Negative = more depleted under 4-OHT (candidate knockout
sensitises cells to tamoxifen, i.e. cells need that gene to tolerate the
drug). Positive = enriched under 4-OHT (candidate knockout confers a
growth advantage under tamoxifen). This sign must be carried through
unchanged in every downstream table, figure, and model output.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from scipy import stats
from scipy.optimize import brentq
from scipy.special import digamma, polygamma
from statsmodels.stats.multitest import multipletests

logger = logging.getLogger(__name__)

REPLICATES: tuple[str, ...] = ("R1", "R2")
CONTROL_T0 = "Control T0"
ARM_E2 = "E2"
ARM_E2_4OHT = "E2 + 4-OHT"
# Order matters only for readability of downstream column lists.
REQUIRED_ARMS: tuple[str, ...] = (CONTROL_T0, ARM_E2, ARM_E2_4OHT)
MODELLED_ARMS: tuple[str, ...] = (ARM_E2, ARM_E2_4OHT)
MIN_GUIDES_PER_GENE = 3
LOG2_PSEUDOCOUNT = 1.0


@dataclass(frozen=True)
class FilterStep:
    """Row/gene bookkeeping for one filtering step, for the audit log."""

    name: str
    rows_in: int
    rows_out: int
    genes_in: int
    genes_out: int

    @property
    def rows_lost(self) -> int:
        return self.rows_in - self.rows_out

    @property
    def genes_lost(self) -> int:
        return self.genes_in - self.genes_out

    def log(self) -> None:
        logger.info(
            "filter[%s]: rows %d -> %d (lost %d); genes %d -> %d (lost %d)",
            self.name,
            self.rows_in,
            self.rows_out,
            self.rows_lost,
            self.genes_in,
            self.genes_out,
            self.genes_lost,
        )

    def as_record(self) -> dict[str, object]:
        return {
            "step": self.name,
            "rows_in": self.rows_in,
            "rows_out": self.rows_out,
            "rows_lost": self.rows_lost,
            "genes_in": self.genes_in,
            "genes_out": self.genes_out,
            "genes_lost": self.genes_lost,
        }


def _clean_condition(name: str) -> str:
    """Strip the stray literal quote left by unescaping Data S1's TSV.

    The source file encodes the names using doubled-quote CSV-style
    escaping, e.g. six quote characters wrapping T0. pandas' default TSV
    reader honours that quoting and unescapes it, leaving a name that
    still carries a literal embedded quote character before and after T0.
    That embedded quote is stripped here.
    """
    return name.replace('"', "").strip()


def _validate_header(columns: list[str]) -> None:
    """Check the flattened header is unambiguous before anything reads it.

    Guards against silently misreading the two-row TSV header: exactly one
    ``sgRNA`` and one ``Gene`` id column, no duplicate flattened sample
    names (which would mean two source columns collapsed onto one), and
    every required arm/replicate combination present exactly once.
    """
    for id_col in ("sgRNA", "Gene"):
        n = columns.count(id_col)
        if n != 1:
            raise ValueError(f"expected exactly one {id_col!r} id column, found {n}")

    sample_cols = [c for c in columns if c not in ("sgRNA", "Gene")]
    duplicates = {c for c in sample_cols if sample_cols.count(c) > 1}
    if duplicates:
        raise ValueError(f"duplicate flattened sample columns: {sorted(duplicates)}")

    required = [f"{arm}__{rep}" for arm in REQUIRED_ARMS for rep in REPLICATES]
    missing = [c for c in required if c not in sample_cols]
    if missing:
        raise ValueError(f"required arm/replicate columns missing from header: {missing}")


def _validate_values(df: pd.DataFrame, sample_cols: list[str]) -> None:
    """Check id and count columns are usable before any arithmetic runs on them.

    Non-finite or negative counts would silently corrupt the zero-at-T0
    filter (``NaN == 0`` is False, so a missing count would pass QC as if
    it were nonzero) and any downstream log/model arithmetic. Null or
    duplicate identifiers would silently vanish from ``groupby("Gene")``
    (which drops NaN keys by default) or make sgRNA-based guide counts
    ambiguous. All of these are checked once, here, rather than trusted
    implicitly downstream.
    """
    if df["sgRNA"].isna().any() or df["Gene"].isna().any() or (df["Gene"].astype(str).str.strip() == "").any():
        raise ValueError("null or blank sgRNA/Gene identifiers in Data S1")
    if df["sgRNA"].duplicated().any():
        raise ValueError("duplicate sgRNA identifiers in Data S1")

    counts = df[sample_cols].apply(pd.to_numeric, errors="coerce")
    if counts.isna().to_numpy().any() or not np.isfinite(counts.to_numpy()).all():
        raise ValueError("non-numeric or non-finite counts in Data S1 sample columns")
    if (counts.to_numpy() < 0).any():
        raise ValueError("negative counts in Data S1 sample columns")


def load_normalised_counts(path: str | Path) -> pd.DataFrame:
    """Load Data S1 into a tidy wide frame with clean column names.

    Data S1 is tab-separated with two header rows: row 0 gives the
    replicate (R1/R2), row 1 gives the condition. Columns are flattened to
    ``"<condition>__<replicate>"`` (e.g. ``"E2 + 4-OHT__R1"``), with the
    stray quote in the control arm names removed (see ``_clean_condition``).
    The header and the id/count columns are validated (see
    ``_validate_header`` and ``_validate_values``) so a malformed or
    reshuffled source file fails loudly instead of silently mis-aligning.
    """
    path = Path(path)
    counts = pd.read_csv(path, sep="\t", header=[0, 1])

    flat_columns = []
    for level0, level1 in counts.columns:
        if level1 in ("sgRNA", "Gene"):
            flat_columns.append(level1)
        else:
            flat_columns.append(f"{_clean_condition(level1)}__{level0}")
    counts.columns = flat_columns
    _validate_header(flat_columns)

    sample_cols = [c for c in flat_columns if c not in ("sgRNA", "Gene")]
    _validate_values(counts, sample_cols)

    logger.info("load_normalised_counts: read %d rows from %s", len(counts), path)
    return counts


def select_arms(df: pd.DataFrame) -> pd.DataFrame:
    """Keep only the id columns and the T0/E2/E2+4-OHT sample columns.

    db-cAMP, ICI (fulvestrant), their combination arms, and Control T4w are
    dropped entirely, per PREANALYSIS.md §3 ("Only E2 and E2+4-OHT arms are
    used"). Control T0 is retained only as the zero-count QC reference for
    ``filter_zero_at_t0`` — it is not part of the modelled contrast. This
    step selects columns, not rows, so no rows are lost here.
    """
    sample_cols = [f"{arm}__{rep}" for arm in REQUIRED_ARMS for rep in REPLICATES]
    keep = ["sgRNA", "Gene"] + sample_cols
    missing = [c for c in keep if c not in df.columns]
    if missing:
        raise KeyError(f"expected columns missing after cleaning: {missing}")

    selected = df[keep].copy()
    logger.info(
        "select_arms: kept arms %s (%d sample columns of %d input data columns); "
        "db-cAMP, ICI, and Control T4w arms excluded",
        REQUIRED_ARMS,
        len(sample_cols),
        df.shape[1] - 2,
    )
    return selected


def add_log2_normalised(
    df: pd.DataFrame, sample_cols: list[str], pseudocount: float = LOG2_PSEUDOCOUNT
) -> pd.DataFrame:
    """Add per-sample log2(normalised count + pseudocount) columns.

    No further rescaling (e.g. counts-per-million) is applied here: Data
    S1's counts are already median-ratio normalised across samples by
    MAGeCK-VISPR (see module docstring), so a second total-count
    normalisation on top of an already-normalised matrix is unjustified.
    """
    out = df.copy()
    for col in sample_cols:
        out[f"log2norm_{col}"] = np.log2(df[col] + pseudocount)
    return out


def filter_zero_at_t0(df: pd.DataFrame) -> tuple[pd.DataFrame, FilterStep]:
    """Drop guides with a zero count in either Control T0 replicate.

    PREANALYSIS.md §3 (quoted verbatim; that document's own wording, not
    edited here): "any guide with zero raw counts in the T0 (plasmid/day-0)
    reference is excluded — a zero baseline makes any fold-change or
    interaction estimate for that guide undefined/unstable."
    """
    t0_cols = [f"{CONTROL_T0}__{rep}" for rep in REPLICATES]
    zero_at_t0 = (df[t0_cols] == 0).any(axis=1)
    kept = df.loc[~zero_at_t0].copy()

    step = FilterStep(
        name="zero_count_at_t0",
        rows_in=len(df),
        rows_out=len(kept),
        genes_in=df["Gene"].nunique(),
        genes_out=kept["Gene"].nunique(),
    )
    step.log()
    return kept, step


def filter_low_representation_genes(
    df: pd.DataFrame, min_guides: int = MIN_GUIDES_PER_GENE
) -> tuple[pd.DataFrame, FilterStep]:
    """Drop genes with fewer than ``min_guides`` surviving guides.

    PREANALYSIS.md §3: "genes represented by fewer than 3 surviving guides
    (after the above filters) are excluded from gene-level aggregation —
    insufficient within-gene replication to estimate a stable interaction
    term."
    """
    guide_counts = df.groupby("Gene")["sgRNA"].nunique()
    keep_genes = guide_counts[guide_counts >= min_guides].index
    kept = df.loc[df["Gene"].isin(keep_genes)].copy()

    step = FilterStep(
        name=f"min_{min_guides}_guides_per_gene",
        rows_in=len(df),
        rows_out=len(kept),
        genes_in=df["Gene"].nunique(),
        genes_out=kept["Gene"].nunique(),
    )
    step.log()
    return kept, step


def _to_long(df: pd.DataFrame) -> pd.DataFrame:
    """Reshape to one row per (gene, guide, treatment, replicate) observation."""
    parts = []
    for arm, treatment in ((ARM_E2, 0), (ARM_E2_4OHT, 1)):
        for rep in REPLICATES:
            col = f"log2norm_{arm}__{rep}"
            parts.append(
                pd.DataFrame(
                    {
                        "Gene": df["Gene"].to_numpy(),
                        "sgRNA": df["sgRNA"].to_numpy(),
                        "treatment": treatment,
                        "log2norm": df[col].to_numpy(),
                    }
                )
            )
    return pd.concat(parts, ignore_index=True)


def fit_gene_treatment_effects(long_df: pd.DataFrame) -> pd.DataFrame:
    """Fit a guide-blocked linear model per gene: log2norm ~ guide + treatment.

    Returns one row per gene with the raw (pre-moderation) treatment
    coefficient, its OLS residual variance, residual degrees of freedom,
    and the (X'X)^-1 diagonal entry for the treatment coefficient (used to
    convert a moderated variance into a moderated standard error later).

    Blocking on guide removes between-guide baseline differences from the
    residual before the treatment effect's variance is estimated, so the
    treatment coefficient is not a plain before/after difference of means:
    guide identity is fit as a nuisance parameter in the same model.
    """
    rows = []
    for gene, g in long_df.groupby("Gene", sort=False):
        guides = g["sgRNA"].to_numpy()
        uniq_guides = np.unique(guides)
        n_guides = len(uniq_guides)
        guide_index = {guide_id: i for i, guide_id in enumerate(uniq_guides)}

        n_obs = len(g)
        design = np.zeros((n_obs, n_guides + 1))
        treatment = g["treatment"].to_numpy()
        for row_i, guide_id in enumerate(guides):
            design[row_i, guide_index[guide_id]] = 1.0
        design[:, -1] = treatment

        y = g["log2norm"].to_numpy()
        beta_hat, _, rank, _ = np.linalg.lstsq(design, y, rcond=None)
        residual_df = n_obs - rank
        if residual_df <= 0:
            logger.warning(
                "fit_gene_treatment_effects: gene %s has no residual df (n_obs=%d, rank=%d); skipped",
                gene,
                n_obs,
                rank,
            )
            continue

        residuals = y - design @ beta_hat
        s2 = float(np.sum(residuals**2) / residual_df)
        xtx_inv = np.linalg.pinv(design.T @ design)
        c_tt = float(xtx_inv[-1, -1])

        rows.append(
            {
                "gene": gene,
                "n_guides": n_guides,
                "beta": float(beta_hat[-1]),
                "s2": s2,
                "c_tt": c_tt,
                "df_resid": residual_df,
            }
        )

    return pd.DataFrame(rows)


def _invert_trigamma(target: float) -> float:
    """Solve trigamma(x) = target for x > 0 (trigamma is monotone decreasing)."""
    return brentq(lambda x: polygamma(1, x) - target, 1e-8, 1e8)


def squeeze_variance(s2: np.ndarray, df_resid: np.ndarray) -> tuple[np.ndarray, float, float]:
    """Empirical-Bayes shrinkage of per-gene residual variances.

    Moment estimator of the prior degrees of freedom (d0) and prior
    variance (s0^2) from Smyth (2004), "Linear Models and Empirical Bayes
    Methods for Assessing Differential Expression in Microarray
    Experiments" — the method underlying ``limma::squeezeVar``/``eBayes``.
    Borrows strength across genes so that per-gene variance estimates
    (each from a handful of guides/replicates) are stabilised toward a
    common prior before being used for significance testing.
    """
    s2 = np.asarray(s2, dtype=float)
    df_resid = np.asarray(df_resid, dtype=float)

    if len(s2) < 2:
        raise ValueError("squeeze_variance needs at least 2 gene fits to estimate a prior")
    if not (np.isfinite(s2).all() and np.isfinite(df_resid).all()):
        raise ValueError("non-finite residual variance or degrees of freedom in gene fits")
    if (s2 <= 0).any():
        raise ValueError(
            "non-positive residual variance in gene fits (exact-zero variance is not "
            "handled by this moment estimator)"
        )
    if (df_resid <= 0).any():
        raise ValueError("non-positive residual degrees of freedom in gene fits")

    e = np.log(s2) - digamma(df_resid / 2) + np.log(df_resid / 2)
    ebar = float(np.mean(e))
    n_genes = len(e)
    sample_var_e = float(np.sum((e - ebar) ** 2) / (n_genes - 1))
    mean_trigamma_df = float(np.mean(polygamma(1, df_resid / 2)))
    target = sample_var_e - mean_trigamma_df

    if target <= 0:
        # No detectable between-gene variance heterogeneity beyond sampling
        # noise: fully pool toward a single common variance.
        d0 = np.inf
        s0_sq = float(np.exp(ebar))
        s2_post = np.full_like(s2, s0_sq)
    else:
        half_d0 = _invert_trigamma(target)
        d0 = 2 * half_d0
        s0_sq = float(np.exp(ebar + digamma(d0 / 2) - np.log(d0 / 2)))
        s2_post = (d0 * s0_sq + df_resid * s2) / (d0 + df_resid)

    logger.info("squeeze_variance: prior df d0=%.3f, prior variance s0^2=%.5f", d0, s0_sq)
    return s2_post, d0, s0_sq


def compute_effect_table(gene_fits: pd.DataFrame) -> pd.DataFrame:
    """Moderate variances, compute moderated t-stats, p-values, and BH-FDR."""
    s2_post, d0, _ = squeeze_variance(
        gene_fits["s2"].to_numpy(), gene_fits["df_resid"].to_numpy()
    )
    se_mod = np.sqrt(s2_post * gene_fits["c_tt"].to_numpy())
    beta = gene_fits["beta"].to_numpy()
    t_mod = beta / se_mod
    df_mod = gene_fits["df_resid"].to_numpy() + d0
    p_value = 2 * stats.t.sf(np.abs(t_mod), df=df_mod)
    _, fdr, _, _ = multipletests(p_value, method="fdr_bh")

    out = gene_fits[["gene", "n_guides"]].copy()
    out["effect_size"] = beta
    out["se"] = se_mod
    out["p_value"] = p_value
    out["fdr"] = fdr
    return out


def _load_config(config_path: str | Path) -> dict:
    with open(config_path) as f:
        return yaml.safe_load(f)


def build_labels(config_path: str | Path = "config/config.yaml") -> pd.DataFrame:
    """Run the full label pipeline and write parquet + TSV outputs.

    Writes ``data/processed/labels.parquet`` (per-gene effect table) and
    ``results/tables/label_summary.tsv`` (row/gene counts at every
    filtering step), using paths from ``config_path``.
    """
    config = _load_config(config_path)
    raw_path = Path(config["data"]["raw"]["hany_data_s1"])
    processed_dir = Path(config["data"].get("processed_dir", "data/processed"))
    results_tables_dir = Path(config["data"].get("results_tables_dir", "results/tables"))

    counts = load_normalised_counts(raw_path)
    selected = select_arms(counts)

    modelled_sample_cols = [f"{arm}__{rep}" for arm in MODELLED_ARMS for rep in REPLICATES]
    logged = add_log2_normalised(selected, modelled_sample_cols)

    filtered_t0, step_t0 = filter_zero_at_t0(logged)
    filtered_genes, step_genes = filter_low_representation_genes(filtered_t0)

    long_df = _to_long(filtered_genes)
    gene_fits = fit_gene_treatment_effects(long_df)
    effect_table = compute_effect_table(gene_fits)

    step_fit = FilterStep(
        name="gene_model_fit",
        rows_in=len(filtered_genes),
        rows_out=len(filtered_genes[filtered_genes["Gene"].isin(effect_table["gene"])]),
        genes_in=filtered_genes["Gene"].nunique(),
        genes_out=len(effect_table),
    )
    step_fit.log()
    if step_fit.genes_lost:
        logger.warning(
            "build_labels: %d genes entered modelling but were not fit (no residual "
            "degrees of freedom) — see prior 'fit_gene_treatment_effects' warnings",
            step_fit.genes_lost,
        )

    processed_dir.mkdir(parents=True, exist_ok=True)
    labels_path = processed_dir / "labels.parquet"
    effect_table.to_parquet(labels_path, index=False)
    logger.info("build_labels: wrote %d gene rows to %s", len(effect_table), labels_path)

    summary_rows = [
        {
            "step": "normalised_input",
            "rows_in": None,
            "rows_out": len(counts),
            "rows_lost": None,
            "genes_in": None,
            "genes_out": counts["Gene"].nunique(),
        },
        step_t0.as_record(),
        step_genes.as_record(),
        step_fit.as_record(),
    ]
    summary = pd.DataFrame(summary_rows)
    results_tables_dir.mkdir(parents=True, exist_ok=True)
    summary_path = results_tables_dir / "label_summary.tsv"
    summary.to_csv(summary_path, sep="\t", index=False)
    logger.info("build_labels: wrote filtering summary to %s", summary_path)

    return effect_table


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    build_labels()
