"""Independent-validation Part 2+3: TCGA-BRCA cohort construction and
candidate-gene expression comparisons.

Reads the already-downloaded TCGA-BRCA expression/clinical files via
src/independent_validation_tcga_data.py. Computes, for USP34/VEZF1/EML5/
CITED2: the overall expression distribution, ER+ vs ER- (primary clinical
IHC receptor status), PAM50 Luminal A vs Luminal B, and tumor vs adjacent
normal (paired where a matched-patient pair exists, unpaired otherwise,
reported as separate labeled rows per the preregistered spec). All
group comparisons use Welch's t-test on log2(TPM+1) values (already
log-scale, so a t-test on this scale is standard) with Cohen's d effect
size and a 95% CI on the mean difference; multiple-testing correction
(Benjamini-Hochberg) is applied across all rows in this table jointly.
Differential expression alone is never described as a mechanism.
"""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

from src.independent_validation_tcga_data import build_cohort_table, load_config, load_expression

logger = logging.getLogger(__name__)

CANDIDATES = ["USP34", "VEZF1", "EML5", "CITED2"]
OUT_TABLE = Path("results/tables/independent_validation/TCGA_candidate_expression.tsv")


def _welch_row(candidate: str, comparison: str, group_a_name: str, group_b_name: str, a: pd.Series, b: pd.Series, paired: bool = False, notes: str = "") -> dict:
    a = a.dropna()
    b = b.dropna()
    n_a, n_b = len(a), len(b)
    if n_a < 2 or n_b < 2:
        return dict(
            candidate=candidate, comparison=comparison, group_a=group_a_name, group_b=group_b_name,
            n_a=n_a, n_b=n_b, mean_a=np.nan, mean_b=np.nan, mean_diff=np.nan, ci_low=np.nan, ci_high=np.nan,
            cohens_d=np.nan, p_value=np.nan, test_used="none", paired=paired,
            notes=(notes + "; " if notes else "") + "insufficient N for a defensible test",
        )
    if paired:
        diff = a.values - b.values
        t_stat, p = stats.ttest_rel(a.values, b.values)
        mean_diff = float(np.mean(diff))
        se = stats.sem(diff)
        ci = stats.t.interval(0.95, len(diff) - 1, loc=mean_diff, scale=se) if se > 0 else (mean_diff, mean_diff)
        pooled_sd = float(np.std(diff, ddof=1))
        d = mean_diff / pooled_sd if pooled_sd > 0 else np.nan
        test_used = "paired t-test"
    else:
        t_stat, p = stats.ttest_ind(a.values, b.values, equal_var=False)
        mean_diff = float(a.mean() - b.mean())
        se = np.sqrt(a.var(ddof=1) / n_a + b.var(ddof=1) / n_b)
        dof = (a.var(ddof=1) / n_a + b.var(ddof=1) / n_b) ** 2 / (
            (a.var(ddof=1) / n_a) ** 2 / (n_a - 1) + (b.var(ddof=1) / n_b) ** 2 / (n_b - 1)
        )
        ci = stats.t.interval(0.95, dof, loc=mean_diff, scale=se) if se > 0 else (mean_diff, mean_diff)
        pooled_sd = np.sqrt(((n_a - 1) * a.var(ddof=1) + (n_b - 1) * b.var(ddof=1)) / (n_a + n_b - 2))
        d = mean_diff / pooled_sd if pooled_sd > 0 else np.nan
        test_used = "Welch's t-test"
    return dict(
        candidate=candidate, comparison=comparison, group_a=group_a_name, group_b=group_b_name,
        n_a=n_a, n_b=n_b, mean_a=float(a.mean()), mean_b=float(b.mean()), mean_diff=mean_diff,
        ci_low=float(ci[0]), ci_high=float(ci[1]), cohens_d=float(d), p_value=float(p),
        test_used=test_used, paired=paired, notes=notes,
    )


def build_expression_table(cfg: dict) -> pd.DataFrame:
    expr = load_expression(cfg, genes=CANDIDATES)
    cohort = build_cohort_table(cfg)
    df = cohort.join(expr)

    rows = []
    for candidate in CANDIDATES:
        vals = df.loc[df["is_primary_tumor"], candidate]
        rows.append(dict(
            candidate=candidate, comparison="distribution_all_primary_tumors", group_a="all primary tumors",
            group_b="", n_a=int(vals.notna().sum()), n_b=np.nan, mean_a=float(vals.mean()), mean_b=np.nan,
            mean_diff=np.nan, ci_low=np.nan, ci_high=np.nan, cohens_d=np.nan, p_value=np.nan,
            test_used="descriptive only", paired=False,
            notes=f"median={vals.median():.3f}, sd={vals.std():.3f}, IQR=[{vals.quantile(.25):.3f},{vals.quantile(.75):.3f}]",
        ))

        primary = df.loc[df["is_primary_tumor"]]
        er_pos = primary.loc[primary["ER_STATUS"] == "Positive", candidate]
        er_neg = primary.loc[primary["ER_STATUS"] == "Negative", candidate]
        rows.append(_welch_row(candidate, "ER+ vs ER- (clinical IHC)", "ER+", "ER-", er_pos, er_neg))

        luma = primary.loc[primary["PAM50_SUBTYPE"] == "Luminal A", candidate]
        lumb = primary.loc[primary["PAM50_SUBTYPE"] == "Luminal B", candidate]
        rows.append(_welch_row(candidate, "PAM50 Luminal A vs Luminal B", "Luminal A", "Luminal B", luma, lumb, notes="PAM50 calls are a molecular-subtype proxy, not clinical ER-IHC status -- reported separately from the ER+/ER- row above, never substituted for it"))

        tumor = df.loc[df["is_primary_tumor"], ["patient_barcode", candidate]].set_index("patient_barcode")[candidate]
        normal = df.loc[df["is_normal"], ["patient_barcode", candidate]].set_index("patient_barcode")[candidate]
        paired_patients = tumor.index.intersection(normal.index)
        paired_patients = paired_patients[~paired_patients.duplicated()] if isinstance(paired_patients, pd.Index) else paired_patients
        tumor_p = tumor.loc[tumor.index.isin(paired_patients)]
        tumor_p = tumor_p[~tumor_p.index.duplicated(keep="first")].reindex(paired_patients)
        normal_p = normal[~normal.index.duplicated(keep="first")].reindex(paired_patients)
        rows.append(_welch_row(candidate, "tumor_vs_normal_PAIRED", "primary tumor", "matched solid tissue normal", tumor_p, normal_p, paired=True))

        unmatched_tumor = tumor[~tumor.index.duplicated(keep="first")]
        unmatched_normal = normal[~normal.index.duplicated(keep="first")]
        rows.append(_welch_row(candidate, "tumor_vs_normal_UNPAIRED_descriptive", "all primary tumors", "all solid tissue normal", unmatched_tumor, unmatched_normal, notes="unpaired/descriptive -- redundant with the PAIRED row above, deliberately excluded from the FDR-correction family below to avoid diluting it with a non-independent duplicate test"))

    out = pd.DataFrame(rows)
    # the UNPAIRED_descriptive row is a redundant re-statement of the PAIRED
    # row on the same candidate/comparison and must not enter the same
    # correction family as it (would double-count that hypothesis test)
    testable = out["p_value"].notna() & (out["comparison"] != "tumor_vs_normal_UNPAIRED_descriptive")
    out["fdr"] = np.nan
    out.loc[testable, "fdr"] = multipletests(out.loc[testable, "p_value"], method="fdr_bh")[1]
    logger.info("build_expression_table: %d rows (%d in the FDR-correction family, joint BH-FDR)", len(out), int(testable.sum()))
    return out


def run(config_path: str = "config/config.yaml", out_table: Path = OUT_TABLE) -> pd.DataFrame:
    cfg = load_config(config_path)
    out = build_expression_table(cfg)
    out_table.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(out_table, sep="\t", index=False)
    logger.info("wrote %s (%d rows)", out_table, len(out))
    return out


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run()
