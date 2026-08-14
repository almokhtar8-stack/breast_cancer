"""Independent-validation Part 6: TCGA-BRCA clinical (overall-survival) Cox
models for the four frozen candidates.

Expression enters each model as a per-cohort standardized (z-scored)
continuous covariate -- never dichotomized as the primary analysis, per
the preregistered spec. Two models per candidate/cohort: (A) univariable,
(B) adjusted for AGE and AJCC pathologic stage (mapped to an ordinal I-IV
scale; "Stage X" / missing stage is dropped from the adjusted model only,
never imputed). Repeated in the full primary-tumor cohort and, if the
event count still clears the preregistered minimum, the ER+ subset.
Proportional-hazards assumptions are checked (lifelines' scaled
Schoenfeld-residual test) and reported, not silently ignored when they
fail. FDR correction (Benjamini-Hochberg) is applied across the four
candidates within each model/cohort stratum. A survival association is
never described as a tamoxifen-response result, and a significant HR is
never described as prognostic independence beyond what the model
actually supports.
"""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from lifelines.statistics import proportional_hazard_test
from statsmodels.stats.multitest import multipletests

from src.independent_validation_tcga_data import build_cohort_table, load_config, load_expression

logger = logging.getLogger(__name__)

OUT_TABLE = Path("results/tables/independent_validation/TCGA_candidate_clinical.tsv")

STAGE_ORDINAL = {
    "STAGE I": 1, "STAGE IA": 1, "STAGE IB": 1,
    "STAGE II": 2, "STAGE IIA": 2, "STAGE IIB": 2,
    "STAGE III": 3, "STAGE IIIA": 3, "STAGE IIIB": 3, "STAGE IIIC": 3,
    "STAGE IV": 4,
}


def _build_survival_frame(cfg: dict) -> pd.DataFrame:
    cohort = build_cohort_table(cfg)
    expr = load_expression(cfg, genes=cfg["independent_validation"]["candidates"])
    df = cohort.loc[cohort["is_primary_tumor"]].join(expr)

    df["event"] = np.select([df["OS_STATUS"] == "1:DECEASED", df["OS_STATUS"] == "0:LIVING"], [1, 0], default=np.nan)
    df["duration"] = pd.to_numeric(df["OS_MONTHS"], errors="coerce")
    df["age_numeric"] = pd.to_numeric(df["AGE"], errors="coerce")
    df["stage_ordinal"] = df["AJCC_PATHOLOGIC_TUMOR_STAGE"].map(STAGE_ORDINAL)
    n_before = len(df)
    df = df.loc[df["event"].notna() & df["duration"].notna() & (df["duration"] > 0)]
    logger.info("_build_survival_frame: %d primary tumors -> %d with usable OS_STATUS/OS_MONTHS (%d dropped: missing/invalid follow-up)", n_before, len(df), n_before - len(df))
    return df


def _fit_one(df: pd.DataFrame, candidate: str, adjusted: bool, cohort_label: str, min_events: int) -> dict:
    cols = [candidate, "duration", "event"] + (["age_numeric", "stage_ordinal"] if adjusted else [])
    sub = df[cols].dropna().copy()
    n = len(sub)
    n_events = int(sub["event"].sum())
    base = dict(candidate=candidate, cohort=cohort_label, model=("adjusted_age_stage" if adjusted else "univariable"), n=n, n_events=n_events)
    if n_events < min_events:
        return {**base, "hr_per_sd": np.nan, "ci_low": np.nan, "ci_high": np.nan, "p_value": np.nan, "ph_assumption_p": np.nan, "notes": f"fewer than {min_events} events -- Cox model not fit (underpowered)"}

    sub[candidate] = (sub[candidate] - sub[candidate].mean()) / sub[candidate].std(ddof=0)
    cph = CoxPHFitter()
    try:
        cph.fit(sub, duration_col="duration", event_col="event")
    except Exception as e:
        return {**base, "hr_per_sd": np.nan, "ci_low": np.nan, "ci_high": np.nan, "p_value": np.nan, "ph_assumption_p": np.nan, "notes": f"Cox fit failed: {e}"}

    row = cph.summary.loc[candidate]
    try:
        ph_test = proportional_hazard_test(cph, sub, time_transform="rank")
        ph_p = float(ph_test.summary.loc[candidate, "p"])
        ph_note = "PH assumption holds (p>0.05)" if ph_p > 0.05 else "PH assumption VIOLATED (p<=0.05) -- HR interpreted as a time-averaged effect only"
    except Exception as e:
        ph_p, ph_note = np.nan, f"PH test failed: {e}"

    return {
        **base,
        "hr_per_sd": float(np.exp(row["coef"])),
        "ci_low": float(np.exp(row["coef lower 95%"])),
        "ci_high": float(np.exp(row["coef upper 95%"])),
        "p_value": float(row["p"]),
        "ph_assumption_p": ph_p,
        "notes": ph_note,
    }


def build_clinical_table(cfg: dict) -> pd.DataFrame:
    df = _build_survival_frame(cfg)
    candidates = cfg["independent_validation"]["candidates"]
    min_events = cfg["independent_validation"]["tcga"]["thresholds"]["min_events_cox"]

    cohorts = {"all_primary_tumors": df, "ER_positive": df.loc[df["ER_STATUS"] == "Positive"]}
    rows = []
    for cohort_label, cdf in cohorts.items():
        for candidate in candidates:
            rows.append(_fit_one(cdf, candidate, adjusted=False, cohort_label=cohort_label, min_events=min_events))
            rows.append(_fit_one(cdf, candidate, adjusted=True, cohort_label=cohort_label, min_events=min_events))

    out = pd.DataFrame(rows)
    out["fdr"] = np.nan
    for (cohort_label, model), grp in out.groupby(["cohort", "model"]):
        testable = grp["p_value"].notna()
        if testable.sum() > 0:
            out.loc[grp.index[testable], "fdr"] = multipletests(grp.loc[testable, "p_value"], method="fdr_bh")[1]
    logger.info("build_clinical_table: %d rows", len(out))
    return out


def run(config_path: str = "config/config.yaml", out_table: Path = OUT_TABLE) -> pd.DataFrame:
    cfg = load_config(config_path)
    out = build_clinical_table(cfg)
    out_table.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(out_table, sep="\t", index=False)
    logger.info("wrote %s (%d rows)", out_table, len(out))
    return out


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run()
