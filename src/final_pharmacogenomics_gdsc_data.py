"""GDSC/CancerRxGene data loaders for the final USP34/VEZF1 pharmacogenomics
phase. Reads only already-downloaded local files (GDSC Release 8.5 raw
files + this project's already-verified DepMap Public 26Q1 Model.csv/
expression) -- no network calls at runtime, per project convention.

Provenance: see config/config.yaml's final_pharmacogenomics.gdsc section
and /ibex/scratch/aljaroaa/tamoxifen-data/gdsc/PROVENANCE.txt.

Response-metric direction (independently verified against the official
GDSC_Fitted_Data_Description.pdf, Sanger, v1.0.0, 21 Sep 2017, not
assumed): LOWER LN_IC50 and LOWER AUC both mean MORE SENSITIVE; HIGHER
means MORE RESISTANT.

This is NOT new candidate discovery -- USP34 and VEZF1 are the only two
genes analyzed, per the frozen lead/backup translational conclusion.
"""

from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd
import yaml

from src.independent_validation_depmap_data import load_expression, load_model

logger = logging.getLogger(__name__)

GENES = ["USP34", "VEZF1"]
MIN_N_FULL_BREAST = 15
MIN_N_ER_LUMINAL_EXPLORATORY = 3  # reported as exploratory, never as validation, per instruction


def load_config(config_path: str | Path = "config/config.yaml") -> dict:
    with open(config_path) as f:
        return yaml.safe_load(f)


def _gdsc_dir(cfg: dict) -> Path:
    return Path(cfg["data"]["raw"][cfg["final_pharmacogenomics"]["gdsc"]["raw_dir_key"]])


def load_gdsc_fitted_response(cfg: dict) -> pd.DataFrame:
    """Combined GDSC1 + GDSC2 fitted dose-response table, DATASET-tagged."""
    raw = cfg["final_pharmacogenomics"]["gdsc"]["raw"]
    d = _gdsc_dir(cfg)
    g1 = pd.read_excel(d / raw["gdsc1_fitted_dose_response"])
    g2 = pd.read_excel(d / raw["gdsc2_fitted_dose_response"])
    g1["DATASET"], g2["DATASET"] = "GDSC1", "GDSC2"
    out = pd.concat([g1, g2], ignore_index=True)
    logger.info(
        "load_gdsc_fitted_response: GDSC1=%d rows (%d lines, %d drugs), GDSC2=%d rows (%d lines, %d drugs)",
        len(g1), g1["SANGER_MODEL_ID"].nunique(), g1["DRUG_ID"].nunique(),
        len(g2), g2["SANGER_MODEL_ID"].nunique(), g2["DRUG_ID"].nunique(),
    )
    return out


def load_gdsc_compounds(cfg: dict) -> pd.DataFrame:
    raw = cfg["final_pharmacogenomics"]["gdsc"]["raw"]
    return pd.read_csv(_gdsc_dir(cfg) / raw["screened_compounds"])


def load_gdsc_cell_line_details(cfg: dict) -> dict[str, pd.DataFrame]:
    raw = cfg["final_pharmacogenomics"]["gdsc"]["raw"]
    return pd.read_excel(_gdsc_dir(cfg) / raw["cell_lines_details"], sheet_name=None)


def build_breast_expression_joined(cfg: dict) -> pd.DataFrame:
    """GDSC breast-cancer (TCGA_DESC=='BRCA') fitted-response rows, joined
    to DepMap Public 26Q1 USP34/VEZF1 expression via SangerModelID -- an
    exact ID-based join (both files carry this column natively), never
    fuzzy cell-line-name matching. is_breast/is_er_luminal come directly
    from this project's own already-verified DepMap Model.csv
    classification (independent_validation_depmap_data.load_model),
    reused unchanged, not redefined here.
    """
    gdsc = load_gdsc_fitted_response(cfg)
    brca = gdsc[gdsc["TCGA_DESC"] == "BRCA"].copy()
    n_before = brca["SANGER_MODEL_ID"].nunique()

    model = load_model(cfg, "26Q1")
    expr = load_expression(cfg, "26Q1", GENES)
    expr_joined = expr.join(model[["SangerModelID", "COSMICID", "is_breast", "is_er_luminal"]], how="inner")
    expr_lookup = expr_joined.reset_index().set_index("SangerModelID")

    out = brca.merge(
        expr_lookup[GENES + ["is_breast", "is_er_luminal", "ModelID"]],
        left_on="SANGER_MODEL_ID", right_index=True, how="inner",
    )
    n_after = out["SANGER_MODEL_ID"].nunique()
    n_er_luminal = out.drop_duplicates("SANGER_MODEL_ID")["is_er_luminal"].sum()
    logger.info(
        "build_breast_expression_joined: %d GDSC breast lines -> %d matched to DepMap 26Q1 expression "
        "(%d lost to no SangerModelID match); %d of those are ER+/luminal by the project's own DepMap classification",
        n_before, n_after, n_before - n_after, n_er_luminal,
    )
    assert (out["is_breast"]).all(), "join must only ever include DepMap-confirmed breast lines"
    return out
