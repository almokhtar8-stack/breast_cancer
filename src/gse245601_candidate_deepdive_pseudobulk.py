"""GSE245601 candidate deep-dive Phases 4, 6, 7: patient-level pseudobulk
aggregation from the per-cell table. The "library size" for any cell
grouping is the sum of each cell's own `nCount_RNA` (total UMI count
across all genes for that cell, already frozen) -- verified to reproduce
the frozen Track A `total_library_size` column exactly (per-sample sums
match to the integer), so aggregating raw candidate-gene counts over any
subset of cells and normalizing by the summed `nCount_RNA` is the same
valid count-based pseudobulk representation the frozen Track A/B
pipeline itself uses, generalized to arbitrary patient x condition x
malignancy groupings (including combinations, e.g. individual patients'
malignant compartment, that the frozen Track B eligibility rule excludes
from formal inference but which are still valid to show descriptively).

This module produces DESCRIPTIVE pseudobulk points only. The frozen
edgeR log2FC/FDR (Track A across all 10 patients; Track B across the
3 eligible patients only) remain the sole inferential statistics used
anywhere in this deep-dive -- no new p-value is computed from any
grouping built here.
"""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from src.gse245601_candidate_deepdive_data import GENES, load_per_cell_table

logger = logging.getLogger(__name__)


def _load_config(config_path: str | Path) -> dict:
    with open(config_path) as f:
        return yaml.safe_load(f)


def build_group_pseudobulk(per_cell: pd.DataFrame, genes: list[str], group_cols: list[str]) -> pd.DataFrame:
    """One row per group x gene: n_cells, raw UMI sum for that gene,
    summed nCount_RNA (pseudobulk library size), log2(CPM+1)."""
    rows = []
    for keys, grp in per_cell.groupby(group_cols, observed=True):
        keys = keys if isinstance(keys, tuple) else (keys,)
        lib_size = int(grp["nCount_RNA"].sum())
        n_cells = len(grp)
        for gene in genes:
            raw = int(grp[f"{gene}_raw_count"].sum())
            norm = np.log2(raw / lib_size * 1e6 + 1) if lib_size > 0 else float("nan")
            row = dict(zip(group_cols, keys))
            row.update({"gene": gene, "n_cells": n_cells, "raw_umi_gene": raw, "library_size": lib_size, "normalized_expression": norm})
            rows.append(row)
    out = pd.DataFrame(rows)
    logger.info("build_group_pseudobulk(%s): %d group x gene rows", group_cols, len(out))
    return out


def build_patient_all_epithelial_pseudobulk(per_cell: pd.DataFrame, genes: list[str]) -> pd.DataFrame:
    pb = build_group_pseudobulk(per_cell, genes, ["patient", "condition"])
    return _add_direction(pb, genes)


def build_patient_malignant_pseudobulk(per_cell: pd.DataFrame, genes: list[str], min_cells: int) -> pd.DataFrame:
    """A patient x condition combination with ZERO malignant cells (e.g.
    Tumor_04 x Control) is a real, informative data point -- it must
    appear as `n_cells=0`, `low_cell_count_warning=True`, not silently
    vanish from `groupby` and later show up as NaN after pivoting (which
    would make the missing-arm case indistinguishable from a
    join/pivot bug rather than "zero malignant cells sampled"). Fixed
    after the Phase 25 Codex review found this exact gap."""
    malignant = per_cell.loc[per_cell["malignancy_status"] == "malignant"]
    pb = build_group_pseudobulk(malignant, genes, ["patient", "condition"])

    all_patients = sorted(per_cell["patient"].unique())
    all_conditions = sorted(per_cell["condition"].unique())
    full_index = pd.MultiIndex.from_product([all_patients, all_conditions, genes], names=["patient", "condition", "gene"])
    pb = pb.set_index(["patient", "condition", "gene"]).reindex(full_index).reset_index()
    for col in ("n_cells", "raw_umi_gene", "library_size"):
        pb[col] = pb[col].fillna(0).astype(int)
    # normalized_expression stays NaN when library_size==0 (nothing to normalize -- never fabricated as 0)

    pb["low_cell_count_warning"] = pb["n_cells"] < min_cells
    return _add_direction(pb, genes)


def _add_direction(pb: pd.DataFrame, genes: list[str]) -> pd.DataFrame:
    """Pivots Control/Tamoxifen onto one row per patient x gene and adds
    the exact-normalized-value direction (increase/decrease/equal) --
    never a re-run inferential test, purely descriptive."""
    # dropna=False: pivot_table's default silently drops a "condition" column that is
    # entirely NaN across the input (e.g. a gene/subset where NO patient has data for
    # that arm) -- rare with the full 10-patient cohort but a real edge case for any
    # smaller subset, and MUST NOT silently vanish the Control/Tamoxifen column itself
    wide = pb.pivot_table(index=[c for c in pb.columns if c not in ("condition", "n_cells", "raw_umi_gene", "library_size", "normalized_expression", "low_cell_count_warning")], columns="condition", values=["n_cells", "raw_umi_gene", "library_size", "normalized_expression"] + (["low_cell_count_warning"] if "low_cell_count_warning" in pb.columns else []), dropna=False)
    wide.columns = [f"{a}_{b}" for a, b in wide.columns]
    wide = wide.reset_index()
    if "normalized_expression_Control" in wide.columns and "normalized_expression_Tamoxifen" in wide.columns:
        delta = wide["normalized_expression_Tamoxifen"] - wide["normalized_expression_Control"]
        wide["log_fold_change_descriptive"] = delta
        # np.select's `default` fires for NaN deltas too (NaN>0 and NaN<0 are both False) --
        # an arm with 0 cells (delta=NaN) must be labeled "not_comparable", never "equal"
        # (Phase 25 Codex review, second pass)
        wide["direction_control_to_tam"] = np.select([delta > 0, delta < 0, delta == 0], ["increase", "decrease", "equal"], default="not_comparable")
    return wide


def build_malignancy_condition_patient_summary(per_cell: pd.DataFrame, genes: list[str], include_pooled_all_epithelial: bool = True) -> pd.DataFrame:
    """Phases 7 and 9 (same table serves both -- Phase 9's five metrics
    are a superset of Phase 7's): per patient x malignancy class x
    condition, pseudobulk (count-based) plus cell-level descriptive
    prevalence/intensity statistics. `malignancy_status` also includes a
    pooled `all_epithelial` pseudo-class (malignant + non-malignant
    combined) when `include_pooled_all_epithelial=True`, so the same
    table can answer "prevalence vs intensity" at the whole-epithelial
    level without a second near-duplicate table."""
    frames = [per_cell]
    if include_pooled_all_epithelial:
        pooled = per_cell.copy()
        pooled["malignancy_status"] = "all_epithelial"
        frames.append(pooled)
    combined = pd.concat(frames, ignore_index=True)

    rows = []
    for (patient, malignancy, condition), grp in combined.groupby(["patient", "malignancy_status", "condition"], observed=True):
        lib_size = int(grp["nCount_RNA"].sum())
        n_cells = len(grp)
        for gene in genes:
            raw_col = grp[f"{gene}_raw_count"]
            log_norm_col = grp[f"{gene}_log_norm"]
            positive = raw_col > 0
            raw_sum = int(raw_col.sum())
            rows.append(
                {
                    "gene": gene, "patient": patient, "malignancy_status": malignancy, "condition": condition,
                    "n_cells": n_cells, "pseudobulk_raw_umi": raw_sum, "pseudobulk_library_size": lib_size,
                    "pseudobulk_normalized_expression": np.log2(raw_sum / lib_size * 1e6 + 1) if lib_size > 0 else float("nan"),
                    "fraction_expressing": float(positive.mean()) if n_cells > 0 else float("nan"),
                    "mean_normalized_per_cell": float(log_norm_col.mean()) if n_cells > 0 else float("nan"),
                    "median_normalized_per_cell": float(log_norm_col.median()) if n_cells > 0 else float("nan"),
                    "mean_normalized_positive_cells_only": float(log_norm_col.loc[positive].mean()) if positive.sum() > 0 else float("nan"),
                    "median_normalized_positive_cells_only": float(log_norm_col.loc[positive].median()) if positive.sum() > 0 else float("nan"),
                }
            )
    out = pd.DataFrame(rows)

    # a patient x malignancy-class x condition combination with ZERO cells (e.g. Tumor_04's
    # malignant Control arm) is a real, informative data point -- it must appear as n_cells=0,
    # not silently vanish from groupby (same class of gap the Phase 25 Codex review found in
    # build_patient_malignant_pseudobulk, fixed here too for consistency)
    all_malignancy_classes = sorted(combined["malignancy_status"].unique())
    full_index = pd.MultiIndex.from_product([sorted(combined["patient"].unique()), all_malignancy_classes, sorted(combined["condition"].unique()), genes], names=["patient", "malignancy_status", "condition", "gene"])
    out = out.set_index(["patient", "malignancy_status", "condition", "gene"]).reindex(full_index).reset_index()
    out["n_cells"] = out["n_cells"].fillna(0).astype(int)
    for col in ("pseudobulk_raw_umi", "pseudobulk_library_size"):
        out[col] = out[col].fillna(0)
    # all normalized/prevalence/intensity columns correctly stay NaN when n_cells==0 -- never fabricated

    logger.info("build_malignancy_condition_patient_summary: %d rows", len(out))
    return out


def run_pseudobulk(config_path: str | Path = "config/config.yaml") -> dict[str, pd.DataFrame]:
    config = _load_config(config_path)
    cfg = config["gse245601_candidate_deepdive"]
    out_dir = Path(cfg["output"]["tables_dir"])
    out_dir.mkdir(parents=True, exist_ok=True)
    min_cells = cfg["malignant_cell_min_count"]

    per_cell = load_per_cell_table(config)

    all_epi = build_patient_all_epithelial_pseudobulk(per_cell, GENES)
    all_epi.to_csv(out_dir / "patient_all_epithelial_pseudobulk.tsv", sep="\t", index=False)

    malignant = build_patient_malignant_pseudobulk(per_cell, GENES, min_cells)
    malignant.to_csv(out_dir / "patient_malignant_pseudobulk.tsv", sep="\t", index=False)

    malignancy_condition = build_malignancy_condition_patient_summary(per_cell, GENES)
    malignancy_condition.to_csv(out_dir / "malignancy_condition_patient_summary.tsv", sep="\t", index=False)
    # Phase 9 wants this exact filename; same table (see docstring) -- avoids a near-duplicate builder
    malignancy_condition.to_csv(out_dir / "expression_prevalence_intensity.tsv", sep="\t", index=False)

    return {"all_epithelial": all_epi, "malignant": malignant, "malignancy_condition": malignancy_condition}


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run_pseudobulk()
