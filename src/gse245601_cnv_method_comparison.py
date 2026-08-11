"""Quantitative InferCNV-vs-CopyKAT comparison for the 20 primary GSE245601
tumor samples (Tumor_01_Control .. Tumor_10_Tamoxifen), audit-only.

This module reruns nothing: it reads exactly three already-frozen,
committed tables --

- ``results/tables/gse245601_malignant_cell_labels.tsv`` (InferCNV-derived
  per-epithelial-cell ``primary_malignancy_label``, written by
  ``scripts/analysis/gse245601_05_infercnv_malignant.R``),
- ``results/tables/gse245601_copykat_sensitivity_labels.tsv`` (CopyKAT
  ``sensitivity_malignancy_label``, written by
  ``scripts/analysis/gse245601_06_copykat_sensitivity.R``),
- ``results/tables/gse245601_malignant_summary_per_sample.tsv`` (per-sample
  epithelial-cell totals, same script) --

and recombines them into two quantitative comparison tables. No InferCNV or
CopyKAT run is triggered anywhere in this module, no malignancy label,
parameter, or eligibility rule is changed, and no candidate gene (13
candidates + PAICS) is read, joined, or referenced anywhere.

Cells are matched strictly by ``sample_id`` + epithelial-cell barcode
(``cell_id``, already sample-prefixed and therefore globally unique -- see
``docs/CNV_METHOD_AUDIT.md`` Section 5.2). An inner join is used
deliberately: any epithelial cell present in one table but not the other
would indicate a provenance problem and must not be silently treated as a
comparison denominator; ``match_cells`` logs and raises if that ever
happens rather than padding/dropping silently.

Denominator convention (documented explicitly because two different counts
are both meaningful):

- ``n_matched`` -- epithelial cells with a label from *both* methods,
  regardless of whether CopyKAT's call was ``not_defined``. Equal to
  ``total_epithelial_cells`` for all 20 samples in the current frozen data
  (verified by ``match_cells``).
- ``n_compared`` -- the subset of ``n_matched`` cells where CopyKAT
  returned a defined call (``malignant`` or ``non-malignant epithelial``,
  i.e. excludes ``not_defined``). This is "cells successfully compared by
  both methods" throughout this module's output, and is exactly the
  denominator ``scripts/analysis/gse245601_06_copykat_sensitivity.R``
  itself used for the previously reported ~54.1% overall concordance
  figure -- so every percentage in ``cnv_method_comparison_by_sample.tsv``
  shares one consistent, reproducible denominator per sample, and the four
  2x2 count columns sum exactly to ``n_compared``.

No cutoff for "bad" or "outlier" samples is computed or applied anywhere in
this module -- it produces counts and percentages only. Identifying
outliers is left to descriptive inspection of the output tables.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402
import yaml  # noqa: E402

logger = logging.getLogger(__name__)

# Okabe-Ito colorblind-safe pair (Okabe & Ito 2008): blue for InferCNV,
# vermillion for CopyKAT. Fixed assignment, never cycled or re-derived.
COLOR_INFERCNV = "#0072B2"
COLOR_COPYKAT = "#D55E00"
COLOR_CONNECTOR = "#B0B0B0"

LABEL_MALIGNANT = "malignant"
LABEL_NONMALIGNANT = "non-malignant epithelial"
LABEL_NOT_DEFINED = "not_defined"

INFERCNV_LABELS = (LABEL_MALIGNANT, LABEL_NONMALIGNANT)
COPYKAT_LABELS = (LABEL_MALIGNANT, LABEL_NONMALIGNANT, LABEL_NOT_DEFINED)

REQUIRED_INFERCNV_COLUMNS = ("cell_id", "sample_id", "patient", "condition", "primary_malignancy_label")
REQUIRED_COPYKAT_COLUMNS = ("cell_id", "sample_id", "sensitivity_malignancy_label")
REQUIRED_SUMMARY_COLUMNS = ("sample_id", "n_epithelial_cells", "status")


def _load_config(config_path: str | Path) -> dict:
    with open(config_path) as f:
        return yaml.safe_load(f)


@dataclass(frozen=True)
class CnvMethodComparisonConfig:
    """Resolved, config-driven paths. No hardcoded paths."""

    infercnv_labels_tsv: Path
    infercnv_sample_summary_tsv: Path
    copykat_labels_tsv: Path
    expected_n_tumor_samples: int
    by_sample_tsv: Path
    by_tumor_tsv: Path
    figure_png: Path
    figure_pdf: Path

    @classmethod
    def from_config(cls, config: dict) -> "CnvMethodComparisonConfig":
        cmc = config["gse245601_cnv_method_comparison"]
        inputs = cmc["inputs"]
        out = cmc["output"]
        return cls(
            infercnv_labels_tsv=Path(inputs["infercnv_labels_tsv"]),
            infercnv_sample_summary_tsv=Path(inputs["infercnv_sample_summary_tsv"]),
            copykat_labels_tsv=Path(inputs["copykat_labels_tsv"]),
            expected_n_tumor_samples=int(config["gse245601"]["expected_n_tumor_samples"]),
            by_sample_tsv=Path(out["by_sample_tsv"]),
            by_tumor_tsv=Path(out["by_tumor_tsv"]),
            figure_png=Path(out["figure_png"]),
            figure_pdf=Path(out["figure_pdf"]),
        )


def load_infercnv_labels(path: str | Path) -> pd.DataFrame:
    """Load the frozen per-cell InferCNV-derived malignancy labels
    (read-only; InferCNV/the thresholding logic are not re-run here)."""
    df = pd.read_csv(path, sep="\t")
    missing = [c for c in REQUIRED_INFERCNV_COLUMNS if c not in df.columns]
    if missing:
        raise ValueError(f"InferCNV labels table missing required columns: {missing}")
    if df["cell_id"].isna().any() or df["sample_id"].isna().any():
        raise ValueError("InferCNV labels table contains null cell_id or sample_id values")
    if df["cell_id"].duplicated().any():
        raise ValueError("InferCNV labels table contains duplicate cell_id values")
    bad = set(df["primary_malignancy_label"].unique()) - set(INFERCNV_LABELS)
    if bad:
        raise ValueError(f"InferCNV labels table contains unexpected primary_malignancy_label values: {bad}")
    logger.info("load_infercnv_labels: read %d epithelial cells across %d samples", len(df), df["sample_id"].nunique())
    return df


def load_copykat_labels(path: str | Path) -> pd.DataFrame:
    """Load the frozen per-cell CopyKAT sensitivity labels (read-only;
    CopyKAT is not re-run here). ``sensitivity_malignancy_label`` values
    are exactly the ones ``gse245601_06_copykat_sensitivity.R`` already
    derived from CopyKAT's own ``copykat.pred`` (aneuploid->malignant,
    diploid->non-malignant epithelial, not.defined->not_defined) -- no
    relabeling happens in this module."""
    df = pd.read_csv(path, sep="\t")
    missing = [c for c in REQUIRED_COPYKAT_COLUMNS if c not in df.columns]
    if missing:
        raise ValueError(f"CopyKAT labels table missing required columns: {missing}")
    if df["cell_id"].isna().any() or df["sample_id"].isna().any():
        raise ValueError("CopyKAT labels table contains null cell_id or sample_id values")
    if df["cell_id"].duplicated().any():
        raise ValueError("CopyKAT labels table contains duplicate cell_id values")
    bad = set(df["sensitivity_malignancy_label"].unique()) - set(COPYKAT_LABELS)
    if bad:
        raise ValueError(f"CopyKAT labels table contains unexpected sensitivity_malignancy_label values: {bad}")
    logger.info("load_copykat_labels: read %d epithelial cells across %d samples", len(df), df["sample_id"].nunique())
    return df


def load_infercnv_sample_summary(path: str | Path) -> pd.DataFrame:
    """Load the frozen per-sample epithelial-cell totals (the
    ``total_epithelial_cells`` denominator context column)."""
    df = pd.read_csv(path, sep="\t")
    missing = [c for c in REQUIRED_SUMMARY_COLUMNS if c not in df.columns]
    if missing:
        raise ValueError(f"InferCNV sample summary table missing required columns: {missing}")
    if df["sample_id"].isna().any():
        raise ValueError("InferCNV sample summary table contains null sample_id values")
    if df["sample_id"].duplicated().any():
        raise ValueError("InferCNV sample summary table contains duplicate sample_id values")
    logger.info("load_infercnv_sample_summary: read %d samples", len(df))
    return df


def match_cells(infercnv_df: pd.DataFrame, copykat_df: pd.DataFrame) -> pd.DataFrame:
    """Strict sample_id + cell_id cell matching. The actual join key is
    ``cell_id`` alone (each loader already rejects duplicate cell_id
    values, and cell_id is sample-prefixed and therefore globally unique
    -- see ``docs/CNV_METHOD_AUDIT.md`` Section 5.2), but ``sample_id``
    agreement between the two tables is independently verified below
    rather than assumed, so the effective match is sample_id+cell_id, not
    cell_id alone. Any InferCNV epithelial cell without a matching CopyKAT
    row, or vice versa, is a data-provenance problem and is logged then
    raised on -- never silently dropped or padded."""
    infercnv_only = set(infercnv_df["cell_id"]) - set(copykat_df["cell_id"])
    copykat_only = set(copykat_df["cell_id"]) - set(infercnv_df["cell_id"])
    if infercnv_only or copykat_only:
        raise ValueError(
            f"cell_id sets do not match between InferCNV and CopyKAT label tables: "
            f"{len(infercnv_only)} cell(s) only in InferCNV, {len(copykat_only)} cell(s) only in CopyKAT"
        )

    merged = infercnv_df.merge(
        copykat_df[["cell_id", "sample_id", "sensitivity_malignancy_label"]],
        on="cell_id",
        how="inner",
        validate="one_to_one",
        suffixes=("", "_copykat"),
    )
    mismatched_sample = merged["sample_id"] != merged["sample_id_copykat"]
    if mismatched_sample.any():
        raise ValueError(
            f"{int(mismatched_sample.sum())} cell_id(s) map to a different sample_id in the CopyKAT table "
            f"than in the InferCNV table -- either a barcode collision across samples or a sample-label "
            f"inconsistency between the two tables"
        )
    merged = merged.drop(columns=["sample_id_copykat"])

    logger.info(
        "match_cells: %d InferCNV epithelial cells in, %d CopyKAT epithelial cells in, "
        "%d matched by sample_id+cell_id (0 lost on either side)",
        len(infercnv_df),
        len(copykat_df),
        len(merged),
    )
    return merged


VALID_CONDITIONS = ("Control", "Tamoxifen")


def _condition_pct_columns(sample_id: str) -> tuple[str, str]:
    patient, condition = sample_id.rsplit("_", 1)
    if condition not in VALID_CONDITIONS:
        raise ValueError(f"sample_id {sample_id!r} has unexpected condition {condition!r}, expected one of {VALID_CONDITIONS}")
    return patient, condition


def build_by_sample_table(matched_df: pd.DataFrame, summary_df: pd.DataFrame) -> pd.DataFrame:
    """One row per sample (Tumor_01_Control .. Tumor_10_Tamoxifen): counts
    and percentages of InferCNV malignant calls, CopyKAT aneuploid calls,
    and their agreement, all computed on the shared ``n_compared``
    denominator (CopyKAT-defined, matched cells -- see module docstring).
    ``total_epithelial_cells`` is reported alongside as context; it is not
    used as a percentage denominator here because CopyKAT could not return
    a defined call for every one of those cells."""
    rows = []
    for sample_id, sub in matched_df.groupby("sample_id"):
        n_matched = len(sub)
        defined = sub.loc[sub["sensitivity_malignancy_label"] != LABEL_NOT_DEFINED]
        n_compared = len(defined)
        n_not_defined = n_matched - n_compared

        ic = defined["primary_malignancy_label"]
        ck = defined["sensitivity_malignancy_label"]

        both_malignant = int(((ic == LABEL_MALIGNANT) & (ck == LABEL_MALIGNANT)).sum())
        both_nonmalignant = int(((ic == LABEL_NONMALIGNANT) & (ck == LABEL_NONMALIGNANT)).sum())
        ic_mal_ck_dip = int(((ic == LABEL_MALIGNANT) & (ck == LABEL_NONMALIGNANT)).sum())
        ic_non_ck_mal = int(((ic == LABEL_NONMALIGNANT) & (ck == LABEL_MALIGNANT)).sum())

        if both_malignant + both_nonmalignant + ic_mal_ck_dip + ic_non_ck_mal != n_compared:
            raise ValueError(f"{sample_id}: 2x2 breakdown does not sum to n_compared -- internal accounting error")

        infercnv_malignant_count = both_malignant + ic_mal_ck_dip
        copykat_aneuploid_count = both_malignant + ic_non_ck_mal
        agreement_count = both_malignant + both_nonmalignant

        total_epithelial = int(summary_df.loc[summary_df["sample_id"] == sample_id, "n_epithelial_cells"].iloc[0])
        if total_epithelial != n_matched:
            raise ValueError(
                f"{sample_id}: total_epithelial_cells ({total_epithelial}) != cells matched to a "
                f"CopyKAT row ({n_matched}) -- every InferCNV epithelial cell was expected to have a "
                f"CopyKAT row in the current frozen data"
            )

        patient, condition = _condition_pct_columns(sample_id)
        pct = lambda n: (100.0 * n / n_compared) if n_compared > 0 else float("nan")  # noqa: E731

        rows.append(
            {
                "sample_id": sample_id,
                "patient": patient,
                "condition": condition,
                "total_epithelial_cells": total_epithelial,
                "n_compared": n_compared,
                "infercnv_malignant_count": infercnv_malignant_count,
                "infercnv_malignant_pct": pct(infercnv_malignant_count),
                "copykat_aneuploid_count": copykat_aneuploid_count,
                "copykat_aneuploid_pct": pct(copykat_aneuploid_count),
                "agreement_count": agreement_count,
                "agreement_pct": pct(agreement_count),
                "both_malignant_aneuploid": both_malignant,
                "both_nonmalignant_diploid": both_nonmalignant,
                "infercnv_malignant_copykat_diploid": ic_mal_ck_dip,
                "infercnv_nonmalignant_copykat_aneuploid": ic_non_ck_mal,
                "copykat_not_defined_count": n_not_defined,
            }
        )

    out = pd.DataFrame(rows).sort_values("sample_id").reset_index(drop=True)
    n_compared_total = int(out["n_compared"].sum())
    if n_compared_total > 0:
        logger.info(
            "build_by_sample_table: %d samples; pooled agreement = %d/%d (%.4f%%)",
            len(out),
            out["agreement_count"].sum(),
            n_compared_total,
            100.0 * out["agreement_count"].sum() / n_compared_total,
        )
    else:
        logger.info("build_by_sample_table: %d samples; no cells were comparable in any sample", len(out))
    return out


def build_by_tumor_table(by_sample_df: pd.DataFrame) -> pd.DataFrame:
    """One row per tumor (Tumor_01 .. Tumor_10), pairing the Control and
    Tamoxifen sample rows from ``build_by_sample_table``. Differences are
    Tamoxifen minus Control (a positive value means the Tamoxifen-arm
    sample had a higher percentage than the Control-arm sample from the
    same patient)."""
    df = by_sample_df.copy()
    control = df.loc[df["condition"] == "Control"].set_index("patient")
    tamoxifen = df.loc[df["condition"] == "Tamoxifen"].set_index("patient")

    if control.index.duplicated().any() or tamoxifen.index.duplicated().any():
        raise ValueError(
            f"duplicate patient rows within a condition: "
            f"Control duplicates={sorted(set(control.index[control.index.duplicated()]))}, "
            f"Tamoxifen duplicates={sorted(set(tamoxifen.index[tamoxifen.index.duplicated()]))}"
        )

    missing_control = set(df["patient"]) - set(control.index)
    missing_tamoxifen = set(df["patient"]) - set(tamoxifen.index)
    if missing_control or missing_tamoxifen:
        raise ValueError(
            f"not every tumor has both a Control and a Tamoxifen row: "
            f"missing_control={missing_control}, missing_tamoxifen={missing_tamoxifen}"
        )

    rows = []
    for patient in sorted(control.index):
        c = control.loc[patient]
        t = tamoxifen.loc[patient]
        rows.append(
            {
                "tumor": patient,
                "control_n_compared": int(c["n_compared"]),
                "tamoxifen_n_compared": int(t["n_compared"]),
                "control_infercnv_malignant_pct": c["infercnv_malignant_pct"],
                "tamoxifen_infercnv_malignant_pct": t["infercnv_malignant_pct"],
                "infercnv_malignant_pct_diff_tam_minus_control": t["infercnv_malignant_pct"] - c["infercnv_malignant_pct"],
                "control_copykat_aneuploid_pct": c["copykat_aneuploid_pct"],
                "tamoxifen_copykat_aneuploid_pct": t["copykat_aneuploid_pct"],
                "copykat_aneuploid_pct_diff_tam_minus_control": t["copykat_aneuploid_pct"] - c["copykat_aneuploid_pct"],
                "control_agreement_pct": c["agreement_pct"],
                "tamoxifen_agreement_pct": t["agreement_pct"],
            }
        )

    out = pd.DataFrame(rows).sort_values("tumor").reset_index(drop=True)
    logger.info("build_by_tumor_table: %d tumors paired (Control + Tamoxifen each)", len(out))
    return out


def verify_overall_agreement(by_sample_df: pd.DataFrame, expected_pct: float, tolerance_pct: float = 0.05) -> float:
    """Recompute the pooled overall agreement (sum(agreement_count) /
    sum(n_compared)) from this module's own matched cells and confirm it
    reproduces the previously reported ~54.1% figure from
    ``gse245601_malignancy_concordance.tsv`` to within ``tolerance_pct``
    percentage points. Raises if it does not -- this check exists
    specifically to catch a silent denominator/matching divergence between
    the two modules."""
    n_compared_total = by_sample_df["n_compared"].sum()
    if n_compared_total == 0:
        raise ValueError("cannot verify overall agreement: no cells were comparable in any sample (n_compared sums to 0)")
    overall_pct = 100.0 * by_sample_df["agreement_count"].sum() / n_compared_total
    if abs(overall_pct - expected_pct) > tolerance_pct:
        raise ValueError(
            f"recomputed overall agreement {overall_pct:.4f}% does not match the previously reported "
            f"{expected_pct:.4f}% within {tolerance_pct} percentage points"
        )
    logger.info(
        "verify_overall_agreement: recomputed overall agreement %.4f%% matches previously reported %.4f%%",
        overall_pct,
        expected_pct,
    )
    return overall_pct


def build_figure(by_sample_df: pd.DataFrame, png_path: str | Path, pdf_path: str | Path) -> None:
    """One dumbbell (paired-dot) plot: InferCNV malignant % and CopyKAT
    aneuploid % for all 20 samples, generated only from
    ``build_by_sample_table``'s output -- no hand-typed numbers. Samples
    are ordered exactly as in the by-sample table (Tumor_01_Control ..
    Tumor_10_Tamoxifen); the horizontal gap between the two dots on a row
    is the disagreement this figure exists to surface."""
    df = by_sample_df.sort_values("sample_id", ascending=False).reset_index(drop=True)
    n = len(df)
    y = range(n)

    fig, ax = plt.subplots(figsize=(7.5, 0.32 * n + 1.5))

    for yi, row in zip(y, df.itertuples(index=False)):
        ax.plot(
            [row.infercnv_malignant_pct, row.copykat_aneuploid_pct],
            [yi, yi],
            color=COLOR_CONNECTOR,
            linewidth=1,
            zorder=1,
        )
    ax.scatter(
        df["infercnv_malignant_pct"], y, color=COLOR_INFERCNV, s=32, zorder=2, label="InferCNV malignant %"
    )
    ax.scatter(
        df["copykat_aneuploid_pct"], y, color=COLOR_COPYKAT, s=32, zorder=2, label="CopyKAT aneuploid %"
    )

    ax.set_yticks(list(y))
    ax.set_yticklabels(df["sample_id"], fontsize=8)
    ax.set_xlabel("% of compared epithelial cells")
    ax.set_xlim(-2, 102)
    ax.set_title("InferCNV malignant % vs CopyKAT aneuploid %, per sample (n_compared denominator)", fontsize=10, pad=28)
    ax.grid(axis="x", color="#E0E0E0", linewidth=0.6, zorder=0)
    for spine in ("top", "right", "left"):
        ax.spines[spine].set_visible(False)
    ax.tick_params(axis="y", length=0)
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=2, frameon=False, fontsize=8)

    fig.tight_layout()
    Path(png_path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(png_path, dpi=200, bbox_inches="tight")
    fig.savefig(pdf_path, bbox_inches="tight")
    plt.close(fig)
    logger.info("build_figure: wrote %s and %s (%d samples)", png_path, pdf_path, n)


def write_outputs(by_sample_df: pd.DataFrame, by_tumor_df: pd.DataFrame, cfg: CnvMethodComparisonConfig) -> None:
    for path in (cfg.by_sample_tsv, cfg.by_tumor_tsv):
        path.parent.mkdir(parents=True, exist_ok=True)
    by_sample_df.to_csv(cfg.by_sample_tsv, sep="\t", index=False)
    by_tumor_df.to_csv(cfg.by_tumor_tsv, sep="\t", index=False)
    logger.info(
        "write_outputs: wrote %s (%d samples) and %s (%d tumors)",
        cfg.by_sample_tsv,
        len(by_sample_df),
        cfg.by_tumor_tsv,
        len(by_tumor_df),
    )


def run_cnv_method_comparison(
    config_path: str | Path = "config/config.yaml",
    expected_overall_agreement_pct: float | None = 54.08,
) -> dict[str, pd.DataFrame]:
    """Run the full load / strict-match / by-sample / by-tumor / verify /
    write pipeline. Pass ``expected_overall_agreement_pct=None`` to skip
    the reproduction check (used only by tests with synthetic data that
    are not expected to reproduce the real-data figure)."""
    config = _load_config(config_path)
    cfg = CnvMethodComparisonConfig.from_config(config)

    infercnv_df = load_infercnv_labels(cfg.infercnv_labels_tsv)
    copykat_df = load_copykat_labels(cfg.copykat_labels_tsv)
    summary_df = load_infercnv_sample_summary(cfg.infercnv_sample_summary_tsv)

    n_samples = infercnv_df["sample_id"].nunique()
    if n_samples != cfg.expected_n_tumor_samples:
        raise ValueError(f"expected exactly {cfg.expected_n_tumor_samples} samples, found {n_samples}")

    matched_df = match_cells(infercnv_df, copykat_df)
    by_sample_df = build_by_sample_table(matched_df, summary_df)
    by_tumor_df = build_by_tumor_table(by_sample_df)

    if expected_overall_agreement_pct is not None:
        verify_overall_agreement(by_sample_df, expected_overall_agreement_pct)

    write_outputs(by_sample_df, by_tumor_df, cfg)
    build_figure(by_sample_df, cfg.figure_png, cfg.figure_pdf)

    return {"by_sample": by_sample_df, "by_tumor": by_tumor_df}


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run_cnv_method_comparison()
