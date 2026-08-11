"""Matched, same-cell InferCNV-vs-CopyKAT heatmap comparison for the 6
priority GSE245601 tumors (docs/CNV_METHOD_AUDIT.md Point 1 follow-up,
visualization step).

This module reruns neither method and changes no malignancy label. It
reads:

- the compact per-sample gene-by-cell matrices written by
  ``scripts/analysis/gse245601_10_extract_cnv_comparison_matrices.R``
  (itself read-only against the frozen ``run.final.infercnv_obj`` and
  CopyKAT ``*_raw_results_gene_by_cell.txt`` outputs),
- the frozen per-cell label tables (``gse245601_malignant_cell_labels.tsv``,
  ``gse245601_copykat_sensitivity_labels.tsv``), and
- the already-verified ``cnv_method_comparison_by_sample.tsv`` (Point 1
  output) for the summary-table fractions, rather than recomputing them a
  second way --

and produces, per selected sample, a side-by-side heatmap figure (same
genes as columns, same cells as rows, in the same order in both panels)
with three annotation strips (InferCNV call, CopyKAT call, and a derived
agreement status), plus one tumor-grouped contact sheet and a short
machine-readable summary table.

InferCNV's ``@expr.data`` (reference-centered at 1.0) and CopyKAT's
gene-level ``raw_results`` (log-ratio-like, centered at 0.0) are DIFFERENT
NATIVE SCALES from two different algorithms -- they are plotted with
separate colormaps and separate colorbars, never a shared scale, and
nothing in this module implies the two heatmaps are numerically the same
kind of object. Chromosome-scale pattern concordance and classification
concordance are visually separable (the two heatmap panels vs. the three
annotation strips) but this module draws no automated conclusion about
which pattern explains a given disagreement -- that is left to visual
inspection.
"""

from __future__ import annotations

import logging
import subprocess
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.patches as mpatches  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import yaml  # noqa: E402
from matplotlib.colors import ListedColormap, TwoSlopeNorm  # noqa: E402

from src.gse245601_cnv_method_comparison import (  # noqa: E402
    LABEL_MALIGNANT,
    LABEL_NONMALIGNANT,
    LABEL_NOT_DEFINED,
    load_copykat_labels,
    load_infercnv_labels,
)

logger = logging.getLogger(__name__)

AGREE_TUMOR_LIKE = "agree_tumor_like"
AGREE_NORMAL_LIKE = "agree_normal_like"
DISAGREE = "disagree"
COPYKAT_NOT_DEFINED_STATUS = "copykat_not_defined"

AGREEMENT_STATUSES = (AGREE_TUMOR_LIKE, AGREE_NORMAL_LIKE, DISAGREE, COPYKAT_NOT_DEFINED_STATUS)

# Fixed, semantically-consistent categorical colors: "malignant-like" is
# always the same red regardless of which method or strip, so visual
# agreement between strips is immediately readable as matching color.
COLOR_MALIGNANT_LIKE = "#B2182B"
COLOR_NORMAL_LIKE = "#2166AC"
COLOR_NOT_DEFINED = "#969696"
COLOR_DISAGREE = "#F1A340"

METHOD_LABEL_COLORS: dict[str, str] = {
    LABEL_MALIGNANT: COLOR_MALIGNANT_LIKE,
    LABEL_NONMALIGNANT: COLOR_NORMAL_LIKE,
    LABEL_NOT_DEFINED: COLOR_NOT_DEFINED,
}
AGREEMENT_COLORS: dict[str, str] = {
    AGREE_TUMOR_LIKE: COLOR_MALIGNANT_LIKE,
    AGREE_NORMAL_LIKE: COLOR_NORMAL_LIKE,
    DISAGREE: COLOR_DISAGREE,
    COPYKAT_NOT_DEFINED_STATUS: COLOR_NOT_DEFINED,
}

INFERCNV_CENTER = 1.0
COPYKAT_CENTER = 0.0


def _load_config(config_path: str | Path) -> dict:
    with open(config_path) as f:
        return yaml.safe_load(f)


@dataclass(frozen=True)
class CnvHeatmapComparisonConfig:
    """Resolved, config-driven paths. No hardcoded paths."""

    selected_samples: tuple[str, ...]
    infercnv_labels_tsv: Path
    copykat_labels_tsv: Path
    method_comparison_by_sample_tsv: Path
    working_dir: Path
    extraction_r_script: Path
    micromamba_binary: str
    micromamba_env: str
    figures_dir: Path
    tables_dir: Path
    summary_tsv: Path
    gene_exclusion_tsv: Path
    contact_sheet_png: Path

    @classmethod
    def from_config(cls, config: dict) -> "CnvHeatmapComparisonConfig":
        cfg = config["gse245601_cnv_heatmap_comparison"]
        inputs = cfg["inputs"]
        out = cfg["output"]
        return cls(
            selected_samples=tuple(cfg["selected_samples"]),
            infercnv_labels_tsv=Path(inputs["infercnv_labels_tsv"]),
            copykat_labels_tsv=Path(inputs["copykat_labels_tsv"]),
            method_comparison_by_sample_tsv=Path(inputs["method_comparison_by_sample_tsv"]),
            working_dir=Path(cfg["working_dir"]),
            extraction_r_script=Path(cfg["extraction_r_script"]),
            micromamba_binary=cfg["micromamba_binary"],
            micromamba_env=cfg["micromamba_env"],
            figures_dir=Path(out["figures_dir"]),
            tables_dir=Path(out["tables_dir"]),
            summary_tsv=Path(out["summary_tsv"]),
            gene_exclusion_tsv=Path(out["gene_exclusion_tsv"]),
            contact_sheet_png=Path(out["contact_sheet_png"]),
        )


def run_r_extraction(cfg: CnvHeatmapComparisonConfig) -> None:
    """Invoke the read-only R extraction script (requires the sc245601
    micromamba env for the infercnv/data.table/yaml packages). Neither
    InferCNV nor CopyKAT is rerun -- this only reshapes their already-saved
    outputs."""
    cmd = [cfg.micromamba_binary, "run", "-n", cfg.micromamba_env, "Rscript", str(cfg.extraction_r_script)]
    logger.info("run_r_extraction: %s", " ".join(cmd))
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"R extraction script failed (exit {result.returncode}):\n{result.stdout}\n{result.stderr}")
    logger.info("run_r_extraction: completed successfully")


def compute_agreement_status(infercnv_label: str, copykat_label: str) -> str:
    """Three-way (four-category) per-cell agreement status. CopyKAT
    ``not_defined`` is its own category -- it is neither agreement nor
    disagreement, since InferCNV made a call CopyKAT could not confirm or
    contradict. Callers are expected to have already validated label
    vocabulary (``load_infercnv_labels``/``load_copykat_labels`` already
    reject unexpected values) -- this function still checks its own inputs
    rather than silently classifying an unrecognized label as
    ``disagree``."""
    if infercnv_label not in (LABEL_MALIGNANT, LABEL_NONMALIGNANT):
        raise ValueError(f"unexpected infercnv_label: {infercnv_label!r}")
    if copykat_label not in (LABEL_MALIGNANT, LABEL_NONMALIGNANT, LABEL_NOT_DEFINED):
        raise ValueError(f"unexpected copykat_label: {copykat_label!r}")
    if copykat_label == LABEL_NOT_DEFINED:
        return COPYKAT_NOT_DEFINED_STATUS
    if infercnv_label == LABEL_MALIGNANT and copykat_label == LABEL_MALIGNANT:
        return AGREE_TUMOR_LIKE
    if infercnv_label == LABEL_NONMALIGNANT and copykat_label == LABEL_NONMALIGNANT:
        return AGREE_NORMAL_LIKE
    return DISAGREE


def build_annotation_table(sample_id: str, infercnv_labels_df: pd.DataFrame, copykat_labels_df: pd.DataFrame) -> pd.DataFrame:
    """Per-cell annotation table for one sample: InferCNV label, CopyKAT
    label, derived agreement status, and cnv_score (used only for row
    ordering). Cells are matched by cell_id via an inner join validated to
    be one-to-one and to lose no cell on either side."""
    ic = infercnv_labels_df.loc[infercnv_labels_df["sample_id"] == sample_id, ["cell_id", "primary_malignancy_label", "cnv_score"]]
    ck = copykat_labels_df.loc[copykat_labels_df["sample_id"] == sample_id, ["cell_id", "sensitivity_malignancy_label"]]

    ic_only = set(ic["cell_id"]) - set(ck["cell_id"])
    ck_only = set(ck["cell_id"]) - set(ic["cell_id"])
    if ic_only or ck_only:
        raise ValueError(
            f"{sample_id}: cell_id sets do not match between InferCNV and CopyKAT label tables: "
            f"{len(ic_only)} cell(s) only in InferCNV, {len(ck_only)} cell(s) only in CopyKAT"
        )

    merged = ic.merge(ck, on="cell_id", how="inner", validate="one_to_one")
    merged = merged.rename(
        columns={"primary_malignancy_label": "infercnv_label", "sensitivity_malignancy_label": "copykat_label"}
    )
    merged["agreement_status"] = [
        compute_agreement_status(ic_lbl, ck_lbl) for ic_lbl, ck_lbl in zip(merged["infercnv_label"], merged["copykat_label"])
    ]
    logger.info("build_annotation_table: %s -- %d matched cells", sample_id, len(merged))
    return merged


def build_cell_order(annotation_df: pd.DataFrame) -> list[str]:
    """Row order shared by both heatmap panels and all three annotation
    strips: InferCNV malignant cells first (highest cnv_score first), then
    InferCNV non-malignant cells (highest cnv_score first). Uses only
    already-frozen InferCNV fields -- no new classification."""
    rank = annotation_df["infercnv_label"].map({LABEL_MALIGNANT: 0, LABEL_NONMALIGNANT: 1})
    if rank.isna().any():
        raise ValueError("build_cell_order: unexpected infercnv_label value(s) present")
    ordered = annotation_df.assign(_rank=rank).sort_values(["_rank", "cnv_score"], ascending=[True, False])
    return ordered["cell_id"].tolist()


def load_compact_matrices(sample_id: str, cfg: CnvHeatmapComparisonConfig) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Load the R-extraction script's per-sample outputs: InferCNV matrix,
    CopyKAT matrix (both genes x epithelial cells, gene as index), and the
    shared gene order (chromosome-sorted). Validates the gene axis is
    unique and identical across all three sources -- the R script is
    expected to guarantee this, but it is checked here rather than
    silently trusted across the R/Python boundary."""
    sample_dir = cfg.working_dir / sample_id
    ic = pd.read_csv(sample_dir / "infercnv_matrix.tsv.gz", sep="\t", index_col="gene")
    ck = pd.read_csv(sample_dir / "copykat_matrix.tsv.gz", sep="\t", index_col="gene")
    gene_order = pd.read_csv(sample_dir / "gene_order.tsv", sep="\t")

    if ic.index.duplicated().any():
        raise ValueError(f"{sample_id}: infercnv_matrix.tsv.gz has duplicate gene rows")
    if ck.index.duplicated().any():
        raise ValueError(f"{sample_id}: copykat_matrix.tsv.gz has duplicate gene rows")
    if gene_order["gene"].duplicated().any():
        raise ValueError(f"{sample_id}: gene_order.tsv has duplicate gene rows")
    if set(ic.index) != set(gene_order["gene"]) or set(ck.index) != set(gene_order["gene"]):
        raise ValueError(f"{sample_id}: infercnv/copykat/gene_order gene sets do not match exactly")

    logger.info(
        "load_compact_matrices: %s -- infercnv %s, copykat %s, %d genes in shared order",
        sample_id,
        ic.shape,
        ck.shape,
        len(gene_order),
    )
    return ic, ck, gene_order


def verify_matched_columns(infercnv_matrix: pd.DataFrame, copykat_matrix: pd.DataFrame, expected_cell_ids: list[str]) -> None:
    """Cell-matching correctness check: both matrices must have exactly
    ``expected_cell_ids`` as their column set (order-independent at this
    check; the caller reindexes to a shared order afterwards), and neither
    matrix nor the expected list may contain a duplicate cell_id -- set
    equality alone would not catch a duplicated column. Raises rather than
    silently reindexing over a mismatch."""
    if pd.Index(expected_cell_ids).duplicated().any():
        raise ValueError("expected_cell_ids contains duplicate cell_id value(s)")
    if infercnv_matrix.columns.duplicated().any():
        raise ValueError("InferCNV matrix has duplicate cell_id column(s)")
    if copykat_matrix.columns.duplicated().any():
        raise ValueError("CopyKAT matrix has duplicate cell_id column(s)")

    expected = set(expected_cell_ids)
    ic_cols = set(infercnv_matrix.columns)
    ck_cols = set(copykat_matrix.columns)
    if ic_cols != expected:
        raise ValueError(
            f"InferCNV matrix columns do not match the expected cell set: "
            f"{len(ic_cols - expected)} extra, {len(expected - ic_cols)} missing"
        )
    if ck_cols != expected:
        raise ValueError(
            f"CopyKAT matrix columns do not match the expected cell set: "
            f"{len(ck_cols - expected)} extra, {len(expected - ck_cols)} missing"
        )


def _chromosome_boundaries(gene_order: pd.DataFrame) -> list[tuple[int, int, int]]:
    """Returns (chr, start_col_index, end_col_index) for each contiguous
    chromosome block in ``gene_order``'s row order."""
    chrs = gene_order["chr"].to_numpy()
    boundaries = []
    start = 0
    for i in range(1, len(chrs) + 1):
        if i == len(chrs) or chrs[i] != chrs[start]:
            boundaries.append((int(chrs[start]), start, i))
            start = i
    return boundaries


def _draw_annotation_strip(ax, statuses: list[str], color_map: dict[str, str], title: str) -> None:
    codes = {status: i for i, status in enumerate(color_map)}
    arr = np.array([[codes[s]] for s in statuses])
    cmap = ListedColormap([color_map[s] for s in color_map])
    ax.imshow(arr, aspect="auto", cmap=cmap, vmin=0, vmax=len(color_map) - 1, interpolation="none")
    ax.set_xticks([0])
    ax.set_xticklabels([title], fontsize=6, rotation=90)
    ax.tick_params(axis="x", length=0)
    ax.set_yticks([])


def _draw_heatmap_panel(ax, matrix: pd.DataFrame, cell_order: list[str], gene_order: pd.DataFrame, center: float, title: str) -> None:
    arr = matrix.reindex(index=gene_order["gene"], columns=cell_order).to_numpy()
    if np.isnan(arr).any():
        raise ValueError(f"{title}: NaN values after reindexing to the shared gene/cell order -- matrix/order mismatch")
    arr = arr.T  # -> (n_cells, n_genes) so rows align with the annotation strips
    clip = np.percentile(np.abs(arr - center), 99)
    clip = max(clip, 1e-6)
    norm = TwoSlopeNorm(vcenter=center, vmin=center - clip, vmax=center + clip)
    im = ax.imshow(arr, aspect="auto", cmap="RdBu_r", norm=norm, interpolation="none")
    ax.set_yticks([])
    for chr_num, col_start, col_end in _chromosome_boundaries(gene_order):
        if col_start > 0:
            ax.axvline(col_start - 0.5, color="black", linewidth=0.3, alpha=0.5)
    tick_pos = [(s + e) / 2 - 0.5 for _, s, e in _chromosome_boundaries(gene_order)]
    tick_lbl = [str(c) for c, _, _ in _chromosome_boundaries(gene_order)]
    ax.set_xticks(tick_pos)
    ax.set_xticklabels(tick_lbl, fontsize=5, rotation=90)
    ax.set_xlabel("chromosome", fontsize=7, labelpad=8)
    ax.set_title(title, fontsize=8)
    return im


def build_sample_figure(
    sample_id: str,
    infercnv_matrix: pd.DataFrame,
    copykat_matrix: pd.DataFrame,
    gene_order: pd.DataFrame,
    annotation_df: pd.DataFrame,
    out_path: str | Path,
) -> None:
    """Side-by-side InferCNV/CopyKAT heatmap for one sample: shared gene
    columns (chromosome-sorted), shared cell rows (InferCNV-label-then-
    score order), three left-side annotation strips (InferCNV call,
    CopyKAT call, agreement status). InferCNV and CopyKAT are plotted on
    separate, independently-scaled colormaps (1.0-centered vs 0.0-centered
    native scales) -- never a shared colorbar."""
    cell_order = build_cell_order(annotation_df)
    verify_matched_columns(infercnv_matrix, copykat_matrix, cell_order)
    annotation_df = annotation_df.set_index("cell_id").loc[cell_order]

    n_cells = len(cell_order)
    n_genes = len(gene_order)
    height = min(max(3.5, 0.0022 * n_cells + 2.0), 14.0)
    top_frac = min(0.92, max(0.72, 1 - 0.9 / height))
    bottom_frac = max(0.08, min(0.25, 0.8 / height))

    fig = plt.figure(figsize=(13, height))
    gs = fig.add_gridspec(1, 5, width_ratios=[0.05, 0.05, 0.05, 1, 1], wspace=0.08, top=top_frac, bottom=bottom_frac)

    ax_agree = fig.add_subplot(gs[0, 0])
    ax_ic_lbl = fig.add_subplot(gs[0, 1])
    ax_ck_lbl = fig.add_subplot(gs[0, 2])
    ax_ic = fig.add_subplot(gs[0, 3])
    ax_ck = fig.add_subplot(gs[0, 4], sharey=ax_ic)

    _draw_annotation_strip(ax_agree, annotation_df["agreement_status"].tolist(), AGREEMENT_COLORS, "agreement")
    _draw_annotation_strip(ax_ic_lbl, annotation_df["infercnv_label"].tolist(), METHOD_LABEL_COLORS, "InferCNV")
    _draw_annotation_strip(ax_ck_lbl, annotation_df["copykat_label"].tolist(), METHOD_LABEL_COLORS, "CopyKAT")

    im_ic = _draw_heatmap_panel(ax_ic, infercnv_matrix, cell_order, gene_order, INFERCNV_CENTER, "InferCNV relative expression (centered at 1.0)")
    im_ck = _draw_heatmap_panel(ax_ck, copykat_matrix, cell_order, gene_order, COPYKAT_CENTER, "CopyKAT relative signal (centered at 0.0)")
    plt.setp(ax_ck.get_yticklabels(), visible=False)

    fig.colorbar(im_ic, ax=ax_ic, fraction=0.02, pad=0.02)
    fig.colorbar(im_ck, ax=ax_ck, fraction=0.02, pad=0.02)

    counts = annotation_df["agreement_status"].value_counts()
    legend_handles = [
        mpatches.Patch(color=color, label=status.replace("_", " ")) for status, color in AGREEMENT_COLORS.items()
    ]
    fig.legend(
        handles=legend_handles, loc="lower center", ncol=4, frameon=False, fontsize=7, bbox_to_anchor=(0.5, top_frac + 0.02)
    )

    fig.text(
        0.5,
        0.99,
        f"{sample_id} -- n={n_cells} matched epithelial cells, {n_genes} shared genes "
        f"(agree tumor-like={counts.get(AGREE_TUMOR_LIKE, 0)}, agree normal-like={counts.get(AGREE_NORMAL_LIKE, 0)}, "
        f"disagree={counts.get(DISAGREE, 0)}, CopyKAT not.defined={counts.get(COPYKAT_NOT_DEFINED_STATUS, 0)})",
        fontsize=9,
        ha="center",
        va="top",
    )

    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info("build_sample_figure: wrote %s", out_path)


def build_contact_sheet(selected_samples: tuple[str, ...], figures_dir: Path, out_path: str | Path) -> None:
    """Tumor-grouped contact sheet: one row per tumor, Control and
    Tamoxifen columns adjacent, built by tiling the already-rendered
    per-sample PNGs (guarantees the contact sheet always matches the
    individual figures exactly)."""
    if len(selected_samples) % 2 != 0:
        raise ValueError("build_contact_sheet: expected an even number of samples (Control/Tamoxifen pairs)")
    pairs = [selected_samples[i : i + 2] for i in range(0, len(selected_samples), 2)]
    seen_patients: set[str] = set()
    for control_id, tam_id in pairs:
        control_patient = control_id.rsplit("_", 1)[0]
        tam_patient, tam_condition = tam_id.rsplit("_", 1)
        if control_id.rsplit("_", 1)[1] != "Control" or tam_condition != "Tamoxifen" or control_patient != tam_patient:
            raise ValueError(f"build_contact_sheet: expected adjacent (Control, Tamoxifen) pair for one tumor, got ({control_id}, {tam_id})")
        if control_patient in seen_patients:
            raise ValueError(f"build_contact_sheet: tumor {control_patient} appears in more than one pair")
        seen_patients.add(control_patient)

    n_rows = len(pairs)
    fig, axes = plt.subplots(n_rows, 2, figsize=(20, 4.2 * n_rows), squeeze=False)
    for row, (control_id, tam_id) in enumerate(pairs):
        for col, sample_id in enumerate((control_id, tam_id)):
            img = plt.imread(figures_dir / f"{sample_id}_heatmap_comparison.png")
            ax = axes[row, col]
            ax.imshow(img, aspect="auto")  # fill the grid cell uniformly; per-sample figures keep their native aspect ratio on disk
            ax.set_title(sample_id, fontsize=9)
            ax.axis("off")

    fig.tight_layout()
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info("build_contact_sheet: wrote %s (%d samples)", out_path, len(selected_samples))


def build_summary_table(cfg: CnvHeatmapComparisonConfig) -> pd.DataFrame:
    """Short machine-readable per-sample summary, filtered from the
    already-verified ``cnv_method_comparison_by_sample.tsv`` (Point 1
    output) rather than recomputed a second way, plus the shared-gene
    count from the extraction step's own exclusion report."""
    by_sample = pd.read_csv(cfg.method_comparison_by_sample_tsv, sep="\t")
    selected = by_sample.loc[by_sample["sample_id"].isin(cfg.selected_samples)].copy()
    missing = set(cfg.selected_samples) - set(selected["sample_id"])
    if missing:
        raise ValueError(f"{cfg.method_comparison_by_sample_tsv}: missing selected sample(s): {missing}")

    exclusion = pd.read_csv(cfg.gene_exclusion_tsv, sep="\t")[["sample_id", "n_shared_genes_final"]]
    selected = selected.merge(exclusion, on="sample_id", how="left", validate="one_to_one")
    if selected["n_shared_genes_final"].isna().any():
        raise ValueError("build_summary_table: gene exclusion report is missing a selected sample")

    out = pd.DataFrame(
        {
            "sample": selected["sample_id"],
            "n_matched_cells": selected["total_epithelial_cells"].astype(int),
            "n_shared_genes": selected["n_shared_genes_final"].astype(int),
            "infercnv_malignant_fraction": selected["infercnv_malignant_pct"] / 100.0,
            "copykat_aneuploid_fraction": selected["copykat_aneuploid_pct"] / 100.0,
            "agreement_fraction": selected["agreement_pct"] / 100.0,
        }
    )
    out["sample"] = pd.Categorical(out["sample"], categories=list(cfg.selected_samples), ordered=True)
    out = out.sort_values("sample").reset_index(drop=True)
    out["sample"] = out["sample"].astype(str)

    Path(cfg.summary_tsv).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(cfg.summary_tsv, sep="\t", index=False)
    logger.info("build_summary_table: wrote %s (%d samples)", cfg.summary_tsv, len(out))
    return out


def run_cnv_heatmap_comparison(
    config_path: str | Path = "config/config.yaml",
    run_r_extraction_step: bool = True,
    samples: tuple[str, ...] | None = None,
) -> dict[str, object]:
    """Full pipeline: (optionally) run the R extraction step, then build
    every selected sample's heatmap figure, the contact sheet, and the
    summary table. ``samples`` restricts which selected samples get
    figures built in this call (used by tests to run only the smallest
    sample pair) -- it must be a subset of the configured
    ``selected_samples``; the summary table and contact sheet still cover
    the full configured set unless ``samples`` is also used to subset
    those explicitly by the caller."""
    config = _load_config(config_path)
    cfg = CnvHeatmapComparisonConfig.from_config(config)

    if samples is not None:
        unknown = set(samples) - set(cfg.selected_samples)
        if unknown:
            raise ValueError(f"samples not in configured selected_samples: {unknown}")
    samples_to_build = samples if samples is not None else cfg.selected_samples

    if run_r_extraction_step:
        run_r_extraction(cfg)

    infercnv_labels_df = load_infercnv_labels(cfg.infercnv_labels_tsv)
    copykat_labels_df = load_copykat_labels(cfg.copykat_labels_tsv)

    figure_paths = {}
    for sample_id in samples_to_build:
        annotation_df = build_annotation_table(sample_id, infercnv_labels_df, copykat_labels_df)
        ic_mat, ck_mat, gene_order = load_compact_matrices(sample_id, cfg)
        out_path = cfg.figures_dir / f"{sample_id}_heatmap_comparison.png"
        build_sample_figure(sample_id, ic_mat, ck_mat, gene_order, annotation_df, out_path)
        figure_paths[sample_id] = out_path

    result: dict[str, object] = {"figure_paths": figure_paths}

    if samples is None:
        build_contact_sheet(cfg.selected_samples, cfg.figures_dir, cfg.contact_sheet_png)
        result["summary"] = build_summary_table(cfg)

    return result


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run_cnv_heatmap_comparison(run_r_extraction_step=False)
