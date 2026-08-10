"""NEBULA final poster figures -- one stable, independently reproducible
plotting module for the five visually-approved poster figures.

Approved sources (frozen appearance, carried over verbatim from the
experimental iteration history that led here):

- Figure 1 (CRISPR landscape): approved v6 design.
- Figure 2 (PCA):               approved v6 design.
- Figure 3 (volcano):           approved v6 design.
- Figure 4 (expression heatmap): approved v8 design (the v7 five-region
  layout, with v8's shared row/column geometry fix so block separators
  and group boundaries are numerically identical across every region).
- Figure 5 (USP34 panel):        approved v8 design.

This module makes no scientific change of any kind: no CRISPR/limma/PCA
recomputation, no candidate reclassification, no new statistic, no changed
threshold, no changed gene/sample ordering. It reads exclusively the
frozen, committed plot-input tables already written by
``src.candidate_evidence_summary`` and ``src.gse118713_qc``
(``results/tables/nebula_plot_inputs/*.tsv``), the checksum-pinned
filtered GSE118713 TPM matrix, and the frozen sample metadata table --
the same frozen inputs every experimental nebula_plots_v2..v8 module used.

Independence: this module does NOT import nebula_plots_v2 through
nebula_plots_v8 (those experimental modules are being retired). All
loader, QC, geometry, and plotting logic needed from them has been copied
here unchanged and is re-verified byte-for-byte identical to their output
in ``tests/test_nebula_plots_final.py``.
"""

from __future__ import annotations

import hashlib
import logging
import math
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from matplotlib.text import Text
from matplotlib.transforms import Bbox

from src.candidate_evidence_summary import (
    BENCHMARK_GENE_SYMBOL,
    BENCHMARK_LABEL,
    DIRECTION_SENSITISING,
    DIRECTION_TOLERANCE,
    EVIDENCE_CLASS_NO_SIGNIFICANT_RNA,
    EVIDENCE_CLASS_PRIMARY,
    EVIDENCE_CLASS_RNA_UNAVAILABLE,
    EVIDENCE_CLASS_SECONDARY,
)

logger = logging.getLogger(__name__)

_DARK_TEXT = "#0B0B12"
_MUTED_TEXT = "#4B5563"
_UNAVAILABLE_GREY = "#D1D5DB"

_SAMPLE_GROUP_ORDER: tuple[str, ...] = ("MCF7", "TAMR", "FASR")

# Four clearly distinct evidence-class accent colors -- approved v6/v8
# palette, verified pairwise distinct and distinct from the heatmap's own
# blue/orange scale.
_EVIDENCE_CLASS_ACCENT = {
    EVIDENCE_CLASS_PRIMARY: "#6E44FF",  # deep purple
    EVIDENCE_CLASS_SECONDARY: "#0F9B8E",  # teal
    EVIDENCE_CLASS_NO_SIGNIFICANT_RNA: "#A16207",  # amber
    EVIDENCE_CLASS_RNA_UNAVAILABLE: "#94A3B8",  # slate grey
}

# Approved v8 plain-language, wrapped labels.
_EVIDENCE_CLASS_PLAIN_LABEL = {
    EVIDENCE_CLASS_PRIMARY: "TAMR vs MCF7\nsupport",
    EVIDENCE_CLASS_SECONDARY: "TAMR vs FASR\ncontext only",
    EVIDENCE_CLASS_NO_SIGNIFICANT_RNA: "No significant\nbulk-RNA support",
    EVIDENCE_CLASS_RNA_UNAVAILABLE: "RNA\nunavailable",
}

FIG1_TITLE = "CRISPR screen identifies 28 candidate 4-OHT response genes"
FIG2_TITLE = "PCA separates MCF7, TAMR and FASR expression profiles"

FIG4_TITLE = "Expression patterns of sensitising CRISPR candidates"
FIG4_SUBTITLE = "GSE118713 · parental and endocrine-resistant MCF7 states"
FIG4_FOOTER = "Relative expression shown as row z-score of log2(TPM+1); visualization only."
FIG4_USP34_ANNOTATION = "↑ TAMR vs MCF7\nFDR {fdr}"
FIG4_USP17L29_ANNOTATION = "Expression unavailable\nfiltered before DE"

FIG5_TITLE = "USP34 links sensitising CRISPR evidence with resistance-associated expression"
FIG5_TAMR_VS_MCF7_LABEL = "TAMR vs MCF7\nFDR {fdr}"
FIG5_TAMR_VS_FASR_LABEL = "TAMR vs FASR\nn.s. · FDR {fdr}"
FIG5_CRISPR_BOX = "CRISPR knockout\neffect {effect} · FDR {fdr}"
FIG5_FOOTER = "Expression is elevated vs parental MCF7 but not different from FASR; experimental validation required."


def _fmt_effect(x: float) -> str:
    return f"{x:+.2f}"


def _fmt_fdr(x: float) -> str:
    return f"{x:.3g}"


def _hex_to_rgb(hex_color: str) -> tuple[int, int, int]:
    hex_color = hex_color.lstrip("#")
    return tuple(int(hex_color[i : i + 2], 16) for i in (0, 2, 4))


def min_pairwise_color_distance(colors: list[str]) -> float:
    """Smallest Euclidean RGB distance between any pair of the given hex
    colors -- a quick, dependency-free sanity check that a small
    categorical palette is not accidentally near-duplicated."""
    rgbs = [_hex_to_rgb(c) for c in colors]
    distances = []
    for i in range(len(rgbs)):
        for j in range(i + 1, len(rgbs)):
            d = math.sqrt(sum((a - b) ** 2 for a, b in zip(rgbs[i], rgbs[j])))
            distances.append(d)
    return min(distances) if distances else float("inf")


def _load_config(config_path: str | Path) -> dict:
    with open(config_path) as f:
        return yaml.safe_load(f)


def _sha256_file(path: str | Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


@dataclass(frozen=True)
class NebulaPlotsFinalConfig:
    """Resolved, config-driven paths and style. No hardcoded paths.

    All frozen-data-source paths are read from the committed,
    version-independent ``nebula_plots`` config section (written by v1 and
    never modified since) plus other pre-existing committed sections
    (``crispr_gse118713_integration``, ``gse118713_phase2b``,
    ``gse118713``) -- no dependency on any nebula_plots_v2..v8 config
    section.
    """

    fig1_input_tsv: Path
    fig1_paics_inset_tsv: Path
    fig2_input_tsv: Path
    fig3_input_tsv: Path
    fig5_summary_input_tsv: Path
    fig5_expression_input_tsv: Path
    rna_significance_fdr: float
    expected_n_hits: int
    expected_n_tolerance: int

    fig4_evidence_input_tsv: Path
    filtered_gene_tpm_tsv_gz: Path
    frozen_filtered_gene_tpm_sha256: str
    sample_metadata_tsv: Path

    expected_n_sensitising: int
    expected_n_samples: int
    expected_primary_gene: str
    panel_background: str
    palette: dict[str, str]
    group_colors: dict[str, str]
    direction_colors: dict[str, str]

    output_dir: Path
    contact_sheet_png: Path
    manifest_tsv: Path
    fig1_png: Path
    fig1_pdf: Path
    fig2_png: Path
    fig2_pdf: Path
    fig3_png: Path
    fig3_pdf: Path
    fig4_png: Path
    fig4_pdf: Path
    fig5_png: Path
    fig5_pdf: Path

    @classmethod
    def from_config(cls, config: dict) -> "NebulaPlotsFinalConfig":
        np_cfg = config["nebula_plots"]
        figs_v1 = np_cfg["figures"]
        style = np_cfg["style"]
        npqc = np_cfg["qc"]
        integ = config["crispr_gse118713_integration"]
        filtering = config["gse118713_phase2b"]["filtering"]
        gse = config["gse118713"]
        final_cfg = config["nebula_plots_final"]
        figs_final = final_cfg["figures"]

        return cls(
            fig1_input_tsv=Path(figs_v1["crispr_landscape"]["plot_input_tsv"]),
            fig1_paics_inset_tsv=Path(figs_v1["crispr_landscape"]["paics_inset_tsv"]),
            fig2_input_tsv=Path(figs_v1["pca"]["plot_input_tsv"]),
            fig3_input_tsv=Path(figs_v1["volcano"]["plot_input_tsv"]),
            fig5_summary_input_tsv=Path(figs_v1["usp34_panel"]["summary_plot_input_tsv"]),
            fig5_expression_input_tsv=Path(figs_v1["usp34_panel"]["expression_plot_input_tsv"]),
            rna_significance_fdr=float(np_cfg["rna_significance_fdr"]),
            expected_n_hits=int(integ["expected_n_hits"]),
            expected_n_tolerance=int(npqc["expected_n_tolerance"]),
            fig4_evidence_input_tsv=Path(figs_v1["evidence_matrix"]["plot_input_tsv"]),
            filtered_gene_tpm_tsv_gz=Path(filtering["filtered_gene_tpm_tsv"]),
            frozen_filtered_gene_tpm_sha256=str(filtering["frozen_filtered_gene_tpm_sha256"]),
            sample_metadata_tsv=Path(gse["output"]["sample_metadata_tsv"]),
            expected_n_sensitising=int(npqc["expected_n_sensitising"]),
            expected_n_samples=int(npqc["expected_n_samples"]),
            expected_primary_gene=str(npqc["expected_primary_gene"]),
            panel_background=str(style["panel_background"]),
            palette=dict(style["palette"]),
            group_colors=dict(style["group_colors"]),
            direction_colors=dict(style["direction_colors"]),
            output_dir=Path(final_cfg["output_dir"]),
            contact_sheet_png=Path(final_cfg["contact_sheet_png"]),
            manifest_tsv=Path(final_cfg["manifest_tsv"]),
            fig1_png=Path(figs_final["crispr_landscape"]["png"]),
            fig1_pdf=Path(figs_final["crispr_landscape"]["pdf"]),
            fig2_png=Path(figs_final["pca"]["png"]),
            fig2_pdf=Path(figs_final["pca"]["pdf"]),
            fig3_png=Path(figs_final["volcano"]["png"]),
            fig3_pdf=Path(figs_final["volcano"]["pdf"]),
            fig4_png=Path(figs_final["evidence_matrix"]["png"]),
            fig4_pdf=Path(figs_final["evidence_matrix"]["pdf"]),
            fig5_png=Path(figs_final["usp34_panel"]["png"]),
            fig5_pdf=Path(figs_final["usp34_panel"]["pdf"]),
        )


# --------------------------------------------------------------------------
# Loading frozen plot inputs (read-only; no recomputation) -- copied
# unchanged from nebula_plots_v2 / nebula_plots_v3.
# --------------------------------------------------------------------------


def load_fig1_input(cfg: NebulaPlotsFinalConfig) -> pd.DataFrame:
    df = pd.read_csv(cfg.fig1_input_tsv, sep="\t")
    if len(df) != cfg.expected_n_hits:
        raise ValueError(f"expected {cfg.expected_n_hits} Gate-1 hits in fig1 input, found {len(df)}")
    if BENCHMARK_GENE_SYMBOL in set(df["gene_symbol"]):
        raise ValueError(f"{BENCHMARK_GENE_SYMBOL} must never appear among the {cfg.expected_n_hits} Gate-1 hits")
    logger.info("load_fig1_input: read %d rows", len(df))
    return df


def load_fig1_paics_inset(cfg: NebulaPlotsFinalConfig) -> pd.DataFrame:
    df = pd.read_csv(cfg.fig1_paics_inset_tsv, sep="\t")
    if len(df) != 1 or df.iloc[0]["gene_symbol"] != BENCHMARK_GENE_SYMBOL:
        raise ValueError("expected exactly one PAICS benchmark row")
    if df.iloc[0]["benchmark_label"] != BENCHMARK_LABEL:
        raise ValueError(f"PAICS benchmark row missing expected label {BENCHMARK_LABEL!r}")
    logger.info("load_fig1_paics_inset: confirmed PAICS is labelled %r", BENCHMARK_LABEL)
    return df


def load_fig2_input(cfg: NebulaPlotsFinalConfig) -> pd.DataFrame:
    df = pd.read_csv(cfg.fig2_input_tsv, sep="\t")
    n_samples = df["sample_id"].nunique()
    if n_samples != cfg.expected_n_samples:
        raise ValueError(f"expected {cfg.expected_n_samples} PCA samples, found {n_samples}")
    expected_per_group = cfg.expected_n_samples // len(cfg.group_colors)
    group_counts = df.drop_duplicates("sample_id")["group"].value_counts().to_dict()
    for group in cfg.group_colors:
        if group_counts.get(group, 0) != expected_per_group:
            raise ValueError(
                f"expected {expected_per_group} PCA samples for group {group!r}, found {group_counts.get(group, 0)}"
            )
    logger.info("load_fig2_input: read %d samples", n_samples)
    return df


def load_fig3_input(cfg: NebulaPlotsFinalConfig) -> pd.DataFrame:
    df = pd.read_csv(cfg.fig3_input_tsv, sep="\t")
    if df.empty:
        raise ValueError("fig3 volcano input is empty")
    logger.info("load_fig3_input: read %d background genes", len(df))
    return df


def load_fig5_summary_input(cfg: NebulaPlotsFinalConfig) -> pd.DataFrame:
    df = pd.read_csv(cfg.fig5_summary_input_tsv, sep="\t")
    if len(df) != 1 or df.iloc[0]["gene_symbol"] != cfg.expected_primary_gene:
        raise ValueError(f"expected exactly one {cfg.expected_primary_gene} row in fig5 summary input")
    return df


def load_fig5_expression_input(cfg: NebulaPlotsFinalConfig) -> pd.DataFrame:
    df = pd.read_csv(cfg.fig5_expression_input_tsv, sep="\t")
    if df["sample_id"].nunique() != cfg.expected_n_samples:
        raise ValueError(f"expected {cfg.expected_n_samples} samples in fig5 expression input, found {df['sample_id'].nunique()}")
    return df


def load_fig4_evidence_input(cfg: NebulaPlotsFinalConfig) -> pd.DataFrame:
    df = pd.read_csv(cfg.fig4_evidence_input_tsv, sep="\t")
    if len(df) != cfg.expected_n_sensitising:
        raise ValueError(f"expected {cfg.expected_n_sensitising} sensitising genes, found {len(df)}")
    if BENCHMARK_GENE_SYMBOL in set(df["gene_symbol"]):
        raise ValueError(f"{BENCHMARK_GENE_SYMBOL} must never appear among the sensitising genes")
    primary_rows = df.loc[df["evidence_class"] == EVIDENCE_CLASS_PRIMARY, "gene_symbol"].tolist()
    if primary_rows != [cfg.expected_primary_gene]:
        raise ValueError(f"expected PRIMARY_RESISTANCE_SUPPORT to be exactly [{cfg.expected_primary_gene!r}], found {primary_rows!r}")
    logger.info("load_fig4_evidence_input: read %d rows, sole PRIMARY=%s", len(df), primary_rows[0])
    return df


def load_filtered_tpm(cfg: NebulaPlotsFinalConfig) -> pd.DataFrame:
    actual_sha256 = _sha256_file(cfg.filtered_gene_tpm_tsv_gz)
    if actual_sha256 != cfg.frozen_filtered_gene_tpm_sha256:
        raise ValueError(
            f"filtered GSE118713 TPM matrix checksum mismatch: expected "
            f"{cfg.frozen_filtered_gene_tpm_sha256}, got {actual_sha256}"
        )
    df = pd.read_csv(cfg.filtered_gene_tpm_tsv_gz, sep="\t")
    logger.info("load_filtered_tpm: read %d genes (checksum verified)", len(df))
    return df


def load_sample_metadata(cfg: NebulaPlotsFinalConfig) -> pd.DataFrame:
    df = pd.read_csv(cfg.sample_metadata_tsv, sep="\t")
    if df["sample_id"].nunique() != cfg.expected_n_samples:
        raise ValueError(f"expected {cfg.expected_n_samples} samples in sample metadata, found {df['sample_id'].nunique()}")
    unexpected_groups = set(df["group"]) - set(_SAMPLE_GROUP_ORDER)
    if unexpected_groups:
        raise ValueError(f"sample metadata contains unexpected group labels: {unexpected_groups}")
    logger.info("load_sample_metadata: read %d samples", len(df))
    return df


def ordered_sample_columns(sample_meta_df: pd.DataFrame) -> list[str]:
    """Fixed column order: MCF7 Rep1-3, TAMR Rep1-3, FASR Rep1-3 -- no
    clustering, sorted only by the declared group order then replicate
    number."""
    ordered = []
    for group in _SAMPLE_GROUP_ORDER:
        sub = sample_meta_df.loc[sample_meta_df["group"] == group].sort_values("replicate")
        ordered.extend(sub["sample_id"].tolist())
    return ordered


def build_fig4_heatmap_inputs(
    evidence_df: pd.DataFrame, filtered_tpm_df: pd.DataFrame, sample_meta_df: pd.DataFrame
) -> tuple[pd.DataFrame, pd.Series, pd.Series]:
    """Return (z_score_df, available_mask, evidence_class), all indexed by
    gene_symbol in ``evidence_df``'s existing order (no re-sort).

    ``z_score_df`` holds the row-wise z-score of log2(TPM+1) across the 9
    samples, for display only -- never fed back into candidate selection
    or significance testing. A gene absent from the frozen, filtered TPM
    matrix (removed by the pre-registered expression filter) gets an
    all-NaN row and ``available=False``; its expression is never
    fabricated.
    """
    sample_columns = ordered_sample_columns(sample_meta_df)
    genes = evidence_df["gene_symbol"].tolist()

    tpm_by_gene = filtered_tpm_df.set_index("gene_symbol")
    z_rows = {}
    available = {}
    for gene in genes:
        if gene not in tpm_by_gene.index:
            z_rows[gene] = pd.Series(np.nan, index=sample_columns)
            available[gene] = False
            continue
        tpm_row = tpm_by_gene.loc[gene, sample_columns].to_numpy(dtype=float)
        if not np.isfinite(tpm_row).all() or (tpm_row < 0).any():
            raise ValueError(f"{gene}: non-finite or negative TPM values in the frozen filtered matrix")
        log2_row = np.log2(tpm_row + 1.0)
        std = log2_row.std(ddof=0)
        if std == 0:
            z_row = np.zeros_like(log2_row)
        else:
            z_row = (log2_row - log2_row.mean()) / std
        z_rows[gene] = pd.Series(z_row, index=sample_columns)
        available[gene] = True

    z_score_df = pd.DataFrame(z_rows).T.loc[genes]
    available_mask = pd.Series(available).loc[genes]
    evidence_class = evidence_df.set_index("gene_symbol")["evidence_class"].loc[genes]

    n_unavailable = int((~available_mask).sum())
    logger.info(
        "build_fig4_heatmap_inputs: %d genes with available expression, %d unavailable (%s)",
        int(available_mask.sum()),
        n_unavailable,
        ", ".join(available_mask.index[~available_mask]) if n_unavailable else "none",
    )
    return z_score_df, available_mask, evidence_class


def contiguous_evidence_blocks(evidence_class: pd.Series) -> list[tuple[int, int, str]]:
    """Collapse ``evidence_class`` (already in the frozen, fixed gene
    order) into contiguous (start, end, class) row blocks for direct
    on-figure labelling.

    Raises if the same class appears in more than one non-adjacent run,
    since the block-label design assumes -- and visually asserts -- that
    each class occupies one contiguous band of rows.
    """
    classes = evidence_class.tolist()
    if not classes:
        return []
    blocks: list[tuple[int, int, str]] = []
    start = 0
    for i in range(1, len(classes) + 1):
        if i == len(classes) or classes[i] != classes[start]:
            blocks.append((start, i, classes[start]))
            start = i
    seen = set()
    for _, _, cls in blocks:
        if cls in seen:
            raise ValueError(f"evidence_class {cls!r} is not block-contiguous in the frozen gene order")
        seen.add(cls)
    return blocks


def run_qc_checks_fig123(
    fig1_df: pd.DataFrame,
    fig2_df: pd.DataFrame,
    fig4_df: pd.DataFrame,
    paics_df: pd.DataFrame,
    cfg: NebulaPlotsFinalConfig,
) -> None:
    """Programmatic QC gate for Figures 1-3, run before any figure is
    written. Pure structural re-verification of the frozen plot inputs --
    no new computation, no new threshold."""
    if len(fig1_df) != cfg.expected_n_hits:
        raise ValueError(f"expected {cfg.expected_n_hits} Gate-1 hits, found {len(fig1_df)}")
    if BENCHMARK_GENE_SYMBOL in set(fig1_df["gene_symbol"]):
        raise ValueError(f"{BENCHMARK_GENE_SYMBOL} must never be counted among the Gate-1 hits")
    if fig2_df["sample_id"].nunique() != cfg.expected_n_samples:
        raise ValueError(f"expected {cfg.expected_n_samples} PCA samples, found {fig2_df['sample_id'].nunique()}")
    if len(fig4_df) != cfg.expected_n_sensitising:
        raise ValueError(f"expected {cfg.expected_n_sensitising} sensitising genes, found {len(fig4_df)}")
    primary_rows = fig4_df.loc[fig4_df["evidence_class"] == EVIDENCE_CLASS_PRIMARY, "gene_symbol"].tolist()
    if primary_rows != [cfg.expected_primary_gene]:
        raise ValueError(f"expected PRIMARY_RESISTANCE_SUPPORT to be exactly [{cfg.expected_primary_gene!r}], found {primary_rows!r}")
    if paics_df.iloc[0]["gene_symbol"] != BENCHMARK_GENE_SYMBOL or paics_df.iloc[0]["benchmark_label"] != BENCHMARK_LABEL:
        raise ValueError("PAICS benchmark row is malformed or mislabelled")
    logger.info(
        "run_qc_checks_fig123: PASSED (n_hits=%d, n_pca_samples=%d, n_sensitising=%d, primary=%s, PAICS benchmark-only)",
        len(fig1_df), fig2_df["sample_id"].nunique(), len(fig4_df), primary_rows[0],
    )


def run_qc_checks_fig4(
    evidence_df: pd.DataFrame,
    sample_meta_df: pd.DataFrame,
    z_score_df: pd.DataFrame,
    available_mask: pd.Series,
    cfg: NebulaPlotsFinalConfig,
) -> None:
    """Programmatic QC gate for Figure 4's heatmap data, run before the
    figure is written."""
    if len(evidence_df) != cfg.expected_n_sensitising:
        raise ValueError(f"expected {cfg.expected_n_sensitising} sensitising genes, found {len(evidence_df)}")
    if BENCHMARK_GENE_SYMBOL in set(evidence_df["gene_symbol"]):
        raise ValueError(f"{BENCHMARK_GENE_SYMBOL} must never be counted among the sensitising genes")
    primary_rows = evidence_df.loc[evidence_df["evidence_class"] == EVIDENCE_CLASS_PRIMARY, "gene_symbol"].tolist()
    if primary_rows != [cfg.expected_primary_gene]:
        raise ValueError(f"expected PRIMARY_RESISTANCE_SUPPORT to be exactly [{cfg.expected_primary_gene!r}], found {primary_rows!r}")

    if sample_meta_df["sample_id"].nunique() != cfg.expected_n_samples:
        raise ValueError(f"expected {cfg.expected_n_samples} samples, found {sample_meta_df['sample_id'].nunique()}")
    expected_per_group = cfg.expected_n_samples // len(_SAMPLE_GROUP_ORDER)
    group_counts = sample_meta_df["group"].value_counts().to_dict()
    for group in _SAMPLE_GROUP_ORDER:
        if group_counts.get(group, 0) != expected_per_group:
            raise ValueError(f"expected {expected_per_group} samples for group {group!r}, found {group_counts.get(group, 0)}")

    if list(z_score_df.columns) != ordered_sample_columns(sample_meta_df):
        raise ValueError("heatmap columns are not in the declared MCF7/TAMR/FASR replicate order")
    if list(z_score_df.index) != list(evidence_df["gene_symbol"]):
        raise ValueError("heatmap rows do not match the frozen evidence table's existing gene order")

    for gene, is_available in available_mask.items():
        row_notna = z_score_df.loc[gene].notna()
        if is_available and not row_notna.all():
            raise ValueError(f"{gene}: marked available but has missing per-sample z-score values")
        if not is_available and row_notna.any():
            raise ValueError(f"{gene}: marked unavailable but has non-NaN z-score values")

    logger.info(
        "run_qc_checks_fig4: PASSED (n_sensitising=%d, n_samples=%d, 3/group, primary=%s, %d genes available, %d unavailable)",
        len(evidence_df), sample_meta_df["sample_id"].nunique(), primary_rows[0],
        int(available_mask.sum()), int((~available_mask).sum()),
    )


# --------------------------------------------------------------------------
# Shared style helpers
# --------------------------------------------------------------------------


def _new_axes(cfg: NebulaPlotsFinalConfig, figsize: tuple[float, float]) -> tuple[plt.Figure, plt.Axes]:
    fig, ax = plt.subplots(figsize=figsize)
    fig.patch.set_facecolor(cfg.panel_background)
    ax.set_facecolor(cfg.panel_background)
    for spine in ax.spines.values():
        spine.set_color(cfg.palette["neutral_grey"])
        spine.set_linewidth(0.8)
    ax.tick_params(colors=_DARK_TEXT, labelsize=9)
    ax.xaxis.label.set_color(_DARK_TEXT)
    ax.yaxis.label.set_color(_DARK_TEXT)
    return fig, ax


def _title_block(ax: plt.Axes, title: str, subtitle: str | None = None) -> None:
    ax.text(0.0, 1.12, title, transform=ax.transAxes, ha="left", va="bottom", fontsize=14, fontweight="bold", color=_DARK_TEXT)
    if subtitle:
        ax.text(0.0, 1.045, subtitle, transform=ax.transAxes, ha="left", va="bottom", fontsize=9.5, color=_MUTED_TEXT)


def _save_figure(fig: plt.Figure, png_path: Path, pdf_path: Path) -> None:
    for path in (png_path, pdf_path):
        path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(png_path, dpi=300, facecolor=fig.get_facecolor(), bbox_inches="tight")
    fig.savefig(pdf_path, metadata={"CreationDate": None}, bbox_inches="tight")
    plt.close(fig)
    logger.info("_save_figure: wrote %s, %s", png_path, pdf_path)


# --------------------------------------------------------------------------
# Shared row/column geometry for Figure 4 -- approved v8 alignment fix.
# --------------------------------------------------------------------------


def fig4_row_geometry(n_genes: int) -> dict[str, object]:
    """One authoritative row-coordinate system for Figure 4, matching the
    convention ``imshow`` uses internally: row ``i`` is centered at
    data-coordinate ``i`` and spans ``[i - 0.5, i + 0.5]``. Every
    row-indexed axis (evidence labels, gene names, heatmap, row
    annotations) uses this exact ``ylim``."""
    centers = np.arange(n_genes, dtype=float)
    edges = np.arange(n_genes + 1, dtype=float) - 0.5
    return {"centers": centers, "edges": edges, "ylim": (float(edges[-1]), float(edges[0]))}


def fig4_col_geometry(n_samples: int) -> dict[str, object]:
    """Column analogue of ``fig4_row_geometry`` -- keeps the top
    MCF7/TAMR/FASR group strip's column boundaries numerically identical
    to the heatmap's own column boundaries."""
    centers = np.arange(n_samples, dtype=float)
    edges = np.arange(n_samples + 1, dtype=float) - 0.5
    return {"centers": centers, "edges": edges, "xlim": (float(edges[0]), float(edges[-1]))}


def fig4_block_boundaries(blocks: list[tuple[int, int, str]], row_edges: np.ndarray) -> list[float]:
    """Internal (non-outer) y-positions where block separator lines must
    be drawn, derived from the same ``row_edges`` used by every axis."""
    return [float(row_edges[end]) for _start, end, _cls in blocks[:-1]]


def fig4_sample_group_boundaries(n_samples: int, col_edges: np.ndarray, group_size: int = 3) -> list[float]:
    """Column-boundary analogue of ``fig4_block_boundaries``, for the
    MCF7/TAMR/FASR group divisions (every ``group_size`` columns)."""
    return [float(col_edges[i]) for i in range(group_size, n_samples, group_size)]


# --------------------------------------------------------------------------
# Rendered-bbox overlap checking (approved v7/v8 layout-verification tool).
# --------------------------------------------------------------------------


def _text_bbox(text: Text, renderer) -> Bbox:
    return text.get_window_extent(renderer=renderer)


def _bboxes_overlap(a: Bbox, b: Bbox) -> bool:
    return a.x0 < b.x1 and a.x1 > b.x0 and a.y0 < b.y1 and a.y1 > b.y0


def check_fig4_layout(fig: plt.Figure, artifacts: dict) -> list[str]:
    """Verify, using real renderer bounding boxes, that no text region in
    Figure 4 visually enters another, and that every row-based axis shares
    the identical y-limit and identical separator positions."""
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    violations: list[str] = []

    evidence_boxes = [(t.get_text(), _text_bbox(t, renderer)) for t in artifacts["evidence_label_texts"]]
    gene_boxes = [(t.get_text(), _text_bbox(t, renderer)) for t in artifacts["gene_name_texts"]]
    usp34_boxes = [_text_bbox(t, renderer) for t in artifacts["usp34_annotation_texts"]]
    usp17l29_boxes = [_text_bbox(t, renderer) for t in artifacts["usp17l29_annotation_texts"]]
    heatmap_bbox = artifacts["heatmap_ax"].get_window_extent(renderer=renderer)
    colorbar_bbox = artifacts["colorbar_ax"].get_window_extent(renderer=renderer)
    fig_bbox = fig.get_window_extent(renderer=renderer)

    for e_label, e_box in evidence_boxes:
        for g_label, g_box in gene_boxes:
            if _bboxes_overlap(e_box, g_box):
                violations.append(f"evidence label {e_label!r} overlaps gene name {g_label!r}")

    for box in usp34_boxes:
        if _bboxes_overlap(box, heatmap_bbox):
            violations.append("USP34 annotation overlaps heatmap")
        if _bboxes_overlap(box, colorbar_bbox):
            violations.append("USP34 annotation overlaps colorbar")

    for box in usp17l29_boxes:
        if _bboxes_overlap(box, heatmap_bbox):
            violations.append("USP17L29 annotation overlaps heatmap")
        if _bboxes_overlap(box, colorbar_bbox):
            violations.append("USP17L29 annotation overlaps colorbar")

    all_texts = evidence_boxes + gene_boxes
    all_boxes = [b for _l, b in all_texts] + usp34_boxes + usp17l29_boxes
    for box in all_boxes:
        if box.x0 < fig_bbox.x0 - 1 or box.x1 > fig_bbox.x1 + 1 or box.y0 < fig_bbox.y0 - 1 or box.y1 > fig_bbox.y1 + 1:
            violations.append("a text box extends outside the figure canvas")

    for i in range(len(evidence_boxes)):
        for j in range(i + 1, len(evidence_boxes)):
            if _bboxes_overlap(evidence_boxes[i][1], evidence_boxes[j][1]):
                violations.append(f"evidence label {evidence_boxes[i][0]!r} overlaps {evidence_boxes[j][0]!r}")

    row_ylims = artifacts["row_ylims"]
    reference = next(iter(row_ylims.values()))
    for name, ylim in row_ylims.items():
        if ylim != reference:
            violations.append(f"axis {name!r} ylim {ylim} differs from reference {reference}")

    separator_positions = artifacts["separator_positions"]
    reference_positions = next(iter(separator_positions.values()))
    for name, positions in separator_positions.items():
        if positions != reference_positions:
            violations.append(f"separator positions in {name!r} ({positions}) differ from reference ({reference_positions})")

    col_boundaries = artifacts["col_boundaries"]
    if col_boundaries["group_strip"] != col_boundaries["heat"]:
        violations.append(
            f"group-strip column separators ({col_boundaries['group_strip']}) differ from heatmap column separators ({col_boundaries['heat']})"
        )

    return violations


def check_fig5_layout(fig: plt.Figure, artifacts: dict) -> list[str]:
    """Verify, using real renderer bounding boxes, that Figure 5's
    annotations do not collide with each other, the title, or the data."""
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    violations: list[str] = []

    mcf7_tamr_box = artifacts["mcf7_tamr_bracket_text"].get_window_extent(renderer=renderer)
    tamr_fasr_box = artifacts["tamr_fasr_bracket_text"].get_window_extent(renderer=renderer)
    crispr_box = artifacts["crispr_box_text"].get_window_extent(renderer=renderer)
    title_box = artifacts["title_text"].get_window_extent(renderer=renderer)
    footer_box = artifacts["footer_text"].get_window_extent(renderer=renderer)
    fig_bbox = fig.get_window_extent(renderer=renderer)

    if _bboxes_overlap(mcf7_tamr_box, tamr_fasr_box):
        violations.append("TAMR-vs-MCF7 bracket overlaps TAMR-vs-FASR bracket")
    if _bboxes_overlap(mcf7_tamr_box, title_box):
        violations.append("TAMR-vs-MCF7 bracket overlaps title")
    if _bboxes_overlap(tamr_fasr_box, title_box):
        violations.append("TAMR-vs-FASR bracket overlaps title")
    if _bboxes_overlap(crispr_box, mcf7_tamr_box):
        violations.append("CRISPR box overlaps TAMR-vs-MCF7 bracket")
    if _bboxes_overlap(crispr_box, tamr_fasr_box):
        violations.append("CRISPR box overlaps TAMR-vs-FASR bracket")

    for name, box in [
        ("mcf7_tamr_bracket", mcf7_tamr_box), ("tamr_fasr_bracket", tamr_fasr_box),
        ("crispr_box", crispr_box), ("title", title_box), ("footer", footer_box),
    ]:
        if box.x0 < fig_bbox.x0 - 1 or box.x1 > fig_bbox.x1 + 1 or box.y0 < fig_bbox.y0 - 1 or box.y1 > fig_bbox.y1 + 1:
            violations.append(f"{name} extends outside the figure canvas")

    return violations


# --------------------------------------------------------------------------
# Figure 1: CRISPR landscape -- approved v6 design
# --------------------------------------------------------------------------


def plot_fig1(fig1_df: pd.DataFrame, paics_row: pd.Series, cfg: NebulaPlotsFinalConfig, primary_gene: str) -> plt.Figure:
    df = fig1_df.sort_values(["crispr_effect_size", "gene_symbol"]).reset_index(drop=True)
    fig, ax = _new_axes(cfg, (10, 9.5))

    y = np.arange(len(df))
    colors = [cfg.direction_colors[d] for d in df["crispr_direction"]]
    ax.hlines(y, 0, df["crispr_effect_size"], color=colors, linewidth=3, zorder=2)
    sizes = 25 + 160 * np.clip(-np.log10(df["crispr_fdr"].to_numpy(dtype=float)) / 7.0, 0, 1)
    ax.scatter(df["crispr_effect_size"], y, s=sizes, color=colors, edgecolor="white", linewidth=0.6, zorder=3)
    ax.axvline(0, color=_DARK_TEXT, linewidth=1.0, zorder=1)

    ax.set_yticks(y)
    labels = ax.set_yticklabels(df["gene_symbol"], fontsize=9)
    for tick_label, gene in zip(labels, df["gene_symbol"]):
        if gene == primary_gene:
            tick_label.set_fontweight("bold")
            tick_label.set_color(cfg.direction_colors[DIRECTION_SENSITISING])
            tick_label.set_fontsize(11)
        else:
            tick_label.set_color(_DARK_TEXT)
    ax.tick_params(axis="y", length=0)
    ax.set_xlabel("CRISPR effect size")
    ax.set_ylim(-1.3, len(df))
    ax.invert_yaxis()

    _title_block(ax, FIG1_TITLE)

    ax.text(0.25, 1.015, "←  Sensitising knockout", transform=ax.transAxes, ha="center", va="bottom", fontsize=10.5, fontweight="bold", color=cfg.direction_colors[DIRECTION_SENSITISING])
    ax.text(0.75, 1.015, "Tolerance-associated knockout  →", transform=ax.transAxes, ha="center", va="bottom", fontsize=10.5, fontweight="bold", color=cfg.direction_colors[DIRECTION_TOLERANCE])

    ax.text(0.015, 0.02, "Larger dots = lower CRISPR FDR", transform=ax.transAxes, ha="left", va="bottom", fontsize=8, color=_MUTED_TEXT, style="italic")

    paics_text = (
        f"PAICS — published benchmark\n"
        f"effect {_fmt_effect(paics_row['crispr_effect_size'])} · FDR {_fmt_fdr(paics_row['crispr_fdr'])} · not a Gate-1 hit"
    )
    ax.text(0.985, 0.985, paics_text, transform=ax.transAxes, ha="right", va="top", fontsize=7.5, color=_MUTED_TEXT, style="italic")

    fig.tight_layout()
    return fig


# --------------------------------------------------------------------------
# Figure 2: PCA -- approved v6 design
# --------------------------------------------------------------------------


def plot_fig2(fig2_df: pd.DataFrame, cfg: NebulaPlotsFinalConfig) -> plt.Figure:
    fig, ax = _new_axes(cfg, (6.5, 5.5))
    centroids = {}
    for group, sub in fig2_df.groupby("group", sort=False):
        ax.scatter(sub["PC1"], sub["PC2"], color=cfg.group_colors[group], s=130, edgecolor="white", linewidth=1.0, zorder=3)
        centroids[group] = (sub["PC1"].mean(), sub["PC2"].mean())

    ax.axhline(0, color=cfg.palette["neutral_grey"], linewidth=0.5, zorder=1, alpha=0.6)
    ax.axvline(0, color=cfg.palette["neutral_grey"], linewidth=0.5, zorder=1, alpha=0.6)
    ax.margins(0.28)

    offsets = {"MCF7": (0, 24), "TAMR": (0, -28), "FASR": (0, 24)}
    for group, (cx, cy) in centroids.items():
        dx, dy = offsets.get(group, (0, 20))
        ax.annotate(
            group, xy=(cx, cy), xytext=(dx, dy), textcoords="offset points", ha="center", va="center",
            fontsize=13, fontweight="bold", color=cfg.group_colors[group],
        )

    ax.set_xlabel(f"PC1 ({fig2_df['pc1_variance_explained_pct'].iloc[0]:.1f}%)")
    ax.set_ylabel(f"PC2 ({fig2_df['pc2_variance_explained_pct'].iloc[0]:.1f}%)")
    _title_block(ax, FIG2_TITLE, "GSE118713 · n=3 per group")
    fig.tight_layout()
    return fig


# --------------------------------------------------------------------------
# Figure 3: TAMR_vs_MCF7 volcano -- approved v6 design
# --------------------------------------------------------------------------


def plot_fig3(fig3_df: pd.DataFrame, cfg: NebulaPlotsFinalConfig, primary_gene: str) -> plt.Figure:
    fig, ax = _new_axes(cfg, (8, 6.5))

    background = fig3_df.loc[~fig3_df["is_gate1_hit"]]
    ax.scatter(background["log2fc"], -np.log10(background["fdr"]), s=6, color=cfg.palette["neutral_grey"], alpha=0.35, linewidth=0, zorder=1, label="Other tested genes")

    hits = fig3_df.loc[fig3_df["is_gate1_hit"]]
    for direction, sub in hits.groupby("crispr_direction"):
        label = "Sensitising knockout" if direction == DIRECTION_SENSITISING else "Tolerance-associated knockout"
        ax.scatter(sub["log2fc"], -np.log10(sub["fdr"]), s=40, color=cfg.direction_colors[direction], edgecolor="white", linewidth=0.5, zorder=2, label=label)

    threshold_y = -np.log10(cfg.rna_significance_fdr)
    ax.axhline(threshold_y, color=cfg.palette["cosmic_magenta"], linestyle="--", linewidth=1.0, zorder=1)
    ax.text(0.99, threshold_y, f"FDR = {cfg.rna_significance_fdr:g}", transform=ax.get_yaxis_transform(), color=cfg.palette["cosmic_magenta"], fontsize=7.5, va="bottom", ha="right")

    usp34_row = fig3_df.loc[fig3_df["gene_symbol"] == primary_gene]
    if not usp34_row.empty:
        r = usp34_row.iloc[0]
        x0, y0 = float(r["log2fc"]), float(-np.log10(r["fdr"]))
        ax.scatter([x0], [y0], s=140, facecolor="none", edgecolor=cfg.direction_colors[DIRECTION_SENSITISING], linewidth=2.2, zorder=4)
        callout = f"{primary_gene}\nlog2FC {_fmt_effect(r['log2fc'])}\nFDR {_fmt_fdr(r['fdr'])}"
        ax.annotate(
            callout, xy=(x0, y0), xytext=(28, 22), textcoords="offset points", ha="left", va="bottom",
            fontsize=9.5, fontweight="bold", color=cfg.direction_colors[DIRECTION_SENSITISING],
            arrowprops=dict(arrowstyle="-", color=cfg.direction_colors[DIRECTION_SENSITISING], lw=1.0),
        )

    ax.set_xlabel("log2FC (TAMR vs MCF7)")
    ax.set_ylabel("-log10(FDR)")
    _title_block(ax, "Expression changes associated with tamoxifen resistance", "TAMR vs parental MCF7")

    handles, labels_ = ax.get_legend_handles_labels()
    seen = set()
    dedup_handles, dedup_labels = [], []
    for h, l in zip(handles, labels_):
        if l not in seen:
            dedup_handles.append(h)
            dedup_labels.append(l)
            seen.add(l)
    leg = ax.legend(dedup_handles, dedup_labels, loc="upper left", fontsize=8.5, frameon=False)
    for text in leg.get_texts():
        text.set_color(_DARK_TEXT)
    fig.tight_layout()
    return fig


# --------------------------------------------------------------------------
# Figure 4: candidate expression heatmap -- approved v8 design
# --------------------------------------------------------------------------


def plot_fig4(
    z_score_df: pd.DataFrame,
    available_mask: pd.Series,
    evidence_class: pd.Series,
    sample_meta_df: pd.DataFrame,
    fig5_summary_df: pd.DataFrame,
    cfg: NebulaPlotsFinalConfig,
    primary_gene: str,
) -> tuple[plt.Figure, dict]:
    genes = list(z_score_df.index)
    samples = list(z_score_df.columns)
    n_genes = len(genes)
    n_samples = len(samples)
    sample_group = sample_meta_df.set_index("sample_id")["group"]
    blocks = contiguous_evidence_blocks(evidence_class)

    row_geom = fig4_row_geometry(n_genes)
    col_geom = fig4_col_geometry(n_samples)
    row_boundaries = fig4_block_boundaries(blocks, row_geom["edges"])
    col_boundaries = fig4_sample_group_boundaries(n_samples, col_geom["edges"])

    fig = plt.figure(figsize=(16.5, 0.62 * n_genes + 3.0))
    fig.patch.set_facecolor(cfg.panel_background)
    gs = fig.add_gridspec(
        2, 5,
        width_ratios=[2.5, 1.5, n_samples, 2.3, 0.6],
        height_ratios=[1.0, n_genes],
        wspace=0.08, hspace=0.05,
    )
    ax_blank_a = fig.add_subplot(gs[0, 0])
    ax_blank_b = fig.add_subplot(gs[0, 1])
    ax_group = fig.add_subplot(gs[0, 2])
    ax_blank_d = fig.add_subplot(gs[0, 3])
    ax_cbar_header = fig.add_subplot(gs[0, 4])

    ax_evidence = fig.add_subplot(gs[1, 0])
    ax_genes = fig.add_subplot(gs[1, 1])
    ax_heat = fig.add_subplot(gs[1, 2])
    ax_annot = fig.add_subplot(gs[1, 3])
    ax_cbar = fig.add_subplot(gs[1, 4])

    for ax in (ax_blank_a, ax_blank_b, ax_group, ax_blank_d, ax_cbar_header, ax_evidence, ax_genes, ax_annot):
        ax.set_facecolor(cfg.panel_background)
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)
    for ax in (ax_blank_a, ax_blank_b, ax_blank_d, ax_cbar_header, ax_evidence, ax_genes, ax_annot):
        ax.axis("off")

    row_ylim = row_geom["ylim"]
    ax_evidence.set_xlim(0, 1)
    ax_evidence.set_ylim(*row_ylim)
    ax_genes.set_xlim(0, 1)
    ax_genes.set_ylim(*row_ylim)
    ax_annot.set_xlim(0, 1)
    ax_annot.set_ylim(*row_ylim)

    # --- Region A: evidence-class labels (swatch + wrapped label) --------
    evidence_label_texts: list[Text] = []
    evidence_separator_lines: list[Line2D] = []
    row_edges = row_geom["edges"]
    for start, end, cls in blocks:
        center = (row_edges[start] + row_edges[end]) / 2
        block_height = row_edges[end] - row_edges[start]
        swatch_pad = min(0.15, block_height * 0.15)
        ax_evidence.add_patch(
            Rectangle((0.04, row_edges[start] + swatch_pad), 0.09, block_height - 2 * swatch_pad, color=_EVIDENCE_CLASS_ACCENT[cls], linewidth=0)
        )
        text = ax_evidence.text(
            0.22, center, _EVIDENCE_CLASS_PLAIN_LABEL[cls], ha="left", va="center",
            fontsize=10, fontweight="bold", color=_DARK_TEXT, linespacing=1.35,
        )
        evidence_label_texts.append(text)
    for y in row_boundaries:
        line = ax_evidence.axhline(y, color=cfg.palette["neutral_grey"], linewidth=0.7, alpha=0.6)
        evidence_separator_lines.append(line)

    # --- Region B: gene names (own dedicated column) ----------------------
    gene_name_texts: list[Text] = []
    genes_separator_lines: list[Line2D] = []
    row_centers = row_geom["centers"]
    for i, gene in enumerate(genes):
        y = row_centers[i]
        if gene == primary_gene:
            text = ax_genes.text(0.05, y, gene, ha="left", va="center", fontsize=13, fontweight="bold", color=cfg.direction_colors[DIRECTION_SENSITISING])
        elif gene == "USP17L29":
            text = ax_genes.text(0.05, y, gene, ha="left", va="center", fontsize=10.5, color=cfg.palette["neutral_grey"])
        else:
            text = ax_genes.text(0.05, y, gene, ha="left", va="center", fontsize=10.5, color=_DARK_TEXT)
        gene_name_texts.append(text)
    for y in row_boundaries:
        line = ax_genes.axhline(y, color=cfg.palette["neutral_grey"], linewidth=0.7, alpha=0.6)
        genes_separator_lines.append(line)

    # --- Region C: heatmap (blue low / cream mid / orange-red high) ------
    cmap = LinearSegmentedColormap.from_list("nebula_final_expression", ["#2F6FED", "#FBF4E8", "#F2650B"])
    z_values = np.ma.masked_invalid(z_score_df.to_numpy(dtype=float))
    vmax = float(np.nanmax(np.abs(z_score_df.to_numpy(dtype=float)))) if available_mask.any() else 1.0
    vmax = min(vmax, 2.5) if vmax > 0 else 1.0
    col_edges = col_geom["edges"]
    heat_extent = (float(col_edges[0]), float(col_edges[-1]), float(row_edges[-1]), float(row_edges[0]))
    im = ax_heat.imshow(
        z_values, aspect="auto", cmap=cmap.with_extremes(bad=_UNAVAILABLE_GREY), vmin=-vmax, vmax=vmax,
        interpolation="none", extent=heat_extent,
    )
    ax_heat.set_xlim(*col_geom["xlim"])
    ax_heat.set_ylim(*row_ylim)

    ax_heat.set_xticks(col_geom["centers"])
    ax_heat.set_xticklabels([str(i % 3 + 1) for i in range(n_samples)], fontsize=8, color=_MUTED_TEXT)
    ax_heat.set_yticks([])
    ax_heat.tick_params(axis="x", length=0, pad=3)
    for spine in ax_heat.spines.values():
        spine.set_visible(False)

    heat_col_separator_lines: list[Line2D] = []
    for x in col_boundaries:
        line = ax_heat.axvline(x, color=cfg.panel_background, linewidth=3, zorder=5)
        heat_col_separator_lines.append(line)
    heat_row_separator_lines: list[Line2D] = []
    for y in row_boundaries:
        line = ax_heat.axhline(y, color=cfg.panel_background, linewidth=4, zorder=5)
        heat_row_separator_lines.append(line)

    if primary_gene in genes:
        idx = genes.index(primary_gene)
        ax_heat.add_patch(
            Rectangle((col_edges[0], row_edges[idx]), col_edges[-1] - col_edges[0], row_edges[idx + 1] - row_edges[idx],
                      fill=False, edgecolor=cfg.direction_colors[DIRECTION_SENSITISING], linewidth=2.8, zorder=6)
        )

    # --- Top strip: sample groups -- shares col_geom with the heatmap ------
    ax_group.set_facecolor(cfg.panel_background)
    ax_group.set_xticks([])
    ax_group.set_yticks([])
    for spine in ax_group.spines.values():
        spine.set_visible(False)
    ax_group.set_xlim(*col_geom["xlim"])
    ax_group.set_ylim(0, 1)
    for i, sample in enumerate(samples):
        ax_group.add_patch(Rectangle((col_edges[i], 0), col_edges[i + 1] - col_edges[i], 1, color=cfg.group_colors[sample_group[sample]], linewidth=0))
    group_separator_lines: list[Line2D] = []
    for x in col_boundaries:
        line = ax_group.axvline(x, color=cfg.panel_background, linewidth=3, zorder=5)
        group_separator_lines.append(line)
    for group in _SAMPLE_GROUP_ORDER:
        idxs = [i for i, s in enumerate(samples) if sample_group[s] == group]
        center = (col_edges[min(idxs)] + col_edges[max(idxs) + 1]) / 2
        ax_group.text(center, 1.5, group, ha="center", va="bottom", fontsize=13, fontweight="bold", color=_DARK_TEXT)

    # --- Region D: row-specific annotations (own dedicated column) --------
    usp34_annotation_texts: list[Text] = []
    usp17l29_annotation_texts: list[Text] = []
    annot_separator_lines: list[Line2D] = []
    if primary_gene in genes and len(fig5_summary_df) == 1:
        idx = genes.index(primary_gene)
        row = fig5_summary_df.iloc[0]
        text = ax_annot.text(
            0.5, row_centers[idx], FIG4_USP34_ANNOTATION.format(fdr=_fmt_fdr(row["tamr_vs_mcf7_fdr"])), ha="center", va="center",
            fontsize=9, fontweight="bold", color=cfg.direction_colors[DIRECTION_SENSITISING], linespacing=1.4,
        )
        usp34_annotation_texts.append(text)
    if "USP17L29" in genes:
        idx = genes.index("USP17L29")
        text = ax_annot.text(
            0.5, row_centers[idx], FIG4_USP17L29_ANNOTATION, ha="center", va="center",
            fontsize=9, color=cfg.palette["neutral_grey"], style="italic", linespacing=1.4,
        )
        usp17l29_annotation_texts.append(text)
    for y in row_boundaries:
        line = ax_annot.axhline(y, color=cfg.palette["neutral_grey"], linewidth=0.7, alpha=0.6)
        annot_separator_lines.append(line)

    # --- Region E: colorbar (own dedicated header + colorbar column) ------
    ax_cbar_header.set_xlim(0, 1)
    ax_cbar_header.set_ylim(0, 1)
    ax_cbar_header.text(0.5, 0.97, "RELATIVE\nEXPRESSION", ha="center", va="top", fontsize=8, fontweight="bold", color=_DARK_TEXT, linespacing=1.25)
    ax_cbar_header.text(0.5, 0.32, "HIGH", ha="center", va="center", fontsize=8, color="#F2650B", fontweight="bold")
    ax_cbar_header.text(0.5, 0.10, "LOW", ha="center", va="center", fontsize=8, color="#2F6FED", fontweight="bold")

    cbar = fig.colorbar(im, cax=ax_cbar)
    cbar.ax.tick_params(labelsize=7.5, colors=_MUTED_TEXT)
    cbar.outline.set_visible(False)
    cbar.set_label("Row z-score of\nlog2(TPM+1)", fontsize=7.5, color=_MUTED_TEXT)

    fig.text(0.015, 0.995, FIG4_TITLE, ha="left", va="top", fontsize=16, fontweight="bold", color=_DARK_TEXT)
    fig.text(0.015, 0.965, FIG4_SUBTITLE, ha="left", va="top", fontsize=10.5, color=_MUTED_TEXT)
    fig.text(0.015, 0.01, FIG4_FOOTER, ha="left", va="bottom", fontsize=7.5, color=_MUTED_TEXT, style="italic")

    artifacts = {
        "evidence_label_texts": evidence_label_texts,
        "gene_name_texts": gene_name_texts,
        "usp34_annotation_texts": usp34_annotation_texts,
        "usp17l29_annotation_texts": usp17l29_annotation_texts,
        "heatmap_ax": ax_heat,
        "colorbar_ax": ax_cbar,
        "row_ylims": {
            "evidence": ax_evidence.get_ylim(),
            "genes": ax_genes.get_ylim(),
            "heat": ax_heat.get_ylim(),
            "annot": ax_annot.get_ylim(),
        },
        "separator_positions": {
            "evidence": [round(line.get_ydata()[0], 9) for line in evidence_separator_lines],
            "genes": [round(line.get_ydata()[0], 9) for line in genes_separator_lines],
            "heat": [round(line.get_ydata()[0], 9) for line in heat_row_separator_lines],
            "annot": [round(line.get_ydata()[0], 9) for line in annot_separator_lines],
        },
        "col_boundaries": {
            "group_strip": [round(line.get_xdata()[0], 9) for line in group_separator_lines],
            "heat": [round(line.get_xdata()[0], 9) for line in heat_col_separator_lines],
        },
        "row_geometry": row_geom,
        "col_geometry": col_geom,
    }
    return fig, artifacts


# --------------------------------------------------------------------------
# Figure 5: USP34 expression -- approved v8 design
# --------------------------------------------------------------------------


def plot_fig5(summary_df: pd.DataFrame, expression_df: pd.DataFrame, cfg: NebulaPlotsFinalConfig, primary_gene: str) -> tuple[plt.Figure, dict]:
    row = summary_df.iloc[0]
    fig, ax = plt.subplots(figsize=(9, 5.6))
    fig.patch.set_facecolor(cfg.panel_background)
    ax.set_facecolor(cfg.panel_background)
    for spine in ax.spines.values():
        spine.set_color(cfg.palette["neutral_grey"])
        spine.set_linewidth(0.8)
    ax.tick_params(colors=_DARK_TEXT, labelsize=9)

    groups = list(_SAMPLE_GROUP_ORDER)
    for gi, group in enumerate(groups):
        sub = expression_df.loc[expression_df["group"] == group]
        jitter = np.linspace(-0.12, 0.12, len(sub))
        ax.scatter(gi + jitter, sub["log2_tpm_plus1"], color=cfg.group_colors[group], s=95, edgecolor="white", linewidth=0.9, zorder=3)
        mean_val = sub["log2_tpm_plus1"].mean()
        ax.hlines(mean_val, gi - 0.22, gi + 0.22, color=cfg.group_colors[group], linewidth=2.6, zorder=2)

    ax.set_xticks(range(len(groups)))
    ax.set_xticklabels(["MCF7", "TAMR", "FASR (context)"], fontsize=10.5)
    ax.set_ylabel("log2(TPM+1)", color=_DARK_TEXT)

    data_min = expression_df["log2_tpm_plus1"].min()
    data_max = expression_df["log2_tpm_plus1"].max()
    data_span = data_max - data_min
    ax.set_ylim(data_min - 0.12 * data_span, data_max + 0.42 * data_span)

    y_annot_1 = data_max + 0.10 * data_span
    ax.plot([0, 0, 1, 1], [y_annot_1, y_annot_1 + 0.02 * data_span, y_annot_1 + 0.02 * data_span, y_annot_1], color=_DARK_TEXT, linewidth=1.1)
    mcf7_tamr_text = ax.text(
        0.5, y_annot_1 + 0.035 * data_span, FIG5_TAMR_VS_MCF7_LABEL.format(fdr=_fmt_fdr(row["tamr_vs_mcf7_fdr"])),
        ha="center", va="bottom", fontsize=9.5, fontweight="bold", color=_DARK_TEXT, linespacing=1.3,
    )

    y_annot_2 = data_max + 0.27 * data_span
    ax.plot([1, 1, 2, 2], [y_annot_2, y_annot_2 + 0.02 * data_span, y_annot_2 + 0.02 * data_span, y_annot_2], color=_MUTED_TEXT, linewidth=0.9, linestyle="--")
    tamr_fasr_text = ax.text(
        1.5, y_annot_2 + 0.035 * data_span, FIG5_TAMR_VS_FASR_LABEL.format(fdr=_fmt_fdr(row["tamr_vs_fasr_fdr"])),
        ha="center", va="bottom", fontsize=8.5, color=_MUTED_TEXT, linespacing=1.3,
    )

    crispr_text = ax.text(
        0.985, 0.03, FIG5_CRISPR_BOX.format(effect=_fmt_effect(row["crispr_effect_size"]), fdr=_fmt_fdr(row["crispr_fdr"])),
        transform=ax.transAxes, ha="right", va="bottom", fontsize=9, color=_DARK_TEXT, linespacing=1.3,
        bbox=dict(boxstyle="round,pad=0.4", facecolor="white", edgecolor=cfg.palette["neutral_grey"], linewidth=0.8),
    )

    title_text = fig.text(0.02, 0.97, FIG5_TITLE, ha="left", va="top", fontsize=13.5, fontweight="bold", color=_DARK_TEXT)
    footer_text = fig.text(
        0.02, 0.015, FIG5_FOOTER, ha="left", va="bottom", fontsize=7.5, color=cfg.palette["neutral_grey"], style="italic",
    )

    fig.tight_layout(rect=(0, 0.08, 1, 0.87))

    artifacts = {
        "mcf7_tamr_bracket_text": mcf7_tamr_text,
        "tamr_fasr_bracket_text": tamr_fasr_text,
        "crispr_box_text": crispr_text,
        "title_text": title_text,
        "footer_text": footer_text,
    }
    return fig, artifacts


# --------------------------------------------------------------------------
# Contact sheet (all 5 final figures -- visual review, not a scientific output)
# --------------------------------------------------------------------------


def build_contact_sheet(cfg: NebulaPlotsFinalConfig) -> None:
    entries = [
        (cfg.fig1_png, "1. CRISPR landscape"),
        (cfg.fig2_png, "2. PCA"),
        (cfg.fig3_png, "3. Volcano"),
        (cfg.fig4_png, "4. Expression heatmap"),
        (cfg.fig5_png, "5. USP34"),
    ]
    images = [plt.imread(str(p)) for p, _ in entries]
    ratios = [im.shape[0] / im.shape[1] for im in images]
    total_width = 9.0
    fig = plt.figure(figsize=(total_width, sum(r * total_width for r in ratios) * 0.55 + 1.2))
    fig.patch.set_facecolor("white")
    gs = fig.add_gridspec(len(images), 1, height_ratios=ratios, hspace=0.25)
    for i, ((_path, label), im) in enumerate(zip(entries, images)):
        ax = fig.add_subplot(gs[i, 0])
        ax.imshow(im)
        ax.axis("off")
        ax.set_title(label, fontsize=11, loc="left", color=_DARK_TEXT, fontweight="bold")
    cfg.contact_sheet_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(cfg.contact_sheet_png, dpi=150, facecolor="white", bbox_inches="tight")
    plt.close(fig)
    logger.info("build_contact_sheet: wrote %s", cfg.contact_sheet_png)


# --------------------------------------------------------------------------
# Manifest
# --------------------------------------------------------------------------


def build_manifest(cfg: NebulaPlotsFinalConfig) -> pd.DataFrame:
    """Write ``figure_manifest.tsv``: one row per final PNG, recording its
    approved source version, pixel dimensions, and SHA-256 -- so the final
    output set is independently auditable without re-running matplotlib."""
    from PIL import Image

    rows = [
        {
            "figure": "fig1_crispr_landscape",
            "filename": cfg.fig1_png.name,
            "approved_source_version": "v6",
            "scientific_source": "results/tables/nebula_plot_inputs/fig1_crispr_hit_landscape_input.tsv",
            "purpose": "CRISPR effect-size landscape of the 28 Gate-1 hits (13 sensitising / 15 tolerance-associated), PAICS shown as a benchmark-only annotation.",
            "path": cfg.fig1_png,
        },
        {
            "figure": "fig2_pca",
            "filename": cfg.fig2_png.name,
            "approved_source_version": "v6",
            "scientific_source": "results/tables/nebula_plot_inputs/fig2_gse118713_pca_input.tsv",
            "purpose": "PCA of GSE118713 (9 samples, 3/group) showing MCF7/TAMR/FASR separation.",
            "path": cfg.fig2_png,
        },
        {
            "figure": "fig3_volcano",
            "filename": cfg.fig3_png.name,
            "approved_source_version": "v6",
            "scientific_source": "results/tables/nebula_plot_inputs/fig3_tamr_vs_mcf7_volcano_input.tsv",
            "purpose": "TAMR-vs-MCF7 volcano plot with USP34 callout.",
            "path": cfg.fig3_png,
        },
        {
            "figure": "fig4_candidate_expression_heatmap",
            "filename": cfg.fig4_png.name,
            "approved_source_version": "v8",
            "scientific_source": "results/tables/nebula_plot_inputs/fig4_sensitising_evidence_matrix_input.tsv + checksum-pinned filtered GSE118713 TPM matrix",
            "purpose": "Row-z-score expression heatmap of the 13 sensitising candidates across 9 samples, with plain-language evidence-class blocks.",
            "path": cfg.fig4_png,
        },
        {
            "figure": "fig5_usp34_expression",
            "filename": cfg.fig5_png.name,
            "approved_source_version": "v8",
            "scientific_source": "results/tables/nebula_plot_inputs/fig5_usp34_summary_input.tsv, fig5_usp34_expression_input.tsv",
            "purpose": "USP34 CRISPR + expression evidence panel (TAMR vs MCF7 significant, TAMR vs FASR not significant).",
            "path": cfg.fig5_png,
        },
    ]

    records = []
    for row in rows:
        image = Image.open(row["path"])
        width_px, height_px = image.size
        records.append(
            {
                "figure": row["figure"],
                "filename": row["filename"],
                "approved_source_version": row["approved_source_version"],
                "width_px": width_px,
                "height_px": height_px,
                "sha256": _sha256_file(row["path"]),
                "scientific_source": row["scientific_source"],
                "purpose": row["purpose"],
            }
        )
    manifest_df = pd.DataFrame(records)
    cfg.manifest_tsv.parent.mkdir(parents=True, exist_ok=True)
    manifest_df.to_csv(cfg.manifest_tsv, sep="\t", index=False)
    logger.info("build_manifest: wrote %s (%d rows)", cfg.manifest_tsv, len(manifest_df))
    return manifest_df


# --------------------------------------------------------------------------
# Orchestration
# --------------------------------------------------------------------------


def run_nebula_plots_final(config_path: str | Path = "config/config.yaml") -> dict[str, object]:
    config = _load_config(config_path)
    cfg = NebulaPlotsFinalConfig.from_config(config)

    fig1_df = load_fig1_input(cfg)
    paics_df = load_fig1_paics_inset(cfg)
    fig2_df = load_fig2_input(cfg)
    fig3_df = load_fig3_input(cfg)
    fig4_evidence_df = load_fig4_evidence_input(cfg)
    filtered_tpm_df = load_filtered_tpm(cfg)
    sample_meta_df = load_sample_metadata(cfg)
    fig5_summary_df = load_fig5_summary_input(cfg)
    fig5_expression_df = load_fig5_expression_input(cfg)

    run_qc_checks_fig123(fig1_df, fig2_df, fig4_evidence_df, paics_df, cfg)
    z_score_df, available_mask, evidence_class = build_fig4_heatmap_inputs(fig4_evidence_df, filtered_tpm_df, sample_meta_df)
    run_qc_checks_fig4(fig4_evidence_df, sample_meta_df, z_score_df, available_mask, cfg)
    contiguous_evidence_blocks(evidence_class)

    primary_gene = cfg.expected_primary_gene
    cfg.output_dir.mkdir(parents=True, exist_ok=True)

    fig1 = plot_fig1(fig1_df, paics_df.iloc[0], cfg, primary_gene)
    _save_figure(fig1, cfg.fig1_png, cfg.fig1_pdf)

    fig2 = plot_fig2(fig2_df, cfg)
    _save_figure(fig2, cfg.fig2_png, cfg.fig2_pdf)

    fig3 = plot_fig3(fig3_df, cfg, primary_gene)
    _save_figure(fig3, cfg.fig3_png, cfg.fig3_pdf)

    fig4, fig4_artifacts = plot_fig4(z_score_df, available_mask, evidence_class, sample_meta_df, fig5_summary_df, cfg, primary_gene)
    fig4_violations = check_fig4_layout(fig4, fig4_artifacts)
    if fig4_violations:
        raise ValueError(f"Figure 4 layout violations: {fig4_violations}")
    _save_figure(fig4, cfg.fig4_png, cfg.fig4_pdf)

    fig5, fig5_artifacts = plot_fig5(fig5_summary_df, fig5_expression_df, cfg, primary_gene)
    fig5_violations = check_fig5_layout(fig5, fig5_artifacts)
    if fig5_violations:
        raise ValueError(f"Figure 5 layout violations: {fig5_violations}")
    _save_figure(fig5, cfg.fig5_png, cfg.fig5_pdf)

    build_contact_sheet(cfg)
    manifest_df = build_manifest(cfg)

    logger.info("run_nebula_plots_final: wrote 5 final poster figures (PNG + PDF), the contact sheet, and the manifest")

    return {
        "fig1_input": fig1_df,
        "fig1_paics_inset": paics_df,
        "fig2_input": fig2_df,
        "fig3_input": fig3_df,
        "fig4_evidence_input": fig4_evidence_df,
        "z_score_df": z_score_df,
        "available_mask": available_mask,
        "evidence_class": evidence_class,
        "fig5_summary_input": fig5_summary_df,
        "fig5_expression_input": fig5_expression_df,
        "manifest": manifest_df,
    }


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run_nebula_plots_final()
