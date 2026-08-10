"""NEBULA poster figures for the candidate-prioritisation phase.

Plotting only. This module performs no CRISPR refit, no limma refit, no
PCA recomputation, and no candidate reclassification -- it reads exclusively
the following frozen, committed result tables and renders poster-ready
figures plus the exact plot-input TSVs used to draw them:

- ``results/tables/candidate_evidence_summary.tsv`` (28 Gate-1 hits, config
  ``candidate_evidence_summary.output.evidence_summary_tsv``), produced by
  ``src.candidate_evidence_summary.run_candidate_evidence_summary``.
- ``results/tables/candidate_shortlist.tsv``,
  ``candidate_sensitisation_candidates.tsv``, ``candidate_tolerance_hits.tsv``,
  ``candidate_paics_benchmark.tsv`` -- same module, same run.
- ``results/tables/gse118713_pca_coordinates.tsv`` (config
  ``gse118713_phase2b.qc.pca_coordinates_tsv``), produced by
  ``src.gse118713_qc.run_sample_qc``.
- ``results/tables/gse118713_differential_expression.tsv.gz`` (config
  ``gse118713_phase2b.limma.differential_expression_tsv_gz``), produced by
  ``scripts/analysis/gse118713_limma.R``.
- the checksum-pinned filtered GSE118713 TPM matrix (config
  ``gse118713_phase2b.filtering.filtered_gene_tpm_tsv``, verified against
  ``frozen_filtered_gene_tpm_sha256``) and
  ``results/tables/gse118713_sample_metadata.tsv``.

PAICS (``src.candidate_evidence_summary.BENCHMARK_GENE_SYMBOL``) is a
published benchmark, not one of the 28 Gate-1 hits. It is never merged into
the 28-hit landscape (Figure 1) or the 13-gene sensitising evidence matrix
(Figure 4); Figure 1 shows it only in a clearly separate, dashed-border
inset panel explicitly labelled "Published PAICS benchmark -- not a Gate-1
hit."

TAMR_vs_MCF7 is treated as the PRIMARY resistance-expression contrast
throughout (Figures 3-5); TAMR_vs_FASR is shown only as SECONDARY/contextual
evidence and is never described or annotated as tamoxifen-specific.

No composite/weighted score is computed anywhere in this module. Every
visual channel (bar length, dot size, dot color, marker fill) maps to one
existing, already-reviewed numeric column; multi-column evidence (Figure 4)
is shown as a transparent small-multiples matrix, not blended into a single
number.
"""

from __future__ import annotations

import hashlib
import logging
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.patches import Rectangle

from src.candidate_evidence_summary import (
    BENCHMARK_GENE_SYMBOL,
    BENCHMARK_LABEL,
    DIRECTION_SENSITISING,
    DIRECTION_TOLERANCE,
    EVIDENCE_CLASS_NO_SIGNIFICANT_RNA,
    EVIDENCE_CLASS_ORDER,
    EVIDENCE_CLASS_PRIMARY,
    EVIDENCE_CLASS_RNA_UNAVAILABLE,
    EVIDENCE_CLASS_SECONDARY,
)

logger = logging.getLogger(__name__)

PAICS_INSET_ANNOTATION = "Published PAICS benchmark — not a Gate-1 hit."
EVIDENCE_CLASS_SORT_ORDER: dict[str, int] = {cls: i for i, cls in enumerate(EVIDENCE_CLASS_ORDER)}

_DARK_TEXT = "#0B0B12"
_LIGHT_TEXT = "#F8FAFC"


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
class NebulaPlotsConfig:
    """Resolved, config-driven paths and style. No hardcoded paths."""

    evidence_summary_tsv: Path
    shortlist_tsv: Path
    sensitisation_candidates_tsv: Path
    tolerance_hits_tsv: Path
    paics_benchmark_tsv: Path
    pca_coordinates_tsv: Path
    differential_expression_tsv_gz: Path
    filtered_gene_tpm_tsv_gz: Path
    frozen_filtered_gene_tpm_sha256: str
    sample_metadata_tsv: Path
    rna_significance_fdr: float

    expected_n_hits: int
    expected_n_sensitising: int
    expected_n_tolerance: int
    expected_n_samples: int
    expected_primary_gene: str

    poster_background: str
    panel_background: str
    palette: dict[str, str]
    group_colors: dict[str, str]
    direction_colors: dict[str, str]

    output_dir: Path
    plot_input_dir: Path

    fig1_png: Path
    fig1_png_transparent: Path
    fig1_pdf: Path
    fig1_input_tsv: Path
    fig1_paics_inset_tsv: Path

    fig2_png: Path
    fig2_png_transparent: Path
    fig2_pdf: Path
    fig2_input_tsv: Path

    fig3_png: Path
    fig3_png_transparent: Path
    fig3_pdf: Path
    fig3_input_tsv: Path

    fig4_png: Path
    fig4_png_transparent: Path
    fig4_pdf: Path
    fig4_input_tsv: Path

    fig5_png: Path
    fig5_png_transparent: Path
    fig5_pdf: Path
    fig5_summary_input_tsv: Path
    fig5_expression_input_tsv: Path

    @classmethod
    def from_config(cls, config: dict) -> "NebulaPlotsConfig":
        ces_out = config["candidate_evidence_summary"]["output"]
        integ = config["crispr_gse118713_integration"]
        qc = config["gse118713_phase2b"]["qc"]
        limma = config["gse118713_phase2b"]["limma"]
        filtering = config["gse118713_phase2b"]["filtering"]
        gse = config["gse118713"]
        np_cfg = config["nebula_plots"]
        style = np_cfg["style"]
        figs = np_cfg["figures"]
        npqc = np_cfg["qc"]

        return cls(
            evidence_summary_tsv=Path(ces_out["evidence_summary_tsv"]),
            shortlist_tsv=Path(ces_out["shortlist_tsv"]),
            sensitisation_candidates_tsv=Path(ces_out["sensitisation_candidates_tsv"]),
            tolerance_hits_tsv=Path(ces_out["tolerance_hits_tsv"]),
            paics_benchmark_tsv=Path(ces_out["paics_benchmark_tsv"]),
            pca_coordinates_tsv=Path(qc["pca_coordinates_tsv"]),
            differential_expression_tsv_gz=Path(limma["differential_expression_tsv_gz"]),
            filtered_gene_tpm_tsv_gz=Path(filtering["filtered_gene_tpm_tsv"]),
            frozen_filtered_gene_tpm_sha256=str(filtering["frozen_filtered_gene_tpm_sha256"]),
            sample_metadata_tsv=Path(gse["output"]["sample_metadata_tsv"]),
            rna_significance_fdr=float(np_cfg["rna_significance_fdr"]),
            expected_n_hits=int(integ["expected_n_hits"]),
            expected_n_sensitising=int(npqc["expected_n_sensitising"]),
            expected_n_tolerance=int(npqc["expected_n_tolerance"]),
            expected_n_samples=int(npqc["expected_n_samples"]),
            expected_primary_gene=str(npqc["expected_primary_gene"]),
            poster_background=str(style["poster_background"]),
            panel_background=str(style["panel_background"]),
            palette=dict(style["palette"]),
            group_colors=dict(style["group_colors"]),
            direction_colors=dict(style["direction_colors"]),
            output_dir=Path(np_cfg["output_dir"]),
            plot_input_dir=Path(np_cfg["plot_input_dir"]),
            fig1_png=Path(figs["crispr_landscape"]["png"]),
            fig1_png_transparent=Path(figs["crispr_landscape"]["png_transparent"]),
            fig1_pdf=Path(figs["crispr_landscape"]["pdf"]),
            fig1_input_tsv=Path(figs["crispr_landscape"]["plot_input_tsv"]),
            fig1_paics_inset_tsv=Path(figs["crispr_landscape"]["paics_inset_tsv"]),
            fig2_png=Path(figs["pca"]["png"]),
            fig2_png_transparent=Path(figs["pca"]["png_transparent"]),
            fig2_pdf=Path(figs["pca"]["pdf"]),
            fig2_input_tsv=Path(figs["pca"]["plot_input_tsv"]),
            fig3_png=Path(figs["volcano"]["png"]),
            fig3_png_transparent=Path(figs["volcano"]["png_transparent"]),
            fig3_pdf=Path(figs["volcano"]["pdf"]),
            fig3_input_tsv=Path(figs["volcano"]["plot_input_tsv"]),
            fig4_png=Path(figs["evidence_matrix"]["png"]),
            fig4_png_transparent=Path(figs["evidence_matrix"]["png_transparent"]),
            fig4_pdf=Path(figs["evidence_matrix"]["pdf"]),
            fig4_input_tsv=Path(figs["evidence_matrix"]["plot_input_tsv"]),
            fig5_png=Path(figs["usp34_panel"]["png"]),
            fig5_png_transparent=Path(figs["usp34_panel"]["png_transparent"]),
            fig5_pdf=Path(figs["usp34_panel"]["pdf"]),
            fig5_summary_input_tsv=Path(figs["usp34_panel"]["summary_plot_input_tsv"]),
            fig5_expression_input_tsv=Path(figs["usp34_panel"]["expression_plot_input_tsv"]),
        )


# --------------------------------------------------------------------------
# Loading frozen sources (read-only; no recomputation)
# --------------------------------------------------------------------------


def load_evidence_summary(cfg: NebulaPlotsConfig) -> pd.DataFrame:
    """Load the frozen 28-gene evidence table. Raises if the row count is
    wrong or if PAICS is ever found among the 28 hits."""
    df = pd.read_csv(cfg.evidence_summary_tsv, sep="\t")
    if len(df) != cfg.expected_n_hits:
        raise ValueError(f"expected {cfg.expected_n_hits} Gate-1 hits in evidence summary, found {len(df)}")
    if BENCHMARK_GENE_SYMBOL in set(df["gene_symbol"]):
        raise ValueError(f"{BENCHMARK_GENE_SYMBOL} must never appear among the 28 Gate-1 hits")
    logger.info("load_evidence_summary: read %d rows", len(df))
    return df


def load_shortlist(cfg: NebulaPlotsConfig) -> pd.DataFrame:
    df = pd.read_csv(cfg.shortlist_tsv, sep="\t")
    if list(df["gene_symbol"]) != [cfg.expected_primary_gene]:
        raise ValueError(
            f"expected the primary shortlist to be exactly [{cfg.expected_primary_gene!r}], "
            f"found {list(df['gene_symbol'])!r}"
        )
    logger.info("load_shortlist: %s confirmed as sole primary candidate", cfg.expected_primary_gene)
    return df


def load_sensitising_candidates(cfg: NebulaPlotsConfig) -> pd.DataFrame:
    df = pd.read_csv(cfg.sensitisation_candidates_tsv, sep="\t")
    if len(df) != cfg.expected_n_sensitising:
        raise ValueError(f"expected {cfg.expected_n_sensitising} sensitising genes, found {len(df)}")
    logger.info("load_sensitising_candidates: read %d rows", len(df))
    return df


def load_tolerance_hits(cfg: NebulaPlotsConfig) -> pd.DataFrame:
    df = pd.read_csv(cfg.tolerance_hits_tsv, sep="\t")
    if len(df) != cfg.expected_n_tolerance:
        raise ValueError(f"expected {cfg.expected_n_tolerance} tolerance-associated genes, found {len(df)}")
    logger.info("load_tolerance_hits: read %d rows", len(df))
    return df


def load_paics_benchmark(cfg: NebulaPlotsConfig) -> pd.DataFrame:
    df = pd.read_csv(cfg.paics_benchmark_tsv, sep="\t")
    if len(df) != 1 or df.iloc[0]["gene_symbol"] != BENCHMARK_GENE_SYMBOL:
        raise ValueError("expected exactly one PAICS benchmark row")
    if df.iloc[0]["benchmark_label"] != BENCHMARK_LABEL:
        raise ValueError(f"PAICS benchmark row missing expected label {BENCHMARK_LABEL!r}")
    logger.info("load_paics_benchmark: confirmed PAICS is labelled %r, separate from Gate-1 hits", BENCHMARK_LABEL)
    return df


def load_pca_coordinates(cfg: NebulaPlotsConfig) -> pd.DataFrame:
    df = pd.read_csv(cfg.pca_coordinates_tsv, sep="\t")
    n_samples = df["sample_id"].nunique()
    if n_samples != cfg.expected_n_samples:
        raise ValueError(f"expected {cfg.expected_n_samples} PCA samples, found {n_samples}")
    # Exactly one PC1 and one PC2 coordinate per sample -- catches
    # duplicated sample IDs and missing components in the same check.
    for pc in ("PC1", "PC2"):
        pc_rows = df.loc[df["pc"] == pc]
        counts = pc_rows["sample_id"].value_counts()
        if len(counts) != n_samples or (counts != 1).any():
            raise ValueError(f"expected exactly one {pc} coordinate per sample ({n_samples} samples)")
        if not np.isfinite(pc_rows["coordinate"].to_numpy(dtype=float)).all():
            raise ValueError(f"PCA {pc} coordinates contain non-finite values")

    # A sample's PC1 row and PC2 row must agree on which group it belongs
    # to (e.g. PC1 saying TAMR while PC2 says MCF7 would silently corrupt
    # Figure 2's coloring/legend).
    pc1_groups = df.loc[df["pc"] == "PC1"].set_index("sample_id")["group"]
    pc2_groups = df.loc[df["pc"] == "PC2"].set_index("sample_id")["group"]
    mismatched = pc1_groups.index[pc1_groups != pc2_groups.reindex(pc1_groups.index)]
    if len(mismatched) > 0:
        raise ValueError(f"PCA group label mismatch between PC1 and PC2 rows for samples: {list(mismatched)}")

    logger.info("load_pca_coordinates: read %d samples", n_samples)
    return df


def load_mcf7_volcano_source(cfg: NebulaPlotsConfig) -> pd.DataFrame:
    """Full TAMR_vs_MCF7 differential-expression background (the primary
    resistance-expression contrast), unfiltered by gene identity."""
    df = pd.read_csv(cfg.differential_expression_tsv_gz, sep="\t")
    mcf7 = df.loc[df["contrast"] == "TAMR_vs_MCF7"].reset_index(drop=True)
    if mcf7.empty:
        raise ValueError("no TAMR_vs_MCF7 rows found in the frozen differential-expression table")
    logger.info("load_mcf7_volcano_source: read %d TAMR_vs_MCF7 background genes", len(mcf7))
    return mcf7


def load_filtered_tpm(cfg: NebulaPlotsConfig) -> pd.DataFrame:
    actual_sha256 = _sha256_file(cfg.filtered_gene_tpm_tsv_gz)
    if actual_sha256 != cfg.frozen_filtered_gene_tpm_sha256:
        raise ValueError(
            f"filtered GSE118713 TPM matrix checksum mismatch: expected "
            f"{cfg.frozen_filtered_gene_tpm_sha256}, got {actual_sha256}"
        )
    df = pd.read_csv(cfg.filtered_gene_tpm_tsv_gz, sep="\t")
    logger.info("load_filtered_tpm: read %d genes (checksum verified)", len(df))
    return df


def load_sample_metadata(cfg: NebulaPlotsConfig) -> pd.DataFrame:
    df = pd.read_csv(cfg.sample_metadata_tsv, sep="\t")
    if df["sample_id"].nunique() != cfg.expected_n_samples:
        raise ValueError(f"expected {cfg.expected_n_samples} samples in sample metadata, found {df['sample_id'].nunique()}")
    logger.info("load_sample_metadata: read %d samples", len(df))
    return df


# --------------------------------------------------------------------------
# Shared style helpers
# --------------------------------------------------------------------------


def _new_axes(cfg: NebulaPlotsConfig, figsize: tuple[float, float]) -> tuple[plt.Figure, plt.Axes]:
    fig, ax = plt.subplots(figsize=figsize)
    fig.patch.set_facecolor(cfg.poster_background)
    ax.set_facecolor(cfg.panel_background)
    for spine in ax.spines.values():
        spine.set_color(cfg.palette["neutral_grey"])
        spine.set_linewidth(0.8)
    # Ticks and axis labels sit just outside the panel, over the dark poster
    # background -- they need light text, not dark, to stay readable. Only
    # text drawn *inside* the white panel area should use dark text.
    ax.tick_params(colors=_LIGHT_TEXT, labelsize=8)
    ax.xaxis.label.set_color(_LIGHT_TEXT)
    ax.yaxis.label.set_color(_LIGHT_TEXT)
    return fig, ax


def _style_title(ax: plt.Axes, title: str, subtitle: str | None = None) -> None:
    ax.set_title(title, color=_LIGHT_TEXT, fontsize=12, fontweight="bold", pad=12)
    if subtitle:
        ax.text(
            0.5,
            1.06,
            subtitle,
            transform=ax.transAxes,
            ha="center",
            va="bottom",
            fontsize=8,
            color="#B799FF",
        )


def _diverging_cmap(negative_hex: str, positive_hex: str) -> LinearSegmentedColormap:
    return LinearSegmentedColormap.from_list("nebula_diverging", [negative_hex, "#FFFFFF", positive_hex])


def _save_figure(fig: plt.Figure, png_path: Path, png_transparent_path: Path, pdf_path: Path) -> None:
    for path in (png_path, png_transparent_path, pdf_path):
        path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(png_path, dpi=300, facecolor=fig.get_facecolor(), bbox_inches="tight")
    fig.savefig(png_transparent_path, dpi=300, transparent=True, bbox_inches="tight")
    fig.savefig(pdf_path, metadata={"CreationDate": None}, bbox_inches="tight")
    plt.close(fig)
    logger.info("_save_figure: wrote %s, %s, %s", png_path, png_transparent_path, pdf_path)


def _place_labels(
    ax: plt.Axes,
    points: list[tuple[float, float, str]],
    color_fn,
    fontsize: float = 7.5,
    base_offset: float = 12.0,
) -> None:
    """Minimal greedy label placer: spiral each label out along an
    increasing angle/radius until it clears previously placed labels in
    display space (clearance scales with the label's approximate rendered
    width), so a small set of selective labels stays readable without a
    third-party layout dependency."""
    placed_boxes: list[tuple[np.ndarray, np.ndarray]] = []  # (center, half_extent)
    for x, y, label in points:
        text = label.rstrip("*")
        is_primary = label.endswith("*")
        size = fontsize + (1.5 if is_primary else 0)
        half_extent = np.array([max(len(text), 1) * size * 0.34, size * 0.9])

        disp = np.array(ax.transData.transform((x, y)))
        offset = np.array([base_offset, base_offset])
        tries = 0
        max_tries = 24
        while tries < max_tries:
            candidate_center = disp + offset
            collides = any(
                abs(candidate_center[0] - c[0]) < (half_extent[0] + h[0])
                and abs(candidate_center[1] - c[1]) < (half_extent[1] + h[1])
                for c, h in placed_boxes
            )
            if not collides:
                break
            angle = tries * (np.pi / 4)
            radius = base_offset * (1.4 + 0.55 * tries)
            offset = radius * np.array([np.cos(angle), np.sin(angle)]) + np.array([base_offset, base_offset])
            tries += 1
        placed_boxes.append((disp + offset, half_extent))
        annotation = ax.annotate(
            text,
            xy=(x, y),
            xytext=offset,
            textcoords="offset points",
            fontsize=size,
            fontweight="bold" if is_primary else "normal",
            color=color_fn(text),
            arrowprops=dict(arrowstyle="-", color="#94A3B8", lw=0.6, alpha=0.8),
            annotation_clip=False,
        )
        annotation.set_clip_on(False)


# --------------------------------------------------------------------------
# Figure 1: CRISPR hit landscape (28 Gate-1 hits, PAICS in a separate inset)
# --------------------------------------------------------------------------


def build_fig1_inputs(evidence_df: pd.DataFrame, paics_df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Plot-input table for the 28-hit landscape (sorted by the existing
    ``crispr_effect_size`` column, ascending -- no new rank/score is
    introduced), plus a strictly separate one-row PAICS inset table."""
    fig1_df = evidence_df[["gene_symbol", "crispr_effect_size", "crispr_fdr", "crispr_direction"]].copy()
    fig1_df = fig1_df.sort_values(["crispr_effect_size", "gene_symbol"]).reset_index(drop=True)

    paics_inset_df = paics_df[
        ["gene_symbol", "crispr_effect_size", "crispr_fdr", "crispr_direction", "benchmark_label"]
    ].copy()
    return fig1_df, paics_inset_df


def plot_fig1_crispr_landscape(
    fig1_df: pd.DataFrame, paics_inset_df: pd.DataFrame, cfg: NebulaPlotsConfig, primary_gene: str
) -> plt.Figure:
    fig = plt.figure(figsize=(9, 10))
    fig.patch.set_facecolor(cfg.poster_background)
    gs = fig.add_gridspec(1, 4, width_ratios=[3, 3, 0.15, 1.1], wspace=0.05)
    ax = fig.add_subplot(gs[0, 0:2])
    ax_inset = fig.add_subplot(gs[0, 3])

    ax.set_facecolor(cfg.panel_background)
    for spine in ax.spines.values():
        spine.set_color(cfg.palette["neutral_grey"])

    y = np.arange(len(fig1_df))
    colors = [cfg.direction_colors[d] for d in fig1_df["crispr_direction"]]
    ax.hlines(y, 0, fig1_df["crispr_effect_size"], color=colors, linewidth=3, zorder=2)

    sizes = 25 + 220 * np.clip(-np.log10(fig1_df["crispr_fdr"].to_numpy(dtype=float)) / 7.0, 0, 1)
    ax.scatter(fig1_df["crispr_effect_size"], y, s=sizes, color=colors, edgecolor="white", linewidth=0.6, zorder=3)

    ax.axvline(0, color=_DARK_TEXT, linewidth=1.0, zorder=1)
    ax.set_yticks(y)
    # Gene-name labels sit outside the panel, over the dark poster
    # background -- light text by default, except USP34 which gets its own
    # bright accent color for emphasis (still readable on dark).
    labels = ax.set_yticklabels(fig1_df["gene_symbol"], fontsize=7.5)
    for tick_label, gene in zip(labels, fig1_df["gene_symbol"]):
        if gene == primary_gene:
            tick_label.set_fontweight("bold")
            tick_label.set_color(cfg.direction_colors[DIRECTION_SENSITISING])
            tick_label.set_fontsize(9.5)
        else:
            tick_label.set_color(_LIGHT_TEXT)
    ax.tick_params(axis="x", colors=_LIGHT_TEXT, labelsize=8)
    ax.tick_params(axis="y", length=0)
    ax.set_xlabel("CRISPR effect_size (E2+4-OHT minus E2)", color=_LIGHT_TEXT, fontsize=9)
    ax.set_ylim(-1, len(fig1_df) + 1.6)
    ax.invert_yaxis()
    _style_title(ax, "Gate-1 CRISPR hit landscape (n=28)", "Ranked by CRISPR effect_size — no composite score")

    legend_handles = [
        plt.Line2D([0], [0], marker="o", color="none", markerfacecolor=cfg.direction_colors[DIRECTION_SENSITISING], markersize=9, label="sensitising_knockout (negative effect)"),
        plt.Line2D([0], [0], marker="o", color="none", markerfacecolor=cfg.direction_colors[DIRECTION_TOLERANCE], markersize=9, label="tolerance_associated_knockout (positive effect)"),
    ]
    leg = ax.legend(
        handles=legend_handles,
        loc="lower center",
        bbox_to_anchor=(0.5, -0.14),
        ncol=1,
        fontsize=7.5,
        frameon=True,
        facecolor="white",
        edgecolor=cfg.palette["neutral_grey"],
    )
    leg.set_clip_on(False)  # placed outside the axes via bbox_to_anchor -- must not be clipped to the axes box
    for text in leg.get_texts():
        text.set_color(_DARK_TEXT)

    ref_fdrs = [0.1, 0.01, 0.001]
    ref_sizes = [25 + 220 * np.clip(-np.log10(f) / 7.0, 0, 1) for f in ref_fdrs]
    size_handles = [
        plt.Line2D([0], [0], marker="o", color="none", markerfacecolor=cfg.palette["neutral_grey"], markeredgecolor="white", markersize=np.sqrt(s) / 1.6, label=f"FDR={f:g}")
        for f, s in zip(ref_fdrs, ref_sizes)
    ]
    size_leg = ax.legend(handles=size_handles, loc="upper right", fontsize=6.5, title="Dot size ∝ -log10(CRISPR FDR)", title_fontsize=6.5, frameon=True, facecolor="white", edgecolor=cfg.palette["neutral_grey"])
    for text in size_leg.get_texts():
        text.set_color(_DARK_TEXT)
    size_leg.get_title().set_color(_DARK_TEXT)
    ax.add_artist(leg)

    ax_inset.set_facecolor(cfg.panel_background)
    for spine in ax_inset.spines.values():
        spine.set_color(cfg.palette["cosmic_magenta"])
        spine.set_linestyle((0, (4, 3)))
        spine.set_linewidth(1.4)
    ax_inset.set_xticks([])
    ax_inset.set_yticks([])
    row = paics_inset_df.iloc[0]
    ax_inset.scatter([0.5], [0.55], s=140, color=cfg.palette["cosmic_magenta"], edgecolor="white", linewidth=0.8, zorder=3)
    ax_inset.text(0.5, 0.42, row["gene_symbol"], ha="center", fontsize=10, fontweight="bold", color=_DARK_TEXT)
    ax_inset.text(
        0.5,
        0.30,
        f"effect_size={row['crispr_effect_size']:.3f}\nFDR={row['crispr_fdr']:.3f}",
        ha="center",
        fontsize=7.5,
        color=_DARK_TEXT,
    )
    ax_inset.text(
        0.5,
        0.90,
        PAICS_INSET_ANNOTATION,
        ha="center",
        va="top",
        fontsize=7,
        fontweight="bold",
        wrap=True,
        color=cfg.palette["cosmic_magenta"],
        transform=ax_inset.transAxes,
    )
    ax_inset.set_xlim(0, 1)
    ax_inset.set_ylim(0, 1)

    return fig


# --------------------------------------------------------------------------
# Figure 2: GSE118713 PCA
# --------------------------------------------------------------------------


def build_fig2_inputs(pca_df: pd.DataFrame) -> pd.DataFrame:
    wide = pca_df.pivot(index=["sample_id", "group"], columns="pc", values="coordinate").reset_index()
    variance = pca_df.drop_duplicates("pc").set_index("pc")["variance_explained_fraction"]
    wide = wide[["sample_id", "group", "PC1", "PC2"]].copy()
    wide["pc1_variance_explained_pct"] = variance["PC1"] * 100
    wide["pc2_variance_explained_pct"] = variance["PC2"] * 100
    return wide.sort_values(["group", "sample_id"]).reset_index(drop=True)


def plot_fig2_pca(fig2_df: pd.DataFrame, cfg: NebulaPlotsConfig) -> plt.Figure:
    fig, ax = _new_axes(cfg, (6.5, 5.5))
    for group, sub in fig2_df.groupby("group", sort=False):
        ax.scatter(
            sub["PC1"],
            sub["PC2"],
            label=group,
            color=cfg.group_colors[group],
            s=110,
            edgecolor="white",
            linewidth=0.8,
            zorder=3,
        )
    ax.margins(0.22)  # leave room for sample-id labels near the panel edges
    label_points = [(float(r["PC1"]), float(r["PC2"]), r["sample_id"]) for _, r in fig2_df.iterrows()]
    _place_labels(ax, label_points, lambda _label: _DARK_TEXT, fontsize=7.0, base_offset=10.0)
    ax.axhline(0, color=cfg.palette["neutral_grey"], linewidth=0.6, zorder=1)
    ax.axvline(0, color=cfg.palette["neutral_grey"], linewidth=0.6, zorder=1)
    ax.set_xlabel(f"PC1 ({fig2_df['pc1_variance_explained_pct'].iloc[0]:.1f}%)")
    ax.set_ylabel(f"PC2 ({fig2_df['pc2_variance_explained_pct'].iloc[0]:.1f}%)")
    _style_title(ax, "GSE118713: sample PCA (n=9)", "MCF7 / TAMR / FASR, n=3 replicates per group — no ellipse fit at this n")
    leg = ax.legend(loc="best", fontsize=8, frameon=True, facecolor="white", edgecolor=cfg.palette["neutral_grey"])
    for text in leg.get_texts():
        text.set_color(_DARK_TEXT)
    fig.tight_layout()
    return fig


# --------------------------------------------------------------------------
# Figure 3: TAMR_vs_MCF7 volcano (primary resistance contrast)
# --------------------------------------------------------------------------


def build_fig3_inputs(volcano_background_df: pd.DataFrame, evidence_df: pd.DataFrame) -> pd.DataFrame:
    """Background genes (all TAMR_vs_MCF7-tested genes) annotated with
    CRISPR-hit membership/direction where applicable -- a left join, no
    filtering of the background cloud.

    FDR is validated here, before any ``-log10`` transform is taken for the
    volcano y-axis: a zero, negative, out-of-range, or non-finite FDR would
    otherwise silently produce an infinite or nonsensical -log10(FDR) value.
    """
    hits = evidence_df[["gene_symbol", "crispr_direction", "evidence_class"]].copy()
    merged = volcano_background_df.merge(hits, on="gene_symbol", how="left")
    merged["is_gate1_hit"] = merged["crispr_direction"].notna()

    fdr_values = merged["fdr"].to_numpy(dtype=float)
    if not np.isfinite(fdr_values).all():
        raise ValueError("TAMR_vs_MCF7 volcano input contains non-finite FDR values")
    if not ((fdr_values > 0) & (fdr_values <= 1)).all():
        raise ValueError("TAMR_vs_MCF7 volcano input contains FDR values outside (0, 1]")

    return merged


def plot_fig3_volcano(
    fig3_df: pd.DataFrame,
    cfg: NebulaPlotsConfig,
    primary_gene: str,
    secondary_genes: list[str],
) -> plt.Figure:
    fig, ax = _new_axes(cfg, (8, 6.5))

    background = fig3_df.loc[~fig3_df["is_gate1_hit"]]
    ax.scatter(
        background["log2fc"],
        -np.log10(background["fdr"]),
        s=6,
        color=cfg.palette["neutral_grey"],
        alpha=0.35,
        linewidth=0,
        zorder=1,
        label="Other tested genes",
    )

    hits = fig3_df.loc[fig3_df["is_gate1_hit"]]
    for direction, sub in hits.groupby("crispr_direction"):
        ax.scatter(
            sub["log2fc"],
            -np.log10(sub["fdr"]),
            s=42,
            color=cfg.direction_colors[direction],
            edgecolor="white",
            linewidth=0.5,
            zorder=2,
            label=direction,
        )

    threshold_y = -np.log10(cfg.rna_significance_fdr)
    ax.axhline(threshold_y, color=cfg.palette["cosmic_magenta"], linestyle="--", linewidth=1.0, zorder=1)
    ax.text(
        ax.get_xlim()[1] if ax.get_xlim()[1] else 1,
        threshold_y,
        f" FDR = {cfg.rna_significance_fdr:g}",
        color=cfg.palette["cosmic_magenta"],
        fontsize=7,
        va="bottom",
        ha="right",
    )

    # Give downward-nudged labels room to stay inside the white panel
    # instead of spilling onto the dark poster background below y=0.
    ax.set_ylim(bottom=-0.8)

    # Secondary-context genes get a "†" marker in their label so they
    # cannot be misread as significant in *this* (TAMR_vs_MCF7) contrast --
    # their evidence-class significance comes from TAMR_vs_FASR instead
    # (see the caption added below).
    label_points: list[tuple[float, float, str]] = []
    usp34_row = fig3_df.loc[fig3_df["gene_symbol"] == primary_gene]
    if not usp34_row.empty:
        r = usp34_row.iloc[0]
        label_points.append((float(r["log2fc"]), float(-np.log10(r["fdr"])), f"{primary_gene}*"))
    for gene in secondary_genes:
        row = fig3_df.loc[fig3_df["gene_symbol"] == gene]
        if not row.empty:
            r = row.iloc[0]
            label_points.append((float(r["log2fc"]), float(-np.log10(r["fdr"])), f"{gene}†"))

    def _color_for(gene: str) -> str:
        if gene == primary_gene:
            return cfg.direction_colors[DIRECTION_SENSITISING]
        return _DARK_TEXT

    _place_labels(ax, label_points, _color_for)

    if secondary_genes:
        ax.text(
            0.02,
            0.02,
            "† Secondary-context CRISPR candidates: labelled from their TAMR_vs_FASR evidence,\n"
            "not necessarily significant in this TAMR_vs_MCF7 volcano.",
            transform=ax.transAxes,
            fontsize=6.5,
            color=_DARK_TEXT,
            va="bottom",
            ha="left",
            style="italic",
        )

    ax.set_xlabel("log2FC (TAMR vs MCF7)")
    ax.set_ylabel("-log10(FDR)")
    _style_title(ax, "TAMR_vs_MCF7 volcano — PRIMARY resistance-expression contrast", "TAMR_vs_FASR is shown separately as secondary context only (Figure 4)")

    handles, labels_ = ax.get_legend_handles_labels()
    seen = set()
    dedup_handles, dedup_labels = [], []
    for h, l in zip(handles, labels_):
        if l not in seen:
            dedup_handles.append(h)
            dedup_labels.append(l)
            seen.add(l)
    leg = ax.legend(dedup_handles, dedup_labels, loc="upper left", fontsize=7.5, frameon=True, facecolor="white", edgecolor=cfg.palette["neutral_grey"])
    for text in leg.get_texts():
        text.set_color(_DARK_TEXT)
    fig.tight_layout()
    return fig


# --------------------------------------------------------------------------
# Figure 4: integrated sensitising-candidate evidence matrix (13 genes)
# --------------------------------------------------------------------------


def build_fig4_inputs(sensitising_df: pd.DataFrame) -> pd.DataFrame:
    """Deterministic order: evidence class (PRIMARY > SECONDARY >
    NO_SIGNIFICANT_RNA_SUPPORT > RNA_UNAVAILABLE), then gene symbol
    alphabetically within class -- no hidden score."""
    df = sensitising_df.copy()
    df["_class_sort"] = df["evidence_class"].map(EVIDENCE_CLASS_SORT_ORDER)
    df = df.sort_values(["_class_sort", "gene_symbol"]).drop(columns="_class_sort").reset_index(drop=True)
    return df


def plot_fig4_evidence_matrix(fig4_df: pd.DataFrame, cfg: NebulaPlotsConfig, primary_gene: str) -> plt.Figure:
    n = len(fig4_df)
    fig, axes = plt.subplots(1, 4, figsize=(13, 0.55 * n + 2.2), gridspec_kw={"width_ratios": [1.1, 1.4, 1.4, 1.0]})
    fig.patch.set_facecolor(cfg.poster_background)

    y = np.arange(n)
    cmap_mcf7 = _diverging_cmap(cfg.palette["soft_lilac"], cfg.palette["rose_pink"])
    cmap_fasr = _diverging_cmap(cfg.palette["soft_lilac"], cfg.palette["rose_pink"])

    # Column A: CRISPR effect_size + FDR
    ax = axes[0]
    ax.set_facecolor(cfg.panel_background)
    vmax = fig4_df["crispr_effect_size"].abs().max()
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)
    cmap_crispr = _diverging_cmap(cfg.palette["nebula_purple"], cfg.palette["rose_pink"])
    for i, (_, r) in enumerate(fig4_df.iterrows()):
        ax.add_patch(Rectangle((0, i - 0.45), 1, 0.9, color=cmap_crispr(norm(r["crispr_effect_size"])), zorder=1))
        ax.text(0.5, i, f"{r['crispr_effect_size']:.2f}", ha="center", va="center", fontsize=7.5, color=_DARK_TEXT, zorder=2)
        ax.text(0.5, i + 0.32, f"FDR={r['crispr_fdr']:.2g}", ha="center", va="center", fontsize=6, color=_DARK_TEXT, zorder=2)
    ax.set_xlim(0, 1)
    ax.set_ylim(-0.6, n - 0.4)
    ax.set_xticks([])
    ax.set_yticks(y)
    # Gene-name labels sit outside the panel, over the dark poster
    # background -- light text by default, USP34 gets the accent color.
    yticklabels = ax.set_yticklabels(fig4_df["gene_symbol"], fontsize=8)
    for label, gene in zip(yticklabels, fig4_df["gene_symbol"]):
        if gene == primary_gene:
            label.set_fontweight("bold")
            label.set_fontsize(9.5)
            label.set_color(cfg.direction_colors[DIRECTION_SENSITISING])
        else:
            label.set_color(_LIGHT_TEXT)
    ax.invert_yaxis()
    ax.set_title("CRISPR effect_size\n(+ FDR)", fontsize=8, color=_LIGHT_TEXT)
    for spine in ax.spines.values():
        spine.set_visible(False)

    # Column B: TAMR_vs_MCF7 (PRIMARY)
    ax = axes[1]
    ax.set_facecolor(cfg.panel_background)
    vmax_b = fig4_df["tamr_vs_mcf7_log2fc"].abs().max()
    norm_b = TwoSlopeNorm(vmin=-vmax_b, vcenter=0, vmax=vmax_b) if vmax_b > 0 else TwoSlopeNorm(vmin=-1, vcenter=0, vmax=1)
    for i, (_, r) in enumerate(fig4_df.iterrows()):
        if pd.isna(r["tamr_vs_mcf7_log2fc"]):
            ax.add_patch(Rectangle((0, i - 0.45), 1, 0.9, color="#D1D5DB", zorder=1))
            ax.text(0.5, i, "N/A", ha="center", va="center", fontsize=7.5, color=_DARK_TEXT, zorder=2)
            continue
        ax.add_patch(Rectangle((0, i - 0.45), 1, 0.9, color=cmap_mcf7(norm_b(r["tamr_vs_mcf7_log2fc"])), zorder=1))
        sig = "*" if r["tamr_vs_mcf7_fdr"] < cfg.rna_significance_fdr else ""
        ax.text(0.5, i, f"{r['tamr_vs_mcf7_log2fc']:.2f}{sig}", ha="center", va="center", fontsize=7.5, color=_DARK_TEXT, fontweight="bold" if sig else "normal", zorder=2)
        ax.text(0.5, i + 0.32, f"FDR={r['tamr_vs_mcf7_fdr']:.2g}", ha="center", va="center", fontsize=6, color=_DARK_TEXT, zorder=2)
    ax.set_xlim(0, 1)
    ax.set_ylim(-0.6, n - 0.4)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.invert_yaxis()
    ax.set_title("TAMR_vs_MCF7 log2FC\nPRIMARY evidence (* FDR<%.2f)" % cfg.rna_significance_fdr, fontsize=8, color=_LIGHT_TEXT)
    for spine in ax.spines.values():
        spine.set_visible(False)

    # Column C: TAMR_vs_FASR (SECONDARY/contextual)
    ax = axes[2]
    ax.set_facecolor(cfg.panel_background)
    vmax_c = fig4_df["tamr_vs_fasr_log2fc"].abs().max()
    norm_c = TwoSlopeNorm(vmin=-vmax_c, vcenter=0, vmax=vmax_c) if vmax_c > 0 else TwoSlopeNorm(vmin=-1, vcenter=0, vmax=1)
    for i, (_, r) in enumerate(fig4_df.iterrows()):
        if pd.isna(r["tamr_vs_fasr_log2fc"]):
            ax.add_patch(Rectangle((0, i - 0.45), 1, 0.9, color="#D1D5DB", zorder=1))
            ax.text(0.5, i, "N/A", ha="center", va="center", fontsize=7.5, color=_DARK_TEXT, zorder=2)
            continue
        face = cmap_fasr(norm_c(r["tamr_vs_fasr_log2fc"]))
        face = (face[0], face[1], face[2], face[3] * 0.65)
        ax.add_patch(Rectangle((0, i - 0.45), 1, 0.9, facecolor=face, zorder=1))
        sig = "*" if r["tamr_vs_fasr_fdr"] < cfg.rna_significance_fdr else ""
        ax.text(0.5, i, f"{r['tamr_vs_fasr_log2fc']:.2f}{sig}", ha="center", va="center", fontsize=7.5, color=_DARK_TEXT, zorder=2)
        ax.text(0.5, i + 0.32, f"FDR={r['tamr_vs_fasr_fdr']:.2g}", ha="center", va="center", fontsize=6, color=_DARK_TEXT, zorder=2)
    ax.set_xlim(0, 1)
    ax.set_ylim(-0.6, n - 0.4)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.invert_yaxis()
    ax.set_title("TAMR_vs_FASR log2FC\nSECONDARY/context only (* FDR<%.2f)" % cfg.rna_significance_fdr, fontsize=8, color=_LIGHT_TEXT)
    for spine in ax.spines.values():
        spine.set_visible(False)

    # Column D: evidence class
    ax = axes[3]
    ax.set_facecolor(cfg.panel_background)
    class_colors = {
        EVIDENCE_CLASS_PRIMARY: cfg.direction_colors[DIRECTION_SENSITISING],
        EVIDENCE_CLASS_SECONDARY: cfg.palette["soft_lilac"],
        EVIDENCE_CLASS_NO_SIGNIFICANT_RNA: cfg.palette["neutral_grey"],
        EVIDENCE_CLASS_RNA_UNAVAILABLE: "#D1D5DB",
    }
    class_short = {
        EVIDENCE_CLASS_PRIMARY: "PRIMARY",
        EVIDENCE_CLASS_SECONDARY: "SECONDARY",
        EVIDENCE_CLASS_NO_SIGNIFICANT_RNA: "NO SIG.\nRNA",
        EVIDENCE_CLASS_RNA_UNAVAILABLE: "RNA\nUNAVAIL.",
    }
    for i, (_, r) in enumerate(fig4_df.iterrows()):
        color = class_colors[r["evidence_class"]]
        ax.add_patch(Rectangle((0, i - 0.45), 1, 0.9, color=color, zorder=1))
        text_color = "white" if r["evidence_class"] in (EVIDENCE_CLASS_PRIMARY,) else _DARK_TEXT
        ax.text(0.5, i, class_short[r["evidence_class"]], ha="center", va="center", fontsize=6.5, color=text_color, fontweight="bold", zorder=2)
    ax.set_xlim(0, 1)
    ax.set_ylim(-0.6, n - 0.4)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.invert_yaxis()
    ax.set_title("Evidence class", fontsize=8, color=_LIGHT_TEXT)
    for spine in ax.spines.values():
        spine.set_visible(False)

    if primary_gene in list(fig4_df["gene_symbol"]):
        idx = list(fig4_df["gene_symbol"]).index(primary_gene)
        for ax in axes:
            ax.add_patch(
                Rectangle(
                    (0, idx - 0.48),
                    1,
                    0.96,
                    fill=False,
                    edgecolor=cfg.direction_colors[DIRECTION_SENSITISING],
                    linewidth=2.2,
                    zorder=4,
                )
            )

    fig.suptitle(
        "Sensitising CRISPR hits (n=13): CRISPR + resistance-expression evidence, by class",
        color=_LIGHT_TEXT,
        fontsize=12,
        fontweight="bold",
        y=1.02,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.99))
    return fig


# --------------------------------------------------------------------------
# Figure 5: USP34 focused panel
# --------------------------------------------------------------------------


def build_fig5_summary_input(evidence_df: pd.DataFrame, primary_gene: str) -> pd.DataFrame:
    row = evidence_df.loc[evidence_df["gene_symbol"] == primary_gene]
    if len(row) != 1:
        raise ValueError(f"expected exactly one {primary_gene} row in the evidence summary, found {len(row)}")
    return row[
        [
            "gene_symbol",
            "crispr_effect_size",
            "crispr_fdr",
            "tamr_vs_mcf7_log2fc",
            "tamr_vs_mcf7_fdr",
            "tamr_vs_fasr_log2fc",
            "tamr_vs_fasr_fdr",
            "mcf7_baseline_log2_tpm_plus1",
        ]
    ].reset_index(drop=True)


def build_fig5_expression_input(filtered_tpm_df: pd.DataFrame, sample_meta_df: pd.DataFrame, primary_gene: str) -> pd.DataFrame:
    gene_row = filtered_tpm_df.loc[filtered_tpm_df["gene_symbol"] == primary_gene]
    if len(gene_row) != 1:
        raise ValueError(f"expected exactly one {primary_gene} row in the filtered TPM matrix, found {len(gene_row)}")
    gene_row = gene_row.iloc[0]

    records = []
    for _, meta_row in sample_meta_df.iterrows():
        sample_id = meta_row["sample_id"]
        if sample_id not in gene_row.index:
            raise ValueError(f"sample {sample_id!r} not found in the filtered TPM matrix columns")
        tpm = float(gene_row[sample_id])
        if not np.isfinite(tpm):
            raise ValueError(f"{primary_gene} TPM value for sample {sample_id!r} is not finite: {tpm!r}")
        if tpm < 0:
            raise ValueError(f"{primary_gene} TPM value for sample {sample_id!r} is negative: {tpm!r}")
        records.append(
            {
                "sample_id": sample_id,
                "group": meta_row["group"],
                "replicate": meta_row["replicate"],
                "tpm": tpm,
                "log2_tpm_plus1": float(np.log2(tpm + 1.0)),
            }
        )
    return pd.DataFrame(records).sort_values(["group", "replicate"]).reset_index(drop=True)


def plot_fig5_usp34_panel(
    summary_df: pd.DataFrame, expression_df: pd.DataFrame, cfg: NebulaPlotsConfig, primary_gene: str
) -> plt.Figure:
    fig = plt.figure(figsize=(9, 6))
    fig.patch.set_facecolor(cfg.poster_background)
    gs = fig.add_gridspec(1, 2, width_ratios=[1.3, 1])
    ax_card = fig.add_subplot(gs[0, 0])
    ax_expr = fig.add_subplot(gs[0, 1])

    row = summary_df.iloc[0]

    ax_card.set_facecolor(cfg.panel_background)
    ax_card.set_xticks([])
    ax_card.set_yticks([])
    ax_card.set_xlim(0, 1)
    ax_card.set_ylim(0, 1)
    for spine in ax_card.spines.values():
        spine.set_color(cfg.direction_colors[DIRECTION_SENSITISING])
        spine.set_linewidth(1.6)

    ax_card.text(0.5, 0.94, primary_gene, ha="center", va="top", fontsize=22, fontweight="bold", color=cfg.direction_colors[DIRECTION_SENSITISING])
    ax_card.text(0.5, 0.82, "Primary sensitisation candidate", ha="center", va="top", fontsize=9, color=_DARK_TEXT)

    lines = [
        f"CRISPR effect_size = {row['crispr_effect_size']:.3f}",
        f"CRISPR FDR = {row['crispr_fdr']:.4f}",
        "",
        "TAMR_vs_MCF7 (PRIMARY evidence)",
        f"  log2FC = {row['tamr_vs_mcf7_log2fc']:+.3f}   FDR = {row['tamr_vs_mcf7_fdr']:.5f}",
        "",
        "TAMR_vs_FASR (secondary context)",
        f"  log2FC = {row['tamr_vs_fasr_log2fc']:+.3f}   FDR = {row['tamr_vs_fasr_fdr']:.3f}",
        "",
        f"MCF7 baseline log2(TPM+1) = {row['mcf7_baseline_log2_tpm_plus1']:.3f}",
    ]
    ax_card.text(0.08, 0.70, "\n".join(lines), ha="left", va="top", fontsize=9, color=_DARK_TEXT, family="monospace")

    caption = (
        "USP34 combines a sensitising CRISPR knockout effect with\n"
        "significantly elevated expression in tamoxifen-resistant versus\n"
        "parental MCF7 cells."
    )
    ax_card.text(0.5, 0.20, caption, ha="center", va="top", fontsize=8.5, color=_DARK_TEXT, fontweight="bold")
    caveat = "Association and functional-screen evidence; not proof of a\nresistance mechanism or therapeutic efficacy."
    ax_card.text(0.5, 0.06, caveat, ha="center", va="top", fontsize=7, color=cfg.palette["neutral_grey"], style="italic")

    ax_expr.set_facecolor(cfg.panel_background)
    groups = ["MCF7", "TAMR", "FASR"]
    for gi, group in enumerate(groups):
        sub = expression_df.loc[expression_df["group"] == group]
        jitter = np.linspace(-0.12, 0.12, len(sub))
        ax_expr.scatter(
            gi + jitter,
            sub["log2_tpm_plus1"],
            color=cfg.group_colors[group],
            s=70,
            edgecolor="white",
            linewidth=0.7,
            zorder=3,
        )
        mean_val = sub["log2_tpm_plus1"].mean()
        ax_expr.hlines(mean_val, gi - 0.2, gi + 0.2, color=cfg.group_colors[group], linewidth=2.2, zorder=2)
    ax_expr.set_xticks(range(len(groups)))
    ax_expr.set_xticklabels(groups)
    ax_expr.set_ylabel("log2(TPM+1)")
    for spine in ax_expr.spines.values():
        spine.set_color(cfg.palette["neutral_grey"])
    ax_expr.tick_params(colors=_LIGHT_TEXT, labelsize=8)
    ax_expr.yaxis.label.set_color(_LIGHT_TEXT)
    _style_title(ax_expr, f"{primary_gene} expression", "actual GSE118713 replicate values (n=3/group)")

    fig.tight_layout()
    return fig


# --------------------------------------------------------------------------
# QC
# --------------------------------------------------------------------------


def _validate_direction_matches_sign(df: pd.DataFrame, context: str) -> None:
    """CRISPR direction labels must agree with the sign of the underlying
    effect_size: negative -> sensitising_knockout, positive ->
    tolerance_associated_knockout. Catches a corrupted/mislabelled row
    before it can be colored or grouped inconsistently in any figure.

    ``crispr_effect_size`` must be finite -- a NaN/+-inf value cannot be
    signed, so it is rejected outright rather than silently passing the
    sign comparisons below. It must also be non-zero: these candidate-facing
    tables are restricted to already-classified sensitising/tolerance genes,
    so an exact-zero (indeterminate) effect here indicates upstream
    corruption, not a legitimate third case.
    """
    effect = df["crispr_effect_size"].to_numpy(dtype=float)
    if not np.isfinite(effect).all():
        bad_genes = df.loc[~np.isfinite(effect), "gene_symbol"].tolist()
        raise ValueError(f"{context}: crispr_effect_size is not finite for {bad_genes}")
    if (effect == 0).any():
        bad_genes = df.loc[effect == 0, "gene_symbol"].tolist()
        raise ValueError(
            f"{context}: crispr_effect_size is exactly zero for {bad_genes} -- candidate tables must "
            "contain only non-zero, already-classified sensitising/tolerance effects"
        )

    direction = df["crispr_direction"].to_numpy()
    negative_mismatch = (effect < 0) & (direction != DIRECTION_SENSITISING)
    positive_mismatch = (effect > 0) & (direction != DIRECTION_TOLERANCE)
    bad = negative_mismatch | positive_mismatch
    if bad.any():
        bad_genes = df.loc[bad, "gene_symbol"].tolist()
        raise ValueError(f"{context}: crispr_direction does not match crispr_effect_size sign for {bad_genes}")


def run_qc_checks(
    evidence_df: pd.DataFrame,
    shortlist_df: pd.DataFrame,
    sensitising_df: pd.DataFrame,
    tolerance_df: pd.DataFrame,
    paics_df: pd.DataFrame,
    pca_df: pd.DataFrame,
    sample_meta_df: pd.DataFrame,
    cfg: NebulaPlotsConfig,
) -> None:
    """Programmatic QC gate run before any figure is written. Every check
    is a direct structural comparison against frozen, already-loaded
    tables -- no new computation."""
    if len(evidence_df) != cfg.expected_n_hits:
        raise ValueError(f"expected {cfg.expected_n_hits} Gate-1 hits, found {len(evidence_df)}")
    if BENCHMARK_GENE_SYMBOL in set(evidence_df["gene_symbol"]):
        raise ValueError(f"{BENCHMARK_GENE_SYMBOL} must never be counted among the {cfg.expected_n_hits} Gate-1 hits")
    for name, df in (
        ("sensitising_df", sensitising_df),
        ("tolerance_df", tolerance_df),
        ("shortlist_df", shortlist_df),
    ):
        if BENCHMARK_GENE_SYMBOL in set(df["gene_symbol"]):
            raise ValueError(f"{BENCHMARK_GENE_SYMBOL} must never appear in {name} -- it is benchmark-only")

    _validate_direction_matches_sign(evidence_df, "evidence_df")
    _validate_direction_matches_sign(sensitising_df, "sensitising_df")
    _validate_direction_matches_sign(tolerance_df, "tolerance_df")

    n_sensitising = int((evidence_df["crispr_direction"] == DIRECTION_SENSITISING).sum())
    n_tolerance = int((evidence_df["crispr_direction"] == DIRECTION_TOLERANCE).sum())
    if n_sensitising != cfg.expected_n_sensitising:
        raise ValueError(f"expected {cfg.expected_n_sensitising} sensitising_knockout genes, found {n_sensitising}")
    if n_tolerance != cfg.expected_n_tolerance:
        raise ValueError(f"expected {cfg.expected_n_tolerance} tolerance_associated_knockout genes, found {n_tolerance}")
    if len(sensitising_df) != cfg.expected_n_sensitising:
        raise ValueError(f"sensitising candidates table has {len(sensitising_df)} rows, expected {cfg.expected_n_sensitising}")
    if len(tolerance_df) != cfg.expected_n_tolerance:
        raise ValueError(f"tolerance hits table has {len(tolerance_df)} rows, expected {cfg.expected_n_tolerance}")

    primary_rows = sensitising_df.loc[sensitising_df["evidence_class"] == EVIDENCE_CLASS_PRIMARY, "gene_symbol"].tolist()
    if primary_rows != [cfg.expected_primary_gene]:
        raise ValueError(f"expected PRIMARY_RESISTANCE_SUPPORT to be exactly [{cfg.expected_primary_gene!r}], found {primary_rows!r}")
    if list(shortlist_df["gene_symbol"]) != [cfg.expected_primary_gene]:
        raise ValueError(f"shortlist is not exactly [{cfg.expected_primary_gene!r}]: {list(shortlist_df['gene_symbol'])!r}")

    if paics_df.iloc[0]["gene_symbol"] != BENCHMARK_GENE_SYMBOL or paics_df.iloc[0]["benchmark_label"] != BENCHMARK_LABEL:
        raise ValueError("PAICS benchmark row is malformed or mislabelled")

    n_pca_samples = pca_df["sample_id"].nunique()
    if n_pca_samples != cfg.expected_n_samples:
        raise ValueError(f"expected {cfg.expected_n_samples} PCA samples, found {n_pca_samples}")
    n_meta_samples = sample_meta_df["sample_id"].nunique()
    if n_meta_samples != cfg.expected_n_samples:
        raise ValueError(f"expected {cfg.expected_n_samples} samples in sample metadata, found {n_meta_samples}")

    # Exactly the expected number of samples per group (e.g. 3 MCF7 / 3
    # TAMR / 3 FASR for n=9), and no unexpected group labels.
    sample_groups = pca_df.drop_duplicates("sample_id").set_index("sample_id")["group"]
    expected_per_group = cfg.expected_n_samples // len(cfg.group_colors)
    group_counts = sample_groups.value_counts().to_dict()
    for group in cfg.group_colors:
        if group_counts.get(group, 0) != expected_per_group:
            raise ValueError(
                f"expected {expected_per_group} PCA samples for group {group!r}, found {group_counts.get(group, 0)}"
            )
    unexpected_groups = set(group_counts) - set(cfg.group_colors)
    if unexpected_groups:
        raise ValueError(f"PCA contains unexpected group labels: {unexpected_groups}")

    # PCA and sample-metadata must describe exactly the same sample set.
    pca_samples = set(pca_df["sample_id"])
    meta_samples = set(sample_meta_df["sample_id"])
    if pca_samples != meta_samples:
        raise ValueError(
            f"PCA sample set does not match sample metadata: only in PCA={pca_samples - meta_samples}, "
            f"only in metadata={meta_samples - pca_samples}"
        )

    logger.info(
        "run_qc_checks: PASSED (n_hits=%d, sensitising=%d, tolerance=%d, primary=%s, "
        "PAICS excluded from hits, n_pca_samples=%d)",
        len(evidence_df),
        n_sensitising,
        n_tolerance,
        primary_rows[0],
        n_pca_samples,
    )


# --------------------------------------------------------------------------
# Orchestration
# --------------------------------------------------------------------------


def run_nebula_plots(config_path: str | Path = "config/config.yaml") -> dict[str, object]:
    config = _load_config(config_path)
    cfg = NebulaPlotsConfig.from_config(config)

    evidence_df = load_evidence_summary(cfg)
    shortlist_df = load_shortlist(cfg)
    sensitising_df = load_sensitising_candidates(cfg)
    tolerance_df = load_tolerance_hits(cfg)
    paics_df = load_paics_benchmark(cfg)
    pca_df = load_pca_coordinates(cfg)
    sample_meta_df = load_sample_metadata(cfg)

    run_qc_checks(evidence_df, shortlist_df, sensitising_df, tolerance_df, paics_df, pca_df, sample_meta_df, cfg)

    primary_gene = cfg.expected_primary_gene
    secondary_genes = sorted(
        sensitising_df.loc[sensitising_df["evidence_class"] == EVIDENCE_CLASS_SECONDARY, "gene_symbol"].tolist()
    )

    for path in (cfg.output_dir, cfg.plot_input_dir):
        path.mkdir(parents=True, exist_ok=True)

    fig1_df, paics_inset_df = build_fig1_inputs(evidence_df, paics_df)
    fig1_df.to_csv(cfg.fig1_input_tsv, sep="\t", index=False)
    paics_inset_df.to_csv(cfg.fig1_paics_inset_tsv, sep="\t", index=False)
    fig1 = plot_fig1_crispr_landscape(fig1_df, paics_inset_df, cfg, primary_gene)
    _save_figure(fig1, cfg.fig1_png, cfg.fig1_png_transparent, cfg.fig1_pdf)

    fig2_df = build_fig2_inputs(pca_df)
    fig2_df.to_csv(cfg.fig2_input_tsv, sep="\t", index=False)
    fig2 = plot_fig2_pca(fig2_df, cfg)
    _save_figure(fig2, cfg.fig2_png, cfg.fig2_png_transparent, cfg.fig2_pdf)

    volcano_background_df = load_mcf7_volcano_source(cfg)
    fig3_df = build_fig3_inputs(volcano_background_df, evidence_df)
    fig3_df.to_csv(cfg.fig3_input_tsv, sep="\t", index=False)
    fig3 = plot_fig3_volcano(fig3_df, cfg, primary_gene, secondary_genes)
    _save_figure(fig3, cfg.fig3_png, cfg.fig3_png_transparent, cfg.fig3_pdf)

    fig4_df = build_fig4_inputs(sensitising_df)
    fig4_df.to_csv(cfg.fig4_input_tsv, sep="\t", index=False)
    fig4 = plot_fig4_evidence_matrix(fig4_df, cfg, primary_gene)
    _save_figure(fig4, cfg.fig4_png, cfg.fig4_png_transparent, cfg.fig4_pdf)

    filtered_tpm_df = load_filtered_tpm(cfg)
    fig5_summary_df = build_fig5_summary_input(evidence_df, primary_gene)
    fig5_summary_df.to_csv(cfg.fig5_summary_input_tsv, sep="\t", index=False)
    fig5_expression_df = build_fig5_expression_input(filtered_tpm_df, sample_meta_df, primary_gene)
    fig5_expression_df.to_csv(cfg.fig5_expression_input_tsv, sep="\t", index=False)
    fig5 = plot_fig5_usp34_panel(fig5_summary_df, fig5_expression_df, cfg, primary_gene)
    _save_figure(fig5, cfg.fig5_png, cfg.fig5_png_transparent, cfg.fig5_pdf)

    logger.info("run_nebula_plots: wrote 5 figures (PNG + transparent PNG + PDF) and their plot-input TSVs")

    return {
        "evidence_summary": evidence_df,
        "shortlist": shortlist_df,
        "sensitising": sensitising_df,
        "tolerance": tolerance_df,
        "paics_benchmark": paics_df,
        "fig1_input": fig1_df,
        "fig1_paics_inset": paics_inset_df,
        "fig2_input": fig2_df,
        "fig3_input": fig3_df,
        "fig4_input": fig4_df,
        "fig5_summary_input": fig5_summary_df,
        "fig5_expression_input": fig5_expression_df,
    }


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run_nebula_plots()
