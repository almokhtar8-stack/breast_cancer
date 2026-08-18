"""Candidate volcano figure, v1 -- a PROPOSED replacement for poster Figure 2.

post_freeze_exploratory. This module computes NO new statistics: every plotted
value is read from a frozen genome-wide differential-expression table. It
exists because the current Figure 2 (a within-dataset z-score heatmap) puts
significance nowhere on an axis, so a gene at FDR 0.49 can be misread as
strongly regulated. A volcano puts -log10(FDR) on the y-axis, where the fact
that only 2 of 16 candidate gene-dataset combinations reach FDR 0.05 is
visible without reading a single number.

Data sources (paths from ``config/config.yaml`` ``cross_dataset_genomewide.inputs``,
all frozen before this figure was conceived):

  A. GSE118713  ``gse118713_de_tsv``            limma, TAMR_vs_MCF7 contrast.
       The UNREDACTED table is used so KDM1A has a point to plot; the
       redacted table still carries the retired 2026-08-10 blinding and
       omits KDM1A entirely. USP34/VEZF1/TLK2 rows are identical in both.
  B. GSE111151  ``gse111151_de_tsv``            edgeR QL, resistant vs parental.
  C. GSE240112  ``gse240112_tumor_cell_tsv``    edgeR QL pseudobulk, recurrent
       vs primary, TUMOUR-CELL track -- the primary analysis per
       ``evidence_long.tsv`` (the all-epithelial track is sensitivity-only)
       and the track that reproduces the frozen candidate FDRs.
  D. GSE245601  ``gse245601_track_a_tsv``       edgeR QL pseudobulk, 12 h
       tamoxifen vs control, Track A (all epithelial). Track B is n=3
       strict-malignant exploratory and is deliberately NOT used.

Before anything is plotted, every candidate value is verified against the
frozen reference numbers from ``results/tables/cross_dataset_genomewide/
evidence_long.tsv`` / ``candidate_evidence_summary.tsv`` (the gate below);
a mismatch raises rather than plots.

Determinism: no network access, no randomness. PDF bytes are reproducible via
``SOURCE_DATE_EPOCH``; SVG element ids are pinned with ``svg.hashsalt``.
"""

from __future__ import annotations

import hashlib
import logging
import os
from dataclasses import dataclass, field
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml

# Gene identity colours are taken from the current Figure 2 renderer, not
# chosen here, so candidate identity stays consistent across the poster.
from src.poster_hero_heatmap_v6 import CONTEXT_TITLE, FOCUS_COLORS, FOCUS_FOUR

logger = logging.getLogger(__name__)

CONFIG_PATH = Path("config/config.yaml")
OUT_DIR = Path("results/figures/poster_candidate_volcano_v1")
STATUS_LABEL = "post_freeze_exploratory"

CANDIDATES: tuple[str, ...] = tuple(FOCUS_FOUR)  # KDM1A, TLK2, USP34, VEZF1

# The single significance threshold used anywhere in this figure. It is the
# project's pre-declared FDR threshold (PREANALYSIS.md); nothing here tunes it.
SIG_FDR: float = 0.05

# Exactly this many candidate points must be significant across all four
# panels (USP34 in GSE118713, VEZF1 in GSE240112) -- asserted, not assumed.
EXPECTED_SIGNIFICANT: int = 2

PANEL_ORDER: tuple[str, ...] = ("GSE118713", "GSE111151", "GSE240112", "GSE245601")

# Config keys under cross_dataset_genomewide.inputs, one per panel.
CONFIG_KEYS: dict[str, str] = {
    "GSE118713": "gse118713_de_tsv",
    "GSE111151": "gse111151_de_tsv",
    "GSE240112": "gse240112_tumor_cell_tsv",
    "GSE245601": "gse245601_track_a_tsv",
}

GENE_COLUMNS: dict[str, str] = {
    "GSE118713": "gene_symbol",
    "GSE111151": "gene_name",
    "GSE240112": "gene",
    "GSE245601": "gene",
}

# ---------------------------------------------------------------------------
# Verification gate: frozen reference values (from evidence_long.tsv /
# candidate_evidence_summary.tsv). Values are kept as STRINGS because the
# comparison tolerance depends on the precision they were quoted at: the
# brief's 1e-3, widened to half a unit in the last quoted decimal place for
# values quoted at only 2 dp (e.g. "0.91" means [0.905, 0.915]).
# These strings gate the plot; they are never themselves plotted.
# ---------------------------------------------------------------------------
REFERENCE_FDR: dict[str, dict[str, str]] = {
    "KDM1A": {"GSE118713": "0.494", "GSE111151": "0.91", "GSE240112": "0.59", "GSE245601": "0.37"},
    "TLK2": {"GSE118713": "0.675", "GSE111151": "0.84", "GSE240112": "0.88", "GSE245601": "0.88"},
    "USP34": {"GSE118713": "0.0073", "GSE111151": "0.63", "GSE240112": "0.23", "GSE245601": "0.90"},
    "VEZF1": {"GSE118713": "0.238", "GSE111151": "0.61", "GSE240112": "0.0195", "GSE245601": "0.89"},
}
REFERENCE_LOG2FC: dict[str, dict[str, str]] = {
    "USP34": {"GSE118713": "0.590"},
    "VEZF1": {"GSE118713": "0.427"},
}


def reference_tolerance(quoted: str) -> float:
    """Tolerance for one quoted reference value: the larger of 1e-3 (the
    brief's tolerance) and half a unit in the last quoted decimal place
    (so '0.91' accepts anything that rounds to 0.91)."""
    decimals = len(quoted.split(".")[1]) if "." in quoted else 0
    return max(1e-3, 0.5 * 10.0 ** (-decimals))


@dataclass
class Panel:
    """One volcano panel: a frozen DE table restricted to nothing at all --
    every gene in the table is plotted; only GSE118713 is filtered, to its
    single TAMR_vs_MCF7 contrast."""

    accession: str
    title: str
    source_path: str
    df: pd.DataFrame  # columns: gene, log2fc, fdr
    rows_in: int  # rows read from disk
    n_genes: int  # rows plotted (== len(df); nothing else is dropped)
    contrast: str | None = None
    candidate_rows: pd.DataFrame = field(default_factory=pd.DataFrame)


def _load_config(config_path: Path = CONFIG_PATH) -> dict:
    with open(config_path) as f:
        return yaml.safe_load(f)


def _select_tamr_vs_mcf7(df: pd.DataFrame) -> tuple[pd.DataFrame, str]:
    """Determine the TAMR-vs-MCF7 contrast label FROM THE FILE (not guessed):
    the unique contrast whose '_vs_'-split endpoints are exactly
    {TAMR, MCF7}."""
    labels = sorted(df["contrast"].unique())
    matches = [c for c in labels if set(c.split("_vs_")) == {"TAMR", "MCF7"}]
    if len(matches) != 1:
        raise ValueError(f"expected exactly one TAMR/MCF7 contrast among {labels}, got {matches}")
    contrast = matches[0]
    out = df[df["contrast"] == contrast]
    logger.info("GSE118713 contrast filter %r: rows in=%d, out=%d, lost=%d",
                contrast, len(df), len(out), len(df) - len(out))
    return out, contrast


def load_panels(config_path: Path = CONFIG_PATH) -> list[Panel]:
    """Load the four frozen genome-wide DE tables. No gene is dropped in any
    panel; the only filter anywhere is GSE118713's contrast selection, and
    its row accounting is logged."""
    cfg = _load_config(config_path)["cross_dataset_genomewide"]["inputs"]
    panels: list[Panel] = []
    for accession in PANEL_ORDER:
        path = cfg[CONFIG_KEYS[accession]]
        raw = pd.read_csv(path, sep="\t")
        rows_in = len(raw)
        contrast: str | None = None
        if accession == "GSE118713":
            raw, contrast = _select_tamr_vs_mcf7(raw)
        gene_col = GENE_COLUMNS[accession]
        df = raw[[gene_col, "log2fc", "fdr"]].rename(columns={gene_col: "gene"}).reset_index(drop=True)
        if df["fdr"].isna().any() or (df["fdr"] <= 0).any():
            raise ValueError(f"{accession}: FDR values must be finite and > 0 for -log10(FDR)")
        cand = df[df["gene"].isin(CANDIDATES)].copy()
        counts = cand["gene"].value_counts()
        if sorted(counts.index) != sorted(CANDIDATES) or not (counts == 1).all():
            raise ValueError(f"{accession}: expected exactly one row per candidate, got {counts.to_dict()}")
        logger.info("%s: rows in=%d, plotted=%d, lost=%d (candidates present: %d/4)",
                    accession, rows_in, len(df), rows_in - len(df) - (0 if contrast is None else 0),
                    len(cand))
        panels.append(Panel(
            accession=accession,
            title=CONTEXT_TITLE[accession],
            source_path=str(path),
            df=df,
            rows_in=rows_in,
            n_genes=len(df),
            contrast=contrast,
            candidate_rows=cand.set_index("gene"),
        ))
    return panels


def verify_against_frozen(panels: list[Panel]) -> pd.DataFrame:
    """The gate: every candidate FDR (and, where quoted, log2FC) must match
    the frozen reference numbers. Raises on ANY mismatch -- a figure that
    disagrees with the frozen tables is worse than no figure."""
    rows: list[dict] = []
    mismatches: list[str] = []
    for panel in panels:
        for gene in CANDIDATES:
            got_fdr = float(panel.candidate_rows.loc[gene, "fdr"])
            got_lfc = float(panel.candidate_rows.loc[gene, "log2fc"])
            quoted = REFERENCE_FDR[gene][panel.accession]
            tol = reference_tolerance(quoted)
            ok_fdr = abs(got_fdr - float(quoted)) <= tol
            row = {"gene": gene, "dataset": panel.accession,
                   "reference_fdr": quoted, "extracted_fdr": got_fdr,
                   "fdr_tolerance": tol, "fdr_match": ok_fdr,
                   "reference_log2fc": "", "extracted_log2fc": got_lfc,
                   "log2fc_match": ""}
            if not ok_fdr:
                mismatches.append(f"{gene}/{panel.accession} FDR: got {got_fdr}, expected {quoted} +/- {tol}")
            q_lfc = REFERENCE_LOG2FC.get(gene, {}).get(panel.accession)
            if q_lfc is not None:
                ok_lfc = abs(got_lfc - float(q_lfc)) <= 1e-3
                row["reference_log2fc"] = q_lfc
                row["log2fc_match"] = ok_lfc
                if not ok_lfc:
                    mismatches.append(f"{gene}/{panel.accession} log2FC: got {got_lfc}, expected {q_lfc} +/- 1e-3")
            rows.append(row)
    if mismatches:
        raise ValueError("verification gate FAILED; refusing to plot:\n" + "\n".join(mismatches))
    n_sig = sum(
        int(float(p.candidate_rows.loc[g, "fdr"]) < SIG_FDR) for p in panels for g in CANDIDATES
    )
    if n_sig != EXPECTED_SIGNIFICANT:
        raise ValueError(f"expected exactly {EXPECTED_SIGNIFICANT} significant candidate points, found {n_sig}")
    logger.info("verification gate passed: 16/16 values match; %d significant candidate points", n_sig)
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Rendering
# ---------------------------------------------------------------------------
BACKDROP_GREY = "#c9c9c9"
SIG_BACKDROP_GREY = "#9a9a9a"
DGRAY = "#262626"
MGRAY = "#8c8c8c"

# Deterministic per-panel label anchor positions in DATA coordinates
# (x, y, ha), tuned once by eye against the frozen values so the four
# candidate labels never collide with each other, the ticks, or the points.
# A thin leader line connects each label to its point.
LABEL_POS: dict[tuple[str, str], tuple[float, float, str]] = {
    ("GSE118713", "USP34"): (1.5, 2.45, "left"), ("GSE118713", "VEZF1"): (2.0, 1.05, "left"),
    ("GSE118713", "KDM1A"): (2.5, 0.50, "left"), ("GSE118713", "TLK2"): (-2.0, 0.50, "right"),
    ("GSE111151", "USP34"): (2.8, 0.60, "left"), ("GSE111151", "VEZF1"): (-3.0, 0.70, "right"),
    ("GSE111151", "KDM1A"): (-3.0, 0.20, "right"), ("GSE111151", "TLK2"): (2.8, 0.25, "left"),
    ("GSE240112", "USP34"): (4.5, 1.00, "left"), ("GSE240112", "VEZF1"): (3.6, 2.30, "left"),
    ("GSE240112", "KDM1A"): (5.2, 0.48, "left"), ("GSE240112", "TLK2"): (-4.5, 0.40, "right"),
    ("GSE245601", "USP34"): (1.1, 0.24, "left"), ("GSE245601", "VEZF1"): (1.1, 0.62, "left"),
    ("GSE245601", "KDM1A"): (-1.4, 0.88, "right"), ("GSE245601", "TLK2"): (-1.4, 0.38, "right"),
}


def build_figure(panels: list[Panel], stub: Path) -> dict[str, Path]:
    """Render the four-panel volcano. Saturated fill is RESERVED for candidate
    points passing FDR < SIG_FDR; non-significant candidates are hollow rings
    in their gene colour. Returns {ext: written path}."""
    # Reproducible outputs: pinned PDF creation date (as scripts/poster/
    # build_all.py does) and pinned SVG element-id salt.
    os.environ.setdefault("SOURCE_DATE_EPOCH", "1600000000")
    plt.rcParams["svg.hashsalt"] = "poster_candidate_volcano_v1"

    y_all = [-np.log10(p.df["fdr"].to_numpy()) for p in panels]
    ymax = max(float(y.max()) for y in y_all) * 1.06
    sig_y = -np.log10(SIG_FDR)

    fig, axes = plt.subplots(1, 4, figsize=(13.3, 7.5), sharey=True)
    fig.subplots_adjust(left=0.055, right=0.985, top=0.755, bottom=0.09, wspace=0.10)

    n_sig_drawn = 0
    for ax, panel, y in zip(axes, panels, y_all):
        lfc = panel.df["log2fc"].to_numpy()
        sig_mask = panel.df["fdr"].to_numpy() < SIG_FDR
        cand_mask = panel.df["gene"].isin(CANDIDATES).to_numpy()
        back = ~sig_mask & ~cand_mask
        mid = sig_mask & ~cand_mask
        # back-to-front: genome-wide backdrop, then FDR<0.05 genes, then candidates
        ax.scatter(lfc[back], y[back], s=4, c=BACKDROP_GREY, alpha=0.35,
                   linewidths=0, zorder=1, rasterized=False)
        ax.scatter(lfc[mid], y[mid], s=5, c=SIG_BACKDROP_GREY, alpha=0.45,
                   linewidths=0, zorder=2, rasterized=False)

        xmax = float(np.abs(lfc).max()) * 1.06
        ax.set_xlim(-xmax, xmax)
        ax.set_ylim(0, ymax)
        ax.axhline(sig_y, color=MGRAY, linestyle="--", linewidth=0.9, zorder=3)
        ax.axvline(0.0, color=MGRAY, linestyle="-", linewidth=0.7, alpha=0.6, zorder=3)

        for gene in CANDIDATES:
            gx = float(panel.candidate_rows.loc[gene, "log2fc"])
            gfdr = float(panel.candidate_rows.loc[gene, "fdr"])
            gy = -np.log10(gfdr)
            colour = FOCUS_COLORS[gene]
            significant = gfdr < SIG_FDR
            if significant:
                # the ONLY saturated fills in the whole figure
                ax.scatter([gx], [gy], s=110, c=colour, edgecolors="white",
                           linewidths=1.2, zorder=5)
                n_sig_drawn += 1
            else:
                ax.scatter([gx], [gy], s=95, facecolors="none", edgecolors=colour,
                           linewidths=1.8, zorder=4)
            lx, ly, ha = LABEL_POS[(panel.accession, gene)]
            label = f"{gene}\nFDR {gfdr:.4f}" if significant else gene
            ax.annotate(label, (gx, gy), xytext=(lx, ly),
                        ha=ha, va="center",
                        fontsize=9.5, color=colour if significant else DGRAY,
                        fontweight="bold" if significant else "normal", zorder=6,
                        arrowprops=dict(arrowstyle="-", color=MGRAY, lw=0.6,
                                        alpha=0.65, shrinkA=2, shrinkB=4))

        # Biological description first, accession second, smaller and lighter
        # (the Figure 2 convention). Long titles wrap at their em dash.
        title_text = panel.title.replace(" — ", "\n")
        ax.text(0.5, 1.10, title_text, transform=ax.transAxes, ha="center",
                va="bottom", fontsize=10.8, color=DGRAY, fontweight="bold")
        ax.text(0.5, 1.035, panel.accession, transform=ax.transAxes,
                ha="center", va="bottom", fontsize=9.0, color=MGRAY)
        ax.set_xlabel("log2 fold change", fontsize=10, color=DGRAY)
        ax.tick_params(labelsize=8.5, colors=DGRAY)
        for spine in ("top", "right"):
            ax.spines[spine].set_visible(False)

    if n_sig_drawn != EXPECTED_SIGNIFICANT:
        raise ValueError(f"drew {n_sig_drawn} filled candidate points, expected {EXPECTED_SIGNIFICANT}")

    axes[0].set_ylabel(r"$-\log_{10}$(FDR)  (BH-adjusted, not raw p)", fontsize=10, color=DGRAY)
    axes[0].text(axes[0].get_xlim()[0] * 0.97, sig_y + 0.09, "FDR 0.05",
                 fontsize=8.0, color=MGRAY, va="bottom", ha="left")

    fig.suptitle("Candidate significance across the four transcriptomic contexts",
                 x=0.055, y=0.980, ha="left", fontsize=17, fontweight="bold", color=DGRAY)
    fig.text(0.055, 0.935,
             "2 of 16 candidate gene-dataset combinations reach FDR 0.05 (filled points); each candidate is "
             "supported by at most one dataset. Hollow rings are candidates not reaching FDR 0.05.\nPanel 4 "
             "measures acute 12 h tamoxifen response, not resistance. Panel 3's wide fold-change range "
             "reflects pseudobulk extremes at low expression; no point is trimmed from view in any panel.",
             fontsize=9.0, color="#555555", va="top")
    # legend: one hollow + one filled exemplar in neutral dark grey, so the
    # legend explains the ENCODING without implying any one gene
    from matplotlib.lines import Line2D
    handles = [
        Line2D([], [], marker="o", linestyle="none", markersize=9,
               markerfacecolor="none", markeredgecolor=DGRAY, markeredgewidth=1.6,
               label="candidate, FDR ≥ 0.05"),
        Line2D([], [], marker="o", linestyle="none", markersize=9.5,
               markerfacecolor=DGRAY, markeredgecolor="white",
               label="candidate, FDR < 0.05"),
        Line2D([], [], marker="o", linestyle="none", markersize=5,
               markerfacecolor=SIG_BACKDROP_GREY, markeredgecolor="none",
               label="other genes, FDR < 0.05"),
        Line2D([], [], marker="o", linestyle="none", markersize=4.5,
               markerfacecolor=BACKDROP_GREY, markeredgecolor="none",
               label="all genes"),
    ]
    fig.legend(handles=handles, loc="lower left", bbox_to_anchor=(0.052, 0.858),
               ncol=4, frameon=False, fontsize=8.6, handletextpad=0.4, columnspacing=1.2)
    fig.text(0.985, 0.012, STATUS_LABEL + " — candidate figure, not the frozen poster Figure 2",
             ha="right", va="bottom", fontsize=7.0, color=MGRAY, style="italic")

    stub.parent.mkdir(parents=True, exist_ok=True)
    meta_common = {"Title": "Candidate volcano v1 (proposed Figure 2 replacement)"}
    written: dict[str, Path] = {}
    for ext, meta in (
        ("png", {**meta_common, "Description": STATUS_LABEL}),
        ("pdf", {**meta_common, "Subject": STATUS_LABEL, "Keywords": STATUS_LABEL}),
        ("svg", {**meta_common, "Description": STATUS_LABEL}),
    ):
        path = stub.with_suffix(f".{ext}")
        fig.savefig(path, facecolor="white", dpi=300, metadata=meta)
        written[ext] = path
    plt.close(fig)
    return written


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_manifests(panels: list[Panel], verification: pd.DataFrame,
                    written: dict[str, Path], out_dir: Path) -> None:
    """Manifest recording sources, the GSE118713 contrast, the candidate
    values plotted, and SHA-256 per output -- following the labelling
    convention of poster/figure_manifest.tsv (which is NOT edited)."""
    manifest = pd.DataFrame([{
        "figure_name": "poster_candidate_volcano_v1",
        "scientific_question": "How many of the 16 candidate gene-dataset combinations are actually significant?",
        "analysis_status": STATUS_LABEL,
        "post_freeze": "yes",
        "render_module": "src/poster_candidate_volcano_v1.py",
        "wrapper_script": "scripts/poster/02b_candidate_volcano.py",
        "gse118713_contrast": next(p.contrast for p in panels if p.accession == "GSE118713"),
        "significance_threshold_fdr": SIG_FDR,
        "n_significant_candidate_points": EXPECTED_SIGNIFICANT,
        **{f"source_{p.accession.lower()}": p.source_path for p in panels},
        **{f"n_genes_{p.accession.lower()}": p.n_genes for p in panels},
        **{f"sha256_{ext}": _sha256(path) for ext, path in written.items()},
    }])
    manifest.to_csv(out_dir / "volcano_manifest.tsv", sep="\t", index=False)

    values = pd.concat(
        [p.candidate_rows.reset_index().assign(dataset=p.accession,
                                               significant=lambda d: d["fdr"] < SIG_FDR)
         for p in panels], ignore_index=True,
    )[["dataset", "gene", "log2fc", "fdr", "significant"]]
    values.insert(0, "analysis_status", STATUS_LABEL)
    values.to_csv(out_dir / "candidate_values_plotted.tsv", sep="\t", index=False)

    verification = verification.copy()
    verification.insert(0, "analysis_status", STATUS_LABEL)
    verification.to_csv(out_dir / "verification_against_frozen.tsv", sep="\t", index=False)
    logger.info("wrote manifest, candidate values and verification table to %s", out_dir)


def main(out_dir: Path = OUT_DIR) -> dict[str, Path]:
    panels = load_panels()
    verification = verify_against_frozen(panels)
    written = build_figure(panels, out_dir / "candidate_volcano_v1")
    write_manifests(panels, verification, written, out_dir)
    return written


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    main()
