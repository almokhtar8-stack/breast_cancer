"""ONE poster-grade hero heatmap, v4: a REAL sample-level heatmap.

v1-v3 collapsed every dataset to a group-mean row, which reads as a
summary table rather than a genomics heatmap. v4 replaces that entirely:
every row is one real biological observation (a replicate, a cell-line
background, a tumour, or a patient/condition pseudobulk) drawn straight
from the same frozen per-sample loaders already used and tested in
`poster_exploration_v2_data` -- no new discovery, no new statistics, no
altered differential-expression result.

The only transform applied is a disclosed, standard, deterministic
visualization-only z-score: for each dataset and each gene independently,

    z = (sample_value - dataset/gene mean) / dataset/gene SD

computed across exactly the biological observations shown for that
dataset. This never touches a p-value, FDR, effect size, or the frozen
log2FC values used in v1-v3; it only re-expresses each gene's own
dataset-level distribution as unitless standard deviations so that one
shared diverging colorbar can honestly span all four dataset blocks
(z-scores are dimensionless, unlike raw TPM/CPM/log2CPM, which are not
comparable across platforms).

Row order is fixed by dataset -> biological condition -> sample/model;
rows are never hierarchically clustered, and gene columns keep the fixed
KDM1A / TLK2 / USP34 / VEZF1 order. See
results/reports/poster_hero_heatmap_v4/NOTE.md for the full source
inventory and design rationale.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.cm import ScalarMappable
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.patches import Rectangle

from src import poster_exploration_v2_data as ed
from src import poster_story_v1_data as sv1

logger = logging.getLogger(__name__)

OUT_DIR = Path("results/figures/poster_hero_heatmap_v4")

FOCUS_FOUR = sv1.FOCUS_FOUR  # ["KDM1A", "TLK2", "USP34", "VEZF1"]
FOCUS_COLORS = sv1.FOCUS_COLORS
DGRAY = "#262626"
MGRAY = "#8c8c8c"

DIVERGING_CMAP = LinearSegmentedColormap.from_list(
    "poster_diverging", ["#2E6C8E", "#7FAFC4", "#F3EEE4", "#E3A180", "#C1543A"], N=256,
)

# Restrained biological-state annotation palette. Deliberately a
# different, muted hue family from both the diverging heatmap cmap and
# the gene-header colors, so it reads as a state tag, not a data value.
STATE_COLORS = {
    "baseline": "#d9d9d9",   # MCF7 / parental / primary tumour / control
    "resistant": "#c2a24f",  # TAMR / TamR-derivative sublines (GSE118713, GSE111151)
    "recurrent": "#7c6f9e",  # recurrent tumour (GSE240112)
    "acute_tam": "#4f8fa6",  # 12h ex vivo tamoxifen arm (GSE245601)
}
STATE_LABELS = {
    "baseline": "Baseline / untreated",
    "resistant": "Resistant (TAMR)",
    "recurrent": "Recurrent tumour",
    "acute_tam": "Acute tamoxifen (12h)",
}


@dataclass
class Row:
    dataset: str
    label: str
    state: str
    subgroup: int  # index used only to add a small extra visual gap between biological sub-blocks
    values: dict = field(default_factory=dict)  # gene -> raw value (pre-z), for the z-score computation
    pair_anchor: str | None = None  # label of the row this one is paired to (bracket target), or None


def _zscore_block(rows: list[Row], genes: list[str]) -> None:
    """Gene-wise z-score across exactly the observations in `rows`
    (one dataset block). Mutates each Row's `values` dict in place,
    replacing raw values with z-scores. Documented, visualization-only;
    never touches a frozen DE statistic."""
    for gene in genes:
        raw = np.array([r.values[gene] for r in rows], dtype=float)
        mean, sd = raw.mean(), raw.std(ddof=1)
        if sd == 0 or len(raw) < 3:
            logger.warning("z-score undefined for %s (N=%d, sd=%.4f) -- rendered as neutral 0", gene, len(raw), sd)
            z = np.zeros_like(raw)
        else:
            z = (raw - mean) / sd
        for r, zi in zip(rows, z):
            r.values[gene] = float(zi)


def _build_gse118713_rows() -> list[Row]:
    """MCF7 vs TAMR replicates only. FASR (a second, independent
    fulvestrant-resistant derivative in the same series) is excluded from
    this hero figure because the poster is specifically about tamoxifen
    sensitisation, not fulvestrant resistance -- see NOTE.md."""
    df = ed.load_gse118713_focus_gene_samples()
    df = df[df["condition"].isin(["MCF7", "TAMR"])]
    rows: list[Row] = []
    for condition, state, subgroup in [("MCF7", "baseline", 0), ("TAMR", "resistant", 1)]:
        samples = sorted(df[df["condition"] == condition]["sample"].unique())
        for i, sample in enumerate(samples, start=1):
            label = f"{condition}-{i}"
            sub = df[df["sample"] == sample]
            values = {g: float(sub[sub["gene_symbol"] == g]["tpm"].iloc[0]) for g in FOCUS_FOUR}
            rows.append(Row(dataset="GSE118713", label=label, state=state, subgroup=subgroup, values=values))
    _zscore_block(rows, FOCUS_FOUR)
    return rows


def _build_gse111151_rows() -> list[Row]:
    """Real per-line parental and TamR-derivative samples, in the actual
    blocked design (4 parental lines, 7 independently-derived TamR
    sublines total -- not a 1:1 pairing). Parental-to-derivative
    relationships are shown with a bracket exactly as encoded in
    `paired_parental_sample_id`; lines are never collapsed to one row."""
    df = ed.load_gse111151_focus_gene_samples()
    line_order = list(dict.fromkeys(df.sort_values("gsm")["parental_line"]))
    rows: list[Row] = []
    for subgroup, line in enumerate(line_order):
        parental = df[(df["parental_line"] == line) & (df["status"] == "parental")]
        parental_sample_id = parental["sample_id"].iloc[0]
        values = {g: float(parental[parental["gene_symbol"] == g]["log2cpm"].iloc[0]) for g in FOCUS_FOUR}
        rows.append(Row(dataset="GSE111151", label=parental_sample_id, state="baseline",
                         subgroup=subgroup, values=values))
        derivatives = df[(df["parental_line"] == line) & (df["status"] == "resistant")]
        for derivative_id in sorted(derivatives["derivative_id"].unique()):
            der = derivatives[derivatives["derivative_id"] == derivative_id]
            der_sample_id = der["sample_id"].iloc[0]
            values = {g: float(der[der["gene_symbol"] == g]["log2cpm"].iloc[0]) for g in FOCUS_FOUR}
            rows.append(Row(dataset="GSE111151", label=der_sample_id, state="resistant",
                             subgroup=subgroup, values=values, pair_anchor=parental_sample_id))
    _zscore_block(rows, FOCUS_FOUR)
    return rows


def _build_gse240112_rows() -> list[Row]:
    """Real per-tumour pseudobulk, 3 primary + 3 recurrent -- UNPAIRED
    (different patients/biobanks). No pairing bracket is ever drawn here;
    the two conditions are separated by a group divider and label only."""
    df = ed.load_gse240112_focus_gene_tumours()
    rows: list[Row] = []
    for group, state, subgroup in [("PT", "baseline", 0), ("RT", "recurrent", 1)]:
        sub = df[df["group"] == group]
        sample_ids = sorted(sub["sample_id"].unique(), key=lambda s: int(s[2:]))
        display_prefix = "Primary" if group == "PT" else "Recurrent"
        for sample_id in sample_ids:
            n = sample_id[2:]
            values = {g: float(sub[(sub["sample_id"] == sample_id) & (sub["gene"] == g)]["log2cpm"].iloc[0])
                      for g in FOCUS_FOUR}
            rows.append(Row(dataset="GSE240112", label=f"{display_prefix}-{n}", state=state,
                             subgroup=subgroup, values=values))
    _zscore_block(rows, FOCUS_FOUR)
    return rows


def _build_gse245601_rows() -> list[Row]:
    """Real per-patient, patient-matched Control vs 12h-ex-vivo-Tamoxifen
    pseudobulk (the 3 patients passing the project's own pre-declared
    pseudobulk eligibility filter). Genuinely paired within patient, so a
    bracket is drawn; explicitly labelled as acute ex vivo, not
    resistance."""
    df = ed.load_gse245601_paired_focus_genes()
    rows: list[Row] = []
    for subgroup, patient in enumerate(sorted(df["patient"].unique())):
        short = "T" + patient.split("_")[1]
        ctrl_label = f"{short} Ctrl"
        tam_label = f"{short} TAM"
        ctrl = df[(df["patient"] == patient) & (df["condition"] == "Control")]
        values = {g: float(ctrl[ctrl["gene"] == g]["log2_expr"].iloc[0]) for g in FOCUS_FOUR}
        rows.append(Row(dataset="GSE245601", label=ctrl_label, state="baseline", subgroup=subgroup, values=values))
        tam = df[(df["patient"] == patient) & (df["condition"] == "Tamoxifen")]
        values = {g: float(tam[tam["gene"] == g]["log2_expr"].iloc[0]) for g in FOCUS_FOUR}
        rows.append(Row(dataset="GSE245601", label=tam_label, state="acute_tam", subgroup=subgroup,
                         values=values, pair_anchor=ctrl_label))
    _zscore_block(rows, FOCUS_FOUR)
    return rows


DATASET_BUILDERS = [
    ("GSE118713", _build_gse118713_rows),
    ("GSE111151", _build_gse111151_rows),
    ("GSE240112", _build_gse240112_rows),
    ("GSE245601", _build_gse245601_rows),
]


def build_hero_heatmap_v4(stub: Path) -> None:
    blocks: list[tuple[str, list[Row]]] = [(name, builder()) for name, builder in DATASET_BUILDERS]
    all_rows = [r for _, rows in blocks for r in rows]
    zmax = max(abs(r.values[g]) for r in all_rows for g in FOCUS_FOUR)
    div_norm = TwoSlopeNorm(vcenter=0, vmin=-zmax, vmax=zmax)

    row_h = 0.34
    row_gap = 0.03
    subgroup_gap = 0.10
    block_gap = 0.55
    col_w = 1.0
    n_genes = len(FOCUS_FOUR)

    n_rows = len(all_rows)
    fig_h = 3.0 + n_rows * (row_h + row_gap) * 1.05
    fig, ax = plt.subplots(figsize=(9.4, fig_h), dpi=300)

    y = 0.0
    row_y: dict[tuple[str, str], float] = {}
    block_spans: list[tuple[str, float, float]] = []
    divider_ys: list[float] = []
    group_labels: list[tuple[float, str]] = []

    for dataset, rows in blocks:
        block_top = y
        prev_subgroup = None
        for r in rows:
            if prev_subgroup is not None and r.subgroup != prev_subgroup:
                y -= subgroup_gap
                if dataset == "GSE240112":
                    divider_ys.append(y + subgroup_gap / 2)
            row_y[(dataset, r.label)] = y
            prev_subgroup = r.subgroup
            y -= (row_h + row_gap)
        block_bottom = y + row_gap
        block_spans.append((dataset, block_top, block_bottom))
        y -= block_gap

    if any(d == "GSE240112" for d, _ in blocks):
        pt_ys = [row_y[("GSE240112", r.label)] for r in dict(blocks)["GSE240112"] if r.state == "baseline"]
        rt_ys = [row_y[("GSE240112", r.label)] for r in dict(blocks)["GSE240112"] if r.state == "recurrent"]
        group_labels.append((float(np.mean(pt_ys)), "PRIMARY"))
        group_labels.append((float(np.mean(rt_ys)), "RECURRENT"))

    top_y = 0.0 + row_h * 0.5 + 0.55
    bottom_y = y + block_gap - row_h * 0.5

    ann_x0, ann_x1 = -0.62, -0.50
    label_x = -0.72
    bracket_x0, bracket_x1 = -1.72, -1.58
    dataset_x = -2.85
    grouplabel_x = -2.05

    for dataset, rows in blocks:
        for r in rows:
            yy = row_y[(dataset, r.label)]
            for j, gene in enumerate(FOCUS_FOUR):
                color = DIVERGING_CMAP(div_norm(r.values[gene]))
                ax.add_patch(Rectangle((j * col_w - 0.46, yy - row_h / 2), 0.92, row_h,
                                        facecolor=color, edgecolor="white", linewidth=0.9, zorder=2))
            ax.add_patch(Rectangle((ann_x0, yy - row_h / 2), ann_x1 - ann_x0, row_h,
                                    facecolor=STATE_COLORS[r.state], edgecolor="white", linewidth=0.6, zorder=2))
            ax.text(label_x, yy, r.label, ha="right", va="center", fontsize=7.3, color=DGRAY)

        # pairing brackets -- only where pairing is real (GSE111151 parent
        # -> derivative(s), GSE245601 patient-matched control/tamoxifen)
        if dataset in ("GSE111151", "GSE245601"):
            anchors: dict[str, list[str]] = {}
            for r in rows:
                if r.pair_anchor is not None:
                    anchors.setdefault(r.pair_anchor, []).append(r.label)
            for anchor_label, child_labels in anchors.items():
                y_anchor = row_y[(dataset, anchor_label)]
                child_ys = [row_y[(dataset, c)] for c in child_labels]
                ax.plot([bracket_x0, bracket_x0], [min(child_ys), y_anchor], color=MGRAY, linewidth=1.1,
                        solid_capstyle="round", zorder=3)
                ax.plot([bracket_x0, bracket_x1], [y_anchor, y_anchor], color=MGRAY, linewidth=1.1,
                        solid_capstyle="round", zorder=3)
                for cy in child_ys:
                    ax.plot([bracket_x0, bracket_x1], [cy, cy], color=MGRAY, linewidth=1.1,
                            solid_capstyle="round", zorder=3)

    for dataset, top, bottom in block_spans:
        mid = (top + bottom) / 2
        ax.text(dataset_x, mid, dataset, ha="left", va="center", fontsize=15.5, fontweight="bold", color=DGRAY)
        ax.plot([dataset_x - 0.15, (n_genes - 1) * col_w + 0.46], [top + row_h / 2 + 0.14] * 2,
                color="#eaeaea", linewidth=1.0, zorder=1)

    for yy, text in group_labels:
        ax.text(grouplabel_x, yy, text, ha="right", va="center", fontsize=7.6, fontweight="bold",
                 color=MGRAY, rotation=90)
    for yy in divider_ys:
        ax.plot([-0.46 - 0.02, (n_genes - 1) * col_w + 0.46], [yy, yy], color="#e3e3e3", linewidth=0.8,
                linestyle=(0, (2, 2)), zorder=1)

    header_y = top_y - 0.12
    for j, gene in enumerate(FOCUS_FOUR):
        ax.text(j * col_w, header_y, gene, ha="center", va="bottom", fontsize=19, fontweight="bold",
                 color=FOCUS_COLORS.get(gene, DGRAY))

    ax.set_xlim(-3.35, (n_genes - 1) * col_w + 0.62)
    ax.set_ylim(bottom_y - 0.05, top_y + 0.42)
    ax.axis("off")

    div_sm = ScalarMappable(norm=div_norm, cmap=DIVERGING_CMAP)
    div_sm.set_array([])
    cax = fig.add_axes([0.16, -0.012, 0.30, 0.011])
    cb = fig.colorbar(div_sm, cax=cax, orientation="horizontal")
    cb.set_ticks([-zmax, 0, zmax])
    cb.set_ticklabels(["low", "0", "high"])
    cb.ax.tick_params(labelsize=8.2, length=0, color=DGRAY, labelcolor=DGRAY)
    cb.outline.set_visible(False)
    fig.text(0.16, -0.020, "Gene-wise z-score, standardized within each dataset (Low ← 0 → High)",
              fontsize=8.4, color=DGRAY, va="top")

    legend_x = 0.60
    for i, state in enumerate(["baseline", "resistant", "recurrent", "acute_tam"]):
        sy = -0.006 - i * 0.0095
        fig.add_artist(plt.Rectangle((legend_x, sy), 0.016, 0.007, transform=fig.transFigure,
                                      facecolor=STATE_COLORS[state], edgecolor="none"))
        fig.text(legend_x + 0.024, sy + 0.0035, STATE_LABELS[state], fontsize=7.6, color=DGRAY, va="center")

    fig.text(0.0, 1.012, "Candidate expression states across resistance, recurrence and acute tamoxifen response",
              fontsize=18.5, fontweight="bold", color=DGRAY, ha="left", va="bottom")
    fig.text(0.0, 1.001,
              "Color shows within-dataset, gene-wise standardized expression; rows are real biological observations.",
              fontsize=10.3, color="#555555", ha="left", va="bottom")

    stub.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(stub.with_suffix(".png"), facecolor="white", bbox_inches="tight", dpi=300)
    fig.savefig(stub.with_suffix(".pdf"), facecolor="white", bbox_inches="tight")
    fig.savefig(stub.with_suffix(".svg"), facecolor="white", bbox_inches="tight")
    plt.close(fig)
    logger.info("wrote %s.png/.pdf/.svg (%d real biological rows)", stub, n_rows)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    build_hero_heatmap_v4(OUT_DIR / "HERO_sample_level_heatmap_v4")
