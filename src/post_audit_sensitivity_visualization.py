"""POST-AUDIT SENSITIVITY ANALYSIS -- diagnostic figures only.

These are diagnostic figures to understand the sensitivity analysis, NOT a
poster redesign. All 4 figures read exclusively from
src.post_audit_sensitivity_data (which itself reads only already-frozen
project tables plus the unmodified original freeze code) -- no new
analysis, no recomputation, no hidden weighting.
"""

from __future__ import annotations

import logging
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

from src import post_audit_sensitivity_data as pad

logger = logging.getLogger(__name__)

FIGURES = Path("results/figures/post_audit")

DGRAY, GRAY, LGRAY = "#333333", "#9a9a9a", "#d8d8d8"
FOCUS_COLORS = {"KDM1A": "#D55E00", "TLK2": "#CC79A7", "USP34": "#0072B2", "VEZF1": "#E69F00"}
GOOD, BAD = "#009E73", "#b0392f"


def _save(fig, stub: Path, vector: bool = True) -> None:
    stub.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(stub.with_suffix(".png"), facecolor="white", bbox_inches="tight", dpi=300)
    if vector:
        fig.savefig(stub.with_suffix(".pdf"), facecolor="white", bbox_inches="tight")
    plt.close(fig)
    logger.info("wrote %s", stub)


def _clean(ax):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


# ---------------------------------------------------------------------------
# Figure A -- full Hany CRISPR landscape
# ---------------------------------------------------------------------------

def build_figure_a_crispr_landscape(stub: Path) -> None:
    df = pad.load_genomewide_crispr().sort_values("effect_size").reset_index(drop=True)
    df["rank"] = np.arange(1, len(df) + 1)
    gate1 = df["fdr"] < pad.CRISPR_GATE1_FDR
    sens13 = set(pad.load_significant_sensitising_hits()["gene"])

    fig, ax = plt.subplots(figsize=(12, 8), dpi=300)
    ax.scatter(df.loc[~gate1, "rank"], df.loc[~gate1, "effect_size"], s=4, color=GRAY, alpha=0.2,
               linewidth=0, rasterized=True, zorder=1)
    other_sens = df[df["gene"].isin(sens13 - set(pad.FOCUS_FOUR))]
    ax.scatter(other_sens["rank"], other_sens["effect_size"], s=26, color=DGRAY, alpha=0.75, edgecolor="white",
               linewidth=0.4, zorder=3, label=f"other significant sensitising hits (n={len(other_sens)})")
    tol = df[gate1 & (df["effect_size"] > 0)]
    ax.scatter(tol["rank"], tol["effect_size"], s=14, color=GRAY, alpha=0.5, zorder=2,
               label=f"tolerance-associated Gate-1 hits (n={len(tol)})")

    ax.axhline(0, color=GRAY, linewidth=0.8, zorder=0)

    ordered = sorted(pad.FOCUS_FOUR, key=lambda g: df.loc[df["gene"] == g, "rank"].iloc[0])
    label_ys = np.linspace(1.3, -0.5, len(ordered))
    label_x0 = df["rank"].min() + 90
    for gene, ly in zip(ordered, label_ys):
        row = df[df["gene"] == gene].iloc[0]
        ax.scatter([row["rank"]], [row["effect_size"]], s=220, color=FOCUS_COLORS[gene], edgecolor=DGRAY,
                   linewidth=1.6, zorder=5)
        ax.plot([label_x0, row["rank"]], [ly, row["effect_size"]], color=FOCUS_COLORS[gene], linewidth=1.0, zorder=4)
        ax.text(label_x0 + 30, ly, f"{gene}  (effect rank {int(row['rank'])}/19,103, FDR={row['fdr']:.2g})",
                ha="left", va="center", fontsize=10.5, fontweight="bold", color=FOCUS_COLORS[gene])

    ax.set_xlabel(f"Gene rank, sorted by CRISPR effect size (n={len(df):,} fitted genes)\n"
                  r"$\bf{negative}$ = sensitising knockout under 4-OHT   |   positive = tolerance-associated knockout", fontsize=11)
    ax.set_ylabel("CRISPR effect size (Hany et al. screen)", fontsize=12)
    _clean(ax)
    ax.legend(loc="lower right", fontsize=9.5, frameon=False)
    ax.set_title("A. Genome-wide CRISPR landscape -- all 13 significant sensitising hits shown honestly", fontsize=13, fontweight="bold", loc="left", pad=10)
    _save(fig, stub)


# ---------------------------------------------------------------------------
# Figure B -- sensitising-hit evidence matrix
# ---------------------------------------------------------------------------

def build_figure_b_evidence_matrix(stub: Path) -> None:
    em = pad.build_evidence_matrix().set_index("gene")
    order = pad.load_significant_sensitising_hits()["gene"].tolist()  # already rank-by-effect order
    em = em.loc[order]

    cols = [
        ("crispr_effect", "CRISPR\neffect", "crispr_fdr"),
        ("gse118713_log2fc", "GSE118713\nlog2FC", "gse118713_fdr"),
        ("gse111151_log2fc", "GSE111151\nlog2FC", "gse111151_fdr"),
        ("gse240112_log2fc", "GSE240112\nlog2FC", "gse240112_fdr"),
        ("gse245601_acute_log2fc", "GSE245601\n[ACUTE] log2FC", "gse245601_acute_fdr"),
        ("tcga_tumor_vs_normal_log2diff", "TCGA\ntumor-normal", "tcga_fdr"),
        ("median_chronos_er_luminal", "DepMap\nChronos", None),
    ]

    fig, ax = plt.subplots(figsize=(11, 8.5), dpi=300)
    n_rows, n_cols = len(order), len(cols)
    for j, (col, label, fdr_col) in enumerate(cols):
        vals = em[col]
        vmax = np.nanmax(np.abs(vals)) if vals.notna().any() else 1.0
        for i, gene in enumerate(order):
            v = vals.loc[gene]
            if pd.isna(v):
                ax.text(j, i, "NA", ha="center", va="center", fontsize=8, color=GRAY)
                continue
            frac = 0.5 + 0.5 * (v / vmax) if vmax else 0.5
            frac = min(max(frac, 0.02), 0.98)
            color = plt.cm.RdBu_r(frac)
            sig = fdr_col is not None and pd.notna(em.loc[gene, fdr_col]) and em.loc[gene, fdr_col] < 0.05
            ax.add_patch(plt.Rectangle((j - 0.48, i - 0.48), 0.96, 0.96, facecolor=color,
                                        edgecolor=DGRAY if sig else "white", linewidth=1.8 if sig else 0.5, zorder=2))
            ax.text(j, i, f"{v:.2f}{'*' if sig else ''}", ha="center", va="center", fontsize=8, color=DGRAY, zorder=3)

    ax.set_xlim(-0.6, n_cols - 0.4)
    ax.set_ylim(-0.6, n_rows - 0.4)
    ax.invert_yaxis()
    ax.set_xticks(range(n_cols))
    ax.set_xticklabels([c[1] for c in cols], fontsize=9)
    for j, (_, _, _) in enumerate(cols):
        if cols[j][1].startswith("GSE245601"):
            ax.get_xticklabels()[j].set_color(BAD)
    ax.set_yticks(range(n_rows))
    ylabels = [f"{g} *" if g in pad.FOCUS_FOUR else g for g in order]
    ax.set_yticklabels(ylabels, fontsize=9.5)
    for i, g in enumerate(order):
        if g in pad.FOCUS_FOUR:
            ax.get_yticklabels()[i].set_color(FOCUS_COLORS[g])
            ax.get_yticklabels()[i].set_fontweight("bold")
    ax.tick_params(length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_title("B. Sensitising-hit evidence matrix (13 genes x 7 evidence sources, real values, * = FDR<0.05)", fontsize=12.5, fontweight="bold", loc="left", pad=10)
    fig.text(0.01, -0.02, "NA = not assessed in this project (evidence gap, not a null result). GSE245601 [ACUTE] is never counted as resistance corroboration. TCGA/DepMap computed only where project data exists.",
             fontsize=8.5, color=GRAY, style="italic")
    _save(fig, stub)


# ---------------------------------------------------------------------------
# Figure C -- Hany x DepMap for significant sensitising hits
# ---------------------------------------------------------------------------

def build_figure_c_hany_vs_depmap(stub: Path) -> None:
    em = pad.build_evidence_matrix()
    em = em[em["in_depmap"] == True]  # noqa: E712

    fig, ax = plt.subplots(figsize=(10, 8), dpi=300)
    other = em[~em["gene"].isin(pad.FOCUS_FOUR)]
    sizes = 30 + 260 * other["frac_strongly_dependent_er_luminal"].fillna(0)
    ax.scatter(other["crispr_effect"], other["median_chronos_er_luminal"], s=sizes, color=GRAY, alpha=0.5,
               edgecolor="white", linewidth=0.5, zorder=2)
    for _, r in other.iterrows():
        ax.annotate(r["gene"], (r["crispr_effect"], r["median_chronos_er_luminal"]), xytext=(6, 4),
                    textcoords="offset points", fontsize=8, color=GRAY)

    for gene in pad.FOCUS_FOUR:
        hit = em[em["gene"] == gene]
        if len(hit) == 0:
            continue
        r = hit.iloc[0]
        size = 80 + 420 * (r["frac_strongly_dependent_er_luminal"] or 0)
        ax.scatter([r["crispr_effect"]], [r["median_chronos_er_luminal"]], s=size, color=FOCUS_COLORS[gene],
                   edgecolor=DGRAY, linewidth=1.8, zorder=5)
        ax.annotate(f"{gene} ({r['frac_strongly_dependent_er_luminal']*100:.0f}% dep.)",
                    (r["crispr_effect"], r["median_chronos_er_luminal"]), xytext=(10, 8),
                    textcoords="offset points", fontsize=10, color=FOCUS_COLORS[gene], fontweight="bold")

    ax.axvline(0, color=GRAY, linewidth=0.8, zorder=0)
    ax.axhline(0, color=GRAY, linewidth=0.8, zorder=0)
    ax.set_xlabel("Hany CRISPR drug-context effect  (more negative = stronger sensitising KO)", fontsize=11.5)
    ax.set_ylabel("Median DepMap 26Q1 Chronos gene effect, ER+/luminal (n=11)\n(more negative = stronger baseline cancer-cell dependency)", fontsize=11)
    _clean(ax)
    handles = [Line2D([0], [0], marker="o", color="none", markerfacecolor=GRAY, markeredgecolor=DGRAY, markersize=11,
                       label="marker size = % ER+/luminal lines dependency-probability>0.5")]
    ax.legend(handles=handles, loc="lower right", fontsize=8.8, frameon=False)
    ax.set_title("C. Drug-context sensitisation vs baseline dependency (11 of 13 significant sensitising hits with DepMap data)", fontsize=12, fontweight="bold", loc="left", pad=10)
    fig.text(0.01, -0.03, "Axis descriptions are literal, not proven mechanistic classes. USP17L29 and ICK have no DepMap 26Q1 gene-effect column and are omitted.", fontsize=8.3, color=GRAY, style="italic")
    _save(fig, stub)


# ---------------------------------------------------------------------------
# Figure D -- selection-rule stability
# ---------------------------------------------------------------------------

def build_figure_d_rule_stability(stub: Path) -> None:
    srs = pad.build_selection_rule_sensitivity()
    rule_order = ["RULE_0_original_frozen_gate", "RULE_1_crispr_only_no_rna_gate", "RULE_2_chronic_rna_corroboration",
                  "RULE_3_gse111151_specific", "RULE_4_human_evidence_first"]
    rule_labels = ["Rule 0\noriginal gate", "Rule 1\nCRISPR only", "Rule 2\nnon-acute RNA", "Rule 3\nGSE111151 only", "Rule 4\nhuman evidence\nfirst"]

    fig, ax = plt.subplots(figsize=(11, 7), dpi=300)
    max_rank = int(np.nanmax(srs["rank"])) + 1
    for gene in pad.FOCUS_FOUR:
        sub = srs[srs["gene"] == gene].set_index("rule").loc[rule_order]
        ys = [max_rank - r if pd.notna(r) else np.nan for r in sub["rank"]]
        xs = np.arange(len(rule_order))
        valid = ~np.isnan(ys)
        ax.plot(xs[valid], np.array(ys)[valid], color=FOCUS_COLORS[gene], linewidth=2.2, marker="o", markersize=9,
                zorder=3, label=gene)
        for x, y, elig in zip(xs, ys, sub["eligible"]):
            if not elig:
                ax.scatter([x], [max_rank - 0.4], marker="x", s=90, color=FOCUS_COLORS[gene], zorder=4)

    ax.set_xticks(range(len(rule_order)))
    ax.set_xticklabels(rule_labels, fontsize=9.5)
    ax.set_ylabel("rank among eligible genes (higher on plot = better rank)", fontsize=10.5)
    ax.set_yticks([])
    _clean(ax)
    ax.spines["left"].set_visible(False)
    ax.legend(loc="upper left", bbox_to_anchor=(1.0, 1.0), fontsize=10.5, frameon=False)
    ax.set_title("D. Selection-rule stability -- KDM1A/TLK2/USP34/VEZF1 rank across 5 rules ('x' = ineligible)", fontsize=12.5, fontweight="bold", loc="left", pad=10)
    fig.text(0.01, -0.04, "Rule 3's eligible set is empty for these 4 genes (no gene in the 13-gene universe reaches GSE111151 FDR<0.05) -- all 4 show 'x' under Rule 3.",
             fontsize=8.5, color=GRAY, style="italic")
    _save(fig, stub)


def run(figures_dir: Path = FIGURES) -> None:
    figures_dir.mkdir(parents=True, exist_ok=True)
    build_figure_a_crispr_landscape(figures_dir / "A_crispr_landscape")
    build_figure_b_evidence_matrix(figures_dir / "B_sensitising_hit_evidence_matrix")
    build_figure_c_hany_vs_depmap(figures_dir / "C_hany_vs_depmap")
    build_figure_d_rule_stability(figures_dir / "D_selection_rule_stability")


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run()
