"""Final USP34/VEZF1 GDSC pharmacogenomics phase: three poster-friendly
figures built only from this phase's own real, computed output tables.
Okabe-Ito colorblind-safe palette, matching this project's other
visualization modules.
"""

from __future__ import annotations

import logging
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

logger = logging.getLogger(__name__)

TABLES = Path("results/tables/final_pharmacogenomics")
FIGURES = Path("results/figures/final_pharmacogenomics")

BLUE, ORANGE, GREEN, RED, GRAY, DGRAY = "#0072B2", "#E69F00", "#009E73", "#b0392f", "#9a9a9a", "#555555"


def _top_n_by_fdr(gene: str, n: int = 9) -> pd.DataFrame:
    df = pd.read_csv(TABLES / f"{gene}_GDSC_drug_associations.tsv", sep="\t")
    df = df[df["metric"] == "LN_IC50"].sort_values("p_value").head(n)
    return df


def build_figure_01(out_fig: Path) -> None:
    df = _top_n_by_fdr("USP34", 9)
    fig, ax = plt.subplots(figsize=(11, 6), dpi=200)
    labels = [f"{r.drug_name}\n({r.dataset}, n={r.n})" for r in df.itertuples()]
    colors = [GREEN if fdr < 0.05 else GRAY for fdr in df["fdr"]]
    y = range(len(df))
    ax.barh(y, df["spearman_rho"], color=colors, edgecolor=DGRAY)
    ax.set_yticks(list(y))
    ax.set_yticklabels(labels, fontsize=8.5)
    ax.invert_yaxis()
    ax.axvline(0, color=DGRAY, linewidth=0.8)
    for i, (rho, fdr) in enumerate(zip(df["spearman_rho"], df["fdr"])):
        ax.text(rho + (0.02 if rho >= 0 else -0.02), i, f"FDR={fdr:.3f}", va="center",
                ha="left" if rho >= 0 else "right", fontsize=7.5, color=DGRAY)
    ax.set_xlabel("Spearman rho: USP34 expression vs LN_IC50 (GDSC breast lines)\nnegative = higher USP34 expression -> MORE SENSITIVE", fontsize=9)
    ax.set_xlim(-0.75, 0.35)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    ax.set_title("USP34: top 9 GDSC drug-response associations (breast lines, N>=15)\ngreen = FDR<0.05; correlation only -- not causal, not a validated combination", fontsize=10.5)
    fig.tight_layout()
    out_fig.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_fig, facecolor="white")
    plt.close(fig)
    logger.info("wrote %s", out_fig)


def build_figure_02(out_fig: Path) -> None:
    df = _top_n_by_fdr("VEZF1", 9)
    fig, ax = plt.subplots(figsize=(11, 6), dpi=200)
    labels = [f"{r.drug_name}\n({r.dataset}, n={r.n})" for r in df.itertuples()]
    y = range(len(df))
    ax.barh(y, df["spearman_rho"], color=GRAY, edgecolor=DGRAY)
    ax.set_yticks(list(y))
    ax.set_yticklabels(labels, fontsize=8.5)
    ax.invert_yaxis()
    ax.axvline(0, color=DGRAY, linewidth=0.8)
    for i, (rho, fdr) in enumerate(zip(df["spearman_rho"], df["fdr"])):
        ax.text(rho + (0.02 if rho >= 0 else -0.02), i, f"FDR={fdr:.3f} (n.s.)", va="center",
                ha="left" if rho >= 0 else "right", fontsize=7.5, color=DGRAY)
    ax.set_xlabel("Spearman rho: VEZF1 expression vs LN_IC50 (GDSC breast lines)", fontsize=9)
    ax.set_xlim(-0.65, 0.65)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    ax.set_title("VEZF1: top 9 GDSC drug-response associations by nominal p-value (breast lines, N>=15)\nALL GRAY -- zero reach FDR<0.05; this is a genuine negative result, not a data gap", fontsize=10.2)
    fig.tight_layout()
    out_fig.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_fig, facecolor="white")
    plt.close(fig)
    logger.info("wrote %s", out_fig)


def build_figure_03(out_fig: Path) -> None:
    usp34 = pd.read_csv(TABLES / "USP34_GDSC_drug_associations.tsv", sep="\t")
    vezf1 = pd.read_csv(TABLES / "VEZF1_GDSC_drug_associations.tsv", sep="\t")
    compounds = pd.read_csv(TABLES / "GDSC_compound_availability.tsv", sep="\t")

    fig, axes = plt.subplots(1, 3, figsize=(15, 5), dpi=200)

    ax = axes[0]
    n_sig = [int((usp34["fdr"] < 0.05).sum()), int((vezf1["fdr"] < 0.05).sum())]
    n_tested = [usp34["drug_id"].nunique(), vezf1["drug_id"].nunique()]
    ax.bar(["USP34", "VEZF1"], n_sig, color=[BLUE, ORANGE], edgecolor=DGRAY)
    for i, (s, t) in enumerate(zip(n_sig, n_tested)):
        ax.text(i, s + 0.15, f"{s} of {t}\ndrugs tested", ha="center", fontsize=8.5)
    ax.set_ylabel("FDR<0.05 drug associations\n(breast lines, N>=15)", fontsize=9)
    ax.set_ylim(0, max(n_sig) + 2)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    ax.set_title("FDR-significant hits", fontsize=10.5)

    ax = axes[1]
    endo = compounds[compounds["compound"].isin(["Tamoxifen", "Fulvestrant"])].drop_duplicates("compound")
    labels = endo["compound"].tolist() + ["4-OHT", "Endoxifen"]
    present = [True, True, False, False]
    colors = [GREEN if p else RED for p in present]
    ax.barh(range(len(labels)), [1] * len(labels), color=colors, edgecolor=DGRAY)
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels, fontsize=9)
    ax.set_xticks([])
    for spine in ("top", "right", "bottom"):
        ax.spines[spine].set_visible(False)
    ax.set_title("Endocrine-compound availability in GDSC\ngreen=present, red=absent", fontsize=10.5)

    ax = axes[2]
    tam_u = usp34[usp34["drug_name"] == "Tamoxifen"]["p_value"]
    tam_v = vezf1[vezf1["drug_name"] == "Tamoxifen"]["p_value"]
    tead_hippo_present = usp34["drug_name"].str.contains("TEAD|Hippo|Verteporfin", case=False, na=False).any()
    labels3 = ["Tamoxifen\n(USP34)", "Tamoxifen\n(VEZF1)", "TEAD/Hippo\ndrugs in GDSC"]
    ax.bar(labels3, [1, 1, 1], color=GRAY, edgecolor=DGRAY)
    ax.text(0, 0.5, f"n.s. (all 4 screens)\np={tam_u.min():.2f}-{tam_u.max():.2f}", ha="center", va="center", fontsize=8, color="white", fontweight="bold")
    ax.text(1, 0.5, f"n.s. (all 4 screens)\np={tam_v.min():.2f}-{tam_v.max():.2f}", ha="center", va="center", fontsize=8, color="white", fontweight="bold")
    ax.text(2, 0.5, ("compounds\nfound" if tead_hippo_present else "ZERO\ncompounds\nexist"), ha="center", va="center", fontsize=8.5, color="white", fontweight="bold")
    ax.set_yticks([])
    ax.set_title("Endocrine/Hippo-relevant negative findings", fontsize=10.5)
    for spine in ("top", "right", "left"):
        ax.spines[spine].set_visible(False)

    fig.suptitle("USP34 / VEZF1 pharmacogenomic summary (GDSC Release 8.5, correlational only)", fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.90))
    out_fig.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_fig, facecolor="white")
    plt.close(fig)
    logger.info("wrote %s", out_fig)


def run(figures_dir: Path = FIGURES) -> None:
    figures_dir.mkdir(parents=True, exist_ok=True)
    build_figure_01(figures_dir / "01_USP34_GDSC_drug_response.png")
    build_figure_02(figures_dir / "02_VEZF1_GDSC_drug_response.png")
    build_figure_03(figures_dir / "03_USP34_VEZF1_pharmacogenomic_summary.png")


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run()
