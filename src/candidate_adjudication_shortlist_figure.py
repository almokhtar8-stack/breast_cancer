"""Candidate adjudication Phase 27 (final_review item 8): a single
compact figure summarizing the four provisional, non-exclusive shortlists
side by side. No score, no ranking across lists -- each list answers a
different question (Phase 21).

Data source: `results/tables/candidate_adjudication/shortlist_*.tsv`.
"""

from __future__ import annotations

import logging
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402
import yaml  # noqa: E402

logger = logging.getLogger(__name__)

LIST_INFO = [
    ("A_multimodal_therapeutic", "A: Multimodal\ntherapeutic", "#C44E52"),
    ("B_resistance_biomarker", "B: Resistance\nbiomarker/pathway", "#4C72B0"),
    ("C_functional_sensitisation", "C: Functional\nsensitisation", "#8172B2"),
    ("D_human_tumor", "D: Human-tumor", "#55A868"),
]


def _load_config(config_path: str | Path) -> dict:
    with open(config_path) as f:
        return yaml.safe_load(f)


def plot_shortlist_summary(tables_dir: Path, out_path: Path) -> None:
    fig, axes = plt.subplots(1, 4, figsize=(12, 4.5), sharey=False)
    for ax, (fname, label, color) in zip(axes, LIST_INFO):
        df = pd.read_csv(tables_dir / f"shortlist_{fname}.tsv", sep="\t")
        genes = df["gene"].tolist()[::-1]
        ax.barh(range(len(genes)), [1] * len(genes), color=color, alpha=0.85)
        for i, g in enumerate(genes):
            ax.text(0.05, i, g, va="center", ha="left", fontsize=9, color="white", fontweight="bold")
        ax.set_yticks([])
        ax.set_xticks([])
        ax.set_xlim(0, 1)
        ax.set_title(label, fontsize=9.5)
        for spine in ax.spines.values():
            spine.set_visible(False)

    all_genes = set()
    membership = {}
    for fname, label, _ in LIST_INFO:
        df = pd.read_csv(tables_dir / f"shortlist_{fname}.tsv", sep="\t")
        for g in df["gene"]:
            membership.setdefault(g, []).append(label.split(":")[0])
            all_genes.add(g)
    multi = {g: ls for g, ls in membership.items() if len(ls) > 1}
    caption = "Genes on >1 list: " + "; ".join(f"{g} ({'+'.join(ls)})" for g, ls in sorted(multi.items())) if multi else "No gene appears on more than one shortlist."
    fig.suptitle("Four provisional, non-exclusive shortlists (no list ranked above another)", fontsize=11, y=1.03)
    fig.text(0.5, -0.02, caption, ha="center", fontsize=8, wrap=True)
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info("wrote %s (multi-list genes: %s)", out_path, list(multi.keys()))


def run_shortlist_figure(config_path: str | Path = "config/config.yaml") -> None:
    config = _load_config(config_path)
    adj = config["candidate_adjudication"]["output"]
    tables_dir = Path(adj["tables_dir"])
    final_review = Path(adj["final_review_dir"])
    plot_shortlist_summary(tables_dir, final_review / "08_shortlist_summary.png")


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run_shortlist_figure()
