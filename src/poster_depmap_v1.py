"""ONE poster-grade DepMap figure -- answers "are the four focus genes
already required for baseline growth/survival of ER+/luminal breast-cancer
cells, or is their tamoxifen-sensitising CRISPR effect relatively distinct
from baseline dependency?"

Two different experiments, never conflated by this figure:
  - This project's CRISPR screen: does knocking the gene out make
    4-OHT/tamoxifen work BETTER (drug-context sensitisation)?
  - DepMap: does knocking it out already impair cancer-cell fitness under
    BASELINE, unconditioned culture?

Data source: DepMap release 26Q1 (the project's frozen `active_release` in
`config/config.yaml`), read through the project's existing frozen rules --
no new release, no new subtype definition, no new threshold:
  - ER+/luminal inclusion: `independent_validation_depmap_data.load_model()`,
    which applies DepMap's own curated `ModelSubtypeFeatures` rule
    (breast lineage AND field contains "ER+" or "ER,"), unmodified.
  - Chronos gene effect: `CRISPRGeneEffect.csv`; dependency probability:
    `CRISPRGeneDependency.csv`; both via the release's config paths.
  - "Strong baseline dependency" = dependency probability >
    `config.independent_validation.depmap.strong_dependency_probability_threshold`
    (= 0.5), the project's existing frozen criterion -- it is a probability
    cutoff, NOT a Chronos cutoff, and is labelled as such.
  - Missing values: per-gene `dropna()` on each matrix independently,
    exactly as `post_audit_sensitivity_data.load_depmap_summary_for_genes`
    does; a model must have CRISPR data to be evaluable (22 ER+/luminal
    models are annotated, 11 are CRISPR-evaluable).
  - CRISPR ranks: `post_audit_sensitivity_data.load_significant_sensitising_hits()`.

Every plotted number is derived from those sources at build time; none is
hand-typed. Scope note enforced throughout: TLK2 is the strongest baseline
dependency AMONG THE FOUR FOCUS GENES only -- the wider 13-hit universe
contains stronger baseline-dependency genes. Baseline dependency is NOT
"toxicity", is NOT evidence of normal-tissue safety (DepMap is cancer cell
lines), and high dependency is NOT "better".
"""

from __future__ import annotations

import logging
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, Normalize

from src import post_audit_sensitivity_data as pad
from src.independent_validation_depmap_data import load_model, raw_dir

logger = logging.getLogger(__name__)

OUT_DIR = Path("results/figures/poster_depmap_v1")
CACHE_DIR = Path("results/tables/poster_depmap_v1")

FOCUS_FOUR = ["KDM1A", "TLK2", "USP34", "VEZF1"]  # column order
FOCUS_COLORS = {"KDM1A": "#D55E00", "TLK2": "#56B4E9", "USP34": "#0072B2", "VEZF1": "#E69F00"}
RELEASE = "26Q1"

DGRAY = "#262626"
MGRAY = "#8c8c8c"
# Sequential DEPENDENCY palette -- deliberately NOT the expression
# heatmap's diverging red/blue up/down semantics. Colors are listed from
# the NORMALIZED-0 end (= the most negative Chronos, strongest dependency,
# darkest) to the NORMALIZED-1 end (= Chronos 0, no dependency, lightest),
# because Normalize(vmin=most_negative, vmax=0) maps the most negative
# value to 0.0. Getting this order backwards silently inverts the whole
# figure's meaning.
DEP_CMAP = LinearSegmentedColormap.from_list(
    "dependency", ["#3B2E63", "#6A5AA8", "#9C8CBF", "#CFC7DC", "#F7F5F2"]
)


def _gene_col(columns: pd.Index, symbol: str) -> str | None:
    for col in columns:
        if col.split(" (")[0] == symbol:
            return col
    return None


def _read_gene_columns(path: Path, genes: list[str]) -> pd.DataFrame:
    """Reads only the requested gene columns from a large DepMap matrix."""
    header = pd.read_csv(path, index_col=0, nrows=0)
    wanted = {g: _gene_col(header.columns, g) for g in genes}
    missing = [g for g, c in wanted.items() if c is None]
    if missing:
        raise KeyError(f"genes not found in {path.name}: {missing}")
    keep = set(wanted.values())
    df = pd.read_csv(path, index_col=0, usecols=lambda c: c in keep or c.startswith("Unnamed") or c == "ModelID")
    return df.rename(columns={v: k for k, v in wanted.items()})[genes]


def build_cellline_table(release: str = RELEASE) -> pd.DataFrame:
    """Cell-line-level Chronos gene effect + dependency probability for the
    four focus genes across the frozen ER+/luminal subset. Uses the exact
    frozen inclusion rule, threshold and missing-value handling described
    in the module docstring."""
    cfg = pad.load_config()
    rel_cfg = cfg["independent_validation"]["depmap"]["releases"][release]
    threshold = cfg["independent_validation"]["depmap"]["strong_dependency_probability_threshold"]

    model = load_model(cfg, release)
    luminal_ids = model.index[model["is_er_luminal"]]
    logger.info("ER+/luminal models annotated in %s: %d", release, len(luminal_ids))

    effect = _read_gene_columns(raw_dir(cfg, release) / rel_cfg["raw"]["crispr_gene_effect_csv"], FOCUS_FOUR)
    depprob = _read_gene_columns(raw_dir(cfg, release) / rel_cfg["raw"]["crispr_gene_dependency_csv"], FOCUS_FOUR)

    effect_lum = effect.loc[effect.index.intersection(luminal_ids)]
    depprob_lum = depprob.loc[depprob.index.intersection(luminal_ids)]
    evaluable = effect_lum.dropna(how="all").index.intersection(depprob_lum.dropna(how="all").index)
    logger.info("rows in: %d annotated ER+/luminal -> %d CRISPR-evaluable (rows lost: %d, no CRISPR data)",
                len(luminal_ids), len(evaluable), len(luminal_ids) - len(evaluable))

    rows = []
    for model_id in evaluable:
        name = model.loc[model_id, "StrippedCellLineName"]
        for gene in FOCUS_FOUR:
            chronos = effect_lum.loc[model_id, gene]
            prob = depprob_lum.loc[model_id, gene]
            rows.append({
                "model_id": model_id,
                "cell_line": name,
                "subtype_features": model.loc[model_id, "ModelSubtypeFeatures"],
                "gene": gene,
                "chronos_gene_effect": float(chronos),
                "dependency_probability": float(prob),
                "strongly_dependent": bool(prob > threshold),
            })
    table = pd.DataFrame(rows)
    assert len(table) == len(evaluable) * len(FOCUS_FOUR)
    logger.info("cell-line table: %d rows (%d lines x %d genes)", len(table), len(evaluable), len(FOCUS_FOUR))
    return table


def load_cellline_table(release: str = RELEASE, use_cache: bool = True) -> pd.DataFrame:
    """Cached view of `build_cellline_table` -- the DepMap source matrices
    are ~440 MB each, so the derived per-line values are cached to a small
    TSV. The cache is a reshaped extract of frozen DepMap data, not a new
    scientific result; delete it to force a rebuild."""
    cache = CACHE_DIR / f"depmap_{release}_er_luminal_cellline_matrix.tsv"
    if use_cache and cache.exists():
        return pd.read_csv(cache, sep="\t")
    table = build_cellline_table(release)
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    table.to_csv(cache, sep="\t", index=False)
    logger.info("wrote cache %s", cache)
    return table


def dependency_summary(table: pd.DataFrame) -> pd.DataFrame:
    """Per-gene n/N and percentage of ER+/luminal lines meeting the frozen
    strong-dependency criterion, plus median Chronos -- all computed, never
    hand-typed."""
    rows = []
    for gene in FOCUS_FOUR:
        sub = table[table["gene"] == gene]
        n_strong = int(sub["strongly_dependent"].sum())
        rows.append({
            "gene": gene,
            "n_lines": len(sub),
            "n_strongly_dependent": n_strong,
            "pct_strongly_dependent": 100.0 * n_strong / len(sub),
            "median_chronos": float(sub["chronos_gene_effect"].median()),
        })
    return pd.DataFrame(rows)


def crispr_ranks() -> dict[str, tuple[int, int]]:
    """gene -> (rank_by_effect, total hits) from the frozen sensitising-hit
    table -- the same loader the CRISPR discovery figure uses."""
    hits = pad.load_significant_sensitising_hits()
    total = len(hits)
    return {r["gene"]: (int(r["rank_by_effect"]), total)
            for _, r in hits[hits["gene"].isin(FOCUS_FOUR)].iterrows()}


def build_depmap_v1(stub: Path) -> None:
    table = load_cellline_table()
    summary = dependency_summary(table)
    ranks = crispr_ranks()

    matrix = table.pivot(index="cell_line", columns="gene", values="chronos_gene_effect")[FOCUS_FOUR]
    strong = table.pivot(index="cell_line", columns="gene", values="strongly_dependent")[FOCUS_FOUR]
    # Order cell lines by mean dependency across the four genes (strongest
    # at the top) -- a display ordering, computed from the data.
    order = matrix.mean(axis=1).sort_values().index
    matrix = matrix.loc[order]
    strong = strong.loc[order]

    n_lines = len(matrix)
    vmin = float(np.floor(matrix.min().min() * 10) / 10)
    norm = Normalize(vmin=vmin, vmax=0.0)

    fig = plt.figure(figsize=(9.6, 8.8), dpi=300)
    gs = fig.add_gridspec(2, 1, height_ratios=[n_lines, 2.9], hspace=0.09,
                          left=0.20, right=0.90, top=0.845, bottom=0.115)
    ax = fig.add_subplot(gs[0])
    ax_bar = fig.add_subplot(gs[1])

    ax.imshow(matrix.to_numpy(), cmap=DEP_CMAP, norm=norm, aspect="auto")

    # Subtle ring on cells meeting the frozen strong-dependency criterion
    # (a dependency-PROBABILITY rule, not a Chronos cutoff).
    for i in range(n_lines):
        for j in range(len(FOCUS_FOUR)):
            if strong.iat[i, j]:
                ax.scatter([j], [i], s=52, facecolor="none", edgecolor="white", linewidth=1.9, zorder=4)

    ax.set_xticks(range(len(FOCUS_FOUR)))
    ax.set_xticklabels([])
    ax.set_yticks(range(n_lines))
    ax.set_yticklabels(matrix.index, fontsize=11, color=DGRAY)
    ax.tick_params(length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_xlim(-0.5, len(FOCUS_FOUR) - 0.5)

    # Candidate names above the columns, in their poster identity colors,
    # with the frozen CRISPR functional rank as small secondary text.
    for j, gene in enumerate(FOCUS_FOUR):
        rank, total = ranks[gene]
        ax.text(j, -1.02, gene, fontsize=14.5, fontweight="bold", color=FOCUS_COLORS[gene],
                 ha="center", va="bottom")
        ax.text(j, -0.62, f"CRISPR {rank}/{total}", fontsize=8.6, color=MGRAY, ha="center", va="bottom")

    # ---- summary strip: % of lines meeting the frozen criterion ----------
    pct = summary.set_index("gene").loc[FOCUS_FOUR, "pct_strongly_dependent"]
    n_strong = summary.set_index("gene").loc[FOCUS_FOUR, "n_strongly_dependent"]
    ax_bar.axhline(0, color="#d5d2cd", linewidth=1.0, zorder=1)
    ax_bar.bar(range(len(FOCUS_FOUR)), pct.to_numpy(), width=0.55,
                color=[DEP_CMAP(0.10) for _ in pct], zorder=3)
    for j, gene in enumerate(FOCUS_FOUR):
        label = f"{n_strong[gene]}/{n_lines}  ({pct[gene]:.0f}%)"
        ax_bar.text(j, pct[gene] + 5, label, fontsize=10.5, color=DGRAY, ha="center", va="bottom",
                     fontweight="bold" if pct[gene] > 0 else "normal")
    ax_bar.set_xlim(-0.5, len(FOCUS_FOUR) - 0.5)
    ax_bar.set_ylim(0, 100)
    ax_bar.set_xticks([])
    ax_bar.set_yticks([])
    for spine in ax_bar.spines.values():
        spine.set_visible(False)
    ax_bar.set_ylabel("Lines with strong\nbaseline dependency", fontsize=10, color=DGRAY, rotation=0,
                       ha="right", va="center", labelpad=14)

    fig.text(0.035, 0.975, "DepMap separates tamoxifen sensitisation from baseline dependency",
              fontsize=17.5, fontweight="bold", color=DGRAY, ha="left", va="top")
    fig.text(0.035, 0.941,
              f"Chronos gene-effect scores across the frozen ER+/luminal breast-cancer cell-line subset "
              f"(DepMap {RELEASE}, n = {n_lines}).",
              fontsize=10.5, color="#555555", ha="left", va="top")

    # ---- colorbar: more negative = stronger baseline dependency ---------
    # Vertical bar: bottom = most negative Chronos (darkest, strongest
    # dependency), top = 0 (lightest, little/no dependency).
    cax = fig.add_axes([0.905, 0.44, 0.020, 0.30])
    gradient = np.linspace(vmin, 0.0, 256).reshape(-1, 1)
    cax.imshow(gradient, aspect="auto", cmap=DEP_CMAP, norm=norm,
                extent=[0, 1, vmin, 0.0], origin="lower")
    cax.set_xticks([])
    cax.yaxis.tick_right()
    cax.set_yticks([vmin, 0.0])
    cax.set_yticklabels([f"{vmin:.1f}", "0"], fontsize=8.5, color=DGRAY)
    cax.tick_params(axis="y", length=2, pad=2)
    for spine in cax.spines.values():
        spine.set_color("#cccccc")
    cax.set_title("Chronos", fontsize=9, color=MGRAY, pad=6)
    cax.annotate("stronger\nbaseline\ndependency", xy=(3.1, vmin), xycoords=cax.get_yaxis_transform(),
                  fontsize=8.6, color=DGRAY, ha="left", va="bottom")
    cax.annotate("little/no\ndependency", xy=(3.1, 0.0), xycoords=cax.get_yaxis_transform(),
                  fontsize=8.6, color=MGRAY, ha="left", va="top")

    fig.text(0.035, 0.055,
              "More negative Chronos score = stronger baseline dependency.   "
              f"○ = dependency probability > 0.5 (frozen strong-dependency criterion).",
              fontsize=9.6, color=MGRAY, ha="left", va="center")
    fig.text(0.035, 0.028,
              "Baseline dependency in cancer cell lines is a different measurement from tamoxifen-specific "
              "sensitisation, and is not evidence of normal-tissue safety.",
              fontsize=8.8, color=MGRAY, ha="left", va="center", style="italic")

    stub.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(stub.with_suffix(".png"), facecolor="white", bbox_inches="tight", dpi=300)
    fig.savefig(stub.with_suffix(".pdf"), facecolor="white", bbox_inches="tight")
    fig.savefig(stub.with_suffix(".svg"), facecolor="white", bbox_inches="tight")
    plt.close(fig)
    logger.info("wrote %s.png/.pdf/.svg", stub)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    build_depmap_v1(OUT_DIR / "DEPMAP_v1")
