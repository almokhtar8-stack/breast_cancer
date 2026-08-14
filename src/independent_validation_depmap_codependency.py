"""Independent-validation Part 10 (optional): DepMap genome-wide
co-dependency profile for each candidate, restricted to breast cancer
cell lines with a CRISPR screen. Exploratory only -- reports the top 10
positively and negatively correlated genes per candidate and never feeds
back into candidate ranking. Skipped automatically (with a logged reason)
if fewer than 20 breast lines have a CRISPR screen, per the preregistered
minimum-N-per-group threshold.
"""

from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd

from src.independent_validation_depmap_data import CANDIDATE_ENTREZ, load_config, load_model, raw_dir

logger = logging.getLogger(__name__)

CANDIDATES = ["USP34", "VEZF1", "EML5", "CITED2"]
OUT_TABLE = Path("results/tables/independent_validation/DepMap_candidate_codependency.tsv")


def build_codependency_table(cfg: dict, release: str) -> pd.DataFrame | None:
    min_n = cfg["independent_validation"]["tcga"]["thresholds"]["min_n_group"]
    model = load_model(cfg, release)
    breast_ids = set(model.index[model["is_breast"]])

    path = raw_dir(cfg, release) / cfg["independent_validation"]["depmap"]["releases"][release]["raw"]["crispr_gene_effect_csv"]
    full = pd.read_csv(path, index_col=0)
    full = full.loc[full.index.isin(breast_ids)]
    if len(full) < min_n:
        logger.info("build_codependency_table (%s): only %d breast lines with a CRISPR screen (< min_n_group=%d) -- skipping per preregistered threshold", release, len(full), min_n)
        return None
    logger.info("build_codependency_table (%s): %d breast lines with a CRISPR screen", release, len(full))

    rows = []
    for candidate in CANDIDATES:
        col = f"{candidate} ({CANDIDATE_ENTREZ[candidate]})"
        target = full[col]
        # min_periods=min_n_group guards against degenerate near-perfect
        # correlations from genes with only a handful of non-null values
        # across the breast cohort (some genes have very sparse screen
        # coverage) -- without this, a gene screened in only 2-3 of the 53
        # lines can trivially produce r=+-1.0 and dominate the top-10 list
        others = full.drop(columns=[col])
        pairwise_n = others.notna().astype(int).T.dot(target.notna().astype(int))
        others = others.loc[:, pairwise_n[pairwise_n >= min_n].index]
        corr = others.corrwith(target).dropna().sort_values()
        top_neg = corr.head(10)
        top_pos = corr.tail(10)[::-1]
        for rank, (gene, rho) in enumerate(top_pos.items(), start=1):
            rows.append(dict(candidate=candidate, depmap_release=release, direction="positive", rank=rank, gene=gene.split(" (")[0], pearson_r=float(rho)))
        for rank, (gene, rho) in enumerate(top_neg.items(), start=1):
            rows.append(dict(candidate=candidate, depmap_release=release, direction="negative", rank=rank, gene=gene.split(" (")[0], pearson_r=float(rho)))
    out = pd.DataFrame(rows)
    logger.info("build_codependency_table (%s): %d rows (4 candidates x 20 top correlations)", release, len(out))
    return out


def run(config_path: str = "config/config.yaml", out_table: Path = OUT_TABLE, release: str | None = None) -> pd.DataFrame | None:
    cfg = load_config(config_path)
    release = release or cfg["independent_validation"]["depmap"]["active_release"]
    out = build_codependency_table(cfg, release)
    if out is None:
        return None
    out_table.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(out_table, sep="\t", index=False)
    logger.info("wrote %s (%d rows, release=%s)", out_table, len(out), release)
    return out


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run()
