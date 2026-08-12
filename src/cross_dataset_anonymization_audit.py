"""Cross-dataset genome-wide integration, Phase 23: anonymized ranking
audit. Assigns every gene a deterministic but NON-alphabetical anonymized
ID (``GeneNNNNN``, ordered by a fixed random seed, not by gene symbol),
reruns the full percentile + global-ranking pipeline on the anonymized
table, and compares the resulting ranks to the original, named run.

Why non-alphabetical matters: the global ranking's final tie-break is
"gene symbol ascending" (a deterministic, but name-dependent, rule for
otherwise-exactly-tied genes). If genes were anonymized in alphabetical
order, the anonymized tie-break order would trivially reproduce the
original order even if some hidden bias existed elsewhere -- the test
would prove nothing. Anonymizing in a *shuffled* order and still getting
matching ranks (outside of genuine ties, which are expected and
harmless) is real evidence that no gene name influenced any percentile,
coverage, or sort computation upstream of that final tie-break.

Data source: `results/tables/cross_dataset_genomewide/all_genes_cross_dataset_evidence.tsv`
(the wide matrix *before* percentiles, so this rebuilds percentiles from
scratch on anonymized data rather than just re-sorting an already-percentiled table).
"""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from src.cross_dataset_ranking import assign_coverage_tier, build_global_ranking, compute_dataset_percentiles

logger = logging.getLogger(__name__)


def _load_config(config_path: str | Path) -> dict:
    with open(config_path) as f:
        return yaml.safe_load(f)


def build_anonymized_mapping(genes: pd.Series, seed: int = 20260812) -> pd.DataFrame:
    """Deterministic (fixed seed) but non-alphabetical shuffle -- the
    seed is fixed so the mapping is reproducible run-to-run, but the
    resulting Gene001, Gene002, ... assignment does not follow the
    genes' own alphabetical order."""
    unique_genes = sorted(genes.unique())  # sort first so the *input* to the shuffle is deterministic regardless of table row order
    rng = np.random.default_rng(seed)
    shuffled_order = rng.permutation(len(unique_genes))
    anon_ids = [f"Gene{i + 1:05d}" for i in range(len(unique_genes))]
    mapping = pd.DataFrame({"gene": [unique_genes[i] for i in shuffled_order], "anon_id": anon_ids})
    logger.info("build_anonymized_mapping: %d genes mapped, non-alphabetical (seed=%d)", len(mapping), seed)
    return mapping


def run_anonymized_ranking(wide: pd.DataFrame, mapping: pd.DataFrame, min_datasets_testable: int) -> pd.DataFrame:
    anon_wide = wide.merge(mapping, on="gene", how="inner")
    anon_wide = anon_wide.drop(columns=["gene"]).rename(columns={"anon_id": "gene"})
    with_pct = compute_dataset_percentiles(anon_wide)
    with_tier = assign_coverage_tier(with_pct)
    _full, anon_ranked = build_global_ranking(with_tier, min_datasets_testable)
    return anon_ranked


def compare_rankings(original_ranked: pd.DataFrame, anon_ranked: pd.DataFrame, mapping: pd.DataFrame) -> pd.DataFrame:
    """One row per gene present in the original ranking: original_rank,
    anon_rank (mapped back to the real gene name), and whether they
    match. A mismatch is only expected for genes that were exactly tied
    on every sort field except the (name-dependent) final tie-break."""
    anon_to_gene = mapping.set_index("anon_id")["gene"]
    anon_mapped_back = anon_ranked.copy()
    anon_mapped_back["gene"] = anon_mapped_back["gene"].map(anon_to_gene)

    merged = original_ranked[["gene", "global_rank"]].merge(
        anon_mapped_back[["gene", "global_rank"]], on="gene", how="outer", suffixes=("_original", "_anonymized")
    )
    merged["ranks_match"] = merged["global_rank_original"] == merged["global_rank_anonymized"]
    n_mismatch = int((~merged["ranks_match"]).sum())
    logger.info("compare_rankings: %d/%d genes match exactly; %d mismatches (expected only among genuine ties)", merged["ranks_match"].sum(), len(merged), n_mismatch)
    return merged


def run_anonymization_audit(config_path: str | Path = "config/config.yaml") -> dict[str, pd.DataFrame]:
    config = _load_config(config_path)
    cfg = config["cross_dataset_genomewide"]
    out = cfg["output"]
    min_testable = cfg["min_datasets_testable_for_primary_ranking"]
    tables_dir = Path(out["wide_matrix_tsv"]).parent

    wide = pd.read_csv(out["wide_matrix_tsv"], sep="\t")
    original_ranked = pd.read_csv(tables_dir / "global_ranking_eligible.tsv", sep="\t")

    mapping = build_anonymized_mapping(wide["gene"])
    anon_ranked = run_anonymized_ranking(wide, mapping, min_testable)
    comparison = compare_rankings(original_ranked, anon_ranked, mapping)

    mapping.to_csv(out["anonymized_mapping_tsv"], sep="\t", index=False)
    anon_ranked.to_csv(out["anonymized_ranking_tsv"], sep="\t", index=False)
    comparison_path = tables_dir / "anonymization_comparison.tsv"
    comparison.to_csv(comparison_path, sep="\t", index=False)
    logger.info("wrote %s, %s, %s", out["anonymized_mapping_tsv"], out["anonymized_ranking_tsv"], comparison_path)

    return {"mapping": mapping, "anon_ranked": anon_ranked, "comparison": comparison}


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run_anonymization_audit()
