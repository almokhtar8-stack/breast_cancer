"""Candidate adjudication Phases 21-22: four provisional, non-exclusive
shortlists derived by transparent hierarchical adjudication (never a
hidden numeric score). Each dimension (CRISPR strength, resistance-RNA
evidence, human evidence) is compared in *bands* (VERY_STRONG/STRONG/
MODERATE/WEAK/NO_EVIDENCE, from `src.candidate_adjudication_axes`), not
by sorting on raw FDR to the last decimal -- a within-band "tie" is
expected and is exactly where the next criterion in the hierarchy is
allowed to decide the order; this is what keeps the ranking a disclosed
adjudication rather than a disguised composite score.

List A (multimodal therapeutic): sensitising CRISPR direction is a hard
gate; candidates also need >=1 resistance-RNA dataset FDR<0.05 or
significant human-tumor support (otherwise they belong to List C, not A).
Ranked by, in order: CRISPR-strength band -> resistance-evidence band ->
human-evidence band -> absence of major direction contradiction
(resistance_direction_consensus all_*/majority_*/mixed, ascending
contradiction) -> leave-one-out stability band (ROBUST < MODERATELY_STABLE
< DATASET_DEPENDENT) -> gene symbol (deterministic final tie-break only).
Every one of these is an already-computed, disclosed column -- no new
statistic and no numeric composite of them.

List B (resistance biomarker/pathway): the frozen resistance-consensus
leaderboard's own top genes, CRISPR not required.

List C (functional sensitisation): the frozen CRISPR-sensitising
leaderboard's own top genes, RNA not required.

List D (human-tumor): the frozen human-only leaderboard's own top genes.

Data source: the Phase 5-8, 20 adjudication tables (all already computed,
read-only here).
"""

from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd
import yaml

logger = logging.getLogger(__name__)

BAND_RANK = {"VERY_STRONG": 0, "STRONG": 1, "MODERATE": 2, "WEAK": 3, "DISCORDANT": 3, "NO_EVIDENCE": 4}


def _load_config(config_path: str | Path) -> dict:
    with open(config_path) as f:
        return yaml.safe_load(f)


CONTRADICTION_BAND = {"all_up": 0, "all_down": 0, "majority_up": 1, "majority_down": 1, "mixed": 2, "insufficient": 3}
STABILITY_BAND = {"ROBUST": 0, "MODERATELY_STABLE": 1, "DATASET_DEPENDENT": 2}


def build_list_a_multimodal_therapeutic(pool: pd.DataFrame, axes: pd.DataFrame, stability: pd.DataFrame, max_n: int = 5) -> pd.DataFrame:
    """`crispr_direction` is an unconditional sign label (every testable
    gene gets `sensitising_KO` or `tolerance_associated_KO` regardless of
    significance -- see `src.cross_dataset_evidence_tables._crispr_direction`),
    so a "sensitising direction" gate on that column alone is nearly
    meaningless: half of all tested genes carry it by chance. The gate
    used here additionally requires at least WEAK CRISPR evidence (axis A
    != NO_EVIDENCE, i.e. CRISPR FDR<0.25 and top-10% within the screen) --
    a gene with no real CRISPR signal at all (e.g. FDR~0.8) must not enter
    a *therapeutic* shortlist merely because its noisy effect happened to
    be negative."""
    sensitising = pool.loc[pool["crispr_direction"] == "sensitising_KO", "gene"].tolist()
    m = axes.loc[axes["gene"].isin(sensitising) & (axes["axis_a_functional"] != "NO_EVIDENCE")].copy()
    # both GSE245601 tracks included in the human-significance check, matching axis C's own
    # definition (classify_axis_c_human treats either track as acute evidence) -- a Track-B
    # -only-significant gene must not be excluded here just because Track A wasn't checked
    m = m.merge(pool[["gene", "resistance_fdr05_count", "resistance_direction_consensus", "gse240112_tumor_fdr", "gse245601_epi_fdr", "gse245601_malignant_fdr"]], on="gene", how="left")
    m = m.merge(stability[["gene", "stability_label"]], on="gene", how="left")
    human_sig = (m["gse240112_tumor_fdr"] < 0.05) | (m["gse245601_epi_fdr"] < 0.05) | (m["gse245601_malignant_fdr"] < 0.05)
    multimodal_eligible = m.loc[(m["resistance_fdr05_count"].fillna(0) >= 1) | human_sig.fillna(False)].copy()

    multimodal_eligible["_a_band"] = multimodal_eligible["axis_a_functional"].map(BAND_RANK)
    multimodal_eligible["_b_band"] = multimodal_eligible["axis_b_resistance"].map(BAND_RANK)
    multimodal_eligible["_c_band"] = multimodal_eligible["axis_c_human"].map(BAND_RANK)
    multimodal_eligible["_contradiction_band"] = multimodal_eligible["resistance_direction_consensus"].map(CONTRADICTION_BAND).fillna(3)
    multimodal_eligible["_stability_band"] = multimodal_eligible["stability_label"].map(STABILITY_BAND).fillna(3)
    ranked = multimodal_eligible.sort_values(
        by=["_a_band", "_b_band", "_c_band", "_contradiction_band", "_stability_band", "gene"],
        ascending=[True, True, True, True, True, True], kind="mergesort",
    )

    out = ranked[["gene", "axis_a_functional", "axis_b_resistance", "axis_c_human", "resistance_direction_consensus", "stability_label", "global_rank"]].head(max_n).reset_index(drop=True)
    out.insert(0, "shortlist_rank", range(1, len(out) + 1))
    logger.info("build_list_a_multimodal_therapeutic: %d genes (of %d hard-gate-eligible)", len(out), len(multimodal_eligible))
    return out


def build_list_b_resistance_biomarker(top_resistance: pd.DataFrame, max_n: int = 5) -> pd.DataFrame:
    out = top_resistance.head(max_n)[["gene", "resistance_direction_consensus", "resistance_fdr05_count", "crispr_fdr", "global_rank"]].reset_index(drop=True)
    out.insert(0, "shortlist_rank", range(1, len(out) + 1))
    return out


def build_list_c_functional_sensitisation(top_crispr: pd.DataFrame, max_n: int = 5) -> pd.DataFrame:
    out = top_crispr.head(max_n)[["gene", "crispr_fdr", "crispr_evidence_percentile", "resistance_direction_consensus", "global_rank"]].reset_index(drop=True)
    out.insert(0, "shortlist_rank", range(1, len(out) + 1))
    return out


def build_list_d_human_tumor(top_human: pd.DataFrame, wide: pd.DataFrame, max_n: int = 5) -> pd.DataFrame:
    out = top_human.head(max_n).merge(wide[["gene", "gse245601_epi_fdr", "gse240112_tumor_fdr", "crispr_fdr"]], on="gene", how="left")
    out = out[["gene", "n_datasets_fdr05", "gse245601_epi_fdr", "gse240112_tumor_fdr", "crispr_fdr", "rank"]].rename(columns={"rank": "human_only_rank"})
    out = out.head(max_n).reset_index(drop=True)
    out.insert(0, "shortlist_rank", range(1, len(out) + 1))
    return out


def run_shortlists(config_path: str | Path = "config/config.yaml") -> dict[str, pd.DataFrame]:
    config = _load_config(config_path)
    cdx_out = config["cross_dataset_genomewide"]["output"]
    tables_dir = Path(cdx_out["wide_matrix_tsv"]).parent
    adj_tables = Path(config["candidate_adjudication"]["output"]["tables_dir"])

    pool = pd.read_csv(adj_tables / "adjudication_candidate_pool.tsv", sep="\t")
    axes = pd.read_csv(adj_tables / "three_axis_candidate_matrix.tsv", sep="\t")
    top_resistance = pd.read_csv(adj_tables / "top_resistance_genes_exact.tsv", sep="\t")
    top_crispr = pd.read_csv(adj_tables / "top_crispr_sensitising_exact.tsv", sep="\t")
    top_human = pd.read_csv(cdx_out["top20_human_tumor_tsv"], sep="\t")
    wide = pd.read_csv(tables_dir / "all_genes_cross_dataset_evidence_with_ranking.tsv", sep="\t")
    stability = pd.read_csv(cdx_out["ranking_stability_tsv"], sep="\t")

    list_a = build_list_a_multimodal_therapeutic(pool, axes, stability)
    list_b = build_list_b_resistance_biomarker(top_resistance)
    list_c = build_list_c_functional_sensitisation(top_crispr)
    list_d = build_list_d_human_tumor(top_human, wide)

    for name, df in [("A_multimodal_therapeutic", list_a), ("B_resistance_biomarker", list_b), ("C_functional_sensitisation", list_c), ("D_human_tumor", list_d)]:
        df.to_csv(adj_tables / f"shortlist_{name}.tsv", sep="\t", index=False)
    logger.info("run_shortlists: A=%d, B=%d, C=%d, D=%d genes", len(list_a), len(list_b), len(list_c), len(list_d))
    return {"A": list_a, "B": list_b, "C": list_c, "D": list_d}


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run_shortlists()
