"""Integrates the frozen CRISPR + bulk-RNA candidate evidence
(``results/tables/candidate_evidence_summary.tsv``, produced by
``src.candidate_evidence_summary`` -- NOT modified here) with the new
GSE245601 single-cell pseudobulk evidence (Track A, Track B, and the
malignant-vs-nonmalignant context table) into one new, additive table.

No composite/weighted score is computed anywhere in this module, and no
gene is forced to be "positive" across every evidence layer -- this
module only joins already-computed effect sizes, p-values, and FDRs into
one row per candidate, namespaced by layer, plus a small number of
transparent boolean evidence flags (e.g. "does the sign agree with
CRISPR's sensitising direction") that never combine into a single number.
"""

from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd
import yaml

logger = logging.getLogger(__name__)


def _load_config(config_path: str | Path) -> dict:
    with open(config_path) as f:
        return yaml.safe_load(f)


def load_crispr_bulk_evidence(path: str | Path, candidate_genes: list[str]) -> pd.DataFrame:
    """The frozen 28-hit evidence summary, restricted to the 13
    prospectively-tested candidates. Read-only -- this table is never
    written to by this module."""
    df = pd.read_csv(path, sep="\t")
    sub = df.loc[df["gene_symbol"].isin(candidate_genes)].copy()
    missing = set(candidate_genes) - set(sub["gene_symbol"])
    if missing:
        raise ValueError(f"candidate gene(s) missing from the frozen CRISPR/bulk evidence table: {missing}")
    if sub["gene_symbol"].duplicated().any():
        raise ValueError("duplicate gene_symbol rows in the frozen CRISPR/bulk evidence table")
    logger.info("load_crispr_bulk_evidence: %d/%d candidates found", len(sub), len(candidate_genes))
    return sub


def _namespace_sc_track(track_df: pd.DataFrame, prefix: str) -> pd.DataFrame:
    out = track_df[["gene", "tested", "log2fc", "p_value", "candidate_set_bh_fdr", "direction", "n_patients_up", "n_patients_down"]].copy()
    out.columns = ["gene_symbol"] + [f"{prefix}_{c}" for c in out.columns[1:]]
    return out


def build_integrated_table(
    crispr_bulk: pd.DataFrame,
    track_a: pd.DataFrame,
    track_b: pd.DataFrame,
    malignant_context: pd.DataFrame,
    candidate_genes: list[str],
    track_b_n_pairs: int,
) -> pd.DataFrame:
    """One row per candidate, exactly the 13, in the order given. Every
    evidence layer is namespaced (crispr_, bulk_, sc_track_a_,
    sc_track_b_, malignant_context_) and preserved as-is -- no composite
    score, no forced cross-layer agreement."""
    out = pd.DataFrame({"gene_symbol": candidate_genes})

    crispr_cols = crispr_bulk[["gene_symbol", "crispr_effect_size", "crispr_fdr", "crispr_direction"]].rename(
        columns={"crispr_effect_size": "crispr_effect_size", "crispr_fdr": "crispr_fdr", "crispr_direction": "crispr_direction"}
    )
    crispr_cols["crispr_is_sensitising"] = crispr_cols["crispr_direction"] == "sensitising_knockout"
    out = out.merge(crispr_cols.drop(columns="crispr_direction"), on="gene_symbol", how="left", validate="one_to_one")

    bulk_cols = crispr_bulk[["gene_symbol", "tamr_vs_mcf7_log2fc", "tamr_vs_mcf7_fdr", "evidence_class"]].rename(
        columns={"tamr_vs_mcf7_log2fc": "bulk_tamr_vs_mcf7_log2fc", "tamr_vs_mcf7_fdr": "bulk_tamr_vs_mcf7_fdr", "evidence_class": "bulk_evidence_class"}
    )
    out = out.merge(bulk_cols, on="gene_symbol", how="left", validate="one_to_one")

    out = out.merge(_namespace_sc_track(track_a, "sc_track_a"), on="gene_symbol", how="left", validate="one_to_one")
    out = out.merge(_namespace_sc_track(track_b, "sc_track_b"), on="gene_symbol", how="left", validate="one_to_one")
    out["sc_track_b_exploratory_n3"] = True
    out["sc_track_b_n_patient_pairs"] = track_b_n_pairs

    mal_ctx = malignant_context.loc[~malignant_context["is_paics_benchmark"], ["gene", "tested", "mean_delta_malignant_minus_nonmalignant", "p_value", "candidate_set_bh_fdr", "direction", "n_patients"]].copy()
    mal_ctx.columns = [
        "gene_symbol",
        "malignant_context_tested",
        "malignant_context_mean_delta",
        "malignant_context_p_value",
        "malignant_context_bh_fdr",
        "malignant_context_direction",
        "malignant_context_n_patients",
    ]
    out = out.merge(mal_ctx, on="gene_symbol", how="left", validate="one_to_one")

    # Transparent, non-combined evidence flags -- descriptive navigation aids
    # only, never summed into a score and never treated as validation,
    # replication, or statistical evidence on their own. Missingness
    # propagates as pd.NA (a gene not tested in a required layer yields
    # "not evaluable", never a silent False that would look like a checked
    # disagreement.
    def _direction_agrees(tested_a: bool, direction_a: object, tested_b: bool, direction_b: object) -> object:
        if not tested_a or not tested_b or pd.isna(direction_a) or pd.isna(direction_b):
            return pd.NA
        return direction_a == direction_b

    bulk_direction = out["bulk_tamr_vs_mcf7_log2fc"].apply(
        lambda x: pd.NA if pd.isna(x) else ("up" if x > 0 else ("down" if x < 0 else "unchanged"))
    )
    out["flag_crispr_sensitising_and_sc_track_a_same_direction_as_bulk"] = [
        (_direction_agrees(tested_a, dir_a, True, dir_b) if is_sens else pd.NA)
        for is_sens, tested_a, dir_a, dir_b in zip(out["crispr_is_sensitising"], out["sc_track_a_tested"], out["sc_track_a_direction"], bulk_direction)
    ]
    out["flag_sc_track_a_and_track_b_same_direction"] = [
        _direction_agrees(tested_a, dir_a, tested_b, dir_b)
        for tested_a, dir_a, tested_b, dir_b in zip(out["sc_track_a_tested"], out["sc_track_a_direction"], out["sc_track_b_tested"], out["sc_track_b_direction"])
    ]

    logger.info("build_integrated_table: %d candidate rows (expected %d)", len(out), len(candidate_genes))
    return out


def run_candidate_integration(config_path: str | Path = "config/config.yaml") -> pd.DataFrame:
    config = _load_config(config_path)
    pb_cfg = config["gse245601_pseudobulk"]
    ci_cfg = config["gse245601_candidate_integration"]
    candidates_13 = pb_cfg["candidates"]["thirteen"]
    track_b_n_pairs = len(pb_cfg["track_b_eligible_patients"])

    crispr_bulk = load_crispr_bulk_evidence(ci_cfg["inputs"]["evidence_summary_tsv"], candidates_13)
    track_a = pd.read_csv(ci_cfg["output"]["track_a_candidates_tsv"], sep="\t")
    track_b = pd.read_csv(ci_cfg["output"]["track_b_candidates_tsv"], sep="\t")
    malignant_context = pd.read_csv(ci_cfg["output"]["malignant_vs_nonmalignant_candidates_tsv"], sep="\t")

    out = build_integrated_table(crispr_bulk, track_a, track_b, malignant_context, candidates_13, track_b_n_pairs)

    out_path = Path(ci_cfg["output"]["integrated_master_tsv"])
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(out_path, sep="\t", index=False)
    logger.info("run_candidate_integration: wrote %s (%d rows)", out_path, len(out))
    return out


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run_candidate_integration()
