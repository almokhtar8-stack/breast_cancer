"""Extracts the 13 sensitising CRISPR candidates (and, separately, the
PAICS benchmark) from the genome-wide GSE245601 pseudobulk edgeR results
(``scripts/analysis/gse245601_15_pseudobulk_edger.R``), for Track A
(epithelial) and Track B (strict malignant, n=3, exploratory).

BH multiple-testing correction is applied to exactly the candidate genes
present in a track's filtered DE table (genes removed by edgeR's
``filterByExpr`` -- too lowly expressed to test -- are excluded from the
BH family, not treated as non-significant). PAICS is a published
benchmark and is never included in that BH family
(``gse245601_PREANALYSIS.md``, `candidate_evidence_summary.py`'s own
convention). The genome-wide DE table is read for context but never used
to redefine which genes count as candidates.

Patient-direction consistency uses simple per-sample library-size CPM
(log2(CPM+1)) from the raw pseudobulk counts -- not edgeR's internal TMM
values, which are not exported per-sample-pair -- to answer only "how many
patients moved in the same direction," a supplementary, simpler question
than the edgeR model's aggregate coefficient.
"""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from src.gse245601_pseudobulk_qc import compute_log2cpm, load_pseudobulk

logger = logging.getLogger(__name__)

REQUIRED_DE_COLUMNS = ("gene", "log2fc", "avg_log_cpm", "p_value", "fdr")


def _load_config(config_path: str | Path) -> dict:
    with open(config_path) as f:
        return yaml.safe_load(f)


def load_de_table(path: str | Path) -> pd.DataFrame:
    """Load a genome-wide pseudobulk DE table. Read-only."""
    df = pd.read_csv(path, sep="\t")
    missing = [c for c in REQUIRED_DE_COLUMNS if c not in df.columns]
    if missing:
        raise ValueError(f"DE table missing required columns: {missing}")
    logger.info("load_de_table: %d genes tested (%s)", len(df), path)
    return df


def benjamini_hochberg(p_values: pd.Series) -> pd.Series:
    """Standard BH step-up FDR, applied to exactly the p-values passed in
    -- the caller controls the family (e.g. only the candidate genes
    actually present in a track's filtered DE table)."""
    p = p_values.to_numpy(dtype=float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order]
    bh = ranked * n / (np.arange(1, n + 1))
    bh = np.minimum.accumulate(bh[::-1])[::-1]
    bh = np.clip(bh, 0, 1)
    out = np.empty(n)
    out[order] = bh
    return pd.Series(out, index=p_values.index)


def compute_patient_direction(counts: pd.DataFrame, metadata: pd.DataFrame, gene: str) -> pd.DataFrame:
    """Per-patient sign of change (Tamoxifen minus Control) in simple
    log2(CPM+1), for one gene. One row per patient with both arms
    present."""
    sample_ids = metadata["sample_id"].tolist()
    log2cpm = compute_log2cpm(counts, sample_ids)
    if gene not in log2cpm.index:
        return pd.DataFrame(columns=["patient", "control_log2cpm", "tamoxifen_log2cpm", "delta", "direction"])

    gene_vals = log2cpm.loc[gene]
    meta_indexed = metadata.set_index("sample_id")
    rows = []
    for patient in sorted(metadata["patient"].unique()):
        patient_samples = meta_indexed.loc[meta_indexed["patient"] == patient]
        if not {"Control", "Tamoxifen"}.issubset(set(patient_samples["condition"])):
            continue
        control_id = patient_samples.loc[patient_samples["condition"] == "Control"].index[0]
        tam_id = patient_samples.loc[patient_samples["condition"] == "Tamoxifen"].index[0]
        control_val = gene_vals[control_id]
        tam_val = gene_vals[tam_id]
        delta = tam_val - control_val
        rows.append(
            {
                "patient": patient,
                "control_log2cpm": control_val,
                "tamoxifen_log2cpm": tam_val,
                "delta": delta,
                "direction": "up" if delta > 0 else ("down" if delta < 0 else "unchanged"),
            }
        )
    return pd.DataFrame(rows)


def build_candidate_table(
    de_df: pd.DataFrame, counts: pd.DataFrame, metadata: pd.DataFrame, candidate_genes: list[str], track_label: str
) -> pd.DataFrame:
    """One row per candidate gene present in the track's filtered DE
    table (genes removed by filterByExpr are reported as
    ``tested=False``, excluded from the candidate-set BH family, not
    silently dropped from the output). BH FDR is computed only across the
    tested candidates."""
    rows = []
    tested = de_df.loc[de_df["gene"].isin(candidate_genes)].copy()
    tested_genes = set(tested["gene"])
    untested_genes = [g for g in candidate_genes if g not in tested_genes]

    if len(tested) > 0:
        tested["candidate_set_bh_fdr"] = benjamini_hochberg(tested["p_value"])
    else:
        tested["candidate_set_bh_fdr"] = pd.Series(dtype=float)

    for gene in candidate_genes:
        if gene in untested_genes:
            rows.append(
                {
                    "gene": gene,
                    "track": track_label,
                    "tested": False,
                    "reason_not_tested": "filtered out by edgeR filterByExpr (too lowly expressed to test)",
                    "log2fc": np.nan,
                    "avg_log_cpm": np.nan,
                    "p_value": np.nan,
                    "genomewide_fdr": np.nan,
                    "candidate_set_bh_fdr": np.nan,
                    "direction": "not_tested",
                    "n_patients_up": 0,
                    "n_patients_down": 0,
                    "n_patients_tested": 0,
                    "direction_consistent_fraction": np.nan,
                }
            )
            continue

        row = tested.loc[tested["gene"] == gene].iloc[0]
        patient_dir = compute_patient_direction(counts, metadata, gene)
        n_up = int((patient_dir["direction"] == "up").sum())
        n_down = int((patient_dir["direction"] == "down").sum())
        n_tested_patients = len(patient_dir)
        majority = max(n_up, n_down)
        rows.append(
            {
                "gene": gene,
                "track": track_label,
                "tested": True,
                "reason_not_tested": "",
                "log2fc": float(row["log2fc"]),
                "avg_log_cpm": float(row["avg_log_cpm"]),
                "p_value": float(row["p_value"]),
                "genomewide_fdr": float(row["fdr"]),
                "candidate_set_bh_fdr": float(row["candidate_set_bh_fdr"]),
                "direction": "up" if row["log2fc"] > 0 else ("down" if row["log2fc"] < 0 else "unchanged"),
                "n_patients_up": n_up,
                "n_patients_down": n_down,
                "n_patients_tested": n_tested_patients,
                "direction_consistent_fraction": majority / n_tested_patients if n_tested_patients > 0 else np.nan,
            }
        )

    out = pd.DataFrame(rows)
    logger.info(
        "build_candidate_table[%s]: %d/%d candidates tested, %d untested (filtered)",
        track_label,
        int(out["tested"].sum()),
        len(candidate_genes),
        len(untested_genes),
    )
    return out


def build_paics_row(de_df: pd.DataFrame, counts: pd.DataFrame, metadata: pd.DataFrame, paics_gene: str, track_label: str) -> pd.DataFrame:
    """PAICS reported exactly like a candidate row but explicitly labeled
    as a benchmark and never included in any candidate-set BH family."""
    tested = de_df.loc[de_df["gene"] == paics_gene]
    if len(tested) == 0:
        return pd.DataFrame(
            [
                {
                    "gene": paics_gene,
                    "track": track_label,
                    "benchmark_label": "published_benchmark_not_gate1_hit_not_in_13_candidate_bh_family",
                    "tested": False,
                    "log2fc": np.nan,
                    "avg_log_cpm": np.nan,
                    "p_value": np.nan,
                    "genomewide_fdr": np.nan,
                    "direction": "not_tested",
                    "n_patients_up": 0,
                    "n_patients_down": 0,
                    "n_patients_tested": 0,
                }
            ]
        )
    row = tested.iloc[0]
    patient_dir = compute_patient_direction(counts, metadata, paics_gene)
    return pd.DataFrame(
        [
            {
                "gene": paics_gene,
                "track": track_label,
                "benchmark_label": "published_benchmark_not_gate1_hit_not_in_13_candidate_bh_family",
                "tested": True,
                "log2fc": float(row["log2fc"]),
                "avg_log_cpm": float(row["avg_log_cpm"]),
                "p_value": float(row["p_value"]),
                "genomewide_fdr": float(row["fdr"]),
                "direction": "up" if row["log2fc"] > 0 else ("down" if row["log2fc"] < 0 else "unchanged"),
                "n_patients_up": int((patient_dir["direction"] == "up").sum()),
                "n_patients_down": int((patient_dir["direction"] == "down").sum()),
                "n_patients_tested": len(patient_dir),
            }
        ]
    )


def run_candidate_extraction(config_path: str | Path = "config/config.yaml") -> dict[str, pd.DataFrame]:
    config = _load_config(config_path)
    pb_cfg = config["gse245601_pseudobulk"]
    ci_cfg = config["gse245601_candidate_integration"]
    candidates_13 = pb_cfg["candidates"]["thirteen"]
    paics_gene = pb_cfg["candidates"]["paics"]

    track_a_de = load_de_table(pb_cfg["output"]["de"]["track_a_genomewide_tsv"])
    track_b_de = load_de_table(pb_cfg["output"]["de"]["track_b_genomewide_tsv"])
    track_a_counts, track_a_meta = load_pseudobulk(pb_cfg["output"]["track_a"]["counts_tsv"], pb_cfg["output"]["track_a"]["metadata_tsv"])
    track_b_counts, track_b_meta = load_pseudobulk(pb_cfg["output"]["track_b"]["counts_tsv"], pb_cfg["output"]["track_b"]["metadata_tsv"])

    track_a_candidates = build_candidate_table(track_a_de, track_a_counts, track_a_meta, candidates_13, "track_a_epithelial")
    track_b_candidates = build_candidate_table(track_b_de, track_b_counts, track_b_meta, candidates_13, "track_b_malignant")

    paics_a = build_paics_row(track_a_de, track_a_counts, track_a_meta, paics_gene, "track_a_epithelial")
    paics_b = build_paics_row(track_b_de, track_b_counts, track_b_meta, paics_gene, "track_b_malignant")
    paics_out = pd.concat([paics_a, paics_b], ignore_index=True)

    Path(ci_cfg["output"]["track_a_candidates_tsv"]).parent.mkdir(parents=True, exist_ok=True)
    track_a_candidates.to_csv(ci_cfg["output"]["track_a_candidates_tsv"], sep="\t", index=False)
    track_b_candidates.to_csv(ci_cfg["output"]["track_b_candidates_tsv"], sep="\t", index=False)
    paics_out.to_csv(ci_cfg["output"]["paics_singlecell_tsv"], sep="\t", index=False)
    logger.info("run_candidate_extraction: wrote track A/B candidate tables and the PAICS benchmark table")

    return {"track_a_candidates": track_a_candidates, "track_b_candidates": track_b_candidates, "paics": paics_out}


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run_candidate_extraction()
