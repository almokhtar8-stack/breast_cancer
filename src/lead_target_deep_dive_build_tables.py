"""USP34 vs VEZF1 translational deep-dive: builds the thirteen output tables
for this phase from the curated data in src/lead_target_deep_dive_data.py,
combined with the carried-forward, already-verified reference list from
the prior druggability_safety phase (filtered to USP34/VEZF1 rows only)
and a live join against this project's own frozen CRISPR/RNA-seq
resistance-evidence table for the indirect-target cross-check.

USP34 and VEZF1 are the two lead candidates selected for deeper follow-up;
CITED2 and EML5 remain frozen shortlist members, not disproven. No
CRISPR/RNA-seq/TCGA/DepMap computation happens here -- all such values are
read directly from already-frozen output files, never retyped.
"""

from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd

from src.druggability_safety_data import REFERENCES_ROWS as PRIOR_PHASE_REFERENCES_ROWS
from src.lead_target_deep_dive_data import (
    BONE_MARROW_LIABILITY_ROWS,
    EXPERIMENTAL_PLAN_ROWS,
    GENETIC_CONSTRAINT_ROWS,
    HEAD_TO_HEAD_ROWS,
    INDIRECT_TARGET_CROSSCHECK_ROWS,
    MUSCLE_LIABILITY_ROWS,
    NEW_REFERENCES_ROWS,
    TISSUE_EXPRESSION_ROWS,
    TISSUE_LIABILITY_ROWS,
    USP34_DIRECT_TARGETING_ROWS,
    USP34_INDIRECT_TARGET_ROWS,
    VEZF1_DIRECT_TARGETING_ROWS,
    VEZF1_INDIRECT_TARGET_ROWS,
)

logger = logging.getLogger(__name__)

OUT_DIR = Path("results/tables/lead_target_deep_dive")
CROSS_DATASET_TABLE = Path("results/tables/cross_dataset_genomewide/all_genes_cross_dataset_evidence_with_ranking.tsv")


def build_tissue_expression_table() -> pd.DataFrame:
    df = pd.DataFrame(TISSUE_EXPRESSION_ROWS)
    assert set(df["candidate"]) == {"USP34", "VEZF1"}
    logger.info("build_tissue_expression_table: %d rows", len(df))
    return df


def build_tissue_liability_table() -> pd.DataFrame:
    df = pd.DataFrame(TISSUE_LIABILITY_ROWS)
    assert set(df["candidate"]) == {"USP34", "VEZF1"}
    logger.info("build_tissue_liability_table: %d rows", len(df))
    return df


def build_muscle_liability_table() -> pd.DataFrame:
    df = pd.DataFrame(MUSCLE_LIABILITY_ROWS)
    logger.info("build_muscle_liability_table: %d rows", len(df))
    return df


def build_bone_marrow_liability_table() -> pd.DataFrame:
    df = pd.DataFrame(BONE_MARROW_LIABILITY_ROWS)
    assert set(df["candidate"]) == {"USP34", "VEZF1"}
    logger.info("build_bone_marrow_liability_table: %d rows", len(df))
    return df


def build_genetic_constraint_table() -> pd.DataFrame:
    df = pd.DataFrame(GENETIC_CONSTRAINT_ROWS)
    assert set(df["candidate"]) == {"USP34", "VEZF1"}
    logger.info("build_genetic_constraint_table: %d rows", len(df))
    return df


def build_usp34_direct_targeting_table() -> pd.DataFrame:
    df = pd.DataFrame(USP34_DIRECT_TARGETING_ROWS)
    logger.info("build_usp34_direct_targeting_table: %d rows", len(df))
    return df


def build_vezf1_direct_targeting_table() -> pd.DataFrame:
    df = pd.DataFrame(VEZF1_DIRECT_TARGETING_ROWS)
    logger.info("build_vezf1_direct_targeting_table: %d rows", len(df))
    return df


def build_usp34_indirect_targets_table() -> pd.DataFrame:
    df = pd.DataFrame(USP34_INDIRECT_TARGET_ROWS)
    logger.info("build_usp34_indirect_targets_table: %d rows", len(df))
    return df


def build_vezf1_indirect_targets_table() -> pd.DataFrame:
    df = pd.DataFrame(VEZF1_INDIRECT_TARGET_ROWS)
    logger.info("build_vezf1_indirect_targets_table: %d rows", len(df))
    return df


def build_indirect_target_crosscheck_table() -> pd.DataFrame:
    """Live join: re-reads this project's own frozen cross-dataset evidence
    table for each candidate indirect-target gene. ALL numeric evidence
    columns (hany_crispr_effect/fdr/direction, gse118713_fdr,
    gse240112_tumor_fdr) are overwritten from this live read on every
    build -- there is no separate hardcoded copy of these numbers
    anywhere, so they cannot drift out of sync with the source file (a
    prior version of this function appended live_*_check columns
    alongside untouched hardcoded canonical columns, which could show
    stale and live values side by side if the source ever changed; fixed
    per Codex review).
    """
    cross = pd.read_csv(CROSS_DATASET_TABLE, sep="\t").set_index("gene")
    rows = []
    for row in INDIRECT_TARGET_CROSSCHECK_ROWS:
        gene = row["gene_x"]
        out_row = dict(row)
        if gene in cross.index:
            live = cross.loc[gene]
            out_row["hany_crispr_effect"] = float(live["crispr_effect"])
            out_row["hany_crispr_fdr"] = float(live["crispr_fdr"])
            out_row["hany_crispr_direction"] = live["crispr_direction"]
            out_row["gse118713_fdr"] = float(live["gse118713_fdr"]) if pd.notna(live["gse118713_fdr"]) else None
            out_row["gse240112_tumor_fdr"] = float(live["gse240112_tumor_fdr"]) if pd.notna(live["gse240112_tumor_fdr"]) else None
        else:
            out_row["hany_crispr_effect"] = None
            out_row["hany_crispr_fdr"] = None
            out_row["hany_crispr_direction"] = None
            out_row["gse118713_fdr"] = None
            out_row["gse240112_tumor_fdr"] = None
        rows.append(out_row)
    out = pd.DataFrame(rows)
    logger.info("build_indirect_target_crosscheck_table: %d rows, all numeric columns live-read from %s", len(out), CROSS_DATASET_TABLE)
    return out


def build_head_to_head_table() -> pd.DataFrame:
    df = pd.DataFrame(HEAD_TO_HEAD_ROWS)
    logger.info("build_head_to_head_table: %d rows (dimensions)", len(df))
    return df


def build_experimental_plan_table() -> pd.DataFrame:
    df = pd.DataFrame(EXPERIMENTAL_PLAN_ROWS)
    logger.info("build_experimental_plan_table: %d rows", len(df))
    return df


def build_verified_references_table() -> pd.DataFrame:
    prior = pd.DataFrame(PRIOR_PHASE_REFERENCES_ROWS)
    prior = prior[prior["candidate"].isin(["USP34", "VEZF1"])].copy()
    prior["topic"] = prior["candidate"]
    prior["verification_note"] = "carried forward from druggability_safety phase: " + prior["verification_note"]
    prior_out = prior[["topic", "PMID", "title", "journal", "verification_note"]]

    new = pd.DataFrame(NEW_REFERENCES_ROWS)
    new_out = new[["topic", "PMID", "title", "journal", "verification_note"]]

    out = pd.concat([prior_out, new_out], ignore_index=True)
    logger.info("build_verified_references_table: %d rows (%d carried forward, %d new this pass)", len(out), len(prior_out), len(new_out))
    return out


def run(out_dir: Path = OUT_DIR) -> dict[str, pd.DataFrame]:
    out_dir.mkdir(parents=True, exist_ok=True)

    tables = {
        "USP34_VEZF1_full_tissue_expression.tsv": build_tissue_expression_table(),
        "USP34_VEZF1_tissue_liability.tsv": build_tissue_liability_table(),
        "USP34_VEZF1_muscle_liability.tsv": build_muscle_liability_table(),
        "USP34_VEZF1_bone_marrow_liability.tsv": build_bone_marrow_liability_table(),
        "USP34_VEZF1_human_genetic_constraint.tsv": build_genetic_constraint_table(),
        "USP34_direct_targeting.tsv": build_usp34_direct_targeting_table(),
        "USP34_indirect_targets.tsv": build_usp34_indirect_targets_table(),
        "VEZF1_direct_targeting.tsv": build_vezf1_direct_targeting_table(),
        "VEZF1_indirect_targets.tsv": build_vezf1_indirect_targets_table(),
        "indirect_target_project_crosscheck.tsv": build_indirect_target_crosscheck_table(),
        "USP34_VEZF1_head_to_head.tsv": build_head_to_head_table(),
        "USP34_VEZF1_experimental_plan.tsv": build_experimental_plan_table(),
        "verified_references.tsv": build_verified_references_table(),
    }
    for name, df in tables.items():
        path = out_dir / name
        df.to_csv(path, sep="\t", index=False)
        logger.info("wrote %s (%d rows)", path, len(df))
    return tables


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run()
