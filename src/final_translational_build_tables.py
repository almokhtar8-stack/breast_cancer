"""Final USP34/VEZF1 translational + structure phase: builds the seven
output tables from the curated data in src/final_translational_data.py.

USP34 (lead) and VEZF1 (second/backup) rankings are frozen and unchanged.
No CRISPR/RNA-seq/TCGA/DepMap/druggability computation happens here --
this phase is forward-looking experimental design plus one genuinely new
analysis (a real, locally-run structural/pocket analysis of the two
experimental USP34 structures, see final_translational_data.py's module
docstring for exact provenance).
"""

from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd

from src.final_translational_data import (
    DOCKING_DECISION_ROW,
    DOCKING_QUESTIONS_ROWS,
    EXPERIMENTAL_DESIGN_ROWS,
    FINAL_CONCLUSIONS_ROWS,
    NORMAL_CELL_COMPARATOR_ROWS,
    POCKET_ANALYSIS_ROWS,
    STRUCTURE_INVENTORY_ROWS,
    SUCCESS_FAILURE_ROWS,
    TARGETS,
)

logger = logging.getLogger(__name__)

OUT_DIR = Path("results/tables/final_translational")


def build_experimental_design_table() -> pd.DataFrame:
    df = pd.DataFrame(EXPERIMENTAL_DESIGN_ROWS)
    assert set(df["experiment_id"]) == {"EXP-1", "EXP-3", "EXP-5"}
    logger.info("build_experimental_design_table: %d rows", len(df))
    return df


def build_normal_cell_comparators_table() -> pd.DataFrame:
    df = pd.DataFrame(NORMAL_CELL_COMPARATOR_ROWS)
    assert set(df["experiment_id"]) == {"EXP-2A", "EXP-2B", "EXP-4"}
    logger.info("build_normal_cell_comparators_table: %d rows", len(df))
    return df


def build_structure_inventory_table() -> pd.DataFrame:
    df = pd.DataFrame(STRUCTURE_INVENTORY_ROWS)
    logger.info("build_structure_inventory_table: %d rows", len(df))
    return df


def build_pocket_analysis_table() -> pd.DataFrame:
    df = pd.DataFrame(POCKET_ANALYSIS_ROWS)
    logger.info("build_pocket_analysis_table: %d rows", len(df))
    return df


def build_docking_decision_table() -> pd.DataFrame:
    questions = pd.DataFrame(DOCKING_QUESTIONS_ROWS)
    decision_row = dict(question="DECISION", answer=DOCKING_DECISION_ROW["decision"])
    justification_row = dict(question="JUSTIFICATION", answer=DOCKING_DECISION_ROW["justification"])
    alternative_row = dict(question="ALTERNATIVE (structure-based roadmap)", answer=DOCKING_DECISION_ROW["alternative_recommended"])
    out = pd.concat([questions, pd.DataFrame([decision_row, justification_row, alternative_row])], ignore_index=True)
    logger.info("build_docking_decision_table: %d rows (6 questions + decision + justification + alternative)", len(out))
    return out


def build_success_failure_table() -> pd.DataFrame:
    df = pd.DataFrame(SUCCESS_FAILURE_ROWS)
    assert set(df["target"]) == set(TARGETS)
    logger.info("build_success_failure_table: %d rows", len(df))
    return df


def build_final_conclusions_table() -> pd.DataFrame:
    df = pd.DataFrame(FINAL_CONCLUSIONS_ROWS)
    assert set(df["target"]) == set(TARGETS)
    assert set(df["role"]) == {"LEAD TARGET", "SECOND / BACKUP TARGET"}
    logger.info("build_final_conclusions_table: %d rows", len(df))
    return df


def run(out_dir: Path = OUT_DIR) -> dict[str, pd.DataFrame]:
    out_dir.mkdir(parents=True, exist_ok=True)

    tables = {
        "final_experimental_design.tsv": build_experimental_design_table(),
        "final_normal_cell_comparators.tsv": build_normal_cell_comparators_table(),
        "USP34_structure_inventory.tsv": build_structure_inventory_table(),
        "USP34_pocket_analysis.tsv": build_pocket_analysis_table(),
        "USP34_docking_decision.tsv": build_docking_decision_table(),
        "final_target_success_failure_criteria.tsv": build_success_failure_table(),
        "final_translational_conclusions.tsv": build_final_conclusions_table(),
    }
    for name, df in tables.items():
        path = out_dir / name
        df.to_csv(path, sep="\t", index=False)
        logger.info("wrote %s (%d rows)", path, len(df))
    return tables


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run()
