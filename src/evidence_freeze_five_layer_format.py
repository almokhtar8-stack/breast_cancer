"""Evidence freeze: the canonical five-layer / four-RNA-arrow display
format shared by every table and figure in this phase.

Canonical RNA display order (never reordered): GSE118713 | GSE240112 |
GSE111151 || GSE245601. The first three are resistance/recurrence-state
context datasets and are the only three ever summed into a "resistance"
count; GSE245601 is the acute 12h ex vivo tamoxifen response and is
always shown as the fourth, visually-separated arrow -- present in every
display, never folded into the resistance consensus.

CRISPR direction (sensitising_KO / tolerance_associated_KO) and RNA
direction (up/down) are never merged into one symbol or one column: a
gene may legitimately show RNA up in all three resistance datasets while
CRISPR shows tolerance-associated KO -- that is not a contradiction, it
is two different biological questions (perturbation-fitness phenotype vs.
expression association) answered by two different evidence types.
"""

from __future__ import annotations

import logging

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

RESISTANCE_DATASET_ORDER = ["gse118713", "gse240112", "gse111151"]
FULL_RNA_DATASET_ORDER = ["gse118713", "gse240112", "gse111151", "gse245601"]

UP_ARROW = "↑"
DOWN_ARROW = "↓"
FLAT_ARROW = "→"
NA_SYMBOL = "NA"
SIG_MARK = "*"


def arrow_for(log2fc: float, fdr: float) -> str:
    """One arrow, with a trailing '*' iff FDR<0.05. NA (not the digit 0
    and not an unmarked arrow) for an untestable gene -- a missing value
    must never render as if it were a measured null result."""
    if pd.isna(log2fc):
        return NA_SYMBOL
    if log2fc > 0:
        base = UP_ARROW
    elif log2fc < 0:
        base = DOWN_ARROW
    else:
        base = FLAT_ARROW
    if pd.notna(fdr) and fdr < 0.05:
        return base + SIG_MARK
    return base


def resistance_pattern_3(row: pd.Series) -> str:
    """Positions 1-3 ONLY (GSE118713, GSE240112, GSE111151) -- GSE245601
    is never included here; this string is also the one used for any
    resistance-direction calculation."""
    arrows = [arrow_for(row["gse118713_log2fc"], row["gse118713_fdr"]), arrow_for(row["gse240112_log2fc"], row["gse240112_fdr"]), arrow_for(row["gse111151_log2fc"], row["gse111151_fdr"])]
    return " | ".join(arrows)


def acute_direction(row: pd.Series, acute_log2fc_col: str = "gse245601_epi_log2fc", acute_fdr_col: str = "gse245601_epi_fdr") -> str:
    """GSE245601's single summary arrow (Track A / all-epithelial by
    convention, matching the frozen cross-dataset ranking's own
    "representative track" choice) -- always computed and always shown,
    never included in any resistance count."""
    return arrow_for(row[acute_log2fc_col], row[acute_fdr_col])


def full_rna_pattern_4(row: pd.Series, acute_log2fc_col: str = "gse245601_epi_log2fc", acute_fdr_col: str = "gse245601_epi_fdr") -> str:
    """The full, mandatory four-arrow display: resistance_pattern_3,
    a visual '||' divider, then the acute arrow. The divider is not
    cosmetic -- it is the only thing in a plain-text string that marks
    the biological-context boundary between resistance/recurrence state
    and acute 12h response; every consumer of this string (tables,
    figures, docs) must preserve it, never silently drop it."""
    return f"{resistance_pattern_3(row)} || {acute_direction(row, acute_log2fc_col, acute_fdr_col)}"


def resistance_fdr05_count(row: pd.Series) -> int:
    return int(sum(pd.notna(row[f"{d}_fdr"]) and row[f"{d}_fdr"] < 0.05 for d in RESISTANCE_DATASET_ORDER))


def resistance_nominal_p05_count(row: pd.Series) -> int:
    return int(sum(pd.notna(row[f"{d}_p"]) and row[f"{d}_p"] < 0.05 for d in RESISTANCE_DATASET_ORDER))


def resistance_direction_consistency(row: pd.Series) -> str:
    """Direction agreement across the 3 resistance datasets ONLY
    (GSE245601 never enters this calculation) -- same vocabulary as the
    frozen `resistance_direction_consensus` column
    (all_up/all_down/majority_up/majority_down/mixed/insufficient), so
    this is a display-format re-derivation check, not a new statistic."""
    dirs = []
    for d in RESISTANCE_DATASET_ORDER:
        v = row[f"{d}_log2fc"]
        if pd.isna(v):
            continue
        dirs.append("up" if v > 0 else ("down" if v < 0 else None))
    dirs = [d for d in dirs if d is not None]
    if len(dirs) == 0:
        return "insufficient"
    n_up, n_down = dirs.count("up"), dirs.count("down")
    if n_up == len(dirs):
        return "all_up"
    if n_down == len(dirs):
        return "all_down"
    return "majority_up" if n_up > n_down else ("majority_down" if n_down > n_up else "mixed")


LEGEND = {
    UP_ARROW: "higher in resistant/recurrent (or, for GSE245601, tamoxifen-treated) vs. parental/primary/control",
    DOWN_ARROW: "lower in resistant/recurrent (or GSE245601 tamoxifen-treated) vs. parental/primary/control",
    FLAT_ARROW: "log2FC exactly 0",
    NA_SYMBOL: "gene not testable in this dataset (filtered out, not a measured null result)",
    SIG_MARK: "FDR<0.05 in that dataset (appended to the arrow, e.g. '↑*')",
    "||": "divider between resistance/recurrence-context RNA (first 3) and acute 12h tamoxifen response (4th) -- never the same biological state",
}
