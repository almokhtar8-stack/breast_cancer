"""Independent-validation DepMap release comparison: 24Q4 vs 26Q1.

RUNNABLE, COMPLETE: DepMap Public 26Q1's CRISPRGeneEffect.csv, Model.csv,
matched expression file, and CRISPRGeneDependency.csv (added in a
follow-up manual download) were all manually downloaded from the official
portal and verified (see results/reports/independent_validation/
DEPMAP_26Q1_ACCESS_STATUS.md for the access-attempt history and
/ibex/scratch/aljaroaa/tamoxifen-data/depmap/26Q1/PROVENANCE.txt for
verification detail). Probability-dependent fields (strong-dependency
fraction, A-E tier) are now genuine computed values for both releases.

Builds one row per candidate: 24Q4 and 26Q1 side by side (median gene
effect for all-cancer/breast/ER+-luminal, strong-dependency fraction for
all-cancer/breast/ER+-luminal, essentiality-concern tier), the delta
between them, and explicit classification_changed / interpretation_changed
flags -- so a reviewer can see at a glance whether 26Q1 materially changes
any conclusion rather than merely replacing a version number.
`classification_changed` now correctly reads False for all four
candidates: every essentiality-concern tier reproduces exactly between
releases (see the main report's Part 8B).
"""

from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd

from src.independent_validation_depmap import CANDIDATES, build_dependency_table
from src.independent_validation_depmap_data import load_config, load_model

logger = logging.getLogger(__name__)

OUT_TABLE = Path("results/tables/independent_validation/DepMap_24Q4_vs_26Q1_comparison.tsv")


def build_comparison_table(cfg: dict, old_release: str = "24Q4", new_release: str = "26Q1") -> pd.DataFrame:
    old = build_dependency_table(cfg, old_release).set_index("candidate")
    new = build_dependency_table(cfg, new_release).set_index("candidate")

    old_model = load_model(cfg, old_release)
    new_model = load_model(cfg, new_release)

    rows = []
    for candidate in CANDIDATES:
        o, n = old.loc[candidate], new.loc[candidate]
        classification_changed = o["essentiality_concern"] != n["essentiality_concern"]
        rows.append(dict(
            candidate=candidate,
            **{f"{old_release}_all_median": o["median_gene_effect_all_cancer"]},
            **{f"{new_release}_all_median": n["median_gene_effect_all_cancer"]},
            delta_all=n["median_gene_effect_all_cancer"] - o["median_gene_effect_all_cancer"],
            **{f"{old_release}_breast_median": o["median_gene_effect_breast"]},
            **{f"{new_release}_breast_median": n["median_gene_effect_breast"]},
            delta_breast=n["median_gene_effect_breast"] - o["median_gene_effect_breast"],
            **{f"{old_release}_ERluminal_median": o["median_gene_effect_er_luminal"]},
            **{f"{new_release}_ERluminal_median": n["median_gene_effect_er_luminal"]},
            delta_ERluminal=n["median_gene_effect_er_luminal"] - o["median_gene_effect_er_luminal"],
            **{f"{old_release}_strong_dependency_fraction": o["frac_strongly_dependent_all_cancer"]},
            **{f"{new_release}_strong_dependency_fraction": n["frac_strongly_dependent_all_cancer"]},
            **{f"{old_release}_strong_dependency_fraction_breast": o["frac_strongly_dependent_breast"]},
            **{f"{new_release}_strong_dependency_fraction_breast": n["frac_strongly_dependent_breast"]},
            **{f"{old_release}_strong_dependency_fraction_ERluminal": o["frac_strongly_dependent_er_luminal"]},
            **{f"{new_release}_strong_dependency_fraction_ERluminal": n["frac_strongly_dependent_er_luminal"]},
            **{f"{old_release}_classification": o["essentiality_concern"]},
            **{f"{new_release}_classification": n["essentiality_concern"]},
            classification_changed=classification_changed,
            interpretation_changed=classification_changed,  # tracks the A-E tier; a tier match does not by itself rule out a narrative nuance (e.g. VEZF1's ER+/luminal % dropped 36%->27% while staying MODERATE) -- see report Part 8B for the qualitative discussion
        ))
    out = pd.DataFrame(rows)

    n_old_all, n_new_all = len(old_model), len(new_model)
    n_old_breast, n_new_breast = int(old_model["is_breast"].sum()), int(new_model["is_breast"].sum())
    n_old_luminal, n_new_luminal = int(old_model["is_er_luminal"].sum()), int(new_model["is_er_luminal"].sum())
    logger.info(
        "build_comparison_table: models %s=%d -> %s=%d; breast %d -> %d; ER+/luminal %d -> %d",
        old_release, n_old_all, new_release, n_new_all, n_old_breast, n_new_breast, n_old_luminal, n_new_luminal,
    )
    return out


def run(config_path: str = "config/config.yaml", out_table: Path = OUT_TABLE, old_release: str = "24Q4", new_release: str = "26Q1") -> pd.DataFrame:
    cfg = load_config(config_path)
    out = build_comparison_table(cfg, old_release, new_release)
    out_table.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(out_table, sep="\t", index=False)
    logger.info("wrote %s (%d rows)", out_table, len(out))
    return out


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run()
