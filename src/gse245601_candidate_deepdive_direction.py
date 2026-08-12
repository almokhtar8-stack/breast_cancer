"""GSE245601 candidate deep-dive Phase 5: per-gene patient direction
summary (how many of the 10 tumors' all-epithelial pseudobulk increase
vs decrease at 12h Tamoxifen) -- purely a count of the exact normalized
values already computed by `src.gse245601_candidate_deepdive_pseudobulk`,
no new statistic, no arbitrary "approximately flat" band (exact equality
only, which in practice never occurs with continuous normalized values).
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


def build_patient_direction_summary(all_epithelial_pseudobulk: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for gene, grp in all_epithelial_pseudobulk.groupby("gene"):
        valid = grp.dropna(subset=["direction_control_to_tam"])
        valid = valid.loc[valid["direction_control_to_tam"] != "not_comparable"]
        n_increase = int((valid["direction_control_to_tam"] == "increase").sum())
        n_decrease = int((valid["direction_control_to_tam"] == "decrease").sum())
        n_equal = int((valid["direction_control_to_tam"] == "equal").sum())
        rows.append(
            {
                "gene": gene, "n_patients_total": len(valid),
                "n_patients_increase": n_increase, "n_patients_decrease": n_decrease, "n_patients_equal": n_equal,
                "increasing_patients": ",".join(sorted(valid.loc[valid["direction_control_to_tam"] == "increase", "patient"])),
                "decreasing_patients": ",".join(sorted(valid.loc[valid["direction_control_to_tam"] == "decrease", "patient"])),
            }
        )
    out = pd.DataFrame(rows)
    logger.info("build_patient_direction_summary: %s", {r["gene"]: f"{r['n_patients_increase']}up/{r['n_patients_decrease']}down" for _, r in out.iterrows()})
    return out


def run_direction_summary(config_path: str | Path = "config/config.yaml") -> pd.DataFrame:
    config = _load_config(config_path)
    out_dir = Path(config["gse245601_candidate_deepdive"]["output"]["tables_dir"])
    all_epi = pd.read_csv(out_dir / "patient_all_epithelial_pseudobulk.tsv", sep="\t")
    summary = build_patient_direction_summary(all_epi)
    summary.to_csv(out_dir / "patient_direction_summary.tsv", sep="\t", index=False)
    return summary


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run_direction_summary()
