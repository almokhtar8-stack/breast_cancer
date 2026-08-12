"""GSE245601 candidate deep-dive Phases 15-16: malignant-vs-nonmalignant
baseline enrichment (descriptive, patient-paired) and a treatment x
malignancy interaction check.

Phase 16 explicitly requires NOT forcing an interaction model when sample
structure makes it unstable. Only 3/10 tumors (Tumor_02/03/07) have
reliable (>=50-cell) malignant pseudobulk in both arms; a
patient x malignancy x treatment count GLM would have ~6 residual
degrees of freedom on 12 pseudobulk observations for those 3 patients,
or would need to include the other 7 tumors' single-digit/near-zero
-cell malignant "pseudobulk" as if it were equally reliable, which it is
not. Per the task's own instruction, this is reported as NOT
statistically defensible, and the interaction is instead summarized
descriptively (a per-patient difference-of-differences, restricted to
the reliable patients, no p-value attached).
"""

from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd
import yaml

from src.gse245601_candidate_deepdive_data import GENES

logger = logging.getLogger(__name__)

RELIABLE_MALIGNANT_PATIENTS = ["Tumor_02", "Tumor_03", "Tumor_07"]


def _load_config(config_path: str | Path) -> dict:
    with open(config_path) as f:
        return yaml.safe_load(f)


def build_malignant_enrichment(malig_cond: pd.DataFrame, genes: list[str]) -> pd.DataFrame:
    rows = []
    sub = malig_cond.loc[malig_cond["malignancy_status"].isin(["malignant", "non-malignant epithelial"])]
    for gene in genes:
        for condition in ("Control", "Tamoxifen"):
            g = sub.loc[(sub["gene"] == gene) & (sub["condition"] == condition)]
            wide = g.pivot_table(index="patient", columns="malignancy_status", values="pseudobulk_normalized_expression")
            if "malignant" not in wide.columns or "non-malignant epithelial" not in wide.columns:
                continue
            valid = wide.dropna()
            delta = valid["malignant"] - valid["non-malignant epithelial"]
            rows.append(
                {
                    "gene": gene, "condition": condition, "n_patients_compared": len(valid),
                    "n_patients_malignant_higher": int((delta > 0).sum()), "n_patients_nonmalignant_higher": int((delta < 0).sum()),
                    "median_malignant_minus_nonmalignant": float(delta.median()) if len(delta) else float("nan"),
                    "reliable_patients_only_median": float(delta.reindex(RELIABLE_MALIGNANT_PATIENTS).dropna().median()) if delta.reindex(RELIABLE_MALIGNANT_PATIENTS).notna().any() else float("nan"),
                }
            )
    out = pd.DataFrame(rows)
    logger.info("build_malignant_enrichment: %s", {(r["gene"], r["condition"]): round(r["median_malignant_minus_nonmalignant"], 2) for _, r in out.iterrows()})
    return out


def build_interaction_summary(malig_cond: pd.DataFrame, genes: list[str]) -> pd.DataFrame:
    """Descriptive difference-of-differences per gene, restricted to the
    3 reliably-sampled malignant patients: (malignant Tam - malignant
    Ctrl) - (nonmalignant Tam - nonmalignant Ctrl). No p-value: reported
    as NOT statistically defensible (see module docstring)."""
    rows = []
    sub = malig_cond.loc[malig_cond["patient"].isin(RELIABLE_MALIGNANT_PATIENTS) & malig_cond["malignancy_status"].isin(["malignant", "non-malignant epithelial"])]
    for gene in genes:
        g = sub.loc[sub["gene"] == gene]
        wide = g.pivot_table(index="patient", columns=["malignancy_status", "condition"], values="pseudobulk_normalized_expression")
        per_patient = {}
        for patient in RELIABLE_MALIGNANT_PATIENTS:
            if patient not in wide.index:
                continue
            row = wide.loc[patient]
            try:
                mal_delta = row[("malignant", "Tamoxifen")] - row[("malignant", "Control")]
                nonmal_delta = row[("non-malignant epithelial", "Tamoxifen")] - row[("non-malignant epithelial", "Control")]
                per_patient[patient] = mal_delta - nonmal_delta
            except KeyError:
                continue
        rows.append(
            {
                "gene": gene, "n_reliable_patients": len(per_patient),
                "per_patient_interaction_diff": "; ".join(f"{p}={v:.3f}" for p, v in per_patient.items()),
                "median_interaction_diff": pd.Series(per_patient).median() if per_patient else float("nan"),
                "statistically_defensible": False,
                "reason": "only 3/10 tumors have reliable (>=50-cell) malignant pseudobulk in both arms -- a formal count GLM interaction term would have ~6 residual df across 12 pseudobulk observations; reported descriptively only, no p-value",
            }
        )
    out = pd.DataFrame(rows)
    logger.info("build_interaction_summary: descriptive only (not statistically defensible), %d genes", len(out))
    return out


def run_enrichment(config_path: str | Path = "config/config.yaml") -> dict[str, pd.DataFrame]:
    config = _load_config(config_path)
    out_dir = Path(config["gse245601_candidate_deepdive"]["output"]["tables_dir"])
    malig_cond = pd.read_csv(out_dir / "malignancy_condition_patient_summary.tsv", sep="\t")

    enrichment = build_malignant_enrichment(malig_cond, GENES)
    enrichment.to_csv(out_dir / "malignant_enrichment.tsv", sep="\t", index=False)

    interaction = build_interaction_summary(malig_cond, GENES)
    interaction.to_csv(out_dir / "treatment_malignancy_interaction.tsv", sep="\t", index=False)

    return {"enrichment": enrichment, "interaction": interaction}


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run_enrichment()
