"""GSE245601 candidate deep-dive Phases 21-22: final USP34 classification
and the four-gene acute-decomposition summary table. Synthesizes (never
recomputes) every table already built by the other modules in this
package.
"""

from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd
import yaml

from src.gse245601_candidate_deepdive_data import GENES

logger = logging.getLogger(__name__)

USP34_CLASSIFICATION = {
    "gene": "USP34",
    "classification": "H_NO_CLEAR_ACUTE_CHANGE; secondary note: weak, non-significant C_MALIGNANT_SPECIFIC_DECREASE trend",
    "justification": (
        "Frozen Track A (all epithelial): log2FC=-0.033, FDR=0.901 -- provides no evidence of an overall acute "
        "change (a nonsignificant result, not a formal equivalence test). Patient-level direction is an exact "
        "5-up/5-down split (Phase 5), not a majority in either "
        "direction. No single prevalence or intensity metric dominates (Phase 9: both vary heterogeneously "
        "across patients, 6up/4down and 3up/7down respectively, no consistent driver). All 5 sufficiently"
        "-represented epithelial clusters (Phase 13, >=3 tumors with >=10 cells/arm) show a small negative "
        "median descriptive log2FC (-0.008 to -0.250), but within every one of those clusters individual "
        "tumor direction remains frequently mixed (only 3-4 of 7-10 tumors increase per cluster) -- a weak, "
        "consistent central tendency without a single dominant subpopulation driving it. USP34-high treated "
        "cells (Phase 18) are spread across all 10 patients and most clusters, and are NOT disproportionately "
        "malignant (5.8% of high cells vs. 7.7% baseline malignant fraction among Tam cells) -- not "
        "concentrated in one tumor or compartment. Composition-shift correlation is weak (rho~0.24, Phase 20) "
        "-- not the dominant explanation either. "
        "The one partially-consistent secondary signal: among the 3 tumors with reliable (>=50-cell) malignant "
        "pseudobulk (Tumor_02/03/07), all three show a small decrease (Phase 6/8), consistent with the frozen "
        "Track B point estimate (log2FC=-0.181) -- but Track B's own FDR (0.703) means this is not "
        "statistically significant either, and n=3 tumors cannot establish a robust malignant-specific effect."
    ),
}


def _load_config(config_path: str | Path) -> dict:
    with open(config_path) as f:
        return yaml.safe_load(f)


def build_candidate_summary(tables_dir: Path, frozen: pd.DataFrame) -> pd.DataFrame:
    direction = pd.read_csv(tables_dir / "patient_direction_summary.tsv", sep="\t").set_index("gene")
    enrichment = pd.read_csv(tables_dir / "malignant_enrichment.tsv", sep="\t")
    cluster_response = pd.read_csv(tables_dir / "cluster_candidate_response.tsv", sep="\t")
    prevalence = pd.read_csv(tables_dir / "expression_prevalence_intensity.tsv", sep="\t")

    rows = []
    for gene in GENES:
        f = frozen.set_index("gene").loc[gene]
        d = direction.loc[gene]

        mal_tam = enrichment.loc[(enrichment["gene"] == gene) & (enrichment["condition"] == "Tamoxifen")]
        mal_dir = "malignant_higher" if len(mal_tam) and mal_tam["median_malignant_minus_nonmalignant"].iloc[0] > 0 else "nonmalignant_higher_or_similar"

        supported_clusters = cluster_response.loc[(cluster_response["gene"] == gene) & (cluster_response["sufficiently_represented"])]
        n_cluster_up = int((supported_clusters["median_descriptive_log2fc"] > 0).sum())
        n_cluster_down = int((supported_clusters["median_descriptive_log2fc"] < 0).sum())
        cluster_heterogeneity = f"{n_cluster_up} up / {n_cluster_down} down of {len(supported_clusters)} well-supported clusters"

        prev = prevalence.loc[(prevalence["gene"] == gene) & (prevalence["malignancy_status"] == "all_epithelial")]
        prev_wide = prev.pivot_table(index="patient", columns="condition", values=["fraction_expressing", "mean_normalized_positive_cells_only"])
        prev_delta = (prev_wide[("fraction_expressing", "Tamoxifen")] - prev_wide[("fraction_expressing", "Control")]).dropna()
        int_delta = (prev_wide[("mean_normalized_positive_cells_only", "Tamoxifen")] - prev_wide[("mean_normalized_positive_cells_only", "Control")]).dropna()
        prev_summary = f"{int((prev_delta > 0).sum())}up/{int((prev_delta < 0).sum())}down (prevalence)"
        int_summary = f"{int((int_delta > 0).sum())}up/{int((int_delta < 0).sum())}down (intensity)"

        rows.append(
            {
                "gene": gene,
                "overall_acute_direction_track_a": "down" if f["gse245601_epi_log2fc"] < 0 else "up",
                "track_a_log2fc": f["gse245601_epi_log2fc"], "track_a_fdr": f["gse245601_epi_fdr"],
                "track_b_log2fc": f["gse245601_malignant_log2fc"], "track_b_fdr": f["gse245601_malignant_fdr"],
                "patient_consistency": f"{int(d['n_patients_increase'])}up/{int(d['n_patients_decrease'])}down of {int(d['n_patients_total'])}",
                "malignant_vs_nonmalignant_tam": mal_dir,
                "fraction_expressing_response": prev_summary,
                "positive_cell_intensity_response": int_summary,
                "cluster_heterogeneity": cluster_heterogeneity,
                "main_caveat": (
                    "extremely sparse (~1% of cells detectable) -- interpret with caution" if gene == "EML5"
                    else "malignant compartment reliably sampled in only 3/10 tumors" if gene in ("USP34", "CITED2", "VEZF1")
                    else ""
                ),
            }
        )
    out = pd.DataFrame(rows)
    logger.info("build_candidate_summary: %d genes", len(out))
    return out


def run_summary(config_path: str | Path = "config/config.yaml") -> dict[str, pd.DataFrame]:
    config = _load_config(config_path)
    out_dir = Path(config["gse245601_candidate_deepdive"]["output"]["tables_dir"])
    frozen = pd.read_csv("results/tables/evidence_freeze/final_candidate_evidence.tsv", sep="\t")

    summary = build_candidate_summary(out_dir, frozen)
    summary.to_csv(out_dir / "candidate_acute_decomposition_summary.tsv", sep="\t", index=False)

    usp34_class = pd.DataFrame([USP34_CLASSIFICATION])
    usp34_class.to_csv(out_dir / "usp34_final_classification.tsv", sep="\t", index=False)

    return {"summary": summary, "usp34_classification": usp34_class}


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run_summary()
