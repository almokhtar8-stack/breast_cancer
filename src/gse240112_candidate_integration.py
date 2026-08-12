"""GSE240112 Phase 18: additive 4-layer evidence integration table for the
frozen 13 candidates -- CRISPR screen, GSE118713 bulk acquired-resistance
RNA, GSE245601 acute (12h ex vivo) scRNA, and GSE240112 primary-vs-
recurrent (tamoxifen-resistance context) scRNA. This is a NEW table; it
never overwrites any existing frozen master table
(results/tables/candidate_evidence_summary.tsv or
results/tables/gse245601_candidate_integration/integrated_crispr_bulk_singlecell_candidates.tsv).

No weighted composite score is computed -- the four datasets measure
different biology (CRISPR = functional perturbation under 4-OHT;
GSE118713 = acquired-resistance cell-line state; GSE245601 = acute 12h ex
vivo treatment response; GSE240112 = human primary-vs-recurrent/
treatment-history resistance context) and are not assumed to agree in
direction. The "notes" column is transparent text only.

Data sources: Hany et al. CRISPR screen (root PREANALYSIS.md); GSE118713
(Bhat-Nakshatri lab bulk RNA-seq); GSE245601 (Kim/Whitman et al., Clin
Cancer Res 2023, PMID 37747807); GSE240112 (Fang et al., Genome Medicine
2024, PMID 39558215). Versions as already frozen/committed in this repo
plus the GSE240112 tables produced 2026-08-12.
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


def build_integrated_table(
    crispr_bulk: pd.DataFrame,
    track_a: pd.DataFrame,
    track_b: pd.DataFrame,
    gse240112_candidates: pd.DataFrame,
    gse240112_sensitivity: pd.DataFrame,
    candidate_genes: list[str],
) -> pd.DataFrame:
    """One row per frozen candidate, exactly once. Every source table is
    joined by gene; a gene absent from a given source is reported with
    NA in that source's columns, never dropped from the output."""
    crispr_bulk = crispr_bulk.set_index("gene_symbol")
    track_a = track_a.set_index("gene")
    track_b = track_b.set_index("gene")
    gse240112_candidates = gse240112_candidates.set_index("gene")
    gse240112_sensitivity = gse240112_sensitivity.set_index("gene")

    rows = []
    for gene in candidate_genes:
        cb = crispr_bulk.loc[gene] if gene in crispr_bulk.index else None
        ta = track_a.loc[gene] if gene in track_a.index else None
        tb = track_b.loc[gene] if gene in track_b.index else None
        g240 = gse240112_candidates.loc[gene] if gene in gse240112_candidates.index else None
        g240sens = gse240112_sensitivity.loc[gene] if gene in gse240112_sensitivity.index else None

        rows.append(
            {
                "gene": gene,
                "crispr_effect_size": float(cb["crispr_effect_size"]) if cb is not None else float("nan"),
                "crispr_fdr": float(cb["crispr_fdr"]) if cb is not None else float("nan"),
                "crispr_sensitising_direction": True if cb is not None else None,
                "gse118713_tamr_vs_mcf7_log2fc": float(cb["tamr_vs_mcf7_log2fc"]) if cb is not None and pd.notna(cb["tamr_vs_mcf7_log2fc"]) else float("nan"),
                "gse118713_tamr_vs_mcf7_fdr": float(cb["tamr_vs_mcf7_fdr"]) if cb is not None and pd.notna(cb["tamr_vs_mcf7_fdr"]) else float("nan"),
                "gse118713_evidence_class": str(cb["evidence_class"]) if cb is not None else "not_available",
                "gse245601_track_a_epithelial_log2fc": float(ta["log2fc"]) if ta is not None and bool(ta["tested"]) else float("nan"),
                "gse245601_track_a_epithelial_fdr": float(ta["candidate_set_bh_fdr"]) if ta is not None and bool(ta["tested"]) else float("nan"),
                "gse245601_track_b_malignant_log2fc": float(tb["log2fc"]) if tb is not None and bool(tb["tested"]) else float("nan"),
                "gse245601_track_b_malignant_fdr": float(tb["candidate_set_bh_fdr"]) if tb is not None and bool(tb["tested"]) else float("nan"),
                "gse240112_tumor_cell_tested": bool(g240["tested"]) if g240 is not None else False,
                "gse240112_tumor_cell_log2fc": float(g240["log2fc"]) if g240 is not None and bool(g240["tested"]) else float("nan"),
                "gse240112_tumor_cell_p_value": float(g240["p_value"]) if g240 is not None and bool(g240["tested"]) else float("nan"),
                "gse240112_tumor_cell_candidate_bh_fdr": float(g240["candidate_set_bh_fdr"]) if g240 is not None and bool(g240["tested"]) else float("nan"),
                "gse240112_all_epithelial_log2fc": float(g240sens["all_epithelial_log2fc"]) if g240sens is not None else float("nan"),
                "gse240112_all_epithelial_fdr": float(g240sens["all_epithelial_genomewide_fdr"]) if g240sens is not None else float("nan"),
                "gse240112_direction_agreement_tumor_vs_epithelial": g240sens["direction_agreement"] if g240sens is not None else None,
                "notes": (
                    "CRISPR/GSE118713/GSE245601/GSE240112 measure different biology (functional perturbation, "
                    "acquired-resistance cell-line state, acute 12h ex vivo response, and human primary-vs-recurrent "
                    "context respectively) -- no composite score is computed; interpret each column on its own terms."
                ),
            }
        )
    out = pd.DataFrame(rows)
    logger.info("build_integrated_table: %d candidates, %d columns", len(out), len(out.columns))
    return out


def run_candidate_integration(config_path: str | Path = "config/config.yaml") -> pd.DataFrame:
    config = _load_config(config_path)
    cfg = config["gse240112"]
    ces_cfg = config["candidate_evidence_summary"]
    ci_cfg = config["gse245601_candidate_integration"]
    candidates = cfg["candidates"]["thirteen"]

    crispr_bulk = pd.read_csv(ces_cfg["output"]["evidence_summary_tsv"], sep="\t")
    track_a = pd.read_csv(ci_cfg["output"]["track_a_candidates_tsv"], sep="\t")
    track_b = pd.read_csv(ci_cfg["output"]["track_b_candidates_tsv"], sep="\t")
    gse240112_candidates = pd.read_csv(cfg["output"]["candidate_table_tsv"], sep="\t")
    gse240112_sensitivity = pd.read_csv(cfg["output"]["epithelial"]["sensitivity_comparison_tsv"], sep="\t")

    out = build_integrated_table(crispr_bulk, track_a, track_b, gse240112_candidates, gse240112_sensitivity, candidates)
    out_path = Path(cfg["output"]["integration"]["integrated_tsv"])
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(out_path, sep="\t", index=False)
    logger.info("wrote %s", out_path)
    return out


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run_candidate_integration()
