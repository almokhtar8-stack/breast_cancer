"""GSE111151 Phase 10: additive 5-layer evidence integration table for the
frozen 13 candidates -- CRISPR screen, GSE118713 bulk acquired-resistance
RNA, GSE245601 acute (12h ex vivo) scRNA, GSE240112 primary-vs-recurrent
scRNA, and now GSE111151 independent tamoxifen-resistant-cell-line
validation. This is a NEW table; it never overwrites any existing frozen
master table.

No weighted composite score is computed -- the five datasets measure
different biology (CRISPR = functional perturbation under 4-OHT;
GSE118713 = one acquired-resistance cell-line panel's bulk RNA;
GSE245601 = acute 12h ex vivo human tumor response; GSE240112 = human
primary-vs-recurrent/treatment-history context; GSE111151 = a second,
independent panel of 4 isogenic acquired-tamoxifen-resistance cell
lines) and are not assumed to agree in direction.

Data sources: as already frozen/committed in this repo, plus the
GSE111151 tables produced 2026-08-12.
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
    gse111151_candidates: pd.DataFrame,
    gse111151_classification: pd.DataFrame,
    candidate_genes: list[str],
) -> pd.DataFrame:
    """One row per frozen candidate, exactly once. Every source table is
    joined by gene; a gene absent from a given source is reported with
    NA in that source's columns, never dropped from the output."""
    crispr_bulk = crispr_bulk.set_index("gene_symbol")
    track_a = track_a.set_index("gene")
    track_b = track_b.set_index("gene")
    gse240112_candidates = gse240112_candidates.set_index("gene")
    gse111151_candidates = gse111151_candidates.set_index("gene")
    gse111151_classification = gse111151_classification.set_index("gene")

    rows = []
    for gene in candidate_genes:
        cb = crispr_bulk.loc[gene] if gene in crispr_bulk.index else None
        ta = track_a.loc[gene] if gene in track_a.index else None
        tb = track_b.loc[gene] if gene in track_b.index else None
        g240 = gse240112_candidates.loc[gene] if gene in gse240112_candidates.index else None
        g111 = gse111151_candidates.loc[gene] if gene in gse111151_candidates.index else None
        g111cls = gse111151_classification.loc[gene] if gene in gse111151_classification.index else None

        rows.append(
            {
                "gene": gene,
                "crispr_effect_size": float(cb["crispr_effect_size"]) if cb is not None else float("nan"),
                "crispr_fdr": float(cb["crispr_fdr"]) if cb is not None else float("nan"),
                "gse118713_tamr_vs_mcf7_log2fc": float(cb["tamr_vs_mcf7_log2fc"]) if cb is not None and pd.notna(cb["tamr_vs_mcf7_log2fc"]) else float("nan"),
                "gse118713_tamr_vs_mcf7_fdr": float(cb["tamr_vs_mcf7_fdr"]) if cb is not None and pd.notna(cb["tamr_vs_mcf7_fdr"]) else float("nan"),
                "gse245601_track_a_epithelial_log2fc": float(ta["log2fc"]) if ta is not None and bool(ta["tested"]) else float("nan"),
                "gse245601_track_a_epithelial_fdr": float(ta["candidate_set_bh_fdr"]) if ta is not None and bool(ta["tested"]) else float("nan"),
                "gse245601_track_b_malignant_log2fc": float(tb["log2fc"]) if tb is not None and bool(tb["tested"]) else float("nan"),
                "gse245601_track_b_malignant_fdr": float(tb["candidate_set_bh_fdr"]) if tb is not None and bool(tb["tested"]) else float("nan"),
                "gse240112_tumor_cell_log2fc": float(g240["log2fc"]) if g240 is not None and bool(g240["tested"]) else float("nan"),
                "gse240112_tumor_cell_candidate_bh_fdr": float(g240["candidate_set_bh_fdr"]) if g240 is not None and bool(g240["tested"]) else float("nan"),
                "gse111151_tested": bool(g111["tested"]) if g111 is not None else False,
                "gse111151_log2fc": float(g111["log2fc"]) if g111 is not None and bool(g111["tested"]) else float("nan"),
                "gse111151_p_value": float(g111["p_value"]) if g111 is not None and bool(g111["tested"]) else float("nan"),
                "gse111151_candidate_bh_fdr": float(g111["candidate_set_bh_fdr"]) if g111 is not None and bool(g111["tested"]) else float("nan"),
                "gse111151_n_cell_lines_consistent": g111cls["n_cell_lines_consistent"] if g111cls is not None else float("nan"),
                "gse111151_n_cell_lines_with_both_arms": g111cls["n_cell_lines_with_both_arms"] if g111cls is not None else float("nan"),
                "gse111151_interpretation": g111cls["classification"] if g111cls is not None else "not_available",
                "notes": (
                    "CRISPR/GSE118713/GSE245601/GSE240112/GSE111151 measure different biology (functional perturbation, "
                    "one acquired-resistance cell-line panel, acute 12h ex vivo response, human primary-vs-recurrent context, "
                    "and a second independent acquired-resistance cell-line panel respectively) -- no composite score is "
                    "computed; interpret each column on its own terms."
                ),
            }
        )
    out = pd.DataFrame(rows)
    logger.info("build_integrated_table: %d candidates, %d columns", len(out), len(out.columns))
    return out


def run_candidate_integration(config_path: str | Path = "config/config.yaml") -> pd.DataFrame:
    config = _load_config(config_path)
    cfg = config["gse111151"]
    ces_cfg = config["candidate_evidence_summary"]
    ci_cfg = config["gse245601_candidate_integration"]
    g240_cfg = config["gse240112"]
    candidates = cfg["candidates"]["thirteen"]

    crispr_bulk = pd.read_csv(ces_cfg["output"]["evidence_summary_tsv"], sep="\t")
    track_a = pd.read_csv(ci_cfg["output"]["track_a_candidates_tsv"], sep="\t")
    track_b = pd.read_csv(ci_cfg["output"]["track_b_candidates_tsv"], sep="\t")
    gse240112_candidates = pd.read_csv(g240_cfg["output"]["candidate_table_tsv"], sep="\t")
    gse111151_candidates = pd.read_csv(cfg["output"]["candidate_table_tsv"], sep="\t")
    gse111151_classification = pd.read_csv(cfg["output"]["classification_tsv"], sep="\t")

    out = build_integrated_table(crispr_bulk, track_a, track_b, gse240112_candidates, gse111151_candidates, gse111151_classification, candidates)
    out_path = Path(cfg["output"]["integration"]["integrated_tsv"])
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(out_path, sep="\t", index=False)
    logger.info("wrote %s", out_path)
    return out


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run_candidate_integration()
