"""Cross-dataset genome-wide integration, Phases 4-5: the canonical
long-form evidence table (one row per gene x independent dataset -- never
per track) and the wide master matrix (one row per gene, every raw
per-track column preserved for transparency). All downstream summary/
ranking/Top-20 tables are derived programmatically from these two tables,
never hand-constructed.

For a multi-track dataset (GSE245601: Track A/Track B; GSE240112:
tumor-cell/all-epithelial), the long-form table's single row uses that
dataset's designated PRIMARY track as the representative effect/p/fdr
(Track A for GSE245601 -- the larger, all-epithelial population;
tumor-cell for GSE240112 -- the primary biological analysis per
docs/GSE240112_PREANALYSIS.md), with the secondary track's info recorded
in ``robustness_notes`` (direction agreement) rather than a second row --
this is what "one independent-dataset slot" means operationally. The wide
matrix retains both tracks' raw columns in full; percentile computation
(Phase 6-9) reads the wide matrix directly so no secondary-track
information is lost, only prevented from casting a second vote.

Data source: see src.cross_dataset_gene_mapping module docstring.
"""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

logger = logging.getLogger(__name__)


def _load_config(config_path: str | Path) -> dict:
    with open(config_path) as f:
        return yaml.safe_load(f)


def _crispr_direction(effect: float) -> str:
    if pd.isna(effect):
        return "not_applicable"
    if effect < 0:
        return "sensitising_KO"
    if effect > 0:
        return "tolerance_associated_KO"
    return "approximately_neutral"


def build_wide_matrix(
    universe: pd.DataFrame,
    crispr: pd.DataFrame,
    gse118713: pd.DataFrame,
    gse245601_a: pd.DataFrame,
    gse245601_b: pd.DataFrame,
    gse240112_tumor: pd.DataFrame,
    gse240112_epi: pd.DataFrame,
    gse111151: pd.DataFrame,
) -> pd.DataFrame:
    """One row per gene in the full union universe; every dataset's raw
    columns joined in by symbol, left-joined so a gene absent/untestable
    in a dataset gets explicit NA in that dataset's columns rather than
    being dropped from the table."""
    out = universe[["gene"]].copy()

    c = crispr.rename(columns={"symbol": "gene"}).set_index("gene")
    out["crispr_effect"] = out["gene"].map(c["effect"])
    out["crispr_p"] = out["gene"].map(c["p_value"])
    out["crispr_fdr"] = out["gene"].map(c["fdr"])
    out["crispr_testable"] = out["gene"].isin(c.index)
    out["crispr_direction"] = out["crispr_effect"].map(_crispr_direction)

    g118 = gse118713.rename(columns={"symbol": "gene"}).set_index("gene")
    out["gse118713_log2fc"] = out["gene"].map(g118["effect"])
    out["gse118713_p"] = out["gene"].map(g118["p_value"])
    out["gse118713_fdr"] = out["gene"].map(g118["fdr"])
    out["gse118713_testable"] = out["gene"].isin(g118.index)
    out["gse118713_direction"] = np.where(out["gse118713_log2fc"] > 0, "up_in_TAMR", np.where(out["gse118713_log2fc"] < 0, "down_in_TAMR", pd.NA))

    ga = gse245601_a.rename(columns={"symbol": "gene"}).set_index("gene")
    gb = gse245601_b.rename(columns={"symbol": "gene"}).set_index("gene")
    out["gse245601_epi_log2fc"] = out["gene"].map(ga["effect"])
    out["gse245601_epi_p"] = out["gene"].map(ga["p_value"])
    out["gse245601_epi_fdr"] = out["gene"].map(ga["fdr"])
    out["gse245601_malignant_log2fc"] = out["gene"].map(gb["effect"])
    out["gse245601_malignant_p"] = out["gene"].map(gb["p_value"])
    out["gse245601_malignant_fdr"] = out["gene"].map(gb["fdr"])
    out["gse245601_testable"] = out["gene"].isin(ga.index) | out["gene"].isin(gb.index)
    epi_sign = np.sign(out["gse245601_epi_log2fc"])
    mal_sign = np.sign(out["gse245601_malignant_log2fc"])
    out["gse245601_track_direction_agreement"] = np.where(
        out["gse245601_epi_log2fc"].notna() & out["gse245601_malignant_log2fc"].notna(), epi_sign == mal_sign, pd.NA
    )

    gt = gse240112_tumor.rename(columns={"symbol": "gene"}).set_index("gene")
    ge = gse240112_epi.rename(columns={"symbol": "gene"}).set_index("gene")
    out["gse240112_tumor_log2fc"] = out["gene"].map(gt["effect"])
    out["gse240112_tumor_p"] = out["gene"].map(gt["p_value"])
    out["gse240112_tumor_fdr"] = out["gene"].map(gt["fdr"])
    out["gse240112_epi_log2fc"] = out["gene"].map(ge["effect"])
    out["gse240112_epi_p"] = out["gene"].map(ge["p_value"])
    out["gse240112_epi_fdr"] = out["gene"].map(ge["fdr"])
    # Tumor-cell track ONLY -- this must match exactly what contributes the
    # dataset's percentile/vote (Phase 7). Using OR-of-both-tracks here was a
    # real bug (caught by Codex review): it marked genes "testable" purely
    # because the epithelial sensitivity track saw them, even though such a
    # gene gets NO gse240112 percentile at all (compute_dataset_percentiles
    # uses the tumor-cell track exclusively), silently inflating coverage
    # tier and ranking eligibility for genes with no real GSE240112 vote.
    out["gse240112_testable"] = out["gene"].isin(gt.index)
    tumor_sign = np.sign(out["gse240112_tumor_log2fc"])
    epi_sign2 = np.sign(out["gse240112_epi_log2fc"])
    out["gse240112_track_direction_agreement"] = np.where(
        out["gse240112_tumor_log2fc"].notna() & out["gse240112_epi_log2fc"].notna(), tumor_sign == epi_sign2, pd.NA
    )
    # "fragility" flag: tested in the sensitivity (epithelial) track but not in the primary tumor-cell track
    out["gse240112_outlier_fragility"] = out["gse240112_tumor_log2fc"].isna() & out["gse240112_epi_log2fc"].notna()

    g111 = gse111151.rename(columns={"symbol": "gene"}).set_index("gene")
    out["gse111151_log2fc"] = out["gene"].map(g111["effect"])
    out["gse111151_p"] = out["gene"].map(g111["p_value"])
    out["gse111151_fdr"] = out["gene"].map(g111["fdr"])
    out["gse111151_testable"] = out["gene"].isin(g111.index)

    logger.info("build_wide_matrix: %d genes, %d columns", len(out), len(out.columns))
    return out


def build_evidence_long(wide: pd.DataFrame) -> pd.DataFrame:
    """One row per gene x independent dataset (exactly 5 dataset rows
    per gene, never per track). GSE245601 uses Track A (all epithelial)
    and GSE240112 uses the tumor-cell track as each dataset's
    representative primary contribution; the other track's agreement is
    recorded in robustness_notes, not a second row."""
    records = []
    for _, row in wide.iterrows():
        gene = row["gene"]
        records.append(
            {
                "gene": gene, "dataset": "crispr", "biological_context": "functional_perturbation_screen", "modality": "CRISPR_KO_screen",
                "effect": row["crispr_effect"], "p_value": row["crispr_p"], "fdr": row["crispr_fdr"], "testable": row["crispr_testable"],
                "direction_native": row["crispr_direction"], "n_samples_or_blocks": pd.NA,
                "coverage_notes": "single guide-level model fit, no secondary contrast", "robustness_notes": "",
            }
        )
        records.append(
            {
                "gene": gene, "dataset": "gse118713", "biological_context": "acquired_chronic_resistance_cell_line_bulk", "modality": "bulk_RNAseq",
                "effect": row["gse118713_log2fc"], "p_value": row["gse118713_p"], "fdr": row["gse118713_fdr"], "testable": row["gse118713_testable"],
                "direction_native": row["gse118713_direction"], "n_samples_or_blocks": 6,
                "coverage_notes": "primary contrast TAMR_vs_MCF7; TAMR_vs_FASR/FASR_vs_MCF7 are secondary, not a second dataset vote", "robustness_notes": "",
            }
        )
        agree = row["gse245601_track_direction_agreement"]
        records.append(
            {
                "gene": gene, "dataset": "gse245601", "biological_context": "acute_12h_ex_vivo_tamoxifen_response_human_tumor", "modality": "scRNA_pseudobulk",
                "effect": row["gse245601_epi_log2fc"], "p_value": row["gse245601_epi_p"], "fdr": row["gse245601_epi_fdr"], "testable": row["gse245601_testable"],
                "direction_native": "up_after_acute_tamoxifen" if pd.notna(row["gse245601_epi_log2fc"]) and row["gse245601_epi_log2fc"] > 0 else ("down_after_acute_tamoxifen" if pd.notna(row["gse245601_epi_log2fc"]) and row["gse245601_epi_log2fc"] < 0 else pd.NA),
                "n_samples_or_blocks": 10,
                "coverage_notes": "Track A (all epithelial) is the representative track; Track B (strict malignant) is the same libraries, not a second dataset vote",
                "robustness_notes": f"track_direction_agreement={agree}" if pd.notna(agree) else "track_b_untestable_or_track_a_untestable",
            }
        )
        agree2 = row["gse240112_track_direction_agreement"]
        records.append(
            {
                "gene": gene, "dataset": "gse240112", "biological_context": "human_primary_vs_recurrent_tumor_context", "modality": "scRNA_pseudobulk",
                "effect": row["gse240112_tumor_log2fc"], "p_value": row["gse240112_tumor_p"], "fdr": row["gse240112_tumor_fdr"], "testable": row["gse240112_testable"],
                "direction_native": "up_in_recurrent" if pd.notna(row["gse240112_tumor_log2fc"]) and row["gse240112_tumor_log2fc"] > 0 else ("down_in_recurrent" if pd.notna(row["gse240112_tumor_log2fc"]) and row["gse240112_tumor_log2fc"] < 0 else pd.NA),
                "n_samples_or_blocks": 6,
                "coverage_notes": "tumor-cell track is the primary analysis; all-epithelial (same libraries) is sensitivity-only, not a second dataset vote",
                "robustness_notes": f"track_direction_agreement={agree2}" if pd.notna(agree2) else "epithelial_sensitivity_untestable_or_tumor_untestable",
            }
        )
        records.append(
            {
                "gene": gene, "dataset": "gse111151", "biological_context": "acquired_chronic_resistance_second_independent_cell_line_panel", "modality": "bulk_RNAseq",
                "effect": row["gse111151_log2fc"], "p_value": row["gse111151_p"], "fdr": row["gse111151_fdr"], "testable": row["gse111151_testable"],
                "direction_native": "up_in_resistant" if pd.notna(row["gse111151_log2fc"]) and row["gse111151_log2fc"] > 0 else ("down_in_resistant" if pd.notna(row["gse111151_log2fc"]) and row["gse111151_log2fc"] < 0 else pd.NA),
                "n_samples_or_blocks": 4,
                "coverage_notes": "cell-line-blocked model, 4 blocks, unbalanced resistant-subline count per block", "robustness_notes": "",
            }
        )
    out = pd.DataFrame(records)
    logger.info("build_evidence_long: %d rows (%d genes x 5 datasets)", len(out), len(wide))
    return out


def run_evidence_tables(config_path: str | Path = "config/config.yaml") -> dict[str, pd.DataFrame]:
    from src.cross_dataset_gene_mapping import run_gene_mapping

    config = _load_config(config_path)
    cfg = config["cross_dataset_genomewide"]
    out = cfg["output"]

    loaded = run_gene_mapping(config_path)
    wide = build_wide_matrix(
        loaded["universe"], loaded["crispr"], loaded["gse118713"], loaded["gse245601_track_a"], loaded["gse245601_track_b"],
        loaded["gse240112_tumor"], loaded["gse240112_epi"], loaded["gse111151"],
    )
    long_df = build_evidence_long(wide)

    wide.to_csv(out["wide_matrix_tsv"], sep="\t", index=False)
    long_df.to_csv(out["evidence_long_tsv"], sep="\t", index=False)
    logger.info("wrote %s and %s", out["wide_matrix_tsv"], out["evidence_long_tsv"])

    return {"wide": wide, "long": long_df, **loaded}


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run_evidence_tables()
