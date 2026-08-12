"""Evidence freeze Phases 4-6: the full five-layer evidence table, its
compact human-readable summary, and the seven-gene five-layer table --
all built from the frozen candidate-adjudication decision-table gene pool
(~30 genes: the seven MULTIMODAL_STRONG genes, the provisional
therapeutic shortlist, important near-miss multimodal genes, and the top
functional-sensitisation candidates) plus a fresh read of the frozen wide
evidence matrix for fields the decision table does not carry (GSE245601
malignant-track p/fdr, GSE111151 cell-line consistency, GSE240112
same-library robustness flags).

Data source: `results/tables/candidate_adjudication/final_candidate_decision_table.tsv`,
`results/tables/cross_dataset_genomewide/all_genes_cross_dataset_evidence_with_ranking.tsv`,
`results/tables/cross_dataset_genomewide/ranking_stability.tsv` (all
frozen, read-only). Reuses `src.candidate_adjudication_master_table.load_gse111151_cell_line_consistency`
unchanged.
"""

from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd
import yaml

from src.candidate_adjudication_master_table import load_gse111151_cell_line_consistency
from src.evidence_freeze_five_layer_format import (
    acute_direction,
    full_rna_pattern_4,
    resistance_direction_consistency,
    resistance_fdr05_count,
    resistance_nominal_p05_count,
    resistance_pattern_3,
)

logger = logging.getLogger(__name__)


def _load_config(config_path: str | Path) -> dict:
    with open(config_path) as f:
        return yaml.safe_load(f)


def _direction_label(log2fc: float) -> str:
    if pd.isna(log2fc):
        return "not_testable"
    return "up" if log2fc > 0 else ("down" if log2fc < 0 else "unchanged")


def build_candidate_class(row: pd.Series) -> str:
    """A short, factual label for what kind of evidence this gene has --
    never a claim about therapeutic value."""
    if row["gene"] in {"USP34", "VEZF1", "CUX1", "DPP9", "LZTR1", "SOX2", "TFAP2C"}:
        return "MULTIMODAL_STRONG"
    return row.get("evidence_category", "NA")


def build_full_evidence_table(genes: list[str], decision_table: pd.DataFrame, wide: pd.DataFrame, stability: pd.DataFrame, cell_line_consistency: pd.DataFrame) -> pd.DataFrame:
    """Builds the pre-freeze evidence table. Does NOT set freeze_shortlisted
    /freeze_shortlist_rank/eligible_for_freeze -- those are computed later
    by `src.evidence_freeze_shortlist_freeze` from THIS table, then merged
    back in by `annotate_freeze_columns` below. Splitting the steps this
    way (rather than reading `shortlist_A_multimodal_therapeutic.tsv`
    here) is deliberate: every downstream consumer of the frozen shortlist
    (figures, the compact summary, source verification) must read the
    SAME freeze output, never two different pre/post-freeze gene lists
    that could silently diverge (Phase 20 Codex review finding)."""
    w = wide.loc[wide["gene"].isin(genes)].merge(stability[["gene", "stability_label"]], on="gene", how="left")
    w = w.merge(cell_line_consistency, on="gene", how="left")
    dt = decision_table.set_index("gene")

    rows = []
    for _, r in w.iterrows():
        gene = r["gene"]
        d = dt.loc[gene] if gene in dt.index else None
        # src.evidence_freeze_five_layer_format's functions expect generic
        # `gse240112_log2fc`/`gse240112_fdr`/`gse240112_p` keys; the wide
        # table's own column names are `gse240112_tumor_*` (the tumor-cell
        # track is that dataset's one representative vote) -- renamed here,
        # never by renaming the frozen wide table itself
        fmt_row = r.rename({"gse240112_tumor_log2fc": "gse240112_log2fc", "gse240112_tumor_p": "gse240112_p", "gse240112_tumor_fdr": "gse240112_fdr"})

        rows.append(
            {
                "gene": gene,
                "candidate_class": build_candidate_class(pd.Series({"gene": gene, "evidence_category": d["evidence_category"] if d is not None else "NA"})),
                # CRISPR
                "crispr_effect": r["crispr_effect"], "crispr_p": r["crispr_p"], "crispr_fdr": r["crispr_fdr"],
                "crispr_direction": r["crispr_direction"], "crispr_sensitising_yes_no": r["crispr_direction"] == "sensitising_KO",
                # GSE118713
                "gse118713_log2fc": r["gse118713_log2fc"], "gse118713_p": r["gse118713_p"], "gse118713_fdr": r["gse118713_fdr"],
                "gse118713_direction": _direction_label(r["gse118713_log2fc"]), "gse118713_significant": pd.notna(r["gse118713_fdr"]) and r["gse118713_fdr"] < 0.05,
                # GSE240112
                "gse240112_log2fc": r["gse240112_tumor_log2fc"], "gse240112_p": r["gse240112_tumor_p"], "gse240112_fdr": r["gse240112_tumor_fdr"],
                "gse240112_direction": _direction_label(r["gse240112_tumor_log2fc"]), "gse240112_significant": pd.notna(r["gse240112_tumor_fdr"]) and r["gse240112_tumor_fdr"] < 0.05,
                "gse240112_sample_robustness": "epithelial_track_direction_agrees" if r["gse240112_track_direction_agreement"] else "epithelial_track_direction_disagrees",
                "gse240112_outlier_warning": bool(r["gse240112_outlier_fragility"]),
                # GSE111151
                "gse111151_log2fc": r["gse111151_log2fc"], "gse111151_p": r["gse111151_p"], "gse111151_fdr": r["gse111151_fdr"],
                "gse111151_direction": _direction_label(r["gse111151_log2fc"]), "gse111151_significant": pd.notna(r["gse111151_fdr"]) and r["gse111151_fdr"] < 0.05,
                "gse111151_cell_line_consistency": (f"{int(r['gse111151_consistent_cell_lines_n'])}/{int(r['gse111151_consistent_cell_lines_total'])}" if pd.notna(r.get("gse111151_consistent_cell_lines_n")) else "not_computed"),
                # GSE245601 (acute)
                "gse245601_epi_log2fc": r["gse245601_epi_log2fc"], "gse245601_epi_p": r["gse245601_epi_p"], "gse245601_epi_fdr": r["gse245601_epi_fdr"],
                "gse245601_malignant_log2fc": r["gse245601_malignant_log2fc"], "gse245601_malignant_p": r["gse245601_malignant_p"], "gse245601_malignant_fdr": r["gse245601_malignant_fdr"],
                "acute_summary_direction": acute_direction(fmt_row), "acute_significance_summary": ("track_A_significant" if pd.notna(r["gse245601_epi_fdr"]) and r["gse245601_epi_fdr"] < 0.05 else "track_A_not_significant") + (", track_B_significant" if pd.notna(r["gse245601_malignant_fdr"]) and r["gse245601_malignant_fdr"] < 0.05 else ", track_B_not_significant"),
                # summaries
                "resistance_pattern_3": resistance_pattern_3(fmt_row), "full_rna_pattern_4": full_rna_pattern_4(fmt_row),
                "resistance_fdr05_count": resistance_fdr05_count(fmt_row), "resistance_nominal_p05_count": resistance_nominal_p05_count(fmt_row),
                "resistance_direction_consistency": resistance_direction_consistency(fmt_row),
                "human_tumor_support": "significant" if (pd.notna(r["gse240112_tumor_fdr"]) and r["gse240112_tumor_fdr"] < 0.05) or (pd.notna(r["gse245601_epi_fdr"]) and r["gse245601_epi_fdr"] < 0.05) or (pd.notna(r["gse245601_malignant_fdr"]) and r["gse245601_malignant_fdr"] < 0.05) else "not_significant",
                "global_rank": d["global_rank"] if d is not None else None,
                "resistance_rank": d["resistance_rank"] if d is not None else None,
                "crispr_rank": d["crispr_screen_rank"] if d is not None else None,
                "ranking_stability": r.get("stability_label", "NA"),
                # a PURE direction check only -- "eligible" in the fuller sense (real CRISPR
                # evidence + multimodal corroboration) is `eligible_for_freeze`, computed by
                # src.evidence_freeze_shortlist_freeze and merged in by annotate_freeze_columns
                "crispr_direction_supports_inhibition_strategy": r["crispr_direction"] == "sensitising_KO",
                "main_strength": d["main_strength"] if d is not None else "NA",
                "main_limitation": d["main_weakness"] if d is not None else "NA",
            }
        )
    out = pd.DataFrame(rows).sort_values(by=["global_rank"], ascending=True, na_position="last").reset_index(drop=True)
    if out["gene"].duplicated().any():
        raise ValueError("duplicate gene rows in the evidence-freeze full table")
    logger.info("build_full_evidence_table: %d genes x %d columns", len(out), len(out.columns))
    return out


def annotate_freeze_columns(full_table: pd.DataFrame, freeze_manifest: pd.DataFrame, eligibility_audit: pd.DataFrame) -> pd.DataFrame:
    """Merges the independently-computed freeze result
    (`THERAPEUTIC_SHORTLIST_FREEZE.tsv` + `freeze_eligibility_audit.tsv`,
    both built by `src.evidence_freeze_shortlist_freeze` FROM this same
    full table) back into it. This is the only place `freeze_shortlisted`
    is ever set -- figures, the compact summary, and source verification
    must all be built from the table this function returns, so there is
    exactly one shortlist gene list in the entire freeze phase."""
    out = full_table.merge(eligibility_audit[["gene", "eligible_for_freeze", "ineligibility_reason"]], on="gene", how="left")
    rank_by_gene = freeze_manifest.set_index("gene")["freeze_rank"]
    out["freeze_shortlisted"] = out["gene"].isin(freeze_manifest["gene"])
    out["freeze_shortlist_rank"] = out["gene"].map(rank_by_gene)
    out = out.sort_values(by=["freeze_shortlisted", "global_rank"], ascending=[False, True], na_position="last").reset_index(drop=True)
    return out


def build_compact_summary(full_table: pd.DataFrame) -> pd.DataFrame:
    def freeze_status(row: pd.Series) -> str:
        if row.get("freeze_shortlisted"):
            return "FROZEN_THERAPEUTIC_SHORTLIST"
        if row.get("eligible_for_freeze"):
            return "eligible_not_shortlisted"  # passed the real-evidence + multimodal gate but ranked outside the freeze size
        if not row["crispr_direction_supports_inhibition_strategy"]:
            return "not_eligible_wrong_crispr_direction"
        reason = str(row.get("ineligibility_reason", ""))
        if reason.startswith("CRISPR FDR"):
            return "sensitising_sign_but_no_real_crispr_evidence"
        if reason.startswith("no resistance-RNA"):
            return "sensitising_functional_only_no_rna_support"
        return "not_eligible_other"

    out = pd.DataFrame(
        {
            "TherapeuticRank": full_table["freeze_shortlist_rank"],
            "Gene": full_table["gene"],
            "CRISPR_Effect": full_table["crispr_effect"].round(3),
            "CRISPR_FDR": full_table["crispr_fdr"].round(4),
            "CRISPR_Direction": full_table["crispr_direction"],
            "GSE118713_log2FC": full_table["gse118713_log2fc"].round(3),
            "GSE118713_FDR": full_table["gse118713_fdr"].round(4),
            "GSE240112_log2FC": full_table["gse240112_log2fc"].round(3),
            "GSE240112_FDR": full_table["gse240112_fdr"].round(4),
            "GSE111151_log2FC": full_table["gse111151_log2fc"].round(3),
            "GSE111151_FDR": full_table["gse111151_fdr"].round(4),
            "GSE245601_Acute_log2FC": full_table["gse245601_epi_log2fc"].round(3),
            "GSE245601_Acute_FDR": full_table["gse245601_epi_fdr"].round(4),
            "ResistancePattern": full_table["resistance_pattern_3"],
            "FullRNAPattern": full_table["full_rna_pattern_4"],
            "ResistanceSigCount": full_table["resistance_fdr05_count"],
            "HumanSupport": full_table["human_tumor_support"],
            "SampleRobustness": full_table["gse240112_sample_robustness"] + "; GSE111151 " + full_table["gse111151_cell_line_consistency"],
            "Stability": full_table["ranking_stability"],
            "CandidateClass": full_table["candidate_class"],
            "FreezeStatus": full_table.apply(freeze_status, axis=1),
            "MainStrength": full_table["main_strength"],
            "MainLimitation": full_table["main_limitation"],
        }
    )
    out = out.sort_values(by=["TherapeuticRank", "Gene"], na_position="last").reset_index(drop=True)
    logger.info("build_compact_summary: %d rows", len(out))
    return out


def build_multimodal7_five_layer(full_table: pd.DataFrame, seven_genes: list[str]) -> pd.DataFrame:
    sub = full_table.loc[full_table["gene"].isin(seven_genes)].copy()
    missing = set(seven_genes) - set(sub["gene"])
    if len(seven_genes) != 7 or len(sub) != 7 or missing:
        raise ValueError(f"expected exactly 7 MULTIMODAL_STRONG genes, requested {sorted(seven_genes)}, found {sorted(sub['gene'])}, missing {sorted(missing)}")

    def interpretation(row: pd.Series) -> str:
        if row["crispr_direction_supports_inhibition_strategy"]:
            return "sensitising_KO: knockout was relatively depleted under 4-OHT -- eligible for a knockout-sensitisation strategy on functional-direction grounds"
        return "tolerance_associated_KO: knockout was relatively favored under 4-OHT -- NOT eligible for an inhibition/sensitisation strategy despite any RNA association"

    out = pd.DataFrame(
        {
            "Gene": sub["gene"],
            "CRISPR_direction": sub["crispr_direction"],
            "CRISPR_effect_FDR": sub["crispr_effect"].round(3).astype(str) + " / " + sub["crispr_fdr"].round(4).astype(str),
            "GSE118713_direction_effect_FDR": sub["gse118713_direction"] + " / " + sub["gse118713_log2fc"].round(3).astype(str) + " / " + sub["gse118713_fdr"].round(4).astype(str),
            "GSE240112_direction_effect_FDR": sub["gse240112_direction"] + " / " + sub["gse240112_log2fc"].round(3).astype(str) + " / " + sub["gse240112_fdr"].round(4).astype(str),
            "GSE111151_direction_effect_FDR": sub["gse111151_direction"] + " / " + sub["gse111151_log2fc"].round(3).astype(str) + " / " + sub["gse111151_fdr"].round(4).astype(str),
            "GSE245601_acute_direction_effect_FDR": sub["acute_summary_direction"] + " / " + sub["gse245601_epi_log2fc"].round(3).astype(str) + " / " + sub["gse245601_epi_fdr"].round(4).astype(str),
            "ResistancePattern": sub["resistance_pattern_3"],
            "FullRNAPattern": sub["full_rna_pattern_4"],
            "TherapeuticInhibitionEligible": sub["crispr_direction_supports_inhibition_strategy"],
            "MainInterpretation": sub.apply(interpretation, axis=1),
        }
    )
    out = out.sort_values(by=["TherapeuticInhibitionEligible", "Gene"], ascending=[False, True]).reset_index(drop=True)
    return out


def run_evidence_freeze_tables(config_path: str | Path = "config/config.yaml") -> dict[str, pd.DataFrame]:
    """Builds the PRE-freeze evidence table only (freeze_shortlisted not
    yet known). Call `src.evidence_freeze_shortlist_freeze.run_freeze()`
    next on this output, then `finalize_tables_with_freeze()` below to
    produce the final, freeze-annotated tables."""
    config = _load_config(config_path)
    adj_tables = Path(config["candidate_adjudication"]["output"]["tables_dir"])
    cdx_out = config["cross_dataset_genomewide"]["output"]
    wide_dir = Path(cdx_out["wide_matrix_tsv"]).parent
    out_dir = Path(config["evidence_freeze"]["output"]["tables_dir"])
    out_dir.mkdir(parents=True, exist_ok=True)

    decision_table = pd.read_csv(adj_tables / "final_candidate_decision_table.tsv", sep="\t")
    wide = pd.read_csv(wide_dir / "all_genes_cross_dataset_evidence_with_ranking.tsv", sep="\t")
    stability = pd.read_csv(cdx_out["ranking_stability_tsv"], sep="\t")

    genes = decision_table["gene"].tolist()
    cell_line_consistency = load_gse111151_cell_line_consistency(genes, config)

    full_table = build_full_evidence_table(genes, decision_table, wide, stability, cell_line_consistency)
    full_table.to_csv(out_dir / "final_candidate_evidence.tsv", sep="\t", index=False)
    return {"full": full_table}


def finalize_tables_with_freeze(config_path: str | Path = "config/config.yaml") -> dict[str, pd.DataFrame]:
    """Re-reads the pre-freeze `final_candidate_evidence.tsv` plus the
    freeze module's own outputs, merges freeze status in, and (re)writes
    the final `final_candidate_evidence.tsv`, `final_candidate_summary.tsv`,
    and `multimodal7_five_layer_evidence.tsv` -- the versions every figure
    and downstream table must read."""
    config = _load_config(config_path)
    out_dir = Path(config["evidence_freeze"]["output"]["tables_dir"])

    full_table = pd.read_csv(out_dir / "final_candidate_evidence.tsv", sep="\t")
    freeze_manifest = pd.read_csv(out_dir / "THERAPEUTIC_SHORTLIST_FREEZE.tsv", sep="\t")
    eligibility_audit = pd.read_csv(out_dir / "freeze_eligibility_audit.tsv", sep="\t")

    full_table = annotate_freeze_columns(full_table, freeze_manifest, eligibility_audit)
    compact = build_compact_summary(full_table)
    seven_genes = config["candidate_adjudication"]["multimodal7"]["genes"]
    seven_table = build_multimodal7_five_layer(full_table, seven_genes)

    full_table.to_csv(out_dir / "final_candidate_evidence.tsv", sep="\t", index=False)
    compact.to_csv(out_dir / "final_candidate_summary.tsv", sep="\t", index=False)
    seven_table.to_csv(out_dir / "multimodal7_five_layer_evidence.tsv", sep="\t", index=False)
    logger.info("finalize_tables_with_freeze: %d genes shortlisted, tables finalized", int(full_table["freeze_shortlisted"].sum()))
    return {"full": full_table, "compact": compact, "seven": seven_table}


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run_evidence_freeze_tables()
