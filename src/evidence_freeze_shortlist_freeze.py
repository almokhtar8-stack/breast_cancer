"""Evidence freeze Phases 11-14: freezes the therapeutic (knockout
-inhibition/sensitisation) shortlist directly from
`final_candidate_evidence.tsv`, independently of (though consistent
with) the candidate-adjudication phase's `shortlist_A_multimodal_therapeutic.tsv`
-- this is a fresh application of the eligibility gate and ranking
hierarchy to the frozen evidence table, not a copy of that file, so a bug
in the earlier shortlist-building code cannot silently carry into the
freeze.

Eligibility (hard gate, both required):
1. CRISPR direction is `sensitising_KO`.
2. CRISPR evidence is real, not just a meaningless sign: FDR<0.25 (the
   same "at least WEAK axis-A evidence" bound established during
   candidate adjudication -- a gene with CRISPR FDR~0.8 must never enter
   an inhibition shortlist merely because its noisy effect happened to be
   negative; this project's own resistance-leader genes SLC4A10, DMRTA1,
   GREB1, DLX2, GJB2, ACOT4, IL20 all carry a `sensitising_KO` label with
   CRISPR FDR>0.7 and are excluded by this gate).

RNA significance is explicitly NOT part of the eligibility gate (per this
phase's instruction) -- it is used only as ranking criteria 3-5 below, in
line with "a functional target does not have to be transcriptionally
differentially expressed." A gate on "some multimodal corroboration"
(>=1 resistance dataset FDR<0.05 OR significant human-tumor support) is
retained, however, because it is what distinguishes this list from the
project's separate FUNCTIONAL_SENSITISATION list (RNA support explicitly
not required there, per Phase 13's own class B/C definitions) -- without
it, List A and List C would be identical lists under a different name.

Ranking hierarchy (never a hidden numeric score; every criterion is an
already-computed, disclosed, banded comparison -- a tie within a band is
expected and is exactly where the next criterion decides the order):
2. CRISPR-strength band (VERY_STRONG/STRONG/MODERATE/WEAK, from FDR)
3. resistance-state RNA evidence band (from FDR<0.05 count + consensus)
4. resistance direction consistency (all_*/majority_*/mixed)
5. human-tumor evidence band
6. sample/replicate robustness (GSE111151 cell-line consistency)
7. cross-dataset leave-one-out stability
8. gene symbol (deterministic final tie-break only)

Data source: `results/tables/evidence_freeze/final_candidate_evidence.tsv`
(read-only).
"""

from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd
import yaml

logger = logging.getLogger(__name__)

CRISPR_BAND_ORDER = {"VERY_STRONG": 0, "STRONG": 1, "MODERATE": 2, "WEAK": 3}
CONSISTENCY_BAND = {"all_up": 0, "all_down": 0, "majority_up": 1, "majority_down": 1, "mixed": 2, "insufficient": 3}
STABILITY_BAND = {"ROBUST": 0, "MODERATELY_STABLE": 1, "DATASET_DEPENDENT": 2}


def _load_config(config_path: str | Path) -> dict:
    with open(config_path) as f:
        return yaml.safe_load(f)


def _crispr_band(fdr: float) -> str | None:
    if pd.isna(fdr):
        return None
    if fdr < 0.01:
        return "VERY_STRONG"
    if fdr < 0.05:
        return "STRONG"
    if fdr < 0.10:
        return "MODERATE"
    if fdr < 0.25:
        return "WEAK"
    return None


def _resistance_band(fdr05_count: int, consensus: str) -> str:
    if fdr05_count >= 2 and consensus in ("all_up", "all_down"):
        return "VERY_STRONG"
    if fdr05_count >= 1 and consensus in ("all_up", "all_down", "majority_up", "majority_down"):
        return "STRONG"
    return "WEAK"


def _human_band(human_support: str) -> str:
    return "STRONG" if human_support == "significant" else "NO_EVIDENCE"


def _sample_robustness_band(cell_line_consistency: str) -> str:
    """Parses the GSE111151 "n/total" cell-line-consistency string into a
    band. Never a raw numeric sort key -- ties within a band (e.g. 3/4 and
    3/4) are expected and are broken by the next criterion, not by
    comparing 0.75 to 0.7499999."""
    if not isinstance(cell_line_consistency, str) or "/" not in cell_line_consistency:
        return "UNKNOWN"
    n, total = cell_line_consistency.split("/")
    try:
        n, total = int(n), int(total)
    except ValueError:
        return "UNKNOWN"
    if total == 0:
        return "UNKNOWN"
    frac = n / total
    if frac >= 0.99:
        return "HIGH"
    if frac >= 0.5:
        return "MEDIUM"
    return "LOW"


def determine_eligibility(full_table: pd.DataFrame) -> pd.DataFrame:
    out = full_table.copy()
    out["_crispr_band"] = out["crispr_fdr"].map(_crispr_band)
    sensitising = out["crispr_direction"] == "sensitising_KO"
    real_crispr_evidence = out["_crispr_band"].notna()
    some_multimodal_support = (out["resistance_fdr05_count"].fillna(0) >= 1) | (out["human_tumor_support"] == "significant")
    out["eligible_for_freeze"] = sensitising & real_crispr_evidence & some_multimodal_support
    out["ineligibility_reason"] = ""
    out.loc[~sensitising, "ineligibility_reason"] = "CRISPR direction is tolerance_associated_KO, not sensitising_KO"
    out.loc[sensitising & ~real_crispr_evidence, "ineligibility_reason"] = "CRISPR FDR>=0.25 or untestable -- sensitising sign carries no real evidence"
    out.loc[sensitising & real_crispr_evidence & ~some_multimodal_support, "ineligibility_reason"] = "no resistance-RNA dataset FDR<0.05 and no significant human-tumor support (functional-only -- see List C instead)"
    return out


def rank_eligible(eligible: pd.DataFrame) -> pd.DataFrame:
    """Sort order matches the module docstring's numbered hierarchy
    EXACTLY: CRISPR-strength band (2) -> resistance-RNA band (3) ->
    resistance direction consistency (4) -> human-tumor evidence band (5)
    -> sample/replicate robustness band (6) -> leave-one-out stability (7)
    -> gene symbol (8, deterministic final tie-break only). A prior
    version of this function put human evidence ahead of consistency and
    omitted sample robustness entirely, silently diverging from the
    declared rules -- fixed after the Phase 20 Codex review caught it."""
    df = eligible.copy()
    df["_resistance_band"] = df.apply(lambda r: _resistance_band(r["resistance_fdr05_count"], r["resistance_direction_consistency"]), axis=1)
    df["_consistency_band"] = df["resistance_direction_consistency"].map(CONSISTENCY_BAND).fillna(3)
    df["_human_band"] = df["human_tumor_support"].map(_human_band)
    df["_robustness_band"] = df["gse111151_cell_line_consistency"].map(_sample_robustness_band)
    df["_stability_band"] = df["ranking_stability"].map(STABILITY_BAND).fillna(3)

    band_rank = {"VERY_STRONG": 0, "STRONG": 1, "WEAK": 2, "NO_EVIDENCE": 3}
    robustness_rank = {"HIGH": 0, "MEDIUM": 1, "LOW": 2, "UNKNOWN": 3}
    df["_crispr_band_rank"] = df["_crispr_band"].map(CRISPR_BAND_ORDER)
    df["_resistance_band_rank"] = df["_resistance_band"].map(band_rank)
    df["_human_band_rank"] = df["_human_band"].map(band_rank)
    df["_robustness_band_rank"] = df["_robustness_band"].map(robustness_rank)

    ranked = df.sort_values(
        by=["_crispr_band_rank", "_resistance_band_rank", "_consistency_band", "_human_band_rank", "_robustness_band_rank", "_stability_band", "gene"],
        ascending=True, kind="mergesort",
    ).reset_index(drop=True)
    return ranked


def build_freeze(full_table: pd.DataFrame, max_n: int = 5) -> tuple[pd.DataFrame, pd.DataFrame]:
    audited = determine_eligibility(full_table)
    eligible = audited.loc[audited["eligible_for_freeze"]]
    ranked = rank_eligible(eligible)
    frozen = ranked.head(max_n).copy()
    frozen.insert(0, "freeze_rank", range(1, len(frozen) + 1))
    logger.info("build_freeze: %d genes eligible, %d frozen: %s", len(eligible), len(frozen), frozen["gene"].tolist())
    return frozen, audited


def build_freeze_manifest_table(frozen: pd.DataFrame) -> pd.DataFrame:
    def reason(row: pd.Series) -> str:
        parts = [f"CRISPR {row['_crispr_band']} (FDR={row['crispr_fdr']:.3g}, sensitising_KO)"]
        parts.append(f"resistance {row['_resistance_band']} ({int(row['resistance_fdr05_count'])}/3 FDR<0.05, {row['resistance_direction_consistency']})")
        if row["_human_band"] == "STRONG":
            parts.append("human-tumor support significant")
        parts.append(f"stability {row['ranking_stability']}")
        return "; ".join(parts)

    out = pd.DataFrame(
        {
            "freeze_rank": frozen["freeze_rank"], "gene": frozen["gene"],
            "crispr_effect": frozen["crispr_effect"], "crispr_fdr": frozen["crispr_fdr"], "crispr_direction": frozen["crispr_direction"],
            "resistance_pattern": frozen["resistance_pattern_3"], "full_rna_pattern": frozen["full_rna_pattern_4"],
            "resistance_sig_count": frozen["resistance_fdr05_count"], "human_support": frozen["human_tumor_support"],
            "sample_robustness": frozen["gse111151_cell_line_consistency"],
            "main_strength": frozen["main_strength"], "main_limitation": frozen["main_limitation"],
            "freeze_reason": frozen.apply(reason, axis=1),
        }
    )
    return out


def run_freeze(config_path: str | Path = "config/config.yaml") -> dict[str, pd.DataFrame]:
    config = _load_config(config_path)
    out_dir = Path(config["evidence_freeze"]["output"]["tables_dir"])
    full_table = pd.read_csv(out_dir / "final_candidate_evidence.tsv", sep="\t")

    frozen, audited = build_freeze(full_table)
    manifest = build_freeze_manifest_table(frozen)

    manifest.to_csv(out_dir / "THERAPEUTIC_SHORTLIST_FREEZE.tsv", sep="\t", index=False)
    audited[["gene", "eligible_for_freeze", "ineligibility_reason", "_crispr_band"]].to_csv(out_dir / "freeze_eligibility_audit.tsv", sep="\t", index=False)
    return {"frozen": frozen, "manifest": manifest, "audited": audited}


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run_freeze()
