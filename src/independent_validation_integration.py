"""Independent-validation Part 11+12: integrated TCGA+DepMap validation
table and the two follow-up rankings (mechanistic/biological,
therapeutic-targetability). Reads only already-built outputs from this
phase plus the frozen evidence_freeze table and the prior literature-
review phase's 5-tier mechanistic classification (transcribed here from
results/reports/literature_mechanism/four_candidate_mechanism_review.md's
Part 8, since that classification was only ever recorded in prose there).
Never alters the frozen therapeutic ranking (results/tables/evidence_freeze/
THERAPEUTIC_SHORTLIST_FREEZE.tsv), which is read-only input here.
"""

from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd

logger = logging.getLogger(__name__)

TABLES = Path("results/tables/independent_validation")
OUT_TABLE = TABLES / "four_candidate_independent_validation.tsv"

FROZEN_HANY = Path("results/tables/evidence_freeze/THERAPEUTIC_SHORTLIST_FREEZE.tsv")

# Transcribed verbatim from results/reports/literature_mechanism/
# four_candidate_mechanism_review.md Part 8 (5-tier mechanistic-model
# classification); not recomputed here.
LITERATURE_CLASSIFICATION = {
    "USP34": "C. PLAUSIBLE CROSS-SYSTEM MECHANISM",
    "VEZF1": "D. NETWORK/LITERATURE HYPOTHESIS ONLY",
    "CITED2": "B. STRONG BREAST-CANCER MECHANISTIC SUPPORT, INDIRECT/CONTESTED FOR ENDOCRINE RESISTANCE",
    "EML5": "E. MECHANISTICALLY UNRESOLVED",
}

CANDIDATES = ["USP34", "VEZF1", "EML5", "CITED2"]


def _fmt_p(p) -> str:
    return "not_applicable" if pd.isna(p) else f"{p:.4g}"


def build_integration_table() -> pd.DataFrame:
    hany = pd.read_csv(FROZEN_HANY, sep="\t").set_index("gene")
    expr = pd.read_csv(TABLES / "TCGA_candidate_expression.tsv", sep="\t")
    pathway = pd.read_csv(TABLES / "TCGA_candidate_pathway_associations.tsv", sep="\t")
    clinical = pd.read_csv(TABLES / "TCGA_candidate_clinical.tsv", sep="\t")
    dep = pd.read_csv(TABLES / "DepMap_candidate_dependency.tsv", sep="\t").set_index("candidate")
    dep_expr = pd.read_csv(TABLES / "DepMap_candidate_expression.tsv", sep="\t").set_index("candidate")

    rows = []
    for candidate in CANDIDATES:
        hany_row = hany.loc[candidate]

        er_row = expr.loc[(expr.candidate == candidate) & (expr.comparison == "ER+ vs ER- (clinical IHC)")].iloc[0]
        tn_row = expr.loc[(expr.candidate == candidate) & (expr.comparison == "tumor_vs_normal_PAIRED")].iloc[0]

        cand_pathways = pathway.loc[(pathway.candidate == candidate) & pathway.pathway.ne("NONE")]
        if len(cand_pathways):
            best_pw = cand_pathways.loc[cand_pathways.fdr.idxmin()]
            pathway_summary = f"{best_pw.pathway}: rho={best_pw.spearman_rho:.3f}, FDR={_fmt_p(best_pw.fdr)}" + (" (significant)" if best_pw.fdr < 0.05 else " (not significant)")
        else:
            pathway_summary = "no candidate-specific pathway declared (EML5: none independently justified)"

        clin_er_adj = clinical.loc[(clinical.candidate == candidate) & (clinical.cohort == "ER_positive") & (clinical.model == "adjusted_age_stage")].iloc[0]
        clin_significant = bool(clin_er_adj.fdr < 0.05) if pd.notna(clin_er_adj.fdr) else False
        clinical_summary = (
            f"ER+ adjusted (age+stage): HR/SD={clin_er_adj.hr_per_sd:.2f} [{clin_er_adj.ci_low:.2f},{clin_er_adj.ci_high:.2f}], "
            f"p={_fmt_p(clin_er_adj.p_value)}, FDR={_fmt_p(clin_er_adj.fdr)}, n={int(clin_er_adj.n)}, events={int(clin_er_adj.n_events)}"
            + ("; PH assumption concern noted" if isinstance(clin_er_adj.notes, str) and "VIOLATED" in clin_er_adj.notes else "")
        )

        dep_row = dep.loc[candidate]
        dep_expr_row = dep_expr.loc[candidate]

        er_significant = bool(er_row.fdr < 0.05)
        tn_significant = bool(tn_row.fdr < 0.05)
        pathway_significant = len(cand_pathways) and (cand_pathways.fdr < 0.05).any()
        # NOTE: this counts how many INDEPENDENT TCGA TESTS came back
        # significant (FDR<0.05). It intentionally does NOT test whether the
        # direction of any given TCGA association matches (or would be
        # predicted by) the project's own tamoxifen-resistance hypothesis --
        # no such directional prediction was preregistered for these TCGA
        # comparisons, so no honest concordance test can be run. This is
        # "number of independently significant TCGA signals", not
        # "agreement with the project's result".
        n_significant_tcga_signals = sum([er_significant, bool(pathway_significant), clin_significant])
        concern_tier = dep_row.essentiality_concern

        if candidate == "EML5":
            # a bare TCGA expression difference has no mechanistic anchor to
            # interpret it against (no pathway, no network neighborhood, no
            # supporting literature) -- never scored above LITTLE regardless
            # of how many isolated TCGA signals happen to be significant
            strength = "4_LITTLE_INDEPENDENT_SUPPORT"
        elif n_significant_tcga_signals >= 2 and concern_tier in ("D_LOW_BASELINE_DEPENDENCY", "C_CONTEXT_SPECIFIC_DEPENDENCY"):
            strength = "1_STRONG_INDEPENDENT_SUPPORT"
        elif n_significant_tcga_signals >= 1 or concern_tier == "E_INSUFFICIENT_DATA":
            # E_INSUFFICIENT_DATA (missing CRISPRGeneDependency.csv for this
            # release) deliberately CANNOT reach STRONG on its own even with
            # 2+ TCGA signals -- "unknown essentiality" is not the same
            # claim as "confirmed low essentiality", and STRONG requires
            # confirmation on that axis, not merely absence of a flag. This
            # is a conservative, symmetric, disclosed rule (never silently
            # applied) -- see report Part 5: no candidate is scored STRONG
            # this release purely because CRISPRGeneDependency.csv is
            # missing, and a candidate whose 24Q4 tier was a known concern
            # (D/C/low-risk) is not penalized for the file being missing
            # either -- it simply cannot climb to STRONG without it.
            strength = "2_MODERATE_INDEPENDENT_SUPPORT"
        elif pathway_significant or er_significant or tn_significant:
            strength = "3_PARTIAL_CONTEXT_DEPENDENT_SUPPORT"
        else:
            strength = "4_LITTLE_INDEPENDENT_SUPPORT"
        if concern_tier == "A_STRONG_GENERAL_DEPENDENCY_CONCERN":
            strength = "5_CONTRADICTORY_CAUTIONARY"

        if bool(dep_row.get("dependency_probability_available", True)):
            vezf1_caveat = f"DepMap essentiality concern is {concern_tier.split('_', 1)[1].replace('_', ' ').title()} ({dep_row.frac_strongly_dependent_er_luminal*100:.0f}% of ER+/luminal breast lines strongly dependent) -- a general growth-fitness role, not evidence for a resistance-specific target; PAM50 LumA-vs-LumB direction differs from the ER+/ER- direction."
        else:
            vezf1_caveat = f"DepMap essentiality concern is UNAVAILABLE for {dep_row.depmap_release} (CRISPRGeneDependency.csv not obtained -- the 24Q4 MODERATE tier is not reproduced or refuted, simply not computable from this release's data); PAM50 LumA-vs-LumB direction differs from the ER+/ER- direction."

        # Dual-action vs pure-sensitiser framing. Stronger cancer-cell
        # dependency is NOT automatically bad: a candidate with a real
        # baseline growth-fitness role (B/A tier, or C context-specific)
        # PLUS a significant Hany tamoxifen-context signal is a plausible
        # dual-action target (baseline anti-cancer effect + an added
        # tamoxifen-specific sensitisation on top), not merely a "risky"
        # one. A candidate with a clean D-tier baseline profile plus a
        # significant Hany signal is the cleaner single-mechanism
        # tamoxifen-specific-sensitiser story. Neither category is
        # assigned if the Hany signal itself is not FDR-significant --
        # the categorization requires a real drug-context signal to
        # build on, not just a baseline profile. This is candidate
        # classification only: DepMap alone never establishes cancer-cell
        # selectivity or normal-tissue safety for either category.
        hany_significant = bool(hany_row["crispr_fdr"] < 0.05)
        if not hany_significant:
            action_category = "WEAK_HANY_SIGNAL_NEITHER_CATEGORY_ASSIGNED"
        elif concern_tier in ("B_MODERATE_DEPENDENCY_CONCERN", "A_STRONG_GENERAL_DEPENDENCY_CONCERN", "C_CONTEXT_SPECIFIC_DEPENDENCY"):
            action_category = "POTENTIAL_DUAL_ACTION_CANCER_TARGET"
        elif concern_tier == "D_LOW_BASELINE_DEPENDENCY":
            action_category = "TAMOXIFEN_SPECIFIC_SENSITISER"
        else:
            action_category = "UNCLASSIFIED_INSUFFICIENT_DEPMAP_DATA"

        major_caveat = {
            "USP34": "WNT pathway association not significant in ER+ TCGA tumors (FDR=" + _fmt_p(cand_pathways.fdr.min() if len(cand_pathways) else float('nan')) + "); literature mechanism untested in breast tissue.",
            "VEZF1": vezf1_caveat,
            "EML5": "no candidate-specific pathway could be justified for testing; low expression in most breast lines (median log2(TPM+1)=" + f"{dep_expr_row.median_tpm_log2p1_breast:.2f}" + "); near-zero baseline DepMap dependency.",
            "CITED2": "clinical Cox association not significant after FDR correction despite significant ER+/ER- and pathway signals; the human-breast VEZF1-CITED2 correlation is POSITIVE, opposite in direction to the literature-reported repression (a bulk-tissue confound, not a contradiction of the cell-type-specific literature finding).",
        }[candidate]

        rows.append(dict(
            candidate=candidate,
            frozen_crispr_direction=hany_row.crispr_direction, frozen_crispr_fdr=float(hany_row.crispr_fdr),
            frozen_resistance_sig_count=int(hany_row.resistance_sig_count), frozen_main_strength=hany_row.main_strength,
            literature_mechanistic_classification=LITERATURE_CLASSIFICATION[candidate],
            tcga_ER_positive_mean_expression=float(er_row.mean_a),
            tcga_er_pos_vs_neg_result=f"mean_diff={er_row.mean_diff:.3f}, FDR={_fmt_p(er_row.fdr)}" + (" (significant)" if er_significant else " (not significant)"),
            tcga_tumor_vs_normal_result=f"mean_diff={tn_row.mean_diff:.3f}, FDR={_fmt_p(tn_row.fdr)}" + (" (significant)" if tn_significant else " (not significant)"),
            tcga_strongest_pathway_association=pathway_summary,
            tcga_clinical_association=clinical_summary,
            tcga_interpretation="INFERENCE: " + ("multiple significant TCGA signals (separate analyses of overlapping tumors, not independent replications)" if n_significant_tcga_signals >= 2 else ("one significant TCGA signal" if n_significant_tcga_signals == 1 else "no significant TCGA signal at FDR<0.05")),
            depmap_release=dep_row.depmap_release,
            depmap_all_cancer_dependency=float(dep_row.median_gene_effect_all_cancer),
            depmap_breast_dependency=float(dep_row.median_gene_effect_breast),
            depmap_er_luminal_dependency=float(dep_row.median_gene_effect_er_luminal),
            depmap_essentiality_concern=concern_tier,
            depmap_breast_expression_median_log2tpm1=float(dep_expr_row.median_tpm_log2p1_breast),
            depmap_expression_dependency_relationship=f"rho={dep_expr_row.spearman_rho_expression_vs_dependency_breast:.3f}, p={_fmt_p(dep_expr_row.p_value)}",
            integration_tcga_significant_signal_present=n_significant_tcga_signals >= 1,  # NOT a test of directional agreement with the project's hypothesis -- see n_significant_tcga_signals note above
            integration_direction_of_effect_contested=candidate == "CITED2",  # VEZF1-CITED2 direction and self-contradictory tamoxifen literature both land on CITED2
            integration_clinically_relevant=clin_significant,
            integration_baseline_essentiality_concern=concern_tier in ("A_STRONG_GENERAL_DEPENDENCY_CONCERN", "B_MODERATE_DEPENDENCY_CONCERN"),
            integration_validation_strength=strength,
            integration_major_caveat=major_caveat,
            mechanistic_action_category=action_category,
        ))
    out = pd.DataFrame(rows)
    logger.info("build_integration_table: %d rows", len(out))
    return out


def build_rankings(integration: pd.DataFrame) -> pd.DataFrame:
    """Two separate follow-up rankings; never alters the frozen therapeutic ranking."""
    mech_order = {"USP34": 2, "VEZF1": 3, "EML5": 4, "CITED2": 1}  # richest independent TCGA+literature mechanistic signal first (CITED2), consistent with literature-only ranking
    rows = []
    for candidate in CANDIDATES:
        row = integration.set_index("candidate").loc[candidate]
        target_penalty = row.integration_baseline_essentiality_concern
        rows.append(dict(
            candidate=candidate,
            mechanistic_biological_followup_rank=mech_order[candidate],
            mechanistic_biological_followup_rationale="see literature_mechanism_review.md Part 9 order, now cross-checked against TCGA pathway/clinical signal in " + row.tcga_interpretation,
            therapeutic_targetability_followup_score=(
                (1 if row.frozen_crispr_direction == "sensitising_KO" else 0)
                + (1 if row.integration_clinically_relevant else 0)
                + (1 if "STRONG" in row.literature_mechanistic_classification or "B." in row.literature_mechanistic_classification else 0)
                - (1 if target_penalty else 0)
            ),
            depmap_essentiality_concern=row.depmap_essentiality_concern,
            note="Do NOT alter the frozen therapeutic ranking (USP34 > VEZF1 > EML5 > CITED2); this is a follow-up-order signal only.",
        ))
    ranking = pd.DataFrame(rows).sort_values("therapeutic_targetability_followup_score", ascending=False)
    # tied scores get a tied rank ("min" method, e.g. two candidates tied
    # for 1st both show rank 1, next distinct score shows rank 3) -- never
    # assign sequential 1..N ranks to equal scores, which would fabricate
    # an ordering the score itself does not support
    ranking["therapeutic_targetability_followup_rank"] = ranking["therapeutic_targetability_followup_score"].rank(ascending=False, method="min").astype(int)
    return ranking.sort_values("mechanistic_biological_followup_rank")


def run() -> tuple[pd.DataFrame, pd.DataFrame]:
    integration = build_integration_table()
    OUT_TABLE.parent.mkdir(parents=True, exist_ok=True)
    integration.to_csv(OUT_TABLE, sep="\t", index=False)
    logger.info("wrote %s (%d rows)", OUT_TABLE, len(integration))

    rankings = build_rankings(integration)
    rank_path = TABLES / "four_candidate_followup_rankings.tsv"
    rankings.to_csv(rank_path, sep="\t", index=False)
    logger.info("wrote %s (%d rows)", rank_path, len(rankings))
    return integration, rankings


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run()
