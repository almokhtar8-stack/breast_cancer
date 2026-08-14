"""Independent-validation Parts 7-9: DepMap baseline CRISPR dependency,
essentiality-concern classification, and matched expression. Supports
multiple releases explicitly (see independent_validation_depmap_data.py).
Both 24Q4 and 26Q1 currently have CRISPRGeneDependency.csv verified and
in use (26Q1's was added in a follow-up manual download); if a future
release is missing that file, the probability-based "fraction strongly
dependent" metric and its essentiality_concern tier are reported as
unavailable (E_INSUFFICIENT_DATA) for that release only, never estimated
from gene_effect or any other quantity -- see
has_dependency_probability() in independent_validation_depmap_data.py.

DepMap's standard genome-wide CRISPR dependency screen is NOT the same
experiment as the project's own Hany screen: DepMap measures baseline
fitness effect of knockout in unconditioned standard culture across many
cancer cell lines (a general-essentiality readout), while Hany measured a
drug-context interaction (E2+4-OHT vs E2) in one competitive screen. A
concordant sign between the two is never treated as replication of the
same phenomenon; Part 8 exists specifically to keep this distinction
explicit.

"More negative gene-effect score = greater baseline dependency" (Chronos
scale, this release). "Strongly dependent" uses DepMap's own per-line
dependency PROBABILITY (CRISPRGeneDependency.csv) at the 0.5 threshold
pinned in config/config.yaml's independent_validation.depmap section --
the same threshold DepMap's own portal uses to call a line "dependent."
"""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

from src.independent_validation_depmap_data import has_dependency_probability, load_config, load_expression, load_gene_dependency, load_gene_effect, load_model

logger = logging.getLogger(__name__)

CANDIDATES = ["USP34", "VEZF1", "EML5", "CITED2"]
OUT_DEPENDENCY = Path("results/tables/independent_validation/DepMap_candidate_dependency.tsv")
OUT_EXPRESSION = Path("results/tables/independent_validation/DepMap_candidate_expression.tsv")

FROZEN_HANY = Path("results/tables/evidence_freeze/THERAPEUTIC_SHORTLIST_FREEZE.tsv")


def _classify_concern(frac_all: float, frac_breast: float) -> str:
    if pd.isna(frac_all):
        return "E_INSUFFICIENT_DATA"
    if frac_all >= 0.90:
        return "A_STRONG_GENERAL_DEPENDENCY_CONCERN"
    if frac_all >= 0.50:
        return "B_MODERATE_DEPENDENCY_CONCERN"
    if (frac_breast - frac_all) >= 0.20 and frac_breast >= 0.30:
        return "C_CONTEXT_SPECIFIC_DEPENDENCY"
    if frac_all < 0.10:
        return "D_LOW_BASELINE_DEPENDENCY"
    return "B_MODERATE_DEPENDENCY_CONCERN"


def build_dependency_table(cfg: dict, release: str) -> pd.DataFrame:
    model = load_model(cfg, release)
    effect = load_gene_effect(cfg, release, CANDIDATES)
    thresh = cfg["independent_validation"]["depmap"]["strong_dependency_probability_threshold"]
    hany = pd.read_csv(FROZEN_HANY, sep="\t").set_index("gene")

    has_prob = has_dependency_probability(cfg, release)
    dep_prob = load_gene_dependency(cfg, release, CANDIDATES) if has_prob else None
    if not has_prob:
        logger.warning(
            "build_dependency_table (%s): CRISPRGeneDependency.csv not available for this release -- "
            "frac_strongly_dependent_* and essentiality_concern will be UNAVAILABLE (E_INSUFFICIENT_DATA), "
            "not estimated from gene_effect or any other quantity, per explicit instruction not to fabricate probabilities.",
            release,
        )

    breast_ids = model.index[model["is_breast"]]
    luminal_ids = model.index[model["is_er_luminal"]]
    logger.info("build_dependency_table: %d breast lines, %d ER+/luminal breast lines (dependency threshold p>%.2f)", len(breast_ids), len(luminal_ids), thresh)

    rows = []
    for candidate in CANDIDATES:
        e_all = effect[candidate].dropna()
        e_breast = effect.loc[effect.index.isin(breast_ids), candidate].dropna()
        e_luminal = effect.loc[effect.index.isin(luminal_ids), candidate].dropna()

        if has_prob:
            p_all = dep_prob[candidate].dropna()
            p_breast = dep_prob.loc[dep_prob.index.isin(breast_ids), candidate].dropna()
            p_luminal = dep_prob.loc[dep_prob.index.isin(luminal_ids), candidate].dropna()
            frac_all = float((p_all > thresh).mean())
            frac_breast = float((p_breast > thresh).mean()) if len(p_breast) else np.nan
            frac_luminal = float((p_luminal > thresh).mean()) if len(p_luminal) else np.nan
            frac_note = f"{frac_all*100:.0f}% of lines strongly dependent (p>{thresh})."
        else:
            frac_all = frac_breast = frac_luminal = np.nan
            frac_note = "strongly-dependent fraction UNAVAILABLE (CRISPRGeneDependency.csv not obtained for this release -- not estimated from any other quantity)."

        top_dependent = e_breast.sort_values().head(5)
        top_dependent_named = "; ".join(f"{model.loc[i, 'StrippedCellLineName']}={v:.2f}" for i, v in top_dependent.items())

        concern = _classify_concern(frac_all, frac_breast)
        hany_row = hany.loc[candidate]
        rows.append(dict(
            candidate=candidate,
            depmap_release=release,
            median_gene_effect_all_cancer=float(e_all.median()), iqr_all_cancer=f"[{e_all.quantile(.25):.3f},{e_all.quantile(.75):.3f}]", n_all_cancer=len(e_all),
            median_gene_effect_breast=float(e_breast.median()) if len(e_breast) else np.nan, iqr_breast=f"[{e_breast.quantile(.25):.3f},{e_breast.quantile(.75):.3f}]" if len(e_breast) else "", n_breast=len(e_breast),
            median_gene_effect_er_luminal=float(e_luminal.median()) if len(e_luminal) else np.nan, iqr_er_luminal=f"[{e_luminal.quantile(.25):.3f},{e_luminal.quantile(.75):.3f}]" if len(e_luminal) else "", n_er_luminal=len(e_luminal),
            frac_strongly_dependent_all_cancer=frac_all, frac_strongly_dependent_breast=frac_breast, frac_strongly_dependent_er_luminal=frac_luminal,
            dependency_probability_available=has_prob,
            most_dependent_breast_lines=top_dependent_named,
            essentiality_concern=concern,
            hany_crispr_effect=float(hany_row["crispr_effect"]), hany_crispr_fdr=float(hany_row["crispr_fdr"]), hany_crispr_direction=hany_row["crispr_direction"],
            hany_vs_depmap_comparison=(
                f"Hany (drug-context, E2+4-OHT vs E2, one competitive screen): {hany_row['crispr_direction']}, "
                f"effect={hany_row['crispr_effect']:.2f}, FDR={hany_row['crispr_fdr']:.3f}. "
                f"DepMap (baseline, unconditioned standard culture, {len(e_all)} lines): "
                f"median gene effect={e_all.median():.3f}, {frac_note} "
                "These are NOT the same experiment -- a concordant sign is not replication of the drug-context effect, "
                "only evidence that the gene also has some baseline growth-fitness role."
            ),
        ))
    out = pd.DataFrame(rows)
    logger.info("build_dependency_table (%s): %d rows", release, len(out))
    return out


def build_expression_table(cfg: dict, release: str) -> pd.DataFrame:
    model = load_model(cfg, release)
    expr = load_expression(cfg, release, CANDIDATES)
    effect = load_gene_effect(cfg, release, CANDIDATES)
    breast_ids = model.index[model["is_breast"]]
    luminal_ids = model.index[model["is_er_luminal"]]

    rows = []
    for candidate in CANDIDATES:
        e_breast = expr.loc[expr.index.isin(breast_ids), candidate].dropna()
        e_luminal = expr.loc[expr.index.isin(luminal_ids), candidate].dropna()
        joined = expr.loc[expr.index.isin(breast_ids), [candidate]].join(effect.loc[effect.index.isin(breast_ids), [candidate]], lsuffix="_expr", rsuffix="_effect").dropna()
        if len(joined) >= 5:
            rho, p = stats.spearmanr(joined[f"{candidate}_expr"], joined[f"{candidate}_effect"])
        else:
            rho, p = np.nan, np.nan
        rows.append(dict(
            candidate=candidate,
            depmap_release=release,
            n_breast_lines_expressed=len(e_breast), median_tpm_log2p1_breast=float(e_breast.median()) if len(e_breast) else np.nan,
            n_er_luminal_lines_expressed=len(e_luminal), median_tpm_log2p1_er_luminal=float(e_luminal.median()) if len(e_luminal) else np.nan,
            fraction_breast_lines_expressed_gt1=float((e_breast > 1.0).mean()) if len(e_breast) else np.nan,
            spearman_rho_expression_vs_dependency_breast=float(rho), p_value=float(p), n_paired=len(joined),
            notes="secondary/exploratory analysis; correlation between expression and gene-effect across breast lines, not causal",
        ))
    out = pd.DataFrame(rows)
    testable = out["p_value"].notna()
    out["fdr"] = np.nan
    if testable.sum():
        out.loc[testable, "fdr"] = multipletests(out.loc[testable, "p_value"], method="fdr_bh")[1]
    return out


def run(config_path: str = "config/config.yaml", release: str | None = None) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Build and write the primary DepMap dependency/expression tables.

    ``release`` defaults to config's independent_validation.depmap.
    active_release (currently "26Q1", manually verified -- see
    DEPMAP_26Q1_ACCESS_STATUS.md; "24Q4" remains fully available for
    comparison). Pass an explicit release to build a specific release's
    table without touching the active/primary output files (used by the
    comparison module to build both releases' tables side by side).
    """
    cfg = load_config(config_path)
    release = release or cfg["independent_validation"]["depmap"]["active_release"]
    dep = build_dependency_table(cfg, release)
    OUT_DEPENDENCY.parent.mkdir(parents=True, exist_ok=True)
    dep.to_csv(OUT_DEPENDENCY, sep="\t", index=False)
    logger.info("wrote %s (%d rows, release=%s)", OUT_DEPENDENCY, len(dep), release)

    expr = build_expression_table(cfg, release)
    expr.to_csv(OUT_EXPRESSION, sep="\t", index=False)
    logger.info("wrote %s (%d rows, release=%s)", OUT_EXPRESSION, len(expr), release)
    return dep, expr


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run()
