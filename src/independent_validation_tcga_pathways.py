"""Independent-validation Part 4+5: TCGA-BRCA candidate-pathway association
and the VEZF1-CITED2 human-breast-cancer consistency check.

Candidate-pathway pairs are declared in config/config.yaml's
independent_validation.tcga.pathways.candidate_pathways BEFORE this module
is run (EML5 gets none: no independently-justified pathway exists per the
project's own systems-network + literature review, and none is invented
here). Pathway activity is scored per ER+ primary tumor sample by ssGSEA
(gseapy), then Spearman-correlated against that candidate's own log2(TPM+1)
expression within the same ER+ tumors. Association language only --
correlation is never described as "regulates."

Part 5 tests whether VEZF1 and CITED2 expression are correlated in human
TCGA-BRCA tumors, framed explicitly as a "human breast-cancer consistency
check" (not mechanistic validation) for the external-literature finding
that VEZF1 represses CITED2 in a non-cancer system.
"""

from __future__ import annotations

import logging
from pathlib import Path

import gseapy as gp
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

from src.independent_validation_tcga_data import build_cohort_table, load_config, load_expression, load_expression_symbols

logger = logging.getLogger(__name__)

OUT_TABLE = Path("results/tables/independent_validation/TCGA_candidate_pathway_associations.tsv")


def _load_gmt(path: Path) -> dict[str, list[str]]:
    sets = {}
    for line in Path(path).read_text().splitlines():
        parts = line.rstrip("\n").split("\t")
        sets[parts[0]] = parts[2:]
    return sets


def _gene_sets_for_candidates(cfg: dict, needed: set[str]) -> dict[str, list[str]]:
    tcga_cfg = cfg["independent_validation"]["tcga"]["pathways"]
    hallmark = _load_gmt(tcga_cfg["hallmark_gmt"])
    go_bp = _load_gmt(tcga_cfg["go_bp_gmt"])
    out = {}
    for name in needed:
        if name in hallmark:
            out[name] = hallmark[name]
        elif name in go_bp:
            out[name] = go_bp[name]
        else:
            raise KeyError(f"pathway {name} not found in hallmark.gmt or go_bp.gmt")
    return out


def build_pathway_scores(cfg: dict) -> tuple[pd.DataFrame, pd.DataFrame]:
    cand_pathways: dict[str, list[str]] = cfg["independent_validation"]["tcga"]["pathways"]["candidate_pathways"]
    needed = {p for plist in cand_pathways.values() for p in plist}
    gene_sets = _gene_sets_for_candidates(cfg, needed)

    cohort = build_cohort_table(cfg)
    er_pos_tumors = cohort.loc[cohort["is_primary_tumor"] & (cohort["ER_STATUS"] == "Positive")].index.tolist()
    logger.info("build_pathway_scores: %d ER+ primary tumors for pathway scoring", len(er_pos_tumors))

    expr_symbols = load_expression_symbols(cfg, sample_barcodes=er_pos_tumors)  # samples x genes
    ss = gp.ssgsea(
        data=expr_symbols.T,  # genes x samples, as gseapy expects
        gene_sets=gene_sets,
        sample_norm_method="rank",
        min_size=5,
        max_size=1000,
        threads=8,
        outdir=None,
        no_plot=True,
        seed=20260814,
    )
    scores = ss.res2d.pivot_table(index="Name", columns="Term", values="NES", aggfunc="first")
    scores.index.name = "sample_barcode"
    scores = scores.apply(pd.to_numeric)
    logger.info("build_pathway_scores: ssGSEA NES matrix %d samples x %d gene sets", *scores.shape)
    return scores, cohort.loc[er_pos_tumors]


def build_association_table(cfg: dict) -> pd.DataFrame:
    cand_pathways: dict[str, list[str]] = cfg["independent_validation"]["tcga"]["pathways"]["candidate_pathways"]
    all_candidates = cfg["independent_validation"]["candidates"]

    rows = []
    scores = None
    for candidate in all_candidates:
        pathways = cand_pathways.get(candidate, [])
        if not pathways:
            rows.append(dict(
                candidate=candidate, pathway="NONE", spearman_rho=float("nan"), p_value=float("nan"),
                n=0, notes="no independently-justified candidate-specific pathway declared for this gene "
                           "(per config/config.yaml's independent_validation.tcga.pathways.candidate_pathways) "
                           "-- no pathway test invented post hoc",
            ))
            continue
        if scores is None:
            scores, er_cohort = build_pathway_scores(cfg)
            candidate_expr = load_expression(cfg, genes=all_candidates)
        for pathway in pathways:
            joined = scores[[pathway]].join(candidate_expr[[candidate]], how="inner").dropna()
            rho, p = stats.spearmanr(joined[pathway], joined[candidate])
            rows.append(dict(
                candidate=candidate, pathway=pathway, spearman_rho=float(rho), p_value=float(p),
                n=len(joined), notes="ER+ primary tumors only; ssGSEA NES vs log2(TPM+1); association only, not regulation",
            ))

    out = pd.DataFrame(rows)
    testable = out["p_value"].notna()
    out["fdr"] = float("nan")
    out.loc[testable, "fdr"] = multipletests(out.loc[testable, "p_value"], method="fdr_bh")[1]
    return out


def build_vezf1_cited2_check(cfg: dict) -> pd.DataFrame:
    cohort = build_cohort_table(cfg)
    expr = load_expression(cfg, genes=["VEZF1", "CITED2"])
    df = cohort.join(expr)
    primary = df.loc[df["is_primary_tumor"]]

    groups = {
        "all_TCGA-BRCA_primary_tumors": primary,
        "ER+": primary.loc[primary["ER_STATUS"] == "Positive"],
        "Luminal_A": primary.loc[primary["PAM50_SUBTYPE"] == "Luminal A"],
        "Luminal_B": primary.loc[primary["PAM50_SUBTYPE"] == "Luminal B"],
    }
    rows = []
    for name, g in groups.items():
        sub = g[["VEZF1", "CITED2"]].dropna()
        if len(sub) < 3:
            rows.append(dict(group=name, n=len(sub), spearman_rho=float("nan"), p_value=float("nan"), notes="insufficient N"))
            continue
        rho, p = stats.spearmanr(sub["VEZF1"], sub["CITED2"])
        rows.append(dict(group=name, n=len(sub), spearman_rho=float(rho), p_value=float(p), notes="human breast-cancer consistency check ONLY -- not mechanistic validation of the literature-reported repression, and a null or opposite result here does not invalidate that cross-system literature finding; groups are NESTED/overlapping subsets of the same tumors, not independent replications"))
    out = pd.DataFrame(rows)
    testable = out["p_value"].notna()
    out["fdr"] = float("nan")
    if testable.sum():
        out.loc[testable, "fdr"] = multipletests(out.loc[testable, "p_value"], method="fdr_bh")[1]
    return out


def run(config_path: str = "config/config.yaml", out_table: Path = OUT_TABLE) -> pd.DataFrame:
    cfg = load_config(config_path)
    assoc = build_association_table(cfg)
    out_table.parent.mkdir(parents=True, exist_ok=True)
    assoc.to_csv(out_table, sep="\t", index=False)
    logger.info("wrote %s (%d rows)", out_table, len(assoc))

    vezf1_cited2 = build_vezf1_cited2_check(cfg)
    vc_path = out_table.parent / "TCGA_VEZF1_CITED2_consistency_check.tsv"
    vezf1_cited2.to_csv(vc_path, sep="\t", index=False)
    logger.info("wrote %s (%d rows)", vc_path, len(vezf1_cited2))
    return assoc


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run()
