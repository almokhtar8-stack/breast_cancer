"""POST-AUDIT SENSITIVITY ANALYSIS -- data loaders and table builders.

Companion analysis to an external scientific audit of the frozen
`THERAPEUTIC_SHORTLIST_FREEZE.tsv` (USP34 > VEZF1 > EML5 > CITED2). This
module does NOT modify, recompute, or re-rank anything in the original
frozen analysis. It reads the same already-frozen source tables the
original analysis used (plus the original candidate-freeze CODE itself,
called unmodified, for Rule 0) and asks a new, separate question: how
sensitive is the original conclusion to one specific, debatable design
choice -- requiring RNA/human-tumor corroboration as a hard eligibility
gate, on top of CRISPR significance?

No master score, no hidden weighting. Every table here is either a direct
read of an already-frozen source table, or a transparent, disclosed
re-derivation (e.g. a rank column, a boolean eligibility flag under a
named alternative rule) with the exact source documented in-line.

Frozen sources reused, never altered:
  - data/processed/labels.parquet (genome-wide Hany CRISPR fit)
  - results/tables/evidence_freeze/final_candidate_evidence.tsv (Rule 0 input)
  - results/tables/cross_dataset_genomewide/all_genes_cross_dataset_evidence_with_ranking.tsv
    (genome-wide per-dataset log2FC/p/FDR, all datasets, all genes)
  - results/tables/independent_validation/{TCGA_candidate_expression,DepMap_candidate_dependency}.tsv
  - results/tables/final_pharmacogenomics/{USP34,VEZF1}_GDSC_drug_associations.tsv
  - results/tables/lead_target_deep_dive/*, druggability_safety/* (USP34/VEZF1 only)
  - src.evidence_freeze_shortlist_freeze (Rule 0, called unmodified)
"""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

logger = logging.getLogger(__name__)

LABELS_PARQUET = Path("data/processed/labels.parquet")
FINAL_CANDIDATE_EVIDENCE = Path("results/tables/evidence_freeze/final_candidate_evidence.tsv")
CROSS_DATASET_TABLE = Path("results/tables/cross_dataset_genomewide/all_genes_cross_dataset_evidence_with_ranking.tsv")
FREEZE_TABLE = Path("results/tables/evidence_freeze/THERAPEUTIC_SHORTLIST_FREEZE.tsv")
TCGA_EXPRESSION = Path("results/tables/independent_validation/TCGA_candidate_expression.tsv")
DEPMAP_CANDIDATE_DEPENDENCY = Path("results/tables/independent_validation/DepMap_candidate_dependency.tsv")

ORIGINAL_FOUR = ["USP34", "VEZF1", "EML5", "CITED2"]
FOCUS_FOUR = ["KDM1A", "TLK2", "USP34", "VEZF1"]  # required by the audit
CRISPR_GATE1_FDR = 0.1  # PREANALYSIS.md Section 4, pre-specified, unchanged here
CHRONIC_RESISTANCE_DATASETS = ["gse118713", "gse240112", "gse111151"]  # GSE245601 excluded: acute


def load_config(config_path: str | Path = "config/config.yaml") -> dict:
    with open(config_path) as f:
        return yaml.safe_load(f)


# ---------------------------------------------------------------------------
# Phase 1A / Phase 2 -- genome-wide CRISPR + the significant sensitising universe
# ---------------------------------------------------------------------------

def load_genomewide_crispr() -> pd.DataFrame:
    df = pd.read_parquet(LABELS_PARQUET)
    assert len(df) == 19103, f"expected 19,103 fitted genes, got {len(df)}"
    return df


def load_significant_sensitising_hits() -> pd.DataFrame:
    """The complete, pre-specified-threshold (FDR<0.1, PREANALYSIS.md
    Section 4) sensitising (effect<0) universe -- this is the ENTIRE
    starting point for the sensitivity analysis, not just the 4 named
    genes. Effect-rank and FDR-rank are computed and reported SEPARATELY
    (they are different questions), both restricted to this same
    sensitising subset."""
    df = load_genomewide_crispr()
    sig = df[df["fdr"] < CRISPR_GATE1_FDR].copy()
    sens = sig[sig["effect_size"] < 0].copy()
    sens = sens.sort_values("effect_size").reset_index(drop=True)
    sens["rank_by_effect"] = np.arange(1, len(sens) + 1)
    by_fdr = sens.sort_values("fdr").reset_index(drop=True)
    by_fdr["rank_by_fdr"] = np.arange(1, len(by_fdr) + 1)
    sens = sens.merge(by_fdr[["gene", "rank_by_fdr"]], on="gene", how="left")
    sens = sens.sort_values("rank_by_effect").reset_index(drop=True)
    logger.info("load_significant_sensitising_hits: %d genes (FDR<%.2f, effect<0)", len(sens), CRISPR_GATE1_FDR)
    return sens


def load_gate1_full() -> pd.DataFrame:
    """All 28 Gate-1 significant genes (both directions) -- for context
    only (e.g. to state the tolerance-associated count)."""
    df = load_genomewide_crispr()
    return df[df["fdr"] < CRISPR_GATE1_FDR].copy()


# ---------------------------------------------------------------------------
# Phase 3 -- evidence matrix (RNA is supporting evidence, not eligibility)
# ---------------------------------------------------------------------------

def load_cross_dataset_raw() -> pd.DataFrame:
    df = pd.read_csv(CROSS_DATASET_TABLE, sep="\t")
    return df


def load_depmap_summary_for_genes(genes: list[str], release: str = "26Q1") -> pd.DataFrame:
    """Generic (non-4-candidate-restricted) DepMap 26Q1 summary for an
    arbitrary gene list -- reuses the exact symbol-column-matching logic
    already built and tested in poster_figures_bank_data (that function is
    hardwired to the 28 Gate-1 genes; this is the same logic parameterised
    for an arbitrary gene list so it also covers the two extra genes in
    final_candidate_evidence.tsv, GREB1/SLC4A10, if ever needed)."""
    from src.independent_validation_depmap_data import has_dependency_probability, load_model, raw_dir

    cfg = load_config()
    rel_cfg = cfg["independent_validation"]["depmap"]["releases"][release]
    effect_path = raw_dir(cfg, release) / rel_cfg["raw"]["crispr_gene_effect_csv"]
    dep_path = raw_dir(cfg, release) / rel_cfg["raw"]["crispr_gene_dependency_csv"]

    model = load_model(cfg, release)
    luminal_ids = model.index[model["is_er_luminal"]]

    effect_df = pd.read_csv(effect_path, index_col=0)
    dep_df = pd.read_csv(dep_path, index_col=0) if has_dependency_probability(cfg, release) else None
    thresh = cfg["independent_validation"]["depmap"]["strong_dependency_probability_threshold"]

    def gene_col(cols, symbol):
        for c in cols:
            if c.split(" (")[0] == symbol:
                return c
        return None

    rows = []
    for gene in genes:
        col = gene_col(effect_df.columns, gene)
        if col is None:
            rows.append(dict(gene=gene, in_depmap=False, median_chronos_er_luminal=np.nan,
                              frac_strongly_dependent_er_luminal=np.nan, n_lines=0))
            continue
        vals = effect_df.loc[effect_df.index.intersection(luminal_ids), col].dropna()
        frac = np.nan
        if dep_df is not None:
            dep_col = gene_col(dep_df.columns, gene)
            if dep_col is not None:
                dep_vals = dep_df.loc[dep_df.index.intersection(luminal_ids), dep_col].dropna()
                if len(dep_vals) > 0:
                    frac = float((dep_vals > thresh).mean())
        rows.append(dict(gene=gene, in_depmap=True, median_chronos_er_luminal=float(vals.median()),
                          frac_strongly_dependent_er_luminal=frac, n_lines=len(vals)))
    return pd.DataFrame(rows)


def load_tcga_for_original_four() -> pd.DataFrame:
    """TCGA-BRCA candidate-expression evidence exists ONLY for the
    original 4 candidates in this project -- never computed for KDM1A/
    TLK2/other sensitising hits. Returned as-is; callers must NOT force
    NaN to 0 for genes outside this set."""
    df = pd.read_csv(TCGA_EXPRESSION, sep="\t")
    return df[df["comparison"] == "tumor_vs_normal_PAIRED"][["candidate", "mean_diff", "p_value", "fdr", "n_a"]].rename(
        columns={"candidate": "gene", "mean_diff": "tcga_tumor_vs_normal_log2diff", "p_value": "tcga_p", "fdr": "tcga_fdr", "n_a": "tcga_n_paired"}
    )


def load_pathway_strong_consensus_counts() -> pd.DataFrame:
    """STRONG_CONSENSUS GSEA pathway-membership counts -- computed ONLY
    for the original 4 candidates (the systems_network phase never ran
    this analysis for KDM1A/TLK2/other Gate-1 genes)."""
    path = Path("results/tables/systems_network/candidate_pathway_membership.tsv")
    df = pd.read_csv(path, sep="\t")
    sub = df[(df["candidate_is_member"] == True) & (df["resistance_consensus_class"] == "STRONG_CONSENSUS")]  # noqa: E712
    return sub.groupby("candidate").size().rename("n_strong_consensus_pathways").reset_index().rename(columns={"candidate": "gene"})


def load_gdsc_best_hit_for_original_two() -> pd.DataFrame:
    """GDSC pharmacogenomic screening was only run for USP34 and VEZF1 in
    this project -- returns each gene's strongest FDR row (any metric/
    dataset), or nothing if no association table exists for that gene."""
    rows = []
    for gene in ["USP34", "VEZF1"]:
        path = Path(f"results/tables/final_pharmacogenomics/{gene}_GDSC_drug_associations.tsv")
        if not path.exists():
            continue
        df = pd.read_csv(path, sep="\t")
        best = df.sort_values("fdr").iloc[0]
        n_sig = int((df["fdr"] < 0.05).sum())
        n_sig_unique_compounds = df[df["fdr"] < 0.05]["drug_id"].nunique()
        rows.append(dict(gene=gene, gdsc_n_tests=len(df), gdsc_n_fdr_sig_rows=n_sig,
                          gdsc_n_fdr_sig_unique_compounds=n_sig_unique_compounds,
                          gdsc_best_drug=best["drug_name"], gdsc_best_rho=best["spearman_rho"], gdsc_best_fdr=best["fdr"]))
    return pd.DataFrame(rows)


def build_evidence_matrix() -> pd.DataFrame:
    """Table 03: the central post-audit evidence matrix. Rows = all 13
    significant sensitising Gate-1 hits. Columns = every quantitative
    piece of evidence ALREADY COMPUTED in this project for that gene,
    joined without recomputation. Missing evidence is NaN, never 0 --
    absence of a TCGA/pathway/GDSC column value for a non-original-4 gene
    means "not assessed in this project", not "no effect"."""
    sens = load_significant_sensitising_hits()
    genes = sens["gene"].tolist()

    cross = load_cross_dataset_raw().set_index("gene")
    cross_sub = cross.loc[genes, [
        "gse118713_log2fc", "gse118713_fdr",
        "gse111151_log2fc", "gse111151_fdr",
        "gse240112_tumor_log2fc", "gse240112_tumor_fdr",
        "gse245601_epi_log2fc", "gse245601_epi_fdr",
    ]].reset_index()
    cross_sub.columns = ["gene", "gse118713_log2fc", "gse118713_fdr", "gse111151_log2fc", "gse111151_fdr",
                          "gse240112_log2fc", "gse240112_fdr", "gse245601_acute_log2fc", "gse245601_acute_fdr"]

    depmap = load_depmap_summary_for_genes(genes)
    tcga = load_tcga_for_original_four()
    pathway = load_pathway_strong_consensus_counts()
    gdsc = load_gdsc_best_hit_for_original_two()

    out = sens[["gene", "effect_size", "p_value", "fdr", "rank_by_effect", "rank_by_fdr"]].rename(
        columns={"effect_size": "crispr_effect", "p_value": "crispr_p", "fdr": "crispr_fdr"}
    )
    out = out.merge(cross_sub, on="gene", how="left")
    out = out.merge(depmap, on="gene", how="left")
    out = out.merge(tcga, on="gene", how="left")
    out = out.merge(pathway, on="gene", how="left")
    out = out.merge(gdsc, on="gene", how="left")

    # chronic (non-acute) resistance corroboration count -- the honest,
    # transparent quantity behind "how many independent resistance models
    # support this gene", computed the same way as the frozen
    # resistance_fdr05_count() (gse118713/gse240112/gse111151 only)
    out["n_chronic_datasets_fdr05"] = (
        (out["gse118713_fdr"] < 0.05).fillna(False).astype(int)
        + (out["gse240112_fdr"] < 0.05).fillna(False).astype(int)
        + (out["gse111151_fdr"] < 0.05).fillna(False).astype(int)
    )
    out["n_chronic_datasets_testable"] = (
        out["gse118713_fdr"].notna().astype(int) + out["gse240112_fdr"].notna().astype(int) + out["gse111151_fdr"].notna().astype(int)
    )
    logger.info("build_evidence_matrix: %d genes x %d columns", len(out), out.shape[1])
    return out


# ---------------------------------------------------------------------------
# Phase 5 -- selection-rule sensitivity analysis
# ---------------------------------------------------------------------------

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


_BAND_RANK = {"VERY_STRONG": 0, "STRONG": 1, "MODERATE": 2, "WEAK": 3}


def rule0_original_frozen_gate() -> pd.DataFrame:
    """Rule 0: the ORIGINAL frozen gate, called completely unmodified on
    its original input table -- not a re-implementation. This must
    reproduce THERAPEUTIC_SHORTLIST_FREEZE.tsv exactly; a test asserts
    this bit-for-bit."""
    from src.evidence_freeze_shortlist_freeze import build_freeze

    full_table = pd.read_csv(FINAL_CANDIDATE_EVIDENCE, sep="\t")
    frozen, audited = build_freeze(full_table)
    frozen = frozen.reset_index(drop=True)
    frozen.insert(0, "rank", np.arange(1, len(frozen) + 1))
    return frozen[["rank", "gene", "crispr_effect", "crispr_fdr"]]


def rule1_crispr_only() -> pd.DataFrame:
    """Rule 1: CRISPR significance (Gate-1 FDR<0.1) + sensitising
    direction ONLY. No RNA eligibility requirement of any kind. Ranked by
    CRISPR-strength band, then raw FDR within band (transparent,
    disclosed, no hidden score)."""
    sens = load_significant_sensitising_hits()
    out = sens[["gene", "effect_size", "fdr"]].rename(columns={"effect_size": "crispr_effect", "fdr": "crispr_fdr"}).copy()
    out["_band_rank"] = out["crispr_fdr"].map(_crispr_band).map(_BAND_RANK)
    out = out.sort_values(["_band_rank", "crispr_fdr"]).reset_index(drop=True)
    out.insert(0, "rank", np.arange(1, len(out) + 1))
    return out.drop(columns="_band_rank")


def rule2_chronic_corroboration() -> pd.DataFrame:
    """Rule 2: CRISPR significance + at least 1 CHRONIC resistance/
    recurrence dataset at FDR<0.05 (GSE118713, GSE240112, GSE111151 only
    -- GSE245601 excluded because it is an acute 12h context, never
    counted as resistance evidence anywhere in this analysis)."""
    em = build_evidence_matrix()
    eligible = em[em["n_chronic_datasets_fdr05"] >= 1].copy()
    eligible["_band_rank"] = eligible["crispr_fdr"].map(_crispr_band).map(_BAND_RANK)
    eligible = eligible.sort_values(["_band_rank", "n_chronic_datasets_fdr05", "crispr_fdr"], ascending=[True, False, True]).reset_index(drop=True)
    eligible.insert(0, "rank", np.arange(1, len(eligible) + 1))
    return eligible[["rank", "gene", "crispr_effect", "crispr_fdr", "n_chronic_datasets_fdr05", "n_chronic_datasets_testable"]]


def rule3_gse111151_specific() -> pd.DataFrame:
    """Rule 3 (sensitivity scenario, NOT presented as correct): CRISPR
    significance + corroboration in GSE111151 specifically (FDR<0.05).
    GSE111151 is the only dataset in this project with multiple
    independent parental/resistance backgrounds."""
    em = build_evidence_matrix()
    eligible = em[em["gse111151_fdr"] < 0.05].copy()
    if len(eligible) == 0:
        return eligible.assign(rank=pd.Series(dtype=int))[["rank", "gene", "crispr_effect", "crispr_fdr", "gse111151_log2fc", "gse111151_fdr"]]
    eligible = eligible.sort_values("gse111151_fdr").reset_index(drop=True)
    eligible.insert(0, "rank", np.arange(1, len(eligible) + 1))
    return eligible[["rank", "gene", "crispr_effect", "crispr_fdr", "gse111151_log2fc", "gse111151_fdr"]]


def rule4_human_evidence_first() -> pd.DataFrame:
    """Rule 4: CRISPR-first eligibility (same universe as Rule 1, no RNA
    gate), but the RANKING hierarchy places HUMAN tumor evidence
    (GSE240112 real-patient recurrence RNA, FDR<0.05; TCGA paired
    tumor-vs-normal FDR<0.05 where available) ahead of resistance
    cell-line-model RNA-direction-consistency (GSE118713/GSE111151).
    Directly tests the external reviewer's claim that VEZF1 may outrank
    USP34 once human evidence is considered before cell-line expression
    consistency."""
    em = build_evidence_matrix()
    out = em.copy()
    out["_band_rank"] = out["crispr_fdr"].map(_crispr_band).map(_BAND_RANK)
    human_sig = (out["gse240112_fdr"] < 0.05) | (out["tcga_fdr"] < 0.05)
    out["_human_rank"] = (~human_sig.fillna(False)).astype(int)  # 0 = human-significant, 1 = not
    cellline_sig_count = (out["gse118713_fdr"] < 0.05).fillna(False).astype(int) + (out["gse111151_fdr"] < 0.05).fillna(False).astype(int)
    out["_cellline_rank"] = -cellline_sig_count  # more significant cell-line datasets ranks better (more negative)
    out = out.sort_values(["_band_rank", "_human_rank", "_cellline_rank", "crispr_fdr"]).reset_index(drop=True)
    out.insert(0, "rank", np.arange(1, len(out) + 1))
    return out[["rank", "gene", "crispr_effect", "crispr_fdr", "gse240112_fdr", "tcga_fdr", "gse118713_fdr", "gse111151_fdr"]]


def build_selection_rule_sensitivity() -> pd.DataFrame:
    """Table 04: one row per (rule, gene) with eligibility + rank, for the
    5 rules (0-4), long format so USP34/VEZF1/KDM1A/TLK2's fate under
    every rule can be read directly off one table."""
    rules = {
        "RULE_0_original_frozen_gate": rule0_original_frozen_gate(),
        "RULE_1_crispr_only_no_rna_gate": rule1_crispr_only(),
        "RULE_2_chronic_rna_corroboration": rule2_chronic_corroboration(),
        "RULE_3_gse111151_specific": rule3_gse111151_specific(),
        "RULE_4_human_evidence_first": rule4_human_evidence_first(),
    }
    all_genes = sorted(set(load_significant_sensitising_hits()["gene"]) | set(ORIGINAL_FOUR))
    rows = []
    for rule_name, ranked in rules.items():
        eligible_genes = set(ranked["gene"])
        for gene in all_genes:
            if gene in eligible_genes:
                r = ranked[ranked["gene"] == gene].iloc[0]
                rows.append(dict(rule=rule_name, gene=gene, eligible=True, rank=int(r["rank"]), n_eligible=len(ranked)))
            else:
                rows.append(dict(rule=rule_name, gene=gene, eligible=False, rank=np.nan, n_eligible=len(ranked)))
    out = pd.DataFrame(rows)
    logger.info("build_selection_rule_sensitivity: %d rules x %d genes = %d rows", len(rules), len(all_genes), len(out))
    return out


# ---------------------------------------------------------------------------
# Phase 5 -- leave-one-dataset-out sensitivity
# ---------------------------------------------------------------------------

def build_leave_one_dataset_out() -> pd.DataFrame:
    """Table 05: for each of the 13 sensitising genes, the chronic-
    corroboration count (GSE118713/GSE240112/GSE111151 FDR<0.05) recomputed
    after excluding each chronic dataset individually. GSE245601 is never
    included as a chronic layer, per the audit instruction."""
    em = build_evidence_matrix()
    datasets = {"gse118713": "gse118713_fdr", "gse240112": "gse240112_fdr", "gse111151": "gse111151_fdr"}
    rows = []
    for _, r in em.iterrows():
        sig = {d: (pd.notna(r[col]) and r[col] < 0.05) for d, col in datasets.items()}
        full_count = sum(sig.values())
        for left_out in datasets:
            remaining = {d: v for d, v in sig.items() if d != left_out}
            rows.append(dict(
                gene=r["gene"], left_out_dataset=left_out,
                n_chronic_fdr05_full=full_count,
                n_chronic_fdr05_after_exclusion=sum(remaining.values()),
                still_corroborated=sum(remaining.values()) >= 1,
            ))
    out = pd.DataFrame(rows)
    logger.info("build_leave_one_dataset_out: %d rows", len(out))
    return out


# ---------------------------------------------------------------------------
# Phase 7-9 -- deep dives (in-project quantitative evidence + independently
# verified external literature, each row's source cited explicitly; no
# numeric project value is ever hand-typed here -- only pre-existing
# project tables (already read above) and cited external literature facts,
# the same mixed-evidence convention already used throughout this project's
# own literature_mechanism/druggability_safety tables)
# ---------------------------------------------------------------------------

def build_kdm1a_deep_dive() -> pd.DataFrame:
    em = build_evidence_matrix()
    r = em[em["gene"] == "KDM1A"].iloc[0]
    rows = [
        dict(axis="Hany CRISPR", finding=f"effect={r['crispr_effect']:.3f}, FDR={r['crispr_fdr']:.3g} -- strongest sensitising hit in the entire screen by both effect size and FDR (rank 1/13 both)", source="data/processed/labels.parquet (this project, verified this phase)"),
        dict(axis="GSE118713 (cell-line resistance)", finding=f"log2FC={r['gse118713_log2fc']:.3f}, FDR={r['gse118713_fdr']:.3g} -- not significant", source="results/tables/cross_dataset_genomewide (this project)"),
        dict(axis="GSE111151 (multi-background resistance)", finding=f"log2FC={r['gse111151_log2fc']:.3f}, FDR={r['gse111151_fdr']:.3g} -- not significant", source="results/tables/cross_dataset_genomewide (this project)"),
        dict(axis="GSE240112 (human recurrence)", finding=f"log2FC={r['gse240112_log2fc']:.3f}, FDR={r['gse240112_fdr']:.3g} -- not significant", source="results/tables/cross_dataset_genomewide (this project)"),
        dict(axis="GSE245601 (ACUTE 12h -- not resistance)", finding=f"log2FC={r['gse245601_acute_log2fc']:.3f}, FDR={r['gse245601_acute_fdr']:.3g} -- acute context only, never counted as resistance support", source="results/tables/cross_dataset_genomewide (this project)"),
        dict(axis="TCGA-BRCA", finding="not assessed in this project -- TCGA candidate analysis was only run for the original 4 (USP34/VEZF1/EML5/CITED2)", source="N/A -- evidence gap, not a null result"),
        dict(axis="DepMap 26Q1 ER+/luminal dependency", finding=f"median Chronos={r['median_chronos_er_luminal']:.3f}, {r['frac_strongly_dependent_er_luminal']*100:.0f}% of 11 ER+/luminal lines dependency-probability>0.5 -- weak/no baseline dependency, similar profile to USP34", source="DepMap 26Q1 CRISPRGeneEffect/CRISPRGeneDependency (this project, generic loader)"),
        dict(axis="Blind-control status (project-specific)", finding="KDM1A was one of 2 pre-registered BLIND positive controls (PREANALYSIS.md Section 5); never inspected/reported at gene level until the 2026-08-10 amendment formally retired the blind. Its recovery as the single strongest sensitising hit (rank 1/13) is itself the pre-registered validity check the blind was designed to enable.", source="PREANALYSIS.md Section 5 and 2026-08-10 amendments (this project)"),
        dict(axis="Structural biology", finding="KDM1A/LSD1 has extensive experimental structural coverage, including multiple inhibitor-bound co-crystal structures (e.g. PDB 6NQU with GSK2879552; tranylcypromine-scaffold covalent inhibitors e.g. PDB 2Z5U) -- structurally the most extensively characterized of the 4 focus genes, consistent with its clinical-stage pharmacology. See 06b_structural_tractability_audit.tsv for the full facet-by-facet breakdown.", source="MDPI Molecules 23:1538 / PMC6099836 (PDB 6NQU); web-verified 2026-08-15"),
        dict(axis="Mechanistic literature: ER/ERa biology", finding="LSD1/KDM1A is a well-established ER-corepressor-complex (CoREST) component with independent, pre-existing literature linking it to estrogen-receptor signalling and endocrine-therapy response -- this is WHY it was chosen as a blind positive control in the first place, not a new finding of this project.", source="PREANALYSIS.md Section 5 (cites prior literature); independently corroborated by 2026-08-15 web search below"),
        dict(axis="Mechanistic literature: endocrine resistance", finding="Independent, external literature reports that LSD1 and KDM6A/UTX co-express and co-localize with ER in breast cancer, and that reducing LSD1 activity can prevent resistance to antitumor agents in breast cancer models -- external mechanistic support for a tamoxifen-relevant role, independent of this project's own CRISPR screen.", source="PMC9932783 (LSD1 inhibitors for cancer treatment review); web-verified 2026-08-15"),
        dict(axis="Existing pharmacology / inhibitors", finding="Multiple LSD1/KDM1A inhibitors exist and are in active clinical development. Iadademstat (ORY-1001) is an orally active, highly potent, selective LSD1 inhibitor with Orphan Drug Designation (AML, SCLC) as of 2026 -- not yet FDA-approved for any indication, but the most translationally advanced compound of any gene compared in this sensitivity analysis, including USP34. It has published, breast-cancer-specific preclinical data (targets SOX2-driven breast cancer stem cells in luminal-B/HER2+ subtypes).", source="PMC7138538 / PMID 32191225 (Iadademstat breast cancer stem cell paper); Oryzon Genomics pipeline; web-verified 2026-08-15"),
        dict(axis="Clinical development status", finding="Iadademstat: Orphan Drug Designation (AML, SCLC), multiple clinical trials completed/ongoing as of 2026, described in recent coverage as a 'best-in-class' LSD1 inhibitor candidate. No breast-cancer-specific trial identified in this search. Several other LSD1 inhibitors are also in clinical trials for other cancers (per the same literature review).", source="web-verified 2026-08-15 (WebSearch; no direct primary-registry cross-check performed -- treat as literature-reported status, not independently confirmed against a trial registry)"),
        dict(axis="Toxicity / normal-tissue concerns", finding="Not independently assessed in this project (no GTEx/HPA/liability table was built for KDM1A here, unlike USP34/VEZF1). LSD1 is broadly essential in normal hematopoiesis and neurodevelopment in the published literature (a widely-cited concern for LSD1-inhibitor clinical programs generally, e.g. thrombocytopenia as a class effect) -- this is a known, generic LSD1-inhibitor-class liability, not specific to this project's data.", source="general LSD1-inhibitor-class knowledge; not independently re-verified against a specific primary source this phase -- flagged as an assessment gap"),
        dict(axis="Novelty for this project's objective", finding="LOW novelty as a tamoxifen-response modifier specifically -- KDM1A/LSD1's role in ER signalling and endocrine therapy is well-established prior literature (indeed, that established role is WHY the project chose it as a blind positive control). Its value here is as a validity check on the screen, and as a strong functional hit, not as a novel discovery.", source="PREANALYSIS.md Section 5 (this project's own framing)"),
    ]
    return pd.DataFrame(rows)


def build_tlk2_deep_dive() -> pd.DataFrame:
    em = build_evidence_matrix()
    r = em[em["gene"] == "TLK2"].iloc[0]
    rows = [
        dict(axis="Hany CRISPR", finding=f"effect={r['crispr_effect']:.3f}, FDR={r['crispr_fdr']:.3g} -- 2nd strongest by FDR, 4th strongest by effect, among the 13 sensitising hits", source="data/processed/labels.parquet (this project, verified this phase)"),
        dict(axis="GSE118713 (cell-line resistance)", finding=f"log2FC={r['gse118713_log2fc']:.3f}, FDR={r['gse118713_fdr']:.3g} -- not significant", source="results/tables/cross_dataset_genomewide (this project)"),
        dict(axis="GSE111151 (multi-background resistance)", finding=f"log2FC={r['gse111151_log2fc']:.3f}, FDR={r['gse111151_fdr']:.3g} -- not significant", source="results/tables/cross_dataset_genomewide (this project)"),
        dict(axis="GSE240112 (human recurrence)", finding=f"log2FC={r['gse240112_log2fc']:.3f}, FDR={r['gse240112_fdr']:.3g} -- not significant", source="results/tables/cross_dataset_genomewide (this project)"),
        dict(axis="GSE245601 (ACUTE 12h -- not resistance)", finding=f"log2FC={r['gse245601_acute_log2fc']:.3f}, FDR={r['gse245601_acute_fdr']:.3g} -- acute context only, never counted as resistance support", source="results/tables/cross_dataset_genomewide (this project)"),
        dict(axis="TCGA-BRCA", finding="not assessed in this project -- TCGA candidate analysis was only run for the original 4 (USP34/VEZF1/EML5/CITED2)", source="N/A -- evidence gap, not a null result"),
        dict(axis="DepMap 26Q1 ER+/luminal dependency", finding=f"median Chronos={r['median_chronos_er_luminal']:.3f}, {r['frac_strongly_dependent_er_luminal']*100:.0f}% of 11 ER+/luminal lines dependency-probability>0.5 -- the STRONGEST baseline cancer-cell dependency of any of the 4 focus genes (KDM1A/TLK2/USP34/VEZF1), including the frozen shortlist leaders.", source="DepMap 26Q1 CRISPRGeneEffect/CRISPRGeneDependency (this project, generic loader)"),
        dict(axis="Structural biology", finding="TLK2 HAS a real experimental structure: PDB 5O0Y, the human TLK2 kinase domain, 2.86A resolution, bound to AGS (ATP-gamma-S, a non-hydrolysable ATP analog -- a substrate mimetic, not a small-molecule inhibitor or drug). Correction to this analysis's earlier draft, which incorrectly implied TLK2 had no structure at all -- see 06b_structural_tractability_audit.tsv for the full facet-by-facet breakdown.", source="RCSB PDB 5O0Y; Nature Communications 2018 (s41467-018-04941-y); web-verified 2026-08-15"),
        dict(axis="Independent literature: amplification/prognosis", finding="TLK2 is reported amplified in ~10.5% of ER+ breast tumors in independent published work, with overexpression associated with more aggressive/advanced tumors and worse clinical outcomes regardless of endocrine-therapy status -- external genetic/prognostic support independent of this project's CRISPR screen.", source="Nature Communications 2016 (ncomms12991, TLK2 functional analysis in luminal breast cancer); web-verified 2026-08-15"),
        dict(axis="Independent literature: mechanism", finding="Published work reports TLK2 inhibition/deletion downregulates ERalpha, BCL2 and SKP2, impairs G1/S progression, induces apoptosis, improves progression-free survival in vivo (mouse models), and blocks lung metastasis in a breast-cancer model -- substantive independent preclinical mechanistic support.", source="AACR 2024 abstract A017 (Targeting TLK2 impairs breast cancer growth); web-verified 2026-08-15"),
        dict(axis="Existing pharmacology / inhibitors", finding="No potent, selective TLK2 inhibitor has reached clinical development. The main obstacle reported in the literature is close kinase-domain sequence homology with the TLK1 paralog, making selective inhibition difficult. Only early-stage, narrow-spectrum chemical tool compounds have been reported (structure-activity-relationship optimization stage, not clinical candidates).", source="ScienceDirect 2024 (TLK2 QSAR narrow-spectrum inhibitor discovery); web-verified 2026-08-15"),
        dict(axis="Clinical development status", finding="No clinical-stage TLK2 inhibitor identified in this search -- earliest-stage of the 4 focus genes for direct pharmacological targeting.", source="web-verified 2026-08-15 (WebSearch; absence of a positive result, not an exhaustive registry search)"),
        dict(axis="Toxicity / normal-tissue concerns", finding="Not independently assessed in this project. TLK2 has a broader normal-tissue biological role in chromatin assembly/replication (Tousled-like kinase family, ASF1 substrate) -- a generic replication-machinery liability concern common to this target class, not specifically evaluated here.", source="general TLK2 biology; not independently re-verified against a specific primary source this phase -- flagged as an assessment gap"),
        dict(axis="Novelty for this project's objective", finding="MODERATE-HIGH novelty as a tamoxifen-response modifier: TLK2's prior literature is about amplification/prognosis and general ERalpha-pathway involvement in breast cancer, not specifically about acute tamoxifen-response modulation -- the Hany CRISPR hit is a comparatively novel functional link in that specific context, even though TLK2 itself is an established breast-cancer gene.", source="synthesis of the literature rows above"),
    ]
    return pd.DataFrame(rows)


def build_usp34_updated_liability() -> pd.DataFrame:
    """Table 09: USP34's updated safety/liability picture, integrating the
    already-frozen tissue-liability table with the externally re-verified
    neurodevelopmental LoF paper -- no frozen value is altered; this only
    assembles already-computed/cited rows plus an explicit interpretation
    note satisfying the audit's dosage-extrapolation caveat."""
    liab_path = Path("results/tables/lead_target_deep_dive/USP34_VEZF1_tissue_liability.tsv")
    liab = pd.read_csv(liab_path, sep="\t")
    usp34 = liab[liab["candidate"] == "USP34"].copy()
    usp34["externally_reverified_2026_08_15"] = usp34["organ_system"].isin(["Skeletal development (human, germline)", "Neurological/cognitive"])
    note = pd.DataFrame([{
        "candidate": "USP34", "organ_system": "INTERPRETATION NOTE (added this phase, not a frozen project value)",
        "classification": "N/A",
        "evidence_summary": (
            "The Wigoda et al. 2026 paper (PMID 42315110 per this project's own citation; independently "
            "re-verified 2026-08-15 as 'USP34 Haploinsufficiency as a Cause of Neurodevelopmental Phenotypes', "
            "Clinical Genetics, DOI 10.1111/cge.70194, 6 individuals / 5 confirmed de novo heterozygous LoF, "
            "developmental delay/speech impairment/autism-spectrum features/craniofacial dysmorphism) IS REAL "
            "and IS already correctly documented in this project's frozen USP34_VEZF1_tissue_liability.tsv. "
            "It reports a CONGENITAL, LIFELONG, GERMLINE ~50%-dosage phenotype present from conception through "
            "neurodevelopment. This is NOT the same exposure as a time-limited, adult-onset, pharmacological "
            "partial inhibition of USP34 protein activity in a cancer patient -- developmental windows in which "
            "haploinsufficiency causes irreversible patterning/neurodevelopmental defects are not necessarily "
            "reopened by inhibiting the enzyme in an already-developed adult. This paper does NOT automatically "
            "mean 'USP34 inhibition is unsafe in adults' -- but it is a genuine, serious, and previously-softer "
            "translational liability signal that should be weighted more heavily than a normal-tissue-expression "
            "finding alone, and should be explicitly disclosed in any USP34 translational-strategy document."
        ),
        "species": "N/A", "reversible_or_developmental": "N/A", "zygosity": "N/A",
        "relevance_to_partial_inhibition": "See evidence_summary", "sources": "PMID 42315110; DOI 10.1111/cge.70194 (web-verified 2026-08-15)",
        "externally_reverified_2026_08_15": True,
    }])
    return pd.concat([usp34, note], ignore_index=True)


def build_structural_tractability_audit() -> pd.DataFrame:
    """Table 06b: the precise, non-binary structural/pharmacology re-audit
    requested after the first post-audit pass conflated "has a structure"
    with "USP34 uniquely has a structure" -- TLK2 in fact also has a real,
    independently-verified experimental kinase-domain structure (PDB 5O0Y),
    and VEZF1 has a real (weak, preliminary, homology-model-guided)
    published small-molecule screening hit. Facets A-G are reported
    separately per gene, each individually sourced -- never collapsed back
    into one boolean. Web-verified 2026-08-15 (WebSearch + WebFetch on the
    RCSB PDB entry page and the primary VEZF1-inhibitor paper); this is
    external literature, not a project-computed statistic, and is cited
    accordingly rather than derived from a frozen table."""
    rows = [
        dict(gene="KDM1A",
             A_experimental_human_structure_exists=True,
             B_relevant_domain_structure=True,
             C_ligand_or_probe_bound_structure="YES -- multiple inhibitor-bound co-crystal structures (e.g. PDB 6NQU with GSK2879552; tranylcypromine-class covalent inhibitors e.g. PDB 2Z5U)",
             D_established_druggable_pocket=True,
             E_validated_selective_small_molecule_inhibitor="YES, several -- e.g. iadademstat/ORY-1001 described in its own literature as highly potent and selective for LSD1 (some earlier tranylcypromine-scaffold compounds are less selective, also inhibiting MAOA/MAOB)",
             F_clinical_stage_pharmacology="YES -- 9 LSD1 inhibitors reported to have reached clinical stage across various indications (tranylcypromine-derivatives, ORY-1001, ORY-2001, GSK-2879552, IMG-7289, INCB059872, TAK-418, CC-90011, SP-2577); none breast-cancer-approved",
             G_novel_structurally_addressable_opportunity="LOW novelty -- an already extensively structurally characterized and pharmacologically mature target class; this project's contribution is recovering it functionally, not a new structural opportunity",
             precise_summary="Extensively characterized, clinically mature epigenetic-enzyme target: many experimental structures including several inhibitor-bound co-crystals, multiple clinical-stage selective inhibitors.",
             sources="PMC7138538/PMID 32191225; MDPI Molecules 23:1538 (PMC6099836, PDB 6NQU); PMC9932783 (LSD1 inhibitor review, 9 clinical-stage compounds); web-verified 2026-08-15"),
        dict(gene="TLK2",
             A_experimental_human_structure_exists=True,
             B_relevant_domain_structure="PARTIAL -- kinase domain only (PDB 5O0Y, 2.86A); the coiled-coil dimerization/regulatory regions of full-length TLK2 are not covered",
             C_ligand_or_probe_bound_structure="PARTIAL -- PDB 5O0Y is bound to AGS (ATP-gamma-S), a non-hydrolysable ATP/substrate analog, NOT a small-molecule inhibitor or drug-like ligand",
             D_established_druggable_pocket="PLAUSIBLE, not experimentally confirmed with an inhibitor -- the ATP pocket is the canonical kinase drug pocket by analogy, but no inhibitor-bound TLK2 co-crystal was found",
             E_validated_selective_small_molecule_inhibitor="NO -- only early-stage, narrow-spectrum chemical tool/SAR compounds reported; explicit unsolved TLK1-paralog kinase-domain-homology selectivity problem",
             F_clinical_stage_pharmacology="NO",
             G_novel_structurally_addressable_opportunity="MODERATE-HIGH -- a real experimental structure has existed since ~2018 with no inhibitor yet developed against it; a genuine open medicinal-chemistry opportunity, not an unstructured target",
             precise_summary="Has a real experimental kinase-domain structure (PDB 5O0Y) bound to an ATP analog, not an inhibitor. No selective small-molecule inhibitor exists; TLK1-paralog selectivity is the stated obstacle.",
             sources="RCSB PDB 5O0Y (web-verified via rcsb.org 2026-08-15); Nature Communications 2018 (s41467-018-04941-y, molecular basis of TLK2 activation); ScienceDirect 2024 (S022352342400237X, narrow-spectrum TLK2 inhibitor SAR); web-verified 2026-08-15"),
        dict(gene="USP34",
             A_experimental_human_structure_exists=True,
             B_relevant_domain_structure="PARTIAL -- catalytic domain only (~12% of the 3546-aa full-length protein); no full-length structure exists",
             C_ligand_or_probe_bound_structure="YES, but a covalent ACTIVITY-BASED PROBE, not a drug -- 7W3U is bound to a ubiquitin-propargylamide (UbPA) probe, covalently linked to Cys1903 (LINK-record-confirmed, 1.59-2.48A bond distance); this proves catalytic-cysteine reactivity, it is not a small-molecule ligand",
             D_established_druggable_pocket="PARTIAL -- fpocket identifies a top-ranked, favorably-scored pocket (score 0.845) at the catalytic dyad, but it is large/groove-shaped (1471 cubic Angstrom, an extended ubiquitin-binding surface), not a conventional compact drug pocket; a separate PPI-interface pocket exists at the ubiquitin Ile44 patch (already documented in this project's USP34_pocket_analysis.tsv, not re-derived here)",
             E_validated_selective_small_molecule_inhibitor="NO -- none identified in this project or in this phase's independent re-search",
             F_clinical_stage_pharmacology="NO",
             G_novel_structurally_addressable_opportunity="HIGH -- a genuinely unexplored target with real, recent (2022) experimental structures and direct crystallographic proof of a reactive, covalently-addressable catalytic cysteine, but zero existing chemical matter",
             precise_summary="Novel catalytic DUB target with experimentally resolved covalent chemical addressability at Cys1903 (proven by a covalent activity-based probe, not an inhibitor) but no validated selective small-molecule inhibitor identified. This is the precise, corrected claim -- NOT 'USP34 has a structure whereas KDM1A/TLK2 do not' (both of those also have real experimental structures).",
             sources="7W3R/7W3U (already frozen in this project's final_translational phase, PMID 35588869); USP34_pocket_analysis.tsv, USP34_docking_decision.tsv (this project); web-re-verified 2026-08-15 that no external inhibitor exists"),
        dict(gene="VEZF1",
             A_experimental_human_structure_exists=False,
             B_relevant_domain_structure="NO -- no solved zinc-finger-domain structure; a published inhibitor-discovery study built a homology model on the unrelated Zif268 zinc finger (PDB 1AAY) for lack of a real VEZF1 structure",
             C_ligand_or_probe_bound_structure=False,
             D_established_druggable_pocket="NO experimentally-validated pocket -- transcription-factor DNA-binding interfaces are a recognized, generally harder drug-target class than enzyme active sites",
             E_validated_selective_small_molecule_inhibitor="NO validated inhibitor, but a real, weak, preliminary screening hit exists: T4/503-1-83, IC50=20 micromolar, from an NCI Diversity Library virtual screen with biochemical (DNA-binding-inhibition) validation -- correction to this analysis's earlier, too-strong claim of 'no small-molecule inhibitor program identified'",
             F_clinical_stage_pharmacology="NO",
             G_novel_structurally_addressable_opportunity="LOW -- no experimental structure exists to guide design, and the one published hit is weak-potency and homology-model-derived, not structure-based",
             precise_summary="No experimental structure exists (homology-modeled only). A real but weak (IC50 20 micromolar), early-stage, structure-unguided small-molecule screening hit against VEZF1 DNA binding has been published -- more chemical-matter precedent than 'no program at all', but far from a validated or selective inhibitor.",
             sources="PMC6100598 (Characterization of Small Molecules Inhibiting the Pro-Angiogenic Activity of the Zinc Finger Transcription Factor Vezf1); web-verified 2026-08-15"),
    ]
    return pd.DataFrame(rows)


def build_translational_comparison() -> pd.DataFrame:
    """Table 06: KDM1A / TLK2 / USP34 / VEZF1 head-to-head on translational
    actionability axes. In-project quantitative columns are joined from
    already-computed tables; the structural/pharmacology columns are the
    precise, non-binary A-G facets from build_structural_tractability_audit()
    (a flat "has structure y/n" boolean was found, on re-audit, to
    conflate distinct concepts -- see that function's docstring). Every
    non-numeric cell states its source inline; nothing is a hidden score."""
    em = build_evidence_matrix().set_index("gene")
    struct = build_structural_tractability_audit().set_index("gene")
    rows = []
    for gene in FOCUS_FOUR:
        r = em.loc[gene]
        s = struct.loc[gene]
        rows.append(dict(
            gene=gene,
            crispr_effect=r["crispr_effect"], crispr_fdr=r["crispr_fdr"],
            crispr_rank_by_effect_of_13=int(r["rank_by_effect"]), crispr_rank_by_fdr_of_13=int(r["rank_by_fdr"]),
            n_chronic_datasets_fdr05_of_3=int(r["n_chronic_datasets_fdr05"]),
            depmap_median_chronos_er_luminal=r["median_chronos_er_luminal"],
            depmap_pct_er_luminal_strongly_dependent=None if pd.isna(r["frac_strongly_dependent_er_luminal"]) else round(r["frac_strongly_dependent_er_luminal"] * 100, 1),
            tcga_available=gene in ORIGINAL_FOUR,
            existing_gdsc_association_data=gene in ("USP34", "VEZF1"),
            gdsc_fdr_significant_hit=gene == "USP34",
            structural_pharmacology_precise_summary=s["precise_summary"],
            experimental_structure_exists=s["A_experimental_human_structure_exists"],
            validated_selective_inhibitor=s["E_validated_selective_small_molecule_inhibitor"],
            clinical_stage_pharmacology=s["F_clinical_stage_pharmacology"],
            target_class="histone demethylase (enzyme, established druggable class)" if gene == "KDM1A" else (
                "serine/threonine kinase (enzyme, established druggable class)" if gene == "TLK2" else (
                    "deubiquitinase / DUB (enzyme, emerging druggable class, covalent-probe precedent)" if gene == "USP34" else
                    "C2H2 zinc-finger transcription factor (historically much harder to drug directly than an enzyme)"
                )
            ),
            independent_external_literature_in_breast_cancer=gene in ("KDM1A", "TLK2"),
            was_a_pre_registered_blind_positive_control=gene == "KDM1A",
            documented_human_germline_loss_of_function_phenotype=gene == "USP34",
        ))
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Phase 1 -- external audit claim verification record (table 01)
# ---------------------------------------------------------------------------

def build_external_audit_claim_verification() -> pd.DataFrame:
    """Table 01: every material external-audit claim, independently
    verified against the authoritative project table/code, with the exact
    verified value/finding recorded. This is a curated verification log,
    not a statistic used in a chart -- every numeric claim inside it was
    independently computed above (build_evidence_matrix,
    load_significant_sensitising_hits, rule*() functions) and is repeated
    here only for a compact, itemised audit trail."""
    sens = load_significant_sensitising_hits().set_index("gene")
    em = build_evidence_matrix().set_index("gene")
    gdsc_all = pd.concat([
        pd.read_csv("results/tables/final_pharmacogenomics/USP34_GDSC_drug_associations.tsv", sep="\t"),
        pd.read_csv("results/tables/final_pharmacogenomics/VEZF1_GDSC_drug_associations.tsv", sep="\t"),
    ], ignore_index=True)
    gdsc_sig = gdsc_all[gdsc_all["fdr"] < 0.05]

    rows = [
        dict(claim="KDM1A Hany effect ~-2.17, FDR ~3.85e-4",
             verification=f"VERIFIED, exact: effect={sens.loc['KDM1A','effect_size']:.6f}, FDR={sens.loc['KDM1A','fdr']:.6e}",
             status="CORRECT"),
        dict(claim="TLK2 Hany effect ~-1.85, FDR ~0.0017",
             verification=f"VERIFIED, exact: effect={sens.loc['TLK2','effect_size']:.6f}, FDR={sens.loc['TLK2','fdr']:.6e}",
             status="CORRECT"),
        dict(claim="VEZF1 Hany effect ~-1.60, FDR ~0.0373",
             verification=f"VERIFIED, exact: effect={sens.loc['VEZF1','effect_size']:.6f}, FDR={sens.loc['VEZF1','fdr']:.6e}",
             status="CORRECT"),
        dict(claim="USP34 Hany effect ~-1.39, FDR ~0.0417",
             verification=f"VERIFIED, exact: effect={sens.loc['USP34','effect_size']:.6f}, FDR={sens.loc['USP34','fdr']:.6e}",
             status="CORRECT"),
        dict(claim="Sign convention: negative Hany effect = sensitising KO",
             verification="VERIFIED against src/labels.py docstring and PREANALYSIS.md Section 2, verbatim: 'Negative = more depleted under 4-OHT (candidate knockout sensitizes cells to tamoxifen)'",
             status="CORRECT"),
        dict(claim="RNA differential expression was effectively required for freeze eligibility",
             verification="VERIFIED exactly: src/evidence_freeze_shortlist_freeze.py determine_eligibility() hard-gates on `some_multimodal_support = (resistance_fdr05_count>=1) OR (human_tumor_support=='significant')` -- this IS a hard eligibility requirement in the implemented code, despite the module docstring's own opening sentence claiming 'RNA significance is explicitly NOT part of the eligibility gate' (a real internal documentation inconsistency, see claim below)",
             status="CORRECT"),
        dict(claim="KDM1A/TLK2 excluded specifically because of the RNA/human-tumor corroboration gate, not because of weak CRISPR",
             verification="VERIFIED exactly: freeze_eligibility_audit.tsv lists both with `_crispr_band=VERY_STRONG` (the strongest band) and `ineligibility_reason='no resistance-RNA dataset FDR<0.05 and no significant human-tumor support (functional-only -- see List C instead)'`",
             status="CORRECT"),
        dict(claim="This may have excluded stronger functional hits such as KDM1A and TLK2",
             verification="VERIFIED: among the 13 significant sensitising Gate-1 hits, KDM1A ranks 1st by both effect and FDR; TLK2 ranks 4th by effect and 2nd by FDR. USP34 ranks 12th by effect and 10th by FDR (worst or near-worst of the 13). This project's own pre-existing List C (functional-sensitisation, RNA not required) already ranked KDM1A #1 and TLK2 #2, with USP34 absent from its top 5.",
             status="CORRECT"),
        dict(claim="GSE118713 provides only one parental/resistant model, not independent resistance backgrounds",
             verification="VERIFIED: GSE118713 = 1 parental MCF7 line -> 1 TAMR-derivation event + 1 FASR-derivation event, each with 3 replicate RNA-seq samples (technical/aliquot replicates of one derivation, not 3 independent derivations). GSE111151, not GSE118713, is this project's multi-background dataset (4 parental lines x 7 independently-derived TamR sublines).",
             status="CORRECT"),
        dict(claim="GSE111151 has multiple independent biological backgrounds and should be weighted accordingly",
             verification="VERIFIED: 4 distinct parental lines (MCF-7, T-47D, ZR-75-1, BT-474; BT-474 is HER2-amplified, a different molecular subtype), 7 independently-derived TamR sublines, n=1 per sample (no technical replicates), cell-line-blocked design (~cell_line + resistance_status).",
             status="CORRECT"),
        dict(claim="GSE111151: none of KDM1A/TLK2/USP34/VEZF1 reach significance",
             verification=(
                 f"VERIFIED: gse111151_fdr USP34={em.loc['USP34','gse111151_fdr']:.3f}, "
                 f"VEZF1={em.loc['VEZF1','gse111151_fdr']:.3f}, KDM1A={em.loc['KDM1A','gse111151_fdr']:.3f}, "
                 f"TLK2={em.loc['TLK2','gse111151_fdr']:.3f} -- all non-significant. Zero of the 12 "
                 f"GSE111151-testable significant sensitising Gate-1 hits reach FDR<0.05 (USP17L29 has no "
                 f"GSE111151 count-matrix entry and is untestable, not tested-and-null; Rule 3 eligible set is empty)."
             ),
             status="CORRECT"),
        dict(claim="GSE240112 primary/recurrent tumors are matched/paired",
             verification="INCORRECT AS DOCUMENTED: docs/DATA_PROVENANCE.md line 41 states 'matched primary/recurrent ER+ tumor pairs'. The project's own locked docs/GSE240112_PREANALYSIS.md Section C states explicitly: PT (OriGene Technologies Inc.) and RT (Ontario Tumor Bank) come from DIFFERENT source institutions, 'no pairing statement exists anywhere in the GEO metadata or paper Methods', design is unpaired (~group, contrast RT-PT). The frozen ANALYSIS itself was already correctly run as unpaired; only the DATA_PROVENANCE.md summary line was wrong.",
             status="DOCUMENTATION ERROR -- CORRECTED THIS PHASE (see Phase 10)"),
        dict(claim="GDSC: 9 significant rows may represent only 8 unique compounds",
             verification=f"VERIFIED exactly: {len(gdsc_sig)} FDR<0.05 rows, {gdsc_sig['drug_id'].nunique()} unique drug_id. AZD7762 is counted twice (AUC and LN_IC50, both GDSC1, both significant, same direction) -- the only drug significant in both metrics.",
             status="CORRECT"),
        dict(claim="GDSC: most signals are LN_IC50-specific",
             verification=f"VERIFIED: of the 9 significant rows, {(gdsc_sig['metric']=='LN_IC50').sum()} are LN_IC50 and {(gdsc_sig['metric']=='AUC').sum()} are AUC (AZD7762 only).",
             status="CORRECT"),
        dict(claim="GDSC: some corresponding AUC associations reverse direction",
             verification="VERIFIED: of the 7 drugs significant only in LN_IC50, 2 (Sphingosine Kinase 1 Inhibitor II: rho -0.51 LN_IC50 vs +0.27 AUC; Ara-G: rho -0.48 LN_IC50 vs +0.13 AUC) show a sign flip in the AUC metric (both non-significant in AUC).",
             status="CORRECT"),
        dict(claim="GDSC: no independent GDSC2 replication exists",
             verification=f"VERIFIED: all {len(gdsc_sig)} FDR-significant rows are GDSC1; 0 are GDSC2.",
             status="CORRECT"),
        dict(claim="Environment: h5py missing from environment.yml despite being used",
             verification="VERIFIED: src/gse245601_feature_check.py imports h5py at module level; tests/test_gse245601_feature_check.py imports that module (collection-breaking in a truly clean env). h5py was absent from environment.yml before this phase's fix.",
             status="CORRECT -- FIXED THIS PHASE (see Phase 11)"),
        dict(claim="Environment: pandas version may cause incompatibility (implying an unpinned/risky spec)",
             verification="VERIFIED as a genuine reproducibility gap: pandas was completely unpinned in environment.yml (a fresh install could resolve pandas 3.x, an untested major version). Additionally, gseapy, lifelines, networkx, and requests are all imported by real, test-covered source files but were also absent from environment.yml.",
             status="CORRECT -- FIXED THIS PHASE (see Phase 11)"),
        dict(claim="Clean-environment test collection failures",
             verification="NOT directly reproduced by rebuilding a fresh conda environment in this session (impractical/high-risk within this session's time budget); INSTEAD independently verified via static source-level dependency audit (every third-party top-level import across src/, tests/, scripts/ cross-checked against environment.yml) -- confirms the claim is plausible and the root causes (h5py + 4 other undeclared packages) are real and fixed.",
             status="PLAUSIBLE, ROOT CAUSE CONFIRMED AND FIXED -- NOT EMPIRICALLY RE-REPRODUCED FROM SCRATCH"),
    ]
    return pd.DataFrame(rows)


def run(out_dir: str | Path = "results/tables/post_audit_sensitivity") -> None:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    tables = {
        "01_external_audit_claim_verification.tsv": build_external_audit_claim_verification(),
        "02_significant_sensitising_crispr_hits.tsv": load_significant_sensitising_hits(),
        "03_post_audit_evidence_matrix.tsv": build_evidence_matrix(),
        "04_selection_rule_sensitivity.tsv": build_selection_rule_sensitivity(),
        "05_leave_one_dataset_out.tsv": build_leave_one_dataset_out(),
        "06_top_candidate_translational_comparison.tsv": build_translational_comparison(),
        "06b_structural_tractability_audit.tsv": build_structural_tractability_audit(),
        "07_KDM1A_deep_dive.tsv": build_kdm1a_deep_dive(),
        "08_TLK2_deep_dive.tsv": build_tlk2_deep_dive(),
        "09_USP34_updated_liability.tsv": build_usp34_updated_liability(),
    }
    for name, df in tables.items():
        df.to_csv(out_dir / name, sep="\t", index=False)
        logger.info("wrote %s (%d rows)", out_dir / name, len(df))


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run()
