"""Data loaders for the FINAL poster figure set (results/figures/poster_final/).

This module performs NO new discovery, NO recomputation of any frozen
statistic, and NO re-ranking. It only reshapes already-frozen project
tables -- the same post-audit-sensitivity and figure-bank source tables
already read (unmodified) by ``src/post_audit_sensitivity_data.py`` and
``src/poster_figures_bank_data.py`` -- into the specific shapes the six
final poster figures need. Every function below documents its exact
upstream source; nothing here is hand-typed.

Frozen sources reused, never altered:
  - data/processed/labels.parquet (genome-wide Hany CRISPR fit, 19,103 genes)
  - results/tables/post_audit_sensitivity/02_significant_sensitising_crispr_hits.tsv
    (equivalently, ``post_audit_sensitivity_data.load_significant_sensitising_hits()``)
  - results/tables/post_audit_sensitivity/03_post_audit_evidence_matrix.tsv
    (equivalently, ``post_audit_sensitivity_data.build_evidence_matrix()``)
  - results/tables/post_audit_sensitivity/06b_structural_tractability_audit.tsv
  - results/tables/systems_network/gsea_{dataset}.tsv (per-dataset GSEA)
  - results/tables/independent_validation/TCGA_candidate_expression.tsv
  - DepMap 26Q1 CRISPRGeneEffect.csv (raw matrix, generic symbol-column match)
  - the two USP34 PDB structures (7W3R, 7W3U), already frozen structural
    evidence, re-rendered here only for a different visual composition
"""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pandas as pd

from src import poster_figures_bank_data as pbd
from src import poster_figures_data as pfd
from src import post_audit_sensitivity_data as pad

logger = logging.getLogger(__name__)

load_config = pad.load_config

FOCUS_FOUR = pad.FOCUS_FOUR  # ["KDM1A", "TLK2", "USP34", "VEZF1"] -- reused unchanged
BLIND_CONTROL_RECOVERED = "KDM1A"
BLIND_CONTROL_NOT_RECOVERED = "RCOR1"  # PREANALYSIS.md Section 5: the other pre-registered blind control

# Fixed candidate identity colors -- Okabe-Ito colorblind-safe set, validated
# with the dataviz skill's validate_palette.js (all pairwise CVD checks
# PASS; the WARN on contrast-vs-white-surface for VEZF1/TLK2 is satisfied
# throughout by always pairing these colors with a direct label or legend,
# never a bare unlabeled swatch). USP34/VEZF1 reuse the exact hex values
# already frozen in poster_figures_visualization.GENE_COLORS for visual
# continuity with the earlier figure bank; KDM1A/TLK2 are new, chosen from
# the same 8-color Okabe-Ito family.
FOCUS_COLORS = {
    "USP34": "#0072B2",
    "VEZF1": "#E69F00",
    "KDM1A": "#D55E00",
    "TLK2": "#56B4E9",
}
GRAY = "#9a9a9a"
LGRAY = "#d8d8d8"
DGRAY = "#333333"


# ---------------------------------------------------------------------------
# Figure 1 -- genome-wide CRISPR discovery
# ---------------------------------------------------------------------------

def load_f1_genomewide() -> pd.DataFrame:
    """All 19,103 fitted genes, with a Gate-1 (FDR<0.1) flag and a
    focus-gene flag. Source: ``post_audit_sensitivity_data.load_genomewide_crispr()``
    (data/processed/labels.parquet), unmodified."""
    df = pad.load_genomewide_crispr().copy()
    df["gate1_significant"] = df["fdr"] < pad.CRISPR_GATE1_FDR
    df["is_focus_gene"] = df["gene"].isin(FOCUS_FOUR)
    return df


def load_f1_focus_gene_ranks() -> pd.DataFrame:
    """Focus-gene rows from the frozen 13-gene significant-sensitising
    table (Table 02), with rank_by_effect / rank_by_fdr out of 13 --
    the "USP34 is not the top hit" evidence, read directly, never
    hand-typed. Source: ``post_audit_sensitivity_data.load_significant_sensitising_hits()``."""
    sens = pad.load_significant_sensitising_hits()
    n = len(sens)
    out = sens[sens["gene"].isin(FOCUS_FOUR)][["gene", "effect_size", "fdr", "rank_by_effect", "rank_by_fdr"]].copy()
    out["n_sensitising_hits"] = n
    return out.reset_index(drop=True)


def load_f1_blind_control_row() -> pd.Series | None:
    """RCOR1's genome-wide row (the pre-registered blind control that was
    NOT recovered at the Gate-1 FDR<0.1 threshold, unlike KDM1A) -- an
    honest, real reference point, not a fabricated control. Returns None
    if RCOR1 is absent from the fitted table (should not happen)."""
    df = pad.load_genomewide_crispr()
    row = df[df["gene"] == BLIND_CONTROL_NOT_RECOVERED]
    if len(row) == 0:
        return None
    return row.iloc[0]


# ---------------------------------------------------------------------------
# Figure 2 -- transcriptomic / pathway systems view
# ---------------------------------------------------------------------------

DATASET_CATEGORY = {
    "gse118713": "resistance model (cell line)",
    "gse111151": "resistance model (cell line)",
    "gse240112": "recurrence-associated (human tumour, unpaired)",
    "gse245601": "acute 12h (not resistance)",
}


def load_f2_pathway_matrix() -> pd.DataFrame:
    """Fixed, pre-declared pathway list x 4 datasets, tidy (pathway_label,
    dataset, NES, FDR, dataset_category). Source:
    ``poster_figures_bank_data.load_pathway_landscape()`` (already-frozen
    per-dataset GSEA tables), with a dataset_category column attached here
    to keep resistance-model / recurrence-associated / acute-context
    datasets visually distinguishable -- GSE240112 is never labeled
    "chronic resistance" or grouped with the resistance-model datasets."""
    df = pbd.load_pathway_landscape()
    df = df.copy()
    df["dataset_category"] = df["dataset_key"].map(DATASET_CATEGORY)
    assert df["dataset_category"].notna().all()
    return df


# ---------------------------------------------------------------------------
# Figure 3 -- candidate evidence divergence / integration (4 focus genes)
# ---------------------------------------------------------------------------

def load_f3_evidence_matrix() -> pd.DataFrame:
    """The 4 focus genes' rows from Table 03 (post-audit evidence matrix),
    unmodified. Source: ``post_audit_sensitivity_data.build_evidence_matrix()``."""
    em = pad.build_evidence_matrix()
    out = em[em["gene"].isin(FOCUS_FOUR)].set_index("gene").loc[FOCUS_FOUR].reset_index()
    return out


def _facet_token(value) -> str:
    """First token (YES/NO/PARTIAL/True/False) of a Table 06b facet cell,
    parsed deterministically -- these cells are already fixed, cited
    strings; this only extracts the leading call for a compact plot
    encoding, it does not reinterpret or re-judge the evidence."""
    if isinstance(value, bool):
        return "YES" if value else "NO"
    s = str(value).strip().upper()
    if s.startswith("YES"):
        return "YES"
    if s.startswith("NO"):
        return "NO"
    if s.startswith("PARTIAL") or s.startswith("PLAUSIBLE"):
        return "PARTIAL"
    return "NO"


def load_f3_structural_facets() -> pd.DataFrame:
    """4 focus genes x 4 structural/pharmacology facets (structure exists,
    ligand/probe-bound structure, validated selective inhibitor,
    clinical-stage pharmacology), each encoded YES/PARTIAL/NO from Table
    06b's own cited text -- never re-judged. Source:
    ``post_audit_sensitivity_data.build_structural_tractability_audit()``."""
    struct = pad.build_structural_tractability_audit()
    struct = struct[struct["gene"].isin(FOCUS_FOUR)].set_index("gene").loc[FOCUS_FOUR]
    rows = []
    for gene, r in struct.iterrows():
        rows.append(dict(
            gene=gene,
            structure_exists=_facet_token(r["A_experimental_human_structure_exists"]),
            ligand_or_probe_bound=_facet_token(r["C_ligand_or_probe_bound_structure"]),
            validated_inhibitor=_facet_token(r["E_validated_selective_small_molecule_inhibitor"]),
            clinical_stage=_facet_token(r["F_clinical_stage_pharmacology"]),
        ))
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Figure 4 -- human (TCGA) / DepMap orthogonal validation
# ---------------------------------------------------------------------------

def load_f4_tcga_forest() -> pd.DataFrame:
    """TCGA-BRCA paired tumor-vs-normal log2FC + 95% CI, for the 2 focus
    genes TCGA evidence actually exists for (USP34, VEZF1) -- KDM1A/TLK2
    are NOT in this table because TCGA follow-up in this project was only
    ever run on the original 4-candidate set, a real evidence-depth gap,
    not a null result. Source:
    ``poster_figures_bank_data.load_tcga_expression_forest()``, subset to
    the 2 focus genes with TCGA coverage."""
    forest = pbd.load_tcga_expression_forest()
    sub = forest[(forest["candidate"].isin(FOCUS_FOUR)) & (forest["comparison"] == "tumor_vs_normal_PAIRED")].copy()
    return sub.reset_index(drop=True)


def load_f4_depmap_effect() -> pd.DataFrame:
    """Real per-cell-line Chronos gene-effect values (DepMap 26Q1), the 11
    ER+/luminal screened lines, for all 4 focus genes (including KDM1A/
    TLK2, which the original candidate-only loader in
    ``poster_figures_data`` cannot reach because it is keyed to a
    hardcoded 4-candidate Entrez map) -- generic symbol-column matching,
    the same approach already documented and used in
    ``poster_figures_bank_data._gene_col_generic`` /
    ``post_audit_sensitivity_data.load_depmap_summary_for_genes``, applied
    here to return the full per-line array (needed for a distribution
    plot) rather than only a summary statistic."""
    from src.independent_validation_depmap_data import load_model, raw_dir

    cfg = load_config()
    release = "26Q1"
    rel_cfg = cfg["independent_validation"]["depmap"]["releases"][release]
    effect_path = raw_dir(cfg, release) / rel_cfg["raw"]["crispr_gene_effect_csv"]

    model = load_model(cfg, release)
    luminal_ids = model.index[model["is_er_luminal"]]

    effect_df = pd.read_csv(effect_path, index_col=0)
    col_map = {}
    for gene in FOCUS_FOUR:
        col = pbd._gene_col_generic(effect_df.columns, gene)
        if col is not None:
            col_map[col] = gene
    out = effect_df[list(col_map)].rename(columns=col_map)
    out = out.loc[out.index.intersection(luminal_ids)]
    assert len(out) == 11, f"expected 11 ER+/luminal screened lines, got {len(out)}"
    id_col = out.index.name or "index"
    long = out.reset_index().melt(id_vars=id_col, var_name="gene", value_name="chronos_effect")
    long = long.rename(columns={id_col: "cell_line"})[["cell_line", "gene", "chronos_effect"]]
    logger.info("load_f4_depmap_effect: %d lines x %d genes", len(out), len(col_map))
    return long


# ---------------------------------------------------------------------------
# Figure 5 -- USP34 structure and tractability
# ---------------------------------------------------------------------------

def load_f5_usp34_structural_row() -> pd.Series:
    """USP34's row from Table 06b (structural tractability audit),
    unmodified -- the caption-strip wording for Figure 5 is built directly
    from this row's own cited fields, never re-typed from memory."""
    struct = pad.build_structural_tractability_audit()
    return struct[struct["gene"] == "USP34"].iloc[0]


def usp34_structure_paths() -> dict[str, Path]:
    """Reuses ``poster_figures_data.usp34_structure_paths`` unchanged --
    the same two frozen PDB files (7W3R apo, 7W3U covalent-probe-bound)."""
    return pfd.usp34_structure_paths()


# ---------------------------------------------------------------------------
# Figure 6 -- final candidate logic / experimental strategy
# ---------------------------------------------------------------------------

def load_f6_role_compass() -> pd.DataFrame:
    """All 13 significant sensitising genes (CRISPR strength = -log10 FDR)
    joined to DepMap ER+/luminal dependency fraction -- the real,
    data-grounded scatter behind Figure 6's "role compass" panel, showing
    sensitisation strength is not the same axis as baseline dependency.
    Source: ``post_audit_sensitivity_data.build_evidence_matrix()``
    (already includes both quantities for the full 13-gene universe)."""
    em = pad.build_evidence_matrix()
    out = em[["gene", "crispr_effect", "crispr_fdr", "rank_by_effect", "frac_strongly_dependent_er_luminal"]].copy()
    out["neg_log10_fdr"] = -np.log10(out["crispr_fdr"])
    out["is_focus_gene"] = out["gene"].isin(FOCUS_FOUR)
    return out
