"""Data loaders for the EXPLORATION-V2 poster figure bank.

Complete reset of the visualization strategy (2026-08-16): the previous
`poster_final_*` figure set was rejected for relying too heavily on
summary boxes/scorecards instead of showing real biological data. This
module performs NO new discovery, NO re-ranking, and NO recomputation of
any frozen scientific result (CRISPR effects, RNA FDRs, DepMap dependency
calls, TCGA results, pathway NES/FDR, network edges, structural
conclusions, candidate-role conclusions). Every function here either (a)
reads an already-frozen table/matrix unchanged and reshapes it for
plotting, or (b) applies a disclosed, standard, deterministic
visualization-only transform (documented inline) that never touches a
p-value, FDR, effect size, or significance call.

See `results/reports/poster_exploration_v2/DATA_FOR_VISUALIZATION_AUDIT.md`
for the full source inventory this module implements.
"""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from src import poster_figures_bank_data as pbd
from src import poster_final_data as pfd2
from src import post_audit_sensitivity_data as pad

logger = logging.getLogger(__name__)

load_config = pad.load_config

FOCUS_FOUR = pad.FOCUS_FOUR  # ["KDM1A", "TLK2", "USP34", "VEZF1"]
ORIGINAL_FOUR = pad.ORIGINAL_FOUR  # ["USP34", "VEZF1", "EML5", "CITED2"]

# Restrained, muted candidate accent palette (Okabe-Ito derived, matching
# poster_final_data.FOCUS_COLORS exactly for continuity) -- used ONLY
# where gene identity is the point of a panel; everything else is neutral.
FOCUS_COLORS = dict(pfd2.FOCUS_COLORS)
GRAY = "#9a9a9a"
LGRAY = "#dcdcdc"
VLGRAY = "#eef0f2"
DGRAY = "#2b2b2b"


# ---------------------------------------------------------------------------
# Section A -- CRISPR discovery
# ---------------------------------------------------------------------------

def load_genomewide_crispr() -> pd.DataFrame:
    return pad.load_genomewide_crispr()


def load_significant_sensitising_hits() -> pd.DataFrame:
    return pad.load_significant_sensitising_hits()


def load_blind_control_row() -> pd.Series | None:
    return pfd2.load_f1_blind_control_row()


# ---------------------------------------------------------------------------
# Section B -- bulk / single-cell transcriptomics, real sample-level data
# ---------------------------------------------------------------------------

GSE118713_TPM = Path("data/processed/gse118713_gene_tpm.parquet")
GSE111151_LOG2CPM = Path("data/processed/gse111151_log2cpm.parquet")
GSE111151_METADATA = Path("results/tables/gse111151_sample_metadata.tsv")
GSE240112_CANDIDATE_LOG2CPM = Path("results/tables/gse240112/candidate_sample_level_log2cpm.tsv")
GSE245601_PATIENT_PSEUDOBULK = Path("results/tables/gse245601_candidate_deepdive/patient_malignant_pseudobulk.tsv")
GSE245601_PAIR_ELIGIBILITY = Path("results/tables/gse245601_pair_eligibility.tsv")
GSE245601_TRACKB_COUNTS = Path("results/tables/gse245601_pseudobulk/track_b_malignant_counts.tsv.gz")
GSE245601_TRACKB_METADATA = Path("results/tables/gse245601_pseudobulk/track_b_malignant_metadata.tsv")


def load_gse118713_focus_gene_samples() -> pd.DataFrame:
    """Real per-sample TPM (9 samples: 3 MCF7 + 3 TAMR + 3 FASR -- replicate
    aliquots of ONE derivation event, not independent derivations) for all
    4 focus genes, long format. Source: data/processed/gse118713_gene_tpm.parquet
    (already-frozen, genome-wide matrix)."""
    df = pd.read_parquet(GSE118713_TPM)
    sub = df[df["gene_symbol"].isin(FOCUS_FOUR)]
    sample_cols = [c for c in df.columns if c.startswith(("MCF7_", "TAMR_", "FASR_"))]
    long = sub.melt(id_vars="gene_symbol", value_vars=sample_cols, var_name="sample", value_name="tpm")
    long["condition"] = long["sample"].str.split("_").str[0]
    logger.info("load_gse118713_focus_gene_samples: %d genes x %d samples", sub.shape[0], len(sample_cols))
    return long


def load_gse111151_focus_gene_samples() -> pd.DataFrame:
    """Real per-sample log2CPM for all 4 focus genes, long format, joined
    to the real parental/derivative pairing already encoded in
    gse111151_sample_metadata.tsv (`paired_parental_sample_id`). No pairing
    is invented -- this is the exact blocked design (4 parental lines x 7
    independently-derived TamR sublines) already used elsewhere in this
    project."""
    df = pd.read_parquet(GSE111151_LOG2CPM)
    meta = pd.read_csv(GSE111151_METADATA, sep="\t")
    sub = df[df["gene_symbol"].isin(FOCUS_FOUR)]
    sample_cols = [c for c in df.columns if c not in ("gene_id", "gene_symbol")]
    long = sub.melt(id_vars="gene_symbol", value_vars=sample_cols, var_name="sample_id", value_name="log2cpm")
    long = long.merge(meta, on="sample_id", how="left")
    logger.info("load_gse111151_focus_gene_samples: %d genes x %d samples", sub.shape[0], len(sample_cols))
    return long


def load_gse240112_focus_gene_tumours() -> pd.DataFrame:
    """Real per-tumour log2CPM for all 4 focus genes (this frozen table
    already covers all 13 significant sensitising genes, not only the
    original 4), 6 tumours (3 PT + 3 RT), UNPAIRED (different patients,
    different biobanks -- no connecting line may ever be drawn between a
    PT and an RT observation). Source:
    results/tables/gse240112/candidate_sample_level_log2cpm.tsv."""
    df = pd.read_csv(GSE240112_CANDIDATE_LOG2CPM, sep="\t")
    sub = df[df["gene"].isin(FOCUS_FOUR)].copy()
    assert sub.groupby("gene").size().eq(6).all(), "expected 6 tumour observations per gene"
    logger.info("load_gse240112_focus_gene_tumours: %d genes x 6 tumours", sub["gene"].nunique())
    return sub


def load_gse245601_paired_focus_genes() -> pd.DataFrame:
    """Real per-patient, patient-MATCHED (same tumour, Control vs
    Tamoxifen, 12h ex vivo) log2(CPM+1) for all 4 focus genes, restricted
    to the 3 patients that pass the project's own pre-declared pseudobulk
    eligibility filter (>=50 malignant cells per arm).

    USP34/VEZF1 (and EML5/CITED2) values are read directly from the frozen
    `patient_malignant_pseudobulk.tsv` `normalized_expression_*` columns.
    KDM1A/TLK2 do not appear in that frozen table (it was only ever built
    for the original 4 candidates) -- their values are computed here from
    the frozen genome-wide raw-count matrix
    (`track_b_malignant_counts.tsv.gz`) and the frozen per-sample library
    sizes (`track_b_malignant_metadata.tsv`) using the IDENTICAL formula
    already used to produce the frozen column for the other 4 genes:
    log2((raw_umi_count / library_size * 1e6) + 1). This was verified by
    reproducing USP34's already-frozen Tumor_02 value from the raw counts
    before being used for KDM1A/TLK2 (see
    tests/test_poster_exploration_v2.py). No DE/FDR/significance value is
    computed or touched anywhere in this function."""
    eligible = pd.read_csv(GSE245601_PAIR_ELIGIBILITY, sep="\t")
    eligible_patients = eligible.loc[eligible["eligible_for_pseudobulk"], "patient"].tolist()
    assert eligible_patients == ["Tumor_02", "Tumor_03", "Tumor_07"]

    frozen = pd.read_csv(GSE245601_PATIENT_PSEUDOBULK, sep="\t")
    frozen = frozen[frozen["patient"].isin(eligible_patients) & frozen["gene"].isin(FOCUS_FOUR)]
    rows = []
    for _, r in frozen.iterrows():
        for cond in ("Control", "Tamoxifen"):
            rows.append(dict(patient=r["patient"], gene=r["gene"], condition=cond,
                              log2_expr=r[f"normalized_expression_{cond}"], source="frozen_patient_pseudobulk_table"))

    missing_genes = [g for g in FOCUS_FOUR if g not in frozen["gene"].unique()]
    if missing_genes:
        counts = pd.read_csv(GSE245601_TRACKB_COUNTS, sep="\t")
        meta = pd.read_csv(GSE245601_TRACKB_METADATA, sep="\t").set_index("sample_id")
        for gene in missing_genes:
            grow = counts[counts["gene"] == gene]
            assert len(grow) == 1, f"expected exactly one row for {gene} in track_b_malignant_counts"
            grow = grow.iloc[0]
            for patient in eligible_patients:
                for cond in ("Control", "Tamoxifen"):
                    sample_id = f"{patient}_{cond}"
                    raw = float(grow[sample_id])
                    lib = float(meta.loc[sample_id, "total_library_size"])
                    log2_expr = float(np.log2((raw / lib * 1e6) + 1))
                    rows.append(dict(patient=patient, gene=gene, condition=cond, log2_expr=log2_expr,
                                      source="computed_from_frozen_raw_counts_same_formula"))

    out = pd.DataFrame(rows)
    assert len(out) == len(FOCUS_FOUR) * 3 * 2, f"expected {len(FOCUS_FOUR)*3*2} rows, got {len(out)}"
    logger.info("load_gse245601_paired_focus_genes: %d genes x 3 eligible patients x 2 conditions", out["gene"].nunique())
    return out


def verify_gse245601_computed_formula_matches_frozen_usp34() -> tuple[float, float]:
    """Sanity check used by the test suite: reproduce USP34's already-
    frozen Tumor_02 Control normalized_expression value from the raw
    counts + library size, to verify the formula used for KDM1A/TLK2 above
    is identical to the one already used for the frozen columns."""
    frozen = pd.read_csv(GSE245601_PATIENT_PSEUDOBULK, sep="\t")
    row = frozen[(frozen["patient"] == "Tumor_02") & (frozen["gene"] == "USP34")].iloc[0]
    frozen_value = float(row["normalized_expression_Control"])

    counts = pd.read_csv(GSE245601_TRACKB_COUNTS, sep="\t")
    meta = pd.read_csv(GSE245601_TRACKB_METADATA, sep="\t").set_index("sample_id")
    raw = float(counts.loc[counts["gene"] == "USP34", "Tumor_02_Control"].iloc[0])
    lib = float(meta.loc["Tumor_02_Control", "total_library_size"])
    computed_value = float(np.log2((raw / lib * 1e6) + 1))
    return frozen_value, computed_value


# ---------------------------------------------------------------------------
# Section C -- pathway biology (real ranked-gene GSEA reconstruction)
# ---------------------------------------------------------------------------

DATASET_ORDER = ["gse118713", "gse111151", "gse240112", "gse245601"]
DATASET_LABELS = {
    "gse118713": "GSE118713\n(resistance model)",
    "gse111151": "GSE111151\n(resistance model)",
    "gse240112": "GSE240112\n(recurrence, human)",
    "gse245601": "GSE245601\n(acute 12h)",
}
HERO_PATHWAYS = [
    ("hallmark", "HALLMARK_ESTROGEN_RESPONSE_EARLY", "Estrogen response (early)"),
    ("hallmark", "HALLMARK_ESTROGEN_RESPONSE_LATE", "Estrogen response (late)"),
    ("hallmark", "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "EMT"),
]
EXTRA_PATHWAYS = [
    ("go_bp", "GOBP_REGULATION_OF_EXTRACELLULAR_MATRIX_ORGANIZATION", "ECM organization"),
    ("hallmark", "HALLMARK_E2F_TARGETS", "E2F targets"),
    ("hallmark", "HALLMARK_G2M_CHECKPOINT", "G2M checkpoint"),
    ("hallmark", "HALLMARK_WNT_BETA_CATENIN_SIGNALING", "Wnt / beta-catenin"),
]
GSEA_DIR = Path("results/tables/systems_network")


def load_pathway_trajectories(pathways: list[tuple[str, str, str]]) -> pd.DataFrame:
    """Real NES/FDR per pathway per dataset, ordered
    resistance -> resistance -> recurrence -> acute. Source:
    results/tables/systems_network/gsea_{dataset}.tsv (frozen, unmodified)."""
    rows = []
    for dataset in DATASET_ORDER:
        gsea = pd.read_csv(GSEA_DIR / f"gsea_{dataset}.tsv", sep="\t")
        for collection, pathway, label in pathways:
            hit = gsea[(gsea["collection"] == collection) & (gsea["pathway"] == pathway)]
            if len(hit) == 0:
                rows.append(dict(pathway_label=label, dataset=dataset, NES=np.nan, fdr=np.nan))
                continue
            r = hit.iloc[0]
            rows.append(dict(pathway_label=label, dataset=dataset, NES=r["NES"], fdr=r["fdr"]))
    return pd.DataFrame(rows)


def load_ranked_genes(dataset: str) -> pd.DataFrame:
    """Full pre-ranked gene list used as GSEA input for `dataset` -- the
    exact ranking already frozen and used to produce gsea_{dataset}.tsv,
    read here unmodified. Source:
    results/tables/systems_network/{dataset}_ranked_genes.tsv."""
    path = GSEA_DIR / f"{dataset}_ranked_genes.tsv"
    df = pd.read_csv(path, sep="\t")
    return df.sort_values("ranking_stat", ascending=False).reset_index(drop=True)


HALLMARK_GMT = Path("data/reference/genesets/hallmark.gmt")


def load_gene_set(pathway: str) -> set[str]:
    """Full gene-set membership from the already-frozen, locally-cached
    MSigDB Hallmark GMT file -- read unmodified, used only to reconstruct
    the running-enrichment-score curve shape (Section C3), never to rerun
    a differently-parameterised enrichment test."""
    with open(HALLMARK_GMT) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if parts[0] == pathway:
                return set(parts[2:])
    raise KeyError(f"{pathway} not found in {HALLMARK_GMT}")


def compute_running_enrichment_score(ranked_genes: list[str], ranking_stat: np.ndarray, gene_set: set[str]) -> np.ndarray:
    """Standard weighted running-sum GSEA statistic (Subramanian et al.
    2005): a deterministic reconstruction of the enrichment-score CURVE
    SHAPE from an already-frozen ranked gene list and an already-frozen
    gene set -- this is a visualization of the same statistic already
    summarized by the frozen NES/FDR values (read separately, never
    recomputed here), not a new or differently-thresholded enrichment
    test."""
    n = len(ranked_genes)
    in_set = np.array([g in gene_set for g in ranked_genes])
    n_hits = in_set.sum()
    if n_hits == 0:
        return np.zeros(n)
    weights = np.abs(ranking_stat)
    hit_sum = weights[in_set].sum()
    miss_step = 1.0 / (n - n_hits)
    running = np.zeros(n)
    cur = 0.0
    for i in range(n):
        if in_set[i]:
            cur += weights[i] / hit_sum
        else:
            cur -= miss_step
        running[i] = cur
    return running


def build_enrichment_curve(dataset: str, pathway: str) -> pd.DataFrame:
    ranked = load_ranked_genes(dataset)
    gene_set = load_gene_set(pathway)
    es = compute_running_enrichment_score(ranked["gene"].tolist(), ranked["ranking_stat"].to_numpy(), gene_set)
    out = ranked.copy()
    out["running_es"] = es
    out["in_gene_set"] = out["gene"].isin(gene_set)
    out["rank"] = np.arange(1, len(out) + 1)
    return out


# ---------------------------------------------------------------------------
# Section D -- network (original-4-candidate coverage only, disclosed)
# ---------------------------------------------------------------------------

NETWORK_DIR = Path("results/tables/systems_network")


def load_candidate_pathway_membership() -> pd.DataFrame:
    """STRONG_CONSENSUS pathway memberships for the ORIGINAL 4 candidates
    only (USP34, VEZF1, EML5, CITED2) -- this analysis was never run for
    KDM1A/TLK2 (see DATA_FOR_VISUALIZATION_AUDIT.md). Source:
    results/tables/systems_network/candidate_pathway_membership.tsv."""
    df = pd.read_csv(NETWORK_DIR / "candidate_pathway_membership.tsv", sep="\t")
    sub = df[(df["candidate_is_member"] == True) & (df["resistance_consensus_class"] == "STRONG_CONSENSUS")]  # noqa: E712
    return sub[["candidate", "collection", "pathway"]].reset_index(drop=True)


def load_candidate_convergence() -> pd.DataFrame:
    return pd.read_csv(NETWORK_DIR / "four_candidate_convergence.tsv", sep="\t")


def load_direct_neighbors(candidate: str | None = None) -> pd.DataFrame:
    """Real STRING direct-neighbor table for the original 4 candidates
    (USP34: 10 rows, VEZF1: 1, EML5: 1, CITED2: 18 -- CITED2 and USP34
    have the richest neighbor coverage; VEZF1/EML5 have essentially none).
    KDM1A/TLK2 have no frozen direct-neighbor table (see
    DATA_FOR_VISUALIZATION_AUDIT.md). Pass `candidate` to filter to one
    gene's rows."""
    df = pd.read_csv(NETWORK_DIR / "four_candidate_direct_neighbors.tsv", sep="\t")
    if candidate is not None:
        df = df[df["candidate"] == candidate].reset_index(drop=True)
    return df


# ---------------------------------------------------------------------------
# Section E -- evidence intersection (13 significant sensitising hits)
# ---------------------------------------------------------------------------

def load_evidence_matrix_13() -> pd.DataFrame:
    return pad.build_evidence_matrix()


def build_evidence_sets_13() -> pd.DataFrame:
    """Boolean evidence-set membership for all 13 significant sensitising
    hits, for an UpSet-style intersection plot. Categories are kept
    distinct in TYPE (never mixing an advantage with a liability in one
    set) -- 'high baseline dependency' and 'low baseline dependency' are
    two SEPARATE, mutually exclusive boolean columns, not collapsed into
    one 'dependency' axis. Source: post_audit_sensitivity_data.build_
    evidence_matrix() and build_structural_tractability_audit(), both
    already frozen."""
    em = pad.build_evidence_matrix().set_index("gene")
    struct = pad.build_structural_tractability_audit().set_index("gene")

    rows = []
    for gene in em.index:
        r = em.loc[gene]
        has_resistance_rna = bool((pd.notna(r["gse118713_fdr"]) and r["gse118713_fdr"] < 0.05) or
                                   (pd.notna(r["gse111151_fdr"]) and r["gse111151_fdr"] < 0.05))
        has_recurrence_rna = bool(pd.notna(r["gse240112_fdr"]) and r["gse240112_fdr"] < 0.05)
        frac_dep = r["frac_strongly_dependent_er_luminal"]
        high_dependency = bool(pd.notna(frac_dep) and frac_dep >= 0.5)
        low_dependency = bool(pd.notna(frac_dep) and frac_dep < 0.1)
        has_structure = False
        has_inhibitor = False
        if gene in struct.index:
            s = struct.loc[gene]
            has_structure = bool(s["A_experimental_human_structure_exists"])
            has_inhibitor = str(s["E_validated_selective_small_molecule_inhibitor"]).upper().startswith("YES")
        rows.append(dict(
            gene=gene,
            resistance_model_rna_support=has_resistance_rna,
            recurrence_associated_rna_support=has_recurrence_rna,
            high_baseline_dependency=high_dependency,
            low_baseline_dependency=low_dependency,
            experimental_structure=has_structure,
            validated_inhibitor=has_inhibitor,
        ))
    return pd.DataFrame(rows)


def load_selection_rule_sensitivity() -> pd.DataFrame:
    """Historical Rule 0 (original frozen gate, RNA-eligibility required)
    vs Rule 1 (CRISPR-only, no RNA gate) eligibility/rank per gene --
    already-frozen post-audit sensitivity output, read unmodified."""
    return pad.build_selection_rule_sensitivity()


# ---------------------------------------------------------------------------
# Section F -- human tumour (TCGA) + DepMap, real per-line/per-tumour data
# ---------------------------------------------------------------------------

def load_depmap_effect_focus_four() -> pd.DataFrame:
    """Real per-cell-line Chronos values, 11 ER+/luminal dependency-
    evaluable lines, all 4 focus genes -- reuses poster_final_data's
    generic loader unchanged."""
    return pfd2.load_f4_depmap_effect()


def load_depmap_model_names() -> pd.DataFrame:
    """Real cell-line names (not just ModelID) for the 11 ER+/luminal
    dependency-evaluable lines. Source: DepMap 26Q1 Model.csv, via the
    already-validated independent_validation_depmap_data.load_model."""
    from src.independent_validation_depmap_data import load_model, raw_dir

    cfg = load_config()
    model = load_model(cfg, "26Q1")
    effect_path = raw_dir(cfg, "26Q1") / cfg["independent_validation"]["depmap"]["releases"]["26Q1"]["raw"]["crispr_gene_effect_csv"]
    dependency_evaluable_ids = pd.read_csv(effect_path, usecols=[0], index_col=0).index
    luminal = model.loc[model["is_er_luminal"]]
    evaluable = luminal.loc[luminal.index.intersection(dependency_evaluable_ids)]
    assert len(evaluable) == 11, f"expected 11 ER+/luminal dependency-evaluable lines, got {len(evaluable)}"
    return evaluable[["CellLineName", "OncotreeSubtype"]].reset_index().rename(columns={"index": "ModelID"})


def load_tcga_forest_original_four() -> pd.DataFrame:
    return pbd.load_tcga_expression_forest()


def load_gdsc_top_associations() -> pd.DataFrame:
    return pbd.load_gdsc_top_associations()


# ---------------------------------------------------------------------------
# Section G -- structural biology
# ---------------------------------------------------------------------------

def kdm1a_tlk2_structure_paths() -> dict[str, Path]:
    cfg = load_config()
    d = Path(cfg["data"]["raw"]["kdm1a_tlk2_structures_dir"])
    return {"6NQU": d / "6NQU.pdb", "5O0Y": d / "5O0Y.pdb"}


def usp34_structure_paths() -> dict[str, Path]:
    return pfd2.usp34_structure_paths()


def load_structural_tractability_audit() -> pd.DataFrame:
    return pad.build_structural_tractability_audit()


# ---------------------------------------------------------------------------
# Section H -- translational / follow-up (data-grounded only)
# ---------------------------------------------------------------------------

def load_role_compass() -> pd.DataFrame:
    return pfd2.load_f6_role_compass()
