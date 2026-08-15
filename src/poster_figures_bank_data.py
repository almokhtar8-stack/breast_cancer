"""Data loaders for the poster FIGURE BANK (candidate-figure phase).

Companion module to `poster_figures_data.py`. This module performs no new
discovery, no re-ranking, and no recomputation of any already-frozen
statistic. Where a small deterministic transformation of frozen data is
required purely to make a plot (e.g. selecting a fixed list of biologically
named pathways out of an already-computed GSEA table, or generically
matching DepMap gene-effect columns by symbol instead of a hardcoded
Entrez-ID subset), that transformation is implemented here, explicitly, and
is documented in `results/reports/poster/SCIENTIFIC_FIGURE_AUDIT.md`.

Every plotted number is read directly from an already-frozen project source
table:

- results/tables/nebula_plot_inputs/fig2_gse118713_pca_input.tsv (frozen PCA
  coordinates, an earlier, visually-approved phase's output -- not recomputed)
- results/tables/nebula_plot_inputs/fig3_tamr_vs_mcf7_volcano_input.tsv
  (frozen genome-wide TAMR-vs-MCF7 limma DE, 14,837 genes)
- results/tables/nebula_plot_inputs/fig1_crispr_hit_landscape_input.tsv
  (frozen 28-gene Gate-1 CRISPR hit table)
- results/tables/gse240112_pseudobulk/tumor_cell_genomewide_de.tsv.gz
  (genome-wide edgeR DE, primary vs recurrent tumor pseudobulk)
- results/tables/systems_network/gsea_{dataset}.tsv (per-dataset GSEA,
  hallmark + go_bp collections)
- results/tables/systems_network/candidate_pathway_membership.tsv
- results/tables/independent_validation/TCGA_candidate_expression.tsv,
  TCGA_candidate_clinical.tsv
- results/tables/final_pharmacogenomics/GDSC_top_associations.tsv
- results/tables/lead_target_deep_dive/USP34_VEZF1_tissue_liability.tsv,
  USP34_VEZF1_full_tissue_expression.tsv
- DepMap 26Q1 CRISPRGeneEffect.csv / CRISPRGeneDependency.csv (raw matrices,
  reused via a generic symbol-column match -- see
  `_gene_col_generic()` -- rather than the candidate-only
  `independent_validation_depmap_data._gene_col()`, so the same file can
  answer the 28-gene Gate-1 question this phase asks, still with no new
  network calls and no recomputation of the underlying Chronos scores)
"""

from __future__ import annotations

import gzip
import logging
from pathlib import Path

import pandas as pd

from src.poster_figures_data import CANDIDATES, load_config  # noqa: F401 (re-exported)

logger = logging.getLogger(__name__)

NEBULA_INPUTS = Path("results/tables/nebula_plot_inputs")
GSE118713_PCA = NEBULA_INPUTS / "fig2_gse118713_pca_input.tsv"
GSE118713_VOLCANO = NEBULA_INPUTS / "fig3_tamr_vs_mcf7_volcano_input.tsv"
CRISPR_GATE1 = NEBULA_INPUTS / "fig1_crispr_hit_landscape_input.tsv"

GSE240112_GENOMEWIDE_DE = Path("results/tables/gse240112_pseudobulk/tumor_cell_genomewide_de.tsv.gz")

GSEA_DIR = Path("results/tables/systems_network")
CANDIDATE_PATHWAY_MEMBERSHIP = GSEA_DIR / "candidate_pathway_membership.tsv"

TCGA_EXPRESSION = Path("results/tables/independent_validation/TCGA_candidate_expression.tsv")
TCGA_CLINICAL = Path("results/tables/independent_validation/TCGA_candidate_clinical.tsv")

GDSC_TOP_ASSOCIATIONS = Path("results/tables/final_pharmacogenomics/GDSC_top_associations.tsv")

TISSUE_LIABILITY = Path("results/tables/lead_target_deep_dive/USP34_VEZF1_tissue_liability.tsv")
TISSUE_EXPRESSION = Path("results/tables/lead_target_deep_dive/USP34_VEZF1_full_tissue_expression.tsv")
NORMAL_TISSUE_CONTEXT = Path("results/tables/druggability_safety/candidate_normal_tissue_context.tsv")


# ---------------------------------------------------------------------------
# CRISPR (genome-wide + Gate-1 subset)
# ---------------------------------------------------------------------------

def load_crispr_gate1() -> pd.DataFrame:
    """The frozen 28-gene Gate-1 CRISPR hit table (FDR<0.1), already
    extracted in an earlier phase -- read here unchanged."""
    df = pd.read_csv(CRISPR_GATE1, sep="\t")
    assert len(df) == 28, f"expected 28 Gate-1 hits, got {len(df)}"
    return df


# ---------------------------------------------------------------------------
# GSE118713 (PCA + genome-wide volcano, frozen nebula-phase inputs)
# ---------------------------------------------------------------------------

def load_gse118713_pca() -> pd.DataFrame:
    df = pd.read_csv(GSE118713_PCA, sep="\t")
    assert len(df) == 9, f"expected 9 samples, got {len(df)}"
    return df


def load_gse118713_volcano() -> pd.DataFrame:
    df = pd.read_csv(GSE118713_VOLCANO, sep="\t")
    logger.info("load_gse118713_volcano: %d genes", len(df))
    return df


# ---------------------------------------------------------------------------
# GSE240112 (genome-wide primary-vs-recurrent volcano)
# ---------------------------------------------------------------------------

def load_gse240112_volcano() -> pd.DataFrame:
    df = pd.read_csv(GSE240112_GENOMEWIDE_DE, sep="\t", compression="gzip")
    df["is_candidate"] = df["gene"].isin(CANDIDATES)
    logger.info("load_gse240112_volcano: %d genes", len(df))
    return df


# ---------------------------------------------------------------------------
# Pathway landscape (fixed, biologically-chosen pathway list)
# ---------------------------------------------------------------------------

PATHWAY_LANDSCAPE = [
    ("hallmark", "HALLMARK_ESTROGEN_RESPONSE_EARLY", "Estrogen response (early)"),
    ("hallmark", "HALLMARK_ESTROGEN_RESPONSE_LATE", "Estrogen response (late)"),
    ("hallmark", "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "EMT"),
    ("go_bp", "GOBP_REGULATION_OF_EXTRACELLULAR_MATRIX_ORGANIZATION", "ECM organization"),
    ("hallmark", "HALLMARK_E2F_TARGETS", "E2F targets"),
    ("hallmark", "HALLMARK_G2M_CHECKPOINT", "G2M checkpoint"),
    ("hallmark", "HALLMARK_UV_RESPONSE_DN", "UV response (down)"),
    ("hallmark", "HALLMARK_WNT_BETA_CATENIN_SIGNALING", "Wnt / beta-catenin"),
    ("go_bp", "GOBP_BLOOD_VESSEL_MORPHOGENESIS", "Blood vessel morphogenesis"),
    ("hallmark", "HALLMARK_HEME_METABOLISM", "Heme metabolism"),
]

# Resistance (transcriptomic, comparable NES/FDR) datasets only. GSE245601 is
# an acute 12h response, not resistance -- included but must stay visually
# separate from the 3 resistance datasets. CRISPR pathway context uses a
# different metric entirely (gene-level dependency ranking, not an
# expression fold-change) and is loaded separately, never on the same color
# scale as the RNA panels.
PATHWAY_RNA_DATASETS = {
    "gse118713": "GSE118713",
    "gse240112": "GSE240112",
    "gse111151": "GSE111151",
    "gse245601": "GSE245601 (acute 12h)",
}


def load_pathway_landscape() -> pd.DataFrame:
    """Tidy (pathway_label, dataset, NES, FDR) table for the fixed pathway
    list above, read directly from the per-dataset GSEA tables -- no
    re-running of GSEA, no new pathway list beyond the biologically-motivated
    selection documented in SCIENTIFIC_FIGURE_AUDIT.md section 7."""
    rows = []
    for dataset_key, dataset_label in PATHWAY_RNA_DATASETS.items():
        path = GSEA_DIR / f"gsea_{dataset_key}.tsv"
        gsea = pd.read_csv(path, sep="\t")
        for collection, pathway, label in PATHWAY_LANDSCAPE:
            hit = gsea[(gsea["collection"] == collection) & (gsea["pathway"] == pathway)]
            assert len(hit) <= 1, f"{dataset_key}: expected at most one row for {pathway}, got {len(hit)}"
            if len(hit) == 0:
                # Gene set did not meet this dataset's testable-size filter --
                # a real absence, not fabricated as zero/non-significant.
                rows.append(dict(
                    pathway_label=label, collection=collection, pathway=pathway,
                    dataset_key=dataset_key, dataset_label=dataset_label,
                    NES=float("nan"), fdr=float("nan"),
                ))
                continue
            r = hit.iloc[0]
            rows.append(dict(
                pathway_label=label, collection=collection, pathway=pathway,
                dataset_key=dataset_key, dataset_label=dataset_label,
                NES=r["NES"], fdr=r["fdr"],
            ))
    out = pd.DataFrame(rows)
    logger.info("load_pathway_landscape: %d pathways x %d datasets = %d rows", len(PATHWAY_LANDSCAPE), len(PATHWAY_RNA_DATASETS), len(out))
    return out


def load_pathway_landscape_crispr() -> pd.DataFrame:
    """CRISPR pathway-context NES/FDR for the same fixed pathway list --
    loaded and returned SEPARATELY from load_pathway_landscape() because the
    CRISPR metric (gene-level dependency-ranking enrichment) is not the same
    quantity as an expression-fold-change NES and must never be merged onto
    the same color scale."""
    gsea = pd.read_csv(GSEA_DIR / "gsea_crispr.tsv", sep="\t")
    rows = []
    for collection, pathway, label in PATHWAY_LANDSCAPE:
        hit = gsea[(gsea["collection"] == collection) & (gsea["pathway"] == pathway)]
        assert len(hit) <= 1
        if len(hit) == 0:
            rows.append(dict(pathway_label=label, NES=float("nan"), fdr=float("nan")))
            continue
        r = hit.iloc[0]
        rows.append(dict(pathway_label=label, NES=r["NES"], fdr=r["fdr"]))
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Candidate -> pathway-program mechanism map (curated, not a force-directed graph)
# ---------------------------------------------------------------------------

CANDIDATE_CONVERGENCE = GSEA_DIR / "four_candidate_convergence.tsv"


def load_candidate_convergence() -> pd.DataFrame:
    """Pairwise candidate-candidate shared-pathway convergence counts (all 6
    pairs) -- the real connectivity data behind the mechanism-map figure's
    dashed edges, read directly rather than hand-typed."""
    df = pd.read_csv(CANDIDATE_CONVERGENCE, sep="\t")
    assert len(df) == 6
    return df


def load_candidate_strong_consensus_pathways() -> pd.DataFrame:
    """Per-candidate STRONG_CONSENSUS pathway memberships -- the same table
    used to verify the mechanism-map connectivity in the audit. EML5 is
    expected to return zero rows (real finding, not a bug)."""
    df = pd.read_csv(CANDIDATE_PATHWAY_MEMBERSHIP, sep="\t")
    sub = df[(df["candidate_is_member"] == True) & (df["resistance_consensus_class"] == "STRONG_CONSENSUS")]  # noqa: E712
    logger.info("load_candidate_strong_consensus_pathways: %d rows", len(sub))
    return sub[["candidate", "collection", "pathway"]].reset_index(drop=True)


# ---------------------------------------------------------------------------
# TCGA-BRCA human tumor validation
# ---------------------------------------------------------------------------

TCGA_FOREST_CONTRASTS = ["tumor_vs_normal_PAIRED", "ER+ vs ER- (clinical IHC)"]


def load_tcga_expression_forest() -> pd.DataFrame:
    """Real effect size + 95% CI for the two strongest, most interpretable
    TCGA expression contrasts (paired tumor-vs-normal, ER+ vs ER-), all 4
    candidates -- excludes the purely descriptive distribution row and the
    unpaired-descriptive row (redundant with the paired row, per the
    table's own `notes` column)."""
    df = pd.read_csv(TCGA_EXPRESSION, sep="\t")
    sub = df[df["comparison"].isin(TCGA_FOREST_CONTRASTS)].copy()
    assert len(sub) == 2 * len(CANDIDATES)
    return sub.reset_index(drop=True)


def load_tcga_clinical_er_adjusted() -> pd.DataFrame:
    """ER+ subgroup, age/stage-adjusted Cox HR -- read for a small secondary
    annotation only, never as the figure's headline claim (several rows
    violate the proportional-hazards assumption, flagged in `notes`)."""
    df = pd.read_csv(TCGA_CLINICAL, sep="\t")
    sub = df[(df["cohort"] == "ER_positive") & (df["model"] == "adjusted_age_stage")]
    assert len(sub) == len(CANDIDATES)
    return sub.reset_index(drop=True)


# ---------------------------------------------------------------------------
# DepMap 26Q1: generic (non-candidate-restricted) gene-effect access, for the
# 28 Gate-1 genes -- a small, deterministic, documented widening of the
# existing candidate-only loader's column-matching logic (symbol parsed
# directly off the "SYMBOL (ENTREZID)" column header instead of a hardcoded
# 4-gene Entrez map), reading the exact same already-downloaded raw files.
# ---------------------------------------------------------------------------

def _gene_col_generic(matrix_columns: pd.Index, symbol: str) -> str | None:
    for col in matrix_columns:
        if col.split(" (")[0] == symbol:
            return col
    return None


def load_depmap_gate1_dependency_summary(release: str = "26Q1") -> pd.DataFrame:
    """For the 28 Gate-1 genes (where present in the DepMap matrix): median
    ER+/luminal Chronos gene-effect and fraction of the 11 ER+/luminal lines
    with dependency probability > 0.5 (the same threshold and same n=11 line
    set already used and frozen in independent_validation_depmap.py /
    poster_figures_data.py). Genes absent from the DepMap release (if any)
    are reported with a NaN summary and an explicit `in_depmap` flag rather
    than silently dropped."""
    from src.independent_validation_depmap_data import (
        has_dependency_probability,
        load_model,
        raw_dir,
    )

    cfg = load_config()
    gate1 = load_crispr_gate1()
    genes = gate1["gene_symbol"].tolist()

    rel_cfg = cfg["independent_validation"]["depmap"]["releases"][release]
    effect_path = raw_dir(cfg, release) / rel_cfg["raw"]["crispr_gene_effect_csv"]
    dep_path = raw_dir(cfg, release) / rel_cfg["raw"]["crispr_gene_dependency_csv"]

    model = load_model(cfg, release)
    luminal_ids = model.index[model["is_er_luminal"]]

    effect_df = pd.read_csv(effect_path, index_col=0)
    dep_df = pd.read_csv(dep_path, index_col=0) if has_dependency_probability(cfg, release) else None
    thresh = cfg["independent_validation"]["depmap"]["strong_dependency_probability_threshold"]

    rows = []
    n_missing = 0
    for gene in genes:
        col = _gene_col_generic(effect_df.columns, gene)
        if col is None:
            n_missing += 1
            rows.append(dict(gene=gene, in_depmap=False, median_chronos_er_luminal=float("nan"),
                              frac_strongly_dependent_er_luminal=float("nan"), n_lines=0))
            continue
        vals = effect_df.loc[effect_df.index.intersection(luminal_ids), col].dropna()
        frac = float("nan")
        if dep_df is not None:
            dep_col = _gene_col_generic(dep_df.columns, gene)
            if dep_col is not None:
                dep_vals = dep_df.loc[dep_df.index.intersection(luminal_ids), dep_col].dropna()
                if len(dep_vals) > 0:
                    frac = float((dep_vals > thresh).mean())
        rows.append(dict(gene=gene, in_depmap=True, median_chronos_er_luminal=float(vals.median()),
                          frac_strongly_dependent_er_luminal=frac, n_lines=len(vals)))
    out = pd.DataFrame(rows)
    logger.info("load_depmap_gate1_dependency_summary: %d/%d Gate-1 genes found in DepMap %s (%d missing)",
                (out["in_depmap"]).sum(), len(genes), release, n_missing)
    return out.merge(gate1.rename(columns={"gene_symbol": "gene"}), on="gene", how="left")


# ---------------------------------------------------------------------------
# GDSC pharmacogenomics -- FDR-significant top associations
# ---------------------------------------------------------------------------

def load_gdsc_top_associations() -> pd.DataFrame:
    df = pd.read_csv(GDSC_TOP_ASSOCIATIONS, sep="\t")
    logger.info("load_gdsc_top_associations: %d rows", len(df))
    return df


# ---------------------------------------------------------------------------
# Tissue liability (expression breadth vs documented functional liability,
# kept as two separate encodings -- never merged into a single score)
# ---------------------------------------------------------------------------

def load_tissue_liability() -> pd.DataFrame:
    df = pd.read_csv(TISSUE_LIABILITY, sep="\t")
    logger.info("load_tissue_liability: %d organ-system rows", len(df))
    return df


def load_normal_tissue_context() -> pd.DataFrame:
    df = pd.read_csv(NORMAL_TISSUE_CONTEXT, sep="\t")
    assert set(df["candidate"]) == set(CANDIDATES)
    return df.set_index("candidate").loc[CANDIDATES].reset_index()
