"""Data access for the POSTER-STORY-V1 figure bank.

This module performs NO new discovery, NO recomputation of any frozen
p-value/FDR/effect-size, and NO re-ranking. It reuses the already-frozen
loaders in `poster_exploration_v2_data.py` / `poster_exploration_v3_data.py`
/ `post_audit_sensitivity_data.py` wherever possible, and adds exactly one
new (but still purely reshaping/visualization) data source: the real
malignant-vs-non-malignant pseudobulk comparison for GSE245601, which was
never used in any prior poster phase despite already being frozen. See
`results/reports/poster_story_v1/DATA_AUDIT.md` for the full source
inventory and `STORY_PLAN.md` for why this phase is built the way it is.

Disclosed transforms in this module (never a new statistical test, no
p-value ever computed here):
  - condition-mean aggregation (e.g. mean TPM across 3 MCF7 replicates)
  - log2 fold-change of a mean vs. its own paired reference mean, in the
    SAME dataset (e.g. mean(TAMR) - mean(MCF7)) -- the standard log2FC
    unit, comparable across studies unlike raw TPM/CPM
  - log2(CPM+1) from frozen raw counts + frozen library sizes, IDENTICAL
    formula already verified against a frozen value in v2/v3
"""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pandas as pd

from src import poster_exploration_v2_data as pv2
from src import poster_exploration_v3_data as pv3
from src import post_audit_sensitivity_data as pad

logger = logging.getLogger(__name__)

FOCUS_FOUR = pv3.FOCUS_FOUR
FOCUS_COLORS = pv3.FOCUS_COLORS
GRAY = pv3.GRAY
LGRAY = pv3.LGRAY
DGRAY = pv3.DGRAY

# ---------------------------------------------------------------------------
# CRISPR discovery
# ---------------------------------------------------------------------------
load_significant_sensitising_hits = pv3.load_significant_sensitising_hits
load_genomewide_crispr = pv3.load_genomewide_crispr

# ---------------------------------------------------------------------------
# Cross-dataset log2FC (already-frozen master table -- used for the
# broader 13-hit heatmap alternate)
# ---------------------------------------------------------------------------
load_cross_dataset_raw = pad.load_cross_dataset_raw


def load_13hit_log2fc_matrix() -> pd.DataFrame:
    """The 4 real dataset log2FC columns for all 13 significant
    sensitising hits, read directly from the already-frozen master
    cross-dataset table -- no recomputation. gse245601_epi_log2fc is the
    acute Control-vs-Tamoxifen track (the same column
    post_audit_sensitivity_data renames to gse245601_acute_log2fc)."""
    sens = load_significant_sensitising_hits()
    cross = load_cross_dataset_raw().set_index("gene")
    cols = ["gse118713_log2fc", "gse111151_log2fc", "gse240112_tumor_log2fc", "gse245601_epi_log2fc"]
    out = cross.loc[sens["gene"], cols].reset_index()
    out.columns = ["gene", "GSE118713", "GSE111151", "GSE240112", "GSE245601"]
    return out


# ---------------------------------------------------------------------------
# Hero heatmap -- 4 focus genes, real paired condition means
# ---------------------------------------------------------------------------

def build_hero_heatmap_pairs() -> pd.DataFrame:
    """Real condition-mean paired rows for all 4 focus genes across the 4
    transcriptomic datasets, plus the log2 fold-change of the second
    (non-reference) condition relative to the first (reference) condition
    in that SAME dataset. Every mean is computed from the frozen per-
    sample/per-tumour/per-patient values already loaded by v2/v3 -- no new
    statistical test, no p-value computed here (significance for the
    broader gene set is read separately from the frozen master table, see
    `load_13hit_log2fc_matrix`)."""
    rows = []

    # GSE118713: MCF7 (baseline) vs TAMR (resistant), TPM
    g118713 = pv2.load_gse118713_focus_gene_samples()
    for gene in FOCUS_FOUR:
        sub = g118713[g118713["gene_symbol"] == gene]
        mcf7 = sub.loc[sub["condition"] == "MCF7", "tpm"].mean()
        tamr = sub.loc[sub["condition"] == "TAMR", "tpm"].mean()
        rows.append(dict(dataset="GSE118713", dataset_label="Resistance model (cell line)",
                          ref_label="MCF7 (baseline)", cmp_label="TAMR (resistant)",
                          gene=gene, ref_value=mcf7, cmp_value=tamr,
                          log2fc=float(np.log2(tamr + 1) - np.log2(mcf7 + 1))))

    # GSE111151: Parental (4 backgrounds) vs Resistant (7 sublines), log2CPM
    g111151 = pv3.load_gse111151_focus_gene_samples()
    for gene in FOCUS_FOUR:
        sub = g111151[g111151["gene_symbol"] == gene]
        parental = sub.loc[sub["status"] == "parental", "log2cpm"].mean()
        resistant = sub.loc[sub["status"] == "resistant", "log2cpm"].mean()
        rows.append(dict(dataset="GSE111151", dataset_label="Resistance model (cell line)",
                          ref_label="Parental (4 lines)", cmp_label="Resistant (7 sublines)",
                          gene=gene, ref_value=parental, cmp_value=resistant,
                          log2fc=float(resistant - parental)))

    # GSE240112: Primary (n=3) vs Recurrent (n=3), UNPAIRED, log2CPM
    g240112 = pv3.load_gse240112_focus_gene_tumours()
    for gene in FOCUS_FOUR:
        sub = g240112[g240112["gene"] == gene]
        primary = sub.loc[sub["group"] == "PT", "log2cpm"].mean()
        recurrent = sub.loc[sub["group"] == "RT", "log2cpm"].mean()
        rows.append(dict(dataset="GSE240112", dataset_label="Recurrence (human, unpaired)",
                          ref_label="Primary (n=3)", cmp_label="Recurrent (n=3)",
                          gene=gene, ref_value=primary, cmp_value=recurrent,
                          log2fc=float(recurrent - primary)))

    # GSE245601 acute: Control vs Tamoxifen (12h), patient-matched, log2(CPM+1)
    g245601 = _load_gse245601_acute_means()
    for gene in FOCUS_FOUR:
        control, tamoxifen = g245601[gene]
        rows.append(dict(dataset="GSE245601", dataset_label="Acute 12h (human, patient-matched)",
                          ref_label="Control", cmp_label="Tamoxifen (12h)",
                          gene=gene, ref_value=control, cmp_value=tamoxifen,
                          log2fc=float(tamoxifen - control)))

    return pd.DataFrame(rows)


def _load_gse245601_acute_means() -> dict[str, tuple[float, float]]:
    df = pv2.load_gse245601_paired_focus_genes()
    out = {}
    for gene in FOCUS_FOUR:
        sub = df[df["gene"] == gene]
        control = sub.loc[sub["condition"] == "Control", "log2_expr"].mean()
        tamoxifen = sub.loc[sub["condition"] == "Tamoxifen", "log2_expr"].mean()
        out[gene] = (float(control), float(tamoxifen))
    return out


DATASET_ORDER = ["GSE118713", "GSE111151", "GSE240112", "GSE245601"]

# ---------------------------------------------------------------------------
# Pathway convergence (reused from v3, no change)
# ---------------------------------------------------------------------------
load_pathway_trajectories = pv3.load_pathway_trajectories
HERO_PATHWAYS = pv3.HERO_PATHWAYS

# ---------------------------------------------------------------------------
# Structural / pharmacological comparison (reused from v2/v3, no change)
# ---------------------------------------------------------------------------
kdm1a_tlk2_structure_paths = pv3.kdm1a_tlk2_structure_paths
usp34_structure_paths = pv3.usp34_structure_paths
load_structural_tractability_audit = pv3.load_structural_tractability_audit

# ---------------------------------------------------------------------------
# Disease / clinical context -- recurrence (reused) + malignant vs. non
# malignant (NEW this phase)
# ---------------------------------------------------------------------------
load_gse240112_focus_gene_tumours = pv3.load_gse240112_focus_gene_tumours

MALIGNANT_COUNTS = Path("results/tables/gse245601_pseudobulk/malignant_vs_nonmalignant_counts.tsv.gz")
MALIGNANT_METADATA = Path("results/tables/gse245601_pseudobulk/malignant_vs_nonmalignant_metadata.tsv")
MALIGNANT_CANDIDATES = Path("results/tables/gse245601_candidate_integration/malignant_vs_nonmalignant_candidates.tsv")


def load_malignant_vs_nonmalignant_per_patient() -> pd.DataFrame:
    """Real per-patient log2(CPM+1) values for all 4 focus genes,
    malignant vs. non-malignant epithelial cells (5 patients: Tumor_02,
    03, 07, 09, 10), computed from the frozen genome-wide raw-count
    matrix + frozen per-sample library sizes -- the IDENTICAL formula
    already verified in v2/v3 against a frozen delta value (re-verified
    here against `malignant_vs_nonmalignant_candidates.tsv`'s own
    mean_delta, see tests/test_poster_story_v1.py)."""
    counts = pd.read_csv(MALIGNANT_COUNTS, sep="\t")
    meta = pd.read_csv(MALIGNANT_METADATA, sep="\t").set_index("sample_id")
    patients = sorted(meta["patient"].unique())

    rows = []
    for gene in FOCUS_FOUR:
        grow = counts[counts["gene"] == gene]
        if len(grow) == 0:
            continue
        grow = grow.iloc[0]
        for patient in patients:
            for status in ("malignant", "nonmalignant"):
                sample_id = f"{patient}_{status}"
                if sample_id not in meta.index:
                    continue
                raw = float(grow[sample_id])
                lib = float(meta.loc[sample_id, "total_library_size"])
                log2_expr = float(np.log2((raw / lib * 1e6) + 1))
                rows.append(dict(patient=patient, gene=gene, status=status, log2_expr=log2_expr))
    return pd.DataFrame(rows)


def build_malignant_vs_nonmalignant_paired_delta() -> pd.DataFrame:
    """Real per-patient PAIRED delta (malignant - non-malignant, same
    patient, same tissue sample) for all 4 focus genes -- a disclosed
    subtraction, not a new statistical test (the frozen significance test
    for this comparison lives in `malignant_vs_nonmalignant_candidates.tsv`,
    read separately, never recomputed here)."""
    df = load_malignant_vs_nonmalignant_per_patient()
    rows = []
    for gene in FOCUS_FOUR:
        sub = df[df["gene"] == gene]
        mal = sub[sub["status"] == "malignant"].set_index("patient")["log2_expr"]
        non = sub[sub["status"] == "nonmalignant"].set_index("patient")["log2_expr"]
        common = mal.index.intersection(non.index)
        for patient in common:
            rows.append(dict(gene=gene, patient=patient, delta_log2=float(mal[patient] - non[patient])))
    return pd.DataFrame(rows)


def build_gse240112_recurrence_delta() -> pd.DataFrame:
    """Real per-tumour delta relative to that gene's own primary-group
    mean (GSE240112 is unpaired, so the group mean is the reference --
    same disclosed centering transform already used and verified in v3)."""
    return pv3.build_gse240112_delta_from_primary_mean()


def load_malignant_vs_nonmalignant_frozen_delta() -> pd.DataFrame:
    """The already-frozen per-gene delta/p-value (13 genes), read
    unmodified from `malignant_vs_nonmalignant_candidates.tsv`."""
    return pd.read_csv(MALIGNANT_CANDIDATES, sep="\t")


# ---------------------------------------------------------------------------
# Baseline DepMap dependency (reused from v2/v3, no change)
# ---------------------------------------------------------------------------
load_depmap_effect_focus_four = pv3.load_depmap_effect_focus_four
load_depmap_model_names = pv3.load_depmap_model_names

# ---------------------------------------------------------------------------
# Final synthesis (reused from v3)
# ---------------------------------------------------------------------------
load_rule0_rule1 = pv3.load_rule0_rule1

# ---------------------------------------------------------------------------
# Network backup (real, honest, explicitly backup-only)
# ---------------------------------------------------------------------------
load_direct_neighbors = pv2.load_direct_neighbors
