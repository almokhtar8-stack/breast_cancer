"""Data access for the EXPLORATION-V3 poster-grade figure bank.

Thin wrapper over the already-frozen loaders built in
`src/poster_exploration_v2_data.py` and `src/post_audit_sensitivity_data.py`
-- this module performs NO new discovery, NO recomputation, and NO
re-ranking of any kind. It exists only to give v3's plotting code a single,
clearly-named import surface and to document, in one place, exactly which
v2/post-audit function backs each v3 figure. See
`results/reports/poster_exploration_v2/DATA_FOR_VISUALIZATION_AUDIT.md`
for the full original source-table provenance (unchanged in v3).

Design brief for v3 (2026-08-16): the v2 bank was judged scientifically
sound but visually too close to an internal analysis/report deck (too many
small subplots, too much grey, too much footnote text). v3 rebuilds a
small set of 6 poster-grade hero figures (+ up to 2 alternates) from the
same frozen data, with a stricter visual design system (see
`src/poster_exploration_v3_visualization.py` module docstring).
"""

from __future__ import annotations

import logging

import numpy as np
import pandas as pd

from src import poster_exploration_v2_data as pv2
from src import post_audit_sensitivity_data as pad

logger = logging.getLogger(__name__)

# Exact palette requested for v3: KDM1A=orange, TLK2=light blue, USP34=strong
# blue, VEZF1=gold/yellow. These are the SAME hex values already frozen in
# pv2.FOCUS_COLORS (Okabe-Ito colorblind-safe set, validated with the
# dataviz skill's palette validator in an earlier phase) -- re-exported
# here unchanged, not redefined, so v3 stays pixel-identical in identity
# color to v2/poster_final.
FOCUS_FOUR = pv2.FOCUS_FOUR
FOCUS_COLORS = pv2.FOCUS_COLORS
GRAY = "#9a9a9a"
LGRAY = "#e2e2e2"
DGRAY = "#1f1f1f"

# ---------------------------------------------------------------------------
# Figure 1 -- CRISPR discovery
# ---------------------------------------------------------------------------
load_significant_sensitising_hits = pv2.load_significant_sensitising_hits
load_genomewide_crispr = pv2.load_genomewide_crispr
load_blind_control_row = pv2.load_blind_control_row

# ---------------------------------------------------------------------------
# Figure 2 -- transcriptomic corroboration (GSE111151 + GSE240112)
# ---------------------------------------------------------------------------
load_gse111151_focus_gene_samples = pv2.load_gse111151_focus_gene_samples
load_gse240112_focus_gene_tumours = pv2.load_gse240112_focus_gene_tumours
GSE111151_METADATA = pv2.GSE111151_METADATA


def build_gse111151_delta_from_parental() -> pd.DataFrame:
    """Per-derivative log2 fold-change relative to ITS OWN real parental
    line (real recorded pairing, from `paired_parental_sample_id` --
    never invented). This is a disclosed centering transform for
    visualization only (subtraction, not a new statistical test); no
    p-value/FDR is computed here. Making this transform lets all 4 focus
    genes (very different absolute baseline expression) share one
    common, meaningful y-axis ("log2 fold-change vs. own parental") in a
    single combined poster panel, instead of 4 separate small-multiple
    panels each on its own absolute scale."""
    long = load_gse111151_focus_gene_samples()
    parental = long[long["status"] == "parental"].set_index(["gene_symbol", "parental_line"])["log2cpm"]
    resistant = long[long["status"] == "resistant"].copy()
    resistant["parental_log2cpm"] = resistant.apply(
        lambda r: parental.loc[(r["gene_symbol"], r["parental_line"])], axis=1
    )
    resistant["delta_log2cpm"] = resistant["log2cpm"] - resistant["parental_log2cpm"]
    return resistant[["gene_symbol", "sample_id", "parental_line", "derivative_id", "delta_log2cpm"]].reset_index(drop=True)


def build_gse240112_delta_from_primary_mean() -> pd.DataFrame:
    """Per-tumour log2 fold-change relative to that gene's own PRIMARY-
    group mean (GSE240112 is unpaired -- there is no per-tumour partner
    to subtract, so the group mean is used as the reference point
    instead, a standard, disclosed centering transform, not a new
    statistical test; no p-value/FDR is computed here). Primary tumours
    end up scattered around ~0 by construction (their own group's
    spread); recurrent tumours show real deviation from that reference.
    Makes all 4 focus genes share one common, meaningful y-axis in a
    single combined panel."""
    df = load_gse240112_focus_gene_tumours()
    primary_mean = df[df["group"] == "PT"].groupby("gene")["log2cpm"].mean()
    out = df.copy()
    out["primary_mean_log2cpm"] = out["gene"].map(primary_mean)
    out["delta_log2cpm"] = out["log2cpm"] - out["primary_mean_log2cpm"]
    return out[["gene", "sample_id", "group", "delta_log2cpm"]].reset_index(drop=True)

# ---------------------------------------------------------------------------
# Figure 3 -- pathway convergence
# ---------------------------------------------------------------------------
load_pathway_trajectories = pv2.load_pathway_trajectories
HERO_PATHWAYS = pv2.HERO_PATHWAYS
DATASET_ORDER = pv2.DATASET_ORDER
DATASET_LABELS = pv2.DATASET_LABELS

# ---------------------------------------------------------------------------
# Figure 4 -- post-audit interpretation framework
# ---------------------------------------------------------------------------
load_selection_rule_sensitivity = pad.build_selection_rule_sensitivity
load_evidence_matrix_13 = pv2.load_evidence_matrix_13
ORIGINAL_FOUR = pad.ORIGINAL_FOUR


def load_rule0_rule1() -> tuple[pd.DataFrame, pd.DataFrame]:
    """Rule 0 (original frozen gate, RNA eligibility required) and Rule 1
    (CRISPR-only, no RNA gate), read from the already-frozen
    `build_selection_rule_sensitivity()` long table, reshaped to two wide
    per-gene frames -- the exact same reshape already used in v2's E3
    figure, repeated here for a v3-native call signature."""
    long = load_selection_rule_sensitivity()
    rule0 = long[long["rule"] == "RULE_0_original_frozen_gate"].set_index("gene")
    rule1 = long[long["rule"] == "RULE_1_crispr_only_no_rna_gate"].set_index("gene")
    return rule0, rule1

# ---------------------------------------------------------------------------
# Figure 5 -- human / DepMap context
# ---------------------------------------------------------------------------
load_depmap_effect_focus_four = pv2.load_depmap_effect_focus_four
load_depmap_model_names = pv2.load_depmap_model_names

# ---------------------------------------------------------------------------
# Figure 6 -- structural / pharmacological comparison
# ---------------------------------------------------------------------------
kdm1a_tlk2_structure_paths = pv2.kdm1a_tlk2_structure_paths
usp34_structure_paths = pv2.usp34_structure_paths
load_structural_tractability_audit = pv2.load_structural_tractability_audit
