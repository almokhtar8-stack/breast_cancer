"""Candidate adjudication Phases 7-8: three separate, never-merged
evidence axes (functional/CRISPR, resistance-state RNA, human-tumor), and
a deterministic archetype classification built from those three axes
plus CRISPR direction. No axis is combined into a numeric composite; the
archetype tree is precedence-ordered and documented, first-match-wins,
exactly like `src.cross_dataset_consensus_views.assign_evidence_category`.

Evaluated over the "adjudication candidate pool" -- the union of the
seven MULTIMODAL_STRONG genes, the top-20 resistance-RNA leaders, the
top-20 CRISPR-sensitising leaders, and the near-miss genes from Phases
4-6 -- rather than the full 37,631-gene universe: axis descriptors like
VERY_STRONG..NO_EVIDENCE are only informative for genes with at least
some candidate-level evidence to begin with; scoring the entire universe
this way would mostly report NO_EVIDENCE and add no information.

Data source: the Phase 2-6 adjudication tables plus the frozen
cross-dataset genome-wide tables they were built from (all read-only).
"""

from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd
import yaml

logger = logging.getLogger(__name__)


def _load_config(config_path: str | Path) -> dict:
    with open(config_path) as f:
        return yaml.safe_load(f)


def classify_axis_a_functional(row: pd.Series) -> str:
    """Thresholds derived from the CRISPR screen's own FDR/percentile
    scale: <0.01 and <0.05 are conventional strict/standard significance
    bands; <0.10 is this project's own established Gate-1 threshold
    (config.yaml gate1.fdr_threshold); the WEAK band catches genes in the
    screen's own top decile that miss FDR<0.10 but are not simply noise."""
    fdr = row.get("crispr_fdr")
    pct = row.get("crispr_evidence_percentile")
    if pd.isna(fdr):
        return "NO_EVIDENCE"
    if fdr < 0.01:
        return "VERY_STRONG"
    if fdr < 0.05:
        return "STRONG"
    if fdr < 0.10:
        return "MODERATE"
    if fdr < 0.25 and pd.notna(pct) and pct >= 0.90:
        return "WEAK"
    return "NO_EVIDENCE"


def classify_axis_b_resistance(row: pd.Series) -> str:
    """Uses the three resistance-state datasets' own FDR count and the
    project's existing direction-consensus vocabulary
    (all_up/all_down/majority_up/majority_down/mixed/insufficient) --
    never a new statistic."""
    consensus = row.get("resistance_direction_consensus")
    fdr05 = row.get("resistance_fdr05_count")
    fdr05 = 0 if pd.isna(fdr05) else fdr05
    median_pct = row.get("resistance_median_percentile")
    if consensus == "mixed":
        return "DISCORDANT"
    if consensus in ("insufficient", None) or pd.isna(consensus):
        return "NO_EVIDENCE"
    if fdr05 >= 2 and consensus in ("all_up", "all_down"):
        return "VERY_STRONG"
    if fdr05 >= 1 and consensus in ("all_up", "all_down", "majority_up", "majority_down"):
        return "STRONG"
    if pd.notna(median_pct) and median_pct >= 0.90:
        return "MODERATE"
    if pd.notna(median_pct) and median_pct >= 0.75:
        return "WEAK"
    return "NO_EVIDENCE"


def classify_axis_c_human(row: pd.Series) -> str:
    """Human-tumor evidence strength only -- never implies the same
    mechanism for GSE245601 (acute 12h ex vivo response) and GSE240112
    (primary-vs-recurrent context); `human_evidence_sources` names which
    contributed."""
    epi_fdr = row.get("gse245601_epi_fdr")
    mal_fdr = row.get("gse245601_malignant_fdr")
    gse240112_fdr = row.get("gse240112_tumor_fdr")
    acute_sig = (pd.notna(epi_fdr) and epi_fdr < 0.05) or (pd.notna(mal_fdr) and mal_fdr < 0.05)
    recurrence_sig = pd.notna(gse240112_fdr) and gse240112_fdr < 0.05
    max_pct = max([v for v in (row.get("gse245601_evidence_percentile"), row.get("gse240112_evidence_percentile")) if pd.notna(v)], default=float("nan"))
    if acute_sig and recurrence_sig:
        return "VERY_STRONG"
    if acute_sig or recurrence_sig:
        return "STRONG"
    if pd.notna(max_pct) and max_pct >= 0.90:
        return "MODERATE"
    if pd.notna(max_pct) and max_pct >= 0.75:
        return "WEAK"
    return "NO_EVIDENCE"


def _human_evidence_sources(row: pd.Series) -> str:
    sources = []
    if pd.notna(row.get("gse245601_epi_fdr")) and row["gse245601_epi_fdr"] < 0.05:
        sources.append("GSE245601_acute_epi")
    if pd.notna(row.get("gse245601_malignant_fdr")) and row["gse245601_malignant_fdr"] < 0.05:
        sources.append("GSE245601_acute_malignant")
    if pd.notna(row.get("gse240112_tumor_fdr")) and row["gse240112_tumor_fdr"] < 0.05:
        sources.append("GSE240112_recurrence")
    return "+".join(sources) if sources else "none_FDR<0.05"


def build_candidate_pool(config: dict) -> pd.DataFrame:
    """Union of gene symbols across every adjudication table built so
    far, deduplicated. Returns the pool joined to the full wide-matrix
    row for each gene (one row per unique gene)."""
    tables_dir = Path(config["candidate_adjudication"]["output"]["tables_dir"])
    genes: set[str] = set(config["candidate_adjudication"]["multimodal7"]["genes"])
    for fname in ["top_resistance_genes_exact.tsv", "top_crispr_sensitising_exact.tsv", "multimodal_near_misses.tsv", "multimodal_low_coverage_pattern_genes.tsv"]:
        df = pd.read_csv(tables_dir / fname, sep="\t")
        genes |= set(df["gene"])

    cdx_out = config["cross_dataset_genomewide"]["output"]
    wide_dir = Path(cdx_out["wide_matrix_tsv"]).parent
    wide = pd.read_csv(wide_dir / "all_genes_cross_dataset_evidence_with_ranking.tsv", sep="\t")
    resistance = pd.read_csv(cdx_out["resistance_consensus_tsv"], sep="\t")
    categories = pd.read_csv(wide_dir / "evidence_categories.tsv", sep="\t")
    ranked = pd.read_csv(wide_dir / "global_ranking_eligible.tsv", sep="\t")

    pool = wide.loc[wide["gene"].isin(genes)].merge(resistance, on="gene", how="left").merge(categories, on="gene", how="left").merge(ranked[["gene", "global_rank"]], on="gene", how="left")
    logger.info("build_candidate_pool: %d unique genes across all adjudication leader/near-miss tables", len(pool))
    return pool


def build_three_axis_matrix(pool: pd.DataFrame) -> pd.DataFrame:
    out = pool[["gene", "global_rank", "evidence_category", "crispr_direction"]].copy()
    out["axis_a_functional"] = pool.apply(classify_axis_a_functional, axis=1)
    out["axis_b_resistance"] = pool.apply(classify_axis_b_resistance, axis=1)
    out["axis_c_human"] = pool.apply(classify_axis_c_human, axis=1)
    out["human_evidence_sources"] = pool.apply(_human_evidence_sources, axis=1)
    out = out.sort_values(by=["axis_a_functional", "axis_b_resistance", "gene"], key=lambda c: c.map({"VERY_STRONG": 0, "STRONG": 1, "MODERATE": 2, "WEAK": 3, "DISCORDANT": 3, "NO_EVIDENCE": 4}) if c.name != "gene" else c).reset_index(drop=True)
    logger.info("build_three_axis_matrix: %d genes classified on 3 independent axes", len(out))
    return out


# Listed in actual decision-precedence order (matches assign_archetype's
# if/elif chain exactly -- G is checked second, right after A, not near
# the end; see assign_archetype's own docstring for why).
ARCHETYPE_ORDER = [
    "A_FUNCTIONAL_RESISTANCE_CONVERGENCE",
    "G_CONTEXT_DEPENDENT",
    "B_FUNCTIONAL_HUMAN_CONTEXT",
    "C_FUNCTIONAL_ONLY",
    "D_RESISTANCE_BIOMARKER_PATHWAY",
    "E_HUMAN_RECURRENCE_DOMINANT",
    "F_ACUTE_RESPONSE_DOMINANT",
    "H_LOW_INSUFFICIENT_EVIDENCE",
]

STRONG_LEVELS = {"VERY_STRONG", "STRONG"}


def assign_archetype(row: pd.Series) -> str:
    """Precedence-ordered, first-match-wins (same discipline as
    `assign_evidence_category`): a gene with strong evidence on multiple
    axes is placed in the archetype nearest the top of this list, not the
    one a human might find most narratively appealing.

    G_CONTEXT_DEPENDENT (discordant resistance direction, but real signal
    on another axis) is checked immediately after A -- NOT after B/C/E/F.
    A discordant gene (axis B == "DISCORDANT") can never satisfy A (which
    requires axis B in STRONG_LEVELS, and DISCORDANT is deliberately
    excluded from STRONG_LEVELS), so moving this check up does not change
    A's behavior; checking it any later makes it structurally unreachable,
    since every combination of (functional_ok, human_ok) that could
    satisfy the discordant condition is already fully covered by B, C, E,
    or F individually (caught and regression-tested after the Phase 34
    Codex review found the original ordering made this branch dead code)."""
    a, b, c = row["axis_a_functional"], row["axis_b_resistance"], row["axis_c_human"]
    functional_ok = a in STRONG_LEVELS
    resistance_ok = b in STRONG_LEVELS
    human_ok = c in STRONG_LEVELS
    discordant = b == "DISCORDANT"

    if functional_ok and resistance_ok:
        return "A_FUNCTIONAL_RESISTANCE_CONVERGENCE"
    if discordant and (functional_ok or human_ok):
        return "G_CONTEXT_DEPENDENT"
    if functional_ok and human_ok:
        return "B_FUNCTIONAL_HUMAN_CONTEXT"
    if functional_ok and not resistance_ok and not human_ok:
        return "C_FUNCTIONAL_ONLY"
    if resistance_ok and not functional_ok:
        return "D_RESISTANCE_BIOMARKER_PATHWAY"
    if c == "VERY_STRONG" or (row.get("gse240112_tumor_fdr") is not None and c in STRONG_LEVELS and "GSE240112" in row.get("human_evidence_sources", "") and "GSE245601" not in row.get("human_evidence_sources", "")):
        return "E_HUMAN_RECURRENCE_DOMINANT"
    if c in STRONG_LEVELS and "GSE245601" in row.get("human_evidence_sources", "") and "GSE240112" not in row.get("human_evidence_sources", ""):
        return "F_ACUTE_RESPONSE_DOMINANT"
    return "H_LOW_INSUFFICIENT_EVIDENCE"


def build_archetypes(matrix: pd.DataFrame, pool: pd.DataFrame) -> pd.DataFrame:
    merged = matrix.merge(pool[["gene", "gse240112_tumor_fdr", "human_evidence_sources"]] if "human_evidence_sources" in pool.columns else pool[["gene", "gse240112_tumor_fdr"]], on="gene", how="left", suffixes=("", "_pool"))
    if "human_evidence_sources" not in merged.columns:
        merged["human_evidence_sources"] = matrix["human_evidence_sources"]
    merged["archetype"] = merged.apply(assign_archetype, axis=1)
    out = merged[["gene", "global_rank", "archetype", "axis_a_functional", "axis_b_resistance", "axis_c_human", "evidence_category"]]
    out = out.sort_values(by=["archetype", "global_rank"], na_position="last").reset_index(drop=True)
    logger.info("build_archetypes: %s", out["archetype"].value_counts().to_dict())
    return out


def run_axes(config_path: str | Path = "config/config.yaml") -> dict[str, pd.DataFrame]:
    config = _load_config(config_path)
    pool = build_candidate_pool(config)
    matrix = build_three_axis_matrix(pool)
    archetypes = build_archetypes(matrix, pool)

    out_dir = Path(config["candidate_adjudication"]["output"]["tables_dir"])
    out_dir.mkdir(parents=True, exist_ok=True)
    matrix.to_csv(out_dir / "three_axis_candidate_matrix.tsv", sep="\t", index=False)
    archetypes.to_csv(out_dir / "candidate_archetypes.tsv", sep="\t", index=False)
    pool.to_csv(out_dir / "adjudication_candidate_pool.tsv", sep="\t", index=False)
    return {"pool": pool, "matrix": matrix, "archetypes": archetypes}


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run_axes()
