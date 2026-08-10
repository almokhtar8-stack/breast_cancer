"""Descriptive evidence summary and direction-aware evidence-class hierarchy
for the 28 Gate-1 CRISPR hits, using only the frozen master integration
table, plus a strictly separate published-benchmark check for PAICS.

Source: ``results/tables/crispr_gse118713_master_table.tsv`` (config
``crispr_gse118713_integration.output.master_table_tsv``), produced by
``src.crispr_gse118713_integration.run_integration``. This module runs no
CRISPR refit, no limma refit, and no new statistical test -- it only reads
the frozen table and applies deterministic, pre-declared classification
rules to columns that already exist in it. The PAICS benchmark additionally
reads the frozen full 19,103-gene CRISPR fit (``gate1.labels_path``) and
the frozen GSE118713 differential-expression table
(``gse118713_phase2b.limma.differential_expression_tsv_gz``), read-only.

CRISPR effect-direction convention (verified, not assumed, from
PREANALYSIS.md S2 and the implemented model in
``src.labels.fit_gene_treatment_effects``/``_to_long``: the design matrix
codes ``treatment=0`` for the E2-alone arm and ``treatment=1`` for the
E2+4-OHT arm, so the fitted coefficient is the E2+4-OHT-vs-E2-alone
difference on the guide-blocked log2-normalised-count scale):

- ``effect_size < 0`` (``sensitising_knockout``): the guide/gene is more
  DEPLETED under 4-OHT than under E2 alone -- knockout is associated with
  increased tamoxifen sensitivity in this assay.
- ``effect_size > 0`` (``tolerance_associated_knockout``): the guide/gene
  is more ENRICHED under 4-OHT than under E2 alone -- knockout is
  associated with relative fitness/tolerance under 4-OHT.
- ``effect_size == 0`` (``indeterminate``): no directional call.

Wording is measurement-proximal throughout: a gene is never described as
"causing" sensitivity or resistance, and per PREANALYSIS.md S1's claim
boundary, no direction is described as "restoring sensitivity" or
"reversing resistance" -- the screen was run in a drug-tolerant parental
clone (MCF7-V), not an acquired-resistant derivative.

PRIMARY TRANSLATIONAL QUESTION this module answers: which genes could
potentially be inhibited to enhance tamoxifen sensitivity, and which of
those are supported by expression changes in tamoxifen-resistant cells?
Only ``sensitising_knockout`` (negative-effect) genes are candidates for
the first half of that question -- pharmacological inhibition could
potentially phenocopy the sensitising knockout effect and requires
experimental validation, so only genes whose knockout is itself associated
with increased sensitivity are considered for this question.
``tolerance_associated_knockout`` (positive-effect) genes are retained and
fully documented as tolerance/resistance biology, but are explicitly NOT
inhibition candidates.

Evidence-class hierarchy (deterministic; declared before any gene identity
was inspected; see ``classify_evidence_class``) applies ONLY to
``sensitising_knockout`` genes and gives TAMR_vs_MCF7 explicit PRIMARY
standing over TAMR_vs_FASR, which is SECONDARY/contextual only:

- ``PRIMARY_RESISTANCE_SUPPORT``: negative CRISPR effect, DE available,
  TAMR_vs_MCF7 FDR < ``rna_significance_fdr`` (default 0.05). This is the
  primary resistance-expression evidence class, regardless of the
  TAMR_vs_FASR result.
- ``SECONDARY_CONTEXT_SUPPORT``: negative CRISPR effect, DE available, not
  in PRIMARY_RESISTANCE_SUPPORT, TAMR_vs_FASR FDR < ``rna_significance_fdr``.
  This is contextual endocrine-resistance evidence only -- it is never
  described as tamoxifen-specific evidence, since FASR resistance is to a
  different endocrine agent (fulvestrant).
- ``NO_SIGNIFICANT_RNA_SUPPORT``: negative CRISPR effect, DE available,
  significant in neither contrast.
- ``RNA_UNAVAILABLE``: negative CRISPR effect, DE unavailable (the gene
  was removed by the GSE118713 expression filter before any DE model was
  fit).

Differential expression is association, not causality, and is never
described as such. TAMR_vs_FASR significance is not treated as proof of
tamoxifen specificity, and TAMR_vs_FASR nonsignificance is never used as
evidence of specificity, per PREANALYSIS.md's 2026-08-05 Phase 2B
amendment and ``src.gse118713_tamr_specificity``'s design. No composite,
weighted, or otherwise combined score is computed anywhere in this
module; evidence-class membership and shortlist membership are both pure
boolean classifications on pre-declared thresholds.

PAICS benchmark: PAICS is the principal sensitisation target highlighted
experimentally in Hany et al. 2023 and is used here strictly as a
published, read-only benchmark check (``build_paics_benchmark``). It is
NOT one of the 28 Gate-1 FDR<0.1 hits, is never merged into the master
table or any candidate/shortlist/secondary table, and is never used to
tune ``rna_significance_fdr`` or any other threshold in this module.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

logger = logging.getLogger(__name__)

MAPPING_UNIQUE_FILTERED_OUT = "unique_filtered_out"

DIRECTION_SENSITISING = "sensitising_knockout"
DIRECTION_TOLERANCE = "tolerance_associated_knockout"
DIRECTION_INDETERMINATE = "indeterminate"

DIRECTION_DESCRIPTIONS: dict[str, str] = {
    DIRECTION_SENSITISING: (
        "Knockout is relatively depleted under E2+4-OHT versus E2 alone "
        "(associated with increased tamoxifen sensitivity)."
    ),
    DIRECTION_TOLERANCE: (
        "Knockout is relatively enriched under E2+4-OHT versus E2 alone "
        "(associated with relative fitness/tolerance under 4-OHT)."
    ),
    DIRECTION_INDETERMINATE: "CRISPR effect_size is exactly zero; no directional call.",
}

EVIDENCE_CLASS_PRIMARY = "PRIMARY_RESISTANCE_SUPPORT"
EVIDENCE_CLASS_SECONDARY = "SECONDARY_CONTEXT_SUPPORT"
EVIDENCE_CLASS_NO_SIGNIFICANT_RNA = "NO_SIGNIFICANT_RNA_SUPPORT"
EVIDENCE_CLASS_RNA_UNAVAILABLE = "RNA_UNAVAILABLE"
EVIDENCE_CLASS_NOT_APPLICABLE = "not_applicable_not_a_sensitising_knockout"

# Explicit priority order -- TAMR_vs_MCF7 (PRIMARY) outranks TAMR_vs_FASR
# (SECONDARY) by construction, not by string/lexical accident.
EVIDENCE_CLASS_ORDER: tuple[str, ...] = (
    EVIDENCE_CLASS_PRIMARY,
    EVIDENCE_CLASS_SECONDARY,
    EVIDENCE_CLASS_NO_SIGNIFICANT_RNA,
    EVIDENCE_CLASS_RNA_UNAVAILABLE,
)
EVIDENCE_CLASS_SORT_ORDER: dict[str, int] = {cls: i for i, cls in enumerate(EVIDENCE_CLASS_ORDER)}

EVIDENCE_CLASS_LABELS: dict[str, str] = {
    EVIDENCE_CLASS_PRIMARY: (
        "negative CRISPR effect, DE available, significant TAMR_vs_MCF7 (primary resistance-expression evidence)"
    ),
    EVIDENCE_CLASS_SECONDARY: (
        "negative CRISPR effect, DE available, not significant in TAMR_vs_MCF7, significant TAMR_vs_FASR only "
        "(secondary/contextual endocrine-resistance evidence, not tamoxifen-specific)"
    ),
    EVIDENCE_CLASS_NO_SIGNIFICANT_RNA: "negative CRISPR effect, DE available, no significant RNA support in either contrast",
    EVIDENCE_CLASS_RNA_UNAVAILABLE: "negative CRISPR effect, RNA unavailable (expression-filtered out)",
    EVIDENCE_CLASS_NOT_APPLICABLE: (
        "not applicable -- gene is not a sensitising_knockout, so it is not scored for the "
        "inhibition-target question"
    ),
}

RNA_DIRECTION_ELEVATED = "elevated_in_TAMR"
RNA_DIRECTION_REDUCED = "reduced_in_TAMR"
RNA_DIRECTION_UNCHANGED = "unchanged_in_TAMR"
RNA_DIRECTION_NOT_APPLICABLE = "not_applicable"
RNA_DIRECTION_NA = "na"

BENCHMARK_GENE_SYMBOL = "PAICS"
BENCHMARK_LABEL = "published_benchmark_not_gate1_hit"

REQUIRED_MASTER_COLUMNS: tuple[str, ...] = (
    "gene_symbol",
    "crispr_effect_size",
    "crispr_fdr",
    "mapping_status",
    "tamr_vs_mcf7_log2fc",
    "tamr_vs_mcf7_fdr",
    "tamr_vs_fasr_log2fc",
    "tamr_vs_fasr_fdr",
    "mcf7_baseline_log2_tpm_plus1",
    "gse118713_de_na_reason",
)

REQUIRED_LABELS_COLUMNS: tuple[str, ...] = ("gene", "n_guides", "effect_size", "fdr")
REQUIRED_DE_COLUMNS: tuple[str, ...] = ("gene_id", "gene_symbol", "log2fc", "fdr", "contrast")


def _load_config(config_path: str | Path) -> dict:
    with open(config_path) as f:
        return yaml.safe_load(f)


@dataclass(frozen=True)
class EvidenceSummaryConfig:
    """Resolved, config-driven paths. No hardcoded paths."""

    master_table_tsv: Path
    expected_n_hits: int
    rna_significance_fdr: float
    labels_parquet: Path
    gse118713_de_tsv_gz: Path
    evidence_summary_tsv: Path
    sensitisation_candidates_tsv: Path
    tolerance_hits_tsv: Path
    shortlist_tsv: Path
    secondary_candidates_tsv: Path
    paics_benchmark_tsv: Path
    plot_crispr_effects_tsv: Path
    plot_gse118713_contrasts_tsv: Path
    plot_shortlist_matrix_tsv: Path

    @classmethod
    def from_config(cls, config: dict) -> "EvidenceSummaryConfig":
        integ = config["crispr_gse118713_integration"]
        ces = config["candidate_evidence_summary"]
        out = ces["output"]
        return cls(
            master_table_tsv=Path(integ["output"]["master_table_tsv"]),
            expected_n_hits=int(integ["expected_n_hits"]),
            rna_significance_fdr=float(ces["rna_significance_fdr"]),
            labels_parquet=Path(config["gate1"]["labels_path"]),
            gse118713_de_tsv_gz=Path(config["gse118713_phase2b"]["limma"]["differential_expression_tsv_gz"]),
            evidence_summary_tsv=Path(out["evidence_summary_tsv"]),
            sensitisation_candidates_tsv=Path(out["sensitisation_candidates_tsv"]),
            tolerance_hits_tsv=Path(out["tolerance_hits_tsv"]),
            shortlist_tsv=Path(out["shortlist_tsv"]),
            secondary_candidates_tsv=Path(out["secondary_candidates_tsv"]),
            paics_benchmark_tsv=Path(out["paics_benchmark_tsv"]),
            plot_crispr_effects_tsv=Path(out["plot_crispr_effects_tsv"]),
            plot_gse118713_contrasts_tsv=Path(out["plot_gse118713_contrasts_tsv"]),
            plot_shortlist_matrix_tsv=Path(out["plot_shortlist_matrix_tsv"]),
        )


def load_master_table(cfg: EvidenceSummaryConfig) -> pd.DataFrame:
    """Load and structurally validate the frozen master table.

    Validates only what this module reads -- required columns, row count,
    and DE availability accounting -- it does not re-verify the master
    table's own checksums or upstream provenance, which are the
    integration module's responsibility.
    """
    df = pd.read_csv(cfg.master_table_tsv, sep="\t")
    missing = [c for c in REQUIRED_MASTER_COLUMNS if c not in df.columns]
    if missing:
        raise ValueError(f"master table missing required columns: {missing}")
    if len(df) != cfg.expected_n_hits:
        raise ValueError(f"expected exactly {cfg.expected_n_hits} rows in master table, got {len(df)}")

    duplicate_symbols = df.loc[df["gene_symbol"].duplicated(), "gene_symbol"].tolist()
    if duplicate_symbols:
        raise ValueError(f"master table contains duplicate gene_symbol values: {duplicate_symbols}")

    if not np.isfinite(df["crispr_effect_size"].to_numpy(dtype=float)).all():
        raise ValueError("master table contains non-finite crispr_effect_size values")
    if not np.isfinite(df["crispr_fdr"].to_numpy(dtype=float)).all():
        raise ValueError("master table contains non-finite crispr_fdr values")
    if not df["crispr_fdr"].between(0.0, 1.0).all():
        raise ValueError("master table contains crispr_fdr values outside [0, 1]")

    n_filtered_out = int((df["mapping_status"] == MAPPING_UNIQUE_FILTERED_OUT).sum())
    n_de_available = int(df["gse118713_de_na_reason"].isna().sum())
    if n_de_available + n_filtered_out != len(df):
        raise ValueError(
            f"DE-available ({n_de_available}) + filtered-out ({n_filtered_out}) does not "
            f"account for all {len(df)} rows -- an unexplained DE gap is present"
        )

    de_available = df.loc[df["gse118713_de_na_reason"].isna()]
    rna_value_columns = ("tamr_vs_mcf7_log2fc", "tamr_vs_mcf7_fdr", "tamr_vs_fasr_log2fc", "tamr_vs_fasr_fdr")
    for col in rna_value_columns:
        if not np.isfinite(de_available[col].to_numpy(dtype=float)).all():
            raise ValueError(f"DE-available rows contain non-finite {col} values")
    for col in ("tamr_vs_mcf7_fdr", "tamr_vs_fasr_fdr"):
        if not de_available[col].between(0.0, 1.0).all():
            raise ValueError(f"DE-available rows contain {col} values outside [0, 1]")

    logger.info(
        "load_master_table: read %d rows (%d DE-available, %d filtered out)",
        len(df),
        n_de_available,
        n_filtered_out,
    )
    return df


def classify_crispr_direction(effect_size: float) -> str:
    """Verified sign convention -- see module docstring. Never inferred
    from gene identity or memory."""
    if effect_size < 0:
        return DIRECTION_SENSITISING
    if effect_size > 0:
        return DIRECTION_TOLERANCE
    return DIRECTION_INDETERMINATE


def _log2fc_direction(log2fc: float) -> str:
    if pd.isna(log2fc):
        return "na"
    if log2fc > 0:
        return "up"
    if log2fc < 0:
        return "down"
    return "unchanged"


def classify_rna_direction(log2fc: float) -> str:
    """Elevated/reduced-in-TAMR annotation. Meaningful only within
    PRIMARY_RESISTANCE_SUPPORT (see ``build_evidence_summary``), where it
    reports the TAMR_vs_MCF7 direction for that primary-evidence gene."""
    if pd.isna(log2fc):
        return RNA_DIRECTION_NA
    if log2fc > 0:
        return RNA_DIRECTION_ELEVATED
    if log2fc < 0:
        return RNA_DIRECTION_REDUCED
    return RNA_DIRECTION_UNCHANGED


def classify_evidence_class(row: pd.Series, rna_significance_fdr: float) -> str:
    """Deterministic evidence-class hierarchy for sensitising_knockout
    genes, with TAMR_vs_MCF7 given explicit primary standing over
    TAMR_vs_FASR (see module docstring). Returns
    ``EVIDENCE_CLASS_NOT_APPLICABLE`` for every other direction, since
    tolerance-associated/indeterminate genes are not candidates for the
    inhibition-target question this hierarchy answers."""
    if row["crispr_direction"] != DIRECTION_SENSITISING:
        return EVIDENCE_CLASS_NOT_APPLICABLE
    if not pd.isna(row["gse118713_de_na_reason"]):
        return EVIDENCE_CLASS_RNA_UNAVAILABLE
    mcf7_sig = row["tamr_vs_mcf7_fdr"] < rna_significance_fdr
    if mcf7_sig:
        return EVIDENCE_CLASS_PRIMARY
    fasr_sig = row["tamr_vs_fasr_fdr"] < rna_significance_fdr
    if fasr_sig:
        return EVIDENCE_CLASS_SECONDARY
    return EVIDENCE_CLASS_NO_SIGNIFICANT_RNA


def describe_rna_pattern(row: pd.Series, rna_significance_fdr: float) -> str:
    """Neutral, deterministic, template-based description -- association
    language only, no causal claims, no tamoxifen-specificity claim, and
    TAMR_vs_FASR nonsignificance is never presented as evidence of
    specificity (see module docstring)."""
    if not pd.isna(row["gse118713_de_na_reason"]):
        return "No GSE118713 differential-expression result available (removed by the expression filter)."

    mcf7_sig = row["tamr_vs_mcf7_fdr"] < rna_significance_fdr
    fasr_sig = row["tamr_vs_fasr_fdr"] < rna_significance_fdr
    mcf7_dir = _log2fc_direction(row["tamr_vs_mcf7_log2fc"])
    fasr_dir = _log2fc_direction(row["tamr_vs_fasr_log2fc"])

    if mcf7_sig and fasr_sig:
        return (
            f"Significant TAMR-associated expression difference vs both MCF7 and FASR "
            f"(FDR<{rna_significance_fdr}): {mcf7_dir} vs MCF7, {fasr_dir} vs FASR."
        )
    if mcf7_sig:
        return (
            f"Significant TAMR-associated expression difference vs MCF7 only "
            f"(FDR<{rna_significance_fdr}, {mcf7_dir}); TAMR_vs_FASR not significant at this threshold."
        )
    if fasr_sig:
        return (
            f"Significant TAMR-associated expression difference vs FASR only "
            f"(FDR<{rna_significance_fdr}, {fasr_dir}); TAMR_vs_MCF7 not significant at this threshold."
        )
    return f"No significant TAMR-associated expression difference detected in either contrast at FDR<{rna_significance_fdr}."


def build_evidence_summary(master_df: pd.DataFrame, rna_significance_fdr: float) -> pd.DataFrame:
    """Attach verified direction, the primary/secondary evidence-class
    hierarchy (negative-effect genes only), a primary-RNA-direction
    annotation (primary-evidence genes only -- refers specifically to the
    TAMR_vs_MCF7 contrast, not TAMR_vs_FASR), and a neutral RNA-pattern
    description to every one of the 28 rows. No score, no ranking, no row
    is dropped."""
    out = master_df.copy()
    out["crispr_direction"] = out["crispr_effect_size"].apply(classify_crispr_direction)
    out["crispr_direction_description"] = out["crispr_direction"].map(DIRECTION_DESCRIPTIONS)
    out["evidence_class"] = out.apply(lambda r: classify_evidence_class(r, rna_significance_fdr), axis=1)
    out["evidence_class_label"] = out["evidence_class"].map(EVIDENCE_CLASS_LABELS)
    out["primary_rna_direction"] = out.apply(
        lambda r: classify_rna_direction(r["tamr_vs_mcf7_log2fc"])
        if r["evidence_class"] == EVIDENCE_CLASS_PRIMARY
        else RNA_DIRECTION_NOT_APPLICABLE,
        axis=1,
    )
    out["rna_pattern_description"] = out.apply(lambda r: describe_rna_pattern(r, rna_significance_fdr), axis=1)

    columns = [
        "gene_symbol",
        "crispr_effect_size",
        "crispr_fdr",
        "crispr_direction",
        "crispr_direction_description",
        "tamr_vs_mcf7_log2fc",
        "tamr_vs_mcf7_fdr",
        "tamr_vs_fasr_log2fc",
        "tamr_vs_fasr_fdr",
        "mcf7_baseline_log2_tpm_plus1",
        "gse118713_de_na_reason",
        "evidence_class",
        "evidence_class_label",
        "primary_rna_direction",
        "rna_pattern_description",
    ]
    out = out[columns]

    direction_counts = out["crispr_direction"].value_counts().to_dict()
    class_counts = out.loc[out["crispr_direction"] == DIRECTION_SENSITISING, "evidence_class"].value_counts().to_dict()
    logger.info(
        "build_evidence_summary: %d genes classified; direction counts=%s; "
        "evidence-class counts (negative-effect genes only)=%s",
        len(out),
        direction_counts,
        class_counts,
    )
    return out


def split_by_direction(evidence_df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Two clearly separate biological groups: sensitising_knockout
    (primary inhibition-target candidates) and tolerance_associated_
    knockout (retained/documented, explicitly not inhibition candidates
    for the sensitisation question). Indeterminate genes (effect_size
    exactly zero), if any, appear in neither group. Sensitising genes are
    sorted by the explicit evidence-class priority order (PRIMARY before
    SECONDARY before NO_SIGNIFICANT_RNA before RNA_UNAVAILABLE), not by
    string/lexical accident."""
    sensitising = evidence_df.loc[evidence_df["crispr_direction"] == DIRECTION_SENSITISING].copy()
    sensitising["_class_sort"] = sensitising["evidence_class"].map(EVIDENCE_CLASS_SORT_ORDER)
    sensitising = sensitising.sort_values(["_class_sort", "crispr_fdr", "gene_symbol"]).drop(columns="_class_sort")
    sensitising = sensitising.reset_index(drop=True)

    tolerance = evidence_df.loc[evidence_df["crispr_direction"] == DIRECTION_TOLERANCE].copy()
    tolerance = tolerance.sort_values(["crispr_fdr", "gene_symbol"]).reset_index(drop=True)
    tolerance["not_inhibition_candidate_reason"] = (
        "Positive CRISPR effect (tolerance_associated_knockout): knockout is associated with a "
        "relative fitness/tolerance advantage under 4-OHT, not increased sensitivity, so it is not "
        "a candidate for inhibition to enhance tamoxifen sensitivity. Retained as tolerance/"
        "resistance-biology evidence only."
    )

    logger.info(
        "split_by_direction: %d sensitising_knockout genes, %d tolerance_associated_knockout genes",
        len(sensitising),
        len(tolerance),
    )
    return sensitising, tolerance


def _primary_shortlist_reason(row: pd.Series) -> str:
    return (
        f"Negative CRISPR effect + resistance-associated expression support (CRISPR FDR="
        f"{row['crispr_fdr']:.3g}). Significant TAMR_vs_MCF7 change (primary resistance-expression "
        f"evidence): log2FC={row['tamr_vs_mcf7_log2fc']:.3g}, FDR={row['tamr_vs_mcf7_fdr']:.3g}, "
        f"{row['primary_rna_direction']}. TAMR_vs_FASR (secondary/contextual only): log2FC="
        f"{row['tamr_vs_fasr_log2fc']:.3g}, FDR={row['tamr_vs_fasr_fdr']:.3g}. This is association "
        f"evidence, not a proven resistance mechanism, and does not indicate restored tamoxifen "
        f"sensitivity."
    )


def _secondary_context_note(row: pd.Series) -> str:
    return (
        f"Negative CRISPR effect (sensitising knockout, CRISPR FDR={row['crispr_fdr']:.3g}); no "
        f"significant TAMR_vs_MCF7 evidence (log2FC={row['tamr_vs_mcf7_log2fc']:.3g}, FDR="
        f"{row['tamr_vs_mcf7_fdr']:.3g}). Contextual, not tamoxifen-specific, endocrine-resistance "
        f"evidence only: significant TAMR_vs_FASR change (log2FC={row['tamr_vs_fasr_log2fc']:.3g}, FDR="
        f"{row['tamr_vs_fasr_fdr']:.3g})."
    )


def select_primary_shortlist(sensitising_df: pd.DataFrame) -> pd.DataFrame:
    """The primary next-validation shortlist is exactly the
    PRIMARY_RESISTANCE_SUPPORT genes -- no fallback to a weaker evidence
    class, no padding to a target size. May be empty or a single gene;
    both are reported transparently rather than forced to a size."""
    subset = sensitising_df.loc[sensitising_df["evidence_class"] == EVIDENCE_CLASS_PRIMARY].copy()
    subset = subset.sort_values(["crispr_fdr", "gene_symbol"]).reset_index(drop=True)
    subset["shortlist_reason"] = subset.apply(_primary_shortlist_reason, axis=1)
    logger.info(
        "select_primary_shortlist: %d gene(s) in PRIMARY_RESISTANCE_SUPPORT selected as the primary "
        "shortlist (no padding, no fallback to a weaker evidence class)",
        len(subset),
    )
    return subset


def select_secondary_candidates(sensitising_df: pd.DataFrame) -> pd.DataFrame:
    """SECONDARY_CONTEXT_SUPPORT genes, kept as a separate, unranked
    later-validation set. Sorted only by gene_symbol -- deliberately not
    sorted by FDR or any other evidence-strength column, so the ordering
    does not imply an internal ranking."""
    subset = sensitising_df.loc[sensitising_df["evidence_class"] == EVIDENCE_CLASS_SECONDARY].copy()
    subset = subset.sort_values(["gene_symbol"]).reset_index(drop=True)
    subset["secondary_context_note"] = subset.apply(_secondary_context_note, axis=1)
    logger.info("select_secondary_candidates: %d gene(s) in SECONDARY_CONTEXT_SUPPORT (unranked)", len(subset))
    return subset


def load_full_crispr_labels(path: str | Path) -> pd.DataFrame:
    """Load the frozen full 19,103-gene CRISPR fit (``gate1.labels_path``).
    Used only for the PAICS benchmark check -- never joined into the
    28-hit master table or any candidate table."""
    df = pd.read_parquet(path)
    missing = [c for c in REQUIRED_LABELS_COLUMNS if c not in df.columns]
    if missing:
        raise ValueError(f"full CRISPR labels table missing required columns: {missing}")
    logger.info("load_full_crispr_labels: read %d fitted genes", len(df))
    return df


def load_gse118713_differential_expression(path: str | Path) -> pd.DataFrame:
    """Load the frozen, long-format GSE118713 differential-expression
    table (one row per gene per contrast). Used only for the PAICS
    benchmark check."""
    df = pd.read_csv(path, sep="\t")
    missing = [c for c in REQUIRED_DE_COLUMNS if c not in df.columns]
    if missing:
        raise ValueError(f"GSE118713 differential expression table missing required columns: {missing}")
    logger.info("load_gse118713_differential_expression: read %d rows", len(df))
    return df


def build_paics_benchmark(labels_df: pd.DataFrame, de_df: pd.DataFrame, master_df: pd.DataFrame) -> pd.DataFrame:
    """Read-only published-benchmark check for PAICS (Hany et al. 2023).

    PAICS is deliberately evaluated outside the 28-hit candidate pipeline:
    it is not a Gate-1 FDR<0.1 hit, must never be merged into the master
    table or any candidate/shortlist/secondary table, and must never be
    used to tune ``rna_significance_fdr`` or any other threshold. This
    function raises if PAICS is ever found in the 28-hit master table,
    since that would signal the benchmark and the candidate pipeline have
    become conflated.
    """
    if (master_df["gene_symbol"] == BENCHMARK_GENE_SYMBOL).any():
        raise ValueError(
            f"{BENCHMARK_GENE_SYMBOL} unexpectedly present in the 28-hit master table -- the published "
            "benchmark must remain strictly separate from Gate-1 candidate discovery"
        )

    crispr_rows = labels_df.loc[labels_df["gene"] == BENCHMARK_GENE_SYMBOL]
    if len(crispr_rows) != 1:
        raise ValueError(
            f"expected exactly one {BENCHMARK_GENE_SYMBOL} row in the full CRISPR fit, found {len(crispr_rows)}"
        )
    crispr_row = crispr_rows.iloc[0]

    de_rows = de_df.loc[de_df["gene_symbol"] == BENCHMARK_GENE_SYMBOL]
    mcf7_rows = de_rows.loc[de_rows["contrast"] == "TAMR_vs_MCF7"]
    fasr_rows = de_rows.loc[de_rows["contrast"] == "TAMR_vs_FASR"]
    if len(mcf7_rows) != 1 or len(fasr_rows) != 1:
        raise ValueError(
            f"expected exactly one TAMR_vs_MCF7 and one TAMR_vs_FASR row for {BENCHMARK_GENE_SYMBOL}, "
            f"found {len(mcf7_rows)} and {len(fasr_rows)}"
        )
    mcf7_row = mcf7_rows.iloc[0]
    fasr_row = fasr_rows.iloc[0]

    effect_size = float(crispr_row["effect_size"])
    crispr_fdr = float(crispr_row["fdr"])
    mcf7_log2fc = float(mcf7_row["log2fc"])
    mcf7_fdr = float(mcf7_row["fdr"])
    fasr_log2fc = float(fasr_row["log2fc"])
    fasr_fdr = float(fasr_row["fdr"])

    finite_values = {
        "crispr_effect_size": effect_size,
        "crispr_fdr": crispr_fdr,
        "tamr_vs_mcf7_log2fc": mcf7_log2fc,
        "tamr_vs_mcf7_fdr": mcf7_fdr,
        "tamr_vs_fasr_log2fc": fasr_log2fc,
        "tamr_vs_fasr_fdr": fasr_fdr,
    }
    for name, value in finite_values.items():
        if not np.isfinite(value):
            raise ValueError(f"{BENCHMARK_GENE_SYMBOL} benchmark value {name} is not finite: {value!r}")
    for name, value in (("crispr_fdr", crispr_fdr), ("tamr_vs_mcf7_fdr", mcf7_fdr), ("tamr_vs_fasr_fdr", fasr_fdr)):
        if not (0.0 <= value <= 1.0):
            raise ValueError(f"{BENCHMARK_GENE_SYMBOL} benchmark value {name} is outside [0, 1]: {value!r}")

    record = {
        "gene_symbol": BENCHMARK_GENE_SYMBOL,
        "crispr_effect_size": effect_size,
        "crispr_fdr": crispr_fdr,
        "crispr_n_guides": int(crispr_row["n_guides"]),
        "crispr_direction": classify_crispr_direction(effect_size),
        "tamr_vs_mcf7_log2fc": mcf7_log2fc,
        "tamr_vs_mcf7_fdr": mcf7_fdr,
        "tamr_vs_fasr_log2fc": fasr_log2fc,
        "tamr_vs_fasr_fdr": fasr_fdr,
        "benchmark_label": BENCHMARK_LABEL,
    }
    logger.info(
        "build_paics_benchmark: %s effect_size=%.4g fdr=%.4g direction=%s (published benchmark only, "
        "not a Gate-1 hit, not used to tune thresholds)",
        BENCHMARK_GENE_SYMBOL,
        record["crispr_effect_size"],
        record["crispr_fdr"],
        record["crispr_direction"],
    )
    return pd.DataFrame([record])


def build_plot_inputs(evidence_df: pd.DataFrame, shortlist_df: pd.DataFrame) -> dict[str, pd.DataFrame]:
    """Plain, long/wide tidy tables for later plotting -- no figure is
    produced here. The sensitising/tolerance groups are not repeated here
    since the caller already holds those DataFrames directly (see
    ``run_candidate_evidence_summary``'s return dict)."""
    crispr_effects = evidence_df[
        ["gene_symbol", "crispr_effect_size", "crispr_fdr", "crispr_direction", "evidence_class"]
    ].copy()
    crispr_effects = crispr_effects.sort_values("crispr_fdr").reset_index(drop=True)

    contrasts_rows = []
    for record in evidence_df.itertuples(index=False):
        for contrast, log2fc, fdr in (
            ("TAMR_vs_MCF7", record.tamr_vs_mcf7_log2fc, record.tamr_vs_mcf7_fdr),
            ("TAMR_vs_FASR", record.tamr_vs_fasr_log2fc, record.tamr_vs_fasr_fdr),
        ):
            contrasts_rows.append(
                {
                    "gene_symbol": record.gene_symbol,
                    "contrast": contrast,
                    "log2fc": log2fc,
                    "fdr": fdr,
                    "crispr_direction": record.crispr_direction,
                    "evidence_class": record.evidence_class,
                }
            )
    gse118713_contrasts = pd.DataFrame(contrasts_rows)

    shortlist_matrix = evidence_df.merge(
        shortlist_df[["gene_symbol"]], on="gene_symbol", how="inner"
    ).sort_values(["crispr_fdr", "gene_symbol"]).reset_index(drop=True)

    return {
        "crispr_effects": crispr_effects,
        "gse118713_contrasts": gse118713_contrasts,
        "shortlist_matrix": shortlist_matrix,
    }


def write_outputs(
    evidence_df: pd.DataFrame,
    sensitising_df: pd.DataFrame,
    tolerance_df: pd.DataFrame,
    shortlist_df: pd.DataFrame,
    secondary_df: pd.DataFrame,
    paics_df: pd.DataFrame,
    plot_inputs: dict[str, pd.DataFrame],
    cfg: EvidenceSummaryConfig,
) -> None:
    paths = (
        cfg.evidence_summary_tsv,
        cfg.sensitisation_candidates_tsv,
        cfg.tolerance_hits_tsv,
        cfg.shortlist_tsv,
        cfg.secondary_candidates_tsv,
        cfg.paics_benchmark_tsv,
        cfg.plot_crispr_effects_tsv,
        cfg.plot_gse118713_contrasts_tsv,
        cfg.plot_shortlist_matrix_tsv,
    )
    for path in paths:
        path.parent.mkdir(parents=True, exist_ok=True)

    evidence_df.to_csv(cfg.evidence_summary_tsv, sep="\t", index=False)
    sensitising_df.to_csv(cfg.sensitisation_candidates_tsv, sep="\t", index=False)
    tolerance_df.to_csv(cfg.tolerance_hits_tsv, sep="\t", index=False)
    shortlist_df.to_csv(cfg.shortlist_tsv, sep="\t", index=False)
    secondary_df.to_csv(cfg.secondary_candidates_tsv, sep="\t", index=False)
    paics_df.to_csv(cfg.paics_benchmark_tsv, sep="\t", index=False)
    plot_inputs["crispr_effects"].to_csv(cfg.plot_crispr_effects_tsv, sep="\t", index=False)
    plot_inputs["gse118713_contrasts"].to_csv(cfg.plot_gse118713_contrasts_tsv, sep="\t", index=False)
    plot_inputs["shortlist_matrix"].to_csv(cfg.plot_shortlist_matrix_tsv, sep="\t", index=False)
    logger.info(
        "write_outputs: wrote evidence summary (%d), sensitising group (%d), tolerance group (%d), "
        "primary shortlist (%d), secondary candidates (%d), PAICS benchmark (%d row), and plot-input tables",
        len(evidence_df),
        len(sensitising_df),
        len(tolerance_df),
        len(shortlist_df),
        len(secondary_df),
        len(paics_df),
    )


def run_candidate_evidence_summary(config_path: str | Path = "config/config.yaml") -> dict[str, object]:
    """Run the full descriptive evidence summary / direction split /
    primary-secondary evidence hierarchy / shortlist / PAICS benchmark /
    plot-input pipeline and write all outputs."""
    config = _load_config(config_path)
    cfg = EvidenceSummaryConfig.from_config(config)

    master_df = load_master_table(cfg)
    evidence_df = build_evidence_summary(master_df, cfg.rna_significance_fdr)
    sensitising_df, tolerance_df = split_by_direction(evidence_df)
    shortlist_df = select_primary_shortlist(sensitising_df)
    secondary_df = select_secondary_candidates(sensitising_df)
    plot_inputs = build_plot_inputs(evidence_df, shortlist_df)

    labels_df = load_full_crispr_labels(cfg.labels_parquet)
    de_df = load_gse118713_differential_expression(cfg.gse118713_de_tsv_gz)
    paics_df = build_paics_benchmark(labels_df, de_df, master_df)

    write_outputs(evidence_df, sensitising_df, tolerance_df, shortlist_df, secondary_df, paics_df, plot_inputs, cfg)

    return {
        "evidence_summary": evidence_df,
        "sensitising": sensitising_df,
        "tolerance": tolerance_df,
        "shortlist": shortlist_df,
        "secondary": secondary_df,
        "paics_benchmark": paics_df,
        **plot_inputs,
    }


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run_candidate_evidence_summary()
