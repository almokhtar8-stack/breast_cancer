"""Gate 1 decision, essentiality-contamination check, and direction-sanity
extraction for the Hany 2023 tamoxifen-response labels.

Source: gene-by-treatment interaction effect table produced by
``src.labels.build_labels`` from Hany et al. 2023, Data S1 (DOI:
10.1126/sciadv.add3685), written to ``data/processed/labels.parquet``.

Reference: Hart et al. 2017 Core Essential Genes 2.0 (CEG2), G3 7:2719-2727,
DOI: 10.1534/g3.117.041277, retrieved via
``scripts/download/download_ceg2.py`` (config ``data.reference.ceg2_gene_list``).

Implements PREANALYSIS.md §4 (Gate 1 decision) exactly as preregistered,
plus two supporting checks: an essentiality-contamination screen (is the
bottom of the effect-size ranking dominated by generically essential genes
rather than tamoxifen-specific hits?) and a direction-sanity extraction for
seven named, non-blind genes. No function in this module inspects,
computes on, logs, or writes any screen-derived value (effect size, p-value,
FDR, rank, or set membership) for the two blind positive-control genes
named in PREANALYSIS.md §5. The test suite contains exactly one reference
to those two gene symbols, and it only checks a static config list against
a hardcoded policy -- it never reads a computed value for them.
"""

from __future__ import annotations

import hashlib
import logging
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from scipy.stats import hypergeom

logger = logging.getLogger(__name__)

REQUIRED_LABEL_COLUMNS: tuple[str, ...] = (
    "gene",
    "n_guides",
    "effect_size",
    "se",
    "p_value",
    "fdr",
)
FINITE_REQUIRED_COLUMNS: tuple[str, ...] = ("effect_size", "se", "p_value", "fdr")

BRANCH_CLASSIFIER = "classifier_permitted"
BRANCH_CONTINUOUS = "continuous_modelling_only"
BRANCH_NONE = "no_predictive_model"


def _load_config(config_path: str | Path) -> dict:
    with open(config_path) as f:
        return yaml.safe_load(f)


def _sha256_file(path: str | Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_and_validate_labels(labels_path: str | Path, n_fitted_genes_expected: int) -> pd.DataFrame:
    """Load the frozen label table and validate its structural invariants.

    Any mismatch (missing columns, null/duplicate gene identifiers,
    non-finite required values, or a row count other than
    ``n_fitted_genes_expected``) is treated as an error, not a silent
    continuation.
    """
    labels_path = Path(labels_path)
    df = pd.read_parquet(labels_path)
    logger.info("load_and_validate_labels: read %d rows from %s", len(df), labels_path)

    missing_cols = [c for c in REQUIRED_LABEL_COLUMNS if c not in df.columns]
    if missing_cols:
        raise ValueError(f"labels table missing required columns: {missing_cols}")

    if df["gene"].isna().any():
        raise ValueError("null gene identifiers in labels table")
    if df["gene"].duplicated().any():
        duplicated_count = int(df["gene"].duplicated().sum())
        raise ValueError(f"{duplicated_count} duplicate gene identifiers in labels table")

    finite = np.isfinite(df[list(FINITE_REQUIRED_COLUMNS)].to_numpy(dtype=float))
    if not finite.all():
        n_non_finite = int((~finite).sum())
        raise ValueError(
            f"{n_non_finite} non-finite values across {FINITE_REQUIRED_COLUMNS} in labels table"
        )

    if len(df) != n_fitted_genes_expected:
        raise ValueError(
            f"expected exactly {n_fitted_genes_expected} fitted genes, got {len(df)}"
        )

    logger.info(
        "load_and_validate_labels: validated %d genes (unique, non-null identifiers; "
        "finite required values)",
        len(df),
    )
    return df


def decide_gate1(
    labels_df: pd.DataFrame,
    fdr_threshold: float,
    classifier_min: int,
    continuous_min: int,
) -> dict[str, object]:
    """Apply the PREANALYSIS.md §4 branches to the calculated FDR<threshold count.

    The number of genes passing the FDR threshold is calculated here from
    ``labels_df``, never hardcoded.
    """
    passing = labels_df["fdr"] < fdr_threshold
    n_passing = int(passing.sum())
    logger.info(
        "decide_gate1: %d of %d fitted genes have fdr < %s",
        n_passing,
        len(labels_df),
        fdr_threshold,
    )

    if n_passing >= classifier_min:
        band = f">= {classifier_min}"
        branch = BRANCH_CLASSIFIER
    elif n_passing >= continuous_min:
        band = f"{continuous_min}-{classifier_min - 1}"
        branch = BRANCH_CONTINUOUS
    else:
        band = f"< {continuous_min}"
        branch = BRANCH_NONE

    logger.info("decide_gate1: branch=%s (band=%s)", branch, band)
    return {
        "total_fitted_genes": len(labels_df),
        "fdr_threshold": fdr_threshold,
        "n_passing": n_passing,
        "threshold_band": band,
        "branch_decision": branch,
    }


def write_gate1_decision(
    decision: dict[str, object],
    labels_path: str | Path,
    output_path: str | Path,
) -> pd.DataFrame:
    row = dict(decision)
    row["labels_input_path"] = str(labels_path)
    row["labels_file_sha256"] = _sha256_file(labels_path)

    out = pd.DataFrame([row])
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(output_path, sep="\t", index=False)
    logger.info("write_gate1_decision: wrote decision to %s", output_path)
    return out


def load_ceg2_reference(ceg2_path: str | Path) -> set[str]:
    """Load the frozen CEG2 gene-symbol set (see scripts/download/download_ceg2.py).

    Validates before collapsing to a set: a null/blank symbol would
    silently vanish, and a duplicate symbol would silently collapse,
    without either being counted.
    """
    ceg2_path = Path(ceg2_path)
    ceg2_df = pd.read_csv(ceg2_path, sep="\t")
    if "gene_symbol" not in ceg2_df.columns:
        raise ValueError(f"CEG2 reference file missing 'gene_symbol' column: {ceg2_path}")

    raw_symbols = ceg2_df["gene_symbol"]
    n_rows = len(raw_symbols)
    blank = raw_symbols.isna() | (raw_symbols.astype(str).str.strip() == "")
    if blank.any():
        raise ValueError(f"{int(blank.sum())} null/blank gene symbols in CEG2 reference {ceg2_path}")

    if raw_symbols.duplicated().any():
        n_duplicates = int(raw_symbols.duplicated().sum())
        raise ValueError(
            f"{n_duplicates} duplicate gene symbols in CEG2 reference {ceg2_path}"
        )

    symbols = set(raw_symbols)
    logger.info(
        "load_ceg2_reference: read %d rows, %d unique CEG2 gene symbols from %s",
        n_rows,
        len(symbols),
        ceg2_path,
    )
    return symbols


def _rank_by_effect_size(labels_df: pd.DataFrame) -> pd.DataFrame:
    """Sort ascending by effect_size, breaking ties by gene symbol ascending.

    The secondary key makes the ranking (and therefore the bottom-N tested
    set drawn from it) deterministic when two or more genes share an
    effect_size value.
    """
    return labels_df.sort_values(["effect_size", "gene"], ascending=[True, True])


def essentiality_contamination_check(
    labels_df: pd.DataFrame,
    ceg2_symbols: set[str],
    tested_set_size: int,
) -> dict[str, object]:
    """Rank fitted genes by effect_size ascending and check CEG2 enrichment
    among the most negative ``tested_set_size`` genes.

    The statistical universe is every gene in ``labels_df`` (the fitted-gene
    table), not the full CEG2 set or any external gene catalogue. CEG2 is
    intersected with this universe before computing chance-expected overlap
    and the hypergeometric enrichment p-value. Individual gene identities
    (tested-set members or CEG2 symbols) are never logged or returned here,
    only counts.
    """
    universe_genes = set(labels_df["gene"])
    universe_size = len(universe_genes)

    ceg2_in_universe = ceg2_symbols & universe_genes
    ceg2_absent = ceg2_symbols - universe_genes
    n_ceg2_in_universe = len(ceg2_in_universe)
    n_ceg2_absent = len(ceg2_absent)
    logger.info(
        "essentiality_contamination_check: fitted-gene universe=%d; CEG2 total=%d; "
        "CEG2 in universe=%d; CEG2 absent from universe=%d",
        universe_size,
        len(ceg2_symbols),
        n_ceg2_in_universe,
        n_ceg2_absent,
    )

    ranked = _rank_by_effect_size(labels_df)
    if len(ranked) < tested_set_size:
        raise ValueError(
            f"fitted-gene universe ({len(ranked)}) smaller than requested "
            f"tested_set_size ({tested_set_size})"
        )
    tested_set = ranked.head(tested_set_size)
    tested_set_genes = set(tested_set["gene"])
    if len(tested_set_genes) != tested_set_size:
        raise ValueError(
            f"tested set has {len(tested_set_genes)} unique genes, expected {tested_set_size} "
            "(duplicate gene identifiers should have been rejected at load time)"
        )

    observed_overlap = len(tested_set_genes & ceg2_in_universe)
    expected_overlap = tested_set_size * (n_ceg2_in_universe / universe_size)
    enrichment = observed_overlap / expected_overlap if expected_overlap > 0 else float("nan")

    # One-sided: P(X >= observed_overlap) under the hypergeometric null.
    p_value = float(
        hypergeom.sf(observed_overlap - 1, universe_size, n_ceg2_in_universe, tested_set_size)
    )

    logger.info(
        "essentiality_contamination_check: tested_set_size=%d observed_overlap=%d "
        "expected_overlap=%.4f enrichment=%.4f p_value=%.6g",
        tested_set_size,
        observed_overlap,
        expected_overlap,
        enrichment,
        p_value,
    )

    return {
        "fitted_gene_universe_size": universe_size,
        "ceg2_genes_in_universe": n_ceg2_in_universe,
        "tested_set_size": tested_set_size,
        "observed_ceg2_overlap": observed_overlap,
        "expected_ceg2_overlap": expected_overlap,
        "enrichment_observed_over_expected": enrichment,
        "hypergeometric_p_value_one_sided": p_value,
        "ceg2_symbols_absent_from_universe": n_ceg2_absent,
    }


def write_essentiality_summary(summary: dict[str, object], output_path: str | Path) -> pd.DataFrame:
    out = pd.DataFrame([summary])
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(output_path, sep="\t", index=False)
    logger.info("write_essentiality_summary: wrote aggregate summary to %s", output_path)
    return out


def _observed_direction(effect_size: float) -> str:
    if effect_size < 0:
        return "depleted_under_4OHT"
    if effect_size > 0:
        return "enriched_under_4OHT"
    return "no_change"


def _sign_interpretation(effect_size: float) -> str:
    if effect_size < 0:
        return "sensitizes (knockout needed for tamoxifen tolerance)"
    if effect_size > 0:
        return "growth_advantage (knockout confers advantage under tamoxifen)"
    return "no_directional_call"


def direction_sanity_extraction(labels_df: pd.DataFrame, genes: list[str]) -> pd.DataFrame:
    """Extract the named genes, in the given order, with observed values only.

    No literature-expectation labelling is applied; only the values as
    computed and the documented sign convention (PREANALYSIS.md §2) are
    reported. Fails if any requested gene is absent from the labels table.
    """
    indexed = labels_df.set_index("gene")
    missing = [g for g in genes if g not in indexed.index]
    if missing:
        raise ValueError(f"requested direction-sanity genes absent from labels table: {missing}")

    rows = []
    for gene in genes:
        row = indexed.loc[gene]
        effect_size = float(row["effect_size"])
        se = float(row["se"])
        rows.append(
            {
                "gene": gene,
                "n_guides": int(row["n_guides"]),
                "effect_size": effect_size,
                "se": se,
                "t_moderated": effect_size / se,
                "p_value": float(row["p_value"]),
                "fdr": float(row["fdr"]),
                "observed_direction": _observed_direction(effect_size),
                "sign_interpretation": _sign_interpretation(effect_size),
            }
        )
    logger.info("direction_sanity_extraction: extracted %d requested genes", len(rows))
    return pd.DataFrame(rows)


def write_direction_sanity(direction_df: pd.DataFrame, output_path: str | Path) -> pd.DataFrame:
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    direction_df.to_csv(output_path, sep="\t", index=False)
    logger.info("write_direction_sanity: wrote %d rows to %s", len(direction_df), output_path)
    return direction_df


def run_gate1_checks(config_path: str | Path = "config/config.yaml") -> dict[str, pd.DataFrame]:
    """Run C1 (Gate 1 decision), C2 (essentiality-contamination), and C3
    (direction-sanity) and write their output tables."""
    config = _load_config(config_path)
    gate1_config = config["gate1"]

    labels_path = gate1_config["labels_path"]
    labels_df = load_and_validate_labels(
        labels_path, n_fitted_genes_expected=gate1_config["n_fitted_genes_expected"]
    )

    decision = decide_gate1(
        labels_df,
        fdr_threshold=gate1_config["fdr_threshold"],
        classifier_min=gate1_config["branch_thresholds"]["classifier_min"],
        continuous_min=gate1_config["branch_thresholds"]["continuous_min"],
    )
    decision_df = write_gate1_decision(decision, labels_path, gate1_config["decision_output"])

    ceg2_symbols = load_ceg2_reference(config["data"]["reference"]["ceg2_gene_list"])
    essentiality = essentiality_contamination_check(
        labels_df,
        ceg2_symbols,
        tested_set_size=gate1_config["essentiality"]["tested_set_size"],
    )
    essentiality_df = write_essentiality_summary(
        essentiality, gate1_config["essentiality"]["output"]
    )

    direction_df = direction_sanity_extraction(labels_df, gate1_config["direction_sanity"]["genes"])
    direction_df = write_direction_sanity(direction_df, gate1_config["direction_sanity"]["output"])

    return {
        "decision": decision_df,
        "essentiality": essentiality_df,
        "direction_sanity": direction_df,
    }


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run_gate1_checks()
