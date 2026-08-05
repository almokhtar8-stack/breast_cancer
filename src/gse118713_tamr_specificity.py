"""TAMR-specificity table for GSE118713 Phase 2B.

Source: the per-gene, per-contrast limma differential-expression table
written by ``scripts/analysis/gse118713_limma.R``
(config ``gse118713_phase2b.limma.differential_expression_tsv_gz``).

Joins the three preregistered contrasts (TAMR_vs_MCF7, FASR_vs_MCF7,
TAMR_vs_FASR) by gene. Per PREANALYSIS.md's 2026-08-05 Phase 2B
statistical-plan amendment, TAMR specificity is defined only by the direct
TAMR_vs_FASR contrast -- a nonsignificant FASR_vs_MCF7 result is never
used, on its own or in combination, to declare a gene TAMR-specific. This
module reports the three contrasts' effects/FDRs and the TAMR_vs_MCF7 /
TAMR_vs_FASR direction concordance as descriptive columns only; it
computes no specificity label or ranking.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path

import pandas as pd
import yaml

logger = logging.getLogger(__name__)

REQUIRED_CONTRASTS: tuple[str, ...] = ("TAMR_vs_MCF7", "FASR_vs_MCF7", "TAMR_vs_FASR")


def _load_config(config_path: str | Path) -> dict:
    with open(config_path) as f:
        return yaml.safe_load(f)


@dataclass(frozen=True)
class SpecificityConfig:
    differential_expression_tsv_gz: Path
    output_tsv_gz: Path

    @classmethod
    def from_config(cls, config: dict) -> "SpecificityConfig":
        phase2b = config["gse118713_phase2b"]
        return cls(
            differential_expression_tsv_gz=Path(phase2b["limma"]["differential_expression_tsv_gz"]),
            output_tsv_gz=Path(phase2b["specificity"]["output_tsv_gz"]),
        )


def load_de_table(path: str | Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    required_cols = {"gene_id", "gene_symbol", "log2fc", "fdr", "contrast"}
    missing = required_cols - set(df.columns)
    if missing:
        raise ValueError(f"differential-expression table missing required columns: {sorted(missing)}")

    present_contrasts = set(df["contrast"].unique())
    missing_contrasts = set(REQUIRED_CONTRASTS) - present_contrasts
    if missing_contrasts:
        raise ValueError(f"differential-expression table missing required contrasts: {sorted(missing_contrasts)}")
    logger.info("load_de_table: read %d rows across %d contrasts from %s", len(df), len(present_contrasts), path)
    return df


def build_specificity_table(de_df: pd.DataFrame) -> pd.DataFrame:
    """Join the three contrasts by gene_id and report effects/FDRs side by side.

    ``tamr_vs_fasr_log2fc`` is the direct specificity score. The
    ``same_direction_tamr_vs_mcf7_and_tamr_vs_fasr`` column records only
    sign concordance between TAMR_vs_MCF7 and TAMR_vs_FASR -- it is not a
    specificity verdict and does not use FASR_vs_MCF7 significance.
    """
    pieces = {}
    for contrast in REQUIRED_CONTRASTS:
        sub = de_df.loc[de_df["contrast"] == contrast, ["gene_id", "gene_symbol", "log2fc", "fdr"]]
        if sub["gene_id"].duplicated().any():
            raise ValueError(f"duplicate gene_id within contrast {contrast} in differential-expression table")
        pieces[contrast] = sub.set_index("gene_id")

    gene_sets = {c: set(p.index) for c, p in pieces.items()}
    reference = gene_sets[REQUIRED_CONTRASTS[0]]
    for contrast, genes in gene_sets.items():
        if genes != reference:
            raise ValueError(f"contrast {contrast} gene set does not match {REQUIRED_CONTRASTS[0]}")

    base = pieces[REQUIRED_CONTRASTS[0]][["gene_symbol"]].copy()
    for contrast in REQUIRED_CONTRASTS:
        prefix = contrast.lower()
        base[f"{prefix}_log2fc"] = pieces[contrast]["log2fc"]
        base[f"{prefix}_fdr"] = pieces[contrast]["fdr"]

    base["tamr_vs_fasr_log2fc"] = base["tamr_vs_fasr_log2fc"]
    same_sign = (base["tamr_vs_mcf7_log2fc"] > 0) == (base["tamr_vs_fasr_log2fc"] > 0)
    zero_involved = (base["tamr_vs_mcf7_log2fc"] == 0) | (base["tamr_vs_fasr_log2fc"] == 0)
    base["same_direction_tamr_vs_mcf7_and_tamr_vs_fasr"] = same_sign & ~zero_involved

    out = base.reset_index().rename(columns={"index": "gene_id"})
    out = out.sort_values(["tamr_vs_fasr_fdr", "gene_id"], ascending=[True, True]).reset_index(drop=True)
    logger.info("build_specificity_table: joined %d genes across %d contrasts", len(out), len(REQUIRED_CONTRASTS))
    return out


def write_specificity_table(specificity_df: pd.DataFrame, cfg: SpecificityConfig) -> None:
    cfg.output_tsv_gz.parent.mkdir(parents=True, exist_ok=True)
    # mtime=0 keeps the gzip header byte-identical across reruns of identical data.
    specificity_df.to_csv(
        cfg.output_tsv_gz,
        sep="\t",
        index=False,
        compression={"method": "gzip", "mtime": 0},
    )
    logger.info("write_specificity_table: wrote %d genes to %s", len(specificity_df), cfg.output_tsv_gz)


def run_tamr_specificity(config_path: str | Path = "config/config.yaml") -> pd.DataFrame:
    config = _load_config(config_path)
    cfg = SpecificityConfig.from_config(config)
    de_df = load_de_table(cfg.differential_expression_tsv_gz)
    specificity_df = build_specificity_table(de_df)
    write_specificity_table(specificity_df, cfg)
    return specificity_df


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run_tamr_specificity()
