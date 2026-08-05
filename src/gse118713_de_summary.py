"""Per-contrast differential-expression summary for GSE118713 Phase 2B.

Source: the per-gene, per-contrast limma differential-expression table
written by ``scripts/analysis/gse118713_limma.R``
(config ``gse118713_phase2b.limma.differential_expression_tsv_gz``).

Every number here is aggregated directly from that table -- genes tested,
genes at FDR < 0.05, counts of significant positive/negative effects,
median absolute effect, and minimum FDR, per contrast. Nothing is manually
typed.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

logger = logging.getLogger(__name__)

FDR_THRESHOLD = 0.05


def _load_config(config_path: str | Path) -> dict:
    with open(config_path) as f:
        return yaml.safe_load(f)


@dataclass(frozen=True)
class SummaryConfig:
    differential_expression_tsv_gz: Path
    output_tsv: Path
    contrasts: tuple[str, ...]

    @classmethod
    def from_config(cls, config: dict) -> "SummaryConfig":
        phase2b = config["gse118713_phase2b"]
        return cls(
            differential_expression_tsv_gz=Path(phase2b["limma"]["differential_expression_tsv_gz"]),
            output_tsv=Path(phase2b["summary"]["output_tsv"]),
            contrasts=tuple(phase2b["limma"]["contrasts"]),
        )


def build_de_summary(de_df: pd.DataFrame, contrasts: tuple[str, ...], fdr_threshold: float = FDR_THRESHOLD) -> pd.DataFrame:
    """Aggregate genes-tested, significance counts, and effect summaries per contrast."""
    present = set(de_df["contrast"].unique())
    missing = set(contrasts) - present
    if missing:
        raise ValueError(f"differential-expression table missing contrasts: {sorted(missing)}")

    rows = []
    for contrast in contrasts:
        sub = de_df.loc[de_df["contrast"] == contrast]
        significant = sub.loc[sub["fdr"] < fdr_threshold]
        rows.append(
            {
                "contrast": contrast,
                "genes_tested": int(len(sub)),
                "n_significant_fdr_lt_0_05": int(len(significant)),
                "n_significant_positive": int((significant["log2fc"] > 0).sum()),
                "n_significant_negative": int((significant["log2fc"] < 0).sum()),
                "median_abs_log2fc": float(np.median(sub["log2fc"].abs())),
                "min_fdr": float(sub["fdr"].min()),
            }
        )
    logger.info("build_de_summary: summarized %d contrasts", len(rows))
    return pd.DataFrame(rows)


def write_de_summary(summary_df: pd.DataFrame, cfg: SummaryConfig) -> None:
    cfg.output_tsv.parent.mkdir(parents=True, exist_ok=True)
    summary_df.to_csv(cfg.output_tsv, sep="\t", index=False)
    logger.info("write_de_summary: wrote %s", cfg.output_tsv)


def run_de_summary(config_path: str | Path = "config/config.yaml") -> pd.DataFrame:
    config = _load_config(config_path)
    cfg = SummaryConfig.from_config(config)
    de_df = pd.read_csv(cfg.differential_expression_tsv_gz, sep="\t")
    summary_df = build_de_summary(de_df, cfg.contrasts)
    write_de_summary(summary_df, cfg)
    return summary_df


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run_de_summary()
