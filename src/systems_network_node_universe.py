"""Systems-network phase 14: focused network node universe.

Implements the inclusion rule declared in docs/SYSTEMS_NETWORK_NODE_RULE.md
(categories A-E) BEFORE this module was written to inspect which genes
land in the network -- the rule (fixed top-40 cuts, FDR<0.05 CRISPR,
score>=0.7 STRING) is applied mechanically here, not tuned to hit a target
node count after the fact.
"""

from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd
import yaml

logger = logging.getLogger(__name__)

CANDIDATES = ["USP34", "VEZF1", "EML5", "CITED2"]
TOP_N_LEADING_EDGE = 40
TOP_N_MULTIMODAL = 40


def _load_config(config_path: str | Path) -> dict:
    with open(config_path) as f:
        return yaml.safe_load(f)


def category_A_candidates() -> pd.DataFrame:
    return pd.DataFrame({"gene": CANDIDATES, "node_category": "candidate"})


def category_B_recurrent_leading_edge(tables_dir: Path) -> pd.DataFrame:
    df = pd.read_csv(tables_dir / "recurrent_leading_edge_genes.tsv", sep="\t")
    # Codex review fix (second pass): rank by `max_same_pathway_dataset_count`
    # (the stricter "reproducibly drives the SAME pathway's enrichment
    # signal across datasets" statistic), not the looser
    # `leading_edge_dataset_count` (reachable via a different pathway per
    # dataset) -- the looser statistic previously let genes like DKK1/EGFR/
    # RET/SRC/VEGFA (same-pathway recurrence <=2) outrank genes with
    # genuine same-pathway 3-dataset recurrence. Deterministic tie-breaker:
    # leading_edge_dataset_count, then pathway_count, then gene symbol.
    df = df.sort_values(
        ["max_same_pathway_dataset_count", "leading_edge_dataset_count", "pathway_count", "gene"],
        ascending=[False, False, False, True],
        kind="mergesort",
    ).head(TOP_N_LEADING_EDGE)
    return pd.DataFrame({"gene": df["gene"].values, "node_category": "resistance_leading_edge"})


def category_C_crispr_hits() -> pd.DataFrame:
    crispr = pd.read_parquet("data/processed/labels.parquet")
    hits = crispr.loc[crispr["fdr"] < 0.05].copy()
    hits["node_category"] = hits["effect_size"].apply(lambda e: "crispr_sensitiser" if e < 0 else "crispr_tolerance")
    return hits[["gene", "node_category"]]


def category_D_multimodal_pathway_genes(tables_dir: Path) -> pd.DataFrame:
    mm = pd.read_csv(tables_dir / "multimodal_pathway_convergence.tsv", sep="\t")
    mm_keys = set(zip(mm.loc[mm["convergence_category"] == "MULTIMODAL_PATHWAY", "collection"], mm.loc[mm["convergence_category"] == "MULTIMODAL_PATHWAY", "pathway"]))

    records = []
    for dataset in ["gse118713", "gse240112", "gse111151"]:
        df = pd.read_csv(tables_dir / f"gsea_{dataset}.tsv", sep="\t")
        df = df.loc[df.apply(lambda r: (r["collection"], r["pathway"]) in mm_keys, axis=1)]
        # Codex review fix: only count a gene as leading-edge in THIS
        # dataset if this dataset's own nominal p-value for the pathway is
        # significant (same fix as Phase 9/10 -- MULTIMODAL_PATHWAY status
        # is a property of the pathway across datasets, not a guarantee
        # every individual dataset's own result for it was significant).
        df = df.loc[df["nom_pvalue"] < 0.05]
        for _, row in df.iterrows():
            if pd.isna(row["leading_edge_genes"]):
                continue
            for gene in row["leading_edge_genes"].split(";"):
                records.append(gene)
    counts = pd.Series(records).value_counts()
    # deterministic tie-breaker (Codex review note): value_counts() ties are
    # broken alphabetically before truncating to the top-N, not by
    # insertion/hash order
    counts = counts.sort_index().sort_values(ascending=False, kind="mergesort").head(TOP_N_MULTIMODAL)
    return pd.DataFrame({"gene": counts.index, "node_category": "multimodal_pathway_driver"})


def category_E_candidate_string_partners(path: str | Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    df = df.loc[df["score"] >= 0.7]
    partners = sorted(set(df["preferredName_A"]) | set(df["preferredName_B"]) - set(CANDIDATES))
    partners = [g for g in partners if g not in CANDIDATES]
    return pd.DataFrame({"gene": partners, "node_category": "candidate_string_partner"})


def build_node_universe(tables_dir: Path, string_candidate_partners_path: str | Path) -> pd.DataFrame:
    parts = [
        category_A_candidates(),
        category_B_recurrent_leading_edge(tables_dir),
        category_C_crispr_hits(),
        category_D_multimodal_pathway_genes(tables_dir),
        category_E_candidate_string_partners(string_candidate_partners_path),
    ]
    combined = pd.concat(parts, ignore_index=True)
    grouped = combined.groupby("gene")["node_category"].apply(lambda s: ",".join(sorted(set(s)))).reset_index()
    logger.info("build_node_universe: %d unique genes across %d category rows", len(grouped), len(combined))
    return grouped


def run_node_universe(config_path: str | Path = "config/config.yaml") -> pd.DataFrame:
    config = _load_config(config_path)
    cfg = config["systems_network"]
    tables_dir = Path(cfg["output"]["tables_dir"])
    networks_dir = Path(cfg["output"]["networks_dir"])
    networks_dir.mkdir(parents=True, exist_ok=True)

    out = build_node_universe(tables_dir, "data/reference/interactions/string_candidate_partners.tsv")
    out.to_csv(networks_dir / "node_universe.tsv", sep="\t", index=False)
    return out


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run_node_universe()
