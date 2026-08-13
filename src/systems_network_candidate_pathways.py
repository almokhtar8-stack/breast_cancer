"""Systems-network phase 7: candidate -> pathway mapping for the four frozen
candidates (USP34, VEZF1, EML5, CITED2). Only pathways with actual evidence
are included -- curated gene-set membership, leading-edge membership in a
dataset that showed at least nominal enrichment, or containing a
high-confidence (STRING score>=700) direct interactor of the candidate in a
dataset that showed at least nominal enrichment for that pathway. A pathway
is never included just because its name sounds biologically relevant.
"""

from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd
import yaml

logger = logging.getLogger(__name__)

CANDIDATES = ["USP34", "VEZF1", "EML5", "CITED2"]
RNA_DATASETS = ["gse118713", "gse240112", "gse111151", "gse245601"]
NOMINAL_ENRICHED_P = 0.05


def _load_config(config_path: str | Path) -> dict:
    with open(config_path) as f:
        return yaml.safe_load(f)


def parse_gmt(path: str | Path) -> dict[str, set[str]]:
    sets: dict[str, set[str]] = {}
    with open(path) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            sets[parts[0]] = set(parts[2:])
    return sets


def load_genesets(cfg: dict) -> dict[str, dict[str, set[str]]]:
    return {c: parse_gmt(cfg["genesets"][f"{c}_gmt"]) for c in ("hallmark", "reactome", "go_bp")}


def load_string_partners(path: str | Path, min_score: float = 0.7) -> dict[str, set[str]]:
    """High-confidence (score>=min_score) direct interactors per candidate,
    symmetrized (A->B implies B is a partner of A and vice versa)."""
    df = pd.read_csv(path, sep="\t")
    df = df.loc[df["score"] >= min_score]
    partners: dict[str, set[str]] = {}
    for _, row in df.iterrows():
        a, b = row["preferredName_A"], row["preferredName_B"]
        partners.setdefault(a, set()).add(b)
        partners.setdefault(b, set()).add(a)
    return partners


def build_candidate_pathway_membership(
    candidates: list[str],
    genesets: dict[str, dict[str, set[str]]],
    gsea_tables: dict[str, pd.DataFrame],
    consensus: pd.DataFrame,
    string_partners: dict[str, set[str]],
) -> pd.DataFrame:
    consensus_idx = consensus.set_index(["collection", "pathway"])["consensus_category"]

    # per (dataset, collection, pathway): leading_edge gene set + nominal enrichment flag
    le_lookup: dict[tuple[str, str, str], tuple[set[str], bool]] = {}
    for dataset, df in gsea_tables.items():
        for _, row in df.iterrows():
            genes = set(row["leading_edge_genes"].split(";")) if pd.notna(row["leading_edge_genes"]) else set()
            le_lookup[(dataset, row["collection"], row["pathway"])] = (genes, bool(row["nom_pvalue"] < NOMINAL_ENRICHED_P))

    rows = []
    for candidate in candidates:
        interactors = string_partners.get(candidate, set())
        for collection, sets in genesets.items():
            for pathway, members in sets.items():
                is_member = candidate in members
                le_datasets = []
                any_enriched = False
                for dataset in RNA_DATASETS:
                    genes, enriched = le_lookup.get((dataset, collection, pathway), (set(), False))
                    if enriched:
                        any_enriched = True
                        # Codex review fix: this gate previously counted a
                        # candidate as "leading-edge" whenever it appeared in
                        # the leading-edge gene list at all, even if that
                        # dataset showed no nominal enrichment for the
                        # pathway (leading_edge_genes is populated by gseapy
                        # for every tested pathway, not only enriched ones).
                        # Must require `enriched` too, per this module's own
                        # documented rule.
                        if candidate in genes:
                            le_datasets.append(dataset)
                interactors_in_pathway = sorted(interactors & members)

                include_via_interactor = bool(interactors_in_pathway) and any_enriched
                if not (is_member or le_datasets or include_via_interactor):
                    continue

                row = {
                    "candidate": candidate,
                    "collection": collection,
                    "pathway": pathway,
                    "candidate_is_member": is_member,
                    "candidate_is_leading_edge_datasets": ",".join(le_datasets) if le_datasets else "",
                    "interactor_genes_in_pathway": ";".join(interactors_in_pathway),
                    "resistance_consensus_class": consensus_idx.get((collection, pathway), "not_in_resistance_consensus_table"),
                }
                for dataset in RNA_DATASETS:
                    df = gsea_tables[dataset]
                    match = df.loc[(df["collection"] == collection) & (df["pathway"] == pathway)]
                    if len(match):
                        row[f"{dataset}_NES"] = float(match["NES"].iloc[0])
                        row[f"{dataset}_FDR"] = float(match["fdr"].iloc[0])
                    else:
                        row[f"{dataset}_NES"] = float("nan")
                        row[f"{dataset}_FDR"] = float("nan")
                rows.append(row)

    out = pd.DataFrame(rows)
    logger.info("build_candidate_pathway_membership: %d rows across %d candidates", len(out), len(candidates))
    for candidate in candidates:
        n = (out["candidate"] == candidate).sum()
        logger.info("  %s: %d associated pathways", candidate, n)
    return out


def run_candidate_pathways(config_path: str | Path = "config/config.yaml") -> pd.DataFrame:
    config = _load_config(config_path)
    cfg = config["systems_network"]
    tables_dir = Path(cfg["output"]["tables_dir"])

    genesets = load_genesets(cfg)
    gsea_tables = {d: pd.read_csv(tables_dir / f"gsea_{d}.tsv", sep="\t") for d in RNA_DATASETS}
    consensus = pd.read_csv(tables_dir / "resistance_pathway_consensus.tsv", sep="\t")
    string_partners = load_string_partners("data/reference/interactions/string_candidate_partners.tsv")

    out = build_candidate_pathway_membership(CANDIDATES, genesets, gsea_tables, consensus, string_partners)
    out.to_csv(tables_dir / "candidate_pathway_membership.tsv", sep="\t", index=False)
    return out


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run_candidate_pathways()
