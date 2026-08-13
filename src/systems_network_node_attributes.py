"""Systems-network phase 15: node attribute table for the focused network.

node_class vocabulary follows the task spec (candidate, leading_edge,
crispr_sensitiser, resistance_gene, candidate_partner, multiple) with one
disclosed extension: crispr_tolerance_associated, for CRISPR hits that pass
FDR<0.05 but in the tolerance (not sensitising) direction -- the task's
five-item vocabulary has no exact slot for these and folding them into
"resistance_gene" would blur a genuinely different functional class.
"""

from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd
import yaml

logger = logging.getLogger(__name__)

CANDIDATES = ["USP34", "VEZF1", "EML5", "CITED2"]
RESISTANCE_DATASETS = ["gse118713", "gse240112", "gse111151"]


def _load_config(config_path: str | Path) -> dict:
    with open(config_path) as f:
        return yaml.safe_load(f)


def _node_class(category_str: str) -> str:
    cats = set(category_str.split(","))
    mapping = {
        "candidate": "candidate",
        "resistance_leading_edge": "leading_edge",
        "crispr_sensitiser": "crispr_sensitiser",
        "crispr_tolerance": "crispr_tolerance_associated",
        "multimodal_pathway_driver": "resistance_gene",
        "candidate_string_partner": "candidate_partner",
    }
    mapped = {mapping[c] for c in cats if c in mapping}
    if len(mapped) > 1:
        return "multiple"
    return next(iter(mapped)) if mapped else "unknown"


def _direction(fdr: float, log2fc: float, sig_only: bool = True) -> str:
    if pd.isna(fdr) or pd.isna(log2fc):
        return "not_tested"
    if sig_only and fdr >= 0.05:
        return "ns"
    return "up" if log2fc > 0 else "down"


def build_node_attributes(
    node_universe: pd.DataFrame,
    ranked_tables: dict[str, pd.DataFrame],
    crispr_ranked: pd.DataFrame,
    consensus: pd.DataFrame,
    gsea_tables: dict[str, pd.DataFrame],
    edges: pd.DataFrame,
) -> pd.DataFrame:
    strong_keys = set(zip(consensus.loc[consensus["consensus_category"] == "STRONG_CONSENSUS", "collection"], consensus.loc[consensus["consensus_category"] == "STRONG_CONSENSUS", "pathway"]))

    # gene -> set of (collection, pathway) it's leading-edge in among STRONG_CONSENSUS pathways, any resistance dataset
    gene_pathways: dict[str, set[tuple[str, str]]] = {}
    for dataset in RESISTANCE_DATASETS:
        df = gsea_tables[dataset]
        df = df.loc[df.apply(lambda r: (r["collection"], r["pathway"]) in strong_keys, axis=1)]
        # Codex review fix: same per-dataset-own-significance gate as
        # Phase 9/10/14 -- otherwise `number_consensus_pathways` and
        # `pathways_supported` could be inflated by a dataset in which this
        # specific pathway was not itself significant.
        df = df.loc[df["nom_pvalue"] < 0.05]
        for _, row in df.iterrows():
            if pd.isna(row["leading_edge_genes"]):
                continue
            for gene in row["leading_edge_genes"].split(";"):
                gene_pathways.setdefault(gene, set()).add((row["collection"], row["pathway"]))

    connection_counts: dict[str, int] = {}
    for gene in node_universe["gene"]:
        connected_candidates = set()
        rows = edges.loc[(edges["source_gene"] == gene) | (edges["target_gene"] == gene)]
        for _, r in rows.iterrows():
            other = r["target_gene"] if r["source_gene"] == gene else r["source_gene"]
            if other in CANDIDATES:
                connected_candidates.add(other)
        connection_counts[gene] = len(connected_candidates)

    rc_idx = {d: ranked_tables[d].set_index("gene") for d in RESISTANCE_DATASETS}
    ac_idx = ranked_tables["gse245601"].set_index("gene")
    cr_idx = crispr_ranked.set_index("gene")

    rows = []
    for _, nrow in node_universe.iterrows():
        gene = nrow["gene"]
        node_class = _node_class(nrow["node_category"])

        crispr_row = cr_idx.loc[gene] if gene in cr_idx.index else None
        crispr_effect = float(crispr_row["log2fc"]) if crispr_row is not None else float("nan")
        crispr_fdr = float(crispr_row["fdr"]) if crispr_row is not None else float("nan")
        if pd.isna(crispr_fdr):
            crispr_direction = "not_tested"
        elif crispr_fdr < 0.05:
            crispr_direction = "sensitising_KO" if crispr_effect < 0 else "tolerance_associated_KO"
        else:
            crispr_direction = "nonsignificant_sensitising_direction" if crispr_effect < 0 else "nonsignificant_tolerance_direction"

        ds_vals = {}
        directions = []
        for d in RESISTANCE_DATASETS:
            r = rc_idx[d].loc[gene] if gene in rc_idx[d].index else None
            log2fc = float(r["log2fc"]) if r is not None else float("nan")
            fdr = float(r["fdr"]) if r is not None else float("nan")
            ds_vals[f"{d}_log2fc"] = log2fc
            ds_vals[f"{d}_fdr"] = fdr
            directions.append(_direction(fdr, log2fc))

        ac_row = ac_idx.loc[gene] if gene in ac_idx.index else None
        acute_log2fc = float(ac_row["log2fc"]) if ac_row is not None else float("nan")
        acute_fdr = float(ac_row["fdr"]) if ac_row is not None else float("nan")

        sig_dirs = [d for d in directions if d in ("up", "down")]
        if not sig_dirs:
            resistance_pattern = "no_significant_resistance_dataset"
        elif len(set(sig_dirs)) == 1:
            resistance_pattern = f"consistent_{sig_dirs[0]}_({len(sig_dirs)}/3_significant)"
        else:
            resistance_pattern = f"mixed_({directions.count('up')}up/{directions.count('down')}down_of_3)"

        gene_pw = gene_pathways.get(gene, set())
        human_tumor_support = bool(
            (gene in rc_idx["gse240112"].index and rc_idx["gse240112"].loc[gene, "fdr"] < 0.05)
            or (gene in ac_idx.index and ac_idx.loc[gene, "fdr"] < 0.05)
        )

        row = {
            "gene": gene,
            "node_class": node_class,
            "node_category_raw": nrow["node_category"],
            "crispr_effect": crispr_effect,
            "crispr_fdr": crispr_fdr,
            "crispr_direction": crispr_direction,
            **ds_vals,
            "gse245601_acute_log2fc": acute_log2fc,
            "gse245601_acute_fdr": acute_fdr,
            "resistance_pattern": resistance_pattern,
            "pathways_supported": ";".join(f"{c}:{p}" for c, p in sorted(gene_pw)[:5]),
            "number_consensus_pathways": len(gene_pw),
            "human_tumor_support": human_tumor_support,
            "candidate_connection_count": connection_counts.get(gene, 0),
        }
        rows.append(row)

    out = pd.DataFrame(rows)
    logger.info("build_node_attributes: %d nodes, node_class counts=%s", len(out), out["node_class"].value_counts().to_dict())
    return out


def run_node_attributes(config_path: str | Path = "config/config.yaml") -> pd.DataFrame:
    config = _load_config(config_path)
    cfg = config["systems_network"]
    tables_dir = Path(cfg["output"]["tables_dir"])
    networks_dir = Path(cfg["output"]["networks_dir"])

    node_universe = pd.read_csv(networks_dir / "node_universe.tsv", sep="\t")
    edges = pd.read_csv(networks_dir / "edges.tsv", sep="\t")
    ranked_tables = {
        "gse118713": pd.read_csv(tables_dir / "gse118713_ranked_genes.tsv", sep="\t"),
        "gse240112": pd.read_csv(tables_dir / "gse240112_ranked_genes.tsv", sep="\t"),
        "gse111151": pd.read_csv(tables_dir / "gse111151_ranked_genes.tsv", sep="\t"),
        "gse245601": pd.read_csv(tables_dir / "gse245601_ranked_genes.tsv", sep="\t"),
    }
    crispr_ranked = pd.read_csv(tables_dir / "crispr_ranked_genes.tsv", sep="\t")
    consensus = pd.read_csv(tables_dir / "resistance_pathway_consensus.tsv", sep="\t")
    gsea_tables = {d: pd.read_csv(tables_dir / f"gsea_{d}.tsv", sep="\t") for d in RESISTANCE_DATASETS}

    out = build_node_attributes(node_universe, ranked_tables, crispr_ranked, consensus, gsea_tables, edges)
    out.to_csv(networks_dir / "nodes.tsv", sep="\t", index=False)
    return out


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run_node_attributes()
