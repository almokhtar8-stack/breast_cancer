"""Ad-hoc focused evidence audit for the four USP34 shortest-path bridge
genes (USP9X, RPS27A, UBC, UBB), identified in
results/tables/systems_network/USP34_shortest_paths.tsv.

Reads ONLY already-frozen systems-network outputs -- does not rerun any
upstream ranking/GSEA/consensus phase and does not modify
results/networks/systems_network/ or any evidence-freeze output:
  - results/networks/systems_network/cytoscape/network_nodes.tsv
    (frozen five-layer per-gene evidence display: CRISPR + the three
    resistance datasets + GSE245601 Track A acute)
  - results/tables/systems_network/{gse118713,gse240112,gse111151,
    gse245601}_ranked_genes.tsv (frozen Phase 2 per-dataset p-values, the
    same rows used to build network_nodes.tsv; GSE118713's raw DE file has
    three contrasts -- TAMR_vs_MCF7, FASR_vs_MCF7, TAMR_vs_FASR -- and the
    frozen ranked-genes table already selects TAMR_vs_MCF7 only, the
    resistance contrast used throughout this project; this module reuses
    that frozen selection rather than re-deriving it)
  - results/networks/systems_network/edges.tsv
    (source of the "bridges USP34 to X" network facts)

CRISPR direction convention (frozen from data/processed/labels.parquet,
the hard-rule CRISPR label source): negative effect_size = sensitising_KO,
positive = tolerance_associated_KO. This is a description of screen sign,
never treated as expression direction and never treated as causal proof.

GSE245601 Track A is the acute (12h tamoxifen) pseudobulk result -- it is
reported here as a fifth evidence layer for completeness only and is never
folded into the "resistance dataset" significance/direction-consistency
counts below, per the acute/resistance separation rule used throughout the
systems-network phase.

Classification (conservative, applied uniformly, not tuned per gene):
  A. DATA_SUPPORTED_BRIDGE   -- FDR<0.05 in CRISPR and/or >=1 resistance
                                 dataset for this gene.
  B. PARTIAL_SUPPORT         -- no FDR<0.05 hit anywhere, but >=1 nominal
                                 (p<0.05) hit in CRISPR or a resistance
                                 dataset (GSE245601-only nominal hits do not
                                 qualify, since that is acute, not
                                 resistance, evidence).
  C. NETWORK_ONLY_GENERIC_BRIDGE -- no FDR<0.05 and no qualifying nominal
                                 hit in CRISPR or any resistance dataset.
No tier is described as mechanistic proof; classification reflects only
whether OUR OWN data (CRISPR + resistance transcriptomics) independently
supports the gene, not whether the STRING edge itself is trustworthy.
"""

from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd

logger = logging.getLogger(__name__)

BRIDGE_GENES = ["USP9X", "RPS27A", "UBC", "UBB"]
RESISTANCE_DATASETS = ["gse118713", "gse240112", "gse111151"]
NOMINAL_P = 0.05
FDR_SIG = 0.05

NODES_PATH = Path("results/networks/systems_network/cytoscape/network_nodes.tsv")
EDGES_PATH = Path("results/networks/systems_network/edges.tsv")
RANKED_PATHS = {
    "gse118713": Path("results/tables/systems_network/gse118713_ranked_genes.tsv"),
    "gse240112": Path("results/tables/systems_network/gse240112_ranked_genes.tsv"),
    "gse111151": Path("results/tables/systems_network/gse111151_ranked_genes.tsv"),
    "gse245601": Path("results/tables/systems_network/gse245601_ranked_genes.tsv"),
}
SHORTEST_PATHS_PATH = Path("results/tables/systems_network/USP34_shortest_paths.tsv")

OUT_TABLE = Path("results/tables/systems_network/USP34_bridge_gene_evidence.tsv")
OUT_REPORT = Path("results/reports/systems_network/USP34_bridge_gene_evidence.md")


def _direction(log2fc: float) -> str:
    return "up" if log2fc > 0 else ("down" if log2fc < 0 else "flat")


def _bridges_to(shortest_paths: pd.DataFrame) -> dict[str, str]:
    ok = shortest_paths.loc[shortest_paths["path"] != "NO_PATH_IN_NETWORK"]
    result: dict[str, set[str]] = {g: set() for g in BRIDGE_GENES}
    for _, row in ok.iterrows():
        path_genes = row["path"].split(" -> ")
        if len(path_genes) != 3:
            continue
        _, mid, dest = path_genes
        if mid in result:
            result[mid].add(dest)
    return {g: ",".join(sorted(dests)) for g, dests in result.items()}


def build_bridge_evidence(
    nodes: pd.DataFrame,
    ranked: dict[str, pd.DataFrame],
    shortest_paths: pd.DataFrame,
) -> pd.DataFrame:
    nodes_idx = nodes.set_index("gene")
    bridges_to = _bridges_to(shortest_paths)

    crispr_labels = pd.read_parquet("data/processed/labels.parquet").set_index("gene")

    resistance_flags: dict[str, list[bool]] = {g: [] for g in BRIDGE_GENES}
    resistance_dirs: dict[str, list[str]] = {g: [] for g in BRIDGE_GENES}

    out_rows = []
    for gene in BRIDGE_GENES:
        n = nodes_idx.loc[gene]
        crispr_row = crispr_labels.loc[gene]
        crispr_p = float(crispr_row["p_value"])
        crispr_fdr = float(crispr_row["fdr"])
        crispr_effect = float(crispr_row["effect_size"])
        crispr_sig_fdr = crispr_fdr < FDR_SIG
        crispr_sig_nom = crispr_p < NOMINAL_P

        rec = {
            "gene": gene,
            "bridges_usp34_to": bridges_to.get(gene, ""),
            "crispr_effect": crispr_effect,
            "crispr_p": crispr_p,
            "crispr_fdr": crispr_fdr,
            "crispr_direction": n["crispr_direction"],
            "crispr_significant_fdr05": crispr_sig_fdr,
            "crispr_nominal_p05": crispr_sig_nom,
        }

        any_fdr_sig = crispr_sig_fdr
        any_nominal_resistance_sig = False

        for dataset in RESISTANCE_DATASETS:
            r = ranked[dataset].set_index("gene").loc[gene]
            log2fc = float(r["log2fc"])
            p = float(r["p_value"])
            fdr = float(r["fdr"])
            sig_fdr = fdr < FDR_SIG
            sig_nom = p < NOMINAL_P
            direction = _direction(log2fc)
            rec[f"{dataset}_log2fc"] = log2fc
            rec[f"{dataset}_p"] = p
            rec[f"{dataset}_fdr"] = fdr
            rec[f"{dataset}_direction"] = direction
            rec[f"{dataset}_significant_fdr05"] = sig_fdr
            rec[f"{dataset}_nominal_p05"] = sig_nom
            resistance_flags[gene].append(sig_fdr)
            resistance_dirs[gene].append(direction)
            any_fdr_sig = any_fdr_sig or sig_fdr
            any_nominal_resistance_sig = any_nominal_resistance_sig or sig_nom

        acute = ranked["gse245601"].set_index("gene").loc[gene]
        acute_log2fc = float(acute["log2fc"])
        acute_p = float(acute["p_value"])
        acute_fdr = float(acute["fdr"])
        rec["gse245601_track_a_log2fc"] = acute_log2fc
        rec["gse245601_track_a_p"] = acute_p
        rec["gse245601_track_a_fdr"] = acute_fdr
        rec["gse245601_track_a_direction"] = _direction(acute_log2fc)
        rec["gse245601_track_a_significant_fdr05"] = acute_fdr < FDR_SIG
        rec["gse245601_track_a_nominal_p05"] = acute_p < NOMINAL_P
        rec["gse245601_acute_note"] = "ACUTE 12h response -- not resistance evidence, not counted in classification"

        n_up = resistance_dirs[gene].count("up")
        n_down = resistance_dirs[gene].count("down")
        rec["resistance_direction_consistency"] = "consistent" if (n_up == 0 or n_down == 0) else "mixed"
        rec["resistance_pattern_frozen"] = n["resistance_pattern"]

        if any_fdr_sig:
            classification = "A_DATA_SUPPORTED_BRIDGE"
        elif crispr_sig_nom or any_nominal_resistance_sig:
            classification = "B_PARTIAL_SUPPORT"
        else:
            classification = "C_NETWORK_ONLY_GENERIC_BRIDGE"
        rec["classification"] = classification

        out_rows.append(rec)

    return pd.DataFrame(out_rows)


def run(
    nodes_path: Path = NODES_PATH,
    ranked_paths: dict[str, Path] = None,
    shortest_paths_path: Path = SHORTEST_PATHS_PATH,
    out_table: Path = OUT_TABLE,
) -> pd.DataFrame:
    ranked_paths = ranked_paths or RANKED_PATHS
    nodes = pd.read_csv(nodes_path, sep="\t")
    ranked = {k: pd.read_csv(v, sep="\t") for k, v in ranked_paths.items()}
    shortest_paths = pd.read_csv(shortest_paths_path, sep="\t")

    out = build_bridge_evidence(nodes, ranked, shortest_paths)
    out_table.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(out_table, sep="\t", index=False)
    logger.info(
        "wrote %s (%d genes; classification counts: %s)",
        out_table,
        len(out),
        out["classification"].value_counts().to_dict(),
    )
    return out


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run()
