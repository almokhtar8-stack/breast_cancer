"""Systems-network phases 13, 15, 16: curated interaction edges and node/edge
attribute tables for the focused node universe (docs/SYSTEMS_NETWORK_NODE_RULE.md).

Edge sources (each edge carries its own source/evidence, never merged into
an unlabeled "interaction"):
  - STRING (data/reference/interactions/string_network_universe.tsv,
    required_score>=0.7): physical_PPI if the pair is ALSO present in
    STRING's own `network_type=physical` query result for the same node
    set (data/reference/interactions/string_network_universe_physical.tsv)
    -- STRING's physical subnetwork is curated by STRING itself to exclude
    purely functional-association-only evidence channels (pure
    coexpression/textmining/genomic-context pairs), which is the correct
    authority for this distinction. An earlier version of this module used
    a manual escore>0-or-dscore>0 heuristic instead, which a Codex review
    correctly flagged as too permissive (dscore, STRING's curated-database
    channel, includes non-physical curated associations, e.g. pathway
    co-membership annotations; a nonzero threshold on either channel let
    29 non-physical pairs through). Everything present in the default
    (functional) network query but absent from the physical one is labeled
    functional_association -- never "physical".
  - TRRUST (data/reference/interactions/trrust_human.tsv): regulatory
    (TF_target) edges, directional (source=TF, target=target gene).
  - pathway_co_membership: two node-universe genes that are both GSEA
    leading-edge members of the same STRONG_CONSENSUS or MULTIMODAL_PATHWAY
    Hallmark/Reactome pathway (GO:BP excluded here for the same
    combinatorial-redundancy reason as Phase 8/14), AND that pathway's
    leading edge within the node universe has at most MAX_PATHWAY_GENES_FOR_EDGES
    (10) genes. Large generic pathways (e.g. REACTOME_CELL_CYCLE_MITOTIC, 39
    node-universe genes) would otherwise turn pairwise co-membership into a
    near-complete graph (4,812 edges among 119 nodes at no size cap -- not a
    selective signal, an explicit anti-pattern per Phase 23's "do not make a
    dense unreadable hairball"). Excluding a pathway from edge generation
    does not remove its evidence: it remains fully recorded in
    multimodal_pathway_convergence.tsv and each node's `pathways_supported`
    attribute (Phase 15) -- only the pairwise-edge representation is capped.

No edge is inferred from expression correlation alone anywhere in this
module.
"""

from __future__ import annotations

import itertools
import logging
from pathlib import Path

import pandas as pd
import yaml

logger = logging.getLogger(__name__)

RESISTANCE_DATASETS = ["gse118713", "gse240112", "gse111151"]
MAX_PATHWAY_GENES_FOR_EDGES = 10


def _load_config(config_path: str | Path) -> dict:
    with open(config_path) as f:
        return yaml.safe_load(f)


def build_string_edges(path: str | Path, physical_path: str | Path, node_genes: set[str]) -> pd.DataFrame:
    # STRING's `network` endpoint can occasionally return a small number of
    # rows involving an identifier outside the exact submitted list (STRING's
    # own alias resolution can expand/substitute a queried symbol) -- both
    # endpoints of every returned row are explicitly re-filtered against
    # `node_genes` here, so the downloaded file is read as an unfiltered raw
    # STRING response, never assumed to already match the node universe
    # exactly (Codex review, second pass).
    df = pd.read_csv(path, sep="\t")
    df = df.loc[df["preferredName_A"].isin(node_genes) & df["preferredName_B"].isin(node_genes)]

    physical_df = pd.read_csv(physical_path, sep="\t")
    physical_pairs = set(zip(physical_df["preferredName_A"], physical_df["preferredName_B"])) | set(zip(physical_df["preferredName_B"], physical_df["preferredName_A"]))

    rows = []
    for _, row in df.iterrows():
        is_physical = (row["preferredName_A"], row["preferredName_B"]) in physical_pairs
        interaction_type = "physical_PPI" if is_physical else "functional_association"
        rows.append(
            {
                "source_gene": row["preferredName_A"],
                "target_gene": row["preferredName_B"],
                "interaction_type": interaction_type,
                "database_source": "STRING",
                "confidence": row["score"],
                "pathway": "",
                "evidence_notes": f"combined_score={row['score']:.3f}; escore={row['escore']:.3f}; dscore={row['dscore']:.3f}; STRING network_type=physical membership={is_physical}",
            }
        )
    out = pd.DataFrame(rows)
    logger.info("build_string_edges: %d edges (%d physical_PPI, %d functional_association)", len(out), (out["interaction_type"] == "physical_PPI").sum(), (out["interaction_type"] == "functional_association").sum())
    return out


def build_trrust_edges(path: str | Path, node_genes: set[str]) -> pd.DataFrame:
    """A handful of TRRUST entries are genuine curated TF autoregulation
    (e.g. TP53 activating its own transcription, PMID 22532570) -- real
    biology, but excluded here: a self-loop adds no inter-gene connectivity
    information for this network's purpose (candidate/module bridging) and
    would silently inflate that gene's degree in the Phase 20 hub analysis."""
    df = pd.read_csv(path, sep="\t", header=None, names=["tf", "target", "mode", "pmid"])
    df = df.loc[df["tf"].isin(node_genes) & df["target"].isin(node_genes) & (df["tf"] != df["target"])]
    out = pd.DataFrame(
        {
            "source_gene": df["tf"],
            "target_gene": df["target"],
            "interaction_type": "regulatory",
            "database_source": "TRRUST",
            "confidence": pd.NA,
            "pathway": "",
            "evidence_notes": df.apply(lambda r: f"mode={r['mode']}; PMID={r['pmid']}", axis=1),
        }
    )
    logger.info("build_trrust_edges: %d TF_target edges", len(out))
    return out


def build_pathway_co_membership_edges(tables_dir: Path, node_genes: set[str]) -> pd.DataFrame:
    consensus = pd.read_csv(tables_dir / "resistance_pathway_consensus.tsv", sep="\t")
    mm = pd.read_csv(tables_dir / "multimodal_pathway_convergence.tsv", sep="\t")

    keys = set(zip(consensus.loc[consensus["consensus_category"] == "STRONG_CONSENSUS", "collection"], consensus.loc[consensus["consensus_category"] == "STRONG_CONSENSUS", "pathway"]))
    keys |= set(zip(mm.loc[mm["convergence_category"] == "MULTIMODAL_PATHWAY", "collection"], mm.loc[mm["convergence_category"] == "MULTIMODAL_PATHWAY", "pathway"]))
    keys = {(c, p) for c, p in keys if c in ("hallmark", "reactome")}

    rows = []
    seen_pairs: set[tuple[str, str, str]] = set()
    for dataset in RESISTANCE_DATASETS:
        df = pd.read_csv(tables_dir / f"gsea_{dataset}.tsv", sep="\t")
        df = df.loc[df.apply(lambda r: (r["collection"], r["pathway"]) in keys, axis=1)]
        # Codex review fix: only draw a co-membership edge from a dataset
        # where that dataset's own nominal p-value for the pathway is
        # significant, not merely because the pathway is STRONG_CONSENSUS/
        # MULTIMODAL_PATHWAY overall (which can be driven by the other two
        # datasets alone).
        df = df.loc[df["nom_pvalue"] < 0.05]
        for _, row in df.iterrows():
            if pd.isna(row["leading_edge_genes"]):
                continue
            genes_in_universe = sorted(set(row["leading_edge_genes"].split(";")) & node_genes)
            if len(genes_in_universe) > MAX_PATHWAY_GENES_FOR_EDGES:
                continue
            for a, b in itertools.combinations(genes_in_universe, 2):
                pair_key = tuple(sorted([a, b])) + (row["pathway"],)
                if pair_key in seen_pairs:
                    continue
                seen_pairs.add(pair_key)
                rows.append(
                    {
                        "source_gene": a,
                        "target_gene": b,
                        "interaction_type": "pathway_co_membership",
                        "database_source": f"MSigDB_{row['collection']}",
                        "confidence": pd.NA,
                        "pathway": row["pathway"],
                        "evidence_notes": f"both leading-edge in {dataset}",
                    }
                )
    out = pd.DataFrame(rows)
    logger.info("build_pathway_co_membership_edges: %d edges across %d pathways", len(out), len(keys))
    return out


def build_all_edges(config_path: str | Path = "config/config.yaml") -> pd.DataFrame:
    config = _load_config(config_path)
    cfg = config["systems_network"]
    tables_dir = Path(cfg["output"]["tables_dir"])
    networks_dir = Path(cfg["output"]["networks_dir"])

    node_genes = set(pd.read_csv(networks_dir / "node_universe.tsv", sep="\t")["gene"])

    string_edges = build_string_edges("data/reference/interactions/string_network_universe.tsv", "data/reference/interactions/string_network_universe_physical.tsv", node_genes)
    trrust_edges = build_trrust_edges("data/reference/interactions/trrust_human.tsv", node_genes)
    pathway_edges = build_pathway_co_membership_edges(tables_dir, node_genes)

    all_edges = pd.concat([string_edges, trrust_edges, pathway_edges], ignore_index=True)
    all_edges.to_csv(networks_dir / "edges.tsv", sep="\t", index=False)
    logger.info("build_all_edges: %d total edges (%s)", len(all_edges), all_edges["interaction_type"].value_counts().to_dict())
    return all_edges


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    build_all_edges()
