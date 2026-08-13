"""Systems-network phase 9: recurrent leading-edge genes across the three
resistance datasets, restricted to STRONG_CONSENSUS pathways (significant,
same-direction enrichment in >=2 of GSE118713/GSE240112/GSE111151, Phase 5).

Two distinct recurrence statistics are reported per gene (Codex review,
second pass -- an earlier version reported only the first and it was easy
to over-read as the second):
  - `leading_edge_dataset_count`: number of DISTINCT datasets in which the
    gene is a significant leading-edge member of ANY qualifying pathway
    (can be reached via a different pathway in each dataset).
  - `max_same_pathway_dataset_count`: the stricter statistic -- the
    largest number of datasets for which the gene is a significant
    leading-edge member of the SAME single pathway. This is the one that
    actually supports "reproducibly drives this specific enrichment
    signal" and is what node-universe category B (Phase 14) ranks by.
"""

from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd
import yaml

logger = logging.getLogger(__name__)

RESISTANCE_DATASETS = ["gse118713", "gse240112", "gse111151"]
NOMINAL_ENRICHED_P = 0.05


def _load_config(config_path: str | Path) -> dict:
    with open(config_path) as f:
        return yaml.safe_load(f)


def build_recurrent_leading_edge_genes(
    consensus: pd.DataFrame, gsea_tables: dict[str, pd.DataFrame], crispr_ranked: pd.DataFrame
) -> pd.DataFrame:
    strong_keys = set(zip(consensus.loc[consensus["consensus_category"] == "STRONG_CONSENSUS", "collection"], consensus.loc[consensus["consensus_category"] == "STRONG_CONSENSUS", "pathway"]))
    consensus_dir = consensus.set_index(["collection", "pathway"])["median_NES"]

    records = []
    for dataset in RESISTANCE_DATASETS:
        df = gsea_tables[dataset]
        df = df.loc[df.apply(lambda r: (r["collection"], r["pathway"]) in strong_keys, axis=1)]
        # Codex review fix: a pathway being STRONG_CONSENSUS (significant in
        # >=2 of the 3 datasets, Phase 5) does not mean EVERY dataset's own
        # result for that pathway was itself significant -- gseapy still
        # reports a leading_edge_genes list for the third, non-significant
        # dataset too. Only count a gene as leading-edge-in-this-dataset if
        # that dataset's own nominal p-value for the pathway is <0.05
        # (NOMINAL_ENRICHED_P, same threshold used in Phase 7), otherwise
        # "recurrent across datasets" silently included non-recurrent,
        # nonsignificant provenance.
        df = df.loc[df["nom_pvalue"] < NOMINAL_ENRICHED_P]
        for _, row in df.iterrows():
            if pd.isna(row["leading_edge_genes"]):
                continue
            for gene in row["leading_edge_genes"].split(";"):
                records.append({"gene": gene, "dataset": dataset, "collection": row["collection"], "pathway": row["pathway"], "pathway_NES": row["NES"]})

    long_df = pd.DataFrame(records)
    logger.info("build_recurrent_leading_edge_genes: %d (gene,dataset,pathway) leading-edge rows over %d STRONG_CONSENSUS pathways", len(long_df), len(strong_keys))

    rows = []
    for gene, grp in long_df.groupby("gene"):
        n_datasets = grp["dataset"].nunique()
        n_pathways = grp[["collection", "pathway"]].drop_duplicates().shape[0]
        # Codex review note (second pass): `leading_edge_dataset_count` is a
        # gene-level aggregate across ALL qualifying pathways -- a gene can
        # reach count=3 via three DIFFERENT pathways (one per dataset)
        # without ever being significantly leading-edge for the SAME
        # pathway in more than one or two of them. `max_same_pathway_dataset_count`
        # is the stricter, unambiguous "reproducibly drives this specific
        # enrichment signal" statistic -- the max, over any single
        # (collection, pathway), of how many datasets that gene was
        # leading-edge-and-significant for. Both are reported; node-universe
        # category B (Phase 14) uses the stricter one.
        per_pathway_dataset_count = grp.groupby(["collection", "pathway"])["dataset"].nunique()
        max_same_pathway_dataset_count = int(per_pathway_dataset_count.max())
        n_pos = (grp["pathway_NES"] > 0).sum()
        n_neg = (grp["pathway_NES"] < 0).sum()
        if n_pos > 0 and n_neg > 0:
            direction = "mixed"
        elif n_pos > 0:
            direction = "consistently_up_in_resistance_pathways"
        else:
            direction = "consistently_down_in_resistance_pathways"
        rows.append(
            {
                "gene": gene,
                "leading_edge_dataset_count": n_datasets,
                "max_same_pathway_dataset_count": max_same_pathway_dataset_count,
                "leading_edge_datasets": ",".join(sorted(grp["dataset"].unique())),
                "pathway_count": n_pathways,
                "resistance_direction": direction,
                "example_pathways": ";".join(sorted(grp[["collection", "pathway"]].drop_duplicates().apply(lambda r: f"{r['collection']}:{r['pathway']}", axis=1))[:10]),
            }
        )
    out = pd.DataFrame(rows).merge(crispr_ranked[["gene", "log2fc", "p_value", "fdr"]].rename(columns={"log2fc": "crispr_effect", "p_value": "crispr_p", "fdr": "crispr_fdr"}), on="gene", how="left")
    out["crispr_direction"] = out.apply(_classify_crispr, axis=1)
    out = out.sort_values(["leading_edge_dataset_count", "pathway_count"], ascending=False)
    return out


def _classify_crispr(row: pd.Series) -> str:
    if pd.isna(row["crispr_fdr"]):
        return "not_tested"
    if row["crispr_fdr"] < 0.05:
        return "sensitising_KO" if row["crispr_effect"] < 0 else "tolerance_associated_KO"
    return "nonsignificant_sensitising_direction" if row["crispr_effect"] < 0 else "nonsignificant_tolerance_direction"


def run_leading_edge(config_path: str | Path = "config/config.yaml") -> pd.DataFrame:
    config = _load_config(config_path)
    cfg = config["systems_network"]
    tables_dir = Path(cfg["output"]["tables_dir"])

    consensus = pd.read_csv(tables_dir / "resistance_pathway_consensus.tsv", sep="\t")
    gsea_tables = {d: pd.read_csv(tables_dir / f"gsea_{d}.tsv", sep="\t") for d in RESISTANCE_DATASETS}
    crispr_ranked = pd.read_csv(tables_dir / "crispr_ranked_genes.tsv", sep="\t")

    out = build_recurrent_leading_edge_genes(consensus, gsea_tables, crispr_ranked)
    out.to_csv(tables_dir / "recurrent_leading_edge_genes.tsv", sep="\t", index=False)
    return out


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run_leading_edge()
