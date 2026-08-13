"""Systems-network phase 4 (+11): preranked GSEA for each dataset's signed
genome-wide ranking (src.systems_network_ranking) against Hallmark,
Reactome, and GO:BP (data/reference/genesets/, see
docs/SYSTEMS_NETWORK_GENESETS.md).

Implementation: gseapy.prerank(method="multilevel"), gseapy's own
documentation describes this method as "a faithful port of the fgsea C++
core" -- fgsea itself is unavailable in this project's R environment (see
docs/SYSTEMS_NETWORK_INPUT_AUDIT.md), so this is the closest available
implementation of the same algorithm, not a different statistical method.

Each dataset is run independently against each of the three gene-set
collections (never combined into one FDR null across datasets); results
from the three collections are concatenated into one per-dataset output
table tagged by `collection`.
"""

from __future__ import annotations

import logging
from pathlib import Path

import gseapy
import pandas as pd
import yaml

logger = logging.getLogger(__name__)

COLLECTIONS = ["hallmark", "reactome", "go_bp"]


def _load_config(config_path: str | Path) -> dict:
    with open(config_path) as f:
        return yaml.safe_load(f)


def _parse_tag_pct(tag_pct: str) -> tuple[int, int]:
    n_lead, n_set = tag_pct.split("/")
    return int(n_lead), int(n_set)


def run_gsea_one_collection(rnk: pd.DataFrame, gmt_path: str, collection: str, cfg: dict, seed: int) -> pd.DataFrame:
    result = gseapy.prerank(
        rnk=rnk,
        gene_sets=gmt_path,
        min_size=cfg["genesets"]["min_size"],
        max_size=cfg["genesets"]["max_size"],
        permutation_num=cfg["gsea"]["n_permutations"],
        seed=seed,
        threads=4,
        no_plot=True,
        outdir=None,
        method="multilevel",
    )
    res = result.res2d.copy()
    lead_n, set_n = zip(*res["Tag %"].map(_parse_tag_pct))
    out = pd.DataFrame(
        {
            "pathway": res["Term"],
            "collection": collection,
            "NES": res["NES"].astype(float),
            "nom_pvalue": res["NOM p-val"].astype(float),
            "fdr": res["FDR q-val"].astype(float),
            "gene_set_size_in_data": set_n,
            "n_leading_edge": lead_n,
            "leading_edge_genes": res["Lead_genes"],
        }
    )
    return out


def run_gsea_for_ranking(rnk_path: str | Path, dataset_label: str, cfg: dict) -> pd.DataFrame:
    rnk = pd.read_csv(rnk_path, sep="\t")[["gene", "ranking_stat"]]
    seed = cfg["gsea"]["seed"]
    frames = []
    for collection in COLLECTIONS:
        gmt_path = cfg["genesets"][f"{collection}_gmt"]
        frame = run_gsea_one_collection(rnk, gmt_path, collection, cfg, seed)
        frame["dataset"] = dataset_label
        frames.append(frame)
        logger.info("run_gsea_for_ranking(%s, %s): %d pathways tested", dataset_label, collection, len(frame))
    out = pd.concat(frames, ignore_index=True)
    cols = ["dataset", "collection", "pathway", "NES", "nom_pvalue", "fdr", "gene_set_size_in_data", "n_leading_edge", "leading_edge_genes"]
    return out[cols]


def run_gsea(config_path: str | Path = "config/config.yaml") -> dict[str, pd.DataFrame]:
    config = _load_config(config_path)
    cfg = config["systems_network"]
    tables_dir = Path(cfg["output"]["tables_dir"])

    jobs = {
        "gse118713": (tables_dir / "gse118713_ranked_genes.tsv", "gsea_gse118713.tsv"),
        "gse240112": (tables_dir / "gse240112_ranked_genes.tsv", "gsea_gse240112.tsv"),
        "gse111151": (tables_dir / "gse111151_ranked_genes.tsv", "gsea_gse111151.tsv"),
        "gse245601": (tables_dir / "gse245601_ranked_genes.tsv", "gsea_gse245601.tsv"),
        "gse245601_track_b": (tables_dir / "gse245601_track_b_ranked_genes.tsv", "gsea_gse245601_track_b.tsv"),
        "crispr": (tables_dir / "crispr_ranked_genes.tsv", "gsea_crispr.tsv"),
    }

    results = {}
    for label, (rnk_path, out_name) in jobs.items():
        logger.info("Running GSEA for %s", label)
        res = run_gsea_for_ranking(rnk_path, label, cfg)
        res.to_csv(tables_dir / out_name, sep="\t", index=False)
        results[label] = res
    return results


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run_gsea()
