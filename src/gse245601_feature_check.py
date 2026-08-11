"""GSE245601 candidate-gene feature-space availability check.

Data source: one representative Cell Ranger filtered_feature_bc_matrix.h5
from GSE245601 (10x Genomics GRCh38 reference, 33,538 features, single
"Gene Expression" feature type, single "GRCh38" genome -- confirmed
identical across all 20 primary tumor samples by SHA-256 of the feature ID
vector before this module was written).

This module reads ONLY the feature metadata (Ensembl gene ID + gene
symbol) of the H5 file's `matrix/features` group -- it never reads the
expression matrix (`matrix/data`) itself, consistent with the frozen
project rule that candidate-gene expression must not be inspected during
the preprocessing/malignant-cell-identification phase.
"""

from __future__ import annotations

from pathlib import Path

import h5py
import pandas as pd
import yaml

# Frozen candidate set (13 sensitising CRISPR hits) plus PAICS
# (published benchmark, reported separately, never pooled with the 13).
# Ensembl IDs are the same unique mappings already frozen and validated in
# this project's GSE118713 bulk RNA-seq pipeline
# (results/tables/crispr_gse118713_master_table.tsv /
# candidate_paics_benchmark.tsv), except where noted.
CANDIDATE_ENSEMBL_IDS: dict[str, str] = {
    "USP34": "ENSG00000115464",
    "CTDNEP1": "ENSG00000175826",
    "EIF4ENIF1": "ENSG00000184708",
    "HMGB1": "ENSG00000189403",
    "KDM1A": "ENSG00000004487",
    "PET117": "ENSG00000232838",
    "TADA2B": "ENSG00000173011",
    "VEZF1": "ENSG00000136451",
    "ICK": "ENSG00000112144",
    "SUPT4H1": "ENSG00000213246",
    "TLK2": "ENSG00000146872",
    "TSR3": "ENSG00000007520",
    "USP17L29": "ENSG00000231637",
}
# PAICS Ensembl ID confirmed directly against the GSE245601 Cell Ranger
# GRCh38-2020-A-style reference feature list (this module) -- corrects an
# earlier, incorrect recollection (ENSG00000128646) from the prior
# design-only audit turn, which was never checked against real reference
# data at that stage.
PAICS_ENSEMBL_ID = "ENSG00000128050"


def _load_config(config_path: str | Path) -> dict:
    with open(config_path) as f:
        return yaml.safe_load(f)


def read_feature_table(h5_path: str | Path) -> pd.DataFrame:
    """Read only the feature (gene) metadata from a Cell Ranger
    filtered_feature_bc_matrix.h5 file -- ids, names, feature types, and
    genome. Does not touch the expression matrix."""
    with h5py.File(h5_path, "r") as f:
        ids = [x.decode() for x in f["matrix/features/id"][:]]
        names = [x.decode() for x in f["matrix/features/name"][:]]
        feature_types = [x.decode() for x in f["matrix/features/feature_type"][:]]
        genomes = [x.decode() for x in f["matrix/features/genome"][:]]
    return pd.DataFrame({"ensembl_id": ids, "gene_symbol": names, "feature_type": feature_types, "genome": genomes})


def check_candidate_feature_availability(feature_table: pd.DataFrame) -> pd.DataFrame:
    """For each of the 13 candidates + PAICS, report presence/absence in
    the feature table by symbol and by Ensembl ID, and flag ambiguity
    (more than one matching row). No expression values are read or
    reported here."""
    all_genes = dict(CANDIDATE_ENSEMBL_IDS)
    all_genes["PAICS"] = PAICS_ENSEMBL_ID

    rows = []
    for gene_symbol, expected_ensembl_id in all_genes.items():
        by_symbol = feature_table.loc[feature_table["gene_symbol"] == gene_symbol]
        by_ensembl = feature_table.loc[feature_table["ensembl_id"] == expected_ensembl_id]

        n_symbol_matches = len(by_symbol)
        n_ensembl_matches = len(by_ensembl)

        if n_symbol_matches == 0:
            status = "absent"
        elif n_symbol_matches > 1:
            status = "ambiguous"
        elif n_ensembl_matches != 1:
            status = "ambiguous"
        elif by_symbol.iloc[0]["ensembl_id"] != expected_ensembl_id:
            status = "ambiguous"
        else:
            status = "present"

        rows.append(
            {
                "gene_symbol": gene_symbol,
                "is_paics_benchmark": gene_symbol == "PAICS",
                "expected_ensembl_id": expected_ensembl_id,
                "n_symbol_matches": n_symbol_matches,
                "n_ensembl_id_matches": n_ensembl_matches,
                "matched_ensembl_id": by_symbol.iloc[0]["ensembl_id"] if n_symbol_matches >= 1 else None,
                "feature_type": by_symbol.iloc[0]["feature_type"] if n_symbol_matches >= 1 else None,
                "status": status,
            }
        )
    return pd.DataFrame(rows)


def run_feature_check(config_path: str | Path = "config/config.yaml") -> pd.DataFrame:
    config = _load_config(config_path)
    cfg = config["gse245601"]
    h5_dir = Path(cfg["h5_dir"])
    representative_sample = next(s for s in cfg["samples"] if s["primary_analysis"])
    h5_path = h5_dir / representative_sample["h5"]

    feature_table = read_feature_table(h5_path)
    result = check_candidate_feature_availability(feature_table)

    output_path = Path(cfg["output"]["feature_availability_tsv"])
    output_path.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(output_path, sep="\t", index=False)
    return result


if __name__ == "__main__":
    result = run_feature_check()
    print(result.to_string(index=False))
