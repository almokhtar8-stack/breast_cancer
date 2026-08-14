"""Shared TCGA-BRCA data loaders for the independent-validation phase.

Data source: GDC-harmonized STAR-pipeline RNA-seq (log2(TPM+1)) mirrored by
UCSC Xena's GDC hub, plus cBioPortal clinical annotation from two TCGA-BRCA
studies (brca_tcga_pub for IHC receptor status/PAM50, brca_tcga_pan_can_atlas_2018
for curated survival). Exact URLs, access rationale, and why three sources
are needed are documented in scripts/download/download_tcga_brca.py's
module docstring; this module only reads the already-downloaded raw files
it pins via config/config.yaml's independent_validation.tcga section.

No network calls happen here. Every function logs rows in/out.
"""

from __future__ import annotations

import json
import logging
from pathlib import Path

import pandas as pd
import yaml

logger = logging.getLogger(__name__)

# Official TCGA barcode sample-type codes (positions 14-15 of the barcode),
# per https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/sample-type-codes.
# A fixed, versionless code table -- not something that changes per release.
SAMPLE_TYPE_LABELS = {
    "01": "Primary Solid Tumor",
    "02": "Recurrent Solid Tumor",
    "03": "Primary Blood Derived Cancer - Peripheral Blood",
    "05": "Additional - New Primary",
    "06": "Metastatic",
    "07": "Additional Metastatic",
    "10": "Blood Derived Normal",
    "11": "Solid Tissue Normal",
}


def load_config(config_path: str | Path = "config/config.yaml") -> dict:
    with open(config_path) as f:
        return yaml.safe_load(f)


def _raw_dir(cfg: dict) -> Path:
    return Path(cfg["data"]["raw"]["tcga_brca_dir"])


def load_expression(cfg: dict, genes: list[str] | None = None) -> pd.DataFrame:
    """Return samples x genes log2(TPM+1) matrix, columns = gene symbols.

    If ``genes`` is given, only those Ensembl IDs (looked up via
    data/reference/tcga_candidate_ensembl_ids.tsv) are kept; otherwise the
    full matrix (~60,660 genes) is returned for pathway-score computation.
    """
    raw = _raw_dir(cfg)
    path = raw / cfg["independent_validation"]["tcga"]["raw"]["star_tpm_tsv_gz"]
    df = pd.read_csv(path, sep="\t", index_col=0)
    logger.info("load_expression: read %s (%d genes x %d samples)", path, *df.shape)
    df.index = df.index.str.split(".").str[0]  # drop Ensembl version suffix
    n_before = df.shape[0]
    df = df[~df.index.duplicated(keep="first")]
    if df.shape[0] != n_before:
        logger.info("load_expression: dropped %d duplicate unversioned Ensembl IDs", n_before - df.shape[0])
    df = df.T  # samples x genes
    if genes is not None:
        candidate_ids = pd.read_csv("data/reference/tcga_candidate_ensembl_ids.tsv", sep="\t")
        id_map = dict(zip(candidate_ids["symbol"], candidate_ids["ensembl_gene_id_unversioned"]))
        missing = [g for g in genes if g not in id_map]
        if missing:
            raise ValueError(f"load_expression: no pinned Ensembl ID for {missing}")
        cols = {id_map[g]: g for g in genes}
        absent = [eid for eid in cols if eid not in df.columns]
        if absent:
            raise ValueError(f"load_expression: Ensembl ID(s) {absent} not found in expression matrix")
        df = df[list(cols)].rename(columns=cols)
        logger.info("load_expression: subset to %d requested genes", len(genes))
    return df


def load_expression_symbols(cfg: dict, sample_barcodes: list[str] | None = None) -> pd.DataFrame:
    """Full samples x gene-symbol log2(TPM+1) matrix for pathway scoring.

    Ensembl IDs are mapped to HGNC symbols via the GENCODE v36 probemap
    (the exact annotation GDC's STAR pipeline used for this matrix).
    Genes with no symbol are dropped; genes where multiple Ensembl IDs map
    to the same symbol are collapsed by mean (both counts are logged).
    """
    raw = _raw_dir(cfg)
    expr_path = raw / cfg["independent_validation"]["tcga"]["raw"]["star_tpm_tsv_gz"]
    probemap_path = raw / cfg["independent_validation"]["tcga"]["raw"]["gencode_v36_probemap"]

    df = pd.read_csv(expr_path, sep="\t", index_col=0)
    n_genes_in = df.shape[0]
    if sample_barcodes is not None:
        df = df[sample_barcodes]

    probemap = pd.read_csv(probemap_path, sep="\t")
    id_to_symbol = dict(zip(probemap["id"], probemap["gene"]))
    df["symbol"] = df.index.map(id_to_symbol)
    n_no_symbol = df["symbol"].isna().sum()
    df = df.dropna(subset=["symbol"])
    n_dup_ensembl_per_symbol = df["symbol"].duplicated().sum()
    df = df.groupby("symbol").mean(numeric_only=True)
    logger.info(
        "load_expression_symbols: %d Ensembl rows in -> %d dropped (no symbol) -> "
        "%d Ensembl rows collapsed into duplicate symbols -> %d unique gene symbols out",
        n_genes_in, n_no_symbol, n_dup_ensembl_per_symbol, df.shape[0],
    )
    return df.T  # samples x gene symbols


def barcode_to_patient(sample_barcode: str) -> str:
    return sample_barcode[:12]


def barcode_sample_type(sample_barcode: str) -> str:
    code = sample_barcode.split("-")[3][:2]
    return SAMPLE_TYPE_LABELS.get(code, f"Unknown code {code}")


def load_receptor_pam50(cfg: dict) -> pd.DataFrame:
    """Sample-level ER/PR/HER2/PAM50 calls from cBioPortal brca_tcga_pub.

    Returns one row per patient barcode (brca_tcga_pub sample IDs already
    equal the 12-character patient barcode for this cohort, since it is a
    per-patient discovery-cohort study), with an explicit
    missingness-tracking column per field.
    """
    raw = _raw_dir(cfg)
    path = raw / cfg["independent_validation"]["tcga"]["raw"]["receptor_clinical_sample_json"]
    records = json.loads(path.read_text())
    n_records = len(records)
    df = pd.DataFrame(records)
    wide = df.pivot_table(index="patientId", columns="clinicalAttributeId", values="value", aggfunc="first")
    wide.index.name = "patient_barcode"
    keep = [c for c in ["ER_STATUS", "PR_STATUS", "HER2_STATUS", "PAM50_SUBTYPE", "SAMPLE_TYPE"] if c in wide.columns]
    wide = wide[keep]
    logger.info(
        "load_receptor_pam50: %s -> %d records -> %d patients; missing ER_STATUS for %d/%d",
        path, n_records, len(wide), wide["ER_STATUS"].isna().sum() if "ER_STATUS" in wide else len(wide), len(wide),
    )
    return wide


def load_survival_clinical(cfg: dict) -> pd.DataFrame:
    """Patient-level OS/PFS/DFS/AGE/STAGE from cBioPortal brca_tcga_pan_can_atlas_2018."""
    raw = _raw_dir(cfg)
    path = raw / cfg["independent_validation"]["tcga"]["raw"]["survival_clinical_json"]
    records = json.loads(path.read_text())
    n_records = len(records)
    df = pd.DataFrame(records)
    wide = df.pivot_table(index="patientId", columns="clinicalAttributeId", values="value", aggfunc="first")
    wide.index.name = "patient_barcode"
    logger.info("load_survival_clinical: %s -> %d records -> %d patients", path, n_records, len(wide))
    return wide


def build_cohort_table(cfg: dict) -> pd.DataFrame:
    """One row per TCGA-BRCA patient-x-sample-type combination (barcode,
    patient, sample type, ER/PR/HER2/PAM50, survival/stage/age), with
    missingness reported via row/column NaN counts (never silently
    imputed).

    A small number of patients (11 in this cohort) have more than one
    primary-tumor RNA-seq aliquot (e.g. barcode suffix -01A and -01B from
    the same patient); every other TCGA sample-type category here has at
    most one sample per patient. Keeping every aliquot would silently
    duplicate that patient's outcome/covariates in any per-sample
    comparison or Cox model (pseudoreplication). One aliquot per
    patient-x-sample-type combination is kept, chosen deterministically by
    the lexicographically first sample barcode (which prefers vial "A"
    over "B"/"C" when multiple exist); the rest are dropped and logged.
    """
    expr_samples = load_expression(cfg, genes=["CITED2"]).index  # cheap: 1-gene pull just to get sample list
    cohort = pd.DataFrame({"sample_barcode": expr_samples})
    cohort["patient_barcode"] = cohort["sample_barcode"].map(barcode_to_patient)
    cohort["sample_type"] = cohort["sample_barcode"].map(barcode_sample_type)
    cohort["is_primary_tumor"] = cohort["sample_type"] == "Primary Solid Tumor"
    cohort["is_normal"] = cohort["sample_type"] == "Solid Tissue Normal"
    cohort["is_metastatic"] = cohort["sample_type"] == "Metastatic"

    cohort = cohort.sort_values("sample_barcode")
    n_before_dedup = len(cohort)
    dup_mask = cohort.duplicated(subset=["patient_barcode", "sample_type"], keep="first")
    n_dropped = int(dup_mask.sum())
    if n_dropped:
        logger.info(
            "build_cohort_table: dropping %d duplicate-aliquot sample(s) (same patient + same sample_type; "
            "kept the lexicographically first barcode per group to avoid pseudoreplication): %s",
            n_dropped, cohort.loc[dup_mask, "sample_barcode"].tolist(),
        )
    cohort = cohort.loc[~dup_mask]
    cohort = cohort.set_index("sample_barcode")

    receptor = load_receptor_pam50(cfg)
    survival = load_survival_clinical(cfg)
    cohort = cohort.join(receptor, on="patient_barcode")
    cohort = cohort.join(survival, on="patient_barcode", rsuffix="_survival")

    n_er_missing = cohort.loc[cohort["is_primary_tumor"], "ER_STATUS"].isna().sum() if "ER_STATUS" in cohort else len(cohort)
    n_primary = int(cohort["is_primary_tumor"].sum())
    logger.info(
        "build_cohort_table: %d samples in -> %d dropped as duplicate aliquots -> %d one-per-patient-per-sample-type samples "
        "(%d primary tumor, %d normal, %d metastatic); ER_STATUS missing for %d/%d primary tumors",
        n_before_dedup, n_dropped, len(cohort), n_primary, int(cohort["is_normal"].sum()), int(cohort["is_metastatic"].sum()),
        n_er_missing, n_primary,
    )
    return cohort
