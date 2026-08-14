"""Download the TCGA-BRCA expression matrix and clinical/receptor-status
annotation used by the independent human-cohort validation phase.

Three official sources are pulled, each documented separately because none
of them alone carries everything this project needs (harmonized GDC
expression, IHC-based ER/PR/HER2 + PAM50, and long-follow-up curated
survival each live in a different official release):

1. GENE EXPRESSION -- GDC-harmonized RNA-seq (STAR pipeline), TPM,
   log2(TPM+1), mirrored by UCSC Xena's GDC hub (a standard, widely-cited
   secondary distribution of the same GDC-processed files; used here
   because GDC's own per-case file API would require ~1,200 individual
   file downloads for the same content). Cohort: "GDC TCGA Breast Cancer
   (BRCA)". File: TCGA-BRCA.star_tpm.tsv.gz.
   https://gdc.xenahubs.net/download/TCGA-BRCA.star_tpm.tsv.gz
   Sample type (primary tumor / solid tissue normal / metastatic) is
   never looked up from a separate phenotype file: it is decoded directly
   from the official TCGA barcode sample-type code (characters 14-15,
   e.g. "01"=Primary Solid Tumor, "11"=Solid Tissue Normal, "06"=
   Metastatic; https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/sample-type-codes).

2. CLINICAL RECEPTOR STATUS + PAM50 -- cBioPortal study brca_tcga_pub
   (TCGA, Nature 2012 discovery cohort, PMID 23000897), which carries the
   original IHC-based ER_STATUS / PR_STATUS / HER2_STATUS and PAM50_SUBTYPE
   calls in one consistent clinical file. This is the PRIMARY clinical
   receptor-status source for cohort construction (Part 2).
   https://www.cbioportal.org/api/studies/brca_tcga_pub/clinical-data

3. SURVIVAL + STAGE/AGE -- cBioPortal study brca_tcga_pan_can_atlas_2018
   (TCGA PanCancer Atlas, Cell 2018, PMID 29625048), which carries the
   most current, longest-follow-up curated OS/PFS/DFS endpoints plus AGE
   and AJCC_PATHOLOGIC_TUMOR_STAGE. This is a DIFFERENT TCGA data freeze
   than source 2 above; the two are joined only by shared patient barcode
   identity, and missingness from imperfect overlap is expected and must
   be reported, not silently dropped.
   https://www.cbioportal.org/api/studies/brca_tcga_pan_can_atlas_2018/clinical-data

This script only downloads. It performs no analysis and is never imported
by src/ or tests/ -- run it manually, once, to populate the local raw-data
cache referenced by config/config.yaml's tcga_brca.raw section.
"""

from __future__ import annotations

import argparse
import gzip
import json
import logging
import urllib.request
from pathlib import Path

logger = logging.getLogger(__name__)

XENA_BASE = "https://gdc.xenahubs.net/download"
XENA_FILES = ["TCGA-BRCA.star_tpm.tsv.gz", "TCGA-BRCA.survival.tsv.gz"]
# GENCODE v36 gene probemap (Ensembl gene ID -> HGNC symbol), the exact
# annotation version GDC's STAR pipeline used to build TCGA-BRCA.star_tpm --
# needed to translate the expression matrix to gene symbols for pathway
# scoring without a version mismatch against a different GENCODE release.
XENA_PROBEMAP = "gencode.v36.annotation.gtf.gene.probemap"

CBIOPORTAL_BASE = "https://www.cbioportal.org/api"
CBIOPORTAL_RECEPTOR_STUDY = "brca_tcga_pub"
CBIOPORTAL_SURVIVAL_STUDY = "brca_tcga_pan_can_atlas_2018"


def _fetch_bytes(url: str, timeout: int = 300) -> bytes:
    logger.info("download_tcga_brca: fetching %s", url)
    request = urllib.request.Request(url, headers={"User-Agent": "breast-cancer-project/1.0"})
    with urllib.request.urlopen(request, timeout=timeout) as response:
        return response.read()


def _fetch_json(url: str, timeout: int = 300):
    return json.loads(_fetch_bytes(url, timeout=timeout))


def _download_xena(out_dir: Path) -> None:
    for name in XENA_FILES:
        dest = out_dir / name
        if dest.exists() and dest.stat().st_size > 0:
            logger.info("download_tcga_brca: %s already present, skipping", name)
            continue
        data = _fetch_bytes(f"{XENA_BASE}/{name}")
        with gzip.GzipFile(fileobj=__import__("io").BytesIO(data)) as gz:
            gz.read(64)  # fail fast on a non-gzip / HTML error body before writing
        dest.write_bytes(data)
        logger.info("download_tcga_brca: wrote %s (%d bytes)", dest, len(data))


def _download_cbioportal_clinical(study_id: str, clinical_data_type: str, out_path: Path) -> None:
    if out_path.exists() and out_path.stat().st_size > 0:
        logger.info("download_tcga_brca: %s already present, skipping", out_path.name)
        return
    records = _fetch_json(
        f"{CBIOPORTAL_BASE}/studies/{study_id}/clinical-data"
        f"?clinicalDataType={clinical_data_type}&projection=SUMMARY&pageSize=100000"
    )
    out_path.write_text(json.dumps(records))
    logger.info("download_tcga_brca: wrote %s (%d records)", out_path, len(records))


def _download_probemap(out_dir: Path) -> None:
    dest = out_dir / XENA_PROBEMAP
    if dest.exists() and dest.stat().st_size > 0:
        logger.info("download_tcga_brca: %s already present, skipping", XENA_PROBEMAP)
        return
    data = _fetch_bytes(f"{XENA_BASE}/{XENA_PROBEMAP}")
    if not data.startswith(b"id\tgene\t"):
        raise RuntimeError(f"download_tcga_brca: unexpected probemap content head: {data[:80]!r}")
    dest.write_bytes(data)
    logger.info("download_tcga_brca: wrote %s (%d bytes)", dest, len(data))


def run(out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    _download_xena(out_dir)
    _download_probemap(out_dir)
    # ER_STATUS / PAM50_SUBTYPE / HER2_STATUS / PR_STATUS / SAMPLE_TYPE /
    # TUMOR_STAGE are SAMPLE-level attributes in brca_tcga_pub, not
    # PATIENT-level, despite some overlapping demographic fields (AGE, SEX)
    # being patient-level -- both are pulled so nothing is silently missed.
    _download_cbioportal_clinical(CBIOPORTAL_RECEPTOR_STUDY, "PATIENT", out_dir / "cbioportal_brca_tcga_pub_clinical_patient.json")
    _download_cbioportal_clinical(CBIOPORTAL_RECEPTOR_STUDY, "SAMPLE", out_dir / "cbioportal_brca_tcga_pub_clinical_sample.json")
    _download_cbioportal_clinical(CBIOPORTAL_SURVIVAL_STUDY, "PATIENT", out_dir / "cbioportal_brca_tcga_pan_can_atlas_2018_clinical_patient.json")

    manifest = out_dir / "PROVENANCE.txt"
    manifest.write_text(
        "TCGA-BRCA independent-validation raw data provenance\n"
        "======================================================\n\n"
        "1. Expression: GDC-harmonized STAR-pipeline TPM, log2(TPM+1),\n"
        "   via UCSC Xena GDC hub.\n"
        "   https://gdc.xenahubs.net/download/TCGA-BRCA.star_tpm.tsv.gz\n\n"
        "2. Survival (Xena, cross-check only): \n"
        "   https://gdc.xenahubs.net/download/TCGA-BRCA.survival.tsv.gz\n\n"
        "3. Clinical receptor status + PAM50 (PRIMARY): cBioPortal study\n"
        f"   {CBIOPORTAL_RECEPTOR_STUDY} (TCGA, Nature 2012 discovery cohort).\n"
        f"   {CBIOPORTAL_BASE}/studies/{CBIOPORTAL_RECEPTOR_STUDY}/clinical-data\n\n"
        "4. Survival + stage/age (PRIMARY for Cox models): cBioPortal study\n"
        f"   {CBIOPORTAL_SURVIVAL_STUDY} (TCGA PanCancer Atlas, Cell 2018).\n"
        f"   {CBIOPORTAL_BASE}/studies/{CBIOPORTAL_SURVIVAL_STUDY}/clinical-data\n"
    )
    logger.info("download_tcga_brca: wrote provenance manifest to %s", manifest)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", type=Path, default=Path("/ibex/scratch/aljaroaa/tamoxifen-data/tcga_brca"))
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO)
    run(args.out_dir)


if __name__ == "__main__":
    main()
