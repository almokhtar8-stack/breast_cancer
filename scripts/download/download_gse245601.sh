#!/usr/bin/env bash
# Download the public GEO supplementary archive for GSE245601 (Kim, Whitman
# et al., Clin Cancer Res 2023, PMID 37747807) and extract per-sample Cell
# Ranger filtered_feature_bc_matrix H5 files. No FASTQs are downloaded --
# raw reads are dbGaP-controlled (phs003186.v1.p1) and are not required for
# the gene-expression pseudobulk validation this project performs.
#
# Usage: scripts/download/download_gse245601.sh
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
DEST_DIR="${REPO_ROOT}/data/raw/gse245601"
URL="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE245nnn/GSE245601/suppl/GSE245601_RAW.tar"

mkdir -p "${DEST_DIR}/h5"

if [ ! -f "${DEST_DIR}/GSE245601_RAW.tar" ]; then
  echo "Downloading ${URL} ..."
  curl -sS -o "${DEST_DIR}/GSE245601_RAW.tar" "${URL}"
fi

echo "Extracting per-sample H5 files ..."
tar -xf "${DEST_DIR}/GSE245601_RAW.tar" -C "${DEST_DIR}/h5"

echo "Writing checksums ..."
(cd "${DEST_DIR}" && sha256sum h5/*.h5 > checksums.sha256)

n_files=$(ls "${DEST_DIR}"/h5/*.h5 | wc -l)
echo "Done. ${n_files} H5 files in ${DEST_DIR}/h5 (expected: 26)."
