#!/usr/bin/env bash
# Downloads MSigDB Hallmark, Reactome (C2 CP:REACTOME), and GO Biological
# Process (C5 GO:BP) gene sets (human, gene-symbol GMT format) for the
# systems-network pathway phase. Run once; downstream analysis modules read
# only the local files this script writes, never the network at runtime.
set -euo pipefail

OUT_DIR="data/reference/genesets"
mkdir -p "$OUT_DIR"

RELEASE="2024.1.Hs"
BASE_URL="https://data.broadinstitute.org/gsea-msigdb/msigdb/release/${RELEASE}"

declare -A FILES=(
  ["h.all.v${RELEASE}.symbols.gmt"]="hallmark"
  ["c2.cp.reactome.v${RELEASE}.symbols.gmt"]="reactome"
  ["c5.go.bp.v${RELEASE}.symbols.gmt"]="go_bp"
)

for fname in "${!FILES[@]}"; do
  label="${FILES[$fname]}"
  dest="${OUT_DIR}/${label}.gmt"
  echo "Downloading ${label} (${fname}) -> ${dest}"
  curl -sSf "${BASE_URL}/${fname}" -o "${dest}"
  n_sets=$(wc -l < "${dest}")
  echo "  ${n_sets} gene sets"
done

date -u +"%Y-%m-%dT%H:%M:%SZ" > "${OUT_DIR}/download_timestamp_utc.txt"
echo "MSigDB release ${RELEASE}" > "${OUT_DIR}/version.txt"
echo "Source: ${BASE_URL}" >> "${OUT_DIR}/version.txt"
echo "Done."
