"""Downloads the two real, already-verified USP34 PDB structures for the
poster's Figure 5 (structural targetability) rendering.

Not part of the deterministic src/ pipeline -- one-time network download.
7W3R and 7W3U were already established as the project's frozen structural
evidence in the final_translational phase (see
results/tables/final_translational/USP34_structure_inventory.tsv and
USP34_pocket_analysis.tsv); this script only re-fetches the same public PDB
coordinate files so they can be rendered as a static poster image. It
performs no new structural analysis, pocket detection, or docking.

Source: RCSB PDB (files.rcsb.org), original structure paper PMID 35588869.
"""

from __future__ import annotations

import hashlib
import logging
from pathlib import Path

import requests

logger = logging.getLogger(__name__)

BASE_URL = "https://files.rcsb.org/download/{pdb_id}.pdb"
PDB_IDS = ["7W3R", "7W3U"]

# SHA256 hashes recorded at download time (2026-08-15), for reproducibility
# verification -- PDB coordinate files for a released structure do not
# change, so a mismatch on re-download would indicate a wrong ID or a
# corrupted transfer, not a legitimate upstream update.
EXPECTED_SHA256 = {
    "7W3R.pdb": "838286cb844ffafe108fbf077ca95fe8192f5f2c2969726f75fc99868ecc5769",
    "7W3U.pdb": "deb96ecda82df99242bdded39b24107de4ddd07fd25ac6783eae4b1f8a21d4c3",
}


def download(out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    for pdb_id in PDB_IDS:
        fname = f"{pdb_id}.pdb"
        out_path = out_dir / fname
        resp = requests.get(BASE_URL.format(pdb_id=pdb_id), timeout=60)
        resp.raise_for_status()
        out_path.write_bytes(resp.content)
        actual = hashlib.sha256(resp.content).hexdigest()
        expected = EXPECTED_SHA256.get(fname)
        status = "recorded this run" if expected is None else ("OK" if actual == expected else "MISMATCH")
        logger.info("wrote %s (%d bytes), sha256=%s [%s]", out_path, len(resp.content), actual, status)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    import yaml

    cfg = yaml.safe_load(open("config/config.yaml"))
    download(Path(cfg["data"]["raw"]["usp34_structures_dir"]))
