"""Downloads two real, published PDB structures for KDM1A/LSD1 and TLK2, for
the poster-exploration-v2 structural comparison figures (Section G).

Not part of the deterministic src/ pipeline -- one-time network download,
following the same pattern already established in
scripts/download/download_usp34_structures.py. Both PDB IDs were already
independently verified (web-verified 2026-08-15) and cited in the frozen
results/tables/post_audit_sensitivity/06b_structural_tractability_audit.tsv:

- 6NQU: human LSD1/KDM1A catalytic domain in complex with the inhibitor
  GSK2879552 (MDPI Molecules 23:1538 / PMC6099836).
- 5O0Y: human TLK2 kinase domain bound to AGS (ATP-gamma-S, a
  non-hydrolysable ATP analog, NOT a small-molecule inhibitor or drug)
  (Nature Communications 2018, s41467-018-04941-y).

This script performs no new structural analysis, pocket detection, or
docking -- it only fetches the public coordinate files so they can be
rendered as static exploratory figures.

Source: RCSB PDB (files.rcsb.org).
"""

from __future__ import annotations

import hashlib
import logging
from pathlib import Path

import requests

logger = logging.getLogger(__name__)

BASE_URL = "https://files.rcsb.org/download/{pdb_id}.pdb"
PDB_IDS = ["6NQU", "5O0Y"]


def download(out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    hashes = {}
    for pdb_id in PDB_IDS:
        fname = f"{pdb_id}.pdb"
        out_path = out_dir / fname
        resp = requests.get(BASE_URL.format(pdb_id=pdb_id), timeout=60)
        resp.raise_for_status()
        out_path.write_bytes(resp.content)
        actual = hashlib.sha256(resp.content).hexdigest()
        hashes[fname] = actual
        logger.info("wrote %s (%d bytes), sha256=%s", out_path, len(resp.content), actual)

    provenance = out_dir / "PROVENANCE.txt"
    provenance.write_text(
        "KDM1A (LSD1) and TLK2 PDB structures -- provenance\n"
        "Downloaded: 2026-08-16\n"
        "Source: https://files.rcsb.org/download/{PDB_ID}.pdb\n"
        "Purpose: static exploratory-figure rendering only "
        "(poster-exploration-v2 phase, Section G structural comparison).\n"
        "The structural facts themselves (KDM1A inhibitor-bound co-crystal, "
        "TLK2 kinase domain bound to an ATP analog not an inhibitor) were "
        "already established and independently web-verified in the "
        "post-audit sensitivity phase -- see "
        "results/tables/post_audit_sensitivity/06b_structural_tractability_audit.tsv. "
        "No new structural analysis, pocket detection, or docking was "
        "performed with these files.\n\n"
        "6NQU.pdb -- human LSD1/KDM1A catalytic domain, inhibitor-bound "
        "(GSK2879552), X-ray. PMC6099836.\n"
        f"  sha256={hashes.get('6NQU.pdb', 'unknown')}\n\n"
        "5O0Y.pdb -- human TLK2 kinase domain, bound to AGS (ATP-gamma-S, "
        "a non-hydrolysable ATP analog -- a substrate mimetic, NOT a "
        "small-molecule inhibitor or drug), 2.86 A, X-ray. "
        "Nature Communications 2018 (s41467-018-04941-y).\n"
        f"  sha256={hashes.get('5O0Y.pdb', 'unknown')}\n"
    )
    logger.info("wrote %s", provenance)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    import yaml

    cfg = yaml.safe_load(open("config/config.yaml"))
    download(Path(cfg["data"]["raw"]["kdm1a_tlk2_structures_dir"]))
