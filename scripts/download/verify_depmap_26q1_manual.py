"""Verify manually-downloaded DepMap Public 26Q1 files before they are
used in any analysis.

Run this AFTER manually downloading the 26Q1 release bundle from the
DepMap portal (see results/reports/independent_validation/
DEPMAP_26Q1_ACCESS_STATUS.md for exact instructions) and placing the
files at the path config/config.yaml's data.raw.depmap_26q1_dir points to
(default /ibex/scratch/aljaroaa/tamoxifen-data/depmap/26Q1).

Checks, for each required file: it exists, is non-empty, has the expected
header schema ("ModelID" for Model.csv; "SYMBOL (Entrez)" gene columns for
the three matrices), and that all four candidate genes' expected Entrez
IDs are present in each matrix. Does NOT invent or accept a checksum --
DepMap's own portal download does not expose one to a manual downloader,
so this script instead does a structural/content sanity check and prints
a locally-computed SHA256 for the user to record. Writes a PROVENANCE.txt
manifest on success. Performs no analysis.
"""

from __future__ import annotations

import argparse
import hashlib
import logging
import sys
from datetime import date
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from src.independent_validation_depmap_data import CANDIDATE_ENTREZ  # noqa: E402

logger = logging.getLogger(__name__)

REQUIRED_FILES = {
    "Model.csv": None,
    "CRISPRGeneEffect.csv": None,
    "CRISPRGeneDependency.csv": None,
    "OmicsExpressionProteinCodingGenesTPMLogp1.csv": None,
}


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _verify_model_csv(path: Path) -> str:
    df = pd.read_csv(path, nrows=5)
    if "ModelID" not in df.columns:
        raise ValueError(f"{path.name}: no ModelID column -- is this really Model.csv?")
    full = pd.read_csv(path)
    if full["ModelID"].duplicated().any():
        raise ValueError(f"{path.name}: duplicate ModelID rows")
    n_breast = int((full.get("OncotreeLineage") == "Breast").sum())
    return f"{len(full)} models, {n_breast} breast (OncotreeLineage=='Breast')"


def _verify_matrix_csv(path: Path) -> str:
    header = pd.read_csv(path, nrows=0).columns
    missing = [f"{sym} ({eid})" for sym, eid in CANDIDATE_ENTREZ.items() if f"{sym} ({eid})" not in header]
    if missing:
        raise ValueError(f"{path.name}: missing expected candidate gene column(s): {missing}")
    n_rows = sum(1 for _ in open(path)) - 1
    return f"{n_rows} model rows, {len(header) - 1} gene columns, all 4 candidate genes present"


def run(release_dir: Path) -> None:
    if not release_dir.exists():
        raise FileNotFoundError(f"{release_dir} does not exist -- create it and place the 4 files there first")

    results = {}
    missing = []
    for fname in REQUIRED_FILES:
        path = release_dir / fname
        if not path.exists():
            missing.append(fname)
            continue
        if fname == "Model.csv":
            summary = _verify_model_csv(path)
        else:
            summary = _verify_matrix_csv(path)
        sha = _sha256(path)
        size = path.stat().st_size
        results[fname] = (summary, sha, size)
        logger.info("verify_depmap_26q1_manual: %s OK -- %s (%d bytes, sha256 %s)", fname, summary, size, sha)

    if missing:
        raise FileNotFoundError(f"still missing: {missing} in {release_dir}")

    manifest_lines = [
        "DepMap Public 26Q1 -- manually-downloaded file verification",
        "=" * 60,
        f"Verified on: {date.today().isoformat()}",
        f"Directory: {release_dir}",
        "",
    ]
    for fname, (summary, sha, size) in results.items():
        manifest_lines += [f"{fname}", f"  size: {size} bytes", f"  sha256: {sha}", f"  content check: {summary}", ""]
    (release_dir / "PROVENANCE.txt").write_text("\n".join(manifest_lines))
    logger.info("verify_depmap_26q1_manual: ALL FILES PRESENT AND STRUCTURALLY VALID. Wrote %s", release_dir / "PROVENANCE.txt")
    logger.info("verify_depmap_26q1_manual: you can now set config/config.yaml's independent_validation.depmap.active_release to '26Q1' and rerun the DepMap analysis modules.")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--release-dir", type=Path, default=Path("/ibex/scratch/aljaroaa/tamoxifen-data/depmap/26Q1"))
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO)
    run(args.release_dir)


if __name__ == "__main__":
    main()
