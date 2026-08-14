"""Downloads GDSC / CancerRxGene Release 8.5 (30 Oct 2023) raw files.

Not part of the deterministic src/ pipeline -- this script performs the
one-time network download documented in config.yaml's
final_pharmacogenomics.gdsc section. Run manually; src/ modules only ever
read the already-downloaded local files this script produces.

Files (all hosted at https://cog.sanger.ac.uk/cancerrxgene/GDSC_release8.5/,
confirmed live and ungated, unlike the DepMap portal):
- GDSC1_fitted_dose_response_27Oct23.xlsx
- GDSC2_fitted_dose_response_27Oct23.xlsx
- screened_compounds_rel_8.5.csv
- Cell_Lines_Details.xlsx

Response-metric definitions (LN_IC50, AUC) were independently confirmed
against the official GDSC_Fitted_Data_Description.pdf (Sanger, v1.0.0,
21 Sep 2017) before any analysis was built on them -- see
results/reports/final_pharmacogenomics/USP34_VEZF1_GDSC_review.md Part 2.
"""

from __future__ import annotations

import hashlib
import logging
from pathlib import Path

import requests

logger = logging.getLogger(__name__)

BASE_URL = "https://cog.sanger.ac.uk/cancerrxgene/GDSC_release8.5/"
FILES = [
    "GDSC1_fitted_dose_response_27Oct23.xlsx",
    "GDSC2_fitted_dose_response_27Oct23.xlsx",
    "screened_compounds_rel_8.5.csv",
    "Cell_Lines_Details.xlsx",
]

# SHA256 hashes recorded at download time (2026-08-14), for reproducibility
# verification -- if a re-download ever produces a different hash, the
# release has changed upstream and this analysis must be re-examined.
EXPECTED_SHA256 = {
    "GDSC1_fitted_dose_response_27Oct23.xlsx": "837b0686500fde75179e490de08f034abd9f882d8b0253d637bafe83e156dafd",
    "GDSC2_fitted_dose_response_27Oct23.xlsx": "f950a7027be265f8a7a74220a27fd18cbd368485349bd8c2048e88bb1cd07560",
    "screened_compounds_rel_8.5.csv": "d3cbb8b595980b36bf92d54af1da1dd995d2e2ce624c7cc70aacf2a754dce782",
    "Cell_Lines_Details.xlsx": "5520e779d25469f864de9cc5d42b99bf9632bd346607c0ee6a9ebad7ece17ac3",
}


def download(out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    for fname in FILES:
        out_path = out_dir / fname
        resp = requests.get(BASE_URL + fname, timeout=120)
        resp.raise_for_status()
        out_path.write_bytes(resp.content)
        actual = hashlib.sha256(resp.content).hexdigest()
        expected = EXPECTED_SHA256[fname]
        status = "OK" if actual == expected else "MISMATCH -- release may have changed upstream"
        logger.info("wrote %s (%d bytes), sha256=%s [%s]", out_path, len(resp.content), actual, status)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    download(Path("/ibex/scratch/aljaroaa/tamoxifen-data/gdsc"))
