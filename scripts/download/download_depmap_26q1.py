"""Attempt to download the DepMap Public 26Q1 release; document exactly
what is and is not obtainable via an official, non-interactive channel.

DepMap's announced 2026 release schedule is 26Q1 and 26Q3 (per
https://forum.depmap.org/t/depmap-quarterly-release-notes/3560, checked
2026-08-14) -- there is no 26Q2. 26Q1 was announced 2026-04-01
(https://forum.depmap.org/t/announcing-the-26q1-release/4606).

ACCESS ATTEMPTS (all checked 2026-08-14; see
results/reports/independent_validation/DEPMAP_26Q1_ACCESS_STATUS.md for
the full log):

1. DepMap portal bulk-download API (depmap.org/portal/api/download/all,
   depmap.org/portal/data_page/?tab=allData) -- returns a Cloudflare
   Turnstile "verifying you're a person" challenge page to every
   non-interactive client tested (curl and the WebFetch tool alike).
   BLOCKED, not bypassed.
2. Figshare/Figshare+ search API (api.figshare.com/v2/articles/search) --
   searched with >15 phrasings for "DepMap Public 26Q1", "DepMap 26Q1
   Public", "Model 26Q1", "CRISPRGeneEffect 26Q1", etc. The ONLY 26Q1 item
   the search API returns is Figshare+ article 31660582, "Chronos
   parameters (Public 26Q1)" (DOI 10.6084/m9.figshare.31660582.v2).
   CORRECTED after independent review (2026-08-14): this item's author is
   an individual account ("Yejie Yun", Figshare author id 21559277), NOT
   the institutional "Broad DepMap" account (id 17476659) that published
   the verified DepMap_24Q4_Public/24Q2/23Q4 bundles -- an earlier version
   of this docstring wrongly claimed a same-account match. This item's
   gene_effect.csv (413,544,937 bytes, MD5
   c89e3ec7e2c3682e5c3535172177a1ee) IS downloaded here for reference, but
   is DepMap-LABELED ONLY, not confirmed-official: its uploader identity
   is unverified and it does NOT resolve whether it is the batch-corrected
   or uncorrected release artifact. It is not used in any analysis. The
   same Figshare item does NOT include Model.csv, CRISPRGeneDependency.csv,
   or a matched expression file (only Chronos-internal nuisance parameters: t0_offset,
   guide_efficacy, replicate_efficacy, library_effect -- not downloaded,
   not used by this project).
3. Figshare author-listing page (figshare.com/authors/Broad_DepMap/...)
   -- returns an empty client-rendered SPA shell to a non-interactive
   client; cannot be scraped without executing the page's JavaScript.
4. Google Cloud Storage bucket referenced in the 26Q1 release notes
   (storage.googleapis.com/shared-portal-files/...) -- individual known
   object paths (e.g. the mutation-pipeline PDF cited in the release
   notes) ARE publicly readable, but the bucket does not permit anonymous
   LISTING, so a data-file object key cannot be discovered without
   already knowing it; no data-file path could be found in any published
   document.
5. AnVIL/Terra (anvilproject.org/news/2026/03/03/depmap-data-release) --
   requires NIH Researcher Authentication Service (RAS) login via
   Login.gov/ID.me and a dbGaP accession (phs003444.v3.p1) for
   controlled-access data; not usable non-interactively, and covers raw
   sequencing data rather than the processed public CRISPRGeneEffect/
   Model.csv files this project needs.

RESULT: CRISPR gene-effect data partially obtained (with the corrected/
uncorrected caveat above). Model.csv, CRISPRGeneDependency.csv, and the
matched expression file could NOT be obtained. Per instruction, this
script does NOT fall back to 24Q4, does NOT scrape past the Cloudflare
challenge, and does NOT fabricate the missing files. It downloads only
what is discoverable via a legitimate API without bypassing any access
control, and reports the rest as missing.

This script only downloads/verifies. It performs no analysis. Run it
manually; it is safe to re-run (skips files already present and
verified).
"""

from __future__ import annotations

import argparse
import hashlib
import logging
import urllib.request
from pathlib import Path

logger = logging.getLogger(__name__)

RELEASE = "26Q1"
CHRONOS_PARAMS_FIGSHARE_ARTICLE = 31660582
CHRONOS_PARAMS_DOI = "10.6084/m9.figshare.31660582.v2"

# Only the gene_effect.csv is relevant to this project; the other files in
# this Figshare item are Chronos-internal nuisance parameters, not used.
GENE_EFFECT_URL = "https://ndownloader.figshare.com/files/67214582"
GENE_EFFECT_EXPECTED_MD5 = "c89e3ec7e2c3682e5c3535172177a1ee"
GENE_EFFECT_EXPECTED_SIZE = 413544937

MISSING_FILES = [
    "Model.csv (cell-line metadata -- required to define breast/ER+/luminal cohorts)",
    "CRISPRGeneDependency.csv (per-line dependency probability)",
    "OmicsExpressionProteinCodingGenesTPMLogp1.csv (matched CCLE expression)",
    "The canonical CRISPRGeneEffect.csv itself (to resolve the corrected/uncorrected ambiguity -- see module docstring)",
]


def _md5(path: Path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _download(url: str, dest: Path) -> None:
    logger.info("download_depmap_26q1: fetching %s -> %s", url, dest)
    request = urllib.request.Request(url, headers={"User-Agent": "breast-cancer-project/1.0"})
    with urllib.request.urlopen(request, timeout=600) as response, open(dest, "wb") as out:
        while True:
            chunk = response.read(1 << 20)
            if not chunk:
                break
            out.write(chunk)


def run(out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    dest = out_dir / "gene_effect_chronos_params.csv"
    if not (dest.exists() and dest.stat().st_size == GENE_EFFECT_EXPECTED_SIZE and _md5(dest) == GENE_EFFECT_EXPECTED_MD5):
        _download(GENE_EFFECT_URL, dest)
        actual_size, actual_md5 = dest.stat().st_size, _md5(dest)
        if actual_size != GENE_EFFECT_EXPECTED_SIZE or actual_md5 != GENE_EFFECT_EXPECTED_MD5:
            raise RuntimeError(f"download_depmap_26q1: gene_effect.csv failed verification (size {actual_size}, md5 {actual_md5})")
        logger.info("download_depmap_26q1: gene_effect.csv verified (%d bytes, md5 %s)", actual_size, actual_md5)
    else:
        logger.info("download_depmap_26q1: gene_effect.csv already present and verified, skipping")

    logger.warning(
        "download_depmap_26q1: BLOCKED -- the following required files could NOT be obtained "
        "via any official non-interactive channel: %s. See results/reports/independent_validation/"
        "DEPMAP_26Q1_ACCESS_STATUS.md for manual-download instructions. No analysis will run against "
        "26Q1 until these are placed at %s.",
        "; ".join(MISSING_FILES), out_dir,
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", type=Path, default=Path("/ibex/scratch/aljaroaa/tamoxifen-data/depmap/26Q1"))
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO)
    run(args.out_dir)


if __name__ == "__main__":
    main()
