"""Download the DepMap Public CRISPR-dependency and Omics-expression release.

Source: Broad Institute Cancer Dependency Map (DepMap) Public release,
distributed via the Figshare+ mirror (the DepMap portal's own download API
sits behind a Cloudflare Turnstile browser challenge that a non-interactive
client cannot pass; Figshare+ is DepMap's own official secondary
distribution channel for the same release bundle).

Release actually accessible during this run: DepMap Public 24Q4
(published 2024-12-10), Figshare+ article 27993248,
DOI 10.25452/figshare.plus.27993248.v1,
https://plus.figshare.com/articles/dataset/DepMap_24Q4_Public/27993248

Newer releases (25Q2, 25Q3, 26Q1) were confirmed to exist via the DepMap
community forum's quarterly release-notes thread
(https://forum.depmap.org/t/depmap-quarterly-release-notes/3560) but no
complete, programmatically-downloadable data bundle for any of them could
be located during this run: the DepMap portal's own download endpoint
(depmap.org/portal/api/download/all) returns a Cloudflare Turnstile
verification page rather than data, and no Figshare+/figshare mirror
article for 25Q2, 25Q3, or 26Q1 could be found via the Figshare search API
or web search. This script therefore pins 24Q4 rather than guessing a
release identifier for data that could not actually be fetched.

Downloads only the four files this project's DepMap analysis needs:
Model.csv (cell-line metadata), CRISPRGeneEffect.csv (Chronos gene-effect
scores, batch-corrected -- the primary dependency metric),
CRISPRGeneDependency.csv (per-gene-per-line dependency probability, used
for the "fraction of strongly dependent lines" threshold), and
OmicsExpressionProteinCodingGenesTPMLogp1.csv (matched CCLE expression).
Each download is verified against the file-level MD5 Figshare reports for
that release.

This script only downloads. It performs no analysis and is never imported
by src/ or tests/ -- run it manually, once, to populate the local raw-data
cache referenced by config/config.yaml's depmap.raw section.
"""

from __future__ import annotations

import argparse
import hashlib
import logging
import urllib.request
from pathlib import Path

logger = logging.getLogger(__name__)

DEPMAP_RELEASE = "24Q4"
FIGSHARE_ARTICLE_ID = 27993248
FIGSHARE_DOI = "10.25452/figshare.plus.27993248.v1"
ARTICLE_PUBLISHED_DATE = "2024-12-10"

# file name -> (figshare file-level download URL, expected MD5, byte size)
FILES = {
    "Model.csv": ("https://ndownloader.figshare.com/files/51065297", "675210d17675f3517b0ce39a3c274f16", 645696),
    "CRISPRGeneEffect.csv": ("https://ndownloader.figshare.com/files/51064667", "6edf7ade09b9b34199210b559d4745d3", 428678699),
    "CRISPRGeneDependency.csv": ("https://ndownloader.figshare.com/files/51064631", "0d3bdadf0c59264e39f7fbadf232ccdb", 421115594),
    "OmicsExpressionProteinCodingGenesTPMLogp1.csv": ("https://ndownloader.figshare.com/files/51065489", "71794802b750ce77c422dad0720a40af", 506628654),
    "README.txt": ("https://ndownloader.figshare.com/files/51065795", "54628c15a10b9ff6536db245cfed5231", 43103),
}


def _md5(path: Path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _download(url: str, dest: Path) -> None:
    logger.info("download_depmap: fetching %s -> %s", url, dest)
    request = urllib.request.Request(url, headers={"User-Agent": "breast-cancer-project/1.0"})
    with urllib.request.urlopen(request, timeout=600) as response, open(dest, "wb") as out:
        while True:
            chunk = response.read(1 << 20)
            if not chunk:
                break
            out.write(chunk)


def run(out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    for name, (url, expected_md5, expected_size) in FILES.items():
        dest = out_dir / name
        if dest.exists() and dest.stat().st_size == expected_size and _md5(dest) == expected_md5:
            logger.info("download_depmap: %s already present and verified, skipping", name)
            continue
        _download(url, dest)
        actual_size = dest.stat().st_size
        actual_md5 = _md5(dest)
        if actual_size != expected_size or actual_md5 != expected_md5:
            raise RuntimeError(
                f"download_depmap: {name} failed verification "
                f"(size {actual_size} vs expected {expected_size}, "
                f"md5 {actual_md5} vs expected {expected_md5})"
            )
        logger.info("download_depmap: %s verified (%d bytes, md5 %s)", name, actual_size, actual_md5)

    manifest = out_dir / "PROVENANCE.txt"
    manifest.write_text(
        f"DepMap release: {DEPMAP_RELEASE}\n"
        f"Figshare+ article: {FIGSHARE_ARTICLE_ID}\n"
        f"DOI: {FIGSHARE_DOI}\n"
        f"Article published date: {ARTICLE_PUBLISHED_DATE}\n"
        f"Source URL: https://plus.figshare.com/articles/dataset/DepMap_24Q4_Public/{FIGSHARE_ARTICLE_ID}\n"
        f"Files: {', '.join(FILES)}\n"
    )
    logger.info("download_depmap: wrote provenance manifest to %s", manifest)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", type=Path, default=Path("/ibex/scratch/aljaroaa/tamoxifen-data/depmap/24Q4"))
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO)
    run(args.out_dir)


if __name__ == "__main__":
    main()
