"""Download the Hart et al. 2017 Core Essential Genes 2.0 (CEG2) reference list.

Source: Hart T et al. "Evaluation and Design of Genome-Wide CRISPR/SpCas9
Knockout Screens." G3 (Bethesda). 2017;7:2719-2727.
DOI: 10.1534/g3.117.041277. PMCID: PMC5555476.

Fetches the official journal supplementary-material archive (mirrored by
Europe PMC, which serves the same publisher-deposited files as
academic.oup.com/PMC without the bot-detection challenges that block a
non-interactive client on those hosts), verifies it is a well-formed ZIP,
extracts Table S2 (the 684-gene CEG2 list: HGNC gene symbol and HGNC ID,
one gene per row, no header), and writes it to a local reference file.

The outer Europe PMC archive is re-stamped with a fresh internal timestamp
on every request, so its own SHA256 is not a usable reproducibility check
(see config/resources.yaml's archive_sha256_note). Reproducibility is
instead enforced against two checksums pinned in config/resources.yaml
that verified stable across independent fetches: the extracted
2719TableS2.txt member itself (archive_member_sha256) and the final
formatted reference file (extracted_reference_sha256). A re-run whose
extracted content, gene count, or either checksum drifts from those pinned
values fails loudly rather than silently overwriting the reference file.

This script only downloads and extracts. It performs no analysis and is
never imported by src/ or tests/, and is not run automatically by any
test or downstream analysis step -- run it manually, once, to populate
the local reference file.
"""

from __future__ import annotations

import argparse
import hashlib
import io
import logging
import urllib.request
import zipfile
from pathlib import Path

import yaml

logger = logging.getLogger(__name__)

# Pinned source: Europe PMC's supplementary-files bundle for PMC5555476,
# i.e. the publisher-deposited supplementary archive for this article.
SOURCE_URL = "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC5555476/supplementaryFiles"
ARCHIVE_MEMBER = "2719TableS2.txt"
EXPECTED_GENE_COUNT = 684
CEG2_RESOURCE_DOI = "10.1534/g3.117.041277"


def _sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _download(url: str) -> bytes:
    logger.info("download_ceg2: fetching %s", url)
    request = urllib.request.Request(url, headers={"User-Agent": "breast-cancer-project/1.0"})
    with urllib.request.urlopen(request, timeout=60) as response:
        data = response.read()
    logger.info("download_ceg2: downloaded %d bytes", len(data))
    return data


def _extract_ceg2_table(archive_bytes: bytes) -> bytes:
    """Pull the CEG2 gene-list member out of the supplementary-archive ZIP."""
    with zipfile.ZipFile(io.BytesIO(archive_bytes)) as archive:
        names = archive.namelist()
        if ARCHIVE_MEMBER not in names:
            raise ValueError(
                f"expected member {ARCHIVE_MEMBER!r} not found in supplementary archive; "
                f"archive contains: {sorted(names)}"
            )
        return archive.read(ARCHIVE_MEMBER)


def _validate_gene_list(raw: bytes) -> list[tuple[str, str]]:
    """Parse and validate the CEG2 table: two tab-separated columns, no header."""
    text = raw.decode("ascii")
    lines = [line for line in text.splitlines() if line != ""]
    blank_lines_dropped = text.count("\n") + 1 - len(lines) if text else 0
    logger.info(
        "download_ceg2: %d non-blank rows read from %s (%d blank lines dropped)",
        len(lines),
        ARCHIVE_MEMBER,
        blank_lines_dropped,
    )

    rows: list[tuple[str, str]] = []
    for line in lines:
        fields = line.split("\t")
        if len(fields) != 2:
            raise ValueError(f"expected 2 tab-separated fields, got {len(fields)}: {line!r}")
        symbol, hgnc_id = fields
        if not symbol.strip() or not hgnc_id.strip():
            raise ValueError(f"blank gene symbol or HGNC id in row: {line!r}")
        rows.append((symbol.strip(), hgnc_id.strip()))

    symbols = [symbol for symbol, _ in rows]
    if len(set(symbols)) != len(symbols):
        duplicates = sorted({s for s in symbols if symbols.count(s) > 1})
        raise ValueError(f"duplicate CEG2 gene symbols: {duplicates}")

    if len(rows) != EXPECTED_GENE_COUNT:
        raise ValueError(
            f"expected exactly {EXPECTED_GENE_COUNT} CEG2 genes, got {len(rows)}"
        )

    logger.info("download_ceg2: validated %d unique CEG2 gene symbols", len(rows))
    return rows


def _load_config(config_path: str | Path) -> dict:
    with open(config_path) as f:
        return yaml.safe_load(f)


def _find_ceg2_resource(resources: dict) -> dict:
    matches = [s for s in resources["sources"] if s.get("accession_or_doi") == CEG2_RESOURCE_DOI]
    if len(matches) != 1:
        raise ValueError(
            f"expected exactly one resources.yaml entry with accession_or_doi "
            f"{CEG2_RESOURCE_DOI!r}, found {len(matches)}"
        )
    return matches[0]


def _load_expected_checksums(resources_path: str | Path) -> tuple[str, str]:
    """Read the pinned, stable checksums for the CEG2 archive member and the
    final extracted reference file from config/resources.yaml.

    These are the reproducibility anchors: unlike the outer Europe PMC zip
    (re-stamped per request), both were verified stable across independent
    fetches and are enforced on every run of this script.
    """
    resources = _load_config(resources_path)
    entry = _find_ceg2_resource(resources)
    try:
        archive_member_sha256 = entry["archive_member_sha256"]
        extracted_reference_sha256 = entry["extracted_reference_sha256"]
    except KeyError as exc:
        raise ValueError(
            f"resources.yaml CEG2 entry missing expected checksum field {exc}; "
            "cannot verify reproducibility without a pinned checksum"
        ) from exc
    return archive_member_sha256, extracted_reference_sha256


def _write_bytes_atomically(path: Path, data: bytes) -> None:
    """Write via a same-directory temp file + atomic rename, so a previously
    good file at ``path`` is never left partially overwritten or clobbered
    by content that turns out to fail validation."""
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = path.with_name(path.name + ".tmp")
    tmp_path.write_bytes(data)
    tmp_path.replace(path)


def download_ceg2(
    config_path: str | Path = "config/config.yaml",
    resources_path: str | Path = "config/resources.yaml",
) -> Path:
    config = _load_config(config_path)
    archive_path = Path(config["data"]["raw"]["ceg2_supplementary_archive"])
    reference_path = Path(config["data"]["reference"]["ceg2_gene_list"])

    expected_member_sha256, expected_reference_sha256 = _load_expected_checksums(resources_path)

    archive_bytes = _download(SOURCE_URL)

    # Extract and verify the member entirely in memory before writing
    # anything to disk: a failed download or a changed member must never
    # clobber a previously good archive or reference file on disk.
    table_bytes = _extract_ceg2_table(archive_bytes)
    table_sha256 = _sha256(table_bytes)
    if table_sha256 != expected_member_sha256:
        raise ValueError(
            f"{ARCHIVE_MEMBER} content changed: sha256={table_sha256} does not match "
            f"pinned resources.yaml archive_member_sha256={expected_member_sha256}"
        )
    logger.info(
        "download_ceg2: %s content verified against pinned checksum (sha256=%s)",
        ARCHIVE_MEMBER,
        table_sha256,
    )

    rows = _validate_gene_list(table_bytes)

    reference_text = "gene_symbol\thgnc_id\n" + "".join(
        f"{symbol}\t{hgnc_id}\n" for symbol, hgnc_id in rows
    )
    reference_bytes = reference_text.encode("utf-8")
    reference_sha256 = _sha256(reference_bytes)
    if reference_sha256 != expected_reference_sha256:
        raise ValueError(
            f"extracted CEG2 reference content changed: sha256={reference_sha256} does not "
            f"match pinned resources.yaml extracted_reference_sha256={expected_reference_sha256}"
        )

    # Only write to disk once both pinned checksums have verified: the
    # member content and the derived reference content are both known good.
    archive_sha256 = _sha256(archive_bytes)
    _write_bytes_atomically(archive_path, archive_bytes)
    logger.info(
        "download_ceg2: wrote raw supplementary archive to %s (sha256=%s, not verified -- "
        "the outer archive is re-stamped per request, see resources.yaml)",
        archive_path,
        archive_sha256,
    )

    _write_bytes_atomically(reference_path, reference_bytes)
    logger.info(
        "download_ceg2: wrote %d CEG2 genes to %s, verified against pinned checksum (sha256=%s)",
        len(rows),
        reference_path,
        reference_sha256,
    )
    logger.info(
        "download_ceg2: provenance -- source_url=%s archive_member=%s",
        SOURCE_URL,
        ARCHIVE_MEMBER,
    )
    return reference_path


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="config/config.yaml")
    parser.add_argument("--resources", default="config/resources.yaml")
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO)
    download_ceg2(args.config, args.resources)
