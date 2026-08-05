import hashlib
import io
import zipfile
from pathlib import Path

import pytest
import yaml

from scripts.download import download_ceg2 as download_ceg2_module
from scripts.download.download_ceg2 import (
    ARCHIVE_MEMBER,
    CEG2_RESOURCE_DOI,
    EXPECTED_GENE_COUNT,
    _extract_ceg2_table,
    _find_ceg2_resource,
    _load_expected_checksums,
    _sha256,
    _validate_gene_list,
    download_ceg2,
)


def _make_ceg2_bytes(n_genes: int = EXPECTED_GENE_COUNT, duplicate: bool = False) -> bytes:
    lines = [f"GENE{i}\tHGNC:{i}" for i in range(n_genes)]
    if duplicate and n_genes >= 2:
        lines[1] = lines[0]
    return ("\n".join(lines) + "\n").encode("ascii")


def _make_archive(member_name: str, content: bytes) -> bytes:
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w") as zf:
        zf.writestr(member_name, content)
    return buf.getvalue()


def _expected_reference_bytes(table_bytes: bytes) -> bytes:
    rows = _validate_gene_list(table_bytes)
    text = "gene_symbol\thgnc_id\n" + "".join(f"{s}\t{h}\n" for s, h in rows)
    return text.encode("utf-8")


def test_extract_ceg2_table_reads_expected_member():
    content = _make_ceg2_bytes()
    archive = _make_archive(ARCHIVE_MEMBER, content)
    extracted = _extract_ceg2_table(archive)
    assert extracted == content


def test_extract_ceg2_table_missing_member_raises():
    archive = _make_archive("SomeOtherFile.txt", b"irrelevant")
    with pytest.raises(ValueError, match="not found"):
        _extract_ceg2_table(archive)


def test_validate_gene_list_accepts_well_formed_684():
    rows = _validate_gene_list(_make_ceg2_bytes())
    assert len(rows) == EXPECTED_GENE_COUNT
    assert rows[0] == ("GENE0", "HGNC:0")


def test_validate_gene_list_rejects_wrong_count():
    with pytest.raises(ValueError, match=str(EXPECTED_GENE_COUNT)):
        _validate_gene_list(_make_ceg2_bytes(n_genes=683))


def test_validate_gene_list_rejects_duplicate_symbols():
    with pytest.raises(ValueError, match="duplicate"):
        _validate_gene_list(_make_ceg2_bytes(duplicate=True))


def test_validate_gene_list_rejects_malformed_row():
    bad = b"ONLY_ONE_COLUMN\n" + _make_ceg2_bytes(n_genes=EXPECTED_GENE_COUNT - 1)
    with pytest.raises(ValueError, match="2 tab-separated fields"):
        _validate_gene_list(bad)


# --- resources.yaml checksum lookup --------------------------------------


def _write_resources_yaml(path: Path, **entry_overrides) -> Path:
    entry = {
        "name": "Hart et al. 2017, Table S2 (test)",
        "accession_or_doi": CEG2_RESOURCE_DOI,
        "archive_member_sha256": "member-hash",
        "extracted_reference_sha256": "reference-hash",
    }
    entry.update(entry_overrides)
    path.write_text(yaml.safe_dump({"sources": [entry]}))
    return path


def test_load_expected_checksums_reads_pinned_values(tmp_path: Path):
    resources_path = _write_resources_yaml(
        tmp_path / "resources.yaml",
        archive_member_sha256="abc123",
        extracted_reference_sha256="def456",
    )
    member_sha, reference_sha = _load_expected_checksums(resources_path)
    assert member_sha == "abc123"
    assert reference_sha == "def456"


def test_find_ceg2_resource_rejects_missing_entry(tmp_path: Path):
    path = tmp_path / "resources.yaml"
    path.write_text(yaml.safe_dump({"sources": [{"accession_or_doi": "not-the-doi"}]}))
    with pytest.raises(ValueError, match="found 0"):
        _load_expected_checksums(path)


def test_find_ceg2_resource_rejects_duplicate_entries(tmp_path: Path):
    path = tmp_path / "resources.yaml"
    entry = {"accession_or_doi": CEG2_RESOURCE_DOI}
    path.write_text(yaml.safe_dump({"sources": [entry, dict(entry)]}))
    with pytest.raises(ValueError, match="found 2"):
        _find_ceg2_resource(yaml.safe_load(path.read_text()))


def test_load_expected_checksums_rejects_missing_checksum_field(tmp_path: Path):
    path = tmp_path / "resources.yaml"
    path.write_text(
        yaml.safe_dump(
            {"sources": [{"accession_or_doi": CEG2_RESOURCE_DOI, "archive_member_sha256": "x"}]}
        )
    )
    with pytest.raises(ValueError, match="missing expected checksum field"):
        _load_expected_checksums(path)


# --- end-to-end checksum enforcement in download_ceg2() ------------------


def _write_config_yaml(path: Path, archive_path: Path, reference_path: Path) -> Path:
    path.write_text(
        yaml.safe_dump(
            {
                "data": {
                    "raw": {"ceg2_supplementary_archive": str(archive_path)},
                    "reference": {"ceg2_gene_list": str(reference_path)},
                }
            }
        )
    )
    return path


def test_download_ceg2_succeeds_when_checksums_match(tmp_path: Path, monkeypatch):
    table_bytes = _make_ceg2_bytes()
    archive_bytes = _make_archive(ARCHIVE_MEMBER, table_bytes)
    reference_bytes = _expected_reference_bytes(table_bytes)
    archive_path = tmp_path / "archive.zip"
    reference_path = tmp_path / "reference.tsv"

    config_path = _write_config_yaml(tmp_path / "config.yaml", archive_path, reference_path)
    resources_path = _write_resources_yaml(
        tmp_path / "resources.yaml",
        archive_member_sha256=_sha256(table_bytes),
        extracted_reference_sha256=_sha256(reference_bytes),
    )
    monkeypatch.setattr(download_ceg2_module, "_download", lambda url: archive_bytes)

    result_path = download_ceg2(config_path, resources_path)
    assert result_path.read_bytes() == reference_bytes
    assert archive_path.read_bytes() == archive_bytes
    # No leftover .tmp files from the atomic-write helper.
    assert sorted(p.name for p in tmp_path.glob("*.tmp")) == []


def _seed_sentinel_files(archive_path: Path, reference_path: Path) -> tuple[bytes, bytes]:
    """Pre-populate archive_path/reference_path with distinct sentinel bytes,
    so a failure test can prove they are left byte-for-byte unchanged, not
    merely that a from-empty run didn't create them."""
    archive_sentinel = b"PREVIOUS_GOOD_ARCHIVE_SENTINEL"
    reference_sentinel = b"PREVIOUS_GOOD_REFERENCE_SENTINEL"
    archive_path.parent.mkdir(parents=True, exist_ok=True)
    reference_path.parent.mkdir(parents=True, exist_ok=True)
    archive_path.write_bytes(archive_sentinel)
    reference_path.write_bytes(reference_sentinel)
    return archive_sentinel, reference_sentinel


def test_download_ceg2_fails_loudly_on_changed_member_content(tmp_path: Path, monkeypatch):
    table_bytes = _make_ceg2_bytes()
    archive_bytes = _make_archive(ARCHIVE_MEMBER, table_bytes)
    reference_bytes = _expected_reference_bytes(table_bytes)
    archive_path = tmp_path / "archive.zip"
    reference_path = tmp_path / "reference.tsv"
    archive_sentinel, reference_sentinel = _seed_sentinel_files(archive_path, reference_path)

    config_path = _write_config_yaml(tmp_path / "config.yaml", archive_path, reference_path)
    resources_path = _write_resources_yaml(
        tmp_path / "resources.yaml",
        archive_member_sha256="0" * 64,  # deliberately wrong
        extracted_reference_sha256=_sha256(reference_bytes),
    )
    monkeypatch.setattr(download_ceg2_module, "_download", lambda url: archive_bytes)

    with pytest.raises(ValueError, match="content changed"):
        download_ceg2(config_path, resources_path)
    assert archive_path.read_bytes() == archive_sentinel
    assert reference_path.read_bytes() == reference_sentinel


def test_download_ceg2_fails_loudly_on_changed_reference_checksum(tmp_path: Path, monkeypatch):
    """Even if the pinned member checksum matches, a mismatched pinned
    reference checksum must still block the write (defense in depth: the
    formatting/serialization step could drift independently)."""
    table_bytes = _make_ceg2_bytes()
    archive_bytes = _make_archive(ARCHIVE_MEMBER, table_bytes)
    archive_path = tmp_path / "archive.zip"
    reference_path = tmp_path / "reference.tsv"
    archive_sentinel, reference_sentinel = _seed_sentinel_files(archive_path, reference_path)

    config_path = _write_config_yaml(tmp_path / "config.yaml", archive_path, reference_path)
    resources_path = _write_resources_yaml(
        tmp_path / "resources.yaml",
        archive_member_sha256=_sha256(table_bytes),
        extracted_reference_sha256="0" * 64,  # deliberately wrong
    )
    monkeypatch.setattr(download_ceg2_module, "_download", lambda url: archive_bytes)

    with pytest.raises(ValueError, match="extracted CEG2 reference content changed"):
        download_ceg2(config_path, resources_path)
    assert archive_path.read_bytes() == archive_sentinel
    assert reference_path.read_bytes() == reference_sentinel


def test_download_ceg2_fails_loudly_on_wrong_gene_count(tmp_path: Path, monkeypatch):
    table_bytes = _make_ceg2_bytes(n_genes=683)
    archive_bytes = _make_archive(ARCHIVE_MEMBER, table_bytes)
    archive_path = tmp_path / "archive.zip"
    reference_path = tmp_path / "reference.tsv"
    archive_sentinel, reference_sentinel = _seed_sentinel_files(archive_path, reference_path)

    config_path = _write_config_yaml(tmp_path / "config.yaml", archive_path, reference_path)
    # Pin the member checksum to this (wrong-count) table so the run gets
    # past the member-checksum check and fails on gene-count validation.
    resources_path = _write_resources_yaml(
        tmp_path / "resources.yaml",
        archive_member_sha256=_sha256(table_bytes),
        extracted_reference_sha256="0" * 64,
    )
    monkeypatch.setattr(download_ceg2_module, "_download", lambda url: archive_bytes)

    with pytest.raises(ValueError, match=str(EXPECTED_GENE_COUNT)):
        download_ceg2(config_path, resources_path)
    assert archive_path.read_bytes() == archive_sentinel
    assert reference_path.read_bytes() == reference_sentinel


def test_download_ceg2_fails_loudly_when_member_missing_from_archive(tmp_path: Path, monkeypatch):
    archive_bytes = _make_archive("WrongMember.txt", b"irrelevant")
    archive_path = tmp_path / "archive.zip"
    reference_path = tmp_path / "reference.tsv"
    archive_sentinel, reference_sentinel = _seed_sentinel_files(archive_path, reference_path)

    config_path = _write_config_yaml(tmp_path / "config.yaml", archive_path, reference_path)
    resources_path = _write_resources_yaml(tmp_path / "resources.yaml")
    monkeypatch.setattr(download_ceg2_module, "_download", lambda url: archive_bytes)

    with pytest.raises(ValueError, match="not found"):
        download_ceg2(config_path, resources_path)
    assert archive_path.read_bytes() == archive_sentinel
    assert reference_path.read_bytes() == reference_sentinel


def test_real_resources_yaml_has_pinned_ceg2_checksums():
    """The actual project config must carry both stable checksums, not just
    the (deliberately unstable, and therefore unenforced) outer archive
    checksum."""
    repo_root = Path(__file__).parent.parent
    member_sha, reference_sha = _load_expected_checksums(repo_root / "config" / "resources.yaml")
    assert len(member_sha) == 64
    assert len(reference_sha) == 64


def test_real_ceg2_reference_file_matches_pinned_checksum():
    repo_root = Path(__file__).parent.parent
    _, expected_reference_sha256 = _load_expected_checksums(repo_root / "config" / "resources.yaml")
    reference_path = repo_root / "data" / "reference" / "hart2017_ceg2_684.tsv"
    actual = hashlib.sha256(reference_path.read_bytes()).hexdigest()
    assert actual == expected_reference_sha256
