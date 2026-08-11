"""GSE245601 sample manifest construction and integrity verification.

Data source: GEO GSE245601 (Kim, Whitman et al., Clin Cancer Res
2023;29(23):4894-4907, PMID 37747807), public supplementary archive
GSE245601_RAW.tar (26 Cell Ranger filtered_feature_bc_matrix H5 files, one
per GSM). Sample/patient/condition mapping is read from config.yaml
(``gse245601.samples``), which was constructed by independently
cross-verifying three primary sources (SRA BioSample XML attributes, GEO
GSM records, and the GEO-provided H5 filenames) -- see
docs/gse245601_PREANALYSIS.md. No FASTQs, no expression values, and no
recomputation of any kind are involved here; this module only verifies file
integrity and assembles the sample-level metadata table.
"""

from __future__ import annotations

import hashlib
from pathlib import Path

import pandas as pd
import yaml


def _load_config(config_path: str | Path) -> dict:
    with open(config_path) as f:
        return yaml.safe_load(f)


def _sha256_file(path: str | Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_frozen_checksums(checksums_path: str | Path) -> dict[str, str]:
    """Parse a `sha256sum`-format file into {relative_path: hexdigest}."""
    checksums: dict[str, str] = {}
    with open(checksums_path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            digest, path = line.split(maxsplit=1)
            checksums[Path(path).name] = digest
    return checksums


def build_sample_manifest(config: dict) -> pd.DataFrame:
    """Build the 26-row GSE245601 sample manifest from config, verifying
    every H5 file's SHA-256 against the frozen checksums file and every
    tumor patient's Control/Tamoxifen pairing.

    Raises ValueError on any checksum mismatch, missing file, or broken
    pairing -- rows are never silently dropped.
    """
    cfg = config["gse245601"]
    h5_dir = Path(cfg["h5_dir"])
    checksums = load_frozen_checksums(cfg["checksums_sha256"])

    rows = []
    for sample in cfg["samples"]:
        h5_path = h5_dir / sample["h5"]
        if not h5_path.exists():
            raise ValueError(f"{sample['gsm']}: expected H5 file not found at {h5_path}")

        expected_sha256 = checksums.get(sample["h5"])
        if expected_sha256 is None:
            raise ValueError(f"{sample['gsm']}: no frozen checksum recorded for {sample['h5']}")
        actual_sha256 = _sha256_file(h5_path)
        if actual_sha256 != expected_sha256:
            raise ValueError(
                f"{sample['gsm']}: checksum mismatch for {sample['h5']} "
                f"(expected {expected_sha256}, got {actual_sha256})"
            )

        rows.append(
            {
                "GSM": sample["gsm"],
                "srx": sample["srx"],
                "biosample": sample["biosample"],
                "subject": sample["subject"],
                "sample_type": sample["sample_type"],
                "patient": sample["patient"],
                "condition": sample["condition"],
                "paired_patient": sample["patient"],
                "h5_filename": sample["h5"],
                "file_size_bytes": h5_path.stat().st_size,
                "sha256": actual_sha256,
                "primary_analysis": bool(sample["primary_analysis"]),
            }
        )

    manifest = pd.DataFrame(rows)

    if len(manifest) != cfg["expected_n_samples"]:
        raise ValueError(f"expected {cfg['expected_n_samples']} samples, found {len(manifest)}")

    primary = manifest.loc[manifest["primary_analysis"]]
    n_tumor_patients = primary["patient"].nunique()
    if n_tumor_patients != cfg["expected_n_tumor_patients"]:
        raise ValueError(f"expected {cfg['expected_n_tumor_patients']} primary tumor patients, found {n_tumor_patients}")
    if len(primary) != cfg["expected_n_tumor_samples"]:
        raise ValueError(f"expected {cfg['expected_n_tumor_samples']} primary tumor samples, found {len(primary)}")

    condition_counts = primary.groupby("patient")["condition"].apply(lambda s: sorted(s.tolist()))
    for patient, conditions in condition_counts.items():
        if conditions != ["Control", "Tamoxifen"]:
            raise ValueError(f"{patient}: expected exactly one Control and one Tamoxifen sample, found {conditions}")

    return manifest


def write_sample_manifest(config_path: str | Path = "config/config.yaml") -> pd.DataFrame:
    config = _load_config(config_path)
    manifest = build_sample_manifest(config)
    output_path = Path(config["gse245601"]["output"]["sample_manifest_tsv"])
    output_path.parent.mkdir(parents=True, exist_ok=True)
    manifest.to_csv(output_path, sep="\t", index=False)
    return manifest


if __name__ == "__main__":
    manifest = write_sample_manifest()
    print(manifest.to_string(index=False))
