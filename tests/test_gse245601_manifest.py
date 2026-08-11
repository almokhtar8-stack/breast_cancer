import copy
from pathlib import Path

import pytest
import yaml

from src.gse245601_manifest import build_sample_manifest, load_frozen_checksums

REPO_ROOT = Path(__file__).parent.parent
CONFIG_PATH = REPO_ROOT / "config" / "config.yaml"


def _real_config() -> dict:
    with open(CONFIG_PATH) as f:
        return yaml.safe_load(f)


class TestSampleManifest:
    def test_26_samples_with_correct_checksums(self):
        config = _real_config()
        manifest = build_sample_manifest(config)
        assert len(manifest) == 26

    def test_20_primary_tumor_samples_10_patients(self):
        config = _real_config()
        manifest = build_sample_manifest(config)
        primary = manifest.loc[manifest["primary_analysis"]]
        assert len(primary) == 20
        assert primary["patient"].nunique() == 10
        assert set(primary["patient"]) == {f"Tumor_{i:02d}" for i in range(1, 11)}

    def test_every_primary_patient_has_exactly_one_control_one_tamoxifen(self):
        config = _real_config()
        manifest = build_sample_manifest(config)
        primary = manifest.loc[manifest["primary_analysis"]]
        for patient, sub in primary.groupby("patient"):
            assert sorted(sub["condition"]) == ["Control", "Tamoxifen"]
            assert len(sub) == 2

    def test_normal_and_t47d_excluded_from_primary(self):
        config = _real_config()
        manifest = build_sample_manifest(config)
        non_primary = manifest.loc[~manifest["primary_analysis"]]
        assert set(non_primary["patient"]) == {"Normal_01", "Normal_02", "T47D"}
        assert len(non_primary) == 6

    def test_gsm_accessions_are_unique_and_well_formed(self):
        config = _real_config()
        manifest = build_sample_manifest(config)
        assert manifest["GSM"].is_unique
        assert manifest["GSM"].str.match(r"^GSM\d+$").all()

    def test_checksum_mismatch_raises(self, tmp_path):
        config = copy.deepcopy(_real_config())
        checksums = load_frozen_checksums(config["gse245601"]["checksums_sha256"])
        first_key = next(iter(checksums))
        corrupted = dict(checksums)
        corrupted[first_key] = "0" * 64
        bad_checksums_path = tmp_path / "checksums.sha256"
        with open(bad_checksums_path, "w") as f:
            for name, digest in corrupted.items():
                f.write(f"{digest}  h5/{name}\n")
        config["gse245601"]["checksums_sha256"] = str(bad_checksums_path)
        with pytest.raises(ValueError, match="checksum mismatch"):
            build_sample_manifest(config)

    def test_h5_filenames_match_gsm_and_patient_condition_labels(self):
        config = _real_config()
        manifest = build_sample_manifest(config)
        for _, row in manifest.iterrows():
            assert row["h5_filename"].startswith(row["GSM"])
            assert row["patient"] in row["h5_filename"]
            assert row["condition"] in row["h5_filename"]
