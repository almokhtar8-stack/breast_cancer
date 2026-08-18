"""Tests for the GSE245601 annotation-concordance feasibility gate (Task 3).

The requested concordance matrix was NOT built, because the authors' per-cell
labels are not publicly available. These tests pin that verdict to its evidence
so it cannot quietly become a guess, and assert that the module is structurally
incapable of emitting a concordance result while the gate fails.
"""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import pytest

from src import post_poster_annotation_feasibility as feas

ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = ROOT / "results" / "post_poster" / "annotation_concordance"
PROBE = OUT_DIR / "feasibility_probe.json"

pytestmark = pytest.mark.skipif(
    not (OUT_DIR / "feasibility_verdict.tsv").exists(),
    reason=("run scripts/post_poster_probe_gse245601_annotations.py then "
            "`python -m src.post_poster_annotation_feasibility` first"),
)


@pytest.fixture(scope="module")
def probe() -> dict:
    return json.loads(PROBE.read_text())


@pytest.fixture(scope="module")
def verdict() -> pd.Series:
    return pd.read_csv(OUT_DIR / "feasibility_verdict.tsv", sep="\t").iloc[0]


def test_gate_failed_and_the_analysis_was_stopped(verdict):
    assert verdict["feasibility_gate"] == "FAIL"
    assert verdict["per_cell_labels_publicly_obtainable"] in (False, "False")
    assert "STOP" in verdict["action_taken"]


def test_no_concordance_output_was_produced():
    """The confusion matrix, kappa and per-tumour concordance must not exist."""
    forbidden = ["confusion_matrix", "kappa", "concordance_by_tumor",
                 "per_tumour_concordance", "cohens_kappa"]
    produced = {p.name.lower() for p in OUT_DIR.glob("*")}
    for name in forbidden:
        assert not any(name in f for f in produced), f"unexpected output matching {name}"


def test_all_three_public_sources_were_checked(probe):
    assert set(probe) >= {"geo", "paper_supplement", "author_repo",
                          "controlled_access_route"}
    sources = pd.read_csv(OUT_DIR / "feasibility_sources_checked.tsv", sep="\t")
    assert len(sources) == 4
    assert (sources["evidence"].str.len() > 0).all()


def test_geo_carries_only_count_matrices(probe):
    geo = probe["geo"]
    assert geo["sample_supplementary_file_extensions"] == [".h5"]
    assert geo["n_samples"] == 26
    assert len(geo["series_supplementary_files"]) == 1
    assert "RAW.tar" in geo["series_supplementary_files"][0]


def test_paper_supplement_contains_no_cell_barcodes(probe):
    """The decisive evidence: every supplementary sheet was opened and scanned."""
    supp = probe["paper_supplement"]
    assert supp["n_xlsx_sheets_scanned"] >= 70
    assert supp["total_barcode_like_cells"] == 0
    # gene-level tables, not a 40k-cell annotation
    assert supp["max_sheet_rows"] < 5000


def test_author_repo_is_at_the_commit_pinned_in_preanalysis(probe):
    repo = probe["author_repo"]
    assert repo["pinned_commit"] == feas_pinned_commit()
    assert repo["head_matches_pinned_commit"] is True


def feas_pinned_commit() -> str:
    """The commit recorded in docs/gse245601_PREANALYSIS.md, read from the doc
    rather than duplicated as a literal."""
    text = (ROOT / "docs" / "gse245601_PREANALYSIS.md").read_text()
    import re
    m = re.search(r"`([0-9a-f]{40})`", text)
    assert m, "no pinned commit found in docs/gse245601_PREANALYSIS.md"
    return m.group(1)


def test_author_repo_never_contained_data_files(probe):
    repo = probe["author_repo"]
    assert repo["n_releases"] == 0
    assert repo["n_tags"] == 0
    assert repo["has_gitattributes_lfs_config"] is False
    assert set(repo["file_extensions_ever_committed"]) <= {".ipynb", ".md", ".png", ".r"}


def test_dbgap_is_recorded_as_controlled_access_not_as_a_public_comparator(probe):
    route = probe["controlled_access_route"]
    assert route["accession"].startswith("phs")
    assert route["status"] == "not queried"
    sources = pd.read_csv(OUT_DIR / "feasibility_sources_checked.tsv", sep="\t")
    dbgap = sources[sources["source"].str.contains("dbGaP")]
    assert len(dbgap) == 1
    # controlled access is neither a yes nor a no for public obtainability
    assert dbgap["per_cell_labels_present"].isna().all()
    assert int(pd.read_csv(OUT_DIR / "feasibility_verdict.tsv", sep="\t")
               ["n_public_sources_checked"].iloc[0]) == 3


def test_verdict_states_why_no_proxy_was_attempted(verdict):
    why = verdict["why_no_proxy"].lower()
    assert "our own pipeline" in why or "circular" in why


def test_our_label_inventory_totals_match_the_frozen_label_table():
    inv = pd.read_csv(OUT_DIR / "our_frozen_label_inventory.tsv", sep="\t")
    frozen = pd.read_csv(feas.OUR_LABELS_PATH, sep="\t")
    assert int(inv["n_cells"].sum()) == len(frozen)
    assert set(inv["patient"]) == set(frozen["patient"])
    for label in frozen["primary_malignancy_label"].unique():
        assert int(inv[label].sum()) == int(
            (frozen["primary_malignancy_label"] == label).sum())


def test_track_b_tumours_are_present_in_the_inventory():
    """Tumor_02/03/07 are the three Track B eligible tumours; the report
    comments on them, so they must actually be in the table."""
    inv = pd.read_csv(OUT_DIR / "our_frozen_label_inventory.tsv", sep="\t")
    assert {"Tumor_02", "Tumor_03", "Tumor_07"} <= set(inv["patient"])


def test_every_output_table_is_labelled_post_freeze():
    for path in sorted(OUT_DIR.glob("*.tsv")):
        df = pd.read_csv(path, sep="\t", nrows=5)
        assert feas.POST_FREEZE_LABEL in df.columns, f"{path.name} not labelled"


def test_probe_record_is_labelled_post_freeze(probe):
    assert probe["post_freeze_exploratory"] is True
