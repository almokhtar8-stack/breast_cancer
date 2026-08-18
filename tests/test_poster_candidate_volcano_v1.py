"""Tests for the candidate volcano figure (proposed Figure 2 replacement).

These exercise the module's logic against the frozen DE tables: the
verification gate, the significance classification, the row accounting, and
the rendering pipeline. post_freeze_exploratory."""

from __future__ import annotations

import gzip
from pathlib import Path

import pandas as pd
import pytest

from src.poster_candidate_volcano_v1 import (
    CANDIDATES,
    EXPECTED_SIGNIFICANT,
    PANEL_ORDER,
    REFERENCE_FDR,
    REFERENCE_LOG2FC,
    SIG_FDR,
    STATUS_LABEL,
    build_figure,
    load_panels,
    reference_tolerance,
    verify_against_frozen,
    write_manifests,
)


@pytest.fixture(scope="module")
def panels():
    return load_panels()


@pytest.fixture(scope="module")
def verification(panels):
    return verify_against_frozen(panels)


# ---------------------------------------------------------------------------
# Verification gate
# ---------------------------------------------------------------------------
def test_loaders_return_the_frozen_candidate_values(panels):
    """Every one of the 16 (gene, dataset) FDRs -- and both quoted log2FCs --
    must match the frozen reference numbers within tolerance."""
    for panel in panels:
        for gene in CANDIDATES:
            quoted = REFERENCE_FDR[gene][panel.accession]
            got = float(panel.candidate_rows.loc[gene, "fdr"])
            assert abs(got - float(quoted)) <= reference_tolerance(quoted), (
                f"{gene}/{panel.accession}: extracted FDR {got} != frozen {quoted}"
            )
            q_lfc = REFERENCE_LOG2FC.get(gene, {}).get(panel.accession)
            if q_lfc is not None:
                got_lfc = float(panel.candidate_rows.loc[gene, "log2fc"])
                assert abs(got_lfc - float(q_lfc)) <= 1e-3


def test_gate_passes_and_covers_all_16_combinations(verification):
    assert len(verification) == 16
    assert verification["fdr_match"].all()
    checked_lfc = verification[verification["reference_log2fc"] != ""]
    assert len(checked_lfc) == 2  # USP34 and VEZF1 in GSE118713
    assert checked_lfc["log2fc_match"].all()


def test_gate_fails_loudly_on_a_mismatched_value(panels):
    """A figure that disagrees with the frozen tables is worse than no
    figure: corrupt one candidate FDR and the gate must raise."""
    import copy

    broken = copy.deepcopy(panels)
    broken[0].candidate_rows.loc["USP34", "fdr"] = 0.5
    with pytest.raises(ValueError, match="verification gate FAILED"):
        verify_against_frozen(broken)


def test_reference_tolerance_follows_quoted_precision():
    # brief's 1e-3 where the reference is quoted to >= 3 decimals
    assert reference_tolerance("0.494") == pytest.approx(1e-3)
    assert reference_tolerance("0.0073") == pytest.approx(1e-3)
    # widened to half a unit in the last place for 2-dp quotes: '0.91'
    # accepts exactly the values that round to 0.91
    assert reference_tolerance("0.91") == pytest.approx(0.005)


# ---------------------------------------------------------------------------
# Significance classification
# ---------------------------------------------------------------------------
def test_threshold_is_fdr_005_asserted_not_assumed():
    assert SIG_FDR == 0.05


def test_exactly_two_candidate_points_are_significant(panels):
    sig = [
        (p.accession, g)
        for p in panels
        for g in CANDIDATES
        if float(p.candidate_rows.loc[g, "fdr"]) < SIG_FDR
    ]
    assert len(sig) == EXPECTED_SIGNIFICANT == 2
    assert set(sig) == {("GSE118713", "USP34"), ("GSE240112", "VEZF1")}


def test_each_candidate_is_supported_by_at_most_one_dataset(panels):
    """The caption's claim, checked from the data it is rendered from."""
    for gene in CANDIDATES:
        n = sum(int(float(p.candidate_rows.loc[gene, "fdr"]) < SIG_FDR) for p in panels)
        assert n <= 1


# ---------------------------------------------------------------------------
# Row accounting: nothing dropped
# ---------------------------------------------------------------------------
def test_no_candidate_is_dropped_by_any_filtering_step(panels):
    for panel in panels:
        assert sorted(panel.candidate_rows.index) == sorted(CANDIDATES), (
            f"{panel.accession} lost a candidate"
        )


def test_gene_counts_entering_each_panel_are_recorded(panels):
    """Each panel records rows read and rows plotted; the only permitted
    difference is GSE118713's contrast selection (one of three contrasts)."""
    by_acc = {p.accession: p for p in panels}
    assert list(by_acc) == list(PANEL_ORDER)
    for p in panels:
        assert p.n_genes == len(p.df)
        if p.accession == "GSE118713":
            assert p.contrast == "TAMR_vs_MCF7"  # determined from the file
            assert p.rows_in == 3 * p.n_genes  # three contrasts, one kept
        else:
            assert p.contrast is None
            assert p.rows_in == p.n_genes  # nothing filtered


def test_panel_gene_counts_match_the_frozen_tables(panels):
    expected = {"GSE118713": 14838, "GSE111151": 27418,
                "GSE240112": 18428, "GSE245601": 17987}
    for p in panels:
        assert p.n_genes == expected[p.accession]


# ---------------------------------------------------------------------------
# Rendering
# ---------------------------------------------------------------------------
def test_figure_builds_and_produces_all_three_formats(panels, tmp_path):
    written = build_figure(panels, tmp_path / "volcano")
    assert set(written) == {"png", "pdf", "svg"}
    for path in written.values():
        assert path.exists() and path.stat().st_size > 0


def test_png_metadata_carries_the_post_freeze_label(panels, tmp_path):
    from PIL import Image

    written = build_figure(panels, tmp_path / "volcano")
    info = Image.open(written["png"]).info
    assert info.get("Description") == STATUS_LABEL


def test_manifest_records_sources_contrast_values_and_hashes(panels, verification, tmp_path):
    written = build_figure(panels, tmp_path / "volcano")
    write_manifests(panels, verification, written, tmp_path)

    manifest = pd.read_csv(tmp_path / "volcano_manifest.tsv", sep="\t")
    assert len(manifest) == 1
    row = manifest.iloc[0]
    assert row["analysis_status"] == STATUS_LABEL
    assert row["post_freeze"] == "yes"
    assert row["gse118713_contrast"] == "TAMR_vs_MCF7"
    assert row["significance_threshold_fdr"] == SIG_FDR
    for p in panels:
        assert row[f"source_{p.accession.lower()}"] == p.source_path
        assert row[f"n_genes_{p.accession.lower()}"] == p.n_genes
    for ext in ("png", "pdf", "svg"):
        assert len(row[f"sha256_{ext}"]) == 64

    values = pd.read_csv(tmp_path / "candidate_values_plotted.tsv", sep="\t")
    assert len(values) == 16
    assert (values["analysis_status"] == STATUS_LABEL).all()
    assert values["significant"].sum() == EXPECTED_SIGNIFICANT

    gate = pd.read_csv(tmp_path / "verification_against_frozen.tsv", sep="\t")
    assert len(gate) == 16
    assert gate["fdr_match"].all()


def test_pdf_bytes_are_reproducible(panels, tmp_path):
    """SOURCE_DATE_EPOCH is pinned, so two renders yield identical PDFs."""
    first = build_figure(panels, tmp_path / "one")["pdf"].read_bytes()
    second = build_figure(panels, tmp_path / "two")["pdf"].read_bytes()
    assert first == second


def test_loader_rejects_a_table_with_nonpositive_fdr(tmp_path, monkeypatch):
    """Stop-and-ask rather than silently degrade: an FDR of 0 cannot be drawn
    at -log10 and must raise, not clip."""
    import src.poster_candidate_volcano_v1 as mod

    src_cfg = mod._load_config()["cross_dataset_genomewide"]["inputs"]
    df = pd.read_csv(src_cfg["gse111151_de_tsv"], sep="\t")
    df.loc[0, "fdr"] = 0.0
    bad = tmp_path / "bad.tsv.gz"
    with gzip.open(bad, "wt") as f:
        df.to_csv(f, sep="\t", index=False)

    def fake_load(config_path=mod.CONFIG_PATH):
        cfg = {"cross_dataset_genomewide": {"inputs": dict(src_cfg)}}
        cfg["cross_dataset_genomewide"]["inputs"]["gse111151_de_tsv"] = str(bad)
        return cfg

    monkeypatch.setattr(mod, "_load_config", fake_load)
    with pytest.raises(ValueError, match="finite and > 0"):
        mod.load_panels()
