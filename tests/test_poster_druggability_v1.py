"""Small dedicated test for the structure/druggability figure -- verifies all
structural/pharmacological claims trace to the audited evidence table, that
no probe or substrate analog is presented as a drug, and that VEZF1 is not
given a fabricated structure."""

from __future__ import annotations

import hashlib
from pathlib import Path

from PIL import Image

from src import poster_druggability_v1 as dg


def _flat(path: str) -> str:
    """File text lowercased with whitespace collapsed, so assertions are not
    defeated by where a sentence happens to wrap."""
    return " ".join(Path(path).read_text().lower().split())

FIGURES = Path("results/figures/poster_druggability_v1")
RENDERS = FIGURES / "renders"
SHORTLIST = Path("results/tables/evidence_freeze/THERAPEUTIC_SHORTLIST_FREEZE.tsv")
FROZEN_SHA256 = "b6990d7e20f1191df7ebbca18be7542cc312e6e8489f05af6e4b4d31e66fb9dc"


def test_exactly_four_focus_genes():
    assert dg.FOCUS_FOUR == ["KDM1A", "TLK2", "USP34", "VEZF1"]
    assert set(dg.PANELS) == set(dg.FOCUS_FOUR)


def test_audit_source_is_the_frozen_audited_table():
    assert dg.AUDIT_TSV == Path(
        "results/tables/post_audit_sensitivity/06b_structural_tractability_audit.tsv")
    audit = dg.load_audit()
    assert list(audit.index) == dg.FOCUS_FOUR


def test_structural_ids_recovered_from_audit_not_fabricated():
    """Every PDB ID the figure renders must appear in the audited table's own
    text for that gene."""
    audit = dg.load_audit()
    recovered = dg.pdb_ids_from_audit(audit)
    assert "6NQU" in recovered["KDM1A"]
    assert "5O0Y" in recovered["TLK2"]
    assert "7W3U" in recovered["USP34"]
    for gene, spec in dg.PANELS.items():
        if spec["render"] is None:
            continue
        pdb_shown = spec["render"].split("_")[1].replace(".png", "")
        assert pdb_shown in recovered[gene], f"{pdb_shown} not traceable to audit for {gene}"


def test_rendered_structure_files_exist_and_are_real_images():
    for gene, spec in dg.PANELS.items():
        if spec["render"] is None:
            continue
        path = RENDERS / spec["render"]
        assert path.exists(), f"missing render for {gene}"
        with Image.open(path) as img:
            assert min(img.size) >= 800


def test_vezf1_has_no_structure_render_and_no_fabricated_pdb():
    assert dg.PANELS["VEZF1"]["render"] is None
    assert not list(RENDERS.glob("VEZF1*"))
    audit = dg.load_audit()
    assert str(audit.loc["VEZF1", "A_experimental_human_structure_exists"]).strip().lower() != "true"
    # the Zif268 homology template must never be presented as a VEZF1 structure
    assert "1AAY" not in str(dg.PANELS["VEZF1"])


def test_no_alphafold_or_homology_model_shown_as_experimental():
    text = Path("src/poster_druggability_v1.py").read_text()
    assert "alphafold" not in text.lower() or "no homology or AlphaFold model is\nsubstituted" in text
    assert dg.PANELS["VEZF1"]["ligand_line"] == "Homology model only"


def test_no_docking_or_affinity_prediction_anywhere():
    for path in ("src/poster_druggability_v1.py", "scripts/render_druggability_structures.py"):
        text = Path(path).read_text().lower()
        for banned in ("dock(", "autodock", "vina", "smina", "glide", "binding_affinity",
                       "predict_pose", "minimize"):
            assert banned not in text, f"{banned} found in {path}"


def test_evidence_track_derived_from_audit_columns():
    audit = dg.load_audit()
    track = dg.evidence_track(audit)
    assert track.loc["KDM1A"].tolist() == ["yes", "yes", "yes"]
    assert track.loc["TLK2"].tolist() == ["limited", "limited", "no"]
    assert track.loc["USP34"].tolist() == ["limited", "limited", "no"]
    assert track.loc["VEZF1"].tolist() == ["no", "no", "no"]
    # TLK2 and USP34 must NOT be forced into a strict ordering
    assert track.loc["TLK2"].tolist() == track.loc["USP34"].tolist()


def test_kdm1a_inhibitor_wording_matches_source():
    audit = dg.load_audit()
    assert str(audit.loc["KDM1A", "F_clinical_stage_pharmacology"]).strip().upper().startswith("YES")
    assert "iadademstat" in str(audit.loc["KDM1A", "E_validated_selective_small_molecule_inhibitor"])
    assert "iadademstat" in dg.PANELS["KDM1A"]["detail"]
    assert "selective inhibitor" in dg.PANELS["KDM1A"]["ligand_line"]


def test_tlk2_atp_analog_never_called_a_therapeutic_inhibitor():
    line = dg.PANELS["TLK2"]["ligand_line"].lower()
    assert "not an inhibitor" in line
    assert "no validated selective" in dg.PANELS["TLK2"]["pharm_line"].lower()
    audit = dg.load_audit()
    assert not str(audit.loc["TLK2", "E_validated_selective_small_molecule_inhibitor"]).strip().upper().startswith("YES")


def test_usp34_probe_never_called_a_selective_drug():
    line = dg.PANELS["USP34"]["ligand_line"].lower()
    assert "activity-based probe" in line and "not a drug" in line
    assert "no validated selective inhibitor" in dg.PANELS["USP34"]["pharm_line"].lower()
    note = _flat("results/reports/poster_druggability_v1/NOTE.md")
    assert "tractability hypothesis" in note
    assert "not proof of druggability" in note


def test_no_absolute_druggable_or_undruggable_claims():
    figure_and_note = (Path("src/poster_druggability_v1.py").read_text()
                       + Path("results/reports/poster_druggability_v1/NOTE.md").read_text()).lower()
    for absolute in ("is druggable", "is undruggable", "proven druggable", "not druggable at all"):
        assert absolute not in figure_and_note


def test_note_states_tractability_is_not_efficacy():
    note = _flat("results/reports/poster_druggability_v1/NOTE.md")
    assert "does not imply efficacy" in note
    assert "no docking" in note
    assert "none is breast-cancer-approved" in note


def test_figure_generation_succeeds_and_outputs_exist():
    stub = FIGURES / "DRUGGABILITY_v1"
    dg.build_druggability_v1(stub)
    for ext in ("png", "pdf", "svg"):
        path = stub.with_suffix(f".{ext}")
        assert path.exists()
        assert path.stat().st_size > 0


def test_png_has_a_sane_minimum_resolution():
    with Image.open(FIGURES / "DRUGGABILITY_v1.png") as image:
        width, height = image.size
        assert width >= 3000
        assert width > height


def test_frozen_shortlist_sha256_unchanged():
    digest = hashlib.sha256(SHORTLIST.read_bytes()).hexdigest()
    assert digest == FROZEN_SHA256
