"""Behaviour of the final poster figure set: gates, outputs, determinism.

post_freeze_exploratory. These exercise the logic, not merely that the code
runs: each gate is corrupted deliberately and must fail loudly, plotted values
are compared against their source tables, and byte-reproducibility is checked
by rendering twice."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from src import (
    poster_corroboration_final,
    poster_dependency_final,
    poster_network_final,
    poster_pathway_final,
    poster_screen_final,
    poster_structure_final,
    poster_workflow_final,
)
from src.poster_candidate_volcano_v1 import CANDIDATES, load_panels
from src.poster_final_common import SIG_FDR, VerificationError, verify

MODULES = [poster_workflow_final, poster_screen_final, poster_corroboration_final,
           poster_pathway_final, poster_network_final, poster_dependency_final,
           poster_structure_final]


# --- the gate itself ---------------------------------------------------------
def test_verify_passes_on_matching_values():
    out = verify("t", [("a", 1.0, 1.0, 0.0), ("b", 2.5, 2.5001, 1e-3)])
    assert len(out) == 2 and out["match"].all()


def test_verify_fails_loudly_and_names_every_mismatch():
    with pytest.raises(VerificationError) as exc:
        verify("t", [("good", 1.0, 1.0, 0.0), ("bad1", 1.0, 2.0, 0.0), ("bad2", 3.0, 9.0, 0.0)])
    msg = str(exc.value)
    assert "bad1" in msg and "bad2" in msg and "good" not in msg
    assert "refusing to plot" in msg


# --- every figure's own gate fails on a corrupted value -----------------------
def test_screen_gate_fails_on_a_corrupted_effect_size(monkeypatch):
    hits, n = poster_screen_final.load_screen()
    broken = hits.copy()
    broken.loc[broken["gene"] == "KDM1A", "effect_size"] = -9.9
    with pytest.raises(VerificationError):
        poster_screen_final.gate(broken, n)


def test_screen_gate_fails_if_the_hit_count_changes():
    hits, n = poster_screen_final.load_screen()
    with pytest.raises(VerificationError):
        poster_screen_final.gate(hits.iloc[:12], n)


def test_corroboration_gate_fails_on_a_corrupted_fdr():
    import copy

    panels = copy.deepcopy(load_panels())
    panels[0].candidate_rows.loc["USP34", "fdr"] = 0.5
    with pytest.raises((ValueError, VerificationError)):
        poster_corroboration_final.gate(panels, poster_corroboration_final.load_pooled_evidence())


def test_pathway_gate_fails_on_a_corrupted_enrichment_score():
    long = poster_pathway_final.load_rows().copy()
    mask = (long["pathway_label"] == "Cell adhesion and motility") & (long["dataset"] == "gse245601")
    long.loc[mask, "NES"] = 2.0          # flip the acute value positive
    with pytest.raises(VerificationError):
        poster_pathway_final.gate(long)


def test_dependency_gate_fails_on_a_corrupted_count():
    summary = poster_dependency_final.load_dependency().copy()
    summary.loc[summary["gene"] == "TLK2", "n_strongly_dependent"] = 4
    with pytest.raises(VerificationError):
        poster_dependency_final.gate(summary)


def test_network_gate_fails_if_a_route_bypasses_the_bridge():
    from src.poster_network_mechanism_v2 import build_network

    stats = poster_network_final.analyse(build_network())
    broken = dict(stats)
    broken["paths"] = [["KDM1A", "SOMETHING_ELSE", "UBB", "USP34"]] + stats["paths"][1:]
    with pytest.raises(VerificationError):
        poster_network_final.gate(broken)


def test_structure_gate_fails_if_vezf1_is_given_a_structure():
    audit = poster_structure_final.load_audit().copy()
    audit.loc["VEZF1", "A_experimental_human_structure_exists"] = "True"
    with pytest.raises(VerificationError):
        poster_structure_final.gate(audit)


# --- outputs ------------------------------------------------------------------
@pytest.mark.parametrize("module", MODULES, ids=lambda m: m.FIGURE)
def test_figure_builds_in_all_three_formats(module, tmp_path):
    written, verification = module.main(tmp_path)
    assert set(written) == {"png", "pdf", "svg"}
    for path in written.values():
        assert path.exists() and path.stat().st_size > 0
    assert verification["match"].all()


@pytest.mark.parametrize("module", MODULES, ids=lambda m: m.FIGURE)
def test_png_and_pdf_bytes_are_reproducible(module, tmp_path):
    first, _ = module.main(tmp_path / "one")
    second, _ = module.main(tmp_path / "two")
    assert first["png"].read_bytes() == second["png"].read_bytes(), f"{module.FIGURE} PNG not reproducible"
    assert first["pdf"].read_bytes() == second["pdf"].read_bytes(), f"{module.FIGURE} PDF not reproducible"


# --- specific claims ----------------------------------------------------------
def test_network_layout_is_deterministic():
    from src.poster_network_mechanism_v2 import build_network

    G = build_network()
    a = poster_network_final.component_layout(G)
    b = poster_network_final.component_layout(G)
    assert set(a) == set(b)
    for node in a:
        assert a[node] == pytest.approx(b[node]), f"{node} moved between layouts"


def test_no_collision_solver_is_used_anywhere():
    """No renderer may IMPORT or CALL a label-collision solver.

    Prose is exempt on purpose: `poster_network_final`'s docstring names
    adjustText to record why it is not used, and a grep over raw source would
    fail on that explanation rather than on a real dependency. So this parses
    the module and inspects imports and call targets instead."""
    import ast
    import inspect

    for module in MODULES:
        tree = ast.parse(inspect.getsource(module))
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                for alias in node.names:
                    assert "adjusttext" not in alias.name.lower(), module.FIGURE
            elif isinstance(node, ast.ImportFrom):
                assert "adjusttext" not in (node.module or "").lower(), module.FIGURE
            elif isinstance(node, ast.Call):
                target = getattr(node.func, "id", "") or getattr(node.func, "attr", "")
                assert target != "adjust_text", module.FIGURE


def test_volcano_plots_source_values_with_zero_displacement(tmp_path):
    """Every candidate marker sits at its exact source log2 fold change."""
    import matplotlib.pyplot as plt

    import src.poster_final_common as common

    panels = load_panels()
    captured = {}
    original = common.save

    def spy(fig, stub):
        out = []
        for ax in fig.axes:
            if ax.get_xlabel() != "log$_2$ fold change":
                continue
            out.append([(float(c.get_offsets()[0][0]), float(c.get_offsets()[0][1]))
                        for c in ax.collections if len(c.get_offsets()) == 1])
        captured["panels"] = out
        return original(fig, stub)

    monkey = poster_corroboration_final.save
    poster_corroboration_final.save = spy
    try:
        poster_corroboration_final.main(tmp_path)
    finally:
        poster_corroboration_final.save = monkey
    plt.close("all")

    assert len(captured["panels"]) == 4
    checked = 0
    for drawn, panel in zip(captured["panels"], panels):
        assert len(drawn) == 4
        source = sorted((float(panel.candidate_rows.loc[g, "log2fc"]),
                         float(-np.log10(panel.candidate_rows.loc[g, "fdr"])))
                        for g in CANDIDATES)
        for (dx, dy), (sx, sy) in zip(sorted(drawn), source):
            assert dx == sx, f"{panel.accession}: plotted x {dx!r} != source {sx!r}"
            assert dy == sy
            checked += 1
    assert checked == 16


def test_bar_heights_equal_the_source_counts(tmp_path):
    """The dependency bars are the counts, not a rescaling of them."""
    import matplotlib.pyplot as plt

    import src.poster_final_common as common

    summary = poster_dependency_final.load_dependency().set_index("gene")
    captured = {}
    original = common.save

    def spy(fig, stub):
        ax = fig.axes[0]
        labels = [t.get_text() for t in ax.get_yticklabels()]
        widths = [p.get_width() for p in ax.patches]
        captured["labels"] = labels
        captured["widths"] = widths
        return original(fig, stub)

    monkey = poster_dependency_final.save
    poster_dependency_final.save = spy
    try:
        poster_dependency_final.main(tmp_path)
    finally:
        poster_dependency_final.save = monkey
    plt.close("all")

    n = len(captured["labels"])
    # first n patches are the denominator bars, next n the count bars
    denominators = captured["widths"][:n]
    counts = captured["widths"][n:2 * n]
    assert all(d == poster_dependency_final.N_LINES for d in denominators)
    for gene, drawn in zip(captured["labels"], counts):
        assert drawn == float(summary.loc[gene, "n_strongly_dependent"]), gene


def test_significance_encoding_is_fill_not_colour():
    """Colour identifies the gene; passing/failing is fill vs hollow ring."""
    from src.poster_final_common import significance_marker
    import matplotlib.pyplot as plt
    from src.poster_palette import GENE_COLOURS

    fig, ax = plt.subplots()
    filled = significance_marker(ax, 0, 0, "USP34", True)
    hollow = significance_marker(ax, 1, 1, "USP34", False)
    assert np.allclose(filled.get_facecolor()[0][:3],
                       [int(GENE_COLOURS["USP34"][i:i + 2], 16) / 255 for i in (1, 3, 5)], atol=0.01)
    assert hollow.get_facecolor().size == 0 or hollow.get_facecolor()[0][3] == 0
    assert np.allclose(hollow.get_edgecolor()[0][:3],
                       [int(GENE_COLOURS["USP34"][i:i + 2], 16) / 255 for i in (1, 3, 5)], atol=0.01)
    plt.close(fig)


def test_manifest_records_every_figure_with_hashes_and_captions():
    m = pd.read_csv("results/figures/poster_final/figure_manifest.tsv", sep="\t")
    assert len(m) == 7
    assert (m["analysis_status"] == "post_freeze_exploratory").all()
    for ext in ("png", "pdf", "svg"):
        assert m[f"sha256_{ext}"].str.len().eq(64).all()
        assert m[f"sha256_{ext}"].nunique() == 7
    caps = pd.read_csv("results/figures/poster_final/figure_captions.tsv", sep="\t")
    assert len(caps) == 7
    assert caps["caption_headline"].str.len().gt(20).all()


def test_verification_table_covers_every_figure_and_all_match():
    v = pd.read_csv("results/figures/poster_final/verification_against_frozen.tsv", sep="\t")
    assert v["match"].all()
    assert set(v["figure"].unique()) == {m.FIGURE for m in MODULES}
