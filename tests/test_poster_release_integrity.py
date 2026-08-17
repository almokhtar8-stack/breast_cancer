"""Release-integrity tests for the public poster deliverable.

These do NOT re-render figures. They verify that what is published under
`poster/` is exactly what the canonical `results/figures/` sources contain, that
the manifest is internally consistent and honest about frozen vs post-freeze
status, and that the science-freeze anchors are unmoved.
"""

from __future__ import annotations

import hashlib
import subprocess
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
POSTER = ROOT / "poster"
FINAL = POSTER / "final_figures"
MANIFEST = POSTER / "figure_manifest.tsv"

FROZEN_SHORTLIST = ROOT / "results/tables/evidence_freeze/THERAPEUTIC_SHORTLIST_FREEZE.tsv"
FROZEN_SHORTLIST_SHA256 = "b6990d7e20f1191df7ebbca18be7542cc312e6e8489f05af6e4b4d31e66fb9dc"
FREEZE_TAG = "science-freeze-2026-08-15"
FREEZE_COMMIT_PREFIX = "9a1b777"

EXPECTED_FIGURES = [
    "01_crispr_discovery", "02_candidate_expression", "03_molecular_networks",
    "04_pathway_remodeling", "05_depmap_dependency", "06_structural_tractability",
]
EXTS = ("png", "pdf", "svg")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


@pytest.fixture(scope="module")
def manifest() -> pd.DataFrame:
    assert MANIFEST.exists(), "poster/figure_manifest.tsv is missing"
    return pd.read_csv(MANIFEST, sep="\t")


def test_manifest_has_exactly_the_six_canonical_figures(manifest):
    assert list(manifest["figure_name"]) == EXPECTED_FIGURES
    assert list(manifest["figure_number"]) == [1, 2, 3, 4, 5, 6]


def test_every_published_figure_file_exists(manifest):
    for name in EXPECTED_FIGURES:
        for ext in EXTS:
            path = FINAL / f"{name}.{ext}"
            assert path.exists(), f"missing {path}"
            assert path.stat().st_size > 0


def test_published_hashes_match_the_files_on_disk(manifest):
    """The manifest is the public provenance record -- it must not drift."""
    for row in manifest.itertuples():
        for ext in EXTS:
            declared = getattr(row, f"sha256_{ext}")
            actual = sha256(ROOT / getattr(row, f"poster_{ext}"))
            assert declared == actual, f"{row.figure_name}.{ext}: manifest hash != file hash"


# Figure 03's label placement uses adjustText, whose collision solver is not
# byte-reproducible across renders (documented in
# src/poster_network_mechanism_v4.py: seeding numpy and pinning PYTHONHASHSEED
# do not fix it, and removing the solver reintroduces real label collisions).
# Its GRAPH is deterministic and asserted in test_poster_network_mechanism_v4.py.
NON_BYTE_REPRODUCIBLE_PNG = {"03_molecular_networks"}


def test_published_png_is_byte_identical_to_its_canonical_source(manifest):
    """poster/final_figures must hold COPIES, never independently re-rendered
    output.

    Only PNG is checked byte-for-byte, and that is deliberate: matplotlib PNG
    output is byte-reproducible, PDF output is reproducible only when
    SOURCE_DATE_EPOCH is set (it otherwise embeds a wall-clock /CreationDate),
    and SVG output is never byte-reproducible because matplotlib emits per-run
    element ids. Re-rendering a figure therefore legitimately changes its
    PDF/SVG bytes while the science is unchanged, so requiring PDF/SVG equality
    here would fail purely as a function of render order.
    """
    checked = 0
    for row in manifest.itertuples():
        if row.figure_name in NON_BYTE_REPRODUCIBLE_PNG:
            continue
        stub = row.canonical_source.replace(".{png,pdf,svg}", "")
        src = ROOT / f"{stub}.png"
        assert src.exists(), f"canonical source missing: {src}"
        assert sha256(src) == sha256(ROOT / row.poster_png), (
            f"{row.figure_name}.png differs from its canonical source -- "
            "re-run: python scripts/poster/build_all.py")
        checked += 1
    assert checked == len(EXPECTED_FIGURES) - len(NON_BYTE_REPRODUCIBLE_PNG)


def test_non_reproducible_figure_is_explicitly_documented():
    """The one figure exempted above must say so in its own source, so the
    exemption can never become silent."""
    src = (ROOT / "src/poster_network_mechanism_v4.py").read_text()
    flat = " ".join(src.split())
    assert "NONDETERMINISM" in src
    assert "adjustText" in src
    assert "LABEL POSITIONS ONLY" in src
    poster_readme = " ".join((POSTER / "README.md").read_text().split())
    assert "not byte-reproducible" in poster_readme


def test_published_pdf_and_svg_sources_exist_and_are_non_empty(manifest):
    """PDF/SVG bytes are not reproducible (see the PNG test), so verify the
    canonical vector sources exist and are real files rather than comparing
    hashes across renders."""
    for row in manifest.itertuples():
        stub = row.canonical_source.replace(".{png,pdf,svg}", "")
        for ext in ("pdf", "svg"):
            src = ROOT / f"{stub}.{ext}"
            assert src.exists(), f"canonical source missing: {src}"
            assert src.stat().st_size > 0
            published = ROOT / getattr(row, f"poster_{ext}")
            assert published.exists() and published.stat().st_size > 0


def test_manifest_render_modules_and_wrappers_exist(manifest):
    for row in manifest.itertuples():
        assert (ROOT / row.render_module).exists(), f"missing {row.render_module}"
        assert (ROOT / row.wrapper_script).exists(), f"missing {row.wrapper_script}"


def test_manifest_declares_the_network_layer_as_post_freeze(manifest):
    """The STRING network analysis was run AFTER the science freeze and must
    never be presented as frozen science."""
    indexed = manifest.set_index("figure_name")
    assert indexed.loc["03_molecular_networks", "analysis_status"] == "post_freeze_exploratory"
    assert indexed.loc["03_molecular_networks", "post_freeze"] == "yes"
    # and nothing else claims post-freeze status
    others = manifest[manifest["figure_name"] != "03_molecular_networks"]
    assert set(others["post_freeze"]) == {"no"}
    assert set(others["analysis_status"]) <= {"frozen", "derived_from_frozen"}


def test_analysis_status_values_are_from_the_documented_vocabulary(manifest):
    allowed = {"frozen", "derived_from_frozen", "post_freeze_exploratory"}
    assert set(manifest["analysis_status"]) <= allowed


def test_frozen_shortlist_sha256_unchanged():
    assert sha256(FROZEN_SHORTLIST) == FROZEN_SHORTLIST_SHA256


def test_science_freeze_tag_still_points_at_the_frozen_commit():
    out = subprocess.run(["git", "log", "-1", "--format=%H", FREEZE_TAG],
                         cwd=ROOT, capture_output=True, text=True, check=True)
    assert out.stdout.strip().startswith(FREEZE_COMMIT_PREFIX)


def test_historical_frozen_shortlist_is_documented_as_distinct_from_poster_focus():
    """A reader must not be able to mistake the poster four for the frozen four."""
    readme = (ROOT / "README.md").read_text()
    poster_readme = (POSTER / "README.md").read_text()
    for text in (readme, poster_readme):
        flat = " ".join(text.split())
        assert "EML5" in flat and "CITED2" in flat, "historical shortlist not stated"
        assert "not" in flat.lower()
    assert "not** the historical frozen" in " ".join(readme.split()) or \
           "not the historical frozen" in " ".join(readme.split()).replace("**", "")


def test_public_docs_do_not_hardcode_machine_local_paths():
    """New public-facing navigation docs must stay portable."""
    for rel in ("poster/README.md", "data/README.md", "tests/README.md",
                "docs/analysis_map.md", "results/figures/README.md",
                "results/reports/README.md"):
        text = (ROOT / rel).read_text()
        assert "/ibex/" not in text, f"{rel} leaks a machine-local path"
        assert "/home/" not in text, f"{rel} leaks a home path"


def test_no_docking_in_any_poster_wrapper_or_build():
    for path in sorted((ROOT / "scripts/poster").glob("*.py")):
        text = path.read_text().lower()
        for banned in ("autodock", "vina", "docking", "binding_affinity"):
            assert banned not in text, f"{path.name} references {banned}"
