#!/usr/bin/env python3
"""Rebuilds the six canonical poster figures, copies them into
`poster/final_figures/`, and verifies/writes `poster/figure_manifest.tsv`.

What this DOES:
  - re-renders each figure from its canonical `src/` implementation, using the
    already-frozen / already-processed tables those implementations read;
  - copies each output to its canonical poster filename and verifies the copy
    is byte-identical to its source (SHA-256);
  - writes the figure manifest with fresh hashes.

Determinism note (measured, not assumed):
  - PNG output is byte-reproducible across runs.
  - PDF output is byte-reproducible only when `SOURCE_DATE_EPOCH` is set (it
    otherwise embeds a wall-clock /CreationDate); this script sets it.
  - SVG output is NOT byte-reproducible -- matplotlib embeds per-run element
    ids. SVG hashes in the manifest therefore describe the published file, and
    are not expected to survive a re-render.

What this does NOT do:
  - it does not re-run the upstream scientific pipelines (CRISPR fitting,
    differential expression, GSEA, DepMap ingestion, STRING download,
    PyMOL structure rendering). Those are separate, heavier steps -- see
    `docs/analysis_map.md`. Figure 06 reuses the PyMOL renders already in
    `results/figures/poster_druggability_v1/renders/`.

Usage:
    python scripts/poster/build_all.py            # build + verify + write manifest
    python scripts/poster/build_all.py --check    # verify only, do not re-render
"""

from __future__ import annotations

import argparse
import hashlib
import logging
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

FINAL_DIR = ROOT / "poster" / "final_figures"
MANIFEST = ROOT / "poster" / "figure_manifest.tsv"
EXTS = ("png", "pdf", "svg")

# figure_number, canonical_name, src module, build fn, output stub, and the
# manifest metadata columns.
FIGURES = [
    dict(number=1, name="01_crispr_discovery",
         module="poster_crispr_discovery_v1", func="build_crispr_discovery_main",
         stub="CRISPR_discovery_main",
         question="Which gene knockouts sensitise ER+ breast-cancer cells to tamoxifen?",
         wrapper="scripts/poster/01_crispr_discovery.py",
         primary_source_data="data/processed/labels.parquet (Hany et al. genome-scale screen, 19,103 fitted genes)",
         analysis_status="frozen", post_freeze="no"),
    dict(number=2, name="02_candidate_expression",
         module="poster_hero_heatmap_v6", func="build_hero_heatmap_v6",
         stub="HERO_sample_level_heatmap_v6",
         question="How do the four candidates behave across resistance, recurrence and acute tamoxifen response?",
         wrapper="scripts/poster/02_candidate_expression.py",
         primary_source_data="GSE118713, GSE111151, GSE240112, GSE245601 processed expression (29 biological rows)",
         analysis_status="derived_from_frozen", post_freeze="no"),
    dict(number=3, name="03_molecular_networks",
         module="poster_network_mechanism_v4", func="build_network_mechanism_main",
         stub="NETWORK_mechanism_v4",
         question="What molecular neighbourhood surrounds each candidate, and are those neighbourhoods connected?",
         wrapper="scripts/poster/03_molecular_networks.py",
         primary_source_data="data/reference/interactions/string_v2_level{1,2}_{functional,physical}.tsv (STRING, score >= 0.7, species 9606)",
         analysis_status="post_freeze_exploratory", post_freeze="yes"),
    dict(number=4, name="04_pathway_remodeling",
         module="poster_pathway_v2", func="build_pathway_v2",
         stub="PATHWAY_v2",
         question="Do network-relevant biological programs change in the transcriptomic contexts?",
         wrapper="scripts/poster/04_pathway_remodeling.py",
         primary_source_data="results/tables/systems_network/gsea_{gse118713,gse111151,gse240112,gse245601}.tsv (frozen GSEA)",
         analysis_status="derived_from_frozen", post_freeze="no"),
    dict(number=5, name="05_depmap_dependency",
         module="poster_depmap_v2", func="build_depmap_v2",
         stub="DEPMAP_v2",
         question="Does tamoxifen sensitisation occur in genes cancer cells already depend on at baseline?",
         wrapper="scripts/poster/05_depmap_dependency.py",
         primary_source_data="DepMap 26Q1 CRISPRGeneEffect/CRISPRGeneDependency + frozen CRISPR effects; cached extract in results/tables/poster_depmap_v1/",
         analysis_status="derived_from_frozen", post_freeze="no"),
    dict(number=6, name="06_structural_tractability",
         module="poster_druggability_v1", func="build_druggability_v1",
         stub="DRUGGABILITY_v1",
         question="Can these candidate vulnerabilities realistically be targeted?",
         wrapper="scripts/poster/06_structural_tractability.py",
         primary_source_data="results/tables/post_audit_sensitivity/06b_structural_tractability_audit.tsv + PDB 6NQU/5O0Y/7W3U renders",
         analysis_status="derived_from_frozen", post_freeze="no"),
]

MANIFEST_COLUMNS = [
    "figure_number", "figure_name", "scientific_question", "poster_png", "poster_pdf",
    "poster_svg", "canonical_source", "render_module", "wrapper_script",
    "primary_source_data", "analysis_status", "post_freeze",
    "sha256_png", "sha256_pdf", "sha256_svg",
]

logger = logging.getLogger("build_all")


def rel(path: Path) -> str:
    """Repo-relative POSIX string. The canonical `src/` modules define
    repo-relative OUT_DIRs, so `path` may already be relative."""
    path = Path(path)
    if path.is_absolute():
        path = path.relative_to(ROOT)
    return path.as_posix()


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def render(fig: dict) -> Path:
    module = __import__(f"src.{fig['module']}", fromlist=[fig["func"], "OUT_DIR"])
    stub = module.OUT_DIR / fig["stub"]
    getattr(module, fig["func"])(stub)
    return stub


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true",
                        help="verify existing outputs and manifest without re-rendering")
    args = parser.parse_args()
    logging.basicConfig(level=logging.WARNING, format="%(levelname)s %(message)s")

    # The canonical renderers read and write repo-relative paths, so the build
    # must run from the repository root regardless of the caller's cwd.
    import os
    os.chdir(ROOT)

    # Make PDF output reproducible: without this matplotlib stamps a wall-clock
    # /CreationDate into every PDF, so byte hashes would change on every run.
    os.environ.setdefault("SOURCE_DATE_EPOCH", "1600000000")

    FINAL_DIR.mkdir(parents=True, exist_ok=True)
    rows, failures = [], []

    for fig in FIGURES:
        module = __import__(f"src.{fig['module']}", fromlist=["OUT_DIR"])
        stub = module.OUT_DIR / fig["stub"]

        if args.check:
            missing = [e for e in EXTS if not stub.with_suffix(f".{e}").exists()]
            if missing:
                failures.append(f"figure {fig['number']}: missing source {missing}")
                continue
        else:
            print(f"[{fig['number']}/6] rendering {fig['name']} ...", flush=True)
            stub = render(fig)

        row = {
            "figure_number": fig["number"],
            "figure_name": fig["name"],
            "scientific_question": fig["question"],
            "canonical_source": rel(stub) + ".{png,pdf,svg}",
            "render_module": f"src/{fig['module']}.py",
            "wrapper_script": fig["wrapper"],
            "primary_source_data": fig["primary_source_data"],
            "analysis_status": fig["analysis_status"],
            "post_freeze": fig["post_freeze"],
        }
        for ext in EXTS:
            src = stub.with_suffix(f".{ext}")
            dst = FINAL_DIR / f"{fig['name']}.{ext}"
            if not args.check:
                shutil.copy2(src, dst)
            if not dst.exists():
                failures.append(f"figure {fig['number']}: missing canonical copy {dst.name}")
                continue
            src_hash, dst_hash = sha256(src), sha256(dst)
            if src_hash != dst_hash:
                failures.append(f"figure {fig['number']}: {dst.name} differs from its source")
            row[f"poster_{ext}"] = rel(dst)
            row[f"sha256_{ext}"] = dst_hash
        rows.append(row)

    if failures:
        print("\nFAILED:")
        for f in failures:
            print("  -", f)
        return 1

    if not args.check:
        lines = ["\t".join(MANIFEST_COLUMNS)]
        for row in rows:
            lines.append("\t".join(str(row.get(c, "")) for c in MANIFEST_COLUMNS))
        MANIFEST.write_text("\n".join(lines) + "\n")

    print(f"\nOK: {len(rows)}/6 canonical figures present, all copies byte-identical to source.")
    print(f"     figures  -> {rel(FINAL_DIR)}/")
    print(f"     manifest -> {rel(MANIFEST)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
