"""Render the complete final poster figure set, with manifest and CVD checks.

post_freeze_exploratory. Renders all seven figures, writes the manifest
(source files, plotted values, SHA-256 per output file), runs the
deuteranopia/protanopia simulation on every rendered PNG, and records the
printed point sizes each figure achieves at its intended footprint on the
poster.

Every renderer runs its own verification gate against the frozen tables and
raises rather than substituting, so a failure here means a real disagreement
with the frozen data, not a rendering problem.
"""

from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd

from src import (
    poster_corroboration_final,
    poster_dependency_final,
    poster_network_final,
    poster_pathway_final,
    poster_screen_final,
    poster_structure_final,
    poster_workflow_final,
)
from src.poster_final_common import CAPTIONS, MIN_FIGURE_PT, OUT_DIR, STATUS_LABEL, sha256
from src.poster_palette import MIN_PRINT_PT, palette_cvd_report, print_point_size, simulate_image

logger = logging.getLogger(__name__)

# figure key -> (module, intended placed width on the A0 portrait sheet in mm,
#                smallest font size used anywhere in that figure, source note)
FIGURES = {
    "F1_methods_workflow": (poster_workflow_final, 380, MIN_FIGURE_PT,
                            "frozen CRISPR loaders, DepMap loader, pinned STRING query"),
    "F2_screen_certainty": (poster_screen_final, 380, MIN_FIGURE_PT,
                            "results/tables (frozen CRISPR screen: effect size, FDR)"),
    "F3_candidate_corroboration": (poster_corroboration_final, 760, MIN_FIGURE_PT,
                                   "4 frozen genome-wide DE tables + meta-analysis + power tables"),
    "F4_programme_signal": (poster_pathway_final, 380, MIN_FIGURE_PT,
                            "frozen Hallmark gene-set enrichment tables"),
    "F5_network_connectivity": (poster_network_final, 380, MIN_FIGURE_PT,
                                "data/reference/interactions/string_v2_* (pinned query)"),
    "F6_baseline_dependency": (poster_dependency_final, 380, MIN_FIGURE_PT,
                               "DepMap Public 26Q1 via the frozen loader"),
    "F7_reachability": (poster_structure_final, 380, MIN_FIGURE_PT,
                        "frozen structural audit table + committed PyMOL renders"),
}


def build_all(out_dir: Path = OUT_DIR) -> pd.DataFrame:
    out_dir.mkdir(parents=True, exist_ok=True)
    cvd_dir = out_dir / "cvd_simulation"
    cvd_dir.mkdir(parents=True, exist_ok=True)

    rows, verifications = [], []
    for key, (module, placed_mm, min_pt, source) in FIGURES.items():
        logger.info("rendering %s ...", key)
        written, verification = module.main(out_dir)
        caption = CAPTIONS.get(key, ("", ""))
        verifications.append(verification.assign(figure=key))
        sims = simulate_image(written["png"], cvd_dir / key)
        fig_width_in = _figure_width_inches(written["png"])
        printed = print_point_size(min_pt, fig_width_in, placed_mm)
        rows.append({
            "figure": key,
            "analysis_status": STATUS_LABEL,
            "render_module": f"src/{module.__name__.split('.')[-1]}.py",
            "primary_source": source,
            "placed_width_mm_on_A0_portrait": placed_mm,
            "rendered_width_in": round(fig_width_in, 2),
            "smallest_font_pt_on_figure": min_pt,
            "caption_headline": caption[0],
            "caption_explanation": caption[1],
            "smallest_printed_pt": round(printed, 1),
            "meets_20pt_floor": printed >= MIN_PRINT_PT,
            **{f"sha256_{ext}": sha256(p) for ext, p in written.items()},
            **{f"sha256_cvd_{k}": sha256(p) for k, p in sims.items()},
        })

    manifest = pd.DataFrame(rows)
    # The captions belong to the poster, not the figure; they are recorded here
    # so POSTER_TEXT.md and the figures can never drift apart.
    manifest[["figure", "caption_headline", "caption_explanation"]].to_csv(
        out_dir / "figure_captions.tsv", sep="\t", index=False)
    manifest.to_csv(out_dir / "figure_manifest.tsv", sep="\t", index=False)
    pd.concat(verifications, ignore_index=True).to_csv(
        out_dir / "verification_against_frozen.tsv", sep="\t", index=False)
    cvd = palette_cvd_report()
    cvd.insert(0, "analysis_status", STATUS_LABEL)
    cvd.to_csv(out_dir / "cvd_palette_simulation.tsv", sep="\t", index=False)

    failing = manifest.loc[~manifest["meets_20pt_floor"], "figure"].tolist()
    if failing:
        logger.warning("figures below the 20 pt print floor at their intended size: %s", failing)
    logger.info("built %d figures into %s", len(manifest), out_dir)
    return manifest


def _figure_width_inches(png: Path) -> float:
    from PIL import Image

    with Image.open(png) as im:
        dpi = im.info.get("dpi", (300, 300))[0] or 300
        return im.size[0] / dpi


def main(out_dir: Path = OUT_DIR) -> pd.DataFrame:
    return build_all(out_dir)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    m = main()
    print(m[["figure", "smallest_printed_pt", "meets_20pt_floor"]].to_string(index=False))
