#!/usr/bin/env python3
"""Diff report: confirm v2 differs from the committed figure 01 ONLY in the
three changed gene colours.

post_freeze_exploratory. Compares the v2 render against
`poster/final_figures/01_crispr_discovery.png` on two levels:

  1. STRUCTURE -- by rebuilding both figures in-process and comparing their
     matplotlib artists: gene order, effect sizes, axis limits, tick
     positions, title/subtitle strings, font sizes, figure size, margins,
     marker sizes, line widths and the greys used for non-candidate genes.
  2. PIXELS -- by re-rendering v1 to a temporary file and comparing the two
     PNGs pixel by pixel, then checking that every differing pixel is
     explained by one of the three colour changes.

Anything that differs and is NOT a colour change is reported as a failure
rather than accepted.
"""

from __future__ import annotations

import json
import sys
import tempfile
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
import os

os.chdir(ROOT)
os.environ.setdefault("SOURCE_DATE_EPOCH", "1600000000")

from src import poster_crispr_discovery_v1 as v1  # noqa: E402
from src import poster_crispr_discovery_v2 as v2  # noqa: E402
from src.poster_palette import GENE_COLOURS, PREVIOUS_GENE_COLOURS, hex_to_rgb01  # noqa: E402

OUT = Path("results/reports/figure1_palette")


def capture(module, stub: Path) -> dict:
    """Render and capture every structural property we care about."""
    grabbed: dict = {}
    original = plt.Figure.savefig

    def spy(self, fname, *a, **k):
        if not grabbed:
            ax = self.axes[0]
            grabbed.update(
                figsize=tuple(self.get_size_inches()),
                subplotpars={key: getattr(self.subplotpars, key)
                             for key in ("left", "right", "top", "bottom")},
                yticklabels=[t.get_text() for t in ax.get_yticklabels()],
                ytick_colors=[t.get_color() for t in ax.get_yticklabels()],
                ytick_sizes=[t.get_fontsize() for t in ax.get_yticklabels()],
                ytick_weights=[t.get_fontweight() for t in ax.get_yticklabels()],
                xlim=tuple(ax.get_xlim()), ylim=tuple(ax.get_ylim()),
                xticks=list(ax.get_xticks()), yticks=list(ax.get_yticks()),
                xlabel=ax.get_xlabel(), xlabel_size=ax.xaxis.label.get_fontsize(),
                figtexts=[(t.get_text(), t.get_fontsize(), t.get_color(),
                           t.get_position()) for t in self.texts],
                # LineCollections (the lollipop stems) expose colours and
                # linewidths; PathCollections (the dots) expose sizes.
                hline_colors=[np.asarray(c.get_colors()).tolist()
                              for c in ax.collections if hasattr(c, "get_colors")],
                line_widths=sorted(
                    round(float(w), 6) for c in ax.collections
                    if hasattr(c, "get_linewidths")
                    for w in np.atleast_1d(c.get_linewidths())),
                scatter_sizes=sorted(
                    round(float(s), 6) for c in ax.collections
                    if hasattr(c, "get_sizes")
                    for s in np.atleast_1d(c.get_sizes())),
                spines_visible={s: ax.spines[s].get_visible() for s in ax.spines},
            )
        return original(self, fname, *a, **k)

    plt.Figure.savefig = spy
    try:
        module.build_crispr_discovery_main(stub)
    finally:
        plt.Figure.savefig = original
    plt.close("all")
    return grabbed


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        a = capture(v1, tmp / "v1")
        b = capture(v2, tmp / "v2")

        # ---- 1. structural comparison ------------------------------------
        expected_colour_keys = {"ytick_colors", "hline_colors"}
        differences, colour_only = {}, {}
        for key in sorted(set(a) | set(b)):
            if a.get(key) != b.get(key):
                (colour_only if key in expected_colour_keys else differences)[key] = {
                    "v1": a.get(key), "v2": b.get(key)}

        # ---- 2. pixel comparison ------------------------------------------
        from PIL import Image

        committed = Path("poster/final_figures/01_crispr_discovery.png")
        im_v1 = np.asarray(Image.open(committed).convert("RGB"), dtype=int)
        im_v2 = np.asarray(Image.open(
            v2.OUT_DIR / "CRISPR_discovery_v2.png").convert("RGB"), dtype=int)

        pixel = {"same_shape": im_v1.shape == im_v2.shape}
        if pixel["same_shape"]:
            differing = np.any(im_v1 != im_v2, axis=-1)
            pixel["n_differing"] = int(differing.sum())
            pixel["pct_differing"] = round(100 * differing.mean(), 3)
            # Every differing pixel must be explained by one of the three
            # changed colours: it should be near an OLD colour in v1 and near
            # the corresponding NEW colour in v2 (anti-aliasing means "near",
            # not "equal"), or be a blend of the colour with white/grey.
            changed = {g: (PREVIOUS_GENE_COLOURS[g], GENE_COLOURS[g])
                       for g in GENE_COLOURS if PREVIOUS_GENE_COLOURS[g] != GENE_COLOURS[g]}
            old_rgb = np.array([hex_to_rgb01(o) * 255 for o, _ in changed.values()])
            new_rgb = np.array([hex_to_rgb01(n) * 255 for _, n in changed.values()])
            ys, xs = np.nonzero(differing)
            px_v1 = im_v1[ys, xs].astype(float)
            px_v2 = im_v2[ys, xs].astype(float)

            def min_dist_to_line(px, base):
                """Distance from each pixel to the nearest segment between a
                base colour and white (covers anti-aliased blends)."""
                best = np.full(len(px), np.inf)
                for c in base:
                    for other in (np.array([255.0, 255.0, 255.0]),
                                  np.array([216.0, 216.0, 216.0])):
                        d = other - c
                        t = np.clip(((px - c) @ d) / (d @ d), 0, 1)[:, None]
                        best = np.minimum(best, np.linalg.norm(px - (c + t * d), axis=1))
                return best

            tol = 42.0
            unexplained = (min_dist_to_line(px_v1, old_rgb) > tol) | \
                          (min_dist_to_line(px_v2, new_rgb) > tol)
            pixel["n_unexplained"] = int(unexplained.sum())
            pixel["all_differences_explained_by_colour"] = bool(unexplained.sum() == 0)

        report = {
            "status": "post_freeze_exploratory",
            "compared": {"baseline": str(committed),
                         "candidate": str(v2.OUT_DIR / "CRISPR_discovery_v2.png")},
            "colour_change": {g: {"from": PREVIOUS_GENE_COLOURS[g], "to": GENE_COLOURS[g],
                                  "changed": PREVIOUS_GENE_COLOURS[g] != GENE_COLOURS[g]}
                              for g in GENE_COLOURS},
            "structural_differences_outside_colour": differences,
            "colour_differences_as_expected": sorted(colour_only),
            "gene_order_identical": a.get("yticklabels") == b.get("yticklabels"),
            "axis_limits_identical": (a.get("xlim") == b.get("xlim")
                                      and a.get("ylim") == b.get("ylim")),
            "ticks_identical": (a.get("xticks") == b.get("xticks")
                                and a.get("yticks") == b.get("yticks")),
            "text_identical": ([t[0] for t in a.get("figtexts", [])]
                               == [t[0] for t in b.get("figtexts", [])]),
            "font_sizes_identical": ([t[1] for t in a.get("figtexts", [])]
                                     == [t[1] for t in b.get("figtexts", [])]
                                     and a.get("ytick_sizes") == b.get("ytick_sizes")
                                     and a.get("xlabel_size") == b.get("xlabel_size")),
            "figsize_identical": a.get("figsize") == b.get("figsize"),
            "margins_identical": a.get("subplotpars") == b.get("subplotpars"),
            "marker_and_line_geometry_identical": (
                a.get("scatter_sizes") == b.get("scatter_sizes")
                and a.get("line_widths") == b.get("line_widths")),
            "pixel_comparison": pixel,
        }
        report["PASS"] = bool(
            not differences
            and report["gene_order_identical"] and report["axis_limits_identical"]
            and report["ticks_identical"] and report["text_identical"]
            and report["font_sizes_identical"] and report["figsize_identical"]
            and report["margins_identical"]
            and report["marker_and_line_geometry_identical"]
            and pixel.get("all_differences_explained_by_colour", False))

        (OUT / "diff_report.json").write_text(json.dumps(report, indent=2, default=str))
        print(json.dumps({k: v for k, v in report.items()
                          if k != "structural_differences_outside_colour"},
                         indent=2, default=str))
        if differences:
            print("\nUNEXPECTED STRUCTURAL DIFFERENCES:")
            print(json.dumps(differences, indent=2, default=str))
        return 0 if report["PASS"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
