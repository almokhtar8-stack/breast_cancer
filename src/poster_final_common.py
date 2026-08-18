"""Shared infrastructure for the final poster figure set.

post_freeze_exploratory. Holds the things every `src/poster_*_final.py`
renderer needs and must not re-implement: the output location, reproducible
saving, the significance-encoding primitives, and the verification-gate
helper that every figure is required to call before it draws anything.

The gate exists because it has already caught real errors: during the
volcano work it exposed that two of the source files named in a brief did
not reproduce the frozen candidate values (a redacted table that still
blinded KDM1A, and the wrong GSE240112 track). A figure that disagrees with
the frozen tables is worse than no figure, so `verify()` raises rather than
substituting.
"""

from __future__ import annotations

import hashlib
import logging
import os
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from src.poster_palette import GENE_COLOURS, NEUTRAL, WHITE

logger = logging.getLogger(__name__)

OUT_DIR = Path("results/figures/poster_final")
STATUS_LABEL = "post_freeze_exploratory"
SIG_FDR = 0.05          # candidate-level corroboration threshold
SCREEN_FDR = 0.10       # pre-registered CRISPR gate (PREANALYSIS.md Section 4)

# Provenance class of each poster gene. The four did NOT come from one
# selection rule and the figures must not imply they did.
PROVENANCE: dict[str, str] = {
    "USP34": "frozen_multimodal_rule",
    "VEZF1": "frozen_multimodal_rule",
    "KDM1A": "post_audit_addition",
    "TLK2": "post_audit_addition",
}
PROVENANCE_LABEL: dict[str, str] = {
    "frozen_multimodal_rule": "frozen multimodal rule (screen + ≥1 corroborating dataset)",
    "post_audit_addition": "added post-audit as a strongest-screen hit (no qualifying corroboration)",
}
PROVENANCE_MARKER: dict[str, str] = {"frozen_multimodal_rule": "o", "post_audit_addition": "s"}


class VerificationError(AssertionError):
    """A figure's plotted values disagree with the frozen source table."""


def verify(name: str, checks: list[tuple[str, float, float, float]]) -> "object":
    """The mandatory gate. `checks` is a list of
    (label, observed, expected, tolerance). Raises VerificationError listing
    every mismatch; returns a DataFrame of the passing comparisons."""
    import pandas as pd

    rows, bad = [], []
    for label, obs, exp, tol in checks:
        ok = abs(float(obs) - float(exp)) <= tol
        rows.append({"figure": name, "check": label, "observed": float(obs),
                     "expected": float(exp), "tolerance": tol, "match": ok})
        if not ok:
            bad.append(f"  {label}: got {obs!r}, expected {exp!r} (tol {tol})")
    if bad:
        raise VerificationError(
            f"{name}: verification gate FAILED against the frozen tables; refusing to plot:\n"
            + "\n".join(bad))
    logger.info("%s: verification gate passed (%d/%d values match the frozen tables)",
                name, len(rows), len(rows))
    return pd.DataFrame(rows)


def pin_reproducibility(salt: str) -> None:
    """Pinned PDF creation date (as scripts/poster/build_all.py does) and
    pinned SVG element-id salt, so output bytes are reproducible."""
    os.environ.setdefault("SOURCE_DATE_EPOCH", "1600000000")
    plt.rcParams["svg.hashsalt"] = salt


def save(fig, stub: Path) -> dict[str, Path]:
    """Write PNG, PDF and SVG with the post-freeze label in the metadata."""
    stub.parent.mkdir(parents=True, exist_ok=True)
    base = {"Title": f"{stub.stem} (poster candidate figure)"}
    written: dict[str, Path] = {}
    for ext, meta in (("png", {**base, "Description": STATUS_LABEL}),
                      ("pdf", {**base, "Subject": STATUS_LABEL, "Keywords": STATUS_LABEL}),
                      ("svg", {**base, "Description": STATUS_LABEL})):
        path = stub.with_suffix(f".{ext}")
        fig.savefig(path, facecolor=WHITE, dpi=300, metadata=meta)
        written[ext] = path
    plt.close(fig)
    logger.info("wrote %s.{png,pdf,svg}", stub)
    return written


def sha256(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def significance_marker(ax, x, y, gene: str, passes: bool, size: float = 200.0,
                        lw: float = 2.4, zorder: int = 10, clip_on: bool = True):
    """THE significance encoding, identical in every figure: solid fill when
    the value passes its threshold, hollow ring in the same gene colour when
    it does not. Colour carries gene identity only."""
    colour = GENE_COLOURS[gene]
    if passes:
        return ax.scatter([x], [y], s=size, c=colour, edgecolors=WHITE,
                          linewidths=1.4, zorder=zorder, clip_on=clip_on)
    return ax.scatter([x], [y], s=size, facecolors="none", edgecolors=colour,
                      linewidths=lw, zorder=zorder, clip_on=clip_on)


def style_axes(ax, *, grid_axis: str | None = None, spines=("top", "right")):
    for s in spines:
        ax.spines[s].set_visible(False)
    for s in ax.spines:
        if ax.spines[s].get_visible():
            ax.spines[s].set_color(NEUTRAL["ink_muted"])
    ax.tick_params(colors=NEUTRAL["ink_muted"], length=3)
    for lbl in ax.get_xticklabels() + ax.get_yticklabels():
        lbl.set_color(NEUTRAL["ink_2"])
    if grid_axis:
        ax.set_axisbelow(True)
        ax.grid(axis=grid_axis, color=NEUTRAL["grid"], linewidth=0.8)


# Every visible string on a figure must clear the 20 pt print floor once the
# figure is placed on the sheet. FONT is the single scale used by all seven
# renderers so the floor can be checked in one place and cannot drift.
FONT = {"tick": 21.0, "axis": 23.0, "annot": 21.0, "panel": 25.0,
        "legend": 21.0, "note": 21.0, "big": 30.0}
MIN_FIGURE_PT = min(FONT.values())


def figure_footer(fig, note: str | None = None):
    """Provenance line and status label. Both are figure text and therefore
    obey the same floor as everything else."""
    fig.text(0.995, 0.006, STATUS_LABEL, ha="right", va="bottom",
             fontsize=FONT["note"], color=NEUTRAL["ink_muted"], style="italic")
    if note:
        fig.text(0.006, 0.006, note, ha="left", va="bottom",
                 fontsize=FONT["note"], color=NEUTRAL["ink_2"])


# Captions are POSTER TEXT, not figure text. The poster bank's convention is a
# declarative finding above each figure with a smaller explanatory line beneath,
# and on the poster both are set in the poster's own type (50-60 pt headings,
# 25-30 pt body). Baking them into the figure instead forced them below the
# 20 pt print floor, so `headline()` no longer draws: it RECORDS the caption
# on the figure object, `poster_final_build` collects the recorded captions into
# the manifest, and `POSTER_TEXT.md` is where they are actually typeset.
CAPTIONS: dict[str, tuple[str, str]] = {}


def headline(fig, finding: str, explain: str, *, key: str | None = None, **_ignored):
    """Record the figure's caption. Draws nothing (see the note above)."""
    fig._poster_caption = (finding.replace("\n", " "), explain.replace("\n", " "))
    if key:
        CAPTIONS[key] = fig._poster_caption
    return fig._poster_caption
