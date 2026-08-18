"""The project's single source of truth for poster colour.

post_freeze_exploratory. No figure defines its own colours: every renderer
under `src/poster_*_final.py` imports from here. Enforced by test
(`tests/test_poster_final_palette.py`): no hex literal may appear outside
this module, no CHROME colour may appear inside a figure, and no colour may
be imported from a v1/v2 renderer.

Encoding rules that the colours exist to serve:

  * Colour identifies the GENE and never its rank or its significance.
  * Significance is encoded by FILL: solid = passes the threshold, hollow
    ring in the same gene colour = does not pass. Identical across figures.
  * CHROME is poster layout only (title band, section rules, QR frame) and
    must never appear inside a plotted figure.

Colour-vision safety is not asserted here, it is computed: `cvd.py`-free,
`simulate_cvd()` below applies the Machado, Oliveira & Fernandes (2009)
severity-1.0 matrices in linear sRGB. No colour-vision library exists in
this environment and none is added.
"""

from __future__ import annotations

import itertools

import numpy as np

# --- data-tier colours ------------------------------------------------------
GENE_COLOURS: dict[str, str] = {
    "USP34": "#6A3D9A",
    "KDM1A": "#D55E00",
    "TLK2": "#009E73",
    "VEZF1": "#56B4E9",
}

NEUTRAL: dict[str, str] = {
    "ink": "#262626",
    "ink_2": "#555555",
    "ink_muted": "#8c8c8c",
    "rule": "#b0b0b0",
    "backdrop": "#c9c9c9",
    "grid": "#e6e6e6",
    "tint": "#f0f0f0",
}

# Sequential/diverging ramp for enrichment scores. Cool = suppressed,
# warm = enriched, near-neutral centre.
DIVERGING: list[str] = ["#2E6C8E", "#7FAFC4", "#F3EEE4", "#E3A180", "#C1543A"]

# --- poster layout only; NEVER inside a figure ------------------------------
CHROME: dict[str, str] = {
    "space": "#1E1233",
    "violet_mid": "#2A1A45",
    "violet": "#6A3D9A",
    "magenta": "#E0459B",
    "violet_light": "#B49BE0",
    "near_white": "#F0EAFB",
    "warm_white": "#F3EEE4",
}

WHITE = "#FFFFFF"  # page ground and marker edge only, not an identity colour

# Every colour a figure is permitted to draw with.
FIGURE_ALLOWED: frozenset[str] = frozenset(
    {c.upper() for c in GENE_COLOURS.values()}
    | {c.upper() for c in NEUTRAL.values()}
    | {c.upper() for c in DIVERGING}
    | {WHITE}
)

# --- typography -------------------------------------------------------------
# An A0 sheet is 841 x 1189 mm. A figure rendered `width_in` inches wide and
# placed `placed_mm` wide on the sheet is scaled by (placed_mm/25.4)/width_in;
# a font of `pt` points on the figure therefore prints at pt * that factor.
# The brief's floor is 20 pt at final print size.
MIN_PRINT_PT: float = 20.0


def print_point_size(fig_pt: float, fig_width_in: float, placed_mm: float) -> float:
    """Printed point size of `fig_pt` text when the figure is placed
    `placed_mm` wide on the poster."""
    return fig_pt * (placed_mm / 25.4) / fig_width_in


# --- colour-vision simulation (computed, not asserted) ----------------------
CVD_MATRICES: dict[str, np.ndarray] = {
    "protanopia": np.array([[0.152286, 1.052583, -0.204868],
                            [0.114503, 0.786281, 0.099216],
                            [-0.003882, -0.048116, 1.051998]]),
    "deuteranopia": np.array([[0.367322, 0.860646, -0.227968],
                              [0.280085, 0.672501, 0.047413],
                              [-0.011820, 0.042940, 0.968881]]),
}


def hex_to_rgb01(h: str) -> np.ndarray:
    h = h.lstrip("#")
    return np.array([int(h[i:i + 2], 16) for i in (0, 2, 4)], dtype=float) / 255.0


def _srgb_to_linear(c: np.ndarray) -> np.ndarray:
    return np.where(c <= 0.04045, c / 12.92, ((c + 0.055) / 1.055) ** 2.4)


def _linear_to_srgb(c: np.ndarray) -> np.ndarray:
    c = np.clip(c, 0.0, 1.0)
    return np.where(c <= 0.0031308, 12.92 * c, 1.055 * np.power(c, 1 / 2.4) - 0.055)


def simulate_cvd(rgb01: np.ndarray, kind: str) -> np.ndarray:
    """rgb01 (..., 3) in [0,1] -> simulated sRGB in [0,1]."""
    return _linear_to_srgb(_srgb_to_linear(rgb01) @ CVD_MATRICES[kind].T)


def rgb01_to_lab(rgb01: np.ndarray) -> np.ndarray:
    lin = _srgb_to_linear(rgb01)
    m = np.array([[0.4124564, 0.3575761, 0.1804375],
                  [0.2126729, 0.7151522, 0.0721750],
                  [0.0193339, 0.1191920, 0.9503041]])
    xyz = (lin @ m.T) / np.array([0.95047, 1.0, 1.08883])
    eps, kappa = 216 / 24389, 24389 / 27
    f = np.where(xyz > eps, np.cbrt(xyz), (kappa * xyz + 16) / 116)
    return np.stack([116 * f[..., 1] - 16,
                     500 * (f[..., 0] - f[..., 1]),
                     200 * (f[..., 1] - f[..., 2])], axis=-1)


def palette_cvd_report(colours: dict[str, str] | None = None):
    """Pairwise CIE76 dE and dL* under normal / protanopic / deuteranopic
    vision. Returns a pandas DataFrame."""
    import pandas as pd

    colours = GENE_COLOURS if colours is None else colours
    rows = []
    for kind in ("normal", "protanopia", "deuteranopia"):
        rgb = {g: hex_to_rgb01(h) for g, h in colours.items()}
        if kind != "normal":
            rgb = {g: simulate_cvd(v, kind) for g, v in rgb.items()}
        lab = {g: rgb01_to_lab(v) for g, v in rgb.items()}
        for a, b in itertools.combinations(sorted(colours), 2):
            d = lab[a] - lab[b]
            rows.append({"vision": kind, "gene_a": a, "gene_b": b,
                         "delta_e76": float(np.sqrt((d ** 2).sum())),
                         "delta_L": float(abs(d[0]))})
    return pd.DataFrame(rows)


def simulate_image(src, dst_stub):
    """Write deuteranopia and protanopia simulations of a rendered PNG.
    Returns {kind: path}."""
    from pathlib import Path

    from PIL import Image

    src, dst_stub = Path(src), Path(dst_stub)
    im = np.asarray(Image.open(src).convert("RGB"), dtype=float) / 255.0
    out = {}
    for kind in ("deuteranopia", "protanopia"):
        sim = (simulate_cvd(im, kind) * 255.0).round().astype(np.uint8)
        path = dst_stub.with_name(f"{dst_stub.stem}_{kind}.png")
        Image.fromarray(sim).save(path)
        out[kind] = path
    return out
