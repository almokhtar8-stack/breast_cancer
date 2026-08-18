"""Gene colours for poster figures -- the single future source of truth.

post_freeze_exploratory. Created so that a colour change cannot silently
recolour a committed figure. The previous arrangement defined the four
candidate colours in `src/post_audit_sensitivity_visualization.py` and imported
them into the poster renderers, which meant editing that one dict would have
changed every figure that already exists. That module is deliberately left
untouched; new renderers import from here instead.

WHY THESE FOUR VALUES. The previous set was
``{"KDM1A": "#D55E00", "TLK2": "#CC79A7", "USP34": "#0072B2", "VEZF1": "#E69F00"}``
-- two oranges (#D55E00, #E69F00) and, at small mark sizes, a blue that pairs
visually with the pink-purple (#0072B2 with #CC79A7). The four genes therefore
read as two similar pairs rather than four distinct identities. The set below
gives four clearly separated hues. Three (#D55E00, #009E73, #56B4E9) are from
the Okabe-Ito colour-vision-safe palette; the purple #6A3D9A is a deliberate
addition, because Okabe-Ito contains no true purple.

Purple and light blue are the pair at risk under red-green colour-vision
deficiency: they converge in hue and must separate on LIGHTNESS instead. That
is not asserted here, it is computed -- see `palette_cvd_report()` below, which
applies the Machado, Oliveira & Fernandes (2009) severity-1.0 matrices in
linear sRGB. No colour-vision library exists in this environment and none is
added.
"""

from __future__ import annotations

import itertools

import numpy as np

GENE_COLOURS: dict[str, str] = {
    "KDM1A": "#D55E00",
    "TLK2": "#009E73",
    "USP34": "#6A3D9A",
    "VEZF1": "#56B4E9",
}

# The set this replaces, kept for the record so the change is auditable and so
# a diff report can state exactly which values moved.
PREVIOUS_GENE_COLOURS: dict[str, str] = {
    "KDM1A": "#D55E00",
    "TLK2": "#CC79A7",
    "USP34": "#0072B2",
    "VEZF1": "#E69F00",
}

# ---------------------------------------------------------------------------
# Colour-vision simulation: Machado, Oliveira & Fernandes (2009), IEEE TVCG
# 15(6), severity-1.0 matrices, applied in LINEAR sRGB.
# ---------------------------------------------------------------------------
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
    """rgb01 (..., 3) in [0, 1] -> simulated sRGB in [0, 1]."""
    return _linear_to_srgb(_srgb_to_linear(rgb01) @ CVD_MATRICES[kind].T)


def rgb01_to_lab(rgb01: np.ndarray) -> np.ndarray:
    lin = _srgb_to_linear(rgb01)
    m = np.array([[0.4124564, 0.3575761, 0.1804375],
                  [0.2126729, 0.7151522, 0.0721750],
                  [0.0193339, 0.1191920, 0.9503041]])
    xyz = (lin @ m.T) / np.array([0.95047, 1.0, 1.08883])   # D65 white
    eps, kappa = 216 / 24389, 24389 / 27
    f = np.where(xyz > eps, np.cbrt(xyz), (kappa * xyz + 16) / 116)
    return np.stack([116 * f[..., 1] - 16,
                     500 * (f[..., 0] - f[..., 1]),
                     200 * (f[..., 1] - f[..., 2])], axis=-1)


def palette_cvd_report(colours: dict[str, str] | None = None):
    """Pairwise CIE76 dE and dL* between gene colours under normal,
    protanopic and deuteranopic vision. Returns a pandas DataFrame."""
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
