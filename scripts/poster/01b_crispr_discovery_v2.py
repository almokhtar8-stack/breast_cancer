#!/usr/bin/env python3
"""Figure 01b -- genome-scale CRISPR discovery, v2 (colour change only).

post_freeze_exploratory. Thin entry point: it imports the canonical
implementation in `src/` and calls it unchanged. No scientific logic or value
is duplicated here -- edit `src/poster_crispr_discovery_v2.py` if the figure
needs to change.

This does NOT replace figure 01. `scripts/poster/01_crispr_discovery.py` and
`src/poster_crispr_discovery_v1.py` are untouched and still build the
committed poster figure.

Usage:
    python scripts/poster/01b_crispr_discovery_v2.py
"""

from __future__ import annotations

import logging
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
# The canonical src/ implementations read and write repo-relative paths, so
# run from the repository root regardless of the caller's working directory.
os.chdir(ROOT)

# Pinned so the PDF is byte-reproducible.
os.environ.setdefault("SOURCE_DATE_EPOCH", "1600000000")

from src.poster_crispr_discovery_v2 import build_crispr_discovery_main, OUT_DIR  # noqa: E402


def main() -> Path:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    stub = OUT_DIR / "CRISPR_discovery_v2"
    build_crispr_discovery_main(stub)
    print(f"wrote {stub}.png / .pdf / .svg")
    return stub


if __name__ == "__main__":
    main()
