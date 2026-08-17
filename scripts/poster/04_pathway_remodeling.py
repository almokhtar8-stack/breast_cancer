#!/usr/bin/env python3
"""Figure 04 -- network-relevant pathway remodeling from frozen GSEA results.

Thin entry point: it imports the canonical implementation in `src/` and calls
it unchanged. No scientific logic or value is duplicated here -- edit
`src/poster_pathway_v2.py` if the science needs to change.

Usage:
    python scripts/poster/04_pathway_remodeling.py
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

from src.poster_pathway_v2 import build_pathway_v2, OUT_DIR  # noqa: E402


def main() -> Path:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    stub = OUT_DIR / "PATHWAY_v2"
    build_pathway_v2(stub)
    print(f"wrote {stub}.png / .pdf / .svg")
    return stub


if __name__ == "__main__":
    main()
