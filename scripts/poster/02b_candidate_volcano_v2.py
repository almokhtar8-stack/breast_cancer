#!/usr/bin/env python3
"""Figure 02b (candidate, v2) -- two volcano variants of candidate significance
across the four transcriptomic contexts. post_freeze_exploratory: a PROPOSED
replacement for poster Figure 2, not a swap; the frozen Figure 2, its manifest,
and the v1 candidate figure are untouched.

Thin entry point: it imports the canonical implementation in `src/` and calls
it unchanged. Edit `src/poster_candidate_volcano_v2.py` if the figure needs to
change.

Usage:
    python scripts/poster/02b_candidate_volcano_v2.py
"""

from __future__ import annotations

import logging
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
os.chdir(ROOT)

# Reproducible PDF bytes, as scripts/poster/build_all.py pins it.
os.environ.setdefault("SOURCE_DATE_EPOCH", "1600000000")

from src.poster_candidate_volcano_v2 import main  # noqa: E402


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    for variant, files in main().items():
        for path in files.values():
            print(f"wrote {path}")
