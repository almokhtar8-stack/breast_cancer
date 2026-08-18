#!/usr/bin/env python3
"""Figure 02b (candidate) -- volcano view of candidate significance across the
four transcriptomic contexts. post_freeze_exploratory: a PROPOSED replacement
for poster Figure 2, not a swap; the frozen Figure 2 and its manifest are
untouched.

Thin entry point: it imports the canonical implementation in `src/` and calls
it unchanged. No scientific logic or value is duplicated here -- edit
`src/poster_candidate_volcano_v1.py` if the figure needs to change.

Usage:
    python scripts/poster/02b_candidate_volcano.py
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

# Reproducible PDF bytes, as scripts/poster/build_all.py pins it.
os.environ.setdefault("SOURCE_DATE_EPOCH", "1600000000")

from src.poster_candidate_volcano_v1 import main  # noqa: E402


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    written = main()
    for ext, path in written.items():
        print(f"wrote {path}")
