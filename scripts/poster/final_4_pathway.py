#!/usr/bin/env python3
"""Poster final figure 4 (pathway) -- thin entry point.

Imports the canonical implementation in src/ and calls it unchanged. Edit
src/poster_pathway_final.py if the figure needs to change.
"""

from __future__ import annotations

import logging
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
os.chdir(ROOT)
os.environ.setdefault("SOURCE_DATE_EPOCH", "1600000000")

from src.poster_pathway_final import main  # noqa: E402

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    written, _ = main()
    for path in written.values():
        print(f"wrote {path}")
