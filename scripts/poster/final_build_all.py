#!/usr/bin/env python3
"""Render the complete final poster figure set: 7 figures x 3 formats,
manifest, colour-vision simulations, and the frozen-value verification table.

Usage:
    python scripts/poster/final_build_all.py
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

from src.poster_final_build import main  # noqa: E402

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    manifest = main()
    print(manifest[["figure", "smallest_printed_pt", "meets_20pt_floor"]].to_string(index=False))
