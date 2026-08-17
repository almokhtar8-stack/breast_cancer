#!/usr/bin/env python3
"""Downloads a consistent, two-level STRING (string-db.org) network for the
CURRENT four poster focus genes (KDM1A, TLK2, USP34, VEZF1), for the
post-freeze exploratory network/mechanism figure v2. Frozen candidate
rankings and all prior scientific results are untouched by this script --
it only writes new local interaction tables.

Methodology matches the exact parameters already used by the project's
earlier network build (`scripts/download_string_interactions.py`,
`src/systems_network_build.py`): STRING `interaction_partners` endpoint,
required_score=700 (STRING's own "high confidence" band, 0-1000 scale),
species=9606 (human), network_type=functional for the primary edge/score
set and network_type=physical (queried on the exact same identifier set)
used only to label which functional edges are also physical_PPI -- never
to add new edges the functional query didn't already return.

Two levels, same rule applied identically to all four candidates:

  Level 1: `interaction_partners` for the four candidate genes together
      (one batched call per network_type) -- returns EVERY partner above
      required_score=700 for each gene independently (not a shared/
      competitive pool), so a hub-heavy candidate cannot crowd out a
      sparse one.
  Level 2: `interaction_partners` for the pooled, capped Level-1 partner
      set (top-12 per candidate by score, deterministic tie-break
      alphabetical -- capping happens in the analysis module, this script
      queries every Level-1 partner returned so the analysis module has
      the full picture to cap and filter from) -- reveals real bridges
      between candidates' neighborhoods and connections to canonical
      endocrine/EMT/WNT/cell-cycle genes, filtered for relevance in the
      analysis module, never invented here.

Run manually / from scripts/; `src/poster_network_mechanism_v2.py` reads
only the local TSV files this writes, never the network at render time.
"""

from __future__ import annotations

import csv
import sys
from io import StringIO
from pathlib import Path

import requests

STRING_API_BASE = "https://string-db.org/api/tsv"
CANDIDATES = ["KDM1A", "TLK2", "USP34", "VEZF1"]
SPECIES = 9606
REQUIRED_SCORE = 700
OUT_DIR = Path("data/reference/interactions")


def fetch_interaction_partners(genes: list[str], network_type: str) -> str:
    identifiers = "%0d".join(genes)
    resp = requests.get(
        f"{STRING_API_BASE}/interaction_partners",
        params={"identifiers": identifiers, "species": SPECIES, "required_score": REQUIRED_SCORE, "network_type": network_type},
        timeout=90,
    )
    resp.raise_for_status()
    return resp.text


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    print(f"Level 1: querying STRING interaction_partners for {CANDIDATES} "
          f"(required_score={REQUIRED_SCORE}, species={SPECIES})", file=sys.stderr)
    level1_functional = fetch_interaction_partners(CANDIDATES, "functional")
    (OUT_DIR / "string_v2_level1_functional.tsv").write_text(level1_functional)
    level1_physical = fetch_interaction_partners(CANDIDATES, "physical")
    (OUT_DIR / "string_v2_level1_physical.tsv").write_text(level1_physical)

    rows = list(csv.DictReader(StringIO(level1_functional), delimiter="\t"))
    counts = {c: sum(1 for r in rows if r["preferredName_A"] == c) for c in CANDIDATES}
    print(f"  Level-1 partner counts (uncapped): {counts}", file=sys.stderr)

    level1_partners = sorted({r["preferredName_B"] for r in rows})
    print(f"Level 2: querying STRING interaction_partners for the pooled "
          f"Level-1 partner set (n={len(level1_partners)})", file=sys.stderr)
    level2_functional = fetch_interaction_partners(level1_partners, "functional")
    (OUT_DIR / "string_v2_level2_functional.tsv").write_text(level2_functional)
    level2_physical = fetch_interaction_partners(level1_partners, "physical")
    (OUT_DIR / "string_v2_level2_physical.tsv").write_text(level2_physical)

    l2_rows = list(csv.DictReader(StringIO(level2_functional), delimiter="\t"))
    print(f"  Level-2 raw edge rows: {len(l2_rows)}", file=sys.stderr)
    print("Done. Wrote 4 files to", OUT_DIR, file=sys.stderr)


if __name__ == "__main__":
    main()
