#!/usr/bin/env python3
"""Downloads STRING (string-db.org) interaction-partner data for a list of
gene symbols. Run manually / from scripts/; downstream systems-network
modules read only the local TSV files this writes, never the network at
analysis runtime.

Usage:
    python3 scripts/download_string_interactions.py --genes-file <path with one gene symbol per line> \
        --out <output tsv path> [--required-score 700] [--species 9606]
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import requests

STRING_API_BASE = "https://string-db.org/api/tsv"


def fetch_partners(genes: list[str], species: int, required_score: int, endpoint: str, network_type: str = "functional") -> str:
    identifiers = "%0d".join(genes)
    params = {"identifiers": identifiers, "species": species, "required_score": required_score, "network_type": network_type}
    resp = requests.get(f"{STRING_API_BASE}/{endpoint}", params=params, timeout=60)
    resp.raise_for_status()
    return resp.text


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--genes-file", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--species", type=int, default=9606)
    parser.add_argument("--required-score", type=int, default=700)
    parser.add_argument("--batch-size", type=int, default=50)
    parser.add_argument("--endpoint", choices=["interaction_partners", "network"], default="interaction_partners")
    parser.add_argument("--network-type", choices=["functional", "physical"], default="functional")
    args = parser.parse_args()

    genes = [g.strip() for g in Path(args.genes_file).read_text().splitlines() if g.strip()]
    genes = sorted(set(genes))
    print(f"Querying STRING {args.endpoint} ({args.network_type}) for {len(genes)} genes (required_score={args.required_score})", file=sys.stderr)

    all_lines: list[str] = []
    header: str | None = None
    batch_size = args.batch_size if args.endpoint == "interaction_partners" else len(genes)
    for i in range(0, len(genes), batch_size):
        batch = genes[i : i + batch_size]
        text = fetch_partners(batch, args.species, args.required_score, args.endpoint, args.network_type)
        lines = text.strip().split("\n")
        if not lines or not lines[0]:
            continue
        if header is None:
            header = lines[0]
            all_lines.append(header)
        all_lines.extend(lines[1:])
        print(f"  batch {i // args.batch_size + 1}: {len(lines) - 1} interaction rows", file=sys.stderr)
        time.sleep(1)

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out).write_text("\n".join(all_lines) + "\n")
    print(f"Wrote {len(all_lines) - 1} rows -> {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
