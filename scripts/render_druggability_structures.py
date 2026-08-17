#!/usr/bin/env python3
"""Renders the experimental candidate structures for the druggability poster
figure with PyMOL. Run manually / from scripts/; the figure module reads only
the PNG files this writes.

Structures are the ALREADY-DOWNLOADED, provenance-recorded local PDB files
(see each directory's PROVENANCE.txt); nothing is re-downloaded and no
docking, pocket detection, energy minimisation or pose prediction is
performed -- these are pure static renders of deposited experimental
coordinates.

  KDM1A  6NQU  chain A + ligand KWM  (GSK2879552, a selective inhibitor)
  TLK2   5O0Y  chain A + ligand AGS  (ATP-gamma-S, a substrate analog --
                                      NOT an inhibitor)
  USP34  7W3U  chain A + chain D ubiquitin + ligand AYE (propargylamide
                                      warhead of a covalent ACTIVITY-BASED
                                      PROBE, covalently linked to Cys1903)

VEZF1 is deliberately absent: no experimental structure exists in the
audited evidence, and no homology/AlphaFold model is substituted.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import yaml

RENDER_W, RENDER_H = 1600, 1250

# Neutral, scientifically restrained protein tones (candidate identity colors
# are used only for the small ligand/catalytic-residue highlights, never to
# recolor whole proteins).
PROTEIN_GREY = "0xB9BEC6"
PARTNER_GREY = "0xC7B79E"   # probe carrier (ubiquitin) -- distinct neutral tan so it reads as the probe, not the target protein

SPEC = {
    "KDM1A": dict(pdb="6NQU", dir_key="kdm1a_tlk2_structures_dir", protein_sel="chain A and polymer",
                  ligand_sel="chain A and resn KWM", ligand_color="0xD55E00", extra_polymer=None, catalytic=None),
    "TLK2": dict(pdb="5O0Y", dir_key="kdm1a_tlk2_structures_dir", protein_sel="chain A and polymer",
                 ligand_sel="chain A and resn AGS", ligand_color="0x2E86C6", extra_polymer=None, catalytic=None),
    "USP34": dict(pdb="7W3U", dir_key="usp34_structures_dir", protein_sel="chain A and polymer",
                  ligand_sel="chain A and resn AYE", ligand_color="0x0072B2", extra_polymer="chain D and polymer",
                  catalytic="chain A and resi 1903"),
}


def render(gene: str, cfg: dict, out_dir: Path) -> Path:
    from pymol import cmd

    spec = SPEC[gene]
    pdb_path = Path(cfg["data"]["raw"][spec["dir_key"]]) / f"{spec['pdb']}.pdb"
    if not pdb_path.exists():
        raise FileNotFoundError(f"{pdb_path} missing -- structures must already be downloaded locally")

    cmd.reinitialize()
    cmd.load(str(pdb_path), gene)
    cmd.remove("solvent")
    cmd.remove("resn SO4+GOL+EDO+PEG+ACT+CL+NA+MPD+DMS")

    keep = f"({spec['protein_sel']}) or ({spec['ligand_sel']})"
    if spec["extra_polymer"]:
        keep += f" or ({spec['extra_polymer']})"
    cmd.remove(f"not ({keep})")

    cmd.hide("everything")
    cmd.show("cartoon", spec["protein_sel"])
    cmd.color(PROTEIN_GREY, spec["protein_sel"])
    if spec["extra_polymer"]:
        cmd.show("cartoon", spec["extra_polymer"])
        cmd.color(PARTNER_GREY, spec["extra_polymer"])

    cmd.show("sticks", spec["ligand_sel"])
    cmd.color(spec["ligand_color"], spec["ligand_sel"])
    cmd.set("stick_radius", 0.32, spec["ligand_sel"])

    if spec["catalytic"]:
        cmd.show("sticks", spec["catalytic"])
        cmd.color(spec["ligand_color"], spec["catalytic"])
        cmd.set("stick_radius", 0.26, spec["catalytic"])

    cmd.set("cartoon_transparency", 0.18)
    cmd.set("ray_opaque_background", 0)
    cmd.set("antialias", 2)
    cmd.set("ray_trace_mode", 0)
    cmd.set("specular", 0.15)
    cmd.set("ambient", 0.28)
    cmd.set("direct", 0.55)
    cmd.set("cartoon_fancy_helices", 1)
    cmd.set("cartoon_smooth_loops", 1)
    cmd.bg_color("white")

    cmd.orient(spec["protein_sel"])
    cmd.zoom(keep, 1.0)
    cmd.turn("y", 12)

    out_dir.mkdir(parents=True, exist_ok=True)
    out_png = out_dir / f"{gene}_{spec['pdb']}.png"
    cmd.png(str(out_png), width=RENDER_W, height=RENDER_H, dpi=300, ray=1)
    print(f"{gene}: {spec['pdb']} -> {out_png}", file=sys.stderr)
    return out_png


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config/config.yaml")
    parser.add_argument("--out-dir", default="results/figures/poster_druggability_v1/renders")
    args = parser.parse_args()

    with open(args.config) as f:
        cfg = yaml.safe_load(f)

    out_dir = Path(args.out_dir)
    for gene in SPEC:
        render(gene, cfg, out_dir)
    print("VEZF1: no experimental structure in the audited evidence -- intentionally not rendered",
          file=sys.stderr)


if __name__ == "__main__":
    main()
