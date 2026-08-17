"""PyMOL rendering script for the exploration-v2 structural comparison
figures (Section G). Renders KDM1A (6NQU, inhibitor-bound) and TLK2 (5O0Y,
ATP-analog-bound) with the SAME camera philosophy, background, and
ray-tracing settings already established for USP34 in
src/poster_structure_render_bank.py, so the three can be visually compared
side by side.

Run standalone via PyMOL's command-line interpreter:

    pymol -cq src/poster_exploration_v2_structure.py -- <6nqu_pdb> <5o0y_pdb> <out_dir>

No new structural analysis, pocket detection, or docking is performed --
these are the same two real, already-fetched PDB structures documented in
results/reports/poster_exploration_v2/DATA_FOR_VISUALIZATION_AUDIT.md,
rendered here for visualization only.
"""

from __future__ import annotations

import sys

from pymol import cmd

COLOR_KDM1A = "0xD55E00"
COLOR_TLK2 = "0x56B4E9"
COLOR_LIGAND_INHIBITOR = "0x009E73"
COLOR_LIGAND_ATP = "0xE69F00"
COLOR_SURFACE = "0xD8DCE0"


def _base_scene(obj: str) -> None:
    cmd.bg_color("white")
    cmd.hide("everything")
    cmd.show("cartoon", obj)
    cmd.set("cartoon_transparency", 0.0)
    cmd.set("ray_shadows", 0)
    cmd.set("ambient", 0.4)
    cmd.set("antialias", 2)
    cmd.set("ray_trace_mode", 0)


def render(kdm1a_pdb: str, tlk2_pdb: str, out_dir: str, width: int = 1800, height: int = 1800, dpi: int = 300) -> None:
    cmd.reinitialize()

    # ---- KDM1A / LSD1, 6NQU, inhibitor KWM ----
    cmd.load(kdm1a_pdb, "kdm1a")
    cmd.remove("kdm1a and not chain A")
    cmd.remove("solvent")
    # LSD1's catalytic (SWIRM+AOL) domain is compact and globular; chain A
    # also carries a long, separate Tower-domain coiled-coil helix that
    # normally packs against a CoREST partner (absent here, monomer only)
    # -- with no partner to pack against it projects ~100 A away from the
    # catalytic core. Showing it in full would shrink the catalytic/
    # ligand-bound region (the actual point of this figure) to a speck.
    # The "full" hero view is therefore restricted to the compact globular
    # core around the ligand (a standard, disclosed cropping choice, not a
    # scientific misrepresentation -- the excluded Tower helix is real and
    # is documented here, not hidden); the untrimmed structure is
    # available in the raw PDB file for anyone who wants the full extent.
    cmd.select("kdm1a_core", "byres (kdm1a within 40 of resn KWM)")
    _base_scene("kdm1a_core")
    cmd.color(COLOR_KDM1A, "kdm1a_core and polymer")
    cmd.show("sticks", "kdm1a_core and resn KWM")
    cmd.color(COLOR_LIGAND_INHIBITOR, "kdm1a_core and resn KWM")
    cmd.set("stick_radius", 0.25, "kdm1a_core and resn KWM")
    cmd.orient("kdm1a_core")
    cmd.zoom("kdm1a_core", buffer=5, complete=1)
    cmd.ray(width, height)
    cmd.png(f"{out_dir}/kdm1a_full.png", dpi=dpi)

    cmd.zoom("kdm1a and resn KWM", buffer=8)
    cmd.ray(width, height)
    cmd.png(f"{out_dir}/kdm1a_pocket.png", dpi=dpi)

    # ---- TLK2, 5O0Y, ATP analog AGS ----
    cmd.reinitialize()
    cmd.load(tlk2_pdb, "tlk2")
    cmd.remove("tlk2 and not chain A")
    cmd.remove("solvent")
    _base_scene("tlk2")
    cmd.color(COLOR_TLK2, "tlk2 and polymer")
    cmd.show("sticks", "tlk2 and resn AGS")
    cmd.color(COLOR_LIGAND_ATP, "tlk2 and resn AGS")
    cmd.set("stick_radius", 0.25, "tlk2 and resn AGS")
    cmd.orient("tlk2")
    cmd.zoom("tlk2", buffer=4, complete=1)
    cmd.ray(width, height)
    cmd.png(f"{out_dir}/tlk2_full.png", dpi=dpi)

    cmd.zoom("tlk2 and resn AGS", buffer=8)
    cmd.ray(width, height)
    cmd.png(f"{out_dir}/tlk2_pocket.png", dpi=dpi)


# PyMOL executes this file with __name__ == "pymol", not "__main__".
if __name__ in ("__main__", "pymol"):
    argv = sys.argv[1:]
    if "--" in argv:
        argv = argv[argv.index("--") + 1:]
    if len(argv) >= 3:
        kdm1a_pdb, tlk2_pdb, out_dir = argv[0], argv[1], argv[2]
        render(kdm1a_pdb, tlk2_pdb, out_dir)
