"""PyMOL rendering script for the poster's USP34 structural-biology figure.

Run standalone via PyMOL's command-line interpreter (not a normal Python
import target, since it uses PyMOL's `cmd` API):

    pymol -cq src/poster_structure_render.py -- <apo_pdb> <bound_pdb> <out_dir>

Renders four real, ray-traced panels from the two already-frozen,
already-verified USP34 PDB structures (7W3R apo, 7W3U covalent-probe-bound)
-- no new structural analysis, pocket detection, or docking is performed;
this only visualizes coordinates already established in the
final_translational phase.

Panels:
  A. 7W3R cartoon, catalytic Cys1903/His2164 shown as sticks.
  B. 7W3U cartoon in a matched orientation, ubiquitin-probe chain shown.
  C. Catalytic-cleft close-up on the bound structure (sticks + surface).
  D. Apo-vs-bound catalytic-domain overlay (structural alignment).
"""

from __future__ import annotations

import sys

from pymol import cmd

CATALYTIC_CYS = 1903
CATALYTIC_HIS = 2164

# Okabe-Ito-derived palette matching the rest of the poster set.
COLOR_APO = "0x8C8C8C"
COLOR_BOUND = "0x0072B2"
COLOR_UB = "0x009E73"
COLOR_CYS = "0xb0392f"
COLOR_HIS = "0xE69F00"


def _base_scene(obj: str) -> None:
    cmd.bg_color("white")
    cmd.hide("everything")
    cmd.show("cartoon", obj)
    cmd.set("cartoon_transparency", 0.0)
    cmd.set("ray_shadows", 0)
    cmd.set("ambient", 0.4)
    cmd.set("antialias", 2)
    cmd.set("ray_trace_mode", 0)


def _highlight_catalytic(obj: str, chain: str = "A") -> None:
    cys_sel = f"{obj} and chain {chain} and resi {CATALYTIC_CYS}"
    his_sel = f"{obj} and chain {chain} and resi {CATALYTIC_HIS}"
    cmd.show("sticks", cys_sel)
    cmd.show("sticks", his_sel)
    cmd.color(COLOR_CYS, cys_sel)
    cmd.color(COLOR_HIS, his_sel)
    cmd.set("stick_radius", 0.25, cys_sel)
    cmd.set("stick_radius", 0.25, his_sel)


def render(apo_pdb: str, bound_pdb: str, out_dir: str, width: int = 1600, height: int = 1600, dpi: int = 300) -> None:
    cmd.reinitialize()
    cmd.load(apo_pdb, "apo")
    cmd.load(bound_pdb, "bound")
    cmd.remove("solvent")
    cmd.remove("apo and not chain A")
    cmd.remove("bound and not (chain A+D)")

    # ---- Panel A: apo overall view ----
    _base_scene("apo")
    cmd.color(COLOR_APO, "apo")
    cmd.util.cnc("apo and resi %d+%d" % (CATALYTIC_CYS, CATALYTIC_HIS))
    _highlight_catalytic("apo")
    cmd.orient("apo")
    cmd.zoom("apo", buffer=3)
    cmd.ray(width, height)
    cmd.png(f"{out_dir}/panel_A_apo_overview.png", dpi=dpi)

    # ---- Panel B: bound overall view (chain A + ubiquitin chain D), matched orientation ----
    cmd.create("bound_domain", "bound and chain A")
    cmd.align("bound_domain", "apo")
    cmd.matrix_copy("bound_domain", "bound")
    _base_scene("bound")
    cmd.color(COLOR_BOUND, "bound and chain A")
    cmd.color(COLOR_UB, "bound and chain D")
    cmd.show("cartoon", "bound and chain D")
    _highlight_catalytic("bound", chain="A")
    cmd.show("sticks", "bound and resn AYE")
    cmd.color(COLOR_UB, "bound and resn AYE")
    cmd.orient("apo")  # reuse apo's matched orientation
    cmd.zoom("bound", buffer=3)
    cmd.ray(width, height)
    cmd.png(f"{out_dir}/panel_B_bound_overview.png", dpi=dpi)

    # ---- Panel C: catalytic-cleft close-up on the bound structure ----
    # No native PyMOL labels here (they visually collide with the sticks at
    # this zoom level) -- residue labels are added later in the matplotlib
    # composite step, where label placement can be controlled precisely.
    cmd.hide("everything", "bound and chain D")
    cmd.set("cartoon_transparency", 0.55, "bound and chain A")
    cleft_sel = f"bound and chain A and resi {CATALYTIC_CYS}+{CATALYTIC_HIS}"
    cmd.zoom(cleft_sel, buffer=8)
    cmd.show("sticks", "bound and resn AYE")
    cmd.set("stick_radius", 0.22, "bound and resn AYE")
    cmd.ray(width, height)
    cmd.png(f"{out_dir}/panel_C_cleft_closeup.png", dpi=dpi)
    cmd.set("cartoon_transparency", 0.0, "bound and chain A")

    # ---- Panel D: apo vs bound catalytic-domain overlay ----
    cmd.show("cartoon", "apo")
    cmd.show("cartoon", "bound and chain A")
    cmd.hide("everything", "bound and chain D")
    cmd.hide("sticks", "bound and resn AYE")
    cmd.color(COLOR_APO, "apo")
    cmd.color(COLOR_BOUND, "bound and chain A")
    _highlight_catalytic("apo")
    _highlight_catalytic("bound", chain="A")
    cmd.set("cartoon_transparency", 0.35, "apo")
    cmd.zoom(f"apo and resi {CATALYTIC_CYS}+{CATALYTIC_HIS}", buffer=10)
    cmd.ray(width, height)
    cmd.png(f"{out_dir}/panel_D_apo_vs_bound_overlay.png", dpi=dpi)

    rms = cmd.align("bound_domain", "apo")[0]
    with open(f"{out_dir}/alignment_rmsd.txt", "w") as f:
        f.write(f"apo (7W3R) vs bound (7W3U) chain-A catalytic-domain alignment RMSD: {rms:.3f} A\n")


# PyMOL executes this file with __name__ == "pymol", not "__main__" --
# guard against import-time execution only.
if __name__ in ("__main__", "pymol"):
    argv = sys.argv[1:]
    if "--" in argv:
        argv = argv[argv.index("--") + 1:]
    if len(argv) >= 3:
        apo_pdb, bound_pdb, out_dir = argv[0], argv[1], argv[2]
        render(apo_pdb, bound_pdb, out_dir)
