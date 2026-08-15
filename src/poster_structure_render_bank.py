"""PyMOL rendering script for the FIGURE BANK's two USP34 structural-biology
candidates (Step 12 / 12b of the figure-bank request).

Run standalone via PyMOL's command-line interpreter:

    pymol -cq src/poster_structure_render_bank.py -- <apo_pdb> <bound_pdb> <out_dir>

Renders five real, ray-traced panels from the two already-frozen,
already-verified USP34 PDB structures (7W3R apo, 7W3U covalent-probe-bound)
-- no new structural analysis, pocket detection, or docking is performed.

Candidate 12 (hero + pocket inset):
  hero_main.png   -- 7W3R, translucent surface + cartoon, catalytic cleft
                      patch highlighted, Cys1903/His2164 as sticks
  hero_inset.png  -- 7W3U catalytic-cleft close-up with covalent probe

Candidate 12b (matched apo-vs-bound comparison):
  comp_apo.png     -- 7W3R, large cartoon, catalytic sticks
  comp_bound.png   -- 7W3U, matched orientation, ubiquitin probe visible
  comp_closeup.png -- catalytic-cleft close-up on the bound structure
"""

from __future__ import annotations

import sys

from pymol import cmd

CATALYTIC_CYS = 1903
CATALYTIC_HIS = 2164

COLOR_APO = "0x8C8C8C"
COLOR_BOUND = "0x0072B2"
COLOR_UB = "0x009E73"
COLOR_CYS = "0xb0392f"
COLOR_HIS = "0xE69F00"
COLOR_SURFACE = "0xD8DCE0"
COLOR_CLEFT_PATCH = "0xAFC6E0"


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


def render(apo_pdb: str, bound_pdb: str, out_dir: str, width: int = 2000, height: int = 2000, dpi: int = 300) -> None:
    cmd.reinitialize()
    cmd.load(apo_pdb, "apo")
    cmd.load(bound_pdb, "bound")
    cmd.remove("solvent")
    cmd.remove("apo and not chain A")
    cmd.remove("bound and not (chain A+D)")

    # ---- shared alignment: bound chain A onto apo, reused for every panel ----
    cmd.create("bound_domain", "bound and chain A")
    cmd.align("bound_domain", "apo")
    cmd.matrix_copy("bound_domain", "bound")

    # =====================================================================
    # Candidate 12 -- hero surface + pocket inset
    # =====================================================================

    _base_scene("apo")
    cmd.show("surface", "apo")
    cmd.set("transparency", 0.45, "apo")
    cmd.color(COLOR_SURFACE, "apo")
    cleft_sel = f"apo and byres (resi {CATALYTIC_CYS}+{CATALYTIC_HIS} around 9)"
    cmd.color(COLOR_CLEFT_PATCH, cleft_sel)
    cmd.set("cartoon_transparency", 0.25, "apo")
    _highlight_catalytic("apo")
    cmd.orient("apo")
    cmd.turn("y", 15)
    cmd.zoom("apo", buffer=4)
    cmd.ray(width, height)
    cmd.png(f"{out_dir}/hero_main.png", dpi=dpi)

    _base_scene("bound and chain A")
    cmd.hide("everything", "bound and chain D")
    cmd.set("cartoon_transparency", 0.55, "bound and chain A")
    cmd.color(COLOR_BOUND, "bound and chain A")
    _highlight_catalytic("bound", chain="A")
    cmd.show("sticks", "bound and resn AYE")
    cmd.set("stick_radius", 0.24, "bound and resn AYE")
    cmd.color(COLOR_UB, "bound and resn AYE")
    cleft_close = f"bound and chain A and resi {CATALYTIC_CYS}+{CATALYTIC_HIS}"
    cmd.zoom(cleft_close, buffer=7)
    cmd.ray(width, height)
    cmd.png(f"{out_dir}/hero_inset.png", dpi=dpi)
    cmd.set("cartoon_transparency", 0.0, "bound and chain A")

    # =====================================================================
    # Candidate 12b -- matched apo-vs-bound comparison
    # =====================================================================

    _base_scene("apo")
    cmd.color(COLOR_APO, "apo")
    _highlight_catalytic("apo")
    cmd.orient("apo")
    cmd.zoom("apo", buffer=3)
    cmd.ray(width, height)
    cmd.png(f"{out_dir}/comp_apo.png", dpi=dpi)

    _base_scene("bound")
    cmd.show("cartoon", "bound and chain D")
    cmd.color(COLOR_BOUND, "bound and chain A")
    cmd.color(COLOR_UB, "bound and chain D")
    _highlight_catalytic("bound", chain="A")
    cmd.show("sticks", "bound and resn AYE")
    cmd.set("stick_radius", 0.24, "bound and resn AYE")
    cmd.color(COLOR_UB, "bound and resn AYE")
    cmd.orient("apo")
    cmd.zoom("bound", buffer=3)
    cmd.ray(width, height)
    cmd.png(f"{out_dir}/comp_bound.png", dpi=dpi)

    cmd.hide("everything", "bound and chain D")
    cmd.set("cartoon_transparency", 0.5, "bound and chain A")
    cmd.zoom(cleft_close, buffer=8)
    cmd.ray(width, height)
    cmd.png(f"{out_dir}/comp_closeup.png", dpi=dpi)

    rms = cmd.align("bound_domain", "apo")[0]
    with open(f"{out_dir}/alignment_rmsd.txt", "w") as f:
        f.write(f"apo (7W3R) vs bound (7W3U) chain-A catalytic-domain alignment RMSD: {rms:.3f} A\n")


# PyMOL executes this file with __name__ == "pymol", not "__main__".
if __name__ in ("__main__", "pymol"):
    argv = sys.argv[1:]
    if "--" in argv:
        argv = argv[argv.index("--") + 1:]
    if len(argv) >= 3:
        apo_pdb, bound_pdb, out_dir = argv[0], argv[1], argv[2]
        render(apo_pdb, bound_pdb, out_dir)
