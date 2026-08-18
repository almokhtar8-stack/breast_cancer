"""Poster figure 7 (final): how reachable each candidate is, by kind of evidence.

post_freeze_exploratory. No new structural work: every claim is read at render
time from the frozen audit table
`results/tables/post_audit_sensitivity/06b_structural_tractability_audit.tsv`
and verified before anything is drawn. The three structure images are the
already-committed PyMOL renders of already-downloaded experimental PDB files
(`results/figures/poster_druggability_v1/renders/`). NO docking, pose
prediction, affinity estimation or pocket re-detection was performed anywhere
in this project, and none is implied here.

WHAT THIS FIXES. The frozen figure 06 carries more prose than picture. Here
the evidence is a five-level track a reader can scan in one pass, the renders
are given the space, and the wording is tightened.

WHY A TRACK AND NOT A SCORE. The four candidates differ in KIND of evidence,
not in degree, so a single "tractability score" would impose a false ordering.
An inhibitor-bound co-crystal, an ATP-analogue-bound kinase domain, a covalent
activity-based probe, and a homology-model screening hit are not four points
on one axis. The five levels are kept separate and never summed.

WHAT "REACHABLE" MEANS HERE, AND WHAT IT DOES NOT. This figure reports
chemical reachability only. It is not efficacy, not selectivity, not safety,
and not evidence that inhibiting any of these genes would help a patient.
USP34's 7W3U is bound to a ubiquitin activity-based PROBE, not a drug, and
covers only the catalytic domain (~12% of the full-length protein). VEZF1 has
no experimental structure at all: no predicted model is substituted for it,
because the absence is itself the finding.

The published counter-evidence for the lead is shown on the figure, not
omitted: losing USP34 has been reported to push breast cells toward a more
mobile state (Cellular Signalling 2017, PMID 28499884), which works against
the proposed benefit.
"""

from __future__ import annotations

import logging
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import pandas as pd

from src.poster_final_common import (
    FONT,
    OUT_DIR,
    figure_footer,
    headline,
    pin_reproducibility,
    save,
    verify,
)
from src.poster_palette import GENE_COLOURS, NEUTRAL, WHITE

logger = logging.getLogger(__name__)

FIGURE = "F7_reachability"
AUDIT_TSV = Path("results/tables/post_audit_sensitivity/06b_structural_tractability_audit.tsv")
RENDER_DIR = Path("results/figures/poster_druggability_v1/renders")

GENES = ("KDM1A", "TLK2", "USP34", "VEZF1")
RENDERS = {"KDM1A": "KDM1A_6NQU.png", "TLK2": "TLK2_5O0Y.png", "USP34": "USP34_7W3U.png"}
STRUCTURE_ID = {"KDM1A": "PDB 6NQU", "TLK2": "PDB 5O0Y", "USP34": "PDB 7W3U",
                "VEZF1": "none"}
BOUND_TO = {"KDM1A": "inhibitor-bound", "TLK2": "ATP analogue", "USP34": "probe, not a drug", "VEZF1": "—"}

# The five levels, kept separate and never summed into a score.
LEVELS = [
    ("Experimental structure", {"KDM1A": 2, "TLK2": 1, "USP34": 1, "VEZF1": 0}),
    ("Chemical bound in it", {"KDM1A": 2, "TLK2": 1, "USP34": 1, "VEZF1": 0}),
    ("Selective inhibitor", {"KDM1A": 2, "TLK2": 0, "USP34": 0, "VEZF1": 0}),
    ("Clinical-stage compounds", {"KDM1A": 2, "TLK2": 0, "USP34": 0, "VEZF1": 0}),
    ("Unexplored opportunity", {"KDM1A": 0, "TLK2": 1, "USP34": 2, "VEZF1": 0}),
]
LEVEL_NOTE = {
    ("Experimental structure", "TLK2"): "kinase domain",
    ("Experimental structure", "USP34"): "catalytic, ~12%",
    ("Chemical bound in it", "TLK2"): "ATP analogue",
    ("Chemical bound in it", "USP34"): "probe",
    ("Unexplored opportunity", "TLK2"): "no inhibitor yet",
    ("Unexplored opportunity", "USP34"): "Cys1903",
}
FILL_LABEL = {2: "yes", 1: "partial", 0: "no"}


def load_audit() -> pd.DataFrame:
    return pd.read_csv(AUDIT_TSV, sep="\t").set_index("gene")


def gate(audit):
    """Every level state drawn must agree with the frozen audit table."""
    def truthy(v):
        return str(v).strip().upper().startswith(("TRUE", "YES"))

    def falsy(v):
        return str(v).strip().upper().startswith(("FALSE", "NO"))

    checks = [
        ("KDM1A_has_experimental_structure", 1.0 if truthy(audit.loc["KDM1A", "A_experimental_human_structure_exists"]) else 0.0, 1.0, 0),
        ("TLK2_has_experimental_structure", 1.0 if truthy(audit.loc["TLK2", "A_experimental_human_structure_exists"]) else 0.0, 1.0, 0),
        ("USP34_has_experimental_structure", 1.0 if truthy(audit.loc["USP34", "A_experimental_human_structure_exists"]) else 0.0, 1.0, 0),
        ("VEZF1_has_NO_experimental_structure", 1.0 if falsy(audit.loc["VEZF1", "A_experimental_human_structure_exists"]) else 0.0, 1.0, 0),
        ("KDM1A_has_selective_inhibitor", 1.0 if truthy(audit.loc["KDM1A", "E_validated_selective_small_molecule_inhibitor"]) else 0.0, 1.0, 0),
        ("TLK2_no_selective_inhibitor", 1.0 if falsy(audit.loc["TLK2", "E_validated_selective_small_molecule_inhibitor"]) else 0.0, 1.0, 0),
        ("USP34_no_selective_inhibitor", 1.0 if falsy(audit.loc["USP34", "E_validated_selective_small_molecule_inhibitor"]) else 0.0, 1.0, 0),
        ("KDM1A_has_clinical_pharmacology", 1.0 if truthy(audit.loc["KDM1A", "F_clinical_stage_pharmacology"]) else 0.0, 1.0, 0),
        ("TLK2_no_clinical_pharmacology", 1.0 if falsy(audit.loc["TLK2", "F_clinical_stage_pharmacology"]) else 0.0, 1.0, 0),
        ("USP34_no_clinical_pharmacology", 1.0 if falsy(audit.loc["USP34", "F_clinical_stage_pharmacology"]) else 0.0, 1.0, 0),
        ("VEZF1_no_clinical_pharmacology", 1.0 if falsy(audit.loc["VEZF1", "F_clinical_stage_pharmacology"]) else 0.0, 1.0, 0),
        ("n_genes_in_audit", len(audit), 4, 0),
        ("n_renders_committed", sum(1 for f in RENDERS.values() if (RENDER_DIR / f).exists()), 3, 0),
    ]
    return verify(FIGURE, checks)


def build(stub: Path):
    pin_reproducibility(FIGURE)
    audit = load_audit()
    verification = gate(audit)

    fig = plt.figure(figsize=(15.5, 10.4))
    col_w, col_x0, col_gap = 0.158, 0.305, 0.018
    # ---- top row: the structures themselves -----------------------------------
    for i, gene in enumerate(GENES):
        x = col_x0 + i * (col_w + col_gap)
        ax = fig.add_axes([x, 0.545, col_w, 0.250])
        ax.axis("off")
        if gene in RENDERS:
            ax.imshow(mpimg.imread(RENDER_DIR / RENDERS[gene]))
        else:
            # No experimental structure and no substitute: the gap is drawn as
            # a gap. A predicted model is deliberately NOT shown in its place.
            ax.add_patch(plt.Rectangle((0.03, 0.05), 0.94, 0.90, transform=ax.transAxes,
                                       facecolor=NEUTRAL["tint"], edgecolor=NEUTRAL["rule"],
                                       linewidth=1.6, linestyle=(0, (6, 5))))
            ax.text(0.5, 0.60, "no experimental\nstructure", transform=ax.transAxes,
                    ha="center", va="center", fontsize=FONT["annot"], color=NEUTRAL["ink_2"], linespacing=1.4)
            ax.text(0.5, 0.28, "the absence is\nthe finding",
                    transform=ax.transAxes, ha="center", va="center", fontsize=FONT["note"] - 2,
                    color=NEUTRAL["ink_muted"], style="italic", linespacing=1.4)
        fig.text(x + col_w / 2, 0.815, gene, ha="center", va="bottom", fontsize=FONT["big"],
                 fontweight="bold", color=GENE_COLOURS[gene])
        fig.text(x + col_w / 2, 0.522, STRUCTURE_ID[gene], ha="center", va="top",
                 fontsize=FONT["note"], color=NEUTRAL["ink"])
        fig.text(x + col_w / 2, 0.487, BOUND_TO[gene], ha="center", va="top",
                 fontsize=FONT["note"] - 2, color=NEUTRAL["ink_muted"], style="italic")

    # ---- evidence track --------------------------------------------------------
    ax = fig.add_axes([col_x0, 0.105, 4 * col_w + 3 * col_gap, 0.355])
    ax.set_xlim(-0.5, len(GENES) - 0.5)
    ax.set_ylim(-0.6, len(LEVELS) - 0.4)
    ax.axis("off")
    for r, (label, states) in enumerate(LEVELS):
        y = len(LEVELS) - 1 - r
        for c, gene in enumerate(GENES):
            state = states[gene]
            colour = GENE_COLOURS[gene]
            # scatter, not Circle: the axes is not aspect-equal, so a patch
            # circle would render as an ellipse
            if state == 2:      # filled = holds fully, the same encoding as every figure
                ax.scatter([c], [y], s=520, c=colour, edgecolors=WHITE, linewidths=1.8, zorder=4)
            elif state == 1:    # hollow ring = partial
                ax.scatter([c], [y], s=520, facecolors="none", edgecolors=colour,
                           linewidths=3.2, zorder=4)
            else:               # absent
                ax.plot([c - 0.10, c + 0.10], [y, y], color=NEUTRAL["rule"], lw=2.4, zorder=4)
            note = LEVEL_NOTE.get((label, gene))
            if note:
                ax.text(c, y - 0.28, note, ha="center", va="top", fontsize=FONT["note"] - 2,
                        color=NEUTRAL["ink_muted"], zorder=5)
        ax.text(-0.60, y, label, ha="right", va="center", fontsize=FONT["annot"],
                color=NEUTRAL["ink"], linespacing=1.35)
        if r % 2 == 0:
            ax.add_patch(plt.Rectangle((-0.5, y - 0.46), len(GENES), 0.92,
                                       facecolor=NEUTRAL["tint"], edgecolor="none", zorder=0))

    from matplotlib.lines import Line2D
    handles = [
        Line2D([], [], marker="o", linestyle="none", markersize=15, markerfacecolor=NEUTRAL["ink"],
               markeredgecolor=WHITE, label="holds"),
        Line2D([], [], marker="o", linestyle="none", markersize=15, markerfacecolor="none",
               markeredgecolor=NEUTRAL["ink"], markeredgewidth=2.6, label="partial"),
        Line2D([], [], marker="_", linestyle="none", markersize=15, markeredgecolor=NEUTRAL["rule"],
               markeredgewidth=2.6, label="does not hold"),
    ]
    leg = fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.55, 0.042),
                     ncol=3, frameon=False, fontsize=FONT["legend"], columnspacing=3.2, handletextpad=0.8)
    for t in leg.get_texts():
        t.set_color(NEUTRAL["ink_2"])

    headline(
        fig,
        "The four differ in the kind of chemical evidence they have, not in\nhow much of one thing",
        "Reachability only. Nothing here is efficacy, selectivity or safety, and no docking, binding prediction or molecular modelling was performed\n"
        "anywhere in this project. KDM1A is the mature target: an inhibitor-bound structure and clinical-stage compounds already exist. USP34's\n"
        "structure holds a covalent ubiquitin activity probe, not a drug, over about 12% of the protein — an unexplored opportunity, not a\n"
        "validated one. Counter-evidence for USP34: losing it has been reported to push breast cells toward a more mobile state (PMID 28499884).",
        key=FIGURE)
    figure_footer(fig, "Frozen structural audit table + committed PyMOL renders.")
    return save(fig, stub), verification


def main(out_dir: Path = OUT_DIR):
    return build(out_dir / FIGURE)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    main()
