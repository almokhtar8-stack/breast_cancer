"""ONE poster-grade structure/druggability figure -- answers "can these four
tamoxifen-sensitising candidate vulnerabilities realistically be targeted?"

The four candidates are NOT equally tractable, and that difference is the
result. Every claim comes from the project's already-audited structural /
pharmacology evidence table,
`results/tables/post_audit_sensitivity/06b_structural_tractability_audit.tsv`
(columns A-G plus `precise_summary`/`sources`), read unmodified. Nothing is
re-searched, re-scored, or re-ranked here.

Structural panels are static PyMOL renders of the ALREADY-DOWNLOADED,
provenance-recorded experimental PDB files (see
`scripts/render_druggability_structures.py`): KDM1A 6NQU (inhibitor-bound),
TLK2 5O0Y (ATP-analog-bound), USP34 7W3U (covalent activity-based-probe-
bound). VEZF1 has no experimental structure in the audited evidence and is
deliberately shown WITHOUT one -- no homology or AlphaFold model is
substituted, because that absence is itself the finding.

NO docking, pose prediction, affinity estimation or pocket re-detection is
performed anywhere in this module or its render script.

Evidence-display choice: a compact 3-row evidence track (experimental
structure / direct ligand-or-probe evidence / selective-or-clinical
pharmacology) rather than a single "maturity score". A single axis would
force TLK2 and USP34 into a false ordering -- they differ in KIND, not
degree (TLK2: strong structural class, no selective chemistry; USP34:
direct covalent catalytic-cysteine probe evidence, no selective chemistry).
No numerical score is invented, because the project has no validated
scoring framework for this.
"""

from __future__ import annotations

import logging
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import pandas as pd

logger = logging.getLogger(__name__)

OUT_DIR = Path("results/figures/poster_druggability_v1")
RENDER_DIR = OUT_DIR / "renders"
AUDIT_TSV = Path("results/tables/post_audit_sensitivity/06b_structural_tractability_audit.tsv")

FOCUS_FOUR = ["KDM1A", "TLK2", "USP34", "VEZF1"]
FOCUS_COLORS = {"KDM1A": "#D55E00", "TLK2": "#56B4E9", "USP34": "#0072B2", "VEZF1": "#E69F00"}

DGRAY = "#262626"
MGRAY = "#8c8c8c"
LGRAY = "#c9c9c9"

# Per-candidate poster-facing presentation. Every factual line here is a
# short restatement of the audited table's own wording (verified against it
# by `check_against_audit()` and by the test suite); the PDB IDs are read
# FROM the audit text, not hand-asserted.
PANELS = {
    "KDM1A": dict(
        target_class="Lysine demethylase",
        render="KDM1A_6NQU.png",
        structure_line="Experimental inhibitor-bound structure",
        ligand_line="GSK2879552 (selective inhibitor)",
        pharm_line="Clinical-stage selective inhibitors",
        detail="e.g. iadademstat / ORY-1001",
    ),
    "TLK2": dict(
        target_class="Ser/Thr kinase",
        render="TLK2_5O0Y.png",
        structure_line="Experimental kinase-domain structure",
        ligand_line="ATP-γ-S substrate analog — not an inhibitor",
        pharm_line="No validated selective clinical inhibitor",
        detail="ATP-site structure; TLK1-paralog selectivity unsolved",
    ),
    "USP34": dict(
        target_class="Deubiquitinase",
        render="USP34_7W3U.png",
        structure_line="Experimental catalytic-domain structure",
        ligand_line="Covalent activity-based probe — not a drug",
        pharm_line="No validated selective inhibitor",
        detail="reactive catalytic Cys1903",
    ),
    "VEZF1": dict(
        target_class="Zinc-finger transcription factor",
        render=None,
        structure_line="No experimental structure identified\nin the audited analysis",
        ligand_line="Homology model only",
        pharm_line="No validated direct pharmacology",
        detail="weak published screening hit (IC50 ≈ 20 µM)",
    ),
}

# The 3-row evidence track. Levels are read from the audited table's own
# boolean/text columns -- see `evidence_track()`.
TRACK_ROWS = [
    ("Experimental structure", "structure"),
    ("Direct ligand / probe evidence", "ligand"),
    ("Selective or clinical pharmacology", "pharmacology"),
]


def _load_trimmed_render(path: Path):
    """Loads a PyMOL render and crops it to its non-transparent content, so
    the structure fills its panel instead of floating inside the render's
    transparent margin."""
    import numpy as np

    img = mpimg.imread(path)
    if img.shape[2] == 4:
        opaque = img[:, :, 3] > 0.02
        rows = np.where(opaque.any(axis=1))[0]
        cols = np.where(opaque.any(axis=0))[0]
        if len(rows) and len(cols):
            pad = 6
            r0, r1 = max(rows[0] - pad, 0), min(rows[-1] + pad + 1, img.shape[0])
            c0, c1 = max(cols[0] - pad, 0), min(cols[-1] + pad + 1, img.shape[1])
            img = img[r0:r1, c0:c1]
    return img


def load_audit() -> pd.DataFrame:
    """The frozen audited structural-tractability table, read unmodified."""
    audit = pd.read_csv(AUDIT_TSV, sep="\t").set_index("gene")
    missing = set(FOCUS_FOUR) - set(audit.index)
    assert not missing, f"audited table missing focus genes: {missing}"
    logger.info("audit table: %d genes in -> %d focus genes used", len(audit), len(FOCUS_FOUR))
    return audit.loc[FOCUS_FOUR]


def pdb_ids_from_audit(audit: pd.DataFrame) -> dict[str, list[str]]:
    """PDB IDs recovered by scanning the audited table's own text, so no
    structure ID in this figure is hand-asserted."""
    import re

    text_cols = ["B_relevant_domain_structure", "C_ligand_or_probe_bound_structure", "sources",
                 "precise_summary"]
    out: dict[str, list[str]] = {}
    for gene in FOCUS_FOUR:
        blob = " ".join(str(audit.loc[gene, c]) for c in text_cols)
        ids = re.findall(r"\b(?:PDB\s+)?(\d[A-Z0-9]{3})\b", blob)
        out[gene] = sorted({i for i in ids if not i.isdigit()})
    return out


def evidence_track(audit: pd.DataFrame) -> pd.DataFrame:
    """Three ordinal evidence levels per gene, derived from the audited
    columns: "yes" / "limited" / "no". "limited" is used where the audit
    itself says PARTIAL, or where the bound species is explicitly not a
    drug (an ATP analog / an activity-based probe)."""
    rows = []
    for gene in FOCUS_FOUR:
        a = audit.loc[gene]
        structure = "yes" if bool(a["A_experimental_human_structure_exists"]) else "no"
        if structure == "yes" and str(a["B_relevant_domain_structure"]).strip().startswith("PARTIAL"):
            structure = "limited"

        c_text = str(a["C_ligand_or_probe_bound_structure"])
        if c_text.strip().lower().startswith(("false", "no")):
            ligand = "no"
        elif c_text.strip().startswith("PARTIAL") or "not a drug" in c_text or "NOT a small-molecule inhibitor" in c_text:
            ligand = "limited"
        else:
            ligand = "yes"

        e_yes = str(a["E_validated_selective_small_molecule_inhibitor"]).strip().upper().startswith("YES")
        f_yes = str(a["F_clinical_stage_pharmacology"]).strip().upper().startswith("YES")
        pharmacology = "yes" if (e_yes and f_yes) else ("limited" if (e_yes or f_yes) else "no")

        rows.append({"gene": gene, "structure": structure, "ligand": ligand, "pharmacology": pharmacology})
    track = pd.DataFrame(rows).set_index("gene")
    logger.info("evidence track derived:\n%s", track.to_string())
    return track


def check_against_audit(audit: pd.DataFrame) -> None:
    """Guards the poster wording against the audited source: the figure must
    not claim a validated/clinical inhibitor where the audit says none, and
    must not call the TLK2 ATP analog or the USP34 probe a drug."""
    tlk2_c = str(audit.loc["TLK2", "C_ligand_or_probe_bound_structure"])
    assert "NOT a small-molecule inhibitor" in tlk2_c or "not an inhibitor" in tlk2_c.lower()
    usp34_c = str(audit.loc["USP34", "C_ligand_or_probe_bound_structure"])
    assert "ACTIVITY-BASED PROBE" in usp34_c or "activity-based probe" in usp34_c.lower()
    assert not str(audit.loc["VEZF1", "A_experimental_human_structure_exists"]).strip().lower() == "true"
    assert str(audit.loc["KDM1A", "F_clinical_stage_pharmacology"]).strip().upper().startswith("YES")
    for gene in ("TLK2", "USP34", "VEZF1"):
        assert not str(audit.loc[gene, "F_clinical_stage_pharmacology"]).strip().upper().startswith("YES")


def _draw_marker(ax, x: float, y: float, level: str, color: str) -> None:
    """filled = yes, half = limited, open = no."""
    if level == "yes":
        ax.scatter([x], [y], s=190, color=color, edgecolor="white", linewidth=1.1, zorder=4)
    elif level == "limited":
        ax.scatter([x], [y], s=190, facecolor="white", edgecolor=color, linewidth=2.2, zorder=4)
        ax.scatter([x], [y], s=190, color=color, edgecolor="none", zorder=5,
                    marker=matplotlib.markers.MarkerStyle("o", fillstyle="left"))
    else:
        ax.scatter([x], [y], s=190, facecolor="white", edgecolor=LGRAY, linewidth=1.8, zorder=4)


def build_druggability_v1(stub: Path) -> None:
    audit = load_audit()
    check_against_audit(audit)
    track = evidence_track(audit)
    pdbs = pdb_ids_from_audit(audit)
    logger.info("PDB IDs recovered from audit text: %s", pdbs)

    fig = plt.figure(figsize=(17.5, 9.8), dpi=300)
    gs = fig.add_gridspec(2, 4, height_ratios=[0.58, 0.42], hspace=0.02, wspace=0.03,
                          left=0.152, right=0.985, top=0.782, bottom=0.095)
    structure_axes: list = []
    evidence_axes: list = []

    for j, gene in enumerate(FOCUS_FOUR):
        spec = PANELS[gene]
        color = FOCUS_COLORS[gene]

        # ---- upper: structural visualization --------------------------
        ax = fig.add_subplot(gs[0, j])
        ax.set_xticks([]); ax.set_yticks([])
        for s in ax.spines.values():
            s.set_visible(False)

        structure_axes.append(ax)
        if spec["render"]:
            img = _load_trimmed_render(RENDER_DIR / spec["render"])
            ax.imshow(img)
            ax.set_xlim(0, img.shape[1]); ax.set_ylim(img.shape[0], 0)
            pdb_shown = spec["render"].split("_")[1].replace(".png", "")
            assert pdb_shown in pdbs[gene], f"{pdb_shown} not recoverable from audit text for {gene}"
        else:
            # Deliberate structural absence: restrained zinc-finger domain
            # schematic, visually unlike the experimental-structure panels.
            ax.set_xlim(0, 1); ax.set_ylim(0, 1)
            ax.add_patch(plt.Rectangle((0.06, 0.40), 0.88, 0.14, facecolor="#F3F3F1",
                                        edgecolor=LGRAY, linewidth=1.0, zorder=2))
            for k in range(6):
                cx = 0.145 + k * 0.142
                ax.add_patch(plt.Circle((cx, 0.47), 0.036, facecolor=color, alpha=0.55,
                                         edgecolor="white", linewidth=1.0, zorder=3))
            ax.text(0.5, 0.335, "six tandem zinc fingers (domain schematic)", fontsize=9.5,
                     color=MGRAY, ha="center", va="top")
            ax.text(0.5, 0.68, "No experimental structure", fontsize=13, color=MGRAY,
                     ha="center", va="center", style="italic")


        # ---- lower: compact evidence track ----------------------------
        axe = fig.add_subplot(gs[1, j])
        evidence_axes.append(axe)
        axe.set_xlim(0, 1); axe.set_ylim(0, 1)
        axe.set_xticks([]); axe.set_yticks([])
        for side in ("top", "right", "bottom"):
            axe.spines[side].set_visible(False)
        axe.spines["left"].set_visible(False)

        axe.text(0.5, 0.965, spec["structure_line"], fontsize=10.8, color=DGRAY,
                  ha="center", va="top", linespacing=1.3)
        axe.text(0.5, 0.735, spec["ligand_line"], fontsize=10.2, color=MGRAY,
                  ha="center", va="top")

        for i, (_, key) in enumerate(TRACK_ROWS):
            y = 0.47 - i * 0.155
            _draw_marker(axe, 0.5, y, track.loc[gene, key], color)

        axe.text(0.5, 0.055, spec["detail"], fontsize=9.3, color=MGRAY, ha="center",
                  va="bottom", style="italic")

    # Gene headings, drawn in FIGURE coords from each structure column's own
    # centre so all four align regardless of whether that panel holds an
    # image (imshow shrinks its axes to the image aspect; the VEZF1
    # schematic panel does not).
    fig.canvas.draw()
    for gene, ax in zip(FOCUS_FOUR, structure_axes):
        box = ax.get_position()
        cx = (box.x0 + box.x1) / 2
        fig.text(cx, 0.855, gene, fontsize=25, fontweight="bold",
                  color=FOCUS_COLORS[gene], ha="center", va="bottom")
        fig.text(cx, 0.828, PANELS[gene]["target_class"], fontsize=12,
                  color=DGRAY, ha="center", va="bottom")
        spec = PANELS[gene]
        pdb_tag = (f"PDB {spec['render'].split('_')[1].replace('.png','')}"
                   if spec["render"] else "no PDB entry")
        fig.text(cx, 0.803, pdb_tag, fontsize=10.5, color=MGRAY, ha="center", va="bottom",
                  family="monospace")

    # Shared evidence-track row labels, positioned from the real marker
    # coordinates of the first evidence panel rather than guessed offsets.
    first = evidence_axes[0]
    label_x = first.get_position().x0 - 0.012
    for i, (label, _) in enumerate(TRACK_ROWS):
        y_axes = 0.47 - i * 0.155
        y_fig = first.transAxes.transform((0, y_axes))[1]
        y_fig = fig.transFigure.inverted().transform((0, y_fig))[1]
        fig.text(label_x, y_fig, label, fontsize=10.5, color=DGRAY, ha="right", va="center")

    fig.text(0.035, 0.978, "Structural and pharmacological evidence differentiates candidate tractability",
              fontsize=21, fontweight="bold", color=DGRAY, ha="left", va="top")
    fig.text(0.035, 0.940,
              "Experimental structural evidence and audited pharmacological maturity across the four focus genes.",
              fontsize=12, color="#555555", ha="left", va="top")

    legend_y = 0.052
    fig.text(0.618, legend_y, "Evidence:", fontsize=10, color=DGRAY, ha="left", va="center")
    for k, (lab, lvl) in enumerate([("established", "yes"), ("limited / not a drug", "limited"), ("none", "no")]):
        x = 0.684 + k * 0.104
        ax_tmp = fig.add_axes([x, legend_y - 0.010, 0.012, 0.020])
        ax_tmp.axis("off"); ax_tmp.set_xlim(0, 1); ax_tmp.set_ylim(0, 1)
        _draw_marker(ax_tmp, 0.5, 0.5, lvl, "#6b6b6b")
        fig.text(x + 0.016, legend_y, lab, fontsize=9.6, color=MGRAY, ha="left", va="center")

    fig.text(0.035, 0.014,
              "Structural evidence indicates tractability, not efficacy. No docking or binding prediction was performed; "
              "a bound ATP analog or activity-based probe is not a therapeutic inhibitor.",
              fontsize=9.4, color=MGRAY, ha="left", va="center", style="italic")

    stub.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(stub.with_suffix(".png"), facecolor="white", bbox_inches="tight", dpi=300)
    fig.savefig(stub.with_suffix(".pdf"), facecolor="white", bbox_inches="tight")
    fig.savefig(stub.with_suffix(".svg"), facecolor="white", bbox_inches="tight")
    plt.close(fig)
    logger.info("wrote %s.png/.pdf/.svg", stub)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    build_druggability_v1(OUT_DIR / "DRUGGABILITY_v1")
