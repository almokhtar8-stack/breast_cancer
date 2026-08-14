"""Druggability + normal-tissue/selectivity review: builds the five output
tables for this phase from the curated dossier data in
src/druggability_safety_data.py, joined against real upstream project
output files (never hand-typed) for the therapeutic-window summary.

USP34, VEZF1, EML5, CITED2 -- frozen therapeutic shortlist and ranking,
unchanged. No CRISPR/RNA-seq/TCGA/DepMap computation happens here; DepMap
and Hany-CRISPR values are read directly from the already-frozen/verified
independent-validation and evidence-freeze output files, never retyped.

No network calls happen here (module is deterministic; all curation was
gathered separately via scripts/research passes and is inlined in
druggability_safety_data.py, per this project's data-hygiene rule that
downloads/lookups live outside src/).
"""

from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd

from src.druggability_safety_data import (
    BONE_MUSCULOSKELETAL_ROWS,
    CANDIDATES,
    DRUGGABILITY_ROWS,
    GENETIC_CONSTRAINT_ROWS,
    NORMAL_TISSUE_ROWS,
    REFERENCES_ROWS,
)

logger = logging.getLogger(__name__)

OUT_DIR = Path("results/tables/druggability_safety")

FROZEN_HANY = Path("results/tables/evidence_freeze/THERAPEUTIC_SHORTLIST_FREEZE.tsv")
DEPMAP_DEPENDENCY = Path("results/tables/independent_validation/DepMap_candidate_dependency.tsv")
INTEGRATION_TABLE = Path("results/tables/independent_validation/four_candidate_independent_validation.tsv")
LITERATURE_COMPARISON = Path("results/tables/literature_mechanism/four_candidate_literature_comparison.tsv")

# Qualitative concern TIER only (a curated judgment call, same status as
# druggability_classification) -- the supporting facts quoted alongside
# each tier are pulled live from the normal_tissue/genetic_constraint/bone
# DataFrames inside build_therapeutic_window_summary below, never
# duplicated as a separate hardcoded narrative, so they cannot drift out
# of sync with those tables if a source row is later corrected. Never
# derived from DepMap (a cancer-cell-line resource, out of scope for
# normal-tissue inference per explicit project rule).
NORMAL_TISSUE_CONCERN_TIER = {"USP34": "MODERATE-HIGH", "VEZF1": "MODERATE", "EML5": "LOWER", "CITED2": "MODERATE-HIGH"}


def _normal_tissue_concern_text(candidate: str, tissue_idx: pd.DataFrame, constraint_idx: pd.DataFrame) -> str:
    tier = NORMAL_TISSUE_CONCERN_TIER[candidate]
    spec = tissue_idx.loc[candidate, "hpa_tissue_specificity_category"]
    loeuf = constraint_idx.loc[candidate, "loeuf"]
    pli = constraint_idx.loc[candidate, "pli"]
    mouse = constraint_idx.loc[candidate, "mouse_ko_phenotype_summary"]
    return f"{tier} -- HPA: {spec}. gnomAD LOEUF={loeuf}, pLI={pli}. Mouse KO: {mouse}"


def _bone_musculoskeletal_concern_text(candidate: str, bone_idx: pd.DataFrame) -> str:
    cat = bone_idx.loc[candidate, "bone_concern_category"]
    species = bone_idx.loc[candidate, "evidence_species"]
    role = bone_idx.loc[candidate, "published_bone_role_summary"]
    role_short = role if len(role) <= 260 else role[:260] + "..."
    return f"{cat} (species: {species}) -- {role_short}"


def build_druggability_table() -> pd.DataFrame:
    df = pd.DataFrame(DRUGGABILITY_ROWS)
    assert list(df["candidate"]) == CANDIDATES, "druggability rows must cover all four frozen candidates in order"
    logger.info("build_druggability_table: %d rows", len(df))
    return df


def build_normal_tissue_table() -> pd.DataFrame:
    df = pd.DataFrame(NORMAL_TISSUE_ROWS)
    assert list(df["candidate"]) == CANDIDATES, "normal-tissue rows must cover all four frozen candidates in order"
    logger.info("build_normal_tissue_table: %d rows", len(df))
    return df


def build_genetic_constraint_table() -> pd.DataFrame:
    df = pd.DataFrame(GENETIC_CONSTRAINT_ROWS)
    assert list(df["candidate"]) == CANDIDATES, "genetic-constraint rows must cover all four frozen candidates in order"
    logger.info("build_genetic_constraint_table: %d rows", len(df))
    return df


def build_bone_musculoskeletal_table() -> pd.DataFrame:
    df = pd.DataFrame(BONE_MUSCULOSKELETAL_ROWS)
    assert list(df["candidate"]) == CANDIDATES, "bone rows must cover all four frozen candidates in order"
    logger.info("build_bone_musculoskeletal_table: %d rows", len(df))
    return df


def build_verified_references_table() -> pd.DataFrame:
    df = pd.DataFrame(REFERENCES_ROWS)
    logger.info("build_verified_references_table: %d rows across %d candidates", len(df), df["candidate"].nunique())
    return df


def build_therapeutic_window_summary(
    druggability: pd.DataFrame,
    normal_tissue: pd.DataFrame,
    genetic_constraint: pd.DataFrame,
    bone: pd.DataFrame,
) -> pd.DataFrame:
    """Separate evidence dimensions, joined per candidate -- never collapsed
    into one master score (per explicit project instruction). Tamoxifen/
    DepMap/mechanism/human-tumor columns are read directly from the
    frozen/independent-validation/literature-mechanism output files;
    druggability/normal-tissue/bone/genetic-constraint columns come from
    this phase's own curated tables above.
    """
    hany = pd.read_csv(FROZEN_HANY, sep="\t").set_index("gene")
    depmap = pd.read_csv(DEPMAP_DEPENDENCY, sep="\t").set_index("candidate")
    integration = pd.read_csv(INTEGRATION_TABLE, sep="\t").set_index("candidate")
    lit = pd.read_csv(LITERATURE_COMPARISON, sep="\t").set_index("candidate")

    drug_idx = druggability.set_index("candidate")
    tissue_idx = normal_tissue.set_index("candidate")
    constraint_idx = genetic_constraint.set_index("candidate")
    bone_idx = bone.set_index("candidate")

    rows = []
    for c in CANDIDATES:
        h = hany.loc[c]
        d = depmap.loc[c]
        integ = integration.loc[c]
        rows.append(dict(
            candidate=c,
            functional_tamoxifen_evidence=f"{h['crispr_direction']}, effect={h['crispr_effect']:.3f}, FDR={h['crispr_fdr']:.4f} (frozen Hany screen)",
            cancer_dependency_depmap=f"{d['essentiality_concern']} (ER+/luminal strongly-dependent fraction={d['frac_strongly_dependent_er_luminal']*100:.1f}%, DepMap {d['depmap_release']})" if pd.notna(d["frac_strongly_dependent_er_luminal"]) else f"{d['essentiality_concern']}",
            mechanism_strength=lit.loc[c, "overall_literature_depth"],
            human_tumor_support=f"{integ['tcga_interpretation']} (validation_strength={integ['integration_validation_strength']})",
            direct_druggability=drug_idx.loc[c, "druggability_classification"],
            alternative_modality_feasibility=drug_idx.loc[c, "plausible_alternative_modalities"],
            normal_tissue_concern=_normal_tissue_concern_text(c, tissue_idx, constraint_idx),
            bone_musculoskeletal_concern=_bone_musculoskeletal_concern_text(c, bone_idx),
            genetic_constraint_concern=f"LOEUF={constraint_idx.loc[c, 'loeuf']}, pLI={constraint_idx.loc[c, 'pli']}",
            novelty=f"mechanism_confidence={lit.loc[c, 'mechanism_confidence']}; no existing chemical probe/inhibitor/degrader verified for this candidate" if drug_idx.loc[c, "known_ligands_tools_probes_degraders"].startswith("NOT FOUND") or "NOT FOUND" in drug_idx.loc[c, "known_ligands_tools_probes_degraders"] or "NONE" in drug_idx.loc[c, "known_ligands_tools_probes_degraders"][:20] else f"mechanism_confidence={lit.loc[c, 'mechanism_confidence']}; one weak tool-compound series exists",
            mechanistic_action_category=integ["mechanistic_action_category"],
            notes=f"Raw targetability rank/score are FROZEN and unchanged by this phase (see results/tables/independent_validation/four_candidate_followup_rankings.tsv); this row adds druggability/normal-tissue/bone dimensions alongside, not a replacement ranking.",
        ))
    out = pd.DataFrame(rows)
    logger.info(
        "build_therapeutic_window_summary: %d candidates joined from %s, %s, %s, %s (frozen/upstream) + this phase's own curated tables",
        len(out), FROZEN_HANY, DEPMAP_DEPENDENCY, INTEGRATION_TABLE, LITERATURE_COMPARISON,
    )
    return out


def run(out_dir: Path = OUT_DIR) -> dict[str, pd.DataFrame]:
    out_dir.mkdir(parents=True, exist_ok=True)

    druggability = build_druggability_table()
    normal_tissue = build_normal_tissue_table()
    genetic_constraint = build_genetic_constraint_table()
    bone = build_bone_musculoskeletal_table()
    references = build_verified_references_table()
    window = build_therapeutic_window_summary(druggability, normal_tissue, genetic_constraint, bone)

    tables = {
        "candidate_druggability.tsv": druggability,
        "candidate_normal_tissue_context.tsv": normal_tissue,
        "candidate_genetic_constraint.tsv": genetic_constraint,
        "candidate_bone_musculoskeletal_context.tsv": bone,
        "candidate_therapeutic_window_summary.tsv": window,
        "verified_references.tsv": references,
    }
    for name, df in tables.items():
        path = out_dir / name
        df.to_csv(path, sep="\t", index=False)
        logger.info("wrote %s (%d rows)", path, len(df))
    return tables


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run()
