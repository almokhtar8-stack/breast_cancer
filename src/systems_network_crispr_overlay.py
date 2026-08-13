"""Systems-network phase 10: CRISPR functional overlay onto resistance
pathways, recurrent leading-edge genes, and candidate-related pathways.

Gene population overlaid = union of: (a) every gene that appears in the
GSEA leading edge of a STRONG_CONSENSUS resistance-consensus pathway in any
of the three resistance datasets (Phase 5/9's population, not restricted to
the >=2-dataset "recurrent" subset here -- the full single-dataset set too),
and (b) every gene associated with a frozen candidate in
candidate_pathway_membership.tsv (Phase 7: curated member, leading-edge
gene, or a candidate's direct STRING interactor).

CRISPR classification follows CLAUDE.md/dataset semantics exactly:
negative effect + FDR<0.05 = sensitising_KO; positive effect + FDR<0.05 =
tolerance_associated_KO; FDR>=0.05 is never called functional evidence,
only a nonsignificant directional lean.
"""

from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd
import yaml

logger = logging.getLogger(__name__)

RESISTANCE_DATASETS = ["gse118713", "gse240112", "gse111151"]


def _load_config(config_path: str | Path) -> dict:
    with open(config_path) as f:
        return yaml.safe_load(f)


def classify_crispr(effect: float, fdr: float) -> str:
    if pd.isna(fdr):
        return "not_tested"
    if fdr < 0.05:
        return "sensitising_KO" if effect < 0 else "tolerance_associated_KO"
    return "nonsignificant_sensitising_direction" if effect < 0 else "nonsignificant_tolerance_direction"


def build_resistance_leading_edge_gene_pathways(consensus: pd.DataFrame, gsea_tables: dict[str, pd.DataFrame]) -> pd.DataFrame:
    strong_keys = set(zip(consensus.loc[consensus["consensus_category"] == "STRONG_CONSENSUS", "collection"], consensus.loc[consensus["consensus_category"] == "STRONG_CONSENSUS", "pathway"]))
    records = []
    for dataset in RESISTANCE_DATASETS:
        df = gsea_tables[dataset]
        df = df.loc[df.apply(lambda r: (r["collection"], r["pathway"]) in strong_keys, axis=1)]
        # Codex review fix (same as Phase 9): only count a gene as
        # leading-edge in THIS dataset if this dataset's own nominal
        # p-value for the pathway is significant, not merely because the
        # pathway is STRONG_CONSENSUS overall.
        df = df.loc[df["nom_pvalue"] < 0.05]
        for _, row in df.iterrows():
            if pd.isna(row["leading_edge_genes"]):
                continue
            for gene in row["leading_edge_genes"].split(";"):
                records.append({"gene": gene, "collection": row["collection"], "pathway": row["pathway"]})
    return pd.DataFrame(records)


def build_pathway_crispr_overlay(
    resistance_le: pd.DataFrame, candidate_membership: pd.DataFrame, crispr_ranked: pd.DataFrame
) -> pd.DataFrame:
    resistance_genes = set(resistance_le["gene"])
    candidate_related_genes: set[str] = set(candidate_membership["candidate"].unique())
    for col in ("interactor_genes_in_pathway",):
        for val in candidate_membership[col].dropna():
            candidate_related_genes.update(g for g in val.split(";") if g)
    le_mask = candidate_membership["candidate_is_leading_edge_datasets"].fillna("") != ""
    candidate_related_genes.update(candidate_membership.loc[candidate_membership["candidate_is_member"] | le_mask, "candidate"])

    all_genes = sorted(resistance_genes | candidate_related_genes)

    rows = []
    crispr_idx = crispr_ranked.set_index("gene")
    for gene in all_genes:
        in_resistance_le = gene in resistance_genes
        in_candidate_related = gene in candidate_related_genes
        source = "both" if (in_resistance_le and in_candidate_related) else ("resistance_consensus_leading_edge" if in_resistance_le else "candidate_related")
        n_pathways = resistance_le.loc[resistance_le["gene"] == gene, ["collection", "pathway"]].drop_duplicates().shape[0] if in_resistance_le else 0
        example = ";".join(sorted((resistance_le.loc[resistance_le["gene"] == gene, "collection"] + ":" + resistance_le.loc[resistance_le["gene"] == gene, "pathway"]).unique())[:5]) if in_resistance_le else ""

        if gene in crispr_idx.index:
            effect = float(crispr_idx.loc[gene, "log2fc"])
            p = float(crispr_idx.loc[gene, "p_value"])
            fdr = float(crispr_idx.loc[gene, "fdr"])
        else:
            effect, p, fdr = float("nan"), float("nan"), float("nan")

        rows.append(
            {
                "gene": gene,
                "source": source,
                "n_strong_consensus_pathways": n_pathways,
                "example_pathways": example,
                "crispr_effect": effect,
                "crispr_p": p,
                "crispr_fdr": fdr,
                "crispr_direction": classify_crispr(effect, fdr),
            }
        )
    out = pd.DataFrame(rows).sort_values(["crispr_direction", "crispr_fdr"], na_position="last")
    return out


def run_crispr_overlay(config_path: str | Path = "config/config.yaml") -> pd.DataFrame:
    config = _load_config(config_path)
    cfg = config["systems_network"]
    tables_dir = Path(cfg["output"]["tables_dir"])

    consensus = pd.read_csv(tables_dir / "resistance_pathway_consensus.tsv", sep="\t")
    gsea_tables = {d: pd.read_csv(tables_dir / f"gsea_{d}.tsv", sep="\t") for d in RESISTANCE_DATASETS}
    candidate_membership = pd.read_csv(tables_dir / "candidate_pathway_membership.tsv", sep="\t")
    crispr_ranked = pd.read_csv(tables_dir / "crispr_ranked_genes.tsv", sep="\t")

    resistance_le = build_resistance_leading_edge_gene_pathways(consensus, gsea_tables)
    out = build_pathway_crispr_overlay(resistance_le, candidate_membership, crispr_ranked)
    out.to_csv(tables_dir / "pathway_crispr_overlay.tsv", sep="\t", index=False)

    counts = out["crispr_direction"].value_counts().to_dict()
    logger.info("build_pathway_crispr_overlay: %d genes overlaid, direction counts=%s", len(out), counts)
    return out


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run_crispr_overlay()
