"""Final USP34/VEZF1 GDSC pharmacogenomics phase: builds the eight output
tables from real, deterministic computation over the local GDSC Release
8.5 files and this project's already-verified DepMap Public 26Q1
expression/model data. No candidate discovery; USP34 and VEZF1 only, per
the frozen lead/backup translational conclusion. No frozen upstream
evidence (Hany, TCGA, DepMap dependency, druggability_safety,
lead_target_deep_dive, final_translational) is altered -- this module
only reads a few of those tables read-only for the project cross-check.
"""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

from src.final_pharmacogenomics_interpretation_data import (
    DRUG_INDIRECT_TARGETING_ROWS,
    FINAL_INTERPRETATION_ROWS,
)
from src.final_pharmacogenomics_gdsc_data import (
    GENES,
    MIN_N_ER_LUMINAL_EXPLORATORY,
    MIN_N_FULL_BREAST,
    build_breast_expression_joined,
    load_config,
    load_gdsc_compounds,
)

logger = logging.getLogger(__name__)

OUT_DIR = Path("results/tables/final_pharmacogenomics")
CROSS_DATASET_TABLE = Path("results/tables/cross_dataset_genomewide/all_genes_cross_dataset_evidence_with_ranking.tsv")

# Conservative, curated (not automated) pathway-coherence keyword groups
# for Part 12's "recurrent drug classes/pathways" check -- only used to
# annotate results already computed from real data, never to select which
# drugs are tested.
PATHWAY_COHERENCE_GROUPS = {
    "estrogen_signaling_endocrine": ["Hormone-related"],
    "wnt_beta_catenin": ["WNT signaling"],
    "cell_cycle": ["Cell cycle"],
    "apoptosis_regulation": ["Apoptosis regulation"],
    "genome_integrity_dna_damage": ["Genome integrity"],
    "rtk_signaling": ["RTK signaling"],
    "pi3k_mtor": ["PI3K/MTOR signaling"],
    "chromatin_histone": ["Chromatin histone acetylation", "Chromatin other"],
    "mitosis": ["Mitosis"],
}


def build_data_provenance_table() -> pd.DataFrame:
    cfg = load_config()
    g = cfg["final_pharmacogenomics"]["gdsc"]
    rows = [
        dict(field="release", value=g["release"]),
        dict(field="release_date", value=g["release_date"]),
        dict(field="source_url", value=g["source_url"]),
        dict(field="gdsc1_file", value=g["raw"]["gdsc1_fitted_dose_response"]),
        dict(field="gdsc2_file", value=g["raw"]["gdsc2_fitted_dose_response"]),
        dict(field="compounds_file", value=g["raw"]["screened_compounds"]),
        dict(field="cell_lines_file", value=g["raw"]["cell_lines_details"]),
        dict(field="LN_IC50_definition", value=g["response_metric_definitions"]["LN_IC50"]),
        dict(field="AUC_definition", value=g["response_metric_definitions"]["AUC"]),
        dict(field="direction_source", value=g["response_metric_definitions"]["source_of_definitions"]),
        dict(field="expression_source", value=g["expression_source"]),
        dict(field="cell_line_join_key", value="SangerModelID / COSMIC_ID (exact ID match, both files carry it natively)"),
        dict(field="raw_data_location", value="OUTSIDE git -- see PROVENANCE.txt in the raw data directory (SHA256 recorded)"),
    ]
    return pd.DataFrame(rows)


def build_compound_availability_table() -> pd.DataFrame:
    cfg = load_config()
    comp = load_gdsc_compounds(cfg)
    fitted = None  # computed lazily below only for the found compounds, to get N/screen
    from src.final_pharmacogenomics_gdsc_data import load_gdsc_fitted_response
    fitted = load_gdsc_fitted_response(cfg)

    targets = {
        "Tamoxifen": r"tamoxifen",
        "4-hydroxytamoxifen / 4-OHT": r"4-hydroxytamoxifen|4-OHT|\b4OHT\b",
        "Endoxifen": r"endoxifen",
        "Fulvestrant": r"fulvestrant",
    }
    rows = []
    for label, pattern in targets.items():
        hits = comp[comp["DRUG_NAME"].str.contains(pattern, case=False, na=False, regex=True) | comp["SYNONYMS"].str.contains(pattern, case=False, na=False, regex=True)]
        if len(hits) == 0:
            rows.append(dict(compound=label, present_in_gdsc=False, drug_id=None, drug_name_in_gdsc=None, target=None, pathway=None, screens=None, n_breast_lines_tested=None))
        else:
            for _, h in hits.iterrows():
                sub = fitted[(fitted["DRUG_ID"] == h["DRUG_ID"]) & (fitted["TCGA_DESC"] == "BRCA")]
                screens = sorted(sub["DATASET"].unique().tolist())
                rows.append(dict(
                    compound=label, present_in_gdsc=True, drug_id=int(h["DRUG_ID"]), drug_name_in_gdsc=h["DRUG_NAME"],
                    target=h["TARGET"], pathway=h["TARGET_PATHWAY"], screens=",".join(screens) if screens else "none (not tested in breast lines)",
                    n_breast_lines_tested=int(sub["SANGER_MODEL_ID"].nunique()),
                ))
    return pd.DataFrame(rows)


def _compute_associations(df: pd.DataFrame, gene: str, min_n: int, subset_label: str) -> pd.DataFrame:
    # Grouped strictly by (DRUG_ID, DATASET) -- the true pseudoreplication-safe
    # unit -- never by DRUG_NAME/TARGET/PATHWAY, which must not affect how
    # many statistical tests a given DRUG_ID contributes. Name/target/pathway
    # are looked up per group via .iloc[0] purely for display, after grouping.
    rows = []
    for (drug_id, dataset), g in df.groupby(["DRUG_ID", "DATASET"]):
        drug_name = g["DRUG_NAME"].iloc[0]
        target = g["PUTATIVE_TARGET"].iloc[0]
        pathway = g["PATHWAY_NAME"].iloc[0]
        assert g["DRUG_NAME"].nunique() == 1, f"DRUG_ID {drug_id} maps to >1 DRUG_NAME within {dataset}"
        for metric in ("LN_IC50", "AUC"):
            gg = g.dropna(subset=[gene, metric])
            n = len(gg)
            if n < min_n:
                continue
            rho, p = stats.spearmanr(gg[gene], gg[metric])
            rows.append(dict(
                gene=gene, subset=subset_label, drug_id=int(drug_id), drug_name=drug_name, target=target, pathway=pathway,
                dataset=dataset, metric=metric, n=n, spearman_rho=rho, p_value=p,
                direction="higher_expression_more_sensitive" if rho < 0 else "higher_expression_more_resistant",
            ))
    out = pd.DataFrame(rows)
    if len(out):
        out["fdr"] = np.nan
        for (ds, met), idx in out.groupby(["dataset", "metric"]).groups.items():
            out.loc[idx, "fdr"] = multipletests(out.loc[idx, "p_value"], method="fdr_bh")[1]
    return out.sort_values("p_value").reset_index(drop=True)


def build_gene_association_tables(joined: pd.DataFrame) -> dict[str, pd.DataFrame]:
    out = {}
    for gene in GENES:
        out[gene] = _compute_associations(joined, gene, MIN_N_FULL_BREAST, "all_breast_lines")
        logger.info("build_gene_association_tables: %s -- %d drug-tests (N>=%d breast lines), %d FDR<0.05",
                    gene, len(out[gene]), MIN_N_FULL_BREAST, int((out[gene]["fdr"] < 0.05).sum()))
    return out


def build_top_associations_table(gene_tables: dict[str, pd.DataFrame]) -> pd.DataFrame:
    rows = []
    for gene, df in gene_tables.items():
        sig = df[df["fdr"] < 0.05].copy()
        sig["tier"] = "FDR_SIGNIFICANT"
        top_effect = df.reindex(df["spearman_rho"].abs().sort_values(ascending=False).index).head(5).copy()
        top_effect["tier"] = "TOP_EFFECT_SIZE_NOT_NECESSARILY_SIGNIFICANT"
        rows.append(pd.concat([sig, top_effect], ignore_index=True).drop_duplicates(subset=["drug_id", "dataset", "metric"]))
    out = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()

    def pathway_group(pathway: str) -> str:
        for group, names in PATHWAY_COHERENCE_GROUPS.items():
            if pathway in names:
                return group
        return "other_ungrouped"
    if len(out):
        out["pathway_coherence_group"] = out["pathway"].apply(pathway_group)
    return out.sort_values(["gene", "fdr", "p_value"]).reset_index(drop=True)


def build_er_luminal_subset_table(joined: pd.DataFrame, top: pd.DataFrame) -> pd.DataFrame:
    """Exploratory-only re-test of the drugs flagged in build_top_associations_table,
    restricted to ER+/luminal lines, using the EXACT drug_id (never
    DRUG_NAME alone, which can silently pool multiple distinct GDSC
    re-screening batches of the same compound into one pseudoreplicated
    test -- caught and fixed during this analysis)."""
    erlum = joined[joined["is_er_luminal"]]
    n_lines = int(erlum["SANGER_MODEL_ID"].nunique())
    rows = []
    seen = set()
    for _, r in top.iterrows():
        key = (r["gene"], r["drug_id"], r["dataset"])
        if key in seen:
            continue
        seen.add(key)
        for metric in ("LN_IC50", "AUC"):
            gg = erlum[erlum["DATASET"] == r["dataset"]]
            gg = gg[gg["DRUG_ID"] == r["drug_id"]].dropna(subset=[r["gene"], metric])
            n = len(gg)
            if n < MIN_N_ER_LUMINAL_EXPLORATORY:
                rows.append(dict(gene=r["gene"], drug_name=r["drug_name"], drug_id=r["drug_id"], dataset=r["dataset"], metric=metric,
                                  n_er_luminal_lines=n, n_total_er_luminal_lines_in_panel=n_lines, spearman_rho=None, p_value=None,
                                  exploratory_only=True, note="N too small to compute (below MIN_N_ER_LUMINAL_EXPLORATORY)"))
                continue
            rho, p = stats.spearmanr(gg[r["gene"]], gg[metric])
            rows.append(dict(gene=r["gene"], drug_name=r["drug_name"], drug_id=r["drug_id"], dataset=r["dataset"], metric=metric,
                              n_er_luminal_lines=n, n_total_er_luminal_lines_in_panel=n_lines, spearman_rho=rho, p_value=p,
                              exploratory_only=True, note="EXPLORATORY ONLY -- small N, not independently FDR-corrected within this subset, not a validation of the full-breast-set finding"))
    return pd.DataFrame(rows)


def _direction_from_log2fc(log2fc) -> str | None:
    if pd.isna(log2fc):
        return None
    return "up" if log2fc > 0 else "down"


def build_project_crosscheck_table(top: pd.DataFrame) -> pd.DataFrame:
    cross = pd.read_csv(CROSS_DATASET_TABLE, sep="\t").set_index("gene")
    targets = set()
    for t in top["target"].dropna():
        for g in str(t).split(","):
            targets.add(g.strip())
    rows = []
    for t in sorted(targets):
        if t in cross.index:
            live = cross.loc[t]
            rows.append(dict(
                drug_target_gene=t, present_in_project_crossdataset_table=True,
                hany_crispr_fdr=float(live["crispr_fdr"]), hany_crispr_direction=live["crispr_direction"],
                gse118713_fdr=float(live["gse118713_fdr"]) if pd.notna(live["gse118713_fdr"]) else None,
                gse118713_direction=live.get("gse118713_direction") if pd.notna(live.get("gse118713_direction")) else None,
                gse240112_tumor_fdr=float(live["gse240112_tumor_fdr"]) if pd.notna(live["gse240112_tumor_fdr"]) else None,
                gse240112_tumor_direction=_direction_from_log2fc(live.get("gse240112_tumor_log2fc")),
                gse111151_fdr=float(live["gse111151_fdr"]) if "gse111151_fdr" in cross.columns and pd.notna(live.get("gse111151_fdr")) else None,
                gse111151_direction=_direction_from_log2fc(live.get("gse111151_log2fc")),
                gse245601_epi_fdr=float(live["gse245601_epi_fdr"]) if pd.notna(live["gse245601_epi_fdr"]) else None,
                gse245601_epi_direction=_direction_from_log2fc(live.get("gse245601_epi_log2fc")),
            ))
        else:
            rows.append(dict(drug_target_gene=t, present_in_project_crossdataset_table=False, hany_crispr_fdr=None, hany_crispr_direction=None,
                              gse118713_fdr=None, gse118713_direction=None, gse240112_tumor_fdr=None, gse240112_tumor_direction=None,
                              gse111151_fdr=None, gse111151_direction=None, gse245601_epi_fdr=None, gse245601_epi_direction=None))
    out = pd.DataFrame(rows)
    out["notes"] = ("Frozen shortlist unaltered by this lookup; DepMap dependency not looked up here (that analysis is scoped to the 4 frozen "
                     "candidates only, not drug-target genes). Direction columns report only the SIGN of that dataset's own log2FC/effect for "
                     "the drug-target gene itself (up/down in resistant vs sensitive, or sensitising/tolerance for Hany) -- they say nothing "
                     "about whether that gene's baseline expression drives GDSC drug response; any apparent convergence with the GDSC finding "
                     "above is a correlational echo across independent datasets, not a validated mechanistic link.")
    return out


def build_final_interpretation_table() -> pd.DataFrame:
    return pd.DataFrame(FINAL_INTERPRETATION_ROWS)


def build_indirect_targeting_table() -> pd.DataFrame:
    return pd.DataFrame(DRUG_INDIRECT_TARGETING_ROWS)


def run(out_dir: Path = OUT_DIR) -> dict[str, pd.DataFrame]:
    out_dir.mkdir(parents=True, exist_ok=True)
    cfg = load_config()

    provenance = build_data_provenance_table()
    compound_avail = build_compound_availability_table()
    joined = build_breast_expression_joined(cfg)
    gene_tables = build_gene_association_tables(joined)
    top = build_top_associations_table(gene_tables)
    er_luminal = build_er_luminal_subset_table(joined, top)
    crosscheck = build_project_crosscheck_table(top)

    tables = {
        "GDSC_data_provenance.tsv": provenance,
        "GDSC_compound_availability.tsv": compound_avail,
        "USP34_GDSC_drug_associations.tsv": gene_tables["USP34"],
        "VEZF1_GDSC_drug_associations.tsv": gene_tables["VEZF1"],
        "GDSC_top_associations.tsv": top,
        "GDSC_ER_luminal_subset.tsv": er_luminal,
        "GDSC_project_crosscheck.tsv": crosscheck,
        "GDSC_indirect_targeting_classification.tsv": build_indirect_targeting_table(),
        "GDSC_final_interpretation.tsv": build_final_interpretation_table(),
    }
    for name, df in tables.items():
        path = out_dir / name
        df.to_csv(path, sep="\t", index=False)
        logger.info("wrote %s (%d rows)", path, len(df))
    return tables


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run()
