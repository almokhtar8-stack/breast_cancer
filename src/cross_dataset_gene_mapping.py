"""Cross-dataset genome-wide integration, Phases 2-3: load each of the
five independent datasets' genome-wide result tables, resolve a known,
benign class of duplicate-symbol rows (pseudoautosomal-region ``ENSGR``
Ensembl duplicates, which carry numerically identical statistics (parsed-value equality, not raw-byte equality) to their
canonical ``ENSG`` counterpart), flag any *genuinely* ambiguous
duplicated symbol (two distinct Ensembl genes sharing one symbol) as
unresolved rather than guessing, and build the union gene universe.

No gene list from any prior phase of this project (the 13 sensitising
candidates, the 28 CRISPR-significant hits, PAICS, or any
individually-discussed gene) is used to select, seed, or filter this
universe -- it is the union of every gene symbol appearing as a testable
feature in any of the five datasets' own frozen output files.

Data sources: `data/processed/labels.parquet` (CRISPR screen),
`results/tables/gse118713_differential_expression_unredacted.tsv.gz`,
`results/tables/gse245601_pseudobulk/track_{a,b}_genomewide_de.tsv.gz`,
`results/tables/gse240112_pseudobulk/{tumor_cell,epithelial}_genomewide_de.tsv.gz`,
`results/tables/gse111151/genomewide_de.tsv.gz`. All already frozen and
committed by prior phases of this project; read-only here.
"""

from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd
import yaml

logger = logging.getLogger(__name__)

DATASET_NAMES = ["crispr", "gse118713", "gse245601", "gse240112", "gse111151"]


def _load_config(config_path: str | Path) -> dict:
    with open(config_path) as f:
        return yaml.safe_load(f)


def load_crispr_screen(path: str | Path) -> pd.DataFrame:
    """CRISPR screen: 19,103 genes, symbol-only identifier, no filtering
    beyond the guide-level model fit itself -- every row is testable."""
    df = pd.read_parquet(path)
    out = df.rename(columns={"gene": "symbol", "effect_size": "effect", "p_value": "p_value", "fdr": "fdr"})
    out = out[["symbol", "effect", "p_value", "fdr", "n_guides"]].copy()
    if out["symbol"].duplicated().any():
        raise ValueError("unexpected duplicate gene symbol in the CRISPR screen table")
    logger.info("load_crispr_screen: %d genes, all testable", len(out))
    return out


def resolve_ensgr_duplicates(df: pd.DataFrame, id_col: str, symbol_col: str) -> tuple[pd.DataFrame, list[str]]:
    """Drop the pseudoautosomal-region (``ENSGR``-prefixed) duplicate row
    for any gene where an ``ENSG``-prefixed row with the *same symbol and
    numerically identical statistics (parsed-value equality, not raw-byte equality)* also exists -- these represent the same
    measured locus annotated twice (once on X, once on the Y PAR), not a
    genuine ambiguity. Returns (cleaned_df, list_of_symbols_resolved)."""
    is_ensgr = df[id_col].str.startswith("ENSGR")
    resolved_symbols: list[str] = []
    if not is_ensgr.any():
        return df, resolved_symbols

    keep_mask = pd.Series(True, index=df.index)
    stat_cols = [c for c in df.columns if c not in (id_col, symbol_col)]
    for symbol, group in df.loc[df[symbol_col].isin(df.loc[is_ensgr, symbol_col])].groupby(symbol_col):
        ensgr_rows = group.loc[group[id_col].str.startswith("ENSGR")]
        ensg_rows = group.loc[~group[id_col].str.startswith("ENSGR")]
        if len(ensgr_rows) == 1 and len(ensg_rows) == 1:
            same_stats = (ensgr_rows[stat_cols].to_numpy() == ensg_rows[stat_cols].to_numpy()).all()
            if same_stats:
                keep_mask.loc[ensgr_rows.index] = False
                resolved_symbols.append(symbol)

    out = df.loc[keep_mask].reset_index(drop=True)
    logger.info("resolve_ensgr_duplicates: resolved %d pseudoautosomal-duplicate symbols (dropped the ENSGR row, kept ENSG)", len(resolved_symbols))
    return out, resolved_symbols


def flag_ambiguous_symbols(df: pd.DataFrame, symbol_col: str) -> tuple[pd.DataFrame, list[str]]:
    """After ENSGR resolution, any symbol still appearing on >1 row maps
    to genuinely distinct genes (different Ensembl IDs, different -- or
    at least not verified identical -- statistics). These are never
    collapsed by picking one row arbitrarily: they are returned
    separately as ``ambiguous_symbols`` and excluded from the
    single-row-per-symbol table this function returns."""
    dup_symbols = df.loc[df[symbol_col].duplicated(keep=False), symbol_col].unique().tolist()
    clean = df.loc[~df[symbol_col].isin(dup_symbols)].reset_index(drop=True)
    logger.info("flag_ambiguous_symbols: %d symbols remain genuinely ambiguous (multiple distinct genes), excluded from the single-row table", len(dup_symbols))
    return clean, sorted(dup_symbols)


def load_gse118713(path: str | Path, contrast: str = "TAMR_vs_MCF7") -> tuple[pd.DataFrame, list[str], list[str]]:
    """Returns (clean_df, ensgr_resolved_symbols, ambiguous_symbols).
    Only the primary contrast is loaded -- TAMR_vs_FASR and FASR_vs_MCF7
    are secondary context, not part of this dataset's one independent
    contribution (docs/CROSS_DATASET_GENOMEWIDE_DATA_AUDIT.md)."""
    df = pd.read_csv(path, sep="\t")
    df = df.loc[df["contrast"] == contrast].copy()
    df = df.rename(columns={"gene_symbol": "symbol", "gene_id": "ensembl_id", "log2fc": "effect"})
    df = df[["ensembl_id", "symbol", "effect", "p_value", "fdr"]]
    resolved, ensgr_resolved = resolve_ensgr_duplicates(df, "ensembl_id", "symbol")
    clean, ambiguous = flag_ambiguous_symbols(resolved, "symbol")
    logger.info("load_gse118713[%s]: %d genes clean, %d ambiguous", contrast, len(clean), len(ambiguous))
    return clean, ensgr_resolved, ambiguous


def load_gse245601_track(path: str | Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    out = df.rename(columns={"log2fc": "effect"})[["gene", "effect", "p_value", "fdr"]].rename(columns={"gene": "symbol"})
    if out["symbol"].duplicated().any():
        raise ValueError("unexpected duplicate gene symbol in a GSE245601 track table")
    return out


def load_gse240112_track(path: str | Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    out = df.rename(columns={"log2fc": "effect"})[["gene", "effect", "p_value", "fdr"]].rename(columns={"gene": "symbol"})
    if out["symbol"].duplicated().any():
        raise ValueError("unexpected duplicate gene symbol in a GSE240112 track table")
    return out


def load_gse111151(path: str | Path) -> tuple[pd.DataFrame, list[str], list[str]]:
    df = pd.read_csv(path, sep="\t")
    df = df.rename(columns={"gene_name": "symbol", "gene_id": "ensembl_id", "log2fc": "effect"})
    df = df[["ensembl_id", "symbol", "effect", "p_value", "fdr"]]
    resolved, ensgr_resolved = resolve_ensgr_duplicates(df, "ensembl_id", "symbol")
    clean, ambiguous = flag_ambiguous_symbols(resolved, "symbol")
    logger.info("load_gse111151: %d genes clean, %d ambiguous", len(clean), len(ambiguous))
    return clean, ensgr_resolved, ambiguous


def build_gene_mapping_audit(
    gse118713_ensgr: list[str], gse118713_ambiguous: list[str], gse111151_ensgr: list[str], gse111151_ambiguous: list[str]
) -> pd.DataFrame:
    """One row per (dataset, symbol, resolution) event -- every mapping
    exception is documented explicitly, per the task's Phase 3 rules."""
    rows = []
    for sym in gse118713_ensgr:
        rows.append({"dataset": "gse118713", "symbol": sym, "issue": "pseudoautosomal_ENSGR_duplicate", "resolution": "dropped_ENSGR_row_kept_ENSG_row_identical_statistics", "excluded_from_testable": False})
    for sym in gse118713_ambiguous:
        rows.append({"dataset": "gse118713", "symbol": sym, "issue": "duplicate_symbol_multiple_distinct_ensembl_ids", "resolution": "not_collapsed_excluded_from_single_row_table", "excluded_from_testable": True})
    for sym in gse111151_ensgr:
        rows.append({"dataset": "gse111151", "symbol": sym, "issue": "pseudoautosomal_ENSGR_duplicate", "resolution": "dropped_ENSGR_row_kept_ENSG_row_identical_statistics", "excluded_from_testable": False})
    for sym in gse111151_ambiguous:
        rows.append({"dataset": "gse111151", "symbol": sym, "issue": "duplicate_symbol_multiple_distinct_ensembl_ids", "resolution": "not_collapsed_excluded_from_single_row_table", "excluded_from_testable": True})
    out = pd.DataFrame(rows)
    logger.info("build_gene_mapping_audit: %d mapping-exception rows", len(out))
    return out


def build_full_gene_universe(datasets: dict[str, pd.DataFrame]) -> pd.DataFrame:
    """``datasets``: mapping of dataset_name -> DataFrame with a
    ``symbol`` column (already ambiguity-resolved, one row per gene =
    testable in that dataset). Returns the union of all symbols with
    presence/testable flags per dataset (identical here by construction:
    every symbol present in a loaded table is, by the loaders above,
    testable -- documented in the data audit as the honest limit of what
    the frozen per-dataset outputs distinguish)."""
    all_symbols: set[str] = set()
    testable_sets: dict[str, set[str]] = {}
    for name, df in datasets.items():
        symbols = set(df["symbol"])
        testable_sets[name] = symbols
        all_symbols |= symbols

    rows = []
    for symbol in sorted(all_symbols):
        present_testable = {name: (symbol in testable_sets[name]) for name in datasets}
        n_present = sum(present_testable.values())
        row = {"gene": symbol, "n_datasets_present": n_present, "n_datasets_testable": n_present}
        for name in datasets:
            row[f"{name}_present"] = present_testable[name]
            row[f"{name}_testable"] = present_testable[name]
        rows.append(row)
    out = pd.DataFrame(rows)
    logger.info("build_full_gene_universe: %d unique genes across %d datasets", len(out), len(datasets))
    return out


def run_gene_mapping(config_path: str | Path = "config/config.yaml") -> dict[str, pd.DataFrame]:
    config = _load_config(config_path)
    cfg = config["cross_dataset_genomewide"]
    inputs = cfg["inputs"]
    out = cfg["output"]

    crispr = load_crispr_screen(inputs["crispr_labels_parquet"])
    gse118713, g118713_ensgr, g118713_amb = load_gse118713(inputs["gse118713_de_tsv"])
    gse245601_a = load_gse245601_track(inputs["gse245601_track_a_tsv"])
    gse245601_b = load_gse245601_track(inputs["gse245601_track_b_tsv"])
    gse240112_tumor = load_gse240112_track(inputs["gse240112_tumor_cell_tsv"])
    gse240112_epi = load_gse240112_track(inputs["gse240112_epithelial_tsv"])
    gse111151, g111151_ensgr, g111151_amb = load_gse111151(inputs["gse111151_de_tsv"])

    # each dataset's ONE independent-contribution symbol set for the universe
    # (GSE245601 = union of Track A/B testable symbols since either track
    # testing a gene makes it testable in "the GSE245601 dataset"; same
    # logic for GSE240112's two tracks)
    gse245601_symbols = pd.DataFrame({"symbol": sorted(set(gse245601_a["symbol"]) | set(gse245601_b["symbol"]))})
    gse240112_symbols = pd.DataFrame({"symbol": sorted(set(gse240112_tumor["symbol"]) | set(gse240112_epi["symbol"]))})

    universe_inputs = {
        "crispr": crispr[["symbol"]],
        "gse118713": gse118713[["symbol"]],
        "gse245601": gse245601_symbols,
        "gse240112": gse240112_symbols,
        "gse111151": gse111151[["symbol"]],
    }
    universe = build_full_gene_universe(universe_inputs)
    mapping_audit = build_gene_mapping_audit(g118713_ensgr, g118713_amb, g111151_ensgr, g111151_amb)

    Path(out["full_gene_universe_tsv"]).parent.mkdir(parents=True, exist_ok=True)
    universe.to_csv(out["full_gene_universe_tsv"], sep="\t", index=False)
    mapping_audit.to_csv(out["gene_mapping_audit_tsv"], sep="\t", index=False)
    logger.info("wrote %s and %s", out["full_gene_universe_tsv"], out["gene_mapping_audit_tsv"])

    return {
        "crispr": crispr,
        "gse118713": gse118713,
        "gse245601_track_a": gse245601_a,
        "gse245601_track_b": gse245601_b,
        "gse240112_tumor": gse240112_tumor,
        "gse240112_epi": gse240112_epi,
        "gse111151": gse111151,
        "universe": universe,
        "mapping_audit": mapping_audit,
    }


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run_gene_mapping()
