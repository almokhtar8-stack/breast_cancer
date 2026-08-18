"""post_freeze_exploratory -- feasibility verdict for the GSE245601
annotation-concordance analysis (Task 3).

The requested analysis was a confusion matrix of the authors' published
Epi. Tumor / Epi. Nontumor per-cell labels against our inferCNV malignant /
non-malignant calls, with Cohen's kappa overall and per tumour. That analysis
requires a public barcode-to-label mapping from the authors. This module turns
the record written by
``scripts/post_poster_probe_gse245601_annotations.py`` into a verdict and a
small evidence table.

It is deliberately incapable of producing a concordance matrix. Reconstructing
the authors' labels in order to compare them against our reconstruction of the
authors' labels would be circular -- the comparison would measure agreement
between two runs of our own pipeline, not agreement with the paper -- so if the
labels are not public the correct output is a documented stop, not a proxy.

Deterministic; no network access (the probe does the fetching).
"""

from __future__ import annotations

import json
import logging
from pathlib import Path

import pandas as pd

logger = logging.getLogger(__name__)

ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = ROOT / "results" / "post_poster" / "annotation_concordance"
PROBE_PATH = OUT_DIR / "feasibility_probe.json"

POST_FREEZE_LABEL = "post_freeze_exploratory"

# Our frozen per-cell calls, which WOULD have been one side of the matrix.
OUR_LABELS_PATH = ROOT / "results" / "tables" / "gse245601_malignant_cell_labels.tsv"


def load_probe(path: Path = PROBE_PATH) -> dict:
    if not path.exists():
        raise FileNotFoundError(
            f"{path} not found -- run scripts/post_poster_probe_gse245601_annotations.py first")
    return json.loads(path.read_text())


def evaluate(probe: dict) -> pd.DataFrame:
    """One row per source checked, each with an explicit obtainable yes/no.

    A source counts as supplying the labels only if it exposes a per-cell
    barcode-to-label mapping. Gene-level supplementary tables, figure panels
    and analysis code do not.
    """
    geo = probe["geo"]
    supp = probe["paper_supplement"]
    repo = probe["author_repo"]
    dbgap = probe["controlled_access_route"]

    rows = [
        {
            "source": "GEO series + samples",
            "identifier": geo["series"],
            "what_exists": (f"{len(geo['series_supplementary_files'])} series supplementary "
                            f"file(s); {geo['n_sample_supplementary_files']} sample "
                            f"supplementary file(s), extensions "
                            f"{','.join(geo['sample_supplementary_file_extensions'])}"),
            "per_cell_labels_present": False,
            "evidence": ("only Cell Ranger filtered feature-barcode matrices; a label table "
                         "would be a separate csv/tsv/rds supplementary file and none exists"),
        },
        {
            "source": "Paper supplementary material",
            "identifier": supp["pmc_id"],
            "what_exists": (f"{len(supp['files'])} supplementary file(s); "
                            f"{supp['n_xlsx_sheets_scanned']} spreadsheet sheets scanned, "
                            f"largest {supp['max_sheet_rows']} rows"),
            "per_cell_labels_present": bool(supp["total_barcode_like_cells"]),
            "evidence": (f"{supp['total_barcode_like_cells']} cells matching a 10x barcode "
                         "across every sheet; all sheets are gene-level differential "
                         "expression, enrichment or signature tables"),
        },
        {
            "source": "Authors' code repository",
            "identifier": f"{repo['repo']}@{repo['pinned_commit'][:7]}",
            "what_exists": (f"{repo['n_blobs_at_pinned_commit']} files at the pinned commit; "
                            f"{repo['n_commits_all_refs']} commits over all refs; "
                            f"{repo['n_releases']} releases, {repo['n_tags']} tags; "
                            f"LFS configured: {repo['has_gitattributes_lfs_config']}; "
                            f"extensions ever committed: "
                            f"{','.join(repo['file_extensions_ever_committed'])}"),
            "per_cell_labels_present": False,
            "evidence": (f"analysis code and executed notebooks only; the largest barcode-like "
                         f"match count in any notebook is "
                         f"{repo['max_notebook_barcode_matches']}, which are truncated "
                         "dataframe previews, not an exported label table"),
        },
        {
            "source": "dbGaP (controlled access)",
            "identifier": dbgap["accession"],
            "what_exists": "processed scRNA-seq data under controlled access",
            "per_cell_labels_present": None,
            "evidence": dbgap["note"],
        },
    ]
    df = pd.DataFrame(rows)
    df["probe_date"] = probe["probe_date"]
    df[POST_FREEZE_LABEL] = True
    return df


def verdict(sources: pd.DataFrame) -> dict:
    """Feasibility gate result. Passes only if some PUBLIC source supplies
    per-cell labels."""
    public = sources[sources["source"] != "dbGaP (controlled access)"]
    obtainable = bool(public["per_cell_labels_present"].eq(True).any())
    return {
        "analysis": "GSE245601 published-vs-reconstructed malignant-call concordance",
        "feasibility_gate": "PASS" if obtainable else "FAIL",
        "per_cell_labels_publicly_obtainable": obtainable,
        "action_taken": ("proceed with concordance matrix" if obtainable else
                         "STOP -- no confusion matrix, kappa or per-tumour concordance "
                         "produced; no proxy comparison attempted"),
        "why_no_proxy": ("reconstructing the authors' labels in order to compare against our "
                         "reconstruction of the authors' labels would measure agreement "
                         "between two runs of our own pipeline, not agreement with the paper"),
        "n_public_sources_checked": int(len(public)),
        "probe_date": str(sources["probe_date"].iloc[0]),
        POST_FREEZE_LABEL: True,
    }


def our_label_inventory(path: Path = OUR_LABELS_PATH) -> pd.DataFrame:
    """What our side of the (unbuildable) matrix would have contributed.

    Reported so the report can state the size and shape of our own calls
    without implying any comparison was made. Frozen labels are read only.
    """
    df = pd.read_csv(path, sep="\t")
    label_col = "primary_malignancy_label"
    counts = (df.groupby(["patient", "condition", label_col]).size()
              .unstack(fill_value=0).reset_index())
    label_cols = [c for c in counts.columns if c not in ("patient", "condition")]
    counts["n_cells"] = counts[label_cols].sum(axis=1)
    counts[POST_FREEZE_LABEL] = True
    logger.info("our frozen labels: %d cells across %d tumours; labels present: %s",
                len(df), df["patient"].nunique(), sorted(df[label_col].unique()))
    return counts


def main(out_dir: Path = OUT_DIR) -> None:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    out_dir.mkdir(parents=True, exist_ok=True)

    probe = load_probe()
    sources = evaluate(probe)
    sources.to_csv(out_dir / "feasibility_sources_checked.tsv", sep="\t", index=False)
    logger.info("wrote feasibility_sources_checked.tsv (%d sources)", len(sources))

    v = verdict(sources)
    pd.DataFrame([v]).to_csv(out_dir / "feasibility_verdict.tsv", sep="\t", index=False)
    logger.info("feasibility gate: %s -- %s", v["feasibility_gate"], v["action_taken"])

    inventory = our_label_inventory()
    inventory.to_csv(out_dir / "our_frozen_label_inventory.tsv", sep="\t", index=False)
    logger.info("wrote our_frozen_label_inventory.tsv (%d tumours)", len(inventory))

    if v["feasibility_gate"] == "PASS":  # pragma: no cover - not reachable today
        raise NotImplementedError(
            "Per-cell labels became publicly available. The concordance matrix, Cohen's "
            "kappa and per-tumour breakdown must be implemented before this module can "
            "report a result -- it deliberately does not guess one.")


if __name__ == "__main__":
    main()
