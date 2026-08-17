"""Tests for the POST-AUDIT SENSITIVITY ANALYSIS (external-audit response).

These tests pin the specific, load-bearing facts that the sensitivity
report depends on: the original frozen gate is reproduced bit-for-bit by
Rule 0 (never re-implemented), the significant-sensitising universe is the
full pre-specified set (not an arbitrary subset), KDM1A/TLK2 rank ahead of
USP34 on pure CRISPR strength, and every dataset-separation rule (acute
vs chronic, GSE111151 blocking) is respected.
"""

from __future__ import annotations

import re
from pathlib import Path

import pandas as pd
import pytest

from src import post_audit_sensitivity_data as pad

FIGURES = Path("results/figures/post_audit")
TABLES = Path("results/tables/post_audit_sensitivity")


class TestSignificantSensitisingUniverse:
    def test_13_significant_sensitising_hits_at_the_prespecified_threshold(self):
        sens = pad.load_significant_sensitising_hits()
        assert len(sens) == 13
        assert (sens["fdr"] < 0.1).all()
        assert (sens["effect_size"] < 0).all()

    def test_focus_four_are_all_present_in_the_universe(self):
        sens = set(pad.load_significant_sensitising_hits()["gene"])
        assert set(pad.FOCUS_FOUR).issubset(sens)

    def test_gate1_has_28_hits_13_sensitising_15_tolerance(self):
        gate1 = pad.load_gate1_full()
        assert len(gate1) == 28
        assert (gate1["effect_size"] < 0).sum() == 13
        assert (gate1["effect_size"] > 0).sum() == 15

    def test_kdm1a_ranks_first_by_both_effect_and_fdr(self):
        sens = pad.load_significant_sensitising_hits().set_index("gene")
        assert sens.loc["KDM1A", "rank_by_effect"] == 1
        assert sens.loc["KDM1A", "rank_by_fdr"] == 1

    def test_usp34_ranks_worse_than_kdm1a_and_tlk2_on_both_crispr_metrics(self):
        sens = pad.load_significant_sensitising_hits().set_index("gene")
        for metric in ("rank_by_effect", "rank_by_fdr"):
            assert sens.loc["USP34", metric] > sens.loc["KDM1A", metric]
            assert sens.loc["USP34", metric] > sens.loc["TLK2", metric]


class TestRule0ReproducesTheFrozenShortlistExactly:
    def test_rule0_gene_order_matches_frozen_freeze_table(self):
        rule0 = pad.rule0_original_frozen_gate()
        frozen = pd.read_csv(pad.FREEZE_TABLE, sep="\t")
        assert rule0["gene"].tolist() == frozen.sort_values("freeze_rank")["gene"].tolist() == ["USP34", "VEZF1", "EML5", "CITED2"]

    def test_rule0_crispr_values_match_frozen_freeze_table_exactly(self):
        rule0 = pad.rule0_original_frozen_gate().set_index("gene")
        frozen = pd.read_csv(pad.FREEZE_TABLE, sep="\t").set_index("gene")
        for gene in ["USP34", "VEZF1", "EML5", "CITED2"]:
            assert rule0.loc[gene, "crispr_effect"] == pytest.approx(frozen.loc[gene, "crispr_effect"], abs=1e-9)
            assert rule0.loc[gene, "crispr_fdr"] == pytest.approx(frozen.loc[gene, "crispr_fdr"], abs=1e-9)

    def test_rule0_calls_the_unmodified_original_freeze_code_not_a_reimplementation(self):
        import inspect
        src = inspect.getsource(pad.rule0_original_frozen_gate)
        assert "from src.evidence_freeze_shortlist_freeze import build_freeze" in src


class TestKdm1aAndTlk2ExclusionReason:
    def test_freeze_eligibility_audit_confirms_functional_only_exclusion(self):
        audit = pd.read_csv("results/tables/evidence_freeze/freeze_eligibility_audit.tsv", sep="\t").set_index("gene")
        for gene in ["KDM1A", "TLK2"]:
            assert audit.loc[gene, "eligible_for_freeze"] == False  # noqa: E712
            assert "no resistance-RNA dataset FDR<0.05" in audit.loc[gene, "ineligibility_reason"]
            assert audit.loc[gene, "_crispr_band"] == "VERY_STRONG"


class TestRule2ChronicCorroboration:
    def test_only_usp34_and_vezf1_pass_chronic_rna_gate(self):
        eligible = set(pad.rule2_chronic_corroboration()["gene"])
        assert eligible == {"USP34", "VEZF1"}

    def test_usp34_sole_chronic_hit_is_gse118713(self):
        em = pad.build_evidence_matrix().set_index("gene")
        assert em.loc["USP34", "gse118713_fdr"] < 0.05
        assert em.loc["USP34", "gse240112_fdr"] >= 0.05
        assert em.loc["USP34", "gse111151_fdr"] >= 0.05

    def test_vezf1_sole_chronic_hit_is_gse240112(self):
        em = pad.build_evidence_matrix().set_index("gene")
        assert em.loc["VEZF1", "gse240112_fdr"] < 0.05
        assert em.loc["VEZF1", "gse118713_fdr"] >= 0.05
        assert em.loc["VEZF1", "gse111151_fdr"] >= 0.05


class TestRule3Gse111151Empty:
    def test_zero_of_13_genes_reach_gse111151_fdr05(self):
        eligible = pad.rule3_gse111151_specific()
        assert len(eligible) == 0


class TestRule4HumanEvidenceFirst:
    def test_vezf1_outranks_usp34_when_human_evidence_ordered_before_cellline_rna(self):
        ranked = pad.rule4_human_evidence_first().set_index("gene")
        assert ranked.loc["VEZF1", "rank"] < ranked.loc["USP34", "rank"]

    def test_kdm1a_and_tlk2_still_outrank_both_frozen_candidates_under_rule4(self):
        ranked = pad.rule4_human_evidence_first().set_index("gene")
        for weaker in ("USP34", "VEZF1"):
            assert ranked.loc["KDM1A", "rank"] < ranked.loc[weaker, "rank"]
            assert ranked.loc["TLK2", "rank"] < ranked.loc[weaker, "rank"]


class TestLeaveOneDatasetOut:
    def test_usp34_loses_all_corroboration_without_gse118713(self):
        loo = pad.build_leave_one_dataset_out()
        row = loo[(loo["gene"] == "USP34") & (loo["left_out_dataset"] == "gse118713")].iloc[0]
        assert not row["still_corroborated"]

    def test_vezf1_loses_all_corroboration_without_gse240112(self):
        loo = pad.build_leave_one_dataset_out()
        row = loo[(loo["gene"] == "VEZF1") & (loo["left_out_dataset"] == "gse240112")].iloc[0]
        assert not row["still_corroborated"]

    def test_gse245601_never_appears_as_a_leave_out_dataset(self):
        loo = pad.build_leave_one_dataset_out()
        assert "gse245601" not in set(loo["left_out_dataset"])


class TestNoMasterScore:
    def test_data_module_never_computes_a_weighted_composite(self):
        src = Path("src/post_audit_sensitivity_data.py").read_text()
        forbidden = ["master_score", "composite_score", "weighted_score", "overall_score"]
        for term in forbidden:
            assert term not in src.lower().replace("_", "_")  # simple substring guard

    def test_evidence_matrix_never_fills_missing_evidence_with_zero(self):
        em = pad.build_evidence_matrix()
        # TCGA/pathway/GDSC are only computed for the original 4 -- the
        # other 9 genes must show real NaN, never a filled 0
        non_original = em[~em["gene"].isin(pad.ORIGINAL_FOUR)]
        assert non_original["tcga_fdr"].isna().all()
        assert non_original["n_strong_consensus_pathways"].isna().all()


class TestGdscAudit:
    def test_9_significant_rows_but_8_unique_compounds(self):
        u = pd.read_csv("results/tables/final_pharmacogenomics/USP34_GDSC_drug_associations.tsv", sep="\t")
        v = pd.read_csv("results/tables/final_pharmacogenomics/VEZF1_GDSC_drug_associations.tsv", sep="\t")
        both = pd.concat([u, v])
        sig = both[both["fdr"] < 0.05]
        assert len(sig) == 9
        assert sig["drug_id"].nunique() == 8

    def test_all_significant_rows_are_gdsc1_none_gdsc2(self):
        u = pd.read_csv("results/tables/final_pharmacogenomics/USP34_GDSC_drug_associations.tsv", sep="\t")
        v = pd.read_csv("results/tables/final_pharmacogenomics/VEZF1_GDSC_drug_associations.tsv", sep="\t")
        sig = pd.concat([u, v])
        sig = sig[sig["fdr"] < 0.05]
        assert set(sig["dataset"]) == {"GDSC1"}

    def test_zero_vezf1_significant_gdsc_associations(self):
        v = pd.read_csv("results/tables/final_pharmacogenomics/VEZF1_GDSC_drug_associations.tsv", sep="\t")
        assert (v["fdr"] < 0.05).sum() == 0


class TestDocumentationCorrection:
    def test_gse240112_no_longer_described_as_matched_in_data_provenance(self):
        text = Path("docs/DATA_PROVENANCE.md").read_text()
        gse240112_section = text.split("## GSE240112")[1].split("## GSE111151")[0]
        assert "matched primary/recurrent" not in gse240112_section
        assert "unpaired" in gse240112_section.lower()

    def test_gse240112_no_longer_described_as_matched_in_readme(self):
        # a second, previously-missed occurrence of the same error was
        # found in the README's dataset table during the science-freeze
        # pass -- this guards against a third recurrence
        text = Path("README.md").read_text()
        gse240112_rows = [line for line in text.splitlines() if line.strip().startswith("| GSE240112")]
        assert gse240112_rows, "no GSE240112 row found in README.md dataset table"
        for row in gse240112_rows:
            assert "matched primary" not in row.lower()
            assert "unpaired" in row.lower()

    def test_no_matched_paired_gse240112_language_anywhere_in_markdown_docs(self):
        # repo-wide sweep: the only legitimate "matched" mentions for
        # GSE240112 concern matched scRNA+scATAC assay modalities per
        # sample, never matched/paired PATIENTS -- guard against a patient
        # -pairing claim reappearing anywhere. Operates on whole PARAGRAPHS
        # (blank-line-delimited blocks), not single physical lines -- a
        # prior version of this test only checked single lines and missed
        # a real occurrence that was wrapped onto a markdown continuation
        # line without the literal string "240112" on it (caught by a
        # Codex final-review pass and fixed here).
        offenders = []
        skip_paths = {Path("results/reports/post_audit/POST_AUDIT_SENSITIVITY_REPORT.md"),
                      Path("results/reports/post_audit/SCIENCE_FREEZE_REPORT.md")}  # these quote the historical error as audit history
        for md in list(Path(".").glob("*.md")) + list(Path("docs").glob("*.md")) + list(Path("results/reports").rglob("*.md")):
            if md in skip_paths:
                continue
            text = md.read_text(errors="ignore")
            for para in text.split("\n\n"):
                # markdown tables have no blank line between rows -- each
                # row is its own logical unit and must not be merged with
                # unrelated rows (e.g. a GSE245601 row's "paired primary"
                # wrongly attributed to a GSE240112 row in the same table)
                units = para.split("\n") if any(ln.strip().startswith("|") for ln in para.split("\n")) else [para.replace("\n", " ")]
                for unit in units:
                    low = unit.lower()
                    if "240112" not in low:
                        continue
                    # \b word boundaries so "unpaired primary" does not
                    # false-positive on the embedded substring "paired primary"
                    if re.search(r"\bmatched primary\b|\bmatched recurrence\b|\bprimary/recurrent pairs\b"
                                 r"|\bpaired primary\b|\brecurrent tumou?r pairs\b", low):
                        offenders.append(f"{md}: {unit.strip()[:200]}")
        assert offenders == [], f"found patient-pairing language for GSE240112: {offenders}"


class TestEnvironmentSpecCompleteness:
    def test_environment_yml_declares_every_imported_third_party_package(self):
        import ast
        import sys

        text = Path("environment.yml").read_text()
        declared = set()
        for line in text.splitlines():
            line = line.strip().lstrip("-").strip()
            if not line or line.startswith("#"):
                continue
            name = line.split("=")[0].split(">")[0].split("<")[0].strip()
            declared.add(name)
        # a small alias map for conda-package-name vs import-name mismatches
        alias = {"pyyaml": "yaml", "bioconductor-edger": None, "bioconductor-limma": None,
                  "r-base": None, "r-statmod": None, "scikit-learn": "sklearn", "pip": None,
                  "python": None, "pymol-open-source": "pymol", "biopython": "Bio"}
        declared_imports = set()
        for d in declared:
            declared_imports.add(alias.get(d, d))
        declared_imports.discard(None)

        imports = set()
        for py in list(Path("src").rglob("*.py")) + list(Path("tests").rglob("*.py")) + list(Path("scripts").rglob("*.py")):
            try:
                tree = ast.parse(py.read_text())
            except SyntaxError:
                continue
            for node in ast.walk(tree):
                if isinstance(node, ast.Import):
                    for n in node.names:
                        imports.add(n.name.split(".")[0])
                elif isinstance(node, ast.ImportFrom) and node.level == 0 and node.module:
                    imports.add(node.module.split(".")[0])

        # Use Python's own authoritative stdlib list rather than a
        # hand-maintained set: the previous literal set silently omitted
        # modules (e.g. `csv`), which made this test demand a conda package
        # for a stdlib module the moment any file imported one.
        stdlib_ish = set(sys.stdlib_module_names) | {"src", "scripts"}
        third_party = imports - stdlib_ish
        missing = sorted(third_party - declared_imports - {"PIL"})  # PIL: transitive via matplotlib, not directly declared upstream either
        assert missing == [], f"third-party imports missing from environment.yml: {missing}"

    def test_pandas_is_pinned_below_major_version_3(self):
        text = Path("environment.yml").read_text()
        assert "pandas>=2,<3" in text or "pandas <3" in text or "pandas<3" in text


class TestFiguresAndTablesExist:
    @pytest.mark.parametrize("stem", ["A_crispr_landscape", "B_sensitising_hit_evidence_matrix", "C_hany_vs_depmap", "D_selection_rule_stability"])
    def test_figure_exists(self, stem):
        assert (FIGURES / f"{stem}.png").exists()

    @pytest.mark.parametrize("n", range(1, 10))
    def test_table_exists(self, n):
        matches = list(TABLES.glob(f"{n:02d}_*.tsv"))
        assert len(matches) == 1, f"expected exactly one table starting with {n:02d}_, found {matches}"
