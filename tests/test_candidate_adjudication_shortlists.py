import ast
import inspect
from pathlib import Path

import pandas as pd

import src.candidate_adjudication_axes as axes_module
import src.candidate_adjudication_decision_table as decision_module
import src.candidate_adjudication_leaders as leaders_module
import src.candidate_adjudication_near_misses as near_misses_module
import src.candidate_adjudication_shortlists as shortlists_module
from src.candidate_adjudication_shortlists import build_list_a_multimodal_therapeutic

REPO_ROOT = Path(__file__).parent.parent

FORBIDDEN_TERMS = ["druggab", "inhibitor", "pathway_enrich", "tractab", "toxicity", "essentiality_score", "literature_score"]

ADJUDICATION_LOGIC_MODULES = [axes_module, decision_module, leaders_module, near_misses_module, shortlists_module]

KNOWN_GENE_SYMBOLS = {
    "USP34", "VEZF1", "CTDNEP1", "EIF4ENIF1", "HMGB1", "KDM1A", "PET117", "TADA2B", "ICK", "SUPT4H1", "TLK2", "TSR3", "USP17L29", "PAICS",
    "CUX1", "DPP9", "LZTR1", "SOX2", "TFAP2C",  # the other five MULTIMODAL_STRONG genes -- flagged as missing from this set by the Phase 34 Codex review
}


def _executable_code_only(module) -> str:
    """Source with module/class/function docstrings stripped, so a
    docstring's *disclaimer* ("no druggability data used") doesn't
    falsely trip a check that only cares about actual identifiers/
    columns/variables in the executable code."""
    tree = ast.parse(inspect.getsource(module))
    for node in ast.walk(tree):
        if isinstance(node, (ast.Module, ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef)):
            if (node.body and isinstance(node.body[0], ast.Expr) and isinstance(node.body[0].value, ast.Constant) and isinstance(node.body[0].value.value, str)):
                node.body[0].value.value = ""
    return ast.unparse(tree)


class TestNoForbiddenVariablesInShortlistLogic:
    def test_no_druggability_or_pathway_terms_in_executable_code(self):
        for module in ADJUDICATION_LOGIC_MODULES:
            code_text = _executable_code_only(module).lower()
            for term in FORBIDDEN_TERMS:
                assert term not in code_text, f"{module.__name__} contains forbidden term '{term}' in executable code (not just a docstring)"


class TestNoHardcodedFavoredGene:
    def test_no_gene_symbol_literal_in_ranking_logic_functions(self):
        """The functions that decide *which* genes make a shortlist must
        not reference any specific gene symbol -- only the CLI/config
        wiring (candidates: [...]) is allowed to name genes."""
        ranking_functions = [
            shortlists_module.build_list_a_multimodal_therapeutic,
            shortlists_module.build_list_b_resistance_biomarker,
            shortlists_module.build_list_c_functional_sensitisation,
            shortlists_module.build_list_d_human_tumor,
            axes_module.assign_archetype,
            axes_module.classify_axis_a_functional,
            axes_module.classify_axis_b_resistance,
            axes_module.classify_axis_c_human,
        ]
        for fn in ranking_functions:
            src_text = inspect.getsource(fn)
            for gene in KNOWN_GENE_SYMBOLS:
                assert gene not in src_text, f"{fn.__name__} references gene symbol {gene}"


class TestListADeterminism:
    def _synthetic_inputs(self):
        pool = pd.DataFrame(
            {
                "gene": ["G1", "G2", "G3"],
                "crispr_direction": ["sensitising_KO", "sensitising_KO", "tolerance_associated_KO"],
                "resistance_fdr05_count": [1, 0, 2],
                "resistance_direction_consensus": ["majority_up", "insufficient", "all_down"],
                "gse240112_tumor_fdr": [0.5, 0.5, 0.01],
                "gse245601_epi_fdr": [0.5, 0.5, 0.5],
                "gse245601_malignant_fdr": [0.5, 0.5, 0.5],
            }
        )
        axes = pd.DataFrame(
            {
                "gene": ["G1", "G2", "G3"],
                "axis_a_functional": ["STRONG", "NO_EVIDENCE", "STRONG"],
                "axis_b_resistance": ["STRONG", "NO_EVIDENCE", "VERY_STRONG"],
                "axis_c_human": ["NO_EVIDENCE", "NO_EVIDENCE", "STRONG"],
                "global_rank": [100, 200, 5],
            }
        )
        stability = pd.DataFrame({"gene": ["G1", "G2", "G3"], "stability_label": ["MODERATELY_STABLE", "DATASET_DEPENDENT", "ROBUST"]})
        return pool, axes, stability

    def test_repeated_calls_produce_identical_output(self):
        pool, axes, stability = self._synthetic_inputs()
        out1 = build_list_a_multimodal_therapeutic(pool, axes, stability)
        out2 = build_list_a_multimodal_therapeutic(pool, axes, stability)
        pd.testing.assert_frame_equal(out1, out2)

    def test_gene_with_no_crispr_evidence_excluded_despite_sensitising_label(self):
        """G2 has crispr_direction='sensitising_KO' (an unconditional sign
        label) but axis_a=NO_EVIDENCE (no real CRISPR signal) and no
        resistance/human support -- must not enter the therapeutic
        shortlist on the strength of a meaningless sign alone."""
        pool, axes, stability = self._synthetic_inputs()
        out = build_list_a_multimodal_therapeutic(pool, axes, stability)
        assert "G2" not in set(out["gene"])

    def test_tolerance_direction_gene_never_enters_list_a(self):
        pool, axes, stability = self._synthetic_inputs()
        out = build_list_a_multimodal_therapeutic(pool, axes, stability)
        assert "G3" not in set(out["gene"])

    def test_contradiction_and_stability_break_ties_after_axes(self):
        """Two genes tied on all three axis bands: the one with a fully
        concordant resistance direction (all_up) and ROBUST stability must
        rank above one with a discordant (mixed) direction and
        DATASET_DEPENDENT stability."""
        pool = pd.DataFrame(
            {
                "gene": ["TIED_GOOD", "TIED_BAD"],
                "crispr_direction": ["sensitising_KO", "sensitising_KO"],
                "resistance_fdr05_count": [1, 1],
                "resistance_direction_consensus": ["all_up", "mixed"],
                "gse240112_tumor_fdr": [0.5, 0.5],
                "gse245601_epi_fdr": [0.5, 0.5],
                "gse245601_malignant_fdr": [0.5, 0.5],
            }
        )
        axes = pd.DataFrame(
            {
                "gene": ["TIED_GOOD", "TIED_BAD"],
                "axis_a_functional": ["STRONG", "STRONG"],
                "axis_b_resistance": ["STRONG", "STRONG"],
                "axis_c_human": ["NO_EVIDENCE", "NO_EVIDENCE"],
                "global_rank": [500, 10],  # deliberately opposite of the expected order, to prove global_rank isn't the tiebreak
            }
        )
        stability = pd.DataFrame({"gene": ["TIED_GOOD", "TIED_BAD"], "stability_label": ["ROBUST", "DATASET_DEPENDENT"]})
        out = build_list_a_multimodal_therapeutic(pool, axes, stability)
        assert out["gene"].tolist() == ["TIED_GOOD", "TIED_BAD"]

    def test_stability_alone_breaks_a_tie_when_contradiction_also_ties(self):
        """Isolates the stability criterion: both genes have IDENTICAL
        resistance_direction_consensus (so the contradiction band ties
        too), and only leave-one-out stability differs -- proves stability
        is consulted as its own, independent tiebreak key, not merely
        riding along with the contradiction comparison."""
        pool = pd.DataFrame(
            {
                "gene": ["TIED_ROBUST", "TIED_DEPENDENT"],
                "crispr_direction": ["sensitising_KO", "sensitising_KO"],
                "resistance_fdr05_count": [1, 1],
                "resistance_direction_consensus": ["majority_up", "majority_up"],
                "gse240112_tumor_fdr": [0.5, 0.5],
                "gse245601_epi_fdr": [0.5, 0.5],
                "gse245601_malignant_fdr": [0.5, 0.5],
            }
        )
        axes = pd.DataFrame(
            {
                "gene": ["TIED_ROBUST", "TIED_DEPENDENT"],
                "axis_a_functional": ["STRONG", "STRONG"],
                "axis_b_resistance": ["STRONG", "STRONG"],
                "axis_c_human": ["NO_EVIDENCE", "NO_EVIDENCE"],
                "global_rank": [500, 10],
            }
        )
        stability = pd.DataFrame({"gene": ["TIED_ROBUST", "TIED_DEPENDENT"], "stability_label": ["ROBUST", "DATASET_DEPENDENT"]})
        out = build_list_a_multimodal_therapeutic(pool, axes, stability)
        assert out["gene"].tolist() == ["TIED_ROBUST", "TIED_DEPENDENT"]
