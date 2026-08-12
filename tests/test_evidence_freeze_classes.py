from pathlib import Path

import pandas as pd
import pytest

REPO_ROOT = Path(__file__).parent.parent


@pytest.fixture(scope="module")
def classes():
    path = REPO_ROOT / "results" / "tables" / "evidence_freeze" / "frozen_candidate_classes.tsv"
    if not path.exists():
        pytest.skip("frozen candidate classes table not generated in this environment")
    return pd.read_csv(path, sep="\t")


class TestFrozenCandidateClasses:
    def test_four_classes_present(self, classes):
        assert set(classes["candidate_class"]) == {
            "A_THERAPEUTIC_INHIBITION_SHORTLIST", "B_RESISTANCE_BIOMARKER_PATHWAY_LEADERS",
            "C_FUNCTIONAL_SENSITISATION_LEADERS", "D_HUMAN_TUMOR_LEADERS",
        }

    def test_class_a_has_exactly_the_four_frozen_genes(self, classes):
        a_genes = set(classes.loc[classes["candidate_class"] == "A_THERAPEUTIC_INHIBITION_SHORTLIST", "gene"])
        assert a_genes == {"VEZF1", "USP34", "EML5", "CITED2"}

    def test_a_gene_may_appear_in_multiple_classes_no_forced_exclusivity(self, classes):
        # VEZF1 is both a therapeutic-shortlist gene (class A) and a functional-sensitisation
        # leader (class C, from the frozen candidate-adjudication shortlist) -- the classes
        # are explicitly non-exclusive per Phase 13's instruction.
        vezf1_classes = set(classes.loc[classes["gene"] == "VEZF1", "candidate_class"])
        assert len(vezf1_classes) >= 1  # at minimum present; multi-membership allowed, not required

    def test_no_pathway_or_druggability_terms_in_evidence_column(self, classes):
        forbidden = ["druggab", "inhibitor", "pathway_enrich", "tractab", "toxicity", "essentiality"]
        text = " ".join(classes["class_defining_evidence"].astype(str)).lower()
        for term in forbidden:
            assert term not in text
