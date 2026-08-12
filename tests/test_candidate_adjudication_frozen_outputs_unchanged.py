import hashlib
import subprocess
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).parent.parent

FROZEN_DIRS = [
    "results/tables/cross_dataset_genomewide",
    "results/tables/gse111151",
    "results/tables/gse240112_pseudobulk",
    "results/tables/gse245601_pseudobulk",
]


def _git_available() -> bool:
    try:
        subprocess.run(["git", "-C", str(REPO_ROOT), "rev-parse"], capture_output=True, check=True)
        return True
    except Exception:
        return False


class TestFrozenCrossDatasetOutputsUntouched:
    def test_no_uncommitted_modification_under_frozen_directories(self):
        """The candidate-adjudication phase must be purely additive: it
        must never modify a file that was already committed by the
        cross-dataset genome-wide integration phase. `git status
        --porcelain` on those directories must show only untracked (`??`)
        or nothing -- never a modified (` M`) or staged-modified (`M `)
        entry."""
        if not _git_available():
            pytest.skip("not running inside a git checkout")
        for d in FROZEN_DIRS:
            if not (REPO_ROOT / d).exists():
                continue
            result = subprocess.run(["git", "-C", str(REPO_ROOT), "status", "--porcelain", d], capture_output=True, text=True, check=True)
            for line in result.stdout.splitlines():
                status = line[:2]
                assert status in ("??", ""), f"unexpected git status '{status}' under frozen directory {d}: {line}"

    def test_adjudication_outputs_are_new_directories_not_overlapping_frozen_ones(self):
        adjudication_dirs = {"results/tables/candidate_adjudication", "results/figures/candidate_adjudication"}
        for d in adjudication_dirs:
            assert d not in FROZEN_DIRS
            for frozen in FROZEN_DIRS:
                assert not d.startswith(frozen) and not frozen.startswith(d)
