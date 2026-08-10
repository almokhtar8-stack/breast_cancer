# CLAUDE.md

Project-specific rules for working in this repository.

## Commits

No co-author trailers, no tool attribution in commit messages, code comments,
docstrings, or documentation.

## Data hygiene

- Every module logs rows in, rows out, and rows lost at every filter and
  join. Never drop rows silently.
- Every module is deterministic, with no network calls at runtime.
  Downloads live in `scripts/`.
- Every module has type hints and a docstring naming its data source and
  version.
- Every module has a pytest exercising its logic, not merely that it runs.

## Hard rules

- Nothing derived from the CRISPR screen may enter the feature table. The
  screen supplies labels only; features come from transcriptomics and
  public annotation.
- RCOR1 and KDM1A blinding was retired 2026-08-10 (PREANALYSIS.md amendment
  2026-08-10): the predictive model became optional future work rather than
  the primary poster deliverable, so the two genes may now be inspected and
  reported like any other CRISPR hit.
- Thresholds are declared in `PREANALYSIS.md` before results are seen.
- Paths come from `config/config.yaml`, never hardcoded.
- Figures are generated from output files, never with hand-typed numbers.
